// makeMCweightsUnfolding.cpp
//
// Computes pp->pPb reweighting maps so that pp MC can be used as the response
// matrix when unfolding a pPb (or Pbp) zT = pT^{D0}/pT^{jet} distribution.
//
// Three maps are produced:
//   hW_det   : w(jet_pt_det, jet_eta_det) = pPb_Reco / pp_Reco  [Reco-driven, detector-level]
//   hW_mc2d  : w(jet_pt_mc,  jet_eta_mc)  = pPb_MC   / pp_MC   [MC truth-level, 2D]
//   hW_mc3d  : w(jet_pt_mc,  jet_eta_mc, d0_z_mc)              [MC truth-level, 3D, if branch exists]
//
// All TTrees are expected to be named "Response" and to contain (at minimum):
//   jet_pt_det, jet_eta_det, jet_pt_mc, jet_eta_mc
// and optionally: d0_z_det, d0_z_mc
//
// Trees/branches (float):
//   jet_pt_det, jet_eta_det, jet_pt_mc, jet_eta_mc, d0_z_det, d0_z_mc, d0_eta_det, d0_eta_mc, and more
//
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH1D.h>
#include <vector>
#include <TMath.h>
#include <iostream>
#include <string>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>
#include <ctime>
#include <algorithm>
#include <cmath>

// ---------------------------------------------------------------------------
// Safe bin-content ratio: returns fallback when denominator is zero or tiny
// ---------------------------------------------------------------------------
static double SafeDiv(double num, double den, double fallback = 1.0)
{
    if (den <= 0.0 || !std::isfinite(num) || !std::isfinite(den)) return fallback;
    double r = num / den;
    return std::isfinite(r) ? r : fallback;
}

// ---------------------------------------------------------------------------
// Clamp a weight to [0, wMax]; bad values -> 0
// ---------------------------------------------------------------------------
static double Clamp(double w, double wMax = 5.0)
{
    if (!std::isfinite(w) || w < 0.0) return 0.0;
    return w > wMax ? wMax : w;
}

// ---------------------------------------------------------------------------
// Build a direct-ratio weight map from two (already-normalised) histograms
// Bins with fewer than minCount raw entries in either sample fall back to 1.
// ---------------------------------------------------------------------------
static TH2D* MakeWeightMap2D(const TH2D* hPP_norm, const TH2D* hPpB_norm,
                              const TH2D* hPP_raw,  const TH2D* hPpB_raw,
                              const char* name, const char* title,
                              double wMax = 5.0, double minCount = 5.0)
{
    TH2D* hW = (TH2D*)hPpB_norm->Clone(name);
    hW->SetTitle(title);
    hW->Reset();
    for (int ix = 1; ix <= hW->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= hW->GetNbinsY(); ++iy) {
            double ppRaw  = hPP_raw  ? hPP_raw ->GetBinContent(ix, iy) : 1e9;
            double pPbRaw = hPpB_raw ? hPpB_raw->GetBinContent(ix, iy) : 1e9;
            double w;
            if (ppRaw < minCount || pPbRaw < minCount) {
                w = 1.0; // not enough statistics — leave unweighted
            } else {
                w = Clamp(SafeDiv(hPpB_norm->GetBinContent(ix, iy),
                                  hPP_norm ->GetBinContent(ix, iy), 1.0), wMax);
            }
            hW->SetBinContent(ix, iy, w);
        }
    }
    return hW;
}

static TH3D* MakeWeightMap3D(const TH3D* hPP_norm, const TH3D* hPpB_norm,
                              const TH3D* hPP_raw,  const TH3D* hPpB_raw,
                              const char* name, const char* title,
                              double wMax = 5.0, double minCount = 3.0)
{
    TH3D* hW = (TH3D*)hPpB_norm->Clone(name);
    hW->SetTitle(title);
    hW->Reset();
    for (int ix = 1; ix <= hW->GetNbinsX(); ++ix)
        for (int iy = 1; iy <= hW->GetNbinsY(); ++iy)
            for (int iz = 1; iz <= hW->GetNbinsZ(); ++iz) {
                double ppRaw  = hPP_raw  ? hPP_raw ->GetBinContent(ix, iy, iz) : 1e9;
                double pPbRaw = hPpB_raw ? hPpB_raw->GetBinContent(ix, iy, iz) : 1e9;
                double w;
                if (ppRaw < minCount || pPbRaw < minCount) {
                    // fall back to the 2D projection ratio along pt×eta
                    w = 1.0;
                } else {
                    w = Clamp(SafeDiv(hPpB_norm->GetBinContent(ix, iy, iz),
                                      hPP_norm ->GetBinContent(ix, iy, iz), 1.0), wMax);
                }
                hW->SetBinContent(ix, iy, iz, w);
            }
    return hW;
}

// ---------------------------------------------------------------------------
// Compute chi2/ndf between two normalised 1D histograms (symmetric errors)
// ---------------------------------------------------------------------------
static void PrintChi2(const TH1D* hA, const TH1D* hB, const char* label)
{
    double chi2 = 0.0; int ndf = 0;
    int nb = std::min(hA->GetNbinsX(), hB->GetNbinsX());
    for (int ib = 1; ib <= nb; ++ib) {
        double va = hA->GetBinContent(ib), vb = hB->GetBinContent(ib);
        double ea = hA->GetBinError(ib),   eb = hB->GetBinError(ib);
        double e2 = ea*ea + eb*eb;
        if (e2 <= 0) continue;
        chi2 += (va - vb)*(va - vb) / e2;
        ++ndf;
    }
    std::cout << "[chi2] " << label
              << "  chi2=" << chi2 << "  ndf=" << ndf
              << "  chi2/ndf=" << (ndf > 0 ? chi2/ndf : -1.0) << "\n";
}

// ---------------------------------------------------------------------------
// Draw an overlay panel (top) + ratio panel (bottom) and write/save
// ---------------------------------------------------------------------------
static void OverlayRatioCanvas(
    const std::vector<TH1D*>& hists,
    const std::vector<const char*>& labels,
    const std::vector<int>& colors,
    TH1D* hNum, TH1D* hDen,           // for the ratio panel (hNum/hDen)
    const char* ratioYtitle,
    const char* canvasName,
    const std::string& saveDir,
    TFile* fout)
{
    TCanvas* c = new TCanvas(canvasName, canvasName, 900, 800);
    TPad* pTop = new TPad("pTop","",0,0.30,1,1.00);
    TPad* pBot = new TPad("pBot","",0,0.00,1,0.30);
    pTop->SetBottomMargin(0.02); pBot->SetTopMargin(0.02); pBot->SetBottomMargin(0.35);
    pTop->Draw(); pBot->Draw();
    pTop->cd();

    double ymax = 0.0;
    for (auto* h : hists)
        for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
            ymax = std::max(ymax, h->GetBinContent(ib) + h->GetBinError(ib));

    bool first = true;
    for (size_t i = 0; i < hists.size(); ++i) {
        hists[i]->SetLineColor(colors[i]); hists[i]->SetMarkerColor(colors[i]);
        hists[i]->SetMarkerStyle(20 + (int)i);
        if (ymax > 0) hists[i]->SetMaximum(1.3 * ymax);
        hists[i]->SetMinimum(0);
        hists[i]->Draw(first ? "E1" : "E1 SAME");
        first = false;
    }
    TLegend* leg = new TLegend(0.55, 0.65, 0.92, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    for (size_t i = 0; i < hists.size(); ++i)
        leg->AddEntry(hists[i], labels[i], "lep");
    leg->Draw();

    pBot->cd();
    TH1D* hRat = (TH1D*)hNum->Clone(Form("%s_ratio", canvasName));
    hRat->Divide(hDen);
    hRat->SetLineColor(kBlack); hRat->SetMarkerStyle(20);
    hRat->SetMinimum(0.5); hRat->SetMaximum(1.5);
    hRat->GetYaxis()->SetTitle(ratioYtitle);
    hRat->GetYaxis()->SetTitleSize(0.12); hRat->GetYaxis()->SetTitleOffset(0.4);
    hRat->GetYaxis()->SetLabelSize(0.10);
    hRat->GetXaxis()->SetLabelSize(0.10); hRat->GetXaxis()->SetTitleSize(0.11);
    hRat->Draw("E1");
    TLine ln(hRat->GetXaxis()->GetXmin(), 1.0, hRat->GetXaxis()->GetXmax(), 1.0);
    ln.SetLineStyle(2); ln.SetLineColor(kGray+2); ln.Draw();

    if (fout) { fout->cd(); hRat->Write(); c->Write(); }
    c->SaveAs((saveDir + "/" + canvasName + ".png").c_str());
    c->SaveAs((saveDir + "/" + canvasName + ".pdf").c_str());
    delete hRat; delete leg; delete c;
}

// ===========================================================================
void makeMCweightsUnfolding(const char* ppRecoFile, const char* ppMCFile,
                             const char* pPbRecoFile, const char* pPbMCFile,
                             const char* outFile = "weights.root")
{
    // -----------------------------------------------------------------------
    // Phase space / binning
    // -----------------------------------------------------------------------
    const int    nPtBins  = 20;
    const double ptMin    =  5.0, ptMax  = 55.0;
    const int    nEtaBins = 12;
    const double etaMin   =  2.5, etaMax =  4.0;
    const int    nD0Bins  = 20;   // for d0_z_mc axis in 3D map
    const double d0Min    =  0.0, d0Max  =  1.0;  // zT = pT^D0/pT^jet in [0,1]
    const double wMax     =  5.0; // maximum allowed weight
    const double minCount2D = 5.0; // min raw entries per bin for direct ratio (2D)
    const double minCount3D = 3.0; // min raw entries per bin for direct ratio (3D)

    // -----------------------------------------------------------------------
    // Open files and get trees
    // -----------------------------------------------------------------------
    TFile* fPPReco  = TFile::Open(ppRecoFile,  "READ");
    TFile* fPPmc    = TFile::Open(ppMCFile,    "READ");
    TFile* fPpBReco = TFile::Open(pPbRecoFile, "READ");
    TFile* fPpBmc   = TFile::Open(pPbMCFile,   "READ");

    auto CheckFile = [](TFile* f, const char* name) -> bool {
        if (!f || f->IsZombie()) { std::cerr << "Cannot open " << name << "\n"; return false; }
        return true;
    };
    if (!CheckFile(fPPReco, ppRecoFile) || !CheckFile(fPPmc,    ppMCFile) ||
        !CheckFile(fPpBReco,pPbRecoFile)|| !CheckFile(fPpBmc,   pPbMCFile)) return;

    // Try several common tree names
    auto GetTree = [](TFile* f, const char* label) -> TTree* {
        for (const char* n : {"Response", "T", "tree", "Tree"}) {
            TTree* t = dynamic_cast<TTree*>(f->Get(n));
            if (t) { std::cout << "  " << label << ": tree='" << n << "'\n"; return t; }
        }
        std::cerr << "Cannot find TTree in " << label << "\n"; return nullptr;
    };
    TTree* tPPReco  = GetTree(fPPReco,  "pp Reco");
    TTree* tPPmc    = GetTree(fPPmc,    "pp MC");
    TTree* tPpBReco = GetTree(fPpBReco, "pPb Reco");
    TTree* tPpBmc   = GetTree(fPpBmc,   "pPb MC");
    if (!tPPReco || !tPPmc || !tPpBReco || !tPpBmc) return;

    // -----------------------------------------------------------------------
    // Create output directory (dated) and output ROOT file
    // -----------------------------------------------------------------------
    {
        std::string s(outFile);
        size_t sl = s.rfind('/'); if (sl != std::string::npos) s = s.substr(sl+1);
        size_t dt = s.rfind('.'); if (dt != std::string::npos) s = s.substr(0, dt);
        char buf[32]{};
        time_t now = time(nullptr); struct tm* lt = localtime(&now);
        if (lt) strftime(buf, sizeof(buf), "%Y-%m-%d", lt); else snprintf(buf,sizeof(buf),"date");
        // store as globals for use below
    }
    // Build outDir from outFile stem + date
    std::string outStem(outFile);
    { size_t p = outStem.rfind('/'); if (p!=std::string::npos) outStem=outStem.substr(p+1);
      size_t d = outStem.rfind('.'); if (d!=std::string::npos) outStem=outStem.substr(0,d); }
    char dateBuf[32]{};
    { time_t now=time(nullptr); struct tm* lt=localtime(&now);
      if(lt) strftime(dateBuf,sizeof(dateBuf),"%Y-%m-%d",lt);
      else   snprintf(dateBuf,sizeof(dateBuf),"date"); }
    std::string outDir = outStem + "_plots_" + dateBuf;
    gSystem->mkdir(outDir.c_str(), /*recursive=*/true);
    std::string rootPath = outDir + "/" + outStem + ".root";

    TFile* fout = TFile::Open(rootPath.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) { std::cerr << "Cannot create " << rootPath << "\n"; return; }
    std::cout << "Output: " << rootPath << "\n";

    // -----------------------------------------------------------------------
    // STEP 1 — Fill 2D detector-level distributions from pp Reco and pPb Reco
    //          -> Reco-driven weight map hW_det
    // -----------------------------------------------------------------------
    TH2D* hPPReco_det_raw  = new TH2D("hPPReco_det_raw",
        "pp Reco det;jet p_{T}^{det} (GeV/c);jet #eta^{det}",
        nPtBins,ptMin,ptMax, nEtaBins,etaMin,etaMax);
    TH2D* hPpBReco_det_raw = new TH2D("hPpBReco_det_raw",
        "pPb Reco det;jet p_{T}^{det} (GeV/c);jet #eta^{det}",
        nPtBins,ptMin,ptMax, nEtaBins,etaMin,etaMax);
    hPPReco_det_raw->Sumw2(); hPpBReco_det_raw->Sumw2();

    {
        TTreeReader r(tPPReco);
        TTreeReaderValue<float> pt(r,"jet_pt_det"), eta(r,"jet_eta_det");
        while (r.Next()) hPPReco_det_raw->Fill(*pt, *eta);
    }
    {
        TTreeReader r(tPpBReco);
        TTreeReaderValue<float> pt(r,"jet_pt_det"), eta(r,"jet_eta_det");
        while (r.Next()) hPpBReco_det_raw->Fill(*pt, *eta);
    }

    // Normalised copies for ratio
    TH2D* hPPReco_det  = (TH2D*)hPPReco_det_raw ->Clone("hPPReco_det");
    TH2D* hPpBReco_det = (TH2D*)hPpBReco_det_raw->Clone("hPpBReco_det");
    if (hPPReco_det ->Integral()>0) hPPReco_det ->Scale(1.0/hPPReco_det ->Integral());
    if (hPpBReco_det->Integral()>0) hPpBReco_det->Scale(1.0/hPpBReco_det->Integral());

    TH2D* hW_det = MakeWeightMap2D(hPPReco_det, hPpBReco_det,
                                    hPPReco_det_raw, hPpBReco_det_raw,
                                    "hW_det",
                                    "Reco-driven det-level weight w(p_{T}^{det},#eta^{det}) = pPb/pp;"
                                    "jet p_{T}^{det} (GeV/c);jet #eta^{det}",
                                    wMax, minCount2D);

    std::cout << "\n[hW_det] Reco-driven detector-level weight map built.\n";

    // -----------------------------------------------------------------------
    // STEP 2 — Fill 2D truth-level distributions from pp MC and pPb MC
    //          -> MC-based weight map hW_mc2d
    // -----------------------------------------------------------------------
    TH2D* hPPmc_mc_raw  = new TH2D("hPPmc_mc_raw",
        "pp MC truth;jet p_{T}^{mc} (GeV/c);jet #eta^{mc}",
        nPtBins,ptMin,ptMax, nEtaBins,etaMin,etaMax);
    TH2D* hPpBmc_mc_raw = new TH2D("hPpBmc_mc_raw",
        "pPb MC truth;jet p_{T}^{mc} (GeV/c);jet #eta^{mc}",
        nPtBins,ptMin,ptMax, nEtaBins,etaMin,etaMax);
    hPPmc_mc_raw->Sumw2(); hPpBmc_mc_raw->Sumw2();

    {
        TTreeReader r(tPPmc);
        TTreeReaderValue<float> pt(r,"jet_pt_mc"), eta(r,"jet_eta_mc");
        while (r.Next()) hPPmc_mc_raw->Fill(*pt, *eta);
    }
    {
        TTreeReader r(tPpBmc);
        TTreeReaderValue<float> pt(r,"jet_pt_mc"), eta(r,"jet_eta_mc");
        while (r.Next()) hPpBmc_mc_raw->Fill(*pt, *eta);
    }

    TH2D* hPPmc_mc  = (TH2D*)hPPmc_mc_raw ->Clone("hPPmc_mc");
    TH2D* hPpBmc_mc = (TH2D*)hPpBmc_mc_raw->Clone("hPpBmc_mc");
    if (hPPmc_mc ->Integral()>0) hPPmc_mc ->Scale(1.0/hPPmc_mc ->Integral());
    if (hPpBmc_mc->Integral()>0) hPpBmc_mc->Scale(1.0/hPpBmc_mc->Integral());

    TH2D* hW_mc2d = MakeWeightMap2D(hPPmc_mc, hPpBmc_mc,
                                     hPPmc_mc_raw, hPpBmc_mc_raw,
                                     "hW_mc2d",
                                     "MC truth-level weight w(p_{T}^{mc},#eta^{mc}) = pPb_MC/pp_MC;"
                                     "jet p_{T}^{mc} (GeV/c);jet #eta^{mc}",
                                     wMax, minCount2D);

    std::cout << "[hW_mc2d] MC truth-level 2D weight map built.\n";

    // -----------------------------------------------------------------------
    // STEP 3 — Optionally build a 3D weight map including d0_z_mc
    //          (zT observable at truth level)
    // -----------------------------------------------------------------------
    TH3D* hW_mc3d = nullptr;
    const bool have3D = tPPmc->GetBranch("d0_z_mc") && tPpBmc->GetBranch("d0_z_mc");

    if (have3D) {
        TH3D* hPPmc_3d_raw  = new TH3D("hPPmc_3d_raw",
            "pp MC (pt,#eta,z_{T}^{mc});jet p_{T}^{mc} (GeV/c);jet #eta^{mc};z_{T}^{mc}",
            nPtBins,ptMin,ptMax, nEtaBins,etaMin,etaMax, nD0Bins,d0Min,d0Max);
        TH3D* hPpBmc_3d_raw = new TH3D("hPpBmc_3d_raw",
            "pPb MC (pt,#eta,z_{T}^{mc});jet p_{T}^{mc} (GeV/c);jet #eta^{mc};z_{T}^{mc}",
            nPtBins,ptMin,ptMax, nEtaBins,etaMin,etaMax, nD0Bins,d0Min,d0Max);
        hPPmc_3d_raw->Sumw2(); hPpBmc_3d_raw->Sumw2();

        {
            TTreeReader r(tPPmc);
            TTreeReaderValue<float> pt(r,"jet_pt_mc"), eta(r,"jet_eta_mc"), d0z(r,"d0_z_mc");
            while (r.Next()) hPPmc_3d_raw->Fill(*pt, *eta, *d0z);
        }
        {
            TTreeReader r(tPpBmc);
            TTreeReaderValue<float> pt(r,"jet_pt_mc"), eta(r,"jet_eta_mc"), d0z(r,"d0_z_mc");
            while (r.Next()) hPpBmc_3d_raw->Fill(*pt, *eta, *d0z);
        }

        TH3D* hPPmc_3d  = (TH3D*)hPPmc_3d_raw ->Clone("hPPmc_3d");
        TH3D* hPpBmc_3d = (TH3D*)hPpBmc_3d_raw->Clone("hPpBmc_3d");
        if (hPPmc_3d ->Integral()>0) hPPmc_3d ->Scale(1.0/hPPmc_3d ->Integral());
        if (hPpBmc_3d->Integral()>0) hPpBmc_3d->Scale(1.0/hPpBmc_3d->Integral());

        hW_mc3d = MakeWeightMap3D(hPPmc_3d, hPpBmc_3d,
                                   hPPmc_3d_raw, hPpBmc_3d_raw,
                                   "hW_mc3d",
                                   "MC 3D weight w(p_{T}^{mc},#eta^{mc},z_{T}^{mc}) = pPb_MC/pp_MC;"
                                   "jet p_{T}^{mc} (GeV/c);jet #eta^{mc};z_{T}^{mc}",
                                   wMax, minCount3D);

        // For low-stat 3D bins, fall back to the 2D weight
        for (int ix = 1; ix <= hW_mc3d->GetNbinsX(); ++ix)
            for (int iy = 1; iy <= hW_mc3d->GetNbinsY(); ++iy)
                for (int iz = 1; iz <= hW_mc3d->GetNbinsZ(); ++iz) {
                    double raw_pp  = hPPmc_3d_raw ->GetBinContent(ix,iy,iz);
                    double raw_pPb = hPpBmc_3d_raw->GetBinContent(ix,iy,iz);
                    if (raw_pp < minCount3D || raw_pPb < minCount3D) {
                        // use 2D fallback
                        hW_mc3d->SetBinContent(ix, iy, iz, hW_mc2d->GetBinContent(ix, iy));
                    }
                }

        std::cout << "[hW_mc3d] MC truth-level 3D weight map built (with 2D fallback for sparse bins).\n";

        fout->cd();
        hPPmc_3d_raw->Write(); hPpBmc_3d_raw->Write();
        hPPmc_3d    ->Write(); hPpBmc_3d    ->Write();
        hW_mc3d     ->Write();
        delete hPPmc_3d; delete hPpBmc_3d;
        delete hPPmc_3d_raw; delete hPpBmc_3d_raw;
    } else {
        std::cout << "[hW_mc3d] branch d0_z_mc not found — skipping 3D weight map.\n";
    }

    // -----------------------------------------------------------------------
    // STEP 4 — Closure tests: apply each weight map to pp MC and compare
    //          against pPb MC (mc-level) or pPb Reco (det-level)
    //
    // Closure A: w_det applied to pp MC det  vs  pPb Reco det  (jet pT)
    // Closure B: w_mc2d applied to pp MC mc  vs  pPb MC mc     (jet pT)
    // -----------------------------------------------------------------------
    const int nPtPlot = 40; const double ptPlotMin = 0, ptPlotMax = 60;

    // -- Closure A --
    TH1D* hA_ppmc_unw = new TH1D("hA_ppmc_detPt_unw",
        "pp MC det (unweighted);jet p_{T}^{det} (GeV/c);Norm. entries", nPtPlot,ptPlotMin,ptPlotMax);
    TH1D* hA_ppmc_w   = new TH1D("hA_ppmc_detPt_w",
        "pp MC det (Reco weighted);jet p_{T}^{det} (GeV/c);Norm. entries", nPtPlot,ptPlotMin,ptPlotMax);
    TH1D* hA_pPbReco  = new TH1D("hA_pPbReco_detPt",
        "pPb Reco det;jet p_{T}^{det} (GeV/c);Norm. entries", nPtPlot,ptPlotMin,ptPlotMax);
    hA_ppmc_unw->Sumw2(); hA_ppmc_w->Sumw2(); hA_pPbReco->Sumw2();

    {
        TTreeReader r(tPPmc);
        TTreeReaderValue<float> pt(r,"jet_pt_det"), eta(r,"jet_eta_det");
        while (r.Next()) {
            hA_ppmc_unw->Fill(*pt);
            int bx = hW_det->GetXaxis()->FindBin(*pt);
            int by = hW_det->GetYaxis()->FindBin(*eta);
            double w = (bx>=1&&bx<=hW_det->GetNbinsX()&&by>=1&&by<=hW_det->GetNbinsY())
                       ? hW_det->GetBinContent(bx,by) : 1.0;
            if (w>0) hA_ppmc_w->Fill(*pt, w);
        }
    }
    { TTreeReader r(tPpBReco); TTreeReaderValue<float> pt(r,"jet_pt_det");
      while(r.Next()) hA_pPbReco->Fill(*pt); }

    // -- Closure B --
    TH1D* hB_ppmc_unw = new TH1D("hB_ppmc_mcPt_unw",
        "pp MC truth (unweighted);jet p_{T}^{mc} (GeV/c);Norm. entries", nPtPlot,ptPlotMin,ptPlotMax);
    TH1D* hB_ppmc_w   = new TH1D("hB_ppmc_mcPt_w",
        "pp MC truth (mc weighted);jet p_{T}^{mc} (GeV/c);Norm. entries", nPtPlot,ptPlotMin,ptPlotMax);
    TH1D* hB_pPbmc    = new TH1D("hB_pPbmc_mcPt",
        "pPb MC truth;jet p_{T}^{mc} (GeV/c);Norm. entries", nPtPlot,ptPlotMin,ptPlotMax);
    hB_ppmc_unw->Sumw2(); hB_ppmc_w->Sumw2(); hB_pPbmc->Sumw2();

    {
        TTreeReader r(tPPmc);
        TTreeReaderValue<float> pt(r,"jet_pt_mc"), eta(r,"jet_eta_mc");
        // optional 3D branch
        float d0z_val = 0.0f;
        TBranch* bd0 = tPPmc->GetBranch("d0_z_mc");
        if (bd0) tPPmc->SetBranchAddress("d0_z_mc",&d0z_val);
        while (r.Next()) {
            if (bd0) tPPmc->GetEntry(r.GetCurrentEntry());
            hB_ppmc_unw->Fill(*pt);
            double w = 1.0;
            if (hW_mc3d) {
                int bx = hW_mc3d->GetXaxis()->FindBin(*pt);
                int by = hW_mc3d->GetYaxis()->FindBin(*eta);
                int bz = hW_mc3d->GetZaxis()->FindBin(d0z_val);
                if (bx>=1&&bx<=hW_mc3d->GetNbinsX()&&by>=1&&by<=hW_mc3d->GetNbinsY()&&bz>=1&&bz<=hW_mc3d->GetNbinsZ())
                    w = hW_mc3d->GetBinContent(bx,by,bz);
            } else {
                int bx = hW_mc2d->GetXaxis()->FindBin(*pt);
                int by = hW_mc2d->GetYaxis()->FindBin(*eta);
                if (bx>=1&&bx<=hW_mc2d->GetNbinsX()&&by>=1&&by<=hW_mc2d->GetNbinsY())
                    w = hW_mc2d->GetBinContent(bx,by);
            }
            if (w>0) hB_ppmc_w->Fill(*pt, w);
        }
    }
    { TTreeReader r(tPpBmc); TTreeReaderValue<float> pt(r,"jet_pt_mc");
      while(r.Next()) hB_pPbmc->Fill(*pt); }

    // Normalise all closure histograms to shape
    auto Norm1D = [](TH1D* h){ if(h->Integral()>0) h->Scale(1.0/h->Integral()); };
    Norm1D(hA_ppmc_unw); Norm1D(hA_ppmc_w); Norm1D(hA_pPbReco);
    Norm1D(hB_ppmc_unw); Norm1D(hB_ppmc_w); Norm1D(hB_pPbmc);

    PrintChi2(hA_ppmc_w, hA_pPbReco, "Closure A: pp MC (Reco-wt det) vs pPb Reco");
    PrintChi2(hB_ppmc_w, hB_pPbmc,   "Closure B: pp MC (mc-wt truth) vs pPb MC");

    // -----------------------------------------------------------------------
    // STEP 5 — zT (d0_z_det) closure in jet-pT slices
    //          pp MC (weighted) vs pPb Reco, and pp MC (unweighted), in bins of jet pT
    // -----------------------------------------------------------------------
    const bool haveZTdet = tPPmc->GetBranch("d0_z_det") && tPpBReco->GetBranch("d0_z_det");
    if (!haveZTdet)
        std::cout << "[zT slices] branch d0_z_det missing in one of the trees — skipping.\n";

    std::vector<std::pair<double,double>> ptSlices = {{5,10},{10,15},{15,20},{20,30}};
    const int nZTbins = 20; const double zTmin=0.0, zTmax=1.2;

    std::vector<TH1D*> vh_ppmc_unw, vh_ppmc_w, vh_pPbReco;
    if (haveZTdet) {
        for (auto& slice : ptSlices) {
            double lo = slice.first, hi = slice.second;
            auto mkh = [&](const char* tag, const char* ttl) {
                TH1D* h = new TH1D(Form("hzT_%s_%.0f_%.0f",tag,lo,hi),
                                   Form("%s [%.0f,%.0f];z_{T};Norm. entries",ttl,lo,hi),
                                   nZTbins,zTmin,zTmax);
                h->Sumw2(); return h;
            };
            vh_ppmc_unw.push_back(mkh("ppmc_unw","pp MC (unw)"));
            vh_ppmc_w  .push_back(mkh("ppmc_w",  "pp MC (wt)" ));
            vh_pPbReco .push_back(mkh("pPbReco",  "pPb Reco"  ));
        }

        // Fill pp MC slices
        {
            TTreeReader r(tPPmc);
            TTreeReaderValue<float> pt_det(r,"jet_pt_det"), eta_det(r,"jet_eta_det");
            TTreeReaderValue<float> zT(r,"d0_z_det");
            while (r.Next()) {
                for (size_t is = 0; is < ptSlices.size(); ++is) {
                    double lo=ptSlices[is].first, hi=ptSlices[is].second;
                    if (*pt_det>=lo && *pt_det<hi) {
                        vh_ppmc_unw[is]->Fill(*zT);
                        int bx = hW_det->GetXaxis()->FindBin(*pt_det);
                        int by = hW_det->GetYaxis()->FindBin(*eta_det);
                        double w = (bx>=1&&bx<=hW_det->GetNbinsX()&&by>=1&&by<=hW_det->GetNbinsY())
                                   ? hW_det->GetBinContent(bx,by) : 1.0;
                        if (w>0) vh_ppmc_w[is]->Fill(*zT, w);
                        break;
                    }
                }
            }
        }
        // Fill pPb Reco slices
        {
            TTreeReader r(tPpBReco);
            TTreeReaderValue<float> pt_det(r,"jet_pt_det");
            TTreeReaderValue<float> zT(r,"d0_z_det");
            while (r.Next()) {
                for (size_t is = 0; is < ptSlices.size(); ++is) {
                    double lo=ptSlices[is].first, hi=ptSlices[is].second;
                    if (*pt_det>=lo && *pt_det<hi) { vh_pPbReco[is]->Fill(*zT); break; }
                }
            }
        }

        for (size_t is = 0; is < ptSlices.size(); ++is) {
            Norm1D(vh_ppmc_unw[is]); Norm1D(vh_ppmc_w[is]); Norm1D(vh_pPbReco[is]);
            PrintChi2(vh_ppmc_w[is], vh_pPbReco[is],
                      Form("zT closure [%.0f-%.0f GeV/c]: pp MC (wt) vs pPb Reco",
                           ptSlices[is].first, ptSlices[is].second));
        }
    }

    // -----------------------------------------------------------------------
    // STEP 6 — Write histograms and diagnostic plots
    // -----------------------------------------------------------------------
    fout->cd();

    // Weight maps
    hPPReco_det_raw ->Write(); hPpBReco_det_raw ->Write();
    hPPReco_det     ->Write(); hPpBReco_det     ->Write();
    hW_det          ->Write();

    hPPmc_mc_raw    ->Write(); hPpBmc_mc_raw    ->Write();
    hPPmc_mc        ->Write(); hPpBmc_mc        ->Write();
    hW_mc2d         ->Write();

    // Closure histograms
    hA_ppmc_unw->Write(); hA_ppmc_w->Write(); hA_pPbReco->Write();
    hB_ppmc_unw->Write(); hB_ppmc_w->Write(); hB_pPbmc  ->Write();

    if (haveZTdet)
        for (size_t is = 0; is < ptSlices.size(); ++is) {
            vh_ppmc_unw[is]->Write(); vh_ppmc_w[is]->Write(); vh_pPbReco[is]->Write();
        }

    // --- Closure A canvas ---
    OverlayRatioCanvas(
        {hA_ppmc_w, hA_pPbReco, hA_ppmc_unw},
        {"pp MC (Reco wt)", "pPb Reco", "pp MC (unw)"},
        {kBlue, kRed, kGreen+2},
        hA_pPbReco, hA_ppmc_w, "pPb / pp MC (wt)",
        "closure_det_ptjet", outDir, fout);

    // --- Closure B canvas ---
    OverlayRatioCanvas(
        {hB_ppmc_w, hB_pPbmc, hB_ppmc_unw},
        {"pp MC (mc wt)", "pPb MC", "pp MC (unw)"},
        {kBlue, kRed, kGreen+2},
        hB_pPbmc, hB_ppmc_w, "pPb MC / pp MC (wt)",
        "closure_mc_ptjet", outDir, fout);

    // --- zT slice canvases ---
    if (haveZTdet) {
        for (size_t is = 0; is < ptSlices.size(); ++is) {
            OverlayRatioCanvas(
                {vh_ppmc_w[is], vh_pPbReco[is], vh_ppmc_unw[is]},
                {"pp MC (wt)", "pPb Reco", "pp MC (unw)"},
                {kBlue, kRed, kGreen+2},
                vh_pPbReco[is], vh_ppmc_w[is], "pPb / pp (wt)",
                Form("zT_slice_%.0f_%.0f", ptSlices[is].first, ptSlices[is].second),
                outDir, fout);
        }
    }

    fout->Close();
    fPPReco->Close(); fPPmc->Close(); fPpBReco->Close(); fPpBmc->Close();

    std::cout << "\nDone.\n"
              << "  Output ROOT : " << rootPath << "\n"
              << "  Plots dir   : " << outDir   << "\n\n"
              << "Weight map names:\n"
              << "  hW_det   -- 2D Reco-driven  w(jet_pt_det, jet_eta_det)\n"
              << "  hW_mc2d  -- 2D MC truth-lev w(jet_pt_mc,  jet_eta_mc)\n"
              << "  hW_mc3d  -- 3D MC truth-lev w(jet_pt_mc,  jet_eta_mc, d0_z_mc)  [if branch present]\n";
}