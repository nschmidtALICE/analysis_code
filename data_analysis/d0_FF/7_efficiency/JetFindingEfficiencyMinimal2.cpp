#include <iostream>
#include <vector>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"
#include "TLatex.h"
#include "TLegend.h"
#include <filesystem>
#include <cmath>

// Minimal jet finding efficiency calculator
// Usage: build with root-config flags, e.g.:
// g++ -std=c++17 JetFindingEfficiencyMinimal.cpp $(root-config --cflags --libs) -o JetFindingEfficiencyMinimal
// ./JetFindingEfficiencyMinimal input.root output_prefix

static double gJetRadius = 0.5;
static double gMinJetPt = 5.0;
static double gMaxJetPt = 60.0;
static double gMinJetEta = 2.5;
static double gMaxJetEta = 4.0;
static double gMinD0Pt = 1.0;
static double gMaxD0Pt = 50.0;
static double gMinD0Eta = 2.0;
static double gMaxD0Eta = 4.5;

inline double DeltaPhi(double a, double b) {
    double d = a - b;
    while (d > M_PI) d -= 2*M_PI;
    while (d <= -M_PI) d += 2*M_PI;
    return d;
}
inline double DeltaR(double eta1, double phi1, double eta2, double phi2) {
    double dphi = DeltaPhi(phi1, phi2);
    double deta = eta1 - eta2;
    return std::sqrt(deta*deta + dphi*dphi);
}

inline bool PassesJetSelectionMC(const std::vector<float>* mc_jet_pt, const std::vector<float>* mc_jet_eta, int idx) {
    if (!mc_jet_pt || !mc_jet_eta) return false;
    if (idx < 0 || idx >= (int)mc_jet_pt->size()) return false;
    double pt = mc_jet_pt->at(idx);
    double eta = mc_jet_eta->at(idx);
    if (pt < gMinJetPt || pt > gMaxJetPt) return false;
    if (eta < gMinJetEta || eta > gMaxJetEta) return false;
    return true;
}
inline bool PassesJetSelectionReco(const std::vector<float>* jet_pt, const std::vector<float>* jet_eta, int idx) {
    if (!jet_pt || !jet_eta) return false;
    if (idx < 0 || idx >= (int)jet_pt->size()) return false;
    double pt = jet_pt->at(idx);
    double eta = jet_eta->at(idx);
    if (pt < gMinJetPt || pt > gMaxJetPt) return false;
    if (eta < gMinJetEta || eta > gMaxJetEta) return false;
    return true;
}
inline bool PassesD0Selection(double d0_pt, double d0_eta) {
    if (d0_pt < gMinD0Pt || d0_pt > gMaxD0Pt) return false;
    if (d0_eta < gMinD0Eta || d0_eta > gMaxD0Eta) return false;
    return true;
}

int JetFindingEfficiencyMinimal2() {
    // std::string inputName = "/media/niviths/SSD2/lhcb_analysis_SSD/20250728_pPb_MC_output/54/54.root";
    std::string inputName = "/media/niviths/SSD2/lhcb_analysis_SSD/20250728_pPb_MC_output/20250728_pPb_MC_output.root";
    // std::string outPrefix = "jetFindingEffMinimal_MBonly";
    std::string outPrefix = "jetFindingEffMinimal_fullMC";

    // create output directory for all plots and root file (append today's date)
    char _datebuf[32];
    {
        time_t _t = time(nullptr);
        struct tm _tm = *localtime(&_t);
        strftime(_datebuf, sizeof(_datebuf), "%Y-%m-%d", &_tm);
    }
    std::string outDir = outPrefix + "_outputs_" + std::string(_datebuf);
    try { std::filesystem::create_directories(outDir); } catch (...) {
        std::cerr << "Warning: failed to create output directory '" << outDir << "' - will attempt to write files to current directory\n"; outDir = "."; }
    // secondary folder for PNG copies
    std::string pngDir = outDir + "/png";
    try { std::filesystem::create_directories(pngDir); } catch (...) {
        std::cerr << "Warning: failed to create png subdirectory '" << pngDir << "'\n"; }

    // helper to save canvas to primary png path and duplicate into pngDir
    auto savePNG = [&](TCanvas &canv, const std::string &fullPath){
        canv.SaveAs(fullPath.c_str());
        try {
            std::filesystem::path p(fullPath);
            std::filesystem::path dest = std::filesystem::path(pngDir) / p.filename();
            canv.SaveAs(dest.string().c_str());
        } catch (...) { /* ignore */ }
    };

    std::cout << "Input: " << inputName << "\n";
    std::cout << "Output prefix: " << outPrefix << "\n";

    TFile *f = TFile::Open(inputName.c_str(), "READ");
    if (!f || f->IsZombie()) { std::cerr << "Failed to open input file\n"; return 1; }
    TTree *t = (TTree*)f->Get("d0jets");
    if (!t) { std::cerr << "d0jets tree not found\n"; return 1; }

    // Branches (pointers)
    std::vector<float>* jet_pt = nullptr; std::vector<float>* jet_eta = nullptr; std::vector<float>* jet_phi = nullptr;
    std::vector<float>* d0_pt = nullptr; std::vector<float>* d0_eta = nullptr; std::vector<int>* d0_jet_idx = nullptr; std::vector<float>* d0_jet_dr = nullptr;
    std::vector<float>* mc_jet_pt = nullptr; std::vector<float>* mc_jet_eta = nullptr; std::vector<float>* mc_jet_phi = nullptr;
    std::vector<float>* mc_d0_pt = nullptr; std::vector<float>* mc_d0_eta = nullptr; std::vector<int>* mc_d0_matched = nullptr; std::vector<int>* mc_d0_jet_idx = nullptr; std::vector<float>* mc_d0_jet_dr = nullptr;

    t->SetBranchAddress("jet_pt", &jet_pt);
    t->SetBranchAddress("jet_eta", &jet_eta);
    t->SetBranchAddress("jet_phi", &jet_phi);

    t->SetBranchAddress("d0_pt", &d0_pt);
    t->SetBranchAddress("d0_eta", &d0_eta);
    t->SetBranchAddress("d0_jet_idx", &d0_jet_idx);
    t->SetBranchAddress("d0_jet_dr", &d0_jet_dr);

    t->SetBranchAddress("mc_jet_pt", &mc_jet_pt);
    t->SetBranchAddress("mc_jet_eta", &mc_jet_eta);
    t->SetBranchAddress("mc_jet_phi", &mc_jet_phi);

    t->SetBranchAddress("mc_d0_pt", &mc_d0_pt);
    t->SetBranchAddress("mc_d0_eta", &mc_d0_eta);
    t->SetBranchAddress("mc_d0_matched", &mc_d0_matched);
    t->SetBranchAddress("mc_d0_jet_idx", &mc_d0_jet_idx);
    t->SetBranchAddress("mc_d0_jet_dr", &mc_d0_jet_dr);

    // Binning: simple fixed binning
    std::vector<double> jetPtBins;
    for (double pt = 0.0; pt <= 60.0; pt += 2.5)
        jetPtBins.push_back(pt);
    std::vector<double> jetEtaBins;
    for (double eta = 2.0; eta <= 4.5; eta += 0.2)
        jetEtaBins.push_back(eta);

    // Create zT vs jetPt numerator/denominator
    const int nZTbins = 20;
    const double zTmin = 0.0;
    const double zTmax = 1.0;
    TH2F *h_num_zT = new TH2F("h_num_zT", "Numerator (zT vs jetPt);z_{T};Jet p_{T}", nZTbins, zTmin, zTmax, jetPtBins.size()-1, &jetPtBins[0]);
    TH2F *h_den_zT = new TH2F("h_den_zT", "Denominator (zT vs jetPt);z_{T};Jet p_{T}", nZTbins, zTmin, zTmax, jetPtBins.size()-1, &jetPtBins[0]);

    Long64_t nEntries = t->GetEntries();
    // nEntries/=10;
    std::cout << "Entries: " << nEntries << "\n";
    Long64_t printInterval = nEntries / 20;
    if (printInterval < 1) printInterval = 1;

    for (Long64_t ie = 0; ie < nEntries; ++ie) {
        t->GetEntry(ie);
        if (ie % printInterval == 0) {
            int pct = (int)std::ceil(100.0 * ie / (double)std::max<Long64_t>(1, nEntries));
            std::cout << "Processing entry " << ie << "/" << nEntries << " (" << pct << "%)" << std::endl;
        }
        const size_t nMCJets = mc_jet_pt ? mc_jet_pt->size() : 0;
        const size_t nRecoJets = jet_pt ? jet_pt->size() : 0;

        // prepare MC per-jet max D0 (generator) info
        std::vector<double> mc_maxD0(nMCJets, -1.0);
        for (size_t i=0;i<(mc_d0_jet_idx?mc_d0_jet_idx->size():0);++i) {
            int mcj = mc_d0_jet_idx->at(i);
            if (mcj<0 || mcj>=(int)nMCJets) continue;
            if (!mc_d0_jet_dr) continue;
            if (mc_d0_jet_dr->at(i) < 0 || mc_d0_jet_dr->at(i) >= gJetRadius) continue;
            double pt = mc_d0_pt ? mc_d0_pt->at(i) : -1.0;
            double eta = mc_d0_eta ? mc_d0_eta->at(i) : -999.0;
            if (!PassesD0Selection(pt, eta)) continue;
            if (pt > mc_maxD0[mcj]) mc_maxD0[mcj] = pt;
        }

        // prepare reco per-jet max D0
        std::vector<double> reco_maxD0(nRecoJets, -1.0);
        for (size_t i=0;i<(d0_jet_idx?d0_jet_idx->size():0);++i) {
            int rj = d0_jet_idx->at(i);
            if (rj<0 || rj>=(int)nRecoJets) continue;
            if (!d0_jet_dr) continue;
            if (d0_jet_dr->at(i) < 0 || d0_jet_dr->at(i) >= gJetRadius) continue;
            double pt = d0_pt ? d0_pt->at(i) : -1.0;
            double eta = d0_eta ? d0_eta->at(i) : -999.0;
            if (!PassesD0Selection(pt, eta)) continue;
            if (pt > reco_maxD0[rj]) reco_maxD0[rj] = pt;
        }

        // denominator: MC jets with a generator D0 in radius and passing selection
        // numerator (filled here): if any reco jet matches this MC jet (deltaR <= gJetRadius)
        for (size_t mcj=0; mcj<nMCJets; ++mcj) {
            if (mc_maxD0[mcj] <= 0) continue;
            if (!PassesJetSelectionMC(mc_jet_pt, mc_jet_eta, mcj)) continue;
            double mcPt = mc_jet_pt->at(mcj);
            double mcEta = mc_jet_eta->at(mcj);
            double mcPhi = mc_jet_phi ? mc_jet_phi->at(mcj) : 0.0;
            // fill denominator zT vs jetPt
            double zT = mc_maxD0[mcj] / mcPt;
            if (zT > 0 && zT < zTmax) h_den_zT->Fill(zT, mcPt);

            // now search for a reconstructed jet that matches this MC jet within deltaR and contains a reconstructed D0
            bool matched = false;
            for (size_t rj=0; rj<nRecoJets; ++rj) {
                // Do not require the reco jet to contain a reconstructed D0 here;
                // we only require a reconstructed jet that matches the MC jet by geometry.
                if (!PassesJetSelectionReco(jet_pt, jet_eta, rj)) continue;
                double rEta = jet_eta ? jet_eta->at(rj) : 0.0;
                double rPhi = jet_phi ? jet_phi->at(rj) : 0.0;
                double dr = DeltaR(mcEta, mcPhi, rEta, rPhi);
                if (dr <= gJetRadius) { matched = true; break; }
            }
            if (matched) {
                if (zT > 0 && zT < zTmax) h_num_zT->Fill(zT, mcPt);
            }
        }

    }

    // create efficiency map (zT vs jetPt)
    TH2F *h_eff = (TH2F*)h_num_zT->Clone("jet_finding_efficiency_zT");
    h_eff->Divide(h_den_zT);
    h_eff->SetTitle("Jet Finding Efficiency vs z_{T};z_{T};Jet p_{T} [GeV]");

    // --- Closure test: apply efficiency to MC denominator to predict reco yield ---
    TH2F *h_pred = (TH2F*)h_den_zT->Clone("predicted_reco_zT");
    h_pred->SetTitle("Predicted reco (den * eff);z_{T};Jet p_{T}");
    TH2F *h_diff = (TH2F*)h_num_zT->Clone("diff_reco_minus_pred");
    h_diff->SetTitle("Reco - Predicted;z_{T};Jet p_{T}");

    double totalNum = 0.0, totalPred = 0.0;
    double chi2 = 0.0; int nBinsUsed = 0;
    int nx = h_eff->GetNbinsX();
    int ny = h_eff->GetNbinsY();
    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            double den = h_den_zT->GetBinContent(ix, iy);
            double eff = h_eff->GetBinContent(ix, iy);
            double pred = den * eff;
            double num = h_num_zT->GetBinContent(ix, iy);
            h_pred->SetBinContent(ix, iy, pred);
            h_diff->SetBinContent(ix, iy, num - pred);
            totalNum += num;
            totalPred += pred;
            double var = pred + 1e-9; // avoid zero
            chi2 += (num - pred) * (num - pred) / var;
            if (den > 0) ++nBinsUsed;
        }
    }

    double relDiffPercent = (totalNum > 0) ? 100.0 * (totalPred - totalNum) / totalNum : 0.0;
    double chi2perBin = (nBinsUsed>0) ? chi2 / nBinsUsed : 0.0;
    std::cout << "Closure test: total reco (data) = " << totalNum << ", total predicted = " << totalPred
              << ", relative diff (%) = " << relDiffPercent << ", chi2/nbins = " << chi2perBin << std::endl;


    // draw and save
    TCanvas c("c_eff", "Jet finding efficiency (zT vs pT)", 800,600);
    c.SetRightMargin(0.2);
    h_eff->SetStats(0);
    h_eff->GetZaxis()->SetTitle("Efficiency");
    h_eff->Draw("COLZ");
    TLatex tex; tex.SetNDC(); tex.SetTextSize(0.04);
    tex.DrawLatex(0.15,0.92, "Jet Finding Efficiency vs z_{T}");
    std::string png = outDir + "/" + outPrefix + "_eff_zT.png";
    std::string pdf = outDir + "/" + outPrefix + "_eff_zT.pdf";
    savePNG(c, png);
    c.SaveAs(pdf.c_str());

    std::cout << "Wrote: " << png << " " << pdf << "\n";
    // save predicted and residual maps for a quick closure diagnostic
    TCanvas c2("c_pred", "Predicted reco (zT vs pT)", 800,600);
    c2.SetRightMargin(0.2);
    h_pred->SetStats(0);
    h_pred->GetZaxis()->SetTitle("Counts");
    h_pred->Draw("COLZ");
    std::string png_pred = outDir + "/" + outPrefix + "_pred_zT.png";
    std::string pdf_pred = outDir + "/" + outPrefix + "_pred_zT.pdf";
    savePNG(c2, png_pred);
    c2.SaveAs(pdf_pred.c_str());

    TCanvas c3("c_res", "Reco - Predicted (zT vs pT)", 800,600);
    c3.SetRightMargin(0.2);
    h_diff->SetStats(0);
    h_diff->GetZaxis()->SetTitle("Counts (data - pred)");
    h_diff->Draw("COLZ");
    std::string png_res = outDir + "/" + outPrefix + "_res_zT.png";
    std::string pdf_res = outDir + "/" + outPrefix + "_res_zT.pdf";
    savePNG(c3, png_res);
    c3.SaveAs(pdf_res.c_str());

    std::cout << "Wrote: " << png_pred << " " << pdf_pred << " and " << png_res << " " << pdf_res << "\n";

    // also write histograms to a ROOT file in the output directory
    std::string rootOut = outDir + "/" + outPrefix + ".root";
    TFile *fout = TFile::Open(rootOut.c_str(), "RECREATE");
    if (fout && !fout->IsZombie()) {
        h_num_zT->Write();
        h_den_zT->Write();
        h_eff->Write();
        h_pred->Write();
        h_diff->Write();
        fout->Close();
        delete fout;
        std::cout << "Wrote root file: " << rootOut << "\n";
    } else {
        std::cerr << "Failed to create output root file '" << rootOut << "'\n";
    }

    // --- per-jetPt-bin efficiency and closure checks ---
    std::vector<double> startPt = {5, 10, 15, 20, 30};
    std::vector<double> endPt = {10, 15, 20, 30, 50};
    if (startPt.size() != endPt.size()) {
        std::cerr << "Pt bin vectors size mismatch\n";
    } else {
        // reopen root file in UPDATE to append per-bin projections
        TFile *fout_upd = TFile::Open(rootOut.c_str(), "UPDATE");
        for (size_t ib=0; ib<startPt.size(); ++ib) {
            double s = startPt[ib];
            double e = endPt[ib];
            int firsty = h_den_zT->GetYaxis()->FindBin(s + 1e-6);
            int lasty = h_den_zT->GetYaxis()->FindBin(e - 1e-6);
            if (lasty < firsty) lasty = firsty;

            TH1D *den_proj = h_den_zT->ProjectionX(Form("den_pt_%.0f_%.0f", s, e), firsty, lasty);
            TH1D *num_proj = h_num_zT->ProjectionX(Form("num_pt_%.0f_%.0f", s, e), firsty, lasty);
            TH1D *pred_proj = h_pred->ProjectionX(Form("pred_pt_%.0f_%.0f", s, e), firsty, lasty);
            TH1D *diff_proj = h_diff->ProjectionX(Form("diff_pt_%.0f_%.0f", s, e), firsty, lasty);

            TH1D *eff_proj = (TH1D*)num_proj->Clone(Form("eff_pt_%.0f_%.0f", s, e));
            eff_proj->SetTitle(Form("Jet finding efficiency zT (pt %.0f-%.0f);z_{T};Efficiency", s, e));
            //do binomial uncertainty on division
            eff_proj->Divide(eff_proj, den_proj, 1.0, 1.0, "B");

            // compute per-bin closure summary
            double tnum = num_proj->Integral();
            double tpred = pred_proj->Integral();
            double chi2b = 0.0; int nbinsb = 0;
            for (int bx=1; bx<=eff_proj->GetNbinsX(); ++bx) {
                double nval = num_proj->GetBinContent(bx);
                double pval = pred_proj->GetBinContent(bx);
                double var = pval==0 ? pval + 1e-9 : pval;
                chi2b += (nval - pval)*(nval - pval)/var;
                if (den_proj->GetBinContent(bx) > 0) ++nbinsb;
            }
            double rel = (tnum>0) ? 100.0 * (tpred - tnum) / tnum : 0.0;
            double chi2per = (nbinsb>0) ? chi2b/nbinsb : 0.0;
            std::cout << Form("Pt bin [%.0f,%.0f): total reco=%.1f, total pred=%.1f, rel diff(%%)=%.2f, chi2/nbins=%.3f", s, e, tnum, tpred, rel, chi2per) << std::endl;

            // draw efficiency vs zT and fit three simple models
            TCanvas cbin(Form("c_eff_pt_%.0f_%.0f", s, e), Form("Efficiency pt %.0f-%.0f", s, e), 700,500);
            eff_proj->SetStats(0);
            eff_proj->SetLineWidth(2);
            eff_proj->SetMarkerStyle(20);
            eff_proj->SetMarkerColor(kBlack);
            eff_proj->SetLineColor(kBlack);
            eff_proj->Draw("pe, E");

            // --- define fit functions ---
            // logistic / sigmoid: [0]/(1+exp(-[1]*(x-[2])))
            TF1 *f_log = new TF1(Form("f_log_pt%.0f_%.0f", s, e), "[0]/(1+exp(-[1]*(x-[2])))", zTmin, zTmax);
            f_log->SetParameters(0.85, 6.0, 0.6);
            f_log->SetParLimits(0, 0.0, 1.0);
            f_log->SetParLimits(1, 0.0, 100.0);
            f_log->SetParLimits(2, zTmin, zTmax);

            // saturating exponential: [0]+[1]*(1-exp(-[2]*x))
            TF1 *f_exp = new TF1(Form("f_exp_pt%.0f_%.0f", s, e), "[0]+[1]*(1-exp(-[2]*x))", zTmin, zTmax);
            f_exp->SetParameters(0.05, 0.8, 2.0);
            f_exp->SetParLimits(0, 0.0, 1.0);
            f_exp->SetParLimits(1, 0.0, 1.0);
            f_exp->SetParLimits(2, 1e-3, 100.0);

            // quadratic polynomial as a fallback
            TF1 *f_pol2 = new TF1(Form("f_pol2_pt%.0f_%.0f", s, e), "pol2", zTmin, zTmax);

            // Fit with weights (use errors in eff_proj)
            eff_proj->Fit(f_log, "R W Q 0");
            eff_proj->Fit(f_exp, "R W Q 0");
            eff_proj->Fit(f_pol2, "R W Q 0");

            // style and draw fits
            f_log->SetLineColor(kRed); f_log->SetLineWidth(2);
            f_exp->SetLineColor(kGreen+2); f_exp->SetLineWidth(2); f_exp->SetLineStyle(2);
            f_pol2->SetLineColor(kBlue+2); f_pol2->SetLineWidth(3); f_pol2->SetLineStyle(4);
            f_log->Draw("SAME");
            f_exp->Draw("SAME");
            f_pol2->Draw("SAME");

            // legend for fits
            TLegend legf(0.6,0.15,0.92,0.38);
            legf.AddEntry(eff_proj, "Efficiency (data)", "l");
            legf.AddEntry(f_log, "Logistic", "l");
            legf.AddEntry(f_exp, "Saturating exp.", "l");
            legf.AddEntry(f_pol2, "Pol2", "l");
            legf.SetBorderSize(0);
            legf.Draw();

            // print fit summaries
            auto printFit = [&](TF1* f, const char* name){
                double chi2f = f->GetChisquare();
                double ndf = f->GetNDF();
                std::cout << Form("Fit %s pt[%.0f,%.0f]: chi2=%.2f, ndf=%.0f, chi2/ndf=%.3f", name, s, e, chi2f, ndf, (ndf>0?chi2f/ndf:0.0)) << std::endl;
                for (int ip=0; ip<f->GetNpar(); ++ip) {
                    std::cout << Form("  %s: p%d = %.5g +/- %.5g", name, ip, f->GetParameter(ip), f->GetParError(ip)) << std::endl;
                }
            };
            printFit(f_log, "Logistic");
            printFit(f_exp, "Exp");
            printFit(f_pol2, "Pol2");

            std::string pngb = outDir + "/" + outPrefix + Form("_eff_zT_pt%.0f-%.0f.png", s, e);
            std::string pdfb = outDir + "/" + outPrefix + Form("_eff_zT_pt%.0f-%.0f.pdf", s, e);
            savePNG(cbin, pngb); cbin.SaveAs(pdfb.c_str());

            // draw closure overlay: num vs pred
            TCanvas ccl(Form("c_closure_pt_%.0f_%.0f", s, e), Form("Closure pt %.0f-%.0f", s, e), 700,500);
            num_proj->SetLineColor(kBlue); num_proj->SetLineWidth(2);
            pred_proj->SetLineColor(kRed); pred_proj->SetLineWidth(2); pred_proj->SetLineStyle(2);
            num_proj->Draw("HIST");
            pred_proj->Draw("HIST SAME");
            TLegend leg(0.6,0.7,0.88,0.88);
            leg.AddEntry(num_proj, "Reco (data)", "l");
            leg.AddEntry(pred_proj, "Predicted (den*eff)", "l");
            leg.Draw();
            std::string pngc = outDir + "/" + outPrefix + Form("_closure_zT_pt%.0f-%.0f.png", s, e);
            std::string pdfc = outDir + "/" + outPrefix + Form("_closure_zT_pt%.0f-%.0f.pdf", s, e);
            savePNG(ccl, pngc); ccl.SaveAs(pdfc.c_str());

            // residual
            TCanvas cres(Form("c_res_pt_%.0f_%.0f", s, e), Form("Residual pt %.0f-%.0f", s, e), 700,500);
            diff_proj->SetLineColor(kBlack); diff_proj->SetLineWidth(2);
            diff_proj->Draw("HIST");
            std::string pngr = outDir + "/" + outPrefix + Form("_res_zT_pt%.0f-%.0f.png", s, e);
            std::string pdfr = outDir + "/" + outPrefix + Form("_res_zT_pt%.0f-%.0f.pdf", s, e);
            savePNG(cres, pngr); cres.SaveAs(pdfr.c_str());

            // write these projections and fits into root file
            if (fout_upd && !fout_upd->IsZombie()) {
                den_proj->Write(); num_proj->Write(); pred_proj->Write(); diff_proj->Write(); eff_proj->Write();
                f_log->Write(); f_exp->Write(); f_pol2->Write();
            }

            // cleanup projections and fit objects
            delete den_proj; delete num_proj; delete pred_proj; delete diff_proj; delete eff_proj;
            delete f_log; delete f_exp; delete f_pol2;
        }
        if (fout_upd && !fout_upd->IsZombie()) { fout_upd->Close(); delete fout_upd; }
    }

    delete h_num_zT; delete h_den_zT; delete h_eff; delete h_pred; delete h_diff;
    f->Close(); delete f;
    return 0;
}
