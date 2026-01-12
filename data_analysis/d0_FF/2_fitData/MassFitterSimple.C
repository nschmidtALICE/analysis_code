
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <memory>
#include <filesystem>
#include <iomanip>
#include <chrono>
#include <thread>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGaxis.h"
#include "TF1.h"
#include "TPad.h"
#include "TSystem.h"

// RooFit includes
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooFit.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooFitResult.h"
#include "RooHist.h"
#include "RooFormulaVar.h"
#include "RooBukinPdf.h"
#include "RooStats/SPlot.h"
#include <fstream>
#include <sstream>

// Include external fitter and plotter headers
#include "PlotHelpers.h"

//==============================================================================
void MassFitterSimple(
    TString inputFile = "/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/59_FF_Pbp_DATA_filtered.root",
    // TString inputFile = "/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/57_FF_pPb_DATA_filtered.root",
    bool isMC = false, bool doPlots = true)
{
    // Simple mass fitting function implementation
    std::cout << "Running simple mass fitter on file: " << inputFile.Data() << " with isMC = " << isMC << std::endl;

    // define the rapidity and jet pT bins
    std::vector<double> yBins = {2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0};
    std::vector<double> jetPtBins = {5, 10, 15, 20, 30, 50};

    // Open input as TTree (FragmNtuple) or TChain from a list
    TTree *tree = nullptr;
    std::unique_ptr<TFile> inFile;
    std::unique_ptr<TChain> chain;

    auto isListPath = [&](const TString &p)
    { return p.EndsWith(".txt") || p.EndsWith(".list"); };
    if (isListPath(inputFile))
    {
        chain = std::make_unique<TChain>("FragmNtuple");
        std::ifstream lst(inputFile.Data());
        if (!lst.is_open())
        {
            std::cerr << "[ERROR] Cannot open list file: " << inputFile << std::endl;
            return;
        }
        std::string line;
        int added = 0;
        while (std::getline(lst, line))
        {
            // trim
            auto ltrim = [](std::string &s)
            { s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char c)
                                              { return !std::isspace(c); })); };
            auto rtrim = [](std::string &s)
            { s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char c)
                                   { return !std::isspace(c); })
                          .base(),
                      s.end()); };
            ltrim(line);
            rtrim(line);
            if (line.empty() || line[0] == '#')
                continue;
            int n = chain->Add(line.c_str());
            if (n > 0)
            {
                ++added;
            }
            else
            {
                std::cerr << "[WARN] Failed to add to chain: " << line << std::endl;
            }
        }
        if (added == 0)
        {
            std::cerr << "[ERROR] No files added to chain from list: " << inputFile << std::endl;
            return;
        }
        tree = chain.get();
    }
    else
    {
        inFile = std::make_unique<TFile>(inputFile, "READ");
        if (!inFile || inFile->IsZombie())
        {
            std::cerr << "[ERROR] Cannot open input file: " << inputFile << std::endl;
            return;
        }
        tree = dynamic_cast<TTree *>(inFile->Get("FragmNtuple"));
        if (!tree)
        {
            std::cerr << "[ERROR] Could not find TTree 'FragmNtuple' in input file" << std::endl;
            return;
        }
    }

    // Check required branches
    if (!tree->GetBranch("tagMass") || !tree->GetBranch("tagY"))
    {
        std::cerr << "[ERROR] Required branches 'tagMass' or 'tagY' are missing in input tree" << std::endl;
        return;
    }
    bool haveJetPt = (tree->GetBranch("jetPt") != nullptr);
    if (!haveJetPt)
    {
        std::cerr << "[WARN] Branch 'jetPt' not found — jet pT binning will be ignored" << std::endl;
    }
    // IP chi2 branches
    bool haveLogIP = (tree->GetBranch("log_tag_ipchi2") != nullptr);
    bool haveIP = (tree->GetBranch("tag_ip_chi2") != nullptr);
    if (!haveLogIP && !haveIP)
    {
        std::cerr << "[WARN] IP chi2 branches ('log_tag_ipchi2' or 'tag_ip_chi2') not found — IP chi2 fits will be skipped" << std::endl;
    }

    // Define RooFit variables
    const int kZTBins = 50;            // zT histogram bins
    const double kZTMin = 0.0;         // zT histogram min
    const double kZTMax = 1.0;         // zT histogram max
    RooRealVar r_mass("tagMass", "D^{0} mass [GeV]", 1.815, 1.925);
    RooRealVar r_y("tagY", "y(D^{0})", yBins.front(), yBins.back());
    RooRealVar r_jetPt("jetPt", "p_{T}^{jet} [GeV]", jetPtBins.front(), jetPtBins.back());
    // Optional zT observable (for weighted histograms). Use wide range to avoid accidental cuts.
    bool haveZ = (tree->GetBranch("tagZ") != nullptr);
    RooRealVar r_zT("tagZ", "z_{T}", -0.1, 1.5);
    RooArgSet obs(r_mass, r_y);
    if (haveJetPt)
        obs.add(r_jetPt);
    if (haveZ)
        obs.add(r_zT);
    else
        std::cerr << "[WARN] Branch 'tagZ' not found — weighted zT histograms will be skipped" << std::endl;

    // Prepare output directory for plots using helper
    // Build a meaningful output directory: <YYYYMMDD>_D0_FF_<DATA|MC>
    auto now = std::chrono::system_clock::now();
    std::time_t now_t = std::chrono::system_clock::to_time_t(now);
    std::tm tm_now;
    localtime_r(&now_t, &tm_now);
    std::ostringstream datebuf;
    datebuf << std::put_time(&tm_now, "%Y%m%d");
    std::string dataLabel = isMC ? "MC" : "DATA";
    // Check inputFile name for pPb/Pbp and append to output directory if present
    std::string inName = inputFile.Data();
    std::string runLabel = "";
    if (inName.find("pPb") != std::string::npos)
        runLabel = "_pPb";
    else if (inName.find("Pbp") != std::string::npos)
        runLabel = "_Pbp";
    std::string baseOutDir = datebuf.str() + std::string("_D0_FF_") + dataLabel + runLabel;
    ensureDirectoryExists(baseOutDir);
    debugLog(std::string("Output directory: ") + baseOutDir);

    // Prepare an output ROOT file to store distributions and fit curves
    std::string fitsNameBase = isMC ? "DistributionsAndFits_D0_MC" : "DistributionsAndFits_D0_DATA";
    std::string fitsPath = baseOutDir + "/" + fitsNameBase + runLabel + ".root";
    std::unique_ptr<TFile> fitFile(TFile::Open(fitsPath.c_str(), "RECREATE"));
    if (!fitFile || fitFile->IsZombie())
    {
        std::cerr << "[ERROR] Cannot create output ROOT file: " << fitsPath << std::endl;
        return;
    }

    // Prepare an output ROOT file to store sPlot weights
    std::string sWNameBase = isMC ? "SWeights_D0_MC" : "SWeights_D0_DATA";
    std::string sWPath = baseOutDir + "/" + sWNameBase + runLabel + ".root";
    std::unique_ptr<TFile> sWFile(TFile::Open(sWPath.c_str(), "RECREATE"));
    if (!sWFile || sWFile->IsZombie())
    {
        std::cerr << "[ERROR] Cannot create sWeights ROOT file: " << sWPath << std::endl;
        return;
    }

    // Reduce RooFit verbosity for speed/readability
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    // Ensure RooStats library is available for SPlot
    gSystem->Load("libRooStats");

    // ------------------------------------------------------------
    // Pre-scan: estimate total number of fits to be performed
    // Mass fit always if bin has entries; IP chi2 fit only if IP branches exist
    // Always count fits for bins with entries (no resume/skip based on existing PNGs)
    int totalFitsPlanned = 0;
    const bool haveIPBranches = (haveLogIP || haveIP);
    for (size_t iy = 0; iy + 1 < yBins.size(); ++iy)
    {
        double yLo = yBins[iy];
        double yHi = yBins[iy + 1];
        size_t jStart = 0, jEnd = (haveJetPt ? jetPtBins.size() - 1 : 1);
        for (size_t ij = jStart; ij < jEnd; ++ij)
        {
            double jLo = haveJetPt ? jetPtBins[ij] : jetPtBins.front();
            double jHi = haveJetPt ? jetPtBins[ij + 1] : jetPtBins.back();

            std::ostringstream cutPS;
            cutPS << std::fixed << std::setprecision(6);
            cutPS << "tagY>=" << yLo << " && tagY<" << yHi;
            if (haveJetPt)
                cutPS << " && jetPt>=" << jLo << " && jetPt<" << jHi;

            Long64_t nSel = tree->GetEntries(cutPS.str().c_str());
            if (nSel <= 0)
                continue;

            // Count mass fit
            ++totalFitsPlanned;
            // Count IP fit if branches exist
            if (haveIPBranches)
                ++totalFitsPlanned;
        }
    }
    int fitsDone = 0;

    // Loop over bins and create RooDataSets
    for (size_t iy = 0; iy + 1 < yBins.size(); ++iy)
    {
        double yLo = yBins[iy];
        double yHi = yBins[iy + 1];

        // either loop jet pT bins if available, else single inclusive bin
        size_t jStart = 0, jEnd = (haveJetPt ? jetPtBins.size() - 1 : 1);
        for (size_t ij = jStart; ij < jEnd; ++ij)
        {
            double jLo = haveJetPt ? jetPtBins[ij] : jetPtBins.front();
            double jHi = haveJetPt ? jetPtBins[ij + 1] : jetPtBins.back();

            // Build cut string
            std::ostringstream cut;
            cut << std::fixed << std::setprecision(6);
            cut << "tagY>=" << yLo << " && tagY<" << yHi;
            if (haveJetPt)
            {
                cut << " && jetPt>=" << jLo << " && jetPt<" << jHi;
            }

            // Name the dataset
            auto mkname = [&](const char *base)
            {
                std::ostringstream nm; nm<<base<<"_y"<<iy<<"_"<<yLo<<"to"<<yHi; if(haveJetPt) nm<<"_j"<<ij<<"_"<<jLo<<"to"<<jHi; return nm.str(); };
            std::string dsName = mkname("data");
            std::string dsTitle = mkname("D0");

            // Create dataset from tree with selection
            RooDataSet ds(dsName.c_str(), dsTitle.c_str(), tree, obs, cut.str().c_str());
            std::cout << "Created dataset " << dsName << " with entries: " << ds.numEntries() << std::endl;

            // No workspace import: operate on the transient dataset directly

            // Skip empty datasets
            if (ds.numEntries() == 0)
            {
                std::cerr << "[INFO] Skipping empty dataset for bin y:[" << yLo << "," << yHi << "]"
                          << (haveJetPt ? " jetPt:[" : " ") << (haveJetPt ? std::to_string(jLo) : "")
                          << (haveJetPt ? ", " : "") << (haveJetPt ? std::to_string(jHi) : "") << std::endl;
                continue;
            }

            // If requested, perform mass fit and plots. Support resume: save under per-jetPt directories and skip if plot exists.
            std::string outDir = baseOutDir;
            if (haveJetPt)
            {
                std::ostringstream jd;
                jd << baseOutDir << "/j" << ij << "_" << jLo << "to" << jHi;
                outDir = jd.str();
                ensureDirectoryExists(outDir);
            }

            // No workspace saving: keep runtime lean

            if (!doPlots)
            {
                continue;
            }

            // ------------------ Mass fit model: double Gaussian signal + exponential background (extended) ------------------
            // Parameters (unique per bin)
            std::ostringstream tag;
            tag << "_y" << iy << "_j" << (haveJetPt ? (int)ij : -1);
            std::string suf = tag.str();

            RooRealVar mean(("mean" + suf).c_str(), "mean", 1.864, 1.82, 1.90);
            RooRealVar sigma1(("sigma1" + suf).c_str(), "sigma1", 0.010, 0.001, 0.050);
            RooRealVar sigma2(("sigma2" + suf).c_str(), "sigma2", 0.020, 0.002, 0.080);
            RooRealVar fracG(("fracG" + suf).c_str(), "fracG", 0.6, 0.0, 1.0);

            RooGaussian g1(("g1" + suf).c_str(), "G1", r_mass, mean, sigma1);
            RooGaussian g2(("g2" + suf).c_str(), "G2", r_mass, mean, sigma2);
            RooAddPdf sig(("sig" + suf).c_str(), "Signal", RooArgList(g1, g2), RooArgList(fracG));

            RooRealVar tau(("tau" + suf).c_str(), "tau", -2.0, -50.0, -0.001);
            RooExponential bkg(("bkg" + suf).c_str(), "Background", r_mass, tau);

            double nTot = std::max(1, (int)ds.numEntries());
            RooRealVar nSig(("nSig" + suf).c_str(), "N_{sig}", 0.6 * nTot, 0.0, 2.0 * nTot);
            RooRealVar nBkg(("nBkg" + suf).c_str(), "N_{bkg}", 0.1 * nTot, 0.0, 2.0 * nTot);

            RooAddPdf model(("model" + suf).c_str(), "Model", RooArgList(sig, bkg), RooArgList(nSig, nBkg));

            // Resume/skip if plot already exists
            std::ostringstream onBase;
            onBase << outDir << "/massfit_y" << iy << "_" << yLo << "to" << yHi;
            if (haveJetPt)
                onBase << "_j" << ij << "_" << jLo << "to" << jHi;
            std::string massPng = onBase.str() + ".png";

            // Always perform mass fit (ensures ROOT outputs reflect fit results)
            RooFitResult *fitRes = model.fitTo(ds, RooFit::Save(true), RooFit::Extended(true), RooFit::PrintLevel(-1), RooFit::PrintEvalErrors(-1));

            // Plot with pulls: create a two-pad canvas
            std::unique_ptr<TCanvas> c(new TCanvas(("c_" + dsName).c_str(), ("Mass fit " + dsTitle).c_str(), 800, 700));
            c->Divide(1, 2, 0, 0);
            TPad *pTop = (TPad *)c->cd(1);
            pTop->SetPad(0, 0.28, 1, 1);
            setupPadMargins(pTop, 0.12, 0.02, 0.01, 0.01);
            TPad *pBot = (TPad *)c->cd(2);
            pBot->SetPad(0, 0, 1, 0.28);
            setupPadMargins(pBot, 0.12, 0.35, 0.01, 0.05);

            // Top: data and fit
            pTop->cd();
            RooPlot *frame = r_mass.frame(RooFit::Bins(100));
            // Name objects so we can build a legend
            std::string nmData = std::string("data") + suf;
            std::string nmModel = std::string("modelCurve") + suf;
            std::string nmBkg = std::string("bkgCurve") + suf;
            std::string nmSig = std::string("sigCurve") + suf;
            ds.plotOn(frame, RooFit::MarkerSize(0.8), RooFit::Name(nmData.c_str()));
            model.plotOn(frame, RooFit::LineColor(kRed + 1), RooFit::LineWidth(2), RooFit::Name(nmModel.c_str()));
            // background component overlay
            model.plotOn(frame, RooFit::Components(bkg), RooFit::LineStyle(kDashed), RooFit::LineColor(kGray + 2), RooFit::Name(nmBkg.c_str()));
            // signal component overlay
            model.plotOn(frame, RooFit::Components(sig), RooFit::LineStyle(kDashed), RooFit::LineColor(kBlue + 1), RooFit::Name(nmSig.c_str()));
            frame->SetTitle("");
            frame->GetXaxis()->SetLabelOffset(999); // hide x labels on top pad
            frame->GetXaxis()->SetTitleOffset(999);
            frame->Draw();

            // Top-left bin label (rapidity and jet pT)
            {
                TPaveText *pt = new TPaveText(0.15, 0.78, 0.48, 0.92, "NDC");
                pt->SetFillColor(0);
                pt->SetFillStyle(0);
                pt->SetBorderSize(0);
                pt->SetTextAlign(13); // left-top
                pt->SetTextFont(42);
                pt->SetTextSize(0.040);
                std::ostringstream lblY;
                lblY << std::fixed << std::setprecision(1) << yLo << " < y < " << yHi;
                pt->AddText(lblY.str().c_str());
                if (haveJetPt)
                {
                    std::ostringstream lblJ;
                    lblJ << std::fixed << std::setprecision(0) << jLo << " < p_{T}^{jet} < " << jHi << " GeV";
                    pt->AddText(lblJ.str().c_str());
                }
                else
                {
                    pt->AddText("jet pT: inclusive");
                }
                pt->Draw();
            }

            // Legend with entries and fit results
            int nfloat = (fitRes ? fitRes->floatParsFinal().getSize() : 6);
            double chi2ndf = frame->chiSquare(nfloat);
            TLegend *leg = new TLegend(0.62, 0.55, 0.92, 0.93);
            setupLegend(leg, 0.035, 0.15, 42, 0, 0, 0);
            TObject *objData = frame->findObject(nmData.c_str());
            TObject *objModel = frame->findObject(nmModel.c_str());
            TObject *objBkg = frame->findObject(nmBkg.c_str());
            if (objData)
                leg->AddEntry(objData, "Data", "lep");
            if (objModel)
                leg->AddEntry(objModel, "Fit", "l");
            if (objBkg)
                leg->AddEntry(objBkg, "Bkg (exp)", "l");
            auto makeLine = [&](const std::string &label, const RooRealVar &v)
            {
                std::ostringstream os; os<<label<<" = "<<std::fixed<<std::setprecision(4)<<v.getVal()<<" +/- "<<v.getError(); return os.str(); };
            leg->AddEntry((TObject *)nullptr, makeLine("mean", mean).c_str(), "");
            leg->AddEntry((TObject *)nullptr, makeLine("sigma1", sigma1).c_str(), "");
            leg->AddEntry((TObject *)nullptr, makeLine("sigma2", sigma2).c_str(), "");
            leg->AddEntry((TObject *)nullptr, makeLine("fracG", fracG).c_str(), "");
            leg->AddEntry((TObject *)nullptr, makeLine("tau", tau).c_str(), "");
            {
                std::ostringstream os;
                os << "Nsig = " << std::fixed << std::setprecision(1) << nSig.getVal() << " +/- " << nSig.getError();
                leg->AddEntry((TObject *)nullptr, os.str().c_str(), "");
            }
            {
                std::ostringstream os;
                os << "Nbkg = " << std::fixed << std::setprecision(1) << nBkg.getVal() << " +/- " << nBkg.getError();
                leg->AddEntry((TObject *)nullptr, os.str().c_str(), "");
            }
            {
                std::ostringstream os;
                os << "chi2/ndf = " << std::setprecision(3) << chi2ndf;
                leg->AddEntry((TObject *)nullptr, os.str().c_str(), "");
            }
            leg->Draw();

            // Bottom: pull distribution
            pBot->cd();
            RooHist *hpull = frame->pullHist();
            RooPlot *pullFrame = r_mass.frame(RooFit::Bins(100));
            pullFrame->addPlotable(hpull, "B");
            pullFrame->SetTitle("");
            pullFrame->GetYaxis()->SetTitle("Pull");
            pullFrame->GetYaxis()->SetTitleSize(0.10);
            pullFrame->GetYaxis()->SetTitleOffset(0.4);
            pullFrame->GetYaxis()->SetLabelSize(0.09);
            pullFrame->GetXaxis()->SetTitleSize(0.12);
            pullFrame->GetXaxis()->SetLabelSize(0.10);
            pullFrame->GetXaxis()->SetTitle("m_{D^{0}} [GeV]");
            pullFrame->Draw();

            // Save figure (PNG) always (overwrite if exists)
            c->SaveAs(massPng.c_str());
            // Progress status: always increment after mass fit
            if (totalFitsPlanned)
            {
                ++fitsDone;
                std::cout << "[Status] Fit " << fitsDone << "/" << totalFitsPlanned
                          << " done: Mass y=[" << yLo << "," << yHi << "]"
                          << (haveJetPt ? (std::string(" jetPt=[") + std::to_string(jLo) + "," + std::to_string(jHi) + "]") : std::string(" inclusive"))
                          << ". Remaining: " << (totalFitsPlanned - fitsDone) << std::endl;
            }

            // Save mass distribution and fit curves to ROOT file
            r_mass.setBins(100);
            TH1 *hMass = ds.createHistogram(("hMass" + suf).c_str(), r_mass);
            if (hMass)
            {
                fitFile->cd();
                hMass->Write();
            }
            // Write curves if present
            if (auto objM = frame->findObject(nmModel.c_str()))
            {
                fitFile->cd();
                objM->Write(("curveMassModel" + suf).c_str());
            }
            if (auto objB = frame->findObject(nmBkg.c_str()))
            {
                fitFile->cd();
                objB->Write(("curveMassBkg" + suf).c_str());
            }
            if (auto objS = frame->findObject(nmSig.c_str()))
            {
                fitFile->cd();
                objS->Write(("curveMassSig" + suf).c_str());
            }
            // Also store canvas for convenience
            fitFile->cd();
            c->Write(("cMass" + suf).c_str());
            // std::string pdfName = on.str(); pdfName.replace(pdfName.end()-3, pdfName.end(), "pdf");
            // c->SaveAs(pdfName.c_str());

            // -------- sPlot for mass fit: store per-event weights (signal/background) --------
            try
            {
                RooStats::SPlot sData((std::string("sData") + suf).c_str(), "sPlot mass", ds, &model, RooArgList(nSig, nBkg));
                sWFile->cd();
                std::string tname = std::string("sWeights_mass") + suf;
                TTree t(tname.c_str(), ("sWeights for mass " + dsTitle).c_str());
                double outY = 0.0, outJet = 0.0, outM = 0.0, wSig = 0.0, wBkg = 0.0;
                t.Branch("tagY", &outY, "tagY/D");
                if (haveJetPt)
                    t.Branch("jetPt", &outJet, "jetPt/D");
                t.Branch("tagMass", &outM, "tagMass/D");
                t.Branch("w_sig", &wSig, "w_sig/D");
                t.Branch("w_bkg", &wBkg, "w_bkg/D");
                const std::string wSigName = std::string(nSig.GetName()) + "_sw";
                const std::string wBkgName = std::string(nBkg.GetName()) + "_sw";
                for (int i = 0; i < ds.numEntries(); ++i)
                {
                    const RooArgSet *row = ds.get(i);
                    outY = row->getRealValue("tagY");
                    if (haveJetPt)
                        outJet = row->getRealValue("jetPt");
                    outM = row->getRealValue("tagMass");
                    wSig = row->getRealValue(wSigName.c_str());
                    wBkg = row->getRealValue(wBkgName.c_str());
                    t.Fill();
                }
                t.Write();

                // Also fill weighted zT histograms directly (no separate macro needed)
                if (haveZ)
                {
                    std::unique_ptr<TH1D> hZTsig(new TH1D(("hZT_sig" + suf).c_str(), ("zT signal " + dsTitle).c_str(), kZTBins, kZTMin, kZTMax));
                    std::unique_ptr<TH1D> hZTbkg(new TH1D(("hZT_bkg" + suf).c_str(), ("zT background " + dsTitle).c_str(), kZTBins, kZTMin, kZTMax));
                    std::unique_ptr<TH1D> hZTRaw(new TH1D(("hZT_raw" + suf).c_str(), ("zT raw " + dsTitle).c_str(), kZTBins, kZTMin, kZTMax));
                    for (int i = 0; i < ds.numEntries(); ++i)
                    {
                        const RooArgSet *row = ds.get(i);
                        double z = row->getRealValue("tagZ");
                        double ws = row->getRealValue(wSigName.c_str());
                        double wb = row->getRealValue(wBkgName.c_str());
                        hZTRaw->Fill(z);          // unweighted
                        hZTsig->Fill(z, ws);      // signal-weighted
                        hZTbkg->Fill(z, wb);      // background-weighted (may be useful later)
                    }
                    // Save to the distributions/fit file
                    fitFile->cd();
                    hZTsig->Write();
                    hZTbkg->Write();
                    hZTRaw->Write();
                }
            }
            catch (std::exception &e)
            {
                std::cerr << "[WARN] SPlot(mass) failed for bin y=[" << yLo << "," << yHi << "]"
                          << (haveJetPt ? (std::string(" jetPt=[") + std::to_string(jLo) + "," + std::to_string(jHi) + "]") : std::string(" inclusive"))
                          << ": " << e.what() << std::endl;
            }
            if (fitRes)
            {
                delete fitRes;
                fitRes = nullptr;
            }

            // ------------------ IP chi2 fit (two Bukin components) ------------------
            if (doPlots && (haveLogIP || haveIP))
            {
                // Build dataset for IP chi2 using the same selection.
                // Prefer the existing log_tag_ipchi2 branch; otherwise fit raw tag_ip_chi2.
                RooArgSet ipObs;
                std::unique_ptr<RooRealVar> r_ip_raw;
                std::unique_ptr<RooRealVar> r_logip;

                if (haveLogIP)
                {
                    // Match range to Fitter.C: [-3, 5]
                    r_logip = std::make_unique<RooRealVar>("log_tag_ipchi2", "log_{10}(IP #chi^{2})", -3.0, 5.0);
                    ipObs.add(*r_logip);
                }
                else if (haveIP)
                {
                    r_ip_raw = std::make_unique<RooRealVar>("tag_ip_chi2", "IP #chi^{2}", 0.0, 1e4);
                    ipObs.add(*r_ip_raw);
                }
                // Include variables referenced in the cut expression
                ipObs.add(r_y);
                if (haveJetPt)
                    ipObs.add(r_jetPt);
                if (haveZ)
                    ipObs.add(r_zT);

                std::string ipdsName = mkname("data_ipchi2");
                RooDataSet ipds(ipdsName.c_str(), ipdsName.c_str(), tree, ipObs, cut.str().c_str());
                debugLog(std::string("Created IP chi2 dataset ") + ipdsName + ", entries: " + std::to_string(ipds.numEntries()));

                // Output naming and resume for IP
                std::ostringstream ipBase;
                ipBase << outDir << "/ipchi2fit_y" << iy << "_" << yLo << "to" << yHi;
                if (haveJetPt)
                    ipBase << "_j" << ij << "_" << jLo << "to" << jHi;
                std::string ipPng = ipBase.str() + ".png";
                // Choose the variable to fit and its frame
                RooRealVar *xVar = nullptr;
                if (haveLogIP)
                {
                    xVar = r_logip.get();
                }
                else if (haveIP)
                {
                    xVar = r_ip_raw.get();
                }

                if (xVar)
                {
                    // Define two Bukin components: prompt and nonprompt
                    // Prompt-like (match Fitter.C start values and limits)
                    RooRealVar xpP(("xpP" + suf).c_str(), "x_{p}^{P}", 0.0, -0.5, 1.0);
                    RooRealVar sigP(("sigP" + suf).c_str(), "#sigma_{P}", 0.7, 0.4, 2.0);
                    RooRealVar xiP(("xiP" + suf).c_str(), "#xi_{P}", 0.0, -0.5, 0.5);
                    RooRealVar rho1P(("rho1P" + suf).c_str(), "#rho_{1}^{P}", -0.08, -0.2, -0.01);
                    RooRealVar rho2P(("rho2P" + suf).c_str(), "#rho_{2}^{P}", 0.01, 0.0001, 0.98);

                    RooBukinPdf bukP(("bukP" + suf).c_str(), "Prompt Bukin", *xVar, xpP, sigP, xiP, rho1P, rho2P);

                    // Nonprompt-like (match Fitter.C start values and limits)
                    RooRealVar xpNP(("xpNP" + suf).c_str(), "x_{p}^{NP}", 1.9, 1.6, 3.0);
                    RooRealVar sigNP(("sigNP" + suf).c_str(), "#sigma_{NP}", 0.4, 0.3, 0.7);
                    RooRealVar xiNP(("xiNP" + suf).c_str(), "#xi_{NP}", 0.1, 0.05, 0.5);
                    RooRealVar rho1NP(("rho1NP" + suf).c_str(), "#rho_{1}^{NP}", -0.10, -0.95, -0.05);
                    RooRealVar rho2NP(("rho2NP" + suf).c_str(), "#rho_{2}^{NP}", 0.20, 0.01, 0.98);

                    RooBukinPdf bukNP(("bukNP" + suf).c_str(), "Nonprompt Bukin", *xVar, xpNP, sigNP, xiNP, rho1NP, rho2NP);

                    // Yields
                    double nTotIP = std::max(1, (int)ipds.numEntries());
                    RooRealVar nP(("nPrompt" + suf).c_str(), "N_{P}", 0.8 * nTotIP, 0.0, 2.0 * nTotIP);
                    RooRealVar nNP(("nNonPr" + suf).c_str(), "N_{NP}", 0.2 * nTotIP, 0.0, 2.0 * nTotIP);

                    RooAddPdf modelIP(("modelIP" + suf).c_str(), "IP Model", RooArgList(bukP, bukNP), RooArgList(nP, nNP));

                    // Fit
                    RooFitResult *fitIP = modelIP.fitTo(ipds, RooFit::Save(true), RooFit::Extended(true), RooFit::PrintLevel(-1), RooFit::PrintEvalErrors(-1));

                    // Plot: two-pad canvas with pulls
                    std::unique_ptr<TCanvas> cIP(new TCanvas(("cip_" + ipdsName).c_str(), ("IPchi2 fit " + dsTitle).c_str(), 800, 700));
                    cIP->Divide(1, 2, 0, 0);
                    TPad *pTopIP = (TPad *)cIP->cd(1);
                    pTopIP->SetPad(0, 0.28, 1, 1);
                    setupPadMargins(pTopIP, 0.12, 0.02, 0.01, 0.01);
                    pTopIP->SetLogy(true);
                    TPad *pBotIP = (TPad *)cIP->cd(2);
                    pBotIP->SetPad(0, 0, 1, 0.28);
                    setupPadMargins(pBotIP, 0.12, 0.35, 0.01, 0.05);

                    pTopIP->cd();
                    RooPlot *frIP = xVar->frame(RooFit::Bins(100));
                    std::string nmIPData = std::string("ipdata") + suf;
                    std::string nmIPMod = std::string("ipmodel") + suf;
                    std::string nmIPP = std::string("ipprompt") + suf;
                    std::string nmIPNP = std::string("ipnonprompt") + suf;
                    ipds.plotOn(frIP, RooFit::MarkerSize(0.8), RooFit::Name(nmIPData.c_str()));
                    modelIP.plotOn(frIP, RooFit::LineColor(kBlue + 1), RooFit::LineWidth(2), RooFit::Name(nmIPMod.c_str()));
                    modelIP.plotOn(frIP, RooFit::Components(bukP), RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen + 2), RooFit::Name(nmIPP.c_str()));
                    modelIP.plotOn(frIP, RooFit::Components(bukNP), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed + 1), RooFit::Name(nmIPNP.c_str()));
                    frIP->SetTitle("");
                    frIP->SetMinimum(0.5); // ensure positive min for log-y
                    frIP->GetXaxis()->SetLabelOffset(999);
                    frIP->GetXaxis()->SetTitleOffset(999);
                    frIP->Draw();

                    // Top-left bin label (rapidity and jet pT) for IP plot
                    {
                        TPaveText *ptIP = new TPaveText(0.15, 0.78, 0.48, 0.92, "NDC");
                        ptIP->SetFillColor(0);
                        ptIP->SetFillStyle(0);
                        ptIP->SetBorderSize(0);
                        ptIP->SetTextAlign(13); // left-top
                        ptIP->SetTextFont(42);
                        ptIP->SetTextSize(0.040);

                        std::ostringstream lblY;
                        lblY << std::fixed << std::setprecision(1) << yLo << " < #it{y} < " << yHi;
                        ptIP->AddText(lblY.str().c_str());
                        if (haveJetPt)
                        {
                            std::ostringstream lblJ;
                            lblJ << std::fixed << std::setprecision(0) << jLo << " < p_{T}^{jet} < " << jHi << " GeV";
                            ptIP->AddText(lblJ.str().c_str());
                        }
                        else
                        {
                            ptIP->AddText("jet pT: inclusive");
                        }
                        ptIP->Draw();
                    }

                    // Set appropriate y-axis range for log scale
                    double max_val = frIP->GetMaximum();
                    frIP->SetAxisRange(0.9, max_val * 5, "Y"); // Min value > 0 for log scale

                    // Legend
                    int nfloatIP = (fitIP ? fitIP->floatParsFinal().getSize() : 10);
                    double chi2ndfIP = frIP->chiSquare(nfloatIP);
                    TLegend *legIP = new TLegend(0.58, 0.55, 0.92, 0.93);
                    setupLegend(legIP, 0.035, 0.15, 42, 0, 0, 0);
                    TObject *oIPD = frIP->findObject(nmIPData.c_str());
                    TObject *oIPM = frIP->findObject(nmIPMod.c_str());
                    TObject *oIPP = frIP->findObject(nmIPP.c_str());
                    TObject *oIPNP = frIP->findObject(nmIPNP.c_str());
                    if (oIPD)
                        legIP->AddEntry(oIPD, "Data", "lep");
                    if (oIPM)
                        legIP->AddEntry(oIPM, "Total fit", "l");
                    if (oIPP)
                        legIP->AddEntry(oIPP, "Prompt (Bukin)", "l");
                    if (oIPNP)
                        legIP->AddEntry(oIPNP, "Nonprompt (Bukin)", "l");
                    {
                        std::ostringstream os;
                        os << "N_{P} = " << std::fixed << std::setprecision(1) << nP.getVal() << " +/- " << nP.getError();
                        legIP->AddEntry((TObject *)nullptr, os.str().c_str(), "");
                    }
                    {
                        std::ostringstream os;
                        os << "N_{NP} = " << std::fixed << std::setprecision(1) << nNP.getVal() << " +/- " << nNP.getError();
                        legIP->AddEntry((TObject *)nullptr, os.str().c_str(), "");
                    }
                    // Prompt fraction and uncertainty via error propagation
                    if ((nP.getVal() + nNP.getVal()) > 0)
                    {
                        RooFormulaVar fPrompt((std::string("fPrompt") + suf).c_str(), "@0/(@0+@1)", RooArgList(nP, nNP));
                        double fP = fPrompt.getVal();
                        double fPerr = (fitIP ? fPrompt.getPropagatedError(*fitIP) : 0.0);
                        std::ostringstream os;
                        os << "f_{prompt} = " << std::fixed << std::setprecision(3) << fP << " +/- " << fPerr;
                        legIP->AddEntry((TObject *)nullptr, os.str().c_str(), "");
                    }
                    {
                        std::ostringstream os;
                        os << "chi2/ndf = " << std::setprecision(3) << chi2ndfIP;
                        legIP->AddEntry((TObject *)nullptr, os.str().c_str(), "");
                    }
                    legIP->Draw();

                    // Pulls
                    pBotIP->cd();
                    RooHist *hIPpull = frIP->pullHist();
                    RooPlot *frPullIP = xVar->frame(RooFit::Bins(100));
                    frPullIP->addPlotable(hIPpull, "P");
                    frPullIP->getAttMarker()->SetMarkerSize(0.5);
                    frPullIP->getAttLine()->SetLineWidth(1);
                    frPullIP->SetTitle("");
                    frPullIP->GetYaxis()->SetTitle("Pull");
                    frPullIP->GetYaxis()->SetTitleSize(0.10);
                    frPullIP->GetYaxis()->SetTitleOffset(0.4);
                    frPullIP->GetYaxis()->SetLabelSize(0.09);
                    frPullIP->GetXaxis()->SetTitleSize(0.12);
                    frPullIP->GetXaxis()->SetLabelSize(0.10);
                    frPullIP->GetXaxis()->SetTitle(haveLogIP ? "log_{10}(IP #chi^{2})" : "IP #chi^{2}");
                    frPullIP->Draw();

                    // Save figure (PNG) always (overwrite if exists)
                    cIP->SaveAs(ipPng.c_str());

                    // Save IP distribution and fit curves to ROOT file
                    xVar->setBins(100);
                    TH1 *hIP = ipds.createHistogram(("hIP" + suf).c_str(), *xVar);
                    if (hIP)
                    {
                        fitFile->cd();
                        hIP->Write();
                    }
                    if (auto oTot = frIP->findObject(nmIPMod.c_str()))
                    {
                        fitFile->cd();
                        oTot->Write(("curveIPTot" + suf).c_str());
                    }
                    if (auto oP = frIP->findObject(nmIPP.c_str()))
                    {
                        fitFile->cd();
                        oP->Write(("curveIPPrompt" + suf).c_str());
                    }
                    if (auto oNP = frIP->findObject(nmIPNP.c_str()))
                    {
                        fitFile->cd();
                        oNP->Write(("curveIPNonPrompt" + suf).c_str());
                    }
                    fitFile->cd();
                    cIP->Write(("cIP" + suf).c_str());
                    // std::string ipPdf = onIP.str(); ipPdf.replace(ipPdf.end()-3, ipPdf.end(), "pdf"); cIP->SaveAs(ipPdf.c_str());

                    // -------- sPlot for IP fit: store per-event weights (prompt/nonprompt) --------
                    try
                    {
                        RooStats::SPlot sIP((std::string("sDataIP") + suf).c_str(), "sPlot IP", ipds, &modelIP, RooArgList(nP, nNP));
                        sWFile->cd();
                        std::string tnameIP = std::string("sWeights_ipchi2") + suf;
                        TTree tIP(tnameIP.c_str(), ("sWeights for IP " + dsTitle).c_str());
                        double outYip = 0.0, outJetIP = 0.0, outX = 0.0, wP = 0.0, wNP = 0.0;
                        tIP.Branch("tagY", &outYip, "tagY/D");
                        if (haveJetPt)
                            tIP.Branch("jetPt", &outJetIP, "jetPt/D");
                        tIP.Branch(xVar->GetName(), &outX, (std::string(xVar->GetName()) + "/D").c_str());
                        tIP.Branch("w_prompt", &wP, "w_prompt/D");
                        tIP.Branch("w_nonprompt", &wNP, "w_nonprompt/D");
                        const std::string wPName = std::string(nP.GetName()) + "_sw";
                        const std::string wNPName = std::string(nNP.GetName()) + "_sw";
                        for (int i = 0; i < ipds.numEntries(); ++i)
                        {
                            const RooArgSet *row = ipds.get(i);
                            outYip = row->getRealValue("tagY");
                            if (haveJetPt)
                                outJetIP = row->getRealValue("jetPt");
                            outX = row->getRealValue(xVar->GetName());
                            wP = row->getRealValue(wPName.c_str());
                            wNP = row->getRealValue(wNPName.c_str());
                            tIP.Fill();
                        }
                        tIP.Write();

                        // Also fill weighted zT histograms for IP prompt/nonprompt
                        if (haveZ)
                        {
                            std::unique_ptr<TH1D> hZTprompt(new TH1D(("hZT_prompt" + suf).c_str(), ("zT prompt " + dsTitle).c_str(), kZTBins, kZTMin, kZTMax));
                            std::unique_ptr<TH1D> hZTnonprompt(new TH1D(("hZT_nonprompt" + suf).c_str(), ("zT nonprompt " + dsTitle).c_str(), kZTBins, kZTMin, kZTMax));
                            for (int i = 0; i < ipds.numEntries(); ++i)
                            {
                                const RooArgSet *row = ipds.get(i);
                                double z = row->getRealValue("tagZ");
                                double wp = row->getRealValue(wPName.c_str());
                                double wnp = row->getRealValue(wNPName.c_str());
                                hZTprompt->Fill(z, wp);
                                hZTnonprompt->Fill(z, wnp);
                            }
                            fitFile->cd();
                            hZTprompt->Write();
                            hZTnonprompt->Write();

                            // Build combined signal*prompt weighted histogram and comparison canvas
                            TH1D *hRawSaved = (TH1D*)fitFile->Get(("hZT_raw" + suf).c_str());
                            TH1D *hSigSaved = (TH1D*)fitFile->Get(("hZT_sig" + suf).c_str());
                            std::unique_ptr<TH1D> hZTSigPrompt(new TH1D(("hZT_sigPrompt" + suf).c_str(), ("zT signal*prompt " + dsTitle).c_str(), kZTBins, kZTMin, kZTMax));
                            int nUse = std::min(ds.numEntries(), ipds.numEntries());
                            const std::string wSigName2 = std::string(nSig.GetName()) + "_sw";
                            for (int i = 0; i < nUse; ++i)
                            {
                                const RooArgSet *rowM = ds.get(i);
                                const RooArgSet *rowIP2 = ipds.get(i);
                                double z = rowM->getRealValue("tagZ");
                                double wS = rowM->getRealValue(wSigName2.c_str());
                                double wP2 = rowIP2->getRealValue(wPName.c_str());
                                hZTSigPrompt->Fill(z, wS * wP2);
                            }
                            fitFile->cd();
                            hZTSigPrompt->Write();

                            auto normClone = [](TH1D *h){ if(!h) return (TH1D*)nullptr; TH1D* c=(TH1D*)h->Clone((std::string(h->GetName())+"_norm").c_str()); double integ=c->Integral(); if(integ>0) c->Scale(1.0/integ); return c; };
                            TH1D *hRawN = normClone(hRawSaved);
                            TH1D *hSigN = normClone(hSigSaved);
                            TH1D *hPromptN = normClone(hZTprompt.get());
                            TH1D *hSigPromptN = normClone(hZTSigPrompt.get());
                            if (hRawN){ hRawN->SetLineColor(kBlack); hRawN->SetLineWidth(2); }
                            if (hSigN){ hSigN->SetLineColor(kBlue+1); hSigN->SetLineWidth(2); }
                            if (hPromptN){ hPromptN->SetLineColor(kGreen+2); hPromptN->SetLineWidth(2); }
                            if (hSigPromptN){ hSigPromptN->SetLineColor(kMagenta+1); hSigPromptN->SetLineWidth(2); }
                            std::ostringstream cmpName;
                            cmpName << outDir << "/zTcompare_y" << iy << "_" << yLo << "to" << yHi;
                            if (haveJetPt) cmpName << "_j" << ij << "_" << jLo << "to" << jHi;
                            cmpName << ".png";
                            std::unique_ptr<TCanvas> cCmp(new TCanvas(("cZTcompare" + suf).c_str(), ("zT comparison " + dsTitle).c_str(), 800, 600));
                            cCmp->cd();
                            if (hRawN){ hRawN->GetXaxis()->SetTitle("z_{T}"); hRawN->GetYaxis()->SetTitle("Normalized entries"); hRawN->Draw("hist"); }
                            if (hSigN) hSigN->Draw("hist same");
                            if (hPromptN) hPromptN->Draw("hist same");
                            if (hSigPromptN) hSigPromptN->Draw("hist same");
                            TLegend *legCmp = new TLegend(0.58,0.60,0.90,0.90);
                            legCmp->SetBorderSize(0); legCmp->SetFillStyle(0);
                            if (hRawN) legCmp->AddEntry(hRawN, "Raw", "l");
                            if (hSigN) legCmp->AddEntry(hSigN, "Signal", "l");
                            if (hPromptN) legCmp->AddEntry(hPromptN, "Prompt", "l");
                            if (hSigPromptN) legCmp->AddEntry(hSigPromptN, "Signal*Prompt", "l");
                            legCmp->Draw();
                            cCmp->SaveAs(cmpName.str().c_str());
                            fitFile->cd();
                            if (hRawN) hRawN->Write();
                            if (hSigN) hSigN->Write();
                            if (hPromptN) hPromptN->Write();
                            if (hSigPromptN) hSigPromptN->Write();
                            cCmp->Write(("cZTcompare" + suf).c_str());
                        }
                    }
                    catch (std::exception &e)
                    {
                        std::cerr << "[WARN] SPlot(IP) failed for bin y=[" << yLo << "," << yHi << "]"
                                  << (haveJetPt ? (std::string(" jetPt=[") + std::to_string(jLo) + "," + std::to_string(jHi) + "]") : std::string(" inclusive"))
                                  << ": " << e.what() << std::endl;
                    }

                    if (fitIP)
                    {
                        delete fitIP;
                        fitIP = nullptr;
                    }
                    // Progress status for IP fit: always increment
                    if (totalFitsPlanned)
                    {
                        ++fitsDone;
                        std::cout << "[Status] Fit " << fitsDone << "/" << totalFitsPlanned
                                  << " done: IPchi2 y=[" << yLo << "," << yHi << "]"
                                  << (haveJetPt ? (std::string(" jetPt=[") + std::to_string(jLo) + "," + std::to_string(jHi) + "]") : std::string(" inclusive"))
                                  << ". Remaining: " << (totalFitsPlanned - fitsDone) << std::endl;
                    }
                }
            }
        }
    }
    // Finalize ROOT file
    if (fitFile)
    {
        fitFile->Write(nullptr, TObject::kOverwrite);
        fitFile->Close();
    }
    if (sWFile)
    {
        sWFile->Write(nullptr, TObject::kOverwrite);
        sWFile->Close();
    }
    std::cout << "Finished. Plots written under: " << baseOutDir << " and ROOT saved to: " << fitsPath
              << "; sWeights saved to: " << sWPath << std::endl;
}