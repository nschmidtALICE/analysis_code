#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2.h>
#include <TTree.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLine.h>
#include <TColor.h>
#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TF1.h>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <cstdio>

// Map underflow/overflow into the boundary bins for 2D histograms (helper)
static void MapUnderOverflowToEdges(TH2 *h)
{
    if (!h)
        return;
    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();
    TH2 *tmp = (TH2 *)h->Clone(Form("%s_uomap_tmp", h->GetName()));
    tmp->Reset();
    for (int i = 0; i <= nx + 1; ++i)
    {
        for (int j = 0; j <= ny + 1; ++j)
        {
            double val = h->GetBinContent(i, j);
            double err = h->GetBinError(i, j);
            int im = (i < 1) ? 1 : ((i > nx) ? nx : i);
            int jm = (j < 1) ? 1 : ((j > ny) ? ny : j);
            double old = tmp->GetBinContent(im, jm);
            double oldErr = tmp->GetBinError(im, jm);
            tmp->SetBinContent(im, jm, old + val);
            tmp->SetBinError(im, jm, std::sqrt(oldErr * oldErr + err * err));
        }
    }
    for (int i = 1; i <= nx; ++i)
    {
        for (int j = 1; j <= ny; ++j)
        {
            h->SetBinContent(i, j, tmp->GetBinContent(i, j));
            h->SetBinError(i, j, tmp->GetBinError(i, j));
        }
    }
    for (int j = 0; j <= ny + 1; ++j)
    {
        h->SetBinContent(0, j, 0.0);
        h->SetBinError(0, j, 0.0);
        h->SetBinContent(nx + 1, j, 0.0);
        h->SetBinError(nx + 1, j, 0.0);
    }
    for (int i = 0; i <= nx + 1; ++i)
    {
        h->SetBinContent(i, 0, 0.0);
        h->SetBinError(i, 0, 0.0);
        h->SetBinContent(i, ny + 1, 0.0);
        h->SetBinError(i, ny + 1, 0.0);
    }
    delete tmp;
}

// Simple Chi2 calculator for 2D histograms
static double ComputeChi2_2D(const TH2 *hData, const TH2 *hModel, int &ndf)
{
    ndf = 0;
    if (!hData || !hModel)
        return 0.0;
    int nx = hData->GetNbinsX();
    int ny = hData->GetNbinsY();
    double chi2 = 0.0;
    for (int ix = 1; ix <= nx; ++ix)
    {
        for (int iy = 1; iy <= ny; ++iy)
        {
            double d = hData->GetBinContent(ix, iy);
            double e = hData->GetBinError(ix, iy);
            double m = hModel->GetBinContent(ix, iy);
            double err = e;
            // fallback: if no error stored, assume Poisson sqrt(d)
            if (err <= 0)
                err = (d > 0) ? std::sqrt(d) : 1.0;
            if (err > 0)
            {
                double r = (d - m) / err;
                chi2 += r * r;
                ++ndf;
            }
        }
    }
    return chi2;
}

static void chooseGrid(int n, int &ncols, int &nrows)
{
    ncols = (int)std::ceil(std::sqrt((double)n));
    nrows = (int)std::ceil((double)n / (double)ncols);
}

// Clone and normalize a TH1D to unit area (accounts for bin width). Returns a heap-allocated clone.
static TH1D *CloneNormalized(const TH1D *h, const char *newName)
{
    if (!h)
        return nullptr;
    TH1D *c = (TH1D *)h->Clone(newName);
    if (!c)
        return nullptr;
    c->SetDirectory(nullptr);
    double integral = c->Integral();
    if (integral > 0.0)
        c->Scale(1.0 / integral);
    return c;
}

// are now fixed to their previous default = true behavior to reduce complexity.
void unfold_2D(
    const std::string &outfile = "unfolded_output.root",
    int nIter = 4,
    const std::vector<std::string> &jetPtBins = {"5_10", "10_15", "15_20", "20_30", "30_50"},
    const std::vector<int> &yBins = {0, 1, 2, 3, 4, 5, 6, 7},
    bool isClosure = false,
    bool verbose = true)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    // select input file depending on closure flag
    const std::string &infileResponse = isClosure
                                            ? "/media/niviths/SSD2/lhcb_analysis_SSD/mc_merge_pPb_Pbp/20250728_pPb_MC_output_response.root"
                                            // : "/media/niviths/SSD2/lhcb_analysis_SSD/mc_merge_pPb_Pbp/response_merged.root";
                                            : "/media/niviths/SSD2/lhcb_analysis_SSD/mc_merge_MBonly/MBresponse.root";
    // (Color palette helpers moved to file scope.)
    // Naming scheme
    std::string unfoldedName = "unfolded_zT";
    std::vector<double> yBinBorders = {2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0};
    // Jet pT bin borders with underflow/overflow bins for response matrix
    std::vector<double> jetPtBinBorders = {2, 5, 10, 15, 20, 30, 50, 70};
    std::vector<double> zBins = {
        0.0, 0.05, 0.1, 0.15, 0.2,
        0.25, 0.3, 0.35, 0.4, 0.45,
        0.5, 0.55, 0.6, 0.65, 0.7,
        0.75, 0.8, 0.85, 0.9, 0.95,
        1.0};

    // Pattern for per-jet measured input files; %s will be replaced with jetPt string
    std::string infilePattern = isClosure ? "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_MC/TagZHistograms_%s.root" : "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2025-10-08_pPb/TagZHistograms_%s.root";
    // std::string infilePattern = isClosure ? "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_MC/TagZHistograms_%s.root" : "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA/TagZHistograms_%s.root";

    TFile *finResponse = TFile::Open(infileResponse.c_str());
    if (!finResponse || finResponse->IsZombie())
    {
        std::cerr << "Error: Cannot open response file " << infileResponse << std::endl;
        return;
    }
    else if (verbose)
    {
        std::cout << "Successfully opened response file: " << infileResponse << std::endl;
    }

    // Prepare output directory and file (dropped into dated folder)
    // current date string YYYY-MM-DD
    time_t t = time(nullptr);
    struct tm *lt = localtime(&t);
    char dateBuf[32] = {0};
    if (lt)
    {
        strftime(dateBuf, sizeof(dateBuf), "%Y-%m-%d", lt);
    }
    else
    {
        snprintf(dateBuf, sizeof(dateBuf), "unknown-date");
    }

    std::string outDir = isClosure ? Form("unfolded2D_zT_closure_%s", dateBuf) : Form("unfolded2D_zT_%s", dateBuf);
    // create directory (recursively) if needed
    std::string mkdirCmd = std::string("mkdir -p ") + outDir;
    int mkret = system(mkdirCmd.c_str());
    (void)mkret;

    std::string outPath = outDir + "/" + outfile;
    TFile *fout = TFile::Open(outPath.c_str(), "RECREATE");
    if (!fout || fout->IsZombie())
    {
        std::cerr << "Error: Cannot open output file " << outPath << std::endl;
        return;
    }
    if (verbose)
        std::cout << "Output directory: " << outDir << "\n";
    std::string treeName = "Response";
    // Get response tree
    TTree *tree = (TTree *)finResponse->Get(treeName.c_str());
    if (!tree)
    {
        std::cerr << "Warning: Missing response tree: " << treeName << std::endl;
        return;
    }

    // loop over yBins
    for (const auto &yBin : yBins)
    {
        // Compact: create a persistent 2D histogram with x = jetPt bin borders and y = z bins
        TH2D *h2_z_vs_jetpt = nullptr;
        if (jetPtBinBorders.size() >= 2 && zBins.size() >= 2)
        {
            h2_z_vs_jetpt = new TH2D(Form("h2_z_vs_jetpt_%d", yBin),
                                     "z_{T} vs Jet p_{T} bin;Jet p_{T} (GeV/c);z_{T}",
                                     (int)jetPtBinBorders.size() - 1, &jetPtBinBorders[0],
                                     (int)zBins.size() - 1, &zBins[0]);
            // Keep the histogram in memory for later filling/inspection; also write now so it exists in file
            if (verbose)
                std::cout << "Created persistent 2D histogram h2_z_vs_jetpt (x bins=" << (jetPtBinBorders.size() - 1) << ", y bins=" << (zBins.size() - 1) << ")\n";
        }
        for (const auto &jetPt : jetPtBins)
        {
            TString infileJet = Form(infilePattern.c_str(), jetPt.c_str());
            TFile *fin = TFile::Open(infileJet.Data());
            if (!fin || fin->IsZombie())
            {
                std::cerr << "Warning: Cannot open input file for jetPt " << jetPt << ": " << infileJet << std::endl;
                if (fin)
                {
                    fin->Close();
                    delete fin;
                }
                continue;
            }

            int bestIterIndexForThisY = 1; // will be updated after closure chi2 evaluation
            TString histName = Form("promptSignalTagZHist_FullyWeighted_%s_bin%d", jetPt.c_str(), yBin);

            if (verbose)
            {
                std::cout << "\n---\nProcessing jetPt bin: " << jetPt
                          << ", yBin: " << yBin << std::endl;
                std::cout << "Histogram name: " << histName << std::endl;
            }

            // Get measured histogram
            TH1D *hMeasured = (TH1D *)fin->Get(histName.Data());
            if (!hMeasured)
            {
                std::cerr << "Warning: Missing measured histogram: " << histName << std::endl;
                fin->Close();
                delete fin;
                continue;
            }
            if (verbose)
                std::cout << "Loaded measured histogram with " << hMeasured->GetEntries() << " entries." << std::endl;
            // Compact fill: deposit the measured 1D z-histogram into the persistent 2D map
            // x-coordinate: numeric jetPt center inferred from jetPt string (e.g. "5_10" -> 7.5)
            if (h2_z_vs_jetpt)
            {
                // parse jetPt string to numeric center
                double jetPtCenter = -1.0;
                size_t us = jetPt.find('_');
                if (us != std::string::npos)
                {
                    try
                    {
                        double low = std::stod(jetPt.substr(0, us));
                        double high = std::stod(jetPt.substr(us + 1));
                        jetPtCenter = 0.5 * (low + high);
                    }
                    catch (...)
                    {
                        jetPtCenter = -1.0;
                    }
                }
                // fallback: if parsing failed, try to place in middle of available jetPt bin borders
                if (jetPtCenter < 0 && jetPtBinBorders.size() >= 2)
                {
                    // use center of first/last as rough fallback
                    jetPtCenter = 0.5 * (jetPtBinBorders.front() + jetPtBinBorders.back());
                }

                // Set the 2D bin content/error directly per z-bin for this jet pT slice
                int nBinsZ = hMeasured->GetNbinsX();
                int ix = h2_z_vs_jetpt->GetXaxis()->FindBin(jetPtCenter);
                int nx = h2_z_vs_jetpt->GetNbinsX();
                int ny = h2_z_vs_jetpt->GetNbinsY();
                if (ix < 1)
                    ix = 1;
                if (ix > nx)
                    ix = nx;

                if (verbose)
                {
                    std::cout << "Mapping jetPt " << jetPt << " (center=" << jetPtCenter << ") to x-bin " << ix
                              << " [" << h2_z_vs_jetpt->GetXaxis()->GetBinLowEdge(ix) << "-"
                              << h2_z_vs_jetpt->GetXaxis()->GetBinUpEdge(ix) << "]" << std::endl;
                    std::cout << "  Measured histogram integral: " << hMeasured->Integral() << std::endl;
                }

                for (int ib = 1; ib <= nBinsZ; ++ib)
                {
                    double zVal = hMeasured->GetBinCenter(ib);
                    int iy = h2_z_vs_jetpt->GetYaxis()->FindBin(zVal);
                    if (iy < 1 || iy > ny)
                        continue;
                    double w = hMeasured->GetBinContent(ib);
                    double e = hMeasured->GetBinError(ib);

                    // Check for existing content and accumulate if needed
                    double existingW = h2_z_vs_jetpt->GetBinContent(ix, iy);
                    double existingE = h2_z_vs_jetpt->GetBinError(ix, iy);
                    if (existingW > 0 && verbose)
                    {
                        std::cout << "Warning: Overwriting bin (" << ix << "," << iy << ") content " << existingW << " with " << w << std::endl;
                    }

                    h2_z_vs_jetpt->SetBinContent(ix, iy, w);
                    h2_z_vs_jetpt->SetBinError(ix, iy, e);
                }
            }
            // Close per-jet input after use to avoid too many open files
            fin->Close();
            delete fin;
        }

        if (h2_z_vs_jetpt)
        {
            if (verbose)
            {
                std::cout << "=== MEASURED HISTOGRAM SUMMARY ===" << std::endl;
                std::cout << "Total measured h2_z_vs_jetpt integral: " << h2_z_vs_jetpt->Integral() << std::endl;
                std::cout << "Measured histogram dimensions: " << h2_z_vs_jetpt->GetNbinsX() << " x " << h2_z_vs_jetpt->GetNbinsY() << std::endl;
                for (int ix = 1; ix <= h2_z_vs_jetpt->GetNbinsX(); ++ix)
                {
                    double binIntegral = 0;
                    for (int iy = 1; iy <= h2_z_vs_jetpt->GetNbinsY(); ++iy)
                    {
                        binIntegral += h2_z_vs_jetpt->GetBinContent(ix, iy);
                    }
                    std::cout << "  x-bin " << ix << " [" << h2_z_vs_jetpt->GetXaxis()->GetBinLowEdge(ix)
                              << "-" << h2_z_vs_jetpt->GetXaxis()->GetBinUpEdge(ix) << "]: " << binIntegral << std::endl;
                }
            }
            // Build 2D response for this y-bin from the response tree
            // Templates inherit the exact binning from the measured map
            TH2D *hMeasTmpl = (TH2D *)h2_z_vs_jetpt->Clone(Form("meas2Dtmpl_%d", yBin));
            hMeasTmpl->Reset();
            hMeasTmpl->SetTitle("Reco (measured) template 2D");
            TH2D *hTruthTmpl = (TH2D *)h2_z_vs_jetpt->Clone(Form("truth2Dtmpl_%d", yBin));
            hTruthTmpl->Reset();
            hTruthTmpl->SetTitle("Truth template 2D");

            RooUnfoldResponse response2D(hMeasTmpl, hTruthTmpl);

            // Variables for response tree branches
            float d0_z_det = 0, d0_z_mc = 0;
            float jet_pt_det = 0, jet_pt_mc = 0;
            float d0_y_det = 0, d0_y_mc = 0; // These might be d0_eta_det/d0_eta_mc
            float jet_nconst_det = 0, jet_nconst_mc = 0;
            float jet_dr = 0;
            float evtWeight = 1.0f;

            tree->SetBranchAddress("d0_z_det", &d0_z_det);
            tree->SetBranchAddress("d0_z_mc", &d0_z_mc);
            tree->SetBranchAddress("jet_pt_det", &jet_pt_det);
            tree->SetBranchAddress("jet_pt_mc", &jet_pt_mc);

            // Try both possible branch names for eta/rapidity
            if (tree->GetBranch("d0_eta_det"))
            {
                tree->SetBranchAddress("d0_eta_det", &d0_y_det);
                tree->SetBranchAddress("d0_eta_mc", &d0_y_mc);
                if (verbose)
                    std::cout << "Using d0_eta_det/d0_eta_mc branches" << std::endl;
            }
            else if (tree->GetBranch("d0_y_det"))
            {
                tree->SetBranchAddress("d0_y_det", &d0_y_det);
                tree->SetBranchAddress("d0_y_mc", &d0_y_mc);
                if (verbose)
                    std::cout << "Using d0_y_det/d0_y_mc branches" << std::endl;
            }
            else
            {
                std::cerr << "Warning: No d0_eta_det or d0_y_det branch found!" << std::endl;
            }

            tree->SetBranchAddress("jet_nconst_det", &jet_nconst_det);
            tree->SetBranchAddress("jet_nconst_mc", &jet_nconst_mc);
            tree->SetBranchAddress("jet_dr", &jet_dr);
            if (tree->GetBranch("weight"))
                tree->SetBranchAddress("weight", &evtWeight);

            // Determine eta bin borders from yBin index
            int yBinIdx = yBin;
            if (yBinIdx < 0 || yBinIdx + 1 >= (int)yBinBorders.size())
            {
                std::cerr << "Warning: yBin index out of range for borders (" << yBinIdx << ")\n";
                // Still write the measured map for visibility
                fout->cd();
                h2_z_vs_jetpt->SetDirectory(fout);
                h2_z_vs_jetpt->Write();
                delete hMeasTmpl;
                delete hTruthTmpl;
                continue;
            }

            // Define eta range for this yBin - make available for entire yBin processing
            double eta_min = yBinBorders[yBinIdx];
            double eta_max = yBinBorders[yBinIdx + 1];

            // === CRITICAL: SCALE MATCHING ===
            // The response matrix must be built using the same statistical weight as the measured data
            double measuredTotal = h2_z_vs_jetpt->Integral();
            if (verbose)
                std::cout << "Measured 2D histogram integral: " << measuredTotal << std::endl;

            Long64_t nEntries = tree->GetEntries();
            int nFilled = 0;
            for (Long64_t i = 0; i < nEntries; ++i)
            {
                tree->GetEntry(i);
                double w = (double)evtWeight;
                bool passDetEta = (d0_y_det >= eta_min && d0_y_det < eta_max);
                bool passMcEta = (d0_y_mc >= eta_min && d0_y_mc < eta_max);

                // Fill truth template from MC data (like unfold_new does)
                if (passMcEta)
                {
                    hTruthTmpl->Fill(jet_pt_mc, d0_z_mc, w);
                }

                if (passDetEta && passMcEta)
                {
                    response2D.Fill(jet_pt_det, d0_z_det, jet_pt_mc, d0_z_mc, w);
                    nFilled++;
                }
                else if (passMcEta && !passDetEta)
                {
                    response2D.Miss(jet_pt_mc, d0_z_mc, w);
                }
                else if (!passMcEta && passDetEta)
                {
                    response2D.Fake(jet_pt_det, d0_z_det, w);
                }
            }

            if (verbose)
            {
                std::cout << "Response building: " << nFilled << " Fill events from " << nEntries << " total entries" << std::endl;
                std::cout << "Truth template integral: " << hTruthTmpl->Integral() << std::endl;
            }

            // Validate response matrix has sufficient statistics (like unfold_new)
            if (nFilled < 50)
            {
                std::cerr << "Warning: Insufficient response statistics (" << nFilled << " events) for yBin " << yBin << std::endl;
                delete hMeasTmpl;
                delete hTruthTmpl;
                continue;
            }

            // Map under/overflow into boundary bins (CRITICAL - missing in original)
            MapUnderOverflowToEdges(h2_z_vs_jetpt); // Apply to measured histogram
            MapUnderOverflowToEdges(hTruthTmpl);    // Apply to truth template

            // // Optionally map under/overflow of the response matrix to the edge bins
            // if (TH2 *hRespInit = response2D.Hresponse())
            // {
            //     MapUnderOverflowToEdges(hRespInit);
            // }

            // Diagnostic: check for x-bins with no measured data
            // (Note: we only print; we do not modify the response here.)
            if (TH2 *hResp = response2D.Hresponse())
            {
                int nRespBins = hResp->GetNbinsX() * hResp->GetNbinsY();
                for (int ix = 1; ix <= h2_z_vs_jetpt->GetNbinsX(); ++ix)
                {
                    bool hasMeasuredData = false;
                    for (int iy = 1; iy <= h2_z_vs_jetpt->GetNbinsY(); ++iy)
                    {
                        if (h2_z_vs_jetpt->GetBinContent(ix, iy) > 0)
                        {
                            hasMeasuredData = true;
                            break;
                        }
                    }
                    // If no measured data in this jet pT slice, zero corresponding response elements
                    if (!hasMeasuredData && verbose)
                    {
                        std::cout << "No measured data in x-bin " << ix << " ["
                                  << h2_z_vs_jetpt->GetXaxis()->GetBinLowEdge(ix) << "-"
                                  << h2_z_vs_jetpt->GetXaxis()->GetBinUpEdge(ix) << "]" << std::endl;
                    }
                }
            }

            // Debug output to check data consistency
            if (verbose)
            {
                std::cout << "=== DEBUG INFO for yBin " << yBin << " ===" << std::endl;
                std::cout << "Measured h2_z_vs_jetpt integral: " << h2_z_vs_jetpt->Integral() << std::endl;
                if (TH2 *hResp = response2D.Hresponse())
                {
                    std::cout << "Response matrix integral: " << hResp->Integral() << std::endl;
                    std::cout << "Response matrix dimensions: " << hResp->GetNbinsX() << " x " << hResp->GetNbinsY() << std::endl;
                }
                if (TH2 *hTruth = dynamic_cast<TH2 *>(response2D.Htruth()))
                {
                    std::cout << "Truth histogram integral: " << hTruth->Integral() << std::endl;
                    std::cout << "Truth dimensions: " << hTruth->GetNbinsX() << " x " << hTruth->GetNbinsY() << std::endl;
                }
                if (TH2 *hMeas = dynamic_cast<TH2 *>(response2D.Hmeasured()))
                {
                    std::cout << "Response measured histogram integral: " << hMeas->Integral() << std::endl;
                    std::cout << "Response measured dimensions: " << hMeas->GetNbinsX() << " x " << hMeas->GetNbinsY() << std::endl;
                }
            }

            // Plot initial measured and truth 2D distributions side-by-side for this y-bin
            {
                TH2 *truth2D = dynamic_cast<TH2 *>(response2D.Htruth());
                int ncols = truth2D ? 2 : 1;
                int nrows = 1;
                TCanvas *cInit = new TCanvas(Form("c_initial2D_meas_truth_y%d", yBin), "Initial 2D: Measured and Truth", 640 * ncols, 540 * nrows);
                cInit->Divide(ncols, nrows);
                // Harmonize color scale if both exist
                double zmax = h2_z_vs_jetpt->GetMaximum();
                if (truth2D)
                    zmax = std::max(zmax, truth2D->GetMaximum());
                // Panel 1: measured
                cInit->cd(1);
                gPad->SetRightMargin(0.12);
                gPad->SetLeftMargin(0.12);
                gPad->SetLogz(1);
                gPad->SetTickx();
                gPad->SetTicky();
                h2_z_vs_jetpt->GetZaxis()->SetRangeUser(0.0, zmax * 1.05);
                h2_z_vs_jetpt->SetTitle(Form("Measured %.1f<#it{y}<%.1f;#it{p}_{T}^{jet} [GeV];z_{T}", eta_min, eta_max));
                h2_z_vs_jetpt->Draw("COLZ");
                // Panel 2: truth (if available)
                if (truth2D)
                {
                    double zmax2 = truth2D->GetMaximum();
                    cInit->cd(2);
                    gPad->SetRightMargin(0.12);
                    gPad->SetLeftMargin(0.12);
                    gPad->SetLogz(1);
                    gPad->SetTickx();
                    gPad->SetTicky();
                    truth2D->GetZaxis()->SetRangeUser(0.0, zmax2 * 1.05);
                    truth2D->SetTitle(Form("Truth (MC) %.1f<#it{y}<%.1f;#it{p}_{T}^{jet} [GeV];z_{T}", eta_min, eta_max));
                    truth2D->Draw("COLZ");
                }
                cInit->SaveAs(Form("%s/initial2D_meas_truth_y%d.png", outDir.c_str(), yBin));
                delete cInit;
            }

            // === SCALE MATCHING (like unfold_new.cpp) ===
            // Ensure measured histogram and response measured are consistently scaled
            TH2D *h2_z_vs_jetpt_forUnfold = h2_z_vs_jetpt;
            TH2D *h2_z_vs_jetpt_forUnfoldOwned = nullptr;
            {
                TH1 *hRespMeasNow = response2D.Hmeasured();
                if (hRespMeasNow && measuredTotal > 0)
                {
                    double respMeasIntegral = hRespMeasNow->Integral();
                    if (verbose)
                    {
                        std::cout << "Scale matching: measured=" << measuredTotal << ", response measured=" << respMeasIntegral << std::endl;
                    }
                    if (std::abs(measuredTotal - respMeasIntegral) > 1e-6)
                    {
                        // Create scaled copy for unfolding
                        double scale = respMeasIntegral / measuredTotal;
                        h2_z_vs_jetpt_forUnfoldOwned = (TH2D *)h2_z_vs_jetpt->Clone(Form("h2_scaled_y%d", yBin));
                        h2_z_vs_jetpt_forUnfoldOwned->Scale(scale);
                        h2_z_vs_jetpt_forUnfold = h2_z_vs_jetpt_forUnfoldOwned;
                        if (verbose)
                        {
                            std::cout << "Applied scale factor: " << scale << std::endl;
                        }
                    }
                }
            }

            // Iteration scan: unfold at iter=1..nIter, compute chi2(measured, refolded)
            std::vector<TH2D *> unfoldedPerIter;
            std::vector<TH2D *> refoldedPerIter;
            std::vector<double> redChi2PerIter; // aligns with refoldedPerIter
            unfoldedPerIter.reserve(nIter);
            refoldedPerIter.reserve(nIter);
            redChi2PerIter.reserve(nIter);

            int bestIter = 1;
            double bestRedChi2 = 1e99;
            for (int iter = 1; iter <= nIter; ++iter)
            {
                RooUnfoldBayes unfoldIter(&response2D, h2_z_vs_jetpt_forUnfold, iter);
                TH2D *hReco = (TH2D *)unfoldIter.Hreco(RooUnfold::kCovToy);
                if (!hReco)
                    continue;
                hReco->SetName(Form("unfolded_zT2D_bin%d_iter%d", yBin, iter));
                unfoldedPerIter.push_back((TH2D *)hReco->Clone());

                // Refold: apply response to unfolded truth estimate
                TH1 *hRefoldBase = response2D.ApplyToTruth(hReco);  // returns TH1* (base for TH2D)
                TH2D *hRefoldD = dynamic_cast<TH2D *>(hRefoldBase); // Direct cast, no clone
                if (hRefoldD)
                {
                    hRefoldD->SetName(Form("refolded2D_bin%d_iter%d", yBin, iter));
                    refoldedPerIter.push_back((TH2D *)hRefoldD->Clone()); // Clone after name set

                    // Debug output
                    if (verbose && iter == 1)
                    {
                        std::cout << "Iter " << iter << " - Measured integral: " << h2_z_vs_jetpt->Integral()
                                  << ", Unfolded integral: " << hReco->Integral()
                                  << ", Refolded integral: " << hRefoldD->Integral() << std::endl;
                    }

                    int ndf = 0;
                    // Compare against the histogram actually used for unfolding (may be scaled)
                    double chi2 = ComputeChi2_2D(h2_z_vs_jetpt_forUnfold, hRefoldD, ndf);
                    double red = (ndf > 0) ? (chi2 / ndf) : 1e99;
                    redChi2PerIter.push_back(red);
                    if (verbose)
                        std::cout << "yBin " << yBin << ": iter=" << iter << ", chi2/ndf=" << red << std::endl;
                    if (red < bestRedChi2)
                    {
                        bestRedChi2 = red;
                    }
                    bestIter = iter;
                }
                if (hRefoldBase)
                    delete hRefoldBase;
            }

            // Write measured and response objects
            fout->cd();
            // Save both original measured and the histogram used for unfolding (if scaled)
            h2_z_vs_jetpt->SetDirectory(fout);
            h2_z_vs_jetpt->Write();
            if (h2_z_vs_jetpt_forUnfold != h2_z_vs_jetpt && h2_z_vs_jetpt_forUnfold)
            {
                TH2D *h2_unf_copy = (TH2D *)h2_z_vs_jetpt_forUnfold->Clone(Form("measured2D_forUnfold_y%d", yBin));
                h2_unf_copy->SetDirectory(fout);
                h2_unf_copy->Write();
                delete h2_unf_copy;
            }
            response2D.Write(Form("response2D_bin%d", yBin));

            // Write all unfolded/refolded iterations and pick best
            for (auto *h : unfoldedPerIter)
                if (h)
                    h->Write();
            for (auto *h : refoldedPerIter)
                if (h)
                    h->Write();

            // Multi-panel: Unfolded 2D per iteration
            if (!unfoldedPerIter.empty())
            {
                int ncols = 1, nrows = 1;
                chooseGrid((int)unfoldedPerIter.size(), ncols, nrows);
                TCanvas *cUnf = new TCanvas(Form("c_unfold2D_iters_y%d", yBin), "Unfolded 2D (per iter)", 600 * ncols, 500 * nrows);
                cUnf->Divide(ncols, nrows);
                for (size_t i = 0; i < unfoldedPerIter.size(); ++i)
                {
                    cUnf->cd((int)i + 1);
                    gPad->SetRightMargin(0.12);
                    gPad->SetLeftMargin(0.08);
                    gPad->SetBottomMargin(0.08);
                    gPad->SetTopMargin(0.01);
                    gPad->SetTickx();
                    gPad->SetTicky();
                    if (unfoldedPerIter[i])
                    {
                        unfoldedPerIter[i]->SetTitle("");
                        gPad->SetLogz(1);
                        unfoldedPerIter[i]->Draw("COLZ");
                        TLatex latex;
                        latex.SetNDC();
                        latex.SetTextSize(0.045);
                        latex.DrawLatex(0.15, 0.90, Form("Iteration %zu", i + 1));
                        if (i == 0)
                        {
                            latex.DrawLatex(0.15, 0.83, Form("%.1f<#it{y}<%.1f", eta_min, eta_max));
                            latex.DrawLatex(0.15, 0.76, Form("Unfolded distribution"));
                        }
                    }
                }
                cUnf->SaveAs(Form("%s/unfold2D_iters_y%d.png", outDir.c_str(), yBin));
            }

            // Multi-panel: Refolded 2D per iteration
            if (!refoldedPerIter.empty())
            {
                int ncols = 1, nrows = 1;
                chooseGrid((int)refoldedPerIter.size(), ncols, nrows);
                TCanvas *cRef = new TCanvas(Form("c_refold2D_iters_y%d", yBin), "Refolded 2D (per iter)", 600 * ncols, 500 * nrows);
                cRef->Divide(ncols, nrows);
                for (size_t i = 0; i < refoldedPerIter.size(); ++i)
                {
                    cRef->cd((int)i + 1);
                    gPad->SetRightMargin(0.12);
                    gPad->SetLeftMargin(0.08);
                    gPad->SetBottomMargin(0.08);
                    gPad->SetTopMargin(0.01);
                    gPad->SetTickx();
                    gPad->SetTicky();
                    gPad->SetLogz(1);
                    if (refoldedPerIter[i])
                    {
                        refoldedPerIter[i]->SetTitle("");
                        refoldedPerIter[i]->Draw("COLZ");
                        TLatex latex;
                        latex.SetNDC();
                        latex.SetTextSize(0.045);
                        latex.DrawLatex(0.15, 0.90, Form("Iteration %zu", i + 1));
                        if (i == 0)
                            latex.DrawLatex(0.15, 0.83, Form("%.1f<#it{y}<%.1f", eta_min, eta_max));
                        if (i == 0)
                            latex.DrawLatex(0.15, 0.76, Form("Refolded distribution"));
                    }
                }
                cRef->SaveAs(Form("%s/refold2D_iters_y%d.png", outDir.c_str(), yBin));
            }

            // Multi-panel: Ratio (refolded / measured) per iteration
            if (!refoldedPerIter.empty())
            {
                std::vector<TH2D *> ratioPerIter;
                ratioPerIter.reserve(refoldedPerIter.size());
                for (size_t i = 0; i < refoldedPerIter.size(); ++i)
                {
                    if (!refoldedPerIter[i])
                    {
                        ratioPerIter.push_back(nullptr);
                        continue;
                    }
                    TH2D *hRatio = (TH2D *)refoldedPerIter[i]->Clone(Form("ratio_refold_over_meas_y%d_iter%zu", yBin, i + 1));
                    // Divide refolded by measured (binomial errors)
                    // Use the histogram fed to unfolding for a consistent normalization in the ratio
                    hRatio->Divide(refoldedPerIter[i], h2_z_vs_jetpt_forUnfold ? h2_z_vs_jetpt_forUnfold : h2_z_vs_jetpt, 1.0, 1.0, "B");
                    ratioPerIter.push_back(hRatio);
                    // write ratio histogram
                    hRatio->Write();
                }
                int ncols = 1, nrows = 1;
                chooseGrid((int)ratioPerIter.size(), ncols, nrows);
                TCanvas *cRatio = new TCanvas(Form("c_refold2D_ratio_iters_y%d", yBin), "Refold/Measured ratio (per iter)", 600 * ncols, 500 * nrows);
                cRatio->Divide(ncols, nrows);
                for (size_t i = 0; i < ratioPerIter.size(); ++i)
                {
                    cRatio->cd((int)i + 1);
                    gPad->SetRightMargin(0.12);
                    gPad->SetLeftMargin(0.08);
                    gPad->SetBottomMargin(0.08);
                    gPad->SetTopMargin(0.01);
                    gPad->SetTickx();
                    gPad->SetTicky();
                    if (ratioPerIter[i])
                    {
                        ratioPerIter[i]->SetTitle("");
                        ratioPerIter[i]->GetZaxis()->SetRangeUser(0.0, 2.0);
                        ratioPerIter[i]->GetZaxis()->SetTitle("Refold / Measured");
                        ratioPerIter[i]->Draw("COLZ");
                        TLatex latex;
                        latex.SetNDC();
                        latex.SetTextSize(0.045);
                        latex.DrawLatex(0.15, 0.90, Form("Iteration %zu", i + 1));
                        if (i == 0)
                        {
                            latex.DrawLatex(0.15, 0.83, Form("%.1f<#it{y}<%.1f", eta_min, eta_max));
                            latex.DrawLatex(0.15, 0.76, Form("Refolded / Measured"));
                        }
                    }
                }
                cRatio->SaveAs(Form("%s/refold2D_ratio_iters_y%d.png", outDir.c_str(), yBin));
                // cleanup ratio clones
                for (auto *h : ratioPerIter)
                    delete h;
            }

            // Save best iteration as a named clone for convenience
            TH2D *bestUnf = nullptr;
            if (bestIter >= 1 && bestIter <= (int)unfoldedPerIter.size())
            {
                bestUnf = (TH2D *)unfoldedPerIter[bestIter - 1]->Clone(Form("unfolded2D_best_y%d_iter%d", yBin, bestIter));
                if (bestUnf)
                    bestUnf->Write();
            }
            // Also save the best refold/measurement ratio for quick inspection
            if (bestIter >= 1 && bestIter <= (int)refoldedPerIter.size() && refoldedPerIter[bestIter - 1])
            {
                TH2D *bestRatio = (TH2D *)refoldedPerIter[bestIter - 1]->Clone(Form("ratio_refold_over_meas_best_y%d_iter%d", yBin, bestIter));
                bestRatio->Divide(refoldedPerIter[bestIter - 1], h2_z_vs_jetpt_forUnfold ? h2_z_vs_jetpt_forUnfold : h2_z_vs_jetpt, 1.0, 1.0, "B");
                bestRatio->Write();
                delete bestRatio;
            }

            // =========================
            // 1D projections and plots
            // =========================
            // Project along zT (Y) integrated over all jet pT bins
            int nxBins = h2_z_vs_jetpt->GetNbinsX();
            // Project the histogram used for unfolding (scaled) for consistent refold comparisons
            TH2D *h2_meas_for_proj = h2_z_vs_jetpt_forUnfold ? h2_z_vs_jetpt_forUnfold : h2_z_vs_jetpt;
            TH1D *measProjAll = h2_meas_for_proj->ProjectionY(Form("measured1D_z_allpt_y%d", yBin), 1, nxBins);
            if (measProjAll)
            {
                measProjAll->SetTitle(Form("Measured z_{T} (all p_{T}) yBin %d", yBin));
                measProjAll->GetXaxis()->SetTitle("z_{T}");
                measProjAll->GetYaxis()->SetTitle("Counts");
                measProjAll->Write();
            }

            // Multi-panel: per-iteration overlay of measured vs refolded 1D z projection
            if (!refoldedPerIter.empty())
            {
                int ncols = 1, nrows = 1;
                chooseGrid((int)refoldedPerIter.size(), ncols, nrows);
                TCanvas *cProjMR = new TCanvas(Form("c_projY_meas_refold_iters_y%d", yBin), "Measured vs Refolded (1D z, per iter)", 700 * ncols, 560 * nrows);
                cProjMR->Divide(ncols, nrows);
                std::vector<TH1D *> refProjToDelete;
                std::vector<TH1D *> normClonesToDelete;
                refProjToDelete.reserve(refoldedPerIter.size());
                normClonesToDelete.reserve(2 * refoldedPerIter.size());
                for (size_t i = 0; i < refoldedPerIter.size(); ++i)
                {
                    cProjMR->cd((int)i + 1);
                    gPad->SetRightMargin(0.01);
                    gPad->SetLeftMargin(0.1);
                    gPad->SetTopMargin(0.01);
                    gPad->SetBottomMargin(0.08);
                    gPad->SetTickx();
                    gPad->SetTicky();
                    gPad->SetTickx();
                    gPad->SetTicky();
                    TH1D *refProj = refoldedPerIter[i] ? refoldedPerIter[i]->ProjectionY(Form("refold1D_z_allpt_y%d_iter%zu", yBin, i + 1), 1, nxBins) : nullptr;
                    refProjToDelete.push_back(refProj);
                    // Draw normalized clones for comparison; keep originals written unmodified
                    TH1D *measNorm = CloneNormalized(measProjAll, Form("measured1D_z_allpt_y%d_norm_iter%zu", yBin, i + 1));
                    TH1D *refNorm = CloneNormalized(refProj, Form("refold1D_z_allpt_y%d_iter%zu_norm", yBin, i + 1));
                    if (measNorm)
                    {
                        measNorm->SetLineColor(kBlack);
                        measNorm->SetMarkerColor(kBlack);
                        measNorm->SetMarkerStyle(20);
                        measNorm->GetXaxis()->SetTitle("z_{T}");
                        measNorm->GetYaxis()->SetTitle("1/N dN/dz_{T}");
                        measNorm->GetYaxis()->SetRangeUser(0.0, 1.4 * measNorm->GetMaximum());
                        measNorm->SetTitle("");
                        measNorm->Draw("E1");
                        normClonesToDelete.push_back(measNorm);
                    }
                    if (refNorm)
                    {
                        refNorm->SetLineColor(kBlue + 1);
                        refNorm->SetLineWidth(2);
                        refNorm->Draw("HIST SAME");
                        normClonesToDelete.push_back(refNorm);
                    }
                    // Legend and labels per pad
                    TLegend *leg = new TLegend(0.55, 0.72, 0.88, 0.88);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    if (measNorm)
                        leg->AddEntry(measNorm, "Measured (norm)", "lep");
                    if (refNorm)
                        leg->AddEntry(refNorm, Form("Refolded (iter %zu, norm)", i + 1), "l");
                    leg->Draw();
                    TLatex latex;
                    latex.SetNDC();
                    latex.SetTextSize(0.045);
                    if (i < redChi2PerIter.size())
                        latex.DrawLatex(0.15, 0.90, Form("iter %zu, #chi^{2}/ndf=%.3g", i + 1, redChi2PerIter[i]));
                    latex.DrawLatex(0.15, 0.83, Form("%.1f<#eta<%.1f", eta_min, eta_max));
                }
                cProjMR->SaveAs(Form("%s/projY_meas_vs_refold_iters_y%d.png", outDir.c_str(), yBin));
                for (auto *h : refProjToDelete)
                    if (h)
                        delete h;
                for (auto *h : normClonesToDelete)
                    if (h)
                        delete h;
                delete cProjMR;
            }

            // Per-jet pT slice projections (best iteration only)
            if (bestUnf && bestIter >= 1 && bestIter <= (int)refoldedPerIter.size() && refoldedPerIter[bestIter - 1])
            {
                TH2D *bestRefold2D = refoldedPerIter[bestIter - 1];

                // Skip first and last bins in projY plots
                int nPlotBins = nxBins - 2; // Exclude first and last
                if (nPlotBins < 1)
                    nPlotBins = 1; // Safety check
                int ncols = 1, nrows = 1;
                chooseGrid(nPlotBins, ncols, nrows);

                // Overlay measured vs refolded for each x-slice (normalized for comparison)
                TCanvas *cSliceMR = new TCanvas(Form("c_projY_meas_vs_refold_best_slices_y%d", yBin), "Measured vs Refolded (best) per jet p_{T} bin", 1400 * ncols, 1000 * nrows);
                cSliceMR->Divide(ncols, nrows);
                std::vector<TH1D *> mSlices, rSlices;
                std::vector<TH1D *> normClonesToDelete2;
                mSlices.reserve(nPlotBins);
                rSlices.reserve(nPlotBins);
                normClonesToDelete2.reserve(2 * nPlotBins);
                int padIndex = 1;
                for (int ix = 2; ix <= nxBins - 1; ++ix) // Skip first (ix=1) and last (ix=nxBins)
                {
                    cSliceMR->cd(padIndex);
                    gPad->SetRightMargin(0.01);
                    gPad->SetLeftMargin(0.10);
                    gPad->SetTopMargin(0.01);
                    gPad->SetBottomMargin(0.08);
                    gPad->SetTickx();
                    gPad->SetTicky();
                    gPad->SetTickx();
                    gPad->SetTicky();
                    TH1D *mSlice = h2_z_vs_jetpt->ProjectionY(Form("meas1D_z_slice_x%d_y%d", ix, yBin), ix, ix);
                    TH1D *rSlice = bestRefold2D->ProjectionY(Form("refold1D_z_slice_x%d_y%d", ix, yBin), ix, ix);
                    mSlices.push_back(mSlice);
                    rSlices.push_back(rSlice);
                    TH1D *mNorm = CloneNormalized(mSlice, Form("meas1D_z_slice_x%d_y%d_normMR", ix, yBin));
                    TH1D *rNorm = CloneNormalized(rSlice, Form("refold1D_z_slice_x%d_y%d_normMR", ix, yBin));
                    if (mNorm)
                    {
                        mNorm->SetLineColor(kBlack);
                        mNorm->SetMarkerColor(kBlack);
                        mNorm->SetMarkerStyle(20);
                        mNorm->GetXaxis()->SetTitle("z_{T}");
                        mNorm->GetYaxis()->SetTitle("1/N dN/dz_{T}");
                        mNorm->GetYaxis()->SetRangeUser(0.0, 1.6 * mNorm->GetMaximum());
                        mNorm->SetTitle("");
                        mNorm->Draw("E1");
                        normClonesToDelete2.push_back(mNorm);
                    }
                    if (rNorm)
                    {
                        rNorm->SetLineColor(kBlue + 1);
                        rNorm->SetLineWidth(2);
                        rNorm->Draw("HISTE SAME");
                        normClonesToDelete2.push_back(rNorm);
                    }
                    // keep originals written unmodified
                    if (mSlice)
                        mSlice->Write();
                    if (rSlice)
                        rSlice->Write();
                    // Legend and labels
                    TLegend *leg = new TLegend(0.55, 0.72, 0.88, 0.88);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    if (mNorm)
                        leg->AddEntry(mNorm, "Measured (norm)", "lep");
                    if (rNorm)
                        leg->AddEntry(rNorm, Form("Refolded (iter %d, norm)", bestIter), "l");
                    leg->Draw();
                    TLatex latex;
                    latex.SetNDC();
                    latex.SetTextSize(0.045);
                    double ptLow = (ix - 1 < (int)jetPtBinBorders.size()) ? jetPtBinBorders[ix - 1] : -1;
                    double ptHigh = (ix < (int)jetPtBinBorders.size()) ? jetPtBinBorders[ix] : -1;
                    latex.DrawLatex(0.15, 0.90, Form("%.0f<p_{T}^{jet}<%.0f GeV", ptLow, ptHigh));
                    latex.DrawLatex(0.15, 0.83, Form("%.1f<#eta<%.1f", eta_min, eta_max));
                    padIndex++;
                }
                cSliceMR->SaveAs(Form("%s/projY_meas_vs_refold_best_slices_y%d.png", outDir.c_str(), yBin));
                for (auto *h : mSlices)
                    if (h)
                        delete h;
                for (auto *h : rSlices)
                    if (h)
                        delete h;
                for (auto *h : normClonesToDelete2)
                    if (h)
                        delete h;
                delete cSliceMR;

                // Show unfolded (truth) per x-slice for best iteration, and overlay measured and MC truth (normalized)
                TCanvas *cSliceUnf = new TCanvas(Form("c_projY_unfold_best_slices_y%d", yBin), "Unfolded vs Measured & Truth (best) per jet p_{T} bin", 1400 * ncols, 1000 * nrows);
                cSliceUnf->Divide(ncols, nrows);
                std::vector<TH1D *> uSlices, mSlices2, tSlices;
                std::vector<TH1D *> normClonesToDelete3;
                uSlices.reserve(nPlotBins);
                mSlices2.reserve(nPlotBins);
                tSlices.reserve(nPlotBins);
                normClonesToDelete3.reserve(3 * nPlotBins);
                // Access MC truth 2D distribution from the response for projections
                TH2D *truth2D = dynamic_cast<TH2D *>(response2D.Htruth());
                padIndex = 1;
                for (int ix = 2; ix <= nxBins - 1; ++ix) // Skip first (ix=1) and last (ix=nxBins)
                {
                    cSliceUnf->cd(padIndex);
                    gPad->SetRightMargin(0.01);
                    gPad->SetLeftMargin(0.08);
                    gPad->SetTopMargin(0.01);
                    gPad->SetBottomMargin(0.08);
                    gPad->SetTickx();
                    gPad->SetTicky();
                    // Build projections for this jet pT slice
                    TH1D *mSlice2 = h2_z_vs_jetpt->ProjectionY(Form("meas1D_z_slice_x%d_y%d_forUnfold", ix, yBin), ix, ix);
                    TH1D *uSlice = bestUnf->ProjectionY(Form("unfold1D_z_slice_x%d_y%d", ix, yBin), ix, ix);
                    TH1D *tSlice = truth2D ? truth2D->ProjectionY(Form("truth1D_z_slice_x%d_y%d", ix, yBin), ix, ix) : nullptr;
                    mSlices2.push_back(mSlice2);
                    uSlices.push_back(uSlice);
                    tSlices.push_back(tSlice);

                    // Draw normalized clones for comparison (measured points, unfolded red, truth green)
                    TH1D *mNorm2 = CloneNormalized(mSlice2, Form("meas1D_z_slice_x%d_y%d_normUnf", ix, yBin));
                    TH1D *uNorm = CloneNormalized(uSlice, Form("unfold1D_z_slice_x%d_y%d_normUnf", ix, yBin));
                    TH1D *tNorm = CloneNormalized(tSlice, Form("truth1D_z_slice_x%d_y%d_normUnf", ix, yBin));
                    if (mNorm2)
                    {
                        mNorm2->SetLineColor(kBlack);
                        mNorm2->SetMarkerColor(kBlack);
                        mNorm2->SetMarkerStyle(20);
                        mNorm2->GetXaxis()->SetTitle("z_{T}");
                        mNorm2->GetYaxis()->SetTitle("1/N dN/dz_{T}");
                        mNorm2->GetYaxis()->SetRangeUser(0.0, 1.6 * mNorm2->GetMaximum());
                        mNorm2->SetTitle("");
                        mNorm2->Draw("E1");
                        normClonesToDelete3.push_back(mNorm2);
                    }
                    if (uNorm)
                    {
                        uNorm->SetLineColor(kRed + 1);
                        uNorm->SetLineWidth(2);
                        uNorm->Draw("HIST SAME");
                        normClonesToDelete3.push_back(uNorm);
                    }
                    if (tNorm)
                    {
                        tNorm->SetLineColor(kGreen + 2);
                        tNorm->SetLineWidth(2);
                        tNorm->Draw("HIST SAME");
                        normClonesToDelete3.push_back(tNorm);
                    }
                    // keep originals written unmodified
                    if (mSlice2)
                        mSlice2->Write();
                    if (uSlice)
                        uSlice->Write();
                    if (tSlice)
                        tSlice->Write();

                    // Legend and labels
                    TLegend *leg = new TLegend(0.48, 0.70, 0.88, 0.90);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    if (mNorm2)
                        leg->AddEntry(mNorm2, "Measured (norm)", "lep");
                    if (uNorm)
                        leg->AddEntry(uNorm, Form("Unfolded (iter %d, norm)", bestIter), "l");
                    if (tNorm)
                        leg->AddEntry(tNorm, "Truth (MC, norm)", "l");
                    leg->Draw();

                    TLatex latex;
                    latex.SetNDC();
                    latex.SetTextSize(0.045);
                    double ptLow = (ix - 1 < (int)jetPtBinBorders.size()) ? jetPtBinBorders[ix - 1] : -1;
                    double ptHigh = (ix < (int)jetPtBinBorders.size()) ? jetPtBinBorders[ix] : -1;
                    latex.DrawLatex(0.15, 0.90, Form("%.0f<p_{T}^{jet}<%.0f GeV", ptLow, ptHigh));
                    latex.DrawLatex(0.15, 0.83, Form("%.1f<#eta<%.1f", eta_min, eta_max));
                    padIndex++;
                }
                cSliceUnf->SaveAs(Form("%s/projY_unfold_best_slices_y%d.png", outDir.c_str(), yBin));
                for (auto *h : uSlices)
                    if (h)
                        delete h;
                for (auto *h : mSlices2)
                    if (h)
                        delete h;
                for (auto *h : tSlices)
                    if (h)
                        delete h;
                for (auto *h : normClonesToDelete3)
                    if (h)
                        delete h;
                delete cSliceUnf;
            }

            // Cleanup
            for (auto *h : unfoldedPerIter)
                delete h;
            for (auto *h : refoldedPerIter)
                delete h;
            if (measProjAll)
                delete measProjAll;
            if (bestUnf)
                delete bestUnf;

            delete hMeasTmpl;
            delete hTruthTmpl;
            if (h2_z_vs_jetpt_forUnfoldOwned)
                delete h2_z_vs_jetpt_forUnfoldOwned;
        }
    }
}