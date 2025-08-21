#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2.h>
#include <TTree.h>
#include <TLatex.h>
#include <TColor.h>
#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <algorithm>

// Helper: map underflow/overflow into the boundary bins for 1D histograms
static void MapUnderOverflowToEdges(TH1* h) {
    if (!h) return;
    int nb = h->GetNbinsX();
    double under = h->GetBinContent(0);
    double over = h->GetBinContent(nb+1);
    double underErr = h->GetBinError(0);
    double overErr = h->GetBinError(nb+1);
    if (under != 0.0) {
        double old = h->GetBinContent(1);
        double oldErr = h->GetBinError(1);
        h->SetBinContent(1, old + under);
        h->SetBinError(1, std::sqrt(oldErr*oldErr + underErr*underErr));
        h->SetBinContent(0, 0.0);
        h->SetBinError(0, 0.0);
    }
    if (over != 0.0) {
        double old = h->GetBinContent(nb);
        double oldErr = h->GetBinError(nb);
        h->SetBinContent(nb, old + over);
        h->SetBinError(nb, std::sqrt(oldErr*oldErr + overErr*overErr));
        h->SetBinContent(nb+1, 0.0);
        h->SetBinError(nb+1, 0.0);
    }
}

// Helper: map underflow/overflow into the boundary bins for 2D histograms
static void MapUnderOverflowToEdges(TH2* h) {
    if (!h) return;
    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();
    // create a temporary histogram to accumulate mapped contents
    TH2* tmp = (TH2*)h->Clone(Form("%s_uomap_tmp", h->GetName()));
    tmp->Reset();
    for (int i = 0; i <= nx+1; ++i) {
        for (int j = 0; j <= ny+1; ++j) {
            double val = h->GetBinContent(i, j);
            double err = h->GetBinError(i, j);
            int im = (i < 1) ? 1 : ( (i > nx) ? nx : i );
            int jm = (j < 1) ? 1 : ( (j > ny) ? ny : j );
            double old = tmp->GetBinContent(im, jm);
            double oldErr = tmp->GetBinError(im, jm);
            tmp->SetBinContent(im, jm, old + val);
            tmp->SetBinError(im, jm, std::sqrt(oldErr*oldErr + err*err));
        }
    }
    // copy back the mapped core bins and zero out under/overflow bins
    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            h->SetBinContent(i, j, tmp->GetBinContent(i, j));
            h->SetBinError(i, j, tmp->GetBinError(i, j));
        }
    }
    // zero under/overflow bins (edges)
    for (int j = 0; j <= ny+1; ++j) {
        h->SetBinContent(0, j, 0.0); h->SetBinError(0, j, 0.0);
        h->SetBinContent(nx+1, j, 0.0); h->SetBinError(nx+1, j, 0.0);
    }
    for (int i = 0; i <= nx+1; ++i) {
        h->SetBinContent(i, 0, 0.0); h->SetBinError(i, 0, 0.0);
        h->SetBinContent(i, ny+1, 0.0); h->SetBinError(i, ny+1, 0.0);
    }
    delete tmp;
}

void unfold_new(
                const std::string& infileResponse = "/media/niviths/SSD2/lhcb_analysis_SSD/mc_merge_pPb_Pbp/response_merge.root",
                // const std::string& infileResponse = "/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/response.root",
                // const std::string& infileResponse = "/media/niviths/SSD2/lhcb_analysis_SSD/20250728_pPb_MC_output/20250728_pPb_MC_output_response.root",
                const std::string& outfile = "unfolded_output.root",
                int nIter = 5,
                // const std::vector<std::string>& jetPtBins = {"10_15"},
                const std::vector<std::string>& jetPtBins = {"5_10", "10_15", "15_20", "20_30"},//, "30_50"},
                const std::vector<int>& yBins = {0, 1, 2, 3, 4, 5, 6, 7},
                bool mapUnderflow = true)
{
    // Helper: convert HSV to RGB (all in [0,1]) and return as array
    auto hsv_to_rgb = [](double h, double s, double v, double& r, double& g, double& b){
        if (s <= 0.0) { r = g = b = v; return; }
        double hh = h * 6.0;
        if (hh >= 6.0) hh = 0.0;
        int i = (int)hh;
        double ff = hh - i;
        double p = v * (1.0 - s);
        double q = v * (1.0 - (s * ff));
        double t = v * (1.0 - (s * (1.0 - ff)));
        switch(i) {
            case 0: r=v; g=t; b=p; break;
            case 1: r=q; g=v; b=p; break;
            case 2: r=p; g=v; b=t; break;
            case 3: r=p; g=q; b=v; break;
            case 4: r=t; g=p; b=v; break;
            case 5: default: r=v; g=p; b=q; break;
        }
    };

    // Helper: generate N visually-appealing colors across the HSV spectrum
    auto make_palette = [&](int N, double sat=0.75, double val=0.85){
        std::vector<int> cols; cols.reserve(N);
        for (int i = 0; i < N; ++i) {
            // sample hue in the middle of each bin to avoid duplicating endpoints
            double h = (N > 0) ? ( (i + 0.5) / (double)N ) : 0.0; // 0..1, avoids h==1
            double r,g,b; hsv_to_rgb(h, sat, val, r, g, b);
            int ci = TColor::GetFreeColorIndex();
            TColor* tc = new TColor(ci, r, g, b);
            cols.push_back(ci);
        }
        return cols;
    };
    // Naming scheme
    std::string unfoldedName = "unfolded_zT";
    std::vector<double> yBinBorders = {2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0};

    // Pattern for per-jet measured input files; %s will be replaced with jetPt string
    std::string infilePattern = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA/TagZHistograms_%s.root";

    TFile* finResponse = TFile::Open(infileResponse.c_str());
    if (!finResponse || finResponse->IsZombie()) {
        std::cerr << "Error: Cannot open response file " << infileResponse << std::endl;
        return;
    }

    // Prepare output directory and file (dropped into dated folder)
    // current date string YYYY-MM-DD
    time_t t = time(nullptr);
    struct tm* lt = localtime(&t);
    char dateBuf[32] = {0};
    if (lt) {
        strftime(dateBuf, sizeof(dateBuf), "%Y-%m-%d", lt);
    } else {
        snprintf(dateBuf, sizeof(dateBuf), "unknown-date");
    }

    std::string outDir = Form("unfolded_zT_%s", dateBuf);
    // create directory (recursively) if needed
    std::string mkdirCmd = std::string("mkdir -p ") + outDir;
    int mkret = system(mkdirCmd.c_str()); (void)mkret;

    std::string outPath = outDir + "/" + outfile;
    TFile* fout = TFile::Open(outPath.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: Cannot open output file " << outPath << std::endl;
        return;
    }
    std::cout << "Output directory: " << outDir << "\n";
    std::string treeName = "Response";
    // Get response tree
    TTree* tree = (TTree*)finResponse->Get(treeName.c_str());
    if (!tree) {
        std::cerr << "Warning: Missing response tree: " << treeName << std::endl;
        return;
    }
    // Loop over jetPt and rapidity bins
    for (const auto& jetPt : jetPtBins) {
        // Open the measured input file for this jetPt
        TString infileJet = Form(infilePattern.c_str(), jetPt.c_str());
        TFile* fin = TFile::Open(infileJet.Data());
        if (!fin || fin->IsZombie()) {
            std::cerr << "Warning: Cannot open input file for jetPt " << jetPt << ": " << infileJet << std::endl;
            if (fin) fin->Close();
            continue;
        }

        for (const auto& yBin : yBins) {
            // Build histogram and tree names
            TString histName = Form("promptSignalTagZHist_%s_bin%d", jetPt.c_str(), yBin);
            // TString histName = Form("promptSignalTagZHist_FullyWeighted_%s_bin%d", jetPt.c_str(), yBin);
            TString outHistName = Form("%s_%s_%d", unfoldedName.c_str(), jetPt.c_str(), yBin);

            std::cout << "\n---\nProcessing jetPt bin: " << jetPt
                      << ", yBin: " << yBin << std::endl;
            std::cout << "Histogram name: " << histName << std::endl;

            // Get measured histogram
            TH1D* hMeasured = (TH1D*)fin->Get(histName.Data());
            if (!hMeasured) {
                std::cerr << "Warning: Missing measured histogram: " << histName << std::endl;
                continue;
            }
            std::cout << "Loaded measured histogram with " << hMeasured->GetEntries() << " entries." << std::endl;

            // Variables for tree branches
            float d0_z_det, d0_z_mc;
            float jet_pt_det, jet_pt_mc;
            float d0_eta_det, d0_eta_mc;
            float jet_nconst_det, jet_nconst_mc;
            float jet_dr;
            tree->SetBranchAddress("d0_z_det", &d0_z_det);
            tree->SetBranchAddress("d0_z_mc", &d0_z_mc);
            tree->SetBranchAddress("jet_pt_det", &jet_pt_det);
            tree->SetBranchAddress("jet_pt_mc", &jet_pt_mc);
            tree->SetBranchAddress("d0_eta_det", &d0_eta_det);
            tree->SetBranchAddress("d0_eta_mc", &d0_eta_mc);
            tree->SetBranchAddress("jet_nconst_det", &jet_nconst_det);
            tree->SetBranchAddress("jet_nconst_mc", &jet_nconst_mc);
            tree->SetBranchAddress("jet_dr", &jet_dr);

            // Determine eta bin borders from yBin index
            int yBinIdx = yBin;
            if (yBinIdx < 0 || yBinIdx+1 >= (int)yBinBorders.size()) {
                std::cerr << "Warning: yBin index out of range: " << yBin << std::endl;
                continue;
            }
            double eta_min = yBinBorders[yBinIdx];
            double eta_max = yBinBorders[yBinIdx+1];
            std::cout << "Eta bin: [" << eta_min << ", " << eta_max << ")" << std::endl;

            // Create a human-readable y-bin string from the eta borders for output names
            TString yBinStr = Form("%.1f_%.1f", eta_min, eta_max);

            // Prepare measured output name early so it can be used by refolding/plotting
            TString measuredOutName = Form("measured_zT_%s_%s", jetPt.c_str(), yBinStr.Data());

            // Parse jetPt bin borders from jetPt string (e.g., "5_10")
            double jetpt_min = 0, jetpt_max = 0;
            size_t jet_sep = jetPt.find('_');
            if (jet_sep != std::string::npos) {
                jetpt_min = std::stod(jetPt.substr(0, jet_sep));
                jetpt_max = std::stod(jetPt.substr(jet_sep + 1));
            } else {
                jetpt_min = std::stod(jetPt);
                jetpt_max = jetpt_min + 5.0; // default width if not specified
            }
            std::cout << "jetPt bin: [" << jetpt_min << ", " << jetpt_max << ")" << std::endl;

            // Create a per-jetPt subdirectory under the dated output folder
            std::string jetOutDir = outDir + "/" + jetPt;
            std::string mkdirJetCmd = std::string("mkdir -p ") + jetOutDir;
            int mkjet = system(mkdirJetCmd.c_str()); (void)mkjet;

            // Fill response matrix from tree with eta and jetPt matching for both jet and D0
            // We'll build a histogram-based response using templates so binning matches the
            // measured histogram exactly. Clone the measured histogram to serve as the
            // MC-truth template, then reset it and fill it from the tree while we build
            // the response in a single pass.
            int nBins = hMeasured->GetNbinsX();
            // Create MC truth template by cloning the measured histogram's binning
            TH1D* hPrior = (TH1D*)hMeasured->Clone(Form("hPrior_%s_%s", jetPt.c_str(), yBinStr.Data()));
            hPrior->SetTitle("MC Truth Distribution");
            hPrior->Reset();

            // Construct RooUnfoldResponse using the measured (reco) and truth templates
            RooUnfoldResponse response(hMeasured, hPrior);

            // Optionally map under/overflow into edge bins in the freshly-built response
            if (mapUnderflow) {
                TH2* hRespInit = response.Hresponse();
                if (hRespInit) MapUnderOverflowToEdges(hRespInit);
            }

            // === CRITICAL: SCALE MATCHING ===
            // The response matrix must be built using the same statistical weight as the measured data
            // If measured data has different normalization than MC, we need to match scales
            double measuredTotal = hMeasured->Integral();
            std::cout << "Pre-fill measured total: " << measuredTotal << std::endl;

            Long64_t nEntries = tree->GetEntries();
            int nFilled = 0;
            int nMissed = 0; // Events that pass MC cuts but fail detector cuts
            int nFake = 0;   // Events that pass detector cuts but fail MC cuts

            // Optional: detect common per-event weight branch names and hook them
            float evtWeight = 1.0f;
            if (tree->GetBranch("weight")) tree->SetBranchAddress("weight", &evtWeight);
            else if (tree->GetBranch("evtWeight")) tree->SetBranchAddress("evtWeight", &evtWeight);
            else if (tree->GetBranch("eventWeight")) tree->SetBranchAddress("eventWeight", &evtWeight);
            else if (tree->GetBranch("totalWeight")) tree->SetBranchAddress("totalWeight", &evtWeight);

            // Cache selected events to allow fast, deterministic replays for weighted response rebuilding
            struct CachedEvent { double d0_z_mc; double d0_z_det; bool passMC; bool passDet; double jet_dr; double weight; };
            std::vector<CachedEvent> events; events.reserve(std::min<Long64_t>(nEntries, 500000));

            for (Long64_t i = 0; i < nEntries; ++i) {
                tree->GetEntry(i);

                // Check MC level cuts
                bool passMC = (d0_eta_mc >= eta_min && d0_eta_mc < eta_max &&
                              jet_pt_mc >= jetpt_min && jet_pt_mc < jetpt_max &&
                              jet_nconst_mc > 1);

                // Check detector level cuts  
                bool passDet = (d0_eta_det >= eta_min && d0_eta_det < eta_max &&
                               jet_pt_det >= jetpt_min && jet_pt_det < jetpt_max &&
                               jet_nconst_det > 1);

                // check distance between jets
                if (jet_dr > 0.15) {
                    continue;
                }

                // Fill MC truth template for any event that passes MC selection (use evtWeight if present)
                if (passMC) {
                    hPrior->Fill(d0_z_mc, evtWeight);
                }

                if (passMC && passDet) {
                    // Both MC and detector level pass - normal response
                    response.Fill(d0_z_det, d0_z_mc, evtWeight);
                    ++nFilled;
                } 
                else if (passMC && !passDet) {
                    // MC level passes but detector fails - missed event
                    response.Miss(d0_z_mc, evtWeight);
                    ++nMissed;
                } else if (!passMC && passDet) {
                    // Detector level passes but MC fails - fake event
                    response.Fake(d0_z_det, evtWeight);
                    ++nFake;
                }
                // If neither passes, skip the event

                // Cache this event for later replay (only events reaching here have jet_dr <= cut)
                events.push_back(CachedEvent{(double)d0_z_mc, (double)d0_z_det, passMC, passDet, (double)jet_dr, (double)evtWeight});
            }
            
            std::cout << "Response matrix statistics:" << std::endl;
            std::cout << "  - Filled: " << nFilled << " entries" << std::endl;
            std::cout << "  - Missed: " << nMissed << " entries" << std::endl;
            std::cout << "  - Fake: " << nFake << " entries" << std::endl;
            std::cout << "  - Total processed: " << nEntries << std::endl;
            
            // Validate response matrix has sufficient statistics
            if (nFilled < 100) {
                std::cerr << "ERROR: Insufficient statistics in response matrix (" 
                          << nFilled << " entries). Skipping this bin." << std::endl;
                continue;
            }

            double priorEntries = hPrior->GetEntries();
            double priorIntegral = hPrior->Integral();

            std::cout << "MC truth statistics: " << priorEntries << " entries, integral = " << priorIntegral << std::endl;

            if (priorIntegral <= 0) {
                std::cerr << "ERROR: No MC truth events found for this bin! Check bin boundaries." << std::endl;
                delete hPrior;
                continue;
            }

            // === RESPONSE MATRIX VALIDATION ===
            TH2* hRespMatrix = response.Hresponse();
            if (hRespMatrix) {
                std::cout << "\n=== Response Matrix Diagnostics ===" << std::endl;
                
                // Check response matrix properties
                double totalResp = hRespMatrix->Integral();
                std::cout << "  Response matrix total: " << totalResp << std::endl;
                
                // Check row and column sums (efficiency and purity)
                std::vector<double> rowSums(nBins+1, 0.0);
                std::vector<double> colSums(nBins+1, 0.0);
                
                for (int i = 1; i <= nBins; ++i) {
                    for (int j = 1; j <= nBins; ++j) {
                        double val = hRespMatrix->GetBinContent(i, j);
                        rowSums[i] += val;  // reco bin i efficiency
                        colSums[j] += val;  // truth bin j acceptance
                    }
                }
                
                // Print efficiency (fraction of truth events detected)
                std::cout << "  Truth bin acceptance:" << std::endl;
                for (int j = 1; j <= nBins; ++j) {
                    double truthInBin = hPrior->GetBinContent(j);
                    double efficiency = (truthInBin > 0) ? colSums[j] / truthInBin : 0.0;
                    if (j <= 5 || efficiency < 0.1 || efficiency > 1.1) {  // Show first 5 bins and problematic ones
                        std::cout << "    Truth bin " << j << ": " << colSums[j] << " detected / " 
                                  << truthInBin << " truth = " << efficiency*100 << "% eff" << std::endl;
                    }
                }
                
                // Check for empty rows/columns that could cause problems
                int emptyRows = 0, emptyCols = 0;
                for (int i = 1; i <= nBins; ++i) {
                    if (rowSums[i] <= 0) emptyRows++;
                    if (colSums[i] <= 0) emptyCols++;
                }
                if (emptyRows > 0) std::cout << "  WARNING: " << emptyRows << " empty reco bins in response!" << std::endl;
                if (emptyCols > 0) std::cout << "  WARNING: " << emptyCols << " empty truth bins in response!" << std::endl;
            }

            // Quick runtime checks for underflow/overflow entries in the prior and measured
            // histograms. Warn if many entries fall outside the histogram binning.
            int underBin = 0;
            int overBin = nBins + 1;
            double priorUnder = hPrior->GetBinContent(underBin);
            double priorOver = hPrior->GetBinContent(overBin);
            double priorInclusive = hPrior->Integral(underBin, overBin);

            double measuredUnder = hMeasured->GetBinContent(underBin);
            double measuredOver = hMeasured->GetBinContent(overBin);
            double measuredInclusive = hMeasured->Integral(underBin, overBin);

            const int minCount = 10;
            const double fracThresh = 0.01; // 1%

            double priorThresh = ( (double)minCount > fracThresh*priorInclusive ) ? (double)minCount : fracThresh*priorInclusive;
            if ((priorUnder + priorOver) > priorThresh) {
                std::cerr << "WARNING: Significant under/overflow in MC truth histogram for " << jetPt
                          << " bin " << yBinStr << ": under=" << priorUnder << ", over=" << priorOver
                          << ", total(incl)=" << priorInclusive << std::endl;
                std::cerr << "  Consider extending bin ranges or mapping under/overflow into boundary bins." << std::endl;
            }

            double measThresh = ( (double)minCount > fracThresh*measuredInclusive ) ? (double)minCount : fracThresh*measuredInclusive;
            if ((measuredUnder + measuredOver) > measThresh) {
                std::cerr << "WARNING: Significant under/overflow in measured histogram " << histName
                          << " for " << jetPt << " bin " << yBinStr << ": under=" << measuredUnder
                          << ", over=" << measuredOver << ", total(incl)=" << measuredInclusive << std::endl;
                std::cerr << "  Consider extending bin ranges, re-binning, or mapping under/overflow into boundary bins." << std::endl;
            }

            // Unfold using Bayesian method - each iteration builds on the previous
            fout->cd();

            // Collect unfolded histograms so we can plot all iterations together
            std::vector<TH1D*> unfoldedVec;
            // Store raw unfolded histograms (unnormalized) for gen-weight calculations and convergence
            std::vector<TH1D*> unfoldedSavedVec;
            // Collect refolded plotting clones for combined refolding comparison
            std::vector<TH1D*> refoldedPlotVec;

            // Map under/overflow into boundary bins for measured and prior before further processing
            // MapUnderOverflowToEdges(hPrior);
            // MapUnderOverflowToEdges(hMeasured);

            // Snapshot original prior before any reweight/reset
            TH1D* hPriorOriginal = (TH1D*)hPrior->Clone(Form("hPriorOriginal_%s_%s", jetPt.c_str(), yBinStr.Data()));
            if (hPriorOriginal) hPriorOriginal->SetDirectory(nullptr);

            // Compute det-level weights (measured / response.Hmeasured) and rebuild a weighted response using cached events
            TH1* hRespMeasured = response.Hmeasured();
            TH1D* hWeights = nullptr;
            if (hRespMeasured) {
                // Compute det-level weights from normalized shapes (shape ratio) rather than raw counts
                // This avoids huge global normalization mismatches when measured is luminosity-weighted
                hWeights = (TH1D*)hMeasured->Clone(Form("hWeights_%s_%s", jetPt.c_str(), yBinStr.Data()));
                hWeights->SetDirectory(nullptr);

                const double eps = 1e-12;
                double Imeas = hMeasured->Integral();
                double Iresp = hRespMeasured->Integral();

                int nb = hWeights->GetNbinsX();
                std::vector<double> wvec; wvec.reserve(nb);
                for (int b = 1; b <= nb; ++b) {
                    double measBin = hMeasured->GetBinContent(b);
                    double respBin = hRespMeasured->GetBinContent(b);
                    double numer = (Imeas > eps) ? (measBin / Imeas) : 0.0; // normalized measured shape
                    double denom = (Iresp > eps) ? (respBin / Iresp) : 0.0; // normalized response-measured shape
                    double w = 1.0;
                    if (denom > 1e-15) {
                        w = numer / denom;
                    } else if (respBin > eps) {
                        // Fallback to absolute ratio if resp normalized value is zero but absolute resp bin exists
                        w = (measBin > 0) ? (measBin / respBin) : 1.0;
                    } else {
                        // No information in response for this bin: leave weight at 1 (neutral)
                        w = 1.0;
                    }
                    if (!std::isfinite(w) || w <= 0) w = 1.0;
                    wvec.push_back(w);
                }

                // Determine adaptive clipping from percentiles to avoid extreme outliers
                std::vector<double> wsorted = wvec;
                std::sort(wsorted.begin(), wsorted.end());
                double p50 = 1.0;
                double p99 = 1.0;
                if (!wsorted.empty()) {
                    int idx50 = std::max(0, (int)std::floor(0.50 * (wsorted.size()-1)));
                    int idx99 = std::max(0, (int)std::floor(0.99 * (wsorted.size()-1)));
                    p50 = wsorted[idx50];
                    p99 = wsorted[idx99];
                }
                const double globalMaxCap = 5.0; // absolute ceiling
                const double globalMinCap = 0.2; // floor to avoid zeroing
                double adaptiveCap = std::min(globalMaxCap, (p99 > 0 ? p99 : globalMaxCap));

                // Now fill hWeights applying adaptive clipping and floor
                for (int b = 1; b <= nb; ++b) {
                    double w = wvec[b-1];
                    if (w > adaptiveCap) w = adaptiveCap;
                    if (w < globalMinCap) w = globalMinCap;
                    hWeights->SetBinContent(b, w);
                }

                // Debug: print measured vs response integrals and hWeights stats
                {
                    double measI = hMeasured->Integral();
                    double respMeasI = hRespMeasured ? hRespMeasured->Integral() : 0.0;
                    int nbw = hWeights->GetNbinsX();
                    double wmin = 1e9, wmax = -1.0, wsum = 0.0;
                    for (int b = 1; b <= nbw; ++b) {
                        double ww = hWeights->GetBinContent(b);
                        if (ww < wmin) wmin = ww;
                        if (ww > wmax) wmax = ww;
                        wsum += ww;
                    }
                    double wmean = (nbw>0) ? wsum/nbw : 0.0;
                    std::cout << "  [DEBUG] Det-weight summary: measuredI=" << measI << ", respMeasuredI=" << respMeasI
                              << ", hWeights(min,max,mean)=(" << wmin << "," << wmax << "," << wmean << ")" << std::endl;
                    int nbshow = std::min(5, nbw);
                    for (int b = 1; b <= nbshow; ++b) {
                        std::cout << "    hWeights bin " << b << " = " << hWeights->GetBinContent(b)
                                  << " (meas=" << hMeasured->GetBinContent(b) << ", respMeas=" << (hRespMeasured? hRespMeasured->GetBinContent(b) : 0.0) << ")" << std::endl;
                    }
                    std::cout << "    Det-weight percentiles: p99=" << p99 << ", adaptiveCap=" << adaptiveCap << std::endl;
                }

                // Rebuild weighted response using cached events
                TH1D* hPriorRebuild = (TH1D*)hPrior->Clone(Form("hPriorRebuild_detw_%s_%s", jetPt.c_str(), yBinStr.Data()));
                hPriorRebuild->Reset();
                RooUnfoldResponse responseWeighted(hMeasured, hPriorRebuild);
                for (const auto &ev : events) {
                    if (ev.jet_dr > 0.15) continue; // already ensured but keep guard
                    double detw = 1.0;
                    int wb = hWeights->FindBin(ev.d0_z_det);
                    if (wb >= 1 && wb <= hWeights->GetNbinsX()) detw = hWeights->GetBinContent(wb);
                    if (ev.passMC) hPriorRebuild->Fill(ev.d0_z_mc, detw);
                    if (ev.passMC && ev.passDet) responseWeighted.Fill(ev.d0_z_det, ev.d0_z_mc, detw);
                    else if (ev.passMC && !ev.passDet) responseWeighted.Miss(ev.d0_z_mc, detw);
                    else if (!ev.passMC && ev.passDet) responseWeighted.Fake(ev.d0_z_det, detw);
                }
                    // Use weighted response for initial unfolding
                response = responseWeighted;
                if (mapUnderflow) {
                    TH2* hRespDet = response.Hresponse();
                    if (hRespDet) MapUnderOverflowToEdges(hRespDet);
                }
                // Debug: inspect det-weighted response column sums
                {
                    TH2* hRespDet2 = response.Hresponse();
                    if (hRespDet2) {
                        int nx = hRespDet2->GetNbinsX();
                        int ny = hRespDet2->GetNbinsY();
                        double colmin = 1e99, colmax = -1.0;
                        for (int j = 1; j <= ny; ++j) {
                            double colsum = 0.0;
                            for (int i = 1; i <= nx; ++i) colsum += hRespDet2->GetBinContent(i, j);
                            if (colsum < colmin) colmin = colsum;
                            if (colsum > colmax) colmax = colsum;
                        }
                        std::cout << "  [DEBUG] Det-weighted response: total=" << hRespDet2->Integral()
                                  << ", colSum(min,max)=(" << colmin << "," << colmax << ")" << std::endl;
                        int jm = std::max(1, ny/2);
                        std::cout << "    colSum first=" << " ";
                        for (int j : {1, jm, ny}) {
                            double cs=0; for (int i=1;i<=nx;++i) cs += hRespDet2->GetBinContent(i,j);
                            std::cout << " j="<<j<<"->"<<cs;
                        }
                        std::cout << std::endl;
                    }
                }
                hRespMeasured = response.Hmeasured();
            }

            // Perform unfolding with different iteration counts
            for (int iter = 1; iter <= nIter; ++iter) {
                std::cout << "\n--- Bayesian Unfolding: " << iter << " iterations ---" << std::endl;
                
                // For iter>1, implement gen-level reweighting using previous unfolded result
                if (iter > 1 && unfoldedSavedVec.size() >= (size_t)(iter-1)) {
                    TH1D* hPrevUnfold = unfoldedSavedVec[iter-2];
                    if (hPrevUnfold && hPriorOriginal) {
                        // Build gen-weight: unfolded_prev / original_prior
                        TH1D* hGenW = (TH1D*)hPrevUnfold->Clone(Form("hGenW_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), iter));
                        hGenW->SetDirectory(nullptr);
                        // Build gen-weight from normalized shapes (compare PDFs, not absolute counts)
                        double Iprev = hPrevUnfold->Integral();
                        double Iprior = hPriorOriginal->Integral();
                        const double maxGW = 5.0; // tuneable clamp
                        const double minGW = 0.2; // avoid zeroing bins entirely
                        for (int b = 1; b <= hGenW->GetNbinsX(); ++b) {
                            double prevVal = (Iprev > 0) ? (hPrevUnfold->GetBinContent(b) / Iprev) : 0.0;
                            double priorVal = (Iprior > 0) ? (hPriorOriginal->GetBinContent(b) / Iprior) : 0.0;
                            double gw = 1.0;
                            if (priorVal > 0) gw = prevVal / priorVal;
                            if (gw > maxGW) gw = maxGW;
                            if (gw < minGW) gw = minGW;
                            hGenW->SetBinContent(b, gw);
                        }

                        // Rebuild response with combined det*gen weights using cached events
                        TH1D* hPriorGenRebuild = (TH1D*)hPriorOriginal->Clone(Form("hPriorGenRebuild_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), iter));
                        hPriorGenRebuild->Reset();
                        RooUnfoldResponse responseGenWeighted(hMeasured, hPriorGenRebuild);
                        for (const auto &ev : events) {
                            if (ev.jet_dr > 0.15) continue;
                            double detw = 1.0;
                            if (hWeights) {
                                int wb = hWeights->FindBin(ev.d0_z_det);
                                if (wb >= 1 && wb <= hWeights->GetNbinsX()) detw = hWeights->GetBinContent(wb);
                            }
                            int tb = hGenW->FindBin(ev.d0_z_mc);
                            double genw = 1.0;
                            if (tb >= 1 && tb <= hGenW->GetNbinsX()) genw = hGenW->GetBinContent(tb);
                            double weight = detw * genw * ev.weight;
                            if (ev.passMC) hPriorGenRebuild->Fill(ev.d0_z_mc, weight);
                            if (ev.passMC && ev.passDet) responseGenWeighted.Fill(ev.d0_z_det, ev.d0_z_mc, weight);
                            else if (ev.passMC && !ev.passDet) responseGenWeighted.Miss(ev.d0_z_mc, weight);
                            else if (!ev.passMC && ev.passDet) responseGenWeighted.Fake(ev.d0_z_det, weight);
                        }
                        response = responseGenWeighted;
                        if (mapUnderflow) {
                            TH2* hRespGen = response.Hresponse();
                            if (hRespGen) MapUnderOverflowToEdges(hRespGen);
                        }
                        // Debug: print gen-weight and gen-weighted response summaries
                        if (hGenW) {
                            int nbg = hGenW->GetNbinsX();
                            double gmin = 1e9, gmax = -1.0, gsum=0.0;
                            for (int b=1;b<=nbg;++b){ double gv = hGenW->GetBinContent(b); if (gv<gmin) gmin=gv; if (gv>gmax) gmax=gv; gsum+=gv; }
                            std::cout << "  [DEBUG] Gen-weight summary: min="<<gmin<<", max="<<gmax<<", mean="<<(nbg?gsum/nbg:0) << std::endl;
                        }
                        {
                            TH2* hRespGen2 = response.Hresponse();
                            if (hRespGen2) {
                                int nx = hRespGen2->GetNbinsX();
                                int ny = hRespGen2->GetNbinsY();
                                double colmin = 1e99, colmax=-1.0;
                                for (int j=1;j<=ny;++j){ double cs=0; for (int i=1;i<=nx;++i) cs += hRespGen2->GetBinContent(i,j); if (cs<colmin) colmin=cs; if (cs>colmax) colmax=cs; }
                                std::cout << "  [DEBUG] Gen-weighted response: total="<<hRespGen2->Integral()<<", colSum(min,max)=("<<colmin<<","<<colmax<<")"<<std::endl;
                            }
                        }
                        delete hGenW;
                    }
                }

                // RooUnfoldBayes will now use the possibly gen-weighted response (or the det-weighted one)
                RooUnfoldBayes unfold(&response, hMeasured, iter);
                
                // Optional: Set regularization parameter (smoothing)
                // unfold.SetRegParm(1.0);
                
                TH1D* hUnfolded = (TH1D*)unfold.Hreco();
                if (!hUnfolded) {
                    std::cerr << "ERROR: Unfolding failed for " << iter << " iterations" << std::endl;
                    continue;
                }

                // Validate unfolded result
                double unfoldedIntegral = hUnfolded->Integral();
                double measuredIntegral = hMeasured->Integral();

                std::cout << "Results after " << iter << " iterations:" << std::endl;
                std::cout << "  - Measured integral: " << measuredIntegral << std::endl;
                std::cout << "  - Unfolded integral: " << unfoldedIntegral << std::endl;
                std::cout << "  - Ratio (should be ~1): " << unfoldedIntegral/measuredIntegral << std::endl;

                TString iterName = Form("%s_%s_%s_iter%d", unfoldedName.c_str(), jetPt.c_str(), yBinStr.Data(), iter);
                TH1D* hUnfoldedSaved = (TH1D*)hUnfolded->Clone(iterName);
                // Write saved clone to output file
                hUnfoldedSaved->Write();
                // Keep raw (unnormalized) unfolded histograms for later gen-weight calculations
                unfoldedSavedVec.push_back(hUnfoldedSaved);
                std::cout << "  - Stored raw unfolded histogram for iter " << iter << ": " << iterName << " (bins=" << hUnfoldedSaved->GetNbinsX() << ")" << std::endl;

                // --- Refolding test: fold the unfolded (truth) histogram back through the
                // response TH2 to produce a predicted measured (reco) distribution and
                // compare it to the original measured histogram. Store plotting clones for
                // combined plotting after all iterations.
                TH2* hRespForRefold = response.Hresponse();
                if (hRespForRefold) {
                    TH1D* hRefolded = (TH1D*)hMeasured->Clone(Form("refolded_%s_iter%d", jetPt.c_str(), iter));
                    hRefolded->Reset();

                    int n = hMeasured->GetNbinsX();
                    std::cout << "  [DEBUG] Refold setup: measuredI="<<hMeasured->Integral() << ", unfoldedI="<< hUnfoldedSaved->Integral() << ", responseI="<< hRespForRefold->Integral() << std::endl;
                    // For each truth bin j, determine the MC truth population in response
                    // and distribute the unfolded truth content into reco bins using the
                    // corresponding column of the response matrix.
                    for (int j = 1; j <= n; ++j) {
                        // Sum MC truth in bin j from the response (sum over reco bins)
                        double mcTruthInBin = 0.0;
                        for (int i = 1; i <= n; ++i) {
                            mcTruthInBin += hRespForRefold->GetBinContent(i, j);
                        }
                        double unfoldedVal = hUnfoldedSaved->GetBinContent(j);
                        if (j==1 || j==n/2 || j==n) {
                            double scale_dbg = (mcTruthInBin>0) ? (unfoldedVal/mcTruthInBin) : 0.0;
                            std::cout << "    [DEBUG] refold bin j="<<j<<": mcTruthInBin="<<mcTruthInBin<<", unfoldedVal="<<unfoldedVal<<", scale="<<scale_dbg<<std::endl;
                        }
                        if (mcTruthInBin <= 0) continue; // nothing to map for this truth bin
                        double scale = unfoldedVal / mcTruthInBin;
                        // Distribute into reco bins
                        for (int i = 1; i <= n; ++i) {
                            double add = scale * hRespForRefold->GetBinContent(i, j);
                            double old = hRefolded->GetBinContent(i);
                            hRefolded->SetBinContent(i, old + add);
                        }
                    }

                    // Give refolded a descriptive name and write it
                    TString refoldName = Form("refolded_zT_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), iter);
                    hRefolded->SetName(refoldName);
                    hRefolded->SetTitle(Form("Refolded measured zT %s %s iter %d", jetPt.c_str(), yBinStr.Data(), iter));
                    hRefolded->Write();

                    // Compute a simple chi2 between measured and refolded (using measured errors)
                    double chi2 = 0.0;
                    int ndf = 0;
                    for (int i = 1; i <= n; ++i) {
                        double m = hMeasured->GetBinContent(i);
                        double r = hRefolded->GetBinContent(i);
                        double err = hMeasured->GetBinError(i);
                        if (err <= 0) continue;
                        chi2 += (m - r) * (m - r) / (err * err);
                        ++ndf;
                    }
                    double redChi2 = (ndf > 0) ? chi2 / ndf : -1.0;
                    std::cout << "  Refolding test (iter " << iter << ") -- chi2=" << chi2
                              << ", ndf=" << ndf << ", chi2/ndf=" << redChi2 << std::endl;

                    // === DETAILED CLOSURE TEST DIAGNOSTICS ===
                    if (iter == 1 || iter == nIter) {  // Debug first and last iteration
                        std::cout << "  === Detailed closure diagnostics (iter " << iter << ") ===" << std::endl;
                        
                        // Check unfolded vs truth comparison
                        double unfoldTruthChi2 = 0.0;
                        int unfoldTruthNdf = 0;
                        for (int i = 1; i <= n; ++i) {
                            double unf = hUnfoldedSaved->GetBinContent(i);
                            double tru = hPrior->GetBinContent(i);
                            double unfErr = hUnfoldedSaved->GetBinError(i);
                            if (unfErr > 0) {
                                unfoldTruthChi2 += (unf - tru) * (unf - tru) / (unfErr * unfErr);
                                unfoldTruthNdf++;
                            }
                        }
                        double unfoldTruthRedChi2 = (unfoldTruthNdf > 0) ? unfoldTruthChi2 / unfoldTruthNdf : -1.0;
                        std::cout << "    Unfolded vs Truth: chi2=" << unfoldTruthChi2 
                                  << ", ndf=" << unfoldTruthNdf << ", chi2/ndf=" << unfoldTruthRedChi2 << std::endl;
                        
                        // Bin-by-bin comparison for problematic bins
                        std::cout << "    Bin-by-bin closure check (first 10 bins):" << std::endl;
                        for (int i = 1; i <= std::min(10, n); ++i) {
                            double meas = hMeasured->GetBinContent(i);
                            double refo = hRefolded->GetBinContent(i);
                            double unfo = hUnfoldedSaved->GetBinContent(i);
                            double truth = hPrior->GetBinContent(i);
                            double ratio = (meas > 0) ? refo/meas : 0.0;
                            std::cout << "      Bin " << i << ": meas=" << meas << ", refold=" << refo 
                                      << " (ratio=" << ratio << "), unfold=" << unfo << ", truth=" << truth << std::endl;
                        }
                        
                        // Check for systematic biases
                        double measTotal = hMeasured->Integral();
                        double refoldTotal = hRefolded->Integral();
                        double unfoldTotal = hUnfoldedSaved->Integral();
                        double truthTotal = hPrior->Integral();
                        
                        std::cout << "    Total integrals: meas=" << measTotal << ", refold=" << refoldTotal 
                                  << " (ratio=" << refoldTotal/measTotal << ")" << std::endl;
                        std::cout << "                     unfold=" << unfoldTotal << ", truth=" << truthTotal 
                                  << " (ratio=" << unfoldTotal/truthTotal << ")" << std::endl;
                        
                        // Alternative refolding method for cross-check
                        TH1D* hRefoldAlt = (TH1D*)hMeasured->Clone(Form("refold_alt_%s_iter%d", jetPt.c_str(), iter));
                        hRefoldAlt->Reset();
                        
                        // Direct matrix multiplication: R_ij * T_j = M_i
                        for (int i = 1; i <= n; ++i) {
                            double refoldSum = 0.0;
                            for (int j = 1; j <= n; ++j) {
                                double respVal = hRespForRefold->GetBinContent(i, j);
                                double unfoldVal = hUnfoldedSaved->GetBinContent(j);
                                refoldSum += respVal * unfoldVal;
                            }
                            hRefoldAlt->SetBinContent(i, refoldSum);
                        }
                        
                        // Compare alternative refolding method
                        double altChi2 = 0.0;
                        int altNdf = 0;
                        for (int i = 1; i <= n; ++i) {
                            double m = hMeasured->GetBinContent(i);
                            double r = hRefoldAlt->GetBinContent(i);
                            double err = hMeasured->GetBinError(i);
                            if (err > 0) {
                                altChi2 += (m - r) * (m - r) / (err * err);
                                altNdf++;
                            }
                        }
                        double altRedChi2 = (altNdf > 0) ? altChi2 / altNdf : -1.0;
                        std::cout << "    Alternative refolding: chi2=" << altChi2 
                                  << ", ndf=" << altNdf << ", chi2/ndf=" << altRedChi2 << std::endl;
                        
                        delete hRefoldAlt;
                    }

                    // Create a plotting clone (normalized) and store for combined plotting
                    TH1D* hRefoldPlot = (TH1D*)hRefolded->Clone(Form("%s_plot", refoldName.Data()));
                    double ri_plot = hRefoldPlot->Integral(); if (ri_plot > 0) hRefoldPlot->Scale(1.0/ri_plot);
                    // Ensure internal Entries counter matches current integral after scaling
                    hRefoldPlot->SetEntries((Long64_t)hRefoldPlot->Integral());
                    refoldedPlotVec.push_back(hRefoldPlot);

                    // cleanup temporary refolded histogram
                    delete hRefolded;
                } else {
                    std::cerr << "Warning: cannot perform refolding test; response.Hresponse() is null." << std::endl;
                }

                // === UNFOLDING CONVERGENCE CHECK ===
                if (iter > 1) {
                    // Compare with previous iteration to check convergence
                    TH1D* hPrevUnfolded = unfoldedVec.back();  // Previous iteration (already stored)
                    double convChi2 = 0.0;
                    int convNdf = 0;
                    double maxDiff = 0.0;
                    int maxDiffBin = 0;
                    
                    for (int i = 1; i <= hUnfoldedSaved->GetNbinsX(); ++i) {
                        double curr = hUnfoldedSaved->GetBinContent(i);
                        double prev = hPrevUnfolded->GetBinContent(i);
                        double currErr = hUnfoldedSaved->GetBinError(i);
                        
                        if (currErr > 0) {
                            double diff = curr - prev;
                            convChi2 += diff * diff / (currErr * currErr);
                            convNdf++;
                            
                            if (std::abs(diff) > maxDiff) {
                                maxDiff = std::abs(diff);
                                maxDiffBin = i;
                            }
                        }
                    }
                    
                    double convRedChi2 = (convNdf > 0) ? convChi2 / convNdf : -1.0;
                    std::cout << "  Convergence test: chi2=" << convChi2 << ", ndf=" << convNdf 
                              << ", chi2/ndf=" << convRedChi2 << std::endl;
                    std::cout << "    Max change: " << maxDiff << " in bin " << maxDiffBin << std::endl;
                    
                    if (convRedChi2 < 0.1) {
                        std::cout << "  --> Unfolding appears converged (chi2/ndf < 0.1)" << std::endl;
                    }
                }

                // Create a separate plotting clone and normalize it
                TH1D* hUnfoldedPlot = (TH1D*)hUnfoldedSaved->Clone(Form("%s_plot", iterName.Data()));
                double ui = hUnfoldedPlot->Integral();
                if (ui > 0) hUnfoldedPlot->Scale(1.0 / ui);
                // Keep Entries consistent after normalization
                hUnfoldedPlot->SetEntries((Long64_t)hUnfoldedPlot->Integral());
                unfoldedVec.push_back(hUnfoldedPlot);
                std::cout << "  - Saved as: " << iterName << std::endl;
            }
            
            // === FINAL CLOSURE TEST SUMMARY ===
            std::cout << "\n=== CLOSURE TEST SUMMARY ===" << std::endl;
            if (!refoldedPlotVec.empty()) {
                // Find best iteration based on chi2
                double bestChi2 = 1e9;
                int bestIter = 1;
                
                for (int iter = 1; iter <= nIter; ++iter) {
                    // Recalculate chi2 for each iteration (this could be stored from above)
                    TString refoldName = Form("refolded_zT_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), iter);
                    TH1D* hRefoldStored = (TH1D*)fout->Get(refoldName);
                    if (hRefoldStored) {
                        double chi2 = 0.0;
                        int ndf = 0;
                        for (int i = 1; i <= nBins; ++i) {
                            double m = hMeasured->GetBinContent(i);
                            double r = hRefoldStored->GetBinContent(i);
                            double err = hMeasured->GetBinError(i);
                            if (err > 0) {
                                chi2 += (m - r) * (m - r) / (err * err);
                                ndf++;
                            }
                        }
                        double redChi2 = (ndf > 0) ? chi2 / ndf : 1e9;
                        std::cout << "  Iteration " << iter << ": chi2/ndf = " << redChi2;
                        if (redChi2 < bestChi2) {
                            bestChi2 = redChi2;
                            bestIter = iter;
                            std::cout << " <- BEST";
                        }
                        std::cout << std::endl;
                    }
                }
                
                std::cout << "  Recommended iteration: " << bestIter << " (chi2/ndf = " << bestChi2 << ")" << std::endl;
                
                if (bestChi2 > 2.0) {
                    std::cout << "  WARNING: Poor closure (chi2/ndf > 2). Possible issues:" << std::endl;
                    std::cout << "    - Response matrix statistics too low" << std::endl;
                    std::cout << "    - Bin migration not well modeled" << std::endl;
                    std::cout << "    - Systematic differences between measured and MC" << std::endl;
                } else if (bestChi2 > 1.5) {
                    std::cout << "  CAUTION: Moderate closure issues (chi2/ndf > 1.5)" << std::endl;
                } else {
                    std::cout << "  GOOD: Closure test passed (chi2/ndf < 1.5)" << std::endl;
                }
            }

            // Save original measured zT distribution
            TH1D* hMeasuredClone = (TH1D*)hMeasured->Clone(measuredOutName);
            hMeasuredClone->Write();
            // Create plotting clone for measured and normalize
            TH1D* hMeasuredPlot = (TH1D*)hMeasuredClone->Clone(Form("%s_plot", measuredOutName.Data()));
            double mi = hMeasuredPlot->Integral();
            if (mi > 0) hMeasuredPlot->Scale(1.0 / mi);
            // Keep Entries consistent after normalization
            hMeasuredPlot->SetEntries((Long64_t)hMeasuredPlot->Integral());
            
            // Save MC truth distribution for comparison
            TString truthOutName = Form("truth_zT_%s_%s", jetPt.c_str(), yBinStr.Data());
            hPrior->SetName(truthOutName);
            hPrior->Write();
            // Create plotting clone for truth and normalize
            TH1D* hPriorPlot = (TH1D*)hPrior->Clone(Form("%s_plot", truthOutName.Data()));
            double ti = hPriorPlot->Integral();
            if (ti > 0) hPriorPlot->Scale(1.0 / ti);
            // Keep Entries consistent after normalization
            hPriorPlot->SetEntries((Long64_t)hPriorPlot->Integral());
            
            // Save response matrix object
            TString responseOutName = Form("response_%s_%s", jetPt.c_str(), yBinStr.Data());
            // Ensure response under/overflow mapping applied if requested
            if (mapUnderflow) {
                TH2* hRespFinal = response.Hresponse();
                if (hRespFinal) MapUnderOverflowToEdges(hRespFinal);
            }
            response.Write(responseOutName);

            // Plot and save the response matrix (TH2) into the output directory as PNG
            TH2* hResp = response.Hresponse();
            if (hResp) {
                // Clone to give a stable name and write to file
                TH2* hRespClone = (TH2*)hResp->Clone(Form("hResponse_%s_%s", jetPt.c_str(), yBinStr.Data()));
                hRespClone->SetTitle(Form("Response matrix %s bin %s", jetPt.c_str(), yBinStr.Data()));
                hRespClone->Write();

                // Draw on a TCanvas and save as PNG so the image is rendered correctly
                TCanvas* cResp = new TCanvas(Form("c_resp_%s_%s", jetPt.c_str(), yBinStr.Data()), "Response Matrix", 800, 600);
                cResp->cd();
                hRespClone->Draw("COLZ");
                TString pngName = Form("%s/response_%s_%s.png", jetOutDir.c_str(), jetPt.c_str(), yBinStr.Data());
                cResp->SaveAs(pngName.Data());
                delete cResp;
                delete hRespClone;
                std::cout << "  - Saved response image: " << pngName << std::endl;
            } else {
                std::cerr << "Warning: response.Hresponse() returned null for " << histName << std::endl;
            }
            
            // Create a comparison plot: measured, MC truth, and all unfolded iterations
            {
                TCanvas* cComp = new TCanvas(Form("c_comp_%s_%s", jetPt.c_str(), yBinStr.Data()), "zT comparison", 900, 700);
                cComp->SetRightMargin(0.02);
                // Style measured and truth (use plotting clones)
                hMeasuredPlot->SetLineColor(kBlack);
                hMeasuredPlot->SetLineWidth(2);
                hMeasuredPlot->SetMarkerStyle(20);

                hPriorPlot->SetLineColor(kRed);
                hPriorPlot->SetLineStyle(2);
                hPriorPlot->SetLineWidth(2);

                // Find maximum for axis scaling among normalized plots
                double maxVal = hMeasuredPlot->GetMaximum();
                if (hPriorPlot->GetMaximum() > maxVal) maxVal = hPriorPlot->GetMaximum();
                for (size_t ii = 0; ii < unfoldedVec.size(); ++ii) {
                    if (unfoldedVec[ii]->GetMaximum() > maxVal) maxVal = unfoldedVec[ii]->GetMaximum();
                }
                hMeasuredPlot->SetMaximum(maxVal * 1.2);

                // Draw base measured and truth plots
                hMeasuredPlot->GetXaxis()->SetTitle("z_{T}");
                hMeasuredPlot->GetYaxis()->SetTitle("dN/dz_{T} (norm)");
                hMeasuredPlot->Draw("E");
                hPriorPlot->Draw("HISTSAME");

                // Generate a smooth palette for unfolded iterations
                std::vector<int> ucols = make_palette((int)unfoldedVec.size());
                for (size_t ii = 0; ii < unfoldedVec.size(); ++ii) {
                    int col = ucols[ii % ucols.size()];
                    unfoldedVec[ii]->SetLineColor(col);
                    unfoldedVec[ii]->SetLineStyle(1);
                    unfoldedVec[ii]->SetLineWidth(2);
                    unfoldedVec[ii]->Draw("HISTSAME");
                }

                // Legend
                // Move legend to top-left
                TLegend leg(0.14, 0.65, 0.45, 0.88);
                leg.SetBorderSize(0);
                leg.SetTextSize(0.03);
                leg.SetMargin(0.12);
                leg.SetFillStyle(0); leg.SetFillColor(0);
                leg.AddEntry(hMeasuredPlot, "Measured", "lep");
                leg.AddEntry(hPriorPlot, "MC truth", "l");
                for (size_t ii = 0; ii < unfoldedVec.size(); ++ii) {
                    TString lbl = Form("Unfold iter %zu", ii+1);
                    leg.AddEntry(unfoldedVec[ii], lbl, "l");
                }
                leg.Draw();

                TLatex texRef2;
                texRef2.SetNDC();
                texRef2.SetTextSize(0.03);
                texRef2.SetTextAlign(31);
                texRef2.DrawLatex(0.95, 0.85, Form("#font[12]{LHCb} work-in-progress"));
                texRef2.DrawLatex(0.95, 0.81, Form("Jet p_{T}: %s", jetPt.c_str()));
                texRef2.DrawLatex(0.95, 0.77, Form("y-bin: %s", yBinStr.Data()));

                // Save comparison image and write canvas to file
                TString cmpName = Form("%s/compare_zT_%s_%s.png", jetOutDir.c_str(), jetPt.c_str(), yBinStr.Data());
                cComp->SaveAs(cmpName.Data());
                cComp->Write();
                delete cComp;
                std::cout << "  - Saved comparison image: " << cmpName << std::endl;

                // Cleanup plotting clones from memory (saved clones remain in the ROOT file)
                for (auto h : unfoldedVec) delete h;
                delete hMeasuredPlot;
                delete hPriorPlot;
            }

            // Combined refolding comparison: overlay all refolded iterations (single canvas)
            if (!refoldedPlotVec.empty()) {
                TCanvas* cRefAll = new TCanvas(Form("c_refold_all_%s_%s", jetPt.c_str(), yBinStr.Data()),
                                                Form("Refolding all iters %s %s", jetPt.c_str(), yBinStr.Data()), 1000, 900);
                cRefAll->SetRightMargin(0.02);
                cRefAll->cd();

                // Top overlay: normalized measured and all refolded iterations
                TH1D* hMeasuredNorm2 = (TH1D*)hMeasured->Clone(Form("%s_norm_allref", measuredOutName.Data()));
                double mi2 = hMeasuredNorm2->Integral(); if (mi2 > 0) hMeasuredNorm2->Scale(1.0/mi2);
                hMeasuredNorm2->SetLineColor(kBlack); hMeasuredNorm2->SetMarkerStyle(20);
                double maxVal2 = hMeasuredNorm2->GetMaximum();
                hMeasuredNorm2->GetXaxis()->SetTitle("z_{T}");
                hMeasuredNorm2->GetYaxis()->SetTitle("dN/dz_{T} (norm)");
                hMeasuredNorm2->Draw("E");

                std::vector<int> rcols = make_palette((int)refoldedPlotVec.size());
                for (size_t ii = 0; ii < refoldedPlotVec.size(); ++ii) {
                    TH1D* hr = refoldedPlotVec[ii];
                    int col = rcols[ii % rcols.size()];
                    hr->SetLineColor(col); hr->SetLineWidth(2); hr->SetLineStyle(1);
                    if (hr->GetMaximum() > maxVal2) maxVal2 = hr->GetMaximum();
                    hr->Draw("HISTSAME");
                }
                hMeasuredNorm2->SetMaximum(maxVal2 * 1.2);

                // Legend
                TLegend legRef(0.15, 0.65, 0.6, 0.88); legRef.SetBorderSize(0);
                legRef.SetTextSize(0.030);
                legRef.SetMargin(0.12);
                legRef.SetFillStyle(0); legRef.SetFillColor(0);
                legRef.AddEntry(hMeasuredNorm2, "Measured", "lep");
                for (size_t ii = 0; ii < refoldedPlotVec.size(); ++ii) {
                    TString lbl = Form("Refold iter %zu", ii+1);
                    legRef.AddEntry(refoldedPlotVec[ii], lbl, "l");
                }
                legRef.Draw();

                TLatex texRef;
                texRef.SetNDC();
                texRef.SetTextSize(0.03);
                texRef.SetTextAlign(31);
                texRef.DrawLatex(0.95, 0.85, Form("#font[12]{LHCb} work-in-progress"));
                texRef.DrawLatex(0.95, 0.81, Form("Jet p_{T}: %s", jetPt.c_str()));
                texRef.DrawLatex(0.95, 0.77, Form("y-bin: %s", yBinStr.Data()));

                TString allPng = Form("%s/refold_compare_all_%s_%s.png", jetOutDir.c_str(), jetPt.c_str(), yBinStr.Data());
                cRefAll->SaveAs(allPng.Data());
                cRefAll->Write();
                std::cout << "  - Saved combined refolding comparison image: " << allPng << std::endl;

                // Cleanup
                delete hMeasuredNorm2; delete cRefAll;
                for (auto h : refoldedPlotVec) delete h;
            }
            
            std::cout << "\nSaved outputs:" << std::endl;
            std::cout << "  - Measured zT: " << measuredOutName << std::endl;
            std::cout << "  - MC truth zT: " << truthOutName << std::endl;
            std::cout << "  - Response matrix: " << responseOutName << std::endl;
            
            // Cleanup
            delete hPrior;
        }
        // Close per-jet input file
        fin->Close();
        delete fin;
    }

    fout->Close();

    std::cout << "Unfolding complete. Output saved to " << outPath << std::endl;
}

// Example usage:
// unfold_new("input.root", "unfolded_output.root", 4, {"5_10", "10_20"}, {"2.5_3.0", "3.0_3.5"});