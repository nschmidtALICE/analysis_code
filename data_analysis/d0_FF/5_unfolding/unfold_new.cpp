#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2.h>
#include <TTree.h>
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

// -----------------------------------------------------------------------------
// Utility helpers (file-scope)
// -----------------------------------------------------------------------------

// HSV -> RGB helper (all channels in [0,1])
static void HSVtoRGB(double h, double s, double v, double &r, double &g, double &b)
{
    if (s <= 0.0)
    {
        r = g = b = v;
        return;
    }
    double hh = h * 6.0;
    if (hh >= 6.0)
        hh = 0.0;
    int i = (int)hh;
    double ff = hh - i;
    double p = v * (1.0 - s);
    double q = v * (1.0 - (s * ff));
    double t = v * (1.0 - (s * (1.0 - ff)));
    switch (i)
    {
    case 0:
        r = v;
        g = t;
        b = p;
        break;
    case 1:
        r = q;
        g = v;
        b = p;
        break;
    case 2:
        r = p;
        g = v;
        b = t;
        break;
    case 3:
        r = p;
        g = q;
        b = v;
        break;
    case 4:
        r = t;
        g = p;
        b = v;
        break;
    default:
        r = v;
        g = p;
        b = q;
        break;
    }
}

// Generate N distinct colors across HSV hue circle
static std::vector<int> MakePalette(int N, double sat = 0.75, double val = 0.85)
{
    std::vector<int> cols;
    cols.reserve(N);
    for (int i = 0; i < N; ++i)
    {
        double h = (N > 0 ? (i + 0.5) / double(N) : 0.0); // center of hue bins
        double r, g, b;
        HSVtoRGB(h, sat, val, r, g, b);
        int ci = TColor::GetFreeColorIndex();
        new TColor(ci, r, g, b); // keep pointer alive implicitly
        cols.push_back(ci);
    }
    return cols;
}

// -----------------------------------------------------------------------------
// Chi2 helper
// -----------------------------------------------------------------------------
struct Chi2Result
{
    double chi2 = 0.0;
    int ndf = 0;
    double redchi2 = -1.0;
};
static Chi2Result ComputeChi2(const TH1 *hData, const TH1 *hModel)
{
    Chi2Result r;
    if (!hData || !hModel)
        return r;
    int nb = hData->GetNbinsX();
    for (int i = 1; i <= nb; ++i)
    {
        double d = hData->GetBinContent(i);
        double m = hModel->GetBinContent(i);
        double e = hData->GetBinError(i);
        if (e > 0)
        {
            r.chi2 += (d - m) * (d - m) / (e * e);
            ++r.ndf;
        }
    }
    if (r.ndf > 0)
        r.redchi2 = r.chi2 / r.ndf;
    return r;
}

// Helper: map underflow/overflow into the boundary bins for 1D histograms
static void MapUnderOverflowToEdges(TH1 *h)
{
    if (!h)
        return;
    int nb = h->GetNbinsX();
    double under = h->GetBinContent(0);
    double over = h->GetBinContent(nb + 1);
    double underErr = h->GetBinError(0);
    double overErr = h->GetBinError(nb + 1);
    if (under != 0.0)
    {
        double old = h->GetBinContent(1);
        double oldErr = h->GetBinError(1);
        h->SetBinContent(1, old + under);
        h->SetBinError(1, std::sqrt(oldErr * oldErr + underErr * underErr));
        h->SetBinContent(0, 0.0);
        h->SetBinError(0, 0.0);
    }
    if (over != 0.0)
    {
        double old = h->GetBinContent(nb);
        double oldErr = h->GetBinError(nb);
        h->SetBinContent(nb, old + over);
        h->SetBinError(nb, std::sqrt(oldErr * oldErr + overErr * overErr));
        h->SetBinContent(nb + 1, 0.0);
        h->SetBinError(nb + 1, 0.0);
    }
}

// Helper: map underflow/overflow into the boundary bins for 2D histograms
static void MapUnderOverflowToEdges(TH2 *h)
{
    if (!h)
        return;
    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();
    // create a temporary histogram to accumulate mapped contents
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
    // copy back the mapped core bins and zero out under/overflow bins
    for (int i = 1; i <= nx; ++i)
    {
        for (int j = 1; j <= ny; ++j)
        {
            h->SetBinContent(i, j, tmp->GetBinContent(i, j));
            h->SetBinError(i, j, tmp->GetBinError(i, j));
        }
    }
    // zero under/overflow bins (edges)
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

// are now fixed to their previous default = true behavior to reduce complexity.
void unfold_new(
    const std::string &outfile = "unfolded_output.root",
    int nIter = 4,
    const std::vector<std::string> &jetPtBins = {"5_10", "10_15", "15_20", "20_30"},
    const std::vector<int> &yBins = {0, 1, 2, 3, 4, 5, 6, 7},
    bool isClosure = false,
    bool verbose = false)
{
    // select input file depending on closure flag
    const std::string &infileResponse = isClosure
                                            ? "/media/niviths/SSD2/lhcb_analysis_SSD/mc_merge_pPb_Pbp/20250728_pPb_MC_output_response.root"
                                            // : "/media/niviths/SSD2/lhcb_analysis_SSD/mc_merge_pPb_Pbp/response_merged.root";
                                            : "/media/niviths/SSD2/lhcb_analysis_SSD/mc_merge_MBonly/MBresponse.root";
    // (Color palette helpers moved to file scope.)
    // Naming scheme
    std::string unfoldedName = "unfolded_zT";
    std::vector<double> yBinBorders = {2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0};

    // Pattern for per-jet measured input files; %s will be replaced with jetPt string
    std::string infilePattern = isClosure ? "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_MC/TagZHistograms_%s.root" : "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA/TagZHistograms_%s.root";

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

    std::string outDir = isClosure ? Form("unfolded_zT_closure_%s", dateBuf) : Form("unfolded_zT_%s", dateBuf);
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
    // Loop over jetPt and rapidity bins
    for (const auto &jetPt : jetPtBins)
    {
        struct ResidualEntry
        {
            TString yBinStr;
        };
        struct RefoldEntry
        {
            TString yBinStr;
            TH1D *meas;
            TH1D *refold;
            int iter;
        };
        struct CompareEntry
        {
            TString yBinStr;
            TH1D *meas;
            TH1D *truth;
            TH1D *unfolded;
            int iter;
        };
        struct MFMFEntry
        {
            TString yBinStr;
            TH1D *meas;
            TH1D *truth;
            TH1D *missed;
            TH1D *fake;
        };
        std::vector<ResidualEntry> residualPerYBin;
        std::vector<TString> responsePerYBin; // y-bin strings with a valid response
        std::vector<RefoldEntry> refoldPerYBin;
        std::vector<CompareEntry> comparePerYBin;
        std::vector<MFMFEntry> mfmfPerYBin;
        // Open the measured input file for this jetPt
        TString infileJet = Form(infilePattern.c_str(), jetPt.c_str());
        TFile *fin = TFile::Open(infileJet.Data());
        if (!fin || fin->IsZombie())
        {
            std::cerr << "Warning: Cannot open input file for jetPt " << jetPt << ": " << infileJet << std::endl;
            if (fin)
                fin->Close();
            continue;
        }

    // Create per-jetPt output directories (one place only) so multiDir is
    // available after the rapidity loop for multi-panel summaries.
    std::string jetOutDir = outDir + "/" + jetPt;
    std::string mkdirJetCmd = std::string("mkdir -p ") + jetOutDir;
    int mkjet = system(mkdirJetCmd.c_str());
    (void)mkjet;
    std::string diagDir = jetOutDir + "/diagnostics";
    std::string mkdirDiagCmd = std::string("mkdir -p ") + diagDir;
    int mkdiag = system(mkdirDiagCmd.c_str());
    (void)mkdiag;
    std::string multiDir = jetOutDir + "/multi";
    std::string mkdirMultiCmd = std::string("mkdir -p ") + multiDir;
    int mkmulti = system(mkdirMultiCmd.c_str());
    (void)mkmulti;

        for (const auto &yBin : yBins)
        {
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
                continue;
            }
            if (verbose)
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
            if (yBinIdx < 0 || yBinIdx + 1 >= (int)yBinBorders.size())
            {
                std::cerr << "Warning: yBin index out of range: " << yBin << std::endl;
                continue;
            }
            double eta_min = yBinBorders[yBinIdx];
            double eta_max = yBinBorders[yBinIdx + 1];
            if (verbose)
                std::cout << "Eta bin: [" << eta_min << ", " << eta_max << ")" << std::endl;

            // Create a human-readable y-bin string from the eta borders for output names
            TString yBinStr = Form("%.1f_%.1f", eta_min, eta_max);

            // Prepare measured output name early so it can be used by refolding/plotting
            TString measuredOutName = Form("measured_zT_%s_%s", jetPt.c_str(), yBinStr.Data());

            // Parse jetPt bin borders from jetPt string (e.g., "5_10")
            double jetpt_min = 0, jetpt_max = 0;
            size_t jet_sep = jetPt.find('_');
            if (jet_sep != std::string::npos)
            {
                jetpt_min = std::stod(jetPt.substr(0, jet_sep));
                jetpt_max = std::stod(jetPt.substr(jet_sep + 1));
            }
            else
            {
                jetpt_min = std::stod(jetPt);
                jetpt_max = jetpt_min + 5.0; // default width if not specified
            }
            if (verbose)
                std::cout << "jetPt bin: [" << jetpt_min << ", " << jetpt_max << ")" << std::endl;

            // per-jetPt directories (jetOutDir/diagDir/multiDir) are created once per-jetPt
            // above the rapidity loop so they are available here. No action needed.

            // Fill response matrix from tree with eta and jetPt matching for both jet and D0
            // We'll build a histogram-based response using templates so binning matches the
            // measured histogram exactly. Clone the measured histogram to serve as the
            // MC-truth template, then reset it and fill it from the tree while we build
            // the response in a single pass.
            int nBins = hMeasured->GetNbinsX();
            // Create MC truth template by cloning the measured histogram's binning
            TH1D *hPrior = (TH1D *)hMeasured->Clone(Form("hPrior_%s_%s", jetPt.c_str(), yBinStr.Data()));
            hPrior->SetTitle("MC Truth Distribution");
            hPrior->Reset();

            // Additional distributions: missed (truth-level events failing reco) and fake (reco-level events failing truth)
            TH1D *hMissedTruth = (TH1D *)hPrior->Clone(Form("missed_zT_%s_%s", jetPt.c_str(), yBinStr.Data()));
            hMissedTruth->SetTitle("Missed truth events (pass MC, fail det)");
            hMissedTruth->Reset();
            TH1D *hFakeReco = (TH1D *)hMeasured->Clone(Form("fake_zT_%s_%s", jetPt.c_str(), yBinStr.Data()));
            hFakeReco->SetTitle("Fake reco events (fail MC, pass det)");
            hFakeReco->Reset();

            // Construct RooUnfoldResponse using the measured (reco) and truth templates
            RooUnfoldResponse response(hMeasured, hPrior);

            if (verbose)
                std::cout << "[INFO] Minimal mode: raw response, no det/gen reweighting." << std::endl;

            // map under/overflow into edge bins in the freshly-built response
            {
                TH2 *hRespInit = response.Hresponse();
                if (hRespInit)
                    MapUnderOverflowToEdges(hRespInit);
            }

            // === CRITICAL: SCALE MATCHING ===
            // The response matrix must be built using the same statistical weight as the measured data
            // If measured data has different normalization than MC, we need to match scales
            double measuredTotal = hMeasured->Integral();
            if (verbose)
                std::cout << "Pre-fill measured total: " << measuredTotal << std::endl;

            Long64_t nEntries = tree->GetEntries();
            int nFilled = 0;
            int nMissed = 0; // Events that pass MC cuts but fail detector cuts

            // Optional: detect common per-event weight branch names and hook them
            float evtWeight = 1.0f;
            if (tree->GetBranch("weight"))
                tree->SetBranchAddress("weight", &evtWeight);
            else if (tree->GetBranch("evtWeight"))
                tree->SetBranchAddress("evtWeight", &evtWeight);
            else if (tree->GetBranch("eventWeight"))
                tree->SetBranchAddress("eventWeight", &evtWeight);
            else if (tree->GetBranch("totalWeight"))
                tree->SetBranchAddress("totalWeight", &evtWeight);

            for (Long64_t i = 0; i < nEntries; ++i)
            {
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
                if (jet_dr > 0.15)
                {
                    continue;
                }

                // Fill MC truth template for any event that passes MC selection (use evtWeight if present)
                if (passMC)
                {
                    hPrior->Fill(d0_z_mc, evtWeight);
                }

                if (passMC && passDet)
                {
                    // Both MC and detector level pass - normal response
                    response.Fill(d0_z_det, d0_z_mc, evtWeight);
                    ++nFilled;
                }
                else if (passMC && !passDet)
                {
                    // MC level passes but detector fails - missed event
                    response.Miss(d0_z_mc, evtWeight);
                    ++nMissed;
                    hMissedTruth->Fill(d0_z_mc, evtWeight);
                }
                else if (!passMC && passDet)
                {
                    // Detector level passes but MC fails - fake events ignored in response
                    hFakeReco->Fill(d0_z_det, evtWeight);
                }
                // If neither passes, skip the event

                // (Event-level caching removed.)
            }

            if (verbose)
            {
                std::cout << "Response matrix statistics:" << std::endl;
                std::cout << "  - Filled: " << nFilled << " entries" << std::endl;
                std::cout << "  - Missed: " << nMissed << " entries" << std::endl;
                std::cout << "  - Fake: (ignored in response, diagnostics only)" << std::endl;
                std::cout << "  - Total processed: " << nEntries << std::endl;
            }

            // Validate response matrix has sufficient statistics
            if (nFilled < 50)
            {
                std::cerr << "ERROR: Insufficient statistics in response matrix ("
                          << nFilled << " entries). Skipping this bin." << std::endl;
                continue;
            }

            double priorEntries = hPrior->GetEntries();
            double priorIntegral = hPrior->Integral();

            if (verbose)
                std::cout << "MC truth statistics: " << priorEntries << " entries, integral = " << priorIntegral << std::endl;

            if (priorIntegral <= 0)
            {
                std::cerr << "ERROR: No MC truth events found for this bin! Check bin boundaries." << std::endl;
                delete hPrior;
                continue;
            }

            // === Added Diagnostic: response.Htruth vs hPrior consistency ===
            {
                TH1 *hRespTruth = response.Htruth();
                if (hRespTruth)
                {
                    double respTruthInt = hRespTruth->Integral();
                    double priorInt = hPrior->Integral();
                    double intRatio = (priorInt > 0 ? respTruthInt / priorInt : 0.0);
                    double maxRelDiff = 0.0;
                    int maxBin = 0;
                    int nBad = 0;
                    for (int b = 1; b <= hPrior->GetNbinsX(); ++b)
                    {
                        double a = hPrior->GetBinContent(b);
                        double c = hRespTruth->GetBinContent(b);
                        double rel = (a != 0 ? std::abs(a - c) / std::abs(a) : (c != 0 ? 1.0 : 0.0));
                        if (rel > maxRelDiff)
                        {
                            maxRelDiff = rel;
                            maxBin = b;
                        }
                        if (rel > 0.05)
                            ++nBad; // >5% discrepancy
                    }
                    if (verbose)
                    {
                        std::cout << "[CHECK] Htruth vs hPrior: integral(prior)=" << priorInt
                                  << " integral(Htruth)=" << respTruthInt
                                  << " ratio=" << intRatio
                                  << " maxRelDiff=" << maxRelDiff
                                  << " (bin=" << maxBin << ") nBins>5%=" << nBad << std::endl;
                        if (std::abs(1.0 - intRatio) > 0.02 || maxRelDiff > 0.1)
                        {
                            std::cout << "  WARNING: Response truth histogram differs notably from constructed prior (possible fill mismatch or selection divergence)." << std::endl;
                        }
                    }
                }
                else if (verbose)
                {
                    std::cout << "[CHECK] response.Htruth() returned null (no truth comparison)." << std::endl;
                }
            }

            // === Added Diagnostic: measured vs response.Hmeasured scaling (luminosity weighting mismatch) ===
            {
                TH1 *hRespMeasBase = response.Hmeasured();
                if (hRespMeasBase)
                {
                    double respMeasInt = hRespMeasBase->Integral();
                    double measInt = hMeasured->Integral();
                    double scale = (respMeasInt > 0 ? measInt / respMeasInt : 0.0);
                    if (verbose)
                    {
                        std::cout << "[CHECK] Measured vs response.Hmeasured: integral(Measured)=" << measInt
                                  << " integral(Hmeasured)=" << respMeasInt
                                  << " scale(Meas/Hmeas)=" << scale << std::endl;
                        if (scale < 0.5 || scale > 2.0)
                        {
                            std::cout << "  WARNING: Large global scale mismatch (>2x) between measured and response.Hmeasured; data may carry extra luminosity or per-event weights not in response." << std::endl;
                        }
                        int nbPrint = std::min(10, hMeasured->GetNbinsX());
                        std::cout << "  First " << nbPrint << " bin ratios (Meas / Hmeasured):" << std::endl;
                        for (int b = 1; b <= nbPrint; ++b)
                        {
                            double m = hMeasured->GetBinContent(b);
                            double r = hRespMeasBase->GetBinContent(b);
                            double ratio = (r > 0 ? m / r : 0.0);
                            std::cout << "    bin " << b << ": m=" << m << " r=" << r << " ratio=" << ratio << std::endl;
                        }
                    }
                    // Write ratio histogram
                    TH1D *hMeasRespRatio = (TH1D *)hMeasured->Clone(Form("hMeasRespRatio_%s_%s", jetPt.c_str(), yBinStr.Data()));
                    hMeasRespRatio->SetTitle("Measured / response.Hmeasured");
                    for (int b = 1; b <= hMeasRespRatio->GetNbinsX(); ++b)
                    {
                        double r = hRespMeasBase->GetBinContent(b);
                        double m = hMeasured->GetBinContent(b);
                        hMeasRespRatio->SetBinContent(b, (r > 0 ? m / r : 0.0));
                        hMeasRespRatio->SetBinError(b, 0.0);
                    }
                    hMeasRespRatio->Write();
                    // Plot Measured / Hmeasured ratio
                    {
                        TCanvas *cMeasResp = new TCanvas(Form("c_measRespRatio_%s_%s", jetPt.c_str(), yBinStr.Data()),
                                                         "Measured / Response", 800, 600);
                        cMeasResp->cd();
                        hMeasRespRatio->SetLineColor(kGreen + 2);
                        hMeasRespRatio->SetMarkerStyle(20);
                        hMeasRespRatio->SetMarkerColor(kGreen + 2);
                        hMeasRespRatio->GetXaxis()->SetTitle("z_{T}");
                        hMeasRespRatio->GetYaxis()->SetTitle("Measured / Hmeasured");
                        // Auto y-range (clip extremes but focus around 1)
                        double ymin = 1e9, ymax = -1e9;
                        for (int b = 1; b <= hMeasRespRatio->GetNbinsX(); ++b)
                        {
                            double v = hMeasRespRatio->GetBinContent(b);
                            if (v <= 0)
                                continue;
                            if (v < ymin)
                                ymin = v;
                            if (v > ymax)
                                ymax = v;
                        }
                        if (ymin > 0 && ymax > 0)
                        {
                            double span = ymax - ymin;
                            if (span < 0.2)
                                span = 0.2;
                            hMeasRespRatio->GetYaxis()->SetRangeUser(std::max(0.0, ymin - 0.1 * span), ymax + 0.1 * span);
                        }
                        hMeasRespRatio->Draw("EP");
                        TLine l1b(hMeasRespRatio->GetXaxis()->GetXmin(), 1.0, hMeasRespRatio->GetXaxis()->GetXmax(), 1.0);
                        l1b.SetLineStyle(2);
                        l1b.SetLineColor(kGray + 2);
                        l1b.Draw();
                        // Save ratio plot into diagnostics directory
                        TString pngMeasResp = Form("%s/measured_responseRatio_%s_%s.png", diagDir.c_str(), jetPt.c_str(), yBinStr.Data());
                        cMeasResp->SaveAs(pngMeasResp.Data());
                        cMeasResp->Write();
                        delete cMeasResp;
                    }
                }
            }

            // === RESPONSE MATRIX VALIDATION ===
            TH2 *hRespMatrix = response.Hresponse();
            if (hRespMatrix)
            {
                if (verbose)
                    std::cout << "\n=== Response Matrix Diagnostics ===" << std::endl;

                // Check response matrix properties
                double totalResp = hRespMatrix->Integral();
                if (verbose)
                    std::cout << "  Response matrix total: " << totalResp << std::endl;

                // Check row and column sums (efficiency and purity)
                std::vector<double> rowSums(nBins + 1, 0.0);
                std::vector<double> colSums(nBins + 1, 0.0);

                for (int i = 1; i <= nBins; ++i)
                {
                    for (int j = 1; j <= nBins; ++j)
                    {
                        double val = hRespMatrix->GetBinContent(i, j);
                        rowSums[i] += val; // reco bin i efficiency
                        colSums[j] += val; // truth bin j acceptance
                    }
                }

                // Print efficiency (fraction of truth events detected)
                if (verbose)
                    std::cout << "  Truth bin acceptance:" << std::endl;
                for (int j = 1; j <= nBins; ++j)
                {
                    double truthInBin = hPrior->GetBinContent(j);
                    double efficiency = (truthInBin > 0) ? colSums[j] / truthInBin : 0.0;
                    if (j <= 5 || efficiency < 0.1 || efficiency > 1.1)
                    { // Show first 5 bins and problematic ones
                        if (verbose)
                            std::cout << "    Truth bin " << j << ": " << colSums[j] << " detected / "
                                      << truthInBin << " truth = " << efficiency * 100 << "% eff" << std::endl;
                    }
                }

                // Check for empty rows/columns that could cause problems
                int emptyRows = 0, emptyCols = 0;
                for (int i = 1; i <= nBins; ++i)
                {
                    if (rowSums[i] <= 0)
                        emptyRows++;
                    if (colSums[i] <= 0)
                        emptyCols++;
                }
                if (emptyRows > 0 && verbose)
                    std::cout << "  WARNING: " << emptyRows << " empty reco bins in response!" << std::endl;
                if (emptyCols > 0 && verbose)
                    std::cout << "  WARNING: " << emptyCols << " empty truth bins in response!" << std::endl;
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

            double priorThresh = ((double)minCount > fracThresh * priorInclusive) ? (double)minCount : fracThresh * priorInclusive;
            if ((priorUnder + priorOver) > priorThresh)
            {
                std::cerr << "WARNING: Significant under/overflow in MC truth histogram for " << jetPt
                          << " bin " << yBinStr << ": under=" << priorUnder << ", over=" << priorOver
                          << ", total(incl)=" << priorInclusive << std::endl;
                std::cerr << "  Consider extending bin ranges or mapping under/overflow into boundary bins." << std::endl;
            }

            double measThresh = ((double)minCount > fracThresh * measuredInclusive) ? (double)minCount : fracThresh * measuredInclusive;
            if ((measuredUnder + measuredOver) > measThresh)
            {
                std::cerr << "WARNING: Significant under/overflow in measured histogram " << histName
                          << " for " << jetPt << " bin " << yBinStr << ": under=" << measuredUnder
                          << ", over=" << measuredOver << ", total(incl)=" << measuredInclusive << std::endl;
                std::cerr << "  Consider extending bin ranges, re-binning, or mapping under/overflow into boundary bins." << std::endl;
            }

            // Unfold using Bayesian method - each iteration builds on the previous
            fout->cd();

            // Collect unfolded histograms so we can plot all iterations together
            std::vector<TH1D *> unfoldedVec;
            // Store raw unfolded histograms (unnormalized) for gen-weight calculations and convergence
            std::vector<TH1D *> unfoldedSavedVec;
            // Collect refolded plotting clones for combined refolding comparison
            std::vector<TH1D *> refoldedPlotVec;
            // Collect closure-only diagnostics across iterations (overlay later)
            std::vector<TH1D *> unfoldTruthRatioVec;
            std::vector<TH1D *> unfoldTruthResidualVec;
            // Store covariance / correlation matrices per iteration for multi-panel summary
            std::vector<TH2D *> covMats;
            std::vector<TH2D *> corrMats;

            // Map under/overflow into boundary bins for measured and prior before further processing (c)
            // Always map under/overflow for all relevant histograms
            MapUnderOverflowToEdges(hPrior);
            MapUnderOverflowToEdges(hMeasured);
            MapUnderOverflowToEdges(hMissedTruth);
            MapUnderOverflowToEdges(hFakeReco);

            // Phase 1: normalization of measured histogram to response.Hmeasured integral
            TH1D *hMeasuredForUnfold = hMeasured;
            TH1D *hMeasuredForUnfoldOwned = nullptr;
            {
                TH1 *hRespMeasNow = response.Hmeasured();
                if (hRespMeasNow)
                {
                    double respInt = hRespMeasNow->Integral();
                    double measInt = hMeasured->Integral();
                    if (respInt > 0 && measInt > 0)
                    {
                        double scale = respInt / measInt;
                        hMeasuredForUnfoldOwned = (TH1D *)hMeasured->Clone(Form("hMeasuredForUnfold_%s_%s", jetPt.c_str(), yBinStr.Data()));
                        hMeasuredForUnfoldOwned->Scale(scale);
                        hMeasuredForUnfold = hMeasuredForUnfoldOwned;
                        std::cout << "[NORMALIZE] Measured scaled by " << scale << " (respInt=" << respInt << ", measInt=" << measInt << ")" << std::endl;
                    }
                }
            }

            // Perform unfolding with different iteration counts
            for (int iter = 1; iter <= nIter; ++iter)
            {
                if (verbose)
                    std::cout << "\n--- Bayesian Unfolding: " << iter << " iterations ---" << std::endl;

                // For iter>1, implement gen-level reweighting using previous unfolded result
                // Gen-level reweighting removed (always minimal mode)

                // RooUnfoldBayes will now use the possibly gen-weighted response (or the det-weighted one)
                RooUnfoldBayes unfold(&response, hMeasuredForUnfold, iter);

                // Optional: Set regularization parameter (smoothing) - this is by default set to iter
                // unfold.SetRegParm(2);

                TH1D *hUnfolded = (TH1D *)unfold.Hreco();
                if (!hUnfolded)
                {
                    std::cerr << "ERROR: Unfolding failed for " << iter << " iterations" << std::endl;
                    continue;
                }

                // Export covariance and correlation matrices
                {
                    TMatrixD cov = unfold.Ereco();
                    int nb = hUnfolded->GetNbinsX();
                    if (cov.GetNrows() == nb && cov.GetNcols() == nb)
                    {
                        TH2D *hCov = new TH2D(Form("cov_zT_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), iter),
                                              Form("Unfold covariance iter %d;bin;bin", iter), nb, 0.5, nb + 0.5, nb, 0.5, nb + 0.5);
                        TH2D *hCorr = new TH2D(Form("corr_zT_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), iter),
                                               Form("Unfold correlation iter %d;bin;bin", iter), nb, 0.5, nb + 0.5, nb, 0.5, nb + 0.5);
                        for (int i = 0; i < nb; ++i)
                        {
                            double di = cov(i, i) > 0 ? std::sqrt(cov(i, i)) : 0.0;
                            for (int j = 0; j < nb; ++j)
                            {
                                double v = cov(i, j);
                                hCov->SetBinContent(i + 1, j + 1, v);
                                double dj = cov(j, j) > 0 ? std::sqrt(cov(j, j)) : 0.0;
                                double corr = (di > 0 && dj > 0) ? v / (di * dj) : 0.0;
                                hCorr->SetBinContent(i + 1, j + 1, corr);
                            }
                        }
                        hCov->Write();
                        hCorr->Write();
                        covMats.push_back(hCov);
                        corrMats.push_back(hCorr);
                        if (verbose)
                            std::cout << "  - Wrote covariance and correlation matrices for iter " << iter << std::endl;
                    }
                    else if (verbose)
                    {
                        std::cout << "  WARNING: Covariance matrix dimension mismatch (" << cov.GetNrows() << "x" << cov.GetNcols() << ") vs bins=" << hUnfolded->GetNbinsX() << std::endl;
                    }
                }

                // Validate unfolded result
                double unfoldedIntegral = hUnfolded->Integral();
                double measuredIntegral = hMeasuredForUnfold->Integral();

                if (verbose)
                {
                    std::cout << "Results after " << iter << " iterations:" << std::endl;
                    std::cout << "  - Measured integral: " << measuredIntegral << std::endl;
                    std::cout << "  - Unfolded integral: " << unfoldedIntegral << std::endl;
                    std::cout << "  - Ratio (should be ~1): " << unfoldedIntegral / measuredIntegral << std::endl;
                }

                // Minimal-mode detailed diagnostics removed for brevity.

                TString iterName = Form("%s_%s_%s_iter%d", unfoldedName.c_str(), jetPt.c_str(), yBinStr.Data(), iter);
                TH1D *hUnfoldedSaved = (TH1D *)hUnfolded->Clone(iterName);
                // Write saved clone to output file
                hUnfoldedSaved->Write();
                // Keep raw (unnormalized) unfolded histograms for later gen-weight calculations
                unfoldedSavedVec.push_back(hUnfoldedSaved);
                if (verbose)
                    std::cout << "  - Stored raw unfolded histogram for iter " << iter << ": " << iterName << " (bins=" << hUnfoldedSaved->GetNbinsX() << ")" << std::endl;

                // --- Refolding test using ApplyToTruth (b) ---
                TH2 *hRespForRefold = response.Hresponse();
                if (hRespForRefold)
                {
                    TH1D *hRefolded = (TH1D *)response.ApplyToTruth(hUnfoldedSaved); // includes fakes
                    if (!hRefolded)
                    {
                        std::cerr << "ERROR: ApplyToTruth returned null" << std::endl;
                    }
                    else
                    {
                        hRefolded->SetDirectory(fout);
                        TString refoldName = Form("refolded_zT_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), iter);
                        hRefolded->SetName(refoldName);
                        hRefolded->Write();

                        Chi2Result refChi = ComputeChi2(hMeasuredForUnfold, hRefolded);
                        if (verbose)
                            std::cout << "  Refolding (ApplyToTruth) iter " << iter << ": chi2/ndf=" << refChi.redchi2 << std::endl;

                        if ((iter == 1 || iter == nIter) && verbose)
                        {
                            std::cout << "    Integral(measUsed)=" << hMeasuredForUnfold->Integral() << ", refold=" << hRefolded->Integral()
                                      << ", ratio=" << (hRefolded->Integral() / (hMeasuredForUnfold->Integral() > 0 ? hMeasuredForUnfold->Integral() : 1)) << std::endl;
                        }

                        // Residual histogram: (Measured - Refold)/Measured per bin (closure shape)
                        {
                            TH1D *hResidual = (TH1D *)hMeasuredForUnfold->Clone(Form("closureResidual_zT_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), iter));
                            hResidual->SetTitle(Form("Closure residual (M-R)/M iter %d", iter));
                            for (int b = 1; b <= hResidual->GetNbinsX(); ++b)
                            {
                                double m = hMeasuredForUnfold->GetBinContent(b);
                                double r = hRefolded->GetBinContent(b);
                                double em = hMeasuredForUnfold->GetBinError(b);
                                double val = (m > 0 ? (m - r) / m : 0.0);
                                // Approximate statistical uncertainty on residual from measured uncertainty only
                                double err = (m > 0 ? em / m : 0.0);
                                hResidual->SetBinContent(b, val);
                                hResidual->SetBinError(b, err);
                            }
                            // Set basic styling so when loaded later for overlay they have markers
                            hResidual->SetLineColor(kAzure + 1);
                            hResidual->SetMarkerStyle(20);
                            hResidual->Write();
                            // Quick summary of largest absolute residual (first iteration and final)
                            if ((iter == 1 || iter == nIter) && verbose)
                            {
                                double maxAbs = 0;
                                int maxBin = 0;
                                for (int b = 1; b <= hResidual->GetNbinsX(); ++b)
                                {
                                    double a = std::abs(hResidual->GetBinContent(b));
                                    if (a > maxAbs)
                                    {
                                        maxAbs = a;
                                        maxBin = b;
                                    }
                                }
                                std::cout << "    [RESIDUAL] Max |(M-R)/M|=" << maxAbs << " at bin=" << maxBin << std::endl;
                            }
                        }

                        // Additional closure-only diagnostics: unfolded vs truth comparison
                        if (isClosure)
                        {
                            // Ratio: unfolded / MC truth
                            TH1D *hUnfoldTruthRatio = (TH1D *)hUnfoldedSaved->Clone(Form("unfoldTruthRatio_zT_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), iter));
                            hUnfoldTruthRatio->SetTitle(Form("Unfold / Truth ratio iter %d", iter));
                            for (int b = 1; b <= hUnfoldTruthRatio->GetNbinsX(); ++b)
                            {
                                double tru = hPrior->GetBinContent(b);
                                double unf = hUnfoldedSaved->GetBinContent(b);
                                double etru = hPrior->GetBinError(b);
                                double eunf = hUnfoldedSaved->GetBinError(b);
                                double ratio = (tru > 0 ? unf / tru : 0.0);
                                // Propagate relative errors (assuming uncorrelated) if both >0
                                double err = 0.0;
                                if (tru > 0 && unf > 0)
                                {
                                    double relUnf = eunf / unf;
                                    double relTru = etru / tru;
                                    err = ratio * std::sqrt(relUnf * relUnf + relTru * relTru);
                                }
                                hUnfoldTruthRatio->SetBinContent(b, ratio);
                                hUnfoldTruthRatio->SetBinError(b, err);
                            }
                            hUnfoldTruthRatio->Write();
                            unfoldTruthRatioVec.push_back(hUnfoldTruthRatio);

                            // Residual: (Unfolded - Truth)/Truth
                            TH1D *hUnfoldTruthResidual = (TH1D *)hPrior->Clone(Form("unfoldTruthResidual_zT_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), iter));
                            hUnfoldTruthResidual->SetTitle(Form("(Truth - Unfold)/Truth iter %d", iter));
                            for (int b = 1; b <= hUnfoldTruthResidual->GetNbinsX(); ++b)
                            {
                                double tru = hPrior->GetBinContent(b);
                                double unf = hUnfoldedSaved->GetBinContent(b);
                                double etru = hPrior->GetBinError(b);
                                double eunf = hUnfoldedSaved->GetBinError(b);
                                double val = (tru > 0 ? (tru - unf) / tru : 0.0);
                                double err = 0.0;
                                if (tru > 0 && unf > 0)
                                {
                                    // (T-U)/T = 1 - U/T ; propagate error of U/T
                                    double ratio = unf / tru;
                                    double dU = eunf;
                                    double dT = etru;
                                    // variance of U/T ~ (dU/T)^2 + (U*dT/T^2)^2 (neglect covariance)
                                    double varRatio = (dU / tru) * (dU / tru) + (unf * dT / (tru * tru)) * (unf * dT / (tru * tru));
                                    err = std::sqrt(varRatio);
                                }
                                hUnfoldTruthResidual->SetBinContent(b, val);
                                hUnfoldTruthResidual->SetBinError(b, err);
                            }
                            hUnfoldTruthResidual->Write();
                            unfoldTruthResidualVec.push_back(hUnfoldTruthResidual);

                            // Pull distributions (aggregate and vs zT) : (Unfold - Truth)/sigma(Unfold)
                            TH1D *hPullVsZ = (TH1D *)hPrior->Clone(Form("pullVsZ_zT_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), iter));
                            hPullVsZ->SetTitle(Form("Pull vs z_{T} iter %d;z_{T};pull", iter));
                            TH1D *hPullDist = new TH1D(Form("pullDist_zT_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), iter),
                                                       Form("Pull distribution iter %d;pull;(bins)", iter), 60, -6, 6);
                            for (int b = 1; b <= hPrior->GetNbinsX(); ++b)
                            {
                                double tru = hPrior->GetBinContent(b);
                                double unf = hUnfoldedSaved->GetBinContent(b);
                                double eunf = hUnfoldedSaved->GetBinError(b);
                                if (eunf > 0)
                                {
                                    double pull = (unf - tru) / eunf; // sign convention
                                    hPullVsZ->SetBinContent(b, pull);
                                    hPullVsZ->SetBinError(b, 0.0);
                                    hPullDist->Fill(pull);
                                }
                                else
                                {
                                    hPullVsZ->SetBinContent(b, 0.0);
                                    hPullVsZ->SetBinError(b, 0.0);
                                }
                            }
                            hPullVsZ->Write();
                            hPullDist->Write();
                            // Plot pull vs zT
                            TCanvas *cPullVsZ = new TCanvas(Form("c_pullVsZ_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), iter), "Pull vs zT", 800, 600);
                            cPullVsZ->cd();
                            hPullVsZ->SetStats(0);
                            hPullVsZ->GetYaxis()->SetRangeUser(-6, 6);
                            hPullVsZ->Draw("HIST");
                            {
                                double xmin = hPullVsZ->GetXaxis()->GetXmin();
                                double xmax = hPullVsZ->GetXaxis()->GetXmax();
                                TLine l0(xmin, 0, xmax, 0);
                                l0.SetLineStyle(2);
                                l0.SetLineColor(kGray + 2);
                                l0.Draw();
                                TLine l1p(xmin, 1, xmax, 1);
                                l1p.SetLineStyle(3);
                                l1p.SetLineColor(kGray + 1);
                                l1p.Draw();
                                TLine l1m(xmin, -1, xmax, -1);
                                l1m.SetLineStyle(3);
                                l1m.SetLineColor(kGray + 1);
                                l1m.Draw();
                                TLine l2p(xmin, 2, xmax, 2);
                                l2p.SetLineStyle(3);
                                l2p.SetLineColor(kGray + 1);
                                l2p.Draw();
                                TLine l2m(xmin, -2, xmax, -2);
                                l2m.SetLineStyle(3);
                                l2m.SetLineColor(kGray + 1);
                                l2m.Draw();
                            }
                            TString pngPullVsZ = Form("%s/pullVsZ_%s_%s_iter%d.png", diagDir.c_str(), jetPt.c_str(), yBinStr.Data(), iter);
                            cPullVsZ->SaveAs(pngPullVsZ);
                            cPullVsZ->Write();
                            delete cPullVsZ;
                            // Plot pull distribution with Gaussian fit if enough bins
                            TCanvas *cPullDist = new TCanvas(Form("c_pullDist_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), iter), "Pull Distribution", 700, 600);
                            cPullDist->cd();
                            hPullDist->SetStats(0);
                            hPullDist->Draw("HIST E");
                            if (hPullDist->GetEntries() > 5)
                            {
                                TF1 *fG = new TF1(Form("fPullG_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), iter), "gaus", -3, 3);
                                hPullDist->Fit(fG, "QNR"); // quiet, no draw, store
                                fG->SetLineColor(kRed);
                                fG->Draw("SAME");
                                if (verbose)
                                    std::cout << "    Pull fit iter " << iter << ": mean=" << fG->GetParameter(1) << " sigma=" << fG->GetParameter(2) << std::endl;
                            }
                            TString pngPullDist = Form("%s/pullDist_%s_%s_iter%d.png", diagDir.c_str(), jetPt.c_str(), yBinStr.Data(), iter);
                            cPullDist->SaveAs(pngPullDist);
                            cPullDist->Write();
                            delete cPullDist;
                            if (verbose)
                            {
                                double mean = hPullDist->GetMean();
                                double rms = hPullDist->GetRMS();
                                std::cout << "    Pull stats iter " << iter << ": mean=" << mean << " rms=" << rms << std::endl;
                            }
                        }

                        // Normalized plotting clone
                        TH1D *hRefoldPlot = (TH1D *)hRefolded->Clone(Form("%s_plot", refoldName.Data()));
                        double ri_plot = hRefoldPlot->Integral();
                        if (ri_plot > 0)
                            hRefoldPlot->Scale(1.0 / ri_plot);
                        hRefoldPlot->SetEntries((Long64_t)hRefoldPlot->Integral());
                        refoldedPlotVec.push_back(hRefoldPlot);
                        // keep hRefolded (stored in file) but not in memory clones
                    }
                }
                else
                {
                    std::cerr << "Warning: cannot perform refolding test; response.Hresponse() is null." << std::endl;
                }
                // Create a separate plotting clone and normalize it (must be inside loop for scope)
                TH1D *hUnfoldedPlot = (TH1D *)hUnfoldedSaved->Clone(Form("%s_plot", iterName.Data()));
                double ui = hUnfoldedPlot->Integral();
                if (ui > 0)
                    hUnfoldedPlot->Scale(1.0 / ui);
                hUnfoldedPlot->SetEntries((Long64_t)hUnfoldedPlot->Integral());
                unfoldedVec.push_back(hUnfoldedPlot);
                std::cout << "  - Saved as: " << iterName << std::endl;
            }

            // === FINAL CLOSURE TEST SUMMARY ===
            if (verbose)
                std::cout << "\n=== CLOSURE TEST SUMMARY ===" << std::endl;
            if (!refoldedPlotVec.empty())
            {
                // Find best iteration based on chi2
                double bestChi2 = 1e9;
                int bestIter = 1;

                for (int iter = 1; iter <= nIter; ++iter)
                {
                    // Recalculate chi2 for each iteration (this could be stored from above)
                    TString refoldName = Form("refolded_zT_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), iter);
                    TH1D *hRefoldStored = (TH1D *)fout->Get(refoldName);
                    if (hRefoldStored)
                    {
                        double chi2 = 0.0;
                        int ndf = 0;
                        for (int i = 1; i <= nBins; ++i)
                        {
                            double m = hMeasuredForUnfold->GetBinContent(i);
                            double r = hRefoldStored->GetBinContent(i);
                            double err = hMeasuredForUnfold->GetBinError(i);
                            if (err > 0)
                            {
                                chi2 += (m - r) * (m - r) / (err * err);
                                ndf++;
                            }
                        }
                        double redChi2 = (ndf > 0) ? chi2 / ndf : 1e9;
                        if (verbose)
                            std::cout << "  Iteration " << iter << ": chi2/ndf = " << redChi2;
                        if (redChi2 < bestChi2)
                        {
                            bestChi2 = redChi2;
                            bestIter = iter;
                            if (verbose)
                                std::cout << " <- BEST";
                        }
                        if (verbose)
                            std::cout << std::endl;
                    }
                }

                if (verbose)
                    std::cout << "  Recommended iteration: " << bestIter << " (chi2/ndf = " << bestChi2 << ")" << std::endl;
                bestIterIndexForThisY = bestIter;

                if (bestChi2 > 2.0)
                {
                    if (verbose)
                    {
                        std::cout << "  WARNING: Poor closure (chi2/ndf > 2). Possible issues:" << std::endl;
                        std::cout << "    - Response matrix statistics too low" << std::endl;
                        std::cout << "    - Bin migration not well modeled" << std::endl;
                        std::cout << "    - Systematic differences between measured and MC" << std::endl;
                    }
                }
                else if (bestChi2 > 1.5)
                {
                    if (verbose)
                        std::cout << "  CAUTION: Moderate closure issues (chi2/ndf > 1.5)" << std::endl;
                }
                else
                {
                    if (verbose)
                        std::cout << "  GOOD: Closure test passed (chi2/ndf < 1.5)" << std::endl;
                }
                // Overlay of closure residuals across iterations
                {
                    TCanvas *cResAll = new TCanvas(Form("c_closureResidualAll_%s_%s", jetPt.c_str(), yBinStr.Data()),
                                                   "Closure Residuals All Iterations", 950, 750);
                    cResAll->cd();
                    std::vector<int> rcols = MakePalette(nIter, 0.7, 0.85);
                    TLegend legRes(0.15, 0.65, 0.5, 0.88);
                    legRes.SetBorderSize(0);
                    legRes.SetFillStyle(0);
                    legRes.SetTextSize(0.03);
                    bool firstDraw = true;
                    double xmin = 0, xmax = 0;
                    bool rangeSet = false;
                    for (int iter = 1; iter <= nIter; ++iter)
                    {
                        TH1D *hRes = (TH1D *)fout->Get(Form("closureResidual_zT_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), iter));
                        if (!hRes)
                            continue;
                        if (firstDraw)
                        {
                            hRes->SetLineColor(rcols[(iter - 1) % rcols.size()]);
                            hRes->SetMarkerColor(rcols[(iter - 1) % rcols.size()]);
                            hRes->SetMarkerStyle(20);
                            hRes->GetXaxis()->SetTitle("z_{T}");
                            hRes->GetYaxis()->SetTitle("(Measured - Refold)/Measured");
                            hRes->GetYaxis()->SetRangeUser(-1.0, 1.0);
                            hRes->Draw("EP");
                            xmin = hRes->GetXaxis()->GetXmin();
                            xmax = hRes->GetXaxis()->GetXmax();
                            rangeSet = true;
                            firstDraw = false;
                        }
                        else
                        {
                            hRes->SetLineColor(rcols[(iter - 1) % rcols.size()]);
                            hRes->SetMarkerColor(rcols[(iter - 1) % rcols.size()]);
                            hRes->SetMarkerStyle(20);
                            hRes->Draw("EPSAME");
                        }
                        legRes.AddEntry(hRes, Form("Iter %d", iter), "lep");
                    }
                    if (!firstDraw && rangeSet)
                    {
                        TLine l0(xmin, 0, xmax, 0);
                        l0.SetLineStyle(2);
                        l0.SetLineColor(kGray + 2);
                        l0.Draw();
                        TLine lPlus(xmin, 0.2, xmax, 0.2);
                        lPlus.SetLineStyle(3);
                        lPlus.SetLineColor(kGray + 1);
                        lPlus.Draw();
                        TLine lMinus(xmin, -0.2, xmax, -0.2);
                        lMinus.SetLineStyle(3);
                        lMinus.SetLineColor(kGray + 1);
                        lMinus.Draw();
                    }
                    legRes.Draw();
                    TLatex tex;
                    tex.SetNDC();
                    tex.SetTextSize(0.03);
                    tex.SetTextAlign(31);
                    tex.DrawLatex(0.95, 0.87, Form("Jet p_{T}: %s", jetPt.c_str()));
                    tex.DrawLatex(0.95, 0.83, Form("y-bin: %s", yBinStr.Data()));
                    TString pngResAll = Form("%s/closureResidual_all_%s_%s.png", jetOutDir.c_str(), jetPt.c_str(), yBinStr.Data());
                    cResAll->SaveAs(pngResAll.Data());
                    cResAll->Write();
                    delete cResAll;
                    // Capture best iteration residual and refold for rapidity multi-panel summary
                    if (bestIter > 0)
                    {
                        TH1D *hBestRes = (TH1D *)fout->Get(Form("closureResidual_zT_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), bestIter));
                        if (hBestRes)
                        {
                            TH1D *hClone = (TH1D *)hBestRes->Clone(Form("closureResidual_best_%s_%s", jetPt.c_str(), yBinStr.Data()));
                            // Store only y-bin tag; we'll retrieve all iteration residuals later
                            residualPerYBin.push_back({yBinStr});
                            delete hClone; // not needed
                        }
                        // Refold best iteration overlay data
                        TString bestRefName = Form("refolded_zT_%s_%s_iter%d", jetPt.c_str(), yBinStr.Data(), bestIter);
                        TH1D *hBestRef = (TH1D *)fout->Get(bestRefName);
                        if (hBestRef)
                        {
                            TH1D *hMClone = (TH1D *)hMeasuredForUnfold->Clone(Form("measured_norm_bestRef_%s_%s", jetPt.c_str(), yBinStr.Data()));
                            double imc = hMClone->Integral();
                            if (imc > 0)
                                hMClone->Scale(1.0 / imc);
                            hMClone->SetEntries((Long64_t)hMClone->Integral());
                            TH1D *hRefClone = (TH1D *)hBestRef->Clone(Form("refold_norm_bestRef_%s_%s", jetPt.c_str(), yBinStr.Data()));
                            double irc = hRefClone->Integral();
                            if (irc > 0)
                                hRefClone->Scale(1.0 / irc);
                            hRefClone->SetEntries((Long64_t)hRefClone->Integral());
                            refoldPerYBin.push_back({yBinStr, hMClone, hRefClone, bestIter});
                        }
                    }
                }
            }

            // Save original measured zT distribution
            TH1D *hMeasuredClone = (TH1D *)hMeasured->Clone(measuredOutName);
            hMeasuredClone->Write();
            // Create plotting clone for measured and normalize
            TH1D *hMeasuredPlot = (TH1D *)hMeasuredClone->Clone(Form("%s_plot", measuredOutName.Data()));
            double mi = hMeasuredPlot->Integral();
            if (mi > 0)
                hMeasuredPlot->Scale(1.0 / mi);
            // Keep Entries consistent after normalization
            hMeasuredPlot->SetEntries((Long64_t)hMeasuredPlot->Integral());

            // Save MC truth distribution for comparison
            TString truthOutName = Form("truth_zT_%s_%s", jetPt.c_str(), yBinStr.Data());
            hPrior->SetName(truthOutName);
            hPrior->Write();
            // Create plotting clone for truth and normalize
            TH1D *hPriorPlot = (TH1D *)hPrior->Clone(Form("%s_plot", truthOutName.Data()));
            double ti = hPriorPlot->Integral();
            if (ti > 0)
                hPriorPlot->Scale(1.0 / ti);
            // Keep Entries consistent after normalization
            hPriorPlot->SetEntries((Long64_t)hPriorPlot->Integral());

            // Save missed and fake distributions (even if empty) for diagnostics
            hMissedTruth->Write();
            hFakeReco->Write();
            // Normalized plotting clones
            TH1D *hMissedPlot = (TH1D *)hMissedTruth->Clone(Form("%s_plot", hMissedTruth->GetName()));
            double mi2 = hMissedPlot->Integral();
            if (mi2 > 0)
                hMissedPlot->Scale(1.0 / mi2);
            hMissedPlot->SetEntries((Long64_t)hMissedPlot->Integral());
            TH1D *hFakePlot = (TH1D *)hFakeReco->Clone(Form("%s_plot", hFakeReco->GetName()));
            double fi2 = hFakePlot->Integral();
            if (fi2 > 0)
                hFakePlot->Scale(1.0 / fi2);
            hFakePlot->SetEntries((Long64_t)hFakePlot->Integral());

            // Save response matrix object (edges already mapped earlier)
            TString responseOutName = Form("response_%s_%s", jetPt.c_str(), yBinStr.Data());
            response.Write(responseOutName);

            // Plot and save the response matrix (TH2) into the output directory as PNG
            TH2 *hResp = response.Hresponse();
            if (hResp)
            {
                // Clone to give a stable name and write to file
                TH2 *hRespClone = (TH2 *)hResp->Clone(Form("hResponse_%s_%s", jetPt.c_str(), yBinStr.Data()));
                hRespClone->SetTitle(Form("Response matrix %s bin %s", jetPt.c_str(), yBinStr.Data()));
                hRespClone->Write();

                // Draw on a TCanvas and save as PNG so the image is rendered correctly
                TCanvas *cResp = new TCanvas(Form("c_resp_%s_%s", jetPt.c_str(), yBinStr.Data()), "Response Matrix", 800, 600);
                cResp->cd();
                hRespClone->Draw("COLZ");
                TString pngName = Form("%s/response_%s_%s.png", jetOutDir.c_str(), jetPt.c_str(), yBinStr.Data());
                cResp->SaveAs(pngName.Data());
                delete cResp;
                delete hRespClone;
                std::cout << "  - Saved response image: " << pngName << std::endl;
                // Record y bin for later multi-panel summary
                responsePerYBin.push_back(yBinStr);
            }
            else
            {
                std::cerr << "Warning: response.Hresponse() returned null for " << histName << std::endl;
            }

            // Create a comparison plot: measured, MC truth, and all unfolded iterations
            {
                TCanvas *cComp = new TCanvas(Form("c_comp_%s_%s", jetPt.c_str(), yBinStr.Data()), "zT comparison", 900, 700);
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
                if (hPriorPlot->GetMaximum() > maxVal)
                    maxVal = hPriorPlot->GetMaximum();
                for (size_t ii = 0; ii < unfoldedVec.size(); ++ii)
                {
                    if (unfoldedVec[ii]->GetMaximum() > maxVal)
                        maxVal = unfoldedVec[ii]->GetMaximum();
                }
                hMeasuredPlot->SetMaximum(maxVal * 1.2);

                // Draw base measured and truth plots
                hMeasuredPlot->GetXaxis()->SetTitle("z_{T}");
                hMeasuredPlot->GetYaxis()->SetTitle("dN/dz_{T} (norm)");
                hMeasuredPlot->Draw("E");
                hPriorPlot->Draw("HISTSAME");

                // Generate a smooth palette for unfolded iterations
                std::vector<int> ucols = MakePalette((int)unfoldedVec.size());
                for (size_t ii = 0; ii < unfoldedVec.size(); ++ii)
                {
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
                leg.SetFillStyle(0);
                leg.SetFillColor(0);
                leg.AddEntry(hMeasuredPlot, "Measured", "lep");
                leg.AddEntry(hPriorPlot, "MC truth", "l");
                for (size_t ii = 0; ii < unfoldedVec.size(); ++ii)
                {
                    TString lbl = Form("Unfold iter %zu", ii + 1);
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
                if (verbose)
                    std::cout << "  - Saved comparison image: " << cmpName << std::endl;

                // Store best iteration comparison (measured, truth, unfolded best) for rapidity multi-panel
                if (bestIterIndexForThisY > 0 && (size_t)(bestIterIndexForThisY - 1) < unfoldedVec.size())
                {
                    TH1D *cMeas = (TH1D *)hMeasuredPlot->Clone(Form("compareBest_meas_%s_%s", jetPt.c_str(), yBinStr.Data()));
                    TH1D *cTruth = (TH1D *)hPriorPlot->Clone(Form("compareBest_truth_%s_%s", jetPt.c_str(), yBinStr.Data()));
                    TH1D *cUnf = (TH1D *)unfoldedVec[bestIterIndexForThisY - 1]->Clone(Form("compareBest_unfolded_%s_%s", jetPt.c_str(), yBinStr.Data()));
                    comparePerYBin.push_back({yBinStr, cMeas, cTruth, cUnf, bestIterIndexForThisY});
                }

                // Cleanup plotting clones from memory (saved clones remain in the ROOT file)
                for (auto h : unfoldedVec)
                    delete h;
                delete hMeasuredPlot;
                delete hPriorPlot;
                // Produce an auxiliary canvas comparing measured, truth, missed, fake (normalized shapes)
                {
                    TCanvas *cMF = new TCanvas(Form("c_meas_truth_miss_fake_%s_%s", jetPt.c_str(), yBinStr.Data()),
                                               "Measured / Truth / Missed / Fake", 900, 700);
                    cMF->SetRightMargin(0.02);
                    // Style
                    hMeasuredClone->SetLineColor(kBlack);
                    hMeasuredClone->SetMarkerStyle(20);
                    hMeasuredClone->SetMarkerColor(kBlack);
                    hMeasuredClone->SetLineWidth(2);
                    hPrior->SetLineColor(kRed);
                    hPrior->SetLineStyle(2);
                    hPrior->SetLineWidth(2);
                    hMissedPlot->SetLineColor(kBlue + 1);
                    hMissedPlot->SetLineStyle(3);
                    hMissedPlot->SetLineWidth(2);
                    hFakePlot->SetLineColor(kMagenta + 2);
                    hFakePlot->SetLineStyle(4);
                    hFakePlot->SetLineWidth(2);
                    // Determine max among normalized shapes (reuse normalized clones except measuredClone/truth need normalized versions)
                    TH1D *hMeasNormTmp = (TH1D *)hMeasuredClone->Clone("_tmp_measNorm");
                    double imn = hMeasNormTmp->Integral();
                    if (imn > 0)
                        hMeasNormTmp->Scale(1.0 / imn);
                    TH1D *hTruthNormTmp = (TH1D *)hPrior->Clone("_tmp_truthNorm");
                    double itn = hTruthNormTmp->Integral();
                    if (itn > 0)
                        hTruthNormTmp->Scale(1.0 / itn);
                    double ymaxMF = std::max({hMeasNormTmp->GetMaximum(), hTruthNormTmp->GetMaximum(), hMissedPlot->GetMaximum(), hFakePlot->GetMaximum()});
                    hMeasNormTmp->SetMaximum(ymaxMF * 1.25);
                    hMeasNormTmp->GetXaxis()->SetTitle("z_{T}");
                    hMeasNormTmp->GetYaxis()->SetTitle("dN/dz_{T} (norm)");
                    hMeasNormTmp->Draw("HIST");
                    hTruthNormTmp->Draw("HISTSAME");
                    if (hMissedPlot->Integral() > 0)
                        hMissedPlot->Draw("HISTSAME");
                    if (hFakePlot->Integral() > 0)
                        hFakePlot->Draw("HISTSAME");
                    TLegend legMF(0.14, 0.62, 0.48, 0.88);
                    legMF.SetBorderSize(0);
                    legMF.SetFillStyle(0);
                    legMF.SetTextSize(0.03);
                    legMF.AddEntry(hMeasNormTmp, "Measured", "l");
                    legMF.AddEntry(hTruthNormTmp, "MC truth", "l");
                    legMF.AddEntry(hMissedPlot, "Missed (truth only)", "l");
                    legMF.AddEntry(hFakePlot, "Fake (reco only)", "l");
                    legMF.Draw();
                    TLatex texMF;
                    texMF.SetNDC();
                    texMF.SetTextSize(0.03);
                    texMF.SetTextAlign(31);
                    texMF.DrawLatex(0.95, 0.85, Form("Jet p_{T}: %s", jetPt.c_str()));
                    texMF.DrawLatex(0.95, 0.81, Form("y-bin: %s", yBinStr.Data()));
                    TString pngMF = Form("%s/meas_truth_miss_fake_%s_%s.png", jetOutDir.c_str(), jetPt.c_str(), yBinStr.Data());
                    cMF->SaveAs(pngMF.Data());
                    cMF->Write();
                    if (verbose)
                        std::cout << "  - Saved Measured/Truth/Missed/Fake comparison: " << pngMF << std::endl;
                    // Store normalized clones for multi-panel; create persistent normalized versions
                    TH1D *hMeasStore = (TH1D *)hMeasuredClone->Clone(Form("mfmf_meas_%s_%s", jetPt.c_str(), yBinStr.Data()));
                    double imS = hMeasStore->Integral();
                    if (imS > 0)
                        hMeasStore->Scale(1.0 / imS);
                    TH1D *hTruthStore = (TH1D *)hPrior->Clone(Form("mfmf_truth_%s_%s", jetPt.c_str(), yBinStr.Data()));
                    double itS = hTruthStore->Integral();
                    if (itS > 0)
                        hTruthStore->Scale(1.0 / itS);
                    TH1D *hMissStore = (TH1D *)hMissedPlot->Clone(Form("mfmf_missed_%s_%s", jetPt.c_str(), yBinStr.Data()));
                    TH1D *hFakeStore = (TH1D *)hFakePlot->Clone(Form("mfmf_fake_%s_%s", jetPt.c_str(), yBinStr.Data()));
                    mfmfPerYBin.push_back({yBinStr, hMeasStore, hTruthStore, hMissStore, hFakeStore});
                    delete hMeasNormTmp;
                    delete hTruthNormTmp;
                    delete cMF;
                }
                delete hMissedPlot;
                delete hFakePlot; // normalized clones
            }

            // Combined refolding comparison: overlay all refolded iterations (single canvas)
            if (!refoldedPlotVec.empty())
            {
                TCanvas *cRefAll = new TCanvas(Form("c_refold_all_%s_%s", jetPt.c_str(), yBinStr.Data()),
                                               Form("Refolding all iters %s %s", jetPt.c_str(), yBinStr.Data()), 900, 700);
                cRefAll->SetRightMargin(0.02);
                cRefAll->cd();

                // Top overlay: normalized measured and all refolded iterations
                TH1D *hMeasuredNorm2 = (TH1D *)hMeasured->Clone(Form("%s_norm_allref", measuredOutName.Data()));
                double mi2 = hMeasuredNorm2->Integral();
                if (mi2 > 0)
                    hMeasuredNorm2->Scale(1.0 / mi2);
                hMeasuredNorm2->SetLineColor(kBlack);
                hMeasuredNorm2->SetMarkerStyle(20);
                double maxVal2 = hMeasuredNorm2->GetMaximum();
                hMeasuredNorm2->GetXaxis()->SetTitle("z_{T}");
                hMeasuredNorm2->GetYaxis()->SetTitle("dN/dz_{T} (norm)");
                hMeasuredNorm2->Draw("E");

                std::vector<int> rcols = MakePalette((int)refoldedPlotVec.size());
                for (size_t ii = 0; ii < refoldedPlotVec.size(); ++ii)
                {
                    TH1D *hr = refoldedPlotVec[ii];
                    int col = rcols[ii % rcols.size()];
                    hr->SetLineColor(col);
                    hr->SetLineWidth(2);
                    hr->SetLineStyle(1);
                    if (hr->GetMaximum() > maxVal2)
                        maxVal2 = hr->GetMaximum();
                    hr->Draw("HISTSAME");
                }
                hMeasuredNorm2->SetMaximum(maxVal2 * 1.2);

                // Legend
                TLegend legRef(0.15, 0.65, 0.6, 0.88);
                legRef.SetBorderSize(0);
                legRef.SetTextSize(0.030);
                legRef.SetMargin(0.12);
                legRef.SetFillStyle(0);
                legRef.SetFillColor(0);
                legRef.AddEntry(hMeasuredNorm2, "Measured", "lep");
                for (size_t ii = 0; ii < refoldedPlotVec.size(); ++ii)
                {
                    TString lbl = Form("Refold iter %zu", ii + 1);
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
                if (verbose)
                    std::cout << "  - Saved combined refolding comparison image: " << allPng << std::endl;

                // (Removed iteration-based multi-panel refold; replaced by rapidity-bin version after loop)

                // Cleanup
                delete hMeasuredNorm2;
                delete cRefAll;
                for (auto h : refoldedPlotVec)
                    delete h;
            }

            // Multi-panel covariance & correlation summary plots across iterations
            if (!covMats.empty())
            {
                int nPanels = (int)covMats.size();
                int nCols = (nPanels <= 3) ? nPanels : (int)std::ceil(std::sqrt(nPanels));
                if (nCols < 1)
                    nCols = 1;
                if (nCols > 6)
                    nCols = 6; // limit
                int nRows = (int)std::ceil((double)nPanels / nCols);
                int w = std::min(300 * nCols, 2400);
                int h = std::min(300 * nRows, 1800);
                // Covariance range
                double cmin = 1e300, cmax = -1e300;
                for (auto *m : covMats)
                {
                    if (!m)
                        continue;
                    int nx = m->GetNbinsX();
                    int ny = m->GetNbinsY();
                    for (int ix = 1; ix <= nx; ++ix)
                        for (int iy = 1; iy <= ny; ++iy)
                        {
                            double v = m->GetBinContent(ix, iy);
                            if (v < cmin)
                                cmin = v;
                            if (v > cmax)
                                cmax = v;
                        }
                }
                if (cmin == cmax)
                {
                    cmin = 0;
                }
                // Covariance canvas
                TCanvas *cCovAll = new TCanvas(Form("c_cov_all_%s_%s", jetPt.c_str(), yBinStr.Data()), "Covariance All", w, h);
                cCovAll->Divide(nCols, nRows, 0.001, 0.001);
                for (int i = 0; i < nPanels; ++i)
                {
                    cCovAll->cd(i + 1);
                    if (!covMats[i])
                        continue;
                    covMats[i]->SetStats(0);
                    covMats[i]->GetZaxis()->SetRangeUser(cmin, cmax);
                    covMats[i]->Draw("COLZ");
                    TLatex tl;
                    tl.SetNDC();
                    tl.SetTextSize(0.04);
                    tl.DrawLatex(0.02, 0.94, Form("Iter %d", i + 1));
                }
                TString pngCovAll = Form("%s/covariance_all_%s_%s.png", diagDir.c_str(), jetPt.c_str(), yBinStr.Data());
                cCovAll->SaveAs(pngCovAll);
                cCovAll->Write();
                if (verbose)
                    std::cout << "  - Saved multi-panel covariance: " << pngCovAll << std::endl;
                delete cCovAll;
            }
            if (!corrMats.empty())
            {
                int nPanels = (int)corrMats.size();
                int nCols = (nPanels <= 3) ? nPanels : (int)std::ceil(std::sqrt(nPanels));
                if (nCols < 1)
                    nCols = 1;
                if (nCols > 6)
                    nCols = 6;
                int nRows = (int)std::ceil((double)nPanels / nCols);
                int w = std::min(300 * nCols, 2400);
                int h = std::min(300 * nRows, 1800);
                TCanvas *cCorrAll = new TCanvas(Form("c_corr_all_%s_%s", jetPt.c_str(), yBinStr.Data()), "Correlation All", w, h);
                cCorrAll->Divide(nCols, nRows, 0.001, 0.001);
                for (int i = 0; i < nPanels; ++i)
                {
                    cCorrAll->cd(i + 1);
                    if (!corrMats[i])
                        continue;
                    corrMats[i]->SetStats(0);
                    corrMats[i]->GetZaxis()->SetRangeUser(-1, 1);
                    corrMats[i]->Draw("COLZ");
                    TLatex tl;
                    tl.SetNDC();
                    tl.SetTextSize(0.04);
                    tl.DrawLatex(0.02, 0.94, Form("Iter %d", i + 1));
                }
                TString pngCorrAll = Form("%s/correlation_all_%s_%s.png", diagDir.c_str(), jetPt.c_str(), yBinStr.Data());
                cCorrAll->SaveAs(pngCorrAll);
                cCorrAll->Write();
                if (verbose)
                    std::cout << "  - Saved multi-panel correlation: " << pngCorrAll << std::endl;
                delete cCorrAll;
            }

            // Overlay all Unfold/Truth ratios across iterations (closure mode)
            if (isClosure && !unfoldTruthRatioVec.empty())
            {
                TCanvas *cRatioAll = new TCanvas(Form("c_unfoldTruthRatio_all_%s_%s", jetPt.c_str(), yBinStr.Data()),
                                                 "Unfold/Truth Ratio All Iterations", 900, 700);
                cRatioAll->SetRightMargin(0.02);
                cRatioAll->cd();
                std::vector<int> cols = MakePalette((int)unfoldTruthRatioVec.size(), 0.75, 0.85);
                bool first = true;
                double xmin = 0, xmax = 0;
                double ymin = 0.0, ymax = 10.5; // fixed range
                TLegend legR(0.15, 0.65, 0.5, 0.88);
                legR.SetBorderSize(0);
                legR.SetFillStyle(0);
                legR.SetTextSize(0.03);
                for (size_t i = 0; i < unfoldTruthRatioVec.size(); ++i)
                {
                    TH1D *h = unfoldTruthRatioVec[i];
                    if (!h)
                        continue;
                    int col = cols[i % cols.size()];
                    h->SetLineColor(col);
                    h->SetMarkerColor(col);
                    h->SetMarkerStyle(20);
                    if (first)
                    {
                        h->GetXaxis()->SetTitle("z_{T}");
                        h->GetYaxis()->SetTitle("Unfold / Truth");
                        h->GetYaxis()->SetRangeUser(ymin, ymax);
                        h->Draw("EP");
                        xmin = h->GetXaxis()->GetXmin();
                        xmax = h->GetXaxis()->GetXmax();
                        first = false;
                    }
                    else
                    {
                        h->Draw("EPSAME");
                    }
                    legR.AddEntry(h, Form("Iter %zu", i + 1), "lep");
                }
                if (!first)
                {
                    TLine l1(xmin, 1.0, xmax, 1.0);
                    l1.SetLineStyle(2);
                    l1.SetLineColor(kGray + 2);
                    l1.Draw();
                    TLine lUp(xmin, 1.2, xmax, 1.2);
                    lUp.SetLineStyle(3);
                    lUp.SetLineColor(kGray + 1);
                    lUp.Draw();
                    TLine lDn(xmin, 0.8, xmax, 0.8);
                    lDn.SetLineStyle(3);
                    lDn.SetLineColor(kGray + 1);
                    lDn.Draw();
                }
                legR.Draw();
                TLatex tex;
                tex.SetNDC();
                tex.SetTextSize(0.03);
                tex.SetTextAlign(31);
                tex.DrawLatex(0.95, 0.87, Form("Jet p_{T}: %s", jetPt.c_str()));
                tex.DrawLatex(0.95, 0.83, Form("y-bin: %s", yBinStr.Data()));
                TString pngRatioAll = Form("%s/unfoldTruthRatio_all_%s_%s.png", jetOutDir.c_str(), jetPt.c_str(), yBinStr.Data());
                cRatioAll->SaveAs(pngRatioAll.Data());
                cRatioAll->Write();
                delete cRatioAll;
            }

            // Overlay all (Truth-Unfold)/Truth residuals across iterations (closure mode)
            if (isClosure && !unfoldTruthResidualVec.empty())
            {
                TCanvas *cResTAll = new TCanvas(Form("c_unfoldTruthResidual_all_%s_%s", jetPt.c_str(), yBinStr.Data()),
                                                "(Truth-Unfold)/Truth All Iterations", 900, 700);
                cResTAll->SetRightMargin(0.02);
                cResTAll->cd();
                std::vector<int> cols = MakePalette((int)unfoldTruthResidualVec.size(), 0.75, 0.85);
                bool first = true;
                double xmin = 0, xmax = 0;
                double ymin = -10.0, ymax = 2.0;
                TLegend legT(0.15, 0.65, 0.5, 0.88);
                legT.SetBorderSize(0);
                legT.SetFillStyle(0);
                legT.SetTextSize(0.03);
                for (size_t i = 0; i < unfoldTruthResidualVec.size(); ++i)
                {
                    TH1D *h = unfoldTruthResidualVec[i];
                    if (!h)
                        continue;
                    int col = cols[i % cols.size()];
                    h->SetLineColor(col);
                    h->SetMarkerColor(col);
                    h->SetMarkerStyle(20);
                    if (first)
                    {
                        h->GetXaxis()->SetTitle("z_{T}");
                        h->GetYaxis()->SetTitle("(Truth-Unfold)/Truth");
                        h->GetYaxis()->SetRangeUser(ymin, ymax);
                        h->Draw("EP");
                        xmin = h->GetXaxis()->GetXmin();
                        xmax = h->GetXaxis()->GetXmax();
                        first = false;
                    }
                    else
                    {
                        h->Draw("EPSAME");
                    }
                    legT.AddEntry(h, Form("Iter %zu", i + 1), "lep");
                }
                if (!first)
                {
                    TLine l0(xmin, 0.0, xmax, 0.0);
                    l0.SetLineStyle(2);
                    l0.SetLineColor(kGray + 2);
                    l0.Draw();
                    TLine lUp(xmin, 0.2, xmax, 0.2);
                    lUp.SetLineStyle(3);
                    lUp.SetLineColor(kGray + 1);
                    lUp.Draw();
                    TLine lDn(xmin, -0.2, xmax, -0.2);
                    lDn.SetLineStyle(3);
                    lDn.SetLineColor(kGray + 1);
                    lDn.Draw();
                }
                legT.Draw();
                TLatex tex;
                tex.SetNDC();
                tex.SetTextSize(0.03);
                tex.SetTextAlign(31);
                tex.DrawLatex(0.95, 0.87, Form("Jet p_{T}: %s", jetPt.c_str()));
                tex.DrawLatex(0.95, 0.83, Form("y-bin: %s", yBinStr.Data()));
                TString pngResAll = Form("%s/unfoldTruthResidual_all_%s_%s.png", jetOutDir.c_str(), jetPt.c_str(), yBinStr.Data());
                cResTAll->SaveAs(pngResAll.Data());
                cResTAll->Write();
                delete cResTAll;
            }

            // cleanup closure-only vectors (their histograms are written; we keep objects for potential later delete if needed)
            unfoldTruthRatioVec.clear();
            unfoldTruthResidualVec.clear();

            if (verbose)
            {
                std::cout << "\nSaved outputs:" << std::endl;
                std::cout << "  - Measured zT: " << measuredOutName << std::endl;
                std::cout << "  - MC truth zT: " << truthOutName << std::endl;
                std::cout << "  - Response matrix: " << responseOutName << std::endl;
            }

            // Cleanup
            if (hMeasuredForUnfoldOwned)
            {
                delete hMeasuredForUnfoldOwned;
                hMeasuredForUnfoldOwned = nullptr;
            }
            delete hPrior;
        }
        // After finishing all rapidity bins for this jetPt, build multi-panel residual vs y bins
        if (!residualPerYBin.empty())
        {
            int nPanels = (int)residualPerYBin.size();
            int nCols = (nPanels <= 3) ? nPanels : (int)std::ceil(std::sqrt(nPanels));
            if (nCols < 1)
                nCols = 1;
            if (nCols > 6)
                nCols = 6;
            int nRows = (int)std::ceil((double)nPanels / nCols);
            int w = std::min(340 * nCols, 2400);
            int h = std::min(320 * nRows, 1800);
            TCanvas *cResYMulti = new TCanvas(Form("c_closureResidual_multi_%s_allY", jetPt.c_str()), "Closure Residuals (all iters) across y bins", w, h);
            cResYMulti->Divide(nCols, nRows, 0.001, 0.001);
            for (int i = 0; i < nPanels; ++i)
            {
                cResYMulti->cd(i + 1);
                gPad->SetTickx();
                gPad->SetTicky();
                gPad->SetRightMargin(0.02);
                gPad->SetLeftMargin(0.12);
                gPad->SetTopMargin(0.02);
                gPad->SetBottomMargin(0.14);
                TString yStr = residualPerYBin[i].yBinStr;
                std::vector<int> cols = MakePalette(nIter, 0.70, 0.85);
                TLegend legPad(0.12, 0.68, 0.55, 0.90);
                legPad.SetBorderSize(0);
                legPad.SetFillStyle(0);
                legPad.SetTextSize(0.035);
                bool first = true;
                double xmin = 0, xmax = 0;
                bool haveRange = false;
                for (int iter = 1; iter <= nIter; ++iter)
                {
                    TH1D *hRes = (TH1D *)fout->Get(Form("closureResidual_zT_%s_%s_iter%d", jetPt.c_str(), yStr.Data(), iter));
                    if (!hRes)
                        continue;
                    hRes->SetLineColor(cols[(iter - 1) % cols.size()]);
                    hRes->SetMarkerColor(cols[(iter - 1) % cols.size()]);
                    hRes->SetMarkerStyle(20);
                    if (first)
                    {
                        hRes->GetYaxis()->SetRangeUser(-1.0, 1.0);
                        hRes->GetXaxis()->SetTitle("z_{T}");
                        hRes->GetYaxis()->SetTitle("(Measured - Refold)/Measured");
                        hRes->SetTitle("");
                        hRes->Draw("EP");
                        xmin = hRes->GetXaxis()->GetXmin();
                        xmax = hRes->GetXaxis()->GetXmax();
                        haveRange = true;
                        first = false;
                    }
                    else
                    {
                        hRes->Draw("EPSAME");
                    }
                    legPad.AddEntry(hRes, Form("Iter %d", iter), "lep");
                }
                // Draw darker, thicker guide lines so they remain visible over markers/lines
                TLine l0(0, 0.0, 1, 0.0);
                l0.SetLineStyle(2);
                l0.SetLineColor(kBlack);
                l0.SetLineWidth(2);
                l0.Draw();
                TLine lUp(0, 0.2, 1, 0.2);
                lUp.SetLineStyle(3);
                lUp.SetLineColor(kBlack);
                lUp.SetLineWidth(1);
                lUp.Draw();
                TLine lDn(0, -0.2, 1, -0.2);
                lDn.SetLineStyle(3);
                lDn.SetLineColor(kBlack);
                lDn.SetLineWidth(1);
                lDn.Draw();
                gPad->Modified();
                gPad->Update();
                legPad.Draw();
                TLatex tl;
                tl.SetNDC();
                tl.SetTextSize(0.045);
                if (i == 0)
                    tl.DrawLatex(0.20, 0.90, "#font[12]{LHCb} in-progress");
                if (i == 0)
                    tl.DrawLatex(0.20, 0.85, Form("Jet p_{T}: %s", jetPt.c_str()));
                tl.DrawLatex(0.20, (i == 0) ? 0.80 : 0.90, Form("y: %s", yStr.Data()));
            }
            TString pngResY = Form("%s/closureResidual_multi_%s_allY.png", multiDir.c_str(), jetPt.c_str());
            cResYMulti->SaveAs(pngResY);
            TString pdfResY = Form("%s/closureResidual_multi_%s_allY.pdf", multiDir.c_str(), jetPt.c_str());
            cResYMulti->SaveAs(pdfResY);
            cResYMulti->Write();
            if (verbose)
                std::cout << "  - Saved multi-panel closure residual (all iterations) across rapidity bins: " << pngResY << std::endl;
            delete cResYMulti;
        }
        // Build refold comparison multi-panel across rapidity bins (best iteration per y)
        if (!refoldPerYBin.empty())
        {
            int nPanels = (int)refoldPerYBin.size();
            int nCols = (nPanels <= 3) ? nPanels : (int)std::ceil(std::sqrt(nPanels));
            if (nCols < 1)
                nCols = 1;
            if (nCols > 6)
                nCols = 6;
            int nRows = (int)std::ceil((double)nPanels / nCols);
            int w = std::min(340 * nCols, 2400);
            int h = std::min(320 * nRows, 1800);
            TCanvas *cRefYMulti = new TCanvas(Form("c_refold_multi_%s_allY", jetPt.c_str()), "Refold (best iter) across y bins", w, h);
            cRefYMulti->Divide(nCols, nRows, 0.001, 0.001);
            for (int i = 0; i < nPanels; ++i)
            {
                cRefYMulti->cd(i + 1);
                gPad->SetTickx();
                gPad->SetTicky();
                gPad->SetRightMargin(0.02);
                gPad->SetLeftMargin(0.12);
                gPad->SetTopMargin(0.02);
                gPad->SetBottomMargin(0.14);
                TH1D *hM = refoldPerYBin[i].meas;
                TH1D *hR = refoldPerYBin[i].refold;
                if (!hM || !hR)
                    continue;
                // Determine max for consistent pad scaling
                double maxv = std::max(hM->GetMaximum(), hR->GetMaximum());
                hM->SetMaximum(maxv * 1.25);
                hM->SetLineColor(kBlack);
                hM->SetMarkerStyle(20);
                hM->GetXaxis()->SetTitle("z_{T}");
                hM->GetYaxis()->SetTitle("dN/dz_{T} (norm)");
                hM->SetTitle("");
                hM->Draw("E");
                hR->SetLineColor(kRed + 1);
                hR->SetLineWidth(2);
                hR->Draw("HISTSAME");
                TLatex tl;
                tl.SetNDC();
                tl.SetTextSize(0.045);
                if (i == 0)
                    tl.DrawLatex(0.20, 0.90, "#font[12]{LHCb} in-progress");
                if (i == 0)
                    tl.DrawLatex(0.20, 0.85, Form("Jet p_{T}: %s", jetPt.c_str()));
                tl.DrawLatex(0.20, (i == 0) ? 0.80 : 0.90, Form("y: %s (iter %d)", refoldPerYBin[i].yBinStr.Data(), refoldPerYBin[i].iter));
                if (i == 0)
                {
                    TLegend *leg = new TLegend(0.55, 0.70, 0.90, 0.88);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    leg->SetTextSize(0.04);
                    leg->AddEntry(hM, "Measured", "lep");
                    leg->AddEntry(hR, "Refold (best)", "l");
                    leg->Draw();
                }
            }
            TString pngRefY = Form("%s/refold_compare_multi_%s_allY.png", multiDir.c_str(), jetPt.c_str());
            cRefYMulti->SaveAs(pngRefY);
            cRefYMulti->Write();
            if (verbose)
                std::cout << "  - Saved multi-panel refold comparison across rapidity bins: " << pngRefY << std::endl;
            delete cRefYMulti;
            for (auto &e : refoldPerYBin)
            {
                delete e.meas;
                delete e.refold;
            }
        }
        // Build comparison multi-panel across rapidity bins (best iteration per y)
        if (!comparePerYBin.empty())
        {
            int nPanels = (int)comparePerYBin.size();
            int nCols = (nPanels <= 3) ? nPanels : (int)std::ceil(std::sqrt(nPanels));
            if (nCols < 1)
                nCols = 1;
            if (nCols > 6)
                nCols = 6;
            int nRows = (int)std::ceil((double)nPanels / nCols);
            int w = std::min(340 * nCols, 2400);
            int h = std::min(320 * nRows, 1800);
            TCanvas *cCompYMulti = new TCanvas(Form("c_compare_zT_multi_%s_allY", jetPt.c_str()), "Measured/Truth/Unfolded (best iter) across y bins", w, h);
            cCompYMulti->Divide(nCols, nRows, 0.001, 0.001);
            for (int i = 0; i < nPanels; ++i)
            {
                cCompYMulti->cd(i + 1);
                gPad->SetTickx();
                gPad->SetTicky();
                gPad->SetRightMargin(0.02);
                gPad->SetLeftMargin(0.12);
                gPad->SetTopMargin(0.02);
                gPad->SetBottomMargin(0.14);
                TH1D *hM = comparePerYBin[i].meas;
                TH1D *hT = comparePerYBin[i].truth;
                TH1D *hU = comparePerYBin[i].unfolded;
                if (!hM || !hT || !hU)
                    continue;
                double maxv = std::max({hM->GetMaximum(), hT->GetMaximum(), hU->GetMaximum()});
                hM->SetMaximum(maxv * 1.25);
                hM->GetXaxis()->SetTitle("z_{T}");
                hM->GetYaxis()->SetTitle("dN/dz_{T} (norm)");
                hM->SetTitle("");
                hM->Draw("E");
                hT->SetLineStyle(2);
                hT->Draw("HISTSAME");
                hU->SetLineColor(kBlue + 1);
                hU->SetLineWidth(2);
                hU->Draw("HISTSAME");
                TLatex tl;
                tl.SetNDC();
                tl.SetTextSize(0.045);
                if (i == 0)
                    tl.DrawLatex(0.2, 0.9, "#font[12]{LHCb} in-progress");
                if (i == 0)
                    tl.DrawLatex(0.2, 0.85, Form("Jet p_{T}: %s", jetPt.c_str()));
                tl.DrawLatex(0.2, (i == 0) ? 0.80 : 0.90, Form("y: %s", comparePerYBin[i].yBinStr.Data()));
                if (i == 0)
                {
                    TLegend *leg = new TLegend(0.6, 0.75, 0.92, 0.95);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    leg->SetTextSize(0.04);
                    leg->AddEntry(hM, "Measured", "lep");
                    leg->AddEntry(hU, Form("Unfolded (iter %d)", comparePerYBin[i].iter), "l");
                    leg->AddEntry(hT, "MC truth", "l");
                    leg->Draw();
                }
            }
            TString pngCompY = Form("%s/compare_zT_multi_%s_allY.png", multiDir.c_str(), jetPt.c_str());
            cCompYMulti->SaveAs(pngCompY);
            cCompYMulti->Write();
            if (verbose)
                std::cout << "  - Saved multi-panel zT comparison across rapidity bins: " << pngCompY << std::endl;
            delete cCompYMulti;
            for (auto &e : comparePerYBin)
            {
                delete e.meas;
                delete e.truth;
                delete e.unfolded;
            }
        }
        // Build Measured/Truth/Missed/Fake multi-panel across rapidity bins
        if (!mfmfPerYBin.empty())
        {
            int nPanels = (int)mfmfPerYBin.size();
            int nCols = (nPanels <= 3) ? nPanels : (int)std::ceil(std::sqrt(nPanels));
            if (nCols < 1)
                nCols = 1;
            if (nCols > 6)
                nCols = 6;
            int nRows = (int)std::ceil((double)nPanels / nCols);
            int w = std::min(340 * nCols, 2400);
            int h = std::min(320 * nRows, 1800);
            TCanvas *cMFMF = new TCanvas(Form("c_meas_truth_miss_fake_multi_%s_allY", jetPt.c_str()), "Measured/Truth/Missed/Fake across y bins", w, h);
            cMFMF->Divide(nCols, nRows, 0.001, 0.001);
            for (int i = 0; i < nPanels; ++i)
            {
                cMFMF->cd(i + 1);
                gPad->SetTickx();
                gPad->SetTicky();
                gPad->SetRightMargin(0.02);
                gPad->SetLeftMargin(0.12);
                gPad->SetTopMargin(0.02);
                gPad->SetBottomMargin(0.14);
                TH1D *hM = mfmfPerYBin[i].meas;
                TH1D *hT = mfmfPerYBin[i].truth;
                TH1D *hMiss = mfmfPerYBin[i].missed;
                TH1D *hFake = mfmfPerYBin[i].fake;
                if (!hM || !hT || !hMiss || !hFake)
                    continue;
                double ymax = std::max({hM->GetMaximum(), hT->GetMaximum(), hMiss->GetMaximum(), hFake->GetMaximum()});
                hM->SetLineColor(kBlack);
                hM->SetMarkerStyle(20);
                hM->SetMarkerColor(kBlack);
                hM->SetLineWidth(2);
                hT->SetLineColor(kRed);
                hT->SetLineStyle(2);
                hT->SetLineWidth(2);
                hMiss->SetLineColor(kBlue + 1);
                hMiss->SetLineStyle(3);
                hMiss->SetLineWidth(2);
                hFake->SetLineColor(kMagenta + 2);
                hFake->SetLineStyle(4);
                hFake->SetLineWidth(2);
                hM->SetMaximum(ymax * 1.25);
                hM->GetXaxis()->SetTitle("z_{T}");
                hM->GetYaxis()->SetTitle("dN/dz_{T} (norm)");
                hM->SetTitle("");
                hM->Draw("HIST");
                hT->Draw("HISTSAME");
                if (hMiss->Integral() > 0)
                    hMiss->Draw("HISTSAME");
                if (hFake->Integral() > 0)
                    hFake->Draw("HISTSAME");
                TLatex tl;
                tl.SetNDC();
                tl.SetTextSize(0.045);
                if (i == 0)
                    tl.DrawLatex(0.20, 0.90, "#font[12]{LHCb} in-progress");
                if (i == 0)
                    tl.DrawLatex(0.20, 0.85, Form("Jet p_{T}: %s", jetPt.c_str()));
                tl.DrawLatex(0.20, (i == 0) ? 0.80 : 0.90, Form("y: %s", mfmfPerYBin[i].yBinStr.Data()));
                if (i == 0)
                {
                    TLegend *leg = new TLegend(0.50, 0.62, 0.90, 0.88);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    leg->SetTextSize(0.035);
                    leg->AddEntry(hM, "Measured", "l");
                    leg->AddEntry(hT, "MC truth", "l");
                    leg->AddEntry(hMiss, "Missed", "l");
                    leg->AddEntry(hFake, "Fake", "l");
                    leg->Draw();
                }
            }
            TString pngMFMF = Form("%s/meas_truth_miss_fake_multi_%s_allY.png", multiDir.c_str(), jetPt.c_str());
            cMFMF->SaveAs(pngMFMF);
            cMFMF->Write();
            if (verbose)
                std::cout << "  - Saved multi-panel Measured/Truth/Missed/Fake across rapidity bins: " << pngMFMF << std::endl;
            delete cMFMF;
            for (auto &e : mfmfPerYBin)
            {
                delete e.meas;
                delete e.truth;
                delete e.missed;
                delete e.fake;
            }
        }
        // Build response matrix multi-panel across rapidity bins
        if (!responsePerYBin.empty())
        {
            int nPanels = (int)responsePerYBin.size();
            int nCols = (nPanels <= 3) ? nPanels : (int)std::ceil(std::sqrt(nPanels));
            if (nCols < 1)
                nCols = 1;
            if (nCols > 6)
                nCols = 6;
            int nRows = (int)std::ceil((double)nPanels / nCols);
            int w = std::min(340 * nCols, 2400);
            int h = std::min(320 * nRows, 1800);
            // Determine global z-range
            double zMin = 1e300, zMax = -1e300;
            for (auto &yStr : responsePerYBin)
            {
                TH2 *hR = (TH2 *)fout->Get(Form("hResponse_%s_%s", jetPt.c_str(), yStr.Data()));
                if (!hR)
                    continue;
                int nx = hR->GetNbinsX();
                int ny = hR->GetNbinsY();
                for (int ix = 1; ix <= nx; ++ix)
                    for (int iy = 1; iy <= ny; ++iy)
                    {
                        double v = hR->GetBinContent(ix, iy);
                        if (v > zMax)
                            zMax = v;
                        if (v < zMin && v > 0)
                            zMin = v; // ignore zeros for lower bound
                    }
            }
            if (!(zMax > 0) || zMin >= zMax)
            {
                zMin = 0;
            }
            TCanvas *cRespMulti = new TCanvas(Form("c_response_multi_%s_allY", jetPt.c_str()), "Response matrices across y bins", w, h);
            cRespMulti->Divide(nCols, nRows, 0.001, 0.001);
            for (int i = 0; i < nPanels; ++i)
            {
                cRespMulti->cd(i + 1);
                gPad->SetTickx();
                gPad->SetTicky();
                gPad->SetRightMargin(0.12); // leave room for z axis
                gPad->SetLeftMargin(0.10);
                gPad->SetTopMargin(0.02);
                gPad->SetBottomMargin(0.14);
                TString yStr = responsePerYBin[i];
                TH2 *hR = (TH2 *)fout->Get(Form("hResponse_%s_%s", jetPt.c_str(), yStr.Data()));
                if (!hR)
                    continue;
                hR->SetStats(0);
                if (zMin < zMax)
                    hR->GetZaxis()->SetRangeUser(zMin, zMax * 1.0001);
                hR->GetXaxis()->SetTitle("Reco z_{T} bin");
                hR->GetYaxis()->SetTitle("Truth z_{T} bin");
                hR->SetTitle("");
                hR->Draw("COLZ");
                TLatex tl;
                tl.SetNDC();
                tl.SetTextSize(0.045);
                if (i == 0)
                    tl.DrawLatex(0.18, 0.92, "#font[12]{LHCb} in-progress");
                if (i == 0)
                    tl.DrawLatex(0.18, 0.86, Form("Jet p_{T}: %s", jetPt.c_str()));
                tl.DrawLatex(0.18, (i == 0) ? 0.80 : 0.92, Form("y: %s", yStr.Data()));
            }
            TString pngRespMulti = Form("%s/response_multi_%s_allY.png", multiDir.c_str(), jetPt.c_str());
            cRespMulti->SaveAs(pngRespMulti);
            cRespMulti->Write();
            if (verbose)
                std::cout << "  - Saved multi-panel response matrix across rapidity bins: " << pngRespMulti << std::endl;
            delete cRespMulti;
        }
        // Close per-jet input file
        fin->Close();
        delete fin;
    }

    // Close output file within function scope
    fout->Close();
    if (verbose)
        std::cout << "Unfolding complete. Output saved to " << outPath << std::endl;
}

// Example usage:
// unfold_new("input.root", "unfolded_output.root", 4, {"5_10", "10_20"}, {"2.5_3.0", "3.0_3.5"});
