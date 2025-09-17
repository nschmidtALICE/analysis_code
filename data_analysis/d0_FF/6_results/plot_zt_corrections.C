#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <iostream>
#include <algorithm>
#include <TLine.h>
#include <regex>
#include <string>
#include <vector>
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>

// plot_zt_corrections.C
// Macro to re-plot the saved tagZAllEfficiencyCorrections inputs
// Usage (in ROOT):
//   root -l plot_zt_corrections.C
//   plot_zt_corrections("/path/to/tagZAllEfficiencyCorrections_bin0.root","out_prefix")

void plot_zt_corrections(const char* infile = "ALL",
                         const char* outprefix = "5_10_tagZAllEfficiencyCorrections_bin2") {
// void plot_zt_corrections(const char* infile = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA/5_10/tagZAllEfficiencyCorrections_bin2.root",
//                          const char* outprefix = "5_10_tagZAllEfficiencyCorrections_bin2") {
    gStyle->SetOptStat(0);

    // Special batch mode: if infile == "ALL" iterate over pt ranges and rapidity bins
    if (std::string(infile) == "ALL") {
        std::vector<std::string> pt_pairs = {"5_10", "10_15", "15_20", "20_30", "30_50"};
        const std::string base_dir = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA";
        for (const auto &pt : pt_pairs) {
            for (int bin = 0; bin <= 7; ++bin) {
                std::string in = base_dir + "/" + pt + "/tagZAllEfficiencyCorrections_bin" + std::to_string(bin) + ".root";
                std::string out = pt + std::string("_tagZAllEfficiencyCorrections_bin") + std::to_string(bin);
                std::cout << "[ALL MODE] Processing: " << in << " -> " << out << std::endl;
                // call the same macro for each file (will follow normal single-file path)
                plot_zt_corrections(in.c_str(), out.c_str());
            }
        }
        return;
    }

    TFile *f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open file: " << infile << std::endl;
        return;
    }
    std::vector<double> yBins = {2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0};
    TString rapidityString = "";
    //search infile for bin number and set rapidityString accordingly
    try {
        std::string s = infile;
        std::smatch m;
        std::regex re("bin(\\d+)");
        if (std::regex_search(s, m, re) && m.size() > 1) {
            int bin = std::stoi(m[1].str());
            if (bin >= 0 && bin + 1 < (int)yBins.size()) {
                rapidityString = TString::Format("%.1f - %.1f", yBins[bin], yBins[bin+1]);
            } else if (bin >= 0 && bin < (int)yBins.size()) {
                rapidityString = TString::Format("%.1f", yBins[bin]);
            } else {
                rapidityString = TString::Format("bin%d", bin);
            }
        }
    } catch (...) {
        rapidityString = "";
    }

    std::cout << "Using rapidity string: " << rapidityString.Data() << std::endl;
    

    // Try to load graphs
    TGraphErrors *gKaon = dynamic_cast<TGraphErrors*>(f->Get("tagZKaonCorrGraph"));
    TGraphErrors *gPion = dynamic_cast<TGraphErrors*>(f->Get("tagZPionCorrGraph"));
    TGraphErrors *gReco = dynamic_cast<TGraphErrors*>(f->Get("tagZRecoEffCorrGraph"));
    TGraphErrors *gAcc = dynamic_cast<TGraphErrors*>(f->Get("tagZAcceptanceCorrGraph"));
    TGraphErrors *gCombined = dynamic_cast<TGraphErrors*>(f->Get("tagZCombinedCorrGraph"));

    // Also try to load underlying histograms if graphs are missing
    TH1D *hKaon = dynamic_cast<TH1D*>(f->Get("tagZKaonCorrHist"));
    TH1D *hPion = dynamic_cast<TH1D*>(f->Get("tagZPionCorrHist"));
    TH1D *hReco = dynamic_cast<TH1D*>(f->Get("tagZRecoEffCorrHist"));
    TH1D *hAcc = dynamic_cast<TH1D*>(f->Get("tagZAcceptanceCorrHist"));
    TH1D *hCombined = dynamic_cast<TH1D*>(f->Get("tagZCombinedCorrHist"));

    if (!gKaon && !hKaon && !gPion && !hPion && !gCombined && !hCombined) {
        std::cerr << "ERROR: no usable inputs found in file. Expected tagZKaonCorrGraph / tagZKaonCorrHist etc." << std::endl;
        f->Close(); delete f;
        return;
    }

    // Function to convert TH1D -> TGraphErrors for drawing with points+errors
    auto histToGraph = [](TH1D* h) -> TGraphErrors* {
        if (!h) return nullptr;
        int nb = h->GetNbinsX();
        TGraphErrors* g = new TGraphErrors();
        int idx = 0;
        for (int b = 1; b <= nb; ++b) {
            double val = h->GetBinContent(b);
            double err = h->GetBinError(b);
            if (val == 0 && err == 0) continue; // skip empty
            double x = h->GetBinCenter(b);
            double ex = h->GetBinWidth(b)/2.0;
            g->SetPoint(idx, x, val);
            g->SetPointError(idx, ex, err);
            ++idx;
        }
        return g;
    };

    // If graphs are missing but histograms exist, convert
    if (!gKaon && hKaon) gKaon = histToGraph(hKaon);
    if (!gPion && hPion) gPion = histToGraph(hPion);
    if (!gReco && hReco) gReco = histToGraph(hReco);
    if (!gAcc && hAcc) gAcc = histToGraph(hAcc);
    if (!gCombined && hCombined) gCombined = histToGraph(hCombined);

    // Compute y-range
    double ymax = 0.0;
    auto considerGraph = [&](TGraphErrors* g) {
        if (!g) return;
        int n = g->GetN();
        for (int i = 0; i < n; ++i) {
            double y; g->GetPoint(i, *(new double), y);
            ymax = std::max(ymax, y);
        }
    };
    // safer approach to avoid dynamic double leak: iterate via arrays
    auto considerGraphSafe = [&](TGraphErrors* g) {
        if (!g) return;
        int n = g->GetN();
        double xval, yval;
        for (int i = 0; i < n; ++i) {
            g->GetPoint(i, xval, yval);
            ymax = std::max(ymax, yval);
        }
    };
    considerGraphSafe(gKaon);
    considerGraphSafe(gPion);
    considerGraphSafe(gReco);
    considerGraphSafe(gAcc);
    considerGraphSafe(gCombined);
    // also consider histograms
    auto considerHist = [&](TH1D* h) {
        if (!h) return;
        ymax = std::max(ymax, h->GetMaximum());
    };
    considerHist(hKaon); considerHist(hPion); considerHist(hReco); considerHist(hAcc); considerHist(hCombined);

    if (ymax < 0.1) ymax = 1.0; // default
    double ytop = std::min(1.1, ymax * 1.2);
    if (ytop < 0.2) ytop = 0.2;

    // Create canvas matching original layout
    TCanvas *c = new TCanvas(outprefix, "All Efficiency Corrections vs TagZ", 1000, 1000);
    c->SetLeftMargin(0.08);
    c->SetRightMargin(0.01);
    c->SetTopMargin(0.01);
    c->SetBottomMargin(0.08);
    c->SetTickx();
    c->SetTicky();

    // Prepare first graph to draw axes
    TGraphErrors* gBase = gKaon ? gKaon : (gPion ? gPion : (gCombined ? gCombined : nullptr));
    if (!gBase) {
        std::cerr << "ERROR: no baseline graph to set axes." << std::endl;
        f->Close(); delete f;
        return;
    }

    // Create an empty histogram to hold axis ranges
    double xmin = 0.0, xmax = 1.0;
    // attempt to find x-range from any available histogram
    // if (hKaon) { xmin = hKaon->GetXaxis()->GetXmin(); xmax = hKaon->GetXaxis()->GetXmax(); }
    // else if (hPion) { xmin = hPion->GetXaxis()->GetXmin(); xmax = hPion->GetXaxis()->GetXmax(); }
    // else if (gBase && gBase->GetN() > 0) { double x,y; gBase->GetPoint(0,x,y); xmin = std::min(xmin,x); gBase->GetPoint(gBase->GetN()-1,x,y); xmax = std::max(xmax,x); }

    TH1D *hAxis = new TH1D("hAxis",";#it{z}_{T};Efficiency", 100, xmin, xmax);
    hAxis->SetDirectory(nullptr);
    hAxis->SetStats(0);
    hAxis->GetYaxis()->SetRangeUser(0.0, 1.39);
    hAxis->GetYaxis()->SetTitleOffset(1.1);
    hAxis->Draw();

    // draw horizontal guide line at y = 1.0
    TLine *hline = new TLine(xmin, 1.0, xmax, 1.0);
    hline->SetLineStyle(2);
    hline->SetLineColor(kBlack);
    hline->SetLineWidth(1);
    hline->Draw("same");



    // Draw labels: LHCb in-progress, jet pT range, rapidity range
    TString ptRangeString = "";
    try {
        std::string s = infile;
        std::smatch mpt;
        std::regex rept("([0-9]+_[0-9]+)");
        if (std::regex_search(s, mpt, rept) && mpt.size() > 1) {
            std::string prs = mpt[1].str();
            // replace underscore with hyphen
            for (auto &c : prs) if (c == '_') c = '-';
            ptRangeString = TString::Format("jet p_{T}: %s GeV/#it{c}", prs.c_str());
        }
    } catch (...) { ptRangeString = ""; }

    TLatex lab;
    lab.SetNDC();
    lab.SetTextFont(62);
    lab.SetTextSize(0.04);
    lab.DrawLatex(0.12, 0.93, "LHCb (in-progress)");
    lab.SetTextFont(42);
    lab.SetTextSize(0.033);
    if (ptRangeString.Length() > 0) lab.DrawLatex(0.12, 0.89, ptRangeString);
    if (rapidityString.Length() > 0) lab.DrawLatex(0.12, 0.85, TString::Format("rapidity: %s", rapidityString.Data()));

    // Draw graphs in order with original style
    if (gKaon) {
        gKaon->SetMarkerStyle(20);
        gKaon->SetMarkerSize(2.0);
        gKaon->SetMarkerColor(kBlue+2);
        gKaon->SetLineColor(kBlue+2);
        gKaon->Draw("PE same");
    }
    if (gPion) {
        gPion->SetMarkerStyle(21);
        gPion->SetMarkerSize(2.0);
        gPion->SetMarkerColor(kRed+2);
        gPion->SetLineColor(kRed+2);
        gPion->Draw("PE same");
    }
    if (gCombined) {
        gCombined->SetMarkerStyle(22);
        gCombined->SetMarkerSize(2.0);
        gCombined->SetMarkerColor(kMagenta+2);
        gCombined->SetLineColor(kMagenta+2);
        gCombined->Draw("PE same");
    }
    if (gReco && gReco->GetN() > 0) {
        gReco->SetMarkerStyle(29);
        gReco->SetMarkerSize(2.5);
        gReco->SetMarkerColor(kOrange+2);
        gReco->SetLineColor(kOrange+2);
        gReco->Draw("PE same");
    }
    if (gAcc && gAcc->GetN() > 0) {
        gAcc->SetMarkerStyle(47);
        gAcc->SetMarkerSize(2.0);
        gAcc->SetMarkerColor(kAzure+2);
        gAcc->SetLineColor(kAzure+2);
        gAcc->Draw("PE same");
    }
    
    // Legend
    TLegend *leg = new TLegend(0.65, 0.76, 0.95, 0.96);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    if (gKaon) leg->AddEntry(gKaon, "Kaon PID", "pe");
    if (gPion) leg->AddEntry(gPion, "Pion PID", "pe");
    if (gCombined) leg->AddEntry(gCombined, "Combined PID", "pe");
    if (gReco && gReco->GetN()>0) leg->AddEntry(gReco, "Reconstruction", "pe");
    if (gAcc && gAcc->GetN()>0) leg->AddEntry(gAcc, "Acceptance", "pe");
    leg->Draw();

    // TLatex lat;

    // Save into a dated directory: ./plots_YYYYMMDD/
    char datebuf[32] = {0};
    time_t t = time(nullptr);
    struct tm tm = *localtime(&t);
    strftime(datebuf, sizeof(datebuf), "%Y%m%d", &tm);
    std::string outdir = std::string("plots_") + datebuf;
    // create directory if it doesn't exist
    struct stat st;
    if (stat(outdir.c_str(), &st) != 0) {
        if (mkdir(outdir.c_str(), 0755) != 0) {
            std::cerr << "WARNING: could not create output directory: " << outdir << " - will save in current dir instead." << std::endl;
            outdir = ".";
        }
    }

    std::string png = outdir + "/" + std::string(outprefix) + ".png";
    std::string pdf = outdir + "/" + std::string(outprefix) + ".pdf";
    c->SaveAs(png.c_str());
    c->SaveAs(pdf.c_str());

    std::cout << "Saved replot to: " << png << " and " << pdf << std::endl;

    // cleanup
    delete hAxis;
    // only delete converted graphs created from histToGraph (we created none that are tracked separately here)
    f->Close(); delete f;
}
