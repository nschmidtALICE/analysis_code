// JetFindingEfficiencyReplot.cpp
// Standalone re-plotter: loads jetFindingEffMinimal ROOT outputs and recreates
// all summary and per-pt-bin plots (efficiency map, predicted, residual,
// per-pt efficiency with fits, closure overlay, residual). Optionally refits.
//
// Build:
//   g++ -std=c++17 JetFindingEfficiencyReplot.cpp $(root-config --cflags --libs) -o JetFindingEfficiencyReplot
// Usage:
//   ./JetFindingEfficiencyReplot path/to/jetFindingEffMinimal.root [outPrefix] [--no-refit]
// If outPrefix omitted, uses "replot" next to input file name.

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <regex>
#include <ctime>
#include <sstream>
#include <iomanip>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TPaletteAxis.h"
#include "TLine.h"

namespace fs = std::filesystem;

struct PtBin
{
    double lo;
    double hi;
};

static TF1 *MakeLogistic(const std::string &name, double xmin, double xmax)
{
    TF1 *f = new TF1(name.c_str(), "[0]/(1+exp(-[1]*(x-[2])))", xmin, xmax);
    f->SetParameters(0.85, 6.0, 0.6);
    f->SetParLimits(0, 0, 1.0);
    f->SetParLimits(1, 0, 100.0);
    f->SetParLimits(2, xmin, xmax);
    f->SetLineColor(kRed);
    f->SetLineWidth(2);
    return f;
}
static TF1 *MakeExpSat(const std::string &name, double xmin, double xmax)
{
    TF1 *f = new TF1(name.c_str(), "[0]+[1]*(1-exp(-[2]*x))", xmin, xmax);
    f->SetParameters(0.05, 0.8, 2.0);
    f->SetParLimits(0, 0, 1.0);
    f->SetParLimits(1, 0, 1.0);
    f->SetParLimits(2, 1e-3, 100.0);
    f->SetLineColor(kGreen + 2);
    f->SetLineStyle(2);
    f->SetLineWidth(2);
    return f;
}
static TF1 *MakePol2(const std::string &name, double xmin, double xmax)
{
    TF1 *f = new TF1(name.c_str(), "pol2", xmin, xmax);
    f->SetLineColor(kBlue + 2);
    f->SetLineStyle(4);
    f->SetLineWidth(3);
    return f;
}

static void PrintFitSummary(const TF1 *f, const PtBin &pb, const std::string &tag)
{
    if (!f)
        return;
    double chi2 = f->GetChisquare();
    double ndf = f->GetNDF();
    std::cout << "Fit " << tag << " pt[" << pb.lo << "," << pb.hi << "]: chi2=" << chi2 << " ndf=" << ndf
              << " chi2/ndf=" << (ndf > 0 ? chi2 / ndf : 0) << "\n";
    for (int i = 0; i < f->GetNpar(); ++i)
    {
        std::cout << "  " << tag << ": p" << i << " = " << f->GetParameter(i)
                  << " +/- " << f->GetParError(i) << "\n";
    }
}

int JetFindingEfficiencyReplot()
{

    std::string inputPath = "/media/niviths/local/analysis_code/data_analysis/d0_FF/7_efficiency/jetFindingEffMinimal_fullMC_outputs_2025-10-14/jetFindingEffMinimal_fullMC.root";
    // std::string inputPath = "/media/niviths/local/analysis_code/data_analysis/d0_FF/7_efficiency/jetFindingEffMinimal_outputs/jetFindingEffMinimal.root";
    std::string outPrefix = "jetFindingEffMinimal_replot";
    bool doRefit = true;

    TFile *fin = TFile::Open(inputPath.c_str(), "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "Cannot open input ROOT file\n";
        return 2;
    }

    // Build date stamp for versioned output directory (YYYYMMDD)
    std::string dateStamp;
    {
        std::time_t t = std::time(nullptr);
        std::tm tm{};
        localtime_r(&t, &tm);
        std::ostringstream oss;
        oss << std::put_time(&tm, "%Y%m%d");
        dateStamp = oss.str();
    }
    // Output directory (parallel to input or cwd) with date tag
    std::string outDir = outPrefix + "_outputs_" + dateStamp;
    try
    {
        fs::create_directories(outDir);
    }
    catch (...)
    {
        std::cerr << "Warning: cannot create output dir\n";
        outDir = ".";
    }
    std::string pngDir = outDir + "/png";
    try
    {
        fs::create_directories(pngDir);
    }
    catch (...)
    {
    }

    auto savePNG = [&](TCanvas &c, const std::string &path)
    {
        c.SaveAs(path.c_str());
        try
        {
            fs::path p(path);
            fs::path d = fs::path(pngDir) / p.filename();
            c.SaveAs(d.string().c_str());
        }
        catch (...)
        {
        }
    };

    // Retrieve main histograms
    TH2 *h_eff = dynamic_cast<TH2 *>(fin->Get("jet_finding_efficiency_zT"));
    TH2 *h_pred = dynamic_cast<TH2 *>(fin->Get("predicted_reco_zT"));
    TH2 *h_diff = dynamic_cast<TH2 *>(fin->Get("diff_reco_minus_pred"));
    if (!h_eff || !h_pred || !h_diff)
    {
        std::cerr << "Missing one or more core histograms in file." << std::endl;
    }

    // Recreate the three main canvases
    if (h_eff)
    {
        TCanvas c("c_eff_re", "Jet finding efficiency", 800, 600);
        // Set a modest right margin; palette will start exactly at the frame edge (1 - rightMargin)
        c.SetRightMargin(0.11);
        c.SetTopMargin(0.02);
        c.SetBottomMargin(0.08);
        c.SetLeftMargin(0.08);
        h_eff->SetStats(0);
        h_eff->GetZaxis()->SetTitle("Jet Finding Efficiency");
        h_eff->SetTitle("");
        h_eff->Draw("COLZ");
        c.Update();
        if (auto *palette = (TPaletteAxis *)h_eff->GetListOfFunctions()->FindObject("palette"))
        {
            // Place palette just to the right of the frame with a small gap
            double frameRight = 1.0 - gPad->GetRightMargin();
            double padTop = 1.0 - gPad->GetTopMargin();
            double padBottom = gPad->GetBottomMargin();
            double gap = 0.004; // small horizontal gap
            double w = 0.030;   // palette width (fraction of pad)
            palette->SetX1NDC(frameRight + gap);
            palette->SetX2NDC(std::min(0.995, frameRight + gap + w));
            palette->SetY1NDC(padBottom);
            palette->SetY2NDC(padTop);
            gPad->Modified();
            gPad->Update();
        }
        // Add labels
        {
            TLatex t;
            t.SetNDC();
            // LHCb work-in-progress
            t.SetTextFont(62); // bold
            t.SetTextSize(0.04);
            t.DrawLatex(0.12, 0.93, "LHCb work-in-progress");
            // Jet finder info
            t.SetTextFont(42);
            t.SetTextSize(0.04);
            t.DrawLatex(0.12, 0.88, "anti-k_{T}, R=0.5");
        }
        std::string base = outDir + "/" + outPrefix + "_eff_zT";
        savePNG(c, base + ".png");
        c.SaveAs((base + ".pdf").c_str());
    }
    if (h_pred)
    {
        TCanvas c("c_pred_re", "Predicted reco", 800, 600);
        c.SetRightMargin(0.20);
        h_pred->SetStats(0);
        h_pred->GetZaxis()->SetTitle("Counts");
        h_pred->Draw("COLZ");
        std::string base = outDir + "/" + outPrefix + "_pred_zT";
        savePNG(c, base + ".png");
        c.SaveAs((base + ".pdf").c_str());
    }
    if (h_diff)
    {
        TCanvas c("c_res_re", "Reco - Pred", 800, 600);
        c.SetRightMargin(0.20);
        h_diff->SetStats(0);
        h_diff->GetZaxis()->SetTitle("Counts (data - pred)");
        h_diff->Draw("COLZ");
        std::string base = outDir + "/" + outPrefix + "_res_zT";
        savePNG(c, base + ".png");
        c.SaveAs((base + ".pdf").c_str());
    }

    // Pt bin definitions (mirror original)
    std::vector<PtBin> ptBins = {{5, 10}, {10, 15}, {15, 20}, {20, 30}, {30, 50}};
    // Containers to collect per-pt projections for later multi-panel plotting
    std::vector<TH1 *> multi_eff_vec;
    std::vector<TH1 *> multi_num_vec;
    std::vector<TH1 *> multi_pred_vec;
    std::vector<TH1 *> multi_diff_vec;
    // store per-pt histograms that represent the Pol2 fit evaluated at bin centers
    std::vector<TH1 *> multi_pol2_hist_vec;
    // Store per-pt fit functions (cloned) so multi-plot can reuse the same fits
    std::vector<TF1 *> multi_f_log_vec;
    std::vector<TF1 *> multi_f_exp_vec;
    std::vector<TF1 *> multi_f_pol2_vec;
    // Derive numerators / denominators for projections
    TH2 *h_num = dynamic_cast<TH2 *>(fin->Get("h_num_zT"));
    TH2 *h_den = dynamic_cast<TH2 *>(fin->Get("h_den_zT"));
    if (!h_num || !h_den)
    {
        std::cerr << "Warning: missing h_num_zT or h_den_zT; per-pt plots will be skipped." << std::endl;
    }

    // Attempt to load precomputed pt-bin objects if present; else compute projections.
    const double zTmin = 0.0, zTmax = 1.0;
    for (const auto &pb : ptBins)
    {
        if (!h_num || !h_den || !h_pred || !h_diff)
            break;
        int firstY = h_den->GetYaxis()->FindBin(pb.lo + 1e-6);
        int lastY = h_den->GetYaxis()->FindBin(pb.hi - 1e-6);
        if (lastY < firstY)
            lastY = firstY;
        TH1 *den_proj = h_den->ProjectionX(Form("den_pt_%.0f_%.0f_re", pb.lo, pb.hi), firstY, lastY);
        TH1 *num_proj = h_num->ProjectionX(Form("num_pt_%.0f_%.0f_re", pb.lo, pb.hi), firstY, lastY);
        TH1 *pred_proj = h_pred->ProjectionX(Form("pred_pt_%.0f_%.0f_re", pb.lo, pb.hi), firstY, lastY);
        TH1 *diff_proj = h_diff->ProjectionX(Form("diff_pt_%.0f_%.0f_re", pb.lo, pb.hi), firstY, lastY);
        TH1 *eff_proj = dynamic_cast<TH1 *>(num_proj->Clone(Form("eff_pt_%.0f_%.0f_re", pb.lo, pb.hi)));
        eff_proj->SetTitle(Form("Jet finding efficiency zT (pt %.0f-%.0f);z_{T};Efficiency", pb.lo, pb.hi));
        eff_proj->Divide(eff_proj, den_proj, 1.0, 1.0, "B");

        // Remove points with exactly y=1.0 efficiency (likely problematic bins)
        for (int bin = 1; bin <= eff_proj->GetNbinsX(); ++bin) {
            double y = eff_proj->GetBinContent(bin);
            if (std::abs(y - 1.0) < 1e-9) {  // exactly 1.0 within numerical precision
                eff_proj->SetBinContent(bin, 0.0);
                eff_proj->SetBinError(bin, 0.0);
            }
        }

        // Clone and stash per-pt projections for multi-panel plots (owned by us)
        try
        {
            multi_eff_vec.push_back(dynamic_cast<TH1 *>(eff_proj->Clone(Form("multi_eff_pt%.0f_%.0f", pb.lo, pb.hi))));
            multi_num_vec.push_back(dynamic_cast<TH1 *>(num_proj->Clone(Form("multi_num_pt%.0f_%.0f", pb.lo, pb.hi))));
            multi_pred_vec.push_back(dynamic_cast<TH1 *>(pred_proj->Clone(Form("multi_pred_pt%.0f_%.0f", pb.lo, pb.hi))));
            multi_diff_vec.push_back(dynamic_cast<TH1 *>(diff_proj->Clone(Form("multi_diff_pt%.0f_%.0f", pb.lo, pb.hi))));
        }
        catch (...)
        {
            // cloning failed; continue without multi-panel entry
        }

        // Efficiency + fits
        TCanvas cEff(Form("c_eff_pt_%.0f_%.0f_re", pb.lo, pb.hi), Form("Efficiency pt %.0f-%.0f", pb.lo, pb.hi), 700, 500);
        gPad->SetTickx();
        gPad->SetTicky();
        eff_proj->SetStats(0);
        eff_proj->SetLineColor(kBlack);
        eff_proj->SetLineWidth(2);
        eff_proj->SetMarkerStyle(20);
        eff_proj->Draw("pe, E");
        TF1 *f_log = nullptr;
        TF1 *f_exp = nullptr;
        TF1 *f_pol2 = nullptr;
        if (doRefit)
        {
            f_log = MakeLogistic(Form("f_log_pt%.0f_%.0f_re", pb.lo, pb.hi), zTmin, zTmax);
            f_exp = MakeExpSat(Form("f_exp_pt%.0f_%.0f_re", pb.lo, pb.hi), zTmin, zTmax);
            f_pol2 = MakePol2(Form("f_pol2_pt%.0f_%.0f_re", pb.lo, pb.hi), zTmin, zTmax);
            eff_proj->Fit(f_log, "R W Q 0");
            eff_proj->Fit(f_exp, "R W Q 0");
            eff_proj->Fit(f_pol2, "R W Q 0");
            f_log->Draw("SAME");
            f_exp->Draw("SAME");
            f_pol2->Draw("SAME");
            PrintFitSummary(f_log, pb, "Logistic");
            PrintFitSummary(f_exp, pb, "Exp");
            PrintFitSummary(f_pol2, pb, "Pol2");
        }
        else
        {
            // Try to fetch stored fits with original names (without _re) if they exist
            f_log = dynamic_cast<TF1 *>(fin->Get(Form("f_log_pt%.0f_%.0f", pb.lo, pb.hi)));
            f_exp = dynamic_cast<TF1 *>(fin->Get(Form("f_exp_pt%.0f_%.0f", pb.lo, pb.hi)));
            f_pol2 = dynamic_cast<TF1 *>(fin->Get(Form("f_pol2_pt%.0f_%.0f", pb.lo, pb.hi)));
            if (f_log)
                f_log->Draw("SAME");
            if (f_exp)
                f_exp->Draw("SAME");
            if (f_pol2)
                f_pol2->Draw("SAME");
        }

            // create a histogram that holds the Pol2 fit values at the bin centers
            TH1 *pol2_hist = dynamic_cast<TH1 *>(eff_proj->Clone(Form("pol2hist_pt%.0f_%.0f", pb.lo, pb.hi)));
            pol2_hist->SetTitle(Form("Pol2-evaluated (pt %.0f-%.0f);z_{T};Efficiency", pb.lo, pb.hi));
            pol2_hist->SetDirectory(nullptr); // keep owned by us
            // if a Pol2 fit exists fill the bin contents with the fit evaluated at bin centers
            if (f_pol2) {
                for (int bin = 1; bin <= pol2_hist->GetNbinsX(); ++bin) {
                    double x = pol2_hist->GetBinCenter(bin);
                    double y = f_pol2->Eval(x);
                    pol2_hist->SetBinContent(bin, y);
                    pol2_hist->SetBinError(bin, 0.0);
                }
            } else {
                // otherwise clear contents
                for (int bin = 1; bin <= pol2_hist->GetNbinsX(); ++bin) {
                    pol2_hist->SetBinContent(bin, 0.0);
                    pol2_hist->SetBinError(bin, 0.0);
                }
            }
            pol2_hist->SetLineColor(kBlue);
            pol2_hist->SetLineWidth(2);
            multi_pol2_hist_vec.push_back(pol2_hist);
        TLegend legF(0.60, 0.15, 0.92, 0.38);
        legF.AddEntry(eff_proj, "Efficiency (data)", "pe");
        if (f_log)
            legF.AddEntry(f_log, "Logistic", "l");
        if (f_exp)
            legF.AddEntry(f_exp, "Saturating exp.", "l");
        // add pol2-evaluated histogram to legend and draw it (blue line)
        TH1 *pol2_hist_local = nullptr;
        if (!multi_pol2_hist_vec.empty() && multi_pol2_hist_vec.size() >= multi_eff_vec.size()) {
            // if we prefilled multi_pol2_hist_vec earlier, match by size
            pol2_hist_local = multi_pol2_hist_vec.back();
        }
        if (f_pol2 && pol2_hist_local == nullptr) {
            // fallback: create on-the-fly pol2 representation from f_pol2
            pol2_hist_local = dynamic_cast<TH1 *>(eff_proj->Clone(Form("pol2_onfly_pt%.0f_%.0f", pb.lo, pb.hi)));
            pol2_hist_local->SetDirectory(nullptr);
            for (int bin = 1; bin <= pol2_hist_local->GetNbinsX(); ++bin) {
                double x = pol2_hist_local->GetBinCenter(bin);
                double y = f_pol2->Eval(x);
                pol2_hist_local->SetBinContent(bin, y);
                pol2_hist_local->SetBinError(bin, 0.0);
            }
            pol2_hist_local->SetLineColor(kBlue);
            pol2_hist_local->SetLineWidth(2);
        }
        if (pol2_hist_local) {
            // pol2_hist_local->SetMarkerSize(0);
            pol2_hist_local->Draw("HIST SAME");
            legF.AddEntry(pol2_hist_local, "Pol2 (eval)", "l");
        }
        if (f_pol2)
            legF.AddEntry(f_pol2, "Pol2", "l");
        legF.SetBorderSize(0);
        legF.Draw();
        std::string baseEff = outDir + "/" + outPrefix + Form("_eff_zT_pt%.0f-%.0f", pb.lo, pb.hi);
        savePNG(cEff, baseEff + ".png");
        cEff.SaveAs((baseEff + ".pdf").c_str());

        // Closure overlay num vs pred
        TCanvas cClos(Form("c_closure_pt_%.0f_%.0f_re", pb.lo, pb.hi), Form("Closure pt %.0f-%.0f", pb.lo, pb.hi), 700, 500);
        gPad->SetTickx();
        gPad->SetTicky();
        num_proj->SetLineColor(kBlue);
        num_proj->SetLineWidth(2);
        pred_proj->SetLineColor(kRed);
        pred_proj->SetLineWidth(2);
        pred_proj->SetLineStyle(2);
        num_proj->Draw("HIST");
        pred_proj->Draw("HIST SAME");
        TLegend legC(0.20, 0.70, 0.48, 0.88);
        legC.AddEntry(num_proj, "Reco (data)", "l");
        legC.AddEntry(pred_proj, "Predicted (den*eff)", "l");
        legC.Draw();
        std::string baseClos = outDir + "/" + outPrefix + Form("_closure_zT_pt%.0f-%.0f", pb.lo, pb.hi);
        savePNG(cClos, baseClos + ".png");
        cClos.SaveAs((baseClos + ".pdf").c_str());

        // Residual
        TCanvas cRes(Form("c_res_pt_%.0f_%.0f_re", pb.lo, pb.hi), Form("Residual pt %.0f-%.0f", pb.lo, pb.hi), 700, 500);
        gPad->SetTickx();
        gPad->SetTicky();
        diff_proj->SetLineColor(kBlack);
        diff_proj->SetLineWidth(2);
        diff_proj->Draw("HIST");
        std::string baseRes = outDir + "/" + outPrefix + Form("_res_zT_pt%.0f-%.0f", pb.lo, pb.hi);
        savePNG(cRes, baseRes + ".png");
        cRes.SaveAs((baseRes + ".pdf").c_str());

        // Clean up local clones (fits only if refit to avoid double delete for borrowed ones)
        delete den_proj;
        delete num_proj;
        delete pred_proj;
        delete diff_proj;
        delete eff_proj;
        // don't delete pol2_hist_local here if it points into multi_pol2_hist_vec; only delete if we created it on the fly
        if (pol2_hist_local && (multi_pol2_hist_vec.empty() || pol2_hist_local != multi_pol2_hist_vec.back())) {
            delete pol2_hist_local;
        }

        // clone the fits (if present) for multi-panel use; clone before deleting local objects
        try
        {
            if (f_log)
                multi_f_log_vec.push_back(dynamic_cast<TF1 *>(f_log->Clone(Form("multi_f_log_pt%.0f_%.0f", pb.lo, pb.hi))));
            else
                multi_f_log_vec.push_back(nullptr);
            if (f_exp)
                multi_f_exp_vec.push_back(dynamic_cast<TF1 *>(f_exp->Clone(Form("multi_f_exp_pt%.0f_%.0f", pb.lo, pb.hi))));
            else
                multi_f_exp_vec.push_back(nullptr);
            if (f_pol2)
                multi_f_pol2_vec.push_back(dynamic_cast<TF1 *>(f_pol2->Clone(Form("multi_f_pol2_pt%.0f_%.0f", pb.lo, pb.hi))));
            else
                multi_f_pol2_vec.push_back(nullptr);
        }
        catch (...)
        {
            // cloning failed; push nulls to keep indices aligned
            while (multi_f_log_vec.size() < multi_eff_vec.size())
                multi_f_log_vec.push_back(nullptr);
            while (multi_f_exp_vec.size() < multi_eff_vec.size())
                multi_f_exp_vec.push_back(nullptr);
            while (multi_f_pol2_vec.size() < multi_eff_vec.size())
                multi_f_pol2_vec.push_back(nullptr);
        }
    }

    // --- Create multi-panel summaries from collected per-pt projections ---
    if (!multi_eff_vec.empty())
    {
        int nPanels = (int)multi_eff_vec.size();
        int nCols = (nPanels <= 3) ? nPanels : (int)std::ceil(std::sqrt(nPanels));
        if (nCols < 1)
            nCols = 1;
        if (nCols > 6)
            nCols = 6;
        int nRows = (int)std::ceil((double)nPanels / nCols);
        int w = std::min(340 * nCols, 2400);
        int h = std::min(320 * nRows, 1800);
        std::cout << "Multi-panel: " << nPanels << " panels in " << nCols << "x" << nRows << " grid, canvas " << w << "x" << h << "\n";
        TCanvas *cMultiEff = new TCanvas(Form("c_replot_eff_multi_%s", outPrefix.c_str()), "Efficiency per-pt bins", w, h);
        cMultiEff->Divide(nCols, nRows, 0.001, 0.001);
        for (int i = 0; i < nPanels; ++i)
        {
            cMultiEff->cd(i + 1);
            gPad->SetRightMargin(0.01);
            gPad->SetLeftMargin(0.09);
            gPad->SetTopMargin(0.01);
            gPad->SetBottomMargin(0.08);
            gPad->SetTickx();
            gPad->SetTicky();
            TH1 *h = multi_eff_vec[i];
            if (!h)
                continue;
            h->SetStats(0);
            h->SetLineWidth(1);
            h->SetLineColor(kBlack);
            h->SetMarkerStyle(20);
            h->GetYaxis()->SetTitle("jet finding efficiency");
            h->SetTitle("");
            // ensure full vertical range so fit curves are not clipped
            if (h->GetMaximum() < 1.0)
                h->GetYaxis()->SetRangeUser(0.41, 1.49);
            h->Draw("PE");


            // Draw fits for this panel using the cloned TF1s we stored earlier from the per-pt loop.
            TF1 *f_log = nullptr;
            TF1 *f_exp = nullptr;
            TF1 *f_pol2 = nullptr;
            if (i < (int)multi_f_log_vec.size())
                f_log = multi_f_log_vec[i];
            if (i < (int)multi_f_exp_vec.size())
                f_exp = multi_f_exp_vec[i];
            if (i < (int)multi_f_pol2_vec.size())
                f_pol2 = multi_f_pol2_vec[i];
            PtBin pb = {0.0, 0.0};
            if (i < (int)ptBins.size())
                pb = ptBins[i];
            
            // Draw fits first
            // if (f_log)
            // {
            //     f_log->SetNpx(500);
            //     f_log->SetLineWidth(2);
            //     f_log->Draw("SAME");
            //     // PrintFitSummary(f_log, pb, Form("Multi_Logistic_pt%d", i));
            // }
            // if (f_exp)
            // {
            //     f_exp->SetNpx(500);
            //     f_exp->SetLineWidth(2);
            //     f_exp->Draw("SAME");
            //     // PrintFitSummary(f_exp, pb, Form("Multi_Exp_pt%d", i));
            // }
            if (f_pol2)
            {
                f_pol2->SetNpx(500);
                f_pol2->SetLineWidth(2);
                f_pol2->Draw("SAME");
                // PrintFitSummary(f_pol2, pb, Form("Multi_Pol2_pt%d", i));
            }
            h->Draw("PE,same");
            // draw the pol2-evaluated histogram if available
            TH1 *pol2hist = nullptr;
            if (i < (int)multi_pol2_hist_vec.size())
                pol2hist = multi_pol2_hist_vec[i];
            if (pol2hist) {
                pol2hist->SetLineColor(kBlue);
                pol2hist->SetMarkerColor(kBlue);
                pol2hist->SetLineWidth(1);
                // pol2hist->SetMarkerSize(0);
                pol2hist->SetMarkerStyle(24);
                pol2hist->Draw("p SAME");
            }
            
            // Create legend as pointer to ensure proper ownership
            TLegend *legA = new TLegend(0.15, 0.72, 0.42, 0.94);
            legA->SetTextSize(0.032);
            legA->SetBorderSize(0);
            legA->SetFillStyle(0);
            legA->SetFillColor(0);
            legA->AddEntry(h, Form("Efficiency (data, jet-p_{T}: %2.0f - %2.0f GeV/#it{c})", pb.lo, pb.hi), "pe");
            // if (f_log) legA->AddEntry(f_log, "Logistic", "l");
            // if (f_exp) legA->AddEntry(f_exp, "Saturating exp.", "l");
            if (f_pol2) legA->AddEntry(f_pol2, "Pol2", "l");
            if( pol2hist ) legA->AddEntry(pol2hist, "Pol2 (eval)", "p");
            legA->Draw();
            // gPad->Modified();
            // gPad->Update();

            // draw a constant y=1.0 line using a pol0 TF1 so it appears as a function (and can be included in legends)
            TF1 *f_const = new TF1(Form("f_const_pt%d", i), "pol0", zTmin, zTmax);
            f_const->SetParameter(0, 1.0);
            f_const->SetLineStyle(2);
            f_const->SetLineWidth(2);
            f_const->SetLineColor(kGray+2);
            f_const->Draw("SAME");
            // delete f_const;
            // TLatex tl;
            // tl.SetNDC();
            // tl.SetTextSize(0.04);
            // tl.DrawLatex(0.15, 0.94, Form("pt bin %d", i + 1));
        }
        std::string multiEffBase = outDir + "/" + outPrefix + "_replot_eff_zT_multi";
        cMultiEff->SaveAs((multiEffBase + ".png").c_str());
        cMultiEff->SaveAs((multiEffBase + ".pdf").c_str());
        delete cMultiEff;
    }

    if (!multi_num_vec.empty() && multi_num_vec.size() == multi_pred_vec.size())
    {
        int nPanels = (int)multi_num_vec.size();
        int nCols = (nPanels <= 3) ? nPanels : (int)std::ceil(std::sqrt(nPanels));
        if (nCols < 1)
            nCols = 1;
        if (nCols > 6)
            nCols = 6;
        int nRows = (int)std::ceil((double)nPanels / nCols);
        int w = std::min(340 * nCols, 2400);
        int h = std::min(320 * nRows, 1800);
        TCanvas *cMultiClos = new TCanvas(Form("c_replot_closure_multi_%s", outPrefix.c_str()), "Closure num vs pred per-pt", w, h);
        cMultiClos->Divide(nCols, nRows, 0.001, 0.001);
        for (int i = 0; i < nPanels; ++i)
        {
            cMultiClos->cd(i + 1);
            gPad->SetRightMargin(0.02);
            gPad->SetLeftMargin(0.12);
            gPad->SetTopMargin(0.02);
            gPad->SetBottomMargin(0.14);
            gPad->SetTickx();
            gPad->SetTicky();
            TH1 *hnum = multi_num_vec[i];
            TH1 *hpred = multi_pred_vec[i];
            if (!hnum || !hpred)
                continue;
            hnum->SetLineColor(kBlack);
            hnum->SetMarkerStyle(20);
            hnum->SetLineWidth(2);
            hpred->SetLineColor(kRed);
            hpred->SetLineWidth(2);
            hpred->SetLineStyle(2);
            hnum->Draw("HIST");
            hpred->Draw("HISTSAME");
            TLatex tl;
            tl.SetNDC();
            tl.SetTextSize(0.04);
            tl.DrawLatex(0.15, 0.94, Form("pt bin %d", i + 1));
        }
        std::string multiClosBase = outDir + "/" + outPrefix + "_replot_closure_zT_multi";
        cMultiClos->SaveAs((multiClosBase + ".png").c_str());
        cMultiClos->SaveAs((multiClosBase + ".pdf").c_str());
        delete cMultiClos;
    }

    if (!multi_diff_vec.empty())
    {
        int nPanels = (int)multi_diff_vec.size();
        int nCols = (nPanels <= 3) ? nPanels : (int)std::ceil(std::sqrt(nPanels));
        if (nCols < 1)
            nCols = 1;
        if (nCols > 6)
            nCols = 6;
        int nRows = (int)std::ceil((double)nPanels / nCols);
        int w = std::min(340 * nCols, 2400);
        int h = std::min(320 * nRows, 1800);
        TCanvas *cMultiRes = new TCanvas(Form("c_replot_res_multi_%s", outPrefix.c_str()), "Residual per-pt bins", w, h);
        cMultiRes->Divide(nCols, nRows, 0.001, 0.001);
        for (int i = 0; i < nPanels; ++i)
        {
            cMultiRes->cd(i + 1);
            gPad->SetRightMargin(0.02);
            gPad->SetLeftMargin(0.12);
            gPad->SetTopMargin(0.02);
            gPad->SetBottomMargin(0.14);
            gPad->SetTickx();
            gPad->SetTicky();
            TH1 *hd = multi_diff_vec[i];
            if (!hd)
                continue;
            hd->SetLineColor(kBlack);
            hd->SetLineWidth(2);
            hd->Draw("HIST");
            TLatex tl;
            tl.SetNDC();
            tl.SetTextSize(0.04);
            tl.DrawLatex(0.15, 0.94, Form("pt bin %d", i + 1));
        }
        std::string multiResBase = outDir + "/" + outPrefix + "_replot_res_zT_multi";
        cMultiRes->SaveAs((multiResBase + ".png").c_str());
        cMultiRes->SaveAs((multiResBase + ".pdf").c_str());
        delete cMultiRes;
    }

    // cleanup collected multi-panel histograms
    for (auto p : multi_eff_vec)
        delete p;
    for (auto p : multi_pol2_hist_vec)
        delete p;
    for (auto p : multi_num_vec)
        delete p;
    for (auto p : multi_pred_vec)
        delete p;
    for (auto p : multi_diff_vec)
        delete p;
    // delete stored fit clones
    for (auto f : multi_f_log_vec)
        if (f)
            delete f;
    for (auto f : multi_f_exp_vec)
        if (f)
            delete f;
    for (auto f : multi_f_pol2_vec)
        if (f)
            delete f;

    fin->Close();
    delete fin;
    std::cout << "Replotting complete. Outputs in: " << outDir << std::endl;
    return 0;
}
