#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TColor.h>
#include <TLatex.h>
#include <TLine.h>
#include <TKey.h>

#include <cmath>
#include <limits>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>

std::string SanitizeLabel(std::string label)
{
    std::replace(label.begin(), label.end(), '.', 'p');
    std::replace(label.begin(), label.end(), '/', '_');
    return label;
}

static std::vector<std::string> SplitCsv(const std::string &csv)
{
    std::vector<std::string> out;
    std::stringstream ss(csv);
    std::string item;
    while (std::getline(ss, item, ','))
    {
        item.erase(0, item.find_first_not_of(" \t\n\r"));
        item.erase(item.find_last_not_of(" \t\n\r") + 1);
        if (!item.empty())
            out.push_back(item);
    }
    return out;
}

TH1D *MakeRatioHistogram(const TH1 *numerator, const TH1 *denominator, const std::string &name, const std::string &yTitle)
{
    if (!numerator || !denominator)
        return nullptr;
    TH1D *ratio = static_cast<TH1D *>(numerator->Clone(name.c_str()));
    ratio->SetDirectory(nullptr);
    ratio->Reset();
    ratio->SetTitle((name + ";z_{T};" + yTitle).c_str());
    const int nb = numerator->GetNbinsX();
    for (int bin = 1; bin <= nb; ++bin)
    {
        const double num = numerator->GetBinContent(bin);
        const double den = denominator->GetBinContent(bin);
        if (den == 0.0)
            continue;
        ratio->SetBinContent(bin, num / den);
        const double denErr = denominator->GetBinError(bin);
        ratio->SetBinError(bin, denErr > 0.0 ? numerator->GetBinError(bin) / den : 0.0);
    }
    return ratio;
}

static std::vector<double> ExtractAxisEdges(const TAxis *axis)
{
    std::vector<double> edges;
    if (!axis)
        return edges;
    const int nb = axis->GetNbins();
    edges.resize(nb + 1);
    for (int i = 1; i <= nb; ++i)
        edges[i - 1] = axis->GetBinLowEdge(i);
    edges[nb] = axis->GetBinUpEdge(nb);
    return edges;
}

// Usage:
// root -l -b -q 'plot_unfold_results.cpp+("d0_unfolded_zt.root","10_15,15_20,20_30,30_100","eta0,eta1,eta2")'
void plot_unfold_results(const std::string &inputFile = "d0_unfolded_zt.root",
                         //  const std::string &jetBinsCsv = "10_15,15_20,20_30,30_100",
                         const std::string &jetBinsCsv = "10_15,15_20,20_30,30_100",
                         const std::string &rapidityCsv = "eta0,eta1,eta2",
                         const std::string &rapidityLabel = "2.5 < y < 3.0, 3.0 < y < 3.5, 3.5 < y < 4.0")
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    std::vector<std::string> jetBins = SplitCsv(jetBinsCsv);
    std::vector<std::string> rapBins = SplitCsv(rapidityCsv);
    std::vector<std::string> rapLabels = SplitCsv(rapidityLabel);

    // automatically use the provided number of jet and rapidity bins
    int ncol = static_cast<int>(jetBins.size());
    int nrow = static_cast<int>(rapBins.size());
    if (ncol <= 0 || nrow <= 0)
    {
        std::cerr << "Error: need at least one jet bin and one rapidity bin." << std::endl;
        return;
    }

    // canvas sizing: scale with number of pads but cap to reasonable max
    int canvasW = std::min(300 * ncol, 2400);
    int canvasH = std::min(300 * nrow, 1800);

    std::unique_ptr<TFile> f(TFile::Open(inputFile.c_str(), "READ"));
    if (!f || f->IsZombie())
    {
        std::cerr << "Failed to open file: " << inputFile << std::endl;
        return;
    }

    TCanvas *c = new TCanvas("c_unfold_summary", "Unfold summary", canvasW, canvasH);
    c->Divide(ncol, nrow);
    // colors
    const int colMeasured = kBlack;
    const int colUnfolded = kBlue + 1;
    const int colRefolded = kGreen + 2;
    const int colTruth = kRed + 1;

    for (int iy = 0; iy < nrow; ++iy)
    {
        // rapidity: top->bottom increasing, so y index maps to rapBins[iy]
        const std::string &rapLabel = rapBins[iy];
        const std::string displayRap = (rapLabels.size() > static_cast<size_t>(iy) ? rapLabels[iy] : rapLabel);
        for (int ix = 0; ix < ncol; ++ix)
        {
            const std::string &jetLabel = jetBins[ix];
            const int pad = iy * ncol + ix + 1; // ROOT pads number left->right top->bottom
            c->cd(pad)->SetLeftMargin(0.12);
            c->cd(pad)->SetRightMargin(0.01);
            c->cd(pad)->SetTopMargin(0.01);
            c->cd(pad)->SetBottomMargin(0.08);
            gPad->SetTickx(1);
            gPad->SetTicky(1);

            const std::string measuredName = "measured_zT_" + jetLabel + "_" + rapLabel;
            const std::string unfoldedName = "best_unfolded_zT_" + jetLabel + "_" + rapLabel;
            const std::string truthName = "truth_zT_" + jetLabel + "_" + rapLabel;

            TH1 *hmeas = dynamic_cast<TH1 *>(f->Get(measuredName.c_str()));
            TH1 *hunf = dynamic_cast<TH1 *>(f->Get(unfoldedName.c_str()));
            TH1 *htruth = dynamic_cast<TH1 *>(f->Get(truthName.c_str()));

            // find refolded with highest iteration if available
            TH1 *href = nullptr;
            for (int iter = 10; iter >= 1; --iter)
            {
                std::ostringstream rn;
                rn << "refolded_zT_" << jetLabel << "_" << rapLabel << "_iter" << iter;
                TH1 *h = dynamic_cast<TH1 *>(f->Get(rn.str().c_str()));
                if (h)
                {
                    href = h;
                    break;
                }
            }

            if (!hmeas && !hunf && !htruth && !href)
            {
                // empty pad
                continue;
            }

            // determine y-range from measured/truth/unfolded
            double yMax = 0.0;
            auto consider = [&](TH1 *h)
            { if (h) yMax = std::max(yMax, h->GetMaximum()); };
            consider(hmeas);
            consider(hunf);
            consider(htruth);
            consider(href);
            if (yMax <= 0)
                yMax = 1.0;

            // unfolded
            if (hunf)
            {
                hunf->SetTitle(";z_{T};counts");
                hunf->SetLineColor(colUnfolded);
                hunf->SetMarkerColor(colUnfolded);
                hunf->SetMarkerStyle(20);
                hunf->SetMaximum(yMax * 1.3);
                hunf->Draw("E1");
            }
            // draw measured first as points
            if (hmeas)
            {
                hmeas->SetMarkerStyle(24);
                hmeas->SetMarkerColor(colMeasured);
                hmeas->SetLineColor(colMeasured);
                hmeas->SetTitle(";z_{T};counts");
                hmeas->Draw(hunf ? "E1 SAME" : "E1");
            }

            // refolded
            if (href)
            {
                href->SetLineColor(colRefolded);
                href->SetLineStyle(7);
                href->SetMarkerColor(colRefolded);
                href->SetMarkerStyle(25);
                href->Draw("E1 SAME");
            }

            // truth (rescaled for shape comparison if needed)
            if (htruth)
            {
                htruth->SetLineColor(colTruth);
                htruth->SetMarkerColor(colTruth);
                htruth->SetMarkerStyle(21);
                // scale truth to unfolded integral for shape comparison
                if (hunf && hunf->Integral() > 0 && htruth->Integral() > 0)
                {
                    double scale = hunf->Integral() / htruth->Integral();
                    htruth->SetDirectory(nullptr);
                    htruth = static_cast<TH1 *>(htruth->Clone((std::string(htruth->GetName()) + "_scaled").c_str()));
                    htruth->Scale(scale);
                }
                htruth->Draw("E1 SAME");
            }

            // legend
            if (pad == 1) // only draw legend on first pad
            {
                std::cout << "Drawing legend for pad " << pad << std::endl;
                TLegend *legend = new TLegend(0.6, 0.7, 0.95, 0.9);
                legend->SetFillStyle(0);
                legend->SetBorderSize(0);
                if (hmeas)
                    legend->AddEntry(hmeas, "Measured", "lep");
                if (hunf)
                    legend->AddEntry(hunf, "Unfolded", "lep");
                if (href)
                    legend->AddEntry(href, "Refolded", "lep");
                if (htruth)
                    legend->AddEntry(htruth, "Truth (shape)", "lep");
                legend->Draw();
            }

            // title text
            TLatex tex;
            tex.SetNDC();
            tex.SetTextSize(0.04);
            std::ostringstream title;
            title << "jet p_{T}: " << jetLabel << " GeV/c   rapidity: " << displayRap;
            tex.DrawLatex(0.22, 0.92, title.str().c_str());
        }
    }

    c->Modified();
    c->Update();
    // Save output in same directory as the input file
    std::string outDir = ".";
    const size_t slashPos = inputFile.find_last_of("/\\");
    if (slashPos != std::string::npos)
    {
        outDir = inputFile.substr(0, slashPos);
    }
    // detect dataset tag in input filename and append to output filenames
    std::string datasetTag = "";
    if (inputFile.find("pPb") != std::string::npos)
        datasetTag = "_pPb";
    else if (inputFile.find("Pbp") != std::string::npos)
        datasetTag = "_Pbp";

    // human-friendly dataset label (no leading underscore)
    std::string datasetLabel = "";
    if (inputFile.find("pPb") != std::string::npos)
        datasetLabel = "pPb";
    else if (inputFile.find("Pbp") != std::string::npos)
        datasetLabel = "Pbp";

    // helper: draw dataset + experiment + sqrt(s_NN) label on the 3rd pad of a canvas (if present)
    auto drawDatasetLabelOnCanvas = [&](TCanvas *can)
    {
        if (!can)
            return;
        const int npads = ncol * nrow;
        if (npads < 3)
            return; // no third pad available
        can->cd(3);
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.01);
        gPad->SetBottomMargin(0.08);
        TLatex lab;
        lab.SetNDC();
        lab.SetTextFont(42);
        // compute a right-aligned X coordinate using pad right margin and a small offset
        const double padRight = gPad->GetRightMargin();
        const double x = 1.0 - padRight - 0.06;
        double y = 0.875;
        // Right-align text (horizontal:3, vertical:1 => 31 -> right/top)
        lab.SetTextAlign(31);
        // Draw experiment on first line (top-right)
        lab.SetTextSize(0.04);
        lab.DrawLatex(x, y, "LHC#it{b}");
        // Draw dataset + center-of-mass energy on next line (right-aligned)
        lab.SetTextSize(0.04);
        std::string energy = "#it{s}_{NN} = 8.16 TeV";
        if (!datasetLabel.empty())
            energy = datasetLabel + std::string(", ") + energy;
        lab.DrawLatex(x, y - 0.045, energy.c_str());
        can->Modified();
        can->Update();
    };

    // draw label on the main summary canvas (third pad)
    drawDatasetLabelOnCanvas(c);

    // const std::string outPath = outDir + "/unfold_summary" + datasetTag + ".png";
    // c->SaveAs(outPath.c_str());
    // const std::string outPathpdf = outDir + "/unfold_summary" + datasetTag + ".pdf";
    // c->SaveAs(outPathpdf.c_str());

    // Second summary using unfolded result with total uncertainties (if available)
    TCanvas *c2 = new TCanvas("c_unfold_summary_totalErr", "Unfold summary (total errors)", canvasW, canvasH);
    c2->Divide(ncol, nrow);
    for (int iy = 0; iy < nrow; ++iy)
    {
        const std::string &rapLabel = rapBins[iy];
        const std::string displayRap = (rapLabels.size() > static_cast<size_t>(iy) ? rapLabels[iy] : rapLabel);
        for (int ix = 0; ix < ncol; ++ix)
        {
            const std::string &jetLabel = jetBins[ix];
            const int pad = iy * ncol + ix + 1;
            c2->cd(pad)->SetLeftMargin(0.12);
            c2->cd(pad)->SetRightMargin(0.01);
            c2->cd(pad)->SetTopMargin(0.01);
            c2->cd(pad)->SetBottomMargin(0.08);
            gPad->SetTickx(1);
            gPad->SetTicky(1);

            const std::string measuredName = "measured_zT_" + jetLabel + "_" + rapLabel;
            const std::string unfoldedName = "best_unfolded_with_totalErr_" + SanitizeLabel(jetLabel) + "_" + rapLabel;
            const std::string fallbackUnfold = "best_unfolded_zT_" + jetLabel + "_" + rapLabel;
            const std::string truthName = "truth_zT_" + jetLabel + "_" + rapLabel;

            TH1 *hmeas = dynamic_cast<TH1 *>(f->Get(measuredName.c_str()));
            TH1 *hunf = dynamic_cast<TH1 *>(f->Get(unfoldedName.c_str()));
            if (!hunf)
                hunf = dynamic_cast<TH1 *>(f->Get(fallbackUnfold.c_str()));
            TH1 *htruth = dynamic_cast<TH1 *>(f->Get(truthName.c_str()));

            // find refolded with highest iteration if available
            TH1 *href = nullptr;
            for (int iter = 10; iter >= 1; --iter)
            {
                std::ostringstream rn;
                rn << "refolded_zT_" << jetLabel << "_" << rapLabel << "_iter" << iter;
                TH1 *h = dynamic_cast<TH1 *>(f->Get(rn.str().c_str()));
                if (h)
                {
                    href = h;
                    break;
                }
            }

            if (!hmeas && !hunf && !htruth && !href)
            {
                continue;
            }

            double yMax = 0.0;
            auto consider = [&](TH1 *h)
            { if (h) yMax = std::max(yMax, h->GetMaximum()); };
            consider(hmeas);
            consider(hunf);
            consider(htruth);
            consider(href);
            if (yMax <= 0)
                yMax = 1.0;

            if (hunf)
            {
                hunf->SetTitle(";z_{T};counts");
                hunf->SetLineColor(colUnfolded);
                hunf->SetMarkerColor(colUnfolded);
                hunf->SetMarkerStyle(20);
                hunf->SetMaximum(yMax * 1.3);
                hunf->Draw("E1");
            }
            if (hmeas)
            {
                hmeas->SetMarkerStyle(24);
                hmeas->SetMarkerColor(colMeasured);
                hmeas->SetLineColor(colMeasured);
                hmeas->SetTitle(";z_{T};counts");
                hmeas->Draw(hunf ? "E1 SAME" : "E1");
            }
            if (href)
            {
                href->SetLineColor(colRefolded);
                href->SetLineStyle(7);
                href->SetMarkerColor(colRefolded);
                href->SetMarkerStyle(25);
                href->Draw("E1 SAME");
            }
            if (htruth)
            {
                htruth->SetLineColor(colTruth);
                htruth->SetMarkerColor(colTruth);
                htruth->SetMarkerStyle(21);
                if (hunf && hunf->Integral() > 0 && htruth->Integral() > 0)
                {
                    double scale = hunf->Integral() / htruth->Integral();
                    htruth->SetDirectory(nullptr);
                    htruth = static_cast<TH1 *>(htruth->Clone((std::string(htruth->GetName()) + "_scaled").c_str()));
                    htruth->Scale(scale);
                }
                htruth->Draw("E1 SAME");
            }

            if (pad == 1)
            {
                TLegend *legend = new TLegend(0.6, 0.7, 0.95, 0.9);
                legend->SetFillStyle(0);
                legend->SetBorderSize(0);
                if (hmeas)
                    legend->AddEntry(hmeas, "Measured", "lep");
                if (hunf)
                    legend->AddEntry(hunf, "Unfolded (tot err)", "lep");
                if (href)
                    legend->AddEntry(href, "Refolded", "lep");
                if (htruth)
                    legend->AddEntry(htruth, "Truth (shape)", "lep");
                legend->Draw();
            }

            TLatex tex2;
            tex2.SetNDC();
            tex2.SetTextSize(0.04);
            std::ostringstream title2;
            title2 << "jet p_{T}: " << jetLabel << " GeV/c   rapidity: " << displayRap;
            tex2.DrawLatex(0.22, 0.92, title2.str().c_str());
        }
    }
    c2->Modified();
    c2->Update();
    // annotate dataset / experiment / energy on pad 3
    drawDatasetLabelOnCanvas(c2);
    const std::string outPath2 = outDir + "/unfold_summary_with_totalErr" + datasetTag + ".png";
    c2->SaveAs(outPath2.c_str());
    const std::string outPath2pdf = outDir + "/unfold_summary_with_totalErr" + datasetTag + ".pdf";
    c2->SaveAs(outPath2pdf.c_str());

    // Ratio summary canvas (measured / unfolded and measured / refolded)
    TCanvas *cr = new TCanvas("c_unfold_ratio", "Unfold ratios", canvasW, canvasH);
    cr->Divide(ncol, nrow);
    for (int iy = 0; iy < nrow; ++iy)
    {
        const std::string &rapLabel = rapBins[iy];
        const std::string displayRap = (rapLabels.size() > static_cast<size_t>(iy) ? rapLabels[iy] : rapLabel);
        for (int ix = 0; ix < ncol; ++ix)
        {
            const std::string &jetLabel = jetBins[ix];
            const int pad = iy * ncol + ix + 1;
            cr->cd(pad)->SetLeftMargin(0.12);
            cr->cd(pad)->SetRightMargin(0.01);
            cr->cd(pad)->SetTopMargin(0.01);
            cr->cd(pad)->SetBottomMargin(0.08);
            gPad->SetTickx(1);
            gPad->SetTicky(1);

            const std::string measuredName = "measured_zT_" + jetLabel + "_" + rapLabel;
            const std::string unfoldedName = "best_unfolded_zT_" + jetLabel + "_" + rapLabel;

            TH1 *hmeas = dynamic_cast<TH1 *>(f->Get(measuredName.c_str()));
            TH1 *hunf = dynamic_cast<TH1 *>(f->Get(unfoldedName.c_str()));

            // find refolded
            TH1 *href = nullptr;
            for (int iter = 10; iter >= 1; --iter)
            {
                std::ostringstream rn;
                rn << "refolded_zT_" << jetLabel << "_" << rapLabel << "_iter" << iter;
                TH1 *h = dynamic_cast<TH1 *>(f->Get(rn.str().c_str()));
                if (h)
                {
                    href = h;
                    break;
                }
            }

            if (!hmeas)
                continue;

            TH1 *hr1 = nullptr;
            TH1 *hr2 = nullptr;
            if (hunf && hunf->GetNbinsX() == hmeas->GetNbinsX())
            {
                hr1 = static_cast<TH1 *>(hmeas->Clone((std::string("ratio_meas_unf_") + jetLabel + "_" + rapLabel).c_str()));
                hr1->SetDirectory(nullptr);
                hr1->Divide(hunf);
                hr1->SetMarkerStyle(20);
                hr1->SetMarkerColor(kBlack);
                hr1->SetLineColor(kBlack);
            }
            if (href && href->GetNbinsX() == hmeas->GetNbinsX())
            {
                hr2 = static_cast<TH1 *>(hmeas->Clone((std::string("ratio_meas_ref_") + jetLabel + "_" + rapLabel).c_str()));
                hr2->SetDirectory(nullptr);
                hr2->Divide(href);
                hr2->SetMarkerStyle(21);
                hr2->SetMarkerColor(kGreen + 2);
                hr2->SetLineColor(kGreen + 2);
            }

            // draw frame
            TH1 *hframe = nullptr;
            if (hr1)
                hframe = static_cast<TH1 *>(hr1->Clone("frame_ratio"));
            else if (hr2)
                hframe = static_cast<TH1 *>(hr2->Clone("frame_ratio"));
            if (!hframe)
                continue;
            hframe->SetDirectory(nullptr);
            hframe->SetTitle(";z_{T};Measured / Model");
            hframe->GetYaxis()->SetRangeUser(0.5, 1.5);
            hframe->Draw("AXIS");
            if (hr1)
                hr1->Draw("E1 SAME");
            if (hr2)
                hr2->Draw("E1 SAME");
            TLine l;
            l.SetLineStyle(7);
            l.SetLineColor(kBlack);
            l.DrawLine(hframe->GetXaxis()->GetXmin(), 1.0, hframe->GetXaxis()->GetXmax(), 1.0);

            TLatex t;
            t.SetNDC();
            t.SetTextSize(0.04);
            std::ostringstream ttl;
            ttl << "jet p_{T}: " << jetLabel << " GeV/c   rapidity: " << displayRap;
            t.DrawLatex(0.22, 0.92, ttl.str().c_str());

            // legend for ratio plot (draw on first pad only)
            if (pad == 1)
            {
                TLegend *lgr = new TLegend(0.4, 0.75, 0.95, 0.9);
                lgr->SetFillStyle(0);
                lgr->SetBorderSize(0);
                if (hr1)
                    lgr->AddEntry(hr1, "Measured / Unfolded", "lep");
                if (hr2)
                    lgr->AddEntry(hr2, "Measured / Refolded", "lep");
                lgr->Draw();
            }
        }
    }
    cr->Modified();
    cr->Update();
    // annotate dataset / experiment / energy on pad 3
    drawDatasetLabelOnCanvas(cr);
    // const std::string outPathR = outDir + "/unfold_ratio_summary" + datasetTag + ".png";
    // cr->SaveAs(outPathR.c_str());
    // const std::string outPathRpdf = outDir + "/unfold_ratio_summary" + datasetTag + ".pdf";
    // cr->SaveAs(outPathRpdf.c_str());

    // Ratio summary for total-uncertainty unfolded (if present)
    TCanvas *cr2 = new TCanvas("c_unfold_ratio_totalErr", "Unfold ratios (total err)", canvasW, canvasH);
    cr2->Divide(ncol, nrow);
    for (int iy = 0; iy < nrow; ++iy)
    {
        const std::string &rapLabel = rapBins[iy];
        const std::string displayRap = (rapLabels.size() > static_cast<size_t>(iy) ? rapLabels[iy] : rapLabel);
        for (int ix = 0; ix < ncol; ++ix)
        {
            const std::string &jetLabel = jetBins[ix];
            const int pad = iy * ncol + ix + 1;
            cr2->cd(pad)->SetLeftMargin(0.12);
            cr2->cd(pad)->SetRightMargin(0.01);
            cr2->cd(pad)->SetTopMargin(0.01);
            cr2->cd(pad)->SetBottomMargin(0.08);
            gPad->SetTickx(1);
            gPad->SetTicky(1);

            const std::string measuredName = "measured_zT_" + jetLabel + "_" + rapLabel;
            const std::string unfoldedNameTotal = "best_unfolded_with_totalErr_" + SanitizeLabel(jetLabel) + "_" + rapLabel;
            const std::string unfoldedName = "best_unfolded_zT_" + jetLabel + "_" + rapLabel;

            TH1 *hmeas = dynamic_cast<TH1 *>(f->Get(measuredName.c_str()));
            TH1 *hunf = dynamic_cast<TH1 *>(f->Get(unfoldedNameTotal.c_str()));
            if (!hunf)
                hunf = dynamic_cast<TH1 *>(f->Get(unfoldedName.c_str()));

            // find refolded
            TH1 *href = nullptr;
            for (int iter = 10; iter >= 1; --iter)
            {
                std::ostringstream rn;
                rn << "refolded_zT_" << jetLabel << "_" << rapLabel << "_iter" << iter;
                TH1 *h = dynamic_cast<TH1 *>(f->Get(rn.str().c_str()));
                if (h)
                {
                    href = h;
                    break;
                }
            }

            if (!hmeas)
                continue;

            TH1 *hr1 = nullptr;
            TH1 *hr2 = nullptr;
            if (hunf && hunf->GetNbinsX() == hmeas->GetNbinsX())
            {
                hr1 = static_cast<TH1 *>(hmeas->Clone((std::string("ratio_meas_unf_tot_") + jetLabel + "_" + rapLabel).c_str()));
                hr1->SetDirectory(nullptr);
                hr1->Divide(hunf);
                hr1->SetMarkerStyle(20);
                hr1->SetMarkerColor(kBlack);
                hr1->SetLineColor(kBlack);
            }
            if (href && href->GetNbinsX() == hmeas->GetNbinsX())
            {
                hr2 = static_cast<TH1 *>(hmeas->Clone((std::string("ratio_meas_ref_tot_") + jetLabel + "_" + rapLabel).c_str()));
                hr2->SetDirectory(nullptr);
                hr2->Divide(href);
                hr2->SetMarkerStyle(21);
                hr2->SetMarkerColor(kGreen + 2);
                hr2->SetLineColor(kGreen + 2);
            }

            TH1 *hframe = nullptr;
            if (hr1)
                hframe = static_cast<TH1 *>(hr1->Clone("frame_ratio2"));
            else if (hr2)
                hframe = static_cast<TH1 *>(hr2->Clone("frame_ratio2"));
            if (!hframe)
                continue;
            hframe->SetDirectory(nullptr);
            hframe->SetTitle(";z_{T};Measured / Model");
            hframe->GetYaxis()->SetRangeUser(0.5, 1.5);
            hframe->Draw("AXIS");
            if (hr1)
                hr1->Draw("E1 SAME");
            if (hr2)
                hr2->Draw("E1 SAME");
            TLine l;
            l.SetLineStyle(7);
            l.SetLineColor(kBlack);
            l.DrawLine(hframe->GetXaxis()->GetXmin(), 1.0, hframe->GetXaxis()->GetXmax(), 1.0);

            TLatex t;
            t.SetNDC();
            t.SetTextSize(0.04);
            std::ostringstream ttl;
            ttl << "jet p_{T}: " << jetLabel << " GeV/c   rapidity: " << displayRap;
            t.DrawLatex(0.22, 0.92, ttl.str().c_str());

            // legend for ratio plot (draw on first pad only)
            if (pad == 1)
            {
                TLegend *lgrb = new TLegend(0.4, 0.75, 0.95, 0.9);
                lgrb->SetFillStyle(0);
                lgrb->SetBorderSize(0);
                if (hr1)
                    lgrb->AddEntry(hr1, "Measured / Unfolded", "lep");
                if (hr2)
                    lgrb->AddEntry(hr2, "Measured / Refolded", "lep");
                lgrb->Draw();
            }
        }
    }
    cr2->Modified();
    cr2->Update();
    // annotate dataset / experiment / energy on pad 3
    drawDatasetLabelOnCanvas(cr2);
    const std::string outPathR2 = outDir + "/unfold_ratio_summary_with_totalErr" + datasetTag + ".png";
    cr2->SaveAs(outPathR2.c_str());
    const std::string outPathR2pdf = outDir + "/unfold_ratio_summary_with_totalErr" + datasetTag + ".pdf";
    cr2->SaveAs(outPathR2pdf.c_str());

    // Prior-variation summary: overlay rFlat and rW ratio histograms in a multi-panel canvas
    TCanvas *cp = new TCanvas("c_prior_summary", "Prior variations (prior / nominal)", canvasW, canvasH);
    cp->Divide(ncol, nrow);
    // helper: try multiple candidate names to find histogram (robust against label sanitization/order differences)
    auto findHist = [&](const std::vector<std::string> &cands) -> TH1 *
    {
        for (const auto &n : cands)
        {
            TH1 *h = dynamic_cast<TH1 *>(f->Get(n.c_str()));
            if (h)
                return h;
        }
        return nullptr;
    };

    for (int iy = 0; iy < nrow; ++iy)
    {
        const std::string &rapLabel = rapBins[iy];
        const std::string displayRap = (rapLabels.size() > static_cast<size_t>(iy) ? rapLabels[iy] : rapLabel);
        for (int ix = 0; ix < ncol; ++ix)
        {
            const std::string &jetLabel = jetBins[ix];
            const int pad = iy * ncol + ix + 1;
            cp->cd(pad)->SetLeftMargin(0.12);
            cp->cd(pad)->SetRightMargin(0.01);
            cp->cd(pad)->SetTopMargin(0.01);
            cp->cd(pad)->SetBottomMargin(0.08);
            gPad->SetTickx(1);
            gPad->SetTicky(1);

            // candidate name permutations
            std::vector<std::string> flatCands = {
                std::string("ratio_prior_flat_") + SanitizeLabel(jetLabel) + "_" + rapLabel,
                std::string("ratio_prior_flat_") + jetLabel + "_" + rapLabel,
                std::string("ratio_prior_flat_") + SanitizeLabel(jetLabel) + "_" + SanitizeLabel(rapLabel)};
            std::vector<std::string> wCands = {
                std::string("ratio_prior_weighted_") + SanitizeLabel(jetLabel) + "_" + rapLabel,
                std::string("ratio_prior_weighted_") + jetLabel + "_" + rapLabel,
                std::string("ratio_prior_weighted_") + SanitizeLabel(jetLabel) + "_" + SanitizeLabel(rapLabel)};

            // Diagnostic: list attempts and which candidates exist
            std::string usedFlatName;
            for (const auto &n : flatCands)
            {
                TH1 *h = dynamic_cast<TH1 *>(f->Get(n.c_str()));
                // std::cout << "[prior_summary] flat candidate: '" << n << "' -> " << (h ? "FOUND" : "missing") << std::endl;
                if (h && usedFlatName.empty())
                    usedFlatName = n;
            }
            std::string usedWName;
            for (const auto &n : wCands)
            {
                TH1 *h = dynamic_cast<TH1 *>(f->Get(n.c_str()));
                // std::cout << "[prior_summary] weighted candidate: '" << n << "' -> " << (h ? "FOUND" : "missing") << std::endl;
                if (h && usedWName.empty())
                    usedWName = n;
            }

            TH1 *hFlat = usedFlatName.empty() ? nullptr : dynamic_cast<TH1 *>(f->Get(usedFlatName.c_str()));
            TH1 *hW = usedWName.empty() ? nullptr : dynamic_cast<TH1 *>(f->Get(usedWName.c_str()));

            if (!hFlat && !hW)
            {
                std::cout << "[prior_summary] no prior ratios found for " << jetLabel << " / " << rapLabel << std::endl;
                continue;
            }

            // More diagnostics: print basic histogram stats for any found histogram
            auto printHistInfo = [&](TH1 *h, const std::string &label, const std::string &name)
            {
                if (!h)
                    return;
                const int nbins = h->GetNbinsX();
                int nonzero = 0;
                double minv = std::numeric_limits<double>::infinity();
                double maxv = -std::numeric_limits<double>::infinity();
                for (int b = 1; b <= nbins; ++b)
                {
                    double c = h->GetBinContent(b);
                    if (std::isfinite(c))
                    {
                        if (c != 0.0)
                            ++nonzero;
                        minv = std::min(minv, c);
                        maxv = std::max(maxv, c);
                    }
                }
            };

            printHistInfo(hFlat, "flat", usedFlatName.empty() ? std::string("(none)") : usedFlatName);
            printHistInfo(hW, "weighted", usedWName.empty() ? std::string("(none)") : usedWName);

            if (hFlat)
            {
                hFlat->SetLineColor(kRed + 1);
                hFlat->SetMarkerColor(kRed + 1);
                hFlat->SetMarkerStyle(21);
            }
            if (hW)
            {
                hW->SetLineColor(kGreen + 2);
                hW->SetMarkerColor(kGreen + 2);
                hW->SetMarkerStyle(20);
            }

            TH1D *frame = nullptr;
            if (hFlat)
                frame = static_cast<TH1D *>(hFlat->Clone((std::string("frame_prior_") + jetLabel + "_" + rapLabel).c_str()));
            else if (hW)
                frame = static_cast<TH1D *>(hW->Clone((std::string("frame_prior_") + jetLabel + "_" + rapLabel).c_str()));
            if (!frame)
                continue;
            frame->SetDirectory(nullptr);
            frame->SetTitle(";z_{T};prior / nominal");
            frame->GetYaxis()->SetRangeUser(0.5, 1.5);
            frame->Draw("AXIS");

            // compute per-bin systematic band from available prior variations (RMS of deviations from 1.0)
            {
                TH1D *sysBandP = (TH1D *)(frame->Clone((std::string("sys_prior_band_") + jetLabel + "_" + rapLabel).c_str()));
                const int nb = frame->GetNbinsX();
                int nonzeroFlat = 0, nonzeroW = 0, nonzeroSys = 0;
                double sigmaSum = 0.0, sigmaMax = 0.0;
                double sigmaMin = std::numeric_limits<double>::infinity();
                for (int b = 1; b <= nb; ++b)
                {
                    int count = 0;
                    double sumsq = 0.0;
                    if (hFlat)
                    {
                        double dv = hFlat->GetBinContent(b) - 1.0;
                        if (std::isfinite(dv))
                        {
                            sumsq += dv * dv;
                            ++count;
                        }
                        double fc = hFlat->GetBinContent(b);
                        if (std::isfinite(fc) && fc != 0.0)
                            ++nonzeroFlat;
                    }
                    if (hW)
                    {
                        double dv = hW->GetBinContent(b) - 1.0;
                        if (std::isfinite(dv))
                        {
                            sumsq += dv * dv;
                            ++count;
                        }
                        double wc = hW->GetBinContent(b);
                        if (std::isfinite(wc) && wc != 0.0)
                            ++nonzeroW;
                    }
                    double sigma = (count > 0) ? std::sqrt(sumsq / count) : 0.0;
                    if (sysBandP)
                    {
                        sysBandP->SetBinContent(b, 1.0);
                        sysBandP->SetBinError(b, sigma);
                    }
                    if (sigma > 0.0 && std::isfinite(sigma))
                    {
                        ++nonzeroSys;
                        sigmaSum += sigma;
                        sigmaMax = std::max(sigmaMax, sigma);
                        sigmaMin = std::min(sigmaMin, sigma);
                    }
                }
                double sigmaMean = (nonzeroSys > 0) ? sigmaSum / nonzeroSys : 0.0;
                if (sigmaMin == std::numeric_limits<double>::infinity())
                    sigmaMin = 0.0;
                if (sysBandP)
                {
                    sysBandP->SetFillColor(kGray + 2);
                    sysBandP->SetMarkerStyle(1);
                    sysBandP->SetMarkerColor(kGray + 2);
                    sysBandP->SetFillStyle(3844);
                    sysBandP->SetLineColor(kGray + 2);
                    sysBandP->Draw("E3 SAME");
                }
            }

            if (hFlat)
            {
                hFlat->Draw("pe SAME");
            }
            else
            {
                std::cout << "[prior_summary] no flat prior ratio found for " << jetLabel << " / " << rapLabel << std::endl;
            }
            if (hW)
            {
                hW->Draw("pe SAME");
            }
            else
            {
                std::cout << "[prior_summary] no weighted prior ratio found for " << jetLabel << " / " << rapLabel << std::endl;
            }

            if (pad == 1)
            {
                TLegend *legp = new TLegend(0.35, 0.72, 0.88, 0.88);
                legp->SetBorderSize(0);
                legp->SetFillStyle(0);
                if (hFlat)
                    legp->AddEntry(hFlat, "Prior: flat / nominal", "lep");
                if (hW)
                    legp->AddEntry(hW, "Prior: weighted / nominal", "lep");
                legp->Draw();
            }
            TLine *unityLinePrior = new TLine(frame->GetXaxis()->GetXmin(), 1.0, frame->GetXaxis()->GetXmax(), 1.0);
            unityLinePrior->SetLineStyle(7);
            unityLinePrior->SetLineColor(kBlack);
            unityLinePrior->Draw();

            TLatex tp;
            tp.SetNDC();
            tp.SetTextSize(0.04);
            std::ostringstream ttlp;
            ttlp << "jet p_{T}: " << jetLabel << " GeV/c   rapidity: " << displayRap;
            tp.DrawLatex(0.22, 0.92, ttlp.str().c_str());
        }
    }
    cp->Modified();
    cp->Update();
    // annotate dataset / experiment / energy on pad 3
    drawDatasetLabelOnCanvas(cp);
    const std::string outPathP = outDir + "/prior_variation_summary" + datasetTag + ".png";
    cp->SaveAs(outPathP.c_str());
    const std::string outPathPpdf = outDir + "/prior_variation_summary" + datasetTag + ".pdf";
    cp->SaveAs(outPathPpdf.c_str());

    // Global multi-panel comparison of (base-1) / base / (base+1) iterations
    // Determine highest available iteration across all keys once (cap at 7)
    int globalMaxIter = 0;
    for (int it = 7; it >= 1; --it)
    {
        TIter next(f->GetListOfKeys());
        TKey *k;
        while ((k = (TKey *)next()))
        {
            const char *kn = k->GetName();
            if (!kn)
                continue;
            std::string sname(kn);
            std::string suf = std::string("_iter") + std::to_string(it);
            if (sname.find(suf) != std::string::npos)
            {
                globalMaxIter = it;
                break;
            }
        }
        if (globalMaxIter)
            break;
    }
    if (globalMaxIter == 0)
    {
        std::cout << "[plot_unfold] no iter keys found; falling back to 7" << std::endl;
        globalMaxIter = 7;
    }
    int globalBaseIter = std::max(1, globalMaxIter - 1);
    std::cout << "[plot_unfold] globalMaxIter=" << globalMaxIter << " globalBaseIter=" << globalBaseIter << std::endl;

    // Overlay counts: base-1 / base / base+1
    TCanvas *citer = new TCanvas("c_iter_plusminus", "Iteration +/-1 comparison (global base)", canvasW, canvasH);
    citer->Divide(ncol, nrow);
    for (int iy = 0; iy < nrow; ++iy)
    {
        const std::string &rapLabel = rapBins[iy];
        const std::string displayRap = (rapLabels.size() > static_cast<size_t>(iy) ? rapLabels[iy] : rapLabel);
        for (int ix = 0; ix < ncol; ++ix)
        {
            const std::string &jetLabel = jetBins[ix];
            const int pad = iy * ncol + ix + 1;
            citer->cd(pad)->SetLeftMargin(0.12);
            citer->cd(pad)->SetRightMargin(0.01);
            citer->cd(pad)->SetTopMargin(0.01);
            citer->cd(pad)->SetBottomMargin(0.08);
            gPad->SetTickx(1);
            gPad->SetTicky(1);

            auto makeName = [&](int it)
            {
                std::ostringstream s;
                s << "unfolded_zT_" << jetLabel << "_" << rapLabel << "_iter" << it;
                return s.str();
            };

            TH1 *cPrev = (globalBaseIter > 1) ? dynamic_cast<TH1 *>(f->Get(makeName(globalBaseIter - 1).c_str())) : nullptr;
            TH1 *cBase = dynamic_cast<TH1 *>(f->Get(makeName(globalBaseIter).c_str()));
            TH1 *cNext = dynamic_cast<TH1 *>(f->Get(makeName(globalBaseIter + 1).c_str()));

            if (!cPrev && !cBase && !cNext)
                continue;

            double yMax = 0.0;
            auto consider = [&](TH1 *h)
            { if (h) yMax = std::max(yMax, h->GetMaximum()); };
            consider(cPrev);
            consider(cBase);
            consider(cNext);
            if (yMax <= 0)
                yMax = 1.0;

            if (cBase)
            {
                cBase->SetTitle(";z_{T};counts");
                cBase->SetLineColor(kBlue + 1);
                cBase->SetMarkerColor(kBlue + 1);
                cBase->SetMarkerStyle(20);
                cBase->SetMaximum(yMax * 1.3);
                cBase->Draw("E1");
            }
            if (cPrev)
            {
                cPrev->SetLineColor(kRed + 1);
                cPrev->SetMarkerColor(kRed + 1);
                cPrev->SetMarkerStyle(21);
                cPrev->Draw(cBase ? "E1 SAME" : "E1");
            }
            if (cNext)
            {
                cNext->SetLineColor(kGreen + 2);
                cNext->SetMarkerColor(kGreen + 2);
                cNext->SetMarkerStyle(24);
                cNext->Draw((cBase || cPrev) ? "E1 SAME" : "E1");
            }

            if (pad == 1)
            {
                TLegend *leg = new TLegend(0.6, 0.7, 0.95, 0.9);
                leg->SetFillStyle(0);
                leg->SetBorderSize(0);
                if (cBase)
                    leg->AddEntry(cBase, (std::string("iter") + std::to_string(globalBaseIter) + " (default)").c_str(), "lep");
                if (cPrev)
                    leg->AddEntry(cPrev, (std::string("iter") + std::to_string(globalBaseIter - 1)).c_str(), "lep");
                if (cNext)
                    leg->AddEntry(cNext, (std::string("iter") + std::to_string(globalBaseIter + 1)).c_str(), "lep");
                leg->Draw();
            }
            TLatex t;
            t.SetNDC();
            t.SetTextSize(0.04);
            std::ostringstream ttl;
            ttl << "jet p_{T}: " << jetLabel << " GeV/c   rapidity: " << displayRap;
            t.DrawLatex(0.22, 0.92, ttl.str().c_str());
        }
    }
    citer->Modified();
    citer->Update();
    // annotate dataset / experiment / energy on pad 3
    drawDatasetLabelOnCanvas(citer);
    const std::string outPathIter = outDir + "/iter_plusminus_summary" + datasetTag + ".png";
    citer->SaveAs(outPathIter.c_str());
    const std::string outPathIterPdf = outDir + "/iter_plusminus_summary" + datasetTag + ".pdf";
    citer->SaveAs(outPathIterPdf.c_str());

    // Ratio multi-panel: draw (prev / base) and (next / base) where available
    TCanvas *crIter = new TCanvas("c_iter_plusminus_ratio", "Iteration +/-1 ratios (global base)", canvasW, canvasH);
    crIter->Divide(ncol, nrow);
    for (int iy = 0; iy < nrow; ++iy)
    {
        const std::string &rapLabel = rapBins[iy];
        const std::string displayRap = (rapLabels.size() > static_cast<size_t>(iy) ? rapLabels[iy] : rapLabel);
        for (int ix = 0; ix < ncol; ++ix)
        {
            const std::string &jetLabel = jetBins[ix];
            const int pad = iy * ncol + ix + 1;
            crIter->cd(pad)->SetLeftMargin(0.12);
            crIter->cd(pad)->SetRightMargin(0.01);
            crIter->cd(pad)->SetTopMargin(0.01);
            crIter->cd(pad)->SetBottomMargin(0.08);
            gPad->SetTickx(1);
            gPad->SetTicky(1);

            auto makeName = [&](int it)
            {
                std::ostringstream s;
                s << "unfolded_zT_" << jetLabel << "_" << rapLabel << "_iter" << it;
                return s.str();
            };

            TH1 *rPrev = (globalBaseIter > 1) ? dynamic_cast<TH1 *>(f->Get(makeName(globalBaseIter - 1).c_str())) : nullptr;
            TH1 *rBase = dynamic_cast<TH1 *>(f->Get(makeName(globalBaseIter).c_str()));
            TH1 *rNext = dynamic_cast<TH1 *>(f->Get(makeName(globalBaseIter + 1).c_str()));

            if (!rBase)
                continue;

            // make ratios
            rPrev = (rPrev && rPrev->GetNbinsX() == rBase->GetNbinsX()) ? static_cast<TH1 *>(rPrev->Clone((std::string("ratio_iter_prev_") + jetLabel + "_" + rapLabel).c_str())) : nullptr;
            rNext = (rNext && rNext->GetNbinsX() == rBase->GetNbinsX()) ? static_cast<TH1 *>(rNext->Clone((std::string("ratio_iter_next_") + jetLabel + "_" + rapLabel).c_str())) : nullptr;
            if (rPrev)
            {
                rPrev->SetDirectory(nullptr);
                rPrev->Divide(rPrev, rBase, 1, 1, "b");
                rPrev->SetMarkerStyle(21);
                rPrev->SetMarkerColor(kRed + 1);
                rPrev->SetLineColor(kRed + 1);
            }
            if (rNext)
            {
                rNext->SetDirectory(nullptr);
                rNext->Divide(rNext, rBase, 1, 1, "b");
                rNext->SetMarkerStyle(24);
                rNext->SetMarkerColor(kGreen + 2);
                rNext->SetLineColor(kGreen + 2);
            }

            TH1D *frameIter = nullptr;
            if (rPrev)
                frameIter = static_cast<TH1D *>(rPrev->Clone("frame_iter_ratio_global"));
            else if (rNext)
                frameIter = static_cast<TH1D *>(rNext->Clone("frame_iter_ratio_global"));
            if (!frameIter)
                continue;
            frameIter->SetDirectory(nullptr);
            frameIter->SetTitle(";z_{T};iteration / iter");
            frameIter->GetYaxis()->SetRangeUser(0.5, 1.5);
            frameIter->Draw("AXIS");

            // compute and draw per-bin systematic band (RMS of deviations from 1.0 of available +/-1 iters)
            // {
            const int nb = frameIter->GetNbinsX();
            TH1D *sysBandG = (TH1D *)(frameIter->Clone((std::string("sys_iter_band_global_") + jetLabel + "_" + rapLabel).c_str()));
            int nonzeroPrev = 0, nonzeroNext = 0, nonzeroSys = 0;
            double sigmaSum = 0.0, sigmaMax = 0.0;
            double sigmaMin = std::numeric_limits<double>::infinity();
            for (int b = 1; b <= nb; ++b)
            {
                if (rPrev)
                {
                    double pc = rPrev->GetBinContent(b);
                    if (std::isfinite(pc) && pc != 0.0)
                        ++nonzeroPrev;
                }
                if (rNext)
                {
                    double nc = rNext->GetBinContent(b);
                    if (std::isfinite(nc) && nc != 0.0)
                        ++nonzeroNext;
                }

                int count = 0;
                double sumsq = 0.0;
                if (rPrev)
                {
                    double dv = rPrev->GetBinContent(b) - 1.0;
                    if (std::isfinite(dv))
                    {
                        sumsq += dv * dv;
                        ++count;
                    }
                }
                if (rNext)
                {
                    double dv = rNext->GetBinContent(b) - 1.0;
                    if (std::isfinite(dv))
                    {
                        sumsq += dv * dv;
                        ++count;
                    }
                }
                double sigma = (count > 0) ? std::sqrt(sumsq / count) : 0.0;
                if (sysBandG)
                {
                    sysBandG->SetBinContent(b, 1.0);
                    sysBandG->SetBinError(b, sigma);
                }
                if (sigma > 0.0 && std::isfinite(sigma))
                {
                    ++nonzeroSys;
                    sigmaSum += sigma;
                    sigmaMax = std::max(sigmaMax, sigma);
                    sigmaMin = std::min(sigmaMin, sigma);
                }
            }

            if (sysBandG)
            {
                sysBandG->SetFillColor(kGray + 2);
                sysBandG->SetMarkerStyle(1);
                sysBandG->SetMarkerColor(kGray + 2);
                sysBandG->SetFillStyle(3844);
                sysBandG->SetLineColor(kGray + 2);
                sysBandG->Draw("E3 SAME");
            }
            // }

            if (rPrev)
            {
                rPrev->SetLineColor(kRed + 1);
                rPrev->SetMarkerColor(kRed + 1);
                rPrev->Draw("E1 SAME");
            }
            if (rNext)
            {
                rNext->SetLineColor(kGreen + 2);
                rNext->SetMarkerColor(kGreen + 2);
                rNext->Draw("E1 SAME");
            }
            TLine l(frameIter->GetXaxis()->GetXmin(), 1.0, frameIter->GetXaxis()->GetXmax(), 1.0);
            l.SetLineStyle(2);
            l.Draw();

            if (pad == 1)
            {
                TLegend *lgr = new TLegend(0.4, 0.75, 0.95, 0.9);
                lgr->SetFillStyle(0);
                lgr->SetBorderSize(0);
                if (rPrev)
                    lgr->AddEntry(rPrev, (std::string("iter") + std::to_string(globalBaseIter - 1) + " / iter" + std::to_string(globalBaseIter)).c_str(), "lep");
                if (rNext)
                    lgr->AddEntry(rNext, (std::string("iter") + std::to_string(globalBaseIter + 1) + " / iter" + std::to_string(globalBaseIter)).c_str(), "lep");
                if (sysBandG)
                    lgr->AddEntry(sysBandG, "RMS syst.", "f");
                lgr->Draw();
            }

            TLine *unitylineiter = new TLine(frameIter->GetXaxis()->GetXmin(), 1.0, frameIter->GetXaxis()->GetXmax(), 1.0);
            unitylineiter->SetLineStyle(7);
            unitylineiter->SetLineColor(kBlack);
            unitylineiter->Draw();

            TLatex t;
            t.SetNDC();
            t.SetTextSize(0.04);
            std::ostringstream ttl;
            ttl << "jet p_{T}: " << jetLabel << " GeV/c   rapidity: " << displayRap;
            t.DrawLatex(0.22, 0.92, ttl.str().c_str());
        }
    }
    crIter->Modified();
    crIter->Update();
    // annotate dataset / experiment / energy on pad 3
    drawDatasetLabelOnCanvas(crIter);
    const std::string outPathIterR = outDir + "/iter_plusminus_ratio_summary" + datasetTag + ".png";
    crIter->SaveAs(outPathIterR.c_str());
    const std::string outPathIterRpdf = outDir + "/iter_plusminus_ratio_summary" + datasetTag + ".pdf";
    crIter->SaveAs(outPathIterRpdf.c_str());

    // Overlay canvas: measured, response-measured, response-truth, unfolded
    TCanvas *co = new TCanvas("c_response_overlay", "Measured + Response + Unfolded", canvasW, canvasH);
    co->Divide(ncol, nrow);
    for (int iy = 0; iy < nrow; ++iy)
    {
        const std::string &rapLabel = rapBins[iy];
        const std::string displayRap = (rapLabels.size() > static_cast<size_t>(iy) ? rapLabels[iy] : rapLabel);
        for (int ix = 0; ix < ncol; ++ix)
        {
            const std::string &jetLabel = jetBins[ix];
            const int pad = iy * ncol + ix + 1;
            co->cd(pad)->SetLeftMargin(0.12);
            co->cd(pad)->SetRightMargin(0.01);
            co->cd(pad)->SetTopMargin(0.01);
            co->cd(pad)->SetBottomMargin(0.08);
            gPad->SetTickx(1);
            gPad->SetTicky(1);

            const std::string measuredName = "measured_zT_" + jetLabel + "_" + rapLabel;
            const std::string unfoldedNameTotal = "best_unfolded_with_totalErr_" + SanitizeLabel(jetLabel) + "_" + rapLabel;
            const std::string unfoldedName = "best_unfolded_zT_" + jetLabel + "_" + rapLabel;

            TH1 *hmeas = dynamic_cast<TH1 *>(f->Get(measuredName.c_str()));
            TH1 *hunf = dynamic_cast<TH1 *>(f->Get(unfoldedNameTotal.c_str()));
            if (!hunf)
                hunf = dynamic_cast<TH1 *>(f->Get(unfoldedName.c_str()));

            // try several candidate names for response-measured and response-truth
            TH1 *hrespM = nullptr;
            TH1 *hrespT = nullptr;
            std::vector<std::string> respM_cands = {
                "responseMeasured_zT_" + jetLabel + "_" + rapLabel,
                "response_measured_zT_" + jetLabel + "_" + rapLabel,
                "responseMeasured_" + jetLabel + "_" + rapLabel,
                "response_measured_" + jetLabel + "_" + rapLabel};
            std::vector<std::string> respT_cands = {
                "responseTruth_zT_" + jetLabel + "_" + rapLabel,
                "response_truth_zT_" + jetLabel + "_" + rapLabel,
                "responseTruth_" + jetLabel + "_" + rapLabel,
                "response_truth_" + jetLabel + "_" + rapLabel};
            for (auto &n : respM_cands)
            {
                if (!hrespM)
                    hrespM = dynamic_cast<TH1 *>(f->Get(n.c_str()));
            }
            for (auto &n : respT_cands)
            {
                if (!hrespT)
                    hrespT = dynamic_cast<TH1 *>(f->Get(n.c_str()));
            }

            if (!hmeas && !hunf && !hrespM && !hrespT)
                continue;

            // prepare clones/scales so we don't modify originals
            TH1 *hrespM_draw = nullptr;
            TH1 *hrespT_draw = nullptr;
            if (hrespM)
            {
                hrespM_draw = static_cast<TH1 *>(hrespM->Clone((std::string("respM_draw_") + jetLabel + "_" + rapLabel).c_str()));
                hrespM_draw->SetDirectory(nullptr);
                if (hmeas && hrespM_draw->Integral() > 0 && hmeas->Integral() > 0)
                    hrespM_draw->Scale(hmeas->Integral() / hrespM_draw->Integral());
            }
            if (hrespT)
            {
                hrespT_draw = static_cast<TH1 *>(hrespT->Clone((std::string("respT_draw_") + jetLabel + "_" + rapLabel).c_str()));
                hrespT_draw->SetDirectory(nullptr);
                if (hmeas && hrespT_draw->Integral() > 0 && hmeas->Integral() > 0)
                    hrespT_draw->Scale(hmeas->Integral() / hrespT_draw->Integral());
            }

            double yMax = 0.0;
            auto consider = [&](TH1 *h)
            { if (h) yMax = std::max(yMax, h->GetMaximum()); };
            consider(hmeas);
            consider(hunf);
            consider(hrespM_draw);
            consider(hrespT_draw);
            if (yMax <= 0)
                yMax = 1.0;
            TLine *meanLine = nullptr;
            double meanval = 0.0;
            if (hunf)
            {
                hunf->SetLineColor(kBlue + 1);
                hunf->SetMarkerColor(kBlue + 1);
                hunf->SetMarkerStyle(20);
                hunf->SetMaximum(yMax * 1.3);
                hunf->Draw("E1");

                //determine mean value of histogram and draw vertical line at that value as well as add that value to the legend
                meanval = hunf->GetMean();
                meanLine = new TLine(meanval, 0, meanval, hunf->GetMaximum() * .65);
                meanLine->SetLineColor(kBlue + 1);
                meanLine->SetLineWidth(2);
                meanLine->SetLineStyle(2);
                meanLine->Draw();
            }
            if (hmeas)
            {
                hmeas->SetMarkerStyle(24);
                hmeas->SetMarkerColor(kBlack);
                hmeas->SetLineColor(kBlack);
                hmeas->Draw(hunf ? "E1 SAME" : "E1");
            }
            if (hrespM_draw)
            {
                hrespM_draw->SetLineColor(kMagenta + 2);
                hrespM_draw->SetLineWidth(2);
                hrespM_draw->Draw("HIST SAME");
            }
            if (hrespT_draw)
            {
                hrespT_draw->SetLineColor(kCyan + 1);
                hrespT_draw->SetLineWidth(2);
                hrespT_draw->SetLineStyle(7);
                hrespT_draw->Draw("HIST SAME");
            }

            if (pad == 1)
            {
                TLegend *leg = new TLegend(0.5, 0.65, 0.95, 0.9);
                leg->SetFillStyle(0);
                leg->SetBorderSize(0);
                if (hmeas)
                    leg->AddEntry(hmeas, "Measured", "lep");
                if (hrespM_draw)
                    leg->AddEntry(hrespM_draw, "Response (measured)", "l");
                if (hrespT_draw)
                    leg->AddEntry(hrespT_draw, "Response (truth, prior)", "l");
                if (hunf){
                    leg->AddEntry(hunf, "Unfolded", "lep");
                    leg->AddEntry(meanLine, (std::string("mean: ") + std::to_string(meanval)).c_str(), "l");
                }
                leg->Draw();
            } else {
                //only add mean value to legend
                TLegend *leg = new TLegend(0.15, 0.8, 0.6, 0.9);
                leg->SetFillStyle(0);
                leg->SetBorderSize(0);
                if (hunf){
                    // leg->AddEntry(hunf, "Unfolded", "lep");
                    leg->AddEntry(meanLine, (std::string("mean: ") + std::to_string(meanval)).c_str(), "l");
                }
                leg->Draw();
             }

            TLatex t;
            t.SetNDC();
            t.SetTextSize(0.04);
            std::ostringstream ttl;
            ttl << "jet p_{T}: " << jetLabel << " GeV/c   rapidity: " << displayRap;
            t.DrawLatex(0.22, 0.92, ttl.str().c_str());
        }
    }
    co->Modified();
    co->Update();
    // annotate dataset / experiment / energy on pad 3
    drawDatasetLabelOnCanvas(co);
    const std::string outPathO = outDir + "/unfold_response_overlay" + datasetTag + ".png";
    co->SaveAs(outPathO.c_str());
    const std::string outPathOpdf = outDir + "/unfold_response_overlay" + datasetTag + ".pdf";
    co->SaveAs(outPathOpdf.c_str());

    // Multi-panel: measured response matrices only
    TCanvas *cRespMeasured = new TCanvas("c_response_matrices_measured", "Response matrices (measured)", canvasW, canvasH);
    cRespMeasured->Divide(ncol, nrow);
    for (int iy = 0; iy < nrow; ++iy)
    {
        const std::string &rapLabel = rapBins[iy];
        const std::string displayRap = (rapLabels.size() > static_cast<size_t>(iy) ? rapLabels[iy] : rapLabel);
        for (int ix = 0; ix < ncol; ++ix)
        {
            const std::string &jetLabel = jetBins[ix];
            const int pad = iy * ncol + ix + 1;
            cRespMeasured->cd(pad);
            gPad->SetLeftMargin(0.12);
            gPad->SetRightMargin(0.16);
            gPad->SetTopMargin(0.08);
            gPad->SetBottomMargin(0.12);

            // try actual response name used by the unfolder
            std::vector<std::string> respCandidates = { std::string("responseMatrix_") + jetLabel + "_" + rapLabel,
                                                        std::string("responseMatrix_") + SanitizeLabel(jetLabel) + "_" + rapLabel };
            TH2 *hResp = nullptr;
            for (const auto &n : respCandidates)
            {
                if (!hResp)
                    hResp = dynamic_cast<TH2 *>(f->Get(n.c_str()));
            }
            if (!hResp)
                continue;
            TH2 *hdraw = static_cast<TH2 *>(hResp->Clone((std::string("respMeasured_panel_") + jetLabel + "_" + rapLabel).c_str()));
            hdraw->SetDirectory(nullptr);
            hdraw->SetTitle("");
            hdraw->GetXaxis()->SetTitle("reco D^{0}_{z} (det)");
            hdraw->GetYaxis()->SetTitle("truth D^{0}_{z} (mc)");
            hdraw->Draw("COLZ");
            TLatex t; t.SetNDC(); t.SetTextSize(0.035);
            std::ostringstream ttl; ttl << "p_{T}: " << jetLabel << "  y: " << displayRap;
            t.DrawLatex(0.02, 0.96, ttl.str().c_str());
        }
    }
    cRespMeasured->Modified(); cRespMeasured->Update();
    const std::string outRespMeasuredPng = outDir + "/response_matrices_measured" + datasetTag + ".png";
    cRespMeasured->SaveAs(outRespMeasuredPng.c_str());
    const std::string outRespMeasuredPdf = outDir + "/response_matrices_measured" + datasetTag + ".pdf";
    cRespMeasured->SaveAs(outRespMeasuredPdf.c_str());

    // Multi-panel: truth view (transposed + per-column normalized)
    TCanvas *cRespTruth = new TCanvas("c_response_matrices_truth", "Response matrices (truth view)", canvasW, canvasH);
    cRespTruth->Divide(ncol, nrow);
    for (int iy = 0; iy < nrow; ++iy)
    {
        const std::string &rapLabel = rapBins[iy];
        const std::string displayRap = (rapLabels.size() > static_cast<size_t>(iy) ? rapLabels[iy] : rapLabel);
        for (int ix = 0; ix < ncol; ++ix)
        {
            const std::string &jetLabel = jetBins[ix];
            const int pad = iy * ncol + ix + 1;
            cRespTruth->cd(pad);
            gPad->SetLeftMargin(0.12);
            gPad->SetRightMargin(0.16);
            gPad->SetTopMargin(0.08);
            gPad->SetBottomMargin(0.12);

            std::vector<std::string> respCandidates = { std::string("responseMatrix_") + jetLabel + "_" + rapLabel,
                                                        std::string("responseMatrix_") + SanitizeLabel(jetLabel) + "_" + rapLabel };
            TH2D *hResp = nullptr;
            for (const auto &n : respCandidates)
            {
                if (!hResp)
                    hResp = dynamic_cast<TH2D *>(f->Get(n.c_str()));
            }
            if (!hResp)
                continue;

            // transpose: create new TH2D with swapped axes
            const int nx = hResp->GetNbinsX();
            const int ny = hResp->GetNbinsY();
            std::vector<double> xedges = ExtractAxisEdges(hResp->GetXaxis());
            std::vector<double> yedges = ExtractAxisEdges(hResp->GetYaxis());
            TH2D *hTrans = new TH2D((std::string("respTrans_") + jetLabel + "_" + rapLabel).c_str(), "", ny, &yedges[0], nx, &xedges[0]);
            for (int i = 1; i <= nx; ++i)
            {
                for (int j = 1; j <= ny; ++j)
                {
                    double val = hResp->GetBinContent(i, j);
                    hTrans->SetBinContent(j, i, val);
                }
            }

            // normalize each column (x) so sum_y = 1.0 if possible
            for (int ixbin = 1; ixbin <= hTrans->GetNbinsX(); ++ixbin)
            {
                double sum = 0.0;
                for (int iybin = 1; iybin <= hTrans->GetNbinsY(); ++iybin) sum += hTrans->GetBinContent(ixbin, iybin);
                if (sum > 0)
                {
                    for (int iybin = 1; iybin <= hTrans->GetNbinsY(); ++iybin)
                    {
                        double v = hTrans->GetBinContent(ixbin, iybin);
                        hTrans->SetBinContent(ixbin, iybin, v / sum);
                    }
                }
            }

            hTrans->SetDirectory(nullptr);
            hTrans->SetTitle("");
            hTrans->GetXaxis()->SetTitle("truth D^{0}_{z} (mc)");
            hTrans->GetYaxis()->SetTitle("reco D^{0}_{z} (det)");
            hTrans->Draw("COLZ");
            TLatex t; t.SetNDC(); t.SetTextSize(0.035);
            std::ostringstream ttl; ttl << "p_{T}: " << jetLabel << "  y: " << displayRap;
            t.DrawLatex(0.02, 0.96, ttl.str().c_str());
        }
    }
    cRespTruth->Modified(); cRespTruth->Update();
    const std::string outRespTruthPng = outDir + "/response_matrices_truth_view" + datasetTag + ".png";
    cRespTruth->SaveAs(outRespTruthPng.c_str());
    const std::string outRespTruthPdf = outDir + "/response_matrices_truth_view" + datasetTag + ".pdf";
    cRespTruth->SaveAs(outRespTruthPdf.c_str());

    // Closure-test summary: closureUnfolded / pseudoTruth
    TCanvas *cClosure = new TCanvas("c_closure_ratio_summary", "Closure ratio summary", canvasW, canvasH);
    cClosure->Divide(ncol, nrow);
    for (int iy = 0; iy < nrow; ++iy)
    {
        const std::string &rapLabel = rapBins[iy];
        const std::string displayRap = (rapLabels.size() > static_cast<size_t>(iy) ? rapLabels[iy] : rapLabel);
        for (int ix = 0; ix < ncol; ++ix)
        {
            const std::string &jetLabel = jetBins[ix];
            const int pad = iy * ncol + ix + 1;
            cClosure->cd(pad)->SetLeftMargin(0.12);
            cClosure->cd(pad)->SetRightMargin(0.01);
            cClosure->cd(pad)->SetTopMargin(0.01);
            cClosure->cd(pad)->SetBottomMargin(0.08);
            gPad->SetTickx(1);
            gPad->SetTicky(1);

            const std::string closureName = "closureUnfolded_" + jetLabel + "_" + rapLabel;
            const std::string pseudoTruthName = "closurePseudoTruth_" + jetLabel + "_" + rapLabel;

            TH1 *hClosure = dynamic_cast<TH1 *>(f->Get(closureName.c_str()));
            TH1 *hPseudoTruth = dynamic_cast<TH1 *>(f->Get(pseudoTruthName.c_str()));

            if (!hClosure || !hPseudoTruth)
            {
                std::cout << "[closure_summary] missing histogram(s) for jet=" << jetLabel
                          << " rap=" << rapLabel
                          << " closure=" << (hClosure ? "Y" : "N")
                          << " pseudoTruth=" << (hPseudoTruth ? "Y" : "N")
                          << std::endl;
                continue;
            }

            if (hClosure->GetNbinsX() != hPseudoTruth->GetNbinsX())
            {
                std::cout << "[closure_summary] NBINS MISMATCH for jet=" << jetLabel << " rap=" << rapLabel
                          << " closure_nbins=" << hClosure->GetNbinsX() << " pseudo_nbins=" << hPseudoTruth->GetNbinsX() << std::endl;
            }

            TH1D *hRatio(MakeRatioHistogram(hClosure, hPseudoTruth,
                                            std::string("ratio_closure_truth_") + jetLabel + "_" + rapLabel,
                                            "closure / truth"));
            if (!hRatio)
            {
                std::cout << "[closure_summary] failed to build ratio for jet=" << jetLabel
                          << " rap=" << rapLabel << std::endl;
                continue;
            }

            hRatio->SetTitle(";z_{T};closure / truth");
            hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
            hRatio->Draw("AXIS");

            hRatio->SetLineColor(kBlue + 1);
            hRatio->SetMarkerColor(kBlue + 1);
            hRatio->SetMarkerStyle(20);
            hRatio->Draw("E1");

            TLine *unityClosure = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
            unityClosure->SetLineStyle(7);
            unityClosure->SetLineColor(kBlack);
            unityClosure->Draw("same");

            // fit the ratio to a constant, plot it and add to legend
            TF1 *fConst = new TF1((std::string("fConstClosure_") + jetLabel + "_" + rapLabel).c_str(), "[0]", hRatio->GetXaxis()->GetXmin(), hRatio->GetXaxis()->GetXmax());
            // fConst->SetParameter(0.0, 0.9);
            hRatio->Fit(fConst, "Q0M"); // quiet, no drawing, use x-range of histogram
            fConst->SetLineColor(kGreen + 2);
            // fConst->SetLineStyle(7);
            fConst->SetLineWidth(2);
            fConst->Draw("SAME");

            if (pad == 1)
            {
                TLegend *legClosure = new TLegend(0.3, 0.78, 0.95, 0.90);
                legClosure->SetFillStyle(0);
                legClosure->SetBorderSize(0);
                legClosure->AddEntry(hRatio, "Closure unfolded / pseudo truth", "lep");
                // legClosure->AddEntry(fConst, Form("Const fit - %.2f", fConst->GetParameter(0)), "l");
                legClosure->Draw();
            }
            // only report contsant fit value on other pads to avoid overcrowding legend
            // std::cout << "[closure_summary] Const fit - " << fConst->GetParameter(0) << std::endl;
            // put legend in bottom left for const fit
            TLegend *legConst = new TLegend(0.15, 0.12, 0.65, 0.2);
            legConst->SetFillStyle(0);
            legConst->SetBorderSize(0);
            legConst->AddEntry(fConst, Form("Const fit - %.2f", fConst->GetParameter(0)), "l");
            legConst->Draw();

            TLatex t;
            t.SetNDC();
            t.SetTextSize(0.04);
            std::ostringstream ttl;
            ttl << "jet p_{T}: " << jetLabel << " GeV/c   rapidity: " << displayRap;
            t.DrawLatex(0.22, 0.92, ttl.str().c_str());
        }
    }
    cClosure->Modified();
    cClosure->Update();
    // annotate dataset / experiment / energy on pad 3
    drawDatasetLabelOnCanvas(cClosure);
    const std::string outPathClosure = outDir + "/closure_ratio_summary" + datasetTag + ".png";
    cClosure->SaveAs(outPathClosure.c_str());
    const std::string outPathClosurePdf = outDir + "/closure_ratio_summary" + datasetTag + ".pdf";
    cClosure->SaveAs(outPathClosurePdf.c_str());

    // Multi-panel: kinematic efficiency overlays (one panel per rapidity bin)
    int canvasH2 = std::min(350 * 1, 2400);
    int canvasW2 = std::min(350 * nrow, 1800);
    TCanvas *cKin = new TCanvas("c_kinematic_efficiency", "Kinematic efficiency by rapidity", canvasW2, canvasH2);
    // arrange as 1 column x nrow panels (one panel per rapidity bin)
    int kinCols = 1;
    int kinRows = nrow > 0 ? nrow : 1;
    cKin->Divide(kinRows, kinCols);
    //define 5 colors for the jet pT bins to use for the kinematic efficiency histograms
    const int jetColors[] = { kBlue+2, kRed+2, kGreen+2, kMagenta+2, kCyan+2 }; // purple

    for (int iy = 0; iy < nrow; ++iy)
    {
        const std::string &rapLabel = rapBins[iy];
        const std::string displayRap = (rapLabels.size() > static_cast<size_t>(iy) ? rapLabels[iy] : rapLabel);
        const int pad = iy + 1; // since 1 column
        cKin->cd(pad);
        gPad->SetLeftMargin(0.09);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.01);
        gPad->SetBottomMargin(0.08);

        // collect kinematic-eff histograms for all jet bins for this rapidity
        std::vector<TH1 *> kinHists;
        for (size_t ix = 0; ix < jetBins.size(); ++ix)
        {
            const std::string &jetLabel = jetBins[ix];
            std::vector<std::string> candidates = { std::string("kinEff_zT_") + jetLabel + "_" + rapLabel,
                                                    std::string("kinEff_zT_") + SanitizeLabel(jetLabel) + "_" + rapLabel };
            TH1 *h = nullptr;
            for (const auto &n : candidates)
            {
                if (!h)
                    h = dynamic_cast<TH1 *>(f->Get(n.c_str()));
            }
            if (h)
            {
                h->SetDirectory(nullptr);
                kinHists.push_back(h);
            }
        }

        if (kinHists.empty())
        {
            // nothing to draw for this rapidity
            TLatex t; t.SetNDC(); t.SetTextSize(0.04);
            t.DrawLatex(0.2, 0.5, (std::string("No kinematic-eff histograms for rapidity: ") + displayRap).c_str());
            continue;
        }

        // determine y-range
        double yMax = 0.0;
        double yMin = 1e9;
        for (auto h : kinHists)
        {
            yMax = std::max(yMax, h->GetMaximum());
            // consider non-zero minima
            for (int b = 1; b <= h->GetNbinsX(); ++b)
            {
                double v = h->GetBinContent(b);
                if (v > 0) yMin = std::min(yMin, v);
            }
        }
        if (yMax <= 0) yMax = 1.0;
        if (yMin <= 0 || yMin > yMax/1000.0) yMin = 0.0;

        // draw first histogram to create axes
        TH1 *h0 = kinHists.front();
        h0->SetTitle((std::string(";z_{T};kinematic efficiency (") + displayRap + ")").c_str());
        h0->GetYaxis()->SetRangeUser(0, 1.149);
        h0->SetLineColor(jetColors[0]);
        h0->SetMarkerColor(jetColors[0]);
        h0->SetMarkerStyle(20);
        h0->Draw("E1");

        // draw remaining histograms
        for (size_t i = 1; i < kinHists.size(); ++i)
        {
            TH1 *h = kinHists[i];
            int col = 2 + (i % 7); // choose varied colors
            h->SetLineColor(jetColors[i]);
            h->SetMarkerColor(jetColors[i]);
            h->SetMarkerStyle(20 + (i % 4));
            h->Draw("E1 SAME");
        }
        TLine* unityLine = new TLine(h0->GetXaxis()->GetXmin(), 1.0, h0->GetXaxis()->GetXmax(), 1.0);
        unityLine->SetLineStyle(7);
        unityLine->SetLineColor(kGray+2);
        unityLine->Draw();
        // legend: list jet pT labels
        TLegend *leg = new TLegend(0.55, 0.15, 0.95, 0.40);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        //add an entry without an object that just says "Jet pT bins:"
        leg->AddEntry((TObject *)nullptr, "Jet p_{T} bins:", "");
        for (size_t i = 0; i < kinHists.size(); ++i)
        {
            const std::string &jetLabel = jetBins[i];
            //replace the underscores in the jet label with spaces for better display in the legend
            std::string displayJetLabel = jetLabel;
            std::replace(displayJetLabel.begin(), displayJetLabel.end(), '_', '-');
            //add GeV/c unit to the jet label if it doesn't already have it
            if (displayJetLabel.find("GeV") == std::string::npos)
                displayJetLabel += " GeV/c";
            leg->AddEntry(kinHists[i], displayJetLabel.c_str(), "lep");
        }
        leg->Draw();

        TLatex t; t.SetNDC(); t.SetTextSize(0.04);
        t.DrawLatex(0.22, 0.92, (std::string("rapidity: ") + displayRap).c_str());
    }

    cKin->Modified(); cKin->Update();
    const std::string outKinPng = outDir + "/kinematic_efficiency_by_rapidity" + datasetTag + ".png";
    cKin->SaveAs(outKinPng.c_str());
    const std::string outKinPdf = outDir + "/kinematic_efficiency_by_rapidity" + datasetTag + ".pdf";
    cKin->SaveAs(outKinPdf.c_str());
    
    // Multi-panel: correlation coefficient plots per iteration (iter 3,4,5)
    const int plotIters[] = {3, 4, 5};
    for (int iterIdx = 0; iterIdx < 3; ++iterIdx)
    {
        const int iteration = plotIters[iterIdx];
        std::ostringstream cname;
        cname << "c_correlation_iter" << iteration;
        TCanvas *cCorr = new TCanvas(cname.str().c_str(), (std::string("Correlation coefficients iter") + std::to_string(iteration)).c_str(), canvasW, canvasH);
        cCorr->Divide(ncol, nrow);
        for (int iy = 0; iy < nrow; ++iy)
        {
            const std::string &rapLabel = rapBins[iy];
            const std::string displayRap = (rapLabels.size() > static_cast<size_t>(iy) ? rapLabels[iy] : rapLabel);
            for (int ix = 0; ix < ncol; ++ix)
            {
                const std::string &jetLabel = jetBins[ix];
                const int pad = iy * ncol + ix + 1;
                cCorr->cd(pad)->SetLeftMargin(0.12);
                cCorr->cd(pad)->SetRightMargin(0.01);
                cCorr->cd(pad)->SetTopMargin(0.01);
                cCorr->cd(pad)->SetBottomMargin(0.08);
                gPad->SetTickx(1);
                gPad->SetTicky(1);

                // try names with explicit iteration suffix
                std::vector<std::string> candidates = {
                    std::string("correlation_") + SanitizeLabel(jetLabel) + "_" + rapLabel + "_iter" + std::to_string(iteration),
                    std::string("correlation_") + jetLabel + "_" + rapLabel + "_iter" + std::to_string(iteration),
                    std::string("correlation_") + SanitizeLabel(jetLabel) + "_" + SanitizeLabel(rapLabel) + "_iter" + std::to_string(iteration)
                };

                TObject *obj = nullptr;
                for (const auto &n : candidates)
                {
                    if (!obj)
                        obj = f->Get(n.c_str());
                }

                // fallback: scan keys for any name starting with 'correlation_' and containing the iteration suffix
                if (!obj)
                {
                    TIter next(f->GetListOfKeys());
                    TKey *k;
                    const std::string iterSuf = std::string("_iter") + std::to_string(iteration);
                    while ((k = (TKey *)next()))
                    {
                        const char *kn = k->GetName();
                        if (!kn)
                            continue;
                        std::string sname(kn);
                        if (sname.rfind("correlation_", 0) == 0 && sname.find(iterSuf) != std::string::npos)
                        {
                            if (sname.find(jetLabel) != std::string::npos || sname.find(rapLabel) != std::string::npos || sname.find(SanitizeLabel(jetLabel)) != std::string::npos)
                            {
                                obj = f->Get(sname.c_str());
                                break;
                            }
                            if (!obj)
                                obj = f->Get(sname.c_str());
                        }
                    }
                }

                if (!obj)
                    continue;

                TH2 *h2 = dynamic_cast<TH2 *>(obj);
                TH1 *h1 = dynamic_cast<TH1 *>(obj);
                if (h2)
                {
                    TH2 *hdraw = static_cast<TH2 *>(h2->Clone((std::string("corr_draw_iter") + std::to_string(iteration) + "_" + jetLabel + "_" + rapLabel).c_str()));
                    hdraw->SetDirectory(nullptr);
                    hdraw->SetTitle("");
                    hdraw->GetXaxis()->SetTitle("variable 1");
                    hdraw->GetYaxis()->SetTitle("variable 2");
                    hdraw->Draw("COLZ");
                }
                else if (h1)
                {
                    TH1 *hdraw = static_cast<TH1 *>(h1->Clone((std::string("corr_draw_iter") + std::to_string(iteration) + "_" + jetLabel + "_" + rapLabel).c_str()));
                    hdraw->SetDirectory(nullptr);
                    hdraw->SetTitle((std::string(";z_{T};correlation (iter=") + std::to_string(iteration) + ")").c_str());
                    hdraw->SetLineColor(kBlue+1);
                    hdraw->SetMarkerColor(kBlue+1);
                    hdraw->SetMarkerStyle(20);
                    hdraw->Draw("E1");
                    TLine unity(hdraw->GetXaxis()->GetXmin(), 0.0, hdraw->GetXaxis()->GetXmax(), 0.0);
                    unity.SetLineStyle(7);
                    unity.SetLineColor(kGray+2);
                    unity.Draw();
                }

                TLatex t; t.SetNDC(); t.SetTextSize(0.04);
                std::ostringstream ttl; ttl << "jet p_{T}: " << jetLabel << " GeV/c   rapidity: " << displayRap;
                t.DrawLatex(0.22, 0.92, ttl.str().c_str());
            }
        }
        cCorr->Modified(); cCorr->Update();
        const std::string outCorrPng = outDir + (std::string("/correlation_summary_iter") + std::to_string(iteration) + datasetTag + ".png");
        cCorr->SaveAs(outCorrPng.c_str());
        const std::string outCorrPdf = outDir + (std::string("/correlation_summary_iter") + std::to_string(iteration) + datasetTag + ".pdf");
        cCorr->SaveAs(outCorrPdf.c_str());
    }
}
