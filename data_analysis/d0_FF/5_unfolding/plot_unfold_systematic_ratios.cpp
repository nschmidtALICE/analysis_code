#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TF1.h>

#include <algorithm>
#include <cctype>
#include <ctime>
#include <iostream>
#include <memory>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

std::string SystematicSanitizeLabel(std::string label)
{
    std::replace(label.begin(), label.end(), '.', 'p');
    std::replace(label.begin(), label.end(), '/', '_');
    std::replace(label.begin(), label.end(), '\\', '_');
    std::replace(label.begin(), label.end(), ' ', '_');
    std::replace(label.begin(), label.end(), ',', '_');
    return label;
}

static std::string BuildTodayDateTag()
{
    std::time_t now = std::time(nullptr);
    std::tm localTime {};
#if defined(_WIN32)
    localtime_s(&localTime, &now);
#else
    localtime_r(&now, &localTime);
#endif

    char buffer[11] = {0};
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d", &localTime);
    return buffer;
}

static std::string ResolveVariationDirectoryName(const std::string &variationName,
                                                 const std::vector<std::string> &labels)
{
    std::string resolvedName = variationName;
    if (resolvedName.empty())
    {
        if (labels.size() > 1 && !labels[1].empty())
        {
            resolvedName = labels[1];
        }
        else
        {
            resolvedName = "unnamed_variation";
        }
    }

    resolvedName = SystematicSanitizeLabel(resolvedName);
    if (resolvedName.empty())
    {
        resolvedName = "unnamed_variation";
    }
    return resolvedName;
}

static std::vector<std::string> SystematicSplitCsv(const std::string &csv)
{
    std::vector<std::string> out;
    std::stringstream ss(csv);
    std::string item;
    while (std::getline(ss, item, ','))
    {
        const size_t first = item.find_first_not_of(" \t\n\r");
        if (first == std::string::npos)
        {
            continue;
        }
        const size_t last = item.find_last_not_of(" \t\n\r");
        out.push_back(item.substr(first, last - first + 1));
    }
    return out;
}

static std::string SystematicBaseName(const std::string &path)
{
    const size_t slashPos = path.find_last_of("/\\");
    std::string base = (slashPos == std::string::npos) ? path : path.substr(slashPos + 1);
    const size_t dotPos = base.find_last_of('.');
    if (dotPos != std::string::npos)
    {
        base = base.substr(0, dotPos);
    }
    return base;
}

static std::vector<std::string> BuildDefaultSystematicLabels(size_t nFiles)
{
    std::vector<std::string> labels;
    if (nFiles == 0)
    {
        return labels;
    }

    labels.push_back("default");
    for (size_t index = 1; index < nFiles; ++index)
    {
        std::ostringstream label;
        label << "variation" << index;
        labels.push_back(label.str());
    }
    return labels;
}

TH1D *MakeSystematicRatioHistogram(const TH1 *numerator, const TH1 *denominator,
                                   const std::string &name, const std::string &yTitle)
{
    if (!numerator || !denominator)
    {
        return nullptr;
    }

    TH1D *ratio = static_cast<TH1D *>(numerator->Clone(name.c_str()));
    ratio->SetDirectory(nullptr);
    ratio->SetTitle((name + ";z_{T};" + yTitle).c_str());
    ratio->Sumw2();
    ratio->Divide(numerator, denominator, 1.0, 1.0, "B");
    return ratio;
}

static std::string BuildUnfoldedHistogramName(const std::string &histPrefix,
                                              const std::string &jetLabel,
                                              const std::string &rapLabel)
{
    std::string jetKey = jetLabel;
    if (histPrefix.find("with_totalErr") != std::string::npos)
    {
        jetKey = SystematicSanitizeLabel(jetLabel);
    }
    return histPrefix + jetKey + "_" + rapLabel;
}

// Usage examples:
// root -l -b -q 'plot_unfold_systematic_ratios.cpp+("nominal.root,var1.root")'
// root -l -b -q 'plot_unfold_systematic_ratios.cpp+("nominal.root,var1.root,var2.root","default,tracking,prior")'
// root -l -b -q 'plot_unfold_systematic_ratios.cpp+("nominal.root,var1.root","default,tracking","10_15,15_20,20_30,30_100","eta0,eta1,eta2","2.5 < y < 3.0, 3.0 < y < 3.5, 3.5 < y < 4.0","best_unfolded_zT_",0.5,1.5,true,"tracking")'
void plot_unfold_systematic_ratios(
    const std::string &inputFilesCsv = "",
    const std::string &labelsCsv = "",
    const std::string &variationName = "",
    const std::string &jetBinsCsv = "10_15,15_20,20_30",
    // const std::string &jetBinsCsv = "10_15,15_20,20_30,30_100",
    const std::string &rapidityCsv = "eta0,eta1,eta2",
    const std::string &rapidityLabel = "2.5 < y < 3.0, 3.0 < y < 3.5, 3.5 < y < 4.0",
    const std::string &histPrefix = "best_unfolded_zT_",
    double ratioMin = 0.5,
    double ratioMax = 1.5,
    bool useTotalErr = true)
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    if (inputFilesCsv.empty())
    {
        std::cerr << "No input files provided. Pass a CSV list with the default file first, followed by up to two systematic variations." << std::endl;
        return;
    }

    const std::vector<std::string> inputFiles = SystematicSplitCsv(inputFilesCsv);
    if (inputFiles.size() < 2)
    {
        std::cerr << "Need at least two input files: one default and one systematic variation." << std::endl;
        return;
    }
    if (inputFiles.size() > 3)
    {
        std::cerr << "This macro supports at most three input files total (one default plus up to two variations)." << std::endl;
        return;
    }

    std::vector<std::string> labels = SystematicSplitCsv(labelsCsv);
    if (labels.empty())
    {
        labels = BuildDefaultSystematicLabels(inputFiles.size());
    }
    if (labels.size() < inputFiles.size())
    {
        const size_t originalSize = labels.size();
        labels.resize(inputFiles.size());
        for (size_t index = originalSize; index < inputFiles.size(); ++index)
        {
            std::ostringstream label;
            label << "variation" << index;
            labels[index] = label.str();
        }
    }
    if (labels.size() > inputFiles.size())
    {
        labels.resize(inputFiles.size());
    }

    for (size_t index = 0; index < inputFiles.size(); ++index)
    {
        std::cout << "[systematic_ratio] input[" << index << "] label='" << labels[index]
                  << "' file='" << inputFiles[index] << "'" << std::endl;
    }

    const std::string resolvedHistPrefix = useTotalErr ? "best_unfolded_with_totalErr_" : histPrefix;
    const std::string outputModeTag = (resolvedHistPrefix.find("with_totalErr") != std::string::npos) ? "_with_totalErr" : "";
    std::cout << "[systematic_ratio] using histogram prefix '" << resolvedHistPrefix << "'" << std::endl;

    const std::vector<std::string> jetBins = SystematicSplitCsv(jetBinsCsv);
    const std::vector<std::string> rapBins = SystematicSplitCsv(rapidityCsv);
    const std::vector<std::string> rapLabels = SystematicSplitCsv(rapidityLabel);

    const int ncol = static_cast<int>(jetBins.size());
    const int nrow = static_cast<int>(rapBins.size());
    if (ncol <= 0 || nrow <= 0)
    {
        std::cerr << "Need at least one jet bin and one rapidity bin." << std::endl;
        return;
    }

    std::vector<std::unique_ptr<TFile>> files;
    files.reserve(inputFiles.size());
    for (const std::string &fileName : inputFiles)
    {
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
        if (!file || file->IsZombie())
        {
            std::cerr << "Failed to open input file: " << fileName << std::endl;
            return;
        }
        files.push_back(std::move(file));
    }

    const int canvasW = std::min(300 * ncol, 2400);
    const int canvasH = std::min(300 * nrow, 1800);

    std::string outDir = ".";
    const size_t slashPos = inputFiles.front().find_last_of("/\\");
    if (slashPos != std::string::npos)
    {
        outDir = inputFiles.front().substr(0, slashPos);
    }

    const std::string dateTag = BuildTodayDateTag();
    const std::string variationDirName = ResolveVariationDirectoryName(variationName, labels);
    const std::string systematicOutDir = "systematic_studies/" + dateTag + "/" + variationDirName;
    // const std::string systematicOutDir = outDir + "/systematic_studies/" + dateTag + "/" + variationDirName;
    if (gSystem->mkdir(systematicOutDir.c_str(), true) != 0 && gSystem->AccessPathName(systematicOutDir.c_str()))
    {
        std::cerr << "Failed to create output directory: " << systematicOutDir << std::endl;
        return;
    }
    std::cout << "[systematic_ratio] writing outputs to '" << systematicOutDir << "'" << std::endl;

    std::string datasetTag;
    std::string datasetLabel;
    if (inputFiles.front().find("pPb") != std::string::npos)
    {
        datasetTag = "_pPb";
        datasetLabel = "pPb";
    }
    else if (inputFiles.front().find("Pbp") != std::string::npos)
    {
        datasetTag = "_Pbp";
        datasetLabel = "Pbp";
    }

    const std::string outRootPath = systematicOutDir + "/unfold_systematic_ratios" + outputModeTag + datasetTag + ".root";
    std::unique_ptr<TFile> outFile(TFile::Open(outRootPath.c_str(), "RECREATE"));
    if (!outFile || outFile->IsZombie())
    {
        std::cerr << "Warning: failed to create ratio output file: " << outRootPath << std::endl;
        outFile.reset();
    }

    auto drawDatasetLabelOnCanvas = [&](TCanvas *canvas) {
        if (!canvas)
        {
            return;
        }
        const int npads = ncol * nrow;
        if (npads < 3)
        {
            return;
        }

        canvas->cd(3);
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.01);
        gPad->SetBottomMargin(0.08);

        TLatex label;
        label.SetNDC();
        label.SetTextFont(42);
        label.SetTextAlign(31);
        label.SetTextSize(0.04);

        const double x = 1.0 - gPad->GetRightMargin() - 0.06;
        const double y = 0.875;
        label.DrawLatex(x, y, "LHC#it{b}");

        std::string energy = "#it{s}_{NN} = 8.16 TeV";
        if (!datasetLabel.empty())
        {
            energy = datasetLabel + std::string(", ") + energy;
        }
        label.DrawLatex(x, y - 0.045, energy.c_str());
    };

    const std::vector<int> colors = {kRed + 2, kBlue + 2};
    const std::vector<int> markers = {20, 21};

    TCanvas *canvas = new TCanvas("c_unfold_systematic_ratios", "Unfold systematic ratios", canvasW, canvasH);
    canvas->Divide(ncol, nrow);

    for (int iy = 0; iy < nrow; ++iy)
    {
        const std::string &rapKey = rapBins[iy];
        const std::string displayRap = (rapLabels.size() > static_cast<size_t>(iy)) ? rapLabels[iy] : rapKey;

        for (int ix = 0; ix < ncol; ++ix)
        {
            const std::string &jetKey = jetBins[ix];
            const int pad = iy * ncol + ix + 1;
            canvas->cd(pad)->SetLeftMargin(0.12);
            canvas->cd(pad)->SetRightMargin(0.01);
            canvas->cd(pad)->SetTopMargin(0.01);
            canvas->cd(pad)->SetBottomMargin(0.08);
            gPad->SetTickx(1);
            gPad->SetTicky(1);

            const std::string defaultHistName = BuildUnfoldedHistogramName(resolvedHistPrefix, jetKey, rapKey);
            TH1 *hDefault = dynamic_cast<TH1 *>(files.front()->Get(defaultHistName.c_str()));
            if (!hDefault)
            {
                std::cout << "[systematic_ratio] missing default histogram '" << defaultHistName
                          << "' in file '" << inputFiles.front() << "'" << std::endl;
                continue;
            }

            std::vector<TH1D *> ratioHists;
            std::vector<size_t> ratioSources;
            ratioHists.reserve(inputFiles.size() - 1);
            ratioSources.reserve(inputFiles.size() - 1);

            for (size_t index = 1; index < inputFiles.size(); ++index)
            {
                const std::string variationHistName = BuildUnfoldedHistogramName(resolvedHistPrefix, jetKey, rapKey);
                TH1 *hVariation = dynamic_cast<TH1 *>(files[index]->Get(variationHistName.c_str()));
                if (!hVariation)
                {
                    std::cout << "[systematic_ratio] missing variation histogram '" << variationHistName
                              << "' in file '" << inputFiles[index] << "'" << std::endl;
                    continue;
                }
                if (hVariation->GetNbinsX() != hDefault->GetNbinsX())
                {
                    std::cout << "[systematic_ratio] nbins mismatch for jet=" << jetKey
                              << " rap=" << rapKey
                              << " default=" << hDefault->GetNbinsX()
                              << " variation=" << hVariation->GetNbinsX() << std::endl;
                    continue;
                }

                const std::string ratioName = "ratio_" + SystematicSanitizeLabel(labels[index]) + "_"
                                              + SystematicSanitizeLabel(jetKey) + "_" + rapKey;
                TH1D *ratio = MakeSystematicRatioHistogram(
                    hVariation,
                    hDefault,
                    ratioName,
                    labels[index] + " / " + labels[0]);
                if (!ratio)
                {
                    continue;
                }

                ratio->SetLineColor(colors[index - 1]);
                ratio->SetMarkerColor(colors[index - 1]);
                ratio->SetMarkerStyle(markers[index - 1]);
                ratio->SetLineWidth(2);
                ratioHists.push_back(ratio);
                ratioSources.push_back(index);

                if (outFile)
                {
                    outFile->cd();
                    ratio->Write();
                }

                std::cout << "[systematic_ratio] built ratio for " << labels[index]
                          << " / " << labels[0]
                          << " jet=" << jetKey << " rap=" << rapKey
                          << " default_integral=" << hDefault->Integral()
                          << " variation_integral=" << hVariation->Integral() << std::endl;
            }

            if (ratioHists.empty())
            {
                continue;
            }

            TH1D *frame = static_cast<TH1D *>(ratioHists.front()->Clone((std::string("frame_systematic_ratio_") + jetKey + "_" + rapKey).c_str()));
            frame->SetDirectory(nullptr);
            frame->Reset();
            frame->SetTitle(";z_{T};variation / default");
            frame->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
            frame->Draw("AXIS");

            for (TH1D *ratio : ratioHists)
            {
                ratio->Draw("E1 SAME");
            }

            TLine* unity = new TLine(frame->GetXaxis()->GetXmin(), 1.0, frame->GetXaxis()->GetXmax(), 1.0);
            unity->SetLineStyle(7);
            unity->SetLineColor(kBlack);
            unity->Draw();

            // Per-panel constant fits and legend (bottom-left)
            TLegend *padLegend = new TLegend(0.15, 0.10, 0.85, 0.28);
            padLegend->SetFillStyle(0);
            padLegend->SetBorderSize(0);
            padLegend->SetTextSize(0.035);

            for (size_t ri = 0; ri < ratioHists.size(); ++ri)
            {
                TH1D *rhist = ratioHists[ri];
                const size_t srcIndex = ratioSources[ri];
                const std::string safeLabel = SystematicSanitizeLabel(labels[srcIndex]);
                const std::string fitName = std::string("fconst_") + safeLabel + "_pad" + std::to_string(pad);
                TF1 *fconst = new TF1(fitName.c_str(), "pol0", frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax());
                rhist->Fit(fconst, "Q0");
                const double p0 = fconst->GetParameter(0);
                const double p0e = fconst->GetParError(0);
                fconst->SetLineColor(colors[ri]);
                fconst->SetLineWidth(2);
                fconst->Draw("SAME");
                std::ostringstream entry;
                entry << labels[srcIndex] << " / " << labels[0] << ": ";
                entry << std::fixed << std::setprecision(3) << p0 << " #pm " << p0e;

                padLegend->AddEntry(fconst, entry.str().c_str(), "l");
            }
            padLegend->Draw();

            if (pad == 1)
            {
                TLegend *legend = new TLegend(0.35, 0.72, 0.95, 0.90);
                legend->SetFillStyle(0);
                legend->SetBorderSize(0);
                for (size_t index = 0; index < ratioHists.size(); ++index)
                {
                    legend->AddEntry(ratioHists[index], (labels[ratioSources[index]] + " / " + labels[0]).c_str(), "lep");
                }
                legend->Draw();
            }

            TLatex panelLabel;
            panelLabel.SetNDC();
            panelLabel.SetTextFont(42);
            panelLabel.SetTextSize(0.035);
            panelLabel.SetTextAlign(13);
            std::ostringstream jetText;
            jetText << "jet p_{T}: " << jetKey << " GeV/c";
            panelLabel.DrawLatex(0.15, 0.92, jetText.str().c_str());
            std::ostringstream rapText;
            rapText << "rapidity: " << displayRap;
            panelLabel.DrawLatex(0.15, 0.865, rapText.str().c_str());
        }
    }

    canvas->Modified();
    canvas->Update();
    drawDatasetLabelOnCanvas(canvas);

    const std::string outPng = systematicOutDir + "/unfold_systematic_ratio_summary" + outputModeTag + datasetTag + ".png";
    const std::string outPdf = systematicOutDir + "/unfold_systematic_ratio_summary" + outputModeTag + datasetTag + ".pdf";
    canvas->SaveAs(outPng.c_str());
    canvas->SaveAs(outPdf.c_str());

    // --- NEW: Shape comparison canvas: overlay original unfolded distributions (normalized) ---
    TCanvas *c_shape = new TCanvas("c_unfold_shape_compare", "Unfolded shape comparison", canvasW, canvasH);
    c_shape->Divide(ncol, nrow);
    for (int iy = 0; iy < nrow; ++iy)
    {
        const std::string &rapKey = rapBins[iy];
        for (int ix = 0; ix < ncol; ++ix)
        {
            const std::string &jetKey = jetBins[ix];
            const int pad = iy * ncol + ix + 1;
            c_shape->cd(pad)->SetLeftMargin(0.12);
            c_shape->cd(pad)->SetRightMargin(0.01);
            c_shape->cd(pad)->SetTopMargin(0.01);
            c_shape->cd(pad)->SetBottomMargin(0.08);
            gPad->SetTickx(1);
            gPad->SetTicky(1);

            const std::string histName = BuildUnfoldedHistogramName(resolvedHistPrefix, jetKey, rapKey);
            TH1 *hDef = dynamic_cast<TH1 *>(files.front()->Get(histName.c_str()));
            if (!hDef) continue;

            // Prepare frame using a clone of the default histogram (reset to show axes)
            TH1 *frame = static_cast<TH1 *>(hDef->Clone((std::string("frame_shape_") + jetKey + "_" + rapKey).c_str()));
            frame->SetDirectory(nullptr);
            frame->Reset();
            frame->SetTitle((std::string(";z_{T};dN/dz_{T} (absolute)")).c_str());

            // Draw default (black solid)
            TH1 *hdef_clone = static_cast<TH1 *>(hDef->Clone((std::string("hdef_norm_") + jetKey + "_" + rapKey).c_str()));
            hdef_clone->SetDirectory(nullptr);
            // keep absolute scale (do not normalize)
            hdef_clone->SetLineColor(kBlack);
            hdef_clone->SetLineWidth(2);
            hdef_clone->SetMarkerStyle(24);
            hdef_clone->SetMarkerSize(0.9);
            hdef_clone->SetMarkerColor(kBlack);

            // Determine y-axis range from max of default and variations
            double ymax = hdef_clone->GetMaximum();
            for (size_t idx = 1; idx < files.size(); ++idx)
            {
                TH1 *hv_tmp = dynamic_cast<TH1 *>(files[idx]->Get(histName.c_str()));
                if (!hv_tmp) continue;
                double m = hv_tmp->GetMaximum();
                if (m > ymax) ymax = m;
            }
            if (ymax <= 0) ymax = 1.0; // fallback
            frame->GetYaxis()->SetRangeUser(0.0, ymax * 1.15);
            frame->Draw("AXIS");
            hdef_clone->Draw("E1 SAME");

            // Draw variations (normalized) and build legend
            TLegend *lg = new TLegend(0.55, 0.65, 0.95, 0.90);
            lg->SetFillStyle(0);
            lg->SetBorderSize(0);
            lg->AddEntry(hdef_clone, labels[0].c_str(), "lep");
            for (size_t idx = 1; idx < files.size(); ++idx)
            {
                TH1 *hv = dynamic_cast<TH1 *>(files[idx]->Get(histName.c_str()));
                if (!hv) continue;
                TH1 *hvc = static_cast<TH1 *>(hv->Clone((std::string("hvar_norm_") + SystematicSanitizeLabel(labels[idx]) + "_" + jetKey + "_" + rapKey).c_str()));
                hvc->SetDirectory(nullptr);
                // keep absolute scale (do not normalize)
                int col = colors[(idx - 1) % colors.size()];
                hvc->SetLineColor(col);
                hvc->SetMarkerStyle(24);
                hvc->SetMarkerSize(0.9);
                hvc->SetMarkerColor(col);
                hvc->SetLineWidth(2);
                hvc->Draw("E1 SAME");
                lg->AddEntry(hvc, labels[idx].c_str(), "lep");
                if (outFile)
                {
                    outFile->cd();
                    hvc->Write();
                }
            }
            lg->Draw();

            // write default clone as well
            if (outFile)
            {
                outFile->cd();
                hdef_clone->Write();
            }

            TLatex panelLabel;
            panelLabel.SetNDC();
            panelLabel.SetTextFont(42);
            panelLabel.SetTextSize(0.035);
            panelLabel.SetTextAlign(13);
            std::ostringstream jetText;
            jetText << "jet p_{T}: " << jetKey << " GeV/c";
            panelLabel.DrawLatex(0.15, 0.92, jetText.str().c_str());
            std::ostringstream rapText;
            rapText << "rapidity: " << ((rapLabels.size() > static_cast<size_t>(iy)) ? rapLabels[iy] : rapKey);
            panelLabel.DrawLatex(0.15, 0.865, rapText.str().c_str());
        }
    }
    c_shape->Modified();
    c_shape->Update();
    drawDatasetLabelOnCanvas(c_shape);
    const std::string outShapePng = systematicOutDir + "/unfold_shape_comparison" + outputModeTag + datasetTag + ".png";
    const std::string outShapePdf = systematicOutDir + "/unfold_shape_comparison" + outputModeTag + datasetTag + ".pdf";
    c_shape->SaveAs(outShapePng.c_str());
    c_shape->SaveAs(outShapePdf.c_str());

    if (outFile)
    {
        outFile->Write();
        outFile->Close();
        std::cout << "[systematic_ratio] wrote ratio histograms to '" << outRootPath << "'" << std::endl;
    }
    std::cout << "[systematic_ratio] saved summary plots to '" << outPng << "' and '" << outPdf << "'" << std::endl;
}