#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>

void shiftGraphX(TGraphErrors* graph, double shiftValue) {
    if (!graph) {
        std::cerr << "Error: Null TGraphErrors provided." << std::endl;
        return;
    }

    int nPoints = graph->GetN(); // Get the number of points in the graph
    for (int i = 0; i < nPoints; ++i) {
        double x, y;
        graph->GetPoint(i, x, y); // Get the current x and y values
        graph->SetPoint(i, x + shiftValue, y); // Shift the x value by the specified amount
    }
}
void removeNFirstPoints(TGraphErrors* graph, int nPointsToRemove) {
    if (!graph) {
        std::cerr << "Error: Null TGraphErrors provided." << std::endl;
        return;
    }

    int nPoints = graph->GetN(); // Get the number of points in the graph
    if (nPoints <= nPointsToRemove) {
        std::cerr << "Error: Graph has fewer than or equal to the number of points to remove." << std::endl;
        return;
    }

    // Remove the first nPointsToRemove points
    for (int i = 0; i < nPointsToRemove; ++i) {
        graph->RemovePoint(0);
    }
}
// void compareSystems(const char* filepPb = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_pPb/FitParametersUnBinnedD0Y_5_10.root", const char* filePbp = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_Pbp/FitParametersUnBinnedD0Y_5_10.root") {
void compareSystems(const char* filepPb = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_pPb/FitParametersUnBinnedD0Y", const char* filePbp = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_Pbp2/FitParametersUnBinnedD0Y") {

    TString jetMomentumBins[] = {"5_10", "10_15", "15_20", "20_30", "30_50"};
    int nBins = sizeof(jetMomentumBins) / sizeof(jetMomentumBins[0]);

    // Arrays to store all loaded graphs
    std::vector<TGraphErrors*> graphs_pPb;
    std::vector<TGraphErrors*> graphs_Pbp;
    std::vector<TFile*> files_pPb;
    std::vector<TFile*> files_Pbp;

    // Loop through jet momentum bins to load files and graphs
    for (int i = 0; i < nBins; ++i) {
        // Construct file paths
        TString filepPb_full = TString(filepPb) + "_" + jetMomentumBins[i] + ".root";
        TString filePbp_full = TString(filePbp) + "_" + jetMomentumBins[i] + ".root";

        // Open pPb file
        TFile* inputfilepPb = TFile::Open(filepPb_full, "READ");
        if (!inputfilepPb || inputfilepPb->IsZombie()) {
            std::cerr << "Error: Could not open file " << filepPb_full << std::endl;
            // Clean up previously opened files
            for (auto file : files_pPb) if (file) file->Close();
            for (auto file : files_Pbp) if (file) file->Close();
            return;
        }
        files_pPb.push_back(inputfilepPb);

        // Open Pbp file
        TFile* inputfilePbp = TFile::Open(filePbp_full, "READ");
        if (!inputfilePbp || inputfilePbp->IsZombie()) {
            std::cerr << "Error: Could not open file " << filePbp_full << std::endl;
            // Clean up previously opened files
            for (auto file : files_pPb) if (file) file->Close();
            for (auto file : files_Pbp) if (file) file->Close();
            return;
        }
        files_Pbp.push_back(inputfilePbp);

        // Load graph from pPb file
        TGraphErrors* graph_pPb = dynamic_cast<TGraphErrors*>(inputfilepPb->Get("BinTagZMeanF"));
        if (!graph_pPb) {
            std::cerr << "Error: Could not find TGraphErrors 'BinTagZMeanF' in file " << filepPb_full << std::endl;
            // Clean up
            for (auto file : files_pPb) if (file) file->Close();
            for (auto file : files_Pbp) if (file) file->Close();
            return;
        }
        removeNFirstPoints(graph_pPb, 2); // Remove the first 2 points from the graph
        graphs_pPb.push_back(graph_pPb);

        // Load graph from Pbp file
        TGraphErrors* graph_Pbp = dynamic_cast<TGraphErrors*>(inputfilePbp->Get("BinTagZMeanF"));
        if (!graph_Pbp) {
            std::cerr << "Error: Could not find TGraphErrors 'BinTagZMeanF' in file " << filePbp_full << std::endl;
            // Clean up
            for (auto file : files_pPb) if (file) file->Close();
            for (auto file : files_Pbp) if (file) file->Close();
            return;
        }
        removeNFirstPoints(graph_Pbp, 2); // Remove the first 2 points from the graph
        graphs_Pbp.push_back(graph_Pbp);
    }

    // Create a canvas to draw the graphs
    TCanvas* canvas = new TCanvas("canvas", "", 600, 800);
    canvas->SetTickx();
    canvas->SetTicky();
    // canvas->SetLeftMargin(leftMargin);
    canvas->SetRightMargin(0.01);
    canvas->SetTopMargin(0.01);
    // canvas->SetBottomMargin(bottomMargin);

    // TH2D* dummyCanvasHisto = new TH2D("dummyCanvasHisto", "", 100, 1.810, 4.19, 100, 0.35, 0.61);
    TH2D* dummyCanvasHisto = new TH2D("dummyCanvasHisto", "", 100, 1.810, 4.59, 100, 0.0, 0.65);
    dummyCanvasHisto->GetXaxis()->SetTitle("Rapidity (y_{cms})");
    dummyCanvasHisto->GetYaxis()->SetTitle("Mean #it{z}_{T}");
    dummyCanvasHisto->SetStats(0);  // Disable stats box
    dummyCanvasHisto->Draw();

    // Define color gradients for blue (pPb) and red (Pbp)
    Int_t blueColors[] = {kBlue+2, kBlue+1, kBlue, kBlue-3, kBlue-7};  // Dark to light blue
    Int_t redColors[] = {kRed+2, kRed+1, kRed, kRed-3, kRed-7};       // Dark to light red
    Int_t markerStyles[] = {20, 21, 22, 23, 34};  // Different marker styles for each bin
    Int_t markerStylesOpen[] = {24, 25, 26, 27, 28}; // Open markers for Pbp

    // Plot all pPb graphs
    for (int i = 0; i < nBins; ++i) {
        TGraphErrors* graph_pPb = graphs_pPb[i];
        graph_pPb->SetMarkerStyle(markerStyles[i]);
        graph_pPb->SetMarkerColor(blueColors[i]);
        graph_pPb->SetLineColor(blueColors[i]);
        graph_pPb->SetMarkerSize(1.0);
        graph_pPb->Draw("p,same");
    }

    // Plot all Pbp graphs (shifted)
    for (int i = 0; i < nBins; ++i) {
        TGraphErrors* graph_Pbp = graphs_Pbp[i];
        shiftGraphX(graph_Pbp, 0.465); // Shift the Pbp graphs to the right
        graph_Pbp->SetMarkerStyle(markerStylesOpen[i]);
        graph_Pbp->SetMarkerColor(redColors[i]);
        graph_Pbp->SetLineColor(redColors[i]);
        graph_Pbp->SetMarkerSize(1.0);
        graph_Pbp->Draw("p,same");
    }

    // Add a comprehensive legend
    TLegend* legend = new TLegend(0.15, 0.13, 0.55, 0.33);
    legend->SetTextFont(42);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetTextSize(0.03);
    legend->SetNColumns(2);

        legend->AddEntry((TObject*)nullptr, "p-Pb p_{T}^{jet}:", "");
        legend->AddEntry((TObject*)nullptr, "Pb-p p_{T}^{jet}", "");
    
    // Add legend entries for each jet momentum bin
    for (int i = 0; i < nBins; ++i) {
        TString binLabel = jetMomentumBins[i].ReplaceAll("_", "-") + " GeV";
        // TString binLabel = "p_{T}^{jet} " + jetMomentumBins[i].ReplaceAll("_", "-") + " GeV";
        legend->AddEntry(graphs_pPb[i], binLabel, "lp");
        legend->AddEntry(graphs_Pbp[i], binLabel, "lp");
    }
    legend->Draw();

    //add label for LHcb in progress
    TLatex* label = new TLatex();
    label->SetTextSize(0.04);
    label->DrawLatexNDC(0.15, 0.93, "LHCb #bf{in progress}");
    // Add label for center of mass energy
    TLatex* label2 = new TLatex();
    label2->SetTextSize(0.04);
    label2->DrawLatexNDC(0.15, 0.90, "#bf{#sqrt{#it{s}_{NN}} = 8.16 TeV}");

    // Save the canvas as an image
    canvas->SaveAs("BinTagZMeanF_Comparison.pdf");

    // Clean up
    for (auto file : files_pPb) if (file) file->Close();
    for (auto file : files_Pbp) if (file) file->Close();
    delete canvas;
}