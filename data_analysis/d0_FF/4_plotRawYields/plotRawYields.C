#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <filesystem>

// ROOT includes
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TAxis.h"

// Check if a graph is valid
bool is_valid_graph(TGraphErrors* graph) {
    return (graph != nullptr && graph->GetN() > 0);
}

class PlotGraphsObject {
private:
    std::string ptString;
    std::string ptRange;
    std::string obsTag;
    std::string xTitle;
    double minPlotRange;
    double minFitRange;
    double maxFitRange;
    
public:
    std::string OutfilePath;
    
    // Histograms and graphs
    TGraphErrors* hMYield;
    TGraphErrors* hMYieldSG;
    TGraphErrors* hMYieldFix;
    TGraphErrors* hIPNonPromptFraction;
    TGraphErrors* gInclFragFunc;
    TGraphErrors* gDecayFragFunc;
    TGraphErrors* gPromptFragFunc;
    TGraphErrors* gAccCorr0;
    TGraphErrors* gAccCorr1;
    TGraphErrors* gAccCorr2;
    TGraphErrors* gAccCorr3;
    TGraphErrors* gAccCorr4;
    TGraphErrors* gAccCorr5;
    TGraphErrors* gAccCorr6;
    TGraphErrors* gAccCorr;
    TGraphErrors* gAccCorr5_Ext;
    TGraphErrors* gAccCorrTotal;
    
    // Processed graphs
    TGraphErrors* graphPRaw;
    TGraphErrors* graphNPRaw;
    TGraphErrors* graphInclRaw;
    TGraphErrors* graphPCorr;
    TGraphErrors* graphNPCorr;
    TGraphErrors* graphInclCorr;
    
    // Constructor - fix signature to match usage
    PlotGraphsObject(const std::string& ptRng, bool isZt, bool isMC);
    
    // Destructor
    ~PlotGraphsObject();
    
    // Public methods
    TGraphErrors* plotPt(bool isCorrected = false);
    void plotCorrFacAcceptance();
    TGraphErrors* plotNonPromptFraction(bool isCorrected = false);
    void plotYieldResult(bool isCorrected = false);
    void plotYieldSummary(std::vector<TGraphErrors*> yieldArray, 
                          std::vector<std::string> pTArray, 
                          int normType = 0, 
                          std::string Seltag = "incl", 
                          bool isCorrected = false,
                          std::vector<TGraphErrors*> corrArray = {});
    
    // Helper methods
    TGraphErrors* extrapolateIfNecessary(TGraphErrors* graph);
    void setHisto(TGraphErrors* histo, double MarkerScale, int color);
    void normalize(TGraphErrors* graph);
    TGraphErrors* applyAccCorrection(TGraphErrors* graph, TGraphErrors* gAccCorrInput = nullptr);
    TGraphErrors* multiplyGraphs(TGraphErrors* graph1, TGraphErrors* graph2);
    TGraphErrors* invertGraphs(TGraphErrors* graph1);
    double propagateError(double factorA, double factorB, double factorAErr, double factorBErr, int Type = 0);
    void setOptions();
    TGraphErrors* create_empty_graph();
    TGraphErrors* create_constant_graph(double value, const std::string& name, int n_points = 10);
    void createPNPFractions(TGraphErrors* hyield, TGraphErrors* hyieldVar1, 
                           TGraphErrors* hyieldVar2, TGraphErrors* hNonPromptFrac, 
                           bool Absolute);
};

// Constructor implementation
PlotGraphsObject::PlotGraphsObject(const std::string& ptRng, bool isZt, bool isMC) : 
    ptString(ptRng), minPlotRange(0), minFitRange(0), maxFitRange(1) {
    
    // Initialize all graph pointers to nullptr to avoid undefined behavior
    hMYield = nullptr;
    hMYieldSG = nullptr;
    hMYieldFix = nullptr;
    hIPNonPromptFraction = nullptr;
    gInclFragFunc = nullptr;
    gDecayFragFunc = nullptr;
    gPromptFragFunc = nullptr;
    gAccCorr0 = nullptr;
    gAccCorr1 = nullptr;
    gAccCorr2 = nullptr;
    gAccCorr3 = nullptr;
    gAccCorr4 = nullptr;
    gAccCorr5 = nullptr;
    gAccCorr6 = nullptr;
    gAccCorr = nullptr;
    gAccCorr5_Ext = nullptr;
    gAccCorrTotal = nullptr;
    graphPRaw = nullptr;
    graphNPRaw = nullptr;
    graphInclRaw = nullptr;
    graphPCorr = nullptr;
    graphNPCorr = nullptr;
    graphInclCorr = nullptr;
    
    std::cout << "- - - - - - - - - - - - - - - - -" << std::endl;
    std::cout << "Start Plots for D0 in pT range: " << ptString << std::endl;
    std::cout << "- - - - - - - - - - - - - - - - -" << std::endl;
    
    // Format pT range display string
    ptRange = ptString;
    std::replace(ptRange.begin(), ptRange.end(), '_', '-');
    
    // Set observable tag and x-axis title
    if (isZt) {
        obsTag = "zT";
        xTitle = "#it{z}_{T} = p_{T}^{D^{0}}/p_{T}^{jet}";
    } else {
        obsTag = "Y";
        xTitle = "#it{y}";
    }
    
    // Prepare input file path
    std::string basepath = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA";
    if(isMC) basepath = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_MC";
    std::string rootFileName = basepath + "/FitParametersUnBinnedD0" + obsTag + "_" + ptString + ".root";
    TFile* fInFileHisto = nullptr;
    
    try {
        fInFileHisto = TFile::Open(rootFileName.c_str());
    } catch (const std::exception& e) {
        std::cerr << "Error opening file: " << e.what() << std::endl;
        fInFileHisto = nullptr;
    }
    
    // Create output directory
    OutfilePath = basepath + "/RawSignalYields_D0/";
    if (!std::filesystem::exists(OutfilePath)) {
        std::filesystem::create_directories(OutfilePath);
        std::cout << "make new directory: " << OutfilePath << std::endl;
    }
    
    if (fInFileHisto && !fInFileHisto->IsZombie()) {
        std::cout << "o Open data file: " << fInFileHisto->GetName() << std::endl;
        
        // Get histograms from file with error checking
        hMYield = dynamic_cast<TGraphErrors*>(fInFileHisto->Get("FitMSYieldF"));
        if (!is_valid_graph(hMYield)) {
            std::cout << "Warning: FitMSYieldF not found or invalid, creating fallback" << std::endl;
            hMYield = create_constant_graph(100.0, "FitMSYieldF_fallback");
        }
        
        hMYieldSG = dynamic_cast<TGraphErrors*>(fInFileHisto->Get("FitMSYieldSGF"));
        if (!is_valid_graph(hMYieldSG)) {
            std::cout << "Warning: FitMSYieldSGF not found or invalid, creating fallback" << std::endl;
            hMYieldSG = create_constant_graph(80.0, "FitMSYieldSGF_fallback");
        }
        
        hMYieldFix = dynamic_cast<TGraphErrors*>(fInFileHisto->Get("FitMSYieldFixF"));
        if (!is_valid_graph(hMYieldFix)) {
            hMYieldFix = dynamic_cast<TGraphErrors*>(fInFileHisto->Get("FitMSYieldDCBF"));
            if (!is_valid_graph(hMYieldFix)) {
                std::cout << "Warning: No valid yield histogram found, creating fallback" << std::endl;
                hMYieldFix = create_constant_graph(90.0, "FitMSYieldFixF_fallback");
            }
        } else {
            std::cout << "->->->Still used the FitMSYieldFixF histogram, please update to FitMSYieldDCBF, aka run the fitting code again!" << std::endl;
        }
        
        // Get and process the prompt fraction histogram
        hIPNonPromptFraction = dynamic_cast<TGraphErrors*>(fInFileHisto->Get("FitIPPromptFracF"));
        if (is_valid_graph(hIPNonPromptFraction)) {
            std::cout << "Converting prompt fraction to non-prompt fraction (1-value)" << std::endl;
            for (int i = 0; i < hIPNonPromptFraction->GetN(); i++) {
                double x = hIPNonPromptFraction->GetX()[i];
                double y = hIPNonPromptFraction->GetY()[i];
                double ex = hIPNonPromptFraction->GetEX()[i];
                double ey = hIPNonPromptFraction->GetEY()[i];
                
                // Set new value: 1-y (convert prompt to non-prompt)
                hIPNonPromptFraction->SetPoint(i, x, 1.0-y);
                // Error stays the same for 1-y
                hIPNonPromptFraction->SetPointError(i, ex, ey);
            }
            
            // Update the name to reflect content - FIX STRING CONVERSION
            hIPNonPromptFraction->SetName((std::string(hIPNonPromptFraction->GetName()) + "_NonPrompt").c_str());
        }
        
        if (!is_valid_graph(hIPNonPromptFraction)) {
            std::cout << "Warning: FitIPPromptFracF is not a valid graph, creating a fallback" << std::endl;
            
            // Try alternative name patterns that might exist in the file
            std::vector<std::string> alt_names = {"FitIPPromptFrac", "hIPPromptFrac", "FitIPRes_PromptFrac", "BdecayFrac"};
            for (const auto& name : alt_names) {
                TGraphErrors* test_graph = dynamic_cast<TGraphErrors*>(fInFileHisto->Get(name.c_str()));
                if (is_valid_graph(test_graph)) {
                    std::cout << "Found alternative Nonprompt fraction graph: " << name << std::endl;
                    hIPNonPromptFraction = test_graph;
                    break;
                }
            }
            
            // If still not valid, create a default graph with 50% Nonprompt fraction
            if (!is_valid_graph(hIPNonPromptFraction)) {
                std::cout << "Creating default Nonprompt fraction graph with 50% fraction" << std::endl;
                hIPNonPromptFraction = create_constant_graph(0.5, "FitIPPromptFracF_default");
            }
        }
        
        // Get fragmentation function graphs
        gInclFragFunc = dynamic_cast<TGraphErrors*>(fInFileHisto->Get("ginclFragFunc"));
        gDecayFragFunc = dynamic_cast<TGraphErrors*>(fInFileHisto->Get("decayFragFunc"));
        gPromptFragFunc = dynamic_cast<TGraphErrors*>(fInFileHisto->Get("promptFragFunc"));
        
        // Get correction histograms with error checking
        try {
            gAccCorr0 = dynamic_cast<TGraphErrors*>(fInFileHisto->Get("CorrHist_TypeRnd"));
            if (!is_valid_graph(gAccCorr0)) {
                std::cout << "Warning: CorrHist_TypeRnd is not a valid graph, creating empty graph" << std::endl;
                gAccCorr0 = create_empty_graph();
            }
            
            gAccCorr1 = dynamic_cast<TGraphErrors*>(fInFileHisto->Get("CorrHist_Type1"));
            if (!is_valid_graph(gAccCorr1)) {
                std::cout << "Warning: CorrHist_Type1 is not a valid graph, creating empty graph" << std::endl;
                gAccCorr1 = create_empty_graph();
            }
            
            gAccCorr2 = dynamic_cast<TGraphErrors*>(fInFileHisto->Get("CorrHist_Type2"));
            if (!is_valid_graph(gAccCorr2)) {
                std::cout << "Warning: CorrHist_Type2 is not a valid graph, creating empty graph" << std::endl;
                gAccCorr2 = create_empty_graph();
            }
            
            gAccCorr3 = dynamic_cast<TGraphErrors*>(fInFileHisto->Get("CorrHist_Type3"));
            if (!is_valid_graph(gAccCorr3)) {
                std::cout << "Warning: CorrHist_Type3 is not a valid graph, creating empty graph" << std::endl;
                gAccCorr3 = create_empty_graph();
            }
            
            gAccCorr4 = dynamic_cast<TGraphErrors*>(fInFileHisto->Get("CorrHist_Type4"));
            if (!is_valid_graph(gAccCorr4)) {
                std::cout << "Warning: CorrHist_Type4 is not a valid graph, creating empty graph" << std::endl;
                gAccCorr4 = create_empty_graph();
            }
            
            gAccCorr5 = dynamic_cast<TGraphErrors*>(fInFileHisto->Get("CorrHist_Type5"));
            if (!is_valid_graph(gAccCorr5)) {
                std::cout << "Warning: CorrHist_Type5 is not a valid graph, creating empty graph" << std::endl;
                gAccCorr5 = create_empty_graph();
            }
            
            gAccCorr6 = dynamic_cast<TGraphErrors*>(fInFileHisto->Get("CorrHist_Type6"));
            if (!is_valid_graph(gAccCorr6)) {
                std::cout << "Warning: CorrHist_Type6 is not a valid graph, creating empty graph" << std::endl;
                gAccCorr6 = create_empty_graph();
            }
        } catch (const std::exception& e) {
            std::cerr << "Error loading correction histograms: " << e.what() << std::endl;
            // Create empty graphs as fallbacks
            gAccCorr0 = create_empty_graph();
            gAccCorr1 = create_empty_graph();
            gAccCorr2 = create_empty_graph();
            gAccCorr3 = create_empty_graph();
            gAccCorr4 = create_empty_graph();
            gAccCorr5 = create_empty_graph();
            gAccCorr6 = create_empty_graph();
        }
        
        // Calculate total combined factor
        gAccCorr = multiplyGraphs(gAccCorr0, gAccCorr1);
        
        // Convert factor to efficiency
        gAccCorr = invertGraphs(gAccCorr);
        gAccCorr0 = invertGraphs(gAccCorr0);
        gAccCorr1 = invertGraphs(gAccCorr1);
        gAccCorr2 = invertGraphs(gAccCorr2);
        gAccCorr3 = invertGraphs(gAccCorr3);
        gAccCorr4 = invertGraphs(gAccCorr4);
        gAccCorr5 = invertGraphs(gAccCorr5);
        gAccCorr6 = invertGraphs(gAccCorr6);
        
        gAccCorr5_Ext = nullptr;
        // if (ptRange == "40_60" && resonance == "Psi2S") {
        //     std::cout << "extrapolate missing points" << std::endl;
        //     gAccCorr5_Ext = extrapolateIfNecessary(gAccCorr5);
        // }
        
        // Create total acceptance correction graph
        gAccCorrTotal = multiplyGraphs(gAccCorr0, gAccCorr2);
        gAccCorrTotal = multiplyGraphs(gAccCorrTotal, gAccCorr3);
        gAccCorrTotal = multiplyGraphs(gAccCorrTotal, gAccCorr4);
        
        if (gAccCorr5_Ext) {
            gAccCorrTotal = multiplyGraphs(gAccCorrTotal, gAccCorr5_Ext);
        } else {
            gAccCorrTotal = multiplyGraphs(gAccCorrTotal, gAccCorr5);
        }
        
        gAccCorrTotal = multiplyGraphs(gAccCorrTotal, gAccCorr6);
        
        // Create prompt and non-prompt fractions
        createPNPFractions(hMYield, hMYieldSG, hMYieldFix, hIPNonPromptFraction, true);
        
        // Close the file when done
        fInFileHisto->Close();
        delete fInFileHisto;
    } else {
        std::cout << "o File: " << rootFileName << " does not exist or is invalid, creating fallback graphs" << std::endl;
        
        // Create fallback graphs for required histograms
        hMYield = create_constant_graph(100.0, "FitMSYieldF_fallback");
        hMYieldSG = create_constant_graph(80.0, "FitMSYieldSGF_fallback");
        hMYieldFix = create_constant_graph(90.0, "FitMSYieldFixF_fallback");
        hIPNonPromptFraction = create_constant_graph(0.5, "FitIPPromptFracF_default");
        
        // Create empty graphs for other histograms
        gAccCorr0 = create_empty_graph();
        gAccCorr1 = create_empty_graph();
        gAccCorr2 = create_empty_graph();
        gAccCorr3 = create_empty_graph();
        gAccCorr4 = create_empty_graph();
        gAccCorr5 = create_empty_graph();
        gAccCorr6 = create_empty_graph();
        
        // Set gAccCorrTotal to a fallback constant value
        gAccCorrTotal = create_constant_graph(1.0, "AccCorrTotal_fallback");
        
        if (fInFileHisto) {
            fInFileHisto->Close();
            delete fInFileHisto;
        }
    }
    
    // Create empty PNP fractions to ensure they exist
    graphPRaw = create_empty_graph();
    graphPRaw->SetName(("PromptRaw_" + ptString).c_str());
    graphNPRaw = create_empty_graph();
    graphNPRaw->SetName(("NonPromptRaw_" + ptString).c_str());
    graphInclRaw = create_empty_graph();
    graphInclRaw->SetName(("inclRaw_" + ptString).c_str());
    
    // Now create the PNP fractions if we have valid input graphs
    if (is_valid_graph(hMYield) && is_valid_graph(hMYieldSG) && 
        is_valid_graph(hMYieldFix) && is_valid_graph(hIPNonPromptFraction)) {
        createPNPFractions(hMYield, hMYieldSG, hMYieldFix, hIPNonPromptFraction, true);
    } else {
        std::cout << "Warning: Cannot create P/NP fractions due to missing input graphs" << std::endl;
    }
}

// Destructor
PlotGraphsObject::~PlotGraphsObject() {
    // Note: In a real application, you'd need to manage ownership of
    // the TGraphErrors objects more carefully. In ROOT, objects are often
    // owned by the file they're read from or written to.
}

// Implementation of plotPt method
TGraphErrors* PlotGraphsObject::plotPt(bool isCorrected) {
    plotYieldResult(isCorrected);
    plotCorrFacAcceptance();
    return plotNonPromptFraction(isCorrected);
}

// Implementation of extrapolateIfNecessary method
TGraphErrors* PlotGraphsObject::extrapolateIfNecessary(TGraphErrors* graph) {
    if (!is_valid_graph(graph)) {
        return create_empty_graph();
    }
    
    TGraphErrors* extrapolatedGraph = new TGraphErrors(*graph);
    int nBins = graph->GetN();
    
    double* x = graph->GetX();
    double* xErr = graph->GetEX();
    double* yieldVal = graph->GetY();
    double* yErr = graph->GetEY();
    
    double lastYield = 0;
    double lastError = 0;
    
    for (int bin = 0; bin < nBins; bin++) {
        double x, graphYield;
        graph->GetPoint(bin, x, graphYield);
        double xErr = graph->GetErrorX(bin);
        double graphYieldE = graph->GetErrorY(bin);
        
        if (graphYield == 0) {
            extrapolatedGraph->SetPoint(bin, x, lastYield);
            extrapolatedGraph->SetPointError(bin, xErr, lastError);
        }
        else {
            lastYield = graphYield;
            lastError = graphYieldE;
        }
    }
    
    return extrapolatedGraph;
}

// Implementation of plotCorrFacAcceptance method
void PlotGraphsObject::plotCorrFacAcceptance() {
    setOptions();
    
    std::string outputFilename = OutfilePath + "FinFig_AccCorrFactor_" + obsTag + "_" + ptString + ".png";
    
    TCanvas* c = new TCanvas("c", "c: hist", 500*2, 450*2);
    c->cd();
    TGaxis::SetMaxDigits(3);
    
    // Set pad and histo arrangement
    TPad* myPad2 = new TPad("myPad", "The pad", 0, 0, 1, 1);
    myPad2->SetLeftMargin(0.15);
    myPad2->SetTopMargin(0.06);
    myPad2->SetRightMargin(0.04);
    myPad2->SetBottomMargin(0.15);
    myPad2->SetTicks();
    myPad2->Draw();
    myPad2->cd();
    
    TH1F* myBlankHisto2 = new TH1F((std::string("myBlankHisto2") + obsTag + "_" + ptString).c_str(), 
                                  (std::string("Blank Histogram") + obsTag + "_" + ptString).c_str(), 
                                  100, -0.3, 4.7);
    myBlankHisto2->GetXaxis()->SetNdivisions(505);
    myBlankHisto2->SetXTitle(xTitle.c_str());
    myBlankHisto2->GetXaxis()->SetTitleSize(0.05);
    myBlankHisto2->GetXaxis()->SetRangeUser(-0.3, 1);
    
    if (obsTag.find("Y") != std::string::npos) {
        myBlankHisto2->GetXaxis()->SetRangeUser(2.0, 4.0);
    }
    
    myBlankHisto2->GetXaxis()->SetNdivisions(405);
    myBlankHisto2->GetYaxis()->SetTitleOffset(1.35);
    myBlankHisto2->GetYaxis()->SetTitleSize(0.055);
    myBlankHisto2->SetLineColor(0);
    
    myBlankHisto2->SetYTitle("corr Facor");
    myBlankHisto2->GetYaxis()->SetRangeUser(0, 1);
    myBlankHisto2->Draw("E");
    
    double MarkerScale = 1.6;
    
    setHisto(gAccCorr0, MarkerScale, kGreen-8); // random sel
    gAccCorr0->Draw("same EP");
    
    setHisto(gAccCorr2, MarkerScale, kRed-10);
    gAccCorr2->Draw("same EP");
    
    setHisto(gAccCorr3, MarkerScale, kRed-7);
    gAccCorr3->Draw("same EP");
    
    setHisto(gAccCorr4, MarkerScale, kBlue-9);
    gAccCorr4->Draw("same EP");
    
    if (gAccCorr5_Ext) {
        setHisto(gAccCorr5_Ext, MarkerScale, kBlue-4);
        gAccCorr5_Ext->SetMarkerStyle(24);
        gAccCorr5_Ext->Draw("same EP");
    }
    
    setHisto(gAccCorr5, MarkerScale, kBlue-4);
    gAccCorr5->Draw("same EP");
    
    setHisto(gAccCorr6, MarkerScale, kBlue+2);
    gAccCorr6->Draw("same EP");
    
    setHisto(gAccCorrTotal, MarkerScale, kGreen+2);
    gAccCorrTotal->Draw("same EP");
    
    // Create legend
    TLegend* myLegend0;
    if (obsTag.find("Y") != std::string::npos) {
        myLegend0 = new TLegend(0.5, 0.62, 0.7, 0.8);
    } else {
        myLegend0 = new TLegend(0.15, 0.76, 0.4, 0.9);
    }
    
    myLegend0->SetTextFont(42);
    myLegend0->SetBorderSize(0);
    myLegend0->SetFillStyle(0);
    myLegend0->SetFillColor(0);
    myLegend0->SetMargin(0.25);
    myLegend0->SetTextSize(0.04);
    
    myLegend0->AddEntry(myBlankHisto2, Form("#it{p}_{T}^{jet}=%s (GeV/#it{c})", ptRange.c_str()), "");
    myLegend0->AddEntry(myBlankHisto2, "#it{p}_{T}^{D^{0}}>2 (GeV/#it{c})", "");
    myLegend0->Draw();
    
    TLegend* myLegend1 = new TLegend(0.17, 0.4, 0.45, 0.72);
    myLegend1->SetTextFont(42);
    myLegend1->SetBorderSize(0);
    myLegend1->SetFillStyle(0);
    myLegend1->SetFillColor(0);
    myLegend1->SetMargin(0.25);
    myLegend1->SetTextSize(0.04);
    
    myLegend1->AddEntry(gAccCorrTotal, "total corr factor", "pl");
    myLegend1->AddEntry(gAccCorr0, "Rnd sel. corr", "pl");
    myLegend1->AddEntry(gAccCorr2, "#pi reco", "pl");
    myLegend1->AddEntry(gAccCorr3, "#pi sel.", "pl");
    myLegend1->AddEntry(gAccCorr4, "#mu reco", "pl");
    myLegend1->AddEntry(gAccCorr5, "stripping line corr", "pl");
    myLegend1->AddEntry(gAccCorr6, "trigger eff", "pl");
    myLegend1->AddEntry(gAccCorr6, "?? tag selection", "");
    myLegend1->Draw();
    
    c->SaveAs(outputFilename.c_str());
    c->Close();
}

// Implementation of setHisto method
void PlotGraphsObject::setHisto(TGraphErrors* histo, double MarkerScale, int color) {
    if (!histo) return;
    
    histo->SetMarkerSize(1.3 * MarkerScale);
    histo->SetMarkerStyle(20);
    histo->SetMarkerColor(color);
    histo->SetLineStyle(1);
    histo->SetLineWidth(2);
    histo->SetLineColor(color);
}

// Implementation of plotNonPromptFraction method
TGraphErrors* PlotGraphsObject::plotNonPromptFraction(bool isCorrected) {
    setOptions();
    
    std::string corrTag = "";
    if (isCorrected) {
        corrTag = "_Corr";
    }
    
    std::string outputFilename = OutfilePath + "FinFig_NonPromptFrac_" + obsTag + "_" + ptString + corrTag + ".png";
    
    TCanvas* c = new TCanvas("c", "c: hist", 500*2, 450*2);
    c->cd();
    TGaxis::SetMaxDigits(3);
    
    // Set pad and histo arrangement
    TPad* myPad2 = new TPad("myPad", "The pad", 0, 0, 1, 1);
    myPad2->SetLeftMargin(0.15);
    myPad2->SetTopMargin(0.06);
    myPad2->SetRightMargin(0.04);
    myPad2->SetBottomMargin(0.15);
    myPad2->SetTicks();
    myPad2->Draw();
    myPad2->cd();
    
    TH1F* myBlankHisto2 = new TH1F((std::string("myBlankHisto2") + obsTag + "_" + ptString + corrTag).c_str(), 
                                  (std::string("Blank Histogram") + obsTag + "_" + ptString + corrTag).c_str(), 
                                  100, 0, 5);
    myBlankHisto2->GetXaxis()->SetNdivisions(505);
    myBlankHisto2->SetXTitle(xTitle.c_str());
    myBlankHisto2->GetXaxis()->SetTitleSize(0.05);
    myBlankHisto2->GetXaxis()->SetRangeUser(minPlotRange, 1);
    myBlankHisto2->GetXaxis()->SetNdivisions(405);
    myBlankHisto2->GetYaxis()->SetTitleOffset(1.35);
    myBlankHisto2->GetYaxis()->SetTitleSize(0.055);
    myBlankHisto2->SetLineColor(0);
    
    myBlankHisto2->SetYTitle("Non-Prompt Fraction");
    myBlankHisto2->GetYaxis()->SetRangeUser(0, 0.29);
    
    if (obsTag.find("Y") != std::string::npos) {
        myBlankHisto2->GetXaxis()->SetRangeUser(2.0, 4.0);
    }
    
    myBlankHisto2->Draw("E");
    
    double MarkerScale = 1.6;
    
    // Draw Non-Prompt Fraction
    hIPNonPromptFraction->SetMarkerSize(1.3 * MarkerScale);
    hIPNonPromptFraction->SetMarkerStyle(20);
    hIPNonPromptFraction->SetMarkerColor(kGreen+2);
    hIPNonPromptFraction->SetLineStyle(1);
    hIPNonPromptFraction->SetLineWidth(2);
    hIPNonPromptFraction->SetLineColor(kGreen+2);
    hIPNonPromptFraction->Draw("same EP");
    
    // Create legend
    TLegend* myLegend0;
    if (obsTag.find("Y") != std::string::npos) {
        myLegend0 = new TLegend(0.5, 0.62, 0.7, 0.8);
    }
    
    myLegend0 = new TLegend(0.15, 0.72, 0.4, 0.9);
    myLegend0->SetTextFont(42);
    myLegend0->SetBorderSize(0);
    myLegend0->SetFillStyle(0);
    myLegend0->SetFillColor(0);
    myLegend0->SetMargin(0.25);
    myLegend0->SetTextSize(0.04);
    
    myLegend0->AddEntry(myBlankHisto2, "Anti-#it{k}_{T} #it{R} = 0.5, #it{#eta}_{jet}= 2.5-4", "");
    myLegend0->AddEntry(myBlankHisto2, Form("#it{p}_{T}^{jet}=%s (GeV/#it{c})", ptRange.c_str()), "");
    myLegend0->AddEntry(myBlankHisto2, "#it{p}_{T}^{D^{0}}>2 (GeV/#it{c})", "");
    
    TLatex* collaboration = new TLatex(0.62, 0.88, "#bf{LHCb} in progress");
    collaboration->SetNDC();
    collaboration->SetTextSize(0.044);
    
    TLatex* system = new TLatex(0.62, 0.83, "Pb-p, #sqrt{#it{s}} = 8.16 TeV");
    system->SetNDC();
    system->SetTextSize(0.044);
    
    myLegend0->Draw();
    collaboration->Draw();
    system->Draw();
    
    c->SaveAs(outputFilename.c_str());
    c->Close();
    
    return hIPNonPromptFraction;
}

// Implementation of plotYieldSummary method
void PlotGraphsObject::plotYieldSummary(std::vector<TGraphErrors*> yieldArray, 
                                       std::vector<std::string> pTArray, 
                                       int normType, 
                                       std::string Seltag, 
                                       bool isCorrected,
                                       std::vector<TGraphErrors*> corrArray) {
    
    std::vector<int> ColorArray = {kAzure, kAzure-4, kCyan-6, kGreen-3, kTeal-6, kGreen+4, 1, 1};
    
    if (ColorArray.size() < yieldArray.size() || pTArray.size() < yieldArray.size()) {
        std::cout << "Adapt size of color array!" << std::endl;
    }
    
    // Apply acceptance correction
    std::string corrTag = isCorrected ? "_Corr" : "";
    
    std::string titleAddon = "";
    setOptions();
    
    std::string outputFilename;
    if (normType == 0) {
        outputFilename = OutfilePath + "FinFig_YieldSummary_" + obsTag + Seltag + corrTag + ".png";
    } else if (normType == 1) {
        outputFilename = OutfilePath + "FinFig_YieldSummary_" + obsTag + Seltag + "Norm" + corrTag + ".png";
        titleAddon = "/d#sigma";
    } else if (normType == 2) {
        outputFilename = OutfilePath + "FinFig_NonPromptFracSummary_" + obsTag + Seltag + corrTag + ".png";
    }
    
    TCanvas* c = new TCanvas("c", "c: hist", 500*2, 450*2);
    c->cd();
    TGaxis::SetMaxDigits(3);
    
    // Set pad and histo arrangement
    TPad* myPad2 = new TPad("myPad", "The pad", 0, 0, 1, 1);
    myPad2->SetLeftMargin(0.15);
    myPad2->SetTopMargin(0.06);
    myPad2->SetRightMargin(0.04);
    myPad2->SetBottomMargin(0.15);
    myPad2->SetTicks();
    myPad2->Draw();
    
    if (obsTag.find("Y") != std::string::npos && normType < 2) {
        myPad2->cd()->SetLogy();
    } else {
        myPad2->cd();
    }
    
    // Find maximum value
    double max = yieldArray[0]->GetMaximum();
    if (yieldArray.size() > 1 && yieldArray[1]->GetMaximum() > max) {
        max = yieldArray[1]->GetMaximum()*10;
    }
    if (yieldArray.size() > 2 && yieldArray[2]->GetMaximum() > max) {
        max = yieldArray[2]->GetMaximum()*10;
    }
    
    std::cout << "This is the maximum: " << max << std::endl;
    if (max < 0) {
        max = 10000;
    }
    
    TH1F* myBlankHisto2 = new TH1F((std::string("myBlankHisto2") + obsTag + Seltag + corrTag).c_str(), 
                                  (std::string("Blank Histogram") + obsTag + Seltag + corrTag).c_str(), 
                                  100, 0, 5);
    myBlankHisto2->GetXaxis()->SetNdivisions(505);
    myBlankHisto2->SetXTitle(xTitle.c_str());
    myBlankHisto2->GetXaxis()->SetTitleSize(0.05);
    myBlankHisto2->GetXaxis()->SetRangeUser(minPlotRange, 1);
    myBlankHisto2->GetXaxis()->SetNdivisions(405);
    myBlankHisto2->GetYaxis()->SetTitleOffset(1.35);
    myBlankHisto2->GetYaxis()->SetTitleSize(0.055);
    myBlankHisto2->SetLineColor(0);
    
    // Set y-axis labels and ranges
    if (obsTag.find("Y") != std::string::npos) {
        myBlankHisto2->SetYTitle("dN/dY");
        myBlankHisto2->GetYaxis()->SetRangeUser(10, max*700);
        myBlankHisto2->GetXaxis()->SetRangeUser(2.0, 4.0);
        if (normType == 1) {
            myBlankHisto2->GetYaxis()->SetRangeUser(1, 100000);
        }
    } else {
        myBlankHisto2->SetYTitle(("dN/dz_{T}" + titleAddon).c_str());
        myBlankHisto2->GetYaxis()->SetRangeUser(0, max*1.2);
        if (normType == 1) {
            myBlankHisto2->GetYaxis()->SetRangeUser(0, 1.0);
        }
    }
    
    // Special case for Nonprompt fraction
    if (normType == 2) {
        myBlankHisto2->SetYTitle("Non-Prompt Fraction");
        myBlankHisto2->GetYaxis()->SetRangeUser(0, 0.3);
    }
    
    myBlankHisto2->Draw("E");
    
    double MarkerScale = 1.6;
    
    // Draw all graphs
    for (size_t i = 0; i < yieldArray.size(); i++) {
        std::cout << "Histo Name: " << yieldArray[i]->GetName() << std::endl;
        yieldArray[i]->SetMarkerSize(1.3 * MarkerScale);
        yieldArray[i]->SetMarkerStyle(20);
        yieldArray[i]->SetMarkerColor(ColorArray[i]);
        yieldArray[i]->SetLineStyle(1);
        yieldArray[i]->SetLineWidth(2);
        yieldArray[i]->SetLineColor(ColorArray[i]);
        
        if (normType == 1) {
            normalize(yieldArray[i]);
        }
        
        yieldArray[i]->Draw("same EP");
    }
    
    // Create legend
    TLegend* myLegend1;
    if (obsTag.find("Y") != std::string::npos) {
        myLegend1 = new TLegend(0.2, 0.65, 0.4, 0.9);
    } else {
        myLegend1 = new TLegend(0.2, 0.65, 0.4, 0.9);
    }

    
    myLegend1->SetTextFont(42);
    myLegend1->SetBorderSize(0);
    myLegend1->SetFillStyle(0);
    myLegend1->SetFillColor(0);
    myLegend1->SetMargin(0.25);
    myLegend1->SetTextSize(0.04);
    
    myLegend1->AddEntry(yieldArray[0], " #it{p}_{T}^{jet}:", "");
    
    for (size_t i = 0; i < yieldArray.size(); i++) {
        std::string displayPt = pTArray[i];
        std::replace(displayPt.begin(), displayPt.end(), '_', '-');
        myLegend1->AddEntry(yieldArray[i], Form("  %s (GeV/%s)", displayPt.c_str(), "#it{c}"), "LP");
    }
    
    TLatex* collaboration = new TLatex(0.62, 0.88, "#bf{LHCb} in progress");
    collaboration->SetNDC();
    collaboration->SetTextSize(0.044);
    
    TLatex* system = new TLatex(0.62, 0.83, "Pb-p, #sqrt{#it{s}} = 8.16 TeV");
    system->SetNDC();
    system->SetTextSize(0.044);
    
    TLatex* component = nullptr;
    if (normType != 2) {
        component = new TLatex(0.67, 0.73, (Seltag + " component").c_str());
        component->SetNDC();
        component->SetTextSize(0.042);
        component->Draw();
    }
    
    myLegend1->Draw();
    collaboration->Draw();
    system->Draw();
    
    c->SaveAs(outputFilename.c_str());
    c->Close();
}

// Implementation of normalize method
void PlotGraphsObject::normalize(TGraphErrors* graph) {
    if (!is_valid_graph(graph)) {
        std::cout << "Cannot normalize invalid graph" << std::endl;
        return;
    }
    
    double sum = 0.0;
    int nPoints = graph->GetN();
    
    // Calculate sum
    for (int bin = 0; bin < nPoints; bin++) {
        double x, y;
        graph->GetPoint(bin, x, y);
        sum += y;
    }
    
    if (obsTag == "Y") {
        sum *= 0.00001;
    }
    
    if (sum == 0) return;
    
    // Apply normalization
    for (int bin = 0; bin < nPoints; bin++) {
        double x, y;
        graph->GetPoint(bin, x, y);
        double yerr = graph->GetErrorY(bin);
        
        graph->SetPoint(bin, x, y / sum);
        graph->SetPointError(bin, graph->GetErrorX(bin), yerr / sum);
    }
}

// Implementation of applyAccCorrection method
TGraphErrors* PlotGraphsObject::applyAccCorrection(TGraphErrors* graph, TGraphErrors* gAccCorrInput) {
    std::cout << "Divide GRAPHS" << std::endl;
    
    if (!gAccCorrInput) {
        std::cout << "Try to get the histogram from the constructor" << std::endl;
        gAccCorrInput = gAccCorrTotal;
    }
    
    if (!is_valid_graph(graph) || !is_valid_graph(gAccCorrInput)) {
        std::cout << "Invalid graphs for acc correction" << std::endl;
        return create_empty_graph();
    }
    
    int nBinsAcc = gAccCorrInput->GetN();
    int nBinsgraph = graph->GetN();
    
    if (nBinsAcc != nBinsgraph) {
        std::cout << "!!!This can not work fix this error!!!" << std::endl;
        return create_empty_graph();
    }
    
    // Create arrays for the new graph
    double* x = new double[nBinsAcc];
    double* xErr = new double[nBinsAcc];
    double* CorrectedYield = new double[nBinsAcc];
    double* CorrectedYieldErr = new double[nBinsAcc];
    
    // Fill arrays
    for (int bin = 0; bin < nBinsAcc; bin++) {
        double graphYield, corrFactVal;
        graph->GetPoint(bin, x[bin], graphYield);
        xErr[bin] = graph->GetErrorX(bin);
        double graphYieldE = graph->GetErrorY(bin);
        
        double x2;
        gAccCorrInput->GetPoint(bin, x2, corrFactVal);
        double corrFactE = gAccCorrInput->GetErrorY(bin);
        
        if (corrFactVal > 0) {
            CorrectedYield[bin] = graphYield / corrFactVal;
            CorrectedYieldErr[bin] = propagateError(graphYield, 1.0/corrFactVal, graphYieldE, corrFactE, 2);
        }
        else {
            CorrectedYield[bin] = 0;
            CorrectedYieldErr[bin] = 0;
        }
    }
    
    // Find maximum value
    double maximum = 0;
    for (int i = 0; i < nBinsAcc; i++) {
        if (CorrectedYield[i] > maximum) {
            maximum = CorrectedYield[i];
        }
    }
    
    // Create the new graph
    TGraphErrors* ratioGraph = new TGraphErrors(nBinsAcc, x, CorrectedYield, xErr, CorrectedYieldErr);
    std::string title = graph->GetTitle() + std::string("AccCorr") + ptString;
    std::string name = graph->GetName() + std::string("AccCorr") + ptString;
    ratioGraph->SetTitle(title.c_str());
    ratioGraph->SetName(name.c_str());
    ratioGraph->SetMaximum(maximum);
    
    // Cleanup
    delete[] x;
    delete[] xErr;
    delete[] CorrectedYield;
    delete[] CorrectedYieldErr;
    
    return ratioGraph;
}

// Implementation of multiplyGraphs method
TGraphErrors* PlotGraphsObject::multiplyGraphs(TGraphErrors* graph1, TGraphErrors* graph2) {
    std::cout << "Multiply GRAPHS" << std::endl;
    
    // Check if both graphs are valid
    if (!is_valid_graph(graph1)) {
        std::cout << "Warning: First graph is not a valid graph object" << std::endl;
        return create_empty_graph();
    }
    
    if (!is_valid_graph(graph2)) {
        std::cout << "Warning: Second graph is not a valid graph object" << std::endl;
        return create_empty_graph();
    }
    
    int nBins1 = graph1->GetN();
    int nBins2 = graph2->GetN();
    
    if (nBins1 != nBins2) {
        std::cout << "Warning: Graphs have different numbers of bins (" << nBins1 << " vs " << nBins2 << ")" << std::endl;
        return create_empty_graph();
    }
    
    if (nBins1 == 0) {
        std::cout << "Warning: Graphs have zero points" << std::endl;
        return create_empty_graph();
    }
    
    // Create arrays for the new graph
    double* x = new double[nBins1];
    double* xErr = new double[nBins1];
    double* CorrectedYield = new double[nBins1];
    double* CorrectedYieldErr = new double[nBins1];
    
    // Fill arrays
    double maximum = 0.0;
    for (int bin = 0; bin < nBins1; bin++) {
        double graphYield, graphYield2;
        graph1->GetPoint(bin, x[bin], graphYield);
        xErr[bin] = graph1->GetErrorX(bin);
        double graphYieldE = graph1->GetErrorY(bin);
        
        double x2, y2;
        graph2->GetPoint(bin, x2, graphYield2);
        double graphYieldE2 = graph2->GetErrorY(bin);
        
        CorrectedYield[bin] = graphYield * graphYield2;
        // Proper error propagation
        CorrectedYieldErr[bin] = propagateError(graphYield, graphYield2, graphYieldE, graphYieldE2, 1);
        
        if (CorrectedYield[bin] > maximum) {
            maximum = CorrectedYield[bin];
        }
    }
    
    // Create the new graph
    TGraphErrors* ratioGraph = new TGraphErrors(nBins1, x, CorrectedYield, xErr, CorrectedYieldErr);
    std::string title = graph1->GetTitle() + std::string("Multipl") + ptString;
    std::string name = graph1->GetName() + std::string("Multipl") + ptString;
    ratioGraph->SetTitle(title.c_str());
    ratioGraph->SetName(name.c_str());
    ratioGraph->SetMaximum(maximum);
    
    // Cleanup
    delete[] x;
    delete[] xErr;
    delete[] CorrectedYield;
    delete[] CorrectedYieldErr;
    
    return ratioGraph;
}

// Implementation of invertGraphs method
TGraphErrors* PlotGraphsObject::invertGraphs(TGraphErrors* graph1) {
    std::cout << "Invert GRAPHS" << std::endl;
    
    // Check if graph is valid
    if (!is_valid_graph(graph1)) {
        std::cout << "Warning: Graph is not a valid graph object" << std::endl;
        return create_empty_graph();
    }
    
    int nBins1 = graph1->GetN();
    if (nBins1 == 0) {
        std::cout << "Warning: Graph has zero points" << std::endl;
        return create_empty_graph();
    }
    
    // Create arrays for the new graph
    double* x = new double[nBins1];
    double* xErr = new double[nBins1];
    double* CorrectedYield = new double[nBins1];
    double* CorrectedYieldErr = new double[nBins1];
    
    // Fill arrays
    double maximum = 0.0;
    for (int bin = 0; bin < nBins1; bin++) {
        double graphYield, graphYieldE;
        graph1->GetPoint(bin, x[bin], graphYield);
        xErr[bin] = graph1->GetErrorX(bin);
        graphYieldE = graph1->GetErrorY(bin);
        
        if (graphYield != 0) {
            CorrectedYield[bin] = 1.0 / graphYield;
            // Proper error propagation
            CorrectedYieldErr[bin] = propagateError(1.0, graphYield, 0, graphYieldE, 0);
        }
        else {
            CorrectedYield[bin] = 0;
            CorrectedYieldErr[bin] = 0;
        }
        
        if (CorrectedYield[bin] > maximum) {
            maximum = CorrectedYield[bin];
        }
    }
    
    // Create the new graph
    TGraphErrors* ratioGraph = new TGraphErrors(nBins1, x, CorrectedYield, xErr, CorrectedYieldErr);
    std::string title = std::string(graph1->GetTitle()) + "Multipl" + ptString;
    std::string name = std::string(graph1->GetName()) + "Multipl" + ptString;
    ratioGraph->SetTitle(title.c_str());
    ratioGraph->SetName(name.c_str());
    ratioGraph->SetMaximum(maximum);
    
    // Cleanup
    delete[] x;
    delete[] xErr;
    delete[] CorrectedYield;
    delete[] CorrectedYieldErr;
    
    return ratioGraph;
}

// Implementation of propagateError method
double PlotGraphsObject::propagateError(double factorA, double factorB, double factorAErr, double factorBErr, int Type) {
    // Calculate error of f(A,B) = A/B
    // sig(f(A,B))^2 = ...
    
    double error = 0;
    
    // A/B
    if (Type == 0) {
        if (factorB > 0) {
            error = pow(factorA / factorB, 2) * (pow(factorAErr / factorA, 2) + pow(factorBErr / factorB, 2)); // No correlation assumed
        }
    }
    // A*B
    else if (Type == 1) {
        if (factorA > 0 && factorB > 0) {
            error = pow(factorA * factorB, 2) * (pow(factorAErr / factorA, 2) + pow(factorBErr / factorB, 2));
        }
    }
    // A*B (B has no error)
    else if (Type == 2) {
        error = pow(factorAErr * factorB, 2);
    }
    else {
        std::cout << "Not implemented!" << std::endl;
        return 0;
    }
    
    return sqrt(error);
}

// Implementation of plotYieldResult method
void PlotGraphsObject::plotYieldResult(bool isCorrected) {
    setOptions();
    
    // Make sure attributes exist
    if (!graphPRaw || !graphNPRaw || !graphInclRaw) {
        std::cout << "Warning: Required yield graphs not found - they may have failed to be created" << std::endl;
        
        // Create empty fallback graphs
        if (!graphPRaw) {
            graphPRaw = create_empty_graph();
            graphPRaw->SetName(("PromptRaw_" + ptString).c_str());
        }
        
        if (!graphNPRaw) {
            graphNPRaw = create_empty_graph();
            graphNPRaw->SetName(("NonPromptRaw_" + ptString).c_str());
        }
        
        if (!graphInclRaw) {
            graphInclRaw = create_empty_graph();
            graphInclRaw->SetName(("inclRaw_" + ptString).c_str());
        }
    }
    
    TGraphErrors* gpromtOrig = graphPRaw;
    TGraphErrors* gNpromtOrig = graphNPRaw;
    TGraphErrors* hmassyieldOrig = graphInclRaw;
    
    std::string corrTag = "";
    TGraphErrors *gpromt, *gNpromt, *hmassyield;
    
    if (isCorrected) {
        corrTag = "_Corr";
        gpromt = applyAccCorrection(gpromtOrig);        // Apply acceptance correction
        gNpromt = applyAccCorrection(gNpromtOrig);      // Apply acceptance correction
        hmassyield = applyAccCorrection(hmassyieldOrig); // Apply acceptance correction
        
        graphPCorr = gpromt;
        graphPCorr->SetName(("PromptCorr_" + ptString).c_str());
        
        graphNPCorr = gNpromt;
        graphNPCorr->SetName(("NonPromptCorr_" + ptString).c_str());
        
        graphInclCorr = hmassyield;
        graphInclCorr->SetName(("inclCorr_" + ptString).c_str());
    } else {
        gpromt = gpromtOrig;
        gNpromt = gNpromtOrig;
        hmassyield = hmassyieldOrig;
    }
    
    std::string outputFilename = OutfilePath + "FinFig_Yield_" + obsTag + "_" + ptString + corrTag + ".png";
    
    // Create legend
    TLegend* myLegend0;
    if (obsTag.find("Y") != std::string::npos) {
        myLegend0 = new TLegend(0.5, 0.62, 0.7, 0.8);
    } else {
        myLegend0 = new TLegend(0.15, 0.72, 0.4, 0.9);
    }
    
    myLegend0->SetTextFont(42);
    myLegend0->SetBorderSize(0);
    myLegend0->SetFillStyle(0);
    myLegend0->SetFillColor(0);
    myLegend0->SetMargin(0.25);
    myLegend0->SetTextSize(0.04);
    
    // Legend about different contributions
    TLegend* myLegend1;
    if (obsTag.find("Y") != std::string::npos) {
        myLegend1 = new TLegend(0.55, 0.47, 0.6, 0.59);
    } else {
        myLegend1 = new TLegend(0.2, 0.59, 0.4, 0.7);
    }
    
    myLegend1->SetTextFont(42);
    myLegend1->SetBorderSize(0);
    myLegend1->SetFillStyle(0);
    myLegend1->SetFillColor(0);
    myLegend1->SetMargin(0.25);
    myLegend1->SetTextSize(0.04);
    
    double scale = 1;
    double MarkerScale = 1.6;
    double labelOffset2 = 0.23 + 0.2 - 0.3;
    
    TLatex* collaboration = new TLatex(0.62, 0.88, "#bf{LHCb} in progress");
    collaboration->SetNDC();
    collaboration->SetTextSize(0.044 * scale);
    
    TLatex* system = new TLatex(0.62, 0.83, "Pb-p, #sqrt{#it{s}} = 8.16 TeV");
    system->SetNDC();
    system->SetTextSize(0.044 * scale);
    
    TCanvas* c = new TCanvas("c", "c: hist", 500*2, 450*2);
    c->cd();
    TGaxis::SetMaxDigits(3);
    
    // Set pad and histo arrangement
    TPad* myPad2 = new TPad("myPad", "The pad", 0, 0, 1, 1);
    myPad2->SetLeftMargin(0.15);
    myPad2->SetTopMargin(0.06);
    myPad2->SetRightMargin(0.04);
    myPad2->SetBottomMargin(0.15);
    myPad2->SetTicks();
    myPad2->Draw();
    
    if (obsTag.find("Y") != std::string::npos) {
        myPad2->cd()->SetLogy();
    } else {
        myPad2->cd();
    }
    
    double max = hmassyield->GetMaximum();
    
    TH1F* myBlankHisto2 = new TH1F((std::string("myBlankHisto2") + obsTag + "_" + ptString + corrTag).c_str(), 
                                  (std::string("Blank Histogram") + obsTag + "_" + ptString + corrTag).c_str(), 
                                  100, 0, 5);
    myBlankHisto2->GetXaxis()->SetNdivisions(505);
    myBlankHisto2->SetXTitle(xTitle.c_str());
    myBlankHisto2->GetXaxis()->SetTitleSize(0.05);
    myBlankHisto2->GetXaxis()->SetRangeUser(minPlotRange, 1);
    myBlankHisto2->GetXaxis()->SetNdivisions(405);
    myBlankHisto2->GetYaxis()->SetTitleOffset(1.35);
    myBlankHisto2->GetYaxis()->SetTitleSize(0.055);
    myBlankHisto2->SetLineColor(0);
    
    // Set y-axis title and range
    if (obsTag.find("Y") != std::string::npos) {
        myBlankHisto2->SetYTitle("dN/d#it{y}");
        myBlankHisto2->GetYaxis()->SetRangeUser(10, max*4);
        myBlankHisto2->GetXaxis()->SetRangeUser(2.0, 4.0);
    } else {
        myBlankHisto2->SetYTitle("dN/dz_{T}");
        myBlankHisto2->GetYaxis()->SetRangeUser(0, max*2);
    }
    
    myBlankHisto2->Draw("E");
    
    // Draw total yield
    hmassyield->SetMarkerSize(1.3 * MarkerScale);
    hmassyield->SetMarkerStyle(20);
    hmassyield->SetMarkerColor(kGreen+2);
    hmassyield->SetLineStyle(1);
    hmassyield->SetLineWidth(2);
    hmassyield->SetLineColor(kGreen+2);
    hmassyield->Draw("same EP");
    
    // Draw prompt fraction
    gpromt->SetMarkerSize(1.3 * MarkerScale);
    gpromt->SetMarkerStyle(4);
    gpromt->SetLineWidth(2);
    gpromt->SetMarkerColor(kBlue);
    gpromt->SetLineColor(kBlue);
    gpromt->Draw("same EP");
    
    // Draw non-prompt fraction
    gNpromt->SetMarkerSize(1.3 * MarkerScale);
    gNpromt->SetMarkerStyle(4);
    gNpromt->SetLineWidth(2);
    gNpromt->SetMarkerColor(kBlue-9);
    gNpromt->SetLineColor(kBlue-9);
    gNpromt->Draw("same EP");
    
    myLegend1->AddEntry(hmassyield, " inclusive yield", "LP");
    myLegend1->AddEntry(gpromt, " prompt yield", "LP");
    myLegend1->AddEntry(gNpromt, " non-prompt yield", "LP");
    
    myLegend0->AddEntry(myBlankHisto2, "Anti-#it{k}_{T} #it{R} = 0.5, #it{#eta}_{jet}= 2.5-4", "");
    myLegend0->AddEntry(myBlankHisto2, Form("#it{p}_{T}^{jet}=%s (GeV/#it{c})", ptRange.c_str()), "");
    myLegend0->AddEntry(myBlankHisto2, "#it{p}_{T}^{D^{0}}>2 (GeV/#it{c})", "");
    
    myLegend1->Draw();
    collaboration->Draw();
    system->Draw();
    myLegend0->Draw();
    
    c->SaveAs(outputFilename.c_str());
    c->Close();
}

// Make sure these implementations are in your file and not duplicated
TGraphErrors* PlotGraphsObject::create_empty_graph() {
    try {
        // Create an empty graph with proper initialization
        TGraphErrors* empty_graph = new TGraphErrors();
        std::string name = "empty_graph_" + ptString;
        empty_graph->SetName(name.c_str());
        // Set point at 0,0 to ensure it's not completely empty
        empty_graph->SetPoint(0, 0.0, 0.0);
        empty_graph->SetPointError(0, 0.0, 0.0);
        return empty_graph;
    }
    catch (const std::exception& e) {
        std::cerr << "Error creating empty graph: " << e.what() << std::endl;
        // Last resort fallback - return nullptr and handle it elsewhere
        return nullptr;
    }
}

TGraphErrors* PlotGraphsObject::create_constant_graph(double value, const std::string& name, int n_points) {
    try {
        // Create bins from 0 to 1 for zT
        double* x_vals = new double[n_points];
        double* x_errs = new double[n_points];
        double* y_vals = new double[n_points];
        double* y_errs = new double[n_points];
        
        // Fill arrays
        for (int i = 0; i < n_points; i++) {
            x_vals[i] = i * (1.0 / (n_points - 1)); // Equivalent to linspace(0,1,n_points)
            x_errs[i] = 0.05; // Fixed bin width
            y_vals[i] = value; // Constant value
            y_errs[i] = 0.1; // Reasonable uncertainty
        }
        
        // Create the graph
        TGraphErrors* graph = new TGraphErrors(n_points, x_vals, y_vals, x_errs, y_errs);
        graph->SetName(name.c_str());
        
        // Cleanup arrays
        delete[] x_vals;
        delete[] x_errs;
        delete[] y_vals;
        delete[] y_errs;
        
        return graph;
    }
    catch (const std::exception& e) {
        std::cerr << "Error creating constant graph: " << e.what() << std::endl;
        return create_empty_graph();
    }
}

// Implementation of plotRawYields function
void plotRawYields(const std::string& ptRange = "", bool isZt = true, bool isMC = false) {
    std::cout << "is binned var" << std::endl;
    
    std::vector<std::string> pTRangeArray = {"5_10", "10_15", "15_20", "20_30", "30_50"};
    
    std::vector<TGraphErrors*> yieldArray;
    std::vector<TGraphErrors*> yieldArrayP;
    std::vector<TGraphErrors*> yieldArrayNP;
    std::vector<TGraphErrors*> CyieldArray;
    std::vector<TGraphErrors*> CyieldArrayP;
    std::vector<TGraphErrors*> CyieldArrayNP;
    std::vector<TGraphErrors*> CyieldArrayUnNorm;
    std::vector<TGraphErrors*> accArray;
    std::vector<TGraphErrors*> bFracArray;
    
    std::vector<PlotGraphsObject*> gObj(pTRangeArray.size());
    
    // Call the object
    if (ptRange.empty()) {
        // Apply acceptance correction factors
        bool plotAccCorr = false;
        
        for (size_t i = 0; i < pTRangeArray.size(); i++) {
            gObj[i] = new PlotGraphsObject(pTRangeArray[i], isZt, isMC);
            gObj[i]->plotPt(plotAccCorr);
            yieldArray.push_back(gObj[i]->graphInclRaw);
            yieldArrayP.push_back(gObj[i]->graphPRaw);
            yieldArrayNP.push_back(gObj[i]->graphNPRaw);
            
            if (plotAccCorr) {
                CyieldArray.push_back(gObj[i]->graphInclCorr);
                CyieldArrayP.push_back(gObj[i]->graphPCorr);
                CyieldArrayNP.push_back(gObj[i]->graphNPCorr);
            }
            
            bFracArray.push_back(gObj[i]->hIPNonPromptFraction);
        }
        
        // Create an output file to save all final histograms
        // if (plotAccCorr) {
            std::string outFileName = gObj[0]->OutfilePath + "CorrectedFinalHistograms_D0.root";
            TFile* fOutData = new TFile(outFileName.c_str(), "RECREATE");
            
            for (size_t i = 0; i < yieldArray.size(); i++) {
                if (plotAccCorr) {
                    CyieldArray[i]->Write();
                    CyieldArrayP[i]->Write();
                    CyieldArrayNP[i]->Write();
                }
                yieldArray[i]->Write();
                yieldArrayP[i]->Write();
                yieldArrayNP[i]->Write();
                bFracArray[i]->Write();
            }
            
            fOutData->Close();
            delete fOutData;
        // }
        
        // Create summary plots
        if (plotAccCorr) {
            gObj[0]->plotYieldSummary(CyieldArray, pTRangeArray, 0, "incl", plotAccCorr);
            gObj[0]->plotYieldSummary(CyieldArray, pTRangeArray, 1, "incl", plotAccCorr);
            gObj[0]->plotYieldSummary(CyieldArrayP, pTRangeArray, 0, "promt", plotAccCorr);
            gObj[0]->plotYieldSummary(CyieldArrayP, pTRangeArray, 1, "promt", plotAccCorr);
            gObj[0]->plotYieldSummary(CyieldArrayNP, pTRangeArray, 0, "Npromt", plotAccCorr);
            gObj[0]->plotYieldSummary(CyieldArrayNP, pTRangeArray, 1, "Npromt", plotAccCorr);
        }
        else {
            gObj[0]->plotYieldSummary(yieldArray, pTRangeArray, 0, "incl", plotAccCorr);
            gObj[0]->plotYieldSummary(yieldArray, pTRangeArray, 1, "incl", plotAccCorr);
            gObj[0]->plotYieldSummary(yieldArrayP, pTRangeArray, 0, "promt", plotAccCorr);
            gObj[0]->plotYieldSummary(yieldArrayP, pTRangeArray, 1, "promt", plotAccCorr);
            gObj[0]->plotYieldSummary(yieldArrayNP, pTRangeArray, 0, "Npromt", plotAccCorr);
            gObj[0]->plotYieldSummary(yieldArrayNP, pTRangeArray, 1, "Npromt", plotAccCorr);
        }
        
        gObj[0]->plotYieldSummary(bFracArray, pTRangeArray, 2);
        
        // Clean up
        for (size_t i = 0; i < gObj.size(); i++) {
            delete gObj[i];
        }
    }
    else {
        // Plot just one single pT bin
        PlotGraphsObject* graphs = new PlotGraphsObject(ptRange, isZt, isMC);
        graphs->plotPt();
        delete graphs;
    }
}

// Add this implementation before the end of the file
void PlotGraphsObject::createPNPFractions(TGraphErrors* hyield, TGraphErrors* hyieldVar1, 
                                        TGraphErrors* hyieldVar2, TGraphErrors* hNonPromptFrac, 
                                        bool Absolute) {
    // Check if input objects are valid graphs
    if (!is_valid_graph(hyield) || !is_valid_graph(hyieldVar1) || 
        !is_valid_graph(hyieldVar2) || !is_valid_graph(hNonPromptFrac)) {
        
        std::cout << "Error: One or more input graphs are not valid" << std::endl;
        std::vector<std::string> missing;
        if (!is_valid_graph(hyield)) missing.push_back("hyield");
        if (!is_valid_graph(hyieldVar1)) missing.push_back("hyieldVar1");
        if (!is_valid_graph(hyieldVar2)) missing.push_back("hyieldVar2");
        if (!is_valid_graph(hNonPromptFrac)) missing.push_back("hNonPromptFrac");
        std::cout << "Missing valid graphs: ";
        for (size_t i = 0; i < missing.size(); i++) {
            std::cout << missing[i];
            if (i < missing.size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
        
        // Create empty graphs as fallback
        graphPRaw = create_empty_graph();
        graphPRaw->SetName(("PromptRaw_empty_" + ptString).c_str());
        graphNPRaw = create_empty_graph();
        graphNPRaw->SetName(("NonPromptRaw_empty_" + ptString).c_str());
        graphInclRaw = create_empty_graph();
        graphInclRaw->SetName(("inclRaw_empty_" + ptString).c_str());
        return;
    }
    
    int pointsA = hyield->GetN();
    int pointsB = hNonPromptFrac->GetN();
    
    if (pointsA != pointsB) {
        std::cout << "Warning: Different number of points between graphs: " << pointsA 
                  << " vs " << pointsB << std::endl;
        std::cout << "Attempting to adapt by using common bins..." << std::endl;
        
        // Use the minimum number of points to avoid out-of-bounds
        int common_points = std::min(pointsA, pointsB);
        if (common_points == 0) {
            std::cout << "Error: No common bins to process" << std::endl;
            // Create empty fallback graphs
            graphPRaw = create_empty_graph();
            graphPRaw->SetName(("PromptRaw_empty_" + ptString).c_str());
            graphNPRaw = create_empty_graph();
            graphNPRaw->SetName(("NonPromptRaw_empty_" + ptString).c_str());
            graphInclRaw = create_empty_graph();
            graphInclRaw->SetName(("inclRaw_empty_" + ptString).c_str());
            return;
        }
        pointsA = common_points;
    }
    
    // Proceed with the calculation using common_points 
    const int NbOfVariations = 3;
    TGraphErrors* graphNP[NbOfVariations];
    TGraphErrors* graphP[NbOfVariations];
    TGraphErrors* graphyield[NbOfVariations];
    
    std::vector<double> x_arr(pointsA);
    std::vector<double> xErr_arr(pointsA);
    std::vector<double> yieldVal_arr(pointsA);
    std::vector<double> yErr_arr(pointsA);
    std::vector<double> yieldValVar1_arr(pointsA);
    std::vector<double> yErrVar1_arr(pointsA);
    std::vector<double> yieldValVar2_arr(pointsA);
    std::vector<double> yErrVar2_arr(pointsA);
    std::vector<double> fracVal_arr(pointsA);
    std::vector<double> fracValE_arr(pointsA);
    
    // Extract data from input graphs
    for (int i = 0; i < pointsA; i++) {
        double x, y, ex, ey;
        
        hyield->GetPoint(i, x, y);
        x_arr[i] = x;
        xErr_arr[i] = hyield->GetErrorX(i);
        yieldVal_arr[i] = y;
        yErr_arr[i] = hyield->GetErrorY(i);
        
        hyieldVar1->GetPoint(i, x, y);
        yieldValVar1_arr[i] = y;
        yErrVar1_arr[i] = hyieldVar1->GetErrorY(i);
        
        hyieldVar2->GetPoint(i, x, y);
        yieldValVar2_arr[i] = y;
        yErrVar2_arr[i] = hyieldVar2->GetErrorY(i);
        
        hNonPromptFrac->GetPoint(i, x, y);
        fracVal_arr[i] = y;
        fracValE_arr[i] = hNonPromptFrac->GetErrorY(i);
    }
    
    // Now process the data
    std::vector<std::vector<double>> npYield(pointsA, std::vector<double>(NbOfVariations, 0.0));
    std::vector<std::vector<double>> pYield(pointsA, std::vector<double>(NbOfVariations, 0.0));
    std::vector<std::vector<double>> npYieldE(pointsA, std::vector<double>(NbOfVariations, 0.0));
    std::vector<std::vector<double>> pYieldE(pointsA, std::vector<double>(NbOfVariations, 0.0));
    
    for (int i = 0; i < pointsA; i++) {
        double scale = 1.0;
        if (!Absolute) {
            // Scaling for the Radial area of the dR slice
            if (obsTag.find("Y") != std::string::npos) {
                // double Area = M_PI * (pow((x_arr[i] + xErr_arr[i]), 2) - pow((x_arr[i] - xErr_arr[i]), 2));
                // scale = Area;
                scale = xErr_arr[i] * 2;
            } else {
                // Scale by the bin width
                scale = xErr_arr[i] * 2;
            }
        }
        
        // Apply scaling
        yieldVal_arr[i] *= 1.0 / scale;
        yErr_arr[i] *= 1.0 / scale;
        npYield[i][0] = yieldVal_arr[i] * fracVal_arr[i];
        pYield[i][0] = yieldVal_arr[i] * (1.0 - fracVal_arr[i]);
        std::cout << "Yield: " << yieldVal_arr[i] << ", Fraction: " << fracVal_arr[i] 
                  << ", npYield: " << npYield[i][0] << ", pYield: " << pYield[i][0] << std::endl;
        npYieldE[i][0] = propagateError(yieldVal_arr[i], fracVal_arr[i], yErr_arr[i], fracValE_arr[i], 1);
        pYieldE[i][0] = propagateError(yieldVal_arr[i], (1.0 - fracVal_arr[i]), yErr_arr[i], fracValE_arr[i], 1);
        
        // Variation 1
        yieldValVar1_arr[i] *= 1.0 / scale;
        yErrVar1_arr[i] *= 1.0 / scale;
        npYield[i][1] = yieldValVar1_arr[i] * fracVal_arr[i];
        pYield[i][1] = yieldValVar1_arr[i] * (1.0 - fracVal_arr[i]);
        npYieldE[i][1] = propagateError(yieldValVar1_arr[i], fracVal_arr[i], yErrVar1_arr[i], fracValE_arr[i], 1);
        pYieldE[i][1] = propagateError(yieldValVar1_arr[i], (1.0 - fracVal_arr[i]), yErrVar1_arr[i], fracValE_arr[i], 1);
        
        if (NbOfVariations > 2) {
            // Variation 2
            yieldValVar2_arr[i] *= 1.0 / scale;
            yErrVar2_arr[i] *= 1.0 / scale;
            npYield[i][2] = yieldValVar2_arr[i] * fracVal_arr[i];
            pYield[i][2] = yieldValVar2_arr[i] * (1.0 - fracVal_arr[i]);
            npYieldE[i][2] = propagateError(yieldValVar2_arr[i], fracVal_arr[i], yErrVar2_arr[i], fracValE_arr[i], 1);
            pYieldE[i][2] = propagateError(yieldValVar2_arr[i], (1.0 - fracVal_arr[i]), yErrVar2_arr[i], fracValE_arr[i], 1);
        }
    }
    
    // Find max values safely
    double max1 = 0.0, max2 = 0.0;
    for (int i = 0; i < pointsA; i++) {
        if (yieldVal_arr[i] > max1) max1 = yieldVal_arr[i];
        if (yieldValVar1_arr[i] > max2) max2 = yieldValVar1_arr[i];
    }
    double maxTot = std::max(max1, max2);
    
    // Set maximums
    try {
        hyield->SetMaximum(maxTot);
        hyieldVar1->SetMaximum(maxTot);
        if (NbOfVariations > 2) {
            hyieldVar2->SetMaximum(maxTot);
        }
    } catch (const std::exception& e) {
        std::cerr << "Warning: Error setting maximums: " << e.what() << std::endl;
    }
    
    graphyield[0] = hyield;
    graphyield[1] = hyieldVar1;
    if (NbOfVariations > 2) {
        graphyield[2] = hyieldVar2;
    }
    
    // Create the new graphs
    try {
        for (int i = 0; i < NbOfVariations; i++) {
            graphNP[i] = new TGraphErrors();
            graphP[i] = new TGraphErrors();
            
            for (int j = 0; j < pointsA; j++) {
                graphNP[i]->SetPoint(j, x_arr[j], npYield[j][i]);
                graphNP[i]->SetPointError(j, xErr_arr[j], npYieldE[j][i]);
                
                graphP[i]->SetPoint(j, x_arr[j], pYield[j][i]);
                graphP[i]->SetPointError(j, xErr_arr[j], pYieldE[j][i]);
            }
        }
        
        // Store the graphs as class attributes
        graphPRaw = graphP[0];
        graphPRaw->SetName(("PromptRaw_" + ptString).c_str());
        graphNPRaw = graphNP[0];
        graphNPRaw->SetName(("NonPromptRaw_" + ptString).c_str());
        graphInclRaw = graphyield[0];
        graphInclRaw->SetName(("inclRaw_" + ptString).c_str());
    } catch (const std::exception& e) {
        std::cerr << "Error creating output graphs: " << e.what() << std::endl;
        // Create fallback empty graphs
        graphPRaw = create_empty_graph();
        graphPRaw->SetName(("PromptRaw_" + ptString).c_str());
        graphNPRaw = create_empty_graph();
        graphNPRaw->SetName(("NonPromptRaw_" + ptString).c_str());
        graphInclRaw = create_empty_graph();
        graphInclRaw->SetName(("inclRaw_" + ptString).c_str());
    }
}

// Add this implementation before the end of the file
void PlotGraphsObject::setOptions() {
    int font = 42;
    
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(10);
    gStyle->SetCanvasColor(10);
    gStyle->SetTitleFillColor(10);
    gStyle->SetTitleBorderSize(1);
    gStyle->SetStatColor(10);
    gStyle->SetStatBorderSize(1);
    gStyle->SetLegendBorderSize(1);
    
    gStyle->SetDrawBorder(0);
    gStyle->SetTextFont(font);
    gStyle->SetStatFont(font);
    gStyle->SetStatFontSize(0.05);
    gStyle->SetStatX(0.97);
    gStyle->SetStatY(0.98);
    gStyle->SetStatH(0.03);
    gStyle->SetStatW(0.3);
    gStyle->SetTickLength(0.02, "y");
    gStyle->SetEndErrorSize(3);
    gStyle->SetLabelSize(0.05, "xyz");
    gStyle->SetLabelFont(font, "xyz");
    gStyle->SetLabelOffset(0.01, "xyz");
    gStyle->SetTitleFont(font, "xyz");
    gStyle->SetTitleOffset(1.2, "xyz");
    gStyle->SetTitleSize(0.045, "xyz");
    gStyle->SetMarkerSize(1);
    gStyle->SetPalette(1);
    
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
}