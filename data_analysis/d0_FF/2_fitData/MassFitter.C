// MassFitterScript.cpp - OPTIMIZED VERSION
// D0 Meson Analysis Framework for mass and lifetime fits

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
#include "RooWorkspace.h"
#include "RooPlot.h"

// Include external fitter and plotter headers
#include "PlotHelpers.h"
#include "Fitter.C"
#include "Plotter.C"

/**
 * Main class for fitting mass spectra of D0 mesons in different kinematic regions.
 * Handles creation of datasets, fitting, and visualization of results.
 */
class FitSpectraObject {
private:
    // Core configuration parameters
    bool isMC;
    std::pair<double, double> jetPt;
    std::vector<double> zBins;
    bool zTObservable;
    TTree* inputTree;
    std::pair<double, double> sideBandLimits = std::make_pair(1.840, 1.890);
    std::string outfilePath;
    std::string fOutDataName;
    std::string fOutDataNameB;
    
    int nzTBins;
    int numFitItems = 5;
    bool correctionOnly = false;
    bool enableSPlot = true;
    
    // Helper functions to reduce code duplication
    TCanvas* createStandardCanvas(const std::string& name, const std::string& title, 
                                 int width = 800, int height = 600) {
        TCanvas* canvas = new TCanvas(name.c_str(), title.c_str(), width, height);
        setupCanvasMargins(canvas);
        return canvas;
    }
    
    void setupCanvasMargins(TCanvas* canvas) {
        canvas->SetLeftMargin(0.12);
        canvas->SetRightMargin(0.05);
        canvas->SetTopMargin(0.08);
        canvas->SetBottomMargin(0.12);
    }
    
    void setupHistogramStyle(TH1* hist, int color, int markerStyle, const std::string& title = "") {
        hist->SetMarkerColor(color);
        hist->SetLineColor(color);
        hist->SetMarkerStyle(markerStyle);
        hist->SetMarkerSize(1.2);
        if (!title.empty()) hist->SetTitle(title.c_str());
    }
    
    TLegend* createStandardLegend(double x1 = 0.6, double y1 = 0.7, double x2 = 0.9, double y2 = 0.9) {
        TLegend* legend = new TLegend(x1, y1, x2, y2);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        return legend;
    }
    
    // Template for creating parameter graphs (reduces 500+ lines to ~50)
    template<typename T>
    std::vector<TGraphErrors*> createParameterGraphsTemplate(
        const std::string& baseName,
        const std::vector<std::vector<T>>& data,
        const std::vector<double>& xPos,
        const std::vector<double>& xWidth) {
        
        std::vector<TGraphErrors*> graphs(4, nullptr);
        std::vector<std::string> suffixes = {"F", "Start", "LL", "HL"};
        
        for (int type = 0; type < 4; ++type) {
            std::vector<double> yValues(nzTBins), yErrors(nzTBins, 0.0);
            
            for (int i = 0; i < nzTBins; ++i) {
                if (i < static_cast<int>(data.size())) {
                    int dataIndex = (type == 0) ? 0 : type + 1;
                    if (dataIndex < static_cast<int>(data[i].size())) {
                        yValues[i] = data[i][dataIndex];
                        if (type == 0 && data[i].size() > 1) yErrors[i] = data[i][1];
                    }
                }
            }
            
            graphs[type] = new TGraphErrors(nzTBins, xPos.data(), yValues.data(), 
                                           xWidth.data(), yErrors.data());
            graphs[type]->SetName((baseName + suffixes[type]).c_str());
        }
        return graphs;
    }
    
    template<typename T>
    TGraphAsymmErrors* createAsymmGraphTemplate(
        const std::string& name,
        const std::vector<std::vector<T>>& data,
        const std::vector<double>& xPos,
        const std::vector<double>& xWidth) {
        
        std::vector<double> startArray(nzTBins, 0.0);
        std::vector<double> llArray(nzTBins, 0.0);
        std::vector<double> hlArray(nzTBins, 0.0);
        
        for (int i = 0; i < nzTBins; ++i) {
            if (i < static_cast<int>(data.size())) {
                if (2 < static_cast<int>(data[i].size())) startArray[i] = data[i][2];
                if (3 < static_cast<int>(data[i].size())) llArray[i] = data[i][3];
                if (4 < static_cast<int>(data[i].size())) hlArray[i] = data[i][4];
            }
        }
        
        std::vector<double> errLow(nzTBins), errHigh(nzTBins);
        for (int i = 0; i < nzTBins; ++i) {
            errLow[i] = std::max(0.0, startArray[i] - llArray[i]);
            errHigh[i] = std::max(0.0, hlArray[i] - startArray[i]);
        }
        
        TGraphAsymmErrors* graph = new TGraphAsymmErrors(
            nzTBins, xPos.data(), startArray.data(),
            xWidth.data(), xWidth.data(), errLow.data(), errHigh.data());
        graph->SetName(name.c_str());
        return graph;
    }
    
    void writeGraphsToFile(TFile* file, const std::map<std::string, std::vector<TGraphErrors*>>& graphs) {
        for (const auto& [key, graphVector] : graphs) {
            for (const auto& graph : graphVector) {
                if (graph) graph->Write();
            }
        }
    }
    
    void writeAsymmGraphsToFile(TFile* file, const std::map<std::string, TGraphAsymmErrors*>& graphs) {
        for (const auto& [key, graph] : graphs) {
            if (graph) graph->Write();
        }
    }
    
    // Helper function to create and style histogram comparison plots
    void createComparisonPlot(const std::vector<TH1*>& histograms, 
                             const std::vector<std::string>& labels,
                             const std::vector<int>& colors,
                             const std::string& canvasName,
                             const std::string& outputPath) {
        if (histograms.empty()) return;
        
        TCanvas* canvas = createStandardCanvas(canvasName, labels[0], 800, 600);
        
        // Find global y-range
        double maxY = 0;
        for (const auto& hist : histograms) {
            if (hist && hist->GetMaximum() > maxY) maxY = hist->GetMaximum();
        }
        
        // Draw histograms
        for (size_t i = 0; i < histograms.size() && i < labels.size() && i < colors.size(); ++i) {
            if (!histograms[i]) continue;
            
            setupHistogramStyle(histograms[i], colors[i], 20 + i);
            histograms[i]->GetYaxis()->SetRangeUser(0, maxY * 1.2);
            histograms[i]->Draw(i == 0 ? "PE" : "PE same");
        }
        
        // Add legend
        TLegend* legend = createStandardLegend(0.6, 0.7, 0.9, 0.9);
        for (size_t i = 0; i < histograms.size() && i < labels.size(); ++i) {
            if (histograms[i]) legend->AddEntry(histograms[i], labels[i].c_str(), "pe");
        }
        legend->Draw();
        
        canvas->SaveAs(outputPath.c_str());
        delete canvas;
        delete legend;
    }
    
    // Helper function to create efficiency correction plots
    void createEfficiencyPlot(const std::vector<double>& binCenters,
                             const std::vector<double>& binWidths,
                             const std::vector<double>& correctionValues,
                             const std::vector<double>& correctionErrors,
                             const std::string& title,
                             const std::string& outputPath,
                             int color = kBlue) {
        if (binCenters.size() != correctionValues.size()) return;
        
        TCanvas* canvas = createStandardCanvas("effCanvas", title, 800, 600);
        
        TGraphErrors* graph = new TGraphErrors(binCenters.size());
        for (size_t i = 0; i < binCenters.size(); ++i) {
            graph->SetPoint(i, binCenters[i], correctionValues[i]);
            graph->SetPointError(i, binWidths[i], correctionErrors[i]);
        }
        
        graph->SetMarkerColor(color);
        graph->SetLineColor(color);
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1.2);
        graph->SetTitle(title.c_str());
        graph->Draw("APE");
        
        // Add unity line
        if (!binCenters.empty()) {
            TLine* unityLine = new TLine(binCenters.front() - binWidths.front(), 1.0,
                                        binCenters.back() + binWidths.back(), 1.0);
            unityLine->SetLineStyle(2);
            unityLine->SetLineColor(kBlack);
            unityLine->SetLineWidth(2);
            unityLine->Draw("same");
        }
        
        canvas->SaveAs(outputPath.c_str());
        delete canvas;
        delete graph;
    }
    
    // Debug output control
    bool enableDebugOutput = true;
    
    void debugPrint(const std::string& message) {
        if (enableDebugOutput) {
            std::cout << "DEBUG: " << message << std::endl;
        }
    }
    
    // Helper function to validate array sizes and provide error info
    bool validateArraySizes(const std::string& context, int expectedSize) {
        std::vector<std::pair<std::string, size_t>> arrays = {
            {"FitIPRes_PromptYield", FitIPRes_PromptYield.size()},
            {"FitIPRes_NonPromptYield", FitIPRes_NonPromptYield.size()},
            {"FitMRes_SYield", FitMRes_SYield.size()}
        };
        
        bool allValid = true;
        for (const auto& [name, size] : arrays) {
            if (static_cast<int>(size) != expectedSize) {
                std::cerr << "ERROR in " << context << ": " << name 
                         << " size mismatch (expected " << expectedSize 
                         << ", got " << size << ")" << std::endl;
                allValid = false;
            }
        }
        return allValid;
    }
    
    // Result arrays
    std::vector<std::vector<float>> FitMRes_SYield;
    std::vector<std::vector<float>> FitMRes_BYield;
    std::vector<std::vector<float>> FitMRes_Mean;
    std::vector<std::vector<float>> FitMRes_Sig1;
    std::vector<std::vector<float>> FitMRes_deltaSig;
    std::vector<std::vector<float>> FitMRes_Sig2;
    std::vector<std::vector<float>> FitMRes_alpha;
    std::vector<std::vector<float>> FitMRes_n;
    std::vector<std::vector<float>> FitMRes_DGFrac;
    std::vector<std::vector<float>> FitMRes_pol1;
    std::vector<std::vector<float>> FitMRes_pol2;
    std::vector<std::vector<float>> FitMRes_SYieldLim;
    std::vector<std::vector<float>> FitMRes_BYieldLim;
    std::vector<std::vector<float>> FitMRes_SYieldSG;
    std::vector<std::vector<float>> FitMRes_SYieldDCB;
    std::vector<std::vector<float>> Bin_TagZMean;
    std::vector<std::vector<float>> Bin_TagZMean_weighted;
    
    // IP chi2 result arrays
    std::vector<std::vector<float>> FitIPRes_SYield;
    std::vector<std::vector<float>> FitIPRes_XpPrompt;
    std::vector<std::vector<float>> FitIPRes_SigmaPrompt;
    std::vector<std::vector<float>> FitIPRes_XiPrompt;
    std::vector<std::vector<float>> FitIPRes_XpNonprompt;
    std::vector<std::vector<float>> FitIPRes_SigmaNonprompt;
    std::vector<std::vector<float>> FitIPRes_XiNonprompt;
    std::vector<std::vector<float>> FitIPRes_PromptFrac;
    std::vector<std::vector<float>> FitIPRes_PromptYield;
    std::vector<std::vector<float>> FitIPRes_NonPromptYield;
    
    // Histogram arrays
    std::vector<TH1*> ipchi2HistoArray;
    std::vector<TH1*> massHistoArray;
    
public:
    // Constructor 
    FitSpectraObject(
                    const std::pair<double, double>& ptRange, 
                    bool isMc, 
                    const std::vector<double>& zBins, 
                    bool isZtObservable,
                    TTree* tree,
                    bool enableSPlotAnalysis = true);  
    
    // Destructor
    ~FitSpectraObject();
    
    // Methods
    void startFitting();
    void configureFilePaths();
    void initializeResultArrays();
    Fitter* createFitter();
    RooDataSet* prepareMasterDataset(Fitter* fitter);
    void processCorrectionFactors(Fitter* fitter, RooDataSet* dataMaster);
    void extractCorParam(const std::string& idString, int bin, const std::string& type, 
                        std::vector<TCanvas*>& canvas, std::vector<TCanvas*>& canvasNorm, 
                        std::vector<TCanvas*>& canvasDiv, std::vector<double>& corrValue, 
                        std::vector<double>& corrValueErr, std::vector<TLegend*>& legendList,
                        RooRealVar* massVar, RooDataSet* data, RooDataSet* data2 = nullptr, 
                        RooDataSet* data3 = nullptr);
    void createCorrectionFactorGraphs(const std::vector<double>& corrValueKaon, 
                                     const std::vector<double>& corrValueErrKaon,
                                     const std::vector<double>& corrValuePion, 
                                     const std::vector<double>& corrValueErrPion,
                                     const std::vector<double>& corrValueRecoEff, 
                                     const std::vector<double>& corrValueErrRecoEff,
                                     const std::vector<double>& corrValueAcceptance, 
                                     const std::vector<double>& corrValueErrAcceptance,
                                     const std::vector<double>& corrValueCombinedPID, 
                                     const std::vector<double>& corrValueErrCombinedPID,
                                     const std::vector<double>& corrValueCombinedAll, 
                                     const std::vector<double>& corrValueErrCombinedAll);
        void processFitsByBin(Fitter* fitter, RooDataSet* dataMaster, 
                        const std::vector<double>& binCenters, 
                        const std::vector<double>& binWidths);
    void saveResultsToFile(const std::vector<double>& binCenters, 
                          const std::vector<double>& binWidths);
        
    // Add these helper functions before the saveResultsToFile implementation
    std::vector<TGraphErrors*> createParameterGraphs(
        const std::string& key, int nBins, const std::vector<double>& xPos, 
        const std::vector<std::vector<float>>& yPosArr, const std::vector<double>& xWidth);
    
    TGraphAsymmErrors* createGraphsAsymmErr(
        const std::string& key, int nBins, const std::vector<double>& xPos,
        const std::vector<std::vector<float>>& yPosArr, const std::vector<double>& xWidth);
};


/*
 * COMPLETE sPlot WORKFLOW IMPLEMENTATION
 * =====================================
 * 
 * This implementation provides a complete multi-stage sPlot statistical separation workflow:
 * 
 * STAGE 1: Mass Fit with sPlot
 * - Perform Double Gaussian fit to invariant mass distribution
 * - Extract signal and background yields from the fit
 * - Create sPlot weights for signal/background separation
 * - Save sPlot weights to ROOT file for downstream analysis
 * 
 * STAGE 2: Apply Mass sPlot Weights to IP Chi2 Data
 * - Load sPlot weights from mass fit
 * - Create signal-enhanced dataset by applying signal weights to IP chi2 data
 * - This effectively removes background contamination statistically
 * 
 * STAGE 3: IP Chi2 Fit with Weighted Data
 * - Perform Bukin fit to IP chi2 distribution using signal-enhanced data
 * - Extract prompt and nonprompt yields from the fit
 * - The fit provides the model needed for the next sPlot stage
 * 
 * STAGE 4: Create Prompt/Nonprompt sPlot Weights
 * - Use the IP chi2 fit model and yields to create a new sPlot analysis
 * - Generate prompt and nonprompt enhanced datasets
 * - Extract event-by-event sPlot weights for prompt/nonprompt separation
 * - Save these weights to ROOT file for final analysis
 * 
 * FINAL RESULT:
 * - Each event has both mass-based sPlot weights (signal vs background)
 * - AND IP chi2-based sPlot weights (prompt vs nonprompt)
 * - This allows for clean statistical separation of:
 *   1. Signal from background (mass-based)
 *   2. Prompt from nonprompt signal (IP chi2-based)
 * 
 * The implementation follows the proper sPlot formalism using RooStats::SPlot
 * and provides statistically rigorous weights for each separation.
 */

// Implementation of the FitSpectraObject constructor
FitSpectraObject::FitSpectraObject(
    const std::pair<double, double>& ptRange, 
    bool isMc, 
    const std::vector<double>& zBins, 
    bool isZtObservable,
    TTree* tree,
    bool enableSPlotAnalysis)
    : isMC(isMc), jetPt(ptRange), zBins(zBins), zTObservable(isZtObservable), 
      inputTree(tree), nzTBins(zBins.size() - 1), enableSPlot(enableSPlotAnalysis)  // Initialize enableSPlot
{
    std::cout << "Initializing FitSpectraObject with pT range: " 
              << ptRange.first << " - " << ptRange.second << " GeV/c" << std::endl;
    std::cout << "sPlot analysis: " << (enableSPlot ? "ENABLED" : "DISABLED") << std::endl;

    // Configure file paths
    configureFilePaths();
    
    // Print bin information
    if(isZtObservable){
        std::cout << "Using zT observable for fitting" << std::endl;
        std::cout << "Number of zT bins: " << nzTBins << std::endl;
    } else {
        std::cout << "Using Y observable for fitting" << std::endl;
        std::cout << "Number of Y bins: " << nzTBins << std::endl;
    }
    
    // Print bin boundaries
    std::cout << "Bin boundaries: ";
    for (size_t i = 0; i < zBins.size(); ++i) {
        std::cout << zBins[i];
        if (i < zBins.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;
    
    // Verify the tree was passed correctly
    if (inputTree) {
        std::cout << "Input tree loaded with " << inputTree->GetEntries() << " entries" << std::endl;
    } else {
        std::cerr << "Warning: Input tree is null" << std::endl;
    }
}

// Implementation of the destructor
FitSpectraObject::~FitSpectraObject() {
    std::cout << "FitSpectraObject destroying, cleaning up histograms..." << std::endl;
    
    // Clean up mass histograms
    for (size_t i = 0; i < massHistoArray.size(); ++i) {
        if (massHistoArray[i] != nullptr) {
            delete massHistoArray[i];
            massHistoArray[i] = nullptr;
        }
    }
    massHistoArray.clear();
    
    // Clean up IP chi2 histograms
    for (size_t i = 0; i < ipchi2HistoArray.size(); ++i) {
        if (ipchi2HistoArray[i] != nullptr) {
            delete ipchi2HistoArray[i];
            ipchi2HistoArray[i] = nullptr;
        }
    }
    ipchi2HistoArray.clear();
    
    std::cout << "FitSpectraObject destroyed, histograms cleaned up" << std::endl;
}

// Implementation of startFitting method
void FitSpectraObject::startFitting() {
    std::cout << "\n===== STARTING FITTING PROCESS =====" << std::endl;
    
    // Initialize arrays for results
    std::cout << "Initializing result arrays..." << std::endl;
    initializeResultArrays();
    
    // Calculate bin centers and widths
    std::cout << "Calculating bin properties..." << std::endl;
    std::vector<double> zBinsCent(nzTBins, 0.0);
    std::vector<double> zBinsWidth(nzTBins, 0.0);
    
    for (int bin = 0; bin < nzTBins; ++bin) {
        zBinsCent[bin] = 0.5 * (zBins[bin] + zBins[bin + 1]);
        zBinsWidth[bin] = 0.5 * (zBins[bin + 1] - zBins[bin]);
    }
    
    // Create fitter object and prepare datasets
    Fitter* fitter = new Fitter(
        inputTree,
        "D0",           // hardcoded for D0
        nzTBins,
        zTObservable,
        isMC,
        outfilePath,
        true            // update
    );
    
    
    std::cout << "\nPreparing master dataset..." << std::endl;
    RooDataSet* dataMaster = nullptr;
    try {
        dataMaster = prepareMasterDataset(fitter);
        std::cout << "Master dataset created with " << dataMaster->numEntries() << " entries" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "ERROR creating master dataset: " << e.what() << std::endl;
        return;
    }
    
    // Process correction factors if requested
    std::cout << "\nProcessing correction factors..." << std::endl;
    processCorrectionFactors(fitter, dataMaster);
    
    // If running just for correction factors, return early
    if (correctionOnly) {
        std::cout << "Correction-only mode, returning early" << std::endl;
        return;
    }
    
    std::cout << "\n===== STARTING BIN-BY-BIN FITTING =====" << std::endl;
    // Process fit for each bin (now includes sPlot analysis)
    processFitsByBin(fitter, dataMaster, zBinsCent, zBinsWidth);
    
    // Clean up
    delete fitter;
    delete dataMaster;
}


void FitSpectraObject::configureFilePaths() {
    // Add MC/data and observable tags for file naming
    std::string mcTag = isMC ? "MC" : "";
    std::string obsTag = zTObservable ? "zT" : "Y";
    
    // Set up output directories
    std::string outputDir;
    if (isMC) {
        outputDir = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_MC";
    } else {
        outputDir = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA";
    }
    //create output directory if it doesn't exist
    std::filesystem::create_directories(outputDir);

    // Create path strings with proper formatting
    std::stringstream ptRangeStr;
    ptRangeStr << jetPt.first << "_" << jetPt.second;
    
    outfilePath = outputDir + "/" + ptRangeStr.str() + "/";
    fOutDataName = outputDir + "/FitParametersUnBinnedD0" + obsTag + "_" + ptRangeStr.str();
    fOutDataNameB = outputDir + "/BinnedSpectraD0" + obsTag + "_" + ptRangeStr.str();
    
    // Create output directory if it doesn't exist
    std::filesystem::path dirPath(outfilePath);
    if (!std::filesystem::exists(dirPath)) {
        std::filesystem::create_directories(dirPath);
    }
    
    // Print configuration information
    std::cout << "Output directory: " << outputDir << std::endl;
    std::cout << "Output data file: " << fOutDataName << std::endl;
}



void FitSpectraObject::initializeResultArrays() {
    // No need to set numFitItems here as it's now a class member with default value 5
    
    // Create arrays for mass fit results
    FitMRes_SYield.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_BYield.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_Mean.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_Sig1.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_deltaSig.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_Sig2.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_alpha.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_n.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_DGFrac.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_pol1.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_pol2.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_SYieldLim.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_BYieldLim.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_SYieldSG.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_SYieldDCB.resize(nzTBins, std::vector<float>(numFitItems, 0.0));

    Bin_TagZMean.resize(nzTBins, std::vector<float>(1, 0.0));
    Bin_TagZMean_weighted.resize(nzTBins, std::vector<float>(1, 0.0));
    
    // Create arrays for IP chi2 fit results
    FitIPRes_SYield.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitIPRes_XpPrompt.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitIPRes_SigmaPrompt.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitIPRes_XiPrompt.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitIPRes_XpNonprompt.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitIPRes_SigmaNonprompt.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitIPRes_XiNonprompt.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitIPRes_PromptFrac.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitIPRes_PromptYield.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitIPRes_NonPromptYield.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    
    // Create arrays for histograms
    ipchi2HistoArray.resize(nzTBins, nullptr);
    massHistoArray.resize(nzTBins, nullptr);
}

RooDataSet* FitSpectraObject::prepareMasterDataset(Fitter* fitter) {
    // Get mass range from the fitter's dictionary
    auto resSig = fitter->getMassDict("D0");
    auto massRange = resSig.massRange;
    std::pair<double, double> massCut(massRange.first, massRange.second);
    
    // Create fiducial cut string
    std::string fidCut = fitter->fiducialCutString(jetPt, massCut);
    
    // Create master dataset
    RooDataSet* dataMaster = fitter->createDataSet(
        "D0",
        "Inclusive data filtered",
        fidCut,
        true,   // isMass
        -1      // corrVer
    );
    
    return dataMaster;
}



void FitSpectraObject::processCorrectionFactors(Fitter* fitter, RooDataSet* dataMaster) {
    if (!fitter || !dataMaster) {
        std::cerr << "Error: Null fitter or dataset in processCorrectionFactors" << std::endl;
        return;
    }
    
    std::cout << "Processing correction factors for kaon and pion selection efficiency..." << std::endl;
    std::cout << "Dataset entries: " << dataMaster->numEntries() << std::endl;
    
    // Check if efficiency variables are available in the dataset
    const RooArgSet* vars = dataMaster->get(0);
    if (!vars) {
        std::cerr << "Error: Cannot access dataset variables" << std::endl;
        return;
    }
    
    RooAbsReal* kaonEff = dynamic_cast<RooAbsReal*>(vars->find("kaon_efficiency"));
    RooAbsReal* pionEff = dynamic_cast<RooAbsReal*>(vars->find("pion_efficiency"));
    RooAbsReal* recoEff = dynamic_cast<RooAbsReal*>(vars->find("reconstruction_efficiency"));
    RooAbsReal* combinedEff = dynamic_cast<RooAbsReal*>(vars->find("combined_efficiency"));
    RooAbsReal* combinedPID = dynamic_cast<RooAbsReal*>(vars->find("combined_PID_efficiency"));
    RooAbsReal* acceptance = dynamic_cast<RooAbsReal*>(vars->find("acceptance"));
    RooAbsReal* combinedEffAndAcceptance = dynamic_cast<RooAbsReal*>(vars->find("combined_eff_and_acceptance"));
    
    if (!kaonEff || !pionEff) {
        std::cerr << "Error: Required efficiency variables (kaon_efficiency, pion_efficiency) not found in dataset" << std::endl;
        std::cerr << "Available variables: ";
        vars->Print("v");
        return;
    }
    
    if (!recoEff) {
        std::cerr << "Warning: reconstruction_efficiency not found in dataset, will skip reco efficiency correction" << std::endl;
    }
    
    if (!combinedEff) {
        std::cerr << "Warning: combined_efficiency not found in dataset, will skip combined efficiency correction" << std::endl;
    }
    
    if (!combinedPID) {
        std::cerr << "Warning: combined_PID_efficiency not found in dataset, will skip combined PID efficiency correction" << std::endl;
    }
    
    if (!acceptance) {
        std::cerr << "Warning: acceptance not found in dataset, will skip acceptance correction" << std::endl;
    }
    
    if (!combinedEffAndAcceptance) {
        std::cerr << "Warning: combined_eff_and_acceptance not found in dataset, will skip combined all efficiency correction" << std::endl;
    }
    
    // Process correction factors for different versions
    std::vector<int> corrVersions = {1, 2, 3, 4, 5, 6}; // kaon, pion, reco_eff, acceptance, combined_pid, combined_all
    std::vector<std::string> corrTypes = {"kaon", "pion", "reco_eff", "acceptance", "combined_pid", "combined_all"};
    
    // Create canvases for correction parameter extraction
    std::vector<TCanvas*> canvas(6);
    std::vector<TCanvas*> canvasNorm(6);
    std::vector<TCanvas*> canvasDiv(6);
    
    // Initialize correction value arrays
    std::vector<double> corrValueKaon(nzTBins, 0.0);
    std::vector<double> corrValueErrKaon(nzTBins, 0.0);
    std::vector<double> corrValuePion(nzTBins, 0.0);
    std::vector<double> corrValueErrPion(nzTBins, 0.0);
    std::vector<double> corrValueRecoEff(nzTBins, 0.0);
    std::vector<double> corrValueErrRecoEff(nzTBins, 0.0);
    std::vector<double> corrValueAcceptance(nzTBins, 0.0);
    std::vector<double> corrValueErrAcceptance(nzTBins, 0.0);
    std::vector<double> corrValueCombinedPID(nzTBins, 0.0);
    std::vector<double> corrValueErrCombinedPID(nzTBins, 0.0);
    std::vector<double> corrValueCombinedAll(nzTBins, 0.0);
    std::vector<double> corrValueErrCombinedAll(nzTBins, 0.0);
    
    std::vector<std::vector<double>*> corrValueArrays = {&corrValueKaon, &corrValuePion, &corrValueRecoEff, 
                                                        &corrValueAcceptance, &corrValueCombinedPID, &corrValueCombinedAll};
    std::vector<std::vector<double>*> corrValueErrArrays = {&corrValueErrKaon, &corrValueErrPion, &corrValueErrRecoEff,
                                                           &corrValueErrAcceptance, &corrValueErrCombinedPID, &corrValueErrCombinedAll};
    
    // Create legends
    std::vector<TLegend*> legendList(nzTBins);
    for (int i = 0; i < nzTBins; ++i) {
        legendList[i] = new TLegend(0.65, 0.7, 0.9, 0.9);
        legendList[i]->SetTextFont(42);
        legendList[i]->SetBorderSize(0);
        legendList[i]->SetFillStyle(0);
    }
    
    // Create canvases
    for (int i = 0; i < 6; ++i) {
        std::string canvasName = "canvas_" + corrTypes[i];
        std::string canvasNormName = "canvasNorm_" + corrTypes[i];
        std::string canvasDivName = "canvasDiv_" + corrTypes[i];
        
        canvas[i] = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 800);
        canvasNorm[i] = new TCanvas(canvasNormName.c_str(), canvasNormName.c_str(), 1200, 800);
        canvasDiv[i] = new TCanvas(canvasDivName.c_str(), canvasDivName.c_str(), 1200, 800);
        
        // Divide canvases based on number of bins
        int nCols = std::min(nzTBins, 3);
        int nRows = (nzTBins + nCols - 1) / nCols;
        
        canvas[i]->Divide(nCols, nRows);
        canvasNorm[i]->Divide(nCols, nRows);
        canvasDiv[i]->Divide(nCols, nRows);
    }
    
    // Get mass variable for plotting
    auto resSig = fitter->getMassDict("D0");
    auto massRange = resSig.massRange;
    RooRealVar* massVar = new RooRealVar("tagMass", "tagMass", massRange.first, massRange.second);
    
    for (int corrVer : corrVersions) {
        std::cout << "\nProcessing correction version " << corrVer << std::endl;
        
        // Create datasets with weights for each correction version
        std::string datasetName = "data_corrected_v" + std::to_string(corrVer);
        
        // Create fiducial cut string
        std::pair<double, double> massCut(massRange.first, massRange.second);
        std::string fidCut = fitter->fiducialCutString(jetPt, massCut);
        
        // Create weighted dataset for this correction version
        RooDataSet* correctedData = fitter->createDataSet(
            "D0",
            datasetName,
            fidCut,
            true,   // isMass
            corrVer // correction version
        );
        
        if (correctedData) {
            std::cout << "Created corrected dataset v" << corrVer 
                      << " with " << correctedData->numEntries() << " entries" << std::endl;
            
            // Debug: Check the efficiency weights in the corrected dataset
            if (correctedData->numEntries() > 0) {
                std::cout << "Debug: Checking efficiency weights for corrVer " << corrVer << std::endl;
                double totalWeight = 0.0;
                double minWeight = 1e10;
                double maxWeight = 0.0;
                int validWeights = 0;
                
                // Check the efficiency values in the dataset
                for (int i = 0; i < std::min(20, correctedData->numEntries()); ++i) {
                    const RooArgSet* row = correctedData->get(i);
                    double weight = correctedData->weight();
                    
                    // Also check the actual efficiency values
                    RooAbsReal* kaonEff = dynamic_cast<RooAbsReal*>(row->find("kaon_efficiency"));
                    RooAbsReal* pionEff = dynamic_cast<RooAbsReal*>(row->find("pion_efficiency"));
                    RooAbsReal* combinedPID = dynamic_cast<RooAbsReal*>(row->find("combined_PID_efficiency"));
                    
                    if (i < 5) {  // Print first 5 entries
                        std::cout << "  Entry " << i << ": weight=" << weight;
                        if (kaonEff) std::cout << ", kaon_eff=" << kaonEff->getVal();
                        if (pionEff) std::cout << ", pion_eff=" << pionEff->getVal();
                        if (combinedPID) std::cout << ", combined_PID_eff=" << combinedPID->getVal();
                        std::cout << std::endl;
                    }
                    
                    if (weight > 0) {
                        totalWeight += weight;
                        minWeight = std::min(minWeight, weight);
                        maxWeight = std::max(maxWeight, weight);
                        validWeights++;
                    }
                }
                
                if (validWeights > 0) {
                    std::cout << "  Weight statistics (first 20 entries):" << std::endl;
                    std::cout << "    Average: " << totalWeight / validWeights << std::endl;
                    std::cout << "    Min: " << minWeight << std::endl;
                    std::cout << "    Max: " << maxWeight << std::endl;
                    std::cout << "    Valid weights: " << validWeights << std::endl;
                }
            }
            
            // Process each bin for correction parameter extraction
            for (int iBin = 0; iBin < nzTBins; ++iBin) {
                // Create z bin selection string
                std::string zCut;
                if (zTObservable) {
                    zCut = "tagZ > " + std::to_string(zBins[iBin]) + 
                           " && tagZ < " + std::to_string(zBins[iBin + 1]);
                } else {
                    zCut = "tagY > " + std::to_string(zBins[iBin]) + 
                           " && tagY < " + std::to_string(zBins[iBin + 1]);
                }
                
                // Create bin datasets
                RooDataSet* dataBin = static_cast<RooDataSet*>(dataMaster->reduce(zCut.c_str()));
                RooDataSet* correctedBin = static_cast<RooDataSet*>(correctedData->reduce(zCut.c_str()));
                
                if (dataBin && correctedBin) {
                    // Extract correction parameters
                    extractCorParam(
                        "v" + std::to_string(corrVer),
                        iBin,
                        corrTypes[corrVer - 1],
                        canvas,
                        canvasNorm,
                        canvasDiv,
                        *corrValueArrays[corrVer - 1],
                        *corrValueErrArrays[corrVer - 1],
                        legendList,
                        massVar,
                        dataBin,
                        correctedBin
                    );
                }
                
                delete dataBin;
                delete correctedBin;
            }
            
            // Calculate efficiency statistics
            double totalWeight = 0.0;
            double totalEntries = 0.0;
            double minEff = 1.0;
            double maxEff = 0.0;
            
            for (int i = 0; i < correctedData->numEntries(); ++i) {
                const RooArgSet* row = correctedData->get(i);
                double weight = correctedData->weight();
                
                totalWeight += weight;
                totalEntries += 1.0;
                
                if (weight > 0) {
                    minEff = std::min(minEff, weight);
                    maxEff = std::max(maxEff, weight);
                }
            }
            
            double avgEff = totalWeight / totalEntries;
            
            std::cout << "  Average efficiency: " << avgEff << std::endl;
            std::cout << "  Min efficiency: " << minEff << std::endl;
            std::cout << "  Max efficiency: " << maxEff << std::endl;
            std::cout << "  Total weighted entries: " << totalWeight << std::endl;
            // Clean up
            delete correctedData;
        } else {
            std::cerr << "Error: Failed to create corrected dataset for version " << corrVer << std::endl;
        }
    }
    
    // Create and save correction factor graphs
    createCorrectionFactorGraphs(corrValueKaon, corrValueErrKaon, 
                                 corrValuePion, corrValueErrPion,
                                 corrValueRecoEff, corrValueErrRecoEff,
                                 corrValueAcceptance, corrValueErrAcceptance,
                                 corrValueCombinedPID, corrValueErrCombinedPID,
                                 corrValueCombinedAll, corrValueErrCombinedAll);
    
    // Clean up
    delete massVar;
    for (auto& canvas_ptr : canvas) delete canvas_ptr;
    for (auto& canvas_ptr : canvasNorm) delete canvas_ptr;
    for (auto& canvas_ptr : canvasDiv) delete canvas_ptr;
    for (auto& legend_ptr : legendList) delete legend_ptr;
    
    // Calculate combined efficiency correction factor
    std::cout << "\nCalculating combined efficiency correction factors..." << std::endl;
    
    double totalCombinedEff = 0.0;
    int validEntries = 0;
    
    for (int i = 0; i < dataMaster->numEntries(); ++i) {
        const RooArgSet* row = dataMaster->get(i);
        
        RooAbsReal* kaonEffVar = dynamic_cast<RooAbsReal*>(row->find("kaon_efficiency"));
        RooAbsReal* pionEffVar = dynamic_cast<RooAbsReal*>(row->find("pion_efficiency"));
        
        if (kaonEffVar && pionEffVar) {
            double kaonEffVal = kaonEffVar->getVal();
            double pionEffVal = pionEffVar->getVal();
            double combinedEff = kaonEffVal * pionEffVal;
            
            if (combinedEff > 0) {
                totalCombinedEff += combinedEff;
                validEntries++;
            }
        }
    }
    
    if (validEntries > 0) {
        double avgCombinedEff = totalCombinedEff / validEntries;
        std::cout << "Average combined (kaon x pion) efficiency: " << avgCombinedEff << std::endl;
        std::cout << "Valid entries for efficiency calculation: " << validEntries << std::endl;
    } else {
        std::cerr << "Warning: No valid efficiency entries found" << std::endl;
    }
    
    std::cout << "Correction factors processing completed." << std::endl;
}

void FitSpectraObject::extractCorParam(const std::string& idString, int bin, const std::string& type,
                                      std::vector<TCanvas*>& canvas, std::vector<TCanvas*>& canvasNorm,
                                      std::vector<TCanvas*>& canvasDiv, std::vector<double>& corrValue,
                                      std::vector<double>& corrValueErr, std::vector<TLegend*>& legendList,
                                      RooRealVar* massVar, RooDataSet* data, RooDataSet* data2, 
                                      RooDataSet* data3) {
    
    std::cout << "Extracting correction parameters for " << idString << ", bin " << bin << ", type " << type << std::endl;
    
    // Validate inputs
    if (!data || !massVar || bin < 0 || bin >= static_cast<int>(corrValue.size())) {
        std::cerr << "Error: Invalid inputs to extractCorParam" << std::endl;
        return;
    }
    
    // Set ROOT style
    // gROOT->SetStyle("Plain");
    TGaxis::SetMaxDigits(2);
    gStyle->SetOptStat(0);
    
    // Create histograms from datasets with proper binning
    int nBins = 40;  // Match Python version
    std::string histName = "h_MassSignal" + std::to_string(bin) + "_data";
    TH1* h_data = data->createHistogram(histName.c_str(), *massVar, RooFit::Binning(nBins));
    
    if (!h_data) {
        std::cerr << "Error: Failed to create histogram for data" << std::endl;
        return;
    }
    
    // Debug: Check the first dataset
    // std::cout << "Debug: h_data (unweighted) - entries: " << h_data->GetEntries() 
    //           << ", integral: " << h_data->Integral() << std::endl;
    
    std::vector<double> maxList = {h_data->GetMaximum(), 0.0};
    
    TH1* h_data2 = nullptr;
    if (data2) {
        std::string histName2 = "h_MassSignal" + std::to_string(bin) + "_data2";
        h_data2 = data2->createHistogram(histName2.c_str(), *massVar, RooFit::Binning(nBins));
        if (h_data2) {
            maxList.push_back(h_data2->GetMaximum());
            
            // // Debug: Check the second dataset
            // std::cout << "Debug: h_data2 (weighted) - entries: " << h_data2->GetEntries() 
            //           << ", integral: " << h_data2->Integral() << std::endl;
            
            // // Debug: Check if the datasets are actually different
            // std::cout << "Debug: Dataset comparison:" << std::endl;
            // std::cout << "  data (unweighted): " << data->numEntries() << " entries" << std::endl;
            // std::cout << "  data2 (weighted): " << data2->numEntries() << " entries" << std::endl;
            
            // Check sum of weights
            double sumWeights = 0;
            for (int i = 0; i < data2->numEntries(); ++i) {
                data2->get(i);
                sumWeights += data2->weight();
            }
            // std::cout << "  data2 sum of weights: " << sumWeights << std::endl;
            // std::cout << "  data entries (should be same): " << data->numEntries() << std::endl;
        }
    }
    
    TH1* h_data3 = nullptr;
    if (data3) {
        std::string histName3 = "h_MassSignal" + std::to_string(bin) + "_data3";
        h_data3 = data3->createHistogram(histName3.c_str(), *massVar, RooFit::Binning(nBins));
        if (h_data3) {
            maxList.push_back(h_data3->GetMaximum());
        }
    }
    
    // Find the type index for canvas access
    int typeIndex = 0;
    if (type == "kaon") typeIndex = 0;
    else if (type == "pion") typeIndex = 1;
    else if (type == "reco_eff") typeIndex = 2;
    else if (type == "acceptance") typeIndex = 3;
    else if (type == "combined_pid") typeIndex = 4;
    else if (type == "combined_all") typeIndex = 5;
    else {
        std::cerr << "Warning: Unknown correction type: " << type << std::endl;
        typeIndex = 0;  // Default to first type
    }
    
    // Validate canvas arrays
    if (typeIndex >= static_cast<int>(canvas.size()) || typeIndex >= static_cast<int>(canvasNorm.size()) || typeIndex >= static_cast<int>(canvasDiv.size())) {
        std::cerr << "Error: Invalid canvas array access for type " << type << std::endl;
        delete h_data;
        if (h_data2) delete h_data2;
        if (h_data3) delete h_data3;
        return;
    }
    
    // Normalized canvas plotting
    if (canvasNorm[typeIndex]) {
        canvasNorm[typeIndex]->cd(bin + 1);
        TPad* myPad2 = new TPad(("myPad2_" + std::to_string(bin)).c_str(), 
                               ("The pad " + std::to_string(bin)).c_str(), 0, 0, 1, 1);
        myPad2->SetLeftMargin(0.15);  // Adjusted to match Python
        myPad2->SetTopMargin(0.08);
        myPad2->SetRightMargin(0.1);
        myPad2->SetBottomMargin(0.15);
        myPad2->Draw();
        
        // Set histogram properties and draw normalized
        h_data->SetMarkerStyle(20);
        h_data->SetMarkerSize(0.5);
        h_data->SetMarkerColor(1);
        h_data->SetLineColor(kGray + 2);
        h_data->SetTitle(Form("Signal for %s", idString.c_str()));
        h_data->DrawNormalized("E");
        
        if (h_data2) {
            h_data2->SetMarkerStyle(20);
            h_data2->SetMarkerSize(0.5);
            h_data2->SetMarkerColor(kRed + 2);
            h_data2->SetLineColor(kRed + 2);
            h_data2->DrawNormalized("same E");
        }
        
        if (h_data3) {
            h_data3->SetMarkerStyle(20);
            h_data3->SetMarkerSize(0.5);
            h_data3->SetMarkerColor(kBlue + 2);
            h_data3->SetLineColor(kBlue + 2);
            h_data3->DrawNormalized("same E");
        }
        
        canvasNorm[typeIndex]->cd();
    }
    
    // Regular canvas plotting
    if (canvas[typeIndex]) {
        canvas[typeIndex]->cd(bin + 1);
        TPad* myPad = new TPad(("myPad_" + std::to_string(bin)).c_str(), 
                              ("The pad " + std::to_string(bin)).c_str(), 0, 0, 1, 1);
        myPad->SetLeftMargin(0.15);  // Adjusted to match Python
        myPad->SetTopMargin(0.08);
        myPad->SetRightMargin(0.1);
        myPad->SetBottomMargin(0.15);
        myPad->Draw();
        
        double maxVal = *std::max_element(maxList.begin(), maxList.end());
        h_data->GetYaxis()->SetRangeUser(0, maxVal * 1.3);
        h_data->SetMarkerStyle(20);
        h_data->SetMarkerSize(0.5);
        h_data->SetMarkerColor(1);
        h_data->SetLineColor(kGray + 2);
        h_data->DrawCopy("E");
        
        if (h_data2) {
            h_data2->SetMarkerStyle(20);
            h_data2->SetMarkerSize(0.5);
            h_data2->SetMarkerColor(kRed + 2);
            h_data2->SetLineColor(kRed + 2);
            h_data2->DrawCopy("same E");
        }
        
        if (h_data3) {
            h_data3->SetMarkerStyle(20);
            h_data3->SetMarkerSize(0.5);
            h_data3->SetMarkerColor(kBlue + 2);
            h_data3->SetLineColor(kBlue + 2);
            h_data3->DrawCopy("same E");
        }
        
        canvas[typeIndex]->cd();
    }
    
    h_data->GetYaxis()->UnZoom();
    
    // Initialize correction values
    double fitValue = 0.0;
    double fitValueErr = 0.0;
    
    // Ratio canvas for correction factor extraction
    TGaxis::SetMaxDigits(4);
    if (canvasDiv[typeIndex] && h_data2) {
        canvasDiv[typeIndex]->cd(bin + 1);
        
        TPad* myPad3 = new TPad(("myPad3_" + std::to_string(bin)).c_str(), 
                               ("The pad3 " + std::to_string(bin)).c_str(), 0, 0, 1, 1);
        myPad3->SetLeftMargin(0.15);  // Adjusted to match Python
        myPad3->SetTopMargin(0.08);
        myPad3->SetRightMargin(0.1);
        myPad3->SetBottomMargin(0.15);
        myPad3->Draw();
        
        // Create ratio histogram
        TH1* h_ratio = static_cast<TH1*>(h_data2->Clone(("h_ratio_" + std::to_string(bin)).c_str()));
        if (h_ratio) {
            h_ratio->Divide(h_data);
            h_ratio->SetMarkerStyle(20);
            h_ratio->SetMarkerSize(0.5);
            h_ratio->SetMarkerColor(1);
            h_ratio->SetLineColor(kGray + 2);
            h_ratio->DrawCopy("E");
            
            // // Debug: Check the ratio histogram
            // std::cout << "Debug: Ratio histogram for bin " << bin << ", type " << type << std::endl;
            // std::cout << "  h_data entries: " << h_data->GetEntries() << ", integral: " << h_data->Integral() << std::endl;
            // std::cout << "  h_data2 entries: " << h_data2->GetEntries() << ", integral: " << h_data2->Integral() << std::endl;
            // std::cout << "  h_ratio entries: " << h_ratio->GetEntries() << ", integral: " << h_ratio->Integral() << std::endl;
            
            // // Check a few bins of the ratio
            // for (int i = 1; i <= std::min(5, h_ratio->GetNbinsX()); ++i) {
            //     std::cout << "    Bin " << i << ": h_data=" << h_data->GetBinContent(i) 
            //               << ", h_data2=" << h_data2->GetBinContent(i) 
            //               << ", ratio=" << h_ratio->GetBinContent(i) << std::endl;
            // }
            
            // Fit a horizontal line to extract correction factor
            // Use proper bin center calculation matching Python
            double minL = h_data->GetBinCenter(1);
            double maxL = h_data->GetBinCenter(h_data->GetNbinsX());
            TF1* line = new TF1("Line", "[0]", minL, maxL);
            line->SetParameter(0, 1.0);
            
            // Perform fit with proper options
            int fitStatus = h_ratio->Fit("Line", "NQS", "", minL, maxL);
            if (fitStatus == 0) {  // Successful fit
                fitValue = line->GetParameter(0);
                fitValueErr = line->GetParError(0);
                
                line->SetLineColor(kRed);
                line->SetLineWidth(2);
                line->DrawCopy("same");
                
                // Add correction factor to legend if available
                if (bin < static_cast<int>(legendList.size()) && legendList[bin]) {
                    legendList[bin]->AddEntry(h_ratio, 
                                            Form("Correction Constant %.3f±%.3f", fitValue, fitValueErr), 
                                            "p");
                    legendList[bin]->Draw();
                }
            } else {
                std::cerr << "Warning: Fit failed for bin " << bin << ", type " << type << std::endl;
            }
            
            delete line;
            delete h_ratio;
        }
        
        canvasDiv[typeIndex]->cd();
    }
    
    // Store correction values with bounds checking
    if (bin < static_cast<int>(corrValue.size())) {
        corrValue[bin] = fitValue;
    }
    if (bin < static_cast<int>(corrValueErr.size())) {
        corrValueErr[bin] = fitValueErr;
    }
    
    // Save canvases if this is the last bin
    if (bin == (static_cast<int>(corrValue.size()) - 1)) {
        std::string canvasName = outfilePath + "CorrFac_" + idString + "_" + type;
        
        if (canvas[typeIndex]) {
            canvas[typeIndex]->SaveAs((canvasName + "_zT.png").c_str());
        }
        if (canvasDiv[typeIndex]) {
            canvasDiv[typeIndex]->SaveAs((canvasName + "_zTRatio.png").c_str());
        }
        if (canvasNorm[typeIndex]) {
            canvasNorm[typeIndex]->SaveAs((canvasName + "_zTNorm.png").c_str());
        }
        
        std::cout << "Saved correction parameter canvases for " << idString << ", type " << type << std::endl;
    }
    
    // Clean up histograms
    delete h_data;
    if (h_data2) delete h_data2;
    if (h_data3) delete h_data3;
    
    std::cout << "Correction parameter extraction completed for bin " << bin 
              << " with correction factor: " << fitValue << " ± " << fitValueErr << std::endl;
}

void FitSpectraObject::processFitsByBin(Fitter* fitter, RooDataSet* dataMaster, 
                                      const std::vector<double>& binCenters, 
                                      const std::vector<double>& binWidths) {
    if (!fitter || !dataMaster) {
        std::cerr << "Error: Null fitter or dataset in processFitsByBin" << std::endl;
        return;
    }
    
    // Create sPlot output file if enabled
    TFile* splotFile = nullptr;
    std::string splotFileName = outfilePath + "splot_results.root";  // Declare in broader scope
    
    if (enableSPlot) {
        splotFile = new TFile(splotFileName.c_str(), "RECREATE");
        
        if (!splotFile || splotFile->IsZombie()) {
            std::cerr << "Error: Cannot create sPlot output file: " << splotFileName << std::endl;
            if (splotFile) {
                delete splotFile;
                splotFile = nullptr;
            }
            enableSPlot = false;  // Disable sPlot if file creation fails
        } else {
            std::cout << "Created sPlot output file: " << splotFileName << std::endl;
        }
    }
    
    // Get resonance parameters from fitter
    auto resSig = fitter->getMassDict("D0");
    auto sigRegion = resSig.signalRegion;
    
    // Loop over all bins
    for (int iBin = 0; iBin < nzTBins; ++iBin) {
        
        // Create z bin selection string
        std::string zCut;
        if (zTObservable) {
            std::cout << "\n==== Processing bin " << iBin 
                      << " (z: " << zBins[iBin] << "-" << zBins[iBin+1] << ") ====" << std::endl;
            zCut = "tagZ >= " + std::to_string(zBins[iBin]) + 
                   " && tagZ < " + std::to_string(zBins[iBin+1]);
        } else {
            std::cout << "\n==== Processing bin " << iBin 
                      << " (Y: " << zBins[iBin] << "-" << zBins[iBin+1] << ") ====" << std::endl;
            zCut = "tagY >= " + std::to_string(zBins[iBin]) + 
                   " && tagY < " + std::to_string(zBins[iBin+1]);
        }
        
        // Create bin dataset
        RooDataSet* dataBin = static_cast<RooDataSet*>(dataMaster->reduce(zCut.c_str()));
        std::cout << "Bin dataset created with " << dataBin->numEntries() << " entries" << std::endl;
        
        // Skip if too few entries
        if (dataBin->numEntries() < 10) {
            std::cout << "Too few entries in bin, skipping..." << std::endl;
            delete dataBin;
            continue;
        }
        
        // Plot tagZ distribution and calculate mean
        std::cout << "Plotting tagZ distribution for bin " << iBin << "..." << std::endl;
        TH1D* tagZHist = new TH1D(("tagZHist_bin" + std::to_string(iBin)).c_str(), 
                                  ("TagZ Distribution for Bin " + std::to_string(iBin)).c_str(), 
                                  20, 0, 1);
        
        // Fill histogram with tagZ values
        for (int i = 0; i < dataBin->numEntries(); ++i) {
            const RooArgSet* row = dataBin->get(i);
            RooRealVar* tagZVar = dynamic_cast<RooRealVar*>(row->find("tagZ"));
            if (tagZVar) {
                tagZHist->Fill(tagZVar->getVal());
            }
        }
        
                            TH1D* tagZHistSafe = (TH1D*)tagZHist->Clone(("tagZHistSafe_bin" + std::to_string(iBin)).c_str());
                            tagZHistSafe->SetDirectory(nullptr);  // Prevent ROOT from deleting it
                            
        // Calculate mean and standard deviation
        double meanTagZ = tagZHist->GetMean();
        double stdDevTagZ = tagZHist->GetStdDev();
        std::cout << "Mean tagZ: " << meanTagZ << ", StdDev: " << stdDevTagZ << std::endl;
        Bin_TagZMean[iBin][0] = meanTagZ;  // Store mean tagZ for this bin
        // Save histogram to file
        TCanvas* canvas = new TCanvas(("canvas_tagZ_bin" + std::to_string(iBin)).c_str(), 
                                      ("TagZ Distribution for Bin " + std::to_string(iBin)).c_str(), 
                                      800, 600);
                                      // set plotting style for tagZ distribution
        tagZHist->SetLineColor(kBlue);
        tagZHist->SetMarkerColor(kBlue);
        tagZHist->SetLineWidth(2);
        tagZHist->SetMarkerStyle(20);
        tagZHist->SetTitle(("TagZ Distribution for Bin " + std::to_string(iBin)).c_str());
        tagZHist->GetXaxis()->SetTitle("#it{z}_{T}");
        tagZHist->GetYaxis()->SetTitle("Entries");
        
        tagZHist->Draw("pe");
        std::string outputFile = outfilePath + "tagZDistribution_bin" + std::to_string(iBin) + ".png";
        canvas->SaveAs(outputFile.c_str());
        delete canvas;
        
        // Perform mass fit
        std::string zRangeStr = std::to_string(zBins[iBin]) + "-" + std::to_string(zBins[iBin+1]);
        std::cout << "Performing mass fit..." << std::endl;
        
        TH1* massFitHisto = nullptr;
        std::vector<double> fitParams;
        std::vector<double> fitErrors;
        
        // Declare all variables that will be used in debug output outside try block
        RooDataSet* signalEnhancedData = nullptr;
        RooDataSet* bkgLeft = nullptr;
        RooDataSet* bkgRight = nullptr;
        RooDataSet* bkgData = nullptr;
        TH1* ipChi2Histo = nullptr;
        std::vector<double> ipParams, ipErrors;
        RooAbsPdf* ipChi2Model = nullptr;
        RooRealVar* sig_yieldLim = nullptr;
        RooRealVar* prompt_frac = nullptr;
        
        try {
            // Perform Single Gaussian fit
            std::cout << "  Performing Single Gaussian fit..." << std::endl;
            std::tie(massFitHisto, fitParams, fitErrors) = 
                fitter->massFit("D0", dataBin, "SGauss", iBin, zRangeStr);
            
            // Store SG yield for comparison
            FitMRes_SYieldSG[iBin][0] = fitParams[0];
            FitMRes_SYieldSG[iBin][1] = fitErrors[0];
            
            // Perform Double Gaussian fit - THIS IS NOW OUR PRIMARY FIT
            std::cout << "  Performing Double Gaussian fit..." << std::endl;
            std::vector<double> dgParams, dgErrors;
            std::tie(massFitHisto, dgParams, dgErrors) = 
                fitter->massFit("D0", dataBin, "DGauss", iBin, zRangeStr, enableSPlot, splotFile);
            
            // Store histogram safely
            if (massFitHisto) {
                // Create a safe copy to avoid potential pointer invalidation issues
                TH1* massFitHistoCopy = (TH1*)massFitHisto->Clone(("massFitHisto_bin" + std::to_string(iBin)).c_str());
                if (massFitHistoCopy) {
                    massFitHistoCopy->SetDirectory(nullptr);  // Prevent ROOT from managing it
                    massHistoArray[iBin] = massFitHistoCopy;
                    std::cout << "  Stored safe copy of mass fit histogram for bin " << iBin << std::endl;
                } else {
                    std::cout << "  Warning: Failed to create safe copy of mass fit histogram for bin " << iBin << std::endl;
                    massHistoArray[iBin] = nullptr;
                }
            } else {
                massHistoArray[iBin] = nullptr;
            }

            // Ensure sPlot file is properly flushed and synchronized before reading
            if (enableSPlot && splotFile) {
                splotFile->Flush();
                // Force write all pending data to disk
                splotFile->Write("", TObject::kOverwrite);
                std::cout << "  sPlot file flushed and synchronized" << std::endl;
            }
            
            
            // Store fit results in arrays
            if (dgParams.size() >= 12 && dgErrors.size() >= 10) {
                // Store values in result arrays from DGauss fit
                FitMRes_SYield[iBin][0] = dgParams[0];  // Signal yield
                FitMRes_SYield[iBin][1] = dgErrors[0];  // Signal yield error
                
                FitMRes_BYield[iBin][0] = dgParams[1];  // Background yield
                FitMRes_BYield[iBin][1] = dgErrors[1];  // Background yield error
                
                FitMRes_Mean[iBin][0] = dgParams[2];    // Mean mass
                FitMRes_Mean[iBin][1] = dgErrors[2];    // Mean mass error
                
                FitMRes_Sig1[iBin][0] = dgParams[3];    // Sigma 1
                FitMRes_Sig1[iBin][1] = dgErrors[3];    // Sigma 1 error
                
                FitMRes_deltaSig[iBin][0] = dgParams[4]; // Delta sigma
                FitMRes_deltaSig[iBin][1] = dgErrors[4]; // Delta sigma error
                
                FitMRes_Sig2[iBin][0] = dgParams[3] * dgParams[4]; // Sigma 2 = Sigma1 * deltaSig
                FitMRes_Sig2[iBin][1] = 0.0;  // We don't have direct error on Sigma2
                
                // These parameters don't apply to DGauss but we'll set them to 0 for consistency
                FitMRes_alpha[iBin][0] = 0.0;   // Not used in DGauss
                FitMRes_alpha[iBin][1] = 0.0;   // Not used in DGauss
                
                FitMRes_n[iBin][0] = 0.0;       // Not used in DGauss
                FitMRes_n[iBin][1] = 0.0;       // Not used in DGauss
                
                FitMRes_DGFrac[iBin][0] = dgParams[7];  // Gaussian1 fraction (similar to CB fraction)
                FitMRes_DGFrac[iBin][1] = dgErrors[7];  // Gaussian1 fraction error
                
                FitMRes_pol1[iBin][0] = dgParams[8];    // Pol1
                FitMRes_pol1[iBin][1] = dgErrors[8];    // Pol1 error
                
                FitMRes_pol2[iBin][0] = dgParams[9];    // Pol2
                FitMRes_pol2[iBin][1] = dgErrors[9];    // Pol2 error
                
                FitMRes_SYieldLim[iBin][0] = dgParams[10]; // Signal yield in limit region
                FitMRes_BYieldLim[iBin][0] = dgParams[11]; // Background yield in limit region
                
                // Print fit results
                std::cout << "Double Gaussian fit results for bin " << iBin << ":" << std::endl;
                std::cout << "  Signal yield: " << FitMRes_SYield[iBin][0] 
                          << " ± " << FitMRes_SYield[iBin][1] << std::endl;
                std::cout << "  Background yield: " << FitMRes_BYield[iBin][0] 
                          << " ± " << FitMRes_BYield[iBin][1] << std::endl;
                std::cout << "  Mean: " << FitMRes_Mean[iBin][0] 
                          << " ± " << FitMRes_Mean[iBin][1] << std::endl;
                std::cout << "  Sigma1: " << FitMRes_Sig1[iBin][0] 
                          << " ± " << FitMRes_Sig1[iBin][1] << std::endl;
                std::cout << "  Sigma2: " << FitMRes_Sig2[iBin][0] << std::endl;
                std::cout << "  Gaussian1 fraction: " << FitMRes_DGFrac[iBin][0] 
                          << " ± " << FitMRes_DGFrac[iBin][1] << std::endl;
            }
            
            // Update sig_yieldLim parameter with fitted signal yield
            if (FitMRes_SYield[iBin][0] > 0) {
                fitter->updateSigYieldLim("D0", FitMRes_SYield[iBin][0]);
            }
            
            // Stage 1: Create signal-enhanced dataset from mass sPlot weights
            if (enableSPlot && splotFile) {
                std::cout << "  Stage 1: Creating signal-enhanced dataset from mass sPlot weights..." << std::endl;
                
                // Temporarily close the sPlot file to ensure it's properly written
                std::cout << "  Temporarily closing sPlot file for safe reading..." << std::endl;
                splotFile->Close();
                delete splotFile;
                splotFile = nullptr;
                
                // Add small delay to ensure file system operations complete
                std::this_thread::sleep_for(std::chrono::milliseconds(200));
                
                // Create weighted dataset using sPlot weights from mass fit
                signalEnhancedData = fitter->createWeightedDataset(dataBin, splotFileName, iBin, 
                                                                 "sig_sWeight", "signal_enhanced");
                
                // Reopen the sPlot file for continued writing
                splotFile = new TFile(splotFileName.c_str(), "UPDATE");
                if (!splotFile || splotFile->IsZombie()) {
                    std::cerr << "Error: Cannot reopen sPlot file: " << splotFileName << std::endl;
                    enableSPlot = false;  // Disable sPlot for remaining bins
                    splotFile = nullptr;
                } else {
                    std::cout << "  Successfully reopened sPlot file for continued writing" << std::endl;
                }
                
                if (signalEnhancedData) {
                    std::cout << "  Successfully created signal-enhanced dataset with " 
                              << signalEnhancedData->numEntries() << " entries" << std::endl;
                    std::cout << "  Effective entries: " << signalEnhancedData->sumEntries() << std::endl;
                } else {
                    std::cout << "  Warning: Failed to create signal-enhanced dataset, using signal region" << std::endl;
                    signalEnhancedData = static_cast<RooDataSet*>(dataBin->reduce(
                        ("tagMass > " + std::to_string(sigRegion.first) + 
                        " && tagMass < " + std::to_string(sigRegion.second)).c_str()
                    ));
                }
            } else {
                std::cout << "  Using traditional signal region selection..." << std::endl;
                signalEnhancedData = static_cast<RooDataSet*>(dataBin->reduce(

                    ("tagMass > " + std::to_string(sigRegion.first) + 
                    " && tagMass < " + std::to_string(sigRegion.second)).c_str()
                ));
            }
            
            // Stage 2: Fit IP chi2 distribution with signal-enhanced data to get yield variables
            std::cout << "  Stage 2: Fitting IP chi2 distribution with signal-enhanced data..." << std::endl;
            
            // For background modeling in IP chi2, we can use sideband regions
            bkgLeft = static_cast<RooDataSet*>(dataBin->reduce(
                ("tagMass < " + std::to_string(sideBandLimits.first)).c_str()
            ));
            bkgRight = static_cast<RooDataSet*>(dataBin->reduce(
                ("tagMass > " + std::to_string(sideBandLimits.second)).c_str()
            ));
            bkgData = dynamic_cast<RooDataSet*>(bkgLeft->Clone("BKG"));
            bkgData->append(*bkgRight);
            
            // Perform IP chi2 fit with yields for sPlot
            std::tie(ipChi2Histo, ipParams, ipErrors, ipChi2Model, sig_yieldLim, prompt_frac) = 
                fitter->ipchi2FitWithYields("D0", signalEnhancedData, bkgData, "BKGincluded", iBin, zRangeStr, enableSPlot, splotFile);
            
            if(ipChi2Histo) {
                // Create a safe copy to avoid potential pointer invalidation issues
                TH1* ipChi2HistoCopy = (TH1*)ipChi2Histo->Clone(("ipChi2Histo_bin" + std::to_string(iBin)).c_str());
                if (ipChi2HistoCopy) {
                    ipChi2HistoCopy->SetDirectory(nullptr);  // Prevent ROOT from managing it
                    ipchi2HistoArray[iBin] = ipChi2HistoCopy;
                    std::cout << "  Stored safe copy of IP chi2 histogram for bin " << iBin << std::endl;
                } else {
                    std::cout << "  Warning: Failed to create safe copy of IP chi2 histogram for bin " << iBin << std::endl;
                    ipchi2HistoArray[iBin] = nullptr;
                }
            } else {
                std::cerr << "Error: IP chi2 histogram is null for bin " << iBin << std::endl;
                ipchi2HistoArray[iBin] = nullptr;
            }

            // Ensure sPlot file is properly flushed after IP chi2 fit
            if (enableSPlot && splotFile) {
                splotFile->Flush();
                // Force write all pending data to disk
                splotFile->Write("", TObject::kOverwrite);
                std::cout << "  IP chi2 sPlot file flushed and synchronized" << std::endl;
            }
            
            // Stage 3: Process results and create tagZ distribution if sPlot was enabled
            if (enableSPlot && splotFile && ipChi2Model && sig_yieldLim && prompt_frac) {
                std::cout << "  Stage 3: Processing IP chi2 fit results and creating tagZ distribution..." << std::endl;
                
                // Store yield values before any potential file operations that might invalidate pointers
                double totalYield = sig_yieldLim->getVal();
                double promptFraction = prompt_frac->getVal();
                double promptFractionError = prompt_frac->getError();  // Store error as well
                double calculatedPromptYield = totalYield * promptFraction;
                double calculatedNonpromptYield = totalYield * (1 - promptFraction);
                
                // Display yield information
                std::cout << "  Using sig_yieldLim: " << totalYield << std::endl;
                std::cout << "  Using prompt_frac: " << promptFraction << std::endl;
                std::cout << "  Calculated prompt yield: " << calculatedPromptYield << std::endl;
                std::cout << "  Calculated nonprompt yield: " << calculatedNonpromptYield << std::endl;
                
                // Store IP chi2 fit results
                if (ipParams.size() >= 12 && ipErrors.size() >= 12) {
                    FitIPRes_SYield[iBin][0] = totalYield;
                    FitIPRes_PromptFrac[iBin][0] = promptFraction;
                    FitIPRes_PromptFrac[iBin][1] = promptFractionError;  // Use stored error
                    FitIPRes_PromptYield[iBin][0] = calculatedPromptYield;
                    FitIPRes_NonPromptYield[iBin][0] = calculatedNonpromptYield;
                    
                    // Store other IP chi2 fit parameters
                    FitIPRes_XpPrompt[iBin][0] = ipParams[2];
                    FitIPRes_XpPrompt[iBin][1] = ipErrors[2];
                    FitIPRes_SigmaPrompt[iBin][0] = ipParams[3];
                    FitIPRes_SigmaPrompt[iBin][1] = ipErrors[3];
                    FitIPRes_XiPrompt[iBin][0] = ipParams[4];
                    FitIPRes_XiPrompt[iBin][1] = ipErrors[4];
                    
                    FitIPRes_XpNonprompt[iBin][0] = ipParams[7];
                    FitIPRes_XpNonprompt[iBin][1] = ipErrors[7];
                    FitIPRes_SigmaNonprompt[iBin][0] = ipParams[8];
                    FitIPRes_SigmaNonprompt[iBin][1] = ipErrors[8];
                    FitIPRes_XiNonprompt[iBin][0] = ipParams[9];
                    FitIPRes_XiNonprompt[iBin][1] = ipErrors[9];
                }
                
                // Stage 4: Create prompt signal tagZ distribution using combined sPlot weights
                std::cout << "  Stage 4: Creating prompt signal tagZ distribution with combined sPlot weights..." << std::endl;
                    
                    // Create map to store combined weights (for compatibility with existing method)
                    std::map<std::pair<double, double>, std::pair<double, double>> eventWeights;
                    
                    // Read the IP chi2 sPlot weights from the file we just saved
                    std::string ipSplotTreeName = "ipSplotTree_bin" + std::to_string(iBin);
                    TTree* ipSplotTree = (TTree*)splotFile->Get(ipSplotTreeName.c_str());
                    if (ipSplotTree) {
                        double mass_ip, log_ipchi2_ip, prompt_weight_ip, nonprompt_weight_ip;
                        ipSplotTree->SetBranchAddress("mass", &mass_ip);
                        ipSplotTree->SetBranchAddress("log_ipchi2", &log_ipchi2_ip);
                        ipSplotTree->SetBranchAddress("prompt_sWeight", &prompt_weight_ip);
                        ipSplotTree->SetBranchAddress("nonprompt_sWeight", &nonprompt_weight_ip);
                        
                        for (Long64_t i = 0; i < ipSplotTree->GetEntries(); ++i) {
                            ipSplotTree->GetEntry(i);
                            eventWeights[{mass_ip, log_ipchi2_ip}] = {prompt_weight_ip, nonprompt_weight_ip};
                        }
                        
                        std::cout << "    Loaded " << eventWeights.size() << " event weights from IP chi2 sPlot tree" << std::endl;
                        
                        // Use the existing method to create prompt signal tagZ distribution
                        TH1D* promptSignalTagZHist = fitter->createPromptSignalTagZDistribution(
                            dataBin, splotFileName, iBin, eventWeights, "promptSignalTagZ", 20, 0.0, 1.0);
                        
                        if (promptSignalTagZHist) {
                            // Create comparison histograms
                            TH1D* backgroundSubtractedTagZHist = new TH1D(("backgroundSubtractedTagZHist_bin" + std::to_string(iBin)).c_str(), 
                                                                          ("Background-Subtracted TagZ Distribution for Bin " + std::to_string(iBin)).c_str(), 
                                                                          20, 0, 1);
                            
                            // Create background-subtracted distribution using only mass sPlot weights
                            // Ensure file is synchronized before reading
                            if (splotFile) {
                                splotFile->Flush();
                                splotFile->Write("", TObject::kOverwrite);
                            }
                            
                            TFile* splotFileRead = TFile::Open(splotFileName.c_str(), "READ");
                            if (splotFileRead && !splotFileRead->IsZombie()) {
                                TTree* massSplotTree = (TTree*)splotFileRead->Get(("splotTree_bin" + std::to_string(iBin)).c_str());
                                if (massSplotTree) {
                                    double mass_splot, sig_sWeight_mass, bkg_sWeight_mass;
                                    massSplotTree->SetBranchAddress("mass", &mass_splot);
                                    massSplotTree->SetBranchAddress("sig_sWeight", &sig_sWeight_mass);
                                    massSplotTree->SetBranchAddress("bkg_sWeight", &bkg_sWeight_mass);
                                    
                                    std::map<double, double> massSigWeights;
                                    for (Long64_t i = 0; i < massSplotTree->GetEntries(); ++i) {
                                        massSplotTree->GetEntry(i);
                                        massSigWeights[mass_splot] = sig_sWeight_mass;
                                    }
                                    
                                    // Fill background-subtracted histogram
                                    for (int i = 0; i < dataBin->numEntries(); ++i) {
                                        const RooArgSet* row = dataBin->get(i);
                                        RooRealVar* tagZVar = dynamic_cast<RooRealVar*>(row->find("tagZ"));
                                        RooRealVar* massVar = dynamic_cast<RooRealVar*>(row->find("tagMass"));
                                        
                                        if (tagZVar && massVar) {
                                            double tagZ = tagZVar->getVal();
                                            double mass = massVar->getVal();
                                            
                                            auto massSigIt = massSigWeights.find(mass);
                                            if (massSigIt != massSigWeights.end() && massSigIt->second > 0) {
                                                backgroundSubtractedTagZHist->Fill(tagZ, massSigIt->second);
                                            }
                                        }
                                    }
                                }
                                splotFileRead->Close();
                                delete splotFileRead;
                            }
                            
                            // Calculate statistics
                            double promptSignalEntries = promptSignalTagZHist->GetEntries();
                            double promptSignalIntegral = promptSignalTagZHist->Integral();
                            double bkgSubtractedIntegral = backgroundSubtractedTagZHist->Integral();
                            
                            std::cout << "    Prompt signal tagZ histogram:" << std::endl;
                            std::cout << "      Entries: " << promptSignalEntries << std::endl;
                            std::cout << "      Integral (weighted): " << promptSignalIntegral << std::endl;
                            std::cout << "    Background-subtracted tagZ histogram:" << std::endl;
                            std::cout << "      Integral (weighted): " << bkgSubtractedIntegral << std::endl;
                            
                            if (bkgSubtractedIntegral > 0) {
                                double promptFraction = promptSignalIntegral / bkgSubtractedIntegral;
                                std::cout << "      Effective prompt fraction in tagZ: " << promptFraction << std::endl;
                            }

                            // Create safe copies of all histograms to avoid potential RooRealVar reference issues
                            // Check if histograms exist before accessing them
                            if (!promptSignalTagZHist) {
                                std::cerr << "Error: promptSignalTagZHist is null, cannot create safe copy" << std::endl;
                                return;
                            }
                            if (!backgroundSubtractedTagZHist) {
                                std::cerr << "Error: backgroundSubtractedTagZHist is null, cannot create safe copy" << std::endl;
                                return;
                            }
                            if (!tagZHist) {
                                std::cerr << "Error: tagZHist is null, cannot create safe copy" << std::endl;
                                return;
                            }
                                                        
                            TH1D* backgroundSubtractedTagZHistSafe = new TH1D(("backgroundSubtractedTagZHistSafe_bin" + std::to_string(iBin)).c_str(), 
                                                                             ("Signal TagZ (Mass sPlot) - Bin " + std::to_string(iBin)).c_str(), 
                                                                             backgroundSubtractedTagZHist->GetNbinsX(), 
                                                                             backgroundSubtractedTagZHist->GetXaxis()->GetXmin(), 
                                                                             backgroundSubtractedTagZHist->GetXaxis()->GetXmax());
                            
                            TH1D* promptSignalTagZHistSafe = new TH1D(("promptSignalTagZHistSafe_bin" + std::to_string(iBin)).c_str(), 
                                                                     ("Prompt Signal TagZ (Combined sPlot) - Bin " + std::to_string(iBin)).c_str(), 
                                                                     promptSignalTagZHist->GetNbinsX(), 
                                                                     promptSignalTagZHist->GetXaxis()->GetXmin(), 
                                                                     promptSignalTagZHist->GetXaxis()->GetXmax());
                                                        
                            
                            for (int iBinHist = 1; iBinHist <= backgroundSubtractedTagZHist->GetNbinsX(); ++iBinHist) {
                                backgroundSubtractedTagZHistSafe->SetBinContent(iBinHist, backgroundSubtractedTagZHist->GetBinContent(iBinHist));
                                backgroundSubtractedTagZHistSafe->SetBinError(iBinHist, backgroundSubtractedTagZHist->GetBinError(iBinHist));
                            }
                            backgroundSubtractedTagZHistSafe->SetEntries(backgroundSubtractedTagZHist->GetEntries());
                            
                            for (int iBinHist = 1; iBinHist <= promptSignalTagZHist->GetNbinsX(); ++iBinHist) {
                                promptSignalTagZHistSafe->SetBinContent(iBinHist, promptSignalTagZHist->GetBinContent(iBinHist));
                                promptSignalTagZHistSafe->SetBinError(iBinHist, promptSignalTagZHist->GetBinError(iBinHist));
                            }
                            promptSignalTagZHistSafe->SetEntries(promptSignalTagZHist->GetEntries());
                            
                            // Create comprehensive comparison plot
                            TCanvas* promptTagZCanvas = new TCanvas(("promptTagZCanvas_bin" + std::to_string(iBin)).c_str(), 
                                                                  ("Prompt Signal TagZ Analysis - Bin " + std::to_string(iBin)).c_str(), 
                                                                  1200, 800);
                            promptTagZCanvas->Divide(2, 2);

                            // Plot 1: Original tagZ distribution (using safe copy)
                            promptTagZCanvas->cd(1);
                            tagZHistSafe->SetLineColor(kBlack);
                            tagZHistSafe->SetMarkerColor(kBlack);
                            tagZHistSafe->SetMarkerStyle(20);
                            tagZHistSafe->SetTitle("Original TagZ Distribution");
                            tagZHistSafe->GetXaxis()->SetTitle("#it{z}_{T}");
                            tagZHistSafe->GetYaxis()->SetTitle("Entries");
                            tagZHistSafe->Draw("pe");

                            // Plot 2: Background-subtracted tagZ distribution
                            promptTagZCanvas->cd(2);
                            backgroundSubtractedTagZHistSafe->SetLineColor(kBlue);
                            backgroundSubtractedTagZHistSafe->SetMarkerColor(kBlue);
                            backgroundSubtractedTagZHistSafe->SetMarkerStyle(21);
                            backgroundSubtractedTagZHistSafe->SetTitle("Signal TagZ (Mass sPlot)");
                            backgroundSubtractedTagZHistSafe->GetXaxis()->SetTitle("#it{z}_{T}");
                            backgroundSubtractedTagZHistSafe->GetYaxis()->SetTitle("Weighted Entries");
                            backgroundSubtractedTagZHistSafe->Draw("pe");

                            // Plot 3: Prompt signal tagZ distribution
                            promptTagZCanvas->cd(3);
                            promptSignalTagZHistSafe->SetLineColor(kRed);
                            promptSignalTagZHistSafe->SetMarkerColor(kRed);
                            promptSignalTagZHistSafe->SetMarkerStyle(22);
                            promptSignalTagZHistSafe->SetTitle("Prompt Signal TagZ (Combined sPlot)");
                            promptSignalTagZHistSafe->GetXaxis()->SetTitle("#it{z}_{T}");
                            promptSignalTagZHistSafe->GetYaxis()->SetTitle("Weighted Entries");
                            promptSignalTagZHistSafe->Draw("pe");

                            double meanTagZ_weighted = promptSignalTagZHistSafe->GetMean();
        double stdDevTagZ_weighted = promptSignalTagZHistSafe->GetStdDev();
        std::cout << "Mean tagZ: " << meanTagZ_weighted << ", StdDev: " << stdDevTagZ_weighted << std::endl;
        Bin_TagZMean_weighted[iBin][0] = meanTagZ_weighted;  // Store mean tagZ for this bin

                            // Plot 4: Normalized comparison
                            promptTagZCanvas->cd(4);
                            TH1D* tagZHistNorm = (TH1D*)tagZHistSafe->Clone("tagZHistNorm");
                            TH1D* bkgSubNorm = (TH1D*)backgroundSubtractedTagZHistSafe->Clone("bkgSubNorm");
                            TH1D* promptNorm = (TH1D*)promptSignalTagZHistSafe->Clone("promptNorm");
                            // Normalize histograms
                            if (tagZHistNorm->Integral() > 0) tagZHistNorm->Scale(1.0 / tagZHistNorm->Integral());
                            if (bkgSubNorm->Integral() > 0) bkgSubNorm->Scale(1.0 / bkgSubNorm->Integral());
                            if (promptNorm->Integral() > 0) promptNorm->Scale(1.0 / promptNorm->Integral());
                            
                            // Set maximum for better visualization
                            double maxY = std::max({tagZHistNorm->GetMaximum(), bkgSubNorm->GetMaximum(), promptNorm->GetMaximum()});
                            tagZHistNorm->GetYaxis()->SetRangeUser(0, maxY * 1.2);
                            
                            tagZHistNorm->SetTitle("TagZ Distributions Comparison (Normalized)");
                            tagZHistNorm->GetYaxis()->SetTitle("Normalized Entries");
                            tagZHistNorm->Draw("pe");
                            bkgSubNorm->Draw("pe same");
                            promptNorm->Draw("pe same");
                            
                            // Add legend
                            TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
                            legend->AddEntry(tagZHistNorm, "Original", "pe");
                            legend->AddEntry(bkgSubNorm, "Signal (Mass sPlot)", "pe");
                            legend->AddEntry(promptNorm, "Prompt Signal (Combined)", "pe");
                            legend->SetBorderSize(0);
                            legend->SetFillStyle(0);
                            legend->Draw();
                            
                            // Calculate correction factors for tagZ distributions
                            std::cout << "  Calculating efficiency corrections for tagZ distributions..." << std::endl;
                            
                            // Create histograms for correction factors in tagZ bins
                            TH1D* tagZKaonCorrHist = new TH1D(("tagZKaonCorr_bin" + std::to_string(iBin)).c_str(), 
                                                              ("Kaon Efficiency vs TagZ - Bin " + std::to_string(iBin)).c_str(), 
                                                              20, 0, 1);
                            TH1D* tagZPionCorrHist = new TH1D(("tagZPionCorr_bin" + std::to_string(iBin)).c_str(), 
                                                              ("Pion Efficiency vs TagZ - Bin " + std::to_string(iBin)).c_str(), 
                                                              20, 0, 1);
                            TH1D* tagZRecoEffCorrHist = new TH1D(("tagZRecoEffCorr_bin" + std::to_string(iBin)).c_str(),
                                                                ("Reconstruction Efficiency vs TagZ - Bin " + std::to_string(iBin)).c_str(),
                                                                20, 0, 1);
                            TH1D* tagZAcceptanceCorrHist = new TH1D(("tagZAcceptanceCorr_bin" + std::to_string(iBin)).c_str(),
                                                                   ("Acceptance vs TagZ - Bin " + std::to_string(iBin)).c_str(),
                                                                   20, 0, 1);
                            TH1D* tagZCombinedCorrHist = new TH1D(("tagZCombinedCorr_bin" + std::to_string(iBin)).c_str(), 
                                                                  ("Combined PID Efficiency vs TagZ - Bin " + std::to_string(iBin)).c_str(), 
                                                                  20, 0, 1);
                            
                            // Create counter histograms for averaging
                            TH1D* tagZKaonCountHist = new TH1D(("tagZKaonCount_bin" + std::to_string(iBin)).c_str(), 
                                                               "Kaon Count", 20, 0, 1);
                            TH1D* tagZPionCountHist = new TH1D(("tagZPionCount_bin" + std::to_string(iBin)).c_str(), 
                                                               "Pion Count", 20, 0, 1);
                            TH1D* tagZRecoEffCountHist = new TH1D(("tagZRecoEffCount_bin" + std::to_string(iBin)).c_str(),
                                                                 "Reco Eff Count", 20, 0, 1);
                            TH1D* tagZAcceptanceCountHist = new TH1D(("tagZAcceptanceCount_bin" + std::to_string(iBin)).c_str(),
                                                                    "Acceptance Count", 20, 0, 1);
                            TH1D* tagZCombinedCountHist = new TH1D(("tagZCombinedCount_bin" + std::to_string(iBin)).c_str(), 
                                                                   "Combined Count", 20, 0, 1);
                            
                            // Fill correction histograms and create efficiency-weighted distributions by iterating through the bin dataset
                            for (int i = 0; i < dataBin->numEntries(); ++i) {
                                const RooArgSet* row = dataBin->get(i);
                                RooRealVar* tagZVar = dynamic_cast<RooRealVar*>(row->find("tagZ"));
                                RooAbsReal* kaonEff = dynamic_cast<RooAbsReal*>(row->find("kaon_efficiency"));
                                RooAbsReal* pionEff = dynamic_cast<RooAbsReal*>(row->find("pion_efficiency"));
                                RooAbsReal* recoEff = dynamic_cast<RooAbsReal*>(row->find("reconstruction_efficiency"));
                                RooAbsReal* acceptance = dynamic_cast<RooAbsReal*>(row->find("acceptance"));
                                RooAbsReal* combinedPIDEff = dynamic_cast<RooAbsReal*>(row->find("combined_PID_efficiency"));
                                
                                if (tagZVar && kaonEff && pionEff) {
                                    double tagZ = tagZVar->getVal();
                                    double kaonEffVal = kaonEff->getVal();
                                    double pionEffVal = pionEff->getVal();
                                    double recoEffVal = recoEff ? recoEff->getVal() : 1.0;
                                    double acceptanceVal = acceptance ? acceptance->getVal() : 1.0;
                                    double combinedPIDEffVal = combinedPIDEff ? combinedPIDEff->getVal() : kaonEffVal * pionEffVal;
                                    
                                    // Only use events with valid efficiency values
                                    if (kaonEffVal > 0 && kaonEffVal <= 1.0 && 
                                        pionEffVal > 0 && pionEffVal <= 1.0 &&
                                        tagZ >= 0 && tagZ <= 1.0) {
                                        
                                        // Fill efficiency histograms for averaging
                                        tagZKaonCorrHist->Fill(tagZ, kaonEffVal);
                                        tagZKaonCountHist->Fill(tagZ);
                                        
                                        tagZPionCorrHist->Fill(tagZ, pionEffVal);
                                        tagZPionCountHist->Fill(tagZ);
                                        
                                        if (recoEff) {
                                            tagZRecoEffCorrHist->Fill(tagZ, recoEffVal);
                                            tagZRecoEffCountHist->Fill(tagZ);
                                        }
                                        
                                        if (acceptance) {
                                            tagZAcceptanceCorrHist->Fill(tagZ, acceptanceVal);
                                            tagZAcceptanceCountHist->Fill(tagZ);
                                        }
                                        
                                        tagZCombinedCorrHist->Fill(tagZ, combinedPIDEffVal);
                                        tagZCombinedCountHist->Fill(tagZ);
                                    }
                                }
                            }
                            
                            // Calculate average efficiency in each tagZ bin
                            for (int bin = 1; bin <= tagZKaonCorrHist->GetNbinsX(); ++bin) {
                                if (tagZKaonCountHist->GetBinContent(bin) > 0) {
                                    double avgKaonEff = tagZKaonCorrHist->GetBinContent(bin) / tagZKaonCountHist->GetBinContent(bin);
                                    tagZKaonCorrHist->SetBinContent(bin, avgKaonEff);
                                    tagZKaonCorrHist->SetBinError(bin, 0.01); // Placeholder error
                                }
                                
                                if (tagZPionCountHist->GetBinContent(bin) > 0) {
                                    double avgPionEff = tagZPionCorrHist->GetBinContent(bin) / tagZPionCountHist->GetBinContent(bin);
                                    tagZPionCorrHist->SetBinContent(bin, avgPionEff);
                                    tagZPionCorrHist->SetBinError(bin, 0.01); // Placeholder error
                                }
                                
                                if (tagZRecoEffCountHist->GetBinContent(bin) > 0) {
                                    double avgRecoEff = tagZRecoEffCorrHist->GetBinContent(bin) / tagZRecoEffCountHist->GetBinContent(bin);
                                    tagZRecoEffCorrHist->SetBinContent(bin, avgRecoEff);
                                    tagZRecoEffCorrHist->SetBinError(bin, 0.01); // Placeholder error
                                }
                                
                                if (tagZAcceptanceCountHist->GetBinContent(bin) > 0) {
                                    double avgAcceptance = tagZAcceptanceCorrHist->GetBinContent(bin) / tagZAcceptanceCountHist->GetBinContent(bin);
                                    tagZAcceptanceCorrHist->SetBinContent(bin, avgAcceptance);
                                    tagZAcceptanceCorrHist->SetBinError(bin, 0.01); // Placeholder error
                                }
                                
                                if (tagZCombinedCountHist->GetBinContent(bin) > 0) {
                                    double avgCombinedEff = tagZCombinedCorrHist->GetBinContent(bin) / tagZCombinedCountHist->GetBinContent(bin);
                                    tagZCombinedCorrHist->SetBinContent(bin, avgCombinedEff);
                                    tagZCombinedCorrHist->SetBinError(bin, 0.01); // Placeholder error
                                }
                            }
                            
                            // Create TGraphErrors from histograms for easier handling
                            TGraphErrors* tagZKaonCorrGraph = new TGraphErrors();
                            TGraphErrors* tagZPionCorrGraph = new TGraphErrors();
                            TGraphErrors* tagZRecoEffCorrGraph = new TGraphErrors();
                            TGraphErrors* tagZAcceptanceCorrGraph = new TGraphErrors();
                            TGraphErrors* tagZCombinedCorrGraph = new TGraphErrors();
                            
                            int pointIndex = 0;
                            for (int bin = 1; bin <= tagZKaonCorrHist->GetNbinsX(); ++bin) {
                                if (tagZKaonCorrHist->GetBinContent(bin) > 0) {
                                    double binCenter = tagZKaonCorrHist->GetBinCenter(bin);
                                    double binWidth = tagZKaonCorrHist->GetBinWidth(bin) / 2.0;
                                    
                                    tagZKaonCorrGraph->SetPoint(pointIndex, binCenter, tagZKaonCorrHist->GetBinContent(bin));
                                    tagZKaonCorrGraph->SetPointError(pointIndex, binWidth, tagZKaonCorrHist->GetBinError(bin));
                                    
                                    tagZPionCorrGraph->SetPoint(pointIndex, binCenter, tagZPionCorrHist->GetBinContent(bin));
                                    tagZPionCorrGraph->SetPointError(pointIndex, binWidth, tagZPionCorrHist->GetBinError(bin));
                                    
                                    if (tagZRecoEffCorrHist->GetBinContent(bin) > 0) {
                                        tagZRecoEffCorrGraph->SetPoint(pointIndex, binCenter, tagZRecoEffCorrHist->GetBinContent(bin));
                                        tagZRecoEffCorrGraph->SetPointError(pointIndex, binWidth, tagZRecoEffCorrHist->GetBinError(bin));
                                    }
                                    
                                    if (tagZAcceptanceCorrHist->GetBinContent(bin) > 0) {
                                        tagZAcceptanceCorrGraph->SetPoint(pointIndex, binCenter, tagZAcceptanceCorrHist->GetBinContent(bin));
                                        tagZAcceptanceCorrGraph->SetPointError(pointIndex, binWidth, tagZAcceptanceCorrHist->GetBinError(bin));
                                    }
                                    
                                    tagZCombinedCorrGraph->SetPoint(pointIndex, binCenter, tagZCombinedCorrHist->GetBinContent(bin));
                                    tagZCombinedCorrGraph->SetPointError(pointIndex, binWidth, tagZCombinedCorrHist->GetBinError(bin));
                                    
                                    pointIndex++;
                                }
                            }
                            
                            // Set graph properties
                            tagZKaonCorrGraph->SetName(("tagZKaonCorr_bin" + std::to_string(iBin)).c_str());
                            tagZKaonCorrGraph->SetTitle(("Kaon PID Efficiency vs z_{T} - Bin " + std::to_string(iBin)).c_str());
                            tagZKaonCorrGraph->SetMarkerStyle(20);
                            tagZKaonCorrGraph->SetMarkerColor(kBlue);
                            tagZKaonCorrGraph->SetLineColor(kBlue);
                            
                            tagZPionCorrGraph->SetName(("tagZPionCorr_bin" + std::to_string(iBin)).c_str());
                            tagZPionCorrGraph->SetTitle(("Pion PID Efficiency vs z_{T} - Bin " + std::to_string(iBin)).c_str());
                            tagZPionCorrGraph->SetMarkerStyle(21);
                            tagZPionCorrGraph->SetMarkerColor(kRed);
                            tagZPionCorrGraph->SetLineColor(kRed);
                            
                            tagZRecoEffCorrGraph->SetName(("tagZRecoEffCorr_bin" + std::to_string(iBin)).c_str());
                            tagZRecoEffCorrGraph->SetTitle(("Reconstruction Efficiency vs z_{T} - Bin " + std::to_string(iBin)).c_str());
                            tagZRecoEffCorrGraph->SetMarkerStyle(23);
                            tagZRecoEffCorrGraph->SetMarkerColor(kOrange);
                            tagZRecoEffCorrGraph->SetLineColor(kOrange);
                            
                            tagZAcceptanceCorrGraph->SetName(("tagZAcceptanceCorr_bin" + std::to_string(iBin)).c_str());
                            tagZAcceptanceCorrGraph->SetTitle(("Acceptance vs z_{T} - Bin " + std::to_string(iBin)).c_str());
                            tagZAcceptanceCorrGraph->SetMarkerStyle(24);
                            tagZAcceptanceCorrGraph->SetMarkerColor(kMagenta);
                            tagZAcceptanceCorrGraph->SetLineColor(kMagenta);
                            
                            tagZCombinedCorrGraph->SetName(("tagZCombinedCorr_bin" + std::to_string(iBin)).c_str());
                            tagZCombinedCorrGraph->SetTitle(("Combined PID Efficiency vs z_{T} - Bin " + std::to_string(iBin)).c_str());
                            tagZCombinedCorrGraph->SetMarkerStyle(22);
                            tagZCombinedCorrGraph->SetMarkerColor(kGreen + 2);
                            tagZCombinedCorrGraph->SetLineColor(kGreen + 2);
                            
                            // Save correction factor graphs to ROOT file in main directory
                            std::string mainDir;
                            if (isMC) {
                                mainDir = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_MC";
                            } else {
                                mainDir = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA";
                            }
                            
                            // Create pT range string for graph names
                            std::stringstream ptRangeStr;
                            ptRangeStr << jetPt.first << "_" << jetPt.second;
                            std::string ptString = ptRangeStr.str();
                            
                            std::string rootFileName = mainDir + "/TagZCorrectionFactors.root";
                            TFile* tagZCorrFile = new TFile(rootFileName.c_str(), "UPDATE");
                            if (tagZCorrFile && tagZCorrFile->IsOpen()) {
                                // Include pT range and bin in graph names
                                std::string kaonName = "tagZKaonCorrection_" + ptString + "_bin" + std::to_string(iBin);
                                std::string pionName = "tagZPionCorrection_" + ptString + "_bin" + std::to_string(iBin);
                                std::string recoName = "tagZRecoEffCorrection_" + ptString + "_bin" + std::to_string(iBin);
                                std::string acceptanceName = "tagZAcceptanceCorrection_" + ptString + "_bin" + std::to_string(iBin);
                                std::string combinedName = "tagZCombinedCorrection_" + ptString + "_bin" + std::to_string(iBin);
                                
                                tagZKaonCorrGraph->Write(kaonName.c_str());
                                tagZPionCorrGraph->Write(pionName.c_str());
                                if (tagZRecoEffCorrGraph->GetN() > 0) tagZRecoEffCorrGraph->Write(recoName.c_str());
                                if (tagZAcceptanceCorrGraph->GetN() > 0) tagZAcceptanceCorrGraph->Write(acceptanceName.c_str());
                                tagZCombinedCorrGraph->Write(combinedName.c_str());
                                
                                tagZCorrFile->Close();
                                std::cout << "  Saved tagZ correction factor graphs to: " << rootFileName << std::endl;
                                std::cout << "    Graph names: " << kaonName << ", " << pionName << ", " << recoName << ", " << acceptanceName << ", " << combinedName << std::endl;
                            } else {
                                std::cerr << "  Error: Could not create tagZ correction factors ROOT file" << std::endl;
                            }
                            
                            if (tagZCorrFile) {
                                delete tagZCorrFile;
                                tagZCorrFile = nullptr;
                            }
                            
                            // Create plot showing correction factors vs tagZ (original layout)
                            TCanvas* tagZCorrCanvas = new TCanvas(("tagZCorrCanvas_bin" + std::to_string(iBin)).c_str(), 
                                                                 ("TagZ PID Correction Factors - Bin " + std::to_string(iBin)).c_str(), 
                                                                 800, 600);
                            tagZCorrCanvas->SetLeftMargin(0.12);
                            tagZCorrCanvas->SetRightMargin(0.05);
                            tagZCorrCanvas->SetTopMargin(0.08);
                            tagZCorrCanvas->SetBottomMargin(0.12);
                            
                            // Set axis labels and ranges
                            tagZKaonCorrGraph->GetXaxis()->SetTitle("z_{T}");
                            tagZKaonCorrGraph->GetYaxis()->SetTitle("PID Efficiency");
                            tagZKaonCorrGraph->GetYaxis()->SetTitleOffset(1.2);
                            
                            tagZKaonCorrGraph->Draw("APE");
                            tagZPionCorrGraph->Draw("PE same");
                            tagZCombinedCorrGraph->Draw("PE same");
                            
                            // Add legend
                            TLegend* corrLegend = new TLegend(0.7, 0.7, 0.95, 0.9);
                            corrLegend->SetBorderSize(0);
                            corrLegend->SetFillStyle(0);
                            corrLegend->AddEntry(tagZKaonCorrGraph, "Kaon PID", "pe");
                            corrLegend->AddEntry(tagZPionCorrGraph, "Pion PID", "pe");
                            corrLegend->AddEntry(tagZCombinedCorrGraph, "Combined", "pe");
                            corrLegend->Draw();
                            
                            // Save TagZ correction plots
                            std::string tagZCorrOutputFile = outfilePath + "tagZPIDCorrections_bin" + std::to_string(iBin) + ".png";
                            tagZCorrCanvas->SaveAs(tagZCorrOutputFile.c_str());
                            
                            // Clean up canvas and legends for this plot only
                            delete tagZCorrCanvas;
                            delete corrLegend;
                            
                            std::cout << "  TagZ correction factors calculation and saving completed." << std::endl;
                            
                            // Save canvas
                            std::string promptTagZOutputFile = outfilePath + "promptSignalTagZ_bin" + std::to_string(iBin) + ".png";
                            promptTagZCanvas->SaveAs(promptTagZOutputFile.c_str());
                            
                            // Save histograms to sPlot file
                            splotFile->cd();
                            promptSignalTagZHist->Write();
                            backgroundSubtractedTagZHist->Write();
                            
                            std::cout << "  Stage 4 completed: Prompt signal tagZ distribution created and saved" << std::endl;
                            

                            // Create efficiency-weighted prompt signal distributions before cleanup
                            std::cout << "  Creating efficiency-weighted prompt signal distributions..." << std::endl;
                            
                            // Create efficiency-weighted histograms
                            TH1D* promptSignalTagZHist_PIDWeighted = new TH1D(("promptSignalTagZHist_PIDWeighted_bin" + std::to_string(iBin)).c_str(),
                                                                             ("Prompt Signal TagZ (PID Weighted) - Bin " + std::to_string(iBin)).c_str(),
                                                                             20, 0, 1);
                            TH1D* promptSignalTagZHist_RecoWeighted = new TH1D(("promptSignalTagZHist_RecoWeighted_bin" + std::to_string(iBin)).c_str(),
                                                                              ("Prompt Signal TagZ (Reco Weighted) - Bin " + std::to_string(iBin)).c_str(),
                                                                              20, 0, 1);
                            TH1D* promptSignalTagZHist_AcceptanceWeighted = new TH1D(("promptSignalTagZHist_AcceptanceWeighted_bin" + std::to_string(iBin)).c_str(),
                                                                                    ("Prompt Signal TagZ (Acceptance Weighted) - Bin " + std::to_string(iBin)).c_str(),
                                                                                    20, 0, 1);
                            TH1D* promptSignalTagZHist_FullyWeighted = new TH1D(("promptSignalTagZHist_FullyWeighted_bin" + std::to_string(iBin)).c_str(),
                                                                               ("Prompt Signal TagZ (Fully Weighted) - Bin " + std::to_string(iBin)).c_str(),
                                                                               20, 0, 1);
                            
                            // Fill efficiency-weighted histograms from the original prompt signal histogram
                            for (int bin = 1; bin <= promptSignalTagZHistSafe->GetNbinsX(); ++bin) {
                                double binCenter = promptSignalTagZHistSafe->GetBinCenter(bin);
                                double originalWeight = promptSignalTagZHistSafe->GetBinContent(bin);
                                double originalError = promptSignalTagZHistSafe->GetBinError(bin);
                                
                                if (originalWeight > 0) {
                                    // Get average efficiencies for this tagZ bin from the correction histograms
                                    int corrBin = tagZKaonCorrHist->FindBin(binCenter);
                                    double kaonEff = tagZKaonCorrHist->GetBinContent(corrBin);
                                    double pionEff = tagZPionCorrHist->GetBinContent(corrBin);
                                    double recoEff = tagZRecoEffCorrHist->GetBinContent(corrBin);
                                    double acceptance = tagZAcceptanceCorrHist->GetBinContent(corrBin);
                                    double combinedPIDEff = tagZCombinedCorrHist->GetBinContent(corrBin);
                                    
                                    // Apply efficiency corrections (divide by efficiency to correct)
                                    double pidWeightedVal = (combinedPIDEff > 0) ? originalWeight / combinedPIDEff : originalWeight;
                                    double recoWeightedVal = (recoEff > 0) ? originalWeight / recoEff : originalWeight;
                                    double acceptanceWeightedVal = (acceptance > 0) ? originalWeight / acceptance : originalWeight;
                                    double fullyWeightedVal = originalWeight;
                                    if (combinedPIDEff > 0) fullyWeightedVal /= combinedPIDEff;
                                    if (recoEff > 0) fullyWeightedVal /= recoEff;
                                    if (acceptance > 0) fullyWeightedVal /= acceptance;
                                    
                                    // Set bin contents and errors
                                    promptSignalTagZHist_PIDWeighted->SetBinContent(bin, pidWeightedVal);
                                    promptSignalTagZHist_PIDWeighted->SetBinError(bin, originalError);
                                    
                                    promptSignalTagZHist_RecoWeighted->SetBinContent(bin, recoWeightedVal);
                                    promptSignalTagZHist_RecoWeighted->SetBinError(bin, originalError);
                                    
                                    promptSignalTagZHist_AcceptanceWeighted->SetBinContent(bin, acceptanceWeightedVal);
                                    promptSignalTagZHist_AcceptanceWeighted->SetBinError(bin, originalError);
                                    
                                    promptSignalTagZHist_FullyWeighted->SetBinContent(bin, fullyWeightedVal);
                                    promptSignalTagZHist_FullyWeighted->SetBinError(bin, originalError);
                                }
                            }
                            // Save tagZ histograms to ROOT file
                            std::string tagZHistRootFileName = mainDir + "/TagZHistograms_" + ptString + ".root";
                            TFile* tagZHistFile = new TFile(tagZHistRootFileName.c_str(), "UPDATE");
                            if (tagZHistFile && tagZHistFile->IsOpen()) {
                                // Save the safe copies of tagZ histograms with descriptive names
                                std::string tagZHistName = "tagZHist_" + ptString + "_bin" + std::to_string(iBin);
                                std::string bkgSubTagZHistName = "backgroundSubtractedTagZHist_" + ptString + "_bin" + std::to_string(iBin);
                                std::string promptTagZHistName = "promptSignalTagZHist_" + ptString + "_bin" + std::to_string(iBin);
                                std::string recoWeightedTagZHistName = "promptSignalTagZHist_RecoWeighted_" + ptString + "_bin" + std::to_string(iBin);
                                std::string acceptanceWeightedTagZHistName = "promptSignalTagZHist_AcceptanceWeighted_" + ptString + "_bin" + std::to_string(iBin);
                                std::string fullyWeightedTagZHistName = "promptSignalTagZHist_FullyWeighted_" + ptString + "_bin" + std::to_string(iBin);

                                tagZHistSafe->Write(tagZHistName.c_str());
                                backgroundSubtractedTagZHistSafe->Write(bkgSubTagZHistName.c_str());
                                promptSignalTagZHistSafe->Write(promptTagZHistName.c_str());
                                promptSignalTagZHist_RecoWeighted->Write(recoWeightedTagZHistName.c_str());
                                promptSignalTagZHist_AcceptanceWeighted->Write(acceptanceWeightedTagZHistName.c_str());
                                promptSignalTagZHist_FullyWeighted->Write(fullyWeightedTagZHistName.c_str());

                                tagZHistFile->Close();
                                std::cout << "  Saved tagZ histograms to: " << tagZHistRootFileName << std::endl;
                                std::cout << "    Histogram names: " << tagZHistName << ", " << bkgSubTagZHistName << ", " << promptTagZHistName << ", " << recoWeightedTagZHistName << ", " << acceptanceWeightedTagZHistName << ", " << fullyWeightedTagZHistName << std::endl;
                            } else {
                                std::cerr << "  Error: Could not create tagZ histograms ROOT file: " << tagZHistRootFileName << std::endl;
                            }

                            if (tagZHistFile) {
                                delete tagZHistFile;
                                tagZHistFile = nullptr;
                            }
                            
                            // Create comprehensive efficiency correction plot
                            TCanvas* tagZEffCorrCanvas = new TCanvas(("tagZEffCorrCanvas_bin" + std::to_string(iBin)).c_str(), 
                                                                    ("All Efficiency Corrections vs TagZ - Bin " + std::to_string(iBin)).c_str(), 
                                                                    1000, 800);
                            tagZEffCorrCanvas->SetLeftMargin(0.08);
                            tagZEffCorrCanvas->SetRightMargin(0.01);
                            tagZEffCorrCanvas->SetTopMargin(0.08);
                            tagZEffCorrCanvas->SetBottomMargin(0.08);
                            
                            // Set axis labels and ranges
                            tagZKaonCorrGraph->GetXaxis()->SetTitle("z_{T}");
                            tagZKaonCorrGraph->GetYaxis()->SetTitle("Efficiency");
                            tagZKaonCorrGraph->GetYaxis()->SetRangeUser(0.0, 1.1);
                            
                            tagZKaonCorrGraph->Draw("APE");
                            tagZPionCorrGraph->Draw("PE same");
                            if (tagZRecoEffCorrGraph->GetN() > 0) tagZRecoEffCorrGraph->Draw("PE same");
                            if (tagZAcceptanceCorrGraph->GetN() > 0) tagZAcceptanceCorrGraph->Draw("PE same");
                            tagZCombinedCorrGraph->Draw("PE same");
                            
                            // Add comprehensive legend
                            TLegend* effCorrLegend = new TLegend(0.65, 0.65, 0.95, 0.9);
                            effCorrLegend->SetBorderSize(0);
                            effCorrLegend->SetFillStyle(0);
                            effCorrLegend->AddEntry(tagZKaonCorrGraph, "Kaon PID", "pe");
                            effCorrLegend->AddEntry(tagZPionCorrGraph, "Pion PID", "pe");
                            if (tagZRecoEffCorrGraph->GetN() > 0) effCorrLegend->AddEntry(tagZRecoEffCorrGraph, "Reconstruction", "pe");
                            if (tagZAcceptanceCorrGraph->GetN() > 0) effCorrLegend->AddEntry(tagZAcceptanceCorrGraph, "Acceptance", "pe");
                            effCorrLegend->AddEntry(tagZCombinedCorrGraph, "Combined PID", "pe");
                            effCorrLegend->Draw();
                            
                            // Save comprehensive efficiency plot
                            std::string tagZEffCorrOutputFile = outfilePath + "tagZAllEfficiencyCorrections_bin" + std::to_string(iBin) + ".png";
                            tagZEffCorrCanvas->SaveAs(tagZEffCorrOutputFile.c_str());

                            // Also save all inputs used to make this canvas into a ROOT file so they can be re-plotted later
                            std::string tagZEffInputsRoot = outfilePath + "tagZAllEfficiencyCorrections_bin" + std::to_string(iBin) + ".root";
                            TFile* tagZInputsFile = TFile::Open(tagZEffInputsRoot.c_str(), "RECREATE");
                            if (tagZInputsFile && tagZInputsFile->IsOpen()) {
                                // Write graphs (if they contain points)
                                if (tagZKaonCorrGraph && tagZKaonCorrGraph->GetN() > 0) tagZKaonCorrGraph->Write("tagZKaonCorrGraph");
                                if (tagZPionCorrGraph && tagZPionCorrGraph->GetN() > 0) tagZPionCorrGraph->Write("tagZPionCorrGraph");
                                if (tagZRecoEffCorrGraph && tagZRecoEffCorrGraph->GetN() > 0) tagZRecoEffCorrGraph->Write("tagZRecoEffCorrGraph");
                                if (tagZAcceptanceCorrGraph && tagZAcceptanceCorrGraph->GetN() > 0) tagZAcceptanceCorrGraph->Write("tagZAcceptanceCorrGraph");
                                if (tagZCombinedCorrGraph && tagZCombinedCorrGraph->GetN() > 0) tagZCombinedCorrGraph->Write("tagZCombinedCorrGraph");

                                // Also write the underlying binned histograms used to build the graphs (if present)
                                if (tagZKaonCorrHist) tagZKaonCorrHist->Write("tagZKaonCorrHist");
                                if (tagZPionCorrHist) tagZPionCorrHist->Write("tagZPionCorrHist");
                                if (tagZRecoEffCorrHist) tagZRecoEffCorrHist->Write("tagZRecoEffCorrHist");
                                if (tagZAcceptanceCorrHist) tagZAcceptanceCorrHist->Write("tagZAcceptanceCorrHist");
                                if (tagZCombinedCorrHist) tagZCombinedCorrHist->Write("tagZCombinedCorrHist");

                                tagZInputsFile->Close();
                                std::cout << "  Saved tagZ inputs to ROOT file: " << tagZEffInputsRoot << std::endl;
                            } else {
                                std::cerr << "  Error: could not create tagZ inputs ROOT file: " << tagZEffInputsRoot << std::endl;
                            }
                            if (tagZInputsFile) { delete tagZInputsFile; tagZInputsFile = nullptr; }
                            
                            // Create efficiency-weighted prompt signal comparison plot
                            TCanvas* promptSignalWeightedCanvas = new TCanvas(("promptSignalWeightedCanvas_bin" + std::to_string(iBin)).c_str(),
                                                                             ("Efficiency-Weighted Prompt Signal TagZ - Bin " + std::to_string(iBin)).c_str(),
                                                                             1200, 900);
                            promptSignalWeightedCanvas->Divide(2, 3);
                            
                            // Plot 1: Original prompt signal
                            promptSignalWeightedCanvas->cd(1);
                            promptSignalTagZHistSafe->SetLineColor(kBlack);
                            promptSignalTagZHistSafe->SetMarkerColor(kBlack);
                            promptSignalTagZHistSafe->SetMarkerStyle(20);
                            promptSignalTagZHistSafe->SetTitle("Original Prompt Signal");
                            promptSignalTagZHistSafe->GetXaxis()->SetTitle("#it{z}_{T}");
                            promptSignalTagZHistSafe->GetYaxis()->SetTitle("Weighted Entries");
                            promptSignalTagZHistSafe->Draw("pe");
                            
                            // Plot 2: PID-weighted
                            promptSignalWeightedCanvas->cd(2);
                            promptSignalTagZHist_PIDWeighted->SetLineColor(kBlue);
                            promptSignalTagZHist_PIDWeighted->SetMarkerColor(kBlue);
                            promptSignalTagZHist_PIDWeighted->SetMarkerStyle(21);
                            promptSignalTagZHist_PIDWeighted->SetTitle("PID Efficiency Corrected");
                            promptSignalTagZHist_PIDWeighted->GetXaxis()->SetTitle("#it{z}_{T}");
                            promptSignalTagZHist_PIDWeighted->GetYaxis()->SetTitle("Corrected Entries");
                            promptSignalTagZHist_PIDWeighted->Draw("pe");
                            
                            // Plot 3: Reconstruction-weighted
                            promptSignalWeightedCanvas->cd(3);
                            promptSignalTagZHist_RecoWeighted->SetLineColor(kOrange);
                            promptSignalTagZHist_RecoWeighted->SetMarkerColor(kOrange);
                            promptSignalTagZHist_RecoWeighted->SetMarkerStyle(23);
                            promptSignalTagZHist_RecoWeighted->SetTitle("Reconstruction Efficiency Corrected");
                            promptSignalTagZHist_RecoWeighted->GetXaxis()->SetTitle("#it{z}_{T}");
                            promptSignalTagZHist_RecoWeighted->GetYaxis()->SetTitle("Corrected Entries");
                            promptSignalTagZHist_RecoWeighted->Draw("pe");
                            
                            // Plot 4: Acceptance-weighted
                            promptSignalWeightedCanvas->cd(4);
                            promptSignalTagZHist_AcceptanceWeighted->SetLineColor(kMagenta);
                            promptSignalTagZHist_AcceptanceWeighted->SetMarkerColor(kMagenta);
                            promptSignalTagZHist_AcceptanceWeighted->SetMarkerStyle(24);
                            promptSignalTagZHist_AcceptanceWeighted->SetTitle("Acceptance Corrected");
                            promptSignalTagZHist_AcceptanceWeighted->GetXaxis()->SetTitle("#it{z}_{T}");
                            promptSignalTagZHist_AcceptanceWeighted->GetYaxis()->SetTitle("Corrected Entries");
                            promptSignalTagZHist_AcceptanceWeighted->Draw("pe");
                            
                            // Plot 5: Fully weighted
                            promptSignalWeightedCanvas->cd(5);
                            promptSignalTagZHist_FullyWeighted->SetLineColor(kRed);
                            promptSignalTagZHist_FullyWeighted->SetMarkerColor(kRed);
                            promptSignalTagZHist_FullyWeighted->SetMarkerStyle(25);
                            promptSignalTagZHist_FullyWeighted->SetTitle("Fully Efficiency Corrected");
                            promptSignalTagZHist_FullyWeighted->GetXaxis()->SetTitle("#it{z}_{T}");
                            promptSignalTagZHist_FullyWeighted->GetYaxis()->SetTitle("Corrected Entries");
                            promptSignalTagZHist_FullyWeighted->Draw("pe");
                            
                            // Plot 6: Normalized comparison
                            promptSignalWeightedCanvas->cd(6);
                            TH1D* promptOrigNorm = (TH1D*)promptSignalTagZHistSafe->Clone("promptOrigNorm");
                            TH1D* promptPIDNorm = (TH1D*)promptSignalTagZHist_PIDWeighted->Clone("promptPIDNorm");
                            TH1D* promptFullyNorm = (TH1D*)promptSignalTagZHist_FullyWeighted->Clone("promptFullyNorm");
                            
                            if (promptOrigNorm->Integral() > 0) promptOrigNorm->Scale(1.0 / promptOrigNorm->Integral());
                            if (promptPIDNorm->Integral() > 0) promptPIDNorm->Scale(1.0 / promptPIDNorm->Integral());
                            if (promptFullyNorm->Integral() > 0) promptFullyNorm->Scale(1.0 / promptFullyNorm->Integral());
                            
                            double maxYWeighted = std::max({promptOrigNorm->GetMaximum(), promptPIDNorm->GetMaximum(), promptFullyNorm->GetMaximum()});
                            promptOrigNorm->GetYaxis()->SetRangeUser(0, maxYWeighted * 1.2);
                            promptOrigNorm->SetTitle("Efficiency Correction Comparison (Normalized)");
                            promptOrigNorm->GetYaxis()->SetTitle("Normalized Entries");
                            promptOrigNorm->Draw("pe");
                            promptPIDNorm->Draw("pe same");
                            promptFullyNorm->Draw("pe same");
                            
                            TLegend* weightedLegend = new TLegend(0.6, 0.7, 0.9, 0.9);
                            weightedLegend->SetBorderSize(0);
                            weightedLegend->SetFillStyle(0);
                            weightedLegend->AddEntry(promptOrigNorm, "Original", "pe");
                            weightedLegend->AddEntry(promptPIDNorm, "PID Corrected", "pe");
                            weightedLegend->AddEntry(promptFullyNorm, "Fully Corrected", "pe");
                            weightedLegend->Draw();
                            
                            // Save efficiency-weighted prompt signal plot
                            std::string promptWeightedOutputFile = outfilePath + "promptSignalEfficiencyWeighted_bin" + std::to_string(iBin) + ".png";
                            promptSignalWeightedCanvas->SaveAs(promptWeightedOutputFile.c_str());
                            
                            std::cout << "  Efficiency-weighted distributions created and saved" << std::endl;

                            // Also save the inputs used to make the efficiency-weighted plot into a ROOT file for later re-plotting
                            std::string promptInputsRoot = outfilePath + "promptSignalEfficiencyWeighted_bin" + std::to_string(iBin) + ".root";
                            TFile* promptInputsFile = TFile::Open(promptInputsRoot.c_str(), "RECREATE");
                            if (promptInputsFile && promptInputsFile->IsOpen()) {
                                // Write safe original and weighted histograms
                                if (promptSignalTagZHistSafe) promptSignalTagZHistSafe->Write("promptSignalTagZHistSafe");
                                if (backgroundSubtractedTagZHistSafe) backgroundSubtractedTagZHistSafe->Write("backgroundSubtractedTagZHistSafe");
                                if (promptSignalTagZHist_PIDWeighted) promptSignalTagZHist_PIDWeighted->Write("promptSignalTagZHist_PIDWeighted");
                                if (promptSignalTagZHist_RecoWeighted) promptSignalTagZHist_RecoWeighted->Write("promptSignalTagZHist_RecoWeighted");
                                if (promptSignalTagZHist_AcceptanceWeighted) promptSignalTagZHist_AcceptanceWeighted->Write("promptSignalTagZHist_AcceptanceWeighted");
                                if (promptSignalTagZHist_FullyWeighted) promptSignalTagZHist_FullyWeighted->Write("promptSignalTagZHist_FullyWeighted");

                                // Also save the normalized clones used for the comparison plot
                                if (promptOrigNorm) promptOrigNorm->Write("promptOrigNorm");
                                if (promptPIDNorm) promptPIDNorm->Write("promptPIDNorm");
                                if (promptFullyNorm) promptFullyNorm->Write("promptFullyNorm");

                                promptInputsFile->Close();
                                std::cout << "  Saved prompt-signal weighted inputs to ROOT file: " << promptInputsRoot << std::endl;
                            } else {
                                std::cerr << "  Error: could not create prompt inputs ROOT file: " << promptInputsRoot << std::endl;
                            }
                            if (promptInputsFile) { delete promptInputsFile; promptInputsFile = nullptr; }
                            
                            // Clean up efficiency correction objects
                            delete tagZEffCorrCanvas;
                            delete promptSignalWeightedCanvas;
                            delete effCorrLegend;
                            delete weightedLegend;
                            
                            // Clean up efficiency-weighted histograms
                            delete promptSignalTagZHist_PIDWeighted;
                            delete promptSignalTagZHist_RecoWeighted;
                            delete promptSignalTagZHist_AcceptanceWeighted;
                            delete promptSignalTagZHist_FullyWeighted;
                            
                            // Clean up normalized histograms
                            delete promptOrigNorm;
                            delete promptPIDNorm;
                            delete promptFullyNorm;
                            
                            // Clean up correction histograms and graphs
                            delete tagZKaonCorrHist;
                            delete tagZPionCorrHist;
                            delete tagZRecoEffCorrHist;
                            delete tagZAcceptanceCorrHist;
                            delete tagZCombinedCorrHist;
                            delete tagZKaonCountHist;
                            delete tagZPionCountHist;
                            delete tagZRecoEffCountHist;
                            delete tagZAcceptanceCountHist;
                            delete tagZCombinedCountHist;
                            delete tagZKaonCorrGraph;
                            delete tagZPionCorrGraph;
                            delete tagZRecoEffCorrGraph;
                            delete tagZAcceptanceCorrGraph;
                            delete tagZCombinedCorrGraph;
                            
                            // Clean up
                            delete promptTagZCanvas;
                            delete tagZHistNorm;
                            delete bkgSubNorm;
                            delete promptNorm;
                            delete legend;
                            delete backgroundSubtractedTagZHist;
                            delete promptSignalTagZHist;
                            delete tagZHistSafe;  // Clean up the safe copies
                            delete backgroundSubtractedTagZHistSafe;
                            delete promptSignalTagZHistSafe;
                            
                        } else {
                            std::cout << "    Warning: Failed to create prompt signal tagZ distribution" << std::endl;
                        }
                    } else {
                        std::cout << "    Warning: Could not read IP chi2 sPlot tree from file" << std::endl;
                    }
                } else {
                    std::cout << "  Warning: Failed to save IP chi2 sPlot weights to file" << std::endl;
                }
            
            // Add debug output for sPlot workflow verification
            if (enableSPlot) {
                std::cout << "\n=== sPlot Workflow Debug Information ===" << std::endl;
                std::cout << "Bin " << iBin << " sPlot workflow status:" << std::endl;
                std::cout << "  Mass fit completed: " << (massFitHisto ? "YES" : "NO") << std::endl;
                std::cout << "  Signal-enhanced data created: " << (signalEnhancedData ? "YES" : "NO") << std::endl;
                std::cout << "  IP chi2 fit completed: " << (ipChi2Histo ? "YES" : "NO") << std::endl;
                std::cout << "  IP chi2 model available: " << (ipChi2Model ? "YES" : "NO") << std::endl;
                
                // Safe check for yield variables without accessing their values
                bool yieldVarsValid = (sig_yieldLim != nullptr) && (prompt_frac != nullptr);
                std::cout << "  IP chi2 yield variables available: " << (yieldVarsValid ? "YES" : "NO") << std::endl;
                
                // Use stored values from IP chi2 fit results instead of accessing potentially invalid pointers
                if (enableSPlot && splotFile && ipChi2Model && yieldVarsValid) {
                    // These values were already calculated and stored in Stage 3
                    std::cout << "  Final total yield: " << FitIPRes_SYield[iBin][0] << std::endl;
                    std::cout << "  Final prompt fraction: " << FitIPRes_PromptFrac[iBin][0] << std::endl;
                    std::cout << "  Final prompt yield: " << FitIPRes_PromptYield[iBin][0] << std::endl;
                    std::cout << "  Final nonprompt yield: " << FitIPRes_NonPromptYield[iBin][0] << std::endl;
                }
                std::cout << "=======================================" << std::endl;
            }
            
            // Clean up
            // if (enableSPlot && signalEnhancedData != dataBin) {
            //     delete signalEnhancedData;  // Only delete if it's a different dataset
            // }
            // if (bkgLeft) delete bkgLeft;
            // if (bkgRight) delete bkgRight;
            // if (bkgData) delete bkgData;
            
        } catch (const std::exception& e) {
            std::cerr << "ERROR in bin fitting: " << e.what() << std::endl;
            // Emergency cleanup in case of exception
            // These pointers might be null, so we need to check
        }
        
        // Clean up
        // delete dataBin;
        // delete tagZHist;
    }
    
    // Finalize sPlot output file
    if (enableSPlot && splotFile) {
        std::cout << "\nFinalizing sPlot output file..." << std::endl;
        
        // Create summary TTree with bin information
        splotFile->cd();
        TTree* summaryTree = new TTree("splotSummary", "sPlot Analysis Summary");
        std::cout << "Creating summary tree with bin information..." << std::endl;
        int binIndex;
        double binCenter, binWidth;
        int nEntries;
        std::cout << "Adding branches to summary tree..." << std::endl;
        summaryTree->Branch("binIndex", &binIndex, "binIndex/I");
        summaryTree->Branch("binCenter", &binCenter, "binCenter/D");
        summaryTree->Branch("binWidth", &binWidth, "binWidth/D");
        summaryTree->Branch("nEntries", &nEntries, "nEntries/I");
        std::cout << "Filling summary tree with bin data..." << std::endl;
        for (int i = 0; i < nzTBins; ++i) {
            binIndex = i;
            binCenter = binCenters[i];
            binWidth = binWidths[i];
            
            // Count entries in this bin
            std::string zCut;
            if (zTObservable) {
                zCut = "tagZ >= " + std::to_string(zBins[i]) + 
                       " && tagZ < " + std::to_string(zBins[i+1]);
            } else {
                zCut = "tagY >= " + std::to_string(zBins[i]) + 
                       " && tagY < " + std::to_string(zBins[i+1]);
            }
            std::cout << "  Counting entries for bin " << i << ": " << zCut << std::endl;
            RooDataSet* binData = static_cast<RooDataSet*>(dataMaster->reduce(zCut.c_str()));
            nEntries = binData ? binData->numEntries() : 0;
            if (binData) delete binData;
            std::cout << "  Bin " << i << ": Center = " << binCenter 
                      << ", Width = " << binWidth 
                      << ", Entries = " << nEntries << std::endl;
            summaryTree->Fill();
            std::cout << "  Filled summary tree for bin " << i << std::endl;
        }
        
        summaryTree->Write();
        std::cout << "Summary tree written to sPlot file" << std::endl;
        
        // Close the sPlot file properly
        if (splotFile) {
            splotFile->Close();
            delete splotFile;
            splotFile = nullptr;
        }
        
        std::cout << "sPlot analysis completed and saved to file" << std::endl;
    }
    
    // Save results to file
    saveResultsToFile(binCenters, binWidths);
}

// Helper method to save results to ROOT files
void FitSpectraObject::saveResultsToFile(const std::vector<double>& binCenters, 
                                       const std::vector<double>& binWidths) {
    std::cout << "\n===== DEBUG: Starting saveResultsToFile =====" << std::endl;
    std::cout << "DEBUG: nzTBins = " << nzTBins << std::endl;
    std::cout << "DEBUG: binCenters.size() = " << binCenters.size() << std::endl;
    std::cout << "DEBUG: binWidths.size() = " << binWidths.size() << std::endl;
    
    // Check if sizes match
    if (static_cast<int>(binCenters.size()) != nzTBins || static_cast<int>(binWidths.size()) != nzTBins) {
        std::cerr << "ERROR: Size mismatch - nzTBins=" << nzTBins 
                  << ", binCenters.size()=" << binCenters.size() 
                  << ", binWidths.size()=" << binWidths.size() << std::endl;
        return;
    }
    
    std::cout << "DEBUG: Creating arrays for graph creation..." << std::endl;
    
    // Create arrays for graph creation
    std::vector<double> promptYield(nzTBins);
    std::vector<double> nonpromptYield(nzTBins);
    std::vector<double> totalYield(nzTBins);
    std::vector<double> zeros(nzTBins, 0.0);
    
    std::cout << "DEBUG: Checking result array sizes..." << std::endl;
    std::cout << "DEBUG: FitIPRes_PromptYield.size() = " << FitIPRes_PromptYield.size() << std::endl;
    std::cout << "DEBUG: FitIPRes_NonPromptYield.size() = " << FitIPRes_NonPromptYield.size() << std::endl;
    std::cout << "DEBUG: FitMRes_SYield.size() = " << FitMRes_SYield.size() << std::endl;
    
    // Check if result arrays are properly initialized
    if (FitIPRes_PromptYield.size() != static_cast<size_t>(nzTBins)) {
        std::cerr << "ERROR: FitIPRes_PromptYield size mismatch!" << std::endl;
        return;
    }
    if (FitIPRes_NonPromptYield.size() != static_cast<size_t>(nzTBins)) {
        std::cerr << "ERROR: FitIPRes_NonPromptYield size mismatch!" << std::endl;
        return;
    }
    if (FitMRes_SYield.size() != static_cast<size_t>(nzTBins)) {
        std::cerr << "ERROR: FitMRes_SYield size mismatch!" << std::endl;
        return;
    }
    
    std::cout << "DEBUG: Filling yield arrays..." << std::endl;
    for (int i = 0; i < nzTBins; ++i) {
        std::cout << "DEBUG: Processing bin " << i << std::endl;
        
        // Check if sub-arrays are properly sized
        if (FitIPRes_PromptYield[i].size() == 0) {
            std::cerr << "WARNING: FitIPRes_PromptYield[" << i << "] is empty, setting to 0" << std::endl;
            promptYield[i] = 0.0;
        } else {
            promptYield[i] = FitIPRes_PromptYield[i][0];
        }
        
        if (FitIPRes_NonPromptYield[i].size() == 0) {
            std::cerr << "WARNING: FitIPRes_NonPromptYield[" << i << "] is empty, setting to 0" << std::endl;
            nonpromptYield[i] = 0.0;
        } else {
            nonpromptYield[i] = FitIPRes_NonPromptYield[i][0];
        }
        
        if (FitMRes_SYield[i].size() == 0) {
            std::cerr << "WARNING: FitMRes_SYield[" << i << "] is empty, setting to 0" << std::endl;
            totalYield[i] = 0.0;
        } else {
            totalYield[i] = FitMRes_SYield[i][0];
        }
        
        std::cout << "DEBUG: Bin " << i << " yields - prompt: " << promptYield[i] 
                  << ", nonprompt: " << nonpromptYield[i] << ", total: " << totalYield[i] << std::endl;
    }
    
    std::cout << "DEBUG: Creating fragmentation function graphs..." << std::endl;
    
    try {
        TGraphErrors* inclFragFunc = new TGraphErrors(
            nzTBins, binCenters.data(), totalYield.data(),
            binWidths.data(), zeros.data()
        );
        inclFragFunc->SetName("ginclFragFunc");
        std::cout << "DEBUG: Created inclFragFunc graph successfully" << std::endl;
        
        TGraphErrors* promptFragFunc = new TGraphErrors(
            nzTBins, binCenters.data(), promptYield.data(),
            binWidths.data(), zeros.data()
        );
        promptFragFunc->SetName("promptFragFunc");
        std::cout << "DEBUG: Created promptFragFunc graph successfully" << std::endl;
        
        TGraphErrors* nonpromptFragFunc = new TGraphErrors(
            nzTBins, binCenters.data(), nonpromptYield.data(),
            binWidths.data(), zeros.data()
        );
        nonpromptFragFunc->SetName("nonpromptFragFunc");
        std::cout << "DEBUG: Created nonpromptFragFunc graph successfully" << std::endl;
        
        // Create parameter graphs for mass fit using optimized approach
        std::cout << "DEBUG: Creating mass parameter graphs..." << std::endl;
        std::map<std::string, std::vector<TGraphErrors*>> massGraphs;
        std::map<std::string, TGraphAsymmErrors*> massRangeGraphs;
        
        // Define mass parameter mappings to reduce repetitive code
        std::map<std::string, std::pair<std::string, std::vector<std::vector<float>>*>> massParams = {
            {"SYield", {"FitMSYield", &FitMRes_SYield}},
            {"BYield", {"FitMBYield", &FitMRes_BYield}},
            {"Mean", {"FitMMean", &FitMRes_Mean}},
            {"Sig1", {"FitMSig1", &FitMRes_Sig1}},
            {"deltaSig", {"FitMDeltaSig", &FitMRes_deltaSig}},
            {"Sig2", {"FitMSig2", &FitMRes_Sig2}},
            {"alpha", {"FitMAlpha", &FitMRes_alpha}},
            {"n", {"FitMN", &FitMRes_n}},
            {"CBFrac", {"FitMCBFrac", &FitMRes_DGFrac}},
            {"pol1", {"FitMPol1", &FitMRes_pol1}},
            {"pol2", {"FitMPol2", &FitMRes_pol2}},
            {"SYieldLim", {"FitMSYieldLim", &FitMRes_SYieldLim}},
            {"BYieldLim", {"FitMBYieldLim", &FitMRes_BYieldLim}},
            {"SYieldSG", {"FitMSYieldSG", &FitMRes_SYieldSG}},
            {"SYieldDCB", {"FitMSYieldDCB", &FitMRes_SYieldDCB}},
            {"TagZMean", {"BinTagZMean", &Bin_TagZMean}},
            {"TagZMean_weighted", {"BinTagZMean_weighted", &Bin_TagZMean_weighted}}
        };
        
        // Create graphs and range graphs in one loop
        try {
            for (const auto& [key, paramInfo] : massParams) {
                const auto& [graphName, dataPtr] = paramInfo;
                std::cout << "DEBUG: Creating " << key << " graphs..." << std::endl;
                massGraphs[key] = createParameterGraphs(graphName, nzTBins, binCenters, *dataPtr, binWidths);
                
                // Create range graphs for selected parameters
                if (key == "Mean" || key == "Sig1" || key == "deltaSig" || key == "Sig2" || 
                    key == "alpha" || key == "n" || key == "CBFrac") {
                    massRangeGraphs[key + "R"] = createGraphsAsymmErr(graphName + "Range", nzTBins, binCenters, *dataPtr, binWidths);
                }
            }
            std::cout << "DEBUG: Created all mass parameter graphs successfully" << std::endl;
            
        } catch (const std::exception& e) {
            std::cerr << "ERROR: Exception in createParameterGraphs: " << e.what() << std::endl;
            return;
        }        
        // Create parameter graphs for IP chi2 fit using optimized approach
        std::cout << "DEBUG: Creating IP chi2 parameter graphs..." << std::endl;
        std::map<std::string, std::vector<TGraphErrors*>> ipchi2Graphs;
        std::map<std::string, TGraphAsymmErrors*> ipchi2RangeGraphs;
        
        // Define IP chi2 parameter mappings
        std::map<std::string, std::pair<std::string, std::vector<std::vector<float>>*>> ipchi2Params = {
            {"SYield", {"FitIPSYield", &FitIPRes_SYield}},
            {"XpPrompt", {"FitIPXpPrompt", &FitIPRes_XpPrompt}},
            {"SigmaPrompt", {"FitIPSigmaPrompt", &FitIPRes_SigmaPrompt}},
            {"XiPrompt", {"FitIPXiPrompt", &FitIPRes_XiPrompt}},
            {"XpNonprompt", {"FitIPXpNonprompt", &FitIPRes_XpNonprompt}},
            {"SigmaNonprompt", {"FitIPSigmaNonprompt", &FitIPRes_SigmaNonprompt}},
            {"XiNonprompt", {"FitIPXiNonprompt", &FitIPRes_XiNonprompt}},
            {"PromptFrac", {"FitIPPromptFrac", &FitIPRes_PromptFrac}},
            {"PromptYield", {"FitIPPromptYield", &FitIPRes_PromptYield}},
            {"NonPromptYield", {"FitIPNonPromptYield", &FitIPRes_NonPromptYield}}
        };
        
        try {
            for (const auto& [key, paramInfo] : ipchi2Params) {
                const auto& [graphName, dataPtr] = paramInfo;
                std::cout << "DEBUG: Creating IP chi2 " << key << " graphs..." << std::endl;
                ipchi2Graphs[key] = createParameterGraphs(graphName, nzTBins, binCenters, *dataPtr, binWidths);
                
                // Create range graphs for selected parameters
                if (key.find("Xp") != std::string::npos || key.find("Sigma") != std::string::npos || 
                    key.find("Xi") != std::string::npos || key == "PromptFrac") {
                    ipchi2RangeGraphs[key + "R"] = createGraphsAsymmErr(graphName + "Range", nzTBins, binCenters, *dataPtr, binWidths);
                }
            }
            std::cout << "DEBUG: Created all IP chi2 parameter graphs successfully" << std::endl;
            
        } catch (const std::exception& e) {
            std::cerr << "ERROR: Exception in createParameterGraphs for IP chi2: " << e.what() << std::endl;
            return;
        }        
        // Save to parameter file using optimized writing functions
        std::cout << "DEBUG: Creating parameter output file: " << fOutDataName << ".root" << std::endl;
        TFile* fOutData = new TFile((fOutDataName + ".root").c_str(), "RECREATE");
        
        if (!fOutData || fOutData->IsZombie()) {
            std::cerr << "ERROR: Cannot create parameter output file: " << fOutDataName << ".root" << std::endl;
            return;
        }
        
        // Save fragmentation functions
        inclFragFunc->Write();
        promptFragFunc->Write();
        nonpromptFragFunc->Write();
        
        std::cout << "DEBUG: Writing all parameter graphs..." << std::endl;
        writeGraphsToFile(fOutData, massGraphs);
        writeGraphsToFile(fOutData, ipchi2Graphs);
        writeAsymmGraphsToFile(fOutData, massRangeGraphs);
        writeAsymmGraphsToFile(fOutData, ipchi2RangeGraphs);
        
        std::cout << "DEBUG: Closing parameter file..." << std::endl;
        fOutData->Close();
        
        // Save to histogram file
        std::cout << "DEBUG: Creating histogram output file: " << fOutDataNameB << ".root" << std::endl;
        TFile* fOutHisto = new TFile((fOutDataNameB + ".root").c_str(), "RECREATE");
        
        if (!fOutHisto || fOutHisto->IsZombie()) {
            std::cerr << "ERROR: Cannot create histogram output file: " << fOutDataNameB << ".root" << std::endl;
            return;
        }
        
        std::cout << "DEBUG: Writing fragmentation functions to histogram file..." << std::endl;
        // Save fragmentation functions again
        inclFragFunc->Write();
        promptFragFunc->Write();
        nonpromptFragFunc->Write();
        
        std::cout << "DEBUG: Writing histograms..." << std::endl;
        // Save histograms with robust error checking and validation
        for (int i = 0; i < nzTBins; ++i) {
            std::cout << "DEBUG: Processing histograms for bin " << i << std::endl;
            
            // Write IP chi2 histogram with comprehensive safety checks
            try {
                if (i < static_cast<int>(ipchi2HistoArray.size()) && 
                    ipchi2HistoArray[i] != nullptr) {
                    
                    // Additional validation - check if histogram is still valid
                    TH1* ipchi2Hist = ipchi2HistoArray[i];
                    if (ipchi2Hist->GetEntries() >= 0 && ipchi2Hist->GetNbinsX() > 0) {
                        // Create a safe copy to avoid potential pointer issues
                        TH1* ipchi2HistCopy = (TH1*)ipchi2Hist->Clone(("hIPChi2_" + std::to_string(i)).c_str());
                        if (ipchi2HistCopy) {
                            ipchi2HistCopy->SetDirectory(fOutHisto);  // Set directory first
                            ipchi2HistCopy->Write();
                            std::cout << "DEBUG: Wrote IPChi2 histogram for bin " << i << std::endl;
                            delete ipchi2HistCopy;  // Clean up the copy
                        } else {
                            std::cout << "WARNING: Failed to create IPChi2 histogram copy for bin " << i << std::endl;
                        }
                    } else {
                        std::cout << "WARNING: Invalid IPChi2 histogram data for bin " << i << std::endl;
                    }
                } else {
                    std::cout << "WARNING: Null or out-of-bounds IPChi2 histogram for bin " << i << std::endl;
                }
            } catch (const std::exception& e) {
                std::cerr << "ERROR: Exception writing IPChi2 histogram for bin " << i << ": " << e.what() << std::endl;
            } catch (...) {
                std::cerr << "ERROR: Unknown exception writing IPChi2 histogram for bin " << i << std::endl;
            }
            
            // Write mass histogram with comprehensive safety checks
            std::cout << "DEBUG: Writing mass histogram for bin " << i << std::endl;
            try {
                if (i < static_cast<int>(massHistoArray.size()) && 
                    massHistoArray[i] != nullptr) {
                    
                    // Additional validation - check if histogram is still valid
                    TH1* massHist = massHistoArray[i];
                    if (massHist->GetEntries() >= 0 && massHist->GetNbinsX() > 0) {
                        // Create a safe copy to avoid potential pointer issues
                        TH1* massHistCopy = (TH1*)massHist->Clone(("hMassSpectr_" + std::to_string(i)).c_str());
                        if (massHistCopy) {
                            massHistCopy->SetDirectory(fOutHisto);  // Set directory first
                            massHistCopy->Write();
                            std::cout << "DEBUG: Wrote mass histogram for bin " << i << std::endl;
                            delete massHistCopy;  // Clean up the copy
                        } else {
                            std::cout << "WARNING: Failed to create mass histogram copy for bin " << i << std::endl;
                        }
                    } else {
                        std::cout << "WARNING: Invalid mass histogram data for bin " << i << std::endl;
                    }
                } else {
                    std::cout << "WARNING: Null or out-of-bounds mass histogram for bin " << i << std::endl;
                }
            } catch (const std::exception& e) {
                std::cerr << "ERROR: Exception writing mass histogram for bin " << i << ": " << e.what() << std::endl;
            } catch (...) {
                std::cerr << "ERROR: Unknown exception writing mass histogram for bin " << i << std::endl;
            }
        }        
        std::cout << "DEBUG: Closing histogram file..." << std::endl;
        fOutHisto->Close();
                
        std::cout << "DEBUG: saveResultsToFile completed successfully" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "ERROR: Exception in saveResultsToFile: " << e.what() << std::endl;
        return;
    }
    
    std::cout << "===== DEBUG: saveResultsToFile completed =====" << std::endl;
    
    std::cout << "Results saved to:" << std::endl;
    std::cout << "  Parameters: " << fOutDataName << ".root" << std::endl;
    std::cout << "  Histograms: " << fOutDataNameB << ".root" << std::endl;
}


// Add these implementations after the class declaration but before any other function implementations

// Optimized implementation using templates
std::vector<TGraphErrors*> FitSpectraObject::createParameterGraphs(
    const std::string& key, int nBins, const std::vector<double>& xPos, 
    const std::vector<std::vector<float>>& yPosArr, const std::vector<double>& xWidth) {
    return createParameterGraphsTemplate(key, yPosArr, xPos, xWidth);
}

TGraphAsymmErrors* FitSpectraObject::createGraphsAsymmErr(
    const std::string& key, int nBins, const std::vector<double>& xPos,
    const std::vector<std::vector<float>>& yPosArr, const std::vector<double>& xWidth) {
    return createAsymmGraphTemplate(key, yPosArr, xPos, xWidth);
}

void FitSpectraObject::createCorrectionFactorGraphs(const std::vector<double>& corrValueKaon, 
                                                   const std::vector<double>& corrValueErrKaon,
                                                   const std::vector<double>& corrValuePion, 
                                                   const std::vector<double>& corrValueErrPion,
                                                   const std::vector<double>& corrValueRecoEff, 
                                                   const std::vector<double>& corrValueErrRecoEff,
                                                   const std::vector<double>& corrValueAcceptance, 
                                                   const std::vector<double>& corrValueErrAcceptance,
                                                   const std::vector<double>& corrValueCombinedPID, 
                                                   const std::vector<double>& corrValueErrCombinedPID,
                                                   const std::vector<double>& corrValueCombinedAll, 
                                                   const std::vector<double>& corrValueErrCombinedAll) {
    
    std::cout << "Creating correction factor graphs..." << std::endl;
    
    // Validate input - all vectors must have same size
    int nBins = corrValueKaon.size();
    std::vector<const std::vector<double>*> allVectors = {
        &corrValueKaon, &corrValueErrKaon, &corrValuePion, &corrValueErrPion,
        &corrValueRecoEff, &corrValueErrRecoEff, &corrValueAcceptance, &corrValueErrAcceptance,
        &corrValueCombinedPID, &corrValueErrCombinedPID, &corrValueCombinedAll, &corrValueErrCombinedAll
    };
    
    for (const auto* vec : allVectors) {
        if (static_cast<int>(vec->size()) != nBins) {
            std::cerr << "Error: Inconsistent correction factor vector sizes" << std::endl;
            return;
        }
    }
    
    if (nBins == 0) return;
    
    // Prepare data structures
    struct CorrectionData {
        std::string name, title;
        const std::vector<double>* values;
        const std::vector<double>* errors;
        int color, markerStyle;
    };
    
    std::vector<CorrectionData> corrections = {
        {"Kaon", "Kaon PID Efficiency Correction", &corrValueKaon, &corrValueErrKaon, kBlue, 20},
        {"Pion", "Pion PID Efficiency Correction", &corrValuePion, &corrValueErrPion, kRed, 21},
        {"RecoEff", "Reconstruction Efficiency Correction", &corrValueRecoEff, &corrValueErrRecoEff, kGreen+2, 22},
        {"Acceptance", "Acceptance Correction", &corrValueAcceptance, &corrValueErrAcceptance, kMagenta, 23},
        {"CombinedPID", "Combined PID Efficiency Correction", &corrValueCombinedPID, &corrValueErrCombinedPID, kCyan+2, 24},
        {"CombinedAll", "Combined All Corrections", &corrValueCombinedAll, &corrValueErrCombinedAll, kOrange+2, 25}
    };
    
    // Filter valid data points
    std::vector<double> validBinCenters, validBinErrors;
    std::vector<std::vector<double>> validCorrections(corrections.size()), validErrors(corrections.size());
    
    for (int i = 0; i < nBins; ++i) {
        // Check if bin is in valid range and has valid corrections
        bool isValid = (i < nzTBins - 1 && i < static_cast<int>(zBins.size()) - 1);
        if (isValid) {
            for (size_t j = 0; j < corrections.size(); ++j) {
                double val = (*corrections[j].values)[i];
                double err = (*corrections[j].errors)[i];
                if (!(val > 0 && std::isfinite(val) && err >= 0 && std::isfinite(err))) {
                    isValid = false;
                    break;
                }
            }
        }
        
        if (isValid) {
            double binCenter = (zBins[i] + zBins[i + 1]) / 2.0;
            double binError = (zBins[i + 1] - zBins[i]) / 2.0;
            validBinCenters.push_back(binCenter);
            validBinErrors.push_back(binError);
            
            for (size_t j = 0; j < corrections.size(); ++j) {
                validCorrections[j].push_back((*corrections[j].values)[i]);
                validErrors[j].push_back((*corrections[j].errors)[i]);
            }
        }
    }
    
    int nValidBins = validBinCenters.size();
    if (nValidBins == 0) {
        std::cerr << "Error: No valid correction factors found" << std::endl;
        return;
    }
    
    // Create graphs
    std::vector<TGraphErrors*> graphs;
    for (size_t i = 0; i < corrections.size(); ++i) {
        TGraphErrors* graph = new TGraphErrors(nValidBins, validBinCenters.data(), validCorrections[i].data(),
                                              validBinErrors.data(), validErrors[i].data());
        graph->SetMarkerStyle(corrections[i].markerStyle);
        graph->SetMarkerColor(corrections[i].color);
        graph->SetLineColor(corrections[i].color);
        graph->SetMarkerSize(1.2);
        graph->SetTitle(corrections[i].title.c_str());
        graphs.push_back(graph);
    }
    
    // Create combined plot
    TCanvas* canvas = createStandardCanvas("canvasCorrectionFactors", "All Efficiency Correction Factors", 1400, 1000);
    
    // Set axis and range
    std::string xAxisLabel = zTObservable ? "z_{T}" : "y";
    graphs[0]->GetXaxis()->SetTitle(xAxisLabel.c_str());
    graphs[0]->GetYaxis()->SetTitle("Correction Factor");
    graphs[0]->GetYaxis()->SetTitleOffset(1.2);
    
    // Find y-range
    double yMin = 0.8, yMax = 1.2;
    for (int i = 0; i < nValidBins; ++i) {
        for (size_t j = 0; j < corrections.size(); ++j) {
            yMin = std::min(yMin, validCorrections[j][i] - validErrors[j][i]);
            yMax = std::max(yMax, validCorrections[j][i] + validErrors[j][i]);
        }
    }
    graphs[0]->GetYaxis()->SetRangeUser(yMin * 0.9, yMax * 1.1);
    
    // Draw all graphs
    for (size_t i = 0; i < graphs.size(); ++i) {
        graphs[i]->Draw(i == 0 ? "APE" : "PE same");
    }
    
    // Add unity line and legend
    if (nValidBins > 0) {
        TLine* unityLine = new TLine(validBinCenters[0] - validBinErrors[0], 1.0, 
                                    validBinCenters[nValidBins-1] + validBinErrors[nValidBins-1], 1.0);
        unityLine->SetLineStyle(2);
        unityLine->SetLineColor(kBlack);
        unityLine->SetLineWidth(2);
        unityLine->Draw("same");
        
        TLegend* legend = createStandardLegend(0.15, 0.65, 0.45, 0.92);
        for (size_t i = 0; i < graphs.size(); ++i) {
            legend->AddEntry(graphs[i], corrections[i].name.c_str(), "pe");
        }
        legend->AddEntry(unityLine, "Unity", "l");
        legend->Draw();
    }
    
    // Save plots
    std::string combinedFileName = outfilePath + "CorrectionFactors_Combined.png";
    canvas->SaveAs(combinedFileName.c_str());
    
    // Create individual plots
    for (size_t i = 0; i < graphs.size(); ++i) {
        createEfficiencyPlot(validBinCenters, validBinErrors, validCorrections[i], validErrors[i], 
                           corrections[i].title, outfilePath + "CorrectionFactors_" + corrections[i].name + ".png",
                           corrections[i].color);
    }
    
    // Save to ROOT file
    std::string rootFileName = std::string(isMC ? "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_MC" :
                                                "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA") + 
                              "/CorrectionFactors.root";
    TFile* rootFile = new TFile(rootFileName.c_str(), "UPDATE");
    if (rootFile && rootFile->IsOpen()) {
        std::string ptString = std::to_string(static_cast<int>(jetPt.first)) + "_" + std::to_string(static_cast<int>(jetPt.second));
        for (size_t i = 0; i < graphs.size(); ++i) {
            graphs[i]->Write(("graph" + corrections[i].name + "Correction_" + ptString).c_str());
        }
        rootFile->Close();
    }
    delete rootFile;
    
    // Clean up
    delete canvas;
    for (auto* graph : graphs) delete graph;
    
    std::cout << "Correction factor graphs creation completed. Used " << nValidBins 
              << " valid bins out of " << nBins << " total bins." << std::endl;
}

void MassFitter(TString inputFile = "", bool isMC = false, bool isFitSingleBin = false, bool isZtObservable = false, bool enableSPlot = true)
{

    std::string mcTag = isMC ? "MC" : "";
    std::string obsTag = isZtObservable ? "zT" : "Y";
    
    // Open the input ROOT file
    TFile *file = TFile::Open(inputFile);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << inputFile << std::endl;
        return;
    }

    // Get the tree from the file
    TTree *tree = nullptr;
    file->GetObject("FragmNtuple", tree);
    if (!tree) {
        std::cerr << "Error getting tree from file" << std::endl;
        file->Close();
        return;
    }
    
    if (isFitSingleBin) {
        // Single bin configuration
        std::pair<double, double> jetPt(5, 60);
        std::vector<double> zBins = {0.2, 0.5, 0.65, 0.75, 0.85, 0.95, 1.0};  // D0 zT bins
        
        // Create and run fit object
        FitSpectraObject fitter(
            jetPt, isMC, zBins, 
            isZtObservable,
            tree,  // Pass the tree
            enableSPlot  // Pass the sPlot flag
        );
        fitter.startFitting();
    } 
    else {
        // Multi-bin configuration
        // std::vector<double> zBins = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};  // D0 zT bins
        std::vector<double> zBins = {
            0.0, 0.05, 0.1, 0.15, 0.2, 
            0.25, 0.3, 0.35, 0.4, 0.45, 
            0.5, 0.55, 0.6, 0.65, 0.7, 
            0.75, 0.8, 0.85, 0.9, 0.95, 
            1.0}; // D0 zT bins
        // there are 
        //LHCb y bins (rapidity)
        // std::vector<double> yBins = {2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0};
        std::vector<double> yBins = {2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0};
        
        // pT binning for jets
        std::vector<double> startPt = {5, 10, 15, 20, 30};
        std::vector<double> endPt = {10, 15, 20, 30, 50};

        //print jet pt bins 
        std::cout << "Fitting the following jet pT bins:" << std::endl;
        for (size_t i = 0; i < startPt.size(); ++i) {
            std::cout << "  Bin " << i << ": " << startPt[i] 
                      << " to " << endPt[i] << std::endl;
        }
        
        // Process each pT bin
        for (size_t jetBin = 0; jetBin < startPt.size(); ++jetBin) {
            std::pair<double, double> jetPt(startPt[jetBin], endPt[jetBin]);
            
            // Choose bin array based on observable type
            const std::vector<double>& binArray = isZtObservable ? zBins : yBins;
            
            // Create and run fit object
            FitSpectraObject fitter(
                jetPt, isMC, binArray,
                isZtObservable,
                tree,  // Pass the tree
                enableSPlot  // Pass the sPlot flag
            );
            fitter.startFitting();
        }
    }

}