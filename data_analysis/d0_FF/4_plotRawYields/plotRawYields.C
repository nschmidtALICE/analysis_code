#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <iomanip>

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
#include "TKey.h"
#include "TList.h"
#include "TPaveText.h"
#include "TLine.h"

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
    
    // TagZ-dependent correction factors (maps with bin number as key)
    std::map<int, TGraphErrors*> gTagZKaonCorr;
    std::map<int, TGraphErrors*> gTagZPionCorr;
    std::map<int, TGraphErrors*> gTagZCombinedCorr;
    
    // TagZ histograms (maps with bin number as key)
    std::map<int, TH1D*> hTagZRaw;
    std::map<int, TH1D*> hTagZBackgroundSubtracted;
    std::map<int, TH1D*> hTagZPromptSignal;
    
    // Weighted TagZ histograms (maps with bin number as key)
    std::map<int, TH1D*> hTagZAcceptanceWeighted;
    std::map<int, TH1D*> hTagZRecoWeighted;
    std::map<int, TH1D*> hTagZFullyWeighted;
    
    // PID-corrected histograms (maps with bin number as key)
    std::map<int, TH1D*> hTagZKaonCorrected;
    std::map<int, TH1D*> hTagZPionCorrected;
    std::map<int, TH1D*> hTagZCombinedCorrected;
    
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
    void plotCombinedCorrectionDemo();
    void plotTagZWeightedHistogramsDemo();
    TGraphErrors* plotNonPromptFraction(bool isCorrected = false);
    void plotYieldResult(bool isCorrected = false);
    void plotYieldResultWithTagZCorrections(bool isCorrected = false, int targetBin = -1);
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
    TGraphErrors* applyTagZCorrection(TGraphErrors* graph, const std::string& correctionType, int bin);
    TGraphErrors* applyAllCorrections(TGraphErrors* graph, const std::vector<int>& binIndices,
                                      const std::string& tagZCorrectionType,
                                      TGraphErrors* regularCorrection = nullptr);
    TGraphErrors* multiplyGraphs(TGraphErrors* graph1, TGraphErrors* graph2);
    TGraphErrors* invertGraphs(TGraphErrors* graph1);
    double propagateError(double factorA, double factorB, double factorAErr, double factorBErr, int Type = 0);
    void setOptions();
    TGraphErrors* create_empty_graph();
    TGraphErrors* create_constant_graph(double value, const std::string& name, int n_points = 10);
    void createPNPFractions(TGraphErrors* hyield, TGraphErrors* hyieldVar1, 
                           TGraphErrors* hyieldVar2, TGraphErrors* hNonPromptFrac, 
                           bool Absolute);
    void loadTagZCorrectionFactors(const std::string& basepath, bool isMC);
    TGraphErrors* getTagZCorrectionFactor(const std::string& type, int bin);
    
    // Function to apply tagZ-dependent correction to a yield point
    double applyTagZCorrection(double tagZ, double yield, double yieldError, 
                              const std::string& correctionType, int bin, 
                              double& correctedYieldError);
    
    // Function to apply tagZ corrections to a TGraphErrors
    TGraphErrors* applyTagZCorrectionToGraph(TGraphErrors* inputGraph, 
                                             const std::string& correctionType,
                                             const std::vector<int>& binIndices);
    
    // Function to apply PID corrections to tagZ histograms
    TH1D* applyPIDCorrectionToHistogram(TH1D* inputHist, const std::string& correctionType, int bin);
    
    // Function to demonstrate the usage of tagZ corrections
    void demonstrateTagZCorrections();
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
    
    // Initialize tagZ correction factor maps
    gTagZKaonCorr.clear();
    gTagZPionCorr.clear();
    gTagZCombinedCorr.clear();
    
    // Processed graphs
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

    std::string rootFileNameCorrections = basepath + "/CorrectionFactors.root";

    TFile* fInFileHistoCorrections = nullptr;
    
    try {
        fInFileHistoCorrections = TFile::Open(rootFileNameCorrections.c_str());
    } catch (const std::exception& e) {
        std::cerr << "Error opening corrections file: " << e.what() << std::endl;
        fInFileHistoCorrections = nullptr;
    }
    std::cout << "Opened corrections file: " << rootFileNameCorrections << std::endl;
    
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
            // Load the actual correction graphs created by MassFitter.C with pT range in name
            std::string kaonGraphName = "graphKaonCorrection_" + ptString;
            std::string pionGraphName = "graphPionCorrection_" + ptString;
            std::string recoEffGraphName = "graphRecoEffCorrection_" + ptString;
            std::string acceptanceGraphName = "graphAcceptanceCorrection_" + ptString;
            std::string combinedPIDGraphName = "graphCombinedPIDCorrection_" + ptString;
            std::string combinedAllGraphName = "graphCombinedAllCorrection_" + ptString;
            
            gAccCorr0 = dynamic_cast<TGraphErrors*>(fInFileHistoCorrections->Get(kaonGraphName.c_str()));
            if (!is_valid_graph(gAccCorr0)) {
                std::cout << "Warning: " << kaonGraphName << " is not a valid graph, creating empty graph" << std::endl;
                gAccCorr0 = create_empty_graph();
            }
            std::cout << "Loaded kaon correction graph: " << kaonGraphName << std::endl;
            
            gAccCorr1 = dynamic_cast<TGraphErrors*>(fInFileHistoCorrections->Get(pionGraphName.c_str()));
            if (!is_valid_graph(gAccCorr1)) {
                std::cout << "Warning: " << pionGraphName << " is not a valid graph, creating empty graph" << std::endl;
                gAccCorr1 = create_empty_graph();
            }
            std::cout << "Loaded pion correction graph: " << pionGraphName << std::endl;
            
            gAccCorr2 = dynamic_cast<TGraphErrors*>(fInFileHistoCorrections->Get(combinedPIDGraphName.c_str()));
            if (!is_valid_graph(gAccCorr2)) {
                std::cout << "Warning: " << combinedPIDGraphName << " is not a valid graph, creating empty graph" << std::endl;
                gAccCorr2 = create_empty_graph();
            }
            gAccCorr3 = dynamic_cast<TGraphErrors*>(fInFileHistoCorrections->Get(recoEffGraphName.c_str()));
            if (!is_valid_graph(gAccCorr3)) {
                std::cout << "Warning: " << recoEffGraphName << " is not a valid graph, creating empty graph" << std::endl;
                gAccCorr3 = create_empty_graph();
            }
            gAccCorr4 = dynamic_cast<TGraphErrors*>(fInFileHistoCorrections->Get(acceptanceGraphName.c_str()));
            if (!is_valid_graph(gAccCorr4)) {
                std::cout << "Warning: " << acceptanceGraphName << " is not a valid graph, creating empty graph" << std::endl;
                gAccCorr4 = create_empty_graph();
            }
            gAccCorr5 = dynamic_cast<TGraphErrors*>(fInFileHistoCorrections->Get(combinedAllGraphName.c_str()));
            if (!is_valid_graph(gAccCorr5)) {
                std::cout << "Warning: " << combinedAllGraphName << " is not a valid graph, creating empty graph" << std::endl;
                gAccCorr5 = create_empty_graph();
            }
            // Create empty graphs for unused correction factors (backward compatibility)
            // gAccCorr3 = create_empty_graph();
            // gAccCorr4 = create_empty_graph();
            // gAccCorr5 = create_empty_graph();
            gAccCorr6 = create_empty_graph();
            
            // Load tagZ-dependent correction factors
            loadTagZCorrectionFactors(basepath, isMC);
        } catch (const std::exception& e) {
            std::cerr << "Error loading correction histograms: " << e.what() << std::endl;
            // Create empty graphs as fallbacks
            gAccCorr0 = create_empty_graph();  // Kaon correction
            gAccCorr1 = create_empty_graph();  // Pion correction
            gAccCorr2 = create_empty_graph();  // Combined correction
            gAccCorr3 = create_empty_graph();
            gAccCorr4 = create_empty_graph();
            gAccCorr5 = create_empty_graph();
            gAccCorr6 = create_empty_graph();
        }
        
        // Calculate total combined factor using the combined correction graph
        gAccCorr = gAccCorr2;  // Use the combined correction directly
        
        // Convert factor to efficiency - invert all correction graphs
        // gAccCorr = invertGraphs(gAccCorr);
        // gAccCorr0 = invertGraphs(gAccCorr0);  // Kaon
        // gAccCorr1 = invertGraphs(gAccCorr1);  // Pion
        // gAccCorr2 = invertGraphs(gAccCorr2);  // Combined
        // gAccCorr3 = invertGraphs(gAccCorr3);
        // gAccCorr4 = invertGraphs(gAccCorr4);
        // gAccCorr5 = invertGraphs(gAccCorr5);
        // gAccCorr6 = invertGraphs(gAccCorr6);
        
        gAccCorr5_Ext = nullptr;
        // if (ptRange == "40_60" && resonance == "Psi2S") {
        //     std::cout << "extrapolate missing points" << std::endl;
        //     gAccCorr5_Ext = extrapolateIfNecessary(gAccCorr5);
        // }
        
        // For total acceptance correction, use the combined PID correction
        // (since gAccCorr2 is the combined kaon+pion correction)
        gAccCorrTotal = gAccCorr5;  // Start with combined PID correction
        
        // // Multiply with other correction factors if they exist and are not empty
        // if (is_valid_graph(gAccCorr3) && gAccCorr3->GetN() > 0) {
        //     gAccCorrTotal = multiplyGraphs(gAccCorrTotal, gAccCorr3);
        // }
        // if (is_valid_graph(gAccCorr4) && gAccCorr4->GetN() > 0) {
        //     gAccCorrTotal = multiplyGraphs(gAccCorrTotal, gAccCorr4);
        // }
        
        // if (gAccCorr5_Ext && is_valid_graph(gAccCorr5_Ext)) {
        //     gAccCorrTotal = multiplyGraphs(gAccCorrTotal, gAccCorr5_Ext);
        // } else if (is_valid_graph(gAccCorr5) && gAccCorr5->GetN() > 0) {
        //     gAccCorrTotal = multiplyGraphs(gAccCorrTotal, gAccCorr5);
        // }
        
        // if (is_valid_graph(gAccCorr6) && gAccCorr6->GetN() > 0) {
        //     gAccCorrTotal = multiplyGraphs(gAccCorrTotal, gAccCorr6);
        // }
        
        // Create prompt and non-prompt fractions
        createPNPFractions(hMYield, hMYieldSG, hMYieldFix, hIPNonPromptFraction, true);
        
        // Load tagZ-dependent correction factors
        loadTagZCorrectionFactors(basepath, isMC);
        
        // Demonstrate available tagZ corrections (optional - can be commented out for production)
        demonstrateTagZCorrections();
        
        // Close the file when done
        fInFileHisto->Close();
        delete fInFileHisto;
        
        // Close corrections file if it was opened
        if (fInFileHistoCorrections) {
            fInFileHistoCorrections->Close();
            delete fInFileHistoCorrections;
        }
    } else {
        std::cout << "o File: " << rootFileName << " does not exist or is invalid, creating fallback graphs" << std::endl;
        
        // Create fallback graphs for required histograms
        hMYield = create_constant_graph(100.0, "FitMSYieldF_fallback");
        hMYieldSG = create_constant_graph(80.0, "FitMSYieldSGF_fallback");
        hMYieldFix = create_constant_graph(90.0, "FitMSYieldFixF_fallback");
        hIPNonPromptFraction = create_constant_graph(0.5, "FitIPPromptFracF_default");
        
        // Create empty graphs for correction factors
        gAccCorr0 = create_empty_graph();  // Kaon correction
        gAccCorr1 = create_empty_graph();  // Pion correction
        gAccCorr2 = create_empty_graph();  // Combined correction
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
        
        // Close corrections file if it was opened
        if (fInFileHistoCorrections) {
            fInFileHistoCorrections->Close();
            delete fInFileHistoCorrections;
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
    // Clean up tagZ correction factor graphs
    for (auto& pair : gTagZKaonCorr) {
        if (pair.second) {
            delete pair.second;
        }
    }
    gTagZKaonCorr.clear();
    
    for (auto& pair : gTagZPionCorr) {
        if (pair.second) {
            delete pair.second;
        }
    }
    gTagZPionCorr.clear();
    
    for (auto& pair : gTagZCombinedCorr) {
        if (pair.second) {
            delete pair.second;
        }
    }
    gTagZCombinedCorr.clear();
    
    // Clean up tagZ histograms
    for (auto& pair : hTagZRaw) {
        if (pair.second) {
            delete pair.second;
        }
    }
    hTagZRaw.clear();
    
    for (auto& pair : hTagZBackgroundSubtracted) {
        if (pair.second) {
            delete pair.second;
        }
    }
    hTagZBackgroundSubtracted.clear();
    
    for (auto& pair : hTagZPromptSignal) {
        if (pair.second) {
            delete pair.second;
        }
    }
    hTagZPromptSignal.clear();
    
    // Clean up weighted tagZ histograms
    for (auto& pair : hTagZAcceptanceWeighted) {
        if (pair.second) {
            delete pair.second;
        }
    }
    hTagZAcceptanceWeighted.clear();
    
    for (auto& pair : hTagZRecoWeighted) {
        if (pair.second) {
            delete pair.second;
        }
    }
    hTagZRecoWeighted.clear();
    
    for (auto& pair : hTagZFullyWeighted) {
        if (pair.second) {
            delete pair.second;
        }
    }
    hTagZFullyWeighted.clear();
    
    // Clean up PID-corrected histograms
    for (auto& pair : hTagZKaonCorrected) {
        if (pair.second) {
            delete pair.second;
        }
    }
    hTagZKaonCorrected.clear();
    
    for (auto& pair : hTagZPionCorrected) {
        if (pair.second) {
            delete pair.second;
        }
    }
    hTagZPionCorrected.clear();
    
    for (auto& pair : hTagZCombinedCorrected) {
        if (pair.second) {
            delete pair.second;
        }
    }
    hTagZCombinedCorrected.clear();
    
    // Note: In a real application, you'd need to manage ownership of
    // the TGraphErrors objects more carefully. In ROOT, objects are often
    // owned by the file they're read from or written to.
}

// Implementation of plotPt method
TGraphErrors* PlotGraphsObject::plotPt(bool isCorrected) {
    plotYieldResult(isCorrected);
    plotCorrFacAcceptance();
    // plotCombinedCorrectionDemo();
    plotTagZWeightedHistogramsDemo();
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
    std::cout << "Plotting Acceptance Correction Factors for " << obsTag << " in pT range: " << ptString << std::endl;
    std::string outputFilename = OutfilePath + "FinFig_AccCorrFactor_" + obsTag + "_" + ptString + ".png";
    
    TCanvas* c = new TCanvas("c", "c: hist", 500*2, 450*2);
    c->cd();
    TGaxis::SetMaxDigits(3);
    
    // Set pad and histo arrangement
    TPad* myPad2 = new TPad("myPad", "The pad", 0, 0, 1, 1);
    myPad2->SetLeftMargin(0.15);
    myPad2->SetTopMargin(0.01);
    myPad2->SetRightMargin(0.02);
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
    myBlankHisto2->GetYaxis()->SetRangeUser(0, 1.02);
    myBlankHisto2->Draw("E");
    
    double MarkerScale = 1.6;
    //print gAccCorr0 - Kaon PID correction
    gAccCorr0->Print();
    
    setHisto(gAccCorr0, MarkerScale, kBlue); // Kaon PID correction
    gAccCorr0->Draw("same EP");
    
    setHisto(gAccCorr1, MarkerScale, kRed); // Pion PID correction
    gAccCorr1->Draw("same EP");
    
    setHisto(gAccCorr2, MarkerScale, kGreen+2); // Combined PID correction
    gAccCorr2->Draw("same EP");
    
    // Only draw additional correction factors if they are valid and non-empty
    if (is_valid_graph(gAccCorr3) && gAccCorr3->GetN() > 0) {
        setHisto(gAccCorr3, MarkerScale, kRed-7);
        gAccCorr3->Draw("same EP");
    }
    
    if (is_valid_graph(gAccCorr4) && gAccCorr4->GetN() > 0) {
        setHisto(gAccCorr4, MarkerScale, kBlue-9);
        gAccCorr4->Draw("same EP");
    }
    
    if (gAccCorr5_Ext) {
        setHisto(gAccCorr5_Ext, MarkerScale, kBlue-4);
        gAccCorr5_Ext->SetMarkerStyle(24);
        gAccCorr5_Ext->Draw("same EP");
    } else if (is_valid_graph(gAccCorr5) && gAccCorr5->GetN() > 0) {
        setHisto(gAccCorr5, MarkerScale, kBlue-4);
        gAccCorr5->Draw("same EP");
    }
    
    if (is_valid_graph(gAccCorr6) && gAccCorr6->GetN() > 0) {
        setHisto(gAccCorr6, MarkerScale, kBlue+2);
        gAccCorr6->Draw("same EP");
    }
    
    setHisto(gAccCorrTotal, MarkerScale, kMagenta+2);
    gAccCorrTotal->Draw("same EP");
    
    // Create legend
    TLegend* myLegend0 = new TLegend(0.18, 0.35, 0.45, 0.58);
    myLegend0->SetTextFont(42);
    myLegend0->SetBorderSize(0);
    myLegend0->SetFillStyle(0);
    myLegend0->SetFillColor(0);
    myLegend0->SetMargin(0.25);
    myLegend0->SetTextSize(0.04);

    myLegend0->AddEntry(gAccCorr0, "Kaon PID correction", "lep");
    myLegend0->AddEntry(gAccCorr1, "Pion PID correction", "lep");
    myLegend0->AddEntry(gAccCorr2, "Combined PID correction", "lep");
    if (is_valid_graph(gAccCorr3) && gAccCorr3->GetN() > 0)
        myLegend0->AddEntry(gAccCorr3, "Reco efficiency correction", "lep");
    if (is_valid_graph(gAccCorr4) && gAccCorr4->GetN() > 0)
        myLegend0->AddEntry(gAccCorr4, "Acceptance correction", "lep");
    myLegend0->AddEntry(gAccCorrTotal, "Total correction", "lep");

    myLegend0->Draw();
    
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

// Implementation of applyTagZCorrection method
TGraphErrors* PlotGraphsObject::applyTagZCorrection(TGraphErrors* graph, const std::string& correctionType, int bin) {
    std::cout << "Applying tagZ-dependent " << correctionType << " correction for bin " << bin << std::endl;
    
    if (!is_valid_graph(graph)) {
        std::cout << "Invalid input graph for tagZ correction" << std::endl;
        return create_empty_graph();
    }
    
    // Get the appropriate tagZ correction factor graph
    TGraphErrors* tagZCorrGraph = getTagZCorrectionFactor(correctionType, bin);
    if (!tagZCorrGraph || tagZCorrGraph->GetN() == 0) {
        std::cout << "Warning: No tagZ " << correctionType << " correction available for bin " << bin 
                  << ", returning original graph" << std::endl;
        return (TGraphErrors*)graph->Clone();
    }
    
    std::cout << "Found tagZ " << correctionType << " correction graph with " 
              << tagZCorrGraph->GetN() << " points" << std::endl;
    
    int nBinsGraph = graph->GetN();
    
    // Create arrays for the new corrected graph
    double* x = new double[nBinsGraph];
    double* xErr = new double[nBinsGraph];
    double* correctedYield = new double[nBinsGraph];
    double* correctedYieldErr = new double[nBinsGraph];
    
    // Apply correction to each point
    for (int i = 0; i < nBinsGraph; i++) {
        double xVal, yVal;
        graph->GetPoint(i, xVal, yVal);
        xErr[i] = graph->GetErrorX(i);
        double yErr = graph->GetErrorY(i);
        x[i] = xVal;
        
        // Find the correction factor from the tagZ correction graph
        // We'll interpolate if the exact x value is not found
        double corrFactor = 1.0;  // Default to no correction
        double corrFactorErr = 0.0;
        
        // First try to find exact match or interpolate
        bool foundCorrection = false;
        for (int j = 0; j < tagZCorrGraph->GetN(); j++) {
            double tagZX, tagZCorr;
            tagZCorrGraph->GetPoint(j, tagZX, tagZCorr);
            double tagZCorrErr = tagZCorrGraph->GetErrorY(j);
            
            // If we find an exact match or very close match
            if (std::abs(tagZX - xVal) < 0.001) {
                corrFactor = tagZCorr;
                corrFactorErr = tagZCorrErr;
                foundCorrection = true;
                break;
            }
        }
        
        // If no exact match found, use linear interpolation
        if (!foundCorrection && tagZCorrGraph->GetN() >= 2) {
            for (int j = 0; j < tagZCorrGraph->GetN() - 1; j++) {
                double x1, y1, x2, y2;
                tagZCorrGraph->GetPoint(j, x1, y1);
                tagZCorrGraph->GetPoint(j + 1, x2, y2);
                
                if (xVal >= x1 && xVal <= x2) {
                    // Linear interpolation
                    double fraction = (xVal - x1) / (x2 - x1);
                    corrFactor = y1 + fraction * (y2 - y1);
                    
                    // Interpolate errors as well
                    double y1Err = tagZCorrGraph->GetErrorY(j);
                    double y2Err = tagZCorrGraph->GetErrorY(j + 1);
                    corrFactorErr = y1Err + fraction * (y2Err - y1Err);
                    foundCorrection = true;
                    break;
                }
            }
        }
        
        // If still no correction found, extrapolate from nearest point
        if (!foundCorrection && tagZCorrGraph->GetN() > 0) {
            double nearestX, nearestY;
            double minDistance = 1e10;
            int nearestIndex = 0;
            
            for (int j = 0; j < tagZCorrGraph->GetN(); j++) {
                double tagZX, tagZY;
                tagZCorrGraph->GetPoint(j, tagZX, tagZY);
                double distance = std::abs(tagZX - xVal);
                if (distance < minDistance) {
                    minDistance = distance;
                    nearestIndex = j;
                    nearestX = tagZX;
                    nearestY = tagZY;
                }
            }
            
            corrFactor = nearestY;
            corrFactorErr = tagZCorrGraph->GetErrorY(nearestIndex);
            std::cout << "  Warning: Using extrapolation for x=" << xVal 
                      << " from nearest point at x=" << nearestX << std::endl;
        }
        
        // Apply correction factor
        if (corrFactor > 0) {
            correctedYield[i] = yVal / corrFactor;  // Divide by efficiency to get corrected yield
            correctedYieldErr[i] = propagateError(yVal, 1.0/corrFactor, yErr, corrFactorErr, 2);
        } else {
            correctedYield[i] = yVal;  // No correction applied
            correctedYieldErr[i] = yErr;
            std::cout << "  Warning: Invalid correction factor " << corrFactor 
                      << " for point " << i << ", using original value" << std::endl;
        }
    }
    
    // Create the corrected graph
    TGraphErrors* correctedGraph = new TGraphErrors(nBinsGraph, x, correctedYield, xErr, correctedYieldErr);
    std::string title = std::string(graph->GetTitle()) + "_tagZ" + correctionType + "Corr_bin" + std::to_string(bin);
    std::string name = std::string(graph->GetName()) + "_tagZ" + correctionType + "Corr_bin" + std::to_string(bin);
    correctedGraph->SetTitle(title.c_str());
    correctedGraph->SetName(name.c_str());
    
    // Copy style from original graph
    correctedGraph->SetMarkerStyle(graph->GetMarkerStyle());
    correctedGraph->SetMarkerColor(graph->GetMarkerColor());
    correctedGraph->SetLineColor(graph->GetLineColor());
    correctedGraph->SetMarkerSize(graph->GetMarkerSize());
    
    // Cleanup
    delete[] x;
    delete[] xErr;
    delete[] correctedYield;
    delete[] correctedYieldErr;
    
    std::cout << "TagZ " << correctionType << " correction applied successfully" << std::endl;
    return correctedGraph;
}

// Example function showing how to apply tagZ corrections in addition to regular corrections
TGraphErrors* PlotGraphsObject::applyAllCorrections(TGraphErrors* inputGraph, 
                                                   const std::vector<int>& binIndices,
                                                   const std::string& tagZCorrectionType,
                                                   TGraphErrors* regularCorrection) {
    if (!inputGraph) {
        std::cout << "Error: Input graph is null" << std::endl;
        return create_empty_graph();
    }
    
    std::cout << "\nApplying comprehensive corrections to graph: " << inputGraph->GetName() << std::endl;
    
    // Step 1: Apply regular acceptance corrections (if provided)
    TGraphErrors* correctedGraph = inputGraph;
    if (regularCorrection && is_valid_graph(regularCorrection)) {
        std::cout << "Step 1: Applying regular acceptance corrections..." << std::endl;
        correctedGraph = applyAccCorrection(inputGraph, regularCorrection);
    } else {
        std::cout << "Step 1: No regular acceptance corrections to apply" << std::endl;
        // Create a copy of the input graph
        correctedGraph = (TGraphErrors*)inputGraph->Clone((std::string(inputGraph->GetName()) + "_copy").c_str());
    }
    
    // Step 2: Apply tagZ-dependent corrections
    if (gTagZKaonCorr.size() > 0 || gTagZPionCorr.size() > 0 || gTagZCombinedCorr.size() > 0) {
        std::cout << "Step 2: Applying tagZ-dependent " << tagZCorrectionType << " corrections..." << std::endl;
        TGraphErrors* tagZCorrectedGraph = applyTagZCorrectionToGraph(correctedGraph, tagZCorrectionType, binIndices);
        
        // Clean up intermediate graph if it's not the original input
        if (correctedGraph != inputGraph) {
            delete correctedGraph;
        }
        correctedGraph = tagZCorrectedGraph;
    } else {
        std::cout << "Step 2: No tagZ corrections available" << std::endl;
    }
    
    // Set final graph properties
    std::string finalName = std::string(inputGraph->GetName()) + "_" + tagZCorrectionType + "TagZCorr";
    std::string finalTitle = std::string(inputGraph->GetTitle()) + " (All Corrections Applied)";
    correctedGraph->SetName(finalName.c_str());
    correctedGraph->SetTitle(finalTitle.c_str());
    
    std::cout << "All corrections applied successfully!" << std::endl;
    return correctedGraph;
}

// Function to create a demonstration of tagZ correction usage
void PlotGraphsObject::demonstrateTagZCorrections() {
    std::cout << "\n=== TagZ Correction Demonstration ===" << std::endl;
    std::cout << "Available tagZ correction factors:" << std::endl;
    std::cout << "  Kaon corrections: " << gTagZKaonCorr.size() << " bins" << std::endl;
    std::cout << "  Pion corrections: " << gTagZPionCorr.size() << " bins" << std::endl;
    std::cout << "  Combined corrections: " << gTagZCombinedCorr.size() << " bins" << std::endl;
    
    
    // Create additional plots showing tagZ histograms if available
    if (!hTagZRaw.empty() || !hTagZBackgroundSubtracted.empty() || !hTagZPromptSignal.empty()) {
        std::cout << "\nCreating tagZ histogram demonstration plots..." << std::endl;
        std::cout << "Available tagZ histograms:" << std::endl;
        std::cout << "  Raw tagZ histograms: " << hTagZRaw.size() << " bins" << std::endl;
        std::cout << "  Background-subtracted tagZ histograms: " << hTagZBackgroundSubtracted.size() << " bins" << std::endl;
        std::cout << "  Prompt signal tagZ histograms: " << hTagZPromptSignal.size() << " bins" << std::endl;
        
        // Create canvas showing tagZ distributions for first few bins
        TCanvas* tagZHistCanvas = new TCanvas("tagZHistograms", 
                                             ("TagZ Distributions - " + ptString).c_str(), 
                                             1200, 800);
        tagZHistCanvas->Divide(2, 2);
        
        // Plot 1: Raw tagZ distributions for multiple bins
        tagZHistCanvas->cd(1);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);
        
        TLegend* legRaw = new TLegend(0.65, 0.6, 0.9, 0.9);
        legRaw->SetBorderSize(0);
        legRaw->SetFillStyle(0);
        legRaw->SetHeader("Raw TagZ");
        
        bool firstRaw = true;
        int colorIdx = 0;
        int colors[] = {kBlue, kRed, kGreen+2, kMagenta, kCyan, kOrange};
        
        for (const auto& pair : hTagZRaw) {
            if (colorIdx >= 6) break; // Limit to 6 bins for clarity
            int bin = pair.first;
            TH1D* hist = pair.second;
            if (hist) {
                hist->SetLineColor(colors[colorIdx]);
                hist->SetMarkerColor(colors[colorIdx]);
                hist->SetMarkerStyle(20 + colorIdx);
                hist->SetLineWidth(2);
                hist->SetMarkerSize(0.8);
                
                if (firstRaw) {
                    hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum() * 1.9);
                    hist->SetTitle("Raw TagZ Distributions");
                    hist->GetXaxis()->SetTitle("#it{z}_{T}");
                    hist->GetYaxis()->SetTitle("Entries");
                    hist->Draw("PE");
                    firstRaw = false;
                } else {
                    hist->Draw("PE same");
                }
                
                legRaw->AddEntry(hist, Form("Bin %d", bin), "pe");
                colorIdx++;
            }
        }
        legRaw->Draw();
        
        // Plot 2: Background-subtracted tagZ distributions
        tagZHistCanvas->cd(2);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);
        
        TLegend* legBkgSub = new TLegend(0.65, 0.6, 0.9, 0.9);
        legBkgSub->SetBorderSize(0);
        legBkgSub->SetFillStyle(0);
        legBkgSub->SetHeader("Signal (Mass sPlot)");
        
        bool firstBkgSub = true;
        colorIdx = 0;
        
        for (const auto& pair : hTagZBackgroundSubtracted) {
            if (colorIdx >= 6) break;
            int bin = pair.first;
            TH1D* hist = pair.second;
            if (hist) {
                hist->SetLineColor(colors[colorIdx]);
                hist->SetMarkerColor(colors[colorIdx]);
                hist->SetMarkerStyle(21 + colorIdx);
                hist->SetLineWidth(2);
                hist->SetMarkerSize(0.8);
                
                if (firstBkgSub) {
                    hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum() * 1.9);
                    hist->SetTitle("Background-Subtracted TagZ Distributions");
                    hist->GetXaxis()->SetTitle("#it{z}_{T}");
                    hist->GetYaxis()->SetTitle("Weighted Entries");
                    hist->Draw("PE");
                    firstBkgSub = false;
                } else {
                    hist->Draw("PE same");
                }
                
                legBkgSub->AddEntry(hist, Form("Bin %d", bin), "pe");
                colorIdx++;
            }
        }
        legBkgSub->Draw();
        
        // Plot 3: Prompt signal tagZ distributions
        tagZHistCanvas->cd(3);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);
        
        TLegend* legPrompt = new TLegend(0.65, 0.6, 0.9, 0.9);
        legPrompt->SetBorderSize(0);
        legPrompt->SetFillStyle(0);
        legPrompt->SetHeader("Prompt Signal");
        
        bool firstPrompt = true;
        colorIdx = 0;
        
        for (const auto& pair : hTagZPromptSignal) {
            if (colorIdx >= 6) break;
            int bin = pair.first;
            TH1D* hist = pair.second;
            if (hist) {
                hist->SetLineColor(colors[colorIdx]);
                hist->SetMarkerColor(colors[colorIdx]);
                hist->SetMarkerStyle(22 + colorIdx);
                hist->SetLineWidth(2);
                hist->SetMarkerSize(0.8);
                
                if (firstPrompt) {
                    hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum() * 1.9);
                    hist->SetTitle("Prompt Signal TagZ Distributions");
                    hist->GetXaxis()->SetTitle("#it{z}_{T}");
                    hist->GetYaxis()->SetTitle("Weighted Entries");
                    hist->Draw("PE");
                    firstPrompt = false;
                } else {
                    hist->Draw("PE same");
                }
                
                legPrompt->AddEntry(hist, Form("Bin %d", bin), "pe");
                colorIdx++;
            }
        }
        legPrompt->Draw();
        
        // Plot 4: Comparison for a single bin (if available)
        tagZHistCanvas->cd(4);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);
        
        // Find a bin that has all three types of histograms
        int compareBin = -1;
        for (const auto& pair : hTagZRaw) {
            int bin = pair.first;
            if (hTagZBackgroundSubtracted.find(bin) != hTagZBackgroundSubtracted.end() &&
                hTagZPromptSignal.find(bin) != hTagZPromptSignal.end()) {
                compareBin = bin;
                break;
            }
        }
        compareBin+=1;
        std::cout << "Comparing histograms for bin: " << compareBin << std::endl;
        
        if (compareBin >= 0) {
            TH1D* rawHist = hTagZRaw[compareBin];
            TH1D* bkgSubHist = hTagZBackgroundSubtracted[compareBin];
            TH1D* promptHist = hTagZPromptSignal[compareBin];
            
            // Normalize for comparison
            TH1D* rawNorm = (TH1D*)rawHist->Clone("rawNorm");
            TH1D* bkgSubNorm = (TH1D*)bkgSubHist->Clone("bkgSubNorm");
            TH1D* promptNorm = (TH1D*)promptHist->Clone("promptNorm");
            
            if (rawNorm->Integral() > 0) rawNorm->Scale(1.0 / rawNorm->Integral());
            if (bkgSubNorm->Integral() > 0) bkgSubNorm->Scale(1.0 / bkgSubNorm->Integral());
            if (promptNorm->Integral() > 0) promptNorm->Scale(1.0 / promptNorm->Integral());
            
            rawNorm->SetLineColor(kBlack);
            rawNorm->SetMarkerColor(kBlack);
            rawNorm->SetMarkerStyle(20);
            rawNorm->SetLineWidth(2);
            
            bkgSubNorm->SetLineColor(kBlue);
            bkgSubNorm->SetMarkerColor(kBlue);
            bkgSubNorm->SetMarkerStyle(21);
            bkgSubNorm->SetLineWidth(2);
            
            promptNorm->SetLineColor(kRed);
            promptNorm->SetMarkerColor(kRed);
            promptNorm->SetMarkerStyle(22);
            promptNorm->SetLineWidth(2);
            
            rawNorm->SetTitle(Form("TagZ Comparison - Bin %d (Normalized)", compareBin));
            rawNorm->GetXaxis()->SetTitle("#it{z}_{T}");
            rawNorm->GetYaxis()->SetTitle("Normalized Entries");
            rawNorm->Draw("PE");
            bkgSubNorm->Draw("PE same");
            promptNorm->Draw("PE same");
            
            TLegend* legCompare = new TLegend(0.6, 0.7, 0.9, 0.9);
            legCompare->SetBorderSize(0);
            legCompare->SetFillStyle(0);
            legCompare->AddEntry(rawNorm, "Raw", "pe");
            legCompare->AddEntry(bkgSubNorm, "Signal (Mass sPlot)", "pe");
            legCompare->AddEntry(promptNorm, "Prompt Signal", "pe");
            legCompare->Draw();
            
            // delete rawNorm;
            // delete bkgSubNorm;
            // delete promptNorm;
        } else {
            TLatex* noCompareText = new TLatex(0.5, 0.5, "No complete bin data available");
            noCompareText->SetNDC();
            noCompareText->SetTextAlign(22);
            noCompareText->Draw();
        }
        
        // Save the histogram canvas
        std::string histOutputPath = OutfilePath + "TagZHistogramDemo_" + ptString + ".png";
        tagZHistCanvas->SaveAs(histOutputPath.c_str());
        std::cout << "Saved tagZ histogram demonstration plot: " << histOutputPath << std::endl;
        
        delete tagZHistCanvas;
    }
    
    // Create PID correction demonstration plots for prompt signal histograms
    if (!hTagZPromptSignal.empty() && (!gTagZKaonCorr.empty() || !gTagZPionCorr.empty() || !gTagZCombinedCorr.empty())) {
        std::cout << "\nCreating PID correction demonstration plots for prompt signal histograms..." << std::endl;
        
        // Create and store PID-corrected histograms for all available bins as member variables
        std::cout << "Creating PID-corrected histograms for all bins..." << std::endl;
        
        int createdHistograms = 0;
        
        for (const auto& pair : hTagZPromptSignal) {
            int bin = pair.first;
            TH1D* originalHist = pair.second;
            
            if (!originalHist) continue;
            
            std::cout << "  Processing bin " << bin << "..." << std::endl;
            
            // Apply and store kaon corrections if available
            if (gTagZKaonCorr.find(bin) != gTagZKaonCorr.end()) {
                TH1D* kaonCorrectedHist = applyPIDCorrectionToHistogram(originalHist, "kaon", bin);
                if (kaonCorrectedHist) {
                    hTagZKaonCorrected[bin] = kaonCorrectedHist;
                    createdHistograms++;
                    std::cout << "    Created kaon-corrected histogram for bin " << bin << std::endl;
                }
            }
            
            // Apply and store pion corrections if available
            if (gTagZPionCorr.find(bin) != gTagZPionCorr.end()) {
                TH1D* pionCorrectedHist = applyPIDCorrectionToHistogram(originalHist, "pion", bin);
                if (pionCorrectedHist) {
                    hTagZPionCorrected[bin] = pionCorrectedHist;
                    createdHistograms++;
                    std::cout << "    Created pion-corrected histogram for bin " << bin << std::endl;
                }
            }
            
            // Apply and store combined corrections if available
            if (gTagZCombinedCorr.find(bin) != gTagZCombinedCorr.end()) {
                TH1D* combinedCorrectedHist = applyPIDCorrectionToHistogram(originalHist, "combined", bin);
                if (combinedCorrectedHist) {
                    hTagZCombinedCorrected[bin] = combinedCorrectedHist;
                                       createdHistograms++;
                    std::cout << "    Created combined-corrected histogram for bin " << bin << std::endl;
                }
            }
        }
        
        std::cout << "Successfully created " << createdHistograms << " PID-corrected histograms" << std::endl;
        std::cout << "  Kaon corrected: " << hTagZKaonCorrected.size() << " bins" << std::endl;
        std::cout << "  Pion corrected: " << hTagZPionCorrected.size() << " bins" << std::endl;
        std::cout << "  Combined corrected: " << hTagZCombinedCorrected.size() << " bins" << std::endl;
        
        // Create canvas for PID correction effects
        TCanvas* pidCorrCanvas = new TCanvas("pidCorrectionEffects", 
                                           ("PID Correction Effects - " + ptString).c_str(), 
                                           1600, 1200);
        pidCorrCanvas->Divide(2, 3);
        
        // Find a representative bin that has all correction types
        int demonstrationBin = -1;
        for (const auto& pair : hTagZPromptSignal) {
            int bin = pair.first;
            if (gTagZKaonCorr.find(bin) != gTagZKaonCorr.end() &&
                gTagZPionCorr.find(bin) != gTagZPionCorr.end() &&
                gTagZCombinedCorr.find(bin) != gTagZCombinedCorr.end()) {
                demonstrationBin = bin;
                break;
            }
        }
        
        if (demonstrationBin >= 0) {
            TH1D* originalHist = hTagZPromptSignal[demonstrationBin];
            
            // Use the stored corrected histograms from member variables
            TH1D* kaonCorrectedHist = hTagZKaonCorrected[demonstrationBin];
            TH1D* pionCorrectedHist = hTagZPionCorrected[demonstrationBin];
            TH1D* combinedCorrectedHist = hTagZCombinedCorrected[demonstrationBin];
            
            // Plot 1: Original vs Kaon corrected
            pidCorrCanvas->cd(1);
            gPad->SetLeftMargin(0.12);
            gPad->SetBottomMargin(0.12);
            
            kaonCorrectedHist->SetLineColor(kBlue);
            kaonCorrectedHist->SetMarkerColor(kBlue);
            kaonCorrectedHist->SetMarkerStyle(21);
            kaonCorrectedHist->SetLineWidth(2);
            kaonCorrectedHist->SetTitle(Form("Kaon PID Correction - Bin %d", demonstrationBin));
            kaonCorrectedHist->GetXaxis()->SetTitle("#it{z}_{T}");
            kaonCorrectedHist->GetYaxis()->SetTitle("Weighted Entries");
            kaonCorrectedHist->Draw("PE");

            originalHist->SetLineColor(kBlack);
            originalHist->SetMarkerColor(kBlack);
            originalHist->SetMarkerStyle(20);
            originalHist->SetLineWidth(2);
            originalHist->Draw("PE same");
            
            
            TLegend* legKaon = new TLegend(0.6, 0.75, 0.9, 0.9);
            legKaon->SetBorderSize(0);
            legKaon->SetFillStyle(0);
            legKaon->AddEntry(originalHist, "Original", "pe");
            legKaon->AddEntry(kaonCorrectedHist, "Kaon Corrected", "pe");
            legKaon->Draw();
            
            // Plot 2: Original vs Pion corrected
            pidCorrCanvas->cd(2);
            gPad->SetLeftMargin(0.12);
            gPad->SetBottomMargin(0.12);

            pionCorrectedHist->SetTitle(Form("Pion PID Correction - Bin %d", demonstrationBin));
            pionCorrectedHist->SetLineColor(kRed);
            pionCorrectedHist->SetMarkerColor(kRed);
            pionCorrectedHist->SetMarkerStyle(22);
            pionCorrectedHist->SetLineWidth(2);
            pionCorrectedHist->Draw("PE");

            originalHist->Draw("PE same");
            
            TLegend* legPion = new TLegend(0.6, 0.75, 0.9, 0.9);
            legPion->SetBorderSize(0);
            legPion->SetFillStyle(0);
            legPion->AddEntry(originalHist, "Original", "pe");
            legPion->AddEntry(pionCorrectedHist, "Pion Corrected", "pe");
            legPion->Draw();
            
            // Plot 3: Original vs Combined corrected
            pidCorrCanvas->cd(3);
            gPad->SetLeftMargin(0.12);
            gPad->SetBottomMargin(0.12);
            
            combinedCorrectedHist->SetTitle(Form("Combined PID Correction - Bin %d", demonstrationBin));
            combinedCorrectedHist->SetLineColor(kGreen+2);
            combinedCorrectedHist->SetMarkerColor(kGreen+2);
            combinedCorrectedHist->SetMarkerStyle(23);
            combinedCorrectedHist->SetLineWidth(2);
            combinedCorrectedHist->Draw("PE");
            
            originalHist->Draw("PE same");

            TLegend* legCombined = new TLegend(0.6, 0.75, 0.9, 0.9);
            legCombined->SetBorderSize(0);
            legCombined->SetFillStyle(0);
            legCombined->AddEntry(originalHist, "Original", "pe");
            legCombined->AddEntry(combinedCorrectedHist, "Combined Corrected", "pe");
            legCombined->Draw();
            
            // Plot 4: All corrections compared
            pidCorrCanvas->cd(4);
            gPad->SetLeftMargin(0.12);
            gPad->SetBottomMargin(0.12);
            
            originalHist->SetTitle(Form("All PID Corrections Comparison - Bin %d", demonstrationBin));
            combinedCorrectedHist->Draw("PE");
            originalHist->Draw("PE same");
            kaonCorrectedHist->Draw("PE same");
            pionCorrectedHist->Draw("PE same");
            
            TLegend* legAll = new TLegend(0.55, 0.65, 0.9, 0.9);
            legAll->SetBorderSize(0);
            legAll->SetFillStyle(0);
            legAll->AddEntry(originalHist, "Original", "pe");
            legAll->AddEntry(kaonCorrectedHist, "Kaon Corrected", "pe");
            legAll->AddEntry(pionCorrectedHist, "Pion Corrected", "pe");
            legAll->AddEntry(combinedCorrectedHist, "Combined Corrected", "pe");
            legAll->Draw();
            
            // Plot 5: Correction factors as ratios
            pidCorrCanvas->cd(5);
            gPad->SetLeftMargin(0.12);
            gPad->SetBottomMargin(0.12);
            
            // Create ratio histograms
            TH1D* kaonRatio = (TH1D*)kaonCorrectedHist->Clone("kaonRatio");
            TH1D* pionRatio = (TH1D*)pionCorrectedHist->Clone("pionRatio");
            TH1D* combinedRatio = (TH1D*)combinedCorrectedHist->Clone("combinedRatio");
            
            kaonRatio->Divide(originalHist);
            pionRatio->Divide(originalHist);
            combinedRatio->Divide(originalHist);
            
            kaonRatio->SetTitle(Form("Correction Factors (Corrected/Original) - Bin %d", demonstrationBin));
            kaonRatio->GetYaxis()->SetTitle("Correction Factor");
            kaonRatio->GetYaxis()->SetRangeUser(0.7, 4.0);
            kaonRatio->Draw("PE");
            pionRatio->Draw("PE same");
            combinedRatio->Draw("PE same");
            
            // Add horizontal line at 1.0
            TLine* unityLine = new TLine(kaonRatio->GetXaxis()->GetXmin(), 1.0, 
                                        kaonRatio->GetXaxis()->GetXmax(), 1.0);
            unityLine->SetLineStyle(2);
            unityLine->SetLineColor(kGray+2);
            unityLine->Draw();
            
            TLegend* legRatio = new TLegend(0.6, 0.75, 0.9, 0.9);
            legRatio->SetBorderSize(0);
            legRatio->SetFillStyle(0);
            legRatio->AddEntry(kaonRatio, "Kaon Factor", "pe");
            legRatio->AddEntry(pionRatio, "Pion Factor", "pe");
            legRatio->AddEntry(combinedRatio, "Combined Factor", "pe");
            legRatio->Draw();
            
            // Plot 6: Summary statistics
            pidCorrCanvas->cd(6);
            gPad->SetLeftMargin(0.12);
            gPad->SetBottomMargin(0.12);
            
            // Create summary text
            TPaveText* summaryText = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
            summaryText->SetFillColor(kWhite);
            summaryText->SetBorderSize(1);
            summaryText->SetTextAlign(12);
            summaryText->SetTextFont(42);
            
            summaryText->AddText(Form("PID Correction Summary - Bin %d", demonstrationBin));
            summaryText->AddText("");
            summaryText->AddText(Form("Original integral: %.2f", originalHist->Integral()));
            summaryText->AddText(Form("Kaon corrected integral: %.2f (factor: %.3f)", 
                                     kaonCorrectedHist->Integral(),
                                     kaonCorrectedHist->Integral() / originalHist->Integral()));
            summaryText->AddText(Form("Pion corrected integral: %.2f (factor: %.3f)", 
                                     pionCorrectedHist->Integral(),
                                     pionCorrectedHist->Integral() / originalHist->Integral()));
            summaryText->AddText(Form("Combined corrected integral: %.2f (factor: %.3f)", 
                                     combinedCorrectedHist->Integral(),
                                     combinedCorrectedHist->Integral() / originalHist->Integral()));
            summaryText->AddText("");
            summaryText->AddText(Form("Original mean: %.4f", originalHist->GetMean()));
            summaryText->AddText(Form("Combined corrected mean: %.4f", combinedCorrectedHist->GetMean()));
            summaryText->AddText(Form("Mean shift: %.4f", combinedCorrectedHist->GetMean() - originalHist->GetMean()));
            
            summaryText->Draw();
            
            // Save the PID correction canvas
            std::string pidCorrOutputPath = OutfilePath + "PIDCorrectionDemo_" + ptString + ".png";
            pidCorrCanvas->SaveAs(pidCorrOutputPath.c_str());
            std::cout << "Saved PID correction demonstration plot: " << pidCorrOutputPath << std::endl;
            
            // Clean up temporary ratio histograms (corrected histograms are now stored as member variables)
            delete kaonRatio;
            delete pionRatio;
            delete combinedRatio;
            
        } else {
            std::cout << "Warning: No bin found with complete correction factor data for PID demonstration" << std::endl;
        }
        
        delete pidCorrCanvas;
    }
    
    // Show numerical correction factors for each bin
    for (const auto& pair : gTagZCombinedCorr) {
        int bin = pair.first;
        TGraphErrors* graph = pair.second;
        if (graph && graph->GetN() > 0) {
            std::cout << "\nBin " << bin << " combined correction factors:" << std::endl;
            for (int i = 0; i < std::min(5, graph->GetN()); ++i) { // Show first 5 points
                double tagZ, corrFactor;
                graph->GetPoint(i, tagZ, corrFactor);
                double corrError = graph->GetErrorY(i);
                std::cout << "  tagZ=" << std::fixed << std::setprecision(3) << tagZ 
                          << " → correction=" << std::setprecision(4) << corrFactor 
                          << "±" << corrError << std::endl;
            }
            if (graph->GetN() > 5) {
                std::cout << "  ... and " << (graph->GetN() - 5) << " more points" << std::endl;
            }
        }
    }
    
    std::cout << "=== End TagZ Correction Demonstration ===" << std::endl;
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
        // if (normType == 1) {
        //     myBlankHisto2->GetYaxis()->SetRangeUser(1, 100000);
        // }
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
                
                // Save tagZ histograms and PID-corrected histograms for each pT range
                std::cout << "Saving tagZ and PID-corrected histograms for pT range " << i << std::endl;
                
                // Save original tagZ histograms with pT range identifier
                for (const auto& pair : gObj[i]->hTagZRaw) {
                    int bin = pair.first;
                    TH1D* hist = pair.second;
                    if (hist) {
                        std::string name = "tagZRaw_" + pTRangeArray[i] + "_bin" + std::to_string(bin);
                        hist->Write(name.c_str());
                    }
                }
                
                for (const auto& pair : gObj[i]->hTagZBackgroundSubtracted) {
                    int bin = pair.first;
                    TH1D* hist = pair.second;
                    if (hist) {
                        std::string name = "tagZBackgroundSubtracted_" + pTRangeArray[i] + "_bin" + std::to_string(bin);
                        hist->Write(name.c_str());
                    }
                }
                
                for (const auto& pair : gObj[i]->hTagZPromptSignal) {
                    int bin = pair.first;
                    TH1D* hist = pair.second;
                    if (hist) {
                        std::string name = "tagZPromptSignal_" + pTRangeArray[i] + "_bin" + std::to_string(bin);
                        hist->Write(name.c_str());
                    }
                }
                
                // Save PID-corrected histograms with pT range identifier
                for (const auto& pair : gObj[i]->hTagZKaonCorrected) {
                    int bin = pair.first;
                    TH1D* hist = pair.second;
                    if (hist) {
                        std::string name = "tagZKaonCorrected_" + pTRangeArray[i] + "_bin" + std::to_string(bin);
                        hist->Write(name.c_str());
                    }
                }
                
                for (const auto& pair : gObj[i]->hTagZPionCorrected) {
                    int bin = pair.first;
                    TH1D* hist = pair.second;
                    if (hist) {
                        std::string name = "tagZPionCorrected_" + pTRangeArray[i] + "_bin" + std::to_string(bin);
                        hist->Write(name.c_str());
                    }
                }
                
                for (const auto& pair : gObj[i]->hTagZCombinedCorrected) {
                    int bin = pair.first;
                    TH1D* hist = pair.second;
                    if (hist) {
                        std::string name = "tagZCombinedCorrected_" + pTRangeArray[i] + "_bin" + std::to_string(bin);
                        hist->Write(name.c_str());
                    }
                }
            }
            
            std::cout << "Successfully saved all histograms including tagZ and PID-corrected histograms to: " << outFileName << std::endl;
            
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


// Implementation of loadTagZCorrectionFactors function
void PlotGraphsObject::loadTagZCorrectionFactors(const std::string& basepath, bool isMC) {
    std::cout << "Loading tagZ-dependent correction factors..." << std::endl;
    
    // Construct the path to the TagZCorrectionFactors.root file
    std::string tagZCorrectionFile = basepath + "/TagZCorrectionFactors.root";
    
    // Clear existing maps
    gTagZKaonCorr.clear();
    gTagZPionCorr.clear();
    gTagZCombinedCorr.clear();
    
    // Try to open the file
    TFile* tagZFile = nullptr;
    try {
        tagZFile = TFile::Open(tagZCorrectionFile.c_str(), "READ");
        if (!tagZFile || tagZFile->IsZombie()) {
            std::cout << "Warning: TagZ correction factors file not found or corrupted: " << tagZCorrectionFile << std::endl;
            std::cout << "TagZ-dependent corrections will not be applied." << std::endl;
            if (tagZFile) {
                tagZFile->Close();
                delete tagZFile;
            }
            return;
        }
        
        std::cout << "Opened tagZ correction factors file: " << tagZCorrectionFile << std::endl;
        
        // Get list of keys to find all correction factor graphs
        TList* keyList = tagZFile->GetListOfKeys();
        if (!keyList) {
            std::cout << "Warning: No keys found in tagZ correction factors file" << std::endl;
            tagZFile->Close();
            delete tagZFile;
            return;
        }
        
        int loadedGraphs = 0;
        
        // Iterate through all keys to find our correction factor graphs
        TIter keyIter(keyList);
        TKey* key;
        while ((key = (TKey*)keyIter())) {
            std::string keyName = key->GetName();
            
            // Check if this key matches our naming pattern
            if (keyName.find("tagZKaonCorrection_" + ptString + "_bin") == 0 ||
                keyName.find("tagZPionCorrection_" + ptString + "_bin") == 0 ||
                keyName.find("tagZCombinedCorrection_" + ptString + "_bin") == 0) {
                
                // Extract bin number from the key name
                size_t binPos = keyName.find("_bin");
                if (binPos != std::string::npos) {
                    std::string binStr = keyName.substr(binPos + 4);
                    int binNumber = std::stoi(binStr);
                    
                    // Load the graph
                    TGraphErrors* graph = dynamic_cast<TGraphErrors*>(tagZFile->Get(keyName.c_str()));
                    if (graph && graph->GetN() > 0) {
                        // Store in the appropriate map
                        if (keyName.find("tagZKaonCorrection_") == 0) {
                            gTagZKaonCorr[binNumber] = (TGraphErrors*)graph->Clone();
                            std::cout << "  Loaded kaon tagZ correction for bin " << binNumber 
                                      << " with " << graph->GetN() << " points" << std::endl;
                            loadedGraphs++;
                        } else if (keyName.find("tagZPionCorrection_") == 0) {
                            gTagZPionCorr[binNumber] = (TGraphErrors*)graph->Clone();
                            std::cout << "  Loaded pion tagZ correction for bin " << binNumber 
                                      << " with " << graph->GetN() << " points" << std::endl;
                            loadedGraphs++;
                        } else if (keyName.find("tagZCombinedCorrection_") == 0) {
                            gTagZCombinedCorr[binNumber] = (TGraphErrors*)graph->Clone();
                            std::cout << "  Loaded combined tagZ correction for bin " << binNumber 
                                      << " with " << graph->GetN() << " points" << std::endl;
                            loadedGraphs++;
                        }
                    } else {
                        std::cout << "  Warning: Could not load graph " << keyName << " or graph is empty" << std::endl;
                    }
                }
            }
        }
        
        std::cout << "Successfully loaded " << loadedGraphs << " tagZ correction factor graphs" << std::endl;
        std::cout << "  Kaon corrections: " << gTagZKaonCorr.size() << " bins" << std::endl;
        std::cout << "  Pion corrections: " << gTagZPionCorr.size() << " bins" << std::endl;
        std::cout << "  Combined corrections: " << gTagZCombinedCorr.size() << " bins" << std::endl;
        
        tagZFile->Close();
        delete tagZFile;
        
    } catch (const std::exception& e) {
        std::cerr << "Error loading tagZ correction factors: " << e.what() << std::endl;
        if (tagZFile) {
            tagZFile->Close();
            delete tagZFile;
        }
    }
    
    // Load tagZ histograms from TagZHistograms_*.root file
    std::string tagZHistFile = basepath + "/TagZHistograms_" + ptString + ".root";
    
    // Clear existing histogram maps
    hTagZRaw.clear();
    hTagZBackgroundSubtracted.clear();
    hTagZPromptSignal.clear();
    
    TFile* histFile = nullptr;
    try {
        histFile = TFile::Open(tagZHistFile.c_str(), "READ");
        if (!histFile || histFile->IsZombie()) {
            std::cout << "Warning: TagZ histograms file not found: " << tagZHistFile << std::endl;
            std::cout << "TagZ histogram demonstrations will be limited." << std::endl;
            if (histFile) {
                histFile->Close();
                delete histFile;
            }
        } else {
            std::cout << "Loading tagZ histograms from: " << tagZHistFile << std::endl;
            
            // Get list of keys to find all histograms
            TList* histKeyList = histFile->GetListOfKeys();
            if (histKeyList) {
                int loadedHists = 0;
                
                TIter histKeyIter(histKeyList);
                TKey* histKey;
                while ((histKey = (TKey*)histKeyIter())) {
                    std::string histKeyName = histKey->GetName();
                    
                    // Check if this key matches our naming pattern for this pT range
                    if (histKeyName.find("tagZHist_" + ptString + "_bin") == 0 ||
                        histKeyName.find("backgroundSubtractedTagZHist_" + ptString + "_bin") == 0 ||
                        histKeyName.find("promptSignalTagZHist_" + ptString + "_bin") == 0 ||
                        histKeyName.find("promptSignalTagZHist_AcceptanceWeighted_" + ptString + "_bin") == 0 ||
                        histKeyName.find("promptSignalTagZHist_RecoWeighted_" + ptString + "_bin") == 0 ||
                        histKeyName.find("promptSignalTagZHist_FullyWeighted_" + ptString + "_bin") == 0) {
                        
                        // Extract bin number from the key name
                        size_t binPos = histKeyName.find("_bin");
                        if (binPos != std::string::npos) {
                            std::string binStr = histKeyName.substr(binPos + 4);
                            int binNumber = std::stoi(binStr);
                            
                            // Load the histogram
                            TH1D* hist = dynamic_cast<TH1D*>(histFile->Get(histKeyName.c_str()));
                            if (hist) {
                                // Store in the appropriate map
                                if (histKeyName.find("tagZHist_" + ptString + "_bin") == 0) {
                                    hTagZRaw[binNumber] = (TH1D*)hist->Clone();
                                    hTagZRaw[binNumber]->SetDirectory(nullptr);
                                    std::cout << "  Loaded raw tagZ histogram for bin " << binNumber << std::endl;
                                    loadedHists++;
                                } else if (histKeyName.find("backgroundSubtractedTagZHist_" + ptString + "_bin") == 0) {
                                    hTagZBackgroundSubtracted[binNumber] = (TH1D*)hist->Clone();
                                    hTagZBackgroundSubtracted[binNumber]->SetDirectory(nullptr);
                                    std::cout << "  Loaded background-subtracted tagZ histogram for bin " << binNumber << std::endl;
                                    loadedHists++;
                                } else if (histKeyName.find("promptSignalTagZHist_" + ptString + "_bin") == 0) {
                                    hTagZPromptSignal[binNumber] = (TH1D*)hist->Clone();
                                    hTagZPromptSignal[binNumber]->SetDirectory(nullptr);
                                    std::cout << "  Loaded prompt signal tagZ histogram for bin " << binNumber << std::endl;
                                    loadedHists++;
                                } else if (histKeyName.find("promptSignalTagZHist_AcceptanceWeighted_" + ptString + "_bin") == 0) {
                                    hTagZAcceptanceWeighted[binNumber] = (TH1D*)hist->Clone();
                                    hTagZAcceptanceWeighted[binNumber]->SetDirectory(nullptr);
                                    std::cout << "  Loaded acceptance-weighted tagZ histogram for bin " << binNumber << std::endl;
                                    loadedHists++;
                                } else if (histKeyName.find("promptSignalTagZHist_RecoWeighted_" + ptString + "_bin") == 0) {
                                    hTagZRecoWeighted[binNumber] = (TH1D*)hist->Clone();
                                    hTagZRecoWeighted[binNumber]->SetDirectory(nullptr);
                                    std::cout << "  Loaded reco-weighted tagZ histogram for bin " << binNumber << std::endl;
                                    loadedHists++;
                                } else if (histKeyName.find("promptSignalTagZHist_FullyWeighted_" + ptString + "_bin") == 0) {
                                    hTagZFullyWeighted[binNumber] = (TH1D*)hist->Clone();
                                    hTagZFullyWeighted[binNumber]->SetDirectory(nullptr);
                                    std::cout << "  Loaded fully-weighted tagZ histogram for bin " << binNumber << std::endl;
                                    loadedHists++;
                                }
                            }
                        }
                    }
                }
                
                std::cout << "Successfully loaded " << loadedHists << " tagZ histograms" << std::endl;
                std::cout << "  Raw tagZ histograms: " << hTagZRaw.size() << " bins" << std::endl;
                std::cout << "  Background-subtracted tagZ histograms: " << hTagZBackgroundSubtracted.size() << " bins" << std::endl;
                std::cout << "  Prompt signal tagZ histograms: " << hTagZPromptSignal.size() << " bins" << std::endl;
                std::cout << "  Acceptance-weighted tagZ histograms: " << hTagZAcceptanceWeighted.size() << " bins" << std::endl;
                std::cout << "  Reco-weighted tagZ histograms: " << hTagZRecoWeighted.size() << " bins" << std::endl;
                std::cout << "  Fully-weighted tagZ histograms: " << hTagZFullyWeighted.size() << " bins" << std::endl;
            }
            
            histFile->Close();
            delete histFile;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error loading tagZ histograms: " << e.what() << std::endl;
        if (histFile) {
            histFile->Close();
            delete histFile;
        }
    }
}

// Implementation of getTagZCorrectionFactor function
TGraphErrors* PlotGraphsObject::getTagZCorrectionFactor(const std::string& type, int bin) {
    if (type == "kaon") {
        auto it = gTagZKaonCorr.find(bin);
        if (it != gTagZKaonCorr.end()) {
            return it->second;
        }
    } else if (type == "pion") {
        auto it = gTagZPionCorr.find(bin);
        if (it != gTagZPionCorr.end()) {
            return it->second;
        }
    } else if (type == "combined") {
        auto it = gTagZCombinedCorr.find(bin);
        if (it != gTagZCombinedCorr.end()) {
            return it->second;
        }
    }
    
    // Return nullptr if not found
    return nullptr;
}

// Implementation of applyTagZCorrectionToGraph function
TGraphErrors* PlotGraphsObject::applyTagZCorrectionToGraph(TGraphErrors* inputGraph, 
                                                         const std::string& correctionType,
                                                         const std::vector<int>& binIndices) {
    if (!inputGraph || inputGraph->GetN() == 0) {
        std::cout << "Warning: Input graph is null or empty for tagZ correction" << std::endl;
        return create_empty_graph();
    }
    
    if (binIndices.size() != static_cast<size_t>(inputGraph->GetN())) {
        std::cout << "Error: Number of bin indices (" << binIndices.size() 
                  << ") does not match number of points in graph (" << inputGraph->GetN() << ")" << std::endl;
        return create_empty_graph();
    }
    
    // Create output graph
    TGraphErrors* correctedGraph = new TGraphErrors();
    correctedGraph->SetName((std::string(inputGraph->GetName()) + "_tagZCorrected").c_str());
    correctedGraph->SetTitle((std::string(inputGraph->GetTitle()) + " (TagZ " + correctionType + " Corrected)").c_str());
    
    std::cout << "Applying tagZ " << correctionType << " corrections to graph " 
              << inputGraph->GetName() << std::endl;
    
    for (int i = 0; i < inputGraph->GetN(); ++i) {
        double tagZ, yield;
        inputGraph->GetPoint(i, tagZ, yield);
        double yieldError = inputGraph->GetErrorY(i);
        double xError = inputGraph->GetErrorX(i);
        
        int binIndex = binIndices[i];
        
        // Get the correction factor for this bin and tagZ value
        TGraphErrors* corrGraph = getTagZCorrectionFactor(correctionType, binIndex);
        double correctedYield = yield;
        double correctedYieldError = yieldError;
        
        if (corrGraph && corrGraph->GetN() > 0) {
            // Find correction factor using linear interpolation
            double correctionFactor = 1.0;
            double correctionError = 0.0;
            bool foundCorrection = false;
            
            // Linear interpolation between points
            for (int j = 0; j < corrGraph->GetN() - 1; ++j) {
                double x1, y1, x2, y2;
                corrGraph->GetPoint(j, x1, y1);
                corrGraph->GetPoint(j + 1, x2, y2);
                
                if (tagZ >= x1 && tagZ <= x2) {
                    double t = (tagZ - x1) / (x2 - x1);
                    correctionFactor = y1 + t * (y2 - y1);
                    
                    double err1 = corrGraph->GetErrorY(j);
                    double err2 = corrGraph->GetErrorY(j + 1);
                    correctionError = err1 + t * (err2 - err1);
                    
                    foundCorrection = true;
                    break;
                }
            }
            
            // If outside range, use nearest point
            if (!foundCorrection) {
                if (tagZ < corrGraph->GetX()[0]) {
                    double x_dummy;
                    corrGraph->GetPoint(0, x_dummy, correctionFactor);
                    correctionError = corrGraph->GetErrorY(0);
                } else if (tagZ > corrGraph->GetX()[corrGraph->GetN() - 1]) {
                    double x_dummy;
                    corrGraph->GetPoint(corrGraph->GetN() - 1, x_dummy, correctionFactor);
                    correctionError = corrGraph->GetErrorY(corrGraph->GetN() - 1);
                }
                foundCorrection = true;
            }
            
            if (foundCorrection && correctionFactor > 0) {
                // Apply correction (divide by efficiency to get corrected yield
                correctedYield = yield / correctionFactor;
                
                // Propagate errors
                if (yield > 0) {
                    double relYieldError = yieldError / yield;
                    double relCorrError = correctionError / correctionFactor;
                    correctedYieldError = correctedYield * sqrt(relYieldError * relYieldError + 
                                                               relCorrError * relCorrError);
                }
            }
        }
        
        correctedGraph->SetPoint(i, tagZ, correctedYield);
        correctedGraph->SetPointError(i, xError, correctedYieldError);
    }
    
    std::cout << "Applied tagZ " << correctionType << " corrections to " 
              << correctedGraph->GetN() << " points" << std::endl;
    
    return correctedGraph;
}

// Implementation of applyPIDCorrectionToHistogram function
TH1D* PlotGraphsObject::applyPIDCorrectionToHistogram(TH1D* inputHist, const std::string& correctionType, int bin) {
    if (!inputHist) {
        std::cerr << "Error: Input histogram is null" << std::endl;
        return nullptr;
    }
    
    // Get the correction factor graph for this bin
    TGraphErrors* corrGraph = getTagZCorrectionFactor(correctionType, bin);
    if (!corrGraph || corrGraph->GetN() == 0) {
        std::cout << "Warning: No " << correctionType << " correction factors available for bin " << bin 
                  << ". Returning original histogram." << std::endl;
        return (TH1D*)inputHist->Clone((std::string(inputHist->GetName()) + "_" + correctionType + "_corrected").c_str());
    }
    
    // Create corrected histogram
    TH1D* correctedHist = (TH1D*)inputHist->Clone((std::string(inputHist->GetName()) + "_" + correctionType + "_corrected").c_str());
    correctedHist->SetTitle((std::string(inputHist->GetTitle()) + " (" + correctionType + " PID corrected)").c_str());
    correctedHist->SetDirectory(nullptr);
    
    std::cout << "Applying " << correctionType << " PID corrections to histogram " << inputHist->GetName() 
              << " for bin " << bin << std::endl;
    
    int correctedBins = 0;
    
    // Apply corrections bin by bin
    for (int histBin = 1; histBin <= correctedHist->GetNbinsX(); ++histBin) {
        double tagZ = correctedHist->GetBinCenter(histBin);
        double binContent = correctedHist->GetBinContent(histBin);
        double binError = correctedHist->GetBinError(histBin);
        
        if (binContent <= 0) continue; // Skip empty bins
        
        // Find correction factor using linear interpolation
        double correctionFactor = 1.0;
        double correctionError = 0.0;
        bool foundCorrection = false;
        
        // Linear interpolation between points
        for (int j = 0; j < corrGraph->GetN() - 1; ++j) {
            double x1, y1, x2, y2;
            corrGraph->GetPoint(j, x1, y1);
            corrGraph->GetPoint(j + 1, x2, y2);
            
            if (tagZ >= x1 && tagZ <= x2) {
                double t = (tagZ - x1) / (x2 - x1);
                correctionFactor = y1 + t * (y2 - y1);
                
                double err1 = corrGraph->GetErrorY(j);
                double err2 = corrGraph->GetErrorY(j + 1);
                correctionError = err1 + t * (err2 - err1);
                
                foundCorrection = true;
                break;
            }
        }
        
        // If outside range, use nearest point
        if (!foundCorrection) {
            if (tagZ < corrGraph->GetX()[0]) {
                double x_dummy;
                corrGraph->GetPoint(0, x_dummy, correctionFactor);
                correctionError = corrGraph->GetErrorY(0);
            } else if (tagZ > corrGraph->GetX()[corrGraph->GetN() - 1]) {
                double x_dummy;
                corrGraph->GetPoint(corrGraph->GetN() - 1, x_dummy, correctionFactor);
                correctionError = corrGraph->GetErrorY(corrGraph->GetN() - 1);
            }
            foundCorrection = true;
        }
        
        if (foundCorrection && correctionFactor > 0) {
            // Apply correction (divide by efficiency to correct for detector effects)
            double correctedContent = binContent / correctionFactor;
            
            // Propagate errors
            double relContentError = (binContent > 0) ? binError / binContent : 0.0;
            double relCorrError = (correctionFactor > 0) ? correctionError / correctionFactor : 0.0;
            double correctedError = correctedContent * sqrt(relContentError * relContentError + 
                                                           relCorrError * relCorrError);
            
            correctedHist->SetBinContent(histBin, correctedContent);
            correctedHist->SetBinError(histBin, correctedError);
            correctedBins++;
        }
    }
    
    std::cout << "  Applied corrections to " << correctedBins << " bins out of " 
              << correctedHist->GetNbinsX() << " total bins" << std::endl;
    
    return correctedHist;
}

// Implementation of setOptions method
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

// Implementation of plotTagZWeightedHistogramsDemo method
void PlotGraphsObject::plotTagZWeightedHistogramsDemo() {
    setOptions();
    std::cout << "Plotting TagZ Weighted Histograms Demo for " << obsTag << " in pT range: " << ptString << std::endl;
    std::string outputFilename = OutfilePath + "FinFig_TagZWeightedHistogramsDemo_" + obsTag + "_" + ptString + ".png";

    TCanvas* c = new TCanvas("cTagZWeighted", "TagZ Weighted Histograms Demo", 800, 600);
    c->cd();
    TGaxis::SetMaxDigits(3);

    TPad* pad = new TPad("padTagZWeighted", "The pad", 0, 0, 1, 1);
    pad->SetLeftMargin(0.15);
    pad->SetTopMargin(0.05);
    pad->SetRightMargin(0.05);
    pad->SetBottomMargin(0.15);
    pad->SetTicks();
    pad->Draw();
    pad->cd();

    // Plot histograms for bin 0 as demonstration
    int demoBin = 2;
    bool plotted = false;


    if (hTagZFullyWeighted.find(demoBin) != hTagZFullyWeighted.end() && hTagZFullyWeighted[demoBin]) {
        hTagZFullyWeighted[demoBin]->SetLineColor(kGreen+2);
        hTagZFullyWeighted[demoBin]->SetMarkerColor(kGreen+2);
        hTagZFullyWeighted[demoBin]->SetMarkerStyle(20);
        hTagZFullyWeighted[demoBin]->SetMarkerSize(1.2);
        hTagZFullyWeighted[demoBin]->SetLineWidth(2);
        hTagZFullyWeighted[demoBin]->SetXTitle("tag_{Z}");
        hTagZFullyWeighted[demoBin]->SetYTitle("Weighted Counts");
        hTagZFullyWeighted[demoBin]->SetTitle("");
        hTagZFullyWeighted[demoBin]->Draw("pe");
        plotted = true;
    }

    if (hTagZAcceptanceWeighted.find(demoBin) != hTagZAcceptanceWeighted.end() && hTagZAcceptanceWeighted[demoBin]) {
        hTagZAcceptanceWeighted[demoBin]->SetLineColor(kBlue);
        hTagZAcceptanceWeighted[demoBin]->SetMarkerColor(kBlue);
        hTagZAcceptanceWeighted[demoBin]->SetMarkerStyle(20);
        hTagZAcceptanceWeighted[demoBin]->SetMarkerSize(1.2);
        hTagZAcceptanceWeighted[demoBin]->SetLineWidth(2);
        if( plotted) {
            hTagZAcceptanceWeighted[demoBin]->Draw("pe same");
        } else {
            hTagZAcceptanceWeighted[demoBin]->SetXTitle("tag_{Z}");
            hTagZAcceptanceWeighted[demoBin]->SetYTitle("Weighted Counts");
            hTagZAcceptanceWeighted[demoBin]->SetTitle("");
            hTagZAcceptanceWeighted[demoBin]->Draw("pe");
            plotted = true;
        }
    }

    if (hTagZRecoWeighted.find(demoBin) != hTagZRecoWeighted.end() && hTagZRecoWeighted[demoBin]) {
        hTagZRecoWeighted[demoBin]->SetLineColor(kRed);
        hTagZRecoWeighted[demoBin]->SetMarkerColor(kRed);
        hTagZRecoWeighted[demoBin]->SetMarkerStyle(20);
        hTagZRecoWeighted[demoBin]->SetMarkerSize(1.2);
        hTagZRecoWeighted[demoBin]->SetLineWidth(2);
        if (plotted) {
            hTagZRecoWeighted[demoBin]->Draw("pe same");
        } else {
            hTagZRecoWeighted[demoBin]->SetXTitle("tag_{Z}");
            hTagZRecoWeighted[demoBin]->SetYTitle("Weighted Counts");
            hTagZRecoWeighted[demoBin]->SetTitle("");
            hTagZRecoWeighted[demoBin]->Draw("pe");
            plotted = true;
        }
    }


    if (plotted) {
        TLegend* legend = new TLegend(0.55, 0.7, 0.88, 0.88);
        legend->SetTextFont(42);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetFillColor(0);
        legend->SetMargin(0.25);
        legend->SetTextSize(0.04);

        if (hTagZAcceptanceWeighted.find(demoBin) != hTagZAcceptanceWeighted.end() && hTagZAcceptanceWeighted[demoBin])
            legend->AddEntry(hTagZAcceptanceWeighted[demoBin], "Acceptance weighted", "l");
        if (hTagZRecoWeighted.find(demoBin) != hTagZRecoWeighted.end() && hTagZRecoWeighted[demoBin])
            legend->AddEntry(hTagZRecoWeighted[demoBin], "Reco efficiency weighted", "l");
        if (hTagZFullyWeighted.find(demoBin) != hTagZFullyWeighted.end() && hTagZFullyWeighted[demoBin])
            legend->AddEntry(hTagZFullyWeighted[demoBin], "Fully weighted", "l");

        legend->Draw();

        c->SaveAs(outputFilename.c_str());
        std::cout << "Saved TagZ weighted histograms demo plot: " << outputFilename << std::endl;
    } else {
        std::cout << "Warning: No weighted TagZ histograms found for plotting demo" << std::endl;
    }

    delete pad;
    delete c;
}
