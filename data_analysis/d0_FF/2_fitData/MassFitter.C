// MassFitterScript.cpp
// D0 Meson Analysis Framework for mass and lifetime fits

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <memory>
#include <filesystem>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

// RooFit includes
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooFit.h"
#include "RooMsgService.h"
#include "RooWorkspace.h"
#include "RooPlot.h"

// Include external fitter and plotter headers
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
    std::pair<double, double> jetPt;  // Added to store pT range
    std::vector<double> zBins;
    bool zTObservable;
    TTree* inputTree;
    std::pair<double, double> sideBandLimits = std::make_pair(1.840, 1.890);  // D0 mass sideband region
    std::string outfilePath;
    std::string fOutDataName;
    std::string fOutDataNameB;
    
    // Add these member variables
    int nzTBins;                // Number of zT bins
    int numFitItems = 5;        // Number of fit parameters to store
    bool correctionOnly = false;  // Flag for correction-only mode
    
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
    std::vector<std::vector<float>> FitMRes_pol0;
    std::vector<std::vector<float>> FitMRes_pol1;
    std::vector<std::vector<float>> FitMRes_SYieldLim;
    std::vector<std::vector<float>> FitMRes_BYieldLim;
    std::vector<std::vector<float>> FitMRes_SYieldSG;
    std::vector<std::vector<float>> FitMRes_SYieldDCB;
    
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
                    TTree* tree);  
    
    // Destructor
    ~FitSpectraObject();
    
    // Methods
    void startFitting();
    void configureFilePaths();
    void initializeResultArrays();
    Fitter* createFitter();
    RooDataSet* prepareMasterDataset(Fitter* fitter);
    void processCorrectionFactors(Fitter* fitter, RooDataSet* dataMaster);
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


// Implementation of the FitSpectraObject constructor
FitSpectraObject::FitSpectraObject(
    const std::pair<double, double>& ptRange, 
    bool isMc, 
    const std::vector<double>& zBins, 
    bool isZtObservable,
    TTree* tree)
    : isMC(isMc), jetPt(ptRange), zBins(zBins), zTObservable(isZtObservable), 
      inputTree(tree), nzTBins(zBins.size() - 1)  // Initialize nzTBins here
{
    std::cout << "Initializing FitSpectraObject with pT range: " 
              << ptRange.first << " - " << ptRange.second << " GeV/c" << std::endl;

    // Configure file paths
    configureFilePaths();
    
    // Print bin information
    std::cout << "Number of zT bins: " << nzTBins << std::endl;
    
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
    std::cout << "FitSpectraObject destroyed" << std::endl;
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
    // Process fit for each bin
    processFitsByBin(fitter, dataMaster, zBinsCent, zBinsWidth);
    
    // Clean up
    delete fitter;
    delete dataMaster;
}


void FitSpectraObject::configureFilePaths() {
    // Add MC/data and observable tags for file naming
    std::string mcTag = isMC ? "MC" : "";
    std::string obsTag = zTObservable ? "zT" : "dR";
    
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
    FitMRes_pol0.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_pol1.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_SYieldLim.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_BYieldLim.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_SYieldSG.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    FitMRes_SYieldDCB.resize(nzTBins, std::vector<float>(numFitItems, 0.0));
    
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
    // This is a placeholder implementation for the correction factors
    // Typically this would apply efficiency or acceptance corrections to the data
    
    if (!fitter || !dataMaster) {
        std::cerr << "Error: Null fitter or dataset in processCorrectionFactors" << std::endl;
        return;
    }
    
    std::cout << "Processing correction factors..." << std::endl;
    
    // For now, just print dataset info
    std::cout << "Dataset entries: " << dataMaster->numEntries() << std::endl;
    
    // Print the first few entries (for debugging)
    const RooArgSet* vars = dataMaster->get(0);
    if (vars) {
        std::cout << "Dataset variables: ";
        vars->Print("v");
    }
    
    // No actual corrections applied in this simplified implementation
    std::cout << "No corrections applied - this is a placeholder implementation" << std::endl;
    
    // If you need to implement actual corrections, you would:
    // 1. Iterate through dataset entries
    // 2. Apply correction weights based on kinematics
    // 3. Create a new dataset with weights or modify existing weights
}

void FitSpectraObject::processFitsByBin(Fitter* fitter, RooDataSet* dataMaster, 
                                      const std::vector<double>& binCenters, 
                                      const std::vector<double>& binWidths) {
    if (!fitter || !dataMaster) {
        std::cerr << "Error: Null fitter or dataset in processFitsByBin" << std::endl;
        return;
    }
    
    // Get resonance parameters from fitter
    auto resSig = fitter->getMassDict("D0");
    auto sigRegion = resSig.signalRegion;
    
    // Loop over all bins
    for (int iBin = 0; iBin < nzTBins; ++iBin) {
        std::cout << "\n==== Processing bin " << iBin 
                  << " (z: " << zBins[iBin] << "-" << zBins[iBin+1] << ") ====" << std::endl;
        
        // Create z bin selection string
        std::string zCut;
        if (zTObservable) {
            zCut = "tagZ >= " + std::to_string(zBins[iBin]) + 
                   " && tagZ < " + std::to_string(zBins[iBin+1]);
        } else {
            zCut = "tagJetdR >= " + std::to_string(zBins[iBin]) + 
                   " && tagJetdR < " + std::to_string(zBins[iBin+1]);
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
        
        // Perform mass fit
        std::string zRangeStr = std::to_string(zBins[iBin]) + "-" + std::to_string(zBins[iBin+1]);
        std::cout << "Performing mass fit..." << std::endl;
        
        TH1* massFitHisto = nullptr;
        std::vector<double> fitParams;
        std::vector<double> fitErrors;
        
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
                fitter->massFit("D0", dataBin, "DGauss", iBin, zRangeStr);
            
            // Store histogram
            if (massFitHisto) {
                massHistoArray[iBin] = massFitHisto;
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
                
                FitMRes_pol0[iBin][0] = dgParams[8];    // Pol0
                FitMRes_pol0[iBin][1] = dgErrors[8];    // Pol0 error
                
                FitMRes_pol1[iBin][0] = dgParams[9];    // Pol1
                FitMRes_pol1[iBin][1] = dgErrors[9];    // Pol1 error
                
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
            
            // Perform IP chi2 fit to separate prompt and non-prompt
            std::cout << "\nPerforming IP chi2 fit..." << std::endl;
            
            // Prepare signal and background datasets for IP chi2 fit
            RooDataSet* sigData = static_cast<RooDataSet*>(dataBin->reduce(
                ("tagMass > " + std::to_string(sigRegion.first) + 
                " && tagMass < " + std::to_string(sigRegion.second)).c_str()
            ));
            
            // Left sideband
            RooDataSet* bkgLeft = static_cast<RooDataSet*>(dataBin->reduce(
                ("tagMass < " + std::to_string(sideBandLimits.first)).c_str()
            ));
            
            // Right sideband
            RooDataSet* bkgRight = static_cast<RooDataSet*>(dataBin->reduce(
                ("tagMass > " + std::to_string(sideBandLimits.second)).c_str()
            ));
            
            // Combine sidebands
            RooDataSet* bkgData = dynamic_cast<RooDataSet*>(bkgLeft->Clone("BKG"));
            bkgData->append(*bkgRight);
            
            // Perform IP chi2 fit
            TH1* ipChi2Histo = nullptr;
            std::vector<double> ipParams, ipErrors;
            
            std::tie(ipChi2Histo, ipParams, ipErrors) = 
                fitter->ipchi2Fit("D0", sigData, bkgData, "BKGincluded", iBin, zRangeStr);
            
            // Store IP chi2 histogram
            if (ipChi2Histo) {
                ipchi2HistoArray[iBin] = ipChi2Histo;
            }
            
            // Store IP chi2 fit results if we have enough parameters
            if (ipParams.size() >= 8) {
                // Store the IP chi2 fit parameters
                FitIPRes_SYield[iBin][0] = ipParams[0];  // Total signal yield
                FitIPRes_PromptFrac[iBin][0] = ipParams[1];         // Correct position for prompt fraction
                FitIPRes_XpPrompt[iBin][0] = ipParams[2];           // Prompt peak position
                FitIPRes_SigmaPrompt[iBin][0] = ipParams[3];        // Prompt width
                FitIPRes_XiPrompt[iBin][0] = ipParams[4];           // Prompt asymmetry
                // Store rho1_prompt and rho2_prompt if needed
                FitIPRes_XpNonprompt[iBin][0] = ipParams[7];        // Non-prompt peak position
                FitIPRes_SigmaNonprompt[iBin][0] = ipParams[8];     // Non-prompt width
                FitIPRes_XiNonprompt[iBin][0] = ipParams[9];        // Non-prompt asymmetry
                
                // And calculate yields with the correct prompt fraction
                FitIPRes_PromptYield[iBin][0] = ipParams[0] * ipParams[1];       // Using correct prompt_frac
                FitIPRes_NonPromptYield[iBin][0] = ipParams[0] * (1.0 - ipParams[1]); // Using correct fraction
                
                // Store errors if available
                if (ipErrors.size() >= 10) {  // Update: check for enough errors (at least 10)
                    FitIPRes_SYield[iBin][1] = ipErrors[0];          // Signal yield error
                    FitIPRes_PromptFrac[iBin][1] = ipErrors[1];      // Prompt fraction error
                    FitIPRes_XpPrompt[iBin][1] = ipErrors[2];        // Prompt peak position error
                    FitIPRes_SigmaPrompt[iBin][1] = ipErrors[3];     // Prompt width error
                    FitIPRes_XiPrompt[iBin][1] = ipErrors[4];        // Prompt asymmetry error
                    // Skip storing rho1_prompt and rho2_prompt errors (indices 5, 6)
                    FitIPRes_XpNonprompt[iBin][1] = ipErrors[7];     // Non-prompt peak position error
                    FitIPRes_SigmaNonprompt[iBin][1] = ipErrors[8];  // Non-prompt width error
                    FitIPRes_XiNonprompt[iBin][1] = ipErrors[9];     // Non-prompt asymmetry error
                    
                    // Also propagate errors to yield calculations
                    // Using error propagation for product: σ(Y×f) = Y×f × √((σY/Y)² + (σf/f)²)
                    double relErrorSig = ipErrors[0]/std::max(1e-10, ipParams[0]);
                    double relErrorFrac = ipErrors[1]/std::max(1e-10, ipParams[1]);
                    FitIPRes_PromptYield[iBin][1] = ipParams[0] * ipParams[1] * 
                                                  sqrt(pow(relErrorSig, 2) + pow(relErrorFrac, 2));
                    
                    // For non-prompt yield: Y×(1-f)
                    // Error propagation with correlation term: σ(Y×(1-f)) = √((1-f)²σY² + Y²σf²)
                    FitIPRes_NonPromptYield[iBin][1] = sqrt(pow((1.0-ipParams[1])*ipErrors[0], 2) + 
                                                         pow(ipParams[0]*ipErrors[1], 2));
                }
                
                // Print key IP chi2 fit results
                std::cout << "IP chi2 fit results for bin " << iBin << ":" << std::endl;
                std::cout << "  Prompt fraction: " << FitIPRes_PromptFrac[iBin][0] 
                          << " ± " << FitIPRes_PromptFrac[iBin][1] << std::endl;
                std::cout << "  Prompt yield: " << FitIPRes_PromptYield[iBin][0] << std::endl;
                std::cout << "  Non-prompt yield: " << FitIPRes_NonPromptYield[iBin][0] << std::endl;
            }
            
            // Clean up
            delete sigData;
            delete bkgLeft;
            delete bkgRight;
            delete bkgData;
            
        } catch (const std::exception& e) {
            std::cerr << "ERROR in bin fitting: " << e.what() << std::endl;
        }
        
        // Clean up
        delete dataBin;
    }
    
    // Save results to file
    saveResultsToFile(binCenters, binWidths);
}

// Helper method to save results to ROOT files
void FitSpectraObject::saveResultsToFile(const std::vector<double>& binCenters, 
                                       const std::vector<double>& binWidths) {
    std::cout << "\nSaving results to file..." << std::endl;
    
    // Create arrays for graph creation
    std::vector<double> promptYield(nzTBins);
    std::vector<double> nonpromptYield(nzTBins);
    std::vector<double> totalYield(nzTBins);
    std::vector<double> zeros(nzTBins, 0.0);
    
    for (int i = 0; i < nzTBins; ++i) {
        promptYield[i] = FitIPRes_PromptYield[i][0];
        nonpromptYield[i] = FitIPRes_NonPromptYield[i][0];
        totalYield[i] = FitMRes_SYield[i][0];
    }
    
    // Create fragmentation function graphs
    TGraphErrors* inclFragFunc = new TGraphErrors(
        nzTBins, binCenters.data(), totalYield.data(),
        binWidths.data(), zeros.data()
    );
    inclFragFunc->SetName("ginclFragFunc");
    
    TGraphErrors* promptFragFunc = new TGraphErrors(
        nzTBins, binCenters.data(), promptYield.data(),
        binWidths.data(), zeros.data()
    );
    promptFragFunc->SetName("promptFragFunc");
    
    TGraphErrors* nonpromptFragFunc = new TGraphErrors(
        nzTBins, binCenters.data(), nonpromptYield.data(),
        binWidths.data(), zeros.data()
    );
    nonpromptFragFunc->SetName("nonpromptFragFunc");
    
    // Create parameter graphs for mass fit
    std::map<std::string, std::vector<TGraphErrors*>> massGraphs;
    std::map<std::string, TGraphAsymmErrors*> massRangeGraphs;
    
    // Create graphs for each mass parameter
    massGraphs["SYield"] = createParameterGraphs("FitMSYield", nzTBins, binCenters, FitMRes_SYield, binWidths);
    massGraphs["BYield"] = createParameterGraphs("FitMBYield", nzTBins, binCenters, FitMRes_BYield, binWidths);
    massGraphs["Mean"] = createParameterGraphs("FitMMean", nzTBins, binCenters, FitMRes_Mean, binWidths);
    massGraphs["Sig1"] = createParameterGraphs("FitMSig1", nzTBins, binCenters, FitMRes_Sig1, binWidths);
    massGraphs["deltaSig"] = createParameterGraphs("FitMDeltaSig", nzTBins, binCenters, FitMRes_deltaSig, binWidths);
    massGraphs["Sig2"] = createParameterGraphs("FitMSig2", nzTBins, binCenters, FitMRes_Sig2, binWidths);
    massGraphs["alpha"] = createParameterGraphs("FitMAlpha", nzTBins, binCenters, FitMRes_alpha, binWidths);
    massGraphs["n"] = createParameterGraphs("FitMN", nzTBins, binCenters, FitMRes_n, binWidths);
    massGraphs["CBFrac"] = createParameterGraphs("FitMCBFrac", nzTBins, binCenters, FitMRes_DGFrac, binWidths);
    massGraphs["pol0"] = createParameterGraphs("FitMPol0", nzTBins, binCenters, FitMRes_pol0, binWidths);
    massGraphs["pol1"] = createParameterGraphs("FitMPol1", nzTBins, binCenters, FitMRes_pol1, binWidths);
    massGraphs["SYieldLim"] = createParameterGraphs("FitMSYieldLim", nzTBins, binCenters, FitMRes_SYieldLim, binWidths);
    massGraphs["BYieldLim"] = createParameterGraphs("FitMBYieldLim", nzTBins, binCenters, FitMRes_BYieldLim, binWidths);
    massGraphs["SYieldSG"] = createParameterGraphs("FitMSYieldSG", nzTBins, binCenters, FitMRes_SYieldSG, binWidths);
    massGraphs["SYieldDCB"] = createParameterGraphs("FitMSYieldDCB", nzTBins, binCenters, FitMRes_SYieldDCB, binWidths);
    
    // Create range graphs for mass parameters
    massRangeGraphs["MeanR"] = createGraphsAsymmErr("FitMMeanRange", nzTBins, binCenters, FitMRes_Mean, binWidths);
    massRangeGraphs["Sig1R"] = createGraphsAsymmErr("FitMSig1Range", nzTBins, binCenters, FitMRes_Sig1, binWidths);
    massRangeGraphs["deltaSigR"] = createGraphsAsymmErr("FitMDeltaSigRange", nzTBins, binCenters, FitMRes_deltaSig, binWidths);
    massRangeGraphs["Sig2R"] = createGraphsAsymmErr("FitMSig2Range", nzTBins, binCenters, FitMRes_Sig2, binWidths);
    massRangeGraphs["alphaR"] = createGraphsAsymmErr("FitMAlphaRange", nzTBins, binCenters, FitMRes_alpha, binWidths);
    massRangeGraphs["nR"] = createGraphsAsymmErr("FitMNRange", nzTBins, binCenters, FitMRes_n, binWidths);
    massRangeGraphs["CBFracR"] = createGraphsAsymmErr("FitMCBFracRange", nzTBins, binCenters, FitMRes_DGFrac, binWidths);
    
    // Create parameter graphs for IP chi2 fit
    std::map<std::string, std::vector<TGraphErrors*>> ipchi2Graphs;
    std::map<std::string, TGraphAsymmErrors*> ipchi2RangeGraphs;
    
    // Create graphs for each IP chi2 parameter
    ipchi2Graphs["SYield"] = createParameterGraphs("FitIPSYield", nzTBins, binCenters, FitIPRes_SYield, binWidths);
    // Single Bukin parameters for prompt and non-prompt components
    ipchi2Graphs["XpPrompt"] = createParameterGraphs("FitIPXpPrompt", nzTBins, binCenters, FitIPRes_XpPrompt, binWidths);
    ipchi2Graphs["SigmaPrompt"] = createParameterGraphs("FitIPSigmaPrompt", nzTBins, binCenters, FitIPRes_SigmaPrompt, binWidths);
    ipchi2Graphs["XiPrompt"] = createParameterGraphs("FitIPXiPrompt", nzTBins, binCenters, FitIPRes_XiPrompt, binWidths);
    ipchi2Graphs["XpNonprompt"] = createParameterGraphs("FitIPXpNonprompt", nzTBins, binCenters, FitIPRes_XpNonprompt, binWidths);
    ipchi2Graphs["SigmaNonprompt"] = createParameterGraphs("FitIPSigmaNonprompt", nzTBins, binCenters, FitIPRes_SigmaNonprompt, binWidths);
    ipchi2Graphs["XiNonprompt"] = createParameterGraphs("FitIPXiNonprompt", nzTBins, binCenters, FitIPRes_XiNonprompt, binWidths);
    // Fractions and yields
    ipchi2Graphs["PromptFrac"] = createParameterGraphs("FitIPPromptFrac", nzTBins, binCenters, FitIPRes_PromptFrac, binWidths);
    ipchi2Graphs["PromptYield"] = createParameterGraphs("FitIPPromptYield", nzTBins, binCenters, FitIPRes_PromptYield, binWidths);
    ipchi2Graphs["NonPromptYield"] = createParameterGraphs("FitIPNonPromptYield", nzTBins, binCenters, FitIPRes_NonPromptYield, binWidths);
    
    // Create range graphs for IP chi2 parameters
    ipchi2RangeGraphs["XpPromptR"] = createGraphsAsymmErr("FitIPXpPromptRange", nzTBins, binCenters, FitIPRes_XpPrompt, binWidths);
    ipchi2RangeGraphs["SigmaPromptR"] = createGraphsAsymmErr("FitIPSigmaPromptRange", nzTBins, binCenters, FitIPRes_SigmaPrompt, binWidths);
    ipchi2RangeGraphs["XiPromptR"] = createGraphsAsymmErr("FitIPXiPromptRange", nzTBins, binCenters, FitIPRes_XiPrompt, binWidths);
    ipchi2RangeGraphs["XpNonpromptR"] = createGraphsAsymmErr("FitIPXpNonpromptRange", nzTBins, binCenters, FitIPRes_XpNonprompt, binWidths);
    ipchi2RangeGraphs["SigmaNonpromptR"] = createGraphsAsymmErr("FitIPSigmaNonpromptRange", nzTBins, binCenters, FitIPRes_SigmaNonprompt, binWidths);
    ipchi2RangeGraphs["XiNonpromptR"] = createGraphsAsymmErr("FitIPXiNonpromptRange", nzTBins, binCenters, FitIPRes_XiNonprompt, binWidths);
    ipchi2RangeGraphs["PromptFracR"] = createGraphsAsymmErr("FitIPPromptFracRange", nzTBins, binCenters, FitIPRes_PromptFrac, binWidths);
    
    // Save to parameter file
    TFile* fOutData = new TFile((fOutDataName + ".root").c_str(), "RECREATE");
    
    // Save fragmentation functions
    inclFragFunc->Write();
    promptFragFunc->Write();
    nonpromptFragFunc->Write();
    
    // Save mass parameter graphs
    for (const auto& graphSet : massGraphs) {
        for (const auto& graph : graphSet.second) {
            if (graph) graph->Write();
        }
    }
    
    // Save IP chi2 parameter graphs
    for (const auto& graphSet : ipchi2Graphs) {
        for (const auto& graph : graphSet.second) {
            if (graph) graph->Write();
        }
    }
    
    // Save range graphs
    for (const auto& [key, graph] : massRangeGraphs) {
        if (graph) graph->Write();
    }
    
    for (const auto& [key, graph] : ipchi2RangeGraphs) {
        if (graph) graph->Write();
    }
    
    fOutData->Close();
    
    // Save to histogram file
    TFile* fOutHisto = new TFile((fOutDataNameB + ".root").c_str(), "RECREATE");
    
    // Save fragmentation functions again
    inclFragFunc->Write();
    promptFragFunc->Write();
    nonpromptFragFunc->Write();
    
    // Save histograms
    for (int i = 0; i < nzTBins; ++i) {
        if (ipchi2HistoArray[i]) {
            ipchi2HistoArray[i]->SetName(("hIPChi2_" + std::to_string(i)).c_str());
            ipchi2HistoArray[i]->Write();
        }
        if (massHistoArray[i]) {
            massHistoArray[i]->SetName(("hMassSpectr_" + std::to_string(i)).c_str());
            massHistoArray[i]->Write();
        }
    }
    
    fOutHisto->Close();
    
    // Clean up
    delete inclFragFunc;
    delete promptFragFunc;
    delete nonpromptFragFunc;
    
    // Clean up mass graphs
    for (auto& graphSet : massGraphs) {
        for (auto& graph : graphSet.second) {
            delete graph;
        }
    }
    
    // Clean up IP chi2 graphs
    for (auto& graphSet : ipchi2Graphs) {
        for (auto& graph : graphSet.second) {
            delete graph;
        }
    }
    
    // Clean up range graphs
    for (auto& [key, graph] : massRangeGraphs) {
        delete graph;
    }
    
    for (auto& [key, graph] : ipchi2RangeGraphs) {
        delete graph;
    }
    
    std::cout << "Results saved to:" << std::endl;
    std::cout << "  Parameters: " << fOutDataName << ".root" << std::endl;
    std::cout << "  Histograms: " << fOutDataNameB << ".root" << std::endl;
}


// Add these implementations after the class declaration but before any other function implementations

// Implementation of createParameterGraphs
std::vector<TGraphErrors*> FitSpectraObject::createParameterGraphs(
    const std::string& key, int nBins, const std::vector<double>& xPos, 
    const std::vector<std::vector<float>>& yPosArr, const std::vector<double>& xWidth) {
    
    // Initialize vector of TGraphErrors pointers with nullptr
    std::vector<TGraphErrors*> graphs(4, nullptr);
    std::vector<std::string> endings = {"F", "Start", "LL", "HL"};
    
    // Arrays to hold data for graphs
    std::vector<double> tempArray(nBins);
    std::vector<double> tempArrayErr(nBins);
    
    // Process each value type (value, error, min, max)
    for (int i = 0; i < 5; ++i) {
        if (i == 1) { // Skip the error array for index references
            continue;
        }
        
        // Fill temporary arrays with data from yPosArr
        for (int j = 0; j < nBins; ++j) {
            if (j < static_cast<int>(yPosArr.size()) && i < static_cast<int>(yPosArr[j].size())) {
                tempArray[j] = yPosArr[j][i];
                if (i == 0 && 1 < static_cast<int>(yPosArr[j].size())) {
                    tempArrayErr[j] = yPosArr[j][1]; // Use error values from index 1
                } else {
                    tempArrayErr[j] = 0.0;
                }
            } else {
                tempArray[j] = 0.0;
                tempArrayErr[j] = 0.0;
            }
        }
        
        if (i == 0) {
            // Main graph with values and errors
            graphs[0] = new TGraphErrors(
                nBins,
                xPos.data(),
                tempArray.data(),
                xWidth.data(),
                tempArrayErr.data()
            );
            
            // Set min/max if we have values
            if (!tempArray.empty()) {
                double minVal = *std::min_element(tempArray.begin(), tempArray.end());
                double maxVal = *std::max_element(tempArray.begin(), tempArray.end());
                graphs[0]->SetMinimum(minVal);
                graphs[0]->SetMaximum(maxVal);
            }
            
            graphs[0]->SetName((key + endings[0]).c_str());
        } else {
            // Graphs for range limits
            int idx = i - 1;
            if (idx < 4 && idx >= 0) { // Check if index is valid
                // For these graphs, we only want x errors, not y errors
                graphs[idx] = new TGraphErrors(
                    nBins,
                    xPos.data(),
                    tempArray.data(),
                    xWidth.data(),
                    nullptr  // No y errors for these graphs
                );
                
                // Set min/max if we have values
                if (!tempArray.empty()) {
                    double minVal = *std::min_element(tempArray.begin(), tempArray.end());
                    double maxVal = *std::max_element(tempArray.begin(), tempArray.end());
                    graphs[idx]->SetMinimum(minVal);
                    graphs[idx]->SetMaximum(maxVal);
                }
                
                graphs[idx]->SetName((key + endings[idx]).c_str());
            }
        }
    }
    
    return graphs;
}

// Implementation of createGraphsAsymmErr
TGraphAsymmErrors* FitSpectraObject::createGraphsAsymmErr(
    const std::string& key, int nBins, const std::vector<double>& xPos,
    const std::vector<std::vector<float>>& yPosArr, const std::vector<double>& xWidth) {
    
    // Arrays for start values, lower and upper limits
    std::vector<double> startArray(nBins, 0.0);
    std::vector<double> llArray(nBins, 0.0);
    std::vector<double> hlArray(nBins, 0.0);
    
    // Fill arrays from yPosArr
    for (int i = 0; i < nBins; ++i) {
        if (i < static_cast<int>(yPosArr.size())) {
            if (2 < static_cast<int>(yPosArr[i].size())) startArray[i] = yPosArr[i][2]; // Start values at index 2
            if (3 < static_cast<int>(yPosArr[i].size())) llArray[i] = yPosArr[i][3];    // Lower limits at index 3
            if (4 < static_cast<int>(yPosArr[i].size())) hlArray[i] = yPosArr[i][4];    // Upper limits at index 4
        }
    }
    
    // Calculate errors (difference between start value and limits)
    std::vector<double> errLow(nBins, 0.0);
    std::vector<double> errHigh(nBins, 0.0);
    
    for (int i = 0; i < nBins; ++i) {
        errLow[i] = startArray[i] - llArray[i];
        errHigh[i] = hlArray[i] - startArray[i];
        
        // Make sure errors are non-negative
        errLow[i] = std::max(0.0, errLow[i]);
        errHigh[i] = std::max(0.0, errHigh[i]);
    }
    
    // Create asymmetric error graph
    TGraphAsymmErrors* graph = new TGraphAsymmErrors(
        nBins,
        xPos.data(),
        startArray.data(),
        xWidth.data(),  // X low errors
        xWidth.data(),  // X high errors
        errLow.data(),  // Y low errors
        errHigh.data()  // Y high errors
    );
    
    // Set min/max if we have values
    if (!llArray.empty() && !hlArray.empty()) {
        double minVal = *std::min_element(llArray.begin(), llArray.end());
        double maxVal = *std::max_element(hlArray.begin(), hlArray.end());
        graph->SetMinimum(minVal);
        graph->SetMaximum(maxVal);
    }
    
    graph->SetName(key.c_str());
    
    return graph;
}

void MassFitter(TString inputFile = "/media/niviths/SSD2/lhcb_analysis_SSD/20250501_Pbp_data/Pbp_data_filterV1.root", bool isMC = false, bool isFitSingleBin = false, bool isZtObservable = true)
{

    std::string mcTag = isMC ? "MC" : "";
    std::string obsTag = isZtObservable ? "zT" : "dR";


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
            tree  // Pass the tree
        );
        fitter.startFitting();
    } 
    else {
        // Multi-bin configuration
        std::vector<double> zBins = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};  // D0 zT bins
        std::vector<double> rBins = {0, 0.015, 0.03, 0.06, 0.1, 0.2, 0.5};  // D0 R bins
        
        // pT binning for jets
        std::vector<double> startPt = {5, 10, 15, 20, 30};
        std::vector<double> endPt = {10, 15, 20, 30, 50};
        
        // Process each pT bin
        for (size_t jetBin = 0; jetBin < startPt.size(); ++jetBin) {
            std::pair<double, double> jetPt(startPt[jetBin], endPt[jetBin]);
            
            // Choose bin array based on observable type
            const std::vector<double>& binArray = isZtObservable ? zBins : rBins;
            
            // Create and run fit object
            FitSpectraObject fitter(
                jetPt, isMC, binArray,
                isZtObservable,
                tree  // Pass the tree
            );
            fitter.startFitting();
        }
    }

}