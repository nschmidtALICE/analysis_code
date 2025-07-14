// Fitter.h
#ifndef FITTER_H
#define FITTER_H

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TDirectory.h"
#include "TCanvas.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFormulaVar.h"
#include "RooBukinPdf.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooPlot.h"

// Forward declaration for Plotter class
class Plotter;

// Struct to hold parameter values with ranges
struct ParamConfig {
    double value;
    double min;
    double max;
    
    ParamConfig() : value(0), min(0), max(0) {}
    ParamConfig(double v, double mn, double mx) : value(v), min(mn), max(mx) {}
};

// Struct for mass range configurations
struct MassConfig {
    ParamConfig sigma1;
    ParamConfig deltasigma;
    ParamConfig mean;
    ParamConfig n;
    ParamConfig dg_frac;
    ParamConfig pol1;
    ParamConfig pol2;
    std::pair<double, double> massRange;
    ParamConfig sigYield;
    ParamConfig sigYieldLim;
    ParamConfig bkgYield;
    ParamConfig bkgYieldLim;
    std::pair<double, double> signalRegion;
    std::pair<double, double> sbRegion;
};

// Struct for IP chi2 configurations
struct IPChi2Config {
    std::pair<double, double> logIpchi2Range;
    
    // Prompt component
    ParamConfig xpPrompt;
    ParamConfig sigmaPrompt;
    ParamConfig xiPrompt;
    ParamConfig rho1Prompt;
    ParamConfig rho2Prompt;
    
    // Non-prompt component
    ParamConfig xpNonprompt;
    ParamConfig sigmaNonprompt;
    ParamConfig xiNonprompt;
    ParamConfig rho1Nonprompt;
    ParamConfig rho2Nonprompt;
    
    // Fraction
    ParamConfig promptFrac;
    
    // Background
    ParamConfig bkgParam1;
    ParamConfig bkgParam2;
};

class Fitter {
private:
    // File handling
    TFile* tfilePV;
    TFile* fInFileHisto;
    std::string TestFilename;
    std::string outfilePath;
    std::string resonance;
    
    // Configuration
    bool isMC;
    int nBins;
    bool updateStartValues;
    bool isZtObservable;
    
    // B-decay ranges
    std::vector<double> BdecayMinVal;
    std::vector<double> BdecayMaxVal;
    
    // Data structures
    TTree* inTree;
    std::vector<TH1*> hTime;
    std::vector<TH1*> hYield;
    
    // Dictionaries
    std::map<std::string, MassConfig> massDict;
    std::map<std::string, IPChi2Config> ipchi2Dict;
    
    // Initialize dictionaries
    void initDictionary();

public:
    // Modified constructor that takes TTree* directly
    Fitter(TTree* tree, 
           const std::string& resonanceType = "", 
           int numBins = 1,
           bool isZtObservable = true,
           bool isMCData = false, 
           const std::string& outputPath = ".", 
           bool update = false);
    
    // Destructor
    ~Fitter();
    
    // Dictionary methods
    void updateDictionary(RooAbsPdf* signalPdf, RooAbsData* data, const std::string& fitFunc);
    void updateSBfraction(const std::string& resonance, int bin, double ptLimLow);
    void updateSigYield(const std::string& resonance);
    void updateBKGYield(const std::string& resonance);
    
    // Utility methods
    std::string fiducialCutString(const std::pair<double, double>& jetPt, 
                                  const std::pair<double, double>& tagMass);
    
    // Dataset creation
    RooDataSet* createDataSet(const std::string& resonance, 
                              const std::string& name, 
                              const std::string& fidCutString = "", 
                              bool isMass = true, 
                              int bin = 0, 
                              int corrVer = -1, 
                              const std::string& fitVarName = "ipchi2");
    
    // Fitting methods
    std::tuple<TH1*, std::vector<double>, std::vector<double>> 
    massFit(const std::string& resonance, 
            RooDataSet* data, 
            const std::string& fitTypeName = "DCB", 
            int bin = -1, 
            const std::string& zRange = "", 
            bool splot = false, 
            TFile* sFile = nullptr);
    
    std::tuple<TH1*, std::vector<double>, std::vector<double>> 
    ipchi2Fit(const std::string& resonance, 
              RooDataSet* data, 
              RooDataSet* background = nullptr, 
              const std::string& figKey = "All", 
              int bin = 0, 
              const std::string& zRange = "");
    
    // Accessor methods
    MassConfig getMassDict(const std::string& resonance) const {
        auto it = massDict.find(resonance);
        if (it != massDict.end()) {
            return it->second;
        }
        throw std::runtime_error("Resonance not found in mass dictionary: " + resonance);
    }
    
    IPChi2Config getIPChi2Dict(const std::string& resonance) const {
        std::string key = "Sig" + resonance;
        auto it = ipchi2Dict.find(key);
        if (it != ipchi2Dict.end()) {
            return it->second;
        }
        throw std::runtime_error("Resonance not found in IP chi2 dictionary: " + key);
    }
};

#endif // FITTER_H