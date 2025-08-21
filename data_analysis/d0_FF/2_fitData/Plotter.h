
#ifndef PLOTTER_H
#define PLOTTER_H

#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cstdlib>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TEfficiency.h"
#include "TAxis.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TBranch.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooStats/SPlot.h"

class Plotter {
private:
    std::string resonance;
    int fitBin;
    bool isBinned;
    std::string range;
    std::string nameKey;
    std::string obsSelection;
    std::string format;
    std::string basepath;

    // Helper methods
    void setGraphHisto(TH1* histo, const std::string& xTitle, const std::string& yTitle, 
                      double border = 0.1, bool logPlot = false);
    std::string getUniqueId();
    void ensureDirectoryExists(const std::string& path);

public:
    // Constructor
    Plotter(const std::string& resonance, const std::string& basepath = "./", 
           int bin = 0, bool binned = false, const std::string& range = "", 
           const std::string& name = "");
    
    // Destructor
    ~Plotter() {}
    
    // Main plotting methods
    TH1* individualMassFitPlot(RooRealVar* sigYieldParam, RooAbsPdf* extendedPdf, 
                              RooRealVar* massVar, RooDataSet* data, 
                              const std::string& fitTypeName, bool isZtObservable);
    
    TH1* ipchi2FitPlot(const std::string& resonance, RooRealVar* logIpchi2, 
                      RooDataSet* data, RooAbsPdf* totalPdf, 
                      RooAbsPdf* nonpromptPdf = nullptr, 
                      RooAbsPdf* promptPdf = nullptr, 
                      RooAbsPdf* backgroundPdf = nullptr,
                      RooAbsReal* promptYield = nullptr,    // Changed from RooRealVar* to RooAbsReal*
                      RooAbsReal* nonpromptYield = nullptr); // Changed from RooRealVar* to RooAbsReal*
    
    // Other methods from Python implementation
    void plotMultiple(const std::string& name, RooRealVar* lifetimeTag, 
                     bool drawLogDiff, RooDataSet* data1, 
                     RooDataSet* data2 = nullptr, RooDataSet* data3 = nullptr);
                     
    std::pair<double, double> extractCorParamRnd(const std::string& idString, 
                                               RooRealVar* massVar, 
                                               std::vector<TLegend*>& myLegendList, 
                                               RooDataSet* data, 
                                               RooDataSet* data2 = nullptr, 
                                               RooDataSet* data3 = nullptr);
    
    TH1* teffMap2D(TEfficiency* teff, const std::string& effType, const std::string& dataType, 
                  bool closureTF1 = false, bool closureHisto = false);
    
    std::tuple<double, double, double> splotVals(const std::string& resonance, 
                                               RooAbsPdf* extendedPdf, 
                                               RooRealVar* massVar, 
                                               RooDataSet* data, 
                                               TTree* ttree, 
                                               RooRealVar* ns, 
                                               RooRealVar* nb, 
                                               const std::string& fidCutString, 
                                               TFile* sFile);
};

#endif // PLOTTER_H