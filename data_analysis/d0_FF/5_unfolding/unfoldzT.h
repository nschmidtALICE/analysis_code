// UnfoldSpectra.h
#ifndef UNFOLDSPECTRA_H
#define UNFOLDSPECTRA_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <filesystem>
#include <cmath>
#include <algorithm>

// ROOT includes
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TNtuple.h"
#include "TRandom.h"
#include <TPaveText.h>
#include <TROOT.h>
#include <TDirectory.h>
// RooUnfold includes
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"

class UnfoldSpectraClass {
private:
    // General variables
    bool useTagPt;
    std::string inPathRM;
    std::string dictKey;
    std::string figureTag;
    std::string outpathBase;
    bool isPrompt;
    bool applyRMCut;
    std::string inFileNRM;
    std::string inFileNData;
    
    // Binning variables
    std::vector<double> zBinsArrayDet;
    std::vector<double> zBinsArrayTruth;
    std::vector<int> pTBinsArrayTruth;
    std::vector<int> pTBinsArrayDet;
    
    // Output variables
    std::vector<TH2D*> unfoldedArr2D;
    std::vector<std::vector<TH1D*>> unfoldedArrPerBin;
    TH2D* hOriginalPrior;
    TH2D* hCurrentPrior;
    TH2D* externalPrior;
    int errorType;
    
    // Output path
    std::string outPath;
    
    // Input files
    TFile* fResponse;
    TFile* fData;
    
    // Data spectra
    std::vector<TH1D*> measuredSpectraArray;
    TH2D* measuredSpectra2D;

    // Unfolding label
    std::string unfoldLabel;

public:
    // Constructor
    UnfoldSpectraClass(const std::string& promptFlag, const std::string& inFileNameRM, 
                      const std::string& inFileNameData, const std::string& resonance,
                      const std::vector<int>& ptRangeArray);
    
    // Destructor
    ~UnfoldSpectraClass();
    
    // Public methods
    void unfold1D(int ptBin, int regParam, int iteration, int powerLawOffset = 0);
    void unfold2D(int regParam, int iteration, const std::string& tag = "");
    void provideExtPrior(const std::string& tag, const std::string& FileName = "", 
                         const std::string& specialType = "");

private:
    // Private helper methods
    TH2D* prepareRMWeights(int regParam, int iteration);
    double getEventWeight(TH2D* hweights, double zTValue, double pTValue);
    void RefoldingTest(int ptBin, int nIter, TH1D* histo, RooUnfoldResponse* RM1, RooUnfoldResponse* RM2);
    void RefoldingTest2D(int nIter, RooUnfoldResponse* RM1, RooUnfoldResponse* RM2);
    void UnfoldingTest(RooUnfoldResponse* response, int nIter);
    void UnfoldingTest2D(RooUnfoldResponse* response, int nIter);
    void plotPrior2D(int nIter, RooUnfoldResponse* RooUnfoldRM);
    void plotUnfoldingEffect(int ptBin, int regParam);
    void plotUnfoldingEffect2D(int regParam);
    void StabilityTest2D(int nIter, bool plotAll = false);
    void StabilityTest(int ptBin, int nIter, bool useRatioNorm = false);
    void TestRegParam2D();
    void StatTestRM2D(RooUnfoldResponse* response);
    void StatTestRM(RooUnfoldResponse* response, int ptBin);
    std::pair<std::vector<TH1D*>, std::vector<TH1D*>> calculateDispersion(
        const std::vector<TH2D*>& histogramArray2D,
        std::vector<TH1D*> hvariationResult,
        std::vector<TH1D*> hvariationError,
        int nzTBins, int npTBins);
    TH2D* SmearPoints(TH2D* hist, TRandom* fRandom, int factor);
    RooUnfoldResponse* PrepareResponseMatrix2D(int part = 0);
    RooUnfoldResponse* PrepareResponseMatrix3D(int part, TH2D* weightHistogram);
    void getKinEfficiency(int bin);
    TH2D* getResponseMatrix(int part, const std::string& xAxisVar, const std::string& yAxisVar, 
                          bool fineBinning = false, int bin = 0, bool isCut = false, 
                          TFile* externalFile = nullptr);
    std::vector<TH1D*> getRawSpectra(bool mergeTo2D = false);
    TH2D* getRawSpectra2D();
    void deScaleGraph(TGraphErrors* graph);
    TH1D* scaleHistogram(TH1D* histo);
    void plotCorrelationCoefficients(TMatrixD& covarianceMatrix, int i, const std::string& Name);
    void plotJetResponse();
    std::vector<TLine*> drawLines();
    std::string makeOutDir(const std::string& subDirName);
    void plotHist(TObject* hist, const std::string& outputFilename, 
                 const std::string& drawOptions = "", 
                 const std::string& xTitle = "", 
                 const std::string& yTitle = "");
    void saveResult(int regParam, int nIter, const std::string& tag);
    
};

// Function declarations for the main procedures
void unfoldzTDistributionLHCb(int variation);
void unfoldzTAllDistributionLHCb(int variation);

void SaveCanvasQuietly(TCanvas* canvas, const std::string& filename) {
    int oldLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;
    canvas->SaveAs(filename.c_str());
    gErrorIgnoreLevel = oldLevel;
}

#endif // UNFOLDSPECTRA_H