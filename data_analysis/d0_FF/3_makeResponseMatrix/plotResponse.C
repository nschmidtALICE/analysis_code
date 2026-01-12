// Converted from plot_Response_NoVec.py

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sys/stat.h>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TGaxis.h"
#include "TSystem.h"

// RooFit includes
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooFit.h"
#include "RooBinning.h"

class RMplotter {
private:
    // Timer variables
    std::clock_t startTime;
    double currentTime;

    // Format options
    std::string fileFormat;
    bool isMC;
    std::vector<int> ColorArray;
    std::string dictKey;

    // Drawing objects
    TLine* diagonal_line;
    TLine* diagonal_lineDark;

    // Binning
    std::vector<double> pTBinning;
    std::vector<double> zTBinArray;

    // Files and data
    TFile* fInFile;
    TTree* inTreeRM;
    RooDataSet* dataSet;
    std::string output;

    // Member functions - private helpers
    std::pair<std::string, std::string> setupInputOutput();
    RooDataSet* createDataSet();

public:
    // Constructor and destructor
    RMplotter();
    ~RMplotter();

    // Public methods
    RooDataSet* reduceDataSet(double jetPtMin, double jetPtMax, double dRCut = 1.0);
    void plotQAPlots(RooDataSet* dataReduced, RooRealVar& tagMDet, RooRealVar& tagMPart, const std::string& namekey = "");
    void plotRMJet(double ptLow, double ptHigh, double dR = 1.0);
    void plotRMTag(double dR = 1.0);
    void plotRMRapidity(double dR = 1.0);
    
    // Helper methods for plotting
    TPad*** partitionCanvas(TCanvas* C, int Nx, int Ny, double lMargin, double rMargin, double bMargin, double tMargin);
    TH1* setFrame(TH1* h, double xFactor, double yFactor, int i, int Nx, int j, int Ny);
    void setHisto(TH1* Histo, const std::string& Xtitle, const std::string& Ytitle, int big = 0, double border = 0.1);
    void setTH1Histo(TH1* Histo, const std::string& Xtitle, const std::string& Ytitle, int big, double border = 0.1);
};

// Constructor implementation
RMplotter::RMplotter() {
    // Initialize timing
    startTime = std::clock();
    currentTime = 0;
    
    // Set format options
    fileFormat = "pdf";
    isMC = true;
    
    // Initialize color array
    ColorArray = {kViolet+5, kAzure, kAzure-4, kCyan-6, kGreen-3, kTeal-6};
    
    // Set particle type
    dictKey = "D0";
    
    // Create diagonal lines
    diagonal_line = new TLine();
    diagonal_line->SetX1(0);
    diagonal_line->SetX2(1);
    diagonal_line->SetY1(0);
    diagonal_line->SetY2(1);
    diagonal_line->SetLineWidth(4);
    diagonal_line->SetLineStyle(2);
    diagonal_line->SetLineColor(kGray);
    
    diagonal_lineDark = new TLine();
    diagonal_lineDark->SetX1(0);
    diagonal_lineDark->SetX2(1);
    diagonal_lineDark->SetY1(0);
    diagonal_lineDark->SetY2(1);
    diagonal_lineDark->SetLineWidth(4);
    diagonal_lineDark->SetLineStyle(2);
    diagonal_lineDark->SetLineColor(kGray+2);
    
    // Setup binning
    pTBinning = {0, 5, 10, 15, 20, 30, 200};
    // zTBinArray = {0.2, 0.5, 0.65, 0.75, 0.85, 0.95, 1};
    zTBinArray = {
            0.0, 0.05, 0.1, 0.15, 0.2, 
            0.25, 0.3, 0.35, 0.4, 0.45, 
            0.5, 0.55, 0.6, 0.65, 0.7, 
            0.75, 0.8, 0.85, 0.9, 0.95, 
            1.0};
    // zTBinArray = {0.2, 0.5, 0.65, 0.75, 0.85, 0.95, 1};
    
    // Setup input/output
    std::pair<std::string, std::string> io = setupInputOutput();
    
    // Open input file
    fInFile = new TFile(io.first.c_str());
    if (!fInFile || fInFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << io.first << std::endl;
        exit(1);
    }
    inTreeRM = (TTree*)fInFile->Get("Response");
    if (!inTreeRM) {
        std::cerr << "Error: Could not find Response tree in file" << std::endl;
        exit(1);
    }
    
    // Create dataset
    dataSet = createDataSet();
    output = io.second;
}

// Destructor implementation
RMplotter::~RMplotter() {
    // Clean up
    delete diagonal_line;
    delete diagonal_lineDark;
    delete dataSet;
    
    if (fInFile) {
        fInFile->Close();
        delete fInFile;
    }
    
    // Calculate total execution time
    double duration = (std::clock() - startTime) / (double)CLOCKS_PER_SEC;
    std::cout << "Total execution time: " << duration << " seconds" << std::endl;
}

// Setup input and output paths
std::pair<std::string, std::string> RMplotter::setupInputOutput() {
    std::string inputDir = "/media/niviths/local/analysis_code/data_analysis/d0_FF/3_makeResponseMatrix";
    
    // Input file path
    std::string fFileName = "/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/53to56_response.root";
    // std::string fFileName = "/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/response.root";
    // std::string fFileName = "/media/niviths/SSD2/lhcb_analysis_SSD/20250609_merged/1122665_response.root";
    // std::string fFileName = "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_allMC.root";
    // std::string fFileName = "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_allMC.root";
    // std::string fFileName = "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_21_MC_output_D0FF_filterV1.root";
    // std::string fFileName = "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_21_MC_output_D0FF_filterV1.root";
    
    std::cout << "Found file: " << fFileName << std::endl;
    
    // Output directory with current date stamp (YYYYMMDD)
    std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);
    char dateBuf[16];
    std::strftime(dateBuf, sizeof(dateBuf), "%Y%m%d", &tm);
    std::string outputDirFig = inputDir + "/D0_FF_ResMatr_" + dateBuf + "/";
    std::cout << "Using output directory: " << outputDirFig << std::endl;
    
    // Create output directory if it doesn't exist
    struct stat info;
    if (stat(outputDirFig.c_str(), &info) != 0) {
        std::string mkdirCmd = "mkdir -p " + outputDirFig;
        int status = system(mkdirCmd.c_str());
        if (status == 0) {
            std::cout << "Created directory: " << outputDirFig << std::endl;
        } else {
            std::cerr << "Failed to create directory: " << outputDirFig << std::endl;
        }
    }
    
    return std::make_pair(fFileName, outputDirFig);
}

// Create the dataset with cuts
RooDataSet* RMplotter::createDataSet() {
    // Dataset Parameters with updated branch names
    RooRealVar pt_jetDet("jet_pt_det", "Detector-level jet p_{T} (GeV/c)", 0, 500);
    RooRealVar pt_tagDet("d0_pt_det", "Detector-level D0 p_{T} (GeV/c)", 0, 80);
    RooRealVar etaDet("jet_eta_det", "Detector-level jet #eta", 1.5, 5.5);
    RooRealVar phiDet("jet_phi_det", "Detector-level jet #phi", -3.5, 3.5);
    RooRealVar tagMDet("d0_mass_det", "Detector-level D0 mass (GeV/c^{2})", 1.81, 1.935);
    RooRealVar zTDet("d0_z_det", "Detector-level z_{T}", 0, 1.01);
    RooRealVar d0EtaDet("d0_eta_det", "Detector-level D0 #eta", 1.5, 5.5);
    
    RooRealVar pt_jetPart("jet_pt_mc", "Particle-level jet p_{T} (GeV/c)", 0, 500);
    RooRealVar pt_tagPart("d0_pt_mc", "Particle-level D0 p_{T} (GeV/c)", 0, 80);
    RooRealVar etaPart("jet_eta_mc", "Particle-level jet #eta", 1.5, 5.5);
    RooRealVar phiPart("jet_phi_mc", "Particle-level jet #phi", -3.5, 3.5);
    RooRealVar tagMPart("d0_mass_mc", "Particle-level D0 mass (GeV/c^{2})", 1.81, 1.935);
    RooRealVar zTPart("d0_z_mc", "Particle-level z_{T}", 0, 1.01);
    RooRealVar d0EtaPart("d0_eta_mc", "Particle-level D0 #eta", 1.5, 5.5);
    
    RooRealVar nConstDet("jet_nconst_det", "Detector-level jet constituent count", 0, 50);
    RooRealVar nConstPart("jet_nconst_mc", "Particle-level jet constituent count", 0, 50);
    RooRealVar dR("jet_dr", "ΔR between matched jets", -1, 1);
    
    RooArgSet cutVars;
    cutVars.add(pt_jetDet);
    cutVars.add(pt_tagDet);
    cutVars.add(tagMDet);
    cutVars.add(pt_jetPart);
    cutVars.add(pt_tagPart);
    cutVars.add(tagMPart);
    cutVars.add(dR);
    cutVars.add(zTDet);
    cutVars.add(zTPart);
    cutVars.add(nConstDet);
    cutVars.add(nConstPart);
    cutVars.add(d0EtaDet);
    cutVars.add(d0EtaPart);
    
    // Define cuts
    double massCutMin = 1.81;
    double massCutMax = 1.935;
    std::stringstream cutStream;
    cutStream << "d0_mass_mc > " << massCutMin << " && d0_mass_mc < " << massCutMax
              << " && jet_dr > -1 && jet_nconst_det>1";
    std::string cutString = cutStream.str();
    
    std::cout << "Creating dataset with cuts: " << cutString << std::endl;
    
    // Create dataset
    RooDataSet* data = new RooDataSet("Test1Corr", "Test1Corr", inTreeRM, cutVars, cutString.c_str());
    std::cout << "Dataset created with " << data->numEntries() << " entries" << std::endl;
    
    return data;
}

// Apply additional cuts to the dataset
RooDataSet* RMplotter::reduceDataSet(double jetPtMin, double jetPtMax, double dRCut) {
    // Create cut string with updated branch names
    std::stringstream cutStream;
    cutStream << "jet_pt_det > " << jetPtMin << " && jet_pt_det < " << jetPtMax;
    
    // Add dR cut if needed
    if (dRCut < 1.0) {
        cutStream << " && jet_dr < " << dRCut;
    }
    
    std::string cutString = cutStream.str();
    std::cout << "Applying additional cuts: " << cutString << std::endl;
    
    // Create variables needed for dataset operations with updated branch names
    RooRealVar jetPtDet("jet_pt_det", "jet_pt_det", 0, 500);
    RooRealVar dR("jet_dr", "jet_dr", -1, 1);
    RooArgSet vars(jetPtDet, dR);
    
    // Apply cuts
    RooDataSet* dataReduced = (RooDataSet*)dataSet->reduce(cutString.c_str());
    std::cout << "Reduced dataset has " << dataReduced->numEntries() << " entries" << std::endl;
    
    std::cout << "------------------------" << std::endl;
    
    // Create variables for QA plots with updated branch names
    RooRealVar tagMDet("d0_mass_det", "d0_mass_det", 1.81, 1.935);
    RooRealVar tagMPart("d0_mass_mc", "d0_mass_mc", 1.81, 1.935);
    
    // Plot QA plots
    plotQAPlots(dataReduced, tagMDet, tagMPart, "Mass");
    
    return dataReduced;
}

void RMplotter::plotRMJet(double ptLow, double ptHigh, double dR) {
    // Get reduced dataset
    RooDataSet* reducedData = reduceDataSet(ptLow, ptHigh, dR);
    if (!reducedData) {
        std::cerr << "Failed to create reduced dataset" << std::endl;
        return;
    }
    
    // Define variables for jet pT
    RooRealVar pt_jetDet("jet_pt_det", "Detector-level jet p_{T} (GeV/c)", 0, 100);
    RooRealVar pt_jetPart("jet_pt_mc", "Particle-level jet p_{T} (GeV/c)", 0, 100);
    
    // Create 2D histogram from dataset
    TH2D* hh_data = dynamic_cast<TH2D*>(reducedData->createHistogram("hist", pt_jetDet, 
                                               RooFit::Binning(50, 0, 100), 
                                               RooFit::YVar(pt_jetPart, RooFit::Binning(50, 0, 100))));
    
    if (!hh_data) {
        std::cerr << "Failed to create 2D histogram from dataset" << std::endl;
        delete reducedData;
        return;
    }
    
    // Style and normalize the histogram
    hh_data->SetTitle("");
    hh_data->GetXaxis()->SetTitle("Detector-level jet p_{T} (GeV/c)");
    hh_data->GetYaxis()->SetTitle("Particle-level jet p_{T} (GeV/c)");
    hh_data->GetZaxis()->SetTitle("Entries");
    
    // Set style options
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetPadRightMargin(0.15);
    
    // Create canvas
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800*2, 400*2);
    canvas->cd(1);
    
    // Create pad
    TPad* myPad = new TPad("myPad", "The pad", 0, 0, 1, 1);
    myPad->SetLeftMargin(0.7);
    myPad->SetTopMargin(0.07);
    myPad->SetRightMargin(0.1);
    myPad->SetBottomMargin(0.15);
    myPad->Draw();
    myPad->cd();
    
    // Draw histogram
    hh_data->Draw("colz");
    
    // Save canvas to file
    std::string ptTag = "pT" + std::to_string(int(ptLow)) + "-" + std::to_string(int(ptHigh));
    std::string fileName = output + dictKey + "Bin" + ptTag + "." + fileFormat;
    canvas->SaveAs(fileName.c_str());
    
    // Clean up
    delete canvas;
    delete hh_data;
    delete reducedData;
}

void RMplotter::plotQAPlots(RooDataSet* dataReduced, RooRealVar& tagMDet, RooRealVar& tagMPart, const std::string& namekey) {
    std::cout << "**plot dataset**" << std::endl;
    
    // Create histograms from datasets
    TH1* h_dataFullDet = dataSet->createHistogram("h_MassSignalDetFull_data", tagMDet);
    TH1* h_dataRedDet = dataReduced->createHistogram("h_MassSignalDetRed_data", tagMDet);
    
    TH1* h_dataFullPart = dataSet->createHistogram("h_MassSignalPartFull_data", tagMPart);
    TH1* h_dataRedPart = dataReduced->createHistogram("h_MassSignalDetRed_data", tagMPart);
    
    // Create canvas
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800*2, 400*2);
    canvas->Divide(2);
    
    // Left pad
    canvas->cd(1);
    // TPad* myPad1 = new TPad("myPad1", "The pad", 0, 0, 1, 1);
    // myPad1->SetLeftMargin(0.7);
    // myPad1->SetTopMargin(0.07);
    // myPad1->SetRightMargin(0.1);
    // myPad1->SetBottomMargin(0.15);
    // myPad1->Draw();
    // myPad1->cd();
    
    // Draw detector level histograms
    h_dataFullDet->SetLineColor(1);
    h_dataFullDet->Draw("hist");
    h_dataRedDet->SetLineColor(kBlue);
    h_dataRedDet->Draw("hist same");
    
    // Right pad
    canvas->cd(2);
    // TPad* myPad2 = new TPad("myPad2", "The pad", 0, 0, 1, 1);
    // myPad2->SetLeftMargin(0.7);
    // myPad2->SetTopMargin(0.07);
    // myPad2->SetRightMargin(0.1);
    // myPad2->SetBottomMargin(0.15);
    // myPad2->Draw();
    // myPad2->cd();
    
    // Draw particle level histograms
    h_dataFullPart->SetLineColor(1);
    h_dataFullPart->Draw("hist");
    h_dataRedPart->SetLineColor(kBlue);
    h_dataRedPart->Draw("hist same");
    
    // Create legend
    TLegend* myLegend2 = new TLegend(0.2, 0.8, 0.4, 0.9);
    myLegend2->SetTextFont(42);
    myLegend2->SetBorderSize(0);
    myLegend2->SetFillStyle(0);
    myLegend2->SetFillColor(0);
    myLegend2->SetMargin(0.25);
    myLegend2->SetTextSize(0.04);
    myLegend2->AddEntry(h_dataFullPart, "Data Full", "L");
    myLegend2->AddEntry(h_dataRedPart, "Data Reduced", "L");
    myLegend2->Draw();
    
    // Save canvas
    std::string fileName = output + dictKey + "QA" + namekey + "." + fileFormat;
    canvas->SaveAs(fileName.c_str());
    
    // Clean up
    delete canvas;
    delete h_dataFullDet;
    delete h_dataRedDet;
    delete h_dataFullPart;
    delete h_dataRedPart;
    delete myLegend2;
}


void RMplotter::plotRMTag(double dR) {
    // Get reduced dataset for entire jet pT range
    RooDataSet* reducedData = reduceDataSet(0, 200, dR);
    if (!reducedData) {
        std::cerr << "Failed to create reduced dataset for plotRMTag" << std::endl;
        return;
    }
    
    // Define variables for response matrix
    RooRealVar pt_jetDet("jet_pt_det", "jet_pt_det", 0, 200);
    RooRealVar pt_jetPart("jet_pt_mc", "jet_pt_mc", 0, 200);
    RooRealVar zTDet("d0_z_det", "d0_z_det", 0, 1.01);
    RooRealVar zTPart("d0_z_mc", "d0_z_mc", 0, 1.01);
    
    // Create custom binning for the plot
    std::vector<double> arrayDet(zTBinArray);
    std::vector<double> arrayPart;
    arrayPart.push_back(0); // Insert 0 at the beginning
    arrayPart.insert(arrayPart.end(), arrayDet.begin(), arrayDet.end());
    
    // Create ROOT binning objects
    RooBinning tbinsDet(arrayDet.size()-1, &arrayDet[0], "zTBinningDet");
    RooBinning tbinsPart(arrayPart.size()-1, &arrayPart[0], "zTBinningPart");
    
    // Calculate number of bins for generator and detector levels
    int nBinsGenLvlv = pTBinning.size() - 1;
    int nBinsDetLvlv = pTBinning.size() - 1;
    int maxCombinations = nBinsGenLvlv * nBinsDetLvlv;
    
    std::cout << "Max No combinations: " << maxCombinations << std::endl;
    std::cout << "pT Part bins: " << nBinsGenLvlv << std::endl;
    std::cout << "pT Det bins: " << nBinsDetLvlv << std::endl;
    
    // Create arrays to store histograms and metadata
    std::vector<std::vector<TH2F*>> hh_data(nBinsDetLvlv, std::vector<TH2F*>(nBinsGenLvlv, nullptr));
    std::vector<double> maxArray(maxCombinations, 0);
    std::vector<std::string> legGen(nBinsGenLvlv);
    std::vector<std::string> legDet(nBinsDetLvlv);
    std::vector<std::vector<TLatex*>> systemArray(nBinsDetLvlv, std::vector<TLatex*>(nBinsGenLvlv, nullptr));
    
    // Create all histograms for the zT responses
    int histIndex = 0;
    for (int indexGenerator = 0; indexGenerator < nBinsGenLvlv; ++indexGenerator) {
        // Define generator level pT range
        std::string ptRangePart = "jet_pt_mc > " + std::to_string(pTBinning[indexGenerator]) + 
                                " && jet_pt_mc < " + std::to_string(pTBinning[indexGenerator+1]);
        
        // Build legend entry for generator level
        legGen[indexGenerator] = Form("p_{T}^{jet,gen}: %1.0f-%1.0f", 
                                    pTBinning[indexGenerator], pTBinning[indexGenerator+1]);
        
        for (int indexDet = 0; indexDet < nBinsDetLvlv; ++indexDet) {
            // Define detector level pT range
            std::string ptRangeDet = " && jet_pt_det > " + std::to_string(pTBinning[indexDet]) +
                                   " && jet_pt_det < " + std::to_string(pTBinning[indexDet+1]);
            std::string cuts = ptRangePart + ptRangeDet;
            
            // Build legend entry for detector level
            legDet[indexDet] = Form("p_{T}^{jet,det}: %1.0f-%1.0f", 
                                  pTBinning[indexDet], pTBinning[indexDet+1]);
            
            // Get data in this specific bin
            RooDataSet* dataInBin = dynamic_cast<RooDataSet*>(reducedData->reduce(RooArgSet(zTDet, zTPart, pt_jetDet, pt_jetPart), cuts.c_str()));
            if (!dataInBin) {
                std::cerr << "Failed to reduce dataset for bin " << histIndex << std::endl;
                continue;
            }
            
            // Build histogram in the specific cut ranges for Det and Gen pT
            std::string histName = "hist" + std::to_string(histIndex);
            TH2F* hist = dynamic_cast<TH2F*>(dataInBin->createHistogram(histName.c_str(), zTDet, 
                                                RooFit::Binning(tbinsDet), 
                                                RooFit::YVar(zTPart, RooFit::Binning(tbinsPart))));
            
            if (!hist) {
                std::cerr << "Failed to create histogram for bin " << histIndex << std::endl;
                delete dataInBin;
                continue;
            }
            
            hh_data[indexDet][indexGenerator] = hist;
            maxArray[histIndex] = hist->GetMaximum();
            std::cout << "Integral for bin [" << indexDet << ", " << indexGenerator << "]: " 
                      << hist->Integral() << std::endl;
            
            histIndex++;
            delete dataInBin; // Clean up temporary dataset
        }
    }
    
    // Determine a global range for all plots
    double globalMax = 0;
    double globalMin = 1e6;
    
    // Find maximum value
    for (double max : maxArray) {
        if (max > globalMax) {
            globalMax = max;
        }
    }
    
    // Find minimum non-zero value
    for (double value : maxArray) {
        if (value > 0 && value < globalMin) {
            globalMin = value;
        }
    }
    
    // If no positive values found, set a default minimum
    if (globalMin == 1e6) {
        std::cout << "Warning: No positive values found in response matrix. Using default minimum." << std::endl;
        globalMin = 1e-6;
    }
    
    // Define margins
    double lMargin = 0.08;
    double rMargin = 0.03;
    double bMargin = 0.07;
    double tMargin = 0.03;
    
    // Set style
    gStyle->SetOptStat(0);
    
    // Create canvas
    TCanvas* canvas3 = new TCanvas("C3", "C3", 1600, 1600);
    canvas3->SetFillStyle(4000);
    
    // Create partitioned canvas
    TPad*** pad = partitionCanvas(canvas3, nBinsDetLvlv, nBinsGenLvlv, lMargin, rMargin, bMargin, tMargin);
    if (!pad) {
        std::cerr << "Failed to partition canvas" << std::endl;
        delete canvas3;
        delete reducedData;
        return;
    }
    
    // Fill each pad with the corresponding histogram
    int panelNumber = 0;
    for (int i = 0; i < nBinsDetLvlv; ++i) {  // Jet pT Detector
        for (int j = 0; j < nBinsGenLvlv; ++j) {  // Jet pT Generator
            canvas3->cd(0);
            
            // Get the pad
            std::string pname = "pad_" + std::to_string(i) + "_" + std::to_string(j);
            std::cout << pname << std::endl;
            
            if (!pad[i][j]) {
                std::cerr << "Pad " << i << "," << j << " not found" << std::endl;
                continue;
            }
            
            pad[i][j] = static_cast<TPad*>(gROOT->FindObject(pname.c_str()));
            pad[i][j]->Draw();
            pad[i][j]->SetFillStyle(4000);
            pad[i][j]->SetFrameFillStyle(4000);
            pad[i][j]->SetLogz();
            pad[i][j]->cd();
            
            double xFactor = pad[0][0]->GetAbsWNDC() / pad[i][j]->GetAbsWNDC();
            double yFactor = pad[0][0]->GetAbsHNDC() / pad[i][j]->GetAbsHNDC();
            
            // Invert the y axis for the display
            int invJ = nBinsGenLvlv - j - 1;
            
            if (hh_data[i][invJ]) {
                std::cout << "Name: " << hh_data[i][invJ]->GetName() << std::endl;
                hh_data[i][invJ]->SetTitle("");
                
                // Set frame format
                TH1* hFrame = setFrame(hh_data[i][invJ], xFactor, yFactor, i, nBinsDetLvlv, j, nBinsGenLvlv);
                hFrame->GetZaxis()->SetRangeUser(globalMin, globalMax);
                hFrame->GetXaxis()->SetTitle("z_{T}^{Det. lvl}");
                hFrame->GetYaxis()->SetTitle("z_{T}^{Gen. lvl}");
                hFrame->Draw();
                
                hh_data[i][invJ]->Draw("same colz");
                
                double topM = pad[i][j]->GetTopMargin();
                double leftM = pad[i][j]->GetLeftMargin();
                
                // Add a text label describing the bin
                systemArray[i][invJ] = new TLatex(0.16, 0.8, Form("#splitline{%s}{%s}", 
                                                legGen[invJ].c_str(), legDet[i].c_str()));
                systemArray[i][invJ]->SetTextSize(0.09);
                systemArray[i][invJ]->Draw();
                
                // Draw diagonal line
                if (i == invJ) {
                    diagonal_lineDark->Draw();
                } else {
                    diagonal_line->Draw();
                }
            } else {
                std::cout << "Error: could not find histogram at position: " << panelNumber << std::endl;
            }
            
            panelNumber++;
        }
    }
    
    // Save the canvas
    canvas3->cd();
    std::string ptTag = "pT5-70";
    std::string fileName = output + dictKey + "zTBin" + ptTag + "." + fileFormat;
    canvas3->SaveAs(fileName.c_str());
    
    // Clean up
    for (int i = 0; i < nBinsDetLvlv; ++i) {
        for (int j = 0; j < nBinsGenLvlv; ++j) {
            delete hh_data[i][j];
            delete systemArray[i][j];
            delete pad[i][j];
        }
        delete[] pad[i];
    }
    delete[] pad;
    delete canvas3;
    delete reducedData;
}

void RMplotter::plotRMRapidity(double dR) {
    // Get reduced dataset for entire jet pT range
    RooDataSet* reducedData = reduceDataSet(0, 200, dR);
    if (!reducedData) {
        std::cerr << "Failed to create reduced dataset for plotRMRapidity" << std::endl;
        return;
    }
    
    // Define variables for response matrix
    RooRealVar pt_jetDet("jet_pt_det", "jet_pt_det", 0, 200);
    RooRealVar pt_jetPart("jet_pt_mc", "jet_pt_mc", 0, 200);
    RooRealVar d0EtaDet("d0_eta_det", "d0_eta_det", 1.5, 5.5);
    RooRealVar d0EtaPart("d0_eta_mc", "d0_eta_mc", 1.5, 5.5);
    
    // Create custom binning for D0 rapidity in 0.02 intervals
    std::vector<double> etaBinning = {1.50, 1.52, 1.54, 1.56, 1.58, 1.60, 
                                      1.62, 1.64, 1.66, 1.68, 1.70, 
                                      1.72, 1.74, 1.76, 1.78, 1.80, 
                                      1.82, 1.84, 1.86, 1.88, 1.90,
                                      1.92, 1.94, 1.96, 1.98, 2.00,
                                      2.02, 2.04, 2.06, 2.08, 2.10,
                                      2.12, 2.14, 2.16, 2.18, 2.20,
                                      2.22, 2.24, 2.26, 2.28, 2.30,
                                      2.32, 2.34, 2.36, 2.38, 2.40,
                                      2.42, 2.44, 2.46, 2.48, 2.50,
                                      2.52, 2.54, 2.56, 2.58, 2.60,
                                      2.62, 2.64, 2.66, 2.68, 2.70,
                                      2.72, 2.74, 2.76, 2.78, 2.80,
                                      2.82, 2.84, 2.86, 2.88, 2.90,
                                      2.92, 2.94, 2.96, 2.98, 3.00,
                                      3.02, 3.04, 3.06, 3.08, 3.10,
                                      3.12, 3.14, 3.16, 3.18, 3.20,
                                      3.22, 3.24, 3.26, 3.28, 3.30,
                                      3.32, 3.34, 3.36, 3.38, 3.40,
                                      3.42, 3.44, 3.46, 3.48, 3.50,
                                      3.52, 3.54, 3.56, 3.58, 3.60,
                                      3.62, 3.64, 3.66, 3.68, 3.70,
                                      3.72, 3.74, 3.76, 3.78, 3.80,
                                      3.82, 3.84, 3.86, 3.88, 3.90,
                                      3.92, 3.94, 3.96, 3.98, 4.00,
                                      4.02, 4.04, 4.06, 4.08, 4.10,
                                      4.12, 4.14, 4.16, 4.18, 4.20,
                                      4.22, 4.24, 4.26, 4.28, 4.30,
                                      4.32, 4.34, 4.36, 4.38, 4.40,
                                      4.42, 4.44, 4.46, 4.48, 4.50};
                                    
    // Create ROOT binning objects
    RooBinning etaBinsDet(etaBinning.size()-1, &etaBinning[0], "etaBinningDet");
    RooBinning etaBinsPart(etaBinning.size()-1, &etaBinning[0], "etaBinningPart");
    
    // Calculate number of bins for generator and detector levels
    int nBinsGenLvlv = pTBinning.size() - 1;
    int nBinsDetLvlv = pTBinning.size() - 1;
    int maxCombinations = nBinsGenLvlv * nBinsDetLvlv;
    
    std::cout << "Max No combinations for rapidity: " << maxCombinations << std::endl;
    std::cout << "pT Part bins: " << nBinsGenLvlv << std::endl;
    std::cout << "pT Det bins: " << nBinsDetLvlv << std::endl;
    
    // Create arrays to store histograms and metadata
    std::vector<std::vector<TH2F*>> hh_data(nBinsDetLvlv, std::vector<TH2F*>(nBinsGenLvlv, nullptr));
    std::vector<double> maxArray(maxCombinations, 0);
    std::vector<std::string> legGen(nBinsGenLvlv);
    std::vector<std::string> legDet(nBinsDetLvlv);
    std::vector<std::vector<TLatex*>> systemArray(nBinsDetLvlv, std::vector<TLatex*>(nBinsGenLvlv, nullptr));
    
    // Create all histograms for the D0 rapidity responses
    int histIndex = 0;
    for (int indexGenerator = 0; indexGenerator < nBinsGenLvlv; ++indexGenerator) {
        // Define generator level pT range
        std::string ptRangePart = "jet_pt_mc > " + std::to_string(pTBinning[indexGenerator]) + 
                                " && jet_pt_mc < " + std::to_string(pTBinning[indexGenerator+1]);
        
        // Build legend entry for generator level
        legGen[indexGenerator] = Form("p_{T}^{jet,gen}: %1.0f-%1.0f", 
                                    pTBinning[indexGenerator], pTBinning[indexGenerator+1]);
        
        for (int indexDet = 0; indexDet < nBinsDetLvlv; ++indexDet) {
            // Define detector level pT range
            std::string ptRangeDet = " && jet_pt_det > " + std::to_string(pTBinning[indexDet]) +
                                   " && jet_pt_det < " + std::to_string(pTBinning[indexDet+1]);
            std::string cuts = ptRangePart + ptRangeDet;
            
            // Build legend entry for detector level
            legDet[indexDet] = Form("p_{T}^{jet,det}: %1.0f-%1.0f", 
                                  pTBinning[indexDet], pTBinning[indexDet+1]);
            
            // Get data in this specific bin
            RooDataSet* dataInBin = dynamic_cast<RooDataSet*>(reducedData->reduce(RooArgSet(d0EtaDet, d0EtaPart, pt_jetDet, pt_jetPart), cuts.c_str()));
            if (!dataInBin) {
                std::cerr << "Failed to reduce dataset for rapidity bin " << histIndex << std::endl;
                continue;
            }
            
            // Build histogram in the specific cut ranges for Det and Gen pT
            std::string histName = "histRapidity" + std::to_string(histIndex);
            TH2F* hist = dynamic_cast<TH2F*>(dataInBin->createHistogram(histName.c_str(), d0EtaDet, 
                                                RooFit::Binning(etaBinsDet), 
                                                RooFit::YVar(d0EtaPart, RooFit::Binning(etaBinsPart))));
            
            if (!hist) {
                std::cerr << "Failed to create rapidity histogram for bin " << histIndex << std::endl;
                delete dataInBin;
                continue;
            }
            
            hh_data[indexDet][indexGenerator] = hist;
            maxArray[histIndex] = hist->GetMaximum();
            std::cout << "Rapidity integral for bin [" << indexDet << ", " << indexGenerator << "]: " 
                      << hist->Integral() << std::endl;
            
            histIndex++;
            delete dataInBin; // Clean up temporary dataset
        }
    }
    
    // Determine a global range for all plots
    double globalMax = 0;
    double globalMin = 1e6;
    
    // Find maximum value
    for (double max : maxArray) {
        if (max > globalMax) {
            globalMax = max;
        }
    }
    
    // Find minimum non-zero value
    for (double value : maxArray) {
        if (value > 0 && value < globalMin) {
            globalMin = value;
        }
    }
    
    // If no positive values found, set a default minimum
    if (globalMin == 1e6) {
        std::cout << "Warning: No positive values found in rapidity response matrix. Using default minimum." << std::endl;
        globalMin = 1e-6;
    }
    
    // Define margins
    double lMargin = 0.08;
    double rMargin = 0.03;
    double bMargin = 0.07;
    double tMargin = 0.03;
    
    // Set style
    gStyle->SetOptStat(0);
    
    // Create canvas
    TCanvas* canvas4 = new TCanvas("C4", "C4", 1600, 1600);
    canvas4->SetFillStyle(4000);
    
    // Create partitioned canvas
    TPad*** pad = partitionCanvas(canvas4, nBinsDetLvlv, nBinsGenLvlv, lMargin, rMargin, bMargin, tMargin);
    if (!pad) {
        std::cerr << "Failed to partition canvas for rapidity" << std::endl;
        delete canvas4;
        delete reducedData;
        return;
    }
    
    // Fill each pad with the corresponding histogram
    int panelNumber = 0;
    for (int i = 0; i < nBinsDetLvlv; ++i) {  // Jet pT Detector
        for (int j = 0; j < nBinsGenLvlv; ++j) {  // Jet pT Generator
            canvas4->cd(0);
            
            // Get the pad
            std::string pname = "pad_" + std::to_string(i) + "_" + std::to_string(j);
            std::cout << "Rapidity " << pname << std::endl;
            
            if (!pad[i][j]) {
                std::cerr << "Rapidity pad " << i << "," << j << " not found" << std::endl;
                continue;
            }
            
            pad[i][j] = static_cast<TPad*>(gROOT->FindObject(pname.c_str()));
            pad[i][j]->Draw();
            pad[i][j]->SetFillStyle(4000);
            pad[i][j]->SetFrameFillStyle(4000);
            pad[i][j]->SetLogz();
            pad[i][j]->cd();
            
            double xFactor = pad[0][0]->GetAbsWNDC() / pad[i][j]->GetAbsWNDC();
            double yFactor = pad[0][0]->GetAbsHNDC() / pad[i][j]->GetAbsHNDC();
            
            // Invert the y axis for the display
            int invJ = nBinsGenLvlv - j - 1;
            
            if (hh_data[i][invJ]) {
                std::cout << "Rapidity name: " << hh_data[i][invJ]->GetName() << std::endl;
                hh_data[i][invJ]->SetTitle("");
                
                // Set frame format
                TH1* hFrame = setFrame(hh_data[i][invJ], xFactor, yFactor, i, nBinsDetLvlv, j, nBinsGenLvlv);
                hFrame->GetZaxis()->SetRangeUser(globalMin, globalMax);
                hFrame->GetXaxis()->SetTitle("#eta_{D0}^{Det. lvl}");
                hFrame->GetYaxis()->SetTitle("#eta_{D0}^{Gen. lvl}");
                hFrame->Draw();
                
                hh_data[i][invJ]->Draw("same colz");
                
                double topM = pad[i][j]->GetTopMargin();
                double leftM = pad[i][j]->GetLeftMargin();
                
                // Add a text label describing the bin
                systemArray[i][invJ] = new TLatex(0.16, 0.8, Form("#splitline{%s}{%s}", 
                                                legGen[invJ].c_str(), legDet[i].c_str()));
                systemArray[i][invJ]->SetTextSize(0.09);
                systemArray[i][invJ]->Draw();
                
                // Draw diagonal line
                if (i == invJ) {
                    diagonal_lineDark->Draw();
                } else {
                    diagonal_line->Draw();
                }
            } else {
                std::cout << "Error: could not find rapidity histogram at position: " << panelNumber << std::endl;
            }
            
            panelNumber++;
        }
    }
    
    // Save the canvas
    canvas4->cd();
    std::string ptTag = "pT5-70";
    std::string fileName = output + dictKey + "RapidityBin" + ptTag + "." + fileFormat;
    canvas4->SaveAs(fileName.c_str());
    
    // Clean up
    for (int i = 0; i < nBinsDetLvlv; ++i) {
        for (int j = 0; j < nBinsGenLvlv; ++j) {
            delete hh_data[i][j];
            delete systemArray[i][j];
            delete pad[i][j];
        }
        delete[] pad[i];
    }
    delete[] pad;
    delete canvas4;
    delete reducedData;
}

TPad*** RMplotter::partitionCanvas(TCanvas* C, int Nx, int Ny, double lMargin, double rMargin, double bMargin, double tMargin) {
    if (!C) {
        return nullptr;
    }
    
    // Setup Pad layout
    // Vertical space between pads
    double vSpacing = 0.0;
    // Vertical size of a pad
    double vStep = (1.0 - bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
    
    // Horizontal space between pads
    double hSpacing = 0.0;
    // Horizontal size of a pad
    double hStep = (1.0 - lMargin - rMargin - (Nx-1) * hSpacing) / Nx;
    
    // Allocate 2D array of pads
    TPad*** pad = new TPad**[Nx];
    for (int i = 0; i < Nx; ++i) {
        pad[i] = new TPad*[Ny];
        for (int j = 0; j < Ny; ++j) {
            pad[i][j] = nullptr;
        }
    }
    
    // Column loop
    for (int i = 0; i < Nx; ++i) {
        double hposl, hposr, hfactor, hmarl, hmarr;
        
        if (i == 0) {
            hposl = 0.0;
            hposr = lMargin + hStep;
            hfactor = hposr - hposl;
            hmarl = lMargin / hfactor;
            hmarr = 0.0;
        } else if (i == Nx-1) {
            hposl = hposr + hSpacing;
            hposr = hposl + hStep + rMargin;
            hfactor = hposr - hposl;
            hmarl = 0.0;
            hmarr = rMargin / hfactor;
        } else {
            hposl = hposr + hSpacing;
            hposr = hposl + hStep;
            hmarl = 0.0;
            hmarr = 0.0;
        }
        
        // Row loop
        for (int j = 0; j < Ny; ++j) {
            double vposu, vposd, vfactor, vmard, vmaru;
            
            if (j == 0) {
                vposu = 1.0;
                vposd = 1.0 - vStep - tMargin;
                vfactor = vposu - vposd;
                vmard = 0.0;
                vmaru = tMargin / vfactor;
            } else if (j == Ny-1) {
                vposu = vposd - vSpacing;
                vposd = vposu - vStep - bMargin;
                vfactor = vposu - vposd;
                vmard = bMargin / vfactor;
                vmaru = 0.0;
            } else {
                vposu = vposd - vSpacing;
                vposd = vposu - vStep;
                vmard = 0.0;
                vmaru = 0.0;
            }
            
            if (vposd < 0) vposd = 0;
            
            C->cd(0);
            std::string name = "pad_" + std::to_string(i) + "_" + std::to_string(j);
            std::cout << "Creating pad: " << name << std::endl;
            std::cout << "Position: left=" << hposl << ", right=" << hposr 
                     << ", top=" << vposu << ", bottom=" << vposd << std::endl;
                     
            pad[i][j] = new TPad(name.c_str(), "", hposl, vposd, hposr, vposu);
            pad[i][j]->SetLeftMargin(hmarl);
            pad[i][j]->SetRightMargin(hmarr);
            pad[i][j]->SetBottomMargin(vmard);
            pad[i][j]->SetTopMargin(vmaru);
            pad[i][j]->SetFrameBorderMode(0);
            pad[i][j]->SetBorderMode(0);
            pad[i][j]->SetBorderSize(0);
            pad[i][j]->Draw();
        }
    }
    
    return pad;
}

TH1* RMplotter::setFrame(TH1* h, double xFactor, double yFactor, int i, int Nx, int j, int Ny) {
    std::string fname = "frame_" + std::string(h->GetName());
    TH1* hFrame = static_cast<TH1*>(h->Clone(fname.c_str()));
    
    // Empty histogram
    hFrame->Reset();
    
    // Set axis ranges
    hFrame->GetXaxis()->SetRangeUser(0, 1);
    hFrame->GetYaxis()->SetRangeUser(0, 1);
    
    // Format for y axis
    hFrame->GetYaxis()->SetLabelFont(43);
    hFrame->GetYaxis()->SetLabelSize(32);  // was 16
    hFrame->GetYaxis()->SetLabelOffset(0.01);  // was 0.2
    hFrame->GetYaxis()->SetTitleFont(43);
    hFrame->GetYaxis()->SetTitleSize(35);  // was 16
    hFrame->GetYaxis()->SetTitleOffset(5);  // was 5
    hFrame->GetYaxis()->CenterTitle();
    hFrame->GetYaxis()->SetNdivisions(505);
    
    // TICKS Y Axis
    hFrame->GetYaxis()->SetTickLength(xFactor * 0.04 / yFactor);
    
    // Format for x axis
    hFrame->GetXaxis()->SetLabelFont(43);
    hFrame->GetXaxis()->SetLabelSize(32);  // was 16
    hFrame->GetXaxis()->SetLabelOffset(0.01);
    hFrame->GetXaxis()->SetTitleFont(43);
    hFrame->GetXaxis()->SetTitleSize(38);  // was 16
    hFrame->GetXaxis()->SetTitleOffset(5);  // was 5
    hFrame->GetXaxis()->CenterTitle();
    hFrame->GetXaxis()->SetNdivisions(505);
    
    // TICKS X Axis
    hFrame->GetXaxis()->SetTickLength(yFactor * 0.06 / xFactor);
    
    // Z axis format
    hFrame->GetZaxis()->SetLabelFont(43);
    hFrame->GetZaxis()->SetLabelSize(32);
    hFrame->GetZaxis()->SetTitle("");
    
    return hFrame;
}

void RMplotter::setHisto(TH1* Histo, const std::string& Xtitle, const std::string& Ytitle, int big, double border) {
    Histo->SetStats(0);
    Histo->SetTitle("");
    
    if (big == 0) {
        Histo->GetYaxis()->SetTitleOffset(1.4);
        Histo->GetXaxis()->SetTitleOffset(1.4);
        Histo->GetXaxis()->SetLabelSize(0.05);
        Histo->GetYaxis()->SetLabelSize(0.05);
        Histo->GetXaxis()->SetTitleSize(0.045);
        Histo->GetYaxis()->SetTitleSize(0.045);
    } else if (big == 1) {
        Histo->GetYaxis()->SetTitleOffset(1.0);
        Histo->GetXaxis()->SetTitleOffset(0.82);
        Histo->GetYaxis()->SetLabelSize(0.06);
        Histo->GetXaxis()->SetLabelSize(0.06);
        Histo->GetYaxis()->SetTitleSize(0.07);
        Histo->GetXaxis()->SetTitleSize(0.07);
    } else if (big == 2) {
        Histo->GetYaxis()->SetTitleOffset(1.2);
        Histo->GetXaxis()->SetTitleOffset(1.2);
        Histo->GetXaxis()->SetLabelSize(0.05);
        Histo->GetYaxis()->SetLabelSize(0.05);
        Histo->GetXaxis()->SetTitleSize(0.055);
        Histo->GetYaxis()->SetTitleSize(0.055);
    }
    
    Histo->GetXaxis()->SetNdivisions(505);
    Histo->GetYaxis()->SetNdivisions(505);
    
    // Make nice font
    Histo->GetXaxis()->SetLabelFont(42);
    Histo->GetYaxis()->SetLabelFont(42);
    Histo->GetXaxis()->SetTitleFont(42);
    Histo->GetYaxis()->SetTitleFont(42);
    
    if (Xtitle != "") {
        Histo->GetXaxis()->SetTitle(Xtitle.c_str());
    }
    if (Ytitle != "") {
        Histo->GetYaxis()->SetTitle(Ytitle.c_str());
    }
    
    // Check whether this is TH1 or TH2
    if (Histo->InheritsFrom(TH2::Class())) {
        return;
    } else {
        setTH1Histo(Histo, Xtitle, Ytitle, big, border);
    }
}

void RMplotter::setTH1Histo(TH1* Histo, const std::string& Xtitle, const std::string& Ytitle, int big, double border) {
    // Set a larger y-axis range to leave space for the border
    double min = Histo->GetMinimum(0);
    double max = Histo->GetBinContent(Histo->GetMaximumBin());
    double range = max - min;
    
    double maxNew, minNew;
    if (range > 0) {
        maxNew = min + range * (1.0 + 3 * border);
        minNew = max - range * (1.0 + border);
    } else {
        maxNew = max + (-1) * range * (1.0 + 2 * border);
        minNew = min - (-1) * range * (1.0 + border);
    }
    
    minNew = min;  // Override the calculated value
    Histo->GetYaxis()->SetRangeUser(minNew, maxNew);
    
    Histo->SetLineColor(1);
    Histo->SetMarkerColor(1);
    Histo->SetMarkerStyle(20);
    Histo->SetMarkerSize(0.7);
}

void plotResponse() {
    // Enable batch mode
    gROOT->SetBatch(true);
    
    try {
        // Create the plotter object
        RMplotter* plotterProcess = new RMplotter();
        
        // Call plot methods
        /*
        plotterProcess->plotRMJet(5, 10);
        plotterProcess->plotRMJet(10, 15);
        plotterProcess->plotRMJet(15, 20);
        plotterProcess->plotRMJet(20, 30);
        plotterProcess->plotRMJet(30, 40);
        plotterProcess->plotRMJet(40, 70);
        plotterProcess->plotRMJet(5, 70);
        */
        
        // Plot the tag response matrix
        plotterProcess->plotRMTag();
        
        // Plot the rapidity response matrix
        plotterProcess->plotRMRapidity();
        
        // Clean up
        delete plotterProcess;
        
    } catch (const std::exception& e) {
        std::cerr << "Error in plotResponse: " << e.what() << std::endl;
        return;
    }
}
