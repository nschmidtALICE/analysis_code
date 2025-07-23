#ifndef D0RECOEFFICIENCY_H
#define D0RECOEFFICIENCY_H

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TROOT.h>

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <mutex>
#include <atomic>

/**
 * @class D0RecoEfficiency
 * @brief Standalone C++ implementation for calculating D0 reconstruction efficiency
 *
 * This class reads MC ROOT trees and calculates reconstruction efficiency:
 * How well the detector reconstructs true MC D0 particles.
 *
 * The efficiency is calculated in bins of pT and η for systematic studies.
 */
class D0RecoEfficiency
{
private:
    // Input/Output files
    TFile *m_inputFile;
    TFile *m_outputFile;
    TTree *m_tree;

    // Tree variables for reconstructed D0 particles
    Int_t n_d0s;
    std::vector<float> *d0_pt;
    std::vector<float> *d0_eta;
    std::vector<float> *d0_phi;
    std::vector<float> *d0_mass;
    std::vector<float> *d0_px;
    std::vector<float> *d0_py;
    std::vector<float> *d0_pz;
    std::vector<float> *d0_e;

    // D0 daughter variables for quality cuts
    Int_t n_daughters;
    std::vector<int> *dau_pid;
    std::vector<float> *dau_pt;
    std::vector<float> *dau_eta;
    std::vector<float> *dau_phi;
    std::vector<float> *dau_px;
    std::vector<float> *dau_py;
    std::vector<float> *dau_pz;
    std::vector<int> *dau_d0_idx;
    std::vector<float> *dau_pnn_k;     // ProbNNK for kaon PID
    std::vector<float> *dau_pnn_pi;    // ProbNNπ for pion PID
    std::vector<float> *dau_prb_ghost; // Ghost probability
    std::vector<float> *dau_trckChi2;  // Track chi2

    // MC truth variables
    std::vector<int> *mc_d0_pid;
    std::vector<float> *mc_d0_pt;
    std::vector<float> *mc_d0_eta;
    std::vector<float> *mc_d0_phi;
    std::vector<float> *mc_d0_px;
    std::vector<float> *mc_d0_py;
    std::vector<float> *mc_d0_pz;
    std::vector<int> *mc_d0_matched; // Index of matched reco D0 (-1 if no match)
    
    // MC truth variables for daughters
    std::vector<int> *mc_dau_pid;     // True PID of daughter tracks
    std::vector<int> *mc_dau_matched; // Index of matched MC particle (-1 if no match)

    // Binning for efficiency calculations
    std::vector<double> m_ptBins;
    std::vector<double> m_etaBins;
    std::vector<double> m_pBins;

    // Efficiency histograms
    std::map<std::string, TH2F *> m_efficiencyMaps;
    std::map<std::string, TEfficiency *> m_efficiencyObjects;

    // Configuration parameters
    double m_d0MassWindow;        // D0 mass window (MeV)
    double m_minPt;               // Minimum pT (GeV)
    double m_minEta;              // Minimum eta
    double m_maxEta;              // Maximum eta
    double m_kaonPIDCut;          // Minimum ProbNNK
    double m_pionPIDCut;          // Minimum ProbNNπ
    double m_ghostProbCut;        // Maximum ghost probability
    double m_trackChi2Cut;        // Maximum track chi2
    double m_minDaughterMomentum; // Minimum daughter momentum (GeV)

public:
    D0RecoEfficiency(TString inputFileName,
                     TString outputFileName);

    ~D0RecoEfficiency();

    // Configuration setters
    void SetD0MassWindow(double window) { m_d0MassWindow = window; }
    void SetPtRange(double minPt) { m_minPt = minPt; }
    void SetEtaRange(double minEta, double maxEta)
    {
        m_minEta = minEta;
        m_maxEta = maxEta;
    }
    void SetPIDCuts(double kaonCut, double pionCut)
    {
        m_kaonPIDCut = kaonCut;
        m_pionPIDCut = pionCut;
    }
    void SetQualityCuts(double ghostCut, double chi2Cut, double minMom)
    {
        m_ghostProbCut = ghostCut;
        m_trackChi2Cut = chi2Cut;
        m_minDaughterMomentum = minMom;
    }

    // Binning setup
    void SetPtBins(const std::vector<double> &ptBins) { m_ptBins = ptBins; }
    void SetEtaBins(const std::vector<double> &etaBins) { m_etaBins = etaBins; }
    void SetPBins(const std::vector<double> &pBins) { m_pBins = pBins; }

    // Main calculation methods
    bool Initialize();
    void CalculateRecoEfficiency();
    void SaveResults();
    void PlotEfficiency(const std::string &histName = "reco_efficiency");

    // Helper methods
    bool PassesD0Selection(int d0_idx);
    bool PassesPIDSelection(int d0_idx);
    bool PassesQualityCutsOptimized(int d0_idx); // Optimized version
    double CalculateMCMomentum(int mc_idx);

private:
    void InitializeBranches();
    void CreateHistograms();
    void CleanUp();
};

#endif // D0RECOEFFICIENCY_H
