#ifndef TRACKRECOEFFICIENCY_H
#define TRACKRECOEFFICIENCY_H

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
#include <atomic>
#include <mutex>

/**
 * @class TrackRecoEfficiency
 * @brief Standalone C++ implementation for calculating track reconstruction efficiency
 *
 * This class reads MC ROOT trees and calculates track reconstruction efficiency
 * for pions and kaons (excluding PID efficiency - only geometric and tracking efficiency).
 *
 * The efficiency is calculated as:
 * - Numerator: MC truth particles that have a reconstructed track match
 * - Denominator: All MC truth particles in detector acceptance
 */
class TrackRecoEfficiency
{
private:
    // Input/Output files
    TFile *m_inputFile;
    TFile *m_outputFile;
    TTree *m_tree;

    // Tree variables for MC truth particles
    std::vector<int> *mc_particle_pid;
    std::vector<float> *mc_particle_pt;
    std::vector<float> *mc_particle_eta;
    std::vector<float> *mc_particle_phi;
    std::vector<float> *mc_particle_px;
    std::vector<float> *mc_particle_py;
    std::vector<float> *mc_particle_pz;
    std::vector<float> *mc_particle_e;
    std::vector<int> *mc_particle_matched;    // Index of matched reconstructed track (-1 if no match)
    std::vector<int> *mc_particle_fromD0;    // Flag indicating if particle comes from D0 decay

    // Tree variables for reconstructed tracks
    Int_t n_tracks;
    std::vector<float> *track_pt;
    std::vector<float> *track_eta;
    std::vector<float> *track_phi;
    std::vector<float> *track_px;
    std::vector<float> *track_py;
    std::vector<float> *track_pz;
    std::vector<int> *track_pid;              // Reconstructed PID (not used for efficiency, just for validation)
    std::vector<float> *track_prb_ghost;     // Ghost probability
    std::vector<float> *track_chi2;          // Track chi2
    std::vector<int> *track_mc_matched;      // Index of matched MC particle (-1 if no match)

    // Binning for efficiency calculations
    std::vector<double> m_ptBins;
    std::vector<double> m_etaBins;
    std::vector<double> m_pBins;

    // Efficiency histograms and objects
    std::map<std::string, TH2F *> m_efficiencyMaps;
    std::map<std::string, TEfficiency *> m_efficiencyObjects;
    std::map<std::string, TH1F *> m_1DHistograms;

    // Configuration parameters
    double m_minPt;               // Minimum pT (GeV)
    double m_minEta;              // Minimum eta
    double m_maxEta;              // Maximum eta
    double m_maxGhostProb;        // Maximum ghost probability for quality cuts
    double m_maxTrackChi2;        // Maximum track chi2 for quality cuts
    bool m_requireFromD0;         // Only consider particles from D0 decay
    bool m_useParallel;           // Use parallel processing
    int m_nThreads;               // Number of threads

public:
    TrackRecoEfficiency(TString inputFileName, TString outputFileName);
    ~TrackRecoEfficiency();

    // Configuration setters
    void SetPtRange(double minPt) { m_minPt = minPt; }
    void SetEtaRange(double minEta, double maxEta) { m_minEta = minEta; m_maxEta = maxEta; }
    void SetQualityCuts(double maxGhostProb, double maxChi2) { m_maxGhostProb = maxGhostProb; m_maxTrackChi2 = maxChi2; }
    void SetRequireFromD0(bool require) { m_requireFromD0 = require; }
    void SetUseParallel(bool useParallel) { m_useParallel = useParallel; }
    void SetNumThreads(int nThreads) { m_nThreads = nThreads; }

    // Binning setup
    void SetPtBins(const std::vector<double> &ptBins) { m_ptBins = ptBins; }
    void SetEtaBins(const std::vector<double> &etaBins) { m_etaBins = etaBins; }
    void SetPBins(const std::vector<double> &pBins) { m_pBins = pBins; }

    // Main calculation methods
    bool Initialize();
    void CalculateTrackEfficiency();
    void CalculateTrackEfficiencyParallel();
    void SaveResults();
    void PlotEfficiency(const std::string &particleType = "pion");

    // Helper methods
    bool IsInAcceptance(double pt, double eta);
    bool PassesQualityCuts(int track_idx);
    double CalculateMomentum(double px, double py, double pz);
    std::string GetParticleName(int pid);

private:
    void InitializeBranches();
    void CreateHistograms();
    void CleanUp();
};

#endif // TRACKRECOEFFICIENCY_H
