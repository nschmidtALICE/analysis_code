#ifndef D0EFFICIENCY_H
#define D0EFFICIENCY_H

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>

/**
 * @class D0Efficiency
 * @brief Basic C++ implementation for calculating D0 PID and reconstruction efficiencies
 * 
 * This class reads the ROOT tree produced by main_analysis.py and calculates:
 * 1. PID efficiency: How well particle identification works for K and π from D0 decay
 * 2. Reconstruction efficiency: How well the detector reconstructs true D0 particles
 * 
 * The efficiency is calculated in bins of pT and η for systematic studies.
 */
class D0Efficiency {
private:
    // Input/Output files
    TFile* m_inputFile;
    TFile* m_outputFile;
    TTree* m_tree;
    
    // Tree variables (matching main_analysis.py structure)
    Int_t n_d0s;
    std::vector<float>* d0_pt;
    std::vector<float>* d0_eta;
    std::vector<float>* d0_phi;
    std::vector<float>* d0_mass;
    std::vector<float>* d0_px;
    std::vector<float>* d0_py;
    std::vector<float>* d0_pz;
    std::vector<float>* d0_e;
    
    // D0 daughter variables
    Int_t n_daughters;
    std::vector<int>* dau_pid;
    std::vector<float>* dau_pt;
    std::vector<float>* dau_eta;
    std::vector<float>* dau_phi;
    std::vector<float>* dau_px;
    std::vector<float>* dau_py;
    std::vector<float>* dau_pz;
    std::vector<float>* dau_e;
    std::vector<int>* dau_d0_idx;
    std::vector<float>* dau_pnn_k;    // ProbNNK for kaon PID
    std::vector<float>* dau_pnn_pi;   // ProbNNπ for pion PID
    std::vector<float>* dau_prb_ghost; // Ghost probability
    std::vector<float>* dau_trckChi2;  // Track chi2
    
    // MC truth variables (only available for MC)
    std::vector<int>* mc_d0_pid;
    std::vector<float>* mc_d0_pt;
    std::vector<float>* mc_d0_eta;
    std::vector<float>* mc_d0_phi;
    std::vector<float>* mc_d0_mass;
    std::vector<float>* mc_d0_px;
    std::vector<float>* mc_d0_py;
    std::vector<float>* mc_d0_pz;
    std::vector<float>* mc_d0_e;
    std::vector<int>* mc_d0_origin;
    std::vector<int>* mc_d0_matched;
    std::vector<float>* mc_d0_matched_quality;
    
    // MC daughter variables
    Int_t mc_n_daughters;
    std::vector<int>* mc_dau_pid;
    std::vector<float>* mc_dau_pt;
    std::vector<float>* mc_dau_eta;
    std::vector<float>* mc_dau_phi;
    std::vector<int>* mc_dau_d0_idx;
    
    // Binning for efficiency calculations
    std::vector<double> m_ptBins;
    std::vector<double> m_etaBins;
    std::vector<double> m_pBins;
    
    // Efficiency histograms
    std::map<std::string, TH2F*> m_efficiencyMaps;
    std::map<std::string, TH1F*> m_efficiencyProjections;
    std::map<std::string, TF1*> m_efficiencyFunctions;
    std::map<std::string, TEfficiency*> m_efficiencyObjects;
    
    // Configuration
    bool m_isMC;
    std::string m_efficiencyType;
    double m_d0MassWindow;
    double m_minPt;
    double m_minEta;
    double m_maxEta;
    
    // PID cuts
    double m_kaonPIDCut;
    double m_pionPIDCut;
    double m_ghostProbCut;
    double m_trackChi2Cut;
    
    // Daughter particle cuts
    double m_minDaughterMomentum;
    
public:
    /**
     * @brief Constructor
     * @param inputFileName Path to input ROOT file from main_analysis.py
     * @param outputFileName Path to output ROOT file for efficiency results
     * @param efficiencyType Type of efficiency to calculate ("PID" or "reco")
     * @param isMC Whether input file is MC (true) or data (false)
     */
    D0Efficiency(const std::string& inputFileName, 
                 const std::string& outputFileName,
                 const std::string& efficiencyType = "PID",
                 bool isMC = true);
    
    /**
     * @brief Destructor
     */
    ~D0Efficiency();
    
    /**
     * @brief Initialize the analysis
     * @return true if successful, false otherwise
     */
    bool Initialize();
    
    /**
     * @brief Set binning for efficiency calculations
     * @param ptBins Vector of pT bin edges [GeV]
     * @param etaBins Vector of η bin edges
     */
    void SetBinning(const std::vector<double>& ptBins, const std::vector<double>& etaBins);
    
    /**
     * @brief Set binning for efficiency calculations including momentum
     * @param ptBins Vector of pT bin edges [GeV]
     * @param etaBins Vector of η bin edges
     * @param pBins Vector of p bin edges [GeV]
     */
    void SetBinning(const std::vector<double>& ptBins, const std::vector<double>& etaBins, const std::vector<double>& pBins);
    
    /**
     * @brief Set PID selection criteria
     * @param kaonPIDCut Minimum ProbNNK for kaon identification
     * @param pionPIDCut Minimum ProbNNπ for pion identification
     * @param ghostProbCut Maximum ghost probability
     * @param trackChi2Cut Maximum track chi2/ndf
     */
    void SetPIDCuts(double kaonPIDCut = 0.5, double pionPIDCut = 0.5, 
                    double ghostProbCut = 0.3, double trackChi2Cut = 3.0);
    
    /**
     * @brief Set D0 selection criteria
     * @param massWindow D0 mass window around PDG mass [MeV]
     * @param minPt Minimum D0 pT [GeV]
     * @param minEta Minimum |η| for D0
     * @param maxEta Maximum |η| for D0
     */
    void SetD0Selection(double massWindow = 50.0, double minPt = 2.0, double minEta = 2.0, double maxEta = 5.0);
    
    /**
     * @brief Set daughter particle selection criteria
     * @param minMomentum Minimum momentum for D0 daughters [GeV]
     */
    void SetDaughterCuts(double minMomentum = 2.0);
    
    /**
     * @brief Calculate kaon PID efficiency separately
     * Numerator: D0 candidates with kaon passing PID cuts (regardless of pion)
     * Denominator: All D0 candidates passing basic quality cuts
     */
    void CalculateKaonPIDEfficiency();
    
    /**
     * @brief Calculate pion PID efficiency separately
     * Numerator: D0 candidates with pion passing PID cuts (regardless of kaon)
     * Denominator: All D0 candidates passing basic quality cuts
     */
    void CalculatePionPIDEfficiency();

    /**
     * @brief Calculate D0 PID efficiency (both kaon AND pion pass)
     * Numerator: D0 candidates with BOTH daughters passing PID cuts
     * Denominator: All D0 candidates passing basic quality cuts
     */
    void CalculatePIDEfficiency();
    
    /**
     * @brief Calculate combined D0 PID efficiency (both kaon AND pion pass)
     * Numerator: D0 candidates with BOTH daughters passing PID cuts
     * Denominator: All D0 candidates passing basic quality cuts
     */
    void CalculateCombinedPIDEfficiency();
    
    /**
     * @brief Calculate reconstruction efficiency (MC only)
     * Numerator: MC D0 particles that are reconstructed
     * Denominator: All MC D0 particles in acceptance
     */
    void CalculateRecoEfficiency();
    
    /**
     * @brief Parameterize efficiency as a function of pT for each η bin
     * @param efficiencyName Name of the efficiency histogram to parameterize
     * @param functionType Type of function to use ("exp", "pol2", "erf")
     */
    void ParameterizeEfficiency(const std::string& efficiencyName, 
                               const std::string& functionType = "exp");
    
    /**
     * @brief Get efficiency weight for a given D0 candidate
     * @param pt D0 pT [GeV]
     * @param eta D0 η
     * @param efficiencyName Name of the efficiency to use
     * @param useFunction Use parameterized function (true) or histogram (false)
     * @return Efficiency weight (1/efficiency)
     */
    double GetEfficiencyWeight(double pt, double eta, const std::string& efficiencyName, 
                              bool useFunction = true);
    
    /**
     * @brief Get TEfficiency object directly
     * @param efficiencyName Name of the efficiency object to retrieve
     * @return Pointer to TEfficiency object or nullptr if not found
     */
    TEfficiency* GetEfficiencyObject(const std::string& efficiencyName);
    
    /**
     * @brief Perform closure test on MC
     * Applies efficiency weights and checks if distributions are correctly reproduced
     */
    void PerformClosureTest();
    
    /**
     * @brief Save all histograms and results to output file
     */
    void SaveResults();
    
    /**
     * @brief Print efficiency summary
     */
    void PrintSummary();
    
private:
    /**
     * @brief Setup tree branches
     */
    void SetupBranches();
    
    /**
     * @brief Create efficiency histograms
     */
    void CreateHistograms();
    
    /**
     * @brief Check if D0 candidate passes basic selection
     * @param d0_idx Index of D0 candidate
     * @return true if passes selection
     */
    bool PassesD0Selection(int d0_idx);
    
    /**
     * @brief Check if kaon daughter passes PID selection
     * @param d0_idx Index of D0 candidate
     * @return true if kaon daughter passes PID cuts
     */
    bool PassesKaonPIDSelection(int d0_idx);
    
    /**
     * @brief Check if any daughter passes PID selection
     * @param d0_idx Index of D0 candidate
     * @return true if at least one daughter passes PID cuts
     */
    bool PassesPIDSelection(int d0_idx);

    /**
     * @brief Check if pion daughter passes PID selection
     * @param d0_idx Index of D0 candidate
     * @return true if pion daughter passes PID cuts
     */
    bool PassesPionPIDSelection(int d0_idx);
    
    /**
     * @brief Check if D0 daughters pass basic quality cuts
     * @param d0_idx Index of D0 candidate
     * @return true if both daughters pass quality cuts
     */
    bool PassesQualityCuts(int d0_idx);
    
    /**
     * @brief Find MC truth match for reconstructed D0
     * @param d0_idx Index of reconstructed D0
     * @return Index of matched MC D0, or -1 if no match
     */
    int FindMCMatch(int d0_idx);
    
    /**
     * @brief Calculate ΔR between two particles
     * @param eta1 η of first particle
     * @param phi1 φ of first particle
     * @param eta2 η of second particle
     * @param phi2 φ of second particle
     * @return ΔR
     */
    double CalculateDeltaR(double eta1, double phi1, double eta2, double phi2);
    
    /**
     * @brief Calculate total momentum of D0 candidate
     * @param d0_idx Index of D0 candidate
     * @return Total momentum [GeV]
     */
    double CalculateMomentum(int d0_idx);
    
    /**
     * @brief Calculate total momentum of MC D0 particle
     * @param mc_idx Index of MC D0 particle
     * @return Total momentum [GeV]
     */
    double CalculateMCMomentum(int mc_idx);
    
    /**
     * @brief Get bin indices for given pT and η
     * @param pt pT value [GeV]
     * @param eta η value
     * @param ptBin Output: pT bin index
     * @param etaBin Output: η bin index
     * @return true if within binning range
     */
    bool GetBinIndices(double pt, double eta, int& ptBin, int& etaBin);
};

#endif // D0EFFICIENCY_H
