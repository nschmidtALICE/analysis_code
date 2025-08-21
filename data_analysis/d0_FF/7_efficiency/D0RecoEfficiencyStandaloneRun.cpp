#include <iostream>
#include <vector>
#include <map>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TMath.h"
#include "TString.h"

// Standalone D0 Reconstruction Efficiency Calculator (Reco-based approach)
class D0RecoEfficiencyStandalone {
public:
    D0RecoEfficiencyStandalone(TString inputFileName, TString outputFileName);
    ~D0RecoEfficiencyStandalone();
    
    bool Initialize();
    void CalculateRecoEfficiency();
    void SaveResults();
    void PlotEfficiency(const std::string &histName);
    
    // Configuration methods
    void SetD0MassWindow(double window) { m_d0MassWindow = window; }
    void SetPtRange(double minPt) { m_minPt = minPt; }
    void SetEtaRange(double minEta, double maxEta) { m_minEta = minEta; m_maxEta = maxEta; }
    void SetPIDCuts(double kaonCut, double pionCut) { m_kaonPIDCut = kaonCut; m_pionPIDCut = pionCut; }
    void SetPtBins(const std::vector<double>& bins) { m_ptBins = bins; }
    void SetEtaBins(const std::vector<double>& bins) { m_etaBins = bins; }
    void SetPBins(const std::vector<double>& bins) { m_pBins = bins; }

private:
    // File and tree management
    TFile *m_inputFile;
    TFile *m_outputFile;
    TTree *m_tree;
    
    // Configuration parameters
    double m_d0MassWindow;
    double m_minPt;
    double m_minEta;
    double m_maxEta;
    double m_kaonPIDCut;
    double m_pionPIDCut;
    double m_ghostProbCut;
    double m_trackChi2Cut;
    double m_minDaughterMomentum;
    
    // Binning
    std::vector<double> m_ptBins;
    std::vector<double> m_etaBins;
    std::vector<double> m_pBins;
    
    // Histograms and efficiency objects
    std::map<std::string, TH2F*> m_efficiencyMaps;
    std::map<std::string, TEfficiency*> m_efficiencyObjects;
    
    // Tree branches
    Int_t n_d0s;
    std::vector<float> *d0_pt;
    std::vector<float> *d0_eta;
    std::vector<float> *d0_phi;
    std::vector<float> *d0_mass;
    std::vector<float> *d0_px;
    std::vector<float> *d0_py;
    std::vector<float> *d0_pz;
    std::vector<float> *d0_e;
    
    Int_t n_daughters;
    std::vector<int> *dau_pid;
    std::vector<float> *dau_pt;
    std::vector<float> *dau_eta;
    std::vector<float> *dau_phi;
    std::vector<float> *dau_px;
    std::vector<float> *dau_py;
    std::vector<float> *dau_pz;
    std::vector<int> *dau_d0_idx;
    std::vector<float> *dau_pnn_k;
    std::vector<float> *dau_pnn_pi;
    std::vector<float> *dau_prb_ghost;
    std::vector<float> *dau_trckChi2;
    
    std::vector<int> *mc_d0_pid;
    std::vector<float> *mc_d0_pt;
    std::vector<float> *mc_d0_eta;
    std::vector<float> *mc_d0_phi;
    std::vector<float> *mc_d0_px;
    std::vector<float> *mc_d0_py;
    std::vector<float> *mc_d0_pz;
    std::vector<int> *mc_d0_matched;
    std::vector<float> *mc_d0_matched_quality;
    
    std::vector<int> *mc_dau_pid;
    std::vector<int> *mc_dau_matched;
    
    // Helper methods
    void InitializeBranches();
    void CreateHistograms();
    void CleanUp();
    bool PassesD0Selection(int d0_idx);
    bool PassesQualityCuts(int d0_idx);
    double CalculateMCMomentum(int mc_idx);
    int GetMCMatchIndex(int reco_d0_idx);
    int GetMCMatchQuality(int mc_match_idx);
    bool HasValidMCMatch(int reco_d0_idx, int& mc_match_idx);
};

D0RecoEfficiencyStandalone::D0RecoEfficiencyStandalone(TString inputFileName, TString outputFileName)
    : m_inputFile(nullptr), m_outputFile(nullptr), m_tree(nullptr),
      m_d0MassWindow(50.0), m_minPt(2.0), m_minEta(2.0), m_maxEta(5.0),
      m_kaonPIDCut(0.5), m_pionPIDCut(0.5), m_ghostProbCut(0.3), m_trackChi2Cut(3.0),
      m_minDaughterMomentum(2.0)
{
    // Open input file
    m_inputFile = TFile::Open(inputFileName, "READ");
    if (!m_inputFile || m_inputFile->IsZombie()) {
        std::cerr << "Error: Could not open input file " << inputFileName << std::endl;
        return;
    }
    
    // Create output file
    m_outputFile = new TFile(outputFileName, "RECREATE");
    if (!m_outputFile || m_outputFile->IsZombie()) {
        std::cerr << "Error: Could not create output file " << outputFileName << std::endl;
        return;
    }
    
    // Get tree from input file
    m_tree = (TTree*)m_inputFile->Get("d0jets");
    if (!m_tree) {
        std::cerr << "Error: Could not find tree 'd0jets' in input file" << std::endl;
        return;
    }
    
    // Initialize vector pointers
    d0_pt = nullptr; d0_eta = nullptr; d0_phi = nullptr; d0_mass = nullptr;
    d0_px = nullptr; d0_py = nullptr; d0_pz = nullptr; d0_e = nullptr;
    dau_pid = nullptr; dau_pt = nullptr; dau_eta = nullptr; dau_phi = nullptr;
    dau_px = nullptr; dau_py = nullptr; dau_pz = nullptr; dau_d0_idx = nullptr;
    dau_pnn_k = nullptr; dau_pnn_pi = nullptr; dau_prb_ghost = nullptr; dau_trckChi2 = nullptr;
    mc_d0_pid = nullptr; mc_d0_pt = nullptr; mc_d0_eta = nullptr; mc_d0_phi = nullptr;
    mc_d0_px = nullptr; mc_d0_py = nullptr; mc_d0_pz = nullptr; mc_d0_matched = nullptr;
    mc_d0_matched_quality = nullptr;
    mc_dau_pid = nullptr; mc_dau_matched = nullptr;
    
    // Set default binning
    m_ptBins = {2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0};
    m_etaBins = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
    m_pBins = {5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0, 100.0};
}

D0RecoEfficiencyStandalone::~D0RecoEfficiencyStandalone() {
    CleanUp();
}

void D0RecoEfficiencyStandalone::CleanUp() {
    for (auto& pair : m_efficiencyMaps) {
        if (pair.second) delete pair.second;
    }
    m_efficiencyMaps.clear();
    
    for (auto& pair : m_efficiencyObjects) {
        if (pair.second) delete pair.second;
    }
    m_efficiencyObjects.clear();
    
    if (m_outputFile) {
        m_outputFile->Close();
        delete m_outputFile;
        m_outputFile = nullptr;
    }
    
    if (m_inputFile) {
        m_inputFile->Close();
        m_inputFile = nullptr;
    }
}

bool D0RecoEfficiencyStandalone::Initialize() {
    if (!m_inputFile || !m_outputFile || !m_tree) {
        std::cerr << "Error: Files or tree not properly initialized" << std::endl;
        return false;
    }
    
    InitializeBranches();
    CreateHistograms();
    
    std::cout << "D0RecoEfficiencyStandalone initialized successfully" << std::endl;
    std::cout << "Input file: " << m_inputFile->GetName() << std::endl;
    std::cout << "Output file: " << m_outputFile->GetName() << std::endl;
    std::cout << "Tree entries: " << m_tree->GetEntries() << std::endl;
    
    return true;
}

void D0RecoEfficiencyStandalone::InitializeBranches() {
    // Set branch addresses for reconstructed D0 variables
    m_tree->SetBranchAddress("n_d0s", &n_d0s);
    m_tree->SetBranchAddress("d0_pt", &d0_pt);
    m_tree->SetBranchAddress("d0_eta", &d0_eta);
    m_tree->SetBranchAddress("d0_phi", &d0_phi);
    m_tree->SetBranchAddress("d0_mass", &d0_mass);
    m_tree->SetBranchAddress("d0_px", &d0_px);
    m_tree->SetBranchAddress("d0_py", &d0_py);
    m_tree->SetBranchAddress("d0_pz", &d0_pz);
    m_tree->SetBranchAddress("d0_e", &d0_e);
    
    // Set branch addresses for daughter variables
    m_tree->SetBranchAddress("n_daughters", &n_daughters);
    m_tree->SetBranchAddress("dau_pid", &dau_pid);
    m_tree->SetBranchAddress("dau_pt", &dau_pt);
    m_tree->SetBranchAddress("dau_eta", &dau_eta);
    m_tree->SetBranchAddress("dau_phi", &dau_phi);
    m_tree->SetBranchAddress("dau_px", &dau_px);
    m_tree->SetBranchAddress("dau_py", &dau_py);
    m_tree->SetBranchAddress("dau_pz", &dau_pz);
    m_tree->SetBranchAddress("dau_d0_idx", &dau_d0_idx);
    m_tree->SetBranchAddress("dau_pnn_k", &dau_pnn_k);
    m_tree->SetBranchAddress("dau_pnn_pi", &dau_pnn_pi);
    m_tree->SetBranchAddress("dau_prb_ghost", &dau_prb_ghost);
    m_tree->SetBranchAddress("dau_trckChi2", &dau_trckChi2);
    
    // Set branch addresses for MC truth variables
    if (m_tree->FindBranch("mc_d0_pid")) {
        m_tree->SetBranchAddress("mc_d0_pid", &mc_d0_pid);
    }
    m_tree->SetBranchAddress("mc_d0_pt", &mc_d0_pt);
    m_tree->SetBranchAddress("mc_d0_eta", &mc_d0_eta);
    m_tree->SetBranchAddress("mc_d0_phi", &mc_d0_phi);
    m_tree->SetBranchAddress("mc_d0_px", &mc_d0_px);
    m_tree->SetBranchAddress("mc_d0_py", &mc_d0_py);
    m_tree->SetBranchAddress("mc_d0_pz", &mc_d0_pz);
    m_tree->SetBranchAddress("mc_d0_matched", &mc_d0_matched);
    m_tree->SetBranchAddress("mc_d0_matched_quality", &mc_d0_matched_quality);
    
    if (m_tree->FindBranch("mc_dau_pid")) {
        m_tree->SetBranchAddress("mc_dau_pid", &mc_dau_pid);
    }
    if (m_tree->FindBranch("mc_dau_matched")) {
        m_tree->SetBranchAddress("mc_dau_matched", &mc_dau_matched);
    }
}

void D0RecoEfficiencyStandalone::CreateHistograms() {
    // Create reconstruction efficiency histograms (pT vs eta)
    m_efficiencyMaps["reco_numerator"] = new TH2F("reco_numerator", "Reco Efficiency Numerator;p_{T} [GeV];#eta",
                                                  m_ptBins.size() - 1, &m_ptBins[0],
                                                  m_etaBins.size() - 1, &m_etaBins[0]);
    m_efficiencyMaps["reco_denominator"] = new TH2F("reco_denominator", "Reco Efficiency Denominator;p_{T} [GeV];#eta",
                                                    m_ptBins.size() - 1, &m_ptBins[0],
                                                    m_etaBins.size() - 1, &m_etaBins[0]);
    
    // Create TEfficiency object for proper error handling
    m_efficiencyObjects["reco_efficiency"] = new TEfficiency("reco_efficiency", "Reconstruction Efficiency;p_{T} [GeV];#eta",
                                                             m_ptBins.size() - 1, &m_ptBins[0],
                                                             m_etaBins.size() - 1, &m_etaBins[0]);
    
    // Create momentum-based histograms if p bins are defined
    if (!m_pBins.empty()) {
        m_efficiencyMaps["reco_numerator_p"] = new TH2F("reco_numerator_p", "Reco Efficiency Numerator (p);p [GeV];#eta",
                                                        m_pBins.size() - 1, &m_pBins[0],
                                                        m_etaBins.size() - 1, &m_etaBins[0]);
        m_efficiencyMaps["reco_denominator_p"] = new TH2F("reco_denominator_p", "Reco Efficiency Denominator (p);p [GeV];#eta",
                                                          m_pBins.size() - 1, &m_pBins[0],
                                                          m_etaBins.size() - 1, &m_etaBins[0]);
        
        m_efficiencyObjects["reco_efficiency_p"] = new TEfficiency("reco_efficiency_p", "Reconstruction Efficiency (p);p [GeV];#eta",
                                                                   m_pBins.size() - 1, &m_pBins[0],
                                                                   m_etaBins.size() - 1, &m_etaBins[0]);
    }
    
    std::cout << "Created efficiency histograms with:" << std::endl;
    std::cout << "  pT bins: " << m_ptBins.size() - 1 << " (" << m_ptBins.front() << " - " << m_ptBins.back() << " GeV)" << std::endl;
    std::cout << "  eta bins: " << m_etaBins.size() - 1 << " (" << m_etaBins.front() << " - " << m_etaBins.back() << ")" << std::endl;
    if (!m_pBins.empty()) {
        std::cout << "  p bins: " << m_pBins.size() - 1 << " (" << m_pBins.front() << " - " << m_pBins.back() << " GeV)" << std::endl;
    }
}

int D0RecoEfficiencyStandalone::GetMCMatchIndex(int reco_d0_idx) {
    // Search through MC D0 particles to find which one matches this reconstructed D0
    for (size_t mc_idx = 0; mc_idx < mc_d0_matched->size(); ++mc_idx) {
        if (mc_d0_matched->at(mc_idx) == reco_d0_idx) {
            return mc_idx;
        }
    }
    return -1; // No match found
}

//function to check match quality
int D0RecoEfficiencyStandalone::GetMCMatchQuality(int mc_idx) {
    if (mc_idx < 0 || mc_idx >= (int)mc_d0_matched_quality->size())
        return -1;
    
    float matchquality = mc_d0_matched_quality->at(mc_idx);
    if (matchquality < 0.9) {
        std::cerr << "Warning: Low match quality for MC D0 index " << mc_idx << ": " << matchquality << std::endl;
        return -1; // Low quality match
    }

    return matchquality; // Return match quality
}

bool D0RecoEfficiencyStandalone::HasValidMCMatch(int reco_d0_idx, int& mc_match_idx) {
    mc_match_idx = GetMCMatchIndex(reco_d0_idx);
    
    if(GetMCMatchQuality(mc_match_idx) < 0) {
        return false; // Invalid match quality
    }

    if (mc_match_idx < 0)
        return false;
    
    // Check if the matched MC D0 is in acceptance
    double mc_pt = mc_d0_pt->at(mc_match_idx);
    if (mc_pt < m_minPt)
        return false;
    
    double mc_eta = mc_d0_eta->at(mc_match_idx);
    if (mc_eta > m_maxEta || mc_eta < m_minEta)
        return false;
    
    return true;
}

bool D0RecoEfficiencyStandalone::PassesD0Selection(int d0_idx) {
    if (d0_idx < 0 || d0_idx >= n_d0s)
        return false;
    
    double pt = d0_pt->at(d0_idx);
    double eta = d0_eta->at(d0_idx);
    double mass = d0_mass->at(d0_idx);
    
    // D0 PDG mass is 1864.84 MeV
    const double d0_pdg_mass = 1.86484;
    
    // Apply selection cuts
    if (pt < m_minPt) return false;
    if (eta > m_maxEta || eta < m_minEta) return false;
    if (TMath::Abs(mass - d0_pdg_mass) > m_d0MassWindow * 0.001) return false; // Convert MeV to GeV
    
    // Apply quality cuts
    if (!PassesQualityCuts(d0_idx)) return false;
    
    return true;
}

bool D0RecoEfficiencyStandalone::PassesQualityCuts(int d0_idx) {
    // Find daughters of this D0 and check all cuts in one loop
    bool foundTrueKaon = false, foundTruePion = false;
    
    for (size_t dau_idx = 0; dau_idx < dau_d0_idx->size(); ++dau_idx) {
        if (dau_d0_idx->at(dau_idx) != d0_idx)
            continue;
        
        // Check daughter eta range
        double daughter_eta = dau_eta->at(dau_idx);
        if (daughter_eta < m_minEta || daughter_eta > m_maxEta)
            return false;
        
        // Check daughter momentum
        double daughter_px = dau_px->at(dau_idx);
        double daughter_py = dau_py->at(dau_idx);
        double daughter_pz = dau_pz->at(dau_idx);
        double daughter_p = TMath::Sqrt(daughter_px * daughter_px + daughter_py * daughter_py + daughter_pz * daughter_pz);
        if (daughter_p < m_minDaughterMomentum)
            return false;
        
        // Use MC truth PID instead of reconstructed PID
        if (mc_dau_pid && dau_idx < mc_dau_pid->size()) {
            int mc_pid = mc_dau_pid->at(dau_idx);
            if (abs(mc_pid) == 321) { // True kaon
                foundTrueKaon = true;
            }
            else if (abs(mc_pid) == 211) { // True pion
                foundTruePion = true;
            }
        }
    }
    
    // Check that we found both true kaon and true pion
    return foundTrueKaon && foundTruePion;
}

double D0RecoEfficiencyStandalone::CalculateMCMomentum(int mc_idx) {
    if (mc_idx < 0 || mc_idx >= (int)mc_d0_px->size())
        return -1.0;
    
    double px = mc_d0_px->at(mc_idx);
    double py = mc_d0_py->at(mc_idx);
    double pz = mc_d0_pz->at(mc_idx);
    
    return TMath::Sqrt(px * px + py * py + pz * pz);
}

void D0RecoEfficiencyStandalone::CalculateRecoEfficiency() {
    std::cout << "\nCalculating reconstruction efficiency (reco-based approach)..." << std::endl;
    
    // Get efficiency histograms
    TH2F *h_num = m_efficiencyMaps["reco_numerator"];
    TH2F *h_den = m_efficiencyMaps["reco_denominator"];
    TEfficiency *eff_obj = m_efficiencyObjects["reco_efficiency"];
    
    // Get momentum-based histograms (optional)
    TH2F *h_num_p = nullptr;
    TH2F *h_den_p = nullptr;
    TEfficiency *eff_obj_p = nullptr;
    
    if (!m_pBins.empty()) {
        h_num_p = m_efficiencyMaps["reco_numerator_p"];
        h_den_p = m_efficiencyMaps["reco_denominator_p"];
        eff_obj_p = m_efficiencyObjects["reco_efficiency_p"];
    }
    
    if (!h_num || !h_den || !eff_obj) {
        std::cerr << "Error: Efficiency histograms/objects not available" << std::endl;
        return;
    }
    
    // Loop over all events
    Long64_t nEntries = m_tree->GetEntries();
    std::cout << "Processing " << nEntries << " events..." << std::endl;
    
    // Counters for statistics
    int totalMCInAcceptance = 0;
    int totalRecoMatched = 0;
    int totalRecoSelected = 0;
    
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        m_tree->GetEntry(entry);
        
        if (entry % 50000 == 0) {
            std::cout << "Processing entry " << entry << "/" << nEntries 
                     << " (" << (100.0 * entry / nEntries) << "%)" << std::endl;
        }
        
        // === DENOMINATOR: Loop over MC D0 particles ===
        const size_t nMCD0s = mc_d0_pt->size();
        
        for (size_t mc_idx = 0; mc_idx < nMCD0s; ++mc_idx) {
            // Check if MC D0 is in acceptance
            double mc_pt = mc_d0_pt->at(mc_idx);
            if (mc_pt < m_minPt) continue;
            double mc_eta = mc_d0_eta->at(mc_idx);
            if (mc_eta > m_maxEta || mc_eta < m_minEta) continue;

            // Require both MC daughters in acceptance
            bool kaonInAcc = false, pionInAcc = false;
            if (mc_dau_pid && mc_dau_matched) {
                for (size_t dau_idx = 0; dau_idx < mc_dau_pid->size(); ++dau_idx) {
                    // Check if this daughter belongs to this MC D0
                    // If you have a mc_dau_d0_idx branch, use it for association
                    int pid = mc_dau_pid->at(dau_idx);
                    // For eta, if you have mc_dau_eta branch, use it; else calculate from px,py,pz
                    double dau_eta = 0.0;
                    if (dau_eta < m_minEta || dau_eta > m_maxEta) continue;
                    if (abs(pid) == 321) kaonInAcc = true;
                    if (abs(pid) == 211) pionInAcc = true;
                }
            } else {
                // If no MC daughter info, skip this check
                kaonInAcc = pionInAcc = true;
            }
            if (!(kaonInAcc && pionInAcc)) continue;

            // Calculate momentum only once if needed
            double mc_p = (h_den_p) ? CalculateMCMomentum(mc_idx) : 0.0;

            // Fill denominator (MC D0s in acceptance with both daughters in acceptance)
            h_den->Fill(mc_pt, mc_eta);
            if (h_den_p)
                h_den_p->Fill(mc_p, mc_eta);

            totalMCInAcceptance++;
        }
        
        // === NUMERATOR: Loop over reconstructed D0 particles ===
        for (int reco_idx = 0; reco_idx < n_d0s; ++reco_idx) {
            // First check if the reconstructed D0 passes selection cuts
            if (!PassesD0Selection(reco_idx))
                continue;
            
            // Check if this reconstructed D0 has a valid MC match in acceptance
            int mc_match_idx = -1;
            if (!HasValidMCMatch(reco_idx, mc_match_idx))
                continue;
            
            totalRecoMatched++;
            
            // Get MC kinematics for filling histograms
            double mc_pt = mc_d0_pt->at(mc_match_idx);
            double mc_eta = mc_d0_eta->at(mc_match_idx);
            double mc_p = (h_num_p) ? CalculateMCMomentum(mc_match_idx) : 0.0;
            
            // Fill numerator (reconstructed D0s that pass selection and have MC match in acceptance)
            h_num->Fill(mc_pt, mc_eta);
            if (h_num_p)
                h_num_p->Fill(mc_p, mc_eta);
            
            totalRecoSelected++;
        }
        
        // === Fill TEfficiency objects ===
        // We need to loop over MC particles again for TEfficiency
        for (size_t mc_idx = 0; mc_idx < nMCD0s; ++mc_idx) {
            // Check if MC D0 is in acceptance
            double mc_pt = mc_d0_pt->at(mc_idx);
            if (mc_pt < m_minPt) continue;
            
            double mc_eta = mc_d0_eta->at(mc_idx);
            if (mc_eta > m_maxEta || mc_eta < m_minEta) continue;
            
            double mc_p = (eff_obj_p) ? CalculateMCMomentum(mc_idx) : 0.0;
            
            // Check if this MC D0 has a reconstructed match that passes selection
            int matchedRecoIdx = mc_d0_matched->at(mc_idx);
            bool passesReco = false;
            
            if (matchedRecoIdx >= 0 && matchedRecoIdx < n_d0s) {
                if (PassesD0Selection(matchedRecoIdx)) {
                    passesReco = true;
                }
            }
            
            // Fill TEfficiency object
            eff_obj->Fill(passesReco, mc_pt, mc_eta);
            if (eff_obj_p)
                eff_obj_p->Fill(passesReco, mc_p, mc_eta);
        }
    }
    
    // Create efficiency map (pT vs eta)
    TH2F *h_eff = (TH2F *)h_num->Clone("reco_efficiency");
    h_eff->Divide(h_den);
    h_eff->SetTitle("Reconstruction Efficiency (Reco-based);p_{T} [GeV];#eta");
    m_efficiencyMaps["reco_efficiency"] = h_eff;
    
    // Create efficiency map (p vs eta) if available
    if (h_num_p && h_den_p) {
        TH2F *h_eff_p = (TH2F *)h_num_p->Clone("reco_efficiency_p");
        h_eff_p->Divide(h_den_p);
        h_eff_p->SetTitle("Reconstruction Efficiency (Reco-based);p [GeV];#eta");
        m_efficiencyMaps["reco_efficiency_p"] = h_eff_p;
    }
    
    std::cout << "Reconstruction efficiency calculation completed (reco-based approach)" << std::endl;
    
    // Print statistics
    int totalDen = h_den->GetEntries();
    int totalNum = h_num->GetEntries();
    double overallEff = (totalDen > 0) ? (double)totalNum / totalDen : 0.0;
    
    std::cout << "\nStatistics:" << std::endl;
    std::cout << "  Total MC D0s in acceptance (denominator): " << totalDen << std::endl;
    std::cout << "  Total reconstructed D0s with MC match and passing selection (numerator): " << totalNum << std::endl;
    std::cout << "  Overall reconstruction efficiency: " << overallEff * 100.0 << "% (" << totalNum << "/" << totalDen << ")" << std::endl;
    std::cout << "\nDetailed breakdown:" << std::endl;
    std::cout << "  MC D0s processed in acceptance: " << totalMCInAcceptance << std::endl;
    std::cout << "  Reconstructed D0s with valid MC match: " << totalRecoMatched << std::endl;
    std::cout << "  Reconstructed D0s passing all cuts: " << totalRecoSelected << std::endl;
}

void D0RecoEfficiencyStandalone::SaveResults() {
    if (!m_outputFile) {
        std::cerr << "Error: Output file not available" << std::endl;
        return;
    }
    
    std::cout << "\nSaving results to " << m_outputFile->GetName() << "..." << std::endl;
    
    m_outputFile->cd();
    
    // Save all efficiency histograms
    for (auto& pair : m_efficiencyMaps) {
        if (pair.second) {
            pair.second->Write();
            std::cout << "  Saved histogram: " << pair.first << std::endl;
        }
    }
    
    // Save all efficiency objects
    for (auto& pair : m_efficiencyObjects) {
        if (pair.second) {
            pair.second->Write();
            std::cout << "  Saved efficiency object: " << pair.first << std::endl;
        }
    }
    
    m_outputFile->Write();
    std::cout << "Results saved successfully" << std::endl;
}

void D0RecoEfficiencyStandalone::PlotEfficiency(const std::string &histName) {
    TH2F *h_eff = m_efficiencyMaps[histName];
    if (!h_eff) {
        std::cerr << "Error: Efficiency histogram '" << histName << "' not found" << std::endl;
        return;
    }
    
    std::cout << "Plotting efficiency: " << histName << std::endl;
    
    // Set up canvas
    TCanvas *c1 = new TCanvas("c_eff", "Reconstruction Efficiency", 800, 600);
    c1->SetRightMargin(0.15);
    
    // Set up histogram style
    h_eff->SetStats(0);
    h_eff->GetZaxis()->SetTitle("Efficiency");
    h_eff->GetZaxis()->SetRangeUser(0.0, 1.0);
    
    // Draw efficiency map
    h_eff->Draw("COLZ");
    
    // Add text with overall efficiency
    TH2F *h_num = m_efficiencyMaps["reco_numerator"];
    TH2F *h_den = m_efficiencyMaps["reco_denominator"];
    if (h_num && h_den) {
        double overallEff = (double)h_num->GetEntries() / h_den->GetEntries();
        TString effText = TString::Format("Overall Efficiency: %.3f", overallEff);
        
        TLatex *latex = new TLatex();
        latex->SetNDC();
        latex->SetTextSize(0.04);
        latex->DrawLatex(0.15, 0.92, effText);
    }
    
    // Save canvas
    TString plotName = TString::Format("%s.png", histName.c_str());
    c1->SaveAs(plotName);
    std::cout << "Saved plot: " << plotName << std::endl;
    
    // Also save as PDF
    plotName = TString::Format("%s.pdf", histName.c_str());
    c1->SaveAs(plotName);
    
    delete c1;
}

// Main function for the standalone reconstruction-based efficiency calculation
int D0RecoEfficiencyStandaloneRun(TString inputFile = "/media/niviths/SSD2/lhcb_analysis_SSD/20250728_pPb_MC_output/20250728_pPb_MC_output.root", 
                                 TString outputFile = "output_reco_standalone.root",
                                 double massWindow = 50.0, double minPt = 1.0,
                                 double minEta = 2.0, double maxEta = 4.5,
                                 double kaonPIDCut = 0.5, double pionPIDCut = 0.5,
                                 bool makePlots = true) {
    std::cout << "=== D0 Reconstruction Efficiency Calculator (Standalone Reco-based) ===" << std::endl;
    std::cout << "Input file: " << inputFile << std::endl;
    std::cout << "Output file: " << outputFile << std::endl;
    std::cout << std::endl;
    std::cout << "Configuration:" << std::endl;
    std::cout << "  D0 mass window: " << massWindow << " MeV" << std::endl;
    std::cout << "  Minimum pT: " << minPt << " GeV" << std::endl;
    std::cout << "  Eta range: " << minEta << " - " << maxEta << std::endl;
    std::cout << "  PID cuts (K/π): " << kaonPIDCut << " / " << pionPIDCut << std::endl;
    std::cout << "  Generate plots: " << (makePlots ? "Yes" : "No") << std::endl;
    std::cout << std::endl;
    
    std::cout << "Method: Standalone Reconstruction-based approach" << std::endl;
    std::cout << "  - Denominator: All MC D0s in acceptance" << std::endl;
    std::cout << "  - Numerator: Reconstructed D0s passing selection with MC match in acceptance" << std::endl;
    std::cout << std::endl;
    
    try {
        // Create efficiency calculator with standalone reco-based approach
        D0RecoEfficiencyStandalone calculator(inputFile, outputFile);
        
        // Configure parameters
        calculator.SetD0MassWindow(massWindow);
        calculator.SetPtRange(minPt);
        calculator.SetEtaRange(minEta, maxEta);
        calculator.SetPIDCuts(kaonPIDCut, pionPIDCut);
        
        // Set up custom binning
        std::vector<double> customPtBins;
        for (double pt = 2.0; pt < 6.0; pt += 0.25)
            customPtBins.push_back(pt);
        for (double pt = 6.0; pt < 10.0; pt += 0.5)
            customPtBins.push_back(pt);
        for (double pt = 10.0; pt < 20.0; pt += 1.0)
            customPtBins.push_back(pt);
        for (double pt = 20.0; pt < 30.0; pt += 2.5)
            customPtBins.push_back(pt);
        for (double pt = 30.0; pt < 60.0; pt += 5.0)
            customPtBins.push_back(pt);
        if (customPtBins.back() < 60.0)
            customPtBins.push_back(60.0);
        
        std::vector<double> customEtaBins;
        for (double eta = 2.0; eta <= 4.5; eta += 0.125)
            customEtaBins.push_back(eta);
        
        std::vector<double> customPBins;
        for (double p = 5.0; p < 20.0; p += 2.5)
            customPBins.push_back(p);
        for (double p = 20.0; p < 40.0; p += 2.5)
            customPBins.push_back(p);
        for (double p = 40.0; p < 100.0; p += 10.0)
            customPBins.push_back(p);
        if (customPBins.back() < 100.0)
            customPBins.push_back(100.0);
        
        calculator.SetPtBins(customPtBins);
        calculator.SetEtaBins(customEtaBins);
        calculator.SetPBins(customPBins);
        
        // Initialize
        if (!calculator.Initialize()) {
            std::cerr << "Error: Failed to initialize calculator" << std::endl;
            return 1;
        }
        
        // Calculate efficiency
        std::cout << "Starting efficiency calculation..." << std::endl;
        calculator.CalculateRecoEfficiency();
        
        // Save results
        calculator.SaveResults();
        
        // Generate plots if requested
        if (makePlots) {
            std::cout << "\nGenerating efficiency plots..." << std::endl;
            calculator.PlotEfficiency("reco_efficiency");
            calculator.PlotEfficiency("reco_efficiency_p");
        }
        
        std::cout << "\n=== Calculation completed successfully! ===" << std::endl;
    }
    catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
