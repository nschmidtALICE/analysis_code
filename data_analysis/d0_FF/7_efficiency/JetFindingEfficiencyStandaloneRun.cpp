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

// Standalone Jet Finding Efficiency Calculator
class JetFindingEfficiencyStandalone {
public:
    JetFindingEfficiencyStandalone(TString inputFileName, TString outputFileName);
    ~JetFindingEfficiencyStandalone();
    
    bool Initialize();
    void CalculateJetFindingEfficiency();
    void SaveResults();
    void PlotEfficiency(const std::string &histName);
    
    // Configuration methods
    void SetJetRadius(double radius) { m_jetRadius = radius; }
    void SetJetPtRange(double minPt, double maxPt = 1000.0) { m_minJetPt = minPt; m_maxJetPt = maxPt; }
    void SetJetEtaRange(double minEta, double maxEta) { m_minJetEta = minEta; m_maxJetEta = maxEta; }
    void SetD0PtRange(double minPt, double maxPt = 1000.0) { m_minD0Pt = minPt; m_maxD0Pt = maxPt; }
    void SetD0EtaRange(double minEta, double maxEta) { m_minD0Eta = minEta; m_maxD0Eta = maxEta; }
    void SetJetPtBins(const std::vector<double>& bins) { m_jetPtBins = bins; }
    void SetJetEtaBins(const std::vector<double>& bins) { m_jetEtaBins = bins; }
    void SetD0PtBins(const std::vector<double>& bins) { m_d0PtBins = bins; }
    
private:
    // File and tree management
    TFile *m_inputFile;
    TFile *m_outputFile;
    TTree *m_tree;
    
    // Configuration parameters
    double m_jetRadius;
    double m_minJetPt;
    double m_maxJetPt;
    double m_minJetEta;
    double m_maxJetEta;
    double m_minD0Pt;
    double m_maxD0Pt;
    double m_minD0Eta;
    double m_maxD0Eta;
    
    // Binning
    std::vector<double> m_jetPtBins;
    std::vector<double> m_jetEtaBins;
    std::vector<double> m_d0PtBins;
    
    // Histograms and efficiency objects
    std::map<std::string, TH2F*> m_efficiencyMaps;
    std::map<std::string, TEfficiency*> m_efficiencyObjects;
    
    // Tree branches - Reconstructed jets
    std::vector<float> *jet_pt;
    std::vector<float> *jet_eta;
    std::vector<float> *jet_phi;
    std::vector<float> *jet_mass;
    std::vector<int> *jet_n_const;
    std::vector<int> *jet_n_charged;
    std::vector<int> *jet_n_neutral;
    std::vector<int> *jet_n_d0;
    
    // Tree branches - Reconstructed D0s
    std::vector<float> *d0_pt;
    std::vector<float> *d0_eta;
    std::vector<float> *d0_phi;
    std::vector<float> *d0_mass;
    std::vector<float> *d0_pz;
    std::vector<float> *d0_e;
    std::vector<int> *d0_jet_idx;
    
    // Tree branches - MC truth jets (derived from px,py,pz,e)
    std::vector<float> *mc_jet_px;
    std::vector<float> *mc_jet_py;
    std::vector<float> *mc_jet_pz;
    std::vector<float> *mc_jet_e;
    std::vector<int> *mc_jet_n_const;
    std::vector<int> *mc_jet_n_chr;
    std::vector<int> *mc_jet_n_neu;
    
    // Tree branches - MC truth D0s
    std::vector<int> *mc_d0_pid;
    std::vector<float> *mc_d0_pt;
    std::vector<float> *mc_d0_eta;
    std::vector<float> *mc_d0_phi;
    std::vector<float> *mc_d0_mass;
    std::vector<int> *mc_d0_origin;
    std::vector<int> *mc_d0_matched;
    std::vector<int> *mc_d0_jet_idx;
    std::vector<float> *mc_d0_z;
    
    // Helper methods
    void InitializeBranches();
    void CreateHistograms();
    void CleanUp();
    double CalculateDeltaR(double eta1, double phi1, double eta2, double phi2);
    bool HasD0InJetRadius(double jet_eta, double jet_phi, const std::vector<float>* d0_eta_vec, 
                         const std::vector<float>* d0_phi_vec, const std::vector<float>* d0_pt_vec);
    bool PassesJetSelection(int jet_idx, bool isMC = false);
    bool PassesD0Selection(double d0_pt, double d0_eta);
    
    // MC jet kinematics helpers
    double GetMCJetPt(int mc_jet_idx);
    double GetMCJetEta(int mc_jet_idx);
    double GetMCJetPhi(int mc_jet_idx);
    
    // Jet-D0 association using d0_jet_idx
    bool IsD0AssociatedWithJet(int d0_idx, int jet_idx);
    bool IsD0AssociatedWithMCJet(int mc_d0_idx, int mc_jet_idx);
    
    // MC matching helpers (using d0_jet_idx and mc_d0_jet_idx)
    int GetMCJetForD0(int mc_d0_idx);
    int GetRecoJetForD0(int d0_idx);
    int GetMCJetMatchIndex(int reco_d0_idx);
    bool HasMatchingJetPair(int mc_jet_idx, int& reco_jet_idx);
};

JetFindingEfficiencyStandalone::JetFindingEfficiencyStandalone(TString inputFileName, TString outputFileName)
    : m_inputFile(nullptr), m_outputFile(nullptr), m_tree(nullptr),
      m_jetRadius(0.4), m_minJetPt(5.0), m_maxJetPt(1000.0),
      m_minJetEta(2.0), m_maxJetEta(4.5),
      m_minD0Pt(2.0), m_maxD0Pt(1000.0),
      m_minD0Eta(2.0), m_maxD0Eta(4.5)
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
    jet_pt = nullptr; jet_eta = nullptr; jet_phi = nullptr; jet_mass = nullptr;
    jet_n_const = nullptr; jet_n_charged = nullptr; jet_n_neutral = nullptr; jet_n_d0 = nullptr;
    d0_pt = nullptr; d0_eta = nullptr; d0_phi = nullptr; d0_mass = nullptr;
    d0_pz = nullptr; d0_e = nullptr; d0_jet_idx = nullptr;
    mc_jet_px = nullptr; mc_jet_py = nullptr; mc_jet_pz = nullptr; mc_jet_e = nullptr;
    mc_jet_n_const = nullptr; mc_jet_n_chr = nullptr; mc_jet_n_neu = nullptr;
    mc_d0_pid = nullptr; mc_d0_pt = nullptr; mc_d0_eta = nullptr; mc_d0_phi = nullptr;
    mc_d0_mass = nullptr; mc_d0_origin = nullptr; mc_d0_matched = nullptr;
    mc_d0_jet_idx = nullptr; mc_d0_z = nullptr;
    
    // Set jet pT bins from 0 to 60 with 2.5 GeV width
    m_jetPtBins.clear();
    for (double pt = 0.0; pt <= 60.0; pt += 2.5)
        m_jetPtBins.push_back(pt);
    m_jetEtaBins.clear();
    for (double eta = 2.0; eta <= 4.5; eta += 0.2)
        m_jetEtaBins.push_back(eta);
    m_d0PtBins = m_jetPtBins;
}

JetFindingEfficiencyStandalone::~JetFindingEfficiencyStandalone() {
    CleanUp();
}

void JetFindingEfficiencyStandalone::CleanUp() {
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

bool JetFindingEfficiencyStandalone::Initialize() {
    if (!m_inputFile || !m_outputFile || !m_tree) {
        std::cerr << "Error: Files or tree not properly initialized" << std::endl;
        return false;
    }
    
    InitializeBranches();
    CreateHistograms();
    
    std::cout << "JetFindingEfficiencyStandalone initialized successfully" << std::endl;
    std::cout << "Input file: " << m_inputFile->GetName() << std::endl;
    std::cout << "Output file: " << m_outputFile->GetName() << std::endl;
    std::cout << "Tree entries: " << m_tree->GetEntries() << std::endl;
    std::cout << "Jet radius: " << m_jetRadius << std::endl;
    
    return true;
}

void JetFindingEfficiencyStandalone::InitializeBranches() {
    // Set branch addresses for reconstructed jets
    m_tree->SetBranchAddress("jet_pt", &jet_pt);
    m_tree->SetBranchAddress("jet_eta", &jet_eta);
    m_tree->SetBranchAddress("jet_phi", &jet_phi);
    m_tree->SetBranchAddress("jet_mass", &jet_mass);
    m_tree->SetBranchAddress("jet_n_const", &jet_n_const);
    m_tree->SetBranchAddress("jet_n_charged", &jet_n_charged);
    m_tree->SetBranchAddress("jet_n_neutral", &jet_n_neutral);
    m_tree->SetBranchAddress("jet_n_d0", &jet_n_d0);
    
    // Set branch addresses for reconstructed D0s
    m_tree->SetBranchAddress("d0_pt", &d0_pt);
    m_tree->SetBranchAddress("d0_eta", &d0_eta);
    m_tree->SetBranchAddress("d0_phi", &d0_phi);
    m_tree->SetBranchAddress("d0_mass", &d0_mass);
    m_tree->SetBranchAddress("d0_pz", &d0_pz);
    m_tree->SetBranchAddress("d0_e", &d0_e);
    m_tree->SetBranchAddress("d0_jet_idx", &d0_jet_idx);
    
    // Set branch addresses for MC truth jets
    m_tree->SetBranchAddress("mc_jet_px", &mc_jet_px);
    m_tree->SetBranchAddress("mc_jet_py", &mc_jet_py);
    m_tree->SetBranchAddress("mc_jet_pz", &mc_jet_pz);
    m_tree->SetBranchAddress("mc_jet_e", &mc_jet_e);
    m_tree->SetBranchAddress("mc_jet_n_const", &mc_jet_n_const);
    m_tree->SetBranchAddress("mc_jet_n_charged", &mc_jet_n_chr);
    m_tree->SetBranchAddress("mc_jet_n_neutral", &mc_jet_n_neu);
    
    // Set branch addresses for MC truth D0s
    m_tree->SetBranchAddress("mc_d0_pid", &mc_d0_pid);
    m_tree->SetBranchAddress("mc_d0_pt", &mc_d0_pt);
    m_tree->SetBranchAddress("mc_d0_eta", &mc_d0_eta);
    m_tree->SetBranchAddress("mc_d0_phi", &mc_d0_phi);
    m_tree->SetBranchAddress("mc_d0_mass", &mc_d0_mass);
    m_tree->SetBranchAddress("mc_d0_origin", &mc_d0_origin);
    m_tree->SetBranchAddress("mc_d0_matched", &mc_d0_matched);
    m_tree->SetBranchAddress("mc_d0_jet_idx", &mc_d0_jet_idx);
    m_tree->SetBranchAddress("mc_d0_z", &mc_d0_z);
}

void JetFindingEfficiencyStandalone::CreateHistograms() {
    // Create jet finding efficiency histograms (jet pT vs jet eta)
    m_efficiencyMaps["jet_numerator"] = new TH2F("jet_numerator", "Jet Finding Efficiency Numerator;Jet p_{T} [GeV];Jet #eta",
                                                 m_jetPtBins.size() - 1, &m_jetPtBins[0],
                                                 m_jetEtaBins.size() - 1, &m_jetEtaBins[0]);
    m_efficiencyMaps["jet_denominator"] = new TH2F("jet_denominator", "Jet Finding Efficiency Denominator;Jet p_{T} [GeV];Jet #eta",
                                                   m_jetPtBins.size() - 1, &m_jetPtBins[0],
                                                   m_jetEtaBins.size() - 1, &m_jetEtaBins[0]);
    
    // Create TEfficiency object for proper error handling
    m_efficiencyObjects["jet_finding_efficiency"] = new TEfficiency("jet_finding_efficiency", 
                                                                    "Jet Finding Efficiency;Jet p_{T} [GeV];Jet #eta",
                                                                    m_jetPtBins.size() - 1, &m_jetPtBins[0],
                                                                    m_jetEtaBins.size() - 1, &m_jetEtaBins[0]);
    
    // Create D0 pT-based histograms
    m_efficiencyMaps["jet_numerator_d0pt"] = new TH2F("jet_numerator_d0pt", "Jet Finding Efficiency vs D0 pT;D0 p_{T} [GeV];Jet #eta",
                                                      m_d0PtBins.size() - 1, &m_d0PtBins[0],
                                                      m_jetEtaBins.size() - 1, &m_jetEtaBins[0]);
    m_efficiencyMaps["jet_denominator_d0pt"] = new TH2F("jet_denominator_d0pt", "Jet Finding Efficiency vs D0 pT;D0 p_{T} [GeV];Jet #eta",
                                                        m_d0PtBins.size() - 1, &m_d0PtBins[0],
                                                        m_jetEtaBins.size() - 1, &m_jetEtaBins[0]);
    
    m_efficiencyObjects["jet_finding_efficiency_d0pt"] = new TEfficiency("jet_finding_efficiency_d0pt", 
                                                                         "Jet Finding Efficiency vs D0 pT;D0 p_{T} [GeV];Jet #eta",
                                                                         m_d0PtBins.size() - 1, &m_d0PtBins[0],
                                                                         m_jetEtaBins.size() - 1, &m_jetEtaBins[0]);

    // Add zT vs jet pT efficiency histograms
    int nJetPtBins = m_jetPtBins.size() - 1;
    int nZTbins = 50; // You can adjust this binning as needed
    double zTmin = 0.0, zTmax = 1.0;
    m_efficiencyMaps["jet_numerator_zT"] = new TH2F("jet_numerator_zT", "Jet Finding Efficiency Numerator;z_{T};Jet p_{T} [GeV]",
                                                    nZTbins, zTmin, zTmax, nJetPtBins, &m_jetPtBins[0]);
    m_efficiencyMaps["jet_denominator_zT"] = new TH2F("jet_denominator_zT", "Jet Finding Efficiency Denominator;z_{T};Jet p_{T} [GeV]",
                                                      nZTbins, zTmin, zTmax, nJetPtBins, &m_jetPtBins[0]);
    
    std::cout << "Created efficiency histograms with:" << std::endl;
    std::cout << "  Jet pT bins: " << m_jetPtBins.size() - 1 << " (" << m_jetPtBins.front() << " - " << m_jetPtBins.back() << " GeV)" << std::endl;
    std::cout << "  Jet eta bins: " << m_jetEtaBins.size() - 1 << " (" << m_jetEtaBins.front() << " - " << m_jetEtaBins.back() << ")" << std::endl;
    std::cout << "  D0 pT bins: " << m_d0PtBins.size() - 1 << " (" << m_d0PtBins.front() << " - " << m_d0PtBins.back() << " GeV)" << std::endl;
}

double JetFindingEfficiencyStandalone::CalculateDeltaR(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = TMath::Abs(phi1 - phi2);
    if (dphi > TMath::Pi()) dphi = 2 * TMath::Pi() - dphi;
    return TMath::Sqrt(deta * deta + dphi * dphi);
}

bool JetFindingEfficiencyStandalone::HasD0InJetRadius(double jet_eta, double jet_phi, 
                                                     const std::vector<float>* d0_eta_vec, 
                                                     const std::vector<float>* d0_phi_vec,
                                                     const std::vector<float>* d0_pt_vec) {
    if (!d0_eta_vec || !d0_phi_vec || !d0_pt_vec) return false;
    
    for (size_t i = 0; i < d0_eta_vec->size(); ++i) {
        double d0_pt = d0_pt_vec->at(i);
        double d0_eta = d0_eta_vec->at(i);
        double d0_phi = d0_phi_vec->at(i);
        
        // Check if D0 passes selection cuts
        if (!PassesD0Selection(d0_pt, d0_eta)) continue;
        
        // Check if D0 is within jet radius
        double deltaR = CalculateDeltaR(jet_eta, jet_phi, d0_eta, d0_phi);
        if (deltaR < m_jetRadius) {
            return true;
        }
    }
    
    return false;
}

double JetFindingEfficiencyStandalone::GetMCJetPt(int mc_jet_idx) {
    if (mc_jet_idx < 0 || !mc_jet_px || !mc_jet_py || 
        mc_jet_idx >= (int)mc_jet_px->size() || mc_jet_idx >= (int)mc_jet_py->size()) {
        return -1.0;
    }
    
    double px = mc_jet_px->at(mc_jet_idx);
    double py = mc_jet_py->at(mc_jet_idx);
    return TMath::Sqrt(px * px + py * py);
}

double JetFindingEfficiencyStandalone::GetMCJetEta(int mc_jet_idx) {
    if (mc_jet_idx < 0 || !mc_jet_px || !mc_jet_py || !mc_jet_pz || 
        mc_jet_idx >= (int)mc_jet_px->size() || mc_jet_idx >= (int)mc_jet_py->size() || 
        mc_jet_idx >= (int)mc_jet_pz->size()) {
        return -999.0;
    }
    
    double px = mc_jet_px->at(mc_jet_idx);
    double py = mc_jet_py->at(mc_jet_idx);
    double pz = mc_jet_pz->at(mc_jet_idx);
    double pt = TMath::Sqrt(px * px + py * py);
    
    if (pt <= 0) return -999.0;
    
    return TMath::ASinH(pz / pt);
}

double JetFindingEfficiencyStandalone::GetMCJetPhi(int mc_jet_idx) {
    if (mc_jet_idx < 0 || !mc_jet_px || !mc_jet_py || 
        mc_jet_idx >= (int)mc_jet_px->size() || mc_jet_idx >= (int)mc_jet_py->size()) {
        return -999.0;
    }
    
    double px = mc_jet_px->at(mc_jet_idx);
    double py = mc_jet_py->at(mc_jet_idx);
    return TMath::ATan2(py, px);
}

bool JetFindingEfficiencyStandalone::IsD0AssociatedWithJet(int d0_idx, int jet_idx) {
    if (!d0_jet_idx || d0_idx < 0 || d0_idx >= (int)d0_jet_idx->size()) {
        return false;
    }
    return d0_jet_idx->at(d0_idx) == jet_idx;
}

bool JetFindingEfficiencyStandalone::IsD0AssociatedWithMCJet(int mc_d0_idx, int mc_jet_idx) {
    if (!mc_d0_jet_idx || mc_d0_idx < 0 || mc_d0_idx >= (int)mc_d0_jet_idx->size()) {
        return false;
    }
    return mc_d0_jet_idx->at(mc_d0_idx) == mc_jet_idx;
}

int JetFindingEfficiencyStandalone::GetMCJetForD0(int mc_d0_idx) {
    if (!mc_d0_jet_idx || mc_d0_idx < 0 || mc_d0_idx >= (int)mc_d0_jet_idx->size()) {
        return -1;
    }
    return mc_d0_jet_idx->at(mc_d0_idx);
}

int JetFindingEfficiencyStandalone::GetRecoJetForD0(int d0_idx) {
    if (!d0_jet_idx || d0_idx < 0 || d0_idx >= (int)d0_jet_idx->size()) {
        return -1;
    }
    return d0_jet_idx->at(d0_idx);
}

bool JetFindingEfficiencyStandalone::HasMatchingJetPair(int mc_jet_idx, int& reco_jet_idx) {
    // Find a MC D0 associated with this MC jet that is also reconstructed
    if (!mc_d0_jet_idx || !mc_d0_matched) return false;
    
    for (size_t mc_d0_idx = 0; mc_d0_idx < mc_d0_jet_idx->size(); ++mc_d0_idx) {
        // Check if this MC D0 is associated with the MC jet
        if (mc_d0_jet_idx->at(mc_d0_idx) != mc_jet_idx) continue;
        
        // Check if this MC D0 has a reconstructed match
        if (mc_d0_idx >= mc_d0_matched->size()) continue;
        int matched_reco_d0_idx = mc_d0_matched->at(mc_d0_idx);
        if (matched_reco_d0_idx < 0) continue;
        
        // Get the reconstructed jet associated with the matched reco D0
        int matched_reco_jet_idx = GetRecoJetForD0(matched_reco_d0_idx);
        if (matched_reco_jet_idx >= 0) {
            reco_jet_idx = matched_reco_jet_idx;
            return true;
        }
    }
    
    reco_jet_idx = -1;
    return false;
}

bool JetFindingEfficiencyStandalone::PassesJetSelection(int jet_idx, bool isMC) {
    if (isMC) {
        // For MC jets, calculate kinematics from px,py,pz,e
        double pt = GetMCJetPt(jet_idx);
        double eta = GetMCJetEta(jet_idx);
        
        if (pt < 0 || eta < -900) return false; // Invalid kinematics
        if (pt < m_minJetPt || pt > m_maxJetPt) return false;
        if (eta < m_minJetEta || eta > m_maxJetEta) return false;
    } else {
        // For reconstructed jets, use direct branches
        if (!jet_pt || !jet_eta || jet_idx < 0 || 
            jet_idx >= (int)jet_pt->size() || jet_idx >= (int)jet_eta->size()) {
            return false;
        }
        
        double pt = jet_pt->at(jet_idx);
        double eta = jet_eta->at(jet_idx);
        
        if (pt < m_minJetPt || pt > m_maxJetPt) return false;
        if (eta < m_minJetEta || eta > m_maxJetEta) return false;
    }
    
    return true;
}

bool JetFindingEfficiencyStandalone::PassesD0Selection(double d0_pt, double d0_eta) {
    if (d0_pt < m_minD0Pt || d0_pt > m_maxD0Pt) return false;
    if (d0_eta < m_minD0Eta || d0_eta > m_maxD0Eta) return false;
    return true;
}

int JetFindingEfficiencyStandalone::GetMCJetMatchIndex(int reco_jet_idx) {
    // This function is not available since mc_jet_matched doesn't exist in the tree
    // Jet matching should be done through D0 association instead
    return -1; // Always return no match since we don't have direct jet matching info
}

void JetFindingEfficiencyStandalone::CalculateJetFindingEfficiency() {
    std::cout << "\nCalculating jet finding efficiency..." << std::endl;
    
    // Get efficiency histograms
    TH2F *h_num = m_efficiencyMaps["jet_numerator"];
    TH2F *h_den = m_efficiencyMaps["jet_denominator"];
    TEfficiency *eff_obj = m_efficiencyObjects["jet_finding_efficiency"];
    
    TH2F *h_num_d0pt = m_efficiencyMaps["jet_numerator_d0pt"];
    TH2F *h_den_d0pt = m_efficiencyMaps["jet_denominator_d0pt"];
    TEfficiency *eff_obj_d0pt = m_efficiencyObjects["jet_finding_efficiency_d0pt"];
    
    if (!h_num || !h_den || !eff_obj || !h_num_d0pt || !h_den_d0pt || !eff_obj_d0pt) {
        std::cerr << "Error: Efficiency histograms/objects not available" << std::endl;
        return;
    }
    
    // Loop over all events
    Long64_t nEntries = m_tree->GetEntries();
    std::cout << "Processing " << nEntries << " events..." << std::endl;
    
    // Counters for statistics
    int totalMCJetsWithD0 = 0;
    int totalRecoJetsWithD0 = 0;
    int totalMatchedJets = 0;
    
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        m_tree->GetEntry(entry);
        
        if (entry % 50000 == 0) {
            std::cout << "Processing entry " << entry << "/" << nEntries 
                     << " (" << (100.0 * entry / nEntries) << "%)" << std::endl;
        }
        
        // === DENOMINATOR: Loop over MC jets with generator level D0 ===
        const size_t nMCJets = mc_jet_px ? mc_jet_px->size() : 0;
        
        for (size_t mc_jet_idx = 0; mc_jet_idx < nMCJets; ++mc_jet_idx) {
            // Check if MC jet passes selection
            if (!PassesJetSelection(mc_jet_idx, true)) continue;

            double mc_jet_pt = GetMCJetPt(mc_jet_idx);
            double mc_jet_eta = GetMCJetEta(mc_jet_idx);
            double mc_jet_phi = GetMCJetPhi(mc_jet_idx);

            // Skip if kinematics calculation failed
            if (mc_jet_pt < 0 || mc_jet_eta < -900) continue;

            // Check if this MC jet has a generator level D0 in the jet radius
            if (!HasD0InJetRadius(mc_jet_eta, mc_jet_phi, mc_d0_eta, mc_d0_phi, mc_d0_pt)) {
                continue;
            }

            // Find the highest pT D0 in the jet radius for D0-pT based histograms
            double maxD0Pt = -1.0;
            for (size_t d0_idx = 0; d0_idx < mc_d0_pt->size(); ++d0_idx) {
                double d0_pt = mc_d0_pt->at(d0_idx);
                double d0_eta = mc_d0_eta->at(d0_idx);
                double d0_phi = mc_d0_phi->at(d0_idx);

                if (!PassesD0Selection(d0_pt, d0_eta)) continue;

                double deltaR = CalculateDeltaR(mc_jet_eta, mc_jet_phi, d0_eta, d0_phi);
                if (deltaR < m_jetRadius && d0_pt > maxD0Pt) {
                    maxD0Pt = d0_pt;
                }
            }

            // Fill denominators (MC jets with generator level D0 in jet radius)
            h_den->Fill(mc_jet_pt, mc_jet_eta);
            if (maxD0Pt > 0) {
                h_den_d0pt->Fill(maxD0Pt, mc_jet_eta);
                // Fill zT denominator histogram
                double zT = maxD0Pt / mc_jet_pt;
                if (zT > 0 && zT < 2.0) {
                    m_efficiencyMaps["jet_denominator_zT"]->Fill(zT, mc_jet_pt);
                }
            }

            // Fill TEfficiency for this MC jet
            // Check if there is a reconstructed jet with a D0 in the jet radius
            bool passesReco = false;
            if (jet_pt && jet_eta && jet_phi && d0_eta && d0_phi && d0_pt) {
                for (size_t reco_jet_idx = 0; reco_jet_idx < jet_pt->size(); ++reco_jet_idx) {
                    if (!PassesJetSelection(reco_jet_idx, false)) continue;
                    double reco_jet_eta = jet_eta->at(reco_jet_idx);
                    double reco_jet_phi = jet_phi->at(reco_jet_idx);
                    if (HasD0InJetRadius(reco_jet_eta, reco_jet_phi, d0_eta, d0_phi, d0_pt)) {
                        passesReco = true;
                        break;
                    }
                }
            }
            eff_obj->Fill(passesReco, mc_jet_pt, mc_jet_eta);
            if (maxD0Pt > 0) {
                eff_obj_d0pt->Fill(passesReco, maxD0Pt, mc_jet_eta);
            }

            totalMCJetsWithD0++;
        }
        
        // === NUMERATOR: Loop over reconstructed jets with reconstructed D0 ===
        const size_t nRecoJets = jet_pt ? jet_pt->size() : 0;
        for (size_t reco_jet_idx = 0; reco_jet_idx < nRecoJets; ++reco_jet_idx) {
            // Check if reconstructed jet passes selection
            if (!PassesJetSelection(reco_jet_idx, false)) continue;
            
            double reco_jet_pt = jet_pt->at(reco_jet_idx);
            double reco_jet_eta = jet_eta->at(reco_jet_idx);
            double reco_jet_phi = jet_phi->at(reco_jet_idx);
            
            // Check if this reconstructed jet has a reconstructed D0 in the jet radius
            if (!HasD0InJetRadius(reco_jet_eta, reco_jet_phi, d0_eta, d0_phi, d0_pt)) {
                continue;
            }
            
            // Find the MC jet that should correspond to this reco jet using D0 matching
            // Look for a reco D0 in this jet that has an MC match, then get the MC jet
            int mc_match_idx = -1;
            bool foundMatch = false;
            
            // Check each reco D0 to see if it's in this jet and has an MC match
            if (d0_jet_idx) {
                for (size_t d0_idx = 0; d0_idx < d0_jet_idx->size(); ++d0_idx) {
                    if (d0_jet_idx->at(d0_idx) != (int)reco_jet_idx) continue;
                    
                    // This D0 is in our reco jet, see if it has an MC match
                    if (mc_d0_matched) {
                        for (size_t mc_d0_idx = 0; mc_d0_idx < mc_d0_matched->size(); ++mc_d0_idx) {
                            if (mc_d0_matched->at(mc_d0_idx) == (int)d0_idx) {
                                // Found the MC D0 that matches this reco D0
                                mc_match_idx = GetMCJetForD0(mc_d0_idx);
                                if (mc_match_idx >= 0 && PassesJetSelection(mc_match_idx, true)) {
                                    foundMatch = true;
                                    break;
                                }
                            }
                        }
                    }
                    if (foundMatch) break;
                }
            }
            
            if (!foundMatch || mc_match_idx < 0) {
                continue;
            }
            
            // Get MC jet kinematics for filling histograms
            double mc_jet_pt = GetMCJetPt(mc_match_idx);
            double mc_jet_eta = GetMCJetEta(mc_match_idx);
            double mc_jet_phi = GetMCJetPhi(mc_match_idx);
            
            // Verify that the MC jet also has a generator level D0 in jet radius
            if (!HasD0InJetRadius(mc_jet_eta, mc_jet_phi, mc_d0_eta, mc_d0_phi, mc_d0_pt)) {
                continue;
            }
            
            // Find the highest pT D0 in the MC jet radius for D0-pT based histograms
            double maxMCD0Pt = -1.0;
            for (size_t d0_idx = 0; d0_idx < mc_d0_pt->size(); ++d0_idx) {
                double d0_pt = mc_d0_pt->at(d0_idx);
                double d0_eta = mc_d0_eta->at(d0_idx);
                double d0_phi = mc_d0_phi->at(d0_idx);
                
                if (!PassesD0Selection(d0_pt, d0_eta)) continue;
                
                double deltaR = CalculateDeltaR(mc_jet_eta, mc_jet_phi, d0_eta, d0_phi);
                if (deltaR < m_jetRadius && d0_pt > maxMCD0Pt) {
                    maxMCD0Pt = d0_pt;
                }
            }
            
            // Fill numerators (reconstructed jets with reco D0 that have MC match with gen D0)
            h_num->Fill(mc_jet_pt, mc_jet_eta);
            if (maxMCD0Pt > 0) {
                h_num_d0pt->Fill(maxMCD0Pt, mc_jet_eta);
                // Fill zT numerator histogram
                double zT = maxMCD0Pt / mc_jet_pt;
                if (zT > 0 && zT < 2.0) {
                    m_efficiencyMaps["jet_numerator_zT"]->Fill(zT, mc_jet_pt);
                }
            }
            
            totalRecoJetsWithD0++;
            totalMatchedJets++;
        }
    }
    
    // Create efficiency maps
    TH2F *h_eff = (TH2F *)h_num->Clone("jet_finding_efficiency");
    h_eff->Divide(h_den);
    h_eff->SetTitle("Jet Finding Efficiency;Jet p_{T} [GeV];Jet #eta");
    m_efficiencyMaps["jet_finding_efficiency"] = h_eff;
    
    TH2F *h_eff_d0pt = (TH2F *)h_num_d0pt->Clone("jet_finding_efficiency_d0pt");
    h_eff_d0pt->Divide(h_den_d0pt);
    h_eff_d0pt->SetTitle("Jet Finding Efficiency vs D0 pT;D0 p_{T} [GeV];Jet #eta");
    m_efficiencyMaps["jet_finding_efficiency_d0pt"] = h_eff_d0pt;

    // Create zT vs jet pT efficiency map
    TH2F *h_eff_zT = (TH2F *)m_efficiencyMaps["jet_numerator_zT"]->Clone("jet_finding_efficiency_zT");
    h_eff_zT->Divide(m_efficiencyMaps["jet_denominator_zT"]);
    h_eff_zT->SetTitle("Jet Finding Efficiency vs zT;z_{T};Jet p_{T} [GeV]");
    m_efficiencyMaps["jet_finding_efficiency_zT"] = h_eff_zT;
    
    std::cout << "Jet finding efficiency calculation completed" << std::endl;
    
    // Print statistics
    int totalDen = h_den->GetEntries();
    int totalNum = h_num->GetEntries();
    double overallEff = (totalDen > 0) ? (double)totalNum / totalDen : 0.0;
    
    std::cout << "\nStatistics:" << std::endl;
    std::cout << "  Total MC jets with generator D0 in jet radius (denominator): " << totalDen << std::endl;
    std::cout << "  Total reco jets with reco D0 and MC match with gen D0 (numerator): " << totalNum << std::endl;
    std::cout << "  Overall jet finding efficiency: " << overallEff * 100.0 << "% (" << totalNum << "/" << totalDen << ")" << std::endl;
    std::cout << "\nDetailed breakdown:" << std::endl;
    std::cout << "  MC jets with generator D0 processed: " << totalMCJetsWithD0 << std::endl;
    std::cout << "  Reconstructed jets with reco D0: " << totalRecoJetsWithD0 << std::endl;
    std::cout << "  Successfully matched jets: " << totalMatchedJets << std::endl;
}

void JetFindingEfficiencyStandalone::SaveResults() {
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

void JetFindingEfficiencyStandalone::PlotEfficiency(const std::string &histName) {
    TH2F *h_eff = m_efficiencyMaps[histName];
    if (!h_eff) {
        std::cerr << "Error: Efficiency histogram '" << histName << "' not found" << std::endl;
        return;
    }
    
    std::cout << "Plotting efficiency: " << histName << std::endl;
    
    // Set up canvas
    TCanvas *c1 = new TCanvas("c_jet_eff", "Jet Finding Efficiency", 800, 600);
    c1->SetRightMargin(0.15);
    
    // Set up histogram style
    h_eff->SetStats(0);
    h_eff->GetZaxis()->SetTitle("Efficiency");
    h_eff->GetZaxis()->SetRangeUser(0.0, 1.0);
    
    // Draw efficiency map
    h_eff->Draw("COLZ");
    
    // Add text with overall efficiency
    TH2F *h_num = m_efficiencyMaps["jet_numerator"];
    TH2F *h_den = m_efficiencyMaps["jet_denominator"];
    if (h_num && h_den) {
        double overallEff = (double)h_num->GetEntries() / h_den->GetEntries();
        TString effText = TString::Format("Overall Jet Finding Efficiency: %.3f", overallEff);
        
        TLatex *latex = new TLatex();
        latex->SetNDC();
        latex->SetTextSize(0.04);
        latex->DrawLatex(0.15, 0.92, effText);
        
        TString jetRadiusText = TString::Format("Jet Radius: R = %.1f", m_jetRadius);
        latex->DrawLatex(0.15, 0.87, jetRadiusText);
    }
    
    // Add label for zT plot
    if (histName == "jet_finding_efficiency_zT") {
        TLatex *latex = new TLatex();
        latex->SetNDC();
        latex->SetTextSize(0.04);
        latex->DrawLatex(0.15, 0.92, "Jet Finding Efficiency vs z_{T}");
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

// Main function for the standalone jet finding efficiency calculation
int JetFindingEfficiencyStandaloneRun(TString inputFile = "/media/niviths/SSD2/lhcb_analysis_SSD/20250728_pPb_MC_output/20250728_pPb_MC_output.root", 
                                     TString outputFile = "output_jet_finding_efficiency.root",
                                     double jetRadius = 0.4, double minJetPt = 5.0, double maxJetPt = 60.0,
                                     double minJetEta = 2.5, double maxJetEta = 4.0,
                                     double minD0Pt = 1.0, double maxD0Pt = 50.0,
                                     double minD0Eta = 2.0, double maxD0Eta = 4.5,
                                     bool makePlots = true) {
    std::cout << "=== Jet Finding Efficiency Calculator (Standalone) ===" << std::endl;
    std::cout << "Input file: " << inputFile << std::endl;
    std::cout << "Output file: " << outputFile << std::endl;
    std::cout << std::endl;
    std::cout << "Configuration:" << std::endl;
    std::cout << "  Jet radius: " << jetRadius << std::endl;
    std::cout << "  Jet pT range: " << minJetPt << " - " << maxJetPt << " GeV" << std::endl;
    std::cout << "  Jet eta range: " << minJetEta << " - " << maxJetEta << std::endl;
    std::cout << "  D0 pT range: " << minD0Pt << " - " << maxD0Pt << " GeV" << std::endl;
    std::cout << "  D0 eta range: " << minD0Eta << " - " << maxD0Eta << std::endl;
    std::cout << "  Generate plots: " << (makePlots ? "Yes" : "No") << std::endl;
    std::cout << std::endl;
    
    std::cout << "Method: Jet finding efficiency calculation" << std::endl;
    std::cout << "  - Denominator: MC jets with generator level D0 in jet radius" << std::endl;
    std::cout << "  - Numerator: Reconstructed jets with reconstructed D0 and MC match" << std::endl;
    std::cout << std::endl;
    
    try {
        // Create efficiency calculator
        JetFindingEfficiencyStandalone calculator(inputFile, outputFile);
        
        // Configure parameters
        calculator.SetJetRadius(jetRadius);
        calculator.SetJetPtRange(minJetPt, maxJetPt);
        calculator.SetJetEtaRange(minJetEta, maxJetEta);
        calculator.SetD0PtRange(minD0Pt, maxD0Pt);
        calculator.SetD0EtaRange(minD0Eta, maxD0Eta);
        
        // Set up custom binning
        // std::vector<double> customJetPtBins;
        // for (double pt = 5.0; pt < 10.0; pt += 1.0)
        //     customJetPtBins.push_back(pt);
        // for (double pt = 10.0; pt < 20.0; pt += 2.0)
        //     customJetPtBins.push_back(pt);
        // for (double pt = 20.0; pt <= 100.0; pt += 5.0)
        //     customJetPtBins.push_back(pt);
        
        // std::vector<double> customJetEtaBins;
        // for (double eta = 2.0; eta <= 4.5; eta += 0.25)
        //     customJetEtaBins.push_back(eta);
        
        // std::vector<double> customD0PtBins;
        // for (double pt = 2.0; pt < 10.0; pt += 1.0)
        //     customD0PtBins.push_back(pt);
        // for (double pt = 10.0; pt <= 50.0; pt += 2.5)
        //     customD0PtBins.push_back(pt);
        
        // calculator.SetJetPtBins(customJetPtBins);
        // calculator.SetJetEtaBins(customJetEtaBins);
        // calculator.SetD0PtBins(customD0PtBins);
        
        // Initialize
        if (!calculator.Initialize()) {
            std::cerr << "Error: Failed to initialize calculator" << std::endl;
            return 1;
        }
        
        // Calculate efficiency
        std::cout << "Starting jet finding efficiency calculation..." << std::endl;
        calculator.CalculateJetFindingEfficiency();
        
        // Save results
        calculator.SaveResults();
        
        // Generate plots if requested
        if (makePlots) {
            std::cout << "\nGenerating efficiency plots..." << std::endl;
            calculator.PlotEfficiency("jet_finding_efficiency");
            calculator.PlotEfficiency("jet_finding_efficiency_d0pt");
            calculator.PlotEfficiency("jet_finding_efficiency_zT");
        }
        
        std::cout << "\n=== Jet Finding Efficiency calculation completed successfully! ===" << std::endl;
    }
    catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
