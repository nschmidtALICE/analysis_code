
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
#include <sys/stat.h>
#include <sys/types.h>
#include "TLegend.h"

// Standalone Jet Finding Efficiency Calculator
// This class calculates jet finding efficiency using reconstructed and MC jets/D0s from a ROOT TTree.
// It fills histograms and efficiency objects, and provides plotting and saving routines.
class JetFindingEfficiencyStandalone
{
public:
    JetFindingEfficiencyStandalone(TString inputFileName, TString outputFileName);
    ~JetFindingEfficiencyStandalone();

    bool Initialize();
    void CalculateJetFindingEfficiency();
    void SaveResults();
    void PlotEfficiency(const std::string &histName);

    // Configuration methods
    void SetJetRadius(double radius) { m_jetRadius = radius; }
    void SetJetPtRange(double minPt, double maxPt = 1000.0)
    {
        m_minJetPt = minPt;
        m_maxJetPt = maxPt;
    }
    void SetJetEtaRange(double minEta, double maxEta)
    {
        m_minJetEta = minEta;
        m_maxJetEta = maxEta;
    }
    void SetD0PtRange(double minPt, double maxPt = 1000.0)
    {
        m_minD0Pt = minPt;
        m_maxD0Pt = maxPt;
    }
    void SetD0EtaRange(double minEta, double maxEta)
    {
        m_minD0Eta = minEta;
        m_maxD0Eta = maxEta;
    }
    void SetJetPtBins(const std::vector<double> &bins) { m_jetPtBins = bins; }
    void SetJetEtaBins(const std::vector<double> &bins) { m_jetEtaBins = bins; }
    void SetD0PtBins(const std::vector<double> &bins) { m_d0PtBins = bins; }

    // Add getter for efficiency maps
    TH2F *GetEfficiencyHistogram(const std::string &name) const
    {
        auto it = m_efficiencyMaps.find(name);
        if (it != m_efficiencyMaps.end())
            return it->second;
        return nullptr;
    }

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
    std::map<std::string, TH2F *> m_efficiencyMaps;
    std::map<std::string, TEfficiency *> m_efficiencyObjects;

    // Tree branches - Reconstructed jets
    std::vector<float> *jet_pt;
    std::vector<float> *jet_eta;
    std::vector<float> *jet_phi;
    std::vector<int> *jet_n_d0;

    // Tree branches - Reconstructed D0s
    std::vector<float> *d0_pt;
    std::vector<float> *d0_eta;
    std::vector<float> *d0_phi;
    std::vector<float> *d0_mass;
    std::vector<float> *d0_pz;
    std::vector<float> *d0_e;
    std::vector<int> *d0_jet_idx;
    std::vector<float> *d0_jet_dr; // reco D0-jet distance

    // Tree branches - MC truth jets (derived from px,py,pz,e)
    std::vector<float> *mc_jet_pt;
    std::vector<float> *mc_jet_eta;
    std::vector<float> *mc_jet_phi;

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
    std::vector<float> *mc_d0_jet_dr; // MC D0-jet distance

    void InitializeBranches();
    void CreateHistograms();
    void CleanUp();
    double CalculateDeltaR(double eta1, double phi1, double eta2, double phi2);
    bool HasD0InJetRadius(bool isMC, int jetIdx, int d0Idx);
    bool PassesJetSelection(int jet_idx, bool isMC = false);
    bool PassesD0Selection(double d0_pt, double d0_eta);
    int GetMCJetForD0(int mc_d0_idx);
};

JetFindingEfficiencyStandalone::JetFindingEfficiencyStandalone(TString inputFileName, TString outputFileName)
    : m_inputFile(nullptr), m_outputFile(nullptr), m_tree(nullptr),
      m_jetRadius(0.4), m_minJetPt(5.0), m_maxJetPt(1000.0),
      m_minJetEta(2.0), m_maxJetEta(4.5),
      m_minD0Pt(2.0), m_maxD0Pt(1000.0),
      m_minD0Eta(2.0), m_maxD0Eta(4.5)
{
    // Open input ROOT file containing the TTree with jets and D0s
    m_inputFile = TFile::Open(inputFileName, "READ");
    if (!m_inputFile || m_inputFile->IsZombie())
    {
        std::cerr << "Error: Could not open input file " << inputFileName << std::endl;
        return;
    }
    // Create output ROOT file for results
    m_outputFile = new TFile(outputFileName, "RECREATE");
    if (!m_outputFile || m_outputFile->IsZombie())
    {
        std::cerr << "Error: Could not create output file " << outputFileName << std::endl;
        return;
    }
    // Get analysis tree from input file
    m_tree = (TTree *)m_inputFile->Get("d0jets");
    if (!m_tree)
    {
        std::cerr << "Error: Could not find tree 'd0jets' in input file" << std::endl;
        return;
    }
    // Initialize pointers to branch vectors
    jet_pt = nullptr;
    jet_eta = nullptr;
    jet_phi = nullptr;
    jet_n_d0 = nullptr;
    d0_pt = nullptr;
    d0_eta = nullptr;
    d0_phi = nullptr;
    d0_mass = nullptr;
    d0_pz = nullptr;
    d0_e = nullptr;
    d0_jet_idx = nullptr;
    d0_jet_dr = nullptr;
    mc_jet_pt = nullptr;
    mc_jet_eta = nullptr;
    mc_jet_phi = nullptr;
    mc_d0_pid = nullptr;
    mc_d0_pt = nullptr;
    mc_d0_eta = nullptr;
    mc_d0_phi = nullptr;
    mc_d0_mass = nullptr;
    mc_d0_origin = nullptr;
    mc_d0_matched = nullptr;
    mc_d0_jet_idx = nullptr;
    mc_d0_z = nullptr;
    mc_d0_jet_dr = nullptr;
    // Set default binning for jet pT, eta, and D0 pT
    m_jetPtBins.clear();
    for (double pt = 0.0; pt <= 60.0; pt += 2.5)
        m_jetPtBins.push_back(pt);
    m_jetEtaBins.clear();
    for (double eta = 2.0; eta <= 4.5; eta += 0.2)
        m_jetEtaBins.push_back(eta);
    m_d0PtBins = m_jetPtBins;
}

JetFindingEfficiencyStandalone::~JetFindingEfficiencyStandalone()
{
    CleanUp();
}

void JetFindingEfficiencyStandalone::CleanUp()
{
    for (auto &pair : m_efficiencyMaps)
    {
        if (pair.second)
            delete pair.second;
    }
    m_efficiencyMaps.clear();

    for (auto &pair : m_efficiencyObjects)
    {
        if (pair.second)
            delete pair.second;
    }
    m_efficiencyObjects.clear();

    if (m_outputFile)
    {
        m_outputFile->Close();
        delete m_outputFile;
        m_outputFile = nullptr;
    }

    if (m_inputFile)
    {
        m_inputFile->Close();
        m_inputFile = nullptr;
    }
}

bool JetFindingEfficiencyStandalone::Initialize()
{
    if (!m_inputFile || !m_outputFile || !m_tree)
    {
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

void JetFindingEfficiencyStandalone::InitializeBranches()
{
    // Set branch addresses for reconstructed jets
    m_tree->SetBranchAddress("jet_pt", &jet_pt);
    m_tree->SetBranchAddress("jet_eta", &jet_eta);
    m_tree->SetBranchAddress("jet_phi", &jet_phi);
    m_tree->SetBranchAddress("jet_n_d0", &jet_n_d0);
    // Set branch addresses for reconstructed D0s
    m_tree->SetBranchAddress("d0_pt", &d0_pt);
    m_tree->SetBranchAddress("d0_eta", &d0_eta);
    m_tree->SetBranchAddress("d0_phi", &d0_phi);
    m_tree->SetBranchAddress("d0_mass", &d0_mass);
    m_tree->SetBranchAddress("d0_pz", &d0_pz);
    m_tree->SetBranchAddress("d0_e", &d0_e);
    m_tree->SetBranchAddress("d0_jet_idx", &d0_jet_idx);
    m_tree->SetBranchAddress("d0_jet_dr", &d0_jet_dr);
    // Set branch addresses for MC truth jets
    m_tree->SetBranchAddress("mc_jet_pt", &mc_jet_pt);
    m_tree->SetBranchAddress("mc_jet_eta", &mc_jet_eta);
    m_tree->SetBranchAddress("mc_jet_phi", &mc_jet_phi);
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
    m_tree->SetBranchAddress("mc_d0_jet_dr", &mc_d0_jet_dr);
}

void JetFindingEfficiencyStandalone::CreateHistograms()
{
    int nZTbins = 20; // You can adjust this binning as needed
    double zTmin = 0.0, zTmax = 1.0;
    // Create jet finding efficiency histograms (jet pT vs jet eta)
    m_efficiencyMaps["jet_numerator"] = new TH2F("jet_numerator", "Jet Finding Efficiency Numerator;Jet p_{T} [GeV];Jet #eta",
                                                 m_jetPtBins.size() - 1, &m_jetPtBins[0],
                                                 m_jetEtaBins.size() - 1, &m_jetEtaBins[0]);
    m_efficiencyMaps["jet_denominator"] = new TH2F("jet_denominator", "Jet Finding Efficiency Denominator;Jet p_{T} [GeV];Jet #eta",
                                                   m_jetPtBins.size() - 1, &m_jetPtBins[0],
                                                   m_jetEtaBins.size() - 1, &m_jetEtaBins[0]);
    m_efficiencyMaps["reco_zT_vs_jetPt"] = new TH2F("reco_zT_vs_jetPt", "Reco z_{T} vs Jet p_{T};z_{T};Jet p_{T} [GeV]", nZTbins, zTmin, zTmax, m_jetPtBins.size() - 1, &m_jetPtBins[0]);
    // Add MC yield histogram: mcPt vs mczT for generator-level jets with reconstructed tag
    m_efficiencyMaps["mc_genYield_mcPt_mczT"] = new TH2F("mc_genYield_mcPt_mczT", "MC Gen Yield (Tag Reco);z_{T}^{MC};Jet p_{T}^{MC} [GeV]", nZTbins, zTmin, zTmax, m_jetPtBins.size() - 1, &m_jetPtBins[0]);

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

    m_efficiencyMaps["jet_numerator_zT"] = new TH2F("jet_numerator_zT", "Jet Finding Efficiency Numerator;z_{T};Jet p_{T} [GeV]",
                                                    nZTbins, zTmin, zTmax, nJetPtBins, &m_jetPtBins[0]);
    m_efficiencyMaps["jet_denominator_zT"] = new TH2F("jet_denominator_zT", "Jet Finding Efficiency Denominator;z_{T};Jet p_{T} [GeV]",
                                                      nZTbins, zTmin, zTmax, nJetPtBins, &m_jetPtBins[0]);

    std::cout << "Created efficiency histograms with:" << std::endl;
    std::cout << "  Jet pT bins: " << m_jetPtBins.size() - 1 << " (" << m_jetPtBins.front() << " - " << m_jetPtBins.back() << " GeV)" << std::endl;
    std::cout << "  Jet eta bins: " << m_jetEtaBins.size() - 1 << " (" << m_jetEtaBins.front() << " - " << m_jetEtaBins.back() << ")" << std::endl;
    std::cout << "  D0 pT bins: " << m_d0PtBins.size() - 1 << " (" << m_d0PtBins.front() << " - " << m_d0PtBins.back() << " GeV)" << std::endl;
}

// Calculate the angular distance DeltaR between two objects
double JetFindingEfficiencyStandalone::CalculateDeltaR(double eta1, double phi1, double eta2, double phi2)
{
    double deta = eta1 - eta2;
    double dphi = TMath::Abs(phi1 - phi2);
    if (dphi > TMath::Pi())
        dphi = 2 * TMath::Pi() - dphi;
    return TMath::Sqrt(deta * deta + dphi * dphi);
}

// Check if a D0 candidate is within the jet radius of a given jet
bool JetFindingEfficiencyStandalone::HasD0InJetRadius(bool isMC, int jetIdx, int d0Idx)
{
    // Use only d0_jet_dr or mc_d0_jet_dr, no fallback
    if (isMC)
    {
        if (!mc_d0_jet_dr)
            return false;
        if (d0Idx < 0 || d0Idx >= (int)mc_d0_jet_dr->size())
            return false;
        return mc_d0_jet_dr->at(d0Idx) < m_jetRadius && mc_d0_jet_dr->at(d0Idx) >= 0;
    }
    else
    {
        if (!d0_jet_dr)
            return false;
        if (d0Idx < 0 || d0Idx >= (int)d0_jet_dr->size())
            return false;
        return d0_jet_dr->at(d0Idx) < m_jetRadius && d0_jet_dr->at(d0Idx) >= 0;
    }
}

// Get the MC jet index associated with a given MC D0
int JetFindingEfficiencyStandalone::GetMCJetForD0(int mc_d0_idx)
{
    if (!mc_d0_jet_idx || mc_d0_idx < 0 || mc_d0_idx >= (int)mc_d0_jet_idx->size())
    {
        return -1;
    }
    return mc_d0_jet_idx->at(mc_d0_idx);
}

// Apply kinematic selection cuts to jets
bool JetFindingEfficiencyStandalone::PassesJetSelection(int jet_idx, bool isMC)
{
    if (isMC)
    {
        // For MC jets, use direct branches
        if (!mc_jet_pt || !mc_jet_eta || !mc_jet_phi || jet_idx < 0 ||
            jet_idx >= (int)mc_jet_pt->size() || jet_idx >= (int)mc_jet_eta->size() ||
            jet_idx >= (int)mc_jet_phi->size())
        {
            return false;
        }
        double pt = mc_jet_pt->at(jet_idx);
        double eta = mc_jet_eta->at(jet_idx);
        // Basic kinematic and acceptance cuts
        if (pt < 0 || eta < -900)
            return false; // Invalid kinematics
        if (pt < m_minJetPt || pt > m_maxJetPt)
            return false;
        if (eta < m_minJetEta || eta > m_maxJetEta)
            return false;
    }
    else
    {
        // For reconstructed jets, use direct branches
        if (!jet_pt || !jet_eta || jet_idx < 0 ||
            jet_idx >= (int)jet_pt->size() || jet_idx >= (int)jet_eta->size())
        {
            return false;
        }
        double pt = jet_pt->at(jet_idx);
        double eta = jet_eta->at(jet_idx);
        if (pt < m_minJetPt || pt > m_maxJetPt)
            return false;
        if (eta < m_minJetEta || eta > m_maxJetEta)
            return false;
    }
    return true;
}

// Apply kinematic selection cuts to D0 candidates
bool JetFindingEfficiencyStandalone::PassesD0Selection(double d0_pt, double d0_eta)
{
    if (d0_pt < m_minD0Pt || d0_pt > m_maxD0Pt)
        return false;
    if (d0_eta < m_minD0Eta || d0_eta > m_maxD0Eta)
        return false;
    return true;
}

// Main event loop: calculate jet finding efficiency and fill all histograms
void JetFindingEfficiencyStandalone::CalculateJetFindingEfficiency()
{
    std::cout << "\nCalculating jet finding efficiency..." << std::endl;
    // Get efficiency histograms and objects
    TH2F *h_num = m_efficiencyMaps["jet_numerator"];
    TH2F *h_den = m_efficiencyMaps["jet_denominator"];
    TEfficiency *eff_obj = m_efficiencyObjects["jet_finding_efficiency"];
    TH2F *h_num_d0pt = m_efficiencyMaps["jet_numerator_d0pt"];
    TH2F *h_den_d0pt = m_efficiencyMaps["jet_denominator_d0pt"];
    TEfficiency *eff_obj_d0pt = m_efficiencyObjects["jet_finding_efficiency_d0pt"];
    if (!h_num || !h_den || !eff_obj || !h_num_d0pt || !h_den_d0pt || !eff_obj_d0pt)
    {
        std::cerr << "Error: Efficiency histograms/objects not available" << std::endl;
        return;
    }
    // Loop over all events in the TTree
    Long64_t nEntries = m_tree->GetEntries();
    std::cout << "Total entries in tree: " << m_tree->GetEntries() << std::endl;
    std::cout << "Processing " << nEntries << " events..." << std::endl;
    // Counters for statistics
    int totalMCJetsWithD0 = 0;
    int totalRecoJetsWithD0 = 0;
    int totalMatchedJets = 0;
    for (Long64_t entry = 0; entry < nEntries; ++entry)
    {
        m_tree->GetEntry(entry);
        if (entry % 50000 == 0)
        {
            std::cout << "Processing entry " << entry << "/" << nEntries << " (" << (100.0 * entry / nEntries) << "%)" << std::endl;
        }
        const size_t nMCJets = mc_jet_pt ? mc_jet_pt->size() : 0;
        const size_t nRecoJets = jet_pt ? jet_pt->size() : 0;
        // Cache D0 info for MC jets (for denominator and MC gen yield)
        // Each MC jet: pair of (max D0 pT in jet, has reco tag)
        std::vector<std::pair<double, bool>> mcJet_maxD0Pt_hasRecoTag(nMCJets, std::make_pair(-1.0, false));
        for (size_t mc_d0_idx = 0; mc_d0_idx < (mc_d0_jet_idx ? mc_d0_jet_idx->size() : 0); ++mc_d0_idx)
        {
            int mc_jet_idx = mc_d0_jet_idx->at(mc_d0_idx);
            if (mc_jet_idx < 0 || mc_jet_idx >= (int)nMCJets)
                continue;
            if (!HasD0InJetRadius(true, mc_jet_idx, mc_d0_idx))
                continue;
            double d0_pt_val = mc_d0_pt->at(mc_d0_idx);
            double d0_eta_val = mc_d0_eta ? mc_d0_eta->at(mc_d0_idx) : -999.0;
            // apply D0 selection to generator-level D0s used for denominator
            if (!PassesD0Selection(d0_pt_val, d0_eta_val))
                continue;
            if (d0_pt_val > mcJet_maxD0Pt_hasRecoTag[mc_jet_idx].first)
                mcJet_maxD0Pt_hasRecoTag[mc_jet_idx].first = d0_pt_val;
            if (mc_d0_matched->at(mc_d0_idx) >= 0)
                mcJet_maxD0Pt_hasRecoTag[mc_jet_idx].second = true;
        }
        // Cache D0 info for reco jets (for reco zT vs jet pT)
        // Each reco jet: max D0 pT in jet
        std::vector<double> recoJet_maxD0Pt(nRecoJets, -1.0);
        for (size_t d0_idx = 0; d0_idx < (d0_jet_idx ? d0_jet_idx->size() : 0); ++d0_idx)
        {
            int reco_jet_idx = d0_jet_idx->at(d0_idx);
            if (reco_jet_idx < 0 || reco_jet_idx >= (int)nRecoJets)
                continue;
            if (!HasD0InJetRadius(false, reco_jet_idx, d0_idx))
                continue;
            double d0_pt_val = d0_pt->at(d0_idx);
            double d0_eta_val = d0_eta ? d0_eta->at(d0_idx) : -999.0;
            // apply D0 selection to reconstructed D0s used for numerator
            if (!PassesD0Selection(d0_pt_val, d0_eta_val))
                continue;
            if (d0_pt_val > recoJet_maxD0Pt[reco_jet_idx])
                recoJet_maxD0Pt[reco_jet_idx] = d0_pt_val;
        }
        // === DENOMINATOR: Loop over MC jets with generator level D0 ===
        for (size_t mc_jet_idx = 0; mc_jet_idx < nMCJets; ++mc_jet_idx)
        {
            double mcJetPt = mc_jet_pt ? mc_jet_pt->at(mc_jet_idx) : -1.0;
            double mcJetEta = mc_jet_eta ? mc_jet_eta->at(mc_jet_idx) : -999.0;
            double mcJetPhi = mc_jet_phi ? mc_jet_phi->at(mc_jet_idx) : -999.0;
            if (mcJetPt < 0 || mcJetEta < -900)
                continue;
            if (!PassesJetSelection(mc_jet_idx, true))
                continue;
            double maxD0Pt = mcJet_maxD0Pt_hasRecoTag[mc_jet_idx].first;
            bool hasGenD0InMCJetDen = maxD0Pt > 0;
            if (!hasGenD0InMCJetDen)
                continue;
            // Fill denominator histograms
            h_den->Fill(mcJetPt, mcJetEta);
            h_den_d0pt->Fill(maxD0Pt, mcJetEta);
            double zT = maxD0Pt / mcJetPt;
            if (zT > 0 && zT < 2.0)
                m_efficiencyMaps["jet_denominator_zT"]->Fill(zT, mcJetPt);
            // Determine whether THIS MC jet was reconstructed:
            // Look for MC D0s that belong to this MC jet, check if any is matched
            // to a reconstructed D0 that sits in a reconstructed jet which passes selection.
            bool passesReco = false;
            if (mc_d0_matched && mc_d0_jet_idx && d0_jet_idx)
            {
                for (size_t mc_d0_i = 0; mc_d0_i < mc_d0_matched->size(); ++mc_d0_i)
                {
                    // only consider MC D0s belonging to this MC jet
                    int mcJetForThisD0 = GetMCJetForD0((int)mc_d0_i);
                    if (mcJetForThisD0 != (int)mc_jet_idx)
                        continue;
                    int recoD0Idx = mc_d0_matched->at(mc_d0_i);
                    if (recoD0Idx < 0)
                        continue;
                    // check that the matched reco D0 was assigned to a reco jet
                    if (recoD0Idx >= (int)d0_jet_idx->size())
                        continue;
                    int recoJetForMatchedD0 = d0_jet_idx->at(recoD0Idx);
                    if (recoJetForMatchedD0 < 0 || recoJetForMatchedD0 >= (int)recoJet_maxD0Pt.size())
                        continue;
                    if (!PassesJetSelection(recoJetForMatchedD0, false))
                        continue;
                    if (recoJet_maxD0Pt[recoJetForMatchedD0] > 0)
                    {
                        passesReco = true;
                        break;
                    }
                }
            }
            eff_obj->Fill(passesReco, mcJetPt, mcJetEta);
            eff_obj_d0pt->Fill(passesReco, maxD0Pt, mcJetEta);
            totalMCJetsWithD0++;
        }
        // === NUMERATOR: Loop over reconstructed jets with reconstructed D0 ===
        // track which MC jets have already been counted in the numerator for this event
        std::vector<char> mcJetAlreadyCounted(nMCJets, 0);
        for (size_t reco_jet_idx = 0; reco_jet_idx < nRecoJets; ++reco_jet_idx)
        {
            if (!PassesJetSelection(reco_jet_idx, false))
                continue;
            double reco_jet_pt = jet_pt->at(reco_jet_idx);
            double reco_jet_eta = jet_eta->at(reco_jet_idx);
            double reco_jet_phi = jet_phi->at(reco_jet_idx);
            if (recoJet_maxD0Pt[reco_jet_idx] <= 0)
                continue;
            // Find MC match for D0 in this reco jet
            int mc_match_idx = -1;
            bool foundMatch = false;
            if (d0_jet_idx)
            {
                for (size_t d0_idx = 0; d0_idx < d0_jet_idx->size(); ++d0_idx)
                {
                    if (d0_jet_idx->at(d0_idx) != (int)reco_jet_idx)
                        continue;
                    if (mc_d0_matched)
                    {
                        for (size_t mc_d0_idx = 0; mc_d0_idx < mc_d0_matched->size(); ++mc_d0_idx)
                        {
                            if (mc_d0_matched->at(mc_d0_idx) == (int)d0_idx)
                            {
                                mc_match_idx = GetMCJetForD0(mc_d0_idx);
                                if (mc_match_idx >= 0 && PassesJetSelection(mc_match_idx, true))
                                {
                                    foundMatch = true;
                                    break;
                                }
                            }
                        }
                    }
                    if (foundMatch)
                        break;
                }
            }
            if (!foundMatch || mc_match_idx < 0)
                continue;
            double mcJetPt = mc_jet_pt ? mc_jet_pt->at(mc_match_idx) : -1.0;
            double mcJetEta = mc_jet_eta ? mc_jet_eta->at(mc_match_idx) : -999.0;
            double mcJetPhi = mc_jet_phi ? mc_jet_phi->at(mc_match_idx) : -999.0;
            double maxMCD0Pt = mcJet_maxD0Pt_hasRecoTag[mc_match_idx].first;
            bool hasGenD0InMCJet = maxMCD0Pt > 0;
            if (!hasGenD0InMCJet)
                continue;
            // Fill numerator histograms, but avoid counting the same MC jet multiple times in one event
            if (mc_match_idx >= 0 && mc_match_idx < (int)mcJetAlreadyCounted.size() && !mcJetAlreadyCounted[mc_match_idx])
            {
                h_num->Fill(mcJetPt, mcJetEta);
                h_num_d0pt->Fill(maxMCD0Pt, mcJetEta);
                double zT = maxMCD0Pt / mcJetPt;
                if (zT > 0 && zT < 2.0)
                    m_efficiencyMaps["jet_numerator_zT"]->Fill(zT, mcJetPt);
                totalMatchedJets++;
                mcJetAlreadyCounted[mc_match_idx] = 1;
            }
            // always count reco jets with a D0 (for statistics)
            totalRecoJetsWithD0++;
        }
        // Fill reco zT vs jet pT histogram (single loop)
        TH2F *h_reco_zT_vs_jetPt = m_efficiencyMaps["reco_zT_vs_jetPt"];
        for (size_t reco_jet_idx = 0; reco_jet_idx < nRecoJets; ++reco_jet_idx)
        {
            if (!PassesJetSelection(reco_jet_idx, false))
                continue;
            double reco_jet_pt = jet_pt->at(reco_jet_idx);
            double maxRecoD0Pt = recoJet_maxD0Pt[reco_jet_idx];
            if (maxRecoD0Pt > 0 && reco_jet_pt > 0)
            {
                double zT = maxRecoD0Pt / reco_jet_pt;
                if (zT > 0 && zT < 2.0)
                    h_reco_zT_vs_jetPt->Fill(zT, reco_jet_pt);
            }
        }
        // Fill MC gen yield histogram for jets with reconstructed tag (single loop)
        TH2F *h_mc_genYield_mcPt_mczT = m_efficiencyMaps["mc_genYield_mcPt_mczT"];
        for (size_t mc_jet_idx = 0; mc_jet_idx < nMCJets; ++mc_jet_idx)
        {
            double mcJetPt = mc_jet_pt->at(mc_jet_idx);
            double maxMCD0Pt = mcJet_maxD0Pt_hasRecoTag[mc_jet_idx].first;
            bool hasRecoTag = mcJet_maxD0Pt_hasRecoTag[mc_jet_idx].second;
            if (hasRecoTag && maxMCD0Pt > 0 && mcJetPt > 0)
            {
                double mczT = maxMCD0Pt / mcJetPt;
                if (mczT > 0 && mczT < 2.0)
                    h_mc_genYield_mcPt_mczT->Fill(mczT, mcJetPt);
            }
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

void JetFindingEfficiencyStandalone::SaveResults()
{
    if (!m_outputFile)
    {
        std::cerr << "Error: Output file not available" << std::endl;
        return;
    }

    std::cout << "\nSaving results to " << m_outputFile->GetName() << "..." << std::endl;

    m_outputFile->cd();

    // Save all efficiency histograms
    for (auto &pair : m_efficiencyMaps)
    {
        if (pair.second)
        {
            pair.second->Write();
            std::cout << "  Saved histogram: " << pair.first << std::endl;
        }
    }

    // Save all efficiency objects
    for (auto &pair : m_efficiencyObjects)
    {
        if (pair.second)
        {
            pair.second->Write();
            std::cout << "  Saved efficiency object: " << pair.first << std::endl;
        }
    }

    m_outputFile->Write();
    std::cout << "Results saved successfully" << std::endl;
}

void ensureDirectoryExists(const std::string &dirName)
{
    struct stat info;
    if (stat(dirName.c_str(), &info) != 0)
    {
        mkdir(dirName.c_str(), 0777);
    }
}

void JetFindingEfficiencyStandalone::PlotEfficiency(const std::string &histName)
{

    // Plot all jet finding efficiency zT projections together
    if (histName == "jet_finding_efficiency_zT_allProjections")
    {
        TH2F *h_num_zT = m_efficiencyMaps["jet_numerator_zT"];
        TH2F *h_den_zT = m_efficiencyMaps["jet_denominator_zT"];
        if (!h_num_zT || !h_den_zT)
        {
            std::cerr << "Error: Numerator/denominator zT histograms not found" << std::endl;
            return;
        }
        std::vector<double> startPt = {5, 10, 15, 20, 30};
        std::vector<double> endPt = {10, 15, 20, 30, 50};
        std::vector<int> colors = {kRed, kBlue, kGreen + 2, kMagenta, kOrange + 7};
        TCanvas *cAll = new TCanvas("c_jet_eff_allProj", "Jet Finding Efficiency vs zT (All pT Slices)", 800, 600);
        cAll->SetRightMargin(0.15);
        TLegend *leg = new TLegend(0.55, 0.15, 0.85, 0.45);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        bool first = true;
        for (size_t i = 0; i < startPt.size(); ++i)
        {
            double ptMin = startPt[i];
            double ptMax = endPt[i];
            int binMin = h_num_zT->GetYaxis()->FindBin(ptMin);
            int binMax = h_num_zT->GetYaxis()->FindBin(ptMax) - 1;
            TString sliceName = TString::Format("jet_finding_efficiency_zT_pt_%g_%g", ptMin, ptMax);
            TH1D *num_proj = h_num_zT->ProjectionX(sliceName + "_num", binMin, binMax);
            TH1D *den_proj = h_den_zT->ProjectionX(sliceName + "_den", binMin, binMax);
            TH1D *eff_proj = (TH1D *)num_proj->Clone(sliceName);
            eff_proj->SetTitle(TString::Format("Jet Finding Efficiency vs z_{T} (%g < p_{T} < %g GeV)", ptMin, ptMax));
            eff_proj->SetStats(0);
            eff_proj->GetYaxis()->SetTitle("Efficiency");
            eff_proj->GetXaxis()->SetTitle("z_{T}");
            eff_proj->SetMinimum(0.0);
            eff_proj->SetMaximum(1.0);
            eff_proj->Divide(eff_proj, den_proj, 1.0, 1.0, "B");
            eff_proj->SetLineColor(colors[i % colors.size()]);
            eff_proj->SetLineWidth(2);
            if (first)
            {
                eff_proj->Draw("E1");
                first = false;
            }
            else
            {
                eff_proj->Draw("E1 SAME");
            }
            leg->AddEntry(eff_proj, TString::Format("%g < p_{T} < %g GeV", ptMin, ptMax), "l");
            delete num_proj;
            delete den_proj;
        }
        leg->Draw();
        TLatex *latex = new TLatex();
        latex->SetNDC();
        latex->SetTextSize(0.04);
        latex->DrawLatex(0.15, 0.92, "Jet Finding Efficiency vs z_{T} (All p_{T} Slices)");
        latex->DrawLatex(0.15, 0.85, "#font[12]{LHCb} work-in-progress");
        latex->DrawLatex(0.15, 0.80, TString::Format("Anti #it{k}_{T}, R = %.1f", m_jetRadius));
        std::string outDir = "jetFindingEfficiencies";
        ensureDirectoryExists(outDir);
        TString plotName = TString::Format("%s/%s.png", outDir.c_str(), histName.c_str());
        cAll->SaveAs(plotName);
        plotName = TString::Format("%s/%s.pdf", outDir.c_str(), histName.c_str());
        cAll->SaveAs(plotName);
        std::cout << "Saved plot: " << plotName << std::endl;
        delete cAll;
        delete leg;
    }
    else
    {
        TH2F *h_eff = m_efficiencyMaps[histName];
        if (!h_eff)
        {
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
        if (!(histName == "reco_zT_vs_jetPt") && !(histName == "mc_genYield_mcPt_mczT"))
            h_eff->GetZaxis()->SetRangeUser(0.0, 1.0);

        // Draw efficiency map
        h_eff->Draw("COLZ");

        TLatex *latex = new TLatex();
        latex->SetNDC();
        latex->SetTextSize(0.04);
        TString lhcb_label = TString::Format("#font[12]{LHCb} work-in-progress");
        latex->DrawLatex(0.15, 0.85, lhcb_label);
        TString jetRadiusText = TString::Format("Anti #it{k}_{T}, R = %.1f", m_jetRadius);
        latex->DrawLatex(0.15, 0.80, jetRadiusText);
        // Add text with overall efficiency
        TH2F *h_num = m_efficiencyMaps["jet_numerator"];
        TH2F *h_den = m_efficiencyMaps["jet_denominator"];
        if (h_num && h_den && !(histName == "reco_zT_vs_jetPt"))
        {
            double overallEff = (double)h_num->GetEntries() / h_den->GetEntries();
            TString effText = TString::Format("Overall Jet Finding Efficiency: %.3f", overallEff);
            latex->DrawLatex(0.15, 0.75, effText);
        }

        // Add label for zT plot
        if (histName == "jet_finding_efficiency_zT")
        {
            latex->DrawLatex(0.15, 0.92, "Jet Finding Efficiency vs z_{T}");
        }

        // Plot reco zT vs jet pT
        if (histName == "reco_zT_vs_jetPt")
        {
            h_eff->SetTitle("Reco z_{T} vs Jet p_{T}");
            h_eff->GetXaxis()->SetTitle("z_{T}");
            h_eff->GetYaxis()->SetTitle("Jet p_{T} [GeV]");
        }

        // Save canvas
        std::string outDir = "jetFindingEfficiencies";
        ensureDirectoryExists(outDir);
        TString plotName = TString::Format("%s/%s.png", outDir.c_str(), histName.c_str());
        c1->SaveAs(plotName);
        std::cout << "Saved plot: " << plotName << std::endl;
        // Also save as PDF
        plotName = TString::Format("%s/%s.pdf", outDir.c_str(), histName.c_str());
        c1->SaveAs(plotName);
        delete c1;
    }
}

// Main function for the standalone jet finding efficiency calculation
int JetFindingEfficiencyStandaloneRun(TString inputFile = "/media/niviths/SSD2/lhcb_analysis_SSD/20250728_pPb_MC_output/20250728_pPb_MC_output.root",
                                      TString outputFile = "output_jet_finding_efficiency.root",
                                      double jetRadius = 0.4, double minJetPt = 5.0, double maxJetPt = 60.0,
                                      double minJetEta = 2.5, double maxJetEta = 4.0,
                                      double minD0Pt = 1.0, double maxD0Pt = 50.0,
                                      double minD0Eta = 2.0, double maxD0Eta = 4.5,
                                      bool makePlots = true)
{
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
    std::cout << "  - Numerator: Reconstructed jets with reconstructed D0 and MC match" << std::endl;
    std::cout << "  - Denominator: MC jets with generator level D0 in jet radius" << std::endl;
    std::cout << std::endl;

    try
    {
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

        // std::vector<double> customJetEtaBins;
        // for (double eta = 2.0; eta <= 4.5; eta += 0.25)
        //     customJetEtaBins.push_back(eta);

        // std::vector<double> customD0PtBins;
        // for (double pt = 2.0; pt < 10.0; pt += 1.0)
        //     customD0PtBins.push_back(pt);

        // calculator.SetJetPtBins(customJetPtBins);
        // calculator.SetJetEtaBins(customJetEtaBins);
        // calculator.SetD0PtBins(customD0PtBins);

        // Initialize
        if (!calculator.Initialize())
        {
            std::cerr << "Error: Failed to initialize calculator" << std::endl;
            return 1;
        }

        // Calculate efficiency
        std::cout << "Starting jet finding efficiency calculation..." << std::endl;
        calculator.CalculateJetFindingEfficiency();

        // Save results
        calculator.SaveResults();

        // Generate plots if requested
        if (makePlots)
        {
            std::cout << "\nGenerating efficiency plots..." << std::endl;
            calculator.PlotEfficiency("jet_finding_efficiency");
            calculator.PlotEfficiency("jet_finding_efficiency_d0pt");
            calculator.PlotEfficiency("jet_finding_efficiency_zT");
            calculator.PlotEfficiency("reco_zT_vs_jetPt");
            calculator.PlotEfficiency("mc_genYield_mcPt_mczT");
            calculator.PlotEfficiency("jet_finding_efficiency_zT_allProjections");

            // Individual jet pt slice plots (project numerator and denominator, then divide)
            std::vector<double> startPt = {5, 10, 15, 20, 30};
            std::vector<double> endPt = {10, 15, 20, 30, 50};
            std::string outDir = "jetFindingEfficiencies";
            ensureDirectoryExists(outDir);
            TH2F *h_num_zT = calculator.GetEfficiencyHistogram("jet_numerator_zT");
            TH2F *h_den_zT = calculator.GetEfficiencyHistogram("jet_denominator_zT");
            if (h_num_zT && h_den_zT)
            {
                for (size_t i = 0; i < startPt.size(); ++i)
                {
                    double ptMin = startPt[i];
                    double ptMax = endPt[i];
                    int binMin = h_num_zT->GetYaxis()->FindBin(ptMin);
                    int binMax = h_num_zT->GetYaxis()->FindBin(ptMax) - 1;
                    TString sliceName = TString::Format("jet_finding_efficiency_zT_pt_%g_%g", ptMin, ptMax);
                    TH1D *num_proj = h_num_zT->ProjectionX(sliceName + "_num", binMin, binMax);
                    TH1D *den_proj = h_den_zT->ProjectionX(sliceName + "_den", binMin, binMax);
                    TH1D *eff_proj = (TH1D *)num_proj->Clone(sliceName);
                    eff_proj->SetTitle(TString::Format("Jet Finding Efficiency vs z_{T} (%g < p_{T} < %g GeV)", ptMin, ptMax));
                    eff_proj->SetStats(0);
                    eff_proj->GetYaxis()->SetTitle("Efficiency");
                    eff_proj->GetXaxis()->SetTitle("z_{T}");
                    eff_proj->SetMinimum(0.0);
                    eff_proj->SetMaximum(1.0);
                    eff_proj->Divide(eff_proj, den_proj, 1.0, 1.0, "B");
                    TCanvas *cSlice = new TCanvas(sliceName, sliceName, 800, 600);
                    eff_proj->Draw("E1");
                    TString pngName = TString::Format("%s/%s.png", outDir.c_str(), sliceName.Data());
                    TString pdfName = TString::Format("%s/%s.pdf", outDir.c_str(), sliceName.Data());
                    cSlice->SaveAs(pngName);
                    cSlice->SaveAs(pdfName);
                    std::cout << "Saved plot: " << pngName << std::endl;
                    delete cSlice;
                    delete num_proj;
                    delete den_proj;
                    delete eff_proj;
                }
            }
        }

        std::cout << "\n=== Jet Finding Efficiency calculation completed successfully! ===" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
