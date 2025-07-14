#include "D0Efficiency.h"

D0Efficiency::D0Efficiency(const std::string& inputFileName, 
                          const std::string& outputFileName,
                          const std::string& efficiencyType,
                          bool isMC) 
    : m_inputFile(nullptr), m_outputFile(nullptr), m_tree(nullptr),
      m_isMC(isMC), m_efficiencyType(efficiencyType),
      m_d0MassWindow(50.0), m_minPt(2.0), m_minEta(2.0), m_maxEta(5.0),
      m_kaonPIDCut(0.5), m_pionPIDCut(0.5), m_ghostProbCut(0.3), m_trackChi2Cut(3.0),
      m_minDaughterMomentum(2.0) {
    
    // Open input file
    m_inputFile = TFile::Open(inputFileName.c_str(), "READ");
    if (!m_inputFile || m_inputFile->IsZombie()) {
        std::cerr << "Error: Could not open input file " << inputFileName << std::endl;
        return;
    }
    
    // Create output file
    m_outputFile = new TFile(outputFileName.c_str(), "RECREATE");
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
    d0_pt = nullptr;
    d0_eta = nullptr;
    d0_phi = nullptr;
    d0_mass = nullptr;
    d0_px = nullptr;
    d0_py = nullptr;
    d0_pz = nullptr;
    d0_e = nullptr;
    
    dau_pid = nullptr;
    dau_pt = nullptr;
    dau_eta = nullptr;
    dau_phi = nullptr;
    dau_px = nullptr;
    dau_py = nullptr;
    dau_pz = nullptr;
    dau_e = nullptr;
    dau_d0_idx = nullptr;
    dau_pnn_k = nullptr;
    dau_pnn_pi = nullptr;
    dau_prb_ghost = nullptr;
    dau_trckChi2 = nullptr;
    
    // MC truth variables (only initialize if MC)
    if (m_isMC) {
        mc_d0_pid = nullptr;
        mc_d0_pt = nullptr;
        mc_d0_eta = nullptr;
        mc_d0_phi = nullptr;
        mc_d0_mass = nullptr;
        mc_d0_px = nullptr;
        mc_d0_py = nullptr;
        mc_d0_pz = nullptr;
        mc_d0_e = nullptr;
        mc_d0_origin = nullptr;
        mc_d0_matched = nullptr;
        mc_d0_matched_quality = nullptr;
        
        mc_dau_pid = nullptr;
        mc_dau_pt = nullptr;
        mc_dau_eta = nullptr;
        mc_dau_phi = nullptr;
        mc_dau_d0_idx = nullptr;
    }
    
    // Set default binning
    m_ptBins = {2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0};
    m_etaBins = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
    m_pBins = {2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0};
    
    std::cout << "D0Efficiency initialized for " << efficiencyType << " efficiency calculation" << std::endl;
    std::cout << "Input file: " << inputFileName << std::endl;
    std::cout << "Output file: " << outputFileName << std::endl;
    std::cout << "MC mode: " << (m_isMC ? "Yes" : "No") << std::endl;
}

D0Efficiency::~D0Efficiency() {
    if (m_inputFile) {
        m_inputFile->Close();
        delete m_inputFile;
    }
    if (m_outputFile) {
        m_outputFile->Close();
        delete m_outputFile;
    }
    
    // Clean up histograms
    for (auto& pair : m_efficiencyMaps) {
        delete pair.second;
    }
    for (auto& pair : m_efficiencyProjections) {
        delete pair.second;
    }
    for (auto& pair : m_efficiencyFunctions) {
        delete pair.second;
    }
    for (auto& pair : m_efficiencyObjects) {
        delete pair.second;
    }
}

bool D0Efficiency::Initialize() {
    if (!m_tree) {
        std::cerr << "Error: Tree not available" << std::endl;
        return false;
    }
    
    SetupBranches();
    CreateHistograms();
    
    std::cout << "D0Efficiency initialized successfully" << std::endl;
    std::cout << "Tree entries: " << m_tree->GetEntries() << std::endl;
    
    return true;
}

void D0Efficiency::SetBinning(const std::vector<double>& ptBins, const std::vector<double>& etaBins) {
    m_ptBins = ptBins;
    m_etaBins = etaBins;
    
    std::cout << "Binning updated:" << std::endl;
    std::cout << "pT bins: ";
    for (size_t i = 0; i < m_ptBins.size(); ++i) {
        std::cout << m_ptBins[i];
        if (i < m_ptBins.size() - 1) std::cout << ", ";
    }
    std::cout << " GeV" << std::endl;
    
    std::cout << "η bins: ";
    for (size_t i = 0; i < m_etaBins.size(); ++i) {
        std::cout << m_etaBins[i];
        if (i < m_etaBins.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;
}

void D0Efficiency::SetBinning(const std::vector<double>& ptBins, const std::vector<double>& etaBins, const std::vector<double>& pBins) {
    m_ptBins = ptBins;
    m_etaBins = etaBins;
    m_pBins = pBins;
    
    std::cout << "Binning updated:" << std::endl;
    std::cout << "pT bins: ";
    for (size_t i = 0; i < m_ptBins.size(); ++i) {
        std::cout << m_ptBins[i];
        if (i < m_ptBins.size() - 1) std::cout << ", ";
    }
    std::cout << " GeV" << std::endl;
    
    std::cout << "η bins: ";
    for (size_t i = 0; i < m_etaBins.size(); ++i) {
        std::cout << m_etaBins[i];
        if (i < m_etaBins.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;
    
    std::cout << "p bins: ";
    for (size_t i = 0; i < m_pBins.size(); ++i) {
        std::cout << m_pBins[i];
        if (i < m_pBins.size() - 1) std::cout << ", ";
    }
    std::cout << " GeV" << std::endl;
}

void D0Efficiency::SetPIDCuts(double kaonPIDCut, double pionPIDCut, 
                              double ghostProbCut, double trackChi2Cut) {
    m_kaonPIDCut = kaonPIDCut;
    m_pionPIDCut = pionPIDCut;
    m_ghostProbCut = ghostProbCut;
    m_trackChi2Cut = trackChi2Cut;
    
    std::cout << "PID cuts updated:" << std::endl;
    std::cout << "Kaon PID cut: " << m_kaonPIDCut << std::endl;
    std::cout << "Pion PID cut: " << m_pionPIDCut << std::endl;
    std::cout << "Ghost prob cut: " << m_ghostProbCut << std::endl;
    std::cout << "Track chi2 cut: " << m_trackChi2Cut << std::endl;
}

void D0Efficiency::SetD0Selection(double massWindow, double minPt, double minEta, double maxEta) {
    m_d0MassWindow = massWindow;
    m_minPt = minPt;
    m_minEta = minEta;
    m_maxEta = maxEta;
    
    std::cout << "D0 selection updated:" << std::endl;
    std::cout << "Mass window: ±" << m_d0MassWindow << " MeV" << std::endl;
    std::cout << "Minimum pT: " << m_minPt << " GeV" << std::endl;
    std::cout << "Minimum |η|: " << m_minEta << std::endl;
    std::cout << "Maximum |η|: " << m_maxEta << std::endl;
}

void D0Efficiency::SetDaughterCuts(double minMomentum) {
    m_minDaughterMomentum = minMomentum;
    
    std::cout << "Daughter cuts updated:" << std::endl;
    std::cout << "Minimum momentum: " << m_minDaughterMomentum << " GeV" << std::endl;
    std::cout << "Daughter η range: " << m_minEta << " < |η| < " << m_maxEta << std::endl;
}

void D0Efficiency::CalculateKaonPIDEfficiency() {
    if (!m_tree) {
        std::cerr << "Error: Tree not available" << std::endl;
        return;
    }
    
    std::cout << "Calculating kaon PID efficiency..." << std::endl;
    
    // Get efficiency histograms
    TH2F* h_num = m_efficiencyMaps["kaon_PID_numerator"];
    TH2F* h_den = m_efficiencyMaps["kaon_PID_denominator"];
    TEfficiency* eff_obj = m_efficiencyObjects["kaon_PID_efficiency"];
    
    // Get momentum-based histograms (optional)
    TH2F* h_num_p = nullptr;
    TH2F* h_den_p = nullptr;
    TEfficiency* eff_obj_p = nullptr;
    
    if (!m_pBins.empty()) {
        h_num_p = m_efficiencyMaps["kaon_PID_numerator_p"];
        h_den_p = m_efficiencyMaps["kaon_PID_denominator_p"];
        eff_obj_p = m_efficiencyObjects["kaon_PID_efficiency_p"];
    }
    
    if (!h_num || !h_den || !eff_obj) {
        std::cerr << "Error: Kaon efficiency histograms/objects not available" << std::endl;
        return;
    }
    
    // Loop over all events
    Long64_t nEntries = m_tree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        m_tree->GetEntry(entry);
        
        if (entry % 10000 == 0) {
            std::cout << "Processing entry " << entry << "/" << nEntries << std::endl;
        }
        
        // Loop over all D0 candidates in this event
        for (int d0_idx = 0; d0_idx < n_d0s; ++d0_idx) {
            // Check if D0 passes basic selection
            if (!PassesD0Selection(d0_idx)) continue;
            
            // Check if daughters pass quality cuts
            if (!PassesQualityCuts(d0_idx)) continue;
            
            double pt = d0_pt->at(d0_idx);
            double eta = d0_eta->at(d0_idx);
            double p = CalculateMomentum(d0_idx);
            
            // Fill denominator (all D0s passing basic selection)
            h_den->Fill(pt, eta);
            if (h_den_p) h_den_p->Fill(p, eta);
            
            // Check if kaon passes PID cuts for numerator
            bool passesKaonPID = PassesKaonPIDSelection(d0_idx);
            
            // Fill TEfficiency object
            eff_obj->Fill(passesKaonPID, pt, eta);
            if (eff_obj_p) eff_obj_p->Fill(passesKaonPID, p, eta);
            
            // Fill numerator (D0s with kaon passing PID cuts)
            if (passesKaonPID) {
                h_num->Fill(pt, eta);
                if (h_num_p) h_num_p->Fill(p, eta);
            }
        }
    }
    
    // Create efficiency map (pT vs eta)
    TH2F* h_eff = (TH2F*)h_num->Clone("kaon_PID_efficiency");
    h_eff->Divide(h_den);
    h_eff->SetTitle("Kaon PID Efficiency;p_{T} [GeV];#eta");
    m_efficiencyMaps["kaon_PID_efficiency"] = h_eff;
    
    // Create efficiency map (p vs eta) if available
    if (h_num_p && h_den_p) {
        TH2F* h_eff_p = (TH2F*)h_num_p->Clone("kaon_PID_efficiency_p");
        h_eff_p->Divide(h_den_p);
        h_eff_p->SetTitle("Kaon PID Efficiency;p [GeV];#eta");
        m_efficiencyMaps["kaon_PID_efficiency_p"] = h_eff_p;
    }
    
    std::cout << "Kaon PID efficiency calculation completed" << std::endl;
    
    // Print some statistics
    int totalDen = h_den->GetEntries();
    int totalNum = h_num->GetEntries();
    double overallEff = (totalDen > 0) ? (double)totalNum / totalDen : 0.0;
    std::cout << "Overall kaon PID efficiency: " << overallEff << " (" << totalNum << "/" << totalDen << ")" << std::endl;
}

void D0Efficiency::CalculatePionPIDEfficiency() {
    if (!m_tree) {
        std::cerr << "Error: Tree not available" << std::endl;
        return;
    }
    
    std::cout << "Calculating pion PID efficiency..." << std::endl;
    
    // Get efficiency histograms
    TH2F* h_num = m_efficiencyMaps["pion_PID_numerator"];
    TH2F* h_den = m_efficiencyMaps["pion_PID_denominator"];
    TEfficiency* eff_obj = m_efficiencyObjects["pion_PID_efficiency"];
    
    // Get momentum-based histograms (optional)
    TH2F* h_num_p = nullptr;
    TH2F* h_den_p = nullptr;
    TEfficiency* eff_obj_p = nullptr;
    
    if (!m_pBins.empty()) {
        h_num_p = m_efficiencyMaps["pion_PID_numerator_p"];
        h_den_p = m_efficiencyMaps["pion_PID_denominator_p"];
        eff_obj_p = m_efficiencyObjects["pion_PID_efficiency_p"];
    }
    
    if (!h_num || !h_den || !eff_obj) {
        std::cerr << "Error: Pion efficiency histograms/objects not available" << std::endl;
        return;
    }
    
    // Loop over all events
    Long64_t nEntries = m_tree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        m_tree->GetEntry(entry);
        
        if (entry % 10000 == 0) {
            std::cout << "Processing entry " << entry << "/" << nEntries << std::endl;
        }
        
        // Loop over all D0 candidates in this event
        for (int d0_idx = 0; d0_idx < n_d0s; ++d0_idx) {
            // Check if D0 passes basic selection
            if (!PassesD0Selection(d0_idx)) continue;
            
            // Check if daughters pass quality cuts
            if (!PassesQualityCuts(d0_idx)) continue;
            
            double pt = d0_pt->at(d0_idx);
            double eta = d0_eta->at(d0_idx);
            double p = CalculateMomentum(d0_idx);
            
            // Fill denominator (all D0s passing basic selection)
            h_den->Fill(pt, eta);
            if (h_den_p) h_den_p->Fill(p, eta);
            
            // Check if pion passes PID cuts for numerator
            bool passesPionPID = PassesPionPIDSelection(d0_idx);
            
            // Fill TEfficiency object
            eff_obj->Fill(passesPionPID, pt, eta);
            if (eff_obj_p) eff_obj_p->Fill(passesPionPID, p, eta);
            
            // Fill numerator (D0s with pion passing PID cuts)
            if (passesPionPID) {
                h_num->Fill(pt, eta);
                if (h_num_p) h_num_p->Fill(p, eta);
            }
        }
    }
    
    // Create efficiency map (pT vs eta)
    TH2F* h_eff = (TH2F*)h_num->Clone("pion_PID_efficiency");
    h_eff->Divide(h_den);
    h_eff->SetTitle("Pion PID Efficiency;p_{T} [GeV];#eta");
    m_efficiencyMaps["pion_PID_efficiency"] = h_eff;
    
    // Create efficiency map (p vs eta) if available
    if (h_num_p && h_den_p) {
        TH2F* h_eff_p = (TH2F*)h_num_p->Clone("pion_PID_efficiency_p");
        h_eff_p->Divide(h_den_p);
        h_eff_p->SetTitle("Pion PID Efficiency;p [GeV];#eta");
        m_efficiencyMaps["pion_PID_efficiency_p"] = h_eff_p;
    }
    
    std::cout << "Pion PID efficiency calculation completed" << std::endl;
    
    // Print some statistics
    int totalDen = h_den->GetEntries();
    int totalNum = h_num->GetEntries();
    double overallEff = (totalDen > 0) ? (double)totalNum / totalDen : 0.0;
    std::cout << "Overall pion PID efficiency: " << overallEff << " (" << totalNum << "/" << totalDen << ")" << std::endl;
}

void D0Efficiency::CalculateCombinedPIDEfficiency() {
    if (!m_tree) {
        std::cerr << "Error: Tree not available" << std::endl;
        return;
    }
    
    std::cout << "Calculating combined D0 PID efficiency..." << std::endl;
    
    // Get efficiency histograms
    TH2F* h_num = m_efficiencyMaps["combined_PID_numerator"];
    TH2F* h_den = m_efficiencyMaps["combined_PID_denominator"];
    TEfficiency* eff_obj = m_efficiencyObjects["combined_PID_efficiency"];
    
    // Get momentum-based histograms (optional)
    TH2F* h_num_p = nullptr;
    TH2F* h_den_p = nullptr;
    TEfficiency* eff_obj_p = nullptr;
    
    if (!m_pBins.empty()) {
        h_num_p = m_efficiencyMaps["combined_PID_numerator_p"];
        h_den_p = m_efficiencyMaps["combined_PID_denominator_p"];
        eff_obj_p = m_efficiencyObjects["combined_PID_efficiency_p"];
    }
    
    if (!h_num || !h_den || !eff_obj) {
        std::cerr << "Error: Combined efficiency histograms/objects not available" << std::endl;
        return;
    }
    
    // Loop over all events
    Long64_t nEntries = m_tree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        m_tree->GetEntry(entry);
        
        if (entry % 10000 == 0) {
            std::cout << "Processing entry " << entry << "/" << nEntries << std::endl;
        }
        
        // Loop over all D0 candidates in this event
        for (int d0_idx = 0; d0_idx < n_d0s; ++d0_idx) {
            // Check if D0 passes basic selection
            if (!PassesD0Selection(d0_idx)) continue;
            
            // Check if daughters pass quality cuts
            if (!PassesQualityCuts(d0_idx)) continue;
            
            double pt = d0_pt->at(d0_idx);
            double eta = d0_eta->at(d0_idx);
            double p = CalculateMomentum(d0_idx);
            
            // Fill denominator (all D0s passing basic selection)
            h_den->Fill(pt, eta);
            if (h_den_p) h_den_p->Fill(p, eta);
            
            // Check if BOTH daughters pass PID cuts for numerator
            bool passesCombinedPID = PassesPIDSelection(d0_idx);
            
            // Fill TEfficiency object
            eff_obj->Fill(passesCombinedPID, pt, eta);
            if (eff_obj_p) eff_obj_p->Fill(passesCombinedPID, p, eta);
            
            // Fill numerator (D0s with both daughters passing PID cuts)
            if (passesCombinedPID) {
                h_num->Fill(pt, eta);
                if (h_num_p) h_num_p->Fill(p, eta);
            }
        }
    }
    
    // Create efficiency map (pT vs eta)
    TH2F* h_eff = (TH2F*)h_num->Clone("combined_PID_efficiency");
    h_eff->Divide(h_den);
    h_eff->SetTitle("Combined D0 PID Efficiency;p_{T} [GeV];#eta");
    m_efficiencyMaps["combined_PID_efficiency"] = h_eff;
    
    // Create efficiency map (p vs eta) if available
    if (h_num_p && h_den_p) {
        TH2F* h_eff_p = (TH2F*)h_num_p->Clone("combined_PID_efficiency_p");
        h_eff_p->Divide(h_den_p);
        h_eff_p->SetTitle("Combined D0 PID Efficiency;p [GeV];#eta");
        m_efficiencyMaps["combined_PID_efficiency_p"] = h_eff_p;
    }
    
    std::cout << "Combined D0 PID efficiency calculation completed" << std::endl;
    
    // Print some statistics
    int totalDen = h_den->GetEntries();
    int totalNum = h_num->GetEntries();
    double overallEff = (totalDen > 0) ? (double)totalNum / totalDen : 0.0;
    std::cout << "Overall combined D0 PID efficiency: " << overallEff << " (" << totalNum << "/" << totalDen << ")" << std::endl;
}

void D0Efficiency::CalculateRecoEfficiency() {
    if (!m_isMC) {
        std::cerr << "Error: Reconstruction efficiency can only be calculated for MC" << std::endl;
        return;
    }
    
    std::cout << "Calculating reconstruction efficiency..." << std::endl;
    
    // Get efficiency histograms
    TH2F* h_num = m_efficiencyMaps["reco_numerator"];
    TH2F* h_den = m_efficiencyMaps["reco_denominator"];
    TEfficiency* eff_obj = m_efficiencyObjects["reco_efficiency"];
    
    // Get momentum-based histograms (optional)
    TH2F* h_num_p = nullptr;
    TH2F* h_den_p = nullptr;
    TEfficiency* eff_obj_p = nullptr;
    
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
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        m_tree->GetEntry(entry);
        
        if (entry % 10000 == 0) {
            std::cout << "Processing entry " << entry << "/" << nEntries << std::endl;
        }
        
        // Loop over all MC D0 particles
        for (size_t mc_idx = 0; mc_idx < mc_d0_pt->size(); ++mc_idx) {
            // Check if MC D0 is in acceptance
            double mc_pt = mc_d0_pt->at(mc_idx);
            double mc_eta = mc_d0_eta->at(mc_idx);
            double mc_p = CalculateMCMomentum(mc_idx);
            
            if (mc_pt < m_minPt || mc_eta > m_maxEta || mc_eta < m_minEta) continue;
            
            // Fill denominator (all MC D0s in acceptance)
            h_den->Fill(mc_pt, mc_eta);
            if (h_den_p) h_den_p->Fill(mc_p, mc_eta);
            
            // Check if this MC D0 has a reconstructed match
            int matchedRecoIdx = mc_d0_matched->at(mc_idx);
            bool passesReco = false;
            
            if (matchedRecoIdx >= 0 && matchedRecoIdx < n_d0s) {
                // Check if reconstructed D0 passes selection
                if (PassesD0Selection(matchedRecoIdx)) {
                    passesReco = true;
                    h_num->Fill(mc_pt, mc_eta);
                    if (h_num_p) h_num_p->Fill(mc_p, mc_eta);
                }
            }
            
            // Fill TEfficiency object
            eff_obj->Fill(passesReco, mc_pt, mc_eta);
            if (eff_obj_p) eff_obj_p->Fill(passesReco, mc_p, mc_eta);
        }
    }
    
    // Create efficiency map (pT vs eta)
    TH2F* h_eff = (TH2F*)h_num->Clone("reco_efficiency");
    h_eff->Divide(h_den);
    h_eff->SetTitle("Reconstruction Efficiency;p_{T} [GeV];#eta");
    m_efficiencyMaps["reco_efficiency"] = h_eff;
    
    // Create efficiency map (p vs eta) if available
    if (h_num_p && h_den_p) {
        TH2F* h_eff_p = (TH2F*)h_num_p->Clone("reco_efficiency_p");
        h_eff_p->Divide(h_den_p);
        h_eff_p->SetTitle("Reconstruction Efficiency;p [GeV];#eta");
        m_efficiencyMaps["reco_efficiency_p"] = h_eff_p;
    }
    
    std::cout << "Reconstruction efficiency calculation completed" << std::endl;
    
    // Print some statistics
    int totalDen = h_den->GetEntries();
    int totalNum = h_num->GetEntries();
    double overallEff = (totalDen > 0) ? (double)totalNum / totalDen : 0.0;
    std::cout << "Overall reconstruction efficiency: " << overallEff << " (" << totalNum << "/" << totalDen << ")" << std::endl;
}

void D0Efficiency::ParameterizeEfficiency(const std::string& efficiencyName, 
                                         const std::string& functionType) {
    TH2F* h_eff = m_efficiencyMaps[efficiencyName];
    if (!h_eff) {
        std::cerr << "Error: Efficiency histogram '" << efficiencyName << "' not found" << std::endl;
        return;
    }
    
    std::cout << "Parameterizing " << efficiencyName << " efficiency..." << std::endl;
    
    // Create functions for each η bin
    for (int etaBin = 1; etaBin <= h_eff->GetNbinsY(); ++etaBin) {
        // Get projection for this η bin
        TH1D* h_proj = h_eff->ProjectionX(Form("%s_proj_eta%d", efficiencyName.c_str(), etaBin), 
                                          etaBin, etaBin);
        
        if (h_proj->GetEntries() == 0) continue;
        
        // Define function based on type
        TF1* func = nullptr;
        if (functionType == "exp") {
            func = new TF1(Form("%s_func_eta%d", efficiencyName.c_str(), etaBin),
                          "[0] - [1]*TMath::Exp(-[2]*x)", m_ptBins[0], m_ptBins.back());
            func->SetParameters(0.8, 0.3, 0.5);
        } else if (functionType == "pol2") {
            func = new TF1(Form("%s_func_eta%d", efficiencyName.c_str(), etaBin),
                          "[0] + [1]*x + [2]*x*x", m_ptBins[0], m_ptBins.back());
            func->SetParameters(0.5, 0.1, -0.01);
        } else if (functionType == "erf") {
            func = new TF1(Form("%s_func_eta%d", efficiencyName.c_str(), etaBin),
                          "[0]*TMath::Erf([1]*(x-[2]))", m_ptBins[0], m_ptBins.back());
            func->SetParameters(0.8, 0.5, 3.0);
        } else {
            std::cerr << "Error: Unknown function type '" << functionType << "'" << std::endl;
            delete h_proj;
            continue;
        }
        
        // Fit the projection
        h_proj->Fit(func, "QN");
        
        // Store function
        m_efficiencyFunctions[Form("%s_eta%d", efficiencyName.c_str(), etaBin)] = func;
        
        delete h_proj;
    }
    
    std::cout << "Parameterization completed for " << efficiencyName << std::endl;
}

double D0Efficiency::GetEfficiencyWeight(double pt, double eta, const std::string& efficiencyName, 
                                        bool useFunction) {
    if (useFunction) {
        // Find the appropriate η bin
        int etaBin = -1;
        for (size_t i = 0; i < m_etaBins.size() - 1; ++i) {
            if (eta >= m_etaBins[i] && eta < m_etaBins[i + 1]) {
                etaBin = i + 1;
                break;
            }
        }
        
        if (etaBin < 0) return 1.0;
        
        // Get the function for this η bin
        TF1* func = m_efficiencyFunctions[Form("%s_eta%d", efficiencyName.c_str(), etaBin)];
        if (!func) return 1.0;
        
        double efficiency = func->Eval(pt);
        return (efficiency > 0) ? 1.0 / efficiency : 0.0;
    } else {
        // Use histogram
        TH2F* h_eff = m_efficiencyMaps[efficiencyName];
        if (!h_eff) return 1.0;
        
        int ptBin = h_eff->GetXaxis()->FindBin(pt);
        int etaBin = h_eff->GetYaxis()->FindBin(eta);
        
        double efficiency = h_eff->GetBinContent(ptBin, etaBin);
        return (efficiency > 0) ? 1.0 / efficiency : 0.0;
    }
}

TEfficiency* D0Efficiency::GetEfficiencyObject(const std::string& efficiencyName) {
    auto it = m_efficiencyObjects.find(efficiencyName);
    if (it != m_efficiencyObjects.end()) {
        return it->second;
    }
    return nullptr;
}

void D0Efficiency::PerformClosureTest() {
    if (!m_isMC) {
        std::cerr << "Error: Closure test can only be performed on MC" << std::endl;
        return;
    }
    
    std::cout << "Performing closure test..." << std::endl;
    
    // Create histograms for closure test
    TH2F* h_weighted = new TH2F("closure_weighted", "Weighted Distribution;p_{T} [GeV];#eta",
                               m_ptBins.size()-1, &m_ptBins[0], m_etaBins.size()-1, &m_etaBins[0]);
    TH2F* h_truth = new TH2F("closure_truth", "Truth Distribution;p_{T} [GeV];#eta",
                            m_ptBins.size()-1, &m_ptBins[0], m_etaBins.size()-1, &m_etaBins[0]);
    
    // Loop over all events
    Long64_t nEntries = m_tree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        m_tree->GetEntry(entry);
        
        // Fill truth distribution
        for (size_t mc_idx = 0; mc_idx < mc_d0_pt->size(); ++mc_idx) {
            double mc_pt = mc_d0_pt->at(mc_idx);
            double mc_eta = mc_d0_eta->at(mc_idx);

            // if (mc_pt < m_minPt || TMath::Abs(mc_eta) > m_maxEta || TMath::Abs(mc_eta) < m_minEta) continue;
            if (mc_pt < m_minPt || mc_eta > m_maxEta || mc_eta < m_minEta) continue;

            h_truth->Fill(mc_pt, mc_eta);
        }
        
        // Fill weighted distribution
        for (int d0_idx = 0; d0_idx < n_d0s; ++d0_idx) {
            if (!PassesD0Selection(d0_idx)) continue;
            
            double pt = d0_pt->at(d0_idx);
            double eta = d0_eta->at(d0_idx);
            
            // Get efficiency weight
            double weight = GetEfficiencyWeight(pt, eta, m_efficiencyType + "_efficiency");
            
            h_weighted->Fill(pt, eta, weight);
        }
    }
    
    // Create ratio histogram
    TH2F* h_ratio = (TH2F*)h_weighted->Clone("closure_ratio");
    h_ratio->Divide(h_truth);
    h_ratio->SetTitle("Closure Test Ratio (Weighted/Truth);p_{T} [GeV];#eta");
    
    // Store results
    m_efficiencyMaps["closure_weighted"] = h_weighted;
    m_efficiencyMaps["closure_truth"] = h_truth;
    m_efficiencyMaps["closure_ratio"] = h_ratio;
    
    std::cout << "Closure test completed" << std::endl;
    
    // Print statistics
    double meanRatio = h_ratio->GetMean();
    double rmsRatio = h_ratio->GetRMS();
    std::cout << "Closure test ratio: " << meanRatio << " ± " << rmsRatio << std::endl;
}

void D0Efficiency::SaveResults() {
    if (!m_outputFile) {
        std::cerr << "Error: Output file not available" << std::endl;
        return;
    }
    
    m_outputFile->cd();
    
    // Save all efficiency histograms
    for (auto& pair : m_efficiencyMaps) {
        pair.second->Write();
    }
    
    // Save all TEfficiency objects
    for (auto& pair : m_efficiencyObjects) {
        pair.second->Write();
    }
    
    // Save all efficiency functions
    for (auto& pair : m_efficiencyFunctions) {
        pair.second->Write();
    }
    
    // Save projections
    for (auto& pair : m_efficiencyProjections) {
        pair.second->Write();
    }
    
    m_outputFile->Write();
    
    std::cout << "Results saved to " << m_outputFile->GetName() << std::endl;
}

void D0Efficiency::PrintSummary() {
    std::cout << "\n=== D0 Efficiency Summary ===" << std::endl;
    std::cout << "Efficiency type: " << m_efficiencyType << std::endl;
    std::cout << "MC mode: " << (m_isMC ? "Yes" : "No") << std::endl;
    std::cout << "Tree entries processed: " << m_tree->GetEntries() << std::endl;
    
    // Print efficiency statistics
    for (auto& pair : m_efficiencyMaps) {
        if (pair.first.find("efficiency") != std::string::npos) {
            TH2F* h = pair.second;
            std::cout << pair.first << ": " << h->GetMean() << " (mean)" << std::endl;
        }
    }
    
    std::cout << "==============================\n" << std::endl;
}

// Private methods implementation

void D0Efficiency::SetupBranches() {
    // Set up branches for reconstructed quantities
    m_tree->SetBranchAddress("n_d0s", &n_d0s);
    m_tree->SetBranchAddress("d0_pt", &d0_pt);
    m_tree->SetBranchAddress("d0_eta", &d0_eta);
    m_tree->SetBranchAddress("d0_phi", &d0_phi);
    m_tree->SetBranchAddress("d0_mass", &d0_mass);
    m_tree->SetBranchAddress("d0_px", &d0_px);
    m_tree->SetBranchAddress("d0_py", &d0_py);
    m_tree->SetBranchAddress("d0_pz", &d0_pz);
    m_tree->SetBranchAddress("d0_e", &d0_e);
    
    // Daughter branches
    m_tree->SetBranchAddress("n_daughters", &n_daughters);
    m_tree->SetBranchAddress("dau_pid", &dau_pid);
    m_tree->SetBranchAddress("dau_pt", &dau_pt);
    m_tree->SetBranchAddress("dau_eta", &dau_eta);
    m_tree->SetBranchAddress("dau_phi", &dau_phi);
    m_tree->SetBranchAddress("dau_px", &dau_px);
    m_tree->SetBranchAddress("dau_py", &dau_py);
    m_tree->SetBranchAddress("dau_pz", &dau_pz);
    m_tree->SetBranchAddress("dau_e", &dau_e);
    m_tree->SetBranchAddress("dau_d0_idx", &dau_d0_idx);
    m_tree->SetBranchAddress("dau_pnn_k", &dau_pnn_k);
    m_tree->SetBranchAddress("dau_pnn_pi", &dau_pnn_pi);
    m_tree->SetBranchAddress("dau_prb_ghost", &dau_prb_ghost);
    m_tree->SetBranchAddress("dau_trckChi2", &dau_trckChi2);
    
    // MC truth branches (only if MC)
    if (m_isMC) {
        m_tree->SetBranchAddress("mc_d0_pid", &mc_d0_pid);
        m_tree->SetBranchAddress("mc_d0_pt", &mc_d0_pt);
        m_tree->SetBranchAddress("mc_d0_eta", &mc_d0_eta);
        m_tree->SetBranchAddress("mc_d0_phi", &mc_d0_phi);
        m_tree->SetBranchAddress("mc_d0_mass", &mc_d0_mass);
        m_tree->SetBranchAddress("mc_d0_px", &mc_d0_px);
        m_tree->SetBranchAddress("mc_d0_py", &mc_d0_py);
        m_tree->SetBranchAddress("mc_d0_pz", &mc_d0_pz);
        m_tree->SetBranchAddress("mc_d0_e", &mc_d0_e);
        m_tree->SetBranchAddress("mc_d0_origin", &mc_d0_origin);
        m_tree->SetBranchAddress("mc_d0_matched", &mc_d0_matched);
        m_tree->SetBranchAddress("mc_d0_matched_quality", &mc_d0_matched_quality);
        
        m_tree->SetBranchAddress("mc_n_daughters", &mc_n_daughters);
        m_tree->SetBranchAddress("mc_dau_pid", &mc_dau_pid);
        m_tree->SetBranchAddress("mc_dau_pt", &mc_dau_pt);
        m_tree->SetBranchAddress("mc_dau_eta", &mc_dau_eta);
        m_tree->SetBranchAddress("mc_dau_phi", &mc_dau_phi);
        m_tree->SetBranchAddress("mc_dau_d0_idx", &mc_dau_d0_idx);
    }
}

void D0Efficiency::CreateHistograms() {
    // Create efficiency histograms based on efficiency type
    if (m_efficiencyType == "PID") {
        // Kaon PID efficiency histograms (pT vs eta)
        m_efficiencyMaps["kaon_PID_numerator"] = new TH2F("kaon_PID_numerator", "Kaon PID Numerator;p_{T} [GeV];#eta",
                                                         m_ptBins.size()-1, &m_ptBins[0], 
                                                         m_etaBins.size()-1, &m_etaBins[0]);
        m_efficiencyMaps["kaon_PID_denominator"] = new TH2F("kaon_PID_denominator", "Kaon PID Denominator;p_{T} [GeV];#eta",
                                                           m_ptBins.size()-1, &m_ptBins[0], 
                                                           m_etaBins.size()-1, &m_etaBins[0]);
        
        // Pion PID efficiency histograms (pT vs eta)
        m_efficiencyMaps["pion_PID_numerator"] = new TH2F("pion_PID_numerator", "Pion PID Numerator;p_{T} [GeV];#eta",
                                                         m_ptBins.size()-1, &m_ptBins[0], 
                                                         m_etaBins.size()-1, &m_etaBins[0]);
        m_efficiencyMaps["pion_PID_denominator"] = new TH2F("pion_PID_denominator", "Pion PID Denominator;p_{T} [GeV];#eta",
                                                           m_ptBins.size()-1, &m_ptBins[0], 
                                                           m_etaBins.size()-1, &m_etaBins[0]);
        
        // Combined D0 PID efficiency histograms (pT vs eta)
        m_efficiencyMaps["combined_PID_numerator"] = new TH2F("combined_PID_numerator", "Combined D0 PID Numerator;p_{T} [GeV];#eta",
                                                             m_ptBins.size()-1, &m_ptBins[0], 
                                                             m_etaBins.size()-1, &m_etaBins[0]);
        m_efficiencyMaps["combined_PID_denominator"] = new TH2F("combined_PID_denominator", "Combined D0 PID Denominator;p_{T} [GeV];#eta",
                                                               m_ptBins.size()-1, &m_ptBins[0], 
                                                               m_etaBins.size()-1, &m_etaBins[0]);
        
        // Create TEfficiency objects for each particle type (pT vs eta)
        m_efficiencyObjects["kaon_PID_efficiency"] = new TEfficiency("kaon_PID_efficiency", "Kaon PID Efficiency;p_{T} [GeV];#eta",
                                                                    m_ptBins.size()-1, &m_ptBins[0], 
                                                                    m_etaBins.size()-1, &m_etaBins[0]);
        m_efficiencyObjects["pion_PID_efficiency"] = new TEfficiency("pion_PID_efficiency", "Pion PID Efficiency;p_{T} [GeV];#eta",
                                                                    m_ptBins.size()-1, &m_ptBins[0], 
                                                                    m_etaBins.size()-1, &m_etaBins[0]);
        m_efficiencyObjects["combined_PID_efficiency"] = new TEfficiency("combined_PID_efficiency", "Combined D0 PID Efficiency;p_{T} [GeV];#eta",
                                                                        m_ptBins.size()-1, &m_ptBins[0], 
                                                                        m_etaBins.size()-1, &m_etaBins[0]);
        
        // Create momentum-based histograms if momentum binning is provided
        if (!m_pBins.empty()) {
            // Kaon PID efficiency histograms (p vs eta)
            m_efficiencyMaps["kaon_PID_numerator_p"] = new TH2F("kaon_PID_numerator_p", "Kaon PID Numerator;p [GeV];#eta",
                                                               m_pBins.size()-1, &m_pBins[0], 
                                                               m_etaBins.size()-1, &m_etaBins[0]);
            m_efficiencyMaps["kaon_PID_denominator_p"] = new TH2F("kaon_PID_denominator_p", "Kaon PID Denominator;p [GeV];#eta",
                                                                 m_pBins.size()-1, &m_pBins[0], 
                                                                 m_etaBins.size()-1, &m_etaBins[0]);
            
            // Pion PID efficiency histograms (p vs eta)
            m_efficiencyMaps["pion_PID_numerator_p"] = new TH2F("pion_PID_numerator_p", "Pion PID Numerator;p [GeV];#eta",
                                                               m_pBins.size()-1, &m_pBins[0], 
                                                               m_etaBins.size()-1, &m_etaBins[0]);
            m_efficiencyMaps["pion_PID_denominator_p"] = new TH2F("pion_PID_denominator_p", "Pion PID Denominator;p [GeV];#eta",
                                                                 m_pBins.size()-1, &m_pBins[0], 
                                                                 m_etaBins.size()-1, &m_etaBins[0]);
            
            // Combined D0 PID efficiency histograms (p vs eta)
            m_efficiencyMaps["combined_PID_numerator_p"] = new TH2F("combined_PID_numerator_p", "Combined D0 PID Numerator;p [GeV];#eta",
                                                                   m_pBins.size()-1, &m_pBins[0], 
                                                                   m_etaBins.size()-1, &m_etaBins[0]);
            m_efficiencyMaps["combined_PID_denominator_p"] = new TH2F("combined_PID_denominator_p", "Combined D0 PID Denominator;p [GeV];#eta",
                                                                     m_pBins.size()-1, &m_pBins[0], 
                                                                     m_etaBins.size()-1, &m_etaBins[0]);
            
            // Create TEfficiency objects for each particle type (p vs eta)
            m_efficiencyObjects["kaon_PID_efficiency_p"] = new TEfficiency("kaon_PID_efficiency_p", "Kaon PID Efficiency;p [GeV];#eta",
                                                                          m_pBins.size()-1, &m_pBins[0], 
                                                                          m_etaBins.size()-1, &m_etaBins[0]);
            m_efficiencyObjects["pion_PID_efficiency_p"] = new TEfficiency("pion_PID_efficiency_p", "Pion PID Efficiency;p [GeV];#eta",
                                                                          m_pBins.size()-1, &m_pBins[0], 
                                                                          m_etaBins.size()-1, &m_etaBins[0]);
            m_efficiencyObjects["combined_PID_efficiency_p"] = new TEfficiency("combined_PID_efficiency_p", "Combined D0 PID Efficiency;p [GeV];#eta",
                                                                              m_pBins.size()-1, &m_pBins[0], 
                                                                              m_etaBins.size()-1, &m_etaBins[0]);
        }
    } else if (m_efficiencyType == "reco" && m_isMC) {
        // Reconstruction efficiency histograms (pT vs eta)
        m_efficiencyMaps["reco_numerator"] = new TH2F("reco_numerator", "Reco Numerator;p_{T} [GeV];#eta",
                                                     m_ptBins.size()-1, &m_ptBins[0], 
                                                     m_etaBins.size()-1, &m_etaBins[0]);
        m_efficiencyMaps["reco_denominator"] = new TH2F("reco_denominator", "Reco Denominator;p_{T} [GeV];#eta",
                                                       m_ptBins.size()-1, &m_ptBins[0], 
                                                       m_etaBins.size()-1, &m_etaBins[0]);
        
        // Create TEfficiency object for reconstruction efficiency (pT vs eta)
        m_efficiencyObjects["reco_efficiency"] = new TEfficiency("reco_efficiency", "Reconstruction Efficiency;p_{T} [GeV];#eta",
                                                                m_ptBins.size()-1, &m_ptBins[0], 
                                                                m_etaBins.size()-1, &m_etaBins[0]);
        
        // Create momentum-based histograms if momentum binning is provided
        if (!m_pBins.empty()) {
            // Reconstruction efficiency histograms (p vs eta)
            m_efficiencyMaps["reco_numerator_p"] = new TH2F("reco_numerator_p", "Reco Numerator;p [GeV];#eta",
                                                           m_pBins.size()-1, &m_pBins[0], 
                                                           m_etaBins.size()-1, &m_etaBins[0]);
            m_efficiencyMaps["reco_denominator_p"] = new TH2F("reco_denominator_p", "Reco Denominator;p [GeV];#eta",
                                                             m_pBins.size()-1, &m_pBins[0], 
                                                             m_etaBins.size()-1, &m_etaBins[0]);
            
            // Create TEfficiency object for reconstruction efficiency (p vs eta)
            m_efficiencyObjects["reco_efficiency_p"] = new TEfficiency("reco_efficiency_p", "Reconstruction Efficiency;p [GeV];#eta",
                                                                      m_pBins.size()-1, &m_pBins[0], 
                                                                      m_etaBins.size()-1, &m_etaBins[0]);
        }
    }
    
    // Set histogram options
    for (auto& pair : m_efficiencyMaps) {
        pair.second->Sumw2();
    }
}

bool D0Efficiency::PassesD0Selection(int d0_idx) {
    if (d0_idx < 0 || d0_idx >= n_d0s) return false;
    
    double pt = d0_pt->at(d0_idx);
    double eta = d0_eta->at(d0_idx);
    double mass = d0_mass->at(d0_idx);
    
    // D0 PDG mass is 1864.84 MeV
    const double d0_pdg_mass = 1.86484;
    
    // Apply selection cuts
    if (pt < m_minPt) {
        // std::cout << "too low pt: " << pt << std::endl;
        return false;
    }
    if (eta > m_maxEta || eta < m_minEta) {
        // std::cout << "too high eta: " << eta << std::endl;
        return false;
    }
    
    if (TMath::Abs(mass - d0_pdg_mass) > m_d0MassWindow){
        // std::cout << "D0 mass outside window: " << mass << std::endl;
        return false;
    }
    // std::cout << "D0 selection passed for index " << d0_idx << ": "
    //           << "pT = " << pt << ", eta = " << eta << ", mass = " << mass << std::endl;
    return true;
}

bool D0Efficiency::PassesPIDSelection(int d0_idx) {
    // Find daughters of this D0
    int kaon_idx = -1, pion_idx = -1;
    
    for (size_t dau_idx = 0; dau_idx < dau_d0_idx->size(); ++dau_idx) {
        if (dau_d0_idx->at(dau_idx) == d0_idx) {
            int pid = dau_pid->at(dau_idx);
            if (abs(pid) == 321) { // Kaon
                kaon_idx = dau_idx;
            } else if (abs(pid) == 211) { // Pion
                pion_idx = dau_idx;
            }
        }
    }
    
    if (kaon_idx < 0 || pion_idx < 0) return false;
    
    // Check PID cuts
    if (dau_pnn_k->at(kaon_idx) < m_kaonPIDCut) return false;
    if (dau_pnn_pi->at(pion_idx) < m_pionPIDCut) return false;
    
    return true;
}

bool D0Efficiency::PassesQualityCuts(int d0_idx) {
    // Find daughters of this D0
    for (size_t dau_idx = 0; dau_idx < dau_d0_idx->size(); ++dau_idx) {
        if (dau_d0_idx->at(dau_idx) == d0_idx) {
            // Check quality cuts
            if (dau_prb_ghost->at(dau_idx) > m_ghostProbCut) return false;
            if (dau_trckChi2->at(dau_idx) > m_trackChi2Cut) return false;
            
            // Check daughter eta range
            double daughter_eta = dau_eta->at(dau_idx);
            if (daughter_eta < m_minEta || daughter_eta > m_maxEta) return false;
            
            // Check daughter momentum
            double daughter_px = dau_px->at(dau_idx);
            double daughter_py = dau_py->at(dau_idx);
            double daughter_pz = dau_pz->at(dau_idx);
            double daughter_p = TMath::Sqrt(daughter_px*daughter_px + daughter_py*daughter_py + daughter_pz*daughter_pz);
            if (daughter_p < m_minDaughterMomentum) return false;
        }
    }
    
    return true;
}

int D0Efficiency::FindMCMatch(int d0_idx) {
    if (!m_isMC) return -1;
    
    // Simple matching based on mc_d0_matched vector
    for (size_t mc_idx = 0; mc_idx < mc_d0_matched->size(); ++mc_idx) {
        if (mc_d0_matched->at(mc_idx) == d0_idx) {
            return mc_idx;
        }
    }
    
    return -1;
}

double D0Efficiency::CalculateDeltaR(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = phi1 - phi2;
    
    // Ensure dphi is in [-π, π]
    while (dphi > TMath::Pi()) dphi -= 2 * TMath::Pi();
    while (dphi < -TMath::Pi()) dphi += 2 * TMath::Pi();
    
    return TMath::Sqrt(deta * deta + dphi * dphi);
}

double D0Efficiency::CalculateMomentum(int d0_idx) {
    if (d0_idx < 0 || d0_idx >= n_d0s) return -1.0;
    
    double px = d0_px->at(d0_idx);
    double py = d0_py->at(d0_idx);
    double pz = d0_pz->at(d0_idx);
    
    return TMath::Sqrt(px * px + py * py + pz * pz);
}

double D0Efficiency::CalculateMCMomentum(int mc_idx) {
    if (mc_idx < 0 || mc_idx >= (int)mc_d0_px->size()) return -1.0;
    
    double px = mc_d0_px->at(mc_idx);
    double py = mc_d0_py->at(mc_idx);
    double pz = mc_d0_pz->at(mc_idx);
    
    return TMath::Sqrt(px * px + py * py + pz * pz);
}

bool D0Efficiency::GetBinIndices(double pt, double eta, int& ptBin, int& etaBin) {
    ptBin = -1;
    etaBin = -1;
    
    // Find pT bin
    for (size_t i = 0; i < m_ptBins.size() - 1; ++i) {
        if (pt >= m_ptBins[i] && pt < m_ptBins[i + 1]) {
            ptBin = i;
            break;
        }
    }
    
    // Find η bin
    for (size_t i = 0; i < m_etaBins.size() - 1; ++i) {
        if (eta >= m_etaBins[i] && eta < m_etaBins[i + 1]) {
            etaBin = i;
            break;
        }
    }
    
    return (ptBin >= 0 && etaBin >= 0);
}

bool D0Efficiency::PassesKaonPIDSelection(int d0_idx) {
    // Find kaon daughter of this D0
    int kaon_idx = -1;
    
    for (size_t dau_idx = 0; dau_idx < dau_d0_idx->size(); ++dau_idx) {
        if (dau_d0_idx->at(dau_idx) == d0_idx) {
            int pid = dau_pid->at(dau_idx);
            if (abs(pid) == 321) { // Kaon
                kaon_idx = dau_idx;
                break;
            }
        }
    }
    
    if (kaon_idx < 0) return false;
    
    // Check kaon PID cut
    return (dau_pnn_k->at(kaon_idx) >= m_kaonPIDCut);
}

bool D0Efficiency::PassesPionPIDSelection(int d0_idx) {
    // Find pion daughter of this D0
    int pion_idx = -1;
    
    for (size_t dau_idx = 0; dau_idx < dau_d0_idx->size(); ++dau_idx) {
        if (dau_d0_idx->at(dau_idx) == d0_idx) {
            int pid = dau_pid->at(dau_idx);
            if (abs(pid) == 211) { // Pion
                pion_idx = dau_idx;
                break;
            }
        }
    }
    
    if (pion_idx < 0) return false;
    
    // Check pion PID cut
    return (dau_pnn_pi->at(pion_idx) >= m_pionPIDCut);
}