#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <cmath>

void createResponseMatrix(const char *inputFile, const char *outputFile);

void nTupleMaker(const char *inputFile = "", int inputMC = 1)
{
    std::cout << "Starting D0 FF Analysis with minimal ntuple maker" << std::endl;

    // Use default files if empty string provided
    TString fInputFileName;
    if (strlen(inputFile) == 0)
    {
        if (inputMC)
        {
            fInputFileName = "/media/niviths/SSD2/lhcb_analysis_SSD/d0ff_outputs/localRunning_MC_ntuple_test.root";
        }
        else
        {
            fInputFileName = "/media/niviths/SSD2/lhcb_analysis_SSD/d0ff_outputs/localRunning_MeasuredData_ntuple_test.root";
        }
    }
    else
    {
        fInputFileName = inputFile;
    }

    TString fOutputFileName = fInputFileName;
    fOutputFileName.ReplaceAll(".root", "_filtered.root");

    std::cout << "Input file: " << fInputFileName << std::endl;
    std::cout << "Output file: " << fOutputFileName << std::endl;

    // Open input file
    TFile *inputRoot = TFile::Open(fInputFileName, "READ");
    if (!inputRoot || inputRoot->IsZombie())
    {
        std::cerr << "ERROR: Cannot open input file: " << fInputFileName << std::endl;
        return;
    }

    // Get input tree
    TTree *inputTree = (TTree *)inputRoot->Get("d0jets");
    if (!inputTree)
    {
        std::cerr << "ERROR: Cannot find 'd0jets' tree in input file" << std::endl;
        inputRoot->Close();
        return;
    }

    // Create output file
    TFile *outputRoot = TFile::Open(fOutputFileName, "RECREATE");
    if (!outputRoot || outputRoot->IsZombie())
    {
        std::cerr << "ERROR: Cannot create output file: " << fOutputFileName << std::endl;
        inputRoot->Close();
        return;
    }

    // Create output tree
    TTree *outputTree = new TTree("FragmNtuple", "D0 Jet Fragmentation Ntuple");

    // Setup branches for reading input tree
    // Event info
    int evt_num, run_num, n_pvs;

    // Jet info vectors
    std::vector<float> *jet_pt = nullptr;
    std::vector<float> *jet_eta = nullptr;
    std::vector<float> *jet_phi = nullptr;
    std::vector<float> *jet_mass = nullptr;
    std::vector<int> *jet_n_const = nullptr;
    std::vector<int> *jet_n_charged = nullptr;
    std::vector<int> *jet_n_neutral = nullptr;
    std::vector<int> *jet_n_d0 = nullptr;

    // D0 info vectors
    std::vector<float> *d0_pt = nullptr;
    std::vector<float> *d0_eta = nullptr;
    std::vector<float> *d0_phi = nullptr;
    std::vector<float> *d0_mass = nullptr;
    std::vector<float> *d0_vtx_chi2 = nullptr;
    std::vector<float> *d0_ip = nullptr;
    std::vector<float> *d0_ip_chi2 = nullptr;
    std::vector<float> *d0_fd = nullptr;
    std::vector<float> *d0_fd_chi2 = nullptr;
    std::vector<float> *d0_DOCA = nullptr;
    std::vector<int> *d0_jet_idx = nullptr;
    std::vector<int> *d0_in_jet = nullptr;
    std::vector<float> *d0_z = nullptr;
    std::vector<float> *d0_jet_dr = nullptr;

    // D0 daughter info vectors
    std::vector<int> *dau_pid = nullptr;
    std::vector<float> *dau_pt = nullptr;
    std::vector<float> *dau_px = nullptr;
    std::vector<float> *dau_py = nullptr;
    std::vector<float> *dau_pz = nullptr;
    std::vector<float> *dau_eta = nullptr;
    std::vector<float> *dau_e = nullptr;
    std::vector<float> *dau_phi = nullptr;
    std::vector<int> *dau_d0_idx = nullptr;
    std::vector<float> *dau_pnn_k = nullptr;
    std::vector<float> *dau_pnn_pi = nullptr;
    std::vector<float> *dau_prb_ghost = nullptr;
    std::vector<float> *dau_trckChi2 = nullptr;

    // MC truth info (only used if inputMC==1)
    std::vector<int> *mc_d0_pid = nullptr;
    std::vector<float> *mc_d0_pt = nullptr;
    std::vector<int> *mc_d0_origin = nullptr;
    std::vector<int> *mc_d0_matched = nullptr;

    // Set branch addresses for input tree
    inputTree->SetBranchAddress("evt_num", &evt_num);
    inputTree->SetBranchAddress("run_num", &run_num);
    inputTree->SetBranchAddress("n_pvs", &n_pvs);

    inputTree->SetBranchAddress("jet_pt", &jet_pt);
    inputTree->SetBranchAddress("jet_eta", &jet_eta);
    inputTree->SetBranchAddress("jet_phi", &jet_phi);
    inputTree->SetBranchAddress("jet_mass", &jet_mass);
    inputTree->SetBranchAddress("jet_n_const", &jet_n_const);
    inputTree->SetBranchAddress("jet_n_charged", &jet_n_charged);
    inputTree->SetBranchAddress("jet_n_neutral", &jet_n_neutral);
    inputTree->SetBranchAddress("jet_n_d0", &jet_n_d0);

    inputTree->SetBranchAddress("d0_pt", &d0_pt);
    inputTree->SetBranchAddress("d0_eta", &d0_eta);
    inputTree->SetBranchAddress("d0_phi", &d0_phi);
    inputTree->SetBranchAddress("d0_mass", &d0_mass);
    inputTree->SetBranchAddress("d0_vtx_chi2", &d0_vtx_chi2);
    inputTree->SetBranchAddress("d0_ip", &d0_ip);
    inputTree->SetBranchAddress("d0_ip_chi2", &d0_ip_chi2);
    inputTree->SetBranchAddress("d0_fd", &d0_fd);
    inputTree->SetBranchAddress("d0_fd_chi2", &d0_fd_chi2);
    inputTree->SetBranchAddress("d0_DOCA", &d0_DOCA);
    inputTree->SetBranchAddress("d0_jet_idx", &d0_jet_idx);
    inputTree->SetBranchAddress("d0_in_jet", &d0_in_jet);
    inputTree->SetBranchAddress("d0_z", &d0_z);
    inputTree->SetBranchAddress("d0_jet_dr", &d0_jet_dr);

    inputTree->SetBranchAddress("dau_pid", &dau_pid);
    inputTree->SetBranchAddress("dau_pt", &dau_pt);
    inputTree->SetBranchAddress("dau_px", &dau_px);
    inputTree->SetBranchAddress("dau_py", &dau_py);
    inputTree->SetBranchAddress("dau_pz", &dau_pz);
    inputTree->SetBranchAddress("dau_eta", &dau_eta);
    inputTree->SetBranchAddress("dau_e", &dau_e);
    inputTree->SetBranchAddress("dau_phi", &dau_phi);
    inputTree->SetBranchAddress("dau_d0_idx", &dau_d0_idx);
    inputTree->SetBranchAddress("dau_pnn_k", &dau_pnn_k);
    inputTree->SetBranchAddress("dau_pnn_pi", &dau_pnn_pi);
    inputTree->SetBranchAddress("dau_prb_ghost", &dau_prb_ghost);
    inputTree->SetBranchAddress("dau_trckChi2", &dau_trckChi2);

    // Setup MC truth branches if using MC
    if (inputMC)
    {
        inputTree->SetBranchAddress("mc_d0_pid", &mc_d0_pid);
        inputTree->SetBranchAddress("mc_d0_pt", &mc_d0_pt);
        inputTree->SetBranchAddress("mc_d0_origin", &mc_d0_origin);
        inputTree->SetBranchAddress("mc_d0_matched", &mc_d0_matched);
    }

    // Variables for output tree
    float v_tagdR = 0;
    float v_tagMass = 0;
    float v_tagPt = 0;
    float v_tagEta = 0;
    float v_tag_idx_jet = 0;
    float v_tag_decVtxChi2 = 0;
    float v_tag_logdecVtxChi2 = 0;
    float v_jetPt = 0;
    float v_jetEta = 0;
    float v_jetnConst = 0;
    float v_tagZ = 0;
    float v_isPrimary = 0;
    float v_KprobNNK = 0;
    float v_KprobGhost = 0;
    float v_KTrckChi2 = 0;
    float v_piPprobNNpi = 0;
    float v_piPprobGhost = 0;
    float v_piPTrckChi2 = 0;
    float v_decayVtxChi2 = 0;
    float v_Dist1 = 0; // Will use d0_DOCA

    // Create branches for output tree
    outputTree->Branch("tagJetdR", &v_tagdR, "tagJetdR/F");
    outputTree->Branch("tagMass", &v_tagMass, "tagMass/F");
    outputTree->Branch("tagPt", &v_tagPt, "tagPt/F");
    outputTree->Branch("tagEta", &v_tagEta, "tagEta/F");
    outputTree->Branch("tagidxjet", &v_tag_idx_jet, "tagidxjet/F");
    outputTree->Branch("tag_ip_chi2", &v_tag_decVtxChi2, "tag_ip_chi2/F");
    outputTree->Branch("log_tag_ipchi2", &v_tag_logdecVtxChi2, "log_tag_ipchi2/F");
    outputTree->Branch("jetPt", &v_jetPt, "jetPt/F");
    outputTree->Branch("jetEta", &v_jetEta, "jetEta/F");
    outputTree->Branch("jetnConst", &v_jetnConst, "jetnConst/F");
    outputTree->Branch("tagZ", &v_tagZ, "tagZ/F");
    outputTree->Branch("piPprobNNpi", &v_piPprobNNpi, "piPprobNN/F");
    outputTree->Branch("piPprobGhost", &v_piPprobGhost, "piPprobGhost/F");
    outputTree->Branch("piPTrckChi2", &v_piPTrckChi2, "piPTrckChi2/F");
    outputTree->Branch("KprobNNK", &v_KprobNNK, "KprobNNK/F");
    outputTree->Branch("KprobGhost", &v_KprobGhost, "KprobGhost/F");
    outputTree->Branch("KTrckChi2", &v_KTrckChi2, "KTrckChi2/F");
    outputTree->Branch("decayVtxChi2", &v_decayVtxChi2, "decayVtxChi2/F");
    outputTree->Branch("Distance1", &v_Dist1, "Distance1/F");
    outputTree->Branch("isPrimary", &v_isPrimary, "isPrimary/F");

    // Additional MC branches if needed
    if (inputMC)
    {
        outputTree->Branch("isPrimary", &v_isPrimary, "isPrimary/F");
    }

    // Process events
    Long64_t nEntries = inputTree->GetEntries();
    std::cout << "Processing " << nEntries << " events" << std::endl;

    int events_processed = 0;
    int d0_accepted = 0;

    for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++)
    {
        inputTree->GetEntry(iEntry);

        if (iEntry % 1000 == 0)
        {
            std::cout << "Processing event " << iEntry << " of " << nEntries << std::endl;
        }

        events_processed++;

        // Process each D0 in the event
        for (size_t i_d0 = 0; i_d0 < d0_pt->size(); i_d0++)
        {
            // Basic D0 selection cuts
            if ((*d0_pt)[i_d0] < 2.0)
                continue; // Minimum pT
            if ((*d0_eta)[i_d0] < 2.0 || (*d0_eta)[i_d0] > 4.5)
                continue; // Eta acceptance

            // D0 mass window cut
            if (std::abs((*d0_mass)[i_d0] - 1.865) > 0.05)
                continue; // Mass window: ±50 MeV around nominal D0 mass

            // Check if D0 is associated with a jet
            int jet_idx = (*d0_jet_idx)[i_d0];
            if (jet_idx < 0)
                continue; // Require associated jet

            // Jet selection cuts
            if ((*jet_pt)[jet_idx] < 5.0)
                continue; // Minimum jet pT
            if ((*jet_eta)[jet_idx] < 2.5 || (*jet_eta)[jet_idx] > 4.0)
                continue; // Jet eta range

            // D0 vertex quality
            if ((*d0_vtx_chi2)[i_d0] > 10.0)
                continue; // Require good vertex fit

            // Find kaon and pion among daughters
            float kaon_pnn_k = -1;
            float kaon_ghost_prob = -1;
            float kaon_chi2 = -1;
            float pion_pnn_pi = -1;
            float pion_ghost_prob = -1;
            float pion_chi2 = -1;
            TLorentzVector kaon_vec, pion_vec;
            bool kaon_found = false;
            bool pion_found = false;

            for (size_t i_dau = 0; i_dau < dau_pid->size(); i_dau++)
            {
                if ((*dau_d0_idx)[i_dau] != i_d0)
                    continue; // Only daughters of this D0

                if (std::abs((*dau_pid)[i_dau]) == 321)
                { // Kaon
                    kaon_pnn_k = (*dau_pnn_k)[i_dau];
                    kaon_ghost_prob = (*dau_prb_ghost)[i_dau];
                    kaon_chi2 = (*dau_trckChi2)[i_dau];

                    // Create kaon 4-vector
                    kaon_vec.SetPxPyPzE(
                        (*dau_px)[i_dau],
                        (*dau_py)[i_dau],
                        (*dau_pz)[i_dau],
                        (*dau_e)[i_dau]);
                    kaon_found = true;
                }
                else if (std::abs((*dau_pid)[i_dau]) == 211)
                { // Pion
                    pion_pnn_pi = (*dau_pnn_pi)[i_dau];
                    pion_ghost_prob = (*dau_prb_ghost)[i_dau];
                    pion_chi2 = (*dau_trckChi2)[i_dau];

                    // Create pion 4-vector
                    pion_vec.SetPxPyPzE(
                        (*dau_px)[i_dau],
                        (*dau_py)[i_dau],
                        (*dau_pz)[i_dau],
                        (*dau_e)[i_dau]);
                    pion_found = true;
                }
            }

            // Skip if we didn't find both K and π
            if (!kaon_found || !pion_found)
                continue;

            // Apply acceptance cuts to kaon and pion
            double kaon_eta = kaon_vec.PseudoRapidity();
            double pion_eta = pion_vec.PseudoRapidity();
            double kaon_pt = kaon_vec.Pt();
            double pion_pt = pion_vec.Pt();
            double kaon_p = kaon_vec.P();
            double pion_p = pion_vec.P();

            // Eta acceptance cuts
            if (kaon_eta < 2.0 || kaon_eta > 4.5 || pion_eta < 2.0 || pion_eta > 4.5)
                continue;

            // Momentum cuts - minimum pT and total momentum
            if (kaon_pt < 0.25 || pion_pt < 0.25)
                continue; // Minimum pT cut
            if (kaon_p < 2.0 || pion_p < 2.0)
                continue; // Minimum momentum cut

            // PID quality cuts - now applied after the acceptance cuts
            if (kaon_pnn_k < 0.1)
                continue; // Loose kaon ID
            if (pion_pnn_pi < 0.1)
                continue; // Loose pion ID
            if (kaon_ghost_prob > 0.5 || pion_ghost_prob > 0.5)
                continue; // Ghost probability cut

            // Fill output variables
            v_tagdR = (*d0_jet_dr)[i_d0];
            v_tagMass = (*d0_mass)[i_d0];
            v_tagPt = (*d0_pt)[i_d0];
            v_tagEta = (*d0_eta)[i_d0];
            v_tag_idx_jet = jet_idx;
            v_tag_decVtxChi2 = (*d0_ip_chi2)[i_d0];
            v_tag_logdecVtxChi2 = (*d0_ip_chi2)[i_d0] > 0 ? std::log10((*d0_ip_chi2)[i_d0]) : -999;
            v_jetPt = (*jet_pt)[jet_idx];
            v_jetEta = (*jet_eta)[jet_idx];
            v_jetnConst = (*jet_n_const)[jet_idx];
            v_tagZ = (*d0_z)[i_d0];
            v_KprobNNK = kaon_pnn_k;
            v_KprobGhost = kaon_ghost_prob;
            v_KTrckChi2 = kaon_chi2;
            v_piPprobNNpi = pion_pnn_pi;
            v_piPprobGhost = pion_ghost_prob;
            v_piPTrckChi2 = pion_chi2;
            v_decayVtxChi2 = (*d0_vtx_chi2)[i_d0];
            v_Dist1 = (*d0_DOCA)[i_d0];

            // MC specific variables
            if (inputMC && i_d0 < mc_d0_origin->size())
            {
                v_isPrimary = (*mc_d0_origin)[i_d0] == 1 ? 1 : 0; // Is it a prompt D0?
            }
            else
            {
                v_isPrimary = -1; // Unknown for data
            }

            // Fill the output tree
            outputTree->Fill();
            d0_accepted++;
        }
    }

    // Print summary
    std::cout << "Processing complete!" << std::endl;
    std::cout << "Events processed: " << events_processed << std::endl;
    std::cout << "D0 candidates accepted: " << d0_accepted << std::endl;
    std::cout << "Acceptance rate: " << 100.0 * d0_accepted / events_processed << "%" << std::endl;

    // Write and close output file
    outputTree->Write();
    outputRoot->Write();
    outputRoot->Close();
    inputRoot->Close();

    std::cout << "Output saved to: " << fOutputFileName << std::endl;

    // Create response matrix tree if processing MC
    if (inputMC)
    {
        std::cout << "Creating response matrix tree..." << std::endl;
        createResponseMatrix(fInputFileName, fOutputFileName.ReplaceAll("_filtered.root", "_response.root"));
    }
}

// Function to create response matrix
void createResponseMatrix(const char *inputFile, const char *outputFile)
{
    // Open input file
    TFile *inputRoot = TFile::Open(inputFile, "READ");
    if (!inputRoot || inputRoot->IsZombie())
    {
        std::cerr << "ERROR: Cannot open input file for response matrix: " << inputFile << std::endl;
        return;
    }

    // Get input tree
    TTree *inputTree = (TTree *)inputRoot->Get("d0jets");
    if (!inputTree)
    {
        std::cerr << "ERROR: Cannot find 'd0jets' tree in input file" << std::endl;
        inputRoot->Close();
        return;
    }

    // Create output file
    TFile *outputRoot = TFile::Open(outputFile, "RECREATE");
    if (!outputRoot || outputRoot->IsZombie())
    {
        std::cerr << "ERROR: Cannot create response output file: " << outputFile << std::endl;
        inputRoot->Close();
        return;
    }

    // Create output tree
    TTree *responseTree = new TTree("Response", "Jet-D0 Response Matrix");

    // Variables for response matrix tree
    float r_jet_pt_det = 0;
    float r_jet_eta_det = 0;
    float r_jet_phi_det = 0;
    float r_jet_nconst_det = 0;
    float r_d0_pt_det = 0;
    float r_d0_eta_det = 0;
    float r_d0_phi_det = 0;
    float r_d0_mass_det = 0;
    float r_d0_z_det = 0;

    float r_jet_pt_mc = 0;
    float r_jet_eta_mc = 0;
    float r_jet_phi_mc = 0;
    float r_jet_nconst_mc = 0;
    float r_d0_pt_mc = 0;
    float r_d0_eta_mc = 0;
    float r_d0_phi_mc = 0;
    float r_d0_mass_mc = 0;
    float r_d0_z_mc = 0;

    float r_jet_dr = 0;      // dR between matched jets
    float r_d0_dr = 0;       // dR between matched D0s
    int r_d0_is_primary = 0; // Is MC D0 primary
    int r_jet_ntags_det = 0; // Number of D0 tags in det jet
    int r_jet_ntags_mc = 0;  // Number of D0 tags in MC jet

    // Create branches
    responseTree->Branch("jet_pt_det", &r_jet_pt_det, "jet_pt_det/F");
    responseTree->Branch("jet_eta_det", &r_jet_eta_det, "jet_eta_det/F");
    responseTree->Branch("jet_phi_det", &r_jet_phi_det, "jet_phi_det/F");
    responseTree->Branch("jet_nconst_det", &r_jet_nconst_det, "jet_nconst_det/F");
    responseTree->Branch("d0_pt_det", &r_d0_pt_det, "d0_pt_det/F");
    responseTree->Branch("d0_eta_det", &r_d0_eta_det, "d0_eta_det/F");
    responseTree->Branch("d0_phi_det", &r_d0_phi_det, "d0_phi_det/F");
    responseTree->Branch("d0_mass_det", &r_d0_mass_det, "d0_mass_det/F");
    responseTree->Branch("d0_z_det", &r_d0_z_det, "d0_z_det/F");

    responseTree->Branch("jet_pt_mc", &r_jet_pt_mc, "jet_pt_mc/F");
    responseTree->Branch("jet_eta_mc", &r_jet_eta_mc, "jet_eta_mc/F");
    responseTree->Branch("jet_phi_mc", &r_jet_phi_mc, "jet_phi_mc/F");
    responseTree->Branch("jet_nconst_mc", &r_jet_nconst_mc, "jet_nconst_mc/F");
    responseTree->Branch("d0_pt_mc", &r_d0_pt_mc, "d0_pt_mc/F");
    responseTree->Branch("d0_eta_mc", &r_d0_eta_mc, "d0_eta_mc/F");
    responseTree->Branch("d0_phi_mc", &r_d0_phi_mc, "d0_phi_mc/F");
    responseTree->Branch("d0_mass_mc", &r_d0_mass_mc, "d0_mass_mc/F");
    responseTree->Branch("d0_z_mc", &r_d0_z_mc, "d0_z_mc/F");

    responseTree->Branch("jet_dr", &r_jet_dr, "jet_dr/F");
    responseTree->Branch("d0_dr", &r_d0_dr, "d0_dr/F");
    responseTree->Branch("d0_is_primary", &r_d0_is_primary, "d0_is_primary/I");
    responseTree->Branch("jet_ntags_det", &r_jet_ntags_det, "jet_ntags_det/I");
    responseTree->Branch("jet_ntags_mc", &r_jet_ntags_mc, "jet_ntags_mc/I");

    // Set up input branches
    std::vector<float> *jet_pt = nullptr;
    std::vector<float> *jet_eta = nullptr;
    std::vector<float> *jet_phi = nullptr;
    std::vector<int> *jet_n_const = nullptr;
    std::vector<int> *jet_n_d0 = nullptr;

    std::vector<float> *d0_pt = nullptr;
    std::vector<float> *d0_eta = nullptr;
    std::vector<float> *d0_phi = nullptr;
    std::vector<float> *d0_mass = nullptr;
    std::vector<int> *d0_jet_idx = nullptr;
    std::vector<int> *d0_in_jet = nullptr;
    std::vector<float> *d0_z = nullptr;
    std::vector<float> *dau_px = nullptr;
    std::vector<float> *dau_py = nullptr;
    std::vector<float> *dau_pz = nullptr;
    std::vector<float> *dau_e = nullptr;
    std::vector<int> *dau_pid = nullptr;
    std::vector<int> *dau_d0_idx = nullptr;

    // MC truth branches
    std::vector<int> *mc_d0_pid = nullptr;
    std::vector<float> *mc_d0_pt = nullptr;
    std::vector<float> *mc_d0_eta = nullptr;
    std::vector<float> *mc_d0_phi = nullptr;
    std::vector<float> *mc_d0_mass = nullptr;
    std::vector<int> *mc_d0_origin = nullptr;
    std::vector<int> *mc_d0_matched = nullptr;
    std::vector<int> *mc_d0_jet_idx = nullptr;
    std::vector<float> *mc_d0_z = nullptr;
    std::vector<float> *mc_jet_px = nullptr;
    std::vector<float> *mc_jet_py = nullptr;
    std::vector<float> *mc_jet_pz = nullptr;
    std::vector<float> *mc_jet_e = nullptr;
    std::vector<int> *mc_jet_n_const = nullptr;
    std::vector<int> *mc_jet_n_chr = nullptr;
    std::vector<int> *mc_jet_n_neu = nullptr;
    std::vector<int> *mc_dau_pid = nullptr;
    std::vector<float> *mc_dau_px = nullptr;
    std::vector<float> *mc_dau_py = nullptr;
    std::vector<float> *mc_dau_pz = nullptr;
    std::vector<float> *mc_dau_e = nullptr;
    std::vector<int> *mc_dau_d0_idx = nullptr;

    // Set branch addresses
    inputTree->SetBranchAddress("jet_pt", &jet_pt);
    inputTree->SetBranchAddress("jet_eta", &jet_eta);
    inputTree->SetBranchAddress("jet_phi", &jet_phi);
    inputTree->SetBranchAddress("jet_n_const", &jet_n_const);
    inputTree->SetBranchAddress("jet_n_d0", &jet_n_d0);

    inputTree->SetBranchAddress("d0_pt", &d0_pt);
    inputTree->SetBranchAddress("d0_eta", &d0_eta);
    inputTree->SetBranchAddress("d0_phi", &d0_phi);
    inputTree->SetBranchAddress("d0_mass", &d0_mass);
    inputTree->SetBranchAddress("d0_jet_idx", &d0_jet_idx);
    inputTree->SetBranchAddress("d0_in_jet", &d0_in_jet);
    inputTree->SetBranchAddress("d0_z", &d0_z);

    inputTree->SetBranchAddress("dau_px", &dau_px);
    inputTree->SetBranchAddress("dau_py", &dau_py);
    inputTree->SetBranchAddress("dau_pz", &dau_pz);
    inputTree->SetBranchAddress("dau_e", &dau_e);
    inputTree->SetBranchAddress("dau_pid", &dau_pid);
    inputTree->SetBranchAddress("dau_d0_idx", &dau_d0_idx);

    inputTree->SetBranchAddress("mc_d0_pid", &mc_d0_pid);
    inputTree->SetBranchAddress("mc_d0_pt", &mc_d0_pt);
    inputTree->SetBranchAddress("mc_d0_eta", &mc_d0_eta);
    inputTree->SetBranchAddress("mc_d0_phi", &mc_d0_phi);
    inputTree->SetBranchAddress("mc_d0_mass", &mc_d0_mass);
    inputTree->SetBranchAddress("mc_d0_origin", &mc_d0_origin);
    inputTree->SetBranchAddress("mc_d0_matched", &mc_d0_matched);
    inputTree->SetBranchAddress("mc_d0_jet_idx", &mc_d0_jet_idx);
    inputTree->SetBranchAddress("mc_d0_z", &mc_d0_z);
    inputTree->SetBranchAddress("mc_jet_px", &mc_jet_px);
    inputTree->SetBranchAddress("mc_jet_py", &mc_jet_py);
    inputTree->SetBranchAddress("mc_jet_pz", &mc_jet_pz);
    inputTree->SetBranchAddress("mc_jet_e", &mc_jet_e);
    inputTree->SetBranchAddress("mc_jet_n_const", &mc_jet_n_const);
    inputTree->SetBranchAddress("mc_jet_n_charged", &mc_jet_n_chr);
    inputTree->SetBranchAddress("mc_jet_n_neutral", &mc_jet_n_neu);
    inputTree->SetBranchAddress("mc_dau_pid", &mc_dau_pid);
    inputTree->SetBranchAddress("mc_dau_px", &mc_dau_px);
    inputTree->SetBranchAddress("mc_dau_py", &mc_dau_py);
    inputTree->SetBranchAddress("mc_dau_pz", &mc_dau_pz);
    inputTree->SetBranchAddress("mc_dau_e", &mc_dau_e);
    inputTree->SetBranchAddress("mc_dau_d0_idx", &mc_dau_d0_idx);

    // Define a function to calculate dR between two eta,phi points
    auto deltaR = [](float eta1, float phi1, float eta2, float phi2)
    {
        float deta = eta1 - eta2;
        float dphi = phi1 - phi2;
        while (dphi > M_PI)
            dphi -= 2 * M_PI;
        while (dphi < -M_PI)
            dphi += 2 * M_PI;
        return sqrt(deta * deta + dphi * dphi);
    };

    // Process events
    Long64_t nEntries = inputTree->GetEntries();
    int matches_found = 0;

    for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++)
    {
        if (iEntry % 1000 == 0)
        {
            std::cout << "Processing event " << iEntry << " of " << nEntries
                      << " for response matrix" << std::endl;
        }

        inputTree->GetEntry(iEntry);

        // Skip empty events
        if (!mc_d0_pid || mc_d0_pid->empty() || !d0_pt || d0_pt->empty())
        {
            continue;
        }

        // Loop through MC D0s and find their matched reconstructed D0
        for (size_t iMC = 0; iMC < mc_d0_pid->size(); iMC++)
        {
            // Skip MC D0s not matched to a jet
            if (iMC >= mc_d0_jet_idx->size() || (*mc_d0_jet_idx)[iMC] < 0)
            {
                continue;
            }

            // Make sure this MC D0 has matching reconstructed D0
            int matched_d0_idx = -1;
            if (iMC < mc_d0_matched->size())
            {
                matched_d0_idx = (*mc_d0_matched)[iMC];
            }

            // Skip if no match
            if (matched_d0_idx < 0 || matched_d0_idx >= (int)d0_pt->size())
            {
                continue;
            }

            // Get MC D0 jet index
            int mc_jet_idx = (*mc_d0_jet_idx)[iMC];

            // Get reconstructed D0 jet index
            int reco_jet_idx = (*d0_jet_idx)[matched_d0_idx];

            // Skip if either D0 not associated with a jet
            if (mc_jet_idx < 0 || reco_jet_idx < 0)
            {
                continue;
            }

            // MC truth info
            r_d0_pt_mc = (*mc_d0_pt)[iMC];
            r_d0_eta_mc = (*mc_d0_eta)[iMC];
            r_d0_phi_mc = (*mc_d0_phi)[iMC];
            r_d0_mass_mc = (*mc_d0_mass)[iMC];
            r_d0_z_mc = (*mc_d0_z)[iMC];
            r_d0_is_primary = (iMC < mc_d0_origin->size()) ? ((*mc_d0_origin)[iMC] == 1) : 0;

            // Get MC jet vectors
            TLorentzVector mc_jet;
            if (mc_jet_idx >= 0 && mc_jet_idx < static_cast<int>(mc_jet_px->size()))
            {
                mc_jet.SetPxPyPzE(
                    (*mc_jet_px)[mc_jet_idx],
                    (*mc_jet_py)[mc_jet_idx],
                    (*mc_jet_pz)[mc_jet_idx],
                    (*mc_jet_e)[mc_jet_idx]);

                // Use actual MC jet properties
                r_jet_pt_mc = mc_jet.Pt();
                r_jet_eta_mc = mc_jet.Eta();
                r_jet_phi_mc = mc_jet.Phi();

                // Use actual constituent counts if available
                if (mc_jet_n_const && mc_jet_idx < static_cast<int>(mc_jet_n_const->size()))
                {
                    r_jet_nconst_mc = (*mc_jet_n_const)[mc_jet_idx];
                }
                else if (mc_jet_n_chr && mc_jet_n_neu &&
                         mc_jet_idx < static_cast<int>(mc_jet_n_chr->size()) &&
                         mc_jet_idx < static_cast<int>(mc_jet_n_neu->size()))
                {
                    r_jet_nconst_mc = (*mc_jet_n_chr)[mc_jet_idx] + (*mc_jet_n_neu)[mc_jet_idx];
                }
                else
                {
                    r_jet_nconst_mc = 0;
                }
            }
            else
            {
                // Fallback if index is invalid
                r_jet_pt_mc = r_d0_pt_mc / r_d0_z_mc;
                r_jet_eta_mc = r_d0_eta_mc;
                r_jet_phi_mc = r_d0_phi_mc;
                r_jet_nconst_mc = 0;
            }

            // Count D0s in MC jet - this part is correct and can stay
            r_jet_ntags_mc = 0;
            for (size_t id0 = 0; id0 < mc_d0_jet_idx->size(); id0++)
            {
                if ((*mc_d0_jet_idx)[id0] == mc_jet_idx)
                {
                    r_jet_ntags_mc++;
                }
            }

            // Reconstructed info
            r_d0_pt_det = (*d0_pt)[matched_d0_idx];
            r_d0_eta_det = (*d0_eta)[matched_d0_idx];
            r_d0_phi_det = (*d0_phi)[matched_d0_idx];
            r_d0_mass_det = (*d0_mass)[matched_d0_idx];
            r_d0_z_det = (*d0_z)[matched_d0_idx];

            // Reconstructed jet info
            if (reco_jet_idx < (int)jet_pt->size())
            {
                r_jet_pt_det = (*jet_pt)[reco_jet_idx];
                r_jet_eta_det = (*jet_eta)[reco_jet_idx];
                r_jet_phi_det = (*jet_phi)[reco_jet_idx];
                r_jet_nconst_det = (*jet_n_const)[reco_jet_idx];
                r_jet_ntags_det = (*jet_n_d0)[reco_jet_idx];
            }
            else
            {
                continue; // Skip if jet index out of bounds
            }

            // Calculate distances
            r_d0_dr = deltaR(r_d0_eta_mc, r_d0_phi_mc, r_d0_eta_det, r_d0_phi_det);
            r_jet_dr = deltaR(r_jet_eta_mc, r_jet_phi_mc, r_jet_eta_det, r_jet_phi_det);

            // Skip matches with too large distance
            if (r_jet_dr > 0.4)
            {
                continue;
            }

            // Check D0 daughter acceptance
            bool daughtersInAcceptance = false;

            // Get the kaon and pion 4-vectors from dau_px, dau_py, dau_pz, dau_e

            TLorentzVector kaon, pion;
            bool kaon_found = false;
            bool pion_found = false;

            try
            {
                for (size_t i = 0; i < dau_pid->size(); i++)
                {
                    // Only look at daughters of the matched reco D0
                    if ((*dau_d0_idx)[i] != matched_d0_idx)
                        continue;

                    if (std::abs((*dau_pid)[i]) == 321)
                    { // Kaon
                        kaon.SetPxPyPzE(
                            (*dau_px)[i],
                            (*dau_py)[i],
                            (*dau_pz)[i],
                            (*dau_e)[i]);
                        kaon_found = true;
                    }
                    else if (std::abs((*dau_pid)[i]) == 211)
                    { // Pion
                        pion.SetPxPyPzE(
                            (*dau_px)[i],
                            (*dau_py)[i],
                            (*dau_pz)[i],
                            (*dau_e)[i]);
                        pion_found = true;
                    }
                }

                // Apply acceptance cuts if both particles found
                if (kaon_found && pion_found)
                {
                    // Eta acceptance cuts
                    bool eta_ok = (kaon.Eta() > 2.0 && kaon.Eta() < 4.5 &&
                                   pion.Eta() > 2.0 && pion.Eta() < 4.5);

                    // Momentum cuts
                    bool pt_ok = (kaon.Pt() > 0.25 && pion.Pt() > 0.25);
                    bool p_ok = (kaon.P() > 2.0 && pion.P() > 2.0);

                    daughtersInAcceptance = eta_ok && pt_ok && p_ok;
                }

                // Skip if daughters don't fulfill acceptance requirements
                if (!daughtersInAcceptance)
                {
                    continue;
                }
            }
            catch (...)
            {
                std::cout << "Error checking daughter acceptance, skipping" << std::endl;
                continue;
            }

            // Now check MC truth daughter acceptance
            bool mcDaughtersInAcceptance = false;
            TLorentzVector mc_kaon, mc_pion;
            bool mc_kaon_found = false;
            bool mc_pion_found = false;

            try
            {
                // Find MC kaon and pion for this D0
                for (size_t i = 0; i < mc_dau_pid->size(); i++)
                {
                    // Only look at daughters of this MC D0
                    if ((*mc_dau_d0_idx)[i] != iMC)
                        continue;

                    if (std::abs((*mc_dau_pid)[i]) == 321)
                    { // Kaon
                        mc_kaon.SetPxPyPzE(
                            (*mc_dau_px)[i],
                            (*mc_dau_py)[i],
                            (*mc_dau_pz)[i],
                            (*mc_dau_e)[i]);
                        mc_kaon_found = true;
                    }
                    else if (std::abs((*mc_dau_pid)[i]) == 211)
                    { // Pion
                        mc_pion.SetPxPyPzE(
                            (*mc_dau_px)[i],
                            (*mc_dau_py)[i],
                            (*mc_dau_pz)[i],
                            (*mc_dau_e)[i]);
                        mc_pion_found = true;
                    }
                }

                // Apply acceptance cuts if both particles found
                if (mc_kaon_found && mc_pion_found)
                {
                    // Eta acceptance cuts
                    bool eta_ok = (mc_kaon.Eta() > 2.0 && mc_kaon.Eta() < 4.5 &&
                                   mc_pion.Eta() > 2.0 && mc_pion.Eta() < 4.5);

                    // Momentum cuts
                    bool pt_ok = (mc_kaon.Pt() > 0.25 && mc_pion.Pt() > 0.25);
                    bool p_ok = (mc_kaon.P() > 2.0 && mc_pion.P() > 2.0);

                    mcDaughtersInAcceptance = eta_ok && pt_ok && p_ok;
                }

                // Skip if MC daughters don't fulfill acceptance requirements
                if (!mcDaughtersInAcceptance)
                {
                    continue;
                }
            }
            catch (std::exception &e)
            {
                std::cout << "Error checking MC daughter acceptance: " << e.what() << std::endl;
                continue;
            }
            catch (...)
            {
                std::cout << "Unknown error checking MC daughter acceptance, skipping" << std::endl;
                continue;
            }

            // If we've reached here, both reconstructed and MC daughters are in acceptance
            responseTree->Fill();
            matches_found++;
        }
    }

    std::cout << "Found " << matches_found << " matched D0-jet pairs for response matrix" << std::endl;

    // Write and close
    outputRoot->cd();
    responseTree->Write();
    outputRoot->Close();
    inputRoot->Close();

    std::cout << "Response matrix saved to: " << outputFile << std::endl;
}