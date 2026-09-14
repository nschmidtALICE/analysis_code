#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <TH2D.h>
#include <TH1.h>
#include <fstream>
#include <string>
// Progress helper for the response matrix loop (reuse same style)
auto printProgress = [](Long64_t current, Long64_t total)
{
    const int barWidth = 50;
    double fraction = total > 0 ? double(current + 1) / double(total) : 1.0;
    int pos = static_cast<int>(barWidth * fraction);
    std::cout << "\r[";
    for (int i = 0; i < barWidth; ++i)
    {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << static_cast<int>(fraction * 100.0) << "% (" << (current + 1) << "/" << total << ")" << std::flush;
    if (current + 1 == total)
        std::cout << std::endl;
};

auto calcRapidityFromEnergyPz = [](float energy, float pz)
{
    if (energy <= std::fabs(pz))
        return -999.f;

    double numerator = static_cast<double>(energy + pz);
    double denominator = static_cast<double>(energy - pz);
    if (numerator <= 0.0 || denominator <= 0.0)
        return -999.f;

    return static_cast<float>(0.5 * std::log(numerator / denominator));
};

auto calcRapidityFromPtEtaMass = [](float pt, float eta, float mass)
{
    double ptD = static_cast<double>(pt);
    double etaD = static_cast<double>(eta);
    double massD = static_cast<double>(mass);
    double pz = ptD * std::sinh(etaD);
    double momentum = ptD * std::cosh(etaD);
    double energySquared = momentum * momentum + massD * massD;
    if (energySquared <= 0.0)
        return eta;

    double energy = std::sqrt(energySquared);
    return calcRapidityFromEnergyPz(static_cast<float>(energy), static_cast<float>(pz));
};

void createResponseMatrix(const char *inputFile, const char *outputFile, int jetMode, bool useWeights=false, const char* multweightsfile = nullptr, bool useJetPtWeights=false);

// doJet/doJetMode semantics:
// - doJet=false or doJetMode=0: No jet requirements; do not read or write jet-related branches.
// - doJet=true  or doJetMode=1: Require a D0-associated jet; apply jet pT/eta cuts to the associated jet; write jet/D0-jet branches.
// - doJetMode=2: Require at least one jet in the EVENT passing pT/eta cuts (no D0-jet association required); do not write D0-jet branches.
void nTupleMaker(const char *inputFile = "", TString pPbORPbp = "", int inputMC = 1, bool responseOnly = true, bool buildResponse = true, bool doJet = true, int doJetMode = -1, int systematicVariation = 0, bool useMultiplicityWeights = true, bool useJetPtWeights = false)
{
    //print all function arguments for debugging
    std::cout << "nTupleMaker called with arguments:" << std::endl;
    std::cout << "  inputFile: " << inputFile << std::endl;
    std::cout << "  pPbORPbp: " << pPbORPbp << std::endl;
    std::cout << "  inputMC: " << inputMC << std::endl;
    std::cout << "  responseOnly: " << responseOnly << std::endl;
    std::cout << "  buildResponse: " << buildResponse << std::endl;
    std::cout << "  doJet: " << doJet << std::boolalpha << doJet << std::endl;
    std::cout << "  doJetMode: " << doJetMode << std::endl;
    std::cout << "  systematicVariation: " << systematicVariation << std::endl;
    std::cout << "  useMultiplicityWeights: " << useMultiplicityWeights << std::endl;
    std::cout << "  useJetPtWeights: " << useJetPtWeights << std::endl;

    std::cout << "Starting D0 FF Analysis with minimal ntuple maker" << std::endl;
    TString multweightsfile = "/media/niviths/local/analysis_code/data_analysis/d0_FF/1_createTuple/outputs/2026-05-28/multiplicity_weights_Pbp/multiplicity_weights_Pbp_smoothed.root";
    // TString multweightsfile = "/media/niviths/local/analysis_code/data_analysis/d0_FF/1_createTuple/outputs/2026-05-28/multiplicity_weights_pPb/multiplicity_weights_pPb_smoothed.root";
    std::cout << "Using multiplicity weights file: " << multweightsfile << std::endl;
    // Derive a default jet-pt weights file path from the multiplicity weights file
    TString jetptweightsfile = multweightsfile;
    jetptweightsfile.ReplaceAll("_smoothed.root", "_jetpt_weights.root");
    std::cout << "Looking for jet-pt weights file: " << jetptweightsfile << std::endl;

    // Enforce mutual exclusion between weight modes
    if (useMultiplicityWeights && useJetPtWeights) {
        std::cerr << "ERROR: Both multiplicity and jet-pt weights requested. They are mutually exclusive. Aborting." << std::endl;
        return;
    }
    // Use default files if empty string provided
    TString fInputFileName = inputFile;

    // Build output file name robustly for both .root and list inputs
    TString fOutputFileName = fInputFileName;
    bool inputIsRoot = fInputFileName.EndsWith(".root");
    if (inputIsRoot) {
        fOutputFileName.ReplaceAll(".root", "_filtered.root");
    } else {
        fOutputFileName.ReplaceAll(".txt", "_filtered.root");
    }

    // Append systematic variation tag to output file name, e.g. _sysvar0, _sysvar1, ...
    if (systematicVariation >= 0) {
        fOutputFileName.ReplaceAll(".root", Form("_sysvar%d.root", systematicVariation));
    }

    std::cout << "Input file: " << fInputFileName << std::endl;
    std::cout << "Output file: " << fOutputFileName << std::endl;

    // Fast path: only build response matrix (skip filtering & FragmenNtuple production)
    if (responseOnly) {
        if (!inputMC) {
            std::cerr << "ERROR: responseOnly requested but inputMC=0 (need MC for response). Aborting." << std::endl;
            return;
        }
        TString respOut = fInputFileName;
        if (respOut.EndsWith(".root")) {
            respOut.ReplaceAll(".root", "_response.root");
        } else {
            respOut.ReplaceAll(".txt", "_response.root");
        }
        if (systematicVariation >= 0) {
            respOut.ReplaceAll(".root", Form("_sysvar%d.root", systematicVariation));
        }
        std::cout << "[responseOnly] Generating response matrix from original file -> " << respOut << std::endl;
        createResponseMatrix(fInputFileName, respOut, doJetMode, useMultiplicityWeights, multweightsfile, useJetPtWeights);
        std::cout << "[responseOnly] Done." << std::endl;
        return;
    }

    // Resolve jet mode (backward compatible with existing doJet flag)
    int jetMode = doJetMode;
    if (jetMode == -1)
    {
        jetMode = doJet ? 1 : 0;
    }

    // Prepare input as either a single TTree from TFile or a TChain from a file list
    TFile *inputRoot = nullptr;
    TChain *inputChain = nullptr;
    TTree *inputTree = nullptr;

    auto isListPath = [&](const TString &path) {
        return path.EndsWith(".txt");
    };

    if (isListPath(fInputFileName))
    {
        // Build a TChain from a list of ROOT files
        inputChain = new TChain("d0jets");
        std::ifstream listFile(fInputFileName.Data());
        if (!listFile.is_open()) {
            std::cerr << "ERROR: Cannot open input list file: " << fInputFileName << std::endl;
            return;
        }
        std::string line;
        int added = 0;
        while (std::getline(listFile, line)) {
            // Trim whitespace
            if (line.empty()) continue;
            // Skip comments
            if (line[0] == '#') continue;
            TString fpath(line.c_str());
            fpath = fpath.Strip(TString::kBoth, ' ');
            if (fpath.Length() == 0) continue;
            Long64_t nAdded = inputChain->Add(fpath.Data());
            if (nAdded > 0) {
                ++added;
            } else {
                std::cerr << "[WARN] Failed to add file to chain: " << fpath << std::endl;
            }
        }
        if (added == 0 || inputChain->GetListOfFiles()->GetEntries() == 0) {
            std::cerr << "ERROR: TChain is empty — no valid ROOT files were added from list: " << fInputFileName << std::endl;
            delete inputChain;
            return;
        }
        inputTree = inputChain; // TChain is-a TTree
    }
    else
    {
        // Single ROOT file path
        inputRoot = TFile::Open(fInputFileName, "READ");
        if (!inputRoot || inputRoot->IsZombie())
        {
            std::cerr << "ERROR: Cannot open input file: " << fInputFileName << std::endl;
            if (inputRoot) inputRoot->Close();
            return;
        }
        inputTree = (TTree *)inputRoot->Get("d0jets");
        if (!inputTree)
        {
            std::cerr << "ERROR: Cannot find 'd0jets' tree in input file" << std::endl;
            inputRoot->Close();
            return;
        }
    }



    //load efficiency maps for pions and kaons
    TH2D* effMapKaon = nullptr;

    TFile* effFile_kaon;
    if(pPbORPbp == "Pbp"){
        if(systematicVariation==0)
            effFile_kaon = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/8_pidcalib/lol/pidcalib_output_Ap_09_pEta_k/effhists-ApTurbo16-down-K-MC15TuneV1_ProbNNk>0.9&MC15TuneV1_ProbNNghost<0.3&TRCHI2NDOF<4-P.ETA-binning.root");
        else if(systematicVariation==1)
            effFile_kaon = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/8_pidcalib/pidcalib_sysvar/pidcalib_output_Ap_095_pEta_k/effhists-ApTurbo16-down-K-MC15TuneV1_ProbNNk>0.95&MC15TuneV1_ProbNNghost<0.2&TRCHI2NDOF<3-P.ETA.root");
        else if(systematicVariation==2)
            effFile_kaon = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/8_pidcalib/pidcalib_sysvar/pidcalib_output_Ap_085_pEta_k/effhists-ApTurbo16-down-K-MC15TuneV1_ProbNNk>0.85&MC15TuneV1_ProbNNghost<0.4&TRCHI2NDOF<5-P.ETA.root");
        else{
            std::cerr << "Warning: invalid systematic variation index for kaon PID efficiency map: " << systematicVariation << std::endl;
            return;
        }

    } else if (pPbORPbp == "pPb" || pPbORPbp == "pp") { //TODO fix
        if(systematicVariation==0)
            effFile_kaon = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/8_pidcalib/lol/pidcalib_output_pA_09_pEta_k/effhists-pATurbo16-down-K-MC15TuneV1_ProbNNk>0.9&MC15TuneV1_ProbNNghost<0.3&TRCHI2NDOF<4-P.ETA-binning.root");
        else if(systematicVariation==1)
            effFile_kaon = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/8_pidcalib/pidcalib_sysvar/pidcalib_output_pA_095_pEta_k/effhists-pATurbo16-down-K-MC15TuneV1_ProbNNk>0.95&MC15TuneV1_ProbNNghost<0.2&TRCHI2NDOF<3-P.ETA.root");
        else if(systematicVariation==2)
            effFile_kaon = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/8_pidcalib/pidcalib_sysvar/pidcalib_output_pA_085_pEta_k/effhists-pATurbo16-down-K-MC15TuneV1_ProbNNk>0.85&MC15TuneV1_ProbNNghost<0.4&TRCHI2NDOF<5-P.ETA.root");
        else{
            std::cerr << "Warning: invalid systematic variation index for kaon PID efficiency map: " << systematicVariation << std::endl;
            return;
        }
    } else {
        std::cerr << "no correction for pp for kaon PID efficiency map" << std::endl;
        return;
    }
    if (!effFile_kaon || effFile_kaon->IsZombie())
    {
        std::cerr << "Error: Could not open kaon pid efficiency map file" << std::endl;
        return;
    }
    if(systematicVariation==0){
        effMapKaon = dynamic_cast<TH2D*>(effFile_kaon->Get("eff_MC15TuneV1_ProbNNk>0.9&MC15TuneV1_ProbNNghost<0.3&TRCHI2NDOF<4"));
    } else if(systematicVariation==1){
        effMapKaon = dynamic_cast<TH2D*>(effFile_kaon->Get("eff_MC15TuneV1_ProbNNk>0.95&MC15TuneV1_ProbNNghost<0.2&TRCHI2NDOF<3"));
    } else if(systematicVariation==2){
        effMapKaon = dynamic_cast<TH2D*>(effFile_kaon->Get("eff_MC15TuneV1_ProbNNk>0.85&MC15TuneV1_ProbNNghost<0.4&TRCHI2NDOF<5"));
    }

    TH2D* effMapPion = nullptr;
    TFile* effFile_Pion;
    if(pPbORPbp == "Pbp"){
        if(systematicVariation==0)
            effFile_Pion = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/8_pidcalib/lol/pidcalib_output_Ap_09_pEta_pi/effhists-ApTurbo16-down-Pi-MC15TuneV1_ProbNNpi>0.9&MC15TuneV1_ProbNNghost<0.3&TRCHI2NDOF<4-P.ETA-binning.root");
        else if(systematicVariation==1)
            effFile_Pion = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/8_pidcalib/pidcalib_sysvar/pidcalib_output_Ap_095_pEta_pi/effhists-ApTurbo16-down-Pi-MC15TuneV1_ProbNNpi>0.95&MC15TuneV1_ProbNNghost<0.2&TRCHI2NDOF<3-P.ETA.root");
        else if(systematicVariation==2)
            effFile_Pion = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/8_pidcalib/pidcalib_sysvar/pidcalib_output_Ap_085_pEta_pi/effhists-ApTurbo16-down-Pi-MC15TuneV1_ProbNNpi>0.85&MC15TuneV1_ProbNNghost<0.4&TRCHI2NDOF<5-P.ETA.root");
        else{
            std::cerr << "Warning: invalid systematic variation index for pion PID efficiency map: " << systematicVariation << std::endl;
            return;
        }
    } else if(pPbORPbp == "pPb" || pPbORPbp == "pp") { //TODO fix
        if(systematicVariation==0)
            effFile_Pion = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/8_pidcalib/lol/pidcalib_output_pA_09_pEta_pi/effhists-pATurbo16-down-Pi-MC15TuneV1_ProbNNpi>0.9&MC15TuneV1_ProbNNghost<0.3&TRCHI2NDOF<4-P.ETA-binning.root");
        else if(systematicVariation==1)
            effFile_Pion = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/8_pidcalib/pidcalib_sysvar/pidcalib_output_pA_095_pEta_pi/effhists-pATurbo16-down-Pi-MC15TuneV1_ProbNNpi>0.95&MC15TuneV1_ProbNNghost<0.2&TRCHI2NDOF<3-P.ETA.root");
        else if(systematicVariation==2)
            effFile_Pion = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/8_pidcalib/pidcalib_sysvar/pidcalib_output_pA_085_pEta_pi/effhists-pATurbo16-down-Pi-MC15TuneV1_ProbNNpi>0.85&MC15TuneV1_ProbNNghost<0.4&TRCHI2NDOF<5-P.ETA.root");
        else{
            std::cerr << "Warning: invalid systematic variation index for pion PID efficiency map: " << systematicVariation << std::endl;
            return;
        }
    } else {
        std::cerr << "no correction for pp for pion PID efficiency map" << std::endl;
        return;
    }

    if(systematicVariation==0){
        effMapPion = dynamic_cast<TH2D*>(effFile_Pion->Get("eff_MC15TuneV1_ProbNNpi>0.9&MC15TuneV1_ProbNNghost<0.3&TRCHI2NDOF<4"));
    } else if(systematicVariation==1){
        effMapPion = dynamic_cast<TH2D*>(effFile_Pion->Get("eff_MC15TuneV1_ProbNNpi>0.95&MC15TuneV1_ProbNNghost<0.2&TRCHI2NDOF<3"));
    } else if(systematicVariation==2){
        effMapPion = dynamic_cast<TH2D*>(effFile_Pion->Get("eff_MC15TuneV1_ProbNNpi>0.85&MC15TuneV1_ProbNNghost<0.4&TRCHI2NDOF<5"));
    }

    if (!effMapKaon || !effMapPion)
    {
        std::cerr << "Error: Could not find efficiency maps in file" << std::endl;
        effFile_kaon->Close();
        effFile_Pion->Close();
        return;
    }

    TFile* effFile_Reconstruction;
    if(pPbORPbp == "Pbp"){
        effFile_Reconstruction = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/7_efficiency/output_reco_standalone_15_16_Pbp_rerun_2026-09-08/output_reco_standalone_15_16_Pbp_rerun.root");
        // effFile_Reconstruction = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/7_efficiency/output_reco_standalone_15_16_Pbp_2026-05-26/output_reco_standalone_15_16_Pbp.root");
    } else if(pPbORPbp == "pPb"){
        effFile_Reconstruction = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/7_efficiency/output_reco_standalone_11_12_pPb_rerun_2026-09-08/output_reco_standalone_11_12_pPb_rerun.root");
        // effFile_Reconstruction = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/7_efficiency/output_reco_standalone_11_12_pPb_2026-05-26/output_reco_standalone_11_12_pPb.root");
    } else if(pPbORPbp == "pp"){
        effFile_Reconstruction = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/7_efficiency/output_reco_standalone_full_pp_2026-04-30/output_reco_standalone_full_pp.root");
    } else {
        std::cerr << "no correction for pp for reconstruction efficiency map" << std::endl;
        return;
    }
    if(!effFile_Reconstruction || effFile_Reconstruction->IsZombie())
    {
        std::cerr << "Error: Could not open reconstruction efficiency file" << std::endl;
        return;
    }

    TH2F* effMapReco = dynamic_cast<TH2F*>(effFile_Reconstruction->Get("reco_efficiency"));
    if(!effMapReco)
    {
        std::cerr << "Error: Could not find reconstruction efficiency map in file" << std::endl;
        effFile_Reconstruction->Close();
        return;
    }

    TFile* acceptanceFile;
    if(pPbORPbp == "Pbp")
        acceptanceFile = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/9_acceptance/D0AcceptanceMap_Pbp_74_75_2026-04-08/D0AcceptanceMap_Pbp_74_75.root");
    else if (pPbORPbp == "pPb")
        acceptanceFile = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/9_acceptance/D0AcceptanceMap_pPb_54plus73_2026-03-31/D0AcceptanceMap_pPb_54plus73.root");
    else if (pPbORPbp == "pp")
        acceptanceFile = TFile::Open("/media/niviths/local/analysis_code/data_analysis/d0_FF/9_acceptance/D0AcceptanceMap_pp_1_full_2026-04-30/D0AcceptanceMap_pp_1_full.root");
    else {
        std::cerr << "no correction for pp for acceptance map" << std::endl;
        return;
    }
    if (!acceptanceFile || acceptanceFile->IsZombie())
    {
        std::cerr << "Error: Could not open acceptance map file" << std::endl;
        return;
    }
    TH2D* acceptanceMap = dynamic_cast<TH2D*>(acceptanceFile->Get("hAcceptance"));
    if(!acceptanceMap)
    {
        std::cerr << "Error: Could not find acceptance map in file" << std::endl;
        acceptanceFile->Close();
        return;
    }


    // Create output file
    TFile *outputRoot = TFile::Open(fOutputFileName, "RECREATE");
    if (!outputRoot || outputRoot->IsZombie())
    {
        std::cerr << "ERROR: Cannot create output file: " << fOutputFileName << std::endl;
        if (inputRoot) inputRoot->Close();
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
    std::vector<float> *jet_rapidity = nullptr;
    std::vector<float> *jet_phi = nullptr;
    std::vector<float> *jet_mass = nullptr;
    std::vector<int> *jet_n_const = nullptr;
    std::vector<int> *jet_n_charged = nullptr;
    std::vector<int> *jet_n_neutral = nullptr;
    std::vector<int> *jet_n_d0 = nullptr;

    // D0 info vectors
    std::vector<float> *d0_pt = nullptr;
    std::vector<float> *d0_pz = nullptr;
    std::vector<float> *d0_e = nullptr;
    std::vector<float> *d0_eta = nullptr;
    std::vector<float> *d0_rapidity = nullptr;
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

    const bool hasJetRapidityBranch = inputTree->GetBranch("jet_rapidity");
    const bool hasD0RapidityBranch = inputTree->GetBranch("d0_rapidity");

    // Set branch addresses for input tree
    inputTree->SetBranchAddress("evt_num", &evt_num);
    inputTree->SetBranchAddress("run_num", &run_num);
    inputTree->SetBranchAddress("n_pvs", &n_pvs);
    // Try to read event multiplicity branch if present
    int event_multiplicity = 0;
    bool has_event_multiplicity = false;
    if (inputTree->GetBranch("event_multiplicity")) {
        inputTree->SetBranchAddress("event_multiplicity", &event_multiplicity);
        has_event_multiplicity = true;
    }

    // Read jet collections when any jet requirement is enabled (mode 1 or 2)
    if(jetMode != 0)
    {
        inputTree->SetBranchAddress("jet_pt", &jet_pt);
        inputTree->SetBranchAddress("jet_eta", &jet_eta);
        if (hasJetRapidityBranch)
            inputTree->SetBranchAddress("jet_rapidity", &jet_rapidity);
        inputTree->SetBranchAddress("jet_phi", &jet_phi);
        inputTree->SetBranchAddress("jet_mass", &jet_mass);
        inputTree->SetBranchAddress("jet_n_const", &jet_n_const);
        inputTree->SetBranchAddress("jet_n_charged", &jet_n_charged);
        inputTree->SetBranchAddress("jet_n_neutral", &jet_n_neutral);
        inputTree->SetBranchAddress("jet_n_d0", &jet_n_d0);
    }
    inputTree->SetBranchAddress("d0_pt", &d0_pt);
    inputTree->SetBranchAddress("d0_pz", &d0_pz);
    inputTree->SetBranchAddress("d0_e", &d0_e);
    inputTree->SetBranchAddress("d0_eta", &d0_eta);
    if (hasD0RapidityBranch)
        inputTree->SetBranchAddress("d0_rapidity", &d0_rapidity);
    inputTree->SetBranchAddress("d0_phi", &d0_phi);
    inputTree->SetBranchAddress("d0_mass", &d0_mass);
    inputTree->SetBranchAddress("d0_vtx_chi2", &d0_vtx_chi2);
    inputTree->SetBranchAddress("d0_ip", &d0_ip);
    inputTree->SetBranchAddress("d0_ip_chi2", &d0_ip_chi2);
    inputTree->SetBranchAddress("d0_fd", &d0_fd);
    inputTree->SetBranchAddress("d0_fd_chi2", &d0_fd_chi2);
    inputTree->SetBranchAddress("d0_DOCA", &d0_DOCA);
    // Read D0-jet association branches only when we require D0-associated jets (mode 1)
    if(jetMode == 1)
    {
        inputTree->SetBranchAddress("d0_z", &d0_z);
        inputTree->SetBranchAddress("d0_jet_idx", &d0_jet_idx);
        inputTree->SetBranchAddress("d0_in_jet", &d0_in_jet);
        inputTree->SetBranchAddress("d0_jet_dr", &d0_jet_dr);
    }
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

    auto getJetRapidity = [&](size_t idx) -> float
    {
        if (hasJetRapidityBranch && jet_rapidity && idx < jet_rapidity->size())
            return (*jet_rapidity)[idx];
        if (jet_pt && jet_eta && jet_mass && idx < jet_pt->size() && idx < jet_eta->size() && idx < jet_mass->size())
            return calcRapidityFromPtEtaMass((*jet_pt)[idx], (*jet_eta)[idx], (*jet_mass)[idx]);
        if (jet_eta && idx < jet_eta->size())
            return (*jet_eta)[idx];
        return -999.f;
    };

    auto getD0Rapidity = [&](size_t idx) -> float
    {
        if (hasD0RapidityBranch && d0_rapidity && idx < d0_rapidity->size())
            return (*d0_rapidity)[idx];
        if (d0_e && d0_pz && idx < d0_e->size() && idx < d0_pz->size())
            return calcRapidityFromEnergyPz((*d0_e)[idx], (*d0_pz)[idx]);
        return -999.f;
    };

    // Variables for output tree
    float v_tagdR = 0;
    float v_tagMass = 0;
    float v_tagPt = 0;
    float v_tagEta = 0;
    float v_tagY = 0;
    float v_tag_idx_jet = 0;
    float v_tag_decVtxChi2 = 0;
    float v_tag_logdecVtxChi2 = 0;
    float v_jetPt = 0;
    float v_jetEta = 0;
    float v_jetY = 0;
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
    float v_kaon_efficiency = 0;
    float v_pion_efficiency = 0;
    float v_reconstruction_efficiency = 0;
    float v_combined_PID_efficiency = 0; // New variable for combined efficiency
    float v_combined_efficiency = 0; // New variable for combined efficiency
    float v_combined_eff_and_acceptance = 0; // New variable for combined efficiency and acceptance
    float v_acceptance = 0; // New variable for acceptance
    // Event-level jet info for mode 2 (no D0 association)
    int   v_evtNJets = 0;
    float v_evtLeadJetPt = 0;
    float v_evtLeadJetEta = 0;
    float v_evtLeadJetY = 0;
    float v_evtLeadJetPhi = 0;
    float v_evtLeadJetNConst = 0;

    // Create branches for output tree
    outputTree->Branch("tagMass", &v_tagMass, "tagMass/F");
    outputTree->Branch("tagPt", &v_tagPt, "tagPt/F");
    outputTree->Branch("tagEta", &v_tagEta, "tagEta/F");
    outputTree->Branch("tagY", &v_tagY, "tagY/F");
    outputTree->Branch("tag_ip_chi2", &v_tag_decVtxChi2, "tag_ip_chi2/F");
    outputTree->Branch("log_tag_ipchi2", &v_tag_logdecVtxChi2, "log_tag_ipchi2/F");
    if(jetMode == 1)
    {
        outputTree->Branch("tagJetdR", &v_tagdR, "tagJetdR/F");
        outputTree->Branch("tagidxjet", &v_tag_idx_jet, "tagidxjet/F");
        outputTree->Branch("jetPt", &v_jetPt, "jetPt/F");
        outputTree->Branch("jetEta", &v_jetEta, "jetEta/F");
        outputTree->Branch("jetY", &v_jetY, "jetY/F");
        outputTree->Branch("jetnConst", &v_jetnConst, "jetnConst/F");
        outputTree->Branch("tagZ", &v_tagZ, "tagZ/F");
    }
    else if (jetMode == 2)
    {
        outputTree->Branch("evtNJets", &v_evtNJets, "evtNJets/I");
        outputTree->Branch("evtLeadJetPt", &v_evtLeadJetPt, "evtLeadJetPt/F");
        outputTree->Branch("evtLeadJetEta", &v_evtLeadJetEta, "evtLeadJetEta/F");
        outputTree->Branch("evtLeadJetY", &v_evtLeadJetY, "evtLeadJetY/F");
        outputTree->Branch("evtLeadJetPhi", &v_evtLeadJetPhi, "evtLeadJetPhi/F");
        outputTree->Branch("evtLeadJetNConst", &v_evtLeadJetNConst, "evtLeadJetNConst/F");
    }
    outputTree->Branch("piPprobNNpi", &v_piPprobNNpi, "piPprobNN/F");
    outputTree->Branch("piPprobGhost", &v_piPprobGhost, "piPprobGhost/F");
    outputTree->Branch("piPTrckChi2", &v_piPTrckChi2, "piPTrckChi2/F");
    outputTree->Branch("KprobNNK", &v_KprobNNK, "KprobNNK/F");
    outputTree->Branch("KprobGhost", &v_KprobGhost, "KprobGhost/F");
    outputTree->Branch("KTrckChi2", &v_KTrckChi2, "KTrckChi2/F");
    outputTree->Branch("decayVtxChi2", &v_decayVtxChi2, "decayVtxChi2/F");
    outputTree->Branch("Distance1", &v_Dist1, "Distance1/F");
    outputTree->Branch("isPrimary", &v_isPrimary, "isPrimary/F");
    outputTree->Branch("kaon_efficiency", &v_kaon_efficiency, "kaon_efficiency/F");
    outputTree->Branch("pion_efficiency", &v_pion_efficiency, "pion_efficiency/F");
    outputTree->Branch("reconstruction_efficiency", &v_reconstruction_efficiency, "reconstruction_efficiency/F");
    outputTree->Branch("combined_PID_efficiency", &v_combined_PID_efficiency, "combined_PID_efficiency/F");
    outputTree->Branch("combined_efficiency", &v_combined_efficiency, "combined_efficiency/F");
    outputTree->Branch("combined_eff_and_acceptance", &v_combined_eff_and_acceptance, "combined_eff_and_acceptance/F");
    outputTree->Branch("acceptance", &v_acceptance, "acceptance/F");

    // Additional MC branches if needed
    if (inputMC)
    {
        outputTree->Branch("isPrimary", &v_isPrimary, "isPrimary/F");
    }

    // Load multiplicity weights (if available). This allows reweighting pp MC to match pPb.
    std::map<int, double> multiplicityWeights;
    // Load jet-pt weights histogram if present
    TH1 *h_jet_weights = nullptr;
    TFile *jwfile = nullptr;
    if (inputMC && pPbORPbp == "pp") {
        if (useMultiplicityWeights) {
            TFile *wfile = TFile::Open(multweightsfile, "READ");
            if (wfile && !wfile->IsZombie()) {
                TTree *wtree = (TTree*)wfile->Get("multiplicity_weights");
                if (wtree) {
                    Int_t mult = 0; Double_t w = 1.0;
                    wtree->SetBranchAddress("multiplicity", &mult);
                    wtree->SetBranchAddress("weight", &w);
                    for (Long64_t iw = 0; iw < wtree->GetEntries(); ++iw) {
                        wtree->GetEntry(iw);
                        multiplicityWeights[(int)mult] = (double)w;
                    }
                }
                wfile->Close();
            }
        }
        // Load jet-pt weights only when requested
        if (useJetPtWeights) {
            jwfile = TFile::Open(jetptweightsfile, "READ");
            if (jwfile && !jwfile->IsZombie()) {
                h_jet_weights = dynamic_cast<TH1*>(jwfile->Get("h_jet_weights"));
                if (h_jet_weights) {
                    h_jet_weights->SetDirectory(0); // detach from file so we can close it
                    std::cout << "Loaded jet-pt weights histogram from: " << jetptweightsfile << std::endl;
                } else {
                    std::cerr << "[WARN] Could not find 'h_jet_weights' in " << jetptweightsfile << std::endl;
                }
                jwfile->Close();
            }
        }
    }

    // Add event weight branch to output tree (default 1.0)
    float v_event_weight = 1.0f;
    outputTree->Branch("event_weight", &v_event_weight, "event_weight/F");

    // Process events
    Long64_t nEntries = inputTree->GetEntries();
    std::cout << "Processing " << nEntries << " events for filtered output" << std::endl;

    int events_processed = 0;
    int d0_accepted = 0;


    Long64_t updateInterval = nEntries / 200; // ~200 updates (0.5% increments)
    if (updateInterval < 1) updateInterval = 1;

    for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++)
    {
        inputTree->GetEntry(iEntry);

        if (iEntry % updateInterval == 0 || iEntry + 1 == nEntries)
        {
            printProgress(iEntry, nEntries);
        }

        events_processed++;

        if(n_pvs > 1) // Only process events with 1 primary vertex
        //TODO VALIDATE THIS IS CORRECT
        //TODO VALIDATE THIS IS CORRECT
        //TODO VALIDATE THIS IS CORRECT
        {
            continue; // Skip events with more than 1 primary vertex
        }
        // Jet mode 2: require at least one jet in the event that satisfies pT/eta cuts
        // and compute event-level jet info (leading jet + counts)
        int evtNJets_pass = 0;
        int evtLeadIdx = -1;
        float evtLeadPt = -1.f;
        float evtLeadEta = 0.f;
        float evtLeadY = 0.f;
        float evtLeadPhi = 0.f;
        float evtLeadNConst = 0.f;
        if (jetMode == 2)
        {
            if (jet_pt && jet_eta)
            {
                for (size_t ij = 0; ij < jet_pt->size(); ++ij)
                {
                    float jetY = getJetRapidity(ij);
                    if ((*jet_pt)[ij] >= 5.0 && jetY >= 2.5 && jetY <= 4.0)
                    {
                        evtNJets_pass++;
                        if ((*jet_pt)[ij] > evtLeadPt)
                        {
                            evtLeadPt = (*jet_pt)[ij];
                            evtLeadIdx = static_cast<int>(ij);
                        }
                    }
                }
            }
            if (evtNJets_pass == 0)
            {
                continue; // Skip events without a qualifying jet
            }
            if (evtLeadIdx >= 0)
            {
                evtLeadEta = (jet_eta && evtLeadIdx < (int)jet_eta->size()) ? (*jet_eta)[evtLeadIdx] : 0.f;
                evtLeadY = getJetRapidity(evtLeadIdx);
                evtLeadPhi = (jet_phi && evtLeadIdx < (int)jet_phi->size()) ? (*jet_phi)[evtLeadIdx] : 0.f;
                evtLeadNConst = (jet_n_const && evtLeadIdx < (int)jet_n_const->size()) ? (*jet_n_const)[evtLeadIdx] : 0.f;
            }
        }


        // Process each D0 in the event
        for (size_t i_d0 = 0; i_d0 < d0_pt->size(); i_d0++)
        {
            // Basic D0 selection cuts
            if ((*d0_pt)[i_d0] < 1.0) //TODO is change from 2 removing possible bias?
                continue; // Minimum pT
            // if ((*d0_eta)[i_d0] < 2.0 || (*d0_eta)[i_d0] > 4.5)
            //     continue; // Eta acceptance
            float d0Y = getD0Rapidity(i_d0);
            if (d0Y < 2.0 || d0Y > 4.5)
                continue; // Rapidity acceptance

            // D0 mass window cut
            if (std::abs((*d0_mass)[i_d0] - 1.865) > 0.07)
                continue; // Mass window: ±50 MeV around nominal D0 mass

            // Check if D0 is associated with a jet
            int jet_idx = -1;
            if(jetMode == 1){
                jet_idx = (*d0_jet_idx)[i_d0];
                if (jet_idx < 0)
                    continue; // Require associated jet

                // Jet selection cuts for the associated jet
                if ((*jet_pt)[jet_idx] < 5.0)
                    continue; // Minimum jet pT
                float jetY = getJetRapidity(jet_idx);
                if (jetY < 2.5 || jetY > 4.0)
                    continue; // Jet rapidity range
            }
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
                if ((*dau_d0_idx)[i_dau] != static_cast<int>(i_d0))
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
            if(systematicVariation==0){ //default
                if (kaon_pnn_k < 0.9)
                    continue; // Loose kaon ID
                if (pion_pnn_pi < 0.9)
                    continue; // Loose pion ID
                if (kaon_ghost_prob > 0.3 || pion_ghost_prob > 0.3)
                    continue; // Ghost probability cut
                if (kaon_chi2 > 4.0 || pion_chi2 > 4.0)
                    continue; // Track chi2 cut
            } else if (systematicVariation == 1) {
                // Tighter PID cuts for systematic variation
                if (kaon_pnn_k < 0.95)
                    continue; // Tighter kaon ID
                if (pion_pnn_pi < 0.95)
                    continue; // Tighter pion ID
                if (kaon_ghost_prob > 0.2 || pion_ghost_prob > 0.2)
                    continue; // Tighter ghost probability cut
                if (kaon_chi2 > 3.0 || pion_chi2 > 3.0)
                    continue; // Tighter track chi2 cut
            } else if (systematicVariation == 2) {
                // Looser PID cuts for systematic variation
                if (kaon_pnn_k < 0.85)
                    continue; // Looser kaon ID
                if (pion_pnn_pi < 0.85)
                    continue; // Looser pion ID
                if (kaon_ghost_prob > 0.4 || pion_ghost_prob > 0.4)
                    continue; // Looser ghost probability cut
                if (kaon_chi2 > 5.0 || pion_chi2 > 5.0)
                    continue; // Looser track chi2 cut
            }
            // get efficiencies for kaon and pion
            float kaon_efficiency = effMapKaon->GetBinContent(
                effMapKaon->GetXaxis()->FindBin(kaon_p*1000),
                effMapKaon->GetYaxis()->FindBin(kaon_eta));
                if(kaon_efficiency>1) kaon_efficiency = 1;
            float pion_efficiency = effMapPion->GetBinContent(
                effMapPion->GetXaxis()->FindBin(pion_p*1000),
                effMapPion->GetYaxis()->FindBin(pion_eta));
                if(pion_efficiency>1) pion_efficiency = 1;
            float combined_PID_efficiency = kaon_efficiency * pion_efficiency;

            float reconstruction_efficiency = effMapReco->GetBinContent(
                effMapReco->GetXaxis()->FindBin((*d0_pt)[i_d0]),
                effMapReco->GetYaxis()->FindBin((*d0_eta)[i_d0]));
            float combined_efficiency = combined_PID_efficiency * reconstruction_efficiency;

            float acceptance = acceptanceMap->GetBinContent(
                acceptanceMap->GetXaxis()->FindBin((*d0_pt)[i_d0]),
                acceptanceMap->GetYaxis()->FindBin((*d0_eta)[i_d0]));
            float combined_eff_and_acceptance = combined_efficiency * acceptance;
            // Fill output variables
            v_tagMass = (*d0_mass)[i_d0];
            v_tagPt = (*d0_pt)[i_d0];
            v_tagEta = (*d0_eta)[i_d0];
            v_tagY = d0Y;
            v_tag_logdecVtxChi2 = (*d0_ip_chi2)[i_d0] > 0 ? std::log10((*d0_ip_chi2)[i_d0]) : -999;
            v_tag_decVtxChi2 = (*d0_ip_chi2)[i_d0];
            if(jetMode == 1){
                v_tagdR = (*d0_jet_dr)[i_d0];
                v_tag_idx_jet = jet_idx;
                v_jetPt = (*jet_pt)[jet_idx];
                v_jetEta = (*jet_eta)[jet_idx];
                v_jetY = getJetRapidity(jet_idx);
                v_jetnConst = (*jet_n_const)[jet_idx];
                if((*d0_in_jet)[i_d0]==1 && (*d0_z)[i_d0] > 1) // Sanity check for z variable, which should be between 0 and 1 for D0's inside jets. If it fails, set to -1 to indicate an issue.
                    v_tagZ = (*d0_z)[i_d0]/1000;
                else
                    v_tagZ = (*d0_z)[i_d0];
            }
            else if (jetMode == 2)
            {
                // Save event-level jet info (same for all D0 in the event)
                v_evtNJets = evtNJets_pass;
                v_evtLeadJetPt = evtLeadPt > 0 ? evtLeadPt : 0.f;
                v_evtLeadJetEta = evtLeadEta;
                v_evtLeadJetY = evtLeadY;
                v_evtLeadJetPhi = evtLeadPhi;
                v_evtLeadJetNConst = evtLeadNConst;
            }
            v_KprobNNK = kaon_pnn_k;
            v_KprobGhost = kaon_ghost_prob;
            v_KTrckChi2 = kaon_chi2;
            v_piPprobNNpi = pion_pnn_pi;
            v_piPprobGhost = pion_ghost_prob;
            v_piPTrckChi2 = pion_chi2;
            v_decayVtxChi2 = (*d0_vtx_chi2)[i_d0];
            v_Dist1 = (*d0_DOCA)[i_d0];
            v_kaon_efficiency = kaon_efficiency;
            v_pion_efficiency = pion_efficiency;
            v_reconstruction_efficiency = reconstruction_efficiency;
            v_combined_PID_efficiency = combined_PID_efficiency;
            v_combined_efficiency = combined_efficiency;
            v_combined_eff_and_acceptance = combined_eff_and_acceptance;
            v_acceptance = acceptance;

            // MC specific variables
            if (inputMC && i_d0 < mc_d0_origin->size())
            {
                v_isPrimary = (*mc_d0_origin)[i_d0] == 1 ? 1 : 0; // Is it a prompt D0?
            }
            else
            {
                v_isPrimary = -1; // Unknown for data
            }

            // // Calculate rapidity
            // if ((*d0_e)[i_d0] > (*d0_pz)[i_d0]) {
            //     v_tagY = 0.5 * log(((*d0_e)[i_d0] + (*d0_pz)[i_d0]) / ((*d0_e)[i_d0] - (*d0_pz)[i_d0]));
            // } else {
            //     v_tagY = -999; // Assign an invalid value if rapidity cannot be calculated
            // }
            // Fill event weight only for pp MC when weights were loaded
            // Decide which weight to apply (mutually exclusive)
            v_event_weight = 1.0f;
            if (inputMC && pPbORPbp == "pp") {
                if (useMultiplicityWeights && has_event_multiplicity) {
                    auto it = multiplicityWeights.find(event_multiplicity);
                    if (it != multiplicityWeights.end()) v_event_weight = static_cast<float>(it->second);
                    else v_event_weight = 1.0f;
                } else if (useJetPtWeights) {
                    if (h_jet_weights && jetMode == 1 && jet_idx >= 0 && jet_pt && jet_idx < (int)jet_pt->size()) {
                        double jpt = static_cast<double>((*jet_pt)[jet_idx]);
                        int bin = h_jet_weights->FindBin(jpt);
                        double wj = h_jet_weights->GetBinContent(bin);
                        if (wj > 0) v_event_weight = static_cast<float>(wj);
                    }
                }
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
    if (inputRoot) inputRoot->Close();
    if (inputChain) { delete inputChain; inputChain = nullptr; }

    std::cout << "Output saved to: " << fOutputFileName << std::endl;

    // Create response matrix tree if requested and processing MC
        if (inputMC && buildResponse && (jetMode == 1))
        {
            std::cout << "Creating response matrix tree..." << std::endl;
            // Construct response filename from the original input and append systematic tag
            TString respName = fInputFileName;
            if (respName.EndsWith(".root")) {
                respName.ReplaceAll(".root", "_response.root");
            } else {
                respName.ReplaceAll(".txt", "_response.root");
            }
            if (systematicVariation >= 0) {
                respName.ReplaceAll(".root", Form("_sysvar%d.root", systematicVariation));
            }
            std::cout << "useMultiplicityWeights for response matrix: " << (useMultiplicityWeights ? "Yes" : "No") << std::endl;
            std::cout << "useJetPtWeights for response matrix: " << (useJetPtWeights ? "Yes" : "No") << std::endl;
            createResponseMatrix(fInputFileName, respName, jetMode, useMultiplicityWeights, multweightsfile, useJetPtWeights);
        } else if (inputMC && !buildResponse) {
        std::cout << "[INFO] Skipping response matrix creation (buildResponse=false)." << std::endl;
    }
}

// Function to create response matrix
void createResponseMatrix(const char *inputFile, const char *outputFile, int jetMode, bool useWeights, const char* multweightsfile, bool useJetPtWeights)
{
    // Accept either a single ROOT file or a text file list to build a TChain
    TFile *inputRoot = nullptr;
    TChain *inputChain = nullptr;
    TTree *inputTree = nullptr;
    std::cout << "Creating response matrix from input: " << inputFile << std::endl;
    std::cout << "Output response matrix file: " << outputFile << std::endl;
    std::cout << "Jet mode for response matrix: " << jetMode << std::endl;
    std::cout << "Use multiplicity weights in response matrix: " << (useWeights ? "Yes" : "No") << std::endl;
    std::cout << "Use jet-pt weights in response matrix: " << (useJetPtWeights ? "Yes" : "No") << std::endl;
    // useWeights = false; // Force disable weights for response matrix, as they should not be applied to the response itself

    TString inPath(inputFile);
    auto isListPath = [&](const TString &path) {
        return path.EndsWith(".txt");
    };

    if (isListPath(inPath))
    {
        inputChain = new TChain("d0jets");
        std::ifstream listFile(inPath.Data());
        if (!listFile.is_open()) {
            std::cerr << "ERROR: Cannot open input list file for response matrix: " << inputFile << std::endl;
            return;
        }
        std::string line;
        int added = 0;
        while (std::getline(listFile, line)) {
            if (line.empty()) continue;
            if (line[0] == '#') continue;
            TString fpath(line.c_str());
            fpath = fpath.Strip(TString::kBoth, ' ');
            if (fpath.Length() == 0) continue;
            Long64_t nAdded = inputChain->Add(fpath.Data());
            if (nAdded > 0) ++added; else std::cerr << "[WARN] Failed to add file to chain (response): " << fpath << std::endl;
        }
        if (added == 0 || inputChain->GetListOfFiles()->GetEntries() == 0) {
            std::cerr << "ERROR: TChain is empty in response matrix — no valid files added from list: " << inputFile << std::endl;
            delete inputChain;
            return;
        }
        inputTree = inputChain;
    }
    else
    {
        inputRoot = TFile::Open(inputFile, "READ");
        if (!inputRoot || inputRoot->IsZombie())
        {
            std::cerr << "ERROR: Cannot open input file for response matrix: " << inputFile << std::endl;
            if (inputRoot) inputRoot->Close();
            return;
        }
        inputTree = (TTree *)inputRoot->Get("d0jets");
        if (!inputTree)
        {
            std::cerr << "ERROR: Cannot find 'd0jets' tree in input file" << std::endl;
            inputRoot->Close();
            return;
        }
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
    float r_jet_rapidity_det = 0;
    float r_jet_phi_det = 0;
    float r_jet_nconst_det = 0;
    float r_d0_pt_det = 0;
    float r_d0_eta_det = 0;
    float r_d0_rapidity_det = 0;
    float r_d0_y_det = 0;
    float r_d0_phi_det = 0;
    float r_d0_mass_det = 0;
    float r_d0_z_det = 0;

    float r_jet_pt_mc = 0;
    float r_jet_eta_mc = 0;
    float r_jet_rapidity_mc = 0;
    float r_jet_phi_mc = 0;
    float r_jet_nconst_mc = 0;
    float r_d0_pt_mc = 0;
    float r_d0_eta_mc = 0;
    float r_d0_rapidity_mc = 0;
    float r_d0_y_mc = 0;
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
    responseTree->Branch("jet_rapidity_det", &r_jet_rapidity_det, "jet_rapidity_det/F");
    responseTree->Branch("jet_phi_det", &r_jet_phi_det, "jet_phi_det/F");
    responseTree->Branch("jet_nconst_det", &r_jet_nconst_det, "jet_nconst_det/F");
    responseTree->Branch("d0_pt_det", &r_d0_pt_det, "d0_pt_det/F");
    responseTree->Branch("d0_eta_det", &r_d0_eta_det, "d0_eta_det/F");
    responseTree->Branch("d0_rapidity_det", &r_d0_rapidity_det, "d0_rapidity_det/F");
    responseTree->Branch("d0_y_det", &r_d0_y_det, "d0_y_det/F");
    responseTree->Branch("d0_phi_det", &r_d0_phi_det, "d0_phi_det/F");
    responseTree->Branch("d0_mass_det", &r_d0_mass_det, "d0_mass_det/F");
    responseTree->Branch("d0_z_det", &r_d0_z_det, "d0_z_det/F");

    responseTree->Branch("jet_pt_mc", &r_jet_pt_mc, "jet_pt_mc/F");
    responseTree->Branch("jet_eta_mc", &r_jet_eta_mc, "jet_eta_mc/F");
    responseTree->Branch("jet_rapidity_mc", &r_jet_rapidity_mc, "jet_rapidity_mc/F");
    responseTree->Branch("jet_phi_mc", &r_jet_phi_mc, "jet_phi_mc/F");
    responseTree->Branch("jet_nconst_mc", &r_jet_nconst_mc, "jet_nconst_mc/F");
    responseTree->Branch("d0_pt_mc", &r_d0_pt_mc, "d0_pt_mc/F");
    responseTree->Branch("d0_eta_mc", &r_d0_eta_mc, "d0_eta_mc/F");
    responseTree->Branch("d0_rapidity_mc", &r_d0_rapidity_mc, "d0_rapidity_mc/F");
    responseTree->Branch("d0_y_mc", &r_d0_y_mc, "d0_y_mc/F");
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
    std::vector<float> *jet_rapidity = nullptr;
    std::vector<float> *jet_phi = nullptr;
    std::vector<float> *jet_mass = nullptr;
    std::vector<int> *jet_n_const = nullptr;
    std::vector<int> *jet_n_d0 = nullptr;

    std::vector<float> *d0_pt = nullptr;
    std::vector<float> *d0_eta = nullptr;
    std::vector<float> *d0_rapidity = nullptr;
    std::vector<float> *d0_e = nullptr;
    std::vector<float> *d0_pz = nullptr;
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
    std::vector<float> *mc_d0_rapidity = nullptr;
    std::vector<float> *mc_d0_e = nullptr;
    std::vector<float> *mc_d0_pz = nullptr;
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
    // New: direct kinematic branches (preferred)
    std::vector<float> *mc_jet_pt = nullptr;
    std::vector<float> *mc_jet_eta = nullptr;
    std::vector<float> *mc_jet_rapidity = nullptr;
    std::vector<float> *mc_jet_phi = nullptr;
    std::vector<int> *mc_jet_n_const = nullptr;
    std::vector<int> *mc_jet_n_chr = nullptr;
    std::vector<int> *mc_jet_n_neu = nullptr;
    std::vector<int> *mc_dau_pid = nullptr;
    std::vector<float> *mc_dau_px = nullptr;
    std::vector<float> *mc_dau_py = nullptr;
    std::vector<float> *mc_dau_pz = nullptr;
    std::vector<float> *mc_dau_e = nullptr;
    std::vector<int> *mc_dau_d0_idx = nullptr;

    const bool hasRespMcD0RapidityBranch = inputTree->GetBranch("mc_d0_rapidity");
    const bool hasRespMcJetRapidityBranch = inputTree->GetBranch("mc_jet_rapidity");

    const bool hasRespJetRapidityBranch = inputTree->GetBranch("jet_rapidity");
    const bool hasRespD0RapidityBranch = inputTree->GetBranch("d0_rapidity");

    // Set branch addresses
    inputTree->SetBranchAddress("jet_pt", &jet_pt);
    inputTree->SetBranchAddress("jet_eta", &jet_eta);
    if (hasRespJetRapidityBranch)
        inputTree->SetBranchAddress("jet_rapidity", &jet_rapidity);
    inputTree->SetBranchAddress("jet_phi", &jet_phi);
    if (inputTree->GetBranch("jet_mass"))
        inputTree->SetBranchAddress("jet_mass", &jet_mass);
    inputTree->SetBranchAddress("jet_n_const", &jet_n_const);
    inputTree->SetBranchAddress("jet_n_d0", &jet_n_d0);

    inputTree->SetBranchAddress("d0_pt", &d0_pt);
    inputTree->SetBranchAddress("d0_eta", &d0_eta);
    if (hasRespD0RapidityBranch)
        inputTree->SetBranchAddress("d0_rapidity", &d0_rapidity);
    inputTree->SetBranchAddress("d0_e", &d0_e);
    inputTree->SetBranchAddress("d0_pz", &d0_pz);
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
    if (hasRespMcD0RapidityBranch)
        inputTree->SetBranchAddress("mc_d0_rapidity", &mc_d0_rapidity);
    inputTree->SetBranchAddress("mc_d0_e", &mc_d0_e);
    inputTree->SetBranchAddress("mc_d0_pz", &mc_d0_pz);
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
    // Set branch addresses for direct jet kinematics if available
    if (inputTree->GetBranch("mc_jet_pt")) inputTree->SetBranchAddress("mc_jet_pt", &mc_jet_pt);
    if (inputTree->GetBranch("mc_jet_eta")) inputTree->SetBranchAddress("mc_jet_eta", &mc_jet_eta);
    if (hasRespMcJetRapidityBranch)
        inputTree->SetBranchAddress("mc_jet_rapidity", &mc_jet_rapidity);
    if (inputTree->GetBranch("mc_jet_phi")) inputTree->SetBranchAddress("mc_jet_phi", &mc_jet_phi);
    inputTree->SetBranchAddress("mc_jet_n_const", &mc_jet_n_const);
    inputTree->SetBranchAddress("mc_jet_n_charged", &mc_jet_n_chr);
    inputTree->SetBranchAddress("mc_jet_n_neutral", &mc_jet_n_neu);
    inputTree->SetBranchAddress("mc_dau_pid", &mc_dau_pid);
    inputTree->SetBranchAddress("mc_dau_px", &mc_dau_px);
    inputTree->SetBranchAddress("mc_dau_py", &mc_dau_py);
    inputTree->SetBranchAddress("mc_dau_pz", &mc_dau_pz);
    inputTree->SetBranchAddress("mc_dau_e", &mc_dau_e);
    inputTree->SetBranchAddress("mc_dau_d0_idx", &mc_dau_d0_idx);

    auto getRespJetRapidity = [&](size_t idx) -> float
    {
        if (hasRespJetRapidityBranch && jet_rapidity && idx < jet_rapidity->size())
            return (*jet_rapidity)[idx];
        if (jet_pt && jet_eta && jet_mass && idx < jet_pt->size() && idx < jet_eta->size() && idx < jet_mass->size())
            return calcRapidityFromPtEtaMass((*jet_pt)[idx], (*jet_eta)[idx], (*jet_mass)[idx]);
        if (jet_eta && idx < jet_eta->size())
            return (*jet_eta)[idx];
        return -999.f;
    };

    auto getRespD0Rapidity = [&](size_t idx) -> float
    {
        if (hasRespD0RapidityBranch && d0_rapidity && idx < d0_rapidity->size())
            return (*d0_rapidity)[idx];
        if (d0_e && d0_pz && idx < d0_e->size() && idx < d0_pz->size())
            return calcRapidityFromEnergyPz((*d0_e)[idx], (*d0_pz)[idx]);
        return -999.f;
    };

    auto getRespMcD0Rapidity = [&](size_t idx) -> float
    {
        if (hasRespMcD0RapidityBranch && mc_d0_rapidity && idx < mc_d0_rapidity->size())
            return (*mc_d0_rapidity)[idx];
        if (mc_d0_e && mc_d0_pz && idx < mc_d0_e->size() && idx < mc_d0_pz->size())
            return calcRapidityFromEnergyPz((*mc_d0_e)[idx], (*mc_d0_pz)[idx]);
        return -999.f;
    };

    auto getRespMcJetRapidity = [&](size_t idx) -> float
    {
        if (hasRespMcJetRapidityBranch && mc_jet_rapidity && idx < mc_jet_rapidity->size())
            return (*mc_jet_rapidity)[idx];
        if (mc_jet_eta && idx < mc_jet_eta->size())
            return (*mc_jet_eta)[idx];
        return -999.f;
    };

    // Optionally read event multiplicity branch and load multiplicity weights when requested
    int event_multiplicity = 0;
    bool has_event_multiplicity = false;
    std::map<int, double> multiplicityWeights;
    // Jet-pt weights histogram for response (if available)
    TH1 *h_jet_weights = nullptr;
    TFile *jwfile = nullptr;
    if (useWeights && useJetPtWeights) {
        std::cerr << "ERROR: Both multiplicity and jet-pt weights requested for response creation. They are mutually exclusive. Aborting." << std::endl;
        return;
    }
    if (useWeights) {
        std::cout << "Loading multiplicity weights for pp MC..." << std::endl;
        if (inputTree->GetBranch("event_multiplicity")) {
            inputTree->SetBranchAddress("event_multiplicity", &event_multiplicity);
            has_event_multiplicity = true;
        }
        TFile *wfile = TFile::Open(multweightsfile, "READ");
        if (wfile && !wfile->IsZombie()) {
            TTree *wtree = (TTree*)wfile->Get("multiplicity_weights");
            if (wtree) {
                Int_t mult = 0; Double_t w = 1.0;
                wtree->SetBranchAddress("multiplicity", &mult);
                wtree->SetBranchAddress("weight", &w);
                for (Long64_t iw = 0; iw < wtree->GetEntries(); ++iw) {
                    wtree->GetEntry(iw);
                    multiplicityWeights[(int)mult] = (double)w;
                }
                std::cout << "Loaded " << multiplicityWeights.size() << " multiplicity weights." << std::endl;
            }
            wfile->Close();
        }
    }
    // Load jet-pt weights only when requested (and not using multiplicity weights)
    if (useJetPtWeights) {
        if (multweightsfile) {
            TString jpw(multweightsfile);
            jpw.ReplaceAll("_smoothed.root", "_jetpt_weights.root");
            jwfile = TFile::Open(jpw, "READ");
            if (jwfile && !jwfile->IsZombie()) {
                h_jet_weights = dynamic_cast<TH1*>(jwfile->Get("h_jet_weights"));
                if (h_jet_weights) {
                    h_jet_weights->SetDirectory(0);
                    std::cout << "Loaded jet-pt weights histogram from: " << jpw << std::endl;
                } else {
                    std::cerr << "[WARN] Could not find 'h_jet_weights' in " << jpw << std::endl;
                }
                jwfile->Close();
            }
        }
    }

    // Add event weight branch to response tree (always present, default 1.0)
    float r_event_weight = 1.0f;
    responseTree->Branch("event_weight", &r_event_weight, "event_weight/F");

    // Define a function to calculate dR between two y,phi points
    auto deltaR = [](float y1, float phi1, float y2, float phi2)
    {
        float dy = y1 - y2;
        float dphi = phi1 - phi2;
        while (dphi > M_PI)
            dphi -= 2 * M_PI;
        while (dphi < -M_PI)
            dphi += 2 * M_PI;
        return sqrt(dy * dy + dphi * dphi);
    };
    // // Define a function to calculate dR between two eta,phi points
    // auto deltaR = [](float eta1, float phi1, float eta2, float phi2)
    // {
    //     float deta = eta1 - eta2;
    //     float dphi = phi1 - phi2;
    //     while (dphi > M_PI)
    //         dphi -= 2 * M_PI;
    //     while (dphi < -M_PI)
    //         dphi += 2 * M_PI;
    //     return sqrt(deta * deta + dphi * dphi);
    // };

    // Process events
    Long64_t nEntries = inputTree->GetEntries();
    int matches_found = 0;

    Long64_t updateIntervalResp = nEntries / 200; // ~200 updates
    if (updateIntervalResp < 1) updateIntervalResp = 1;

    std::cout << "Processing " << nEntries << " events for response matrix" << std::endl;

    for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++)
    {
        if (iEntry % updateIntervalResp == 0 || iEntry + 1 == nEntries)
        {
            printProgress(iEntry, nEntries);
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
            r_d0_rapidity_mc = getRespMcD0Rapidity(iMC);
            float d0_mc_rapidity = 0;
            if ((*mc_d0_e)[iMC] > (*mc_d0_pz)[iMC]) {
                d0_mc_rapidity = 0.5 * log(((*mc_d0_e)[iMC] + (*mc_d0_pz)[iMC]) / ((*mc_d0_e)[iMC] - (*mc_d0_pz)[iMC]));
            } else {
                d0_mc_rapidity = -999; // Assign an invalid value if rapidity cannot be calculated
            }
            r_d0_y_mc = d0_mc_rapidity;
            r_d0_phi_mc = (*mc_d0_phi)[iMC];
            r_d0_mass_mc = (*mc_d0_mass)[iMC];
            r_d0_z_mc = (*mc_d0_z)[iMC];
            r_d0_is_primary = (iMC < mc_d0_origin->size()) ? ((*mc_d0_origin)[iMC] == 1) : 0;

            // Get MC jet kinematics: prefer direct branches if present; fallback to 4-vector reconstruction
            bool haveDirect = (mc_jet_pt && mc_jet_eta && mc_jet_phi &&
                               mc_jet_idx >= 0 && mc_jet_idx < static_cast<int>(mc_jet_pt->size()) &&
                               mc_jet_idx < static_cast<int>(mc_jet_eta->size()) &&
                               mc_jet_idx < static_cast<int>(mc_jet_phi->size()));
            if (haveDirect)
            {
                r_jet_pt_mc = (*mc_jet_pt)[mc_jet_idx];
                r_jet_eta_mc = (*mc_jet_eta)[mc_jet_idx];
                r_jet_rapidity_mc = getRespMcJetRapidity(mc_jet_idx);
                r_jet_phi_mc = (*mc_jet_phi)[mc_jet_idx];
            }
            else if (mc_jet_idx >= 0 && mc_jet_px && mc_jet_idx < static_cast<int>(mc_jet_px->size()))
            {
                TLorentzVector mc_jet_tmp;
                mc_jet_tmp.SetPxPyPzE(
                    (*mc_jet_px)[mc_jet_idx],
                    (*mc_jet_py)[mc_jet_idx],
                    (*mc_jet_pz)[mc_jet_idx],
                    (*mc_jet_e)[mc_jet_idx]);
                r_jet_pt_mc = mc_jet_tmp.Pt();
                r_jet_eta_mc = mc_jet_tmp.Eta();
                r_jet_rapidity_mc = mc_jet_tmp.Rapidity();
                r_jet_phi_mc = mc_jet_tmp.Phi();
            }
            else
            {
                // Final fallback if index invalid
                r_jet_pt_mc = (r_d0_z_mc != 0 ? r_d0_pt_mc / r_d0_z_mc : 0);
                r_jet_eta_mc = r_d0_eta_mc;
                r_jet_rapidity_mc = r_d0_rapidity_mc;
                r_jet_phi_mc = r_d0_phi_mc;
            }

            // Constituent counts
            if (mc_jet_n_const && mc_jet_idx >= 0 && mc_jet_idx < static_cast<int>(mc_jet_n_const->size()))
            {
                r_jet_nconst_mc = (*mc_jet_n_const)[mc_jet_idx];
            }
            else if (mc_jet_n_chr && mc_jet_n_neu && mc_jet_idx >= 0 &&
                     mc_jet_idx < static_cast<int>(mc_jet_n_chr->size()) &&
                     mc_jet_idx < static_cast<int>(mc_jet_n_neu->size()))
            {
                r_jet_nconst_mc = (*mc_jet_n_chr)[mc_jet_idx] + (*mc_jet_n_neu)[mc_jet_idx];
            }
            else
            {
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
            r_d0_rapidity_det = getRespD0Rapidity(matched_d0_idx);
            r_d0_phi_det = (*d0_phi)[matched_d0_idx];
            r_d0_mass_det = (*d0_mass)[matched_d0_idx];
            if((*d0_in_jet)[matched_d0_idx] == 1 && (*d0_z)[matched_d0_idx] > 1) // Only apply z scaling if D0 is in jet and z is non-zero to avoid division by zero
                r_d0_z_det = (*d0_z)[matched_d0_idx]/1000.0;
            else
                r_d0_z_det = (*d0_z)[matched_d0_idx];
            float d0_det_rapidity = 0;
            if ((*d0_e)[matched_d0_idx] > (*d0_pz)[matched_d0_idx]) {
                d0_det_rapidity = 0.5 * log(((*d0_e)[matched_d0_idx] + (*d0_pz)[matched_d0_idx]) / ((*d0_e)[matched_d0_idx] - (*d0_pz)[matched_d0_idx]));
            } else {
                d0_det_rapidity = -999; // Assign an invalid value if rapidity cannot be calculated
            }
            r_d0_y_det = d0_det_rapidity;

            // Reconstructed jet info
            if (reco_jet_idx < (int)jet_pt->size())
            {
                r_jet_pt_det = (*jet_pt)[reco_jet_idx];
                r_jet_eta_det = (*jet_eta)[reco_jet_idx];
                r_jet_rapidity_det = getRespJetRapidity(reco_jet_idx);
                r_jet_phi_det = (*jet_phi)[reco_jet_idx];
                r_jet_nconst_det = (*jet_n_const)[reco_jet_idx];
                r_jet_ntags_det = (*jet_n_d0)[reco_jet_idx];
            }
            else
            {
                continue; // Skip if jet index out of bounds
            }

            // --- Apply the same analysis-level cuts used in the tuple maker where possible ---
            // D0-level cuts
            if (r_d0_pt_det < 1.0) continue; // pT cut
            // if (r_d0_eta_det < 2.0 || r_d0_eta_det > 4.5) continue; // eta acceptance
            if (r_d0_rapidity_det < 2.0 || r_d0_rapidity_det > 4.5) continue; // rapidity acceptance
            if (std::fabs(r_d0_mass_det - 1.865) > 0.07) continue; // mass window
            // Vertex chi2 cut not available here (branch not read). If present, consider adding it.

            // Jet-level cuts: when response is created only for jetMode==1, require same jet kinematics
            if (jetMode == 1)
            {
                if (r_jet_pt_det < 5.0) continue; // jet pT cut
                // if (r_jet_eta_det < 2.5 || r_jet_eta_det > 4.0) continue; // jet eta range
                if (r_jet_rapidity_det < 2.5 || r_jet_rapidity_det > 4.0) continue; // jet rapidity range
            }

            // Calculate distances
            r_d0_dr = deltaR(r_d0_rapidity_mc, r_d0_phi_mc, r_d0_rapidity_det, r_d0_phi_det);
            r_jet_dr = deltaR(r_jet_rapidity_mc, r_jet_phi_mc, r_jet_rapidity_det, r_jet_phi_det);

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
                    if ((*mc_dau_d0_idx)[i] != static_cast<int>(iMC))
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
            // Fill event weight only when requested (useWeights==true) and multiplicity info exists
            // Apply either multiplicity weights or jet-pt weights (mutually exclusive)
            r_event_weight = 1.0f;
            if (useWeights && has_event_multiplicity && (event_multiplicity<600)) {
                auto itw = multiplicityWeights.find(event_multiplicity);
                if (itw != multiplicityWeights.end()) r_event_weight = static_cast<float>(itw->second);
                else r_event_weight = 1.0f;
            } else if (useJetPtWeights) {
                if (h_jet_weights && r_jet_pt_det > 0) {
                    int bin = h_jet_weights->FindBin(static_cast<double>(r_jet_pt_det));
                    double wj = h_jet_weights->GetBinContent(bin);
                    if (wj > 0) r_event_weight = r_event_weight * static_cast<float>(wj);
                    // std::cout << "Applied jet-pt weight: " << wj << " for jet pT: " << r_jet_pt_det << std::endl;
                }
            }
            responseTree->Fill();
            matches_found++;
        }
    }

    std::cout << "Found " << matches_found << " matched D0-jet pairs for response matrix" << std::endl;

    // Write and close
    outputRoot->cd();
    responseTree->Write();
    outputRoot->Close();
    if (inputRoot) inputRoot->Close();
    if (inputChain) { delete inputChain; inputChain = nullptr; }

    std::cout << "Response matrix saved to: " << outputFile << std::endl;
}