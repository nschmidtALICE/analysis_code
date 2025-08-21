#include <iostream>
#include <vector>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2D.h"
#include "TString.h"
#include "TCanvas.h"

void D0Acceptance(TString inputFile = "/media/niviths/SSD2/lhcb_analysis_SSD/20250728_pPb_MC_output/20250728_pPb_MC_output.root", TString outputFile = "D0AcceptanceMap_pPb.root") {

    // Kinematic selection
    double minPt = 0.0, maxPt = 60.0;
    double minEta = 2.0, maxEta = 4.5;
    std::cout << "Calculating D0 acceptance map..." << std::endl;
    // Daughter acceptance
    double dauMinEta = 2.0, dauMaxEta = 5.0;

    // Open file and tree
    TFile* f = TFile::Open(inputFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open input file " << inputFile << std::endl;
        return;
    }
    TTree* t = (TTree*)f->Get("d0jets");
    if (!t) {
        std::cerr << "Error: Could not find tree 'd0jets' in input file" << std::endl;
        return;
    }

    // MC D0 branches as in nTupleMaker.C
    std::vector<int>* mc_d0_pid = nullptr;
    std::vector<float>* mc_d0_pt = nullptr;
    std::vector<float>* mc_d0_eta = nullptr;
    std::vector<float>* mc_d0_phi = nullptr;
    std::vector<float>* mc_d0_mass = nullptr;
    std::vector<int>* mc_d0_origin = nullptr;
    std::vector<int>* mc_d0_matched = nullptr;
    std::vector<int>* mc_d0_jet_idx = nullptr;
    std::vector<float>* mc_d0_z = nullptr;

    std::vector<int>* mc_dau_pid = nullptr;
    std::vector<float>* mc_dau_px = nullptr;
    std::vector<float>* mc_dau_py = nullptr;
    std::vector<float>* mc_dau_pz = nullptr;
    std::vector<float>* mc_dau_e = nullptr;
    std::vector<int>* mc_dau_d0_idx = nullptr;
    std::vector<float>* mc_dau_eta = nullptr;

    t->SetBranchAddress("mc_d0_pid", &mc_d0_pid);
    t->SetBranchAddress("mc_d0_pt", &mc_d0_pt);
    t->SetBranchAddress("mc_d0_eta", &mc_d0_eta);
    t->SetBranchAddress("mc_d0_phi", &mc_d0_phi);
    t->SetBranchAddress("mc_d0_mass", &mc_d0_mass);
    t->SetBranchAddress("mc_d0_origin", &mc_d0_origin);
    t->SetBranchAddress("mc_d0_matched", &mc_d0_matched);
    t->SetBranchAddress("mc_d0_jet_idx", &mc_d0_jet_idx);
    t->SetBranchAddress("mc_d0_z", &mc_d0_z);

    t->SetBranchAddress("mc_dau_pid", &mc_dau_pid);
    t->SetBranchAddress("mc_dau_px", &mc_dau_px);
    t->SetBranchAddress("mc_dau_py", &mc_dau_py);
    t->SetBranchAddress("mc_dau_pz", &mc_dau_pz);
    t->SetBranchAddress("mc_dau_e", &mc_dau_e);
    t->SetBranchAddress("mc_dau_d0_idx", &mc_dau_d0_idx);
    t->SetBranchAddress("mc_dau_eta", &mc_dau_eta);

    std::cout << "Branches set up successfully." << std::endl;

    // Define binning for pt and eta
    const int nPtBins = 300;
    const int nEtaBins = 25;
    double ptBins[nPtBins+1];
    double etaBins[nEtaBins+1];
    double ptMin = minPt, ptMax = maxPt;
    double etaMin = minEta, etaMax = maxEta;
    for (int i = 0; i <= nPtBins; ++i) ptBins[i] = ptMin + (ptMax-ptMin)*i/nPtBins;
    for (int i = 0; i <= nEtaBins; ++i) etaBins[i] = etaMin + (etaMax-etaMin)*i/nEtaBins;

    // Numerator and denominator histograms
    TH2D* hDen = new TH2D("hDen", "D0 in range;pt [GeV];eta", nPtBins, ptBins, nEtaBins, etaBins);
    TH2D* hNum = new TH2D("hNum", "D0 with both daughters in acceptance;pt [GeV];eta", nPtBins, ptBins, nEtaBins, etaBins);

    Long64_t nEntries = t->GetEntries();
    int nD0InRange = 0;
    int nD0WithDausInAcc = 0;

    std::cout << "Processing " << nEntries << " entries..." << std::endl;
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        t->GetEntry(entry);
        //progress
        if (entry % 50000 == 0) {
            std::cout << "Processing entry " << entry << "/" << nEntries 
                      << " (" << (100.0 * entry / nEntries) << "%)" << std::endl;
        }
        size_t nD0 = mc_d0_pid->size();
        for (size_t i = 0; i < nD0; ++i) {
            // Only D0 or anti-D0
            if (abs(mc_d0_pid->at(i)) != 421) continue;
            float pt = mc_d0_pt->at(i);
            float eta = mc_d0_eta->at(i);
            if (pt < minPt || pt > maxPt) continue;
            if (eta < minEta || eta > maxEta) continue;
            nD0InRange++;
            hDen->Fill(pt, eta);

            // Find daughters of this D0 (as in nTupleMaker.C, using mc_dau_d0_idx)
            bool kaonInAcc = false, pionInAcc = false;
            for (size_t j = 0; j < mc_dau_pid->size(); ++j) {
                if (mc_dau_d0_idx->at(j) != (int)i) continue;
                int dau_pid = mc_dau_pid->at(j);
                float dau_eta = mc_dau_eta->at(j);
                if (abs(dau_pid) == 321 && dau_eta >= dauMinEta && dau_eta <= dauMaxEta) kaonInAcc = true;
                if (abs(dau_pid) == 211 && dau_eta >= dauMinEta && dau_eta <= dauMaxEta) pionInAcc = true;
            }
            if (kaonInAcc && pionInAcc) {
                nD0WithDausInAcc++;
                hNum->Fill(pt, eta);
            }
        }
    }

    double acceptance = (nD0InRange > 0) ? double(nD0WithDausInAcc) / nD0InRange : 0.0;
    std::cout << "D0 in range: " << nD0InRange << std::endl;
    std::cout << "D0 with both daughters in acceptance: " << nD0WithDausInAcc << std::endl;
    std::cout << "LHCb acceptance correction: " << acceptance << std::endl;

    // Create acceptance map
    TH2D* hAcc = (TH2D*)hNum->Clone("hAcceptance");
    hAcc->SetTitle("D0 acceptance map;pt [GeV];eta");
    hAcc->Divide(hDen);



    // Write to output file
    TFile* fout = TFile::Open(outputFile, "RECREATE");
    hDen->Write();
    hNum->Write();
    hAcc->Write();
    fout->Close();

    // Plot the acceptance map
    TCanvas* cAcc = new TCanvas("cAcc", "D0 Acceptance Map", 800, 600);
    hAcc->SetStats(0);
    hAcc->GetZaxis()->SetTitle("Acceptance");
    hAcc->Draw("COLZ");
    cAcc->SetRightMargin(0.15);
    TString pngName = "D0AcceptanceMap.png";
    TString pdfName = "D0AcceptanceMap.pdf";
    cAcc->SaveAs(pngName);
    cAcc->SaveAs(pdfName);
    std::cout << "Saved acceptance map plot: " << pngName << std::endl;
    delete cAcc;

    f->Close();
    return;
}