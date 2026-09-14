#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH1D.h>
#include <vector>
#include <TMath.h>
#include <iostream>
#include <string>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>

void makepPbTrigMCweights(const std::string& inputFileNameMBMC, 
                          const std::string& inputFileNameTrgMC, 
                          const std::string& outputFileName)
{
    std::cout << "Making pPb trigger MC weights..." << std::endl;

    // small helper to print a simple console progress bar
    auto printProgress = [](Long64_t idx, Long64_t total, const std::string &label){
        if(total <= 0) return;
        const int barWidth = 40;
        double fraction = double(idx+1) / double(total);
        if(fraction < 0) fraction = 0;
        if(fraction > 1) fraction = 1;
        int pos = int(barWidth * fraction);
        std::ostringstream oss;
        oss << "\r" << label << " [";
        for(int i=0;i<barWidth;++i){
            if(i < pos) oss << "=";
            else if(i==pos) oss << ">";
            else oss << " ";
        }
        oss << "] " << std::fixed << std::setprecision(1) << (fraction*100.0) << "% (" << (idx+1) << "/" << total << ")";
        std::cout << oss.str() << std::flush;
        if(idx+1 == total) std::cout << std::endl;
    };
    // Support input being either a single .root file or a .txt list of ROOT files
    TFile* fMBMC = nullptr;
    TFile* fTrgMC = nullptr;
    TChain* chainMB = nullptr;
    TChain* chainTrg = nullptr;
    TTree* tMBMC = nullptr;
    TTree* tTrgMC = nullptr;

    auto makeChainFromTxt = [&](const std::string &txtpath, TChain* &chainOut)->bool {
        std::ifstream listf(txtpath);
        if (!listf.is_open()) return false;
        chainOut = new TChain("d0jets");
        std::string line;
        size_t added = 0;
        while (std::getline(listf, line)) {
            // trim
            auto l = line.find_first_not_of(" \t\r\n");
            if (l == std::string::npos) continue;
            auto r = line.find_last_not_of(" \t\r\n");
            std::string part = line.substr(l, r - l + 1);
            if (part.empty()) continue;
            if (part[0] == '#') continue;
            Long64_t nAdded = chainOut->Add(part.c_str());
            if (nAdded > 0) ++added;
            else std::cerr << "[WARN] Failed to add '" << part << "' from " << txtpath << std::endl;
        }
        return (added > 0);
    };

    // MB input
    if (inputFileNameMBMC.size() >= 4 && inputFileNameMBMC.substr(inputFileNameMBMC.size()-4) == ".txt") {
        if (!makeChainFromTxt(inputFileNameMBMC, chainMB)) {
            std::cerr << "Error: Could not build TChain from list " << inputFileNameMBMC << std::endl;
            return;
        }
        tMBMC = chainMB;
    } else {
        fMBMC = TFile::Open(inputFileNameMBMC.c_str(), "READ");
        if (!fMBMC || fMBMC->IsZombie()) {
            std::cerr << "Error: Could not open input file " << inputFileNameMBMC << std::endl;
            return;
        }
        tMBMC = (TTree*)fMBMC->Get("d0jets");
        if (!tMBMC) {
            std::cerr << "Error: Could not find tree 'd0jets' in input file " << inputFileNameMBMC << std::endl;
            if (fMBMC) fMBMC->Close();
            return;
        }
    }

    // Trigger input
    if (inputFileNameTrgMC.size() >= 4 && inputFileNameTrgMC.substr(inputFileNameTrgMC.size()-4) == ".txt") {
        if (!makeChainFromTxt(inputFileNameTrgMC, chainTrg)) {
            std::cerr << "Error: Could not build TChain from list " << inputFileNameTrgMC << std::endl;
            if (fMBMC) fMBMC->Close();
            return;
        }
        tTrgMC = chainTrg;
    } else {
        fTrgMC = TFile::Open(inputFileNameTrgMC.c_str(), "READ");
        if (!fTrgMC || fTrgMC->IsZombie()) {
            std::cerr << "Error: Could not open input file " << inputFileNameTrgMC << std::endl;
            if (fMBMC) fMBMC->Close();
            return;
        }
        tTrgMC = (TTree*)fTrgMC->Get("d0jets");
        if (!tTrgMC) {
            std::cerr << "Error: Could not find tree 'd0jets' in input file " << inputFileNameTrgMC << std::endl;
            if (fTrgMC) fTrgMC->Close();
            if (fMBMC) fMBMC->Close();
            return;
        }
    }

    std::cout << "Successfully opened input files and found trees." << std::endl;

    //loop over trees and fill histograms for pt distributions of D0s in MB MC and Trigger MC, then make ratio to get weights
    Long64_t nEntriesMB = tMBMC->GetEntries();
    Long64_t nEntriesTrg = tTrgMC->GetEntries();
    std::cout << "Number of entries in MB MC: " << nEntriesMB << std::endl;
    std::cout << "Number of entries in Trigger MC: " << nEntriesTrg << std::endl;

    // Define histograms for pt distributions
    TH1D* hPtMB = new TH1D("hPtMB", "p_{T} distribution in MB MC; p_{T} (GeV/c); Counts", 50, 0, 50);
    TH1D* hPtTrg = new TH1D("hPtTrg", "p_{T} distribution in Trigger MC; p_{T} (GeV/c); Counts", 50, 0, 50);
    // reconstructed-level histograms
    TH1D* hRecPtMB = new TH1D("hRecPtMB", "Reconstructed p_{T} distribution in MB MC; p_{T} (GeV/c); Counts", 50, 0, 50);
    TH1D* hRecPtTrg = new TH1D("hRecPtTrg", "Reconstructed p_{T} distribution in Trigger MC; p_{T} (GeV/c); Counts", 50, 0, 50);
    // 2D histograms for d0_pt vs jet_pt (MC and reco)
    TH2D* h2MC_MB = new TH2D("h2MC_MB", "MC d0_pt vs mc_jet_pt MB; d0 p_{T} (GeV/c); mc jet p_{T} (GeV/c)", 50, 0, 50, 40, 0, 80);
    TH2D* h2MC_Trg = new TH2D("h2MC_Trg", "MC d0_pt vs mc_jet_pt Trg; d0 p_{T} (GeV/c); mc jet p_{T} (GeV/c)", 50, 0, 50, 40, 0, 80);
    TH2D* h2MC_Weights = nullptr;

    TH2D* h2Rec_MB = new TH2D("h2Rec_MB", "Reco d0_pt vs jet_pt MB; d0 p_{T} (GeV/c); reco jet p_{T} (GeV/c)", 50, 0, 50, 40, 0, 80);
    TH2D* h2Rec_Trg = new TH2D("h2Rec_Trg", "Reco d0_pt vs jet_pt Trg; d0 p_{T} (GeV/c); reco jet p_{T} (GeV/c)", 50, 0, 50, 40, 0, 80);
    TH2D* h2Rec_Weights = nullptr;

    // Loop over MB MC tree
        // Loop over MB MC tree using TTreeReader
        {
            TTreeReader readerMB(tMBMC);
            TTreeReaderValue<std::vector<float>> mc_d0_pt_MB(readerMB, "mc_d0_pt");
            TTreeReaderValue<std::vector<float>> rec_d0_pt_MB(readerMB, "d0_pt");
            TTreeReaderValue<std::vector<float>> mc_jet_pt_MB(readerMB, "mc_jet_pt");
            TTreeReaderValue<std::vector<float>> rec_jet_pt_MB(readerMB, "jet_pt");
            Long64_t entryMB = 0;
            while (readerMB.Next()) {
                const std::vector<float> &vec = *mc_d0_pt_MB;
                const std::vector<float> &vecRec = *rec_d0_pt_MB;
                const std::vector<float> &vecJetMC = *mc_jet_pt_MB;
                const std::vector<float> &vecJetRec = *rec_jet_pt_MB;
                // Fill 1D and 2D MB histograms. Fill matching by index when possible; otherwise fill all combinations conservatively.
                for (size_t j = 0; j < vec.size(); ++j) {
                    hPtMB->Fill(vec[j]);
                    if (j < vecJetMC.size()) h2MC_MB->Fill(vec[j], vecJetMC[j]);
                }
                for (size_t j = 0; j < vecRec.size(); ++j) {
                    hRecPtMB->Fill(vecRec[j]);
                    if (j < vecJetRec.size()) h2Rec_MB->Fill(vecRec[j], vecJetRec[j]);
                }
                if ((entryMB % 1000) == 0 || entryMB == nEntriesMB-1) printProgress(entryMB, nEntriesMB, "MB MC");
                ++entryMB;
            }
        }
    // Loop over Trigger MC tree
        // Loop over Trigger MC tree using TTreeReader
        {
            TTreeReader readerTrg(tTrgMC);
            TTreeReaderValue<std::vector<float>> mc_d0_pt_Trg(readerTrg, "mc_d0_pt");
            TTreeReaderValue<std::vector<float>> rec_d0_pt_Trg(readerTrg, "d0_pt");
            TTreeReaderValue<std::vector<float>> mc_jet_pt_Trg(readerTrg, "mc_jet_pt");
            TTreeReaderValue<std::vector<float>> rec_jet_pt_Trg(readerTrg, "jet_pt");
            Long64_t entryTrg = 0;
            while (readerTrg.Next()) {
                const std::vector<float> &vec = *mc_d0_pt_Trg;
                const std::vector<float> &vecRec = *rec_d0_pt_Trg;
                const std::vector<float> &vecJetMC = *mc_jet_pt_Trg;
                const std::vector<float> &vecJetRec = *rec_jet_pt_Trg;
                for (size_t j = 0; j < vec.size(); ++j) {
                    hPtTrg->Fill(vec[j]);
                    if (j < vecJetMC.size()) h2MC_Trg->Fill(vec[j], vecJetMC[j]);
                }
                for (size_t j = 0; j < vecRec.size(); ++j) {
                    hRecPtTrg->Fill(vecRec[j]);
                    if (j < vecJetRec.size()) h2Rec_Trg->Fill(vecRec[j], vecJetRec[j]);
                }
                if ((entryTrg % 1000) == 0 || entryTrg == nEntriesTrg-1) printProgress(entryTrg, nEntriesTrg, "Trg MC");
                ++entryTrg;
            }
        }

    //plot the histograms to check
    TCanvas* c1 = new TCanvas("c1", "pT Distributions", 800, 600);
    c1->SetLogy();
    hPtMB->SetLineColor(kBlue);
    hPtMB->Draw();
    hPtTrg->SetLineColor(kRed);
    hPtTrg->Draw("SAME");
    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(hPtMB, "MB MC", "l");
    legend->AddEntry(hPtTrg, "Trigger MC", "l");
    legend->Draw();

    // plot reconstructed-level distributions
    TCanvas* c1Rec = new TCanvas("c1Rec", "Reconstructed pT Distributions", 800, 600);
    c1Rec->SetLogy();
    hRecPtMB->SetLineColor(kBlue);
    hRecPtMB->Draw();
    hRecPtTrg->SetLineColor(kRed);
    hRecPtTrg->Draw("SAME");
    TLegend* legendRec = new TLegend(0.6, 0.7, 0.9, 0.9);
    legendRec->AddEntry(hRecPtMB, "MB Reco", "l");
    legendRec->AddEntry(hRecPtTrg, "Trigger Reco", "l");
    legendRec->Draw();

    // create output directory with current date
    time_t now = time(nullptr);
    struct tm *lt = localtime(&now);
    char datebuf[64];
    strftime(datebuf, sizeof(datebuf), "%Y-%m-%d", lt);
    std::string outDir = std::string("output_") + datebuf;
    // create directory (recursive if needed)
    gSystem->mkdir(outDir.c_str(), true);

    // save canvas and histograms into dated directory
    std::string pngPath = outDir + "/pT_distributions.png";
    c1->SaveAs(pngPath.c_str());



    // Create histogram for weights (ratio of Trigger MC to MB MC)
    TH1D* hWeights = (TH1D*)hPtTrg->Clone("hWeights");
    TH1D* hWeightsDenom = (TH1D*)hPtMB->Clone("hWeightsDenom");
    hWeights->Scale(1.0/hWeights->Integral());
    hWeightsDenom->Scale(1.0/hWeightsDenom->Integral());

    hWeights->SetTitle("p_{T} weights (Trigger MC / MB MC); p_{T} (GeV/c); Weight");
    hWeights->Divide(hWeightsDenom);

    // Create histogram for reconstructed-level weights
    TH1D* hWeightsRec = (TH1D*)hRecPtTrg->Clone("hWeightsRec");
    hWeightsRec->SetTitle("Reconstructed p_{T} weights (Trigger MC / MB MC); p_{T} (GeV/c); Weight");
    hWeightsRec->Divide(hRecPtMB);

    // Create 2D weight maps by dividing Trg / MB
    if (h2MC_MB->Integral() > 0) {
        h2MC_Weights = (TH2D*)h2MC_Trg->Clone("h2MC_Weights");
        h2MC_Weights->SetTitle("MC weight map (Trg / MB); d0 p_{T}; mc jet p_{T}");
        h2MC_Weights->Divide(h2MC_MB);
    }
    if (h2Rec_MB->Integral() > 0) {
        h2Rec_Weights = (TH2D*)h2Rec_Trg->Clone("h2Rec_Weights");
        h2Rec_Weights->SetTitle("Reco weight map (Trg / MB); d0 p_{T}; reco jet p_{T}");
        h2Rec_Weights->Divide(h2Rec_MB);
    }

    //plot weights
    TCanvas* c2 = new TCanvas("c2", "pT Weights", 800, 600);
    hWeights->SetLineColor(kGreen+2);
    hWeights->Draw();
    //fit a const function to the weights to smooth them out (optional)
    TF1* fitFunc = new TF1("fitFunc", "pol0", 8, 50);
    hWeights->Fit(fitFunc, "R");
    fitFunc->SetLineColor(kMagenta);
    fitFunc->Draw("SAME");
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.04);
    double fitVal = fitFunc->GetParameter(0);
    latex->DrawLatex(0.15, 0.72, TString::Format("Fit weight: %.3f", fitVal));
    std::string weightPlotPath = outDir + "/pT_weights.png";
    c2->SaveAs(weightPlotPath.c_str());

    // Fit and draw reconstructed-level weights
    TF1* fitFuncRec = new TF1("fitFuncRec", "pol0", 8, 50);
    hWeightsRec->Fit(fitFuncRec, "R");
    fitFuncRec->SetLineColor(kMagenta);
    // draw reco weights on a separate canvas
    TCanvas* c2Rec = new TCanvas("c2Rec", "Reconstructed pT Weights", 800, 600);
    hWeightsRec->SetLineColor(kGreen+2);
    hWeightsRec->Draw();
    fitFuncRec->Draw("SAME");
    double fitValRec = fitFuncRec->GetParameter(0);
    TLatex *latexRec = new TLatex();
    latexRec->SetNDC();
    latexRec->SetTextSize(0.04);
    latexRec->DrawLatex(0.15, 0.72, TString::Format("Fit weight (reco): %.3f", fitValRec));
    std::string weightRecPlotPath = outDir + "/pT_weights_reco.png";
    c2Rec->SaveAs(weightRecPlotPath.c_str());

    // Draw and save 2D weight maps
    if (h2MC_Weights) {
        TCanvas* c2dMC = new TCanvas("c2dMC", "MC 2D weight map", 900, 700);
        h2MC_Weights->SetStats(0);
        h2MC_Weights->Draw("COLZ");
        std::string mapMCPath = outDir + "/MC_2D_weight_map.png";
        c2dMC->SaveAs(mapMCPath.c_str());
    }
    if (h2Rec_Weights) {
        TCanvas* c2dRec = new TCanvas("c2dRec", "Reco 2D weight map", 900, 700);
        h2Rec_Weights->SetStats(0);
        h2Rec_Weights->Draw("COLZ");
        std::string mapRecPath = outDir + "/Reco_2D_weight_map.png";
        c2dRec->SaveAs(mapRecPath.c_str());
    }

    // Save weights and histograms to output file
    std::string outRootPath = outDir + "/" + outputFileName;
    TFile* fout = TFile::Open(outRootPath.c_str(), "RECREATE");
    if (fout && !fout->IsZombie()) {
        fout->cd();
        hPtMB->Write();
        hPtTrg->Write();
        hRecPtMB->Write();
        hRecPtTrg->Write();
        c1->Write();
        c1Rec->Write();
        hWeights->Write();
        hWeightsRec->Write();
        h2MC_MB->Write();
        h2MC_Trg->Write();
        if (h2MC_Weights) h2MC_Weights->Write();
        h2Rec_MB->Write();
        h2Rec_Trg->Write();
        if (h2Rec_Weights) h2Rec_Weights->Write();
        fout->Close();
    } else {
        std::cerr << "Warning: could not create output ROOT file: " << outRootPath << std::endl;
    }
}