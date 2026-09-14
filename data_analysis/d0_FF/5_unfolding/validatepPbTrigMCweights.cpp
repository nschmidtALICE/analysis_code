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

void validatepPbTrigMCweights(const std::string &inputFileNameMBMC,
                              const std::string &inputFileNameTrgMC,
                              const std::string &inputWeightsFile,
                            //   float weightFitVal = 1. / 13.5436)
                              float weightFitVal = 1. / 44.4518)
{
    std::cout << "Validating pPb trigger MC weights..." << std::endl;
    std::cout << "Input MB MC file: " << inputFileNameMBMC << std::endl;
    std::cout << "Input Trigger MC file: " << inputFileNameTrgMC << std::endl;
    std::cout << "Input weights file: " << inputWeightsFile << std::endl;
    std::cout << "Using weight fit value: " << weightFitVal << std::endl;

    // small helper to print a simple console progress bar
    auto printProgress = [](Long64_t idx, Long64_t total, const std::string &label)
    {
        if (total <= 0)
            return;
        const int barWidth = 40;
        double fraction = double(idx + 1) / double(total);
        if (fraction < 0)
            fraction = 0;
        if (fraction > 1)
            fraction = 1;
        int pos = int(barWidth * fraction);
        std::ostringstream oss;
        oss << "\r" << label << " [";
        for (int i = 0; i < barWidth; ++i)
        {
            if (i < pos)
                oss << "=";
            else if (i == pos)
                oss << ">";
            else
                oss << " ";
        }
        oss << "] " << std::fixed << std::setprecision(1) << (fraction * 100.0) << "% (" << (idx + 1) << "/" << total << ")";
        std::cout << oss.str() << std::flush;
        if (idx + 1 == total)
            std::cout << std::endl;
    };
    TFile *fMBMC = TFile::Open(inputFileNameMBMC.c_str(), "READ");
    TFile *fTrgMC = TFile::Open(inputFileNameTrgMC.c_str(), "READ");
    if (!fMBMC || fMBMC->IsZombie())
    {
        std::cerr << "Error: Could not open input file " << inputFileNameMBMC << std::endl;
        return;
    }
    if (!fTrgMC || fTrgMC->IsZombie())
    {
        std::cerr << "Error: Could not open input file " << inputFileNameTrgMC << std::endl;
        return;
    }

    TTree *tMBMC = (TTree *)fMBMC->Get("d0jets");
    TTree *tTrgMC = (TTree *)fTrgMC->Get("d0jets");
    if (!tMBMC)
    {
        std::cerr << "Error: Could not find tree 'd0jets' in input file " << inputFileNameMBMC << std::endl;
        return;
    }
    if (!tTrgMC)
    {
        std::cerr << "Error: Could not find tree 'd0jets' in input file " << inputFileNameTrgMC << std::endl;
        return;
    }

    std::cout << "Successfully opened input files and found trees." << std::endl;

    // loop over trees and fill histograms for pt distributions of D0s in MB MC and Trigger MC, then make ratio to get weights
    Long64_t nEntriesMB = tMBMC->GetEntries();
    Long64_t nEntriesTrg = tTrgMC->GetEntries();
    std::cout << "Number of entries in MB MC: " << nEntriesMB << std::endl;
    std::cout << "Number of entries in Trigger MC: " << nEntriesTrg << std::endl;

    TFile *fWeights = TFile::Open(inputWeightsFile.c_str(), "READ");
    if (!fWeights || fWeights->IsZombie())
    {
        std::cerr << "Error: Could not open input weights file " << inputWeightsFile << std::endl;
        return;
    }
    TH1D *hWeights = (TH1D *)fWeights->Get("hWeights");
    if (!hWeights)
    {
        std::cerr << "Error: Could not find histogram 'hWeights' in weights file" << std::endl;
        return;
    }
    int nBins = 50;
    // create validation histograms
    TH1D *hD0PtMCMB = new TH1D("hD0PtMCMB", "D^{0} p_{T} distribution in MB MC; p_{T} (GeV/c); Counts", nBins, 0, 50);
    TH1D *hD0PtMCTrg = new TH1D("hD0PtMCTrg", "D^{0} p_{T} distribution in Trigger MC; p_{T} (GeV/c); Counts", nBins, 0, 50);
    TH1D *hD0PtMCTrigWeighted = new TH1D("hD0PtMCTrigWeighted", "D^{0} p_{T} distribution in MB MC weighted by Trigger MC weights; p_{T} (GeV/c); Counts", nBins, 0, 50);

    // zT distributions (true and reco): z = d0_pt / jet_pt
    TH1D *hZTrueMB = new TH1D("hZTrueMB", "z_{true} distribution in MB MC; z; Counts", 100, 0, 1.2);
    TH1D *hZTrueTrg = new TH1D("hZTrueTrg", "z_{true} distribution in Trigger MC; z; Counts", 100, 0, 1.2);
    TH1D *hZTrueTrgWeighted = new TH1D("hZTrueTrgWeighted", "z_{true} in Trigger MC weighted; z; Counts", 100, 0, 1.2);

    TH1D *hZRecoMB = new TH1D("hZRecoMB", "z_{reco} distribution in MB MC; z; Counts", 100, 0, 1.2);
    TH1D *hZRecoTrg = new TH1D("hZRecoTrg", "z_{reco} distribution in Trigger MC; z; Counts", 100, 0, 1.2);
    TH1D *hZRecoTrgWeighted = new TH1D("hZRecoTrgWeighted", "z_{reco} in Trigger MC weighted; z; Counts", 100, 0, 1.2);

    // response matrices (true MC pT vs reconstructed pT)
    TH2D *hRespMB = new TH2D("hRespMB", "MB response matrix; true p_{T} (GeV/c); reco p_{T} (GeV/c)", nBins, 0, 50, nBins, 0, 50);
    TH2D *hRespTrg = new TH2D("hRespTrg", "Trig response matrix; true p_{T} (GeV/c); reco p_{T} (GeV/c)", nBins, 0, 50, nBins, 0, 50);
    TH2D *hRespTrgWeighted = new TH2D("hRespTrgWeighted", "TrigWeighted response matrix; true p_{T} (GeV/c); reco p_{T} (GeV/c)", nBins, 0, 50, nBins, 0, 50);

    TH1D *hJetPtMCMB = new TH1D("hJetPtMCMB", "Jet p_{T} distribution in MB MC; p_{T} (GeV/c); Counts", nBins, 0, 80);
    TH1D *hJetPtMCTrg = new TH1D("hJetPtMCTrg", "Jet p_{T} distribution in Trigger MC; p_{T} (GeV/c); Counts", nBins, 0, 80);
    TH1D *hJetPtMCTrigWeighted = new TH1D("hJetPtMCTrigWeighted", "Jet p_{T} distribution in MB MC weighted by Trigger MC weights; p_{T} (GeV/c); Counts", nBins, 0, 80);

    TH1D *hD0PtRecMB = new TH1D("hD0PtRecMB", "Reconstructed D^{0} p_{T} distribution in MB MC; p_{T} (GeV/c); Counts", nBins, 0, 50);
    TH1D *hD0PtRecTrg = new TH1D("hD0PtRecTrg", "Reconstructed D^{0} p_{T} distribution in Trigger MC; p_{T} (GeV/c); Counts", nBins, 0, 50);
    TH1D *hD0PtRecTrigWeighted = new TH1D("hD0PtRecTrigWeighted", "Reconstructed D^{0} p_{T} distribution in MB MC weighted by Trigger MC weights; p_{T} (GeV/c); Counts", nBins, 0, 50);

    TH1D *hJetPtRecMB = new TH1D("hJetPtRecMB", "Reconstructed Jet p_{T} distribution in MB MC; p_{T} (GeV/c); Counts", nBins, 0, 80);
    TH1D *hJetPtRecTrg = new TH1D("hJetPtRecTrg", "Reconstructed Jet p_{T} distribution in Trigger MC; p_{T} (GeV/c); Counts", nBins, 0, 80);
    TH1D *hJetPtRecTrigWeighted = new TH1D("hJetPtRecTrigWeighted", "Reconstructed Jet p_{T} distribution in MB MC weighted by Trigger MC weights; p_{T} (GeV/c); Counts", nBins, 0, 80);

    // fill histograms for MB MC
    {
        TTreeReader readerMB(tMBMC);
        TTreeReaderValue<std::vector<float>> mc_d0_pt_MB(readerMB, "mc_d0_pt");
        TTreeReaderValue<std::vector<float>> rec_d0_pt_MB(readerMB, "d0_pt");
        TTreeReaderValue<std::vector<float>> mc_jet_pt_MB(readerMB, "mc_jet_pt");
        TTreeReaderValue<std::vector<float>> rec_jet_pt_MB(readerMB, "jet_pt");
        Long64_t entryMB = 0;
        while (readerMB.Next())
        {
            const std::vector<float> &vecD0PtMB = *mc_d0_pt_MB;
            const std::vector<float> &vecRecD0PtMB = *rec_d0_pt_MB;
            const std::vector<float> &vecJetPtMB = *mc_jet_pt_MB;
            const std::vector<float> &vecRecJetPtMB = *rec_jet_pt_MB;
            for (size_t j = 0; j < vecD0PtMB.size(); ++j)
            {
                hD0PtMCMB->Fill(vecD0PtMB[j]);
                if (j < vecJetPtMB.size() && vecJetPtMB[j] > 0.) {
                    double zt = vecD0PtMB[j] / vecJetPtMB[j];
                    hZTrueMB->Fill(zt);
                }
                if (j < vecRecD0PtMB.size())
                    hRespMB->Fill(vecD0PtMB[j], vecRecD0PtMB[j]);
            }
            for (size_t j = 0; j < vecRecD0PtMB.size(); ++j)
            {
                hD0PtRecMB->Fill(vecRecD0PtMB[j]);
                if (j < vecRecJetPtMB.size() && vecRecJetPtMB[j] > 0.) {
                    double zreco = vecRecD0PtMB[j] / vecRecJetPtMB[j];
                    hZRecoMB->Fill(zreco);
                }
            }
            for (size_t j = 0; j < vecJetPtMB.size(); ++j)
            {
                hJetPtMCMB->Fill(vecJetPtMB[j]);
            }
            for (size_t j = 0; j < vecRecJetPtMB.size(); ++j)
            {
                hJetPtRecMB->Fill(vecRecJetPtMB[j]);
            }
            if ((entryMB % 1000) == 0 || entryMB == nEntriesMB - 1)
                printProgress(entryMB, nEntriesMB, "Filling MB MC histograms");
            ++entryMB;
        }
    }

    //fill histograms for Trig MC
    {
        TTreeReader readerTrg(tTrgMC);
        TTreeReaderValue<std::vector<float>> mc_d0_pt_Trg(readerTrg, "mc_d0_pt");
        TTreeReaderValue<std::vector<float>> rec_d0_pt_Trg(readerTrg, "d0_pt");
        TTreeReaderValue<std::vector<float>> mc_jet_pt_Trg(readerTrg, "mc_jet_pt");
        TTreeReaderValue<std::vector<float>> rec_jet_pt_Trg(readerTrg, "jet_pt");
        Long64_t entryTrg = 0;
        while (readerTrg.Next())
        {
            const std::vector<float> &vecD0PtTrg = *mc_d0_pt_Trg;
            const std::vector<float> &vecRecD0PtTrg = *rec_d0_pt_Trg;
            const std::vector<float> &vecJetPtTrg = *mc_jet_pt_Trg;
            const std::vector<float> &vecRecJetPtTrg = *rec_jet_pt_Trg;
            for (size_t j = 0; j < vecD0PtTrg.size(); ++j)
            {
                hD0PtMCTrg->Fill(vecD0PtTrg[j]);
                hD0PtMCTrigWeighted->Fill(vecD0PtTrg[j], weightFitVal);
                if (j < vecJetPtTrg.size() && vecJetPtTrg[j] > 0.) {
                    double zt = vecD0PtTrg[j] / vecJetPtTrg[j];
                    hZTrueTrg->Fill(zt);
                    hZTrueTrgWeighted->Fill(zt, weightFitVal);
                }
                if (j < vecRecD0PtTrg.size()) {
                    hRespTrg->Fill(vecD0PtTrg[j], vecRecD0PtTrg[j]);
                    hRespTrgWeighted->Fill(vecD0PtTrg[j], vecRecD0PtTrg[j], weightFitVal);
                }
            }
            for (size_t j = 0; j < vecRecD0PtTrg.size(); ++j)
            {
                hD0PtRecTrg->Fill(vecRecD0PtTrg[j]);
                hD0PtRecTrigWeighted->Fill(vecRecD0PtTrg[j], weightFitVal);
            }
            for (size_t j = 0; j < vecJetPtTrg.size(); ++j)
            {
                hJetPtMCTrg->Fill(vecJetPtTrg[j]);
                hJetPtMCTrigWeighted->Fill(vecJetPtTrg[j], weightFitVal);
            }
            for (size_t j = 0; j < vecRecJetPtTrg.size(); ++j)
            {
                hJetPtRecTrg->Fill(vecRecJetPtTrg[j]);
                hJetPtRecTrigWeighted->Fill(vecRecJetPtTrg[j], weightFitVal);
                if (j < vecRecJetPtTrg.size() && vecRecJetPtTrg[j] > 0.) {
                    double zreco = vecRecD0PtTrg[j] / vecRecJetPtTrg[j];
                    hZRecoTrg->Fill(zreco);
                    hZRecoTrgWeighted->Fill(zreco, weightFitVal);
                }
            }
        }
    }


    // -- After filling histograms, produce comparison plots and ratios (absolute counts)
    // create output directory with current date
    time_t now = time(nullptr);
    struct tm *lt = localtime(&now);
    char datebuf[64];
    strftime(datebuf, sizeof(datebuf), "%Y-%m-%d", lt);
    std::string outDir = std::string("output_validation_") + datebuf;
    gSystem->mkdir(outDir.c_str(), true);

    // D0 MC-level comparison canvas (absolute counts)
    TCanvas *cD0 = new TCanvas("cD0", "D0 pT Comparison", 900, 700);
    cD0->SetLogy();
    hD0PtMCMB->SetLineColor(kBlue);
    hD0PtMCTrg->SetLineColor(kRed);
    hD0PtMCTrigWeighted->SetLineColor(kMagenta);
    hD0PtMCMB->SetTitle("D0 MC p_{T} comparison (counts)");
    hD0PtMCMB->Draw("HIST");
    hD0PtMCTrg->Draw("HIST SAME");
    hD0PtMCTrigWeighted->Draw("HIST SAME");
    TLegend *legD0 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legD0->AddEntry(hD0PtMCMB, "MB MC", "l");
    legD0->AddEntry(hD0PtMCTrg, "Trig MC", "l");
    legD0->AddEntry(hD0PtMCTrigWeighted, "TrigWeighted", "l");
    legD0->Draw();
    std::string cD0png = outDir + "/D0_pT_comparison_MC_counts.png";
    cD0->SaveAs(cD0png.c_str());

    // D0 reconstructed-level comparison canvas (absolute counts)
    TCanvas *cD0Rec = new TCanvas("cD0Rec", "D0 pT Reco Comparison", 900, 700);
    cD0Rec->SetLogy();
    hD0PtRecMB->SetLineColor(kBlue);
    hD0PtRecTrg->SetLineColor(kRed);
    hD0PtRecTrigWeighted->SetLineColor(kMagenta);
    hD0PtRecMB->SetTitle("D0 Reco p_{T} comparison (counts)");
    hD0PtRecMB->Draw("HIST");
    hD0PtRecTrg->Draw("HIST SAME");
    hD0PtRecTrigWeighted->Draw("HIST SAME");
    TLegend *legD0Rec = new TLegend(0.6, 0.7, 0.9, 0.9);
    legD0Rec->AddEntry(hD0PtRecMB, "MB Reco", "l");
    legD0Rec->AddEntry(hD0PtRecTrg, "Trig Reco", "l");
    legD0Rec->AddEntry(hD0PtRecTrigWeighted, "TrigWeighted Reco", "l");
    legD0Rec->Draw();
    std::string cD0Recpng = outDir + "/D0_pT_comparison_Reco_counts.png";
    cD0Rec->SaveAs(cD0Recpng.c_str());

    // Jet MC-level comparison canvas (absolute counts)
    TCanvas *cJet = new TCanvas("cJet", "Jet pT Comparison", 900, 700);
    cJet->SetLogy();
    hJetPtMCMB->SetLineColor(kBlue);
    hJetPtMCTrg->SetLineColor(kRed);
    hJetPtMCTrigWeighted->SetLineColor(kMagenta);
    hJetPtMCMB->SetTitle("Jet MC p_{T} comparison (counts)");
    hJetPtMCMB->Draw("HIST");
    hJetPtMCTrg->Draw("HIST SAME");
    hJetPtMCTrigWeighted->Draw("HIST SAME");
    TLegend *legJet = new TLegend(0.6, 0.7, 0.9, 0.9);
    legJet->AddEntry(hJetPtMCMB, "MB MC", "l");
    legJet->AddEntry(hJetPtMCTrg, "Trig MC", "l");
    legJet->AddEntry(hJetPtMCTrigWeighted, "TrigWeighted", "l");
    legJet->Draw();
    std::string cJetpng = outDir + "/Jet_pT_comparison_MC_counts.png";
    cJet->SaveAs(cJetpng.c_str());

    // Jet reconstructed-level comparison canvas (absolute counts)
    TCanvas *cJetRec = new TCanvas("cJetRec", "Jet pT Reco Comparison", 900, 700);
    cJetRec->SetLogy();
    hJetPtRecMB->SetLineColor(kBlue);
    hJetPtRecTrg->SetLineColor(kRed);
    hJetPtRecTrigWeighted->SetLineColor(kMagenta);
    hJetPtRecMB->SetTitle("Jet Reco p_{T} comparison (counts)");
    hJetPtRecMB->Draw("HIST");
    hJetPtRecTrg->Draw("HIST SAME");
    hJetPtRecTrigWeighted->Draw("HIST SAME");
    TLegend *legJetRec = new TLegend(0.6, 0.7, 0.9, 0.9);
    legJetRec->AddEntry(hJetPtRecMB, "MB Reco", "l");
    legJetRec->AddEntry(hJetPtRecTrg, "Trig Reco", "l");
    legJetRec->AddEntry(hJetPtRecTrigWeighted, "TrigWeighted Reco", "l");
    legJetRec->Draw();
    std::string cJetRecpng = outDir + "/Jet_pT_comparison_Reco_counts.png";
    cJetRec->SaveAs(cJetRecpng.c_str());

    // Ratio plots: MB / TrigWeighted (use original counts, not normalized)
    TH1D *hD0RatioMC = (TH1D *)hD0PtMCMB->Clone("hD0Ratio_MC_MB_over_TrgW");
    hD0RatioMC->SetTitle("D0 MC MB / TrigWeighted (counts)");
    hD0RatioMC->Divide(hD0PtMCTrigWeighted);
    TH1D *hJetRatioMC = (TH1D *)hJetPtMCMB->Clone("hJetRatio_MC_MB_over_TrgW");
    hJetRatioMC->SetTitle("Jet MC MB / TrigWeighted (counts)");
    hJetRatioMC->Divide(hJetPtMCTrigWeighted);

    TH1D *hD0RatioRec = (TH1D *)hD0PtRecMB->Clone("hD0Ratio_Rec_MB_over_TrgW");
    hD0RatioRec->SetTitle("D0 Reco MB / TrigWeighted (counts)");
    hD0RatioRec->Divide(hD0PtRecTrigWeighted);
    TH1D *hJetRatioRec = (TH1D *)hJetPtRecMB->Clone("hJetRatio_Rec_MB_over_TrgW");
    hJetRatioRec->SetTitle("Jet Reco MB / TrigWeighted (counts)");
    hJetRatioRec->Divide(hJetPtRecTrigWeighted);

    TCanvas *cD0Ratio = new TCanvas("cD0Ratio", "D0 ratio MC MB/TrigWeighted", 900, 700);
    hD0RatioMC->SetLineColor(kBlack);
    hD0RatioMC->GetYaxis()->SetRangeUser(0, 2);
    hD0RatioMC->Draw("EP");
    std::string cD0Rpng = outDir + "/D0_ratio_MC_MB_over_TrgWeighted.png";
    cD0Ratio->SaveAs(cD0Rpng.c_str());

    TCanvas *cD0RatioRec = new TCanvas("cD0RatioRec", "D0 ratio Reco MB/TrigWeighted", 900, 700);
    hD0RatioRec->SetLineColor(kBlack);
    hD0RatioRec->GetYaxis()->SetRangeUser(0, 2);
    hD0RatioRec->Draw("EP");
    std::string cD0RRecpng = outDir + "/D0_ratio_Rec_MB_over_TrgWeighted.png";
    cD0RatioRec->SaveAs(cD0RRecpng.c_str());

    TCanvas *cJetRatio = new TCanvas("cJetRatio", "Jet ratio MC MB/TrigWeighted", 900, 700);
    hJetRatioMC->SetLineColor(kBlack);
    hJetRatioMC->GetYaxis()->SetRangeUser(0, 2);
    hJetRatioMC->Draw("EP");
    std::string cJetRpng = outDir + "/Jet_ratio_MC_MB_over_TrgWeighted.png";
    cJetRatio->SaveAs(cJetRpng.c_str());

    TCanvas *cJetRatioRec = new TCanvas("cJetRatioRec", "Jet ratio Reco MB/TrigWeighted", 900, 700);
    hJetRatioRec->SetLineColor(kBlack);
    hJetRatioRec->GetYaxis()->SetRangeUser(0, 2);
    hJetRatioRec->Draw("EP");
    std::string cJetRRecpng = outDir + "/Jet_ratio_Rec_MB_over_TrgWeighted.png";
    cJetRatioRec->SaveAs(cJetRRecpng.c_str());

    // Build combined MB + TrigWeighted response matrix and ratio
    TH2D *hRespMBPlusTrigW = (TH2D *)hRespMB->Clone("hRespMBPlusTrigW");
    hRespMBPlusTrigW->SetTitle("MB + TrigWeighted response matrix; true p_{T}; reco p_{T}");
    hRespMBPlusTrigW->Add(hRespTrgWeighted);

    TH2D *hRespRatio = (TH2D *)hRespMBPlusTrigW->Clone("hRespRatio_MBPlusTrigW_over_MB");
    hRespRatio->SetTitle("(MB+TrigWeighted)/MB response ratio; true p_{T}; reco p_{T}");
    hRespRatio->Divide(hRespMB);

    // Draw response matrices and ratio
    TCanvas *cResp = new TCanvas("cResp", "Response matrices", 1400, 600);
    cResp->Divide(3, 1);
    cResp->cd(1);
    hRespMB->SetStats(0);
    hRespMB->Draw("COLZ");
    cResp->cd(2);
    hRespTrgWeighted->SetStats(0);
    hRespTrgWeighted->Draw("COLZ");
    cResp->cd(3);
    hRespMBPlusTrigW->SetStats(0);
    hRespMBPlusTrigW->Draw("COLZ");
    std::string cRespPng = outDir + "/response_matrices_MB_TrgW_MBplusTrgW.png";
    cResp->SaveAs(cRespPng.c_str());

    TCanvas *cRespRatio = new TCanvas("cRespRatio", "Response ratio (MB+TrigW)/MB", 900, 700);
    hRespRatio->SetStats(0);
    hRespRatio->Draw("COLZ");
    std::string cRespRatioPng = outDir + "/response_ratio_MBplusTrgW_over_MB.png";
    cRespRatio->SaveAs(cRespRatioPng.c_str());

    // Combine z distributions: MB + TrigWeighted
    TH1D *hZTrueCombined = (TH1D *)hZTrueMB->Clone("hZTrueCombined");
    hZTrueCombined->SetTitle("z_{true} MB + TrigWeighted; z; Counts");
    hZTrueCombined->Add(hZTrueTrgWeighted);

    TH1D *hZRecoCombined = (TH1D *)hZRecoMB->Clone("hZRecoCombined");
    hZRecoCombined->SetTitle("z_{reco} MB + TrigWeighted; z; Counts");
    hZRecoCombined->Add(hZRecoTrgWeighted);

    // Draw z comparisons and ratios
    TCanvas *cZTrue = new TCanvas("cZTrue", "z_true comparison", 900, 700);
    hZTrueMB->SetLineColor(kBlue);
    hZTrueTrgWeighted->SetLineColor(kMagenta);
    hZTrueCombined->SetLineColor(kBlack);
    hZTrueMB->Draw("HIST");
    hZTrueTrgWeighted->Draw("HIST SAME");
    hZTrueCombined->Draw("HIST SAME");
    TLegend *legZ = new TLegend(0.6, 0.7, 0.9, 0.9);
    legZ->AddEntry(hZTrueMB, "MB true", "l");
    legZ->AddEntry(hZTrueTrgWeighted, "TrigWeighted true", "l");
    legZ->AddEntry(hZTrueCombined, "MB + TrigWeighted", "l");
    legZ->Draw();
    std::string cZTruePng = outDir + "/z_true_comparison_MB_TrgW_combined.png";
    cZTrue->SaveAs(cZTruePng.c_str());

    TCanvas *cZReco = new TCanvas("cZReco", "z_reco comparison", 900, 700);
    hZRecoMB->SetLineColor(kBlue);
    hZRecoTrgWeighted->SetLineColor(kMagenta);
    hZRecoCombined->SetLineColor(kBlack);
    hZRecoMB->Draw("HIST");
    hZRecoTrgWeighted->Draw("HIST SAME");
    hZRecoCombined->Draw("HIST SAME");
    TLegend *legZRec = new TLegend(0.6, 0.7, 0.9, 0.9);
    legZRec->AddEntry(hZRecoMB, "MB reco", "l");
    legZRec->AddEntry(hZRecoTrgWeighted, "TrigWeighted reco", "l");
    legZRec->AddEntry(hZRecoCombined, "MB + TrigWeighted", "l");
    legZRec->Draw();
    std::string cZRecoPng = outDir + "/z_reco_comparison_MB_TrgW_combined.png";
    cZReco->SaveAs(cZRecoPng.c_str());

    // Ratio: combined / MB for z distributions
    TH1D *hZTrueRatio = (TH1D *)hZTrueCombined->Clone("hZTrueRatio_combined_over_MB");
    hZTrueRatio->SetTitle("(MB+TrigWeighted)/MB for z_true");
    hZTrueRatio->Divide(hZTrueMB);
    TCanvas *cZTrueRatio = new TCanvas("cZTrueRatio", "z_true ratio", 900, 700);
    hZTrueRatio->Draw("EP");
    std::string cZTrueRatioPng = outDir + "/z_true_ratio_combined_over_MB.png";
    cZTrueRatio->SaveAs(cZTrueRatioPng.c_str());

    TH1D *hZRecoRatio = (TH1D *)hZRecoCombined->Clone("hZRecoRatio_combined_over_MB");
    hZRecoRatio->SetTitle("(MB+TrigWeighted)/MB for z_reco");
    hZRecoRatio->Divide(hZRecoMB);
    TCanvas *cZRecoRatio = new TCanvas("cZRecoRatio", "z_reco ratio", 900, 700);
    hZRecoRatio->Draw("EP");
    std::string cZRecoRatioPng = outDir + "/z_reco_ratio_combined_over_MB.png";
    cZRecoRatio->SaveAs(cZRecoRatioPng.c_str());

    // Write all histograms and canvases to a ROOT file in the output directory
    std::string outRootPath = outDir + "/validate_plots.root";
    TFile *fout = TFile::Open(outRootPath.c_str(), "RECREATE");
    if (fout && !fout->IsZombie())
    {
        fout->cd();
        hD0PtMCMB->Write();
        hD0PtMCTrg->Write();
        hD0PtMCTrigWeighted->Write();
        hZTrueMB->Write();
        hZTrueTrg->Write();
        hZTrueTrgWeighted->Write();
        hZRecoMB->Write();
        hZRecoTrg->Write();
        hZRecoTrgWeighted->Write();
        hJetPtMCMB->Write();
        hJetPtMCTrg->Write();
        hJetPtMCTrigWeighted->Write();
        hD0PtRecMB->Write();
        hD0PtRecTrg->Write();
        hD0PtRecTrigWeighted->Write();
        hJetPtRecMB->Write();
        hJetPtRecTrg->Write();
        hJetPtRecTrigWeighted->Write();
        hD0RatioMC->Write();
        hJetRatioMC->Write();
        hD0RatioRec->Write();
        hJetRatioRec->Write();
        cD0->Write();
        cD0Rec->Write();
        cJet->Write();
        cJetRec->Write();
        cD0Ratio->Write();
        cD0RatioRec->Write();
        cJetRatio->Write();
        cJetRatioRec->Write();
        cZTrue->Write();
        cZReco->Write();
        cZTrueRatio->Write();
        cZRecoRatio->Write();
        fout->Close();
    }
}