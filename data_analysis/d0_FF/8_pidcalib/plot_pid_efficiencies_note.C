// plot_pidcalib.C
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>
#include <TLatex.h>
#include <ctime>
#include <memory>
#include <cstdio>
#include <string>
#include <vector>

// Updated: added optional species and cut parameters to display on the plot
void plot_pid_efficiencies_note(const char* filenameBase = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-09-03_pPb_DGaussDefault",
                   
    const std::vector<std::string> &jetPtBins = {"10_15", "15_20", "20_30", "30_100"},
    const std::vector<int> &yBins = {0, 1, 2},
    std::vector<double> yBinBorders = {2.5, 3.0, 3.5, 4.0}
                  )
{
  // visual settings
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  //we need to loop over the jetPtBins and then the yBins in order to load all input histograms correctly
  // the files are named as follows: /media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-09-03_pPb_DGaussDefault/{jetPtBins}/tagZAllEfficiencyCorrections_bin{yBins}.root

  //define histogram storage for later plotting. we need to save histograms for each jet pt bin and each y bin. histograms will be for kaon, pion and combined PID efficiency
  std::vector<std::vector<std::unique_ptr<TH1D>>> kaonCorrHists(jetPtBins.size());
  std::vector<std::vector<std::unique_ptr<TH1D>>> pionCorrHists(jetPtBins.size());
  std::vector<std::vector<std::unique_ptr<TH1D>>> recoEffCorrHists(jetPtBins.size());
  for (size_t jetIdx = 0; jetIdx < jetPtBins.size(); ++jetIdx) {
    kaonCorrHists[jetIdx].resize(yBins.size());
    pionCorrHists[jetIdx].resize(yBins.size());
    recoEffCorrHists[jetIdx].resize(yBins.size());
  }

  for (size_t jetIdx = 0; jetIdx < jetPtBins.size(); ++jetIdx) {
    for (size_t yIdx = 0; yIdx < yBins.size(); ++yIdx) {
      const TString inputFile = Form("%s/%s/tagZAllEfficiencyCorrections_bin%d.root",
                                     filenameBase,
                                     jetPtBins[jetIdx].c_str(),
                                     yBins[yIdx]);

      std::unique_ptr<TFile> input(TFile::Open(inputFile, "READ"));
      if (!input || input->IsZombie()) {
        std::printf("ERROR: could not open input file '%s'\n", inputFile.Data());
        continue;
      }

      auto loadHist = [&](const char* histName, std::unique_ptr<TH1D>& slot) {
        TH1D* hist = dynamic_cast<TH1D*>(input->Get(histName));
        if (!hist) {
          std::printf("ERROR: histogram '%s' not found in '%s'\n", histName, inputFile.Data());
          return;
        }
        slot.reset(dynamic_cast<TH1D*>(hist->Clone(Form("%s_%s_bin%d_%s",
                                                     histName,
                                                     jetPtBins[jetIdx].c_str(),
                                                     yBins[yIdx],
                                                     histName))));
        if (!slot) {
          std::printf("ERROR: could not clone histogram '%s' from '%s'\n", histName, inputFile.Data());
        } else {
          slot->SetDirectory(nullptr);
        }
      };

      loadHist("tagZKaonCorrHist", kaonCorrHists[jetIdx][yIdx]);
      loadHist("tagZPionCorrHist", pionCorrHists[jetIdx][yIdx]);
      loadHist("tagZRecoEffCorrHist", recoEffCorrHists[jetIdx][yIdx]);
    }
  }

  //create overview plot of PID efficiencies, the number of jetptbin as number of pads in x and the number of yBins as number of pads in y
  // plot kaon, pion and combined PID efficiencies on each pad. put legend on first pad and label each pad with its corresponding jet pt range and y range

  const int nJetBins = static_cast<int>(jetPtBins.size());
  const int nYBins = static_cast<int>(yBins.size());
  TCanvas* cOverview = new TCanvas("c_pid_eff_overview", "PID efficiency overview", 300 * nYBins, 250 * nJetBins);
  cOverview->Divide(nYBins, nJetBins, 0.001, 0.001);

  TLegend* legend = nullptr;

  for (int jetIdx = 0; jetIdx < nJetBins; ++jetIdx) {
    for (int yIdx = 0; yIdx < nYBins; ++yIdx) {
      const int padIdx = jetIdx * nYBins + yIdx + 1;
      cOverview->cd(padIdx);

      gPad->SetTopMargin(0.02);
      gPad->SetRightMargin(0.03);
      gPad->SetLeftMargin(0.10);
      gPad->SetBottomMargin(0.12);
      gPad->SetTickx();
      gPad->SetTicky();

      TH1D* kaonHist = kaonCorrHists[jetIdx][yIdx].get();
      TH1D* pionHist = pionCorrHists[jetIdx][yIdx].get();
      TH1D* recoHist = recoEffCorrHists[jetIdx][yIdx].get();

      if (!kaonHist || !pionHist || !recoHist) {
        TLatex missing;
        missing.SetNDC();
        missing.SetTextSize(0.04);
        missing.DrawLatex(0.18, 0.50, "missing input histogram");
        continue;
      }

      kaonHist->SetLineColor(kAzure + 2);
      kaonHist->SetLineWidth(2);
      kaonHist->SetMarkerStyle(20);
      kaonHist->SetMarkerColor(kAzure + 2);
      kaonHist->SetMarkerSize(1.2);
      kaonHist->SetTitle("");
      pionHist->SetLineColor(kRed + 1);
      pionHist->SetLineWidth(2);
      pionHist->SetMarkerStyle(21);
      pionHist->SetMarkerColor(kRed + 1);
      pionHist->SetMarkerSize(1.2);
      pionHist->SetTitle("");
      recoHist->SetLineColor(kGreen + 2);
      recoHist->SetLineWidth(2);
      recoHist->SetMarkerStyle(22);
      recoHist->SetMarkerColor(kGreen + 2);
      recoHist->SetMarkerSize(1.4);
      recoHist->SetTitle("");

      kaonHist->SetMinimum(0.0);
      kaonHist->SetMaximum(1.05);
      kaonHist->GetYaxis()->SetTitle("Efficiency");
      kaonHist->GetXaxis()->SetTitle("#it{z}_{T}");

      kaonHist->Draw("pe");
      pionHist->Draw("pe same");
      recoHist->Draw("pe same");

      if (!legend && padIdx == 1) {
        legend = new TLegend(0.18, 0.18, 0.58, 0.30);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.04);
        legend->AddEntry(kaonHist, "Kaon correction", "l");
        legend->AddEntry(pionHist, "Pion correction", "l");
        legend->AddEntry(recoHist, "Reco efficiency correction", "l");
        legend->Draw();
      }

      TLatex labels;
      labels.SetNDC();
      labels.SetTextFont(42);
      labels.SetTextSize(0.035);
      const double yLow = (yIdx < static_cast<int>(yBinBorders.size())) ? yBinBorders[yIdx] : 0.0;
      const double yHigh = (yIdx + 1 < static_cast<int>(yBinBorders.size())) ? yBinBorders[yIdx + 1] : 0.0;
      labels.DrawLatex(0.16, 0.92, Form("y: %.1f - %.1f", yLow, yHigh));
      labels.DrawLatex(0.16, 0.84, Form("jet p_{T}: %s GeV", jetPtBins[jetIdx].c_str()));
    }
  }

  cOverview->cd();
  cOverview->Update();

  const TString outBase = Form("%s/pid_efficiencies_overview", filenameBase);
  cOverview->SaveAs(outBase + ".png");
  cOverview->SaveAs(outBase + ".pdf");


}