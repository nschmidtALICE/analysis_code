// plot_pidcalib.C
#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <TLatex.h>
#include <ctime>
#include <memory>
#include <cstdio>

// Updated: added optional species and cut parameters to display on the plot
void plot_pidcalib(const char* filename = "/media/niviths/local/analysis_code/data_analysis/d0_FF/8_pidcalib/pidcalib_output_20250319/pidcalib/pidcalib_output_Ap_09_pEta_pi/effhists-ApTurbo16-down-Pi-MC15TuneV1_ProbNNpi>0.9&MC15TuneV1_ProbNNghost<0.3&TRCHI2NDOF<3-P.ETA.root",
                   const char* histname = "eff_MC15TuneV1_ProbNNpi>0.9&MC15TuneV1_ProbNNghost<0.3&TRCHI2NDOF<3",
                   const char* species = "",             // "pi" or "K" (auto-detected from histname when empty)
                   double probNNcut = -1.0,               // if <0, will attempt to parse or leave blank
                   double probNNghostCut = -1.0,
                   double trchi2ndofCut = -1.0)
{
  // visual settings
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  // decide species: if user passed species use that, else attempt to infer from histname
  bool isPion = false;
  if (TString(species).Length() > 0) {
    TString sp(species);
    if (sp.Contains("pi", TString::kIgnoreCase)) isPion = true;
  } else {
    if (TString(histname).Contains("ProbNNpi")) isPion = true;
  }
  TString part = isPion ? "Pion" : "Kaon";
  TString partPlot = isPion ? "#pi" : "K";

  // open file
  std::unique_ptr<TFile> f(TFile::Open(filename));
  if (!f || f->IsZombie()) {
    std::printf("ERROR: could not open file '%s'\n", filename);
    return;
  }

  // retrieve histogram (works for TH2D/TH2F via TH2 base)
  TH2D* h = dynamic_cast<TH2D*>(f->Get(histname));
  if (!h) {
    std::printf("ERROR: histogram '%s' not found in '%s'. Keys in file:\n", histname, filename);
    f->ls();
    return;
  }

  // place outputs in the same directory as the input file
  std::string fname(filename);
  std::string outdir;
  size_t pos = fname.find_last_of("/\\");
  if (pos == std::string::npos) {
    outdir = "."; // current directory
  } else {
    outdir = fname.substr(0, pos);
  }
  // ensure directory exists (no-op if it does)
  gSystem->mkdir(outdir.c_str(), kTRUE);

  // optional: set axis titles appropriate for PID efficiency (edit as needed)
  // h->SetTitle("PID efficiency; momentum p [GeV/c]; p_{T} [GeV/c]");
  // draw
  TCanvas* c = new TCanvas("c_pid", "PID efficiency (2D)", 900, 700);
  c->cd();
  c->SetRightMargin(0.12); // make room for color palette
  c->SetLeftMargin(0.08);
  c->SetBottomMargin(0.08);
  c->SetTopMargin(0.01); // increased to make room for TLatex labels
  //set tickx and ticky
    c->SetTickx();
    c->SetTicky();
  h->SetTitle("");
  h->GetXaxis()->SetTitle(Form("#it{p}^{%s} [MeV/#it{c}]", partPlot.Data()));
  h->GetYaxis()->SetTitleOffset(0.8);
  h->GetYaxis()->SetTitle(Form("#eta^{%s}", partPlot.Data()));
  h->GetZaxis()->SetTitle(Form("PID efficiency for %s", part.Data()));
  h->GetZaxis()->SetRangeUser(0.0, 1.0); // efficiency range
  h->GetXaxis()->SetRangeUser(2000,100000);
  h->Draw("colz");
  // // Draw a cropped view so the canvas only shows the requested momentum window.
  // std::unique_ptr<TH2D> hview(dynamic_cast<TH2D*>(h->Clone("h_pid_view")));
  // if (!hview) {
  //   std::printf("ERROR: could not clone histogram '%s'\n", histname);
  //   return;
  // }
  // const int firstXBin = hview->GetXaxis()->FindBin(2.0 + 1e-9);
  // const int lastXBin = hview->GetXaxis()->FindBin(100.0 - 1e-9);
  // hview->GetXaxis()->SetRange(firstXBin, lastXBin);
  // hview->Draw("colz");

  // prepare TLatex and draw labels
  TLatex tex;
  tex.SetNDC();
  tex.SetTextFont(62); // bold-ish for first line
  tex.SetTextSize(0.04);

  // LHCb in-progress label (top-left)
  tex.DrawLatexNDC(0.12, 0.94, "#font[12]{LHCb} (in-progress)");

  // particle species label (below LHCb label)
  tex.SetTextFont(42);
  tex.SetTextSize(0.04);
  TString speciesLabel = Form("%s PID efficiency", part.Data());
  tex.DrawLatexNDC(0.12, 0.90, speciesLabel);

  // cuts label (below species)
  tex.SetTextSize(0.033);

  // If user didn't provide numeric cuts, try to extract from histname simple patterns
  if (probNNcut < 0) {
    // attempt simple parse for "ProbNNpi>0.9" or "ProbNNK>0.9"
    TString hn(histname);
    Ssiz_t p1 = hn.Index("ProbNN");
    if (p1 != kNPOS) {
      // find '>' after ProbNN
      Ssiz_t gt = hn.Index('>', p1);
      if (gt != kNPOS) {
        TString num = hn(gt+1, 4);
        probNNcut = num.Atof();
      }
    }
  }
  if (probNNghostCut < 0) {
    TString hn(histname);
    Ssiz_t p2 = hn.Index("ProbNNghost");
    if (p2 != kNPOS) {
      Ssiz_t lt = hn.Index('<', p2);
      if (lt != kNPOS) {
        TString num = hn(lt+1, 4);
        probNNghostCut = num.Atof();
      }
    }
  }
  if (trchi2ndofCut < 0) {
    TString hn(histname);
    Ssiz_t p3 = hn.Index("TRCHI2NDOF");
    if (p3 != kNPOS) {
      Ssiz_t lt = hn.Index('<', p3);
      if (lt != kNPOS) {
        TString num = hn(lt+1, 4);
        trchi2ndofCut = num.Atof();
      }
    }
  }

  TString cutsLabel;
  if (probNNcut > 0) cutsLabel += Form("ProbNN%s > %.2f", (isPion?"pi":"K"), probNNcut);
  if (probNNghostCut >= 0) {
    if (cutsLabel.Length()>0) cutsLabel += ", ";
    cutsLabel += Form("ProbNNghost < %.2f", probNNghostCut);
  }
  if (trchi2ndofCut >= 0) {
    if (cutsLabel.Length()>0) cutsLabel += ", ";
    cutsLabel += Form("TRCHI2NDOF < %.1f", trchi2ndofCut);
  }

  if (cutsLabel.Length()>0) tex.DrawLatexNDC(0.12, 0.86, cutsLabel);

  // nice colorbar sizing (optional)
  // c->Update();

  // save outputs
  TString base = TString(outdir.c_str()) + "/" + TString(histname);
  c->SaveAs(base + ".png");
  c->SaveAs(base + ".pdf");

  std::printf("Saved plots to: %s.png and %s.pdf\n", base.Data(), base.Data());
}