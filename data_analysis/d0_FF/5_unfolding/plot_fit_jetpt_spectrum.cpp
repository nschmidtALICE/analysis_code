// Plot jet-pt spectrum and fit with sum of two exponentials
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLatex.h>
#include <iostream>

void plot_fit_jetpt_spectrum()
{
    const char *inName = "jet_pt_spectrum.root";
    const char *histName = "hJetPt";
    TFile *f = TFile::Open(inName, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: cannot open " << inName << "\n";
        return;
    }

    TH1D *h = dynamic_cast<TH1D*>(f->Get(histName));
    if (!h) {
        std::cerr << "Error: histogram '" << histName << "' not found in " << inName << "\n";
        f->ls();
        f->Close();
        return;
    }

    gStyle->SetOptStat(1110);
    TCanvas *c = new TCanvas("c","Jet pT spectrum",800,600);
    c->SetLogy(0);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.8);
    h->SetLineWidth(1);
    h->Draw("E");

    // Fit function: A*exp(-x/p1) + B*exp(-x/p2)
    double xmin = 0.5; // avoid zero
    double xmax = h->GetXaxis()->GetXmax();
    TF1 *f2 = new TF1("f2","[0]*exp(-x/[1]) + [2]*exp(-x/[3])", xmin, xmax);
    // initial parameter guesses
    double peak = h->GetMaximum();
    f2->SetParameters(peak, 5.0, peak*0.1, 30.0);
    f2->SetParNames("A","p1","B","p2");

    // Perform fit in a reasonable range (exclude very high-pt if noisy)
    double fitMax = std::min(100.0, xmax);
    h->Fit(f2, "R", "", xmin, fitMax);

    // Draw fit overlay
    f2->SetLineColor(kRed);
    f2->SetLineWidth(2);
    f2->Draw("SAME");

    // Legend with fit parameters
    TLegend *leg = new TLegend(0.55,0.65,0.88,0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(h, "Data", "lep");
    leg->AddEntry(f2, "A e^{-x/p1} + B e^{-x/p2}", "l");
    leg->Draw();

    TLatex tex;
    tex.SetNDC();
    tex.SetTextSize(0.03);
    char buf[256];
    snprintf(buf, sizeof(buf), "p1 = %.3g, p2 = %.3g", f2->GetParameter(1), f2->GetParameter(3));
    tex.DrawLatex(0.15,0.85, buf);

    c->SaveAs("jet_pt_spectrum_fit.png");
    c->SaveAs("jet_pt_spectrum_fit.root");

    std::cout << "Fit results:\n";
    std::cout << " A = " << f2->GetParameter(0) << "  p1 = " << f2->GetParameter(1) << "\n";
    std::cout << " B = " << f2->GetParameter(2) << "  p2 = " << f2->GetParameter(3) << "\n";

    f->Close();
}
