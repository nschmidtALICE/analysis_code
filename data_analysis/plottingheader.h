#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TChain.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <TLine.h>
#include <TStyle.h>
#include <TString.h>
#include <TMath.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TLatex.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooCBShape.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooPolynomial.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <TSystem.h>
#include <TPaveText.h>
#include <RooHist.h>


using namespace std;


void makeLHCbStyleMassPlot(TH1F* h_Mass_Part, TString outputdir, TString collsys, TString xLabel = "K^{-}#pi^{+}", TString yLabel = "D^{0}", TString title = "D0") {
    // Load RooFit
    gSystem->Load("libRooFit");
    using namespace RooFit;
    
    // Create RooFit observable and dataset
    double part_mass = 1864.84;
    double masswindow = 60;
    if(title == "D0"){
        part_mass = 1864.84;
    } else if(title == "Ds"){
        part_mass = 1968.49;
    } else if(title == "Lb"){
        part_mass = 5621.0;
    } else if(title == "Lc"){
        part_mass = 2286.46;
    } else {
        cout << "Unknown particle title" << endl;
        return;
    }
    RooRealVar mass("mass", Form("m(%s)", xLabel.Data()), part_mass-masswindow, part_mass+masswindow, "MeV/c^{2}");
    RooDataHist data("data", Form("%s candidates", yLabel.Data()), RooArgList(mass), h_Mass_Part);
    
    // Create fit model: Crystal Ball + Gaussian
    // Crystal Ball parameters
    RooRealVar mean("mean", "Mean", part_mass, part_mass-10, part_mass+10);
    RooRealVar sigma1("sigma1", "Width of Gaussian core", 7, 3, 15);
    RooRealVar alpha("alpha", "Crystal Ball tail parameter", 1.5, 0.5, 5.0);
    RooRealVar n("n", "Crystal Ball power parameter", 2, 0.5, 10.0);
    
    // Gaussian parameters
    RooRealVar sigma2("sigma2", "Width of wide Gaussian", 12, 8, 25);
    RooRealVar frac("frac", "Fraction", 0.7, 0.0, 1.0);
    
    
    // Background parameters
    RooRealVar c0("c0", "c0", 0, -10, 10);
    RooRealVar c1("c1", "c1", -0.1, -10, 10);
    
    // Define functions
    RooCBShape crystalBall("crystalBall", "Crystal Ball", mass, mean, sigma1, alpha, n);
    RooGaussian gauss("gauss", "Gaussian", mass, mean, sigma2);
    RooAddPdf signal("signal", "Signal", RooArgList(crystalBall, gauss), RooArgList(frac));
    
    RooPolynomial background("background", "Background", mass, RooArgList(c0, c1));
    
    // Signal and background yields
    RooRealVar nsig("nsig", "Number of signal events", h_Mass_Part->Integral()*0.3, 0, h_Mass_Part->Integral()*2);
    RooRealVar nbkg("nbkg", "Number of background events", h_Mass_Part->Integral()*0.3, 0, h_Mass_Part->Integral()*2);
    
    // Combined model
    RooAddPdf model("model", "model", RooArgList(signal, background), RooArgList(nsig, nbkg));
    
    // Perform fit
    RooFitResult* fitResult = model.fitTo(data, Save(), Extended(), PrintLevel(-1));

    // create chi2
    // RooAbsReal * chi2_o_ndf = model.createChi2(data, Range("fullRange"), Extended(true), DataError(RooAbsData::Poisson));

    // // mass.setRange("signal",mean.getVal()-2*sigma2.getVal(),mean.getVal()+2*sigma2.getVal());
    // // RooAbsReal* fsigregion_model = model.createIntegral(mass,NormSet(mass),Range("signal")); //The "NormSet(x)" normalizes it to the total number of events to give you the fraction n_signal_region_events/n_total_events
    // // RooAbsReal* fsigregion_bkg = background.createIntegral(mass,NormSet(mass),Range("signal")); 
    // // Double_t nsigevents = fsigregion_model->getVal()*(nsig.getVal()+nbkg.getVal())-fsigregion_bkg->getVal()*nbkg.getVal(); 
    // // Double_t fsig = nsigevents/(fsigregion_model->getVal()*(nsig.getVal()+nbkg.getVal()));

    // std::cout << "Signal region fraction: " << fsig << std::endl;
    // std::cout << "Signal region events: " << nsigevents << std::endl;

    // Create the frame for plotting
    TCanvas* canvas = new TCanvas(Form("c_%s_MASS_%s_fit", title.Data(), collsys.Data()), Form("%s Mass in %s", yLabel.Data(), collsys.Data()), 900, 700);
    
    // Create a TPad for the mass plot
    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.05);
    pad1->Draw();
    pad1->cd();
    
    RooPlot* frame = mass.frame(Title(""));
    data.plotOn(frame, Name("data"));
    model.plotOn(frame, Name("model"));
    model.plotOn(frame, Name("signal"), Components("signal"), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(frame, Name("background"), Components("background"), LineStyle(kDashed), LineColor(kGreen+2));
    
    // Set axis titles and labels
    frame->GetXaxis()->SetTitle("");
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetTitle("Candidates / (1 MeV/c^{2})");
    frame->GetYaxis()->SetTitleOffset(1.2);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetLabelSize(0.05);
    frame->SetMaximum(frame->GetMaximum()*1.2);
    frame->Draw();
    
    // Add LHCb label
    TLatex* lhcbLabel = new TLatex();
    lhcbLabel->SetNDC();
    lhcbLabel->SetTextFont(42);
    lhcbLabel->SetTextSize(0.04);
    if(collsys == "pAr") lhcbLabel->DrawLatex(0.2, 0.85, Form("LHCb #it{p}Ar, %s #rightarrow %s",yLabel.Data(), xLabel.Data()));
    if(collsys == "pHe") lhcbLabel->DrawLatex(0.2, 0.85, Form("LHCb #it{p}He, %s #rightarrow %s",yLabel.Data(), xLabel.Data()));
    if(collsys == "pH") lhcbLabel->DrawLatex(0.2, 0.85, Form("LHCb #it{p}H, %s #rightarrow %s",yLabel.Data(), xLabel.Data()));

    // lhcbLabel->DrawLatex(0.2, 0.85, "LHCb #it{p}Ar, %s #rightarrow %s",yLabel.Data(), xLabel.Data()));
    lhcbLabel->DrawLatex(0.2, 0.8, "#sqrt{s_{NN}} = 113 GeV");

    // // Add collision system info below the lhcb label, snn = 113 GeV
    // TLatex* collisionSystemLabel = new TLatex();
    // collisionSystemLabel->SetNDC();
    // collisionSystemLabel->SetTextFont(42);
    // collisionSystemLabel->SetTextSize(0.04);
    // collisionSystemLabel->DrawLatex(0.2, 0.8, "#sqrt{s_{NN}} = 113 GeV");
    // // collisionSystemLabel->DrawLatex(0.2, 0.75, "L_{int} = 0.5 nb^{-1}");

    
    // Add fit results to the plot
    double meanVal = mean.getVal();
    double meanErr = mean.getError();
    double sigmaVal = sigma1.getVal();
    double sigmaErr = sigma1.getError();
    double nSigVal = nsig.getVal();
    double nSigErr = nsig.getError();
    double nBckgVal = nbkg.getVal();
    double nBckgErr = nbkg.getError();
    // double chi2ndf = frame->chiSquare();
    // double nfitparam = fitResult->floatParsFinal().getSize();
    // std::cout << "Chi2/ndf = " << chi2ndf/nfitparam << std::endl;
    

    // Int_t npar = signal.getParameters(data)->selectByAttrib("Constant",kFALSE)->getSize();
    // Double_t chi2ndf = frame->chiSquare(npar); // calculate chi2/ndf (model and data selected from last plotted, but can be specified)

    TPaveText* pave = new TPaveText(0.65, 0.70, 0.88, 0.90, "NDC");
    pave->SetBorderSize(0);
    pave->SetFillStyle(0);
    pave->SetTextAlign(12);
    pave->SetTextFont(42);
    pave->SetTextSize(0.04);
    // pave->AddText(Form("#chi^{2}/ndf = %.2f", chi2ndf));
    pave->AddText(Form("Mass = %.1f #pm %.1f MeV/c^{2}", meanVal, meanErr));
    pave->AddText(Form("#sigma = %.1f #pm %.1f MeV/c^{2}", sigmaVal, sigmaErr));
    pave->AddText(Form("N(%s) = %.0f #pm %.0f", yLabel.Data(), nSigVal, nSigErr));
    pave->AddText(Form("N(bkg) = %.0f #pm %.0f", nBckgVal, nBckgErr));
    pave->Draw();

    //add a legend for data, sig+bg, signal, and background
    TLegend *leg = new TLegend(0.2, 0.6, 0.4, 0.8);
    //make the legend to not have a frame
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(frame->findObject("data"), "Data", "lep");
    leg->AddEntry(frame->findObject("model"), "Sig+bkg", "l");
    leg->AddEntry(frame->findObject("signal"), "Signal", "l");
    leg->AddEntry(frame->findObject("background"), "Background", "l");
    leg->Draw();
    
    // Go back to the main canvas to create the residuals pad
    canvas->cd();
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.05);
    pad2->Draw();
    pad2->cd();
    

// Make sure the model is plotted again to create a curve that can be compared with data
model.plotOn(frame, Name("model_for_pull"));

// Create pull histogram without debugging output
RooHist* pullHist = frame->pullHist("data", "model_for_pull");

// Create a new frame for the pulls with explicit range
RooPlot* pullFrame = mass.frame(Title(""));
pullFrame->addPlotable((RooPlotable*)pullHist, "P");

// Make points more visible
pullHist->SetMarkerStyle(20);  // Filled circle
pullHist->SetMarkerSize(0.8);  // Larger marker
pullHist->SetMarkerColor(kBlack);
pullHist->SetLineColor(kBlack);

// Set axis titles and labels for the pull plot
pullFrame->GetXaxis()->SetTitle(Form("m(%s) [MeV/c^{2}]", xLabel.Data()));
pullFrame->GetXaxis()->SetTitleSize(0.12);
pullFrame->GetXaxis()->SetLabelSize(0.10);
pullFrame->GetXaxis()->SetTitleOffset(1.0);

pullFrame->GetYaxis()->SetTitle("Pull");
pullFrame->GetYaxis()->SetTitleSize(0.12);
pullFrame->GetYaxis()->SetLabelSize(0.10);
pullFrame->GetYaxis()->SetTitleOffset(0.4);
pullFrame->GetYaxis()->SetNdivisions(505);

// Set more appropriate range for pulls (using smaller range to see detail)
pullFrame->SetMinimum(-3);
pullFrame->SetMaximum(3);
pullFrame->Draw();

// Draw horizontal lines at -2, 0, and 2 for reference
TLine* line0 = new TLine(part_mass-masswindow, 0, part_mass+masswindow, 0);
line0->SetLineStyle(1);
line0->SetLineColor(kRed);
line0->SetLineWidth(2);
line0->Draw("same");

TLine* lineP2 = new TLine(part_mass-masswindow, 2, part_mass+masswindow, 2);
lineP2->SetLineStyle(2);
lineP2->SetLineColor(kRed);
lineP2->Draw("same");

TLine* lineM2 = new TLine(part_mass-masswindow, -2, part_mass+masswindow, -2);
lineM2->SetLineStyle(2);
lineM2->SetLineColor(kRed);
lineM2->Draw("same");
// cout << "chi2_o_ndf: " << chi2_o_ndf->getVal() << endl;

    // Save the canvas
    canvas->SaveAs(Form("%s%s_MASS_%s_fit.pdf", outputdir.Data(), title.Data(), collsys.Data()));
    canvas->SaveAs(Form("%s%s_MASS_%s_fit.png", outputdir.Data(), title.Data(), collsys.Data()));
    
    // Print fit results
    std::cout << "\nFit results for " << title.Data() << " mass in " << collsys.Data() << ":" << std::endl;
    std::cout << "Mean = " << meanVal << " ± " << meanErr << " MeV/c²" << std::endl;
    std::cout << "Sigma = " << sigmaVal << " ± " << sigmaErr << " MeV/c²" << std::endl;
    std::cout << "N(" << title.Data() << ") = " << nSigVal << " ± " << nSigErr << std::endl;
    std::cout << "Fraction in CB = " << frac.getVal() << " ± " << frac.getError() << std::endl;
    std::cout << "Alpha = " << alpha.getVal() << " ± " << alpha.getError() << std::endl;
    std::cout << "n = " << n.getVal() << " ± " << n.getError() << std::endl;
}
