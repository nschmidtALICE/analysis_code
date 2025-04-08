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
#include <TROOT.h>
#include <RooAddPdf.h>
#include <RooPolynomial.h>
#include <RooExponential.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooFit.h>
#include <RooGlobalFunc.h>
#include <RooFFTConvPdf.h>
#include <RooDataSet.h>
#include <RooProdPdf.h>
#include <RooGenericPdf.h>
#include <RooWorkspace.h>
#include <RooCategory.h>
#include <RooGaussModel.h>
#include <RooDecay.h>
#include <RooAddModel.h>
#include <RooBDecay.h>
#include <RooSimultaneous.h>
#include <RooChi2Var.h>
#include <RooHistPdf.h>
#include <RooConstVar.h>
#include <RooFormulaVar.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <TSystem.h>
#include <TPaveText.h>
#include <RooHist.h>
#include <RooStats/SPlot.h>
using namespace std;

std::vector<double> makeLHCbStyleMassPlot(TH1F *h_Mass_Part, TString outputdir, TString collsys = "pAr", TString xLabel = "K^{-}#pi^{+}", TString yLabel = "D^{0}", TString title = "D0")
{

    // Load RooFit
    gSystem->Load("libRooFit");
    using namespace RooFit;

    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    // Create RooFit observable and dataset
    double part_mass = 1864.84;
    double masswindow = 60;
    double nSigRange = 0.3;
    double sigma1_default = 7;
    double mean_spread = 10;
    if (title == "D0")
    {
        part_mass = 1864.84;
        nSigRange = 0.6;
    }
    else if (title == "Ds")
    {
        part_mass = 1968.49;
    }
    else if (title == "Lb")
    {
        part_mass = 5621.0;
    }
    else if (title == "Lc")
    {
        part_mass = 2286.46;
        sigma1_default = 6;
        mean_spread = 5;
    }
    else if (title == "X3872")
    {
        part_mass = 3871.69;
    }
    else
    {
        cout << "Unknown particle title" << endl;
        std::vector<double> resultsErr = {0, 0, 0, 0};
        return resultsErr;
    }
    RooRealVar mass("mass", Form("m(%s)", xLabel.Data()), part_mass - masswindow, part_mass + masswindow, "MeV/c^{2}");
    RooDataHist data("data", Form("%s candidates", yLabel.Data()), RooArgList(mass), h_Mass_Part);

    // Create fit model: Crystal Ball + Gaussian
    // Crystal Ball parameters
    RooRealVar mean("mean", "Mean", part_mass, part_mass - mean_spread, part_mass + mean_spread);
    RooRealVar sigma1("sigma1", "Width of Gaussian core", sigma1_default, 3, 15);
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
    RooRealVar nsig("nsig", "Number of signal events", h_Mass_Part->Integral() * nSigRange, 0, h_Mass_Part->Integral() * 2);
    RooRealVar nbkg("nbkg", "Number of background events", h_Mass_Part->Integral() * 0.3, 0, h_Mass_Part->Integral() * 2);

    // Combined model
    RooAddPdf model("model", "model", RooArgList(signal, background), RooArgList(nsig, nbkg));

    // Perform fit
    RooFitResult *fitResult = model.fitTo(data, Save(), Extended(), PrintLevel(4));

    // create chi2
    // RooAbsReal * chi2_o_ndf = model.createChi2(data, Range("fullRange"), Extended(true), DataError(RooAbsData::Poisson));

    // Create the frame for plotting
    TCanvas *canvas = new TCanvas(Form("c_%s_MASS_%s_fit", title.Data(), collsys.Data()), Form("%s Mass in %s", yLabel.Data(), collsys.Data()), 900, 700);

    // Create a TPad for the mass plot
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.05);
    pad1->Draw();
    pad1->cd();

    RooPlot *frame = mass.frame(Title(""));
    data.plotOn(frame, Name("data"));
    model.plotOn(frame, Name("model"));
    model.plotOn(frame, Name("signal"), Components("signal"), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(frame, Name("background"), Components("background"), LineStyle(kDashed), LineColor(kGreen + 2));

    // Set axis titles and labels
    frame->GetXaxis()->SetTitle("");
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetTitle("Candidates / (1 MeV/c^{2})");
    frame->GetYaxis()->SetTitleOffset(1.2);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetLabelSize(0.05);
    frame->SetMaximum(frame->GetMaximum() * 1.2);
    frame->Draw();

    // Add LHCb label
    TLatex *lhcbLabel = new TLatex();
    lhcbLabel->SetNDC();
    lhcbLabel->SetTextFont(42);
    lhcbLabel->SetTextSize(0.04);
    if (collsys.Contains("pAr"))
        lhcbLabel->DrawLatex(0.2, 0.85, Form("LHCb #it{p}Ar, %s #rightarrow %s", yLabel.Data(), xLabel.Data()));
    if (collsys.Contains("pHe"))
        lhcbLabel->DrawLatex(0.2, 0.85, Form("LHCb #it{p}He, %s #rightarrow %s", yLabel.Data(), xLabel.Data()));
    if (collsys.Contains("pH"))
        lhcbLabel->DrawLatex(0.2, 0.85, Form("LHCb #it{p}H, %s #rightarrow %s", yLabel.Data(), xLabel.Data()));

    lhcbLabel->DrawLatex(0.2, 0.8, "#sqrt{s_{NN}} = 113 GeV");

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

    TPaveText *pave = new TPaveText(0.65, 0.70, 0.88, 0.90, "NDC");
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

    // add a legend for data, sig+bg, signal, and background
    TLegend *leg = new TLegend(0.2, 0.55, 0.4, 0.75);
    // make the legend to not have a frame
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(frame->findObject("data"), "Data", "lep");
    leg->AddEntry(frame->findObject("model"), "Sig+bkg", "l");
    leg->AddEntry(frame->findObject("signal"), "Signal", "l");
    leg->AddEntry(frame->findObject("background"), "Background", "l");
    leg->Draw();

    // Go back to the main canvas to create the residuals pad
    canvas->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.05);
    pad2->Draw();
    pad2->cd();

    // Make sure the model is plotted again to create a curve that can be compared with data
    model.plotOn(frame, Name("model_for_pull"));

    // Create pull histogram without debugging output
    RooHist *pullHist = frame->pullHist("data", "model_for_pull");

    // Create a new frame for the pulls with explicit range
    RooPlot *pullFrame = mass.frame(Title(""));
    pullFrame->addPlotable((RooPlotable *)pullHist, "P");

    // Make points more visible
    pullHist->SetMarkerStyle(20); // Filled circle
    pullHist->SetMarkerSize(0.8); // Larger marker
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
    TLine *line0 = new TLine(part_mass - masswindow, 0, part_mass + masswindow, 0);
    line0->SetLineStyle(1);
    line0->SetLineColor(kRed);
    line0->SetLineWidth(2);
    line0->Draw("same");

    TLine *lineP2 = new TLine(part_mass - masswindow, 2, part_mass + masswindow, 2);
    lineP2->SetLineStyle(2);
    lineP2->SetLineColor(kRed);
    lineP2->Draw("same");

    TLine *lineM2 = new TLine(part_mass - masswindow, -2, part_mass + masswindow, -2);
    lineM2->SetLineStyle(2);
    lineM2->SetLineColor(kRed);
    lineM2->Draw("same");
    // cout << "chi2_o_ndf: " << chi2_o_ndf->getVal() << endl;

    // Save the canvas
    canvas->SaveAs(Form("%s%s_MASS_%s_fit.pdf", outputdir.Data(), title.Data(), collsys.Data()));
    // canvas->SaveAs(Form("%s%s_MASS_%s_fit.png", outputdir.Data(), title.Data(), collsys.Data()));

    // Print fit results
    std::cout << "\nFit results for " << title.Data() << " mass in " << collsys.Data() << ":" << std::endl;
    std::cout << "Mean = " << meanVal << " ± " << meanErr << " MeV/c²" << std::endl;
    std::cout << "Sigma = " << sigmaVal << " ± " << sigmaErr << " MeV/c²" << std::endl;
    std::cout << "N(" << title.Data() << ") = " << nSigVal << " ± " << nSigErr << std::endl;
    std::cout << "Fraction in CB = " << frac.getVal() << " ± " << frac.getError() << std::endl;
    std::cout << "Alpha = " << alpha.getVal() << " ± " << alpha.getError() << std::endl;
    std::cout << "n = " << n.getVal() << " ± " << n.getError() << std::endl;

    std::vector<double> results = {nSigVal, nSigErr, nBckgVal, nBckgErr};
    return results;
}

std::vector<double> makeLHCbStyleMassPlotFromTree(TFile *inputfile, TString treeName, TString massVarName, TString outputdir, TString selection = "",
                                                  TString collsys = "pAr", TString xLabel = "K^{-}#pi^{+}", TString yLabel = "D^{0}",
                                                  TString title = "D0", bool useWeights = false, TString weightVar = "")
{
    // Load RooFit
    gSystem->Load("libRooFit");
    using namespace RooFit;
    // RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    // Check input file
    if (!inputfile || inputfile->IsZombie())
    {
        std::cerr << "Error: Invalid input file!" << std::endl;
        std::vector<double> resultsErr = {0, 0, 0, 0};
        return resultsErr;
    }

    // Get the tree from the file
    TTree *tree = (TTree *)inputfile->Get(treeName);
    if (!tree)
    {
        std::cerr << "Error: Could not find tree '" << treeName << "' in the input file!" << std::endl;
        std::vector<double> resultsErr = {0, 0, 0, 0};
        return resultsErr;
    }

    // Create output directory if it doesn't exist
    gSystem->Exec("mkdir -p " + outputdir);

    // Set particle properties
    double part_mass = 1864.84;
    double masswindow = 60;
    double nSigRange = 0.3;
    double mean_spread = 10;
    double sigma1_default = 7;
    double sigma1_upper = 15;
    double sigma1_lower = 3;
    double alpha_default = 1.5;
    double alpha_lower = 0.5;
    double alpha_upper = 5.0;
    double n_default = 2;
    double n_lower = 0.5;
    double n_upper = 10.0;
    double sigma2_default = 12;
    double sigma2_lower = 8;
    double sigma2_upper = 25;
    double frac_default = 0.7;
    double frac_lower = 0.0;
    double frac_upper = 1.0;
    double c0_default = 0;
    double c0_lower = -10;
    double c0_upper = 10;
    double c1_default = -0.1;
    double c1_lower = -10;
    double c1_upper = 10;

    if (title == "D0")
    {
        part_mass = 1864.84;
        masswindow = 60;
        if (collsys.Contains("pNe"))
        {
            alpha_default = 1.88;
            c0_default = 7.8;
            c1_default = -0.00316;
            frac_default = 0.505;
            n_default = 3.57;
            sigma1_default = 6.86;
            sigma2_default = 10.33;
        }
    }
    else if (title == "Ds")
    {
        part_mass = 1968.49;
    }
    else if (title == "Lb")
    {
        part_mass = 5621.0;
    }
    else if (title == "Lc")
    {
        part_mass = 2286.46;
        sigma1_default = 6;
        mean_spread = 5;
    }
    else if (title == "X3872")
    {
        part_mass = 3871.69;
    }
    else
    {
        cout << "Unknown particle title" << endl;
        std::vector<double> resultsErr = {0, 0, 0, 0};
        return resultsErr;
    }

    // Create RooFit observable
    RooRealVar mass(massVarName, Form("m(%s)", xLabel.Data()), part_mass - masswindow, part_mass + masswindow, "MeV/c^{2}");

    // Add D0_BPVIPCHI2 variable to the observables
    RooRealVar ipchi2("D0_BPVIPCHI2", "D0_BPVIPCHI2", 0, 10000); // Add appropriate range

    // Create a RooDataSet from the tree including D0_BPVIPCHI2
    RooArgSet observables(mass, ipchi2); // Now includes both mass and IPCHI2

    RooDataSet *data = nullptr;

    if (useWeights && !weightVar.IsNull())
    {
        // Create a weight variable if needed
        RooRealVar weight(weightVar, "weight", 0, 1000);
        observables.add(weight);
        data = new RooDataSet("data", Form("%s candidates", yLabel.Data()),
                              tree, observables, selection, weightVar);
    }
    else
    {
        data = new RooDataSet("data", Form("%s candidates", yLabel.Data()),
                              tree, observables, selection);
    }

    if (data->numEntries() == 0)
    {
        std::cerr << "Error: Empty dataset after selection!" << std::endl;
        delete data;
        std::vector<double> resultsErr = {0, 0, 0, 0};
        return resultsErr;
    }

    std::cout << "Created dataset with " << data->numEntries() << " entries" << std::endl;

    // Use approximate count for estimating signal and background
    Long64_t totalEntries = data->numEntries();

    // Create fit model: Crystal Ball + Gaussian
    // Crystal Ball parameters
    RooRealVar mean("mean", "Mean", part_mass, part_mass - mean_spread, part_mass + mean_spread);
    RooRealVar sigma1("sigma1", "Width of Gaussian core", sigma1_default, sigma1_lower, sigma1_upper);
    RooRealVar alpha("alpha", "Crystal Ball tail parameter", alpha_default, alpha_lower, alpha_upper);
    RooRealVar n("n", "Crystal Ball power parameter", n_default, n_lower, n_upper);

    // Gaussian parameters
    RooRealVar sigma2("sigma2", "Width of wide Gaussian", sigma2_default, sigma2_lower, sigma2_upper);
    RooRealVar frac("frac", "Fraction", frac_default, frac_lower, frac_upper);

    // Background parameters
    RooRealVar c0("c0", "c0", c0_default, c0_lower, c0_upper);
    RooRealVar c1("c1", "c1", c1_default, c1_lower, c1_upper);

    // Define functions
    RooCBShape crystalBall("crystalBall", "Crystal Ball", mass, mean, sigma1, alpha, n);
    RooGaussian gauss("gauss", "Gaussian", mass, mean, sigma2);
    RooAddPdf signal("signal", "Signal", RooArgList(crystalBall, gauss), RooArgList(frac));

    RooPolynomial background("background", "Background", mass, RooArgList(c0, c1));

    // Signal and background yields
    RooRealVar nsig("nsig", "Number of signal events", totalEntries * nSigRange, totalEntries * 0.1, totalEntries * 2);
    RooRealVar nbkg("nbkg", "Number of background events", totalEntries * 0.3, 0, totalEntries * 2);

    // Combined model
    RooAddPdf model("model", "model", RooArgList(signal, background), RooArgList(nsig, nbkg));

    cout << "set up model" << endl;
    // Perform fit
    RooFitResult *fitResult = model.fitTo(*data, Save(), Extended(), PrintLevel(1));

    // Create the frame for plotting
    TCanvas *canvas = new TCanvas(Form("c_%s_MASS_%s_fit", title.Data(), collsys.Data()),
                                  Form("%s Mass in %s", yLabel.Data(), collsys.Data()), 900, 700);

    // Create a TPad for the mass plot
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.05);
    pad1->Draw();
    pad1->cd();

    RooPlot *frame = mass.frame(Title(""));
    data->plotOn(frame, Name("data"));
    model.plotOn(frame, Name("model"));
    model.plotOn(frame, Name("signal"), Components("signal"), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(frame, Name("background"), Components("background"), LineStyle(kDashed), LineColor(kGreen + 2));

    // Set axis titles and labels
    frame->GetXaxis()->SetTitle("");
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetTitle("Candidates / (1 MeV/c^{2})");
    frame->GetYaxis()->SetTitleOffset(1.2);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetLabelSize(0.05);
    frame->SetMaximum(frame->GetMaximum() * 1.2);
    frame->Draw();

    // Add LHCb label
    TLatex *lhcbLabel = new TLatex();
    lhcbLabel->SetNDC();
    lhcbLabel->SetTextFont(42);
    lhcbLabel->SetTextSize(0.04);
    if (collsys.Contains("pAr"))
        lhcbLabel->DrawLatex(0.2, 0.85, Form("LHCb #it{p}Ar, %s #rightarrow %s", yLabel.Data(), xLabel.Data()));
    else if (collsys.Contains("pHe"))
        lhcbLabel->DrawLatex(0.2, 0.85, Form("LHCb #it{p}He, %s #rightarrow %s", yLabel.Data(), xLabel.Data()));
    else if (collsys.Contains("pH"))
        lhcbLabel->DrawLatex(0.2, 0.85, Form("LHCb #it{p}H, %s #rightarrow %s", yLabel.Data(), xLabel.Data()));
    else if (collsys.Contains("pNe"))
        lhcbLabel->DrawLatex(0.2, 0.85, Form("LHCb #it{p}Ne, %s #rightarrow %s", yLabel.Data(), xLabel.Data()));

    lhcbLabel->DrawLatex(0.2, 0.8, "#sqrt{s_{NN}} = 113 GeV");

    // Add fit results to the plot
    double meanVal = mean.getVal();
    double meanErr = mean.getError();
    double sigmaVal = sigma1.getVal();
    double sigmaErr = sigma1.getError();
    double nSigVal = nsig.getVal();
    double nSigErr = nsig.getError();
    double nBckgVal = nbkg.getVal();
    double nBckgErr = nbkg.getError();

    TPaveText *pave = new TPaveText(0.65, 0.70, 0.88, 0.90, "NDC");
    pave->SetBorderSize(0);
    pave->SetFillStyle(0);
    pave->SetTextAlign(12);
    pave->SetTextFont(42);
    pave->SetTextSize(0.04);
    pave->AddText(Form("Mass = %.1f #pm %.1f MeV/c^{2}", meanVal, meanErr));
    pave->AddText(Form("#sigma = %.1f #pm %.1f MeV/c^{2}", sigmaVal, sigmaErr));
    pave->AddText(Form("N(%s) = %.0f #pm %.0f", yLabel.Data(), nSigVal, nSigErr));
    pave->AddText(Form("N(bkg) = %.0f #pm %.0f", nBckgVal, nBckgErr));
    pave->Draw();

    // Add a legend for data, sig+bg, signal, and background
    TLegend *leg = new TLegend(0.2, 0.55, 0.4, 0.75);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(frame->findObject("data"), "Data", "lep");
    leg->AddEntry(frame->findObject("model"), "Sig+bkg", "l");
    leg->AddEntry(frame->findObject("signal"), "Signal", "l");
    leg->AddEntry(frame->findObject("background"), "Background", "l");
    leg->Draw();

    // Go back to the main canvas to create the residuals pad
    canvas->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.05);
    pad2->Draw();
    pad2->cd();

    // Make sure the model is plotted again to create a curve that can be compared with data
    model.plotOn(frame, Name("model_for_pull"));

    // Create pull histogram
    RooHist *pullHist = frame->pullHist("data", "model_for_pull");

    // Create a new frame for the pulls with explicit range
    RooPlot *pullFrame = mass.frame(Title(""));
    pullFrame->addPlotable((RooPlotable *)pullHist, "P");

    // Make points more visible
    pullHist->SetMarkerStyle(20); // Filled circle
    pullHist->SetMarkerSize(0.8); // Larger marker
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

    // Set more appropriate range for pulls
    pullFrame->SetMinimum(-3);
    pullFrame->SetMaximum(3);
    pullFrame->Draw();

    // Draw horizontal lines at -2, 0, and 2 for reference
    TLine *line0 = new TLine(part_mass - masswindow, 0, part_mass + masswindow, 0);
    line0->SetLineStyle(1);
    line0->SetLineColor(kRed);
    line0->SetLineWidth(2);
    line0->Draw("same");

    TLine *lineP2 = new TLine(part_mass - masswindow, 2, part_mass + masswindow, 2);
    lineP2->SetLineStyle(2);
    lineP2->SetLineColor(kRed);
    lineP2->Draw("same");

    TLine *lineM2 = new TLine(part_mass - masswindow, -2, part_mass + masswindow, -2);
    lineM2->SetLineStyle(2);
    lineM2->SetLineColor(kRed);
    lineM2->Draw("same");

    // Save the canvas
    canvas->SaveAs(Form("%s%s_MASS_%s_fit.pdf", outputdir.Data(), title.Data(), collsys.Data()));

    // Print fit results
    std::cout << "\nFit results for " << title.Data() << " mass in " << collsys.Data() << ":" << std::endl;
    std::cout << "Mean = " << meanVal << " ± " << meanErr << " MeV/c²" << std::endl;
    std::cout << "Sigma = " << sigmaVal << " ± " << sigmaErr << " MeV/c²" << std::endl;
    std::cout << "N(" << title.Data() << ") = " << nSigVal << " ± " << nSigErr << std::endl;
    std::cout << "Fraction in CB = " << frac.getVal() << " ± " << frac.getError() << std::endl;
    std::cout << "Alpha = " << alpha.getVal() << " ± " << alpha.getError() << std::endl;
    std::cout << "n = " << n.getVal() << " ± " << n.getError() << std::endl;

    // Option to create sPlot weights
    if (useWeights)
    {
        // Create and save sPlot weights
        std::cout << "Creating sPlot weights..." << std::endl;

        // The SPlot constructor requires yield variables
        RooStats::SPlot *sPlot = new RooStats::SPlot("sPlot", "sPlot", *data, &model,
                                                     RooArgList(nsig, nbkg));

        // Save the sPlot tree to a file
        TFile *sPlotFile = new TFile(Form("%s%s_%s_sPlot.root", outputdir.Data(),
                                          title.Data(), collsys.Data()),
                                     "RECREATE");
        // Check if the sPlot creation was successful
        if (sPlot)
        {
            std::cout << "SPlot creation successful!" << std::endl;

            // Create a new tree with sPlot weights
            TTree *sPlotTree = new TTree("sPlotTree", "Tree with sPlot weights");

            // Variables to hold data
            double massVal;
            double ipchi2Val; // Value from dataset
            double sigWeight, bkgWeight;

            // Create branches
            sPlotTree->Branch(massVarName, &massVal, massVarName + "/D");
            sPlotTree->Branch("D0_BPVIPCHI2", &ipchi2Val, "D0_BPVIPCHI2/F");
            sPlotTree->Branch("sigWeight", &sigWeight, "sigWeight/D");
            sPlotTree->Branch("bkgWeight", &bkgWeight, "bkgWeight/D");

            // Use hardcoded sWeight names - these are standard in ROOT
            TString sigWeightName = "nsig_sw";  // Standard naming convention
            TString bkgWeightName = "nbkg_sw"; // Standard naming convention

            std::cout << "Using signal weight name: " << sigWeightName << std::endl;
            std::cout << "Using background weight name: " << bkgWeightName << std::endl;

            // Fill the tree with sPlot weights and IPChi2 values
            for (int i = 0; i < data->numEntries(); i++)
            {
                const RooArgSet *row = data->get(i);
                massVal = row->getRealValue(massVarName);
                ipchi2Val = row->getRealValue("D0_BPVIPCHI2"); // Get IPCHI2 directly from dataset

                try
                {
                    sigWeight = sPlot->GetSWeight(i, sigWeightName.Data());
                    bkgWeight = sPlot->GetSWeight(i, bkgWeightName.Data());
                    sPlotTree->Fill();
                }
                catch (const std::exception &e)
                {
                    std::cerr << "Error getting sWeights for entry " << i << ": " << e.what() << std::endl;
                }
            }

            // Save original dataset variables for reference - FIX HERE
            // Create a non-const copy we can modify
            RooArgSet *vars = (RooArgSet *)data->get()->Clone("variables");
            vars->Write("variables");

            // Copy the original dataset instead of creating a new one
            RooDataSet *dataCopy = (RooDataSet *)data->Clone("sData");
            dataCopy->Write();

            // Save the tree and close file
            sPlotTree->Write();
            sPlotFile->Close();

            std::cout << "sPlot weights saved to " << Form("%s%s_%s_sPlot.root", outputdir.Data(),
                                                           title.Data(), collsys.Data()) << std::endl;

            // Clean up
            delete sPlotFile;
            delete sPlotTree;
            delete vars; // Clean up the cloned variables
        }
        else
        {
            std::cerr << "Error: sPlot creation failed!" << std::endl;
        }

        delete sPlot;
    }

    // Store and return results
    std::vector<double> results = {nSigVal, nSigErr, nBckgVal, nBckgErr, meanVal, meanErr, sigmaVal, sigmaErr,
                                   alpha.getVal(), alpha.getError(), n.getVal(), n.getError()};

    // Cleanup
    delete canvas;
    delete data;

    return results;
}

// Fit a histogram of log(BPVIPCHI2) with Crystal Ball + Restricted Gaussian
std::vector<double> fitLogIPChi2Distribution(TH1F *histo, TString outputdir, TString collsys = "pAr",
    TString particle = "D^{0}", TString particleType = "D0",
    bool isWeighted = false)
{
    using namespace RooFit;
    gSystem->Load("libRooFit");

    // Set output path
    gSystem->Exec("mkdir -p " + outputdir);

    // Create variables
    double minIPchi2 = -5.0;
    double maxIPchi2 = 4.0;
    if (particleType == "D0")
    {
        if (collsys == "pNe")
        {
            if(isWeighted){
                minIPchi2 = -7.0;
                maxIPchi2 = 8.0;
            }
        }
    }
    if (particleType == "Lc")
    {
        minIPchi2 = -3.0;
        maxIPchi2 = 6.0;
    }
    RooRealVar logIPChi2("logIPChi2", "log(#chi^{2}_{IP})", minIPchi2, maxIPchi2);

    // Create dataset from histogram
    RooDataHist data("data", "log IP chi2 data", RooArgList(logIPChi2), histo);
    double maxcounts = histo->GetMaximum() * 1.5;
    // Crystal Ball parameters for the prompt component (mainly for logIPChi2 < 0)
    // create doubles for all base parameters
    double mean_val = 0.55;
    double mean_min = 0.4;
    double mean_max = 0.6;
    double sigma_val = 0.4;
    double sigma_min = 0.1;
    double sigma_max = 0.6;
    double alpha_val = 1.0;
    double alpha_min = 0.1;
    double alpha_max = 1.1;
    double n_val = 130.;
    double n_min = 100.;
    double n_max = 200.;

    // create doubles for gaussian parameters
    double gmean_val = 1.9;
    double gsigma_val = 0.4;
    double gmean_min = 1.7;
    double gmean_max = 2.2;
    double gsigma_min = 0.3;
    double gsigma_max = 0.5;

    // create double for fraction
    double fraction_val = 0.9;
    double fraction_min = 0.8;
    double fraction_max = 1.0;
    if (particleType == "D0")
    {
        if (collsys == "pNe")
        {
            // sigma_max = 0.6;
            alpha_max = 1.5;
            if(isWeighted){
                fraction_min = 0.7;
                mean_val = 1.1;
                mean_min = 0.7;
                mean_max = 2.0;
                sigma_val = 0.9;
                sigma_min = 0.5;
                sigma_max = 1.5;
                n_val = 130;
                n_min = 100.;
                n_max = 200.;

                gmean_val = 3.5;
                gmean_min = 2.0;
                gmean_max = 5.5;
                gsigma_val = 1.0;
                gsigma_min = 0.4;
                gsigma_max = 1.5;
            }
        }
    }
    if (particleType == "Lc")
    {
        mean_val = 1.5;
        mean_min = 1.4;
        mean_max = 1.6;
        if (collsys == "pNe")
        {
            mean_val = 1.6;
            mean_min = 1.5;
            mean_max = 1.8;
        }
        sigma_val = 0.5;
        sigma_min = 0.4;
        sigma_max = 0.6;
        alpha_val = 0.85;
        alpha_min = 0.75;
        alpha_max = 1.0;
        n_val = 135;
        n_min = 120.;
        n_max = 150.;
        gmean_val = 2.9;
        gmean_min = 2.8;
        gmean_max = 3.1;
        gsigma_val = 0.5;
        gsigma_min = 0.4;
        gsigma_max = 0.6;
    }

    RooRealVar mean("mean", "mean", mean_val, mean_min, mean_max);
    RooRealVar sigma("sigma", "sigma", sigma_val, sigma_min, sigma_max);
    RooRealVar alpha("alpha", "alpha", alpha_val, alpha_min, alpha_max);
    RooRealVar n("n", "n", n_val, n_min, n_max);
    RooCBShape crystalBall("crystalBall", "Crystal Ball", logIPChi2, mean, sigma, alpha, n);

    // Gaussian parameters for the nonprompt component (restricted to logIPChi2 > 0)
    RooRealVar gmean("gmean", "Gaussian mean", gmean_val, gmean_min, gmean_max);
    RooRealVar gsigma("gsigma", "Gaussian sigma", gsigma_val, gsigma_min, gsigma_max);
    RooGaussian fullGauss("fullGauss", "Full Gaussian", logIPChi2, gmean, gsigma);

    // // Create a step function to restrict the Gaussian to logIPChi2 > 0
    // RooFormulaVar step("step", "step(logIPChi2)", "logIPChi2 > 0.0", RooArgList(logIPChi2));

    // Create the partial Gaussian (only defined for logIPChi2 > 0)
    RooProdPdf gauss("gauss", "Partial Gaussian (x > 0)", RooArgList(fullGauss));

    // Fraction for combining the components
    RooRealVar fraction("fraction", "CB fraction", fraction_val, fraction_min, fraction_max);

    // Create combined model
    RooAddPdf model("model", "Crystal Ball + Restricted Gaussian", RooArgList(crystalBall, gauss), RooArgList(fraction));

    // Perform the fit
    RooFitResult *fitResult;
    if (isWeighted) {
        // For weighted data from sPlot, use SumW2Error(true)
        fitResult = model.fitTo(data, Save(), PrintLevel(1), SumW2Error(true));
        std::cout << "Fit performed with SumW2Error(true) for weighted data" << std::endl;
    } else {
        // For normal histograms (not weighted)
        fitResult = model.fitTo(data, Save(), PrintLevel(1));
    }
    std::cout << "Fit status for " << collsys << ": " << fitResult->status() << std::endl;

    // Create canvas for plotting
    TCanvas *canvas = new TCanvas("c_fit_logIPChi2_" + collsys, "Log IP Chi2 Fit", 800, 800);
    canvas->SetLeftMargin(0.15);
    canvas->Divide(1, 2);

    // Main panel for the fit
    TPad *mainPad = (TPad *)canvas->GetPad(1);
    mainPad->SetPad(0, 0.3, 1, 1.0);
    mainPad->SetBottomMargin(0.03);
    mainPad->SetLeftMargin(0.15);
    mainPad->SetRightMargin(0.05);
    mainPad->cd();

    // Create frame and plot data
    RooPlot *frame = logIPChi2.frame(Title(""));
    data.plotOn(frame, Name("data"));

    // Plot the model and components
    model.plotOn(frame, Name("fit"), LineColor(kBlue));
    model.plotOn(frame, Components(crystalBall), LineStyle(kDashed), LineColor(kRed), Name("crystalball"));
    model.plotOn(frame, Components(gauss), LineStyle(kDashed), LineColor(kGreen), Name("gauss"));

    // Customize frame
    frame->GetXaxis()->SetTitle("log(#chi^{2}_{IP})");
    frame->GetYaxis()->SetTitle("Candidates / bin");
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->Draw();
    frame->SetMinimum(1);
    frame->SetMaximum(maxcounts);
    mainPad->SetLogy();
    mainPad->Modified();
    mainPad->Update();
    // Draw a vertical line at logIPChi2 = 0 to highlight the Gaussian restriction
    // TLine* zeroline = new TLine(0, 0, 0, frame->GetMaximum());
    // zeroline->SetLineColor(kBlack);
    // zeroline->SetLineStyle(kDotted);
    // zeroline->SetLineWidth(2);
    // zeroline->Draw();

    // Add a legend
    TLegend *legend = new TLegend(0.70, 0.65, 0.92, 0.92);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry("data", "Data", "p");
    legend->AddEntry("fit", "Total Fit", "l");
    legend->AddEntry("crystalball", "Crystal Ball (Prompt)", "l");
    legend->AddEntry("gauss", "Gaussian (Nonprompt)", "l");
    legend->Draw();

    // Add text with fit parameters
    TPaveText *fitInfo = new TPaveText(0.17, 0.65, 0.58, 0.92, "NDC");
    fitInfo->SetBorderSize(0);
    fitInfo->SetFillStyle(0);
    fitInfo->SetTextAlign(12);
    fitInfo->SetTextSize(0.035);
    fitInfo->AddText(Form("%s data in %s", particle.Data(), collsys.Data()));
    fitInfo->AddText(Form("CB Mean = %.3f #pm %.3f", mean.getVal(), mean.getError()));
    fitInfo->AddText(Form("CB Sigma = %.3f #pm %.3f", sigma.getVal(), sigma.getError()));
    fitInfo->AddText(Form("CB Alpha = %.2f #pm %.2f", alpha.getVal(), alpha.getError()));
    fitInfo->AddText(Form("CB n = %.2f #pm %.2f", n.getVal(), n.getError()));
    fitInfo->AddText(Form("Gauss Mean = %.3f #pm %.3f", gmean.getVal(), gmean.getError()));
    fitInfo->AddText(Form("Gauss Sigma = %.3f #pm %.3f", gsigma.getVal(), gsigma.getError()));
    fitInfo->AddText(Form("Prompt Fraction = %.2f #pm %.2f", fraction.getVal(), fraction.getError()));
    fitInfo->Draw();

    // Add LHCb label
    TLatex *lhcbLabel = new TLatex();
    lhcbLabel->SetTextSize(0.05);
    lhcbLabel->SetTextFont(132);
    lhcbLabel->SetNDC();
    lhcbLabel->DrawLatex(0.2, 0.94, "LHCb SMOG2");

    // Bottom panel for pulls
    TPad *pullPad = (TPad *)canvas->GetPad(2);
    pullPad->SetPad(0, 0.0, 1, 0.3);
    pullPad->SetTopMargin(0.03);
    pullPad->SetBottomMargin(0.3);
    pullPad->SetLeftMargin(0.15);
    pullPad->SetRightMargin(0.05);
    pullPad->cd();

    // Create pull plot
    RooHist *pullHist = frame->pullHist();
    RooPlot *pullFrame = logIPChi2.frame(Title(""));
    pullFrame->addPlotable(pullHist, "P");
    pullFrame->SetTitle("");

    // Customize pull frame
    pullFrame->GetXaxis()->SetTitle("log(#chi^{2}_{IP})");
    pullFrame->GetXaxis()->SetLabelSize(0.12);
    pullFrame->GetXaxis()->SetTitleSize(0.12);
    pullFrame->GetXaxis()->SetTitleOffset(0.9);
    pullFrame->GetYaxis()->SetTitle("Pull");
    pullFrame->GetYaxis()->SetLabelSize(0.12);
    pullFrame->GetYaxis()->SetTitleSize(0.12);
    pullFrame->GetYaxis()->SetTitleOffset(0.5);
    pullFrame->GetYaxis()->SetNdivisions(505);
    pullFrame->GetYaxis()->SetRangeUser(-5, 5);
    pullFrame->Draw();

    // Draw reference lines at y=0, -3, 3
    TLine *zeroLine = new TLine(pullFrame->GetXaxis()->GetXmin(), 0, pullFrame->GetXaxis()->GetXmax(), 0);
    zeroLine->SetLineStyle(2);
    zeroLine->Draw();

    // Draw a vertical line at logIPChi2 = 0 in the pull plot as well
    TLine *zeroLineX = new TLine(0, -5, 0, 5);
    zeroLineX->SetLineColor(kBlack);
    zeroLineX->SetLineStyle(kDotted);
    zeroLineX->SetLineWidth(2);
    zeroLineX->Draw();

    // Save the plot
    canvas->SaveAs(outputdir + "/" + particleType + "_" + collsys + "_logIPChi2_fit.pdf");
    canvas->SaveAs(outputdir + "/" + particleType + "_" + collsys + "_logIPChi2_fit.png");

    // Calculate the fraction of nonprompt D0 candidates (IPChi2 > 0 and fitted by Gaussian)
    double totalArea = 1.0; // Normalized PDFs
    double nonpromptFraction = (1 - fraction.getVal());
    double nonpromptFractionErr = fraction.getError();

    // Print fraction information
    std::cout << collsys << " nonprompt D0 fraction: " << nonpromptFraction * 100 << "% ± "
              << nonpromptFractionErr * 100 << "%" << std::endl;

    // Store results in a vector
    std::vector<double> results;
    results.push_back(mean.getVal());        // [0] CB mean
    results.push_back(mean.getError());      // [1] CB mean error
    results.push_back(sigma.getVal());       // [2] CB sigma
    results.push_back(sigma.getError());     // [3] CB sigma error
    results.push_back(alpha.getVal());       // [4] CB alpha
    results.push_back(alpha.getError());     // [5] CB alpha error
    results.push_back(n.getVal());           // [6] CB n
    results.push_back(n.getError());         // [7] CB n error
    results.push_back(gmean.getVal());       // [8] Gauss mean
    results.push_back(gmean.getError());     // [9] Gauss mean error
    results.push_back(gsigma.getVal());      // [10] Gauss sigma
    results.push_back(gsigma.getError());    // [11] Gauss sigma error
    results.push_back(fraction.getVal());    // [12] Prompt fraction
    results.push_back(fraction.getError());  // [13] Prompt fraction error
    results.push_back(nonpromptFraction);    // [14] nonprompt fraction
    results.push_back(nonpromptFractionErr); // [15] nonprompt fraction error

    delete canvas;
    delete fitResult;

    return results;
}