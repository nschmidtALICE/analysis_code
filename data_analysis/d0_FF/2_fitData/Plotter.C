#include "Plotter.h"
#include <sys/stat.h>
#include <ctime>
#include <cstdlib>
#include <random>
#include <sstream>

// Constructor implementation
Plotter::Plotter(const std::string& resonanceType, const std::string& basepath, 
                int bin, bool binned, const std::string& range, const std::string& name)
    : resonance(resonanceType), fitBin(bin), isBinned(binned), 
      range(range), nameKey(name), format("png") 
{
    // Set observation selection based on range
    if (range.find("z") != std::string::npos) {
        obsSelection = "zT";
    } else {
        obsSelection = "dR";
    }
    
    // Set basepath with trailing slash
    if (basepath.back() == '/') {
        this->basepath = basepath;
    } else {
        this->basepath = basepath + "/";
    }
    
    std::cout << "Save all histograms to: " << this->basepath << std::endl;
    
    // Create directories if they don't exist
    ensureDirectoryExists(this->basepath);
    ensureDirectoryExists(this->basepath + "MassFits_" + obsSelection + "/");
}

// Helper method to ensure a directory exists
void Plotter::ensureDirectoryExists(const std::string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0 || !(info.st_mode & S_IFDIR)) {
        std::cout << "Creating directory: " << path << std::endl;
        // Create directory with mode 0755 (owner:rwx, group:r-x, other:r-x)
        #ifdef _WIN32
        int result = mkdir(path.c_str());
        #else
        int result = mkdir(path.c_str(), 0755);
        #endif
        
        if (result != 0) {
            std::cerr << "Failed to create directory: " << path << std::endl;
        }
    }
}

// Generate a unique identifier for plot objects
std::string Plotter::getUniqueId() {
    std::time_t timestamp = std::time(nullptr);
    
    // Create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(10000, 99999);
    int rand_val = dis(gen);
    
    std::stringstream ss;
    ss << timestamp << "_" << rand_val;
    return ss.str();
}

// Helper method to set histogram properties
void Plotter::setGraphHisto(TH1* histo, const std::string& xTitle, const std::string& yTitle, 
                           double border, bool logPlot) {
    if (!histo) return;
    
    gStyle->SetOptStat(0);
    TGaxis::SetMaxDigits(2);
    
    // Set axis ranges
    double min = histo->GetMinimum();
    double max = histo->GetMaximum();
    double range = max - min;
    
    double maxNew, minNew;
    if (range > 0) {
        maxNew = min + range * (1.0 + 3 * border);
        minNew = max - range * (1.0 + border);
    } else {
        maxNew = max + (-1) * range * (1.0 + 2 * border);
        minNew = min - (-1) * range * (1.0 + border);
    }
    
    minNew = min;
    
    // Handle log scale minimum
    if (logPlot && minNew < 1) {
        minNew = 1;
    }
    
    histo->GetYaxis()->SetRangeUser(minNew, maxNew);
    
    // Set axis titles
    if (!xTitle.empty()) {
        histo->GetXaxis()->SetTitle(xTitle.c_str());
    }
    
    if (!yTitle.empty()) {
        histo->GetYaxis()->SetTitle(yTitle.c_str());
    }
}

// Implementation of individualMassFitPlot method
TH1* Plotter::individualMassFitPlot(RooRealVar* sigYieldParam, RooAbsPdf* extendedPdf, 
                                   RooRealVar* massVar, RooDataSet* data, 
                                   const std::string& fitTypeName) {
    std::cout << "\n==== Creating mass fit plot for bin " << fitBin 
              << " with model " << fitTypeName << " ====" << std::endl;
    
    // Generate unique identifier
    std::string uniqueId = getUniqueId();
    
    // Set output directory
    std::string outputDir = basepath + "MassFits_" + obsSelection + "/";
    std::cout << "  Output directory: " << outputDir << std::endl;
    
    // Check data
    std::cout << "  Data entries: " << data->numEntries() << std::endl;
    
    // Extract all the parameters from binning Range
    RooAbsBinning& binFull = massVar->getBinning("fullRange");   // Use reference instead of pointer
    RooAbsBinning& binSig = massVar->getBinning("signalRange");  // Use reference instead of pointer
    std::cout << "  Mass range: " << binFull.lowBound() << "-" << binFull.highBound() << std::endl;
    
    massVar->setRange(binFull.lowBound(), binFull.highBound());
    RooPlot* frame = massVar->frame(RooFit::Bins(40));
    frame->SetTitle(("Signal for " + range).c_str());
    
    std::cout << "  Plotting data..." << std::endl;
    data->plotOn(frame, RooFit::Name("datahistogram"), 
                RooFit::LineWidth(1), RooFit::MarkerSize(0.5), 
                RooFit::MarkerStyle(20));
    
    // Create histogram from data
    std::cout << "  Creating histogram from data..." << std::endl;
    std::string histName = "h_MassSignal" + std::to_string(fitBin) + "_" + 
                          fitTypeName + "_" + uniqueId + "_data";
    TH1* h_data = data->createHistogram(histName.c_str(), *massVar);
    
    std::cout << "  Plotting fit models..." << std::endl;
    try {
        // Plot total fit
        extendedPdf->plotOn(frame, RooFit::LineWidth(2), RooFit::Name("TotalFit"));
        
        // Plot background fit component
        extendedPdf->plotOn(frame, RooFit::Components("bkg_pdf_ext"),
                          RooFit::Name("BackgroundFit"),
                          RooFit::LineStyle(kDashed),
                          RooFit::LineWidth(2));
        
        // Plot signal fit component
        extendedPdf->plotOn(frame, RooFit::Components("sig_pdf_ext"),
                          RooFit::Name("SignalFit"),
                          RooFit::LineColor(kRed),
                          RooFit::LineStyle(kDashed),
                          RooFit::LineWidth(2));
    } catch (std::exception& e) {
        std::cerr << "  ERROR plotting fit components: " << e.what() << std::endl;
    }
    
    // Save output file path
    std::string output_file = outputDir + "Bin" + std::to_string(fitBin) + "_" + 
                             fitTypeName + (isBinned ? "Binned" : "Unbinned") + "." + format;
    std::cout << "  Will save plot to: " << output_file << std::endl;
    
    try {
        // Draw and save plot
        std::cout << "  Creating canvas..." << std::endl;
        TCanvas* canvas = new TCanvas(("canvas_" + uniqueId).c_str(), 
                                    ("canvas_" + uniqueId).c_str(), 800*2, 400*2);
        
        // Set up pad, draw frame
        TPad* pad = new TPad(("myPad_" + uniqueId).c_str(), 
                           ("The pad_" + uniqueId).c_str(), 0, 0, 1, 1);
        pad->SetLeftMargin(0.15);
        pad->SetBottomMargin(0.15);
        pad->Draw();
        pad->cd();
        
        frame->Draw();
        
        std::cout << "  Saving canvas to " << output_file << std::endl;
        canvas->SaveAs(output_file.c_str());
        std::cout << "  Plot saved successfully!" << std::endl;
        
        // Clean up to prevent memory leaks
        std::cout << "  Cleaning up..." << std::endl;
        delete canvas;
        // std::cout << "  Canvas deleted!" << std::endl;
        // delete pad;
        // std::cout << "  pad cleaned up!" << std::endl;
    } catch (std::exception& e) {
        std::cerr << "  ERROR saving plot: " << e.what() << std::endl;
    }
    
    return h_data;
}

// Implementation of ipchi2FitPlot method
TH1* Plotter::ipchi2FitPlot(const std::string& resonance, RooRealVar* logIpchi2, 
                           RooDataSet* data, RooAbsPdf* totalPdf, 
                           RooAbsPdf* nonpromptPdf, RooAbsPdf* promptPdf, 
                           RooAbsPdf* backgroundPdf,
                           RooAbsReal* promptYield, RooAbsReal* nonpromptYield) {
    // Get parameters from the fit
    RooArgSet* params = totalPdf->getParameters(*data);
    double promptFrac = -1;
    
    RooRealVar* promptFracVar = dynamic_cast<RooRealVar*>(params->find("prompt_frac"));
    if (promptFracVar) {
        promptFrac = promptFracVar->getVal();
    }
    
    // Parameters for Bukin functions
    // Prompt component
    double xpPrompt = -1, sigmaPrompt = -1, xiPrompt = -1, rho1Prompt = -1, rho2Prompt = -1;
    
    RooRealVar* var = dynamic_cast<RooRealVar*>(params->find("xp_prompt"));
    if (var) xpPrompt = var->getVal();
    
    var = dynamic_cast<RooRealVar*>(params->find("sigma_prompt"));
    if (var) sigmaPrompt = var->getVal();
    
    var = dynamic_cast<RooRealVar*>(params->find("xi_prompt"));
    if (var) xiPrompt = var->getVal();
    
    var = dynamic_cast<RooRealVar*>(params->find("rho1_prompt"));
    if (var) rho1Prompt = var->getVal();
    
    var = dynamic_cast<RooRealVar*>(params->find("rho2_prompt"));
    if (var) rho2Prompt = var->getVal();
    
    // Non-prompt component
    double xpNonprompt = -1, sigmaNonprompt = -1, xiNonprompt = -1, rho1Nonprompt = -1, rho2Nonprompt = -1;
    
    var = dynamic_cast<RooRealVar*>(params->find("xp_nonprompt"));
    if (var) xpNonprompt = var->getVal();
    
    var = dynamic_cast<RooRealVar*>(params->find("sigma_nonprompt"));
    if (var) sigmaNonprompt = var->getVal();
    
    var = dynamic_cast<RooRealVar*>(params->find("xi_nonprompt"));
    if (var) xiNonprompt = var->getVal();
    
    var = dynamic_cast<RooRealVar*>(params->find("rho1_nonprompt"));
    if (var) rho1Nonprompt = var->getVal();
    
    var = dynamic_cast<RooRealVar*>(params->find("rho2_nonprompt"));
    if (var) rho2Nonprompt = var->getVal();
    
    // Get yields
    double promptVal = 0, nonpromptVal = 0, sigYield = 0, bkgYield = 0;
    
    if (promptYield && nonpromptYield) {
        promptVal = promptYield->getVal();
        nonpromptVal = nonpromptYield->getVal();
        sigYield = promptVal + nonpromptVal;
    } else {
        var = dynamic_cast<RooRealVar*>(params->find("sig_yieldLim"));
        if (var) {
            sigYield = var->getVal();
            promptVal = sigYield * promptFrac;
            nonpromptVal = sigYield * (1 - promptFrac);
        }
    }
    
    if (backgroundPdf) {
        var = dynamic_cast<RooRealVar*>(params->find("bkg_yieldLim"));
        if (var) bkgYield = var->getVal();
    }
    
    // Create frame for plotting
    RooPlot* frame1 = logIpchi2->frame(RooFit::Bins(50));
    frame1->SetTitle(("log(IP Chi2) Distribution for " + range).c_str());
    
    // Plot data
    data->plotOn(frame1, RooFit::Name("datahistogram"), 
               RooFit::LineColor(kGray+2), RooFit::LineWidth(1), 
               RooFit::MarkerSize(0.5), RooFit::MarkerStyle(20));
    
    // Plot models
    totalPdf->plotOn(frame1, RooFit::LineColor(kGreen+1), 
                    RooFit::LineStyle(1), RooFit::LineWidth(2));
    
    // Plot prompt component
    totalPdf->plotOn(frame1, RooFit::Components("prompt_pdf"), 
                    RooFit::LineColor(kBlue), 
                    RooFit::LineStyle(1), RooFit::LineWidth(2));
    
    // Plot non-prompt component
    totalPdf->plotOn(frame1, RooFit::Components("nonprompt_pdf"), 
                    RooFit::LineColor(kRed), 
                    RooFit::LineStyle(2), RooFit::LineWidth(2));
    
    // Plot background if available
    if (backgroundPdf) {
        totalPdf->plotOn(frame1, RooFit::Components("bkg_pdf"), 
                        RooFit::LineColor(kBlack), 
                        RooFit::LineStyle(4), RooFit::LineWidth(2));
    }
    
    // Create legend with fit parameters
    TLegend* myLegend1 = new TLegend(0.55, 0.50, 0.75, 0.88);
    myLegend1->SetTextFont(42);
    myLegend1->SetBorderSize(0);
    myLegend1->SetFillStyle(0);
    myLegend1->SetFillColor(0);
    myLegend1->SetMargin(0.25);
    myLegend1->SetTextSize(0.030);
    
    // Add Bukin function parameters to legend
    myLegend1->AddEntry(frame1, "Prompt Bukin:", "");
    myLegend1->AddEntry(frame1, Form(" Peak: %.3f", xpPrompt), "");
    myLegend1->AddEntry(frame1, Form(" Width: %.3f", sigmaPrompt), "");
    myLegend1->AddEntry(frame1, Form(" Asym: %.3f", xiPrompt), "");
    myLegend1->AddEntry(frame1, Form(" Tails: %.3f, %.3f", rho1Prompt, rho2Prompt), "");
    
    myLegend1->AddEntry(frame1, "Non-prompt Bukin:", "");
    myLegend1->AddEntry(frame1, Form(" Peak: %.3f", xpNonprompt), "");
    myLegend1->AddEntry(frame1, Form(" Width: %.3f", sigmaNonprompt), "");
    myLegend1->AddEntry(frame1, Form(" Asym: %.3f", xiNonprompt), "");
    myLegend1->AddEntry(frame1, Form(" Tails: %.3f, %.3f", rho1Nonprompt, rho2Nonprompt), "");
    
    myLegend1->AddEntry(frame1, Form("Prompt frac: %.1f%%", promptFrac*100), "");
    
    // Create yield legend
    TLegend* myLegend3 = new TLegend(0.25, 0.81, 0.45, 0.91);
    myLegend3->SetTextFont(42);
    myLegend3->SetBorderSize(0);
    myLegend3->SetFillStyle(0);
    myLegend3->SetFillColor(0);
    myLegend3->SetMargin(0.25);
    myLegend3->SetTextSize(0.03);
    
    myLegend3->AddEntry(frame1, Form("Prompt yield: %.2f", promptVal), "");
    myLegend3->AddEntry(frame1, Form("Non-prompt yield: %.2f", nonpromptVal), "");
    if (backgroundPdf) {
        myLegend3->AddEntry(frame1, Form("BKG yield: %.2f", bkgYield), "");
    }
    
    // Create pull distribution
    RooHist* dataHist = frame1->getHist("datahistogram");
    RooCurve* curve1 = dynamic_cast<RooCurve*>(frame1->getObject(1));
    RooHist* hpull = nullptr;
    if (dataHist && curve1) {
        hpull = dataHist->makePullHist(*curve1, true);
    }
    
    RooPlot* frame2 = logIpchi2->frame(RooFit::Title("Pull Distribution"), RooFit::Bins(50));
    if (hpull) {
        frame2->addPlotable(hpull, "P");
        frame2->getAttMarker()->SetMarkerSize(0.5);
        frame2->getAttLine()->SetLineWidth(1);
    }
    
    // Draw fit and pull distribution
    TCanvas* canvasFit = new TCanvas("canvasFit", "canvasFit", 800*2, 400*2);
    canvasFit->Divide(2);
    
    // First pad with log scale for full range
    canvasFit->cd(1);
    canvasFit->GetPad(1)->SetLogy();  // Set y-axis to logarithmic scale
    
    // Set appropriate y-axis range for log scale
    double max_val = frame1->GetMaximum();
    frame1->SetAxisRange(0.9, max_val*5, "Y");  // Min value > 0 for log scale
    
    frame1->Draw();
    myLegend1->Draw();
    myLegend3->Draw();
    
    // Second pad for zoomed region to see asymmetry better
    canvasFit->cd(2);
    // Clone frame for zoomed view
    RooPlot* frame_zoom = logIpchi2->frame(RooFit::Title("Zoomed View"), RooFit::Bins(50));
    data->plotOn(frame_zoom, RooFit::Name("datahistogram_zoom"), 
               RooFit::LineColor(kGray+2), RooFit::LineWidth(1), 
               RooFit::MarkerSize(0.5), RooFit::MarkerStyle(20));
    
    totalPdf->plotOn(frame_zoom, RooFit::LineColor(kGreen+1), 
                    RooFit::LineStyle(1), RooFit::LineWidth(2));
    
    totalPdf->plotOn(frame_zoom, RooFit::Components("prompt_pdf"), 
                    RooFit::LineColor(kBlue), 
                    RooFit::LineStyle(1), RooFit::LineWidth(2));
    
    totalPdf->plotOn(frame_zoom, RooFit::Components("nonprompt_pdf"), 
                    RooFit::LineColor(kRed), 
                    RooFit::LineStyle(2), RooFit::LineWidth(2));
    
    // Focus on central region to see asymmetry
    frame_zoom->SetAxisRange(-1.5, 4, "X");
    frame_zoom->Draw();
    
    // Save plot
    std::string output_dir = basepath + "IPChi2Fits_" + obsSelection + "/";
    ensureDirectoryExists(output_dir);
    
    std::string outputFile = output_dir + "IPChi2Fit_SingleBukin_" + std::to_string(fitBin) + "." + format;
    canvasFit->SaveAs(outputFile.c_str());
    
    // Create pull canvas
    TCanvas* canvasPull = new TCanvas("canvasPull", "canvasPull", 800, 400);
    canvasPull->cd();
    gPad->SetLeftMargin(0.15);
    frame2->GetYaxis()->SetTitleOffset(1.6);
    frame2->Draw();
    
    std::string pullOutputFile = output_dir + "IPChi2Fit_SingleBukin_Pull_" + std::to_string(fitBin) + "." + format;
    canvasPull->SaveAs(pullOutputFile.c_str());
    
    // Create histogram for return
    TH1* h_data = data->createHistogram(
        Form("h_IPChi2_SingleBukin_%d_data", fitBin), *logIpchi2
    );
    
    if (h_data) {
        double min = logIpchi2->getMin();
        double max = logIpchi2->getMax();
        h_data->GetXaxis()->SetRangeUser(min, max);
    }
    
    // Clean up
    delete canvasFit;
    delete canvasPull;
    delete myLegend1;
    delete myLegend3;
    
    return h_data;
}

// Implementation of plotMultiple method
void Plotter::plotMultiple(const std::string& name, RooRealVar* lifetimeTag, 
                          bool drawLogDiff, RooDataSet* data1, 
                          RooDataSet* data2, RooDataSet* data3) {
    // Construct plot frames
    RooPlot* frame1 = lifetimeTag->frame(RooFit::Bins(225));
    RooPlot* frame2 = lifetimeTag->frame(RooFit::Bins(225));
    
    // Create histograms from datasets
    TH1* dummyHist1 = nullptr;
    TH1* dummyHist2 = nullptr;
    TH1* dummyHist3 = nullptr;
    TH1* dummyHist1C = nullptr;
    TH1* dummyHist2C = nullptr;
    TH1* dummyHist3C = nullptr;
    
    double nData = 100;
    std::vector<double> maxList;
    
    if (data1) {
        data1->plotOn(frame1, RooFit::DataError(RooAbsData::SumW2), // Fix to RooAbsData::SumW2
                     RooFit::LineWidth(1), RooFit::MarkerSize(0.6), 
                     RooFit::MarkerStyle(20), RooFit::MarkerColor(kRed));
        
        std::string histName = "d1" + std::to_string(fitBin) + name;
        dummyHist1 = data1->createHistogram(histName.c_str(), *lifetimeTag);
        nData = data1->sumEntries();
        dummyHist1C = static_cast<TH1*>(dummyHist1->Clone((histName + "Clone").c_str()));
        
        maxList.push_back(dummyHist1->GetMaximum());
    }
    
    if (data2) {
        data2->plotOn(frame1, RooFit::DataError(RooAbsData::SumW2), // Fix to RooAbsData::SumW2
                     RooFit::LineWidth(1), RooFit::MarkerSize(0.6), 
                     RooFit::MarkerStyle(21), RooFit::MarkerColor(kBlue));
        
        std::string histName = "d2" + std::to_string(fitBin) + name;
        dummyHist2 = data2->createHistogram(histName.c_str(), *lifetimeTag);
        dummyHist2C = static_cast<TH1*>(dummyHist2->Clone((histName + "Clone").c_str()));
        
        maxList.push_back(dummyHist2->GetMaximum());
    }
    
    if (data3) {
        data3->plotOn(frame1, RooFit::DataError(RooAbsData::SumW2), // Fix to RooAbsData::SumW2
                     RooFit::LineWidth(1), RooFit::MarkerSize(0.6), 
                     RooFit::MarkerStyle(22), RooFit::MarkerColor(kGreen));
        
        std::string histName = "d3" + std::to_string(fitBin) + name;
        dummyHist3 = data3->createHistogram(histName.c_str(), *lifetimeTag);
        dummyHist3C = static_cast<TH1*>(dummyHist3->Clone((histName + "Clone").c_str()));
        
        maxList.push_back(dummyHist3->GetMaximum());
    }
    
    // Create legends
    TLegend* myLegend1 = new TLegend(0.3, 0.7, 0.55, 0.88);
    myLegend1->SetTextFont(42);
    myLegend1->SetBorderSize(0);
    myLegend1->SetFillStyle(0);
    myLegend1->SetFillColor(0);
    myLegend1->SetMargin(0.25);
    myLegend1->SetTextSize(0.04);
    
    TLegend* myLegend2 = new TLegend(0.3, 0.7, 0.55, 0.88);
    myLegend2->SetTextFont(42);
    myLegend2->SetBorderSize(0);
    myLegend2->SetFillStyle(0);
    myLegend2->SetFillColor(0);
    myLegend2->SetMargin(0.25);
    myLegend2->SetTextSize(0.04);
    
    // Create canvas and set up first pad
    TCanvas* canvas = new TCanvas(("canvasMulti" + std::to_string(fitBin)).c_str(), 
                                 ("canvasMulti" + std::to_string(fitBin)).c_str(), 
                                 800*2, 400*2);
    canvas->Divide(2);
    canvas->cd(1)->SetLogy();
    
    TPad* myPad = new TPad("myPad", "The pad", 0, 0, 1, 1);
    myPad->SetLeftMargin(0.7);
    myPad->SetTopMargin(0.07);
    myPad->SetRightMargin(0.1);
    myPad->SetBottomMargin(0.15);
    myPad->Draw();
    
    // Find max y value
    double max = 0;
    if (!maxList.empty()) {
        max = *std::max_element(maxList.begin(), maxList.end());
    }
    std::cout << "oo max in frame: " << max << std::endl;
    frame1->SetAxisRange(1, max*1.1, "Y");
    frame1->Draw("E");
    
    // Add histograms to legend
    if (data1 && dummyHist1C) {
        dummyHist1C->SetLineColor(kRed);
        dummyHist1C->SetMarkerColor(kRed);
        dummyHist1C->SetMarkerSize(0.7);
        dummyHist1C->SetMarkerStyle(20);
        myLegend1->AddEntry(dummyHist1C, 
                          Form("%s: %.0f Evts", data1->GetName(), data1->sumEntries()), 
                          "PL");
    }
    
    if (data2 && dummyHist2C) {
        dummyHist2C->SetLineColor(kBlue);
        dummyHist2C->SetMarkerColor(kBlue);
        dummyHist2C->SetMarkerSize(0.7);
        dummyHist2C->SetMarkerStyle(21);
        myLegend1->AddEntry(dummyHist2C, 
                          Form("%s: %.0f Evts", data2->GetName(), data2->sumEntries()), 
                          "PL");
    }
    
    if (data3 && dummyHist3C) {
        dummyHist3C->SetLineColor(kGreen);
        dummyHist3C->SetMarkerColor(kGreen);
        dummyHist3C->SetMarkerSize(0.7);
        dummyHist3C->SetMarkerStyle(22);
        myLegend1->AddEntry(dummyHist3C, 
                          Form("%s: %.0f Evts", data3->GetName(), data3->sumEntries()), 
                          "PL");
    }
    
    myLegend1->Draw();
    
    // Second pad for ratio/difference
    if (drawLogDiff) {
        canvas->cd(2)->SetLogy();
    } else {
        canvas->cd(2);
    }
    
    // Create the difference/ratio plot
    if (data1 && data2 && dummyHist1C && dummyHist2C) {
        TPad* myPad3 = new TPad("myPad3", "The pad3", 0, 0, 1, 1);
        myPad3->SetLeftMargin(0.8);
        myPad3->SetTopMargin(0.08);
        myPad3->SetRightMargin(0.1);
        myPad3->SetBottomMargin(0.15);
        myPad3->Draw();
        
        std::string yLegend;
        if (drawLogDiff) {
            yLegend = std::string(data1->GetName()) + "-" + std::string(data2->GetName());
        } else {
            double nC1 = dummyHist1C->GetEntries();
            double nC2 = dummyHist2C->GetEntries();
            if (nC1 > 0) {
                dummyHist1C->Scale(1.0/nC1);
            }
            if (nC2 > 0) {
                dummyHist2C->Scale(1.0/nC2);
            }
            yLegend = std::string(data1->GetName()) + "_{norm}-" + std::string(data2->GetName()) + "_{norm}";
        }
        
        setGraphHisto(dummyHist1C, "t [s]", yLegend);
        dummyHist1C->Add(dummyHist2C, -1);
        
        if (drawLogDiff) {
            dummyHist1C->GetYaxis()->SetRangeUser(1, dummyHist1C->GetMaximum()*1.1);
        }
        
        if (data2) {
            dummyHist2C->Draw("E");
        }
        if (data3) {
            dummyHist3C->Draw("E");
        }
        
        dummyHist1C->Draw("E");
        
        double yieldDiff = dummyHist1C->Integral();
        myLegend2->AddEntry(dummyHist2C, Form("yield of diff: %.0f", yieldDiff), "PL");
        myLegend2->Draw();
    }
    
    // Save the plot
    std::string outputPath = basepath + "TimeFits_" + obsSelection + "/Mult_Bin" + 
                            name + "_" + std::to_string(fitBin) + ".pdf";
    canvas->SaveAs(outputPath.c_str());
    
    // Clean up
    delete canvas;
    delete myPad;
    delete myLegend1;
    delete myLegend2;
    
    // Delete histograms if they were created
    if (dummyHist1) delete dummyHist1;
    if (dummyHist2) delete dummyHist2;
    if (dummyHist3) delete dummyHist3;
    if (dummyHist1C) delete dummyHist1C;
    if (dummyHist2C) delete dummyHist2C;
    if (dummyHist3C) delete dummyHist3C;
}

// Implementation of extractCorParamRnd method
std::pair<double, double> Plotter::extractCorParamRnd(const std::string& idString, 
                                                    RooRealVar* massVar, 
                                                    std::vector<TLegend*>& myLegendList, 
                                                    RooDataSet* data, 
                                                    RooDataSet* data2, 
                                                    RooDataSet* data3) {
    // Initialize canvases
    TCanvas* canvas_Rnd = new TCanvas("cMulti_Rnd", "cMultiCorr", 800, 800);
    TCanvas* canvasNorm_Rnd = new TCanvas("cMultiNorm_Rnd", "cMultiCorrNorm", 800, 800);
    TCanvas* canvasDiv_Rnd = new TCanvas("cMultiDiv_Rnd", "cMultiCorrDiv", 800, 800);
    
    TGaxis::SetMaxDigits(2);
    gStyle->SetOptStat(0);
    
    // Create histograms from datasets
    TH1* h_data = data->createHistogram("h_MassSignalRnd_data", *massVar);
    h_data->Rebin(100);
    
    std::vector<double> maxList = {h_data->GetMaximum(), 0.0};
    
    TH1* h_data2 = nullptr;
    if (data2) {
        h_data2 = data2->createHistogram("h_MassSignalRnd_data2", *massVar);
        h_data2->Rebin(100);
        maxList[1] = h_data2->GetMaximum();
    }
    
    TH1* h_data3 = nullptr;
    if (data3) {
        h_data3 = data3->createHistogram("h_MassSignalRnd_data3", *massVar);
        h_data3->Rebin(100);
        maxList.push_back(h_data3->GetMaximum());
    }
    
    // First canvas - normalized
    canvasNorm_Rnd->cd();
    TPad* myPad2 = new TPad("myPad2Rnd", "The padRnd", 0, 0, 1, 1);
    myPad2->SetLeftMargin(0.8);
    myPad2->SetTopMargin(0.08);
    myPad2->SetRightMargin(0.1);
    myPad2->SetBottomMargin(0.15);
    myPad2->Draw();
    
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.5);
    h_data->SetMarkerColor(1);
    h_data->SetLineColor(kGray+2);
    h_data->DrawNormalized("E");
    
    if (data2 && h_data2) {
        h_data2->SetMarkerStyle(20);
        h_data2->SetMarkerSize(0.5);
        h_data2->SetMarkerColor(kRed+2);
        h_data2->SetLineColor(kRed+2);
        h_data2->DrawNormalized("same E");
    }
    
    if (data3 && h_data3) {
        h_data3->SetMarkerStyle(20);
        h_data3->SetMarkerSize(0.5);
        h_data3->SetMarkerColor(kBlue+2);
        h_data3->SetLineColor(kBlue+2);
        h_data3->DrawNormalized("same E");
    }
    
    // Second canvas - raw
    canvas_Rnd->cd();
    TPad* myPad = new TPad("myPadRnd", "The padRnd", 0, 0, 1, 1);
    myPad->SetLeftMargin(0.8);
    myPad->SetTopMargin(0.08);
    myPad->SetRightMargin(0.1);
    myPad->SetBottomMargin(0.15);
    myPad->Draw();
    
    double maxVal = 0;
    for (double val : maxList) {
        if (val > maxVal) maxVal = val;
    }
    
    h_data->GetYaxis()->SetRangeUser(0, maxVal*1.3);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.5);
    h_data->SetMarkerColor(1);
    h_data->SetLineColor(kGray+2);
    h_data->DrawCopy("E");
    
    if (data2 && h_data2) {
        h_data2->SetMarkerStyle(20);
        h_data2->SetMarkerSize(0.5);
        h_data2->SetMarkerColor(kRed+2);
        h_data2->SetLineColor(kRed+2);
        h_data2->DrawCopy("same E");
    }
    
    if (data3 && h_data3) {
        h_data3->SetMarkerStyle(20);
        h_data3->SetMarkerSize(0.5);
        h_data3->SetMarkerColor(kBlue+2);
        h_data3->SetLineColor(kBlue+2);
        h_data3->DrawCopy("same E");
    }
    
    h_data->GetYaxis()->UnZoom();
    
    // Third canvas - ratio
    TGaxis::SetMaxDigits(4);
    canvasDiv_Rnd->cd();
    
    double fitValue = 0.0;
    double fitValueErr = 0.0;
    
    if (data2 && h_data2) {
        TPad* myPad3 = new TPad("myPad3Rnd", "The pad3Rnd", 0, 0, 1, 1);
        myPad3->SetLeftMargin(0.8);
        myPad3->SetTopMargin(0.08);
        myPad3->SetRightMargin(0.1);
        myPad3->SetBottomMargin(0.15);
        myPad3->Draw();
        
        h_data2->Divide(h_data);
        h_data2->SetMarkerStyle(20);
        h_data2->SetMarkerSize(0.5);
        h_data2->SetMarkerColor(1);
        h_data2->SetLineColor(kGray+2);
        h_data2->DrawCopy("E");
        
        // Fit a horizontal line and save the value for the efficiency correction
        double minL = h_data->GetBinCenter(0);
        double maxL = h_data->GetBinCenter(h_data->GetNbinsX());
        TF1* line = new TF1("Line", "[0]", minL, maxL);
        h_data2->Fit("Line", "NQ", "", minL, maxL);
        fitValue = line->GetParameter(0);
        fitValueErr = line->GetParError(0);
        
        line->DrawCopy("same");
        
        if (!myLegendList.empty() && myLegendList[0]) {
            myLegendList[0]->AddEntry(h_data2, 
                                     Form(" Correction Constant %2.3f#pm%2.3f", fitValue, fitValueErr), 
                                     "");
            myLegendList[0]->Draw();
        }
        
        delete line;
    }
    
    // Save canvases
    canvas_Rnd->SaveAs((basepath + "CorrFac" + idString + "_zT.png").c_str());
    canvasDiv_Rnd->SaveAs((basepath + "CorrFac" + idString + "_zTRatio.png").c_str());
    canvasNorm_Rnd->SaveAs((basepath + "CorrFac" + idString + "_zTNorm.png").c_str());
    
    // Clean up
    delete canvas_Rnd;
    delete canvasNorm_Rnd;
    delete canvasDiv_Rnd;
    delete myPad;
    delete myPad2;
    delete h_data;
    if (h_data2) delete h_data2;
    if (h_data3) delete h_data3;
    
    return std::make_pair(fitValue, fitValueErr);
}

// Implementation of teffMap2D method
TH1* Plotter::teffMap2D(TEfficiency* teff, const std::string& effType, const std::string& dataType, 
                       bool closureTF1, bool closureHisto) {
    TCanvas* canvas = new TCanvas("canvas", "canvas", 0, 0, 600, 400);
    
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetNumberContours(100);
    gPad->SetLogx();
    
    // Set title based on efficiency type
    if (effType == "trigger") {
        teff->SetTitle("Trigger Efficiency; sqrt(pt(muon1)*pt(muon2)); eta(tag); eff");
    } else if (effType == "probnnmu") {
        if (!closureTF1 && !closureHisto) {
            teff->SetTitle("ProbNNmu Efficiency; pt(muon); eta(muon); eff");
        } else {
            teff->SetTitle("ProbNNmu Efficiency; pt(tag); eta(tag); eff");
        }
    } else if (effType == "reco") {
        if (closureHisto) {
            teff->SetTitle("Muon reco * IsMuon Efficiency; pt(tag); eta(tag); eff");
        }
    }
    
    teff->Draw("colz");
    canvas->Update();
    
    TH1* effHist = teff->GetPaintedHistogram();
    if (effHist) {
        effHist->SetMinimum(0);
        effHist->SetMaximum(1);
        effHist->Draw("colz");
        canvas->Update();
    }
    
    canvas->SaveAs((basepath + "/eff2DMap_" + effType + "Swap.pdf").c_str());
    
    delete canvas;
    return effHist;
}

// Implementation of splotVals method
std::tuple<double, double, double> Plotter::splotVals(const std::string& resonance, 
                                                    RooAbsPdf* extendedPdf, 
                                                    RooRealVar* massVar, 
                                                    RooDataSet* data, 
                                                    TTree* ttree, 
                                                    RooRealVar* ns, 
                                                    RooRealVar* nb, 
                                                    const std::string& fidCutString, 
                                                    TFile* sFile) {
    RooStats::SPlot* splot = new RooStats::SPlot("splot", "splot", *data, extendedPdf, RooArgList(*ns, *nb));
    
    data->Print("v");
    splot->Print();
    
    // Fix 4: Correct TFile::Open usage - we need to extract the filename from sFile
    std::string outfileName = std::string(sFile->GetName());
    TFile* fvar = TFile::Open(outfileName.c_str(), "recreate");
    
    // Create a new tree with the sweights
    TTree* varTree = nullptr;
    if (fidCutString.empty()) {
        std::stringstream cutStr;
        cutStr << "mass_tag_measured>" << massVar->getMin() 
               << "&&mass_tag_measured<" << massVar->getMax();
        varTree = ttree->CopyTree(cutStr.str().c_str());
    } else {
        std::stringstream cutStr;
        cutStr << "mass_tag_measured>" << massVar->getMin() 
               << "&&mass_tag_measured<" << massVar->getMax() 
               << " && " << fidCutString;
        varTree = ttree->CopyTree(cutStr.str().c_str());
    }
    
    varTree->SetName("ntuple");
    
    // Add branches for sweights
    Double_t sig_sw = 0.0;
    Double_t bkg_sw = 0.0;
    TBranch* br1 = varTree->Branch("sig_sw", &sig_sw, "sig_sw/D");
    TBranch* br2 = varTree->Branch("bkg_sw", &bkg_sw, "bkg_sw/D");
    
    double scale1 = 0.0;
    double scale2 = 0.0;
    
    // Fill the branches with sWeights
    for (Int_t i = 0; i < data->numEntries(); ++i) {
        const RooArgSet* row = data->get(i);
        
        RooAbsReal* sig_sw_var = dynamic_cast<RooAbsReal*>(row->find("sig_yield_sw"));
        RooAbsReal* bkg_sw_var = dynamic_cast<RooAbsReal*>(row->find("bkg_yield_sw"));
        
        if (sig_sw_var && bkg_sw_var) {
            sig_sw = sig_sw_var->getVal();
            bkg_sw = bkg_sw_var->getVal();
            
            br1->Fill();
            br2->Fill();
            
            scale1 += sig_sw;
            scale2 += sig_sw * sig_sw;
        }
    }
    
    std::cout << sFile << std::endl;
    varTree->Write("", TObject::kOverwrite);
    fvar->Close();
    
    delete splot;
    
    double ratio = 0.0;
    if (scale2 != 0.0) {
        ratio = scale1 / scale2;
    }
    
    return std::make_tuple(scale1, scale2, ratio);
}