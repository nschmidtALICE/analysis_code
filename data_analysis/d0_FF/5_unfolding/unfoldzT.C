// UnfoldSpectra.cpp
#include "unfoldzT.h"

// Constructor implementation
UnfoldSpectraClass::UnfoldSpectraClass(const std::string& promptFlag, 
                                     const std::string& inFileNameRM,
                                     const std::string& inFileNameData, 
                                     const std::string& resonance,
                                     const std::vector<int>& ptRangeArray) : 
    useTagPt(false),
    inPathRM("/media/niviths/SSD2/lhcb_analysis_SSD/"),
    dictKey(resonance),
    errorType(RooUnfold::kCovToy) {
    

    figureTag = "D0_" + promptFlag;

    // Set output path
    outpathBase = "/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/OutputMay_TagpT_" + resonance;
    
    // Set prompt/non-prompt flag
    if (promptFlag == "P") {
        isPrompt = true;
        outPath = outpathBase + "/Prompt_0/";
    } else if (promptFlag == "NP") {
        isPrompt = false;
        outPath = outpathBase + "/NonPrompt_0/";
    }
    
    // Set file paths
    applyRMCut = true;
    inFileNRM = inPathRM + inFileNameRM + ".root";
    inFileNData = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA/RawSignalYields_D0/" + inFileNameData + ".root";

    // Set binning for input and output spectrum
    zBinsArrayTruth = {0.2, 0.5, 0.65, 0.75, 0.85, 0.95, 1};
    
    // Set pT bins
    pTBinsArrayTruth = ptRangeArray;
    pTBinsArrayDet = ptRangeArray;
    
    // Initialize arrays for unfolding results
    for (size_t i = 0; i < ptRangeArray.size() - 1; i++) {
        unfoldedArrPerBin.push_back(std::vector<TH1D*>());
    }
    
    // Open input files
    fResponse = new TFile(inFileNRM.c_str());
    fData = new TFile(inFileNData.c_str());
    
    std::cout << "o Open input File with Response Matrix: " << fResponse->GetName() << std::endl;
    std::cout << "o Open input File with Data: " << fData->GetName() << std::endl;
    
    // Get raw spectra
    measuredSpectraArray = getRawSpectra();
    measuredSpectra2D = nullptr;
    measuredSpectra2D = getRawSpectra2D();
    
    // Create output directory
    if (!std::filesystem::exists(outPath)) {
        std::filesystem::create_directories(outPath);
        std::cout << "make new directory: " << outPath << std::endl;
    }
}

void UnfoldSpectraClass::RefoldingTest(int ptBin, int nIter, TH1D* histo, 
                                      RooUnfoldResponse* RM1, RooUnfoldResponse* RM2) {
    // Unfold the measured spectrum with half the RM statistic and refold it
    // back with the other half of the RM statistic
    RooUnfoldBayes unfoldBayes1(RM1, histo, nIter);
    TH1D* unfoldedSpectrum = dynamic_cast<TH1D*>(unfoldBayes1.Hreco(
        static_cast<RooUnfold::ErrorTreatment>(errorType)));
    
    if (!unfoldedSpectrum) {
        std::cerr << "Error: Unfolding failed in RefoldingTest!" << std::endl;
        return;
    }
    
    TH1D* refoldedSpectrum = dynamic_cast<TH1D*>(RM2->ApplyToTruth(unfoldedSpectrum));
    
    if (!refoldedSpectrum) {
        std::cerr << "Error: Refolding failed in RefoldingTest!" << std::endl;
        delete unfoldedSpectrum;
        return;
    }
    
    // Draw the result
    TCanvas* c = new TCanvas("c", "c: pT", 800, 850);
    c->cd();
    c->SetLeftMargin(0.15);
    
    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.05);
    pad1->SetTopMargin(0.05);
    pad1->Draw();
    pad1->cd();
    
    histo->GetYaxis()->SetTitle("counts per bin");
    histo->GetYaxis()->SetTitleSize(0.06);
    histo->GetYaxis()->SetLabelFont(43);
    histo->GetYaxis()->SetLabelSize(20);
    histo->SetNdivisions(505);
    histo->SetLineColor(4);
    histo->SetLineWidth(3);
    histo->SetLineStyle(1);
    histo->Draw("hist E");
    
    refoldedSpectrum->SetLineColor(1);
    refoldedSpectrum->SetLineWidth(2);
    refoldedSpectrum->SetLineStyle(1);
    refoldedSpectrum->DrawCopy("hist E same");
    
    TLegend* leg2 = new TLegend(0.2, 0.8, 0.5, 0.93, "");
    leg2->SetFillColor(10);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.04);
    leg2->AddEntry(refoldedSpectrum, "measured spectrum unfolded (RM1)+refolded (RM2)", "l");
    leg2->AddEntry(histo, "measured spectrum, detLvlv", "l");
    leg2->Draw("same");
    
    c->cd();
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.05);
    pad2->Draw();
    pad2->cd();
    
    TH1D* refoldedSpectrumC = static_cast<TH1D*>(refoldedSpectrum->Clone(
        (std::string(refoldedSpectrum->GetName()) + "Clone").c_str()));
    refoldedSpectrumC->Divide(histo);
    refoldedSpectrumC->SetMarkerStyle(21);
    refoldedSpectrumC->SetMarkerColor(2);
    
    refoldedSpectrumC->GetXaxis()->SetTitleSize(30);
    refoldedSpectrumC->GetXaxis()->SetTitleFont(43);
    refoldedSpectrumC->GetXaxis()->SetTitleOffset(3.0); // was 4
    refoldedSpectrumC->GetXaxis()->SetLabelFont(43);
    refoldedSpectrumC->GetXaxis()->SetLabelSize(20);
    refoldedSpectrumC->GetYaxis()->SetRangeUser(0, 2);
    refoldedSpectrumC->GetYaxis()->SetTitle(("Ratio to orig, nIter=" + std::to_string(nIter)).c_str());
    refoldedSpectrumC->GetYaxis()->SetTitleSize(20);
    refoldedSpectrumC->GetYaxis()->SetTitleFont(43);
    refoldedSpectrumC->GetYaxis()->SetTitleOffset(2.2);
    refoldedSpectrumC->GetYaxis()->SetLabelFont(43);
    refoldedSpectrumC->GetYaxis()->SetLabelSize(20);
    refoldedSpectrumC->GetYaxis()->SetNdivisions(505);
    refoldedSpectrumC->GetXaxis()->SetTitle("z_{T}");
    refoldedSpectrumC->DrawCopy("hist E");
    
    std::vector<TLine*> lines = drawLines();
    for (auto line : lines) {
        line->Draw("same");
    }
    
    c->cd();
    
    std::string fileName = outPath + "RefoldingTestBin" + std::to_string(ptBin) + 
                           "_nIter" + std::to_string(nIter) + "_" + figureTag + ".png";
    c->SaveAs(fileName.c_str());
    
    // Clean up
    delete c;  // This also deletes the pads
    delete refoldedSpectrumC;
    delete unfoldedSpectrum;
    delete refoldedSpectrum;
}

// Implementation of getRawSpectra method
std::vector<TH1D*> UnfoldSpectraClass::getRawSpectra(bool mergeTo2D) {
    std::cout << "\n==== getRawSpectra DEBUG START ====" << std::endl;
    std::cout << "Input parameter: mergeTo2D = " << (mergeTo2D ? "true" : "false") << std::endl;
    
    std::vector<TH1D*> histoArray;
    
    // Check if we already have the spectra array
    if (!measuredSpectraArray.empty()) {
        std::cout << "Using existing measuredSpectraArray with " << measuredSpectraArray.size() 
                  << " histograms" << std::endl;
        histoArray = measuredSpectraArray;
    } else {
        std::cout << "measuredSpectraArray is empty, loading spectra from file" << std::endl;
        
        // Fill zBinsArrayDet if it's empty
        if (zBinsArrayDet.empty() && !mergeTo2D) {
            std::cout << "Warning: zBinsArrayDet is empty!" << std::endl;
        }
        
        std::cout << "Will process " << (pTBinsArrayDet.size() - 1) << " pT bins" << std::endl;
        
        // Open the data file and get the spectra for each pT bin
        for (size_t i = 0; i < pTBinsArrayDet.size() - 1; i++) {
            std::string histname;
            if (isPrompt) {
                // histname = "PromptCorr_" + std::to_string(pTBinsArrayDet[i]) + "_" + 
                histname = "PromptRaw_" + std::to_string(pTBinsArrayDet[i]) + "_" + 
                           std::to_string(pTBinsArrayDet[i+1]);
            } else {
                // histname = "NonPromptCorr_" + std::to_string(pTBinsArrayDet[i]) + "_" + 
                histname = "NonPromptRaw_" + std::to_string(pTBinsArrayDet[i]) + "_" + 
                           std::to_string(pTBinsArrayDet[i+1]);
            }
            
            std::cout << "Processing pT bin " << i << ": " << pTBinsArrayDet[i] << "-" 
                      << pTBinsArrayDet[i+1] << " GeV/c, histogram name: " << histname << std::endl;
            
            // Get the object from file
            TObject* histObj = fData->Get(histname.c_str());
            if (!histObj) {
                std::cerr << "Error: Could not find object " << histname << std::endl;
                continue;
            }
            
            std::cout << "Found object of class: " << histObj->ClassName() << std::endl;
            
            // Check if it's a histogram or a graph
            if (!histObj->InheritsFrom("TH1")) {
                std::cout << "Object is not a histogram, attempting to convert from graph" << std::endl;
                
                // Plot the graph first
                plotHist(histObj, "HmeasuredGraph" + std::to_string(i), "E");
                
                // It's a TGraphErrors, convert to histogram
                TGraphErrors* graph = dynamic_cast<TGraphErrors*>(histObj);
                if (!graph) {
                    std::cerr << "Error: Object is not a TGraphErrors" << std::endl;
                    continue;
                }
                
                std::cout << "Converting TGraphErrors with " << graph->GetN() << " points to histogram" << std::endl;
                
                // Undo scaling by bin width
                std::cout << "Undoing bin width scaling..." << std::endl;
                deScaleGraph(graph);
                
                // Get arrays of x positions and errors
                int nPoints = graph->GetN();
                double* xArray = graph->GetX();
                double* xErrorArray = graph->GetEX();
                double* yArray = graph->GetY();
                double* yErrorArray = graph->GetEY();
                
                std::cout << "Graph points summary:" << std::endl;
                for (int j = 0; j < std::min(nPoints, 5); j++) {
                    std::cout << "  Point " << j << ": x = " << xArray[j] << " ± " << xErrorArray[j]
                              << ", y = " << yArray[j] << " ± " << yErrorArray[j] << std::endl;
                }
                if (nPoints > 5) {
                    std::cout << "  ... (showing only first 5 points)" << std::endl;
                }
                
                // Create array of bin edges
                std::vector<double> binEdges;
                for (int j = 0; j < nPoints; j++) {
                    binEdges.push_back(xArray[j] - xErrorArray[j]);
                }
                // Add final edge
                binEdges.push_back(xArray[nPoints-1] + xErrorArray[nPoints-1]);
                
                std::cout << "Created " << binEdges.size() << " bin edges: ";
                for (size_t j = 0; j < std::min(binEdges.size(), size_t(10)); j++) {
                    std::cout << binEdges[j] << " ";
                }
                if (binEdges.size() > 10) std::cout << "...";
                std::cout << std::endl;
                
                // Create a histogram with the bin edges
                TH1D* histogramTransform = new TH1D((histname + "_hist").c_str(), histname.c_str(), 
                                                  nPoints, binEdges.data());
                
                // Store bin edges for later use
                if (zBinsArrayDet.empty()) {
                    std::cout << "Storing bin edges in zBinsArrayDet" << std::endl;
                    zBinsArrayDet = binEdges;
                }
                
                // Fill the histogram with data from the graph
                for (int j = 0; j < nPoints; j++) {
                    histogramTransform->SetBinContent(j+1, yArray[j]);
                    histogramTransform->SetBinError(j+1, yErrorArray[j]);
                }
                
                std::cout << "Created histogram: " << histogramTransform->GetNbinsX() << " bins, "
                          << "integral = " << histogramTransform->Integral() << ", "
                          << "max = " << histogramTransform->GetMaximum() << std::endl;
                
                // Add to the array
                plotHist(histogramTransform, "HmeasuredHisto" + std::to_string(i), "E");
                histoArray.push_back(histogramTransform);
            } else {
                std::cout << "Object is a histogram, using directly" << std::endl;
                
                // It's already a histogram
                TH1D* histogram = dynamic_cast<TH1D*>(histObj);
                if (!histogram) {
                    std::cerr << "Error: Could not cast to TH1D" << std::endl;
                    continue;
                }
                
                std::cout << "Histogram details: " << histogram->GetNbinsX() << " bins, "
                          << "integral = " << histogram->Integral() << ", "
                          << "max = " << histogram->GetMaximum() << std::endl;
                
                plotHist(histogram, "HmeasuredHisto" + std::to_string(i), "E");
                histogram->Sumw2();
                histoArray.push_back(histogram);
                
                // We need to determine zBinsArrayDet if it's empty
                if (zBinsArrayDet.empty()) {
                    std::cout << "Warning: Need to extract zBinsArrayDet from histogram" << std::endl;
                    int nBins = histogram->GetNbinsX();
                    for (int j = 1; j <= nBins + 1; j++) {
                        zBinsArrayDet.push_back(histogram->GetBinLowEdge(j));
                    }
                    
                    std::cout << "Extracted " << zBinsArrayDet.size() << " bin edges from histogram" << std::endl;
                }
            }
        }
        
        std::cout << "Loaded " << histoArray.size() << " spectra histograms" << std::endl;
        
        // Store the array for later use
        measuredSpectraArray = histoArray;
    }
    
    // If we're not merging to 2D, return the array of 1D histograms
    if (!mergeTo2D) {
        std::cout << "Returning 1D histograms array with " << histoArray.size() << " histograms" << std::endl;
        std::cout << "==== getRawSpectra DEBUG END ====\n" << std::endl;
        return histoArray;
    }
    
    // Otherwise, build and return the 2D histogram
    std::cout << "mergeTo2D is true, returning histogram array for 2D construction" << std::endl;
    std::cout << "==== getRawSpectra DEBUG END ====\n" << std::endl;
    return histoArray;
}

// Implementation of getRawSpectra2D method
TH2D* UnfoldSpectraClass::getRawSpectra2D() {
    std::cout << "\n==== getRawSpectra2D DEBUG START ====" << std::endl;
    
    // If we already have the 2D spectra, return it
    if (measuredSpectra2D) {
        std::cout << "Using existing measuredSpectra2D: " << measuredSpectra2D->GetNbinsX() << "x" 
                  << measuredSpectra2D->GetNbinsY() << " bins, integral = " 
                  << measuredSpectra2D->Integral() << std::endl;
        std::cout << "==== getRawSpectra2D DEBUG END ====\n" << std::endl;
        return measuredSpectra2D;
    }
    
    std::cout << "No existing measuredSpectra2D, creating new one" << std::endl;
    
    // Make sure we have the 1D spectra
    if (measuredSpectraArray.empty()) {
        std::cout << "measuredSpectraArray is empty, calling getRawSpectra(false)" << std::endl;
        measuredSpectraArray = getRawSpectra(false);
    }
    
    std::cout << "Working with " << measuredSpectraArray.size() << " 1D histograms" << std::endl;
    
    // Create temporary vectors of doubles for the pT bin edges
    std::vector<double> pTBinsDouble;
    for (int val : pTBinsArrayDet) {
        pTBinsDouble.push_back(static_cast<double>(val));
    }
    
    std::cout << "Created pTBinsDouble with " << pTBinsDouble.size() << " edges: ";
    for (size_t i = 0; i < std::min(pTBinsDouble.size(), size_t(10)); i++) {
        std::cout << pTBinsDouble[i] << " ";
    }
    if (pTBinsDouble.size() > 10) std::cout << "...";
    std::cout << std::endl;
    
    std::cout << "zBinsArrayDet has " << zBinsArrayDet.size() << " edges: ";
    for (size_t i = 0; i < std::min(zBinsArrayDet.size(), size_t(10)); i++) {
        std::cout << zBinsArrayDet[i] << " ";
    }
    if (zBinsArrayDet.size() > 10) std::cout << "...";
    std::cout << std::endl;
    
    // Create a 2D histogram from all the 1D histograms
    // x-Axis = jet pT bins, y-Axis = zT bins
    std::cout << "Creating 2D histogram with " << (pTBinsDouble.size() - 1) << "x" 
              << (zBinsArrayDet.size() - 1) << " bins" << std::endl;
    
    TH2D* histo2DMeasurement = new TH2D("mesured2D", "mesured2D", 
                                      pTBinsDouble.size() - 1, pTBinsDouble.data(), 
                                      zBinsArrayDet.size() - 1, zBinsArrayDet.data());
    
    // Fill the 2D histogram with data from the 1D histograms
    std::cout << "Filling 2D histogram from 1D histograms..." << std::endl;
    double totalIntegral = 0.0;
    
    for (size_t binx = 1; binx <= pTBinsArrayDet.size() - 1; binx++) {
        if (binx > measuredSpectraArray.size()) {
            std::cerr << "Error: bin index " << binx << " exceeds array size " 
                      << measuredSpectraArray.size() << std::endl;
            continue;
        }
        
        TH1D* hist1D = measuredSpectraArray[binx-1];
        std::cout << "  Filling from 1D histogram for pT bin " << binx << ": " 
                  << "integral = " << hist1D->Integral() << ", nbins = " << hist1D->GetNbinsX() << std::endl;
        
        double binIntegral = 0.0;
        for (int biny = 1; biny <= hist1D->GetNbinsX(); biny++) {
            double content = hist1D->GetBinContent(biny);
            double error = hist1D->GetBinError(biny);
            histo2DMeasurement->SetBinContent(binx, biny, content);
            histo2DMeasurement->SetBinError(binx, biny, error);
            binIntegral += content;
        }
        totalIntegral += binIntegral;
        std::cout << "    bin integral = " << binIntegral << std::endl;
    }
    
    // Set Sumw2 for correct error propagation
    histo2DMeasurement->Sumw2();
    
    std::cout << "Created 2D histogram with integral = " << histo2DMeasurement->Integral() 
              << " (sum of bin contents = " << totalIntegral << ")" << std::endl;
              
    // Check for empty/negative bins
    int emptyBins = 0, negativeBins = 0;
    for (int i = 1; i <= histo2DMeasurement->GetNbinsX(); i++) {
        for (int j = 1; j <= histo2DMeasurement->GetNbinsY(); j++) {
            double content = histo2DMeasurement->GetBinContent(i, j);
            if (content == 0) emptyBins++;
            if (content < 0) negativeBins++;
        }
    }
    std::cout << "2D histogram has " << emptyBins << " empty bins and " 
              << negativeBins << " negative bins" << std::endl;
    
    // Store the 2D histogram for later use
    measuredSpectra2D = histo2DMeasurement;
    
    std::cout << "==== getRawSpectra2D DEBUG END ====\n" << std::endl;
    return histo2DMeasurement;
}
void UnfoldSpectraClass::RefoldingTest2D(int nIter, RooUnfoldResponse* RM1, RooUnfoldResponse* RM2) {
    // Unfold the measured spectrum with half the RM statistic and refold it
    // back with the other half of the RM statistic.
    RooUnfoldBayes unfoldBayes1(RM1, measuredSpectra2D, nIter);
    TH2D* unfoldedSpectrum = dynamic_cast<TH2D*>(unfoldBayes1.Hreco(
        static_cast<RooUnfold::ErrorTreatment>(errorType)));
    
    if (!unfoldedSpectrum) {
        std::cerr << "Error: Unfolding failed in RefoldingTest2D!" << std::endl;
        return;
    }
    
    TH2D* refoldedSpectrum = dynamic_cast<TH2D*>(RM2->ApplyToTruth(unfoldedSpectrum));
    
    if (!refoldedSpectrum) {
        std::cerr << "Error: Refolding failed in RefoldingTest2D!" << std::endl;
        delete unfoldedSpectrum;
        return;
    }
    
    // Draw the result - create a canvas with a pad for each pT bin
    int npTBins = pTBinsArrayTruth.size() - 1;
    TCanvas* c = new TCanvas("c", "c: pT", 800 * npTBins, 850);
    c->Divide(npTBins, 2, 0, 0);
    
    // Create legend
    TLegend* leg2 = new TLegend(0.2, 0.8, 0.5, 0.93, "");
    leg2->SetFillColor(10);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.04);
    
    for (int ptBin = 1; ptBin <= npTBins; ptBin++) {
        // Upper pads for histograms
        c->cd(ptBin);
        TVirtualPad* pad = c->cd(ptBin);
        pad->SetLeftMargin(0.15);
        pad->SetRightMargin(0.05);
        pad->SetTopMargin(0.05);
        pad->SetBottomMargin(0);
        
        // Project histograms to 1D for this pT bin
        TH1D* origPtBinProj = measuredSpectra2D->ProjectionY(
            ("_pyOrig_" + std::to_string(ptBin)).c_str(), ptBin, ptBin);
        TH1D* refoldedPtBinProj = refoldedSpectrum->ProjectionY(
            ("_pyRefolded_" + std::to_string(ptBin)).c_str(), ptBin, ptBin);
        
        origPtBinProj->GetYaxis()->SetTitle("counts per bin");
        origPtBinProj->GetYaxis()->SetTitleSize(0.06);
        origPtBinProj->GetYaxis()->SetLabelFont(43);
        origPtBinProj->GetYaxis()->SetLabelSize(20);
        origPtBinProj->SetNdivisions(505);
        origPtBinProj->SetLineColor(4);
        origPtBinProj->SetLineWidth(3);
        origPtBinProj->SetLineStyle(1);
        origPtBinProj->Draw("hist E");
        
        refoldedPtBinProj->SetLineColor(1);
        refoldedPtBinProj->SetLineWidth(2);
        refoldedPtBinProj->SetLineStyle(1);
        refoldedPtBinProj->Draw("hist E same");
        
        if (ptBin == 1) {
            leg2->AddEntry(refoldedPtBinProj, "measured spectrum unfolded (RM1)+refolded (RM2)", "l");
            leg2->AddEntry(origPtBinProj, "measured spectrum, detLvlv", "l");
            leg2->Draw("same");
        }
        
        // Lower pads for ratios
        c->cd(npTBins + ptBin);
        pad = c->cd(npTBins + ptBin);
        pad->SetLeftMargin(0.15);
        pad->SetRightMargin(0.05);
        pad->SetTopMargin(0);
        pad->SetBottomMargin(0.35);
        
        // Create ratio histogram
        TH1D* refoldedPtBinProjClone = static_cast<TH1D*>(refoldedPtBinProj->Clone(
            (std::string(refoldedPtBinProj->GetName()) + "Clone").c_str()));
        refoldedPtBinProjClone->Divide(origPtBinProj);
        refoldedPtBinProjClone->SetMarkerStyle(21);
        refoldedPtBinProjClone->SetMarkerColor(2);
        
        refoldedPtBinProjClone->GetXaxis()->SetTitleSize(30);
        refoldedPtBinProjClone->GetXaxis()->SetTitleFont(43);
        refoldedPtBinProjClone->GetXaxis()->SetTitleOffset(3.0);
        refoldedPtBinProjClone->GetXaxis()->SetLabelFont(43);
        refoldedPtBinProjClone->GetXaxis()->SetLabelSize(20);
        refoldedPtBinProjClone->GetYaxis()->SetRangeUser(0, 2);
        refoldedPtBinProjClone->GetYaxis()->SetTitle(("Ratio to orig, nIter=" + std::to_string(nIter)).c_str());
        refoldedPtBinProjClone->GetYaxis()->SetTitleSize(20);
        refoldedPtBinProjClone->GetYaxis()->SetTitleFont(43);
        refoldedPtBinProjClone->GetYaxis()->SetTitleOffset(2.2);
        refoldedPtBinProjClone->GetYaxis()->SetLabelFont(43);
        refoldedPtBinProjClone->GetYaxis()->SetLabelSize(20);
        refoldedPtBinProjClone->GetYaxis()->SetNdivisions(505);
        refoldedPtBinProjClone->GetXaxis()->SetTitle("z_{T}");
        refoldedPtBinProjClone->Draw("hist E");
        
        // Draw reference lines
        std::vector<TLine*> lines = drawLines();
        for (TLine* line : lines) {
            TLine* cLine = static_cast<TLine*>(line->Clone("c"));
            cLine->Draw("same");
            delete cLine;
        }
        
        // Clean up projections
        delete refoldedPtBinProjClone;
    }
    
    // Save the canvas
    std::string fileName = outPath + "RefoldingTest2D_nIter" + std::to_string(nIter) + figureTag + ".png";
    c->SaveAs(fileName.c_str());
    
    // Clean up
    delete c;
    delete unfoldedSpectrum;
    delete refoldedSpectrum;
}

void UnfoldSpectraClass::UnfoldingTest2D(RooUnfoldResponse* response, int nIter) {
    // Get the response matrix histograms - use dynamic_cast for proper type conversion
    TH2D* TH2DRM = dynamic_cast<TH2D*>(response->Hresponse());
    TH2D* hSpectrumTruePerBin = dynamic_cast<TH2D*>(response->Htruth());
    TH2D* hSpectrumMCDetPerBin = dynamic_cast<TH2D*>(response->Hmeasured());
    
    if (!TH2DRM || !hSpectrumTruePerBin || !hSpectrumMCDetPerBin) {
        std::cerr << "Error: Failed to cast response matrix histograms to TH2D!" << std::endl;
        return;
    }
    
    // This is a closure test
    // Take MC det lvl spectrum to see if one can reliably unfold to the gen lvl spectrum
    RooUnfoldBayes bayesUnfold(response, hSpectrumMCDetPerBin, nIter);
    TH2D* hSpectrumMCDetPerBinUnfolded = dynamic_cast<TH2D*>(bayesUnfold.Hreco(
        static_cast<RooUnfold::ErrorTreatment>(errorType)));
    
    if (!hSpectrumMCDetPerBinUnfolded) {
        std::cerr << "Error: Unfolding failed in UnfoldingTest2D!" << std::endl;
        return;
    }
    
    // Draw the results
    int npTBins = pTBinsArrayTruth.size() - 1;
    TCanvas* c = new TCanvas("c", "c: pT", 800 * npTBins, 850);
    c->Divide(npTBins, 2, 0, 0);
    
    // Create legend
    TLegend* leg2 = new TLegend(0.2, 0.8, 0.5, 0.93, "");
    leg2->SetFillColor(10);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.04);
    
    for (int ptBin = 1; ptBin <= npTBins; ptBin++) {
        // Upper pads for histograms
        c->cd(ptBin);
        TVirtualPad* pad = c->cd(ptBin);
        pad->SetLeftMargin(0.15);
        pad->SetRightMargin(0.05);
        pad->SetTopMargin(0.05);
        pad->SetBottomMargin(0);
        
        // Project histograms to 1D for this pT bin
        TH1D* MCTruthPtBinProj = hSpectrumTruePerBin->ProjectionY(
            ("_pyOrig_" + std::to_string(ptBin)).c_str(), ptBin, ptBin);
        TH1D* MCMeasuredUnfoldedPtBinProj = hSpectrumMCDetPerBinUnfolded->ProjectionY(
            ("_pyRefolded_" + std::to_string(ptBin)).c_str(), ptBin, ptBin);
        
        MCMeasuredUnfoldedPtBinProj->GetYaxis()->SetTitle("counts per bin");
        MCMeasuredUnfoldedPtBinProj->GetYaxis()->SetTitleSize(0.06);
        MCMeasuredUnfoldedPtBinProj->GetYaxis()->SetLabelFont(43);
        MCMeasuredUnfoldedPtBinProj->GetYaxis()->SetLabelSize(20);
        MCMeasuredUnfoldedPtBinProj->SetNdivisions(505);
        MCMeasuredUnfoldedPtBinProj->SetLineColor(4);
        MCMeasuredUnfoldedPtBinProj->SetLineWidth(3);
        MCMeasuredUnfoldedPtBinProj->SetLineStyle(1);
        MCMeasuredUnfoldedPtBinProj->Draw("hist E");
        
        MCTruthPtBinProj->SetLineColor(1);
        MCTruthPtBinProj->SetLineWidth(2);
        MCTruthPtBinProj->SetLineStyle(1);
        MCTruthPtBinProj->Draw("hist E same");
        
        if (ptBin == 1) {
            leg2->AddEntry(MCTruthPtBinProj, "MC Gen. lvl spectrum", "l");
            leg2->AddEntry(MCMeasuredUnfoldedPtBinProj, "MC Det. lvl spectrum unfolded", "l");
            leg2->Draw("same");
        }
        
        // Lower pads for ratios
        c->cd(npTBins + ptBin);
        pad = c->cd(npTBins + ptBin);
        pad->SetLeftMargin(0.15);
        pad->SetRightMargin(0.05);
        pad->SetTopMargin(0);
        pad->SetBottomMargin(0.35);
        
        // Create ratio histogram
        TH1D* MCMeasuredUnfoldedPtBinProjClone = static_cast<TH1D*>(MCMeasuredUnfoldedPtBinProj->Clone(
            (std::string(MCMeasuredUnfoldedPtBinProj->GetName()) + "Clone").c_str()));
        MCMeasuredUnfoldedPtBinProjClone->Divide(MCTruthPtBinProj);
        MCMeasuredUnfoldedPtBinProjClone->SetMarkerStyle(21);
        MCMeasuredUnfoldedPtBinProjClone->SetMarkerColor(2);
        
        MCMeasuredUnfoldedPtBinProjClone->GetXaxis()->SetTitleSize(30);
        MCMeasuredUnfoldedPtBinProjClone->GetXaxis()->SetTitleFont(43);
        MCMeasuredUnfoldedPtBinProjClone->GetXaxis()->SetTitleOffset(3.0);
        MCMeasuredUnfoldedPtBinProjClone->GetXaxis()->SetLabelFont(43);
        MCMeasuredUnfoldedPtBinProjClone->GetXaxis()->SetLabelSize(20);
        MCMeasuredUnfoldedPtBinProjClone->GetYaxis()->SetRangeUser(0.5, 1.5);
        MCMeasuredUnfoldedPtBinProjClone->GetYaxis()->SetTitle(
            ("Ratio Unfolded to MC Gen lvl, nIter=" + std::to_string(nIter)).c_str());
        MCMeasuredUnfoldedPtBinProjClone->GetYaxis()->SetTitleSize(20);
        MCMeasuredUnfoldedPtBinProjClone->GetYaxis()->SetTitleFont(43);
        MCMeasuredUnfoldedPtBinProjClone->GetYaxis()->SetTitleOffset(2.2);
        MCMeasuredUnfoldedPtBinProjClone->GetYaxis()->SetLabelFont(43);
        MCMeasuredUnfoldedPtBinProjClone->GetYaxis()->SetLabelSize(20);
        MCMeasuredUnfoldedPtBinProjClone->GetYaxis()->SetNdivisions(505);
        MCMeasuredUnfoldedPtBinProjClone->GetXaxis()->SetTitle("z_{T}");
        MCMeasuredUnfoldedPtBinProjClone->Draw("hist E");
        
        // Draw reference lines
        std::vector<TLine*> lines = drawLines();
        for (TLine* line : lines) {
            TLine* cLine = static_cast<TLine*>(line->Clone("c"));
            cLine->Draw("same");
            delete cLine;
        }
    }
    
    // Save the canvas
    std::string fileName = outPath + "UnfoldingClosureTest2D_nIter" + std::to_string(nIter) + figureTag + ".png";
    c->SaveAs(fileName.c_str());
    
    // Clean up
    delete c;
    delete hSpectrumMCDetPerBinUnfolded;
}

void UnfoldSpectraClass::UnfoldingTest(RooUnfoldResponse* response, int nIter) {
    // Get the response matrix - use dynamic_cast for proper type conversion
    TH2D* TH2DRM = dynamic_cast<TH2D*>(response->Hresponse());
    
    if (!TH2DRM) {
        std::cerr << "Error: Failed to cast response matrix histogram to TH2D!" << std::endl;
        return;
    }
    
    // Get the truth-level jet spectrum (matched) from response matrix (already re-binned)
    TH1D* hSpectrumTruePerBin = TH2DRM->ProjectionY("_py", 1, TH2DRM->GetNbinsX());
    hSpectrumTruePerBin->SetName("hSpectrumTruePerBin");
    
    TH1D* hSpectrumMCDetPerBin = TH2DRM->ProjectionX("_px", 1, TH2DRM->GetNbinsY());
    hSpectrumMCDetPerBin->SetName("hSpectrumMCDetPerBin");
    
    // This is a closure test
    // Take MC det lvl spectrum to see if one can reliably unfold to the gen lvl spectrum
    RooUnfoldBayes bayesUnfold(response, hSpectrumMCDetPerBin, nIter);
    TH1D* hSpectrumMCDetPerBinUnfolded = dynamic_cast<TH1D*>(bayesUnfold.Hreco(
        static_cast<RooUnfold::ErrorTreatment>(errorType)));
    
    if (!hSpectrumMCDetPerBinUnfolded) {
        std::cerr << "Error: Unfolding failed in UnfoldingTest!" << std::endl;
        delete hSpectrumTruePerBin;
        delete hSpectrumMCDetPerBin;
        return;
    }
    
    // Draw the result
    TCanvas* c = new TCanvas("c", "c: pT", 800, 850);
    c->cd();
    c->SetLeftMargin(0.15);
    
    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.05);
    pad1->SetTopMargin(0.05);
    pad1->Draw();
    pad1->cd();
    
    hSpectrumMCDetPerBinUnfolded->GetYaxis()->SetTitle("counts per bin");
    hSpectrumMCDetPerBinUnfolded->GetYaxis()->SetTitleSize(0.06);
    hSpectrumMCDetPerBinUnfolded->GetYaxis()->SetLabelFont(43);
    hSpectrumMCDetPerBinUnfolded->GetYaxis()->SetLabelSize(20);
    hSpectrumMCDetPerBinUnfolded->SetNdivisions(505);
    hSpectrumMCDetPerBinUnfolded->SetLineColor(4);
    hSpectrumMCDetPerBinUnfolded->SetLineWidth(3);
    hSpectrumMCDetPerBinUnfolded->SetLineStyle(1);
    hSpectrumMCDetPerBinUnfolded->Draw("hist E");
    
    hSpectrumTruePerBin->SetLineColor(1);
    hSpectrumTruePerBin->SetLineWidth(2);
    hSpectrumTruePerBin->SetLineStyle(1);
    hSpectrumTruePerBin->Draw("hist E same");
    
    TLegend* leg2 = new TLegend(0.2, 0.8, 0.5, 0.93, "");
    leg2->SetFillColor(10);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.04);
    leg2->AddEntry(hSpectrumTruePerBin, "MC Gen. lvl spectrum", "l");
    leg2->AddEntry(hSpectrumMCDetPerBinUnfolded, "MC Det. lvl spectrum unfolded", "l");
    leg2->Draw("same");
    
    c->cd();
    
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.05);
    pad2->Draw();
    pad2->cd();
    
    TH1D* hSpectrumMCDetPerBinUnfoldedClone = static_cast<TH1D*>(hSpectrumMCDetPerBinUnfolded->Clone(
        (std::string(hSpectrumMCDetPerBinUnfolded->GetName()) + "Clone").c_str()));
    hSpectrumMCDetPerBinUnfoldedClone->Divide(hSpectrumTruePerBin);
    hSpectrumMCDetPerBinUnfoldedClone->SetMarkerStyle(21);
    hSpectrumMCDetPerBinUnfoldedClone->SetMarkerColor(2);
    
    hSpectrumMCDetPerBinUnfoldedClone->GetXaxis()->SetTitleSize(30);
    hSpectrumMCDetPerBinUnfoldedClone->GetXaxis()->SetTitleFont(43);
    hSpectrumMCDetPerBinUnfoldedClone->GetXaxis()->SetTitleOffset(3.0);
    hSpectrumMCDetPerBinUnfoldedClone->GetXaxis()->SetLabelFont(43);
    hSpectrumMCDetPerBinUnfoldedClone->GetXaxis()->SetLabelSize(20);
    hSpectrumMCDetPerBinUnfoldedClone->GetYaxis()->SetRangeUser(0, 2);
    hSpectrumMCDetPerBinUnfoldedClone->GetYaxis()->SetTitle(
        ("Ratio Unfolded to MC Gen lvl, nIter=" + std::to_string(nIter)).c_str());
    hSpectrumMCDetPerBinUnfoldedClone->GetYaxis()->SetTitleSize(20);
    hSpectrumMCDetPerBinUnfoldedClone->GetYaxis()->SetTitleFont(43);
    hSpectrumMCDetPerBinUnfoldedClone->GetYaxis()->SetTitleOffset(2.2);
    hSpectrumMCDetPerBinUnfoldedClone->GetYaxis()->SetLabelFont(43);
    hSpectrumMCDetPerBinUnfoldedClone->GetYaxis()->SetLabelSize(20);
    hSpectrumMCDetPerBinUnfoldedClone->GetYaxis()->SetNdivisions(505);
    hSpectrumMCDetPerBinUnfoldedClone->GetXaxis()->SetTitle("z_{T}");
    hSpectrumMCDetPerBinUnfoldedClone->Draw("hist E");
    
    // Draw reference lines
    std::vector<TLine*> lines = drawLines();
    for (TLine* line : lines) {
        TLine* cLine = static_cast<TLine*>(line->Clone("c"));
        cLine->Draw("same");
        delete cLine;
    }
    
    c->cd();
    
    std::string fileName = outPath + "UnfoldingClosureTest_nIter" + std::to_string(nIter) + figureTag + ".png";
    c->SaveAs(fileName.c_str());
    
    // Clean up
    delete c;
    delete hSpectrumMCDetPerBinUnfoldedClone;
    delete hSpectrumTruePerBin;
    delete hSpectrumMCDetPerBin;
    delete hSpectrumMCDetPerBinUnfolded;
}

std::vector<TLine*> UnfoldSpectraClass::drawLines() {
    std::vector<TLine*> lineList;
    
    TLine* line = new TLine(0, 1, 1, 1);
    line->SetLineColor(1);
    line->SetLineStyle(1);
    lineList.push_back(line);
    
    TLine* line2 = new TLine(0, 1.2, 1, 1.2);
    line2->SetLineColor(1);
    line2->SetLineStyle(2);
    lineList.push_back(line2);
    
    TLine* line3 = new TLine(0, 0.8, 1, 0.8);
    line3->SetLineColor(1);
    line3->SetLineStyle(2);
    lineList.push_back(line3);
    
    return lineList;
}

TH1D* UnfoldSpectraClass::scaleHistogram(TH1D* histo) {
    // This function takes a histogram and scales the bin content by bin width
    TH1D* histoScale = static_cast<TH1D*>(histo->Clone((std::string(histo->GetName()) + "Scaled").c_str()));
    histoScale->Reset("ICESM");
    
    for (int i = 1; i <= histo->GetNbinsX(); i++) {
        double num = histo->GetBinContent(i);
        double err = histo->GetBinError(i);
        double den = histo->GetBinWidth(i);
        double value = 0;
        double error = 0;
        
        if (den != 0) {
            value = num / den;
            error = err / den;
        }
        
        histoScale->SetBinContent(i, value);
        histoScale->SetBinError(i, error);
    }
    
    return histoScale;
}

std::string UnfoldSpectraClass::makeOutDir(const std::string& subDirName) {
    std::string subDirPath;
    
    if (outPath.back() != '/') {
        subDirPath = outPath + "/" + subDirName;
    } else {
        subDirPath = outPath + subDirName;
    }
    
    if (!std::filesystem::exists(subDirPath)) {
        std::filesystem::create_directories(subDirPath);
    }
    
    return subDirPath;
}

TH2D* UnfoldSpectraClass::SmearPoints(TH2D* hist, TRandom* fRandom, int factor) {
    TH2D* hnew = static_cast<TH2D*>(hist->Clone("hnew"));
    
    // The factor is an experimental factor.
    // What would happen if we had factor times more statistic
    // The new error with factor would be = sqrt(factor)*oldError
    // The bin counts would be factor*counts
    // So if we scale down to the old counts to have a comparatively RM
    // We have new bin content = content*factor/factor
    //         new error       = error*sqrt(factor)/factor
    
    for (int binx = 1; binx <= hnew->GetNbinsX(); binx++) {
        for (int biny = 1; biny <= hnew->GetNbinsY(); biny++) {
            double cont = hist->GetBinContent(binx, biny);
            double err = hist->GetBinError(binx, biny) * std::sqrt(factor) / factor;
            double newContent = fRandom->Gaus(cont, err); // 1Sigma variation
            
            // Ensure positive content
            if (newContent > 0) {
                hnew->SetBinContent(binx, biny, newContent);
                hnew->SetBinError(binx, biny, std::sqrt(newContent));
            } else {
                hnew->SetBinContent(binx, biny, 0);
                hnew->SetBinError(binx, biny, 0);
            }
        }
    }
    
    return hnew;
}
void UnfoldSpectraClass::saveResult(int regParam, int nIter, const std::string& tag) {
    // Check if regParam is valid
    if (regParam < 0 || regParam >= unfoldedArr2D.size()) {
        std::cerr << "Error: regParam " << regParam << " out of range for unfoldedArr2D (size: " 
                  << unfoldedArr2D.size() << ")" << std::endl;
        return;
    }
    
    // Create output file name
    std::string fileName = outPath + "Spectra_reg" + std::to_string(regParam) + ".root";
    TFile* fOutData = new TFile(fileName.c_str(), "RECREATE");
    
    if (!fOutData || fOutData->IsZombie()) {
        std::cerr << "Error: Could not create output file " << fileName << std::endl;
        if (fOutData) delete fOutData;
        return;
    }
    
    // Create flag based on prompt/non-prompt status
    std::string flag = isPrompt ? "P" : "NP";
    
    // Get histograms
    TH2D* mainUnfolded = unfoldedArr2D[regParam];
    TH2D* original = measuredSpectra2D;
    
    // Process for each pT bin
    for (size_t ptBin = 1; ptBin < pTBinsArrayTruth.size(); ptBin++) {
        // Get histos and save them
        TH1D* mainUnfoldedProj = mainUnfolded->ProjectionY(
            ("_pyMainUnfolded_Bin" + std::to_string(ptBin)).c_str(), ptBin, ptBin);
        TH1D* originalProj = original->ProjectionY(
            ("_pyOrig_Bin" + std::to_string(ptBin)).c_str(), ptBin, ptBin);
        
        // Scale histograms by bin width
        TH1D* mainUnfoldedS = scaleHistogram(mainUnfoldedProj);
        TH1D* originalS = scaleHistogram(originalProj);
        
        // Set names and write the histograms
        std::string origName = "hOrig" + flag + "_" + std::to_string(pTBinsArrayTruth[ptBin-1]) + 
                             "_" + std::to_string(pTBinsArrayTruth[ptBin]);
        originalS->SetName(origName.c_str());
        originalS->Write();
        
        std::string recoName = "hReco" + flag + "_" + std::to_string(pTBinsArrayTruth[ptBin-1]) + 
                             "_" + std::to_string(pTBinsArrayTruth[ptBin]);
        mainUnfoldedS->SetName(recoName.c_str());
        mainUnfoldedS->Write();
        
        // Clean up
        delete mainUnfoldedProj;
        delete originalProj;
        delete mainUnfoldedS;
        delete originalS;
    }
    
    // Process for each zT bin
    for (size_t ztBin = 1; ztBin < zBinsArrayTruth.size(); ztBin++) {
        // Get histos and save them
        TH1D* mainUnfoldedProj = mainUnfolded->ProjectionX(
            ("_pxMainUnfolded_Bin" + std::to_string(ztBin)).c_str(), ztBin, ztBin);
        TH1D* originalProj = original->ProjectionX(
            ("_pxOrig_Bin" + std::to_string(ztBin)).c_str(), ztBin, ztBin);
        
        // Scale histograms by bin width
        TH1D* mainUnfoldedSII = scaleHistogram(mainUnfoldedProj);
        TH1D* originalSII = scaleHistogram(originalProj);
        
        // Set names and write the histograms
        std::string origName = "hpTOrig" + flag + "_" + std::to_string(zBinsArrayTruth[ztBin-1]) + 
                             "_" + std::to_string(zBinsArrayTruth[ztBin]);
        originalSII->SetName(origName.c_str());
        originalSII->Write();
        
        std::string recoName = "hpTReco" + flag + "_" + std::to_string(zBinsArrayTruth[ztBin-1]) + 
                             "_" + std::to_string(zBinsArrayTruth[ztBin]);
        mainUnfoldedSII->SetName(recoName.c_str());
        mainUnfoldedSII->Write();
        
        // Clean up
        delete mainUnfoldedProj;
        delete originalProj;
        delete mainUnfoldedSII;
        delete originalSII;
    }
    
    // Close the file
    fOutData->Close();
    delete fOutData;
}

RooUnfoldResponse* UnfoldSpectraClass::PrepareResponseMatrix2D(int part) {
    // Get the RM from the input file
    TH2D* RM = getResponseMatrix(part, "zTDet", "zTPart");
    plotHist(RM, RM->GetName(), "colz");
    
    // Get truth-level spectrum (matched) from response matrix projection
    TH1D* hGenLevelMatchedUncutPerBin = RM->ProjectionY();
    hGenLevelMatchedUncutPerBin->SetName(("hGenLevelMatchedUncutPerBin_part" + std::to_string(part)).c_str());
    
    TH1D* hDetLevelMatchedUncutPerBin = RM->ProjectionX();
    hDetLevelMatchedUncutPerBin->SetName(("hDetLevelMatchedUncutPerBin" + std::to_string(part)).c_str());
    
    // Plot for QA
    plotHist(hGenLevelMatchedUncutPerBin, "hGenLevelMatchedUncutPerBin" + std::to_string(part), "E");
    plotHist(hDetLevelMatchedUncutPerBin, "hDetLevelMatchedUncutPerBin" + std::to_string(part), "E");
    
    // Plot jet response QA plots
    plotJetResponse();
    
    // Prepare a RooUnfold RM object
    // The first entry is empty since we do not account for fakes
    // For the second entry, we pass in the projection before the RM was cut to a specific range
    RooUnfoldResponse* RooUnfoldRM = new RooUnfoldResponse(
        0, 
        static_cast<TH1D*>(hGenLevelMatchedUncutPerBin->Clone("T1")), 
        static_cast<TH2D*>(RM->Clone("T")), 
        ("hResponseMatrixMain" + unfoldLabel).c_str(), 
        ("hResponseMatrixMain" + unfoldLabel).c_str());
    
    RooUnfoldRM->UseOverflow(false);
    
    return RooUnfoldRM;
}

void UnfoldSpectraClass::plotJetResponse() {
    // Create subdirectory for the response matrix QA
    std::string subdirpath = makeOutDir("Response");
    
    // In C++ version, this method would create plots of jet response characteristics
    // This would include jet energy resolution, jet energy scale shifts, efficiency, etc.
    // For this implementation, we'll leave it as a placeholder
    
    // Note: The full implementation would require more complex analysis of the response matrix
    // or additional input data that's not clear from the Python code
    std::cout << "plotJetResponse: Response QA directory created at " << subdirpath << std::endl;
}

void UnfoldSpectraClass::unfold1D(int ptBin, int regParam, int iteration, int powerLawOffset) {
    unfoldLabel = "1DBayes" + std::to_string(regParam) + "_Round" + std::to_string(iteration);
    
    // Create specific output directory
    std::string outputDirUnfolding = outPath + unfoldLabel;
    if (!std::filesystem::exists(outputDirUnfolding)) {
        std::cout << "directory does not exist. Make directory: " << outputDirUnfolding << std::endl;
        std::filesystem::create_directories(outputDirUnfolding);
    }
    
    // Prepare the response matrices
    RooUnfoldResponse* responseMatrix = PrepareResponseMatrix2D();
    RooUnfoldResponse* responseMatrix1 = PrepareResponseMatrix2D(1);
    RooUnfoldResponse* responseMatrix2 = PrepareResponseMatrix2D(2);
    
    // Call RooUnfold for different regularization parameters
    for (int i = 1; i <= regParam + 2; i++) {
        // Set up the Bayesian unfolding object
        RooUnfoldBayes unfoldBayes(responseMatrix, measuredSpectraArray[ptBin], i);
        
        // Perform the unfolding
        TH1D* hSpectrumUnfoldedPerBin = dynamic_cast<TH1D*>(unfoldBayes.Hreco(
            static_cast<RooUnfold::ErrorTreatment>(errorType)));
        
        if (!hSpectrumUnfoldedPerBin) {
            std::cerr << "Error: Unfolding failed for iteration " << i << std::endl;
            continue;
        }
        
        hSpectrumUnfoldedPerBin->SetName(("UnfoldedSpectraPerBin_nIter" + std::to_string(i)).c_str());
        plotHist(hSpectrumUnfoldedPerBin, "UnfoldedSpectraBin" + std::to_string(ptBin) + 
                "_nIter" + std::to_string(i), "E");
        unfoldedArrPerBin[ptBin].push_back(hSpectrumUnfoldedPerBin);
        
        // Plot Pearson correlation coefficients
        TMatrixD covarianceMatrix = unfoldBayes.Ereco(
            static_cast<RooUnfold::ErrorTreatment>(errorType));
        plotCorrelationCoefficients(covarianceMatrix, i, "CorMatr_Bin" + std::to_string(ptBin));
        
        // Refolding Test
        RefoldingTest(ptBin, i, measuredSpectraArray[ptBin], responseMatrix1, responseMatrix2);
        
        // Unfolding test
        UnfoldingTest(responseMatrix, i);
    }
    
    // Stability tests
    StabilityTest(ptBin, regParam);
    StabilityTest(ptBin, regParam + 1);
    StabilityTest(ptBin, regParam, true);
    
    // Test of the robustness of the RM in terms of statistics
    StatTestRM(responseMatrix, ptBin);
    
    // Plot the result before and after unfolding
    plotUnfoldingEffect(ptBin, regParam);
    
    // Clean up
    delete responseMatrix;
    delete responseMatrix1;
    delete responseMatrix2;
}
void UnfoldSpectraClass::StatTestRM2D(RooUnfoldResponse* response) {
    std::cout << "\n==== StatTestRM2D DEBUG START ====" << std::endl;
    
    // Get response matrix histograms
    TH2D* TH2DRM = dynamic_cast<TH2D*>(response->Hresponse());
    TH2D* TH2DTruth = dynamic_cast<TH2D*>(response->Htruth());
    TH2D* TH2DMeas = dynamic_cast<TH2D*>(response->Hmeasured());
    
    if (!TH2DRM || !TH2DTruth || !TH2DMeas) {
        std::cerr << "ERROR: Failed to cast response matrix histograms!" << std::endl;
        return;
    }
    
    // Get dimensions and setup
    int npTBins = pTBinsArrayTruth.size() - 1;
    int RegulForStatTest = 3;  // Fixed regularization parameter for this test
    int errorType = RooUnfold::kCovariance;  // Only for this test
    TRandom3* fRandom = new TRandom3(0);  // Random generator with fixed seed
    int factorStatistic = 1;  // Test factor for statistics in the RM
    
    std::cout << "Beginning statistical test with " << npTBins << " pT bins" << std::endl;
    std::cout << "Will use regularization parameter " << RegulForStatTest << std::endl;
    
    // Store unfolded histograms from statistical variations
    std::vector<TH2D*> unfoldTesthistograms;
    const int nVariations = 1000;
    
    std::cout << "Will perform " << nVariations << " random variations of response matrix" << std::endl;
    
    // Loop over variations
    int successfulVariations = 0;
    for (int i = 1; i <= nVariations; i++) {
        if (i % 100 == 0) {
            std::cout << "Processing variation " << i << "/" << nVariations << std::endl;
        }
        
        // Smear the response matrix within statistical uncertainties
        TH2D* RMVariation = SmearPoints(TH2DRM, fRandom, factorStatistic);
        
        // Setup empty RooUnfoldResponse object with the varied matrix
        RooUnfoldResponse* RooUnfoldRMVar = new RooUnfoldResponse(
            ("hResponseMatrix3D" + std::to_string(i)).c_str(), 
            ("hResponseMatrix3D" + std::to_string(i)).c_str());
        
        RooUnfoldRMVar->Setup(TH2DMeas, TH2DTruth, RMVariation);
        
        // Perform unfolding with the varied response matrix
        RooUnfoldBayes unfoldVar(RooUnfoldRMVar, measuredSpectra2D, RegulForStatTest);
        TH2D* hunf = dynamic_cast<TH2D*>(unfoldVar.Hreco(
            static_cast<RooUnfold::ErrorTreatment>(errorType)));
            
        if (!hunf) {
            std::cerr << "ERROR: Variation " << i << " unfolding failed!" << std::endl;
            delete RMVariation;
            delete RooUnfoldRMVar;
            continue;
        }
        
        unfoldTesthistograms.push_back(hunf);
        successfulVariations++;
        
        // Clean up temporary objects
        delete RMVariation;
        delete RooUnfoldRMVar;
    }
    
    std::cout << "Completed " << successfulVariations << " successful variations out of " 
              << nVariations << " attempted" << std::endl;
    
    // Create empty histograms to store mean and std dev for each pT bin
    std::vector<TH1D*> hvariationResult;
    std::vector<TH1D*> hvariationError;
    
    // Create template histogram
    TH1D* histo = unfoldTesthistograms[0]->ProjectionY("_pymain", 1, 1);
    histo->Reset("ICESM");
    
    // Initialize result arrays
    for (int ptBin = 0; ptBin < npTBins; ptBin++) {
        hvariationResult.push_back(static_cast<TH1D*>(histo->Clone(
            ("pTBinCont_" + std::to_string(ptBin)).c_str())));
        hvariationError.push_back(static_cast<TH1D*>(histo->Clone(
            ("pTBinErr_" + std::to_string(ptBin)).c_str())));
    }
    
    // Get number of zT bins from the first unfolded histogram
    int nzTBins = unfoldTesthistograms[0]->ProjectionY("_pymain", 1, 1)->GetNbinsX();
    
    // Calculate dispersion of the unfolded results
    std::cout << "Calculating statistical dispersion from " << unfoldTesthistograms.size() 
              << " variations..." << std::endl;
              
    std::tie(hvariationResult, hvariationError) = calculateDispersion(
        unfoldTesthistograms, hvariationResult, hvariationError, nzTBins, npTBins);
    
    // Create canvas for plots
    TCanvas* c = new TCanvas("c", "c: pT", 800 * npTBins, 850);
    c->Divide(npTBins, 2, 0, 0);
    
    // Plot results for each pT bin
    for (int ptBin = 1; ptBin <= npTBins; ptBin++) {
        std::cout << "Processing pT bin " << ptBin << "/" << npTBins << std::endl;
        // Upper pad for showing the mean values
        c->cd(ptBin);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.05);
        gPad->SetBottomMargin(0);
        
        hvariationResult[ptBin-1]->GetYaxis()->SetTitle("counts per bin");
        hvariationResult[ptBin-1]->GetYaxis()->SetTitleSize(0.06);
        hvariationResult[ptBin-1]->GetYaxis()->SetLabelFont(43);
        hvariationResult[ptBin-1]->GetYaxis()->SetLabelSize(20);
        hvariationResult[ptBin-1]->SetNdivisions(505);
        hvariationResult[ptBin-1]->SetLineColor(4);
        hvariationResult[ptBin-1]->SetLineWidth(3);
        hvariationResult[ptBin-1]->SetLineStyle(1);
        hvariationResult[ptBin-1]->SetMarkerStyle(21);
        hvariationResult[ptBin-1]->Draw("hist E");
        
        // Add a label for the pT bin
        TPaveText* ptLabel = new TPaveText(0.15, 0.92, 0.7, 0.99, "NDC");
        ptLabel->SetBorderSize(0);
        ptLabel->SetFillStyle(0);
        ptLabel->SetTextFont(42);
        ptLabel->SetTextSize(0.04);
        ptLabel->AddText(Form("p_{T} bin: %d-%d GeV/c", 
                            pTBinsArrayTruth[ptBin-1], pTBinsArrayTruth[ptBin]));
        ptLabel->Draw();
        
        // Lower pad for showing relative errors
        c->cd(npTBins + ptBin);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0);
        gPad->SetBottomMargin(0.35);
        
        hvariationError[ptBin-1]->GetYaxis()->SetRangeUser(0.5, 1.5);  // Adjust range as needed
        hvariationError[ptBin-1]->GetYaxis()->SetTitle("Relative Uncertainty");
        hvariationError[ptBin-1]->GetXaxis()->SetTitle("z_{T}");
        hvariationError[ptBin-1]->SetLineColor(4);
        hvariationError[ptBin-1]->SetLineWidth(3);
        hvariationError[ptBin-1]->SetLineStyle(1);
        hvariationError[ptBin-1]->SetMarkerStyle(21);
        hvariationError[ptBin-1]->Draw("hist E");
        
        // Draw reference lines
        std::vector<TLine*> lines = drawLines();
        for (auto line : lines) {
            TLine* cLine = static_cast<TLine*>(line->Clone("c"));
            cLine->Draw("same");
            delete cLine;
        }
    }
    
    // Save canvas
    std::string fileName = outPath + "RMStatTest2D_nIter" + std::to_string(RegulForStatTest) + 
                         "_" + figureTag + ".png";
    c->SaveAs(fileName.c_str());
    
    // Clean up
    // delete c;
    // delete fRandom;
    // delete histo;
    
    // Clean up variation histograms
    for (auto hist : unfoldTesthistograms) {
        delete hist;
    }
    unfoldTesthistograms.clear();
    
    // Clean up result histograms
    for (auto hist : hvariationResult) {
        delete hist;
    }
    hvariationResult.clear();
    
    for (auto hist : hvariationError) {
        delete hist;
    }
    hvariationError.clear();
    
    std::cout << "==== StatTestRM2D DEBUG END ====\n" << std::endl;
}

void UnfoldSpectraClass::plotPrior2D(int nIter, RooUnfoldResponse* RooUnfoldRM) {
    // Set ROOT style parameters
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    
    // Set axis title based on particle type
    std::string xTitle = "p_{T}^{D^{0}}/p_{T}^{jet}";
    
    // Clone histograms for working copies
    TH2D* mainUnfolded = static_cast<TH2D*>(unfoldedArr2D[nIter-1]->Clone("TUnfold"));
    TH2D* original = static_cast<TH2D*>(measuredSpectra2D->Clone("TOrig"));
    TH2D* prior = static_cast<TH2D*>(RooUnfoldRM->Htruth()->Clone("TTruth"));
    
    // Define color arrays
    std::vector<int> colorArray = {kRed-9, kGreen-9, kGreen-8, kBlue-9, kGreen, kSpring+8, kSpring};
    std::vector<int> colorArrayFinalFigs = {kAzure, kAzure-4, kCyan-6, kGreen-3, kTeal-6, kGreen+4, 1, 1};
    
    // Arrays to store histogram properties
    std::vector<double> maxArrayOrig;
    std::vector<double> maxArrayUnfold;
    std::vector<TH1D*> hArrayOrig;
    std::vector<TH1D*> hArrayUnfold;
    
    // Create legend
    TLegend* leg1 = new TLegend(0.2, 0.65, 0.5, 0.85, "");
    leg1->SetFillColor(10);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.06);
    
    // Create canvas based on particle type
    int npTBins = pTBinsArrayTruth.size() - 1;
    TCanvas* c1 = new TCanvas("c1", "c1: pT", 800*npTBins, 800);
    c1->Divide(npTBins, 1, 0, 0);
    
    // Process each pT bin
    for (int ptBin = 1; ptBin <= npTBins; ptBin++) {
        c1->cd(ptBin);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.05);
        gPad->SetBottomMargin(0.15);
        
        // Create projections for this pT bin
        TH1D* mainUnfoldedProj = mainUnfolded->ProjectionY(
            ("_pyMainUnfolded_Bin" + std::to_string(ptBin)).c_str(), ptBin, ptBin);
        TH1D* originalProj = original->ProjectionY(
            ("_pyOrig_Bin" + std::to_string(ptBin)).c_str(), ptBin, ptBin);
        TH1D* priorProj = prior->ProjectionY(
            ("_pyPrior_Bin" + std::to_string(ptBin)).c_str(), ptBin, ptBin);
        
        // Scale histograms by bin width
        TH1D* mainUnfoldedS = scaleHistogram(mainUnfoldedProj);
        TH1D* originalS = scaleHistogram(originalProj);
        TH1D* priorS = scaleHistogram(priorProj);
        
        // Scale prior to match unfolded integral
        if (priorS->Integral() > 0) {
            priorS->Scale(mainUnfoldedS->Integral() / priorS->Integral());
        }
        
        // Set histogram display properties
        originalS->SetNdivisions(505);
        
        // Calculate maximum value for y-axis scaling
        std::vector<double> maxList = {originalS->GetMaximum(), priorS->GetMaximum(), mainUnfoldedS->GetMaximum()};
        double maxY = *std::max_element(maxList.begin(), maxList.end());
        
        // Style the original histogram
        originalS->GetYaxis()->SetTitle("dN/dz_{T}");
        originalS->GetYaxis()->SetTitleSize(0.06);
        originalS->GetYaxis()->SetLabelFont(43);
        originalS->GetYaxis()->SetLabelSize(20);
        originalS->GetXaxis()->SetTitle(xTitle.c_str());
        originalS->GetXaxis()->SetTitleSize(0.05);
        originalS->GetXaxis()->SetRangeUser(0, 1);
        originalS->GetXaxis()->SetNdivisions(405);
        originalS->GetYaxis()->SetRangeUser(0, maxY * 1.3);
        originalS->SetLineColor(4);
        originalS->SetMarkerColor(4);
        originalS->SetLineWidth(3);
        originalS->SetMarkerStyle(20);
        originalS->SetLineStyle(1);
        originalS->DrawCopy("hist E");
        
        // Style the unfolded histogram
        mainUnfoldedS->SetLineColor(2);
        mainUnfoldedS->SetMarkerColor(2);
        mainUnfoldedS->SetLineWidth(2);
        mainUnfoldedS->SetMarkerStyle(20);
        mainUnfoldedS->SetLineStyle(1);
        mainUnfoldedS->DrawCopy("same hist E");
        
        // Style the prior histogram
        priorS->SetLineColor(kSpring+3);
        priorS->SetMarkerColor(kSpring+3);
        priorS->SetLineWidth(2);
        priorS->SetLineStyle(1);
        priorS->SetMarkerStyle(20);
        priorS->DrawCopy("hist E same");
        
        // Store histograms and maximum values in arrays
        maxArrayOrig.push_back(originalS->GetMaximum());
        maxArrayUnfold.push_back(mainUnfoldedS->GetMaximum());
        hArrayOrig.push_back(originalS);
        hArrayUnfold.push_back(mainUnfoldedS);
        
        // Create legend entries for the first pT bin only
        if (ptBin == 1) {
            leg1->AddEntry(originalS, "measured", "l");
            leg1->AddEntry(mainUnfoldedS, ("unfolded nIter=" + std::to_string(nIter)).c_str(), "l");
            leg1->AddEntry(priorS, "prior", "l");
        }
        leg1->Draw("same");
        
        // Add pT bin label
        TPaveText* textFit = new TPaveText(0.2, 0.87, 0.7, 0.92, "NDC");
        textFit->SetBorderSize(0);
        textFit->SetFillStyle(0);
        textFit->SetTextSize(0.06);
        textFit->AddText(Form(" #it{p}_{T}^{jet}: %d-%d GeV", 
                            pTBinsArrayTruth[ptBin-1], pTBinsArrayTruth[ptBin]));
        textFit->Draw("same");
        
        // Clean up projections - we keep scaled versions in arrays
        delete mainUnfoldedProj;
        delete originalProj;
        delete priorProj;
    }
    
    // Save the canvas
    c1->cd();
    std::string fileName = outPath + "PriorComparision2D_Yaxis_nIter" + std::to_string(nIter) + figureTag + ".png";
    c1->SaveAs(fileName.c_str());
    
    // Create a second legend for the zT projections
    TLegend* leg2 = new TLegend(0.5, 0.7, 0.8, 0.93, "");
    leg2->SetFillColor(10);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.06);
    
    // Get number of zT bins
    int nzTBins = zBinsArrayTruth.size() - 1;
    
    // Create canvas based on particle type
    TCanvas* c2 = new TCanvas("c2", "c2: pT", 800*3, 850*2);
    c2->Divide(3, 2);
    
    // Process each zT bin
    for (int zTBin = 1; zTBin <= nzTBins; zTBin++) {
        c2->cd(zTBin);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.05);
        gPad->SetBottomMargin(0.15);
        
        // Create projections for this zT bin
        TH1D* mainUnfoldedProj = mainUnfolded->ProjectionX(
            ("_pxMainUnfolded_Bin" + std::to_string(zTBin)).c_str(), zTBin, zTBin);
        TH1D* originalProj = original->ProjectionX(
            ("_pxOrig_Bin" + std::to_string(zTBin)).c_str(), zTBin, zTBin);
        TH1D* priorProj = prior->ProjectionX(
            ("_pxPrior_Bin" + std::to_string(zTBin)).c_str(), zTBin, zTBin);
        
        // Scale histograms by bin width
        TH1D* mainUnfoldedS = scaleHistogram(mainUnfoldedProj);
        TH1D* originalS = scaleHistogram(originalProj);
        TH1D* priorS = scaleHistogram(priorProj);
        
        // Scale prior to match unfolded integral
        if (priorS->Integral() > 0) {
            priorS->Scale(mainUnfoldedS->Integral() / priorS->Integral());
        }
        
        // Calculate maximum value for y-axis scaling
        std::vector<double> maxList = {originalS->GetMaximum(), priorS->GetMaximum(), mainUnfoldedS->GetMaximum()};
        double maxY = *std::max_element(maxList.begin(), maxList.end());
        
        // Style the unfolded histogram
        mainUnfoldedS->GetYaxis()->SetTitle("dN/dp_{T}");
        mainUnfoldedS->GetYaxis()->SetTitleSize(0.06);
        mainUnfoldedS->GetYaxis()->SetLabelFont(43);
        mainUnfoldedS->GetYaxis()->SetLabelSize(20);
        mainUnfoldedS->GetXaxis()->SetTitle("jet p_{T}");
        mainUnfoldedS->GetXaxis()->SetTitleSize(0.05);
        mainUnfoldedS->GetXaxis()->SetRangeUser(0, 1);
        mainUnfoldedS->GetXaxis()->SetNdivisions(405);
        mainUnfoldedS->GetYaxis()->SetRangeUser(0, maxY * 1.3);
        
        mainUnfoldedS->SetNdivisions(505);
        mainUnfoldedS->SetLineColor(2);
        mainUnfoldedS->SetMarkerColor(2);
        mainUnfoldedS->SetLineWidth(2);
        mainUnfoldedS->SetMarkerStyle(20);
        mainUnfoldedS->SetLineStyle(1);
        mainUnfoldedS->Draw("hist E");
        
        // Style the original histogram
        originalS->SetLineColor(4);
        originalS->SetMarkerColor(4);
        originalS->SetLineWidth(3);
        originalS->SetMarkerStyle(20);
        originalS->SetLineStyle(1);
        originalS->Draw("hist E same");
        
        // Style the prior histogram
        priorS->SetLineColor(kSpring+3);
        priorS->SetMarkerColor(kSpring+3);
        priorS->SetLineWidth(2);
        priorS->SetLineStyle(1);
        priorS->SetMarkerStyle(20);
        priorS->Draw("hist E same");
        
        // Create legend entries for the first bin only
        if (zTBin == 1) {
            leg2->AddEntry(originalS, "measured", "l");
            leg2->AddEntry(mainUnfoldedS, ("unfolded nIter=" + std::to_string(nIter)).c_str(), "l");
            leg2->AddEntry(priorS, "prior", "l");
        }
        leg2->Draw("same");
        
        // Add zT bin label
        TPaveText* textFit = new TPaveText(0.4, 0.8, 0.7, 0.85, "NDC");
        textFit->SetBorderSize(0);
        textFit->SetFillStyle(0);
        textFit->SetTextSize(0.06);
        textFit->AddText(Form(" #it{z}_{T}: %.2f-%.2f", 
                            zBinsArrayTruth[zTBin-1], zBinsArrayTruth[zTBin]));
        textFit->Draw("same");
        
        // Clean up projections
        delete mainUnfoldedProj;
        delete originalProj;
        delete priorProj;
        delete mainUnfoldedS;
        delete originalS;
        delete priorS;
    }
    
    // Save the second canvas
    c2->cd();
    std::string fileName2 = outPath + "PriorComparision2D_Xaxis_nIter" + std::to_string(nIter) + figureTag + ".png";
    c2->SaveAs(fileName2.c_str());
    
    // Create a third canvas for summary plot
    TCanvas* c3 = new TCanvas("c3", "c3: hist", 500*2, 450*2);
    c3->cd();
    TGaxis::SetMaxDigits(3);
    
    // Set pad and histogram arrangement
    TPad* myPad3 = new TPad("myPad3", "The pad", 0, 0, 1, 1);
    myPad3->SetLeftMargin(0.15);
    myPad3->SetTopMargin(0.06);
    myPad3->SetRightMargin(0.04);
    myPad3->SetBottomMargin(0.15);
    myPad3->SetTicks();
    myPad3->Draw();
    myPad3->cd();
    
    // Find maximum value
    double maximum = 0;
    if (!maxArrayUnfold.empty()) {
        maximum = *std::max_element(maxArrayUnfold.begin(), maxArrayUnfold.end());
    }
    
    // Create blank histogram for styling
    TH1F* myBlankHisto3 = new TH1F("myBlankHisto3", "Blank Histogram", 20, 0, 1);
    myBlankHisto3->GetXaxis()->SetNdivisions(505);
    myBlankHisto3->SetXTitle(xTitle.c_str());
    myBlankHisto3->GetXaxis()->SetTitleSize(0.05);
    myBlankHisto3->GetXaxis()->SetRangeUser(0, 1);
    myBlankHisto3->GetXaxis()->SetNdivisions(405);
    myBlankHisto3->GetYaxis()->SetTitleOffset(1.35);
    myBlankHisto3->GetYaxis()->SetTitleSize(0.055);
    myBlankHisto3->SetLineColor(0);
    
    // Set y-axis properties
    myBlankHisto3->SetYTitle("dN/dz_{T} unfolded");
    if (zBinsArrayTruth.size() == zBinsArrayDet.size()) {
        myBlankHisto3->GetYaxis()->SetRangeUser(0, 0.45); // When drawn normalized
    } else {
        myBlankHisto3->GetYaxis()->SetRangeUser(0, 0.8); // When drawn normalized
    }
    myBlankHisto3->Draw("E");
    
    // Create legend
    TLegend* myLegend3 = nullptr;
    if (hArrayUnfold.size() == 4) {
        myLegend3 = new TLegend(0.2, 0.7, 0.4, 0.9);
    } else if (hArrayUnfold.size() == 6) {
        myLegend3 = new TLegend(0.2, 0.6, 0.4, 0.9);
    } else {
        myLegend3 = new TLegend(0.2, 0.7, 0.4, 0.9);
    }
    
    myLegend3->SetTextFont(42);
    myLegend3->SetBorderSize(0);
    myLegend3->SetFillStyle(0);
    myLegend3->SetFillColor(0);
    myLegend3->SetMargin(0.25);
    myLegend3->SetTextSize(0.04);
    
    if (!hArrayUnfold.empty()) {
        myLegend3->AddEntry(hArrayUnfold[0], " #it{p}_{T}^{jet}:", "");
    }
    
    // Draw all unfolded histograms normalized
    double markerScale = 1.6;
    for (size_t i = 0; i < hArrayUnfold.size(); i++) {
        hArrayUnfold[i]->SetMarkerSize(1.3 * markerScale);
        hArrayUnfold[i]->SetMarkerStyle(20);
        hArrayUnfold[i]->SetMarkerColor(colorArrayFinalFigs[i % colorArrayFinalFigs.size()]);
        hArrayUnfold[i]->SetLineStyle(1);
        hArrayUnfold[i]->SetLineWidth(2);
        hArrayUnfold[i]->SetLineColor(colorArrayFinalFigs[i % colorArrayFinalFigs.size()]);
        hArrayUnfold[i]->DrawNormalized("same EP");
        
        myLegend3->AddEntry(hArrayUnfold[i], 
            Form("  %d-%d (GeV/%s)", pTBinsArrayTruth[i], pTBinsArrayTruth[i+1], "#it{c}"), "LP");
    }
    
    myLegend3->Draw();
    
    // Save the canvas
    c3->cd();
    std::string fileName3 = outPath + "UnfoldingEffect2DFinal_" + figureTag + ".png";
    c3->SaveAs(fileName3.c_str());
    
    // Clean up
    // delete c1;
    // delete c2;
    // delete c3;
    // delete leg1;
    // delete leg2;
    // delete myPad3;
    // delete myBlankHisto3;
    // delete myLegend3;
    // delete mainUnfolded;
    // delete original;
    // delete prior;
    
    // Clean up histogram arrays
    for (auto hist : hArrayOrig) {
        delete hist;
    }
    hArrayOrig.clear();
    
    for (auto hist : hArrayUnfold) {
        delete hist;
    }
    hArrayUnfold.clear();
}

void UnfoldSpectraClass::plotUnfoldingEffect2D(int nIter) {
    // Set ROOT style parameters
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    
    // Set axis title based on particle type
    std::string xTitle = "p_{T}^{D^{0}}/p_{T}^{jet}";
    
    // Clone histograms for working copies
    TH2D* mainUnfolded = static_cast<TH2D*>(unfoldedArr2D[nIter-1]->Clone("TUnfold"));
    TH2D* original = static_cast<TH2D*>(measuredSpectra2D->Clone("TOrig"));
    
    // Define color arrays for plotting
    std::vector<int> colorArray = {kRed-9, kGreen-9, kGreen-8, kBlue-9, kGreen, kSpring+8, kSpring};
    std::vector<int> colorArrayFinalFigs = {kAzure, kAzure-4, kCyan-6, kGreen-3, kTeal-6, kGreen+4, 1, 1};
    
    // Arrays to store histogram properties
    std::vector<double> maxArrayOrig;
    std::vector<double> maxArrayUnfold;
    std::vector<TH1D*> hArrayOrig;
    std::vector<TH1D*> hArrayUnfold;
    
    // Create legend
    TLegend* leg1 = new TLegend(0.2, 0.7, 0.5, 0.85, "");
    leg1->SetFillColor(10);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.06);
    
    // Create canvas based on particle type
    int npTBins = pTBinsArrayTruth.size() - 1;
    TCanvas* c1= new TCanvas("c1", "c1: pT", 800*npTBins, 800);
    c1->Divide(npTBins, 1);
    
    // Process each pT bin
    for (int ptBin = 1; ptBin <= npTBins; ptBin++) {
        c1->cd(ptBin);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.05);
        gPad->SetBottomMargin(0.15);
        
        // Create projections for this pT bin
        TH1D* mainUnfoldedProj = mainUnfolded->ProjectionY(
            ("_pyMainUnfolded_Bin" + std::to_string(ptBin)).c_str(), ptBin, ptBin);
        TH1D* originalProj = original->ProjectionY(
            ("_pyOrig_Bin" + std::to_string(ptBin)).c_str(), ptBin, ptBin);
        
        // Scale histograms by bin width
        TH1D* mainUnfoldedS = scaleHistogram(mainUnfoldedProj);
        TH1D* originalS = scaleHistogram(originalProj);
        
        // Find maximum value for y-axis scaling
        std::vector<double> maxList = {originalS->GetMaximum(), mainUnfoldedS->GetMaximum()};
        double maxY = *std::max_element(maxList.begin(), maxList.end());
        originalS->GetYaxis()->SetRangeUser(0, maxY * 1.3);
        
        // Style the original histogram
        originalS->GetYaxis()->SetTitle("dN/dz_{T}");
        originalS->GetYaxis()->SetTitleSize(0.06);
        originalS->GetYaxis()->SetLabelFont(43);
        originalS->GetYaxis()->SetLabelSize(20);
        originalS->GetXaxis()->SetTitle(xTitle.c_str());
        originalS->GetXaxis()->SetTitleSize(0.05);
        originalS->GetXaxis()->SetRangeUser(0, 1);
        originalS->GetXaxis()->SetNdivisions(405);
        originalS->SetNdivisions(505);
        originalS->SetLineColor(4);
        originalS->SetLineWidth(3);
        originalS->SetLineStyle(1);
        originalS->DrawCopy("hist E");
        
        // Style the unfolded histogram
        mainUnfoldedS->SetLineColor(2);
        mainUnfoldedS->SetLineWidth(2);
        mainUnfoldedS->SetLineStyle(1);
        mainUnfoldedS->DrawCopy("hist E same");
        
        // Store histograms and maximum values in arrays
        maxArrayOrig.push_back(originalS->GetMaximum());
        maxArrayUnfold.push_back(mainUnfoldedS->GetMaximum());
        hArrayOrig.push_back(originalS);
        hArrayUnfold.push_back(mainUnfoldedS);
        
        // Create legend entries for the first pT bin only
        if (ptBin == 1) {
            leg1->AddEntry(originalS, "measured", "l");
            leg1->AddEntry(mainUnfoldedS, ("unfolded nIter=" + std::to_string(nIter)).c_str(), "l");
        }
        leg1->Draw("same");
        
        // Add pT bin label
        TPaveText* textFit = new TPaveText(0.2, 0.87, 0.7, 0.92, "NDC");
        textFit->SetBorderSize(0);
        textFit->SetFillStyle(0);
        textFit->SetTextSize(0.06);
        textFit->AddText(Form(" #it{p}_{T}^{jet}: %d-%d GeV", 
                            pTBinsArrayTruth[ptBin-1], pTBinsArrayTruth[ptBin]));
        textFit->Draw("same");
        
        // Clean up projections - we keep scaled versions in arrays
        delete mainUnfoldedProj;
        delete originalProj;
    }
    
    // Save the canvas
    std::string fileName = outPath + "UnfoldingEffect2D_nIter" + std::to_string(nIter) + figureTag + ".png";
    c1->SaveAs(fileName.c_str());
    
    // Create second canvas for summary of measured spectra
    TCanvas* c2 = new TCanvas("c2", "c: hist", 1000, 900);
    c2->cd();
    TGaxis::SetMaxDigits(3);
    
    // Set pad and histogram arrangement
    TPad* myPad2 = new TPad("myPad2", "The pad", 0, 0, 1, 1);
    myPad2->SetLeftMargin(0.15);
    myPad2->SetTopMargin(0.06);
    myPad2->SetRightMargin(0.04);
    myPad2->SetBottomMargin(0.15);
    myPad2->SetTicks();
    myPad2->Draw();
    myPad2->cd();
    
    // Create blank histogram for styling
    TH1F* myBlankHisto2 = new TH1F("myBlankHisto2", "Blank Histogram", 20, 0, 1);
    myBlankHisto2->GetXaxis()->SetNdivisions(505);
    myBlankHisto2->SetXTitle(xTitle.c_str());
    myBlankHisto2->GetXaxis()->SetTitleSize(0.05);
    myBlankHisto2->GetXaxis()->SetRangeUser(0, 1);
    myBlankHisto2->GetXaxis()->SetNdivisions(405);
    myBlankHisto2->GetYaxis()->SetTitleOffset(1.35);
    myBlankHisto2->GetYaxis()->SetTitleSize(0.055);
    myBlankHisto2->SetLineColor(0);
    
    // Set y-axis properties
    myBlankHisto2->SetYTitle("dN/dz_{T} measured");
    myBlankHisto2->GetYaxis()->SetRangeUser(0, 0.45); // When drawn normalized
    myBlankHisto2->Draw("E");
    
    // Create legend for measured spectra
    TLegend* myLegend2 = nullptr;
    if (hArrayOrig.size() == 4) {
        myLegend2 = new TLegend(0.2, 0.7, 0.4, 0.9);
    } else if (hArrayOrig.size() == 6) {
        myLegend2 = new TLegend(0.2, 0.6, 0.4, 0.9);
    } else {
        myLegend2 = new TLegend(0.2, 0.7, 0.4, 0.9);
    }
    
    myLegend2->SetTextFont(42);
    myLegend2->SetBorderSize(0);
    myLegend2->SetFillStyle(0);
    myLegend2->SetFillColor(0);
    myLegend2->SetMargin(0.25);
    myLegend2->SetTextSize(0.04);
    myLegend2->AddEntry(hArrayOrig[0], " #it{p}_{T}^{jet}:", "");
    
    // Draw all measured histograms normalized
    double markerScale = 1.6;
    for (size_t i = 0; i < hArrayOrig.size(); i++) {
        hArrayOrig[i]->SetMarkerSize(1.3 * markerScale);
        hArrayOrig[i]->SetMarkerStyle(20);
        hArrayOrig[i]->SetMarkerColor(colorArrayFinalFigs[i % colorArrayFinalFigs.size()]);
        hArrayOrig[i]->SetLineStyle(1);
        hArrayOrig[i]->SetLineWidth(2);
        hArrayOrig[i]->SetLineColor(colorArrayFinalFigs[i % colorArrayFinalFigs.size()]);
        hArrayOrig[i]->DrawNormalized("same EP");
        
        myLegend2->AddEntry(hArrayOrig[i], 
            Form("  %d-%d (GeV/%s)", pTBinsArrayDet[i], pTBinsArrayDet[i+1], "#it{c}"), "LP");
    }
    
    myLegend2->Draw();
    
    // Save the canvas
    std::string fileName2 = outPath + "UnfoldingEffect2DOrig_" + figureTag + ".png";
    c2->SaveAs(fileName2.c_str());
    
    // Create third canvas for summary of unfolded spectra
    TCanvas* c3 = new TCanvas("c3", "c: hist", 1000, 900);
    c3->cd();
    TGaxis::SetMaxDigits(3);
    
    // Set pad and histogram arrangement
    TPad* myPad3 = new TPad("myPad3", "The pad", 0, 0, 1, 1);
    myPad3->SetLeftMargin(0.15);
    myPad3->SetTopMargin(0.06);
    myPad3->SetRightMargin(0.04);
    myPad3->SetBottomMargin(0.15);
    myPad3->SetTicks();
    myPad3->Draw();
    myPad3->cd();
    
    // Create blank histogram for styling
    TH1F* myBlankHisto3 = new TH1F("myBlankHisto3", "Blank Histogram", 20, 0, 1);
    myBlankHisto3->GetXaxis()->SetNdivisions(505);
    myBlankHisto3->SetXTitle(xTitle.c_str());
    myBlankHisto3->GetXaxis()->SetTitleSize(0.05);
    myBlankHisto3->GetXaxis()->SetRangeUser(0, 1);
    myBlankHisto3->GetXaxis()->SetNdivisions(405);
    myBlankHisto3->GetYaxis()->SetTitleOffset(1.35);
    myBlankHisto3->GetYaxis()->SetTitleSize(0.055);
    myBlankHisto3->SetLineColor(0);
    
    // Set y-axis properties
    myBlankHisto3->SetYTitle("dN/dz_{T} unfolded");
    if (zBinsArrayTruth.size() == zBinsArrayDet.size()) {
        myBlankHisto3->GetYaxis()->SetRangeUser(0, 0.45); // When drawn normalized
    } else {
        myBlankHisto3->GetYaxis()->SetRangeUser(0, 0.8); // When drawn normalized
    }
    myBlankHisto3->Draw("E");
    
    // Create legend for unfolded spectra
    TLegend* myLegend3 = nullptr;
    if (hArrayUnfold.size() == 4) {
        myLegend3 = new TLegend(0.2, 0.7, 0.4, 0.9);
    } else if (hArrayUnfold.size() == 6) {
        myLegend3 = new TLegend(0.2, 0.6, 0.4, 0.9);
    } else {
        myLegend3 = new TLegend(0.2, 0.7, 0.4, 0.9);
    }
    
    myLegend3->SetTextFont(42);
    myLegend3->SetBorderSize(0);
    myLegend3->SetFillStyle(0);
    myLegend3->SetFillColor(0);
    myLegend3->SetMargin(0.25);
    myLegend3->SetTextSize(0.04);
    myLegend3->AddEntry(hArrayUnfold[0], " #it{p}_{T}^{jet}:", "");
    
    // Draw all unfolded histograms normalized
    for (size_t i = 0; i < hArrayUnfold.size(); i++) {
        hArrayUnfold[i]->SetMarkerSize(1.3 * markerScale);
        hArrayUnfold[i]->SetMarkerStyle(20);
        hArrayUnfold[i]->SetMarkerColor(colorArrayFinalFigs[i % colorArrayFinalFigs.size()]);
        hArrayUnfold[i]->SetLineStyle(1);
        hArrayUnfold[i]->SetLineWidth(2);
        hArrayUnfold[i]->SetLineColor(colorArrayFinalFigs[i % colorArrayFinalFigs.size()]);
        hArrayUnfold[i]->DrawNormalized("same EP");
        
        myLegend3->AddEntry(hArrayUnfold[i], 
            Form("  %d-%d (GeV/%s)", pTBinsArrayTruth[i], pTBinsArrayTruth[i+1], "#it{c}"), "LP");
    }
    
    myLegend3->Draw();
    
    // Save the canvas
    std::string fileName3 = outPath + "UnfoldingEffect2DFinal_" + figureTag + ".png";
    c3->SaveAs(fileName3.c_str());
    
    // Create fourth canvas for jet_pT > 20 GeV comparison
    int firstBin = -1;
    for (size_t i = 0; i < pTBinsArrayTruth.size(); i++) {
        if (pTBinsArrayTruth[i] == 20) {
            firstBin = i + 1;
            break;
        }
    }
    
    if (firstBin > 0) {
        int lastpTBin = pTBinsArrayTruth.size() - 1;
        
        std::cout << "Project from bin " << firstBin << " to bin " << lastpTBin << std::endl;
        TH1D* mainUnfolded20Proj = mainUnfolded->ProjectionY(
            "_pyMainUnfolded_Bin20", firstBin, lastpTBin);
        TH1D* mainUnfolded20S = scaleHistogram(mainUnfolded20Proj);
        
        TCanvas* c4 = new TCanvas("c4", "c: hist", 1000, 900);
        c4->cd();
        TGaxis::SetMaxDigits(3);
        
        // Set pad and histogram arrangement
        TPad* myPad4 = new TPad("myPad4", "The pad", 0, 0, 1, 1);
        myPad4->SetLeftMargin(0.15);
        myPad4->SetTopMargin(0.06);
        myPad4->SetRightMargin(0.04);
        myPad4->SetBottomMargin(0.15);
        myPad4->SetTicks();
        myPad4->Draw();
        myPad4->cd();
        
        // Create blank histogram for styling
        TH1F* myBlankHisto4 = new TH1F("myBlankHisto4", "Blank Histogram", 20, 0, 1);
        myBlankHisto4->GetXaxis()->SetNdivisions(505);
        myBlankHisto4->SetXTitle(xTitle.c_str());
        myBlankHisto4->GetXaxis()->SetTitleSize(0.05);
        myBlankHisto4->GetXaxis()->SetRangeUser(0, 1);
        myBlankHisto4->GetXaxis()->SetNdivisions(405);
        myBlankHisto4->GetYaxis()->SetTitleOffset(1.35);
        myBlankHisto4->GetYaxis()->SetTitleSize(0.055);
        myBlankHisto4->SetLineColor(0);
        
        // Set y-axis properties
        myBlankHisto4->SetYTitle("dN/dz_{T} unfolded");
        myBlankHisto4->GetYaxis()->SetRangeUser(0, 0.7);
        myBlankHisto4->Draw("E");
        
        // Create legend
        TLegend* myLegend4 = new TLegend(0.2, 0.7, 0.4, 0.9);
        myLegend4->SetTextFont(42);
        myLegend4->SetBorderSize(0);
        myLegend4->SetFillStyle(0);
        myLegend4->SetFillColor(0);
        myLegend4->SetMargin(0.25);
        myLegend4->SetTextSize(0.04);
        myLegend4->AddEntry(hArrayUnfold[0], " #it{p}_{T}^{jet}:", "");
        
        // Style and draw the merged spectrum
        mainUnfolded20S->SetMarkerSize(1.3 * markerScale);
        mainUnfolded20S->SetMarkerStyle(20);
        mainUnfolded20S->SetMarkerColor(colorArrayFinalFigs[0]);
        mainUnfolded20S->SetLineStyle(1);
        mainUnfolded20S->SetLineWidth(2);
        mainUnfolded20S->SetLineColor(colorArrayFinalFigs[0]);
        mainUnfolded20S->DrawNormalized("same EP");
        
        myLegend4->AddEntry(mainUnfolded20S, 
            Form("  %d-%d (GeV/%s)", pTBinsArrayTruth[firstBin-1], pTBinsArrayTruth[lastpTBin], "#it{c}"), "LP");
        myLegend4->Draw();
        
        // Save the canvas
        std::string fileName4 = outPath + "UnfoldingEffect_SpecialJPsiComp_" + figureTag + ".png";
        c4->SaveAs(fileName4.c_str());
    }
    
    // Clean up histogram arrays
    for (auto hist : hArrayOrig) {
        delete hist;
    }
    hArrayOrig.clear();
    
    for (auto hist : hArrayUnfold) {
        delete hist;
    }
    hArrayUnfold.clear();
}

std::pair<std::vector<TH1D*>, std::vector<TH1D*>> 
UnfoldSpectraClass::calculateDispersion(const std::vector<TH2D*>& histogramArray2D,
                                      std::vector<TH1D*> hvariationResult,
                                      std::vector<TH1D*> hvariationError,
                                      int nzTBins, int npTBins) {
    
    std::cout << "\n==== calculateDispersion DEBUG START ====" << std::endl;
    std::cout << "Input parameters:"
              << "\n - histogramArray2D size: " << histogramArray2D.size() 
              << "\n - hvariationResult size: " << hvariationResult.size() 
              << "\n - hvariationError size: " << hvariationError.size()
              << "\n - nzTBins: " << nzTBins 
              << "\n - npTBins: " << npTBins << std::endl;
    
    // Track execution time
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Check input validity
    if (histogramArray2D.empty()) {
        std::cerr << "ERROR: Empty histogram array, cannot calculate dispersion!" << std::endl;
        std::cout << "==== calculateDispersion DEBUG END (early exit) ====\n" << std::endl;
        return std::make_pair(hvariationResult, hvariationError);
    }
    
    if (hvariationResult.size() != npTBins || hvariationError.size() != npTBins) {
        std::cerr << "ERROR: Output histogram arrays don't match npTBins!" 
                  << " hvariationResult: " << hvariationResult.size()
                  << " hvariationError: " << hvariationError.size()
                  << " npTBins: " << npTBins << std::endl;
    }
    
    std::cout << "Processing " << npTBins << " pT bins with " << nzTBins << " zT bins each" << std::endl;
    
    for (int pTBin = 1; pTBin <= npTBins; pTBin++) {
        std::cout << "\n--- Processing pT bin " << pTBin << " (" << pTBinsArrayTruth[pTBin-1] 
                  << "-" << pTBinsArrayTruth[pTBin] << " GeV/c) ---" << std::endl;
        
        // Arrays to store histograms with bin content distributions
        std::vector<TH1D*> hArray1D_v1(nzTBins, nullptr);
        std::vector<TH1D*> hArray1D_v2(nzTBins, nullptr);
        
        std::cout << "Created histogram arrays: hArray1D_v1[" << hArray1D_v1.size() << "] and hArray1D_v2[" 
                  << hArray1D_v2.size() << "]" << std::endl;
        
        // First create canvases to display distributions
        std::cout << "Creating canvas for pT bin " << pTBin << std::endl;
        TCanvas* cA = new TCanvas(("c_" + std::to_string(pTBin)).c_str(), 
                                 "c: pT", 1600, 1600);
        cA->Divide(3, 4);  // Divide canvas into pads for each zT bin
        
        // Process each zT bin
        for (int i_zTBin = 1; i_zTBin <= nzTBins; i_zTBin++) {
            std::cout << "\n  Processing zT bin " << i_zTBin << " for pT bin " << pTBin << std::endl;
            double binCenter = -1;
            
            // First iteration - get rough estimate of mean and width
            std::cout << "  First iteration: collecting values from " << histogramArray2D.size() 
                      << " variations" << std::endl;
            
            std::vector<double> binValues;
            double sumValues = 0;
            int countValues = 0;
            
            for (size_t i_var = 0; i_var < histogramArray2D.size(); i_var++) {
                if (i_var % 250 == 0 && i_var > 0) {
                    std::cout << "    Processed " << i_var << "/" << histogramArray2D.size() 
                              << " variations" << std::endl;
                }
                
                if (!histogramArray2D[i_var]) {
                    std::cerr << "    ERROR: Null histogram at index " << i_var << std::endl;
                    continue;
                }
                
                TH1D* proj = histogramArray2D[i_var]->ProjectionY("_pymain", pTBin, pTBin);
                if (!proj) {
                    std::cerr << "    ERROR: Projection failed for variation " << i_var << std::endl;
                    continue;
                }
                
                double entry = proj->GetBinContent(i_zTBin);
                binCenter = proj->GetBinCenter(i_zTBin);
                binValues.push_back(entry);
                sumValues += entry;
                countValues++;
                                
                // Create histogram for first iteration from first variation
                if (i_var == 0) {
                    // Calculate rough statistics for histogram range
                    double mean = 0;
                    if (!binValues.empty()) {
                        mean = binValues[0]; // Start with first value as initial estimate
                    }
                    
                    std::cout << "    Creating initial histogram with rough range around " << mean << std::endl;
                    
                    // Create histogram with rough range
                    hArray1D_v1[i_zTBin-1] = new TH1D(
                        Form("Var_%zu_zTBin%d", i_var, i_zTBin),
                        Form("Value Distribution for zT bin %d (first pass)", i_zTBin),
                        200, mean * 0.5, mean * 1.5); // Rough range
                        
                    std::cout << "    Created histogram: " << hArray1D_v1[i_zTBin-1]->GetName() 
                              << " with range [" << hArray1D_v1[i_zTBin-1]->GetXaxis()->GetXmin() 
                              << ", " << hArray1D_v1[i_zTBin-1]->GetXaxis()->GetXmax() 
                              << "] and " << hArray1D_v1[i_zTBin-1]->GetNbinsX() << " bins" << std::endl;
                }
                
                // Fill histogram
                if (hArray1D_v1[i_zTBin-1]) {
                    hArray1D_v1[i_zTBin-1]->Fill(entry);
                }
            }
            
            // Calculate rough statistics from collected values
            double meanValue = 0;
            double rmsValue = 0;
            
            if (!binValues.empty()) {
                meanValue = sumValues / countValues;
                
                double sumSquaredDiff = 0;
                for (double val : binValues) {
                    sumSquaredDiff += std::pow(val - meanValue, 2);
                }
                rmsValue = std::sqrt(sumSquaredDiff / countValues);
            }
            
            std::cout << "  First iteration complete. Collected " << binValues.size() 
                      << " values with mean = " << meanValue << ", RMS = " << rmsValue << std::endl;
            
            if (!hArray1D_v1[i_zTBin-1]) {
                std::cerr << "  ERROR: First iteration histogram is null for zT bin " << i_zTBin << std::endl;
                continue;
            }
            
            // Get statistics from first histogram
            int entriesV1 = hArray1D_v1[i_zTBin-1]->GetEntries();
            int maxBin = hArray1D_v1[i_zTBin-1]->GetMaximumBin();
            double maxBinCenter = hArray1D_v1[i_zTBin-1]->GetBinCenter(maxBin);
            double maxBinContent = hArray1D_v1[i_zTBin-1]->GetBinContent(maxBin);
            double histMean = hArray1D_v1[i_zTBin-1]->GetMean();
            double histRMS = hArray1D_v1[i_zTBin-1]->GetRMS();
            
            std::cout << "  First histogram statistics:"
                      << "\n    - Entries: " << entriesV1
                      << "\n    - Maximum bin: " << maxBin << " at x = " << maxBinCenter 
                      << " with content " << maxBinContent
                      << "\n    - Mean: " << histMean
                      << "\n    - RMS: " << histRMS << std::endl;
            
            // Second iteration - improved range based on first histogram
            std::cout << "  Second iteration: Creating histogram with refined range" << std::endl;
            
            for (size_t i_var = 0; i_var < histogramArray2D.size(); i_var++) {
                if (i_var % 250 == 0 && i_var > 0) {
                    std::cout << "    Processed " << i_var << "/" << histogramArray2D.size() 
                              << " variations" << std::endl;
                }
                
                TH1D* proj = histogramArray2D[i_var]->ProjectionY("_pymain", pTBin, pTBin);
                double entry = proj->GetBinContent(i_zTBin);
                
                // Create improved histogram using results from first iteration
                if (i_var == 0) {
                    double mean2 = maxBinCenter; // Use mode instead of mean for better centering
                    double error = histRMS;
                    
                    if (error <= 0) {
                        std::cout << "    WARNING: RMS <= 0, using default error value" << std::endl;
                        error = meanValue * 0.1; // Use 10% of mean as default if RMS is too small
                    }
                    
                    // Create histogram with better range, 4 sigma on each side
                    double rangeMin = mean2 - 4 * error;
                    double rangeMax = mean2 + 4 * error;
                    
                    if (rangeMin < 0 && meanValue > 0) {
                        rangeMin = 0; // Avoid negative range for count histograms
                    }
                    
                    std::cout << "    Creating refined histogram with range [" << rangeMin 
                              << ", " << rangeMax << "] centered at " << mean2 
                              << " with width 4*" << error << std::endl;
                    
                    hArray1D_v2[i_zTBin-1] = new TH1D(
                        Form("Var2_%zu_zTBin%d", i_var, i_zTBin),
                        Form("Value Distribution for zT bin %d (bin center = %.3f)", i_zTBin, binCenter),
                        100, rangeMin, rangeMax);
                }
                
                // Fill second histogram
                if (hArray1D_v2[i_zTBin-1]) {
                    hArray1D_v2[i_zTBin-1]->Fill(entry);
                }
            }
            
            if (!hArray1D_v2[i_zTBin-1]) {
                std::cerr << "  ERROR: Second iteration histogram is null for zT bin " << i_zTBin << std::endl;
                continue;
            }
            
            // Get statistics from second histogram before fitting
            int entriesV2 = hArray1D_v2[i_zTBin-1]->GetEntries();
            double histMeanV2 = hArray1D_v2[i_zTBin-1]->GetMean();
            double histRMSV2 = hArray1D_v2[i_zTBin-1]->GetRMS();
            
            std::cout << "  Second histogram statistics before fitting:"
                      << "\n    - Entries: " << entriesV2
                      << "\n    - Mean: " << histMeanV2
                      << "\n    - RMS: " << histRMSV2 << std::endl;
            
            // Fit Gaussian to distribution if we have enough entries
            std::cout << "  Fitting Gaussian function to distribution..." << std::endl;
            
            double fitMean = histMeanV2;
            double fitSigma = histRMSV2;
            int fitStatus = -1;
            
            if (entriesV2 >= 10) { // Only fit if we have enough entries
                TF1* gaussFunc = new TF1("fb", "gaus(0)", 
                                   hArray1D_v2[i_zTBin-1]->GetBinCenter(1),
                                   hArray1D_v2[i_zTBin-1]->GetBinCenter(
                                       hArray1D_v2[i_zTBin-1]->GetNbinsX()));
                                       
                // Initial parameter guesses improve fit convergence
                gaussFunc->SetParameter(0, hArray1D_v2[i_zTBin-1]->GetMaximum()); // Amplitude
                gaussFunc->SetParameter(1, histMeanV2); // Mean
                gaussFunc->SetParameter(2, histRMSV2);  // Sigma
                
                fitStatus = hArray1D_v2[i_zTBin-1]->Fit(gaussFunc, "QR"); // Quiet fit with restricted range
                
                // Get fit results
                fitMean = gaussFunc->GetParameter(1);
                fitSigma = gaussFunc->GetParameter(2);
                double fitChi2 = gaussFunc->GetChisquare();
                int fitNDF = gaussFunc->GetNDF();
                
                std::cout << "  Gaussian fit results:"
                          << "\n    - Status: " << fitStatus
                          << "\n    - Mean: " << fitMean << " ± " << gaussFunc->GetParError(1)
                          << "\n    - Sigma: " << fitSigma << " ± " << gaussFunc->GetParError(2)
                          << "\n    - χ²/NDF: " << (fitNDF > 0 ? fitChi2/fitNDF : -1) << std::endl;
                
            } else {
                std::cout << "  Too few entries (" << entriesV2 << ") for fitting, using histogram statistics" << std::endl;
            }
            
            // Store the fit results in output histograms
            std::cout << "  Setting results in output histograms for pT bin " << pTBin-1 
                      << ", zT bin " << i_zTBin << std::endl;
            
            if (pTBin-1 < hvariationResult.size()) {
                hvariationResult[pTBin-1]->SetBinContent(i_zTBin, fitMean);
                hvariationResult[pTBin-1]->SetBinError(i_zTBin, fitSigma);
                
                hvariationError[pTBin-1]->SetBinContent(i_zTBin, 1);
                if (fitMean > 0) {
                    double relativeError = fitSigma / fitMean;
                    hvariationError[pTBin-1]->SetBinError(i_zTBin, relativeError);
                    std::cout << "  Relative error = " << relativeError << std::endl;
                } else {
                    std::cout << "  WARNING: fitMean <= 0, cannot calculate relative error" << std::endl;
                    hvariationError[pTBin-1]->SetBinError(i_zTBin, 0);
                }
            } else {
                std::cerr << "  ERROR: pTBin index " << pTBin-1 << " out of range for output histograms" << std::endl;
            }
            
            // Draw the distribution on canvas
            cA->cd(i_zTBin);
            
            if (hArray1D_v2[i_zTBin-1]->GetEntries() > 0) {
                hArray1D_v2[i_zTBin-1]->SetTitle(Form("Distribution for pT bin %d, zT bin %d", pTBin, i_zTBin));
                hArray1D_v2[i_zTBin-1]->GetXaxis()->SetTitle("Bin Content");
                hArray1D_v2[i_zTBin-1]->GetYaxis()->SetTitle("Frequency");
                hArray1D_v2[i_zTBin-1]->Draw();
                
                // Add text with fit results
                TPaveText* fitResults = new TPaveText(0.6, 0.7, 0.9, 0.9, "NDC");
                fitResults->SetFillColor(0);
                fitResults->SetTextSize(0.04);
                fitResults->SetBorderSize(1);
                fitResults->AddText(Form("Mean = %.4g", fitMean));
                fitResults->AddText(Form("Sigma = %.4g", fitSigma));
                if (fitMean > 0) {
                    fitResults->AddText(Form("Rel.Err = %.2f%%", 100*fitSigma/fitMean));
                }
                fitResults->AddText(Form("N = %d", entriesV2));
                fitResults->Draw();
                
                std::cout << "  Distribution drawn on canvas pad " << i_zTBin << std::endl;
            } else {
                std::cout << "  WARNING: Empty histogram, not drawing on canvas" << std::endl;
            }
        }
        
        // Save the canvas with distributions
        std::string fileName = outPath + "RMStatTest2D_TestAttempt_pT" + 
                             std::to_string(pTBin) + "_" + figureTag + ".png";
        std::cout << "Saving canvas to file: " << fileName << std::endl;
        cA->SaveAs(fileName.c_str());
        
        // Clean up histograms
        std::cout << "Cleaning up temporary histograms for pT bin " << pTBin << std::endl;
        for (auto hist : hArray1D_v1) {
            if (hist) delete hist;
        }
        
        for (auto hist : hArray1D_v2) {
            if (hist) delete hist;
        }
    }
    
    // Summarize results
    std::cout << "\nDispersion calculation complete. Summary:" << std::endl;
    for (int pTBin = 0; pTBin < npTBins; pTBin++) {
        std::cout << "pT bin " << pTBin+1 << " (" << pTBinsArrayTruth[pTBin] 
                  << "-" << pTBinsArrayTruth[pTBin+1] << " GeV/c):" << std::endl;
        
        if (pTBin < hvariationResult.size()) {
            TH1D* result = hvariationResult[pTBin];
            TH1D* error = hvariationError[pTBin];
            
            std::cout << "  Bin   zT        Mean      Error     Rel.Error" << std::endl;
            for (int i = 1; i <= result->GetNbinsX(); i++) {
                double mean = result->GetBinContent(i);
                double sigma = result->GetBinError(i);
                double relError = (mean > 0) ? sigma/mean : 0;
                double ztVal = result->GetBinCenter(i);
                
                std::cout << "  " << std::setw(3) << i 
                          << "  " << std::fixed << std::setprecision(3) << std::setw(6) << ztVal
                          << "  " << std::scientific << std::setprecision(4) << std::setw(10) << mean
                          << "  " << std::setw(10) << sigma
                          << "  " << std::fixed << std::setprecision(4) << std::setw(10) << relError
                          << std::endl;
            }
        } else {
            std::cerr << "  ERROR: No result histogram for this pT bin" << std::endl;
        }
    }
    
    // Calculate elapsed time
    auto endTime = std::chrono::high_resolution_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
    std::cout << "Total execution time: " << elapsedTime << " ms (" << (elapsedTime/1000.0) << " s)" << std::endl;
    
    std::cout << "==== calculateDispersion DEBUG END ====\n" << std::endl;
    return std::make_pair(hvariationResult, hvariationError);
}

void UnfoldSpectraClass::StatTestRM(RooUnfoldResponse* response, int ptBin) {
    std::cout << "\n==== StatTestRM DEBUG START ====" << std::endl;
    std::cout << "Testing response matrix statistical variation for pT bin " << ptBin << std::endl;
    
    // Get histogram for this pT bin
    TH1D* histo = measuredSpectraArray[ptBin];
    if (!histo) {
        std::cerr << "ERROR: No measured spectrum available for pT bin " << ptBin << std::endl;
        return;
    }
    
    // Get the 2D response matrix
    TH2D* TH2DRM = dynamic_cast<TH2D*>(response->Hresponse());
    if (!TH2DRM) {
        std::cerr << "ERROR: Failed to cast response matrix to TH2D!" << std::endl;
        return;
    }
    
    // Settings for test
    int RegulForStatTest = 3;  // Fixed regularization parameter
    int errorType = RooUnfold::kCovariance;  // Error type for this test
    TRandom3* fRandom = new TRandom3(0);  // Random generator with fixed seed
    std::vector<TH1D*> unfoldTesthistograms;  // Store unfolded distributions
    
    const int nVariations = 1000;  // Number of random variations to perform
    std::cout << "Will perform " << nVariations << " random variations of response matrix" << std::endl;
    
    // Perform N toy trials with varied response matrix
    int successfulVariations = 0;
    for (int i = 1; i <= nVariations; i++) {
        if (i % 100 == 0) {
            std::cout << "Processing variation " << i << "/" << nVariations << std::endl;
        }
        
        // Create smeared version of the response matrix
        TH2D* RMVariation = SmearPoints(TH2DRM, fRandom, 1);
        
        // Create projections needed for RooUnfoldResponse
        TH1D* hGenLevelMatchedUncutPerBinVar = RMVariation->ProjectionY(
            ("hGenLevelMatchedUncutPerBin_Var" + std::to_string(i)).c_str());
            
        // Create RooUnfoldResponse object from the varied matrix
        RooUnfoldResponse* RooUnfoldRMVar = new RooUnfoldResponse(
            0, hGenLevelMatchedUncutPerBinVar, RMVariation,
            ("hResponseMatrixMain" + unfoldLabel + "_V" + std::to_string(i)).c_str(),
            ("hResponseMatrixMain" + unfoldLabel + "_V" + std::to_string(i)).c_str());
        
        RooUnfoldRMVar->UseOverflow(false);
        
        // Perform the unfolding
        RooUnfoldBayes unfoldVar(RooUnfoldRMVar, histo, RegulForStatTest);
        TH1D* hunf = dynamic_cast<TH1D*>(unfoldVar.Hreco(
            static_cast<RooUnfold::ErrorTreatment>(errorType)));
            
        if (!hunf) {
            std::cerr << "ERROR: Variation " << i << " unfolding failed!" << std::endl;
            delete hGenLevelMatchedUncutPerBinVar;
            delete RMVariation;
            delete RooUnfoldRMVar;
            continue;
        }
        
        unfoldTesthistograms.push_back(hunf);
        successfulVariations++;
        
        // Clean up temporary objects
        delete hGenLevelMatchedUncutPerBinVar;
        delete RMVariation;
        delete RooUnfoldRMVar;
    }
    
    std::cout << "Completed " << successfulVariations << " successful variations out of " 
              << nVariations << " attempted" << std::endl;
    
    // Calculate statistics from the variations
    if (unfoldTesthistograms.empty()) {
        std::cerr << "ERROR: No successful unfolding variations!" << std::endl;
        delete fRandom;
        return;
    }
    
    // Create histograms to store mean and standard deviation
    TH1D* hvariationResult = static_cast<TH1D*>(unfoldTesthistograms[0]->Clone("hvariationResult"));
    TH1D* hvariationError = static_cast<TH1D*>(unfoldTesthistograms[0]->Clone("hvariationError"));
    
    // Reset the histograms
    hvariationResult->Reset("ICESM");
    hvariationError->Reset("ICESM");
    
    // Calculate mean and standard deviation for each bin
    for (int bin = 1; bin <= unfoldTesthistograms[0]->GetNbinsX(); bin++) {
        std::vector<double> binValues;
        
        // Collect values from all variations
        for (auto hist : unfoldTesthistograms) {
            binValues.push_back(hist->GetBinContent(bin));
        }
        
        // Calculate statistics
        double sum = 0.0;
        for (double val : binValues) {
            sum += val;
        }
        double mean = sum / binValues.size();
        
        double variance = 0.0;
        for (double val : binValues) {
            variance += std::pow(val - mean, 2);
        }
        double stdDev = std::sqrt(variance / binValues.size());
        
        // Store in result histograms
        hvariationResult->SetBinContent(bin, mean);
        hvariationResult->SetBinError(bin, stdDev);
        
        hvariationError->SetBinContent(bin, 1);
        if (mean > 0) {
            hvariationError->SetBinError(bin, stdDev / mean);
        }
    }
    
    // Create canvas for plotting results
    TCanvas* c = new TCanvas("c", "c: pT", 800, 850);
    c->Divide(1, 2);
    
    // Upper pad for the mean values
    c->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0);
    
    hvariationResult->GetYaxis()->SetTitle("counts per bin");
    hvariationResult->GetYaxis()->SetTitleSize(0.06);
    hvariationResult->GetYaxis()->SetLabelFont(43);
    hvariationResult->GetYaxis()->SetLabelSize(20);
    hvariationResult->SetNdivisions(505);
    hvariationResult->SetLineColor(4);
    hvariationResult->SetLineWidth(3);
    hvariationResult->SetLineStyle(1);
    hvariationResult->SetMarkerStyle(21);
    hvariationResult->Draw("hist E");
    
    // Lower pad for relative errors
    c->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0.15);
    
    hvariationError->GetYaxis()->SetTitle("Relative Uncertainty");
    hvariationError->GetXaxis()->SetTitle("z_{T}");
    hvariationError->SetLineColor(4);
    hvariationError->SetLineWidth(3);
    hvariationError->SetLineStyle(1);
    hvariationError->SetMarkerStyle(21);
    hvariationError->Draw("hist E");
    
    // Save canvas
    std::string fileName = outPath + "StatTestRM_ptBin" + std::to_string(ptBin) + 
                          "_nIter" + std::to_string(RegulForStatTest) + "_" + figureTag + ".png";
    c->SaveAs(fileName.c_str());
    
    // Clean up
    delete c;
    delete fRandom;
    delete hvariationResult;
    delete hvariationError;
    
    for (auto hist : unfoldTesthistograms) {
        delete hist;
    }
    
    std::cout << "==== StatTestRM DEBUG END ====\n" << std::endl;
}

void UnfoldSpectraClass::StabilityTest2D(int nIter, bool plotAll) {
    std::cout << "\n==== StabilityTest2D DEBUG START ====" << std::endl;
    std::cout << "Input parameters: nIter = " << nIter << ", plotAll = " << (plotAll ? "true" : "false") << std::endl;
    
    // Track execution time
    auto startTime = std::chrono::high_resolution_clock::now();
    
    if (unfoldedArr2D.empty()) {
        std::cerr << "ERROR: unfoldedArr2D array is empty!" << std::endl;
        return;
    }
    
    if (nIter >= unfoldedArr2D.size()) {
        std::cerr << "ERROR: regParam (" << nIter << ") out of range for unfoldedArr2D (size: " 
                  << unfoldedArr2D.size() << ")" << std::endl;
        return;
    }
    
    // Print bin arrays for verification
    std::cout << "pT bin array (Truth): [";
    for (size_t i = 0; i < pTBinsArrayTruth.size(); ++i) {
        std::cout << pTBinsArrayTruth[i] << (i < pTBinsArrayTruth.size()-1 ? ", " : "");
    }
    std::cout << "] (" << (pTBinsArrayTruth.size()-1) << " bins)" << std::endl;
    
    std::cout << "zT bin array (Truth): [";
    for (size_t i = 0; i < zBinsArrayTruth.size(); ++i) {
        std::cout << zBinsArrayTruth[i] << (i < zBinsArrayTruth.size()-1 ? ", " : "");
    }
    std::cout << "] (" << (zBinsArrayTruth.size()-1) << " bins)" << std::endl;
    
    // Make directory for output plots
    std::string outputDirStability = makeOutDir("Stability");
    std::cout << "Created stability output directory: " << outputDirStability << std::endl;
    
    // Get the original unfolded result for reference
    TH2D* unfolded = unfoldedArr2D[nIter];
    if (!unfolded) {
        std::cerr << "ERROR: Null pointer for unfolded histogram at index " << nIter << std::endl;
        return;
    }
    
    std::cout << "Using unfolded histogram from index " << nIter << ": " 
              << unfolded->GetName() << " (" << unfolded->GetNbinsX() << "x" << unfolded->GetNbinsY() 
              << " bins), integral = " << unfolded->Integral() << std::endl;
    
    // Create random number generator
    TRandom3* random = new TRandom3();
    random->SetSeed(0);
    std::cout << "Created TRandom3 with seed 0" << std::endl;
    
    // Create canvas for stability plots
    int npTBins = pTBinsArrayTruth.size() - 1;
    std::cout << "Creating canvas with " << npTBins << " pT bins, total pads = " << (npTBins * 2) << std::endl;
    
    // Check that we have a valid number of bins before creating the canvas
    if (npTBins <= 0) {
        std::cerr << "ERROR: Invalid number of pT bins: " << npTBins << std::endl;
        return;
    }
    
    TCanvas* cStabilityTest2D = new TCanvas("cStabilityTest2D", "cStabilityTest2D: pT", std::max(800, 400 * npTBins), 800);
    if (!cStabilityTest2D) {
        std::cerr << "ERROR: Failed to create canvas" << std::endl;
        return;
    }
    
    // Divide canvas and verify successful division
    cStabilityTest2D->Divide(npTBins+1, 2);
    int totalPads = (npTBins+1) * 2;
    std::cout << "Canvas divided into " << totalPads << " pads (" << npTBins << " columns x 2 rows)" << std::endl;
    
    // Create legend - will be used in first pad
    TLegend* legend = new TLegend(0.55, 0.6, 0.9, 0.85);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.04);
    legend->AddEntry((TObject*)nullptr, " #it{p}_{T}^{jet}:", "");  
    
    // Create arrays to store mean and std dev of ratios
    std::vector<TH1D*> meanHistos;
    std::vector<TH1D*> stdDevHistos;
    
    // Loop through pT bins
    for (int ptBin = 1; ptBin <= npTBins; ptBin++) {
        std::cout << "\nProcessing pT bin " << ptBin << " (" << pTBinsArrayTruth[ptBin-1] 
                  << "-" << pTBinsArrayTruth[ptBin] << " GeV/c)" << std::endl;
        
        // Get projections for this pT bin
        TH1D* unfoldedProj = unfolded->ProjectionY(
            ("_pyUnfolded_Bin" + std::to_string(ptBin)).c_str(), ptBin, ptBin);
            
        std::cout << "Created projection histogram: " << unfoldedProj->GetName() 
                  << " with " << unfoldedProj->GetNbinsX() << " bins, integral = " 
                  << unfoldedProj->Integral() << ", max = " << unfoldedProj->GetMaximum() << std::endl;
        
        // Set up histogram to store mean and std dev for each bin
        TH1D* meanHisto = static_cast<TH1D*>(unfoldedProj->Clone(
            ("meanRatios_ptBin" + std::to_string(ptBin)).c_str()));
        meanHisto->Reset("ICESM");
        std::cout << "Created mean histogram: " << meanHisto->GetName() 
                  << " with " << meanHisto->GetNbinsX() << " bins" << std::endl;
        TH1D* stdDevHisto = static_cast<TH1D*>(unfoldedProj->Clone(
            ("stdDevRatios_ptBin" + std::to_string(ptBin)).c_str()));
        stdDevHisto->Reset("ICESM");
        std::cout << "Created std dev histogram: " << stdDevHisto->GetName() 
                  << " with " << stdDevHisto->GetNbinsX() << " bins" << std::endl;
        // Fill arrays with mean and std dev histos
        meanHistos.push_back(meanHisto);
        stdDevHistos.push_back(stdDevHisto);
        
        // Array to collect ratio histograms for this pT bin
        std::vector<TH1D*> ratioHistos;
        
        // Create upper pad for original histogram
        if (ptBin > npTBins) {
            std::cerr << "ERROR: Attempting to access invalid pad " << ptBin << std::endl;
            continue;
        }
        std::cout << "Processing upper pad " << ptBin << " for pT bin " << ptBin << std::endl;

        cStabilityTest2D->cd(ptBin);

        unfoldedProj->SetLineColor(kBlack);
        unfoldedProj->SetLineWidth(2);
        unfoldedProj->SetMarkerStyle(20);
        unfoldedProj->GetYaxis()->SetTitle("counts");
        unfoldedProj->SetTitle(("p_{T} bin: " + std::to_string(pTBinsArrayTruth[ptBin-1]) + 
                              "-" + std::to_string(pTBinsArrayTruth[ptBin]) + " GeV/#it{c}").c_str());
        unfoldedProj->Draw("E");
        std::cout << "Drawn unfolded projection for pT bin " << ptBin 
                  << ", integral = " << unfoldedProj->Integral() 
                  << ", max = " << unfoldedProj->GetMaximum() << std::endl;
        // Arrays for storing test statistics
        std::vector<double> means(unfoldedProj->GetNbinsX(), 0.0);
        std::vector<double> stdDevs(unfoldedProj->GetNbinsX(), 0.0);
        std::vector<std::vector<double>> ratioValues(unfoldedProj->GetNbinsX());
        
        // Perform N toy trials - fluctuate the measured data
        const int nToys = 20;
        std::cout << "Performing " << nToys << " toy trials for stability test..." << std::endl;
        
        int successfulToys = 0;
        for (int toy = 0; toy < nToys; toy++) {
            std::cout << "  Processing toy " << toy << "..." << std::endl;
            
            // Create smeared data by fluctuating each bin of the measured data
            TH2D* smearedData = static_cast<TH2D*>(measuredSpectra2D->Clone(
                Form("smearedData_ptBin%d_toy%d", ptBin, toy)));
                
            int fluctuatedBins = 0, zeroBins = 0, negBins = 0;
            double originalSum = 0, smearedSum = 0;
            
            for (int binx = 1; binx <= smearedData->GetNbinsX(); binx++) {
                for (int biny = 1; biny <= smearedData->GetNbinsY(); biny++) {
                    double content = smearedData->GetBinContent(binx, biny);
                    double error = smearedData->GetBinError(binx, biny);
                    originalSum += content;
                    
                    double newContent = random->Gaus(content, error);
                    if (newContent < 0) {
                        newContent = 0;
                        zeroBins++;
                    }
                    
                    smearedData->SetBinContent(binx, biny, newContent);
                    smearedData->SetBinError(binx, biny, std::sqrt(newContent));
                    smearedSum += newContent;
                    fluctuatedBins++;
                }
            }
            
            std::cout << "    Smeared data: " << fluctuatedBins << " bins fluctuated, " 
                      << zeroBins << " bins set to zero, original sum = " << originalSum 
                      << ", smeared sum = " << smearedSum 
                      << " (diff: " << (smearedSum - originalSum) << ")" << std::endl;
            
            // Create the response matrix and perform unfolding
            std::cout << "    Creating response matrix..." << std::endl;
            RooUnfoldResponse* toyRM = PrepareResponseMatrix3D(0, nullptr);
            
            if (!toyRM) {
                std::cerr << "    ERROR: Failed to create response matrix for toy " << toy << std::endl;
                continue;
            }
            
            std::cout << "    Performing unfolding with " << nIter << " iterations..." << std::endl;
            RooUnfoldBayes unfoldBayes(toyRM, smearedData, nIter);
            TH2D* toyUnfolded = dynamic_cast<TH2D*>(unfoldBayes.Hreco(
                static_cast<RooUnfold::ErrorTreatment>(errorType)));
                
            if (!toyUnfolded) {
                std::cerr << "    ERROR: Toy unfolding failed for ptBin=" << ptBin << ", toy=" << toy << std::endl;
                continue;
            }
            
            std::cout << "    Unfolding successful. Toy unfolded integral = " << toyUnfolded->Integral() << std::endl;
            successfulToys++;
            
            // Project to 1D for this pT bin
            TH1D* toyUnfoldedProj = toyUnfolded->ProjectionY(
                Form("_pyToyUnfolded_Bin%d_toy%d", ptBin, toy), ptBin, ptBin);
            
            std::cout << "    Projection integral = " << toyUnfoldedProj->Integral() 
                      << ", orig. projection integral = " << unfoldedProj->Integral() << std::endl;
                
            // Calculate ratio to original unfolded
            TH1D* ratioHisto = static_cast<TH1D*>(toyUnfoldedProj->Clone(
                Form("ratioHisto_ptBin%d_toy%d", ptBin, toy)));
            ratioHisto->Divide(unfoldedProj);
            
            // Calculate min, max, mean of ratio
            double ratioMin = 999, ratioMax = -999, ratioSum = 0;
            int ratioBins = 0;
            
            for (int bin = 1; bin <= ratioHisto->GetNbinsX(); bin++) {
                double val = ratioHisto->GetBinContent(bin);
                if (val > 0) {
                    ratioMin = std::min(ratioMin, val);
                    ratioMax = std::max(ratioMax, val);
                    ratioSum += val;
                    ratioBins++;
                }
            }
            
            std::cout << "    Ratio stats: min = " << ratioMin << ", max = " << ratioMax 
                      << ", mean = " << (ratioBins > 0 ? ratioSum/ratioBins : 0) 
                      << " (" << ratioBins << " bins)" << std::endl;
            
            // Set style properties for ratio histo
            int colorIndex = toy % 9 + 1; // Avoid black (0)
            ratioHisto->SetLineColor(colorIndex);
            ratioHisto->SetMarkerColor(colorIndex);
            
            // Draw toy unfolded result on upper pad if requested
            if (plotAll || toy < 5) {
                toyUnfoldedProj->SetLineColor(colorIndex);
                toyUnfoldedProj->SetLineStyle(2);
                toyUnfoldedProj->Draw("hist same");
                
                if (ptBin == 1 && toy == 0) {
                    legend->AddEntry(unfoldedProj, "Original unfolded", "lp");
                    legend->AddEntry(toyUnfoldedProj, "Toy variations", "l");
                }
            }
            
            // Add ratio histogram to collection
            ratioHistos.push_back(ratioHisto);
            
            // Store ratio values for statistics
            for (int bin = 1; bin <= ratioHisto->GetNbinsX(); bin++) {
                ratioValues[bin-1].push_back(ratioHisto->GetBinContent(bin));
            }
            
            // Clean up temporary objects
        }
        
        std::cout << "Completed " << successfulToys << " successful toys out of " << nToys 
                  << " attempted (" << (100.0*successfulToys/nToys) << "%)" << std::endl;
        
        // Draw legend in first pad
        if (ptBin == 1) {
            legend->Draw("same");
        }
        
        // Calculate mean and std dev for each bin
        std::cout << "Calculating statistics for pT bin " << ptBin << "..." << std::endl;
        for (int bin = 1; bin <= unfoldedProj->GetNbinsX(); bin++) {
            // Skip bins with no values
            if (ratioValues[bin-1].empty()) {
                std::cout << "  Bin " << bin << " (zT = " << unfoldedProj->GetBinCenter(bin) 
                          << "): No valid ratio values" << std::endl;
                continue;
            }
            
            // Calculate mean of ratios
            double sum = 0.0;
            for (double val : ratioValues[bin-1]) {
                sum += val;
            }
            double mean = sum / ratioValues[bin-1].size();
            means[bin-1] = mean;
            
            // Calculate std dev of ratios
            double variance = 0.0;
            for (double val : ratioValues[bin-1]) {
                variance += std::pow(val - mean, 2);
            }
            double stdDev = std::sqrt(variance / ratioValues[bin-1].size());
            stdDevs[bin-1] = stdDev;
            
            // Report statistics
            std::cout << "  Bin " << bin << " (zT = " << unfoldedProj->GetBinCenter(bin) 
                      << "): mean = " << mean << ", stdDev = " << stdDev 
                      << " from " << ratioValues[bin-1].size() << " values" << std::endl;
            
            // Set values in histograms
            meanHisto->SetBinContent(bin, mean);
            meanHisto->SetBinError(bin, stdDev / std::sqrt(ratioValues[bin-1].size())); // Standard error
            
            stdDevHisto->SetBinContent(bin, stdDev);
        }
                std::cout << __LINE__ << std::endl;

        // Before accessing lower pad, verify it exists
        int lowerPadIndex = npTBins + ptBin;
        if (lowerPadIndex > totalPads) {
            std::cerr << "ERROR: Attempting to access invalid pad " << lowerPadIndex 
                      << " (total pads: " << totalPads << ")" << std::endl;
            return;
        }
        
        std::cout << "Processing lower pad " << lowerPadIndex << " for pT bin " << ptBin << std::endl;
        TVirtualPad* lowerPad = cStabilityTest2D->cd(lowerPadIndex);
        
        // Verify pad is valid before proceeding
        if (!lowerPad) {
            std::cerr << "ERROR: Failed to get pad " << lowerPadIndex << std::endl;
            return;
        }

        // Set up blank histogram for ratio plot
        TH1D* ratioFrame = static_cast<TH1D*>(unfoldedProj->Clone("ratioFrame"));
        ratioFrame->Reset("ICESM");
        ratioFrame->SetTitle("");
        ratioFrame->GetYaxis()->SetTitle("Ratio to nominal");
        ratioFrame->GetYaxis()->SetRangeUser(0.8, 1.2);
        ratioFrame->Draw("hist");
                std::cout << __LINE__ << std::endl;

        // Draw reference lines
        TLine* line1 = new TLine(ratioFrame->GetXaxis()->GetXmin(), 1.0, 
                                ratioFrame->GetXaxis()->GetXmax(), 1.0);
        line1->SetLineColor(kBlack);
        line1->SetLineWidth(2);
        line1->Draw("same");
                std::cout << __LINE__ << std::endl;

        TLine* line2 = new TLine(ratioFrame->GetXaxis()->GetXmin(), 1.1, 
                                ratioFrame->GetXaxis()->GetXmax(), 1.1);
        line2->SetLineColor(kBlack);
        line2->SetLineStyle(2);
        line2->Draw("same");
        
        TLine* line3 = new TLine(ratioFrame->GetXaxis()->GetXmin(), 0.9, 
                                ratioFrame->GetXaxis()->GetXmax(), 0.9);
        line3->SetLineColor(kBlack);
        line3->SetLineStyle(2);
        line3->Draw("same");
        
        // Draw all ratio histograms
        std::cout << "Drawing " << (plotAll ? "all" : "up to 5") << " ratio histograms..." << std::endl;
        int drawnHistos = 0;
        for (size_t i = 0; i < ratioHistos.size(); i++) {
            if (plotAll || i < 5) {
                ratioHistos[i]->Draw("hist same");
                drawnHistos++;
            }
        }
        std::cout << "Drew " << drawnHistos << " ratio histograms out of " << ratioHistos.size() << std::endl;
        
        // Draw mean and stdDev histograms on top
        meanHisto->SetLineColor(kBlack);
        meanHisto->SetLineWidth(3);
        meanHisto->SetMarkerStyle(20);
        meanHisto->SetMarkerColor(kBlack);
        meanHisto->Draw("E same");
        
        // Clean up ratio histograms
        ratioHistos.clear();
        std::cout << "Cleaned up " << ratioHistos.capacity() << " ratio histograms" << std::endl;        
    }
    
    // Save canvas
    std::string fileName = outputDirStability + "/StabilityTest2D_nIter" + std::to_string(nIter) + ".png";
    std::cout << "Saving plot to: " << fileName << std::endl;
    cStabilityTest2D->SaveAs(fileName.c_str());
        
    // Create canvas for mean and std dev summary
    TCanvas* c2 = new TCanvas("c2", "StabilityTest Summary", 1200, 800);
    c2->Divide(2, 1);
    std::cout << "Created summary canvas with 2x1 pads" << std::endl;
    
    // Plot mean values
    c2->cd(1);
    TLegend* leg2 = new TLegend(0.55, 0.6, 0.9, 0.85);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.04);
    
    std::cout << "Drawing mean histograms..." << std::endl;
    for (int i = 0; i < npTBins; i++) {
        meanHistos[i]->SetLineColor(i+1);
        meanHistos[i]->SetMarkerColor(i+1);
        meanHistos[i]->SetMarkerStyle(20);
        
        if (i == 0) {
            meanHistos[i]->SetTitle("Mean of toy ratios");
            meanHistos[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
            meanHistos[i]->Draw("E");
            std::cout << "  First mean histogram drawn with range [0.8, 1.2]" << std::endl;
        } else {
            meanHistos[i]->Draw("E same");
        }
        
        leg2->AddEntry(meanHistos[i], 
            ("p_{T} " + std::to_string(pTBinsArrayTruth[i]) + "-" + 
             std::to_string(pTBinsArrayTruth[i+1]) + " GeV/#it{c}").c_str(), "lp");
    }
    
    // Draw reference lines
    TLine* line1m = new TLine(meanHistos[0]->GetXaxis()->GetXmin(), 1.0, 
                            meanHistos[0]->GetXaxis()->GetXmax(), 1.0);
    line1m->SetLineColor(kBlack);
    line1m->SetLineWidth(2);
    line1m->Draw("same");
    
    TLine* line2m = new TLine(meanHistos[0]->GetXaxis()->GetXmin(), 1.1, 
                            meanHistos[0]->GetXaxis()->GetXmax(), 1.1);
    line2m->SetLineColor(kBlack);
    line2m->SetLineStyle(2);
    line2m->Draw("same");
    
    TLine* line3m = new TLine(meanHistos[0]->GetXaxis()->GetXmin(), 0.9, 
                            meanHistos[0]->GetXaxis()->GetXmax(), 0.9);
    line3m->SetLineColor(kBlack);
    line3m->SetLineStyle(2);
    line3m->Draw("same");
    
    leg2->Draw("same");
    
    // Plot standard deviation values
    c2->cd(2);
    TLegend* leg3 = new TLegend(0.55, 0.6, 0.9, 0.85);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->SetTextSize(0.04);
    
    std::cout << "Drawing stdDev histograms..." << std::endl;
    for (int i = 0; i < npTBins; i++) {
        stdDevHistos[i]->SetLineColor(i+1);
        stdDevHistos[i]->SetMarkerColor(i+1);
        stdDevHistos[i]->SetMarkerStyle(20);
        
        if (i == 0) {
            stdDevHistos[i]->SetTitle("Standard deviation of toy ratios");
            stdDevHistos[i]->GetYaxis()->SetTitle("Standard Deviation");
            stdDevHistos[i]->GetYaxis()->SetRangeUser(0, 0.2);
            stdDevHistos[i]->Draw("E");
            std::cout << "  First stdDev histogram drawn with range [0, 0.2]" << std::endl;
        } else {
            stdDevHistos[i]->Draw("E same");
        }
        
        leg3->AddEntry(stdDevHistos[i], 
            ("p_{T} " + std::to_string(pTBinsArrayTruth[i]) + "-" + 
             std::to_string(pTBinsArrayTruth[i+1]) + " GeV/#it{c}").c_str(), "lp");
    }
    
    leg3->Draw("same");
    
    // Save summary canvas
    std::string fileNameSummary = outputDirStability + "/StabilityTest2D_Summary_nIter" + 
                                 std::to_string(nIter) + ".png";
    std::cout << "Saving summary plot to: " << fileNameSummary << std::endl;
    c2->SaveAs(fileNameSummary.c_str());
    
    // Calculate total execution time
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
    std::cout << "Total execution time: " << duration << " ms (" << (duration/1000.0) << " s)" << std::endl;
    
    // Clean up
    // delete c;
    // delete c2;
    // delete random;
    
    std::cout << "Cleaning up " << meanHistos.size() << " mean histograms and " 
              << stdDevHistos.size() << " stdDev histograms..." << std::endl;
    
    for (auto hist : meanHistos) {
        delete hist;
    }
    meanHistos.clear();
    
    for (auto hist : stdDevHistos) {
        delete hist;
    }
    stdDevHistos.clear();
    
    std::cout << "==== StabilityTest2D DEBUG END ====\n" << std::endl;
}

void UnfoldSpectraClass::StabilityTest(int ptBin, int nIter, bool useRatioNorm) {
    std::cout << "\n==== StabilityTest DEBUG START ====" << std::endl;
    std::cout << "Testing stability for ptBin=" << ptBin << ", nIter=" << nIter 
              << ", useRatioNorm=" << (useRatioNorm ? "true" : "false") << std::endl;
    
    if (ptBin >= unfoldedArrPerBin.size() || nIter >= unfoldedArrPerBin[ptBin].size()) {
        std::cerr << "Error: Invalid ptBin or nIter!" << std::endl;
        return;
    }
    
    // Get the original unfolded result
    TH1D* unfoldedResult = unfoldedArrPerBin[ptBin][nIter];
    if (!unfoldedResult) {
        std::cerr << "Error: Null pointer for unfolded histogram!" << std::endl;
        return;
    }
    
    // Create random number generator
    TRandom3* random = new TRandom3();
    random->SetSeed(0);
    
    // Make directory for output plots
    std::string outputDirStability = makeOutDir("Stability");
    
    // Create canvas
    TCanvas* c = new TCanvas("c", "Stability Test", 800, 800);
    c->Divide(1, 2);
    
    // Upper pad for original and varied histograms
    c->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    
    // Style the original histogram
    unfoldedResult->SetLineColor(kBlack);
    unfoldedResult->SetLineWidth(2);
    unfoldedResult->SetMarkerStyle(20);
    unfoldedResult->GetYaxis()->SetTitle("counts per bin");
    unfoldedResult->SetTitle(("Stability Test, p_{T} bin: " + std::to_string(pTBinsArrayTruth[ptBin]) + 
                            "-" + std::to_string(pTBinsArrayTruth[ptBin+1]) + " GeV/#it{c}").c_str());
    unfoldedResult->Draw("E");
    
    // Create legend
    TLegend* legend = new TLegend(0.55, 0.6, 0.9, 0.85);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.04);
    legend->AddEntry(unfoldedResult, "Original unfolded", "lp");
    
    // Arrays for storing mean and std dev
    std::vector<double> means(unfoldedResult->GetNbinsX(), 0.0);
    std::vector<double> stdDevs(unfoldedResult->GetNbinsX(), 0.0);
    std::vector<std::vector<double>> ratioValues(unfoldedResult->GetNbinsX());
    std::vector<TH1D*> ratioHistos;
    
    // Perform N toy trials
    const int nToys = 20;
    const bool plotAll = false;
    
    for (int toy = 0; toy < nToys; toy++) {
        // Create smeared data by fluctuating each bin of the measured data
        TH1D* smearedData = static_cast<TH1D*>(measuredSpectraArray[ptBin]->Clone(
            Form("smearedData_ptBin%d_toy%d", ptBin, toy)));
        
        for (int bin = 1; bin <= smearedData->GetNbinsX(); bin++) {
            double content = smearedData->GetBinContent(bin);
            double error = smearedData->GetBinError(bin);
            
            double newContent = random->Gaus(content, error);
            if (newContent < 0) newContent = 0;
            
            smearedData->SetBinContent(bin, newContent);
            smearedData->SetBinError(bin, std::sqrt(newContent));
        }
        
        // Create response matrix and unfold
        RooUnfoldResponse* toyRM = PrepareResponseMatrix2D(0);
        RooUnfoldBayes unfoldBayes(toyRM, smearedData, nIter);
        TH1D* toyUnfolded = dynamic_cast<TH1D*>(unfoldBayes.Hreco(
            static_cast<RooUnfold::ErrorTreatment>(errorType)));
            
        if (!toyUnfolded) {
            std::cerr << "ERROR: Toy unfolding failed for toy=" << toy << std::endl;
            delete smearedData;
            delete toyRM;
            continue;
        }
        
        // Draw toy result if needed
        int colorIndex = toy % 9 + 1; // Avoid black (0)
        toyUnfolded->SetLineColor(colorIndex);
        toyUnfolded->SetLineStyle(2);
        
        if (plotAll || toy < 5) {
            toyUnfolded->Draw("hist same");
            
            if (toy == 0) {
                legend->AddEntry(unfoldedResult, "Original unfolded", "lp");
                legend->AddEntry(toyUnfolded, "Toy variations", "l");
            }
        }
        
        // Calculate ratio to original
        TH1D* ratioHisto = static_cast<TH1D*>(toyUnfolded->Clone(
            Form("ratioHisto_toy%d", toy)));
        
        if (useRatioNorm) {
            // Normalize before dividing to preserve shape comparison
            double toyIntegral = toyUnfolded->Integral();
            double origIntegral = unfoldedResult->Integral();
            
            if (toyIntegral > 0 && origIntegral > 0) {
                toyUnfolded->Scale(origIntegral / toyIntegral);
            }
        }
        
        ratioHisto->Divide(unfoldedResult);
        ratioHisto->SetLineColor(colorIndex);
        ratioHistos.push_back(ratioHisto);
        
        // Store ratio values for statistics
        for (int bin = 1; bin <= ratioHisto->GetNbinsX(); bin++) {
            ratioValues[bin-1].push_back(ratioHisto->GetBinContent(bin));
        }
        
        // Clean up temporary objects
        delete smearedData;
        delete toyRM;
        delete toyUnfolded;
    }
    
    legend->Draw("same");
    
    // Lower pad for ratios
    c->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    
    // Create blank histogram for ratios
    TH1D* ratioFrame = static_cast<TH1D*>(unfoldedResult->Clone("ratioFrame"));
    ratioFrame->Reset("ICESM");
    ratioFrame->SetTitle("");
    ratioFrame->GetYaxis()->SetTitle("Ratio to nominal");
    ratioFrame->GetYaxis()->SetRangeUser(0.8, 1.2);
    ratioFrame->Draw("hist");
    
    // Draw reference lines
    TLine* line1 = new TLine(ratioFrame->GetXaxis()->GetXmin(), 1.0, 
                            ratioFrame->GetXaxis()->GetXmax(), 1.0);
    line1->SetLineColor(kBlack);
    line1->SetLineWidth(2);
    line1->Draw("same");
    
    TLine* line2 = new TLine(ratioFrame->GetXaxis()->GetXmin(), 1.1, 
                            ratioFrame->GetXaxis()->GetXmax(), 1.1);
    line2->SetLineColor(kBlack);
    line2->SetLineStyle(2);
    line2->Draw("same");
    
    TLine* line3 = new TLine(ratioFrame->GetXaxis()->GetXmin(), 0.9, 
                            ratioFrame->GetXaxis()->GetXmax(), 0.9);
    line3->SetLineColor(kBlack);
    line3->SetLineStyle(2);
    line3->Draw("same");
    
    // Calculate mean and std dev for each bin
    TH1D* meanHisto = static_cast<TH1D*>(unfoldedResult->Clone("meanRatios"));
    meanHisto->Reset("ICESM");
    
    TH1D* stdDevHisto = static_cast<TH1D*>(unfoldedResult->Clone("stdDevRatios"));
    stdDevHisto->Reset("ICESM");
    
    for (int bin = 1; bin <= unfoldedResult->GetNbinsX(); bin++) {
        // Skip bins with no values
        if (ratioValues[bin-1].empty()) continue;
        
        // Calculate mean of ratios
        double sum = 0.0;
        for (double val : ratioValues[bin-1]) {
            sum += val;
        }
        double mean = sum / ratioValues[bin-1].size();
        means[bin-1] = mean;
        
        // Calculate std dev of ratios
        double variance = 0.0;
        for (double val : ratioValues[bin-1]) {
            variance += std::pow(val - mean, 2);
        }
        double stdDev = std::sqrt(variance / ratioValues[bin-1].size());
        stdDevs[bin-1] = stdDev;
        
        // Set values in histograms
        meanHisto->SetBinContent(bin, mean);
        meanHisto->SetBinError(bin, stdDev / std::sqrt(ratioValues[bin-1].size())); // Standard error
        
        stdDevHisto->SetBinContent(bin, stdDev);
    }
    
    // Draw all ratio histograms
    for (size_t i = 0; i < ratioHistos.size(); i++) {
        if (plotAll || i < 5) {
            ratioHistos[i]->Draw("hist same");
        }
    }
    
    // Draw mean histogram
    meanHisto->SetLineColor(kBlack);
    meanHisto->SetLineWidth(3);
    meanHisto->SetMarkerStyle(20);
    meanHisto->SetMarkerColor(kBlack);
    meanHisto->Draw("E same");
    
    // Save canvas
    std::string fileName = outputDirStability + "/StabilityTest_ptBin" + std::to_string(ptBin) + 
                         "_nIter" + std::to_string(nIter);
    if (useRatioNorm) {
        fileName += "_normalized";
    }
    fileName += ".png";
    c->SaveAs(fileName.c_str());
    
    // Create second canvas for stdDev histogram
    TCanvas* c2 = new TCanvas("c2", "Standard Deviation", 800, 600);
    c2->SetLeftMargin(0.15);
    c2->SetRightMargin(0.05);
    
    stdDevHisto->SetTitle("Standard Deviation of Toy Ratios");
    stdDevHisto->GetYaxis()->SetTitle("Standard Deviation");
    stdDevHisto->GetYaxis()->SetRangeUser(0, 0.2);
    stdDevHisto->SetLineColor(kBlue);
    stdDevHisto->SetMarkerColor(kBlue);
    stdDevHisto->SetMarkerStyle(20);
    stdDevHisto->Draw("E");
    
    // Save stdDev canvas
    std::string fileName2 = outputDirStability + "/StabilityTest_StdDev_ptBin" + std::to_string(ptBin) + 
                          "_nIter" + std::to_string(nIter);
    if (useRatioNorm) {
        fileName2 += "_normalized";
    }
    fileName2 += ".png";
    c2->SaveAs(fileName2.c_str());
    
    // Clean up
    delete c;
    delete c2;
    delete random;
    
    // FIX: Remove references to non-existent vectors and clean up individual histograms
    std::cout << "Cleaning up histograms..." << std::endl;
    
    // Delete individual histograms
    delete meanHisto;
    delete stdDevHisto;
    delete ratioFrame;
    
    // Delete ratio histograms
    for (auto hist : ratioHistos) {
        delete hist;
    }
    ratioHistos.clear();
    
    std::cout << "==== StabilityTest DEBUG END ====\n" << std::endl;
}

void UnfoldSpectraClass::getKinEfficiency(int bin) {
    // Define colors for plotting different pT bins
    std::vector<int> colorArrayFinalFigs = {kAzure, kAzure-4, kCyan-6, kGreen-3, kTeal-6, kGreen+4, 1, 1};
    std::vector<TH1D*> kinEffArr;
    
    for (size_t bin = 0; bin < pTBinsArrayTruth.size() - 1; bin++) {
        // Get full response without cuts
        TH2D* Full = getResponseMatrix(0, "zTDet", "zTPart", true, bin, false);
        
        // Get response with cuts
        TH2D* Cut = getResponseMatrix(0, "zTDet", "zTPart", true, bin, true);
        
        // Project to 1D and calculate efficiency
        TH1D* projFull = Full->ProjectionY(("_pyFull_" + std::to_string(bin)).c_str(), 1, Full->GetNbinsX());
        TH1D* projCut = Cut->ProjectionY(("_pyCut_" + std::to_string(bin)).c_str(), 1, Cut->GetNbinsX());
        projCut->SetName(("ProjectionBin_" + std::to_string(bin)).c_str());
        projCut->Divide(projFull);
        kinEffArr.push_back(projCut);
        
        // Clean up temporary histograms
        delete Full;
        delete Cut;
        delete projFull;
    }
    
    // Define plot title based on particle type
    std::string xTitle = "p_{T}^{D^{0}}/p_{T}^{jet}";
    
    // Create canvas for plotting efficiencies
    TCanvas* c3 = new TCanvas("c2", "c: hist", 1000, 900);
    c3->cd();
    TGaxis::SetMaxDigits(3);
    
    // Set pad and histogram arrangement
    TPad* myPad3 = new TPad("myPad", "The pad", 0, 0, 1, 1);
    myPad3->SetLeftMargin(0.15);
    myPad3->SetTopMargin(0.06);
    myPad3->SetRightMargin(0.04);
    myPad3->SetBottomMargin(0.15);
    myPad3->SetTicks();
    myPad3->Draw();
    myPad3->cd();
    
    // Create blank histogram for styling
    TH1F* myBlankHisto3 = new TH1F("myBlankHisto3", "Blank Histogram", 20, 0, 1);
    myBlankHisto3->GetXaxis()->SetNdivisions(505);
    myBlankHisto3->SetXTitle(xTitle.c_str());
    myBlankHisto3->GetXaxis()->SetTitleSize(0.05);
    myBlankHisto3->GetXaxis()->SetRangeUser(0, 1);
    myBlankHisto3->GetXaxis()->SetNdivisions(405);
    myBlankHisto3->GetYaxis()->SetTitleOffset(1.35);
    myBlankHisto3->GetYaxis()->SetTitleSize(0.055);
    myBlankHisto3->SetLineColor(0);
    myBlankHisto3->SetYTitle("kinematic efficiency");
    myBlankHisto3->GetYaxis()->SetRangeUser(0.7, 1.3);
    myBlankHisto3->Draw("P");
    
    // Create legend
    TLegend* myLegend3 = new TLegend(0.2, 0.6, 0.4, 0.9);
    myLegend3->SetTextFont(42);
    myLegend3->SetBorderSize(0);
    myLegend3->SetFillStyle(0);
    myLegend3->SetFillColor(0);
    myLegend3->SetMargin(0.25);
    myLegend3->SetTextSize(0.04);
    myLegend3->AddEntry(kinEffArr[0], " #it{p}_{T}^{jet}:", "");
    
    // Draw all efficiency curves
    double markerScale = 1.6;
    for (size_t i = 0; i < kinEffArr.size(); i++) {
        kinEffArr[i]->SetMarkerSize(1.3 * markerScale);
        kinEffArr[i]->SetMarkerStyle(20);
        kinEffArr[i]->SetMarkerColor(colorArrayFinalFigs[i]);
        kinEffArr[i]->SetLineStyle(1);
        kinEffArr[i]->SetLineWidth(2);
        kinEffArr[i]->SetLineColor(colorArrayFinalFigs[i]);
        kinEffArr[i]->Draw("hist P same");
        
        myLegend3->AddEntry(kinEffArr[i], 
            ("  " + std::to_string(pTBinsArrayTruth[i]) + "-" + std::to_string(pTBinsArrayTruth[i+1]) + " (GeV/#it{c})").c_str(), 
            "LP");
    }
    
    myLegend3->Draw();
    
    // Save the canvas
    std::string fileName = outPath + "kinematicEff_" + figureTag + ".png";
    c3->SaveAs(fileName.c_str());
    
    // Clean up
    delete c3;
    for (TH1D* hist : kinEffArr) {
        delete hist;
    }
}



// Implementation of deScaleGraph helper method
void UnfoldSpectraClass::deScaleGraph(TGraphErrors* graph) {
    // This function undoes normalization by bin width
    int nPoints = graph->GetN();
    for (int i = 0; i < nPoints; i++) {
        double x, y;
        graph->GetPoint(i, x, y);
        double ex = graph->GetErrorX(i);
        double ey = graph->GetErrorY(i);
        
        // Scale factor is twice the x error (bin width)
        double scale = ex * 2.0;
        
        // Scale bin content and error
        y *= scale;
        ey *= scale;
        
        // Update the graph
        graph->SetPoint(i, x, y);
        graph->SetPointError(i, ex, ey);
    }
}

// Destructor
UnfoldSpectraClass::~UnfoldSpectraClass() {
    // Clean up allocated memory
    if (fResponse) delete fResponse;
    if (fData) delete fData;
    
    // Clean up histograms
    for (auto hist : unfoldedArr2D) {
        if (hist) delete hist;
    }
    
    for (auto& binArray : unfoldedArrPerBin) {
        for (auto hist : binArray) {
            if (hist) delete hist;
        }
    }
    
    for (auto hist : measuredSpectraArray) {
        if (hist) delete hist;
    }
    
    if (measuredSpectra2D) delete measuredSpectra2D;
    if (hOriginalPrior) delete hOriginalPrior;
    if (hCurrentPrior) delete hCurrentPrior;
    if (externalPrior) delete externalPrior;
}


// Implementation of PrepareResponseMatrix3D method
RooUnfoldResponse* UnfoldSpectraClass::PrepareResponseMatrix3D(int part, TH2D* weightHistogram) {
    
    // Create empty RooUnfoldResponse object
    std::string responseName = "hResponseMatrix3DMain" + unfoldLabel + "_part" + std::to_string(part);
    RooUnfoldResponse* RooUnfoldRM = new RooUnfoldResponse(responseName.c_str(), responseName.c_str());

    // Get 2D histograms for true and detector level
    TH2D* RM2DTrue = getResponseMatrix(part, "jetPtPart", "zTPart", false, 0, false);
    TH2D* RM2DDet = getResponseMatrix(part, "jetPtDet", "zTDet", false, 0, false);
    
    // Setup the RooUnfoldResponse with these histograms (will be overwritten in Fill function)
    RooUnfoldRM->Setup(RM2DDet, RM2DTrue);

    // Define mass cut based on particle type
    std::pair<double, double> MassCut = std::make_pair(1.81, 1.935);
    

    const double jetPtMin = 5.0;

    // Create random number generator
    TRandom* R = new TRandom();
    
    // Get the TTree from the file
    TTree* tTree = static_cast<TTree*>(fResponse->Get("Response"));
    if (!tTree) {
        std::cerr << "Error: Could not find 'Response' TTree in file!" << std::endl;
        delete R;
        return RooUnfoldRM;
    }
    
    // Variables to hold tree data - FIXED: changed to float to match Float_t in branches
    float zTDet, zTPart, jetPtDet, jetPtPart, tagPtDet, tagPtPart;
    float etaDet, nConstDet, nConstPart;
    double weight, b_weight;  // These can stay double as they're not directly from branches
    
    // Set branch addresses with correct types
    tTree->SetBranchAddress("zTDet", &zTDet);
    tTree->SetBranchAddress("zTPart", &zTPart);
    tTree->SetBranchAddress("jetPtDet", &jetPtDet);
    tTree->SetBranchAddress("jetPtPart", &jetPtPart);
    tTree->SetBranchAddress("tagPtDet", &tagPtDet);
    tTree->SetBranchAddress("tagPtPart", &tagPtPart);
    tTree->SetBranchAddress("nConstDet", &nConstDet);
    tTree->SetBranchAddress("nConstPart", &nConstPart);
    tTree->SetBranchAddress("etaDet", &etaDet);
    
    // B-decay weight branch
    if (!isPrompt) {
        tTree->SetBranchAddress("weight", &b_weight);
    } else {
        b_weight = 1.0;
    }
    
    // Loop over all events
    for (int iEvt = 0; iEvt < tTree->GetEntries(); ++iEvt) {
        double evt_weight = 1.0;
        b_weight = 1.0;
        
        tTree->GetEntry(iEvt);
        
        double rndm = R->Uniform(1.0);
        
        // Apply jet pT cut
        if (jetPtDet < jetPtMin) continue;
        
        // Apply response matrix cuts if needed
        bool passCuts = !applyRMCut || 
                        (applyRMCut && nConstDet > 1 && nConstPart > 1 && 
                         etaDet > 2.5 && etaDet < 4.0);
                         
        if (passCuts) {
            // Apply weight if provided
            if (weightHistogram) {
                evt_weight = getEventWeight(weightHistogram, zTPart, jetPtPart);
                weight = b_weight * evt_weight;
            } else {
                weight = b_weight;
            }
            
            // Fill response matrix based on part parameter
            if (part == 0) {
                if (useTagPt) {
                    RooUnfoldRM->Fill(tagPtDet, zTDet, tagPtPart, zTPart, weight);
                } else {
                    RooUnfoldRM->Fill(jetPtDet, zTDet, jetPtPart, zTPart, weight);
                }
            }
            // Fill with half of statistics (for testing purposes)
            else if (part == 1 && rndm > 0.5) {
                if (useTagPt) {
                    RooUnfoldRM->Fill(tagPtDet, zTDet, tagPtPart, zTPart, weight);
                } else {
                    RooUnfoldRM->Fill(jetPtDet, zTDet, jetPtPart, zTPart, weight);
                }
            }
            // Fill with other half of statistics (for testing purposes)
            else if (part == 2 && rndm < 0.5) {
                if (useTagPt) {
                    RooUnfoldRM->Fill(tagPtDet, zTDet, tagPtPart, zTPart, weight);
                } else {
                    RooUnfoldRM->Fill(jetPtDet, zTDet, jetPtPart, zTPart, weight);
                }
            }
        }
    }
    
    // Don't use overflow bins
    RooUnfoldRM->UseOverflow(false);
    
    // Plot histograms for QA
    plotHist(RooUnfoldRM->Htruth(), "h2DRMtrue_RM_" + std::to_string(part), "colz", 
             "p_{T,jet} (GeV/c)", "z_{T}");
    plotHist(RooUnfoldRM->Hmeasured(), "h2DRMdet_RM_" + std::to_string(part), "colz", 
             "p_{T,jet} (GeV/c)", "z_{T}");
    plotHist(RooUnfoldRM->Hresponse(), "h2DRM_RM_" + std::to_string(part), "colz", " ", " ");
    
    // Clean up
    delete R;
    
    return RooUnfoldRM;
}

TH2D* UnfoldSpectraClass::prepareRMWeights(int regParam, int iteration) {
    // This function finds the weights that need to be used while filling the RM
    
    std::cout << "\n==== prepareRMWeights DEBUG ====" << std::endl;
    std::cout << "Input parameters: regParam = " << regParam << ", iteration = " << iteration << std::endl;
    
    TH2D* hweights = nullptr;
    
    if (externalPrior && iteration == 0) {
        std::cout << "Using external prior for weights calculation (first iteration)" << std::endl;
        
        // Check original prior histogram validity
        if (!hOriginalPrior) {
            std::cerr << "ERROR: hOriginalPrior is null!" << std::endl;
            return nullptr;
        }
        
        // Debug histogram properties
        std::cout << "Original prior: " << hOriginalPrior->GetNbinsX() << "x" 
                  << hOriginalPrior->GetNbinsY() << " bins, integral = " 
                  << hOriginalPrior->Integral() << ", max = " << hOriginalPrior->GetMaximum() << std::endl;
        
        std::cout << "External prior: " << externalPrior->GetNbinsX() << "x" 
                  << externalPrior->GetNbinsY() << " bins, integral = " 
                  << externalPrior->Integral() << ", max = " << externalPrior->GetMaximum() << std::endl;
        
        // Build the 2D ratio between original prior and externally provided prior
        hweights = static_cast<TH2D*>(externalPrior->Clone(("hweightsForPrior_Ext" + std::to_string(iteration)).c_str()));
        hweights->Divide(hOriginalPrior);
        
        std::cout << "After division: weights histogram integral = " << hweights->Integral() 
                  << ", max = " << hweights->GetMaximum() << ", min = " << hweights->GetMinimum() << std::endl;
        
        // Scale by maximum bin
        // Alternative option: normalize by integral
        // hweights->Scale(1.0/hweights->Integral());
        double oldMax = hOriginalPrior->GetMaximum();
        double newMax = hweights->GetMaximum();
        
        double scaleFactor = newMax!= 0 ? oldMax / newMax : 1.0; // Avoid division by zero
        if (scaleFactor < 0) {
            std::cerr << "ERROR: Scale factor is negative! oldMax = " << oldMax 
                      << ", newMax = " << newMax << std::endl;
            return nullptr;
        }
        
        std::cout << "Scaling weights by factor: " << scaleFactor 
                  << " (oldMax = " << oldMax << ", newMax = " << newMax << ")" << std::endl;
        
        hweights->Scale(scaleFactor);
        
        std::cout << "After scaling: weights histogram integral = " << hweights->Integral() 
                  << ", max = " << hweights->GetMaximum() << ", min = " << hweights->GetMinimum() << std::endl;
        
        // Check for problematic values
        int zeroBins = 0, negBins = 0, largeBins = 0;
        for (int i = 1; i <= hweights->GetNbinsX(); i++) {
            for (int j = 1; j <= hweights->GetNbinsY(); j++) {
                double content = hweights->GetBinContent(i, j);
                if (content == 0) zeroBins++;
                if (content < 0) negBins++;
                if (content > 10) largeBins++;
            }
        }
        std::cout << "Weight histogram statistics: " << zeroBins << " zero bins, " 
                  << negBins << " negative bins, " << largeBins << " bins with weight > 10" << std::endl;
        
        plotHist(hweights, "hweights_extPrior", "colz", "p_{T,jet} (GeV/c)", "z_{T}");
    } else {
        std::cout << "Using adaptive unfolding with weights from previous iteration" << std::endl;
        
        // Check validity of source histograms
        if (!hOriginalPrior) {
            std::cerr << "ERROR: hOriginalPrior is null!" << std::endl;
            return nullptr;
        }
        
        if (regParam >= unfoldedArr2D.size()) {
            std::cerr << "ERROR: regParam (" << regParam << ") out of range for unfoldedArr2D (size: " 
                      << unfoldedArr2D.size() << ")" << std::endl;
            return nullptr;
        }
        
        if (!unfoldedArr2D[regParam]) {
            std::cerr << "ERROR: unfoldedArr2D[" << regParam << "] is null!" << std::endl;
            return nullptr;
        }
        
        // Debug histogram properties
        std::cout << "Original prior: " << hOriginalPrior->GetNbinsX() << "x" 
                  << hOriginalPrior->GetNbinsY() << " bins, integral = " 
                  << hOriginalPrior->Integral() << ", max = " << hOriginalPrior->GetMaximum() << std::endl;
        
        std::cout << "Unfolded result: " << unfoldedArr2D[regParam]->GetNbinsX() << "x" 
                  << unfoldedArr2D[regParam]->GetNbinsY() << " bins, integral = " 
                  << unfoldedArr2D[regParam]->Integral() << ", max = " << unfoldedArr2D[regParam]->GetMaximum() << std::endl;
        
        // Build the 2D ratio between prior and reconstructed data after unfolding
        hweights = static_cast<TH2D*>(unfoldedArr2D[regParam]->Clone(("hweightsForPrior_" + std::to_string(iteration)).c_str()));
        hweights->Divide(hOriginalPrior);
        
        std::cout << "After division: weights histogram integral = " << hweights->Integral() 
                  << ", max = " << hweights->GetMaximum() << ", min = " << hweights->GetMinimum() << std::endl;
        
        // Scale by maximum bin
        double oldMax = hOriginalPrior->GetMaximum();
        double newMax = hweights->GetMaximum();
        double scaleFactor = oldMax/newMax;
        
        std::cout << "Scaling weights by factor: " << scaleFactor 
                  << " (oldMax = " << oldMax << ", newMax = " << newMax << ")" << std::endl;
        
        hweights->Scale(scaleFactor);
        
        std::cout << "After scaling: weights histogram integral = " << hweights->Integral() 
                  << ", max = " << hweights->GetMaximum() << ", min = " << hweights->GetMinimum() << std::endl;
        
        // Check for problematic values
        int zeroBins = 0, negBins = 0, largeBins = 0;
        for (int i = 1; i <= hweights->GetNbinsX(); i++) {
            for (int j = 1; j <= hweights->GetNbinsY(); j++) {
                double content = hweights->GetBinContent(i, j);
                if (content == 0) zeroBins++;
                if (content < 0) negBins++;
                if (content > 10) largeBins++;
            }
        }
        std::cout << "Weight histogram statistics: " << zeroBins << " zero bins, " 
                  << negBins << " negative bins, " << largeBins << " bins with weight > 10" << std::endl;
        
        plotHist(hweights, "hweights_Iter" + std::to_string(iteration), "colz", "p_{T,jet} (GeV/c)", "z_{T}");
    }
    
    std::cout << "==== prepareRMWeights END ====\n" << std::endl;
    return hweights;
}

// Implementation of getEventWeight method
double UnfoldSpectraClass::getEventWeight(TH2D* hweights, double zTValue, double pTValue) {
    double weight = 1.0;
    
    int xBin = hweights->GetXaxis()->FindBin(pTValue);
    int yBin = hweights->GetYaxis()->FindBin(zTValue);
    weight = hweights->GetBinContent(xBin, yBin);
    
    return weight;
}

// Implementation of getResponseMatrix method (simplified - full implementation would be more complex)
TH2D* UnfoldSpectraClass::getResponseMatrix(int part, const std::string& xAxisVar, const std::string& yAxisVar, 
                                           bool fineBinning, int bin, bool isCut, TFile* externalFile) {
    // Determine mass cut range
    std::pair<double, double> MassCut = std::make_pair(1.81, 1.935);
    
    // Create filter strings
    std::string resonanceCutString = "tagMPart > " + std::to_string(MassCut.first) + 
                                     " && tagMPart < " + std::to_string(MassCut.second);
    std::string jetCuts = "nConstPart>1 && etaDet>2.5 && etaDet<4";
    
    // Create binning arrays
    std::vector<double> xBins, yBins;
    int nBinsX = 0, nBinsY = 0;
    
    // Determine binning based on axis variables
    if (xAxisVar == "jetPtPart") {
        xBins.clear();
        for (int val : pTBinsArrayTruth) {
            xBins.push_back(static_cast<double>(val));
        }
        nBinsX = xBins.size() - 1;
    } else if (xAxisVar == "jetPtDet") {
        xBins.clear();
        for (int val : pTBinsArrayDet) {
            xBins.push_back(static_cast<double>(val));
        }
        nBinsX = xBins.size() - 1;
    } else if (xAxisVar == "zTDet") {
        xBins = zBinsArrayDet;
        nBinsX = xBins.size() - 1;
    }
    
    if (yAxisVar == "zTPart") {
        yBins = zBinsArrayTruth;
        nBinsY = yBins.size() - 1;
    } else if (yAxisVar == "zTDet") {
        yBins = zBinsArrayDet;
        nBinsY = yBins.size() - 1;
    }
    
    // Create histogram
    std::string histName = "RM_" + xAxisVar + "_" + yAxisVar;
    if (part == 0) {
        histName += "Full";
    } else if (part == 1) {
        histName += "Half1";
    } else {
        histName += "Half2";
    }
    
    TH2D* RM = new TH2D(histName.c_str(), histName.c_str(), nBinsX, xBins.data(), nBinsY, yBins.data());
    
    // Fill histogram from TTree (would need to implement the actual filling)
    // This would involve retrieving the TTree and looping through entries
    // with appropriate filters based on the 'part' parameter
    
    // Set Sumw2 for correct error propagation
    RM->Sumw2();
    
    return RM;
}
// Example implementation of one method (continuing with other methods would make this answer too long)
void UnfoldSpectraClass::unfold2D(int regParam, int iteration, const std::string& tag) {
    std::cout << "\n=============== UNFOLD2D DEBUG START ===============" << std::endl;
    std::cout << "Input parameters: regParam = " << regParam << ", iteration = " << iteration << ", tag = '" << tag << "'" << std::endl;

    unfoldLabel = "2DBayes" + std::to_string(regParam) + "_Round" + std::to_string(iteration);
    std::cout << "Setting unfoldLabel = " << unfoldLabel << std::endl;

    // Update output path
    std::string oldPath = outPath;
    if (isPrompt) {
        outPath = outpathBase + "/Prompt" + tag + "_" + std::to_string(iteration) + "/";
    } else {
        outPath = outpathBase + "/NonPrompt" + tag + "_" + std::to_string(iteration) + "/";
    }
    std::cout << "Updated output path from " << oldPath << " to " << outPath << std::endl;
    
    if (!std::filesystem::exists(outPath)) {
        std::cout << "Creating output directory: " << outPath << std::endl;
        std::filesystem::create_directories(outPath);
    }
    
    // Prepare response matrix
    TH2D* hweights = nullptr;
    if (iteration > 0) {
        std::cout << "Iteration > 0, preparing RM weights..." << std::endl;
        hweights = prepareRMWeights(regParam, iteration);
        if (!hweights) {
            std::cerr << "ERROR: Failed to get weights histogram from prepareRMWeights!" << std::endl;
        } else {
            std::cout << "Successfully created weights histogram with integral = " << hweights->Integral() << std::endl;
        }
    } else {
        std::cout << "First iteration (iteration = 0), not using weights for initial RM" << std::endl;
    }
    
    std::cout << "Preparing main response matrix (part 0)..." << std::endl;
    RooUnfoldResponse* responseMatrix3D = PrepareResponseMatrix3D(0, hweights);
    
    if (!responseMatrix3D) {
        std::cerr << "ERROR: Failed to create main response matrix!" << std::endl;
        return;
    }
    
    std::cout << "Response matrix information:" << std::endl;
    std::cout << " - Measured bins: " << responseMatrix3D->Hmeasured()->GetNbinsX() << "x" 
              << responseMatrix3D->Hmeasured()->GetNbinsY() << std::endl;
    std::cout << " - Truth bins: " << responseMatrix3D->Htruth()->GetNbinsX() << "x" 
              << responseMatrix3D->Htruth()->GetNbinsY() << std::endl;
    std::cout << " - Measured integral: " << responseMatrix3D->Hmeasured()->Integral() << std::endl;
    std::cout << " - Truth integral: " << responseMatrix3D->Htruth()->Integral() << std::endl;
    
    if (iteration == 0) {
        std::cout << "First iteration, saving original prior" << std::endl;
        hOriginalPrior = static_cast<TH2D*>(responseMatrix3D->Htruth()->Clone("hOriginalPrior"));
    } else {
        std::cout << "Updating current prior for iteration " << iteration << std::endl;
        hCurrentPrior = static_cast<TH2D*>(responseMatrix3D->Htruth()->Clone("hCurrentPrior"));
    }
    
    std::cout << "Preparing secondary response matrices (parts 1 & 2)..." << std::endl;
    RooUnfoldResponse* responseMatrix3D_Pt1 = PrepareResponseMatrix3D(1, hweights);
    RooUnfoldResponse* responseMatrix3D_Pt2 = PrepareResponseMatrix3D(2, hweights);
    
    if (!responseMatrix3D_Pt1 || !responseMatrix3D_Pt2) {
        std::cerr << "ERROR: Failed to create secondary response matrices!" << std::endl;
        delete responseMatrix3D;
        if (hweights) delete hweights;
        return;
    }
    
    if (iteration == 0 && externalPrior) {
        std::cout << "First iteration with external prior, preparing all RMs again with weights" << std::endl;
        hweights = prepareRMWeights(regParam, iteration);
        if (!hweights) {
            std::cerr << "ERROR: Failed to get weights histogram for external prior!" << std::endl;
        } else {
            std::cout << "Created weights from external prior with integral = " << hweights->Integral() << std::endl;
            
            std::cout << "Recreating all response matrices with external prior weights..." << std::endl;
            delete responseMatrix3D;
            delete responseMatrix3D_Pt1;
            delete responseMatrix3D_Pt2;
            
            responseMatrix3D = PrepareResponseMatrix3D(0, hweights);
            responseMatrix3D_Pt1 = PrepareResponseMatrix3D(1, hweights);
            responseMatrix3D_Pt2 = PrepareResponseMatrix3D(2, hweights);
            
            if (!responseMatrix3D || !responseMatrix3D_Pt1 || !responseMatrix3D_Pt2) {
                std::cerr << "ERROR: Failed to recreate response matrices with external prior!" << std::endl;
                if (hweights) delete hweights;
                return;
            }
        }
    }
    
    // Plot measured data
    std::cout << "Plotting measured data histogram" << std::endl;
    if (measuredSpectra2D) {
        std::cout << "Measured data 2D: " << measuredSpectra2D->GetNbinsX() << "x" 
                  << measuredSpectra2D->GetNbinsY() << " bins, integral = " 
                  << measuredSpectra2D->Integral() << std::endl;
        plotHist(measuredSpectra2D, "hMeasuredData2D", "colz", "p_{T,jet} (GeV/c)", "z_{T}");
    } else {
        std::cerr << "ERROR: measuredSpectra2D is null!" << std::endl;
    }
    
    // Clear previous unfolding results
    std::cout << "Clearing previous unfolding results (" << unfoldedArr2D.size() << " histograms)" << std::endl;
    for (auto hist : unfoldedArr2D) {
        if (hist) delete hist;
    }
    unfoldedArr2D.clear();
    
    // Perform unfolding for different iteration values
    for (int i = 1; i <= regParam + 2; i++) {
        std::cout << "\nPerforming unfolding with " << i << " Bayes iterations..." << std::endl;
        
        // Set up the Bayesian unfolding
        RooUnfoldBayes unfoldBayes2D(responseMatrix3D, measuredSpectra2D, i);
        
        // Perform the unfolding
        std::cout << "Running unfolding with error treatment: " << errorType << std::endl;
        TH2D* h2DUnfoldedPerBin = dynamic_cast<TH2D*>(unfoldBayes2D.Hreco(
            static_cast<RooUnfold::ErrorTreatment>(errorType)));
            
        if (!h2DUnfoldedPerBin) {
            std::cerr << "ERROR: Unfolding failed - returned null histogram for " << i << " iterations" << std::endl;
            continue;
        }
        
        std::string histName = "hUnfoldedData2DPerBin_nIter" + std::to_string(i);
        h2DUnfoldedPerBin->SetName(histName.c_str());
        std::cout << "Unfolded result: " << h2DUnfoldedPerBin->GetNbinsX() << "x" 
                  << h2DUnfoldedPerBin->GetNbinsY() << " bins, integral = " 
                  << h2DUnfoldedPerBin->Integral() << ", max = " << h2DUnfoldedPerBin->GetMaximum() << std::endl;
        
        // Plot unfolded histogram
        plotHist(h2DUnfoldedPerBin, "hUnfoldedData2D_nIter" + std::to_string(i), "colz", 
                "p_{T,jet} (GeV/c)", "z_{T}");
        unfoldedArr2D.push_back(h2DUnfoldedPerBin);
        
        // Plot correlation coefficients
        std::cout << "Getting covariance matrix for iteration " << i << "..." << std::endl;
        TMatrixD covarianceMatrix = unfoldBayes2D.Ereco(
            static_cast<RooUnfold::ErrorTreatment>(errorType));
        std::cout << "Covariance matrix dimensions: " << covarianceMatrix.GetNrows() << " x " 
                  << covarianceMatrix.GetNcols() << std::endl;
                  
        // Check if matrix elements are finite
        bool hasNaNs = false;
        for (int row = 0; row < covarianceMatrix.GetNrows(); row++) {
            for (int col = 0; col < covarianceMatrix.GetNcols(); col++) {
                if (!std::isfinite(covarianceMatrix(row, col))) {
                    hasNaNs = true;
                    break;
                }
            }
            if (hasNaNs) break;
        }
        
        if (hasNaNs) {
            std::cerr << "WARNING: Covariance matrix contains NaN or Inf values!" << std::endl;
        }
        
        std::cout << "Plotting correlation coefficients..." << std::endl;
        plotCorrelationCoefficients(covarianceMatrix, i, "CorMatr2D");
        
        // Perform refolding test
        std::cout << "Performing refolding test for iteration " << i << "..." << std::endl;
        RefoldingTest2D(i, responseMatrix3D_Pt1, responseMatrix3D_Pt2);
    }
    
    // Perform additional tests
    std::cout << "\nPerforming additional unfolding tests..." << std::endl;
    std::cout << "Running UnfoldingTest2D with regParam = " << regParam << std::endl;
    UnfoldingTest2D(responseMatrix3D, regParam);
    
    std::cout << "Running StabilityTest2D with regParam = " << regParam << std::endl;
    StabilityTest2D(regParam);
    
    std::cout << "Running StabilityTest2D with regParam = " << (regParam + 1) << std::endl;
    StabilityTest2D(regParam + 1);
    
    std::cout << "Running TestRegParam2D" << std::endl;
    TestRegParam2D();
    
    if (iteration == 0) {
        std::cout << "First iteration, running StatTestRM2D" << std::endl;
        StatTestRM2D(responseMatrix3D);
    }
    
    // Plot prior vs unfolded result
    std::cout << "Plotting prior vs unfolded result..." << std::endl;
    plotPrior2D(regParam, responseMatrix3D);
    
    // Plot before/after unfolding
    std::cout << "Plotting before/after unfolding comparison..." << std::endl;
    plotUnfoldingEffect2D(regParam);
    
    // Save results
    std::cout << "Saving results for different regularization parameters..." << std::endl;
    std::cout << "Saving result for regParam = " << (regParam - 1) << std::endl;
    saveResult(regParam - 1, iteration, tag);
    std::cout << "Saving result for regParam = " << regParam << std::endl;
    saveResult(regParam, iteration, tag);
    std::cout << "Saving result for regParam = " << (regParam + 1) << std::endl;
    saveResult(regParam + 1, iteration, tag);
    
    // Plot kinematic efficiency
    std::cout << "Plotting kinematic efficiency..." << std::endl;
    getKinEfficiency(0);
    
    // Clean up
    std::cout << "Cleaning up response matrices..." << std::endl;
    delete responseMatrix3D;
    delete responseMatrix3D_Pt1;
    delete responseMatrix3D_Pt2;
    if (hweights) delete hweights;
    
    std::cout << "=============== UNFOLD2D DEBUG END ===============\n" << std::endl;
}

void UnfoldSpectraClass::TestRegParam2D() {
    std::cout << "\n==== TestRegParam2D DEBUG START ====" << std::endl;
    
    // Check if we have unfolded histograms to analyze
    if (unfoldedArr2D.empty()) {
        std::cerr << "ERROR: unfoldedArr2D array is empty, cannot test regularization parameters" << std::endl;
        return;
    }
    
    std::cout << "Testing regularization parameters for " << unfoldedArr2D.size() 
              << " iterations" << std::endl;
    
    // Create output directory for plots
    std::string outputDir = makeOutDir("RegParamTest");
    std::cout << "Created output directory: " << outputDir << std::endl;

    // Define colors for different iterations
    std::vector<int> colorArray = {kRed-9, kGreen-9, kGreen-8, kBlue-9, kGreen, 
                                  kSpring+8, kSpring, kRed, kRed+2, kRed+5};
    
    // Get number of pT bins
    int npTBins = pTBinsArrayTruth.size() - 1;
    std::cout << "Working with " << npTBins << " pT bins" << std::endl;

    // Create arrays of histograms to store relative error as function of zT for all iterations
    std::vector<std::vector<TH1D*>> hArrpTBins(npTBins + 1);
    
    // Create a dummy histogram to clone for storing relative errors
    TH1D* dummyHisto = unfoldedArr2D[0]->ProjectionY("_pymain", 1, 1);
    dummyHisto->Reset("ICESM");
    
    // Initialize histograms for each pT bin and iteration
    std::cout << "Initializing histograms for " << unfoldedArr2D.size() << " iterations and " 
              << (npTBins + 1) << " pT bins" << std::endl;
    
    for (size_t iterStep = 0; iterStep < unfoldedArr2D.size(); iterStep++) {
        for (int pTBin = 0; pTBin <= npTBins; pTBin++) {
            TH1D* hist = static_cast<TH1D*>(dummyHisto->Clone(
                Form("hpTBin%d_iter%zu", pTBin, iterStep)));
            hist->SetTitle(Form("Relative error for pT bin %d, iteration %zu", pTBin, iterStep));
            hist->GetXaxis()->SetTitle("z_{T}");
            hist->GetYaxis()->SetTitle("Relative error");
            hArrpTBins[pTBin].push_back(hist);
        }
    }
    
    delete dummyHisto; // Clean up the dummy histogram
    
    // Evaluate error change for different iterations
    std::cout << "Evaluating relative errors for each iteration and pT bin..." << std::endl;
    
    for (size_t iterStep = 0; iterStep < unfoldedArr2D.size(); iterStep++) {
        std::cout << "Processing iteration " << (iterStep + 1) << "..." << std::endl;
        
        if (!unfoldedArr2D[iterStep]) {
            std::cerr << "ERROR: Null histogram found for iteration " << (iterStep + 1) << std::endl;
            continue;
        }
        
        for (int pTBin = 1; pTBin <= npTBins; pTBin++) {
            // Project the 2D histogram to get the 1D distribution for this pT bin
            TH1D* proj = unfoldedArr2D[iterStep]->ProjectionY(
                Form("_pymain_iter%zu_pT%d", iterStep, pTBin), pTBin, pTBin);
                
            if (!proj) {
                std::cerr << "ERROR: Projection failed for iteration " << (iterStep + 1) 
                          << ", pT bin " << pTBin << std::endl;
                continue;
            }
            
            std::cout << "  pT bin " << pTBin << ": " << proj->GetNbinsX() << " bins, integral = "
                      << proj->Integral() << std::endl;
                
            // Calculate relative error for each zT bin
            for (int zTBin = 1; zTBin <= proj->GetNbinsX(); zTBin++) {
                double content = proj->GetBinContent(zTBin);
                double error = proj->GetBinError(zTBin);
                double relError = 0;
                
                if (content > 0) {
                    relError = error / content;
                    hArrpTBins[pTBin][iterStep]->SetBinContent(zTBin, relError);
                    
                    // Debug output for first few bins
                    if (zTBin <= 3) {
                        std::cout << "    zT bin " << zTBin << ": content = " << content 
                                  << ", error = " << error << ", rel.error = " << relError << std::endl;
                    }
                }
            }
            
            delete proj; // Clean up the projection
        }
    }
    
    // Create canvas for the plots
    TCanvas* c = new TCanvas("c", "Regularization Parameter Test", 800 * npTBins, 800);
    c->Divide(npTBins);
    
    // Create legend
    TLegend* leg = new TLegend(0.2, 0.7, 0.5, 0.93, "");
    leg->SetFillColor(10);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(hArrpTBins[1][0], "nIter=1", "l");
    
    // Plot relative errors for each pT bin
    for (int pTBin = 1; pTBin <= npTBins; pTBin++) {
        c->cd(pTBin);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.05);
        gPad->SetBottomMargin(0.15);
        
        // Draw first iteration
        hArrpTBins[pTBin][0]->GetYaxis()->SetRangeUser(0, 0.5); // Adjust range as needed
        hArrpTBins[pTBin][0]->SetLineColor(kBlue);
        hArrpTBins[pTBin][0]->SetLineWidth(2);
        hArrpTBins[pTBin][0]->SetMarkerStyle(20);
        hArrpTBins[pTBin][0]->Draw("hist");
        
        // Draw other iterations
        for (size_t iterStep = 1; iterStep < unfoldedArr2D.size(); iterStep++) {
            hArrpTBins[pTBin][iterStep]->SetLineColor(colorArray[iterStep % colorArray.size()]);
            hArrpTBins[pTBin][iterStep]->SetMarkerStyle(20 + iterStep);
            hArrpTBins[pTBin][iterStep]->Draw("hist same");
            
            // Add to legend for first pT bin only
            if (pTBin == 1) {
                leg->AddEntry(hArrpTBins[pTBin][iterStep], 
                             Form("nIter=%zu", iterStep + 1), "l");
            }
        }
        
        // Draw title for this pT bin
        TPaveText* paveText = new TPaveText(0.15, 0.92, 0.7, 0.99, "NDC");
        paveText->SetBorderSize(0);
        paveText->SetFillStyle(0);
        paveText->SetTextFont(42);
        paveText->SetTextSize(0.04);
        paveText->AddText(Form("p_{T} bin: %d-%d GeV/c", 
                             pTBinsArrayTruth[pTBin-1], pTBinsArrayTruth[pTBin]));
        paveText->Draw();
        
        // Draw legend for first pT bin only
        if (pTBin == 1) {
            leg->Draw("same");
        }
    }
    
    // Save the canvas
    std::string fileName = outputDir + "/RegParamTest2D_" + figureTag + ".png";
    std::cout << "Saving plot to: " << fileName << std::endl;
    c->SaveAs(fileName.c_str());
    
    // Clean up
    delete c;
    
    // Clean up histograms
    std::cout << "Cleaning up " << hArrpTBins.size() << " histogram arrays" << std::endl;
    for (auto& hArray : hArrpTBins) {
        for (auto hist : hArray) {
            delete hist;
        }
    }
    
    std::cout << "==== TestRegParam2D DEBUG END ====\n" << std::endl;
}

void UnfoldSpectraClass::plotCorrelationCoefficients(TMatrixD& covarianceMatrix, int i, const std::string& Name) {
    int nBinsX = covarianceMatrix.GetNrows();
    int nBinsY = covarianceMatrix.GetNcols();
      
    TH2D* correlationCoefficientMatrix = new TH2D("correlationCoefficientMatrix", 
                                                 "correlationCoefficientMatrix", 
                                                 nBinsX, 0, nBinsX, 
                                                 nBinsY, 0, nBinsY);

    for (int xbin = 0; xbin < nBinsX; xbin++) {
        double varianceX = covarianceMatrix(xbin, xbin);
        double sigmaX = sqrt(varianceX);

        for (int ybin = 0; ybin < nBinsY; ybin++) {
            double varianceY = covarianceMatrix(ybin, ybin);
            double sigmaY = sqrt(varianceY);

            double covXY = covarianceMatrix(xbin, ybin);
            if (sigmaX > 0 && sigmaY > 0) {
                double Cxy = covXY / (sigmaX * sigmaY);
                correlationCoefficientMatrix->SetBinContent(xbin+1, ybin+1, Cxy);
            }
        }
    }

    plotHist(correlationCoefficientMatrix, "hCorrelationCoefficientMatrix" + std::to_string(i), "colz");
    delete correlationCoefficientMatrix;
}

void UnfoldSpectraClass::provideExtPrior(const std::string& tag, const std::string& fileName, const std::string& specialType) {
    // Return if fileName is empty
    if (fileName.empty()) return;
    
    // Construct file path and open the file
    std::string inFileNameExt_RM = inPathRM + fileName + ".root";
    TFile* inFileExt_RM = new TFile(inFileNameExt_RM.c_str());
    
    if (!inFileExt_RM || inFileExt_RM->IsZombie()) {
        std::cerr << "Error: Could not open external prior file: " << inFileNameExt_RM << std::endl;
        if (inFileExt_RM) delete inFileExt_RM;
        return;
    }
    
    // Get the response matrix from the file
    externalPrior = getResponseMatrix(0, "jetPtPart", "zTPart", false, 0, false, inFileExt_RM);
    
    // Provide a flat prior if requested
    if (specialType == "Flat") {
        // Set all bin entries for a flat prior
        for (int binx = 1; binx <= externalPrior->GetNbinsX(); binx++) {
            for (int biny = 1; biny <= externalPrior->GetNbinsY(); biny++) {
                // Get bin widths
                double xWidth = externalPrior->GetXaxis()->GetBinWidth(binx);
                double yWidth = externalPrior->GetYaxis()->GetBinWidth(biny);
                
                // Set bin content to be proportional to area (flat per area)
                externalPrior->SetBinContent(binx, biny, xWidth * yWidth);
            }
        }
    }
    
    // Disconnect histogram from file to prevent memory issues
    if (externalPrior) {
        externalPrior->SetDirectory(0);
    }
    
    // Close the file
    inFileExt_RM->Close();
    delete inFileExt_RM;
    
    // Update output path
    if (isPrompt) {
        outPath = outpathBase + "/Prompt" + tag + "_0/";
    } else {
        outPath = outpathBase + "/NonPrompt" + tag + "_0/";
    }
    
    // Create the output directory if it doesn't exist
    if (!std::filesystem::exists(outPath)) {
        std::filesystem::create_directories(outPath);
    }
    
    // Plot the external prior
    plotHist(externalPrior, "h_extPrior", "colz", "p_{T,jet} (GeV/c)", "z_{T}");
}

void UnfoldSpectraClass::plotHist(TObject* hist, const std::string& outputFilename, 
                               const std::string& drawOptions, 
                               const std::string& xTitle, 
                               const std::string& yTitle) {
    // Set ROOT style parameters
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    
    // Create canvas
    TCanvas* c = new TCanvas("c", "c: hist", 600, 450);
    c->cd();
    c->SetLeftMargin(0.15);
    
    // Handle different types of histograms
    if (hist->InheritsFrom(TH1::Class())) {
        if (hist->InheritsFrom(TH2::Class())) {
            // 2D histograms
            TH2* h2 = dynamic_cast<TH2*>(hist);
            h2->SetNdivisions(505);
            if (!xTitle.empty()) h2->SetXTitle(xTitle.c_str());
            if (!yTitle.empty()) h2->SetYTitle(yTitle.c_str());
            h2->Draw(drawOptions.c_str());
        } else {
            // 1D histograms
            TH1* h1 = dynamic_cast<TH1*>(hist);
            h1->SetNdivisions(505);
            h1->Draw(drawOptions.c_str());
        }
    } else {
        // Other objects (like TGraph)
        TGraph* gr = dynamic_cast<TGraph*>(hist);
        if (gr) {
            TH1* myBlankHisto = gr->GetHistogram();
            if (myBlankHisto) {
                myBlankHisto->SetNdivisions(505);
                myBlankHisto->GetXaxis()->SetRangeUser(0, 1);
                myBlankHisto->GetYaxis()->SetRangeUser(0, gr->GetMaximum() * 1.2);
                myBlankHisto->Draw("E");
                gr->Draw(("same " + drawOptions).c_str());
            } else {
                gr->Draw(drawOptions.c_str());
            }
        }
    }
    
    // Save the canvas to a file
    std::string fileName = outPath + outputFilename + ".png";
    // c->SaveAs(fileName.c_str());
    SaveCanvasQuietly(c, fileName.c_str());
    
    // Clean up
    delete c;
}

// Main function implementations
void unfoldzT(int variation) {
    int maxIter = 4;
    std::string extFilename_P = "";
    std::string extFilename_NP = "";
    std::string priorType = "";
    std::string tag = "";
    
    // Default file names
    std::string FileRM_P = "20250514_Pbp_21_MC_output_D0FF_filterV1";
    std::string FileRM_NP = "20250514_Pbp_21_MC_output_D0FF_filterV1";
    // std::string FileRM_P = "All_PromptD0";
    // std::string FileRM_NP = "16_17_18AllBdecay_ScaledNom_V2";
    
    // Apply variations
    switch(variation) {
        case 1:
            tag = "SWP_prior";
            extFilename_P = FileRM_NP;
            extFilename_NP = FileRM_P;
            maxIter = 1;
            break;
        case 2:
            tag = "Flat_prior";
            extFilename_P = FileRM_NP;
            extFilename_NP = FileRM_P;
            priorType = "Flat";
            maxIter = 1;
            break;
        case 3:
            tag = "moreHyperons";
            FileRM_NP = "16_17_18AllBdecay_Scaled_p20pc_V2";
            maxIter = 1;
            break;
        case 4:
            tag = "lessHyperons";
            FileRM_NP = "16_17_18AllBdecay_Scaled_m20pc_V2";
            maxIter = 1;
            break;
    }
    
    // D0 analysis
    std::vector<int> pTBinArray = {2, 5, 10, 15, 20, 30};
    int regularizationParam = 4;
    // int regularizationParam = 3;
    
    // Create and run unfolding for prompt
    UnfoldSpectraClass unfoldObjectPrompt("P", FileRM_P, "CorrectedFinalHistograms_D0", "D0", pTBinArray);
    unfoldObjectPrompt.provideExtPrior(tag, extFilename_P, priorType);
    
    for (int iter = 0; iter < maxIter; iter++) {
        std::cout << "Running unfold2D for prompt D0, iteration " << iter << std::endl;
        unfoldObjectPrompt.unfold2D(regularizationParam, iter, tag);
        std::cout << "Completed unfold2D for prompt D0, iteration " << iter << std::endl;
    }
    std::cout << "==== Unfolding for prompt D0 completed ====\n" << std::endl;
    std::cout << "==== Unfolding for non-prompt D0 started ====\n" << std::endl;
    
    // Create and run unfolding for non-prompt
    UnfoldSpectraClass unfoldObjectNonPrompt("NP", FileRM_NP, "CorrectedFinalHistograms_D0", "D0", pTBinArray);
    unfoldObjectNonPrompt.provideExtPrior(tag, extFilename_NP, priorType);
    
    for (int iter = 0; iter < maxIter; iter++) {
        std::cout << "Running unfold2D for non-prompt D0, iteration " << iter << std::endl;
        unfoldObjectNonPrompt.unfold2D(regularizationParam, iter, tag);
        std::cout << "Completed unfold2D for non-prompt D0, iteration " << iter << std::endl;
    }
    std::cout << "==== Unfolding for non-prompt D0 completed ====\n" << std::endl;
}