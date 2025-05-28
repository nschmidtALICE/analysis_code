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
    delete c;
    delete fRandom;
    delete histo;
    
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
                
                delete proj; // Clean up projection
                
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
                delete proj;
                
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
                
                delete gaussFunc;
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
        delete cA;
        
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
    int npTBins = pTBinsArrayDet.size() - 1;
    TCanvas* c = new TCanvas("c", "c: pT", 1200, 800);
    c->Divide(npTBins, 2);
    std::cout << "Created canvas divided into " << npTBins << "x2 pads" << std::endl;
    
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
        std::cout << __LINE__ << std::endl;
        stdDevHistos.push_back(stdDevHisto);
        std::cout << __LINE__ << std::endl;
        
        // Array to collect ratio histograms for this pT bin
        std::vector<TH1D*> ratioHistos;
        std::cout << __LINE__ << std::endl;
        
        // Create upper pad for original histogram
        c->cd(ptBin);
        std::cout << __LINE__ << std::endl;
        unfoldedProj->SetLineColor(kBlack);
        std::cout << __LINE__ << std::endl;
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
                delete smearedData;
                continue;
            }
            
            std::cout << "    Performing unfolding with " << nIter << " iterations..." << std::endl;
            RooUnfoldBayes unfoldBayes(toyRM, smearedData, nIter);
            TH2D* toyUnfolded = dynamic_cast<TH2D*>(unfoldBayes.Hreco(
                static_cast<RooUnfold::ErrorTreatment>(errorType)));
                
            if (!toyUnfolded) {
                std::cerr << "    ERROR: Toy unfolding failed for ptBin=" << ptBin << ", toy=" << toy << std::endl;
                delete smearedData;
                delete toyRM;
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
            delete smearedData;
            delete toyRM;
            delete toyUnfolded;
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
        
        // Create lower pad for ratio histos        c->cd(npTBins + ptBin);
        
        // Set up blank histogram for ratio plot
        TH1D* ratioFrame = static_cast<TH1D*>(unfoldedProj->Clone("ratioFrame"));
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
        for (auto hist : ratioHistos) {
            delete hist;
        }
        ratioHistos.clear();
        std::cout << "Cleaned up " << ratioHistos.capacity() << " ratio histograms" << std::endl;
        
        // Clean up original projection
        delete unfoldedProj;
    }
    
    // Save canvas
    std::string fileName = outputDirStability + "/StabilityTest2D_nIter" + std::to_string(nIter) + ".png";
    std::cout << "Saving plot to: " << fileName << std::endl;
    c->SaveAs(fileName.c_str());
    
    // Clean up
    delete c;
    
    // Clean up histograms
        std::cout << "Cleaning up " << meanHistos.size() << " mean histograms and " 
              << stdDevHistos.size() << " stdDev histograms" << std::endl;
    
    // Fix: Use meanHistos and stdDevHistos for cleanup
    for (auto hist : meanHistos) {
        delete hist;
    }
    meanHistos.clear();
    
    for (auto hist : stdDevHistos) {
        delete hist;
    }
    stdDevHistos.clear();
    
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
    c->SaveAs(fileName.c_str());
    
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
    int regularizationParam = 6;
    // int regularizationParam = 3;
    
    // Create and run unfolding for prompt
    UnfoldSpectraClass unfoldObjectPrompt("P", FileRM_P, "CorrectedFinalHistograms_D0", "D0", pTBinArray);
    unfoldObjectPrompt.provideExtPrior(tag, extFilename_P, priorType);
    
    for (int iter = 0; iter < maxIter; iter++) {
        unfoldObjectPrompt.unfold2D(regularizationParam, iter, tag);
    }
    
    // Create and run unfolding for non-prompt
    UnfoldSpectraClass unfoldObjectNonPrompt("NP", FileRM_NP, "CorrectedFinalHistograms_D0", "D0", pTBinArray);
    unfoldObjectNonPrompt.provideExtPrior(tag, extFilename_NP, priorType);
    
    for (int iter = 0; iter < maxIter; iter++) {
        unfoldObjectNonPrompt.unfold2D(regularizationParam, iter, tag);
    }
}