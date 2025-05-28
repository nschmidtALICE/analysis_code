// Update the implementation signature to match

TH1* Plotter::ipchi2FitPlot(const std::string& resonance, RooRealVar* logIpchi2, 
                           RooDataSet* data, RooAbsPdf* totalPdf, 
                           RooAbsPdf* nonpromptPdf, RooAbsPdf* promptPdf, 
                           RooAbsPdf* backgroundPdf,
                           RooAbsReal* promptYield, RooAbsReal* nonpromptYield) {
    // The rest of the implementation stays the same
    // ...

    // When using these parameters, make sure to access them as RooAbsReal
    if (promptYield && nonpromptYield) {
        promptVal = promptYield->getVal();
        nonpromptVal = nonpromptYield->getVal();
        sigYield = promptVal + nonpromptVal;
    }
    
    // ...
}