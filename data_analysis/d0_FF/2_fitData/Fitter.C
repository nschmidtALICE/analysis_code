// Fitter.cpp
#include "Fitter.h"
#include "Plotter.h"

Fitter::Fitter(TTree* tree, const std::string& resonanceType, int numBins,
               bool isMCData, const std::string& outputPath, bool update)
    : tfilePV(nullptr), TestFilename(""), 
      outfilePath(outputPath), isMC(isMCData), nBins(numBins),
      resonance(resonanceType), updateStartValues(update), inTree(tree), fInFileHisto(nullptr)
{
    // Initialize dictionary
    initDictionary();
    
    std::cout << "This is the MC status: " << (isMC ? "true" : "false") << std::endl;
    
    // Add D0 specific B-decay ranges
    BdecayMinVal = {0, 0, 0, 0, 0};
    BdecayMaxVal = {30, 25, 15, 10, 5}; // for zT
    
    // Check array lengths
    if (BdecayMinVal.size() != static_cast<size_t>(nBins) || BdecayMaxVal.size() != static_cast<size_t>(nBins)) {
        std::cout << "Error: Fix limits" << std::endl;
        std::cout << "number of mass bins= " << nBins << ", number of limits for B-decay fraction: " << BdecayMinVal.size() << std::endl;
    }
    
    // Verify the input tree
    if (inTree) {
        std::cout << "Using input tree with " << inTree->GetEntries() << " entries" << std::endl;
        
        // Get the file that owns the tree
        if (inTree->GetCurrentFile()) {
            std::cout << "Tree source file: " << inTree->GetCurrentFile()->GetName() << std::endl;
        }
    } else {
        std::cout << "Warning: Input tree is null" << std::endl;
    }
}

Fitter::~Fitter() {
    // Clean up
    // Note: Don't close the input files here as they might be used elsewhere
}

void Fitter::initDictionary() {
    // Initialize D0 mass parameters
    MassConfig d0Config;
    d0Config.sigma1 = ParamConfig(0.008, 0.005, 0.012);
    d0Config.deltasigma = ParamConfig(1.5, 1.1, 2.5);
    d0Config.mean = ParamConfig(1.865, 1.860, 1.870);
    d0Config.alpha1 = ParamConfig(2, 1, 5);
    d0Config.n = ParamConfig(1, 0.2, 5);
    d0Config.dg_frac = ParamConfig(0.5, 0.0, 0.99999);
    d0Config.pol0 = ParamConfig(100, -5e3, 5e3);
    d0Config.pol1 = ParamConfig(58, -2e2, 5e2);
    d0Config.pol2 = ParamConfig(0, 0, 0);
    d0Config.massRange = std::make_pair(1.81, 1.935);  // 150 MeV range
    d0Config.sigYield = ParamConfig(100, 1, 5e6);
    d0Config.sigYieldLim = ParamConfig(100, 1, 5e6);
    d0Config.bkgYield = ParamConfig(1000, 0, 4e6);
    d0Config.bkgYieldLim = ParamConfig(1000, 0, 4e6);
    d0Config.signalRegion = std::make_pair(1.845, 1.885);  // 40 MeV range
    d0Config.sbRegion = std::make_pair(1.820, 1.840);     // Sideband region
    
    massDict["D0"] = d0Config;
    
    // Initialize D0 IP chi2 parameters
    IPChi2Config d0IPConfig;
    d0IPConfig.logIpchi2Range = std::make_pair(-3, 5);
    
    // Prompt component
    d0IPConfig.xpPrompt = ParamConfig(0.0, -0.5, 1.0);
    d0IPConfig.sigmaPrompt = ParamConfig(0.7, 0.4, 2.0);
    d0IPConfig.xiPrompt = ParamConfig(0.0, -0.5, 0.5);
    d0IPConfig.rho1Prompt = ParamConfig(-0.2, -0.98, -0.01);
    d0IPConfig.rho2Prompt = ParamConfig(0.2, 0.01, 0.98);
    
    // Non-prompt component
    d0IPConfig.xpNonprompt = ParamConfig(2.5, 2.0, 4.0);
    d0IPConfig.sigmaNonprompt = ParamConfig(0.6, 0.4, 0.9);
    d0IPConfig.xiNonprompt = ParamConfig(0.1, 0.05, 0.5);
    d0IPConfig.rho1Nonprompt = ParamConfig(-0.3, -0.8, -0.1);
    d0IPConfig.rho2Nonprompt = ParamConfig(0.2, 0.01, 0.98);
    
    // Fraction
    d0IPConfig.promptFrac = ParamConfig(0.95, 0.8, 0.999);
    
    // Background
    d0IPConfig.bkgParam1 = ParamConfig(0.5, 0, 1);
    d0IPConfig.bkgParam2 = ParamConfig(0.5, 0, 1);
    
    ipchi2Dict["SigD0"] = d0IPConfig;
}

void Fitter::updateDictionary(RooAbsPdf* signalPdf, RooAbsData* data, const std::string& fitFunc) {
    if (!signalPdf || !data || resonance.empty()) return;
    
    auto& res = massDict[resonance];
    
    RooArgSet* paramSet = signalPdf->getParameters(*data);
    if (!paramSet) return;
    
    std::vector<std::string> keyList;
    
    if (fitFunc == "noSig") {
        keyList = {"pol0", "pol1", "pol2"};
    }
    if (fitFunc == "DGauss") {
        keyList = {"mean", "sigma1", "deltasigma", "dg_frac", "sig_yield", 
                   "bkg_yield", "pol0", "pol1", "pol2"};
    }
    else if (fitFunc == "SGauss") {
        keyList = {"mean", "sigma1", "sig_yield", "bkg_yield", "pol0", "pol1", "pol2"};
    }
    else {
        return;
    }
    
    // Update the parameters in our dictionary
    for (const auto& key : keyList) {
        RooRealVar* param = dynamic_cast<RooRealVar*>(paramSet->find(key.c_str()));
        if (param) {
            double newVal = param->getVal();
            
            // Update the appropriate parameter
            if (key == "mean") res.mean.value = newVal;
            else if (key == "sigma1") res.sigma1.value = newVal;
            else if (key == "deltasigma") res.deltasigma.value = newVal;
            else if (key == "alpha1") res.alpha1.value = newVal;
            else if (key == "n") res.n.value = newVal;
            else if (key == "dg_frac") res.dg_frac.value = newVal;
            else if (key == "sig_yield") res.sigYield.value = newVal;
            else if (key == "bkg_yield") res.bkgYield.value = newVal;
            else if (key == "pol0") res.pol0.value = newVal;
            else if (key == "pol1") res.pol1.value = newVal;
            else if (key == "pol2") res.pol2.value = newVal;
        }
    }
    
    delete paramSet;
}

void Fitter::fixNAlphaValue(const std::string& resonance, double ptLimLow) {
    auto& res = massDict[this->resonance];
    
    double startValAlpha = 0.0;
    double minValAlpha = 0.0;
    double maxValAlpha = 0.0;
    double startValN = 0.0;
    double minValN = 0.0;
    double maxValN = 0.0;
    
    if (resonance.find("Psi2S") != std::string::npos) {
        // These values were determined by fitting the combined mass spectrum
        startValAlpha = 2.17;   // 2.461+-0.044 (before Apr. 22)
        minValAlpha = startValAlpha - 3 * 0.007;
        maxValAlpha = startValAlpha + 3 * 0.007;
        startValN = 0.923;  // 0.554+-0.097
        minValN = startValN - 3 * 0.01;
        maxValN = startValN + 3 * 0.01;
    }
    else if (resonance.find("X3872") != std::string::npos) {
        startValAlpha = 2.464; // 1.625 (before Apr. 22)
        minValAlpha = startValAlpha - 3 * 0.025;
        maxValAlpha = startValAlpha + 3 * 0.025;
        startValN = 0.593;     // 3.356 (before Apr. 22)
        minValN = startValN - 3 * 0.033;
        maxValN = startValN + 3 * 0.033;
    }
    else if (resonance.find("D0") != std::string::npos) {
        startValAlpha = 2.5;
        minValAlpha = startValAlpha - 3 * 0.05;
        maxValAlpha = startValAlpha + 3 * 0.05;
        startValN = 0.5;
        minValN = startValN - 3 * 0.05;
        maxValN = startValN + 3 * 0.05;
    }
    
    // Update alpha1 parameter
    res.alpha1.value = startValAlpha;
    res.alpha1.min = minValAlpha;
    res.alpha1.max = maxValAlpha;
    
    // Update n parameter
    res.n.value = startValN;
    res.n.min = minValN;
    res.n.max = maxValN;
    
    // Special case for Psi2S dg_frac
    if (resonance.find("Psi2S") != std::string::npos) {
        res.dg_frac.value = 0.4;
        res.dg_frac.min = 0.0;
        res.dg_frac.max = 0.5;
    }
}

void Fitter::updateSigYield(const std::string& resonance) {
    std::cout << "Update SigYield limits:" << std::endl;
    auto& res = massDict[resonance];
    
    // Signal yield
    std::cout << "Sig yield Before: " << res.sigYield.value << ", " 
              << res.sigYield.min << ", " << res.sigYield.max << std::endl;
    
    res.sigYield.min = res.sigYield.value * 0.8;
    if (res.sigYield.min < 0) res.sigYield.min = 0;
    res.sigYield.max = res.sigYield.value * 1.2;
    
    std::cout << "Sig yield After: " << res.sigYield.value << ", " 
              << res.sigYield.min << ", " << res.sigYield.max << std::endl;
    
    // Signal yield limit
    std::cout << "Sig yield lim Before: " << res.sigYieldLim.value << ", " 
              << res.sigYieldLim.min << ", " << res.sigYieldLim.max << std::endl;
    
    res.sigYieldLim.min = res.sigYieldLim.value;
    res.sigYieldLim.max = res.sigYieldLim.value;
    
    std::cout << "Sig lim yield After: " << res.sigYieldLim.value << ", " 
              << res.sigYieldLim.min << ", " << res.sigYieldLim.max << std::endl;
}

void Fitter::updateBKGYield(const std::string& resonance) {
    std::cout << "Update BKGYield limits:" << std::endl;
    auto& res = massDict[resonance];
    
    // Background yield
    std::cout << "BKG yield Before: " << res.bkgYield.value << ", " 
              << res.bkgYield.min << ", " << res.bkgYield.max << std::endl;
    
    res.bkgYield.min = res.bkgYield.value * 0.8;
    if (res.bkgYield.min < 0) res.bkgYield.min = 0;
    res.bkgYield.max = res.bkgYield.value * 1.2;
    
    std::cout << "BKG yield After: " << res.bkgYield.value << ", " 
              << res.bkgYield.min << ", " << res.bkgYield.max << std::endl;
    
    // Background yield limit
    std::cout << "BKG yield lim Before: " << res.bkgYieldLim.value << ", " 
              << res.bkgYieldLim.min << ", " << res.bkgYieldLim.max << std::endl;
    
    res.bkgYieldLim.min = res.bkgYieldLim.value;
    res.bkgYieldLim.max = res.bkgYieldLim.value;
    
    std::cout << "BKG lim yield After: " << res.bkgYieldLim.value << ", " 
              << res.bkgYieldLim.min << ", " << res.bkgYieldLim.max << std::endl;
}

std::string Fitter::fiducialCutString(const std::pair<double, double>& jetPt, 
                                      const std::pair<double, double>& tagMass) {
    std::vector<std::string> outputCuts;
    
    // Add jet pT cuts
    if (jetPt.first > 0) {
        outputCuts.push_back("jetPt > " + std::to_string(jetPt.first));
    }
    if (jetPt.second > 0) {
        outputCuts.push_back("jetPt < " + std::to_string(jetPt.second));
    }
    
    // Add mass cuts
    if (tagMass.first > 0) {
        outputCuts.push_back("tagMass > " + std::to_string(tagMass.first));
    }
    if (tagMass.second > 0) {
        outputCuts.push_back("tagMass < " + std::to_string(tagMass.second));
    }
    
    // Join all cuts with && operator
    std::string finalCutString;
    for (size_t i = 0; i < outputCuts.size(); ++i) {
        finalCutString += outputCuts[i];
        if (i < outputCuts.size() - 1) {
            finalCutString += " && ";
        }
    }
    
    return finalCutString;
}

RooDataSet* Fitter::createDataSet(const std::string& resonance, const std::string& name, 
                                 const std::string& fidCutString, bool isMass, int bin, 
                                 int corrVer, const std::string& fitVarName) {
    
    std::cout << "Unbinned Histogram Fit" << std::endl;
    
    // Check if inTree exists
    if (!inTree) {
        std::cout << "ERROR: No TTree found for analysis." << std::endl;
        return nullptr;
    }
    
    // Get ranges from dictionaries
    auto res = massDict[resonance];
    auto ipchi2_params = ipchi2Dict["Sig" + resonance];
    auto mass_range = res.massRange;
    
    // Declare variables
    RooRealVar* tagMass = new RooRealVar("tagMass", "tagMass", 
                                        mass_range.first, mass_range.second);
    
    // For IP chi2 fits
    RooRealVar* tag_ipchi2 = new RooRealVar("tag_ip_chi2", "tag_ip_chi2", 0, 10000);
    RooRealVar* log_tag_ipchi2 = new RooRealVar("log_tag_ipchi2", "log_tag_ipchi2", 
                                                ipchi2_params.logIpchi2Range.first, 
                                                ipchi2_params.logIpchi2Range.second);
    
    // Other variables
    RooRealVar* pt_jet = new RooRealVar("jetPt", "jetPt", 0, 200);
    RooRealVar* pt_tag = new RooRealVar("tagPt", "tagPt", 0, 200);
    RooRealVar* nConst = new RooRealVar("jetnConst", "jetnConst", 0.0, 300.0);
    RooRealVar* qValue = new RooRealVar("QValue", "QValue", -2, 0.5);
    RooRealVar* dRValue = new RooRealVar("tagJetdR", "tagJetdR", 0.0, 1.0);
    RooRealVar* tagZ = new RooRealVar("tagZ", "tagZ", 0.0, 1.01);
    
    // Create RooArgSet with all variables
    RooArgSet* cutVars = new RooArgSet();
    cutVars->add(*tagMass);
    cutVars->add(*tag_ipchi2);
    cutVars->add(*log_tag_ipchi2);
    cutVars->add(*pt_jet);
    cutVars->add(*pt_tag);
    cutVars->add(*nConst);
    cutVars->add(*qValue);
    cutVars->add(*tagZ);
    cutVars->add(*dRValue);
    
    // Add distance variables
    RooRealVar* distance1 = new RooRealVar("Distance1", "Distance1", -10, 200);
    RooRealVar* distance2 = new RooRealVar("Distance2", "Distance2", -10, 200);
    RooRealVar* distance3 = new RooRealVar("Distance3", "Distance3", -10, 200);
    cutVars->add(*distance1);
    cutVars->add(*distance2);
    cutVars->add(*distance3);
    
    // Add weights
    RooRealVar* effWeight = new RooRealVar("EffWeight", "EffWeight", 0, 25000);
    RooRealVar* effWeight_0 = new RooRealVar("EffWeight_0", "EffWeight_0", 0, 5);
    RooRealVar* effWeight_1 = new RooRealVar("EffWeight_1", "EffWeight_1", 0, 500);
    RooRealVar* effWeight_2 = new RooRealVar("EffWeight_2", "EffWeight_2", 0, 4);
    RooRealVar* effWeight_3 = new RooRealVar("EffWeight_3", "EffWeight_3", -1, 2);
    RooRealVar* effWeight_4 = new RooRealVar("EffWeight_4", "EffWeight_4", 0, 2000);
    RooRealVar* effWeight_Rnd = new RooRealVar("tagnRnd", "tagnRnd", -1, 25);
    cutVars->add(*effWeight);
    cutVars->add(*effWeight_0);
    cutVars->add(*effWeight_1);
    cutVars->add(*effWeight_2);
    cutVars->add(*effWeight_3);
    cutVars->add(*effWeight_4);
    cutVars->add(*effWeight_Rnd);
    
    // Create final cut string
    std::string finalCutString = fidCutString;
    finalCutString += "&& tagPt > 2";
    finalCutString += "&& jetnConst > 1";
    
    if (corrVer > -2) {
        finalCutString += "&& Distance1 < 0.5 && Distance2 < 0.5 && Distance3 < 0.5";
    }
    
    // Create dataset with appropriate weight
    RooDataSet* data = nullptr;
    const char* weightName = nullptr;
    
    if (corrVer == 0) weightName = effWeight_Rnd->GetName();
    else if (corrVer == 1) weightName = effWeight->GetName();
    else if (corrVer == 2) weightName = effWeight_0->GetName();
    else if (corrVer == 3) weightName = effWeight_1->GetName();
    else if (corrVer == 4) weightName = effWeight_2->GetName();
    else if (corrVer == 5) weightName = effWeight_3->GetName();
    else if (corrVer == 6) weightName = effWeight_4->GetName();
    
    if (weightName) {
        data = new RooDataSet(name.c_str(), name.c_str(), inTree, *cutVars, finalCutString.c_str(), weightName);
    }
    else {
        data = new RooDataSet(name.c_str(), name.c_str(), inTree, *cutVars, finalCutString.c_str());
    }
    
    return data;
}

std::tuple<TH1*, std::vector<double>, std::vector<double>> 
Fitter::massFit(const std::string& resonance, RooDataSet* data, const std::string& fitTypeName, 
               int bin, const std::string& zRange, bool splot, TFile* sFile) {
    std::cout << "\n==== Starting massFit with " << fitTypeName << " model for bin " << bin << " ====" << std::endl;
    
    // Suppress RooFit messages except errors
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    
    // Initialize output arrays
    std::vector<double> parameterArr(12, 0.0);
    std::vector<double> parameterErrArr(10, 0.0);
    TH1* histogram = nullptr;
    
    try {
        // Get mass parameters
        auto res = massDict[resonance];
        std::pair<double, double> fullRange = res.massRange;
        std::pair<double, double> signalRange = res.signalRegion;
        std::pair<double, double> SBRange = res.sbRegion;
        
        std::cout << "  Mass range: (" << fullRange.first << ", " << fullRange.second 
                  << "), Signal region: (" << signalRange.first << ", " << signalRange.second << ")" << std::endl;
        std::cout << "  Data entries: " << data->numEntries() << std::endl;
        
        // Create mass variable
        RooRealVar* mass_tag_measured = new RooRealVar("tagMass", "tagMass", fullRange.first, fullRange.second);
        mass_tag_measured->setRange("fullRange", fullRange.first, fullRange.second);
        mass_tag_measured->setRange("signalRange", signalRange.first, signalRange.second);
        mass_tag_measured->setRange("SBleft", fullRange.first, SBRange.first);
        mass_tag_measured->setRange("SBright", SBRange.second, fullRange.second);
        
        // Signal parameters
        RooRealVar* sigma1 = new RooRealVar("sigma1", "sigma1", 
                                          res.sigma1.value, res.sigma1.min, res.sigma1.max);
        RooRealVar* deltasigma = new RooRealVar("deltasigma", "deltasigma", 
                                              res.deltasigma.value, res.deltasigma.min, res.deltasigma.max);
        RooFormulaVar* sigma2 = new RooFormulaVar("sigma2", "sigma2", "sigma1*deltasigma", 
                                                 RooArgList(*sigma1, *deltasigma));
        RooRealVar* mean = new RooRealVar("mean", "mean", 
                                        res.mean.value, res.mean.min, res.mean.max);
        RooRealVar* alpha1 = new RooRealVar("alpha1", "alpha1", 
                                          res.alpha1.value, res.alpha1.min, res.alpha1.max);
        RooFormulaVar* alpha2 = new RooFormulaVar("alpha2", "alpha2", "-1*alpha1", RooArgList(*alpha1));
        RooRealVar* n = new RooRealVar("n", "n", 
                                      res.n.value, res.n.min, res.n.max);
        RooRealVar* dg_frac = new RooRealVar("dg_frac", "dg_frac", 
                                           res.dg_frac.value, res.dg_frac.min, res.dg_frac.max);
        RooRealVar* sig_yield = new RooRealVar("sig_yield", "sig_yield", 
                                             res.sigYield.value, res.sigYield.min, res.sigYield.max);
        
        // Create signal PDF
        RooAbsPdf* sig_pdf = nullptr;
        if (fitTypeName == "DGauss") {
            // Double Gaussian
            RooGaussian* Gauss1_pdf = new RooGaussian("Sig1_pdf", "Sig1_pdf", 
                                                     *mass_tag_measured, *mean, *sigma1);
            RooGaussian* Gauss2_pdf = new RooGaussian("Sig2_pdf", "Sig2_pdf", 
                                                     *mass_tag_measured, *mean, *sigma2);
            sig_pdf = new RooAddPdf("sig_pdf", "Signal", 
                                    RooArgList(*Gauss2_pdf, *Gauss1_pdf), RooArgList(*dg_frac));
            std::cout << "  Using Double Gaussian PDF" << std::endl;
        }
        else if (fitTypeName == "SGauss" || fitTypeName == "noSig") {
            // Single Gaussian
            RooGaussian* Gauss1_pdf = new RooGaussian("Gauss1_pdf", "Gauss1_pdf", 
                                                     *mass_tag_measured, *mean, *sigma1);
            sig_pdf = Gauss1_pdf;
            std::cout << "  Using Single Gaussian PDF" << std::endl;
        }
        
        // Signal extended PDF
        RooExtendPdf* sig_pdf_ext = nullptr;
        if (sig_pdf) {
            sig_pdf_ext = new RooExtendPdf("sig_pdf_ext", "sig_pdf_ext", *sig_pdf, *sig_yield, "fullRange");
        }
        
        // Background parameters and PDF
        RooRealVar* bkg_yield = new RooRealVar("bkg_yield", "bkg_yield", 
                                             res.bkgYield.value, res.bkgYield.min, res.bkgYield.max);
        RooRealVar* poly0 = new RooRealVar("pol0", "pol0", 
                                         res.pol0.value, res.pol0.min, res.pol0.max);
        RooRealVar* poly1 = new RooRealVar("pol1", "pol1", 
                                         res.pol1.value, res.pol1.min, res.pol1.max);
        
        RooPolynomial* bkg_pdf = new RooPolynomial("bkg_pdf", "bkg_pdf", 
                                                 *mass_tag_measured, RooArgList(*poly0, *poly1));
        RooExtendPdf* bkg_pdf_ext = new RooExtendPdf("bkg_pdf_ext", "bkg_pdf_ext", 
                                                    *bkg_pdf, *bkg_yield, "fullRange");
        
        // Build the combined model
        RooAddPdf* extended_pdf = nullptr;
        if (fitTypeName == "noSig") {
            extended_pdf = new RooAddPdf("model", "model", RooArgList(*bkg_pdf_ext));
        }
        else {
            extended_pdf = new RooAddPdf("model", "model", RooArgList(*sig_pdf_ext, *bkg_pdf_ext));
        }
        
        // Fit the data
        RooFitResult* fit_result = nullptr;
        if (fitTypeName == "noSig") {
            std::cout << "  Fitting in sideband regions only" << std::endl;
            fit_result = extended_pdf->fitTo(*data, RooFit::Save(true), 
                                            RooFit::PrintLevel(1), 
                                            RooFit::Range("SBleft,SBright"));
        }
        else {
            std::cout << "  Fitting in full range" << std::endl;
            fit_result = extended_pdf->fitTo(*data, RooFit::Save(true), 
                                            RooFit::PrintLevel(1), 
                                            RooFit::Range("fullRange"));
        }
        
        // Update dictionary if requested
        if (updateStartValues) {
            updateDictionary(extended_pdf, data, fitTypeName);
        }
        
        // Calculate yields in signal region
        RooAbsReal* integral_fullB = bkg_pdf_ext->createIntegral(*mass_tag_measured, 
                                                               RooFit::NormSet(*mass_tag_measured), 
                                                               RooFit::Range("fullRange"));
        RooAbsReal* integral_sigB = bkg_pdf_ext->createIntegral(*mass_tag_measured, 
                                                              RooFit::NormSet(*mass_tag_measured), 
                                                              RooFit::Range("signalRange"));
        
        RooAbsReal* integral_fullS = nullptr;
        RooAbsReal* integral_sigS = nullptr;
        if (sig_pdf_ext) {
            integral_fullS = sig_pdf_ext->createIntegral(*mass_tag_measured, 
                                                       RooFit::NormSet(*mass_tag_measured), 
                                                       RooFit::Range("fullRange"));
            integral_sigS = sig_pdf_ext->createIntegral(*mass_tag_measured, 
                                                      RooFit::NormSet(*mass_tag_measured), 
                                                      RooFit::Range("signalRange"));
        }
        
        // Get parameters
        RooArgSet* parameters = extended_pdf->getParameters(*data);
        
        // Calculate yields
        double fullRangeSYield = sig_yield->getVal();
        double fullRangeBYield = bkg_yield->getVal();
        
        double SfactorBKG = 0.0;
        double SfactorS = 0.0;
        if (integral_fullB && integral_fullB->getVal() > 0) {
            SfactorBKG = fullRangeBYield / integral_fullB->getVal();
        }
        
        if (integral_fullS && integral_fullS->getVal() > 0) {
            SfactorS = fullRangeSYield / integral_fullS->getVal();
        }
        
        double NevtBKG_SignalRange = SfactorBKG * integral_sigB->getVal();
        double NevtS_SignalRange = 0.0;
        if (integral_sigS) {
            NevtS_SignalRange = SfactorS * integral_sigS->getVal();
        }
        
        // Create plot of the fit using the Plotter class
        Plotter plotter(resonance, outfilePath, bin, false, zRange);
        histogram = plotter.individualMassFitPlot(sig_yield, extended_pdf, mass_tag_measured, data, fitTypeName);
        
        // Extract fit parameters for return
        parameterArr[0] = sig_yield->getVal();
        parameterArr[1] = bkg_yield->getVal();
        parameterArr[2] = mean->getVal();
        parameterArr[3] = sigma1->getVal();
        parameterArr[4] = deltasigma->getVal();
        parameterArr[5] = alpha1->getVal();
        parameterArr[6] = n->getVal();
        parameterArr[7] = dg_frac->getVal();
        parameterArr[8] = poly0->getVal();
        parameterArr[9] = poly1->getVal();
        parameterArr[10] = NevtS_SignalRange;
        parameterArr[11] = NevtBKG_SignalRange;
        
        // Extract parameter errors
        if (fitTypeName != "noSig") {
            parameterErrArr[0] = sig_yield->getError();
        }
        parameterErrArr[1] = bkg_yield->getError();
        parameterErrArr[2] = mean->getError();
        parameterErrArr[3] = sigma1->getError();
        parameterErrArr[4] = deltasigma->getError();
        parameterErrArr[5] = alpha1->getError();
        parameterErrArr[6] = n->getError();
        parameterErrArr[7] = dg_frac->getError();
        parameterErrArr[8] = poly0->getError();
        parameterErrArr[9] = poly1->getError();
        
        std::cout << "==== Fit completed successfully ====" << std::endl;
        
        // Clean up
        delete fit_result;
        delete parameters;
        delete integral_fullB;
        delete integral_sigB;
        if (integral_fullS) delete integral_fullS;
        if (integral_sigS) delete integral_sigS;
        
        // Note: We're not deleting the RooFit objects created with new here
        // to avoid segmentation faults due to ROOT's object ownership model.
        // In a real implementation, proper ownership management would be needed.
        
        return std::make_tuple(histogram, parameterArr, parameterErrArr);
    }
    catch (std::exception& e) {
        std::cerr << "Error in massFit: " << e.what() << std::endl;
        return std::make_tuple(nullptr, parameterArr, parameterErrArr);
    }
}
std::tuple<TH1*, std::vector<double>, std::vector<double>> 
Fitter::ipchi2Fit(const std::string& resonance, RooDataSet* data, RooDataSet* background,
                 const std::string& figKey, int bin, const std::string& zRange) {
    std::cout << "\n==== Starting IP chi2 fit with Bukin function for bin " << bin << " ====" << std::endl;
    
    // Suppress RooFit messages
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    
    // Initialize return values
    std::vector<double> parameterArr(12, 0.0);
    std::vector<double> parameterErrArr(12, 0.0);
    TH1* histogram = nullptr;
    
    try {
        // Get dictionaries with parameters
        auto ipchi2_params = ipchi2Dict["Sig" + resonance];
        auto mass_params = massDict[resonance];
        
        // Print dataset info
        std::cout << "  Signal dataset entries: " << data->numEntries() << std::endl;
        if (background) {
            std::cout << "  Background dataset entries: " << background->numEntries() << std::endl;
        }
        
        // Create RooFit variable for log(IP Chi2)
        RooRealVar* log_ipchi2 = new RooRealVar("log_tag_ipchi2", "log(tag_ip_chi2)", 
                                              ipchi2_params.logIpchi2Range.first, 
                                              ipchi2_params.logIpchi2Range.second);
        
        // Create variables for the model
        RooRealVar* sig_yield = new RooRealVar("sig_yield", "sig_yield", 
                                            mass_params.sigYield.value, 
                                            mass_params.sigYield.min, 
                                            mass_params.sigYield.max);
        
        RooRealVar* sig_yieldLim = new RooRealVar("sig_yieldLim", "sig_yieldLim", 
                                               mass_params.sigYieldLim.value, 
                                               mass_params.sigYieldLim.min, 
                                               mass_params.sigYieldLim.max);
        
        RooRealVar* bkg_yieldLim = new RooRealVar("bkg_yieldLim", "bkg_yieldLim", 
                                               mass_params.bkgYieldLim.value, 
                                               mass_params.bkgYieldLim.min, 
                                               mass_params.bkgYieldLim.max);
        
        RooRealVar* prompt_frac = new RooRealVar("prompt_frac", "prompt_frac", 
                                              ipchi2_params.promptFrac.value, 
                                              ipchi2_params.promptFrac.min, 
                                              ipchi2_params.promptFrac.max);
        
        // Create single Bukin for prompt component
        RooRealVar* xp_prompt = new RooRealVar("xp_prompt", "xp_prompt", 
                                            ipchi2_params.xpPrompt.value, 
                                            ipchi2_params.xpPrompt.min, 
                                            ipchi2_params.xpPrompt.max);
        
        RooRealVar* sigma_prompt = new RooRealVar("sigma_prompt", "sigma_prompt", 
                                               ipchi2_params.sigmaPrompt.value, 
                                               ipchi2_params.sigmaPrompt.min, 
                                               ipchi2_params.sigmaPrompt.max);
        
        RooRealVar* xi_prompt = new RooRealVar("xi_prompt", "xi_prompt", 
                                            ipchi2_params.xiPrompt.value, 
                                            ipchi2_params.xiPrompt.min, 
                                            ipchi2_params.xiPrompt.max);
        
        RooRealVar* rho1_prompt = new RooRealVar("rho1_prompt", "rho1_prompt", 
                                              ipchi2_params.rho1Prompt.value, 
                                              ipchi2_params.rho1Prompt.min, 
                                              ipchi2_params.rho1Prompt.max);
        
        RooRealVar* rho2_prompt = new RooRealVar("rho2_prompt", "rho2_prompt", 
                                              ipchi2_params.rho2Prompt.value, 
                                              ipchi2_params.rho2Prompt.min, 
                                              ipchi2_params.rho2Prompt.max);
        
        // Create the prompt Bukin PDF
        RooBukinPdf* prompt_pdf = new RooBukinPdf("prompt_pdf", "prompt_pdf", 
                                               *log_ipchi2, *xp_prompt, *sigma_prompt, *xi_prompt, 
                                               *rho1_prompt, *rho2_prompt);
        
        // Create non-prompt component using a single Bukin
        RooRealVar* xp_nonprompt = new RooRealVar("xp_nonprompt", "xp_nonprompt", 
                                               ipchi2_params.xpNonprompt.value, 
                                               ipchi2_params.xpNonprompt.min, 
                                               ipchi2_params.xpNonprompt.max);
        
        RooRealVar* sigma_nonprompt = new RooRealVar("sigma_nonprompt", "sigma_nonprompt", 
                                                  ipchi2_params.sigmaNonprompt.value, 
                                                  ipchi2_params.sigmaNonprompt.min, 
                                                  ipchi2_params.sigmaNonprompt.max);
        
        RooRealVar* xi_nonprompt = new RooRealVar("xi_nonprompt", "xi_nonprompt", 
                                               ipchi2_params.xiNonprompt.value, 
                                               ipchi2_params.xiNonprompt.min, 
                                               ipchi2_params.xiNonprompt.max);
        
        RooRealVar* rho1_nonprompt = new RooRealVar("rho1_nonprompt", "rho1_nonprompt", 
                                                 ipchi2_params.rho1Nonprompt.value, 
                                                 ipchi2_params.rho1Nonprompt.min, 
                                                 ipchi2_params.rho1Nonprompt.max);
        
        RooRealVar* rho2_nonprompt = new RooRealVar("rho2_nonprompt", "rho2_nonprompt", 
                                                 ipchi2_params.rho2Nonprompt.value, 
                                                 ipchi2_params.rho2Nonprompt.min, 
                                                 ipchi2_params.rho2Nonprompt.max);
        
        RooBukinPdf* nonprompt_pdf = new RooBukinPdf("nonprompt_pdf", "nonprompt_pdf", 
                                                  *log_ipchi2, *xp_nonprompt, *sigma_nonprompt, *xi_nonprompt, 
                                                  *rho1_nonprompt, *rho2_nonprompt);
        
        // Calculate yields for prompt and non-prompt
        RooFormulaVar* prompt_yield = new RooFormulaVar("prompt_yield", "prompt_yield", 
                                                     "sig_yieldLim*prompt_frac", 
                                                     RooArgList(*sig_yieldLim, *prompt_frac));
        
        RooFormulaVar* nonprompt_yield = new RooFormulaVar("nonprompt_yield", "nonprompt_yield", 
                                                        "sig_yieldLim*(1-prompt_frac)", 
                                                        RooArgList(*sig_yieldLim, *prompt_frac));
        
        // Create the total PDF
        RooAbsPdf* total_pdf = nullptr;
        RooAbsPdf* background_pdf = nullptr;
        
        // Create background PDF if provided
        if (background && background->numEntries() > 0) {
            RooRealVar* bkg_param1 = new RooRealVar("bkg_param1", "bkg_param1", 
                                                 ipchi2_params.bkgParam1.value, 
                                                 ipchi2_params.bkgParam1.min, 
                                                 ipchi2_params.bkgParam1.max);
            
            RooRealVar* bkg_param2 = new RooRealVar("bkg_param2", "bkg_param2", 
                                                 ipchi2_params.bkgParam2.value, 
                                                 ipchi2_params.bkgParam2.min, 
                                                 ipchi2_params.bkgParam2.max);
            
            // Use a polynomial for background
            background_pdf = new RooPolynomial("bkg_pdf", "bkg_pdf", 
                                           *log_ipchi2, RooArgList(*bkg_param1, *bkg_param2));
            
            // Combined model with signal and background
            RooAddPdf* signal_model = new RooAddPdf("signal_model", "signal_model", 
                                                 RooArgList(*prompt_pdf, *nonprompt_pdf), 
                                                 RooArgList(*prompt_yield, *nonprompt_yield));
            
            total_pdf = new RooAddPdf("total_pdf", "total_pdf", 
                                     RooArgList(*signal_model, *background_pdf), 
                                     RooArgList(*sig_yieldLim, *bkg_yieldLim));
        } else {
            // Signal-only model
            total_pdf = new RooAddPdf("total_pdf", "total_pdf", 
                                     RooArgList(*prompt_pdf, *nonprompt_pdf), 
                                     RooArgList(*prompt_yield, *nonprompt_yield));
        }
        
        // Perform the fit
        std::cout << "  Performing IP chi2 fit with Bukin function..." << std::endl;
        RooFitResult* result = total_pdf->fitTo(*data, RooFit::Save(true), RooFit::PrintLevel(1));
        
        // Create plot
        std::cout << "  Creating IP chi2 fit plot..." << std::endl;
        Plotter plotter(resonance, outfilePath, bin, false, zRange, figKey);
        histogram = plotter.ipchi2FitPlot(resonance, log_ipchi2, data, total_pdf, 
                                        nonprompt_pdf, prompt_pdf, background_pdf,
                                        prompt_yield, nonprompt_yield);
        
        // Extract fit parameters for simplified Bukin model
        parameterArr[0] = sig_yield->getVal();
        parameterArr[1] = prompt_frac->getVal();
        
        // Prompt Bukin parameters
        parameterArr[2] = xp_prompt->getVal();
        parameterArr[3] = sigma_prompt->getVal();
        parameterArr[4] = xi_prompt->getVal();
        parameterArr[5] = rho1_prompt->getVal();
        parameterArr[6] = rho2_prompt->getVal();
        
        // Non-prompt Bukin parameters
        parameterArr[7] = xp_nonprompt->getVal();
        parameterArr[8] = sigma_nonprompt->getVal();
        parameterArr[9] = xi_nonprompt->getVal();
        parameterArr[10] = rho1_nonprompt->getVal();
        parameterArr[11] = rho2_nonprompt->getVal();
        
        // Extract errors for key parameters
        parameterErrArr[0] = sig_yield->getError();
        parameterErrArr[1] = prompt_frac->getError();
        
        parameterErrArr[2] = xp_prompt->getError();
        parameterErrArr[3] = sigma_prompt->getError();
        parameterErrArr[4] = xi_prompt->getError();
        
        parameterErrArr[7] = xp_nonprompt->getError();
        parameterErrArr[8] = sigma_nonprompt->getError();
        parameterErrArr[9] = xi_nonprompt->getError();
        
        std::cout << "  IP chi2 fit with simplified Bukin function completed successfully" << std::endl;
        std::cout << "  Prompt fraction: " << prompt_frac->getVal() << " ± " << prompt_frac->getError() << std::endl;
        std::cout << "  Prompt asymmetry parameter: xi=" << xi_prompt->getVal() << std::endl;
        std::cout << "  Non-prompt asymmetry parameter: xi=" << xi_nonprompt->getVal() << std::endl;
        
        // Clean up (in a real implementation, more thorough cleanup would be needed)
        delete result;
        
        // Return results
        return std::make_tuple(histogram, parameterArr, parameterErrArr);
        
    } catch (const std::exception& e) {
        std::cerr << "Error in ipchi2Fit with Bukin function: " << e.what() << std::endl;
        return std::make_tuple(nullptr, parameterArr, parameterErrArr);
    }
}