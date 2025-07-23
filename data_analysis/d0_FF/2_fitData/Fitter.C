// Fitter.cpp
#include "Fitter.h"
#include "Plotter.h"
#include <chrono>
#include <thread>
#include "TSystem.h"
#include "RooStats/SPlot.h"

Fitter::Fitter(TTree *tree, const std::string &resonanceType, int numBins, bool zTObservable,
               bool isMCData, const std::string &outputPath, bool update)
    : tfilePV(nullptr), TestFilename(""),
      outfilePath(outputPath), isMC(isMCData), nBins(numBins), isZtObservable(zTObservable),
      resonance(resonanceType), updateStartValues(update), inTree(tree), fInFileHisto(nullptr)
{
    // Initialize dictionary
    initDictionary();

    std::cout << "This is the MC status: " << (isMC ? "true" : "false") << std::endl;

    // Verify the input tree
    if (inTree)
    {
        std::cout << "Using input tree with " << inTree->GetEntries() << " entries" << std::endl;

        // Get the file that owns the tree
        if (inTree->GetCurrentFile())
        {
            std::cout << "Tree source file: " << inTree->GetCurrentFile()->GetName() << std::endl;
        }
    }
    else
    {
        std::cout << "Warning: Input tree is null" << std::endl;
    }
}

Fitter::~Fitter()
{
    // Clean up
    // Note: Don't close the input files here as they might be used elsewhere
}

void Fitter::initDictionary()
{
    // Initialize D0 mass parameters
    MassConfig d0Config;
    d0Config.sigma1 = ParamConfig(0.008, 0.005, 0.012);
    d0Config.deltasigma = ParamConfig(1.5, 1.1, 2.5);
    d0Config.mean = ParamConfig(1.865, 1.860, 1.870);
    d0Config.n = ParamConfig(1, 0.2, 5);
    d0Config.dg_frac = ParamConfig(0.5, 0.0, 0.99999);
    d0Config.pol1 = ParamConfig(-0.5, -10, 10);
    d0Config.pol2 = ParamConfig(58, -2e2, 5e2);
    // d0Config.pol2 = ParamConfig(58, -2e2, 5e2);
    d0Config.massRange = std::make_pair(1.815, 1.925); // 150 MeV range
    // d0Config.massRange = std::make_pair(1.81, 1.935);  // 150 MeV range
    d0Config.sigYield = ParamConfig(100, 0, 5e6);
    d0Config.sigYieldLim = ParamConfig(100, 0, 5e6);
    d0Config.bkgYield = ParamConfig(1000, 0, 4e6);
    d0Config.bkgYieldLim = ParamConfig(1000, 0, 4e6);
    d0Config.signalRegion = std::make_pair(1.845, 1.885); // 40 MeV range
    d0Config.sbRegion = std::make_pair(1.825, 1.910);     // Sideband region
    // d0Config.sbRegion = std::make_pair(1.820, 1.840);     // Sideband region

    massDict["D0"] = d0Config;

    // Initialize D0 IP chi2 parameters
    IPChi2Config d0IPConfig;
    d0IPConfig.logIpchi2Range = std::make_pair(-3, 5);

    // Prompt component
    d0IPConfig.xpPrompt = ParamConfig(0.0, -0.5, 1.0);
    d0IPConfig.sigmaPrompt = ParamConfig(0.7, 0.4, 2.0);
    d0IPConfig.xiPrompt = ParamConfig(0.0, -0.5, 0.5);
    d0IPConfig.rho1Prompt = ParamConfig(-0.08, -0.2, -0.01);
    d0IPConfig.rho2Prompt = ParamConfig(0.01, 0.0001, 0.98);

    // Non-prompt component
    d0IPConfig.xpNonprompt = ParamConfig(1.9, 1.6, 3.0);
    d0IPConfig.sigmaNonprompt = ParamConfig(0.4, 0.3, 0.7);
    d0IPConfig.xiNonprompt = ParamConfig(0.1, 0.05, 0.5);
    d0IPConfig.rho1Nonprompt = ParamConfig(-0.1, -0.95, -0.05);
    d0IPConfig.rho2Nonprompt = ParamConfig(0.2, 0.01, 0.98);

    // Fraction - starting with a more central value to avoid boundary issues
    d0IPConfig.promptFrac = ParamConfig(0.95, 0.7, 0.999);

    // Background
    d0IPConfig.bkgParam1 = ParamConfig(0.5, 0, 1);
    d0IPConfig.bkgParam2 = ParamConfig(0.5, 0, 1);

    ipchi2Dict["SigD0"] = d0IPConfig;
}

void Fitter::updateDictionary(RooAbsPdf *signalPdf, RooAbsData *data, const std::string &fitFunc)
{
    if (!signalPdf || !data || resonance.empty())
        return;

    auto &res = massDict[resonance];

    RooArgSet *paramSet = signalPdf->getParameters(*data);
    if (!paramSet)
        return;

    std::vector<std::string> keyList;

    if (fitFunc == "noSig")
    {
        keyList = {"pol1", "pol2"};
    }
    if (fitFunc == "DGauss")
    {
        keyList = {"mean", "sigma1", "deltasigma", "dg_frac", "sig_yield",
                   "bkg_yield", "pol1", "pol2"};
    }
    else if (fitFunc == "SGauss")
    {
        keyList = {"mean", "sigma1", "sig_yield", "bkg_yield", "pol1", "pol2"};
    }
    else
    {
        return;
    }

    // Update the parameters in our dictionary
    for (const auto &key : keyList)
    {
        RooRealVar *param = dynamic_cast<RooRealVar *>(paramSet->find(key.c_str()));
        if (param)
        {
            double newVal = param->getVal();

            // Update the appropriate parameter
            if (key == "mean")
                res.mean.value = newVal;
            else if (key == "sigma1")
                res.sigma1.value = newVal;
            else if (key == "deltasigma")
                res.deltasigma.value = newVal;
            else if (key == "n")
                res.n.value = newVal;
            else if (key == "dg_frac")
                res.dg_frac.value = newVal;
            else if (key == "sig_yield")
                res.sigYield.value = newVal;
            else if (key == "bkg_yield")
                res.bkgYield.value = newVal;
            else if (key == "pol1")
                res.pol1.value = newVal;
            else if (key == "pol2")
                res.pol2.value = newVal;
        }
    }

    delete paramSet;
}

void Fitter::updateSigYield(const std::string &resonance)
{
    std::cout << "Update SigYield limits:" << std::endl;
    auto &res = massDict[resonance];

    // Signal yield
    std::cout << "Sig yield Before: " << res.sigYield.value << ", "
              << res.sigYield.min << ", " << res.sigYield.max << std::endl;

    res.sigYield.min = res.sigYield.value * 0.8;
    if (res.sigYield.min < 0)
        res.sigYield.min = 0;
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

void Fitter::updateBKGYield(const std::string &resonance)
{
    std::cout << "Update BKGYield limits:" << std::endl;
    auto &res = massDict[resonance];

    // Background yield
    std::cout << "BKG yield Before: " << res.bkgYield.value << ", "
              << res.bkgYield.min << ", " << res.bkgYield.max << std::endl;

    res.bkgYield.min = res.bkgYield.value * 0.8;
    if (res.bkgYield.min < 0)
        res.bkgYield.min = 0;
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

std::string Fitter::fiducialCutString(const std::pair<double, double> &jetPt,
                                      const std::pair<double, double> &tagMass)
{
    std::vector<std::string> outputCuts;

    // Add jet pT cuts
    if (jetPt.first > 0)
    {
        outputCuts.push_back("jetPt > " + std::to_string(jetPt.first));
    }
    if (jetPt.second > 0)
    {
        outputCuts.push_back("jetPt < " + std::to_string(jetPt.second));
    }

    // Add mass cuts
    if (tagMass.first > 0)
    {
        outputCuts.push_back("tagMass > " + std::to_string(tagMass.first));
    }
    if (tagMass.second > 0)
    {
        outputCuts.push_back("tagMass < " + std::to_string(tagMass.second));
    }

    // Join all cuts with && operator
    std::string finalCutString;
    for (size_t i = 0; i < outputCuts.size(); ++i)
    {
        finalCutString += outputCuts[i];
        if (i < outputCuts.size() - 1)
        {
            finalCutString += " && ";
        }
    }

    return finalCutString;
}

RooDataSet *Fitter::createDataSet(const std::string &resonance, const std::string &name,
                                  const std::string &fidCutString, bool isMass,
                                  int corrVer, const std::string &fitVarName)
{

    std::cout << "Unbinned Histogram Fit" << std::endl;

    // Check if inTree exists
    if (!inTree)
    {
        std::cout << "ERROR: No TTree found for analysis." << std::endl;
        return nullptr;
    }

    // Get ranges from dictionaries
    auto res = massDict[resonance];
    auto ipchi2_params = ipchi2Dict["Sig" + resonance];
    auto mass_range = res.massRange;

    // Declare variables
    RooRealVar *tagMass = new RooRealVar("tagMass", "tagMass",
                                         mass_range.first, mass_range.second);

    // For IP chi2 fits
    RooRealVar *tag_ipchi2 = new RooRealVar("tag_ip_chi2", "tag_ip_chi2", 0, 10000);
    RooRealVar *log_tag_ipchi2 = new RooRealVar("log_tag_ipchi2", "log_tag_ipchi2",
                                                ipchi2_params.logIpchi2Range.first,
                                                ipchi2_params.logIpchi2Range.second);

    // Other variables
    RooRealVar *pt_jet = new RooRealVar("jetPt", "jetPt", 0, 200);
    RooRealVar *pt_tag = new RooRealVar("tagPt", "tagPt", 0, 200);
    RooRealVar *nConst = new RooRealVar("jetnConst", "jetnConst", 0.0, 300.0);
    // RooRealVar* qValue = new RooRealVar("QValue", "QValue", -2, 0.5);
    RooRealVar *dRValue = new RooRealVar("tagJetdR", "tagJetdR", 0.0, 1.0);
    RooRealVar *tagZ = new RooRealVar("tagZ", "tagZ", 0.0, 1.01);
    RooRealVar *tagY = new RooRealVar("tagY", "tagY", 0.0, 5.01);

    // Create RooArgSet with all variables
    RooArgSet *cutVars = new RooArgSet();
    cutVars->add(*tagMass);
    cutVars->add(*tag_ipchi2);
    cutVars->add(*log_tag_ipchi2);
    cutVars->add(*pt_jet);
    cutVars->add(*pt_tag);
    cutVars->add(*nConst);
    // cutVars->add(*qValue);
    cutVars->add(*tagZ);
    cutVars->add(*tagY);
    cutVars->add(*dRValue);

    // Add distance variables
    RooRealVar *distance1 = new RooRealVar("Distance1", "Distance1", -10, 200);
    // RooRealVar* distance2 = new RooRealVar("Distance2", "Distance2", -10, 200);
    // RooRealVar* distance3 = new RooRealVar("Distance3", "Distance3", -10, 200);
    cutVars->add(*distance1);
    // cutVars->add(*distance2);
    // cutVars->add(*distance3);

    std::cout << "Using correction version: " << corrVer << std::endl;

    // Add weights
    RooRealVar *kaon_efficiency = new RooRealVar("kaon_efficiency", "kaon_efficiency", 0, 1);
    RooRealVar *pion_efficiency = new RooRealVar("pion_efficiency", "pion_efficiency", 0, 1);
    RooRealVar *combined_efficiency = new RooRealVar("combined_efficiency", "combined_efficiency", 0, 1);
    // RooRealVar* effWeight = new RooRealVar("EffWeight", "EffWeight", 0, 25000);
    // RooRealVar* effWeight_0 = new RooRealVar("EffWeight_0", "EffWeight_0", 0, 5);
    // RooRealVar* effWeight_1 = new RooRealVar("EffWeight_1", "EffWeight_1", 0, 500);
    // RooRealVar* effWeight_2 = new RooRealVar("EffWeight_2", "EffWeight_2", 0, 4);
    // RooRealVar* effWeight_3 = new RooRealVar("EffWeight_3", "EffWeight_3", -1, 2);
    // RooRealVar* effWeight_4 = new RooRealVar("EffWeight_4", "EffWeight_4", 0, 2000);
    // RooRealVar* effWeight_Rnd = new RooRealVar("tagnRnd", "tagnRnd", -1, 25);
    cutVars->add(*kaon_efficiency);
    cutVars->add(*pion_efficiency);
    cutVars->add(*combined_efficiency);
    // cutVars->add(*effWeight);
    // cutVars->add(*effWeight_0);
    // cutVars->add(*effWeight_1);
    // cutVars->add(*effWeight_2);
    // cutVars->add(*effWeight_3);
    // cutVars->add(*effWeight_4);
    // cutVars->add(*effWeight_Rnd);

    // Create final cut string
    std::string finalCutString = fidCutString;
    finalCutString += "&& tagPt > 2";
    finalCutString += "&& jetnConst > 1";

    if (corrVer > -2)
    {
        finalCutString += "&& Distance1 < 0.5";
        // finalCutString += "&& Distance1 < 0.5 && Distance2 < 0.5 && Distance3 < 0.5";
    }

    // Create dataset with appropriate weight
    RooDataSet *data = nullptr;
    const char *weightName = nullptr;

    std::cout << "Debug: Creating dataset with corrVer = " << corrVer << std::endl;

    if (corrVer == 1)
    {
        weightName = "kaon_efficiency";
        std::cout << "Debug: Using kaon_efficiency as weight" << std::endl;
    }
    else if (corrVer == 2)
    {
        weightName = "pion_efficiency";
        std::cout << "Debug: Using pion_efficiency as weight" << std::endl;
    }
    else if (corrVer == 3)
    {
        weightName = "combined_efficiency";
        std::cout << "Debug: Using combined_efficiency as weight" << std::endl;
    }
    // if (corrVer == 0) weightName = effWeight_Rnd->GetName();
    // else if (corrVer == 1) weightName = effWeight->GetName();
    // else if (corrVer == 2) weightName = effWeight_0->GetName();
    // else if (corrVer == 3) weightName = effWeight_1->GetName();
    // else if (corrVer == 4) weightName = effWeight_2->GetName();
    // else if (corrVer == 5) weightName = effWeight_3->GetName();
    // else if (corrVer == 6) weightName = effWeight_4->GetName();

    if (weightName)
    {
        std::cout << "Debug: Creating weighted dataset with weight: " << weightName << std::endl;
        data = new RooDataSet(name.c_str(), name.c_str(), inTree, *cutVars, finalCutString.c_str(), weightName);
    }
    else
    {
        std::cout << "Debug: Creating unweighted dataset" << std::endl;
        data = new RooDataSet(name.c_str(), name.c_str(), inTree, *cutVars, finalCutString.c_str());
    }

    return data;
}

std::tuple<TH1 *, std::vector<double>, std::vector<double>>
Fitter::massFit(const std::string &resonance, RooDataSet *data, const std::string &fitTypeName,
                int bin, const std::string &zRange, bool splot, TFile *sFile)
{
    std::cout << "\n==== Starting massFit with " << fitTypeName << " model for bin " << bin << " ====" << std::endl;

    // Suppress RooFit messages except errors
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

    // Initialize output arrays
    std::vector<double> parameterArr(12, 0.0);
    std::vector<double> parameterErrArr(10, 0.0);
    TH1 *histogram = nullptr;

    try
    {
        // Get mass parameters
        auto res = massDict[resonance];
        std::pair<double, double> fullRange = res.massRange;
        std::pair<double, double> signalRange = res.signalRegion;
        std::pair<double, double> SBRange = res.sbRegion;

        std::cout << "  Mass range: (" << fullRange.first << ", " << fullRange.second
                  << "), Signal region: (" << signalRange.first << ", " << signalRange.second << ")" << std::endl;
        std::cout << "  Data entries: " << data->numEntries() << std::endl;

        // Create mass variable
        RooRealVar *mass_tag_measured = new RooRealVar("tagMass", "tagMass", fullRange.first, fullRange.second);
        mass_tag_measured->setRange("fullRange", fullRange.first, fullRange.second);
        mass_tag_measured->setRange("signalRange", signalRange.first, signalRange.second);
        mass_tag_measured->setRange("SBleft", fullRange.first, SBRange.first);
        mass_tag_measured->setRange("SBright", SBRange.second, fullRange.second);

        // Signal parameters
        RooRealVar *sigma1 = new RooRealVar("sigma1", "sigma1",
                                            res.sigma1.value, res.sigma1.min, res.sigma1.max);
        RooRealVar *deltasigma = new RooRealVar("deltasigma", "deltasigma",
                                                res.deltasigma.value, res.deltasigma.min, res.deltasigma.max);
        RooFormulaVar *sigma2 = new RooFormulaVar("sigma2", "sigma2", "sigma1*deltasigma",
                                                  RooArgList(*sigma1, *deltasigma));
        RooRealVar *mean = new RooRealVar("mean", "mean",
                                          res.mean.value, res.mean.min, res.mean.max);
        RooRealVar *n = new RooRealVar("n", "n",
                                       res.n.value, res.n.min, res.n.max);
        RooRealVar *dg_frac = new RooRealVar("dg_frac", "dg_frac",
                                             res.dg_frac.value, res.dg_frac.min, res.dg_frac.max);
        RooRealVar *sig_yield = new RooRealVar("sig_yield", "sig_yield",
                                               res.sigYield.value, res.sigYield.min, res.sigYield.max);

        // Create signal PDF
        RooAbsPdf *sig_pdf = nullptr;
        if (fitTypeName == "DGauss")
        {
            // Double Gaussian
            RooGaussian *Gauss1_pdf = new RooGaussian("Sig1_pdf", "Sig1_pdf",
                                                      *mass_tag_measured, *mean, *sigma1);
            RooGaussian *Gauss2_pdf = new RooGaussian("Sig2_pdf", "Sig2_pdf",
                                                      *mass_tag_measured, *mean, *sigma2);
            sig_pdf = new RooAddPdf("sig_pdf", "Signal",
                                    RooArgList(*Gauss2_pdf, *Gauss1_pdf), RooArgList(*dg_frac));
            std::cout << "  Using Double Gaussian PDF" << std::endl;
        }
        else if (fitTypeName == "SGauss" || fitTypeName == "noSig")
        {
            // Single Gaussian
            RooGaussian *Gauss1_pdf = new RooGaussian("Gauss1_pdf", "Gauss1_pdf",
                                                      *mass_tag_measured, *mean, *sigma1);
            sig_pdf = Gauss1_pdf;
            std::cout << "  Using Single Gaussian PDF" << std::endl;
        }

        // Signal extended PDF
        RooExtendPdf *sig_pdf_ext = nullptr;
        if (sig_pdf)
        {
            sig_pdf_ext = new RooExtendPdf("sig_pdf_ext", "sig_pdf_ext", *sig_pdf, *sig_yield, "fullRange");
        }

        // Background parameters and PDF
        RooRealVar *bkg_yield = new RooRealVar("bkg_yield", "bkg_yield",
                                               res.bkgYield.value, res.bkgYield.min, res.bkgYield.max);
        RooRealVar *poly0 = new RooRealVar("pol0", "pol0",
                                           res.pol1.value, res.pol1.min, res.pol1.max);
        RooRealVar *poly1 = new RooRealVar("pol1", "pol1",
                                           res.pol2.value, res.pol2.min, res.pol2.max);

        RooPolynomial *bkg_pdf = new RooPolynomial("bkg_pdf", "bkg_pdf",
                                                   //  *mass_tag_measured, RooArgList(*poly0));
                                                   *mass_tag_measured, RooArgList(*poly0, *poly1));
        RooExtendPdf *bkg_pdf_ext = new RooExtendPdf("bkg_pdf_ext", "bkg_pdf_ext",
                                                     *bkg_pdf, *bkg_yield, "fullRange");

        // Build the combined model
        RooAddPdf *extended_pdf = nullptr;
        if (fitTypeName == "noSig")
        {
            extended_pdf = new RooAddPdf("model", "model", RooArgList(*bkg_pdf_ext));
        }
        else
        {
            extended_pdf = new RooAddPdf("model", "model", RooArgList(*sig_pdf_ext, *bkg_pdf_ext));
        }

        // Fit the data
        RooFitResult *fit_result = nullptr;
        if (fitTypeName == "noSig")
        {
            std::cout << "  Fitting in sideband regions only" << std::endl;
            fit_result = extended_pdf->fitTo(*data, RooFit::Save(true),
                                             RooFit::PrintLevel(0),
                                             RooFit::Range("SBleft,SBright"));
        }
        else
        {
            std::cout << "  Fitting in full range" << std::endl;
            fit_result = extended_pdf->fitTo(*data, RooFit::Save(true),
                                             RooFit::PrintLevel(0),
                                             RooFit::Range("fullRange"));
        }

        // Update dictionary if requested
        if (updateStartValues)
        {
            updateDictionary(extended_pdf, data, fitTypeName);
        }

        // Calculate yields in signal region
        RooAbsReal *integral_fullB = bkg_pdf_ext->createIntegral(*mass_tag_measured,
                                                                 RooFit::NormSet(*mass_tag_measured),
                                                                 RooFit::Range("fullRange"));
        RooAbsReal *integral_sigB = bkg_pdf_ext->createIntegral(*mass_tag_measured,
                                                                RooFit::NormSet(*mass_tag_measured),
                                                                RooFit::Range("signalRange"));

        RooAbsReal *integral_fullS = nullptr;
        RooAbsReal *integral_sigS = nullptr;
        if (sig_pdf_ext)
        {
            integral_fullS = sig_pdf_ext->createIntegral(*mass_tag_measured,
                                                         RooFit::NormSet(*mass_tag_measured),
                                                         RooFit::Range("fullRange"));
            integral_sigS = sig_pdf_ext->createIntegral(*mass_tag_measured,
                                                        RooFit::NormSet(*mass_tag_measured),
                                                        RooFit::Range("signalRange"));
        }

        // Get parameters
        RooArgSet *parameters = extended_pdf->getParameters(*data);

        // Calculate yields
        double fullRangeSYield = sig_yield->getVal();
        double fullRangeBYield = bkg_yield->getVal();

        double SfactorBKG = 0.0;
        double SfactorS = 0.0;
        if (integral_fullB && integral_fullB->getVal() > 0)
        {
            SfactorBKG = fullRangeBYield / integral_fullB->getVal();
        }

        if (integral_fullS && integral_fullS->getVal() > 0)
        {
            SfactorS = fullRangeSYield / integral_fullS->getVal();
        }

        double NevtBKG_SignalRange = SfactorBKG * integral_sigB->getVal();
        double NevtS_SignalRange = 0.0;
        if (integral_sigS)
        {
            NevtS_SignalRange = SfactorS * integral_sigS->getVal();
        }

        // Create plot of the fit using the Plotter class
        Plotter plotter(resonance, outfilePath, bin, false, zRange);
        histogram = plotter.individualMassFitPlot(sig_yield, extended_pdf, mass_tag_measured, data, fitTypeName, isZtObservable);

        // Perform sPlot analysis if requested
        RooDataSet *splotData = nullptr;
        if (splot && fitTypeName != "noSig" && sFile)
        {
            std::cout << "  Creating sPlot weights..." << std::endl;

            // Create sPlot
            RooArgList yieldsList;
            yieldsList.add(*sig_yield);
            yieldsList.add(*bkg_yield);

            // Create sPlot object
            RooStats::SPlot *splotObj = new RooStats::SPlot("splotObj", "splotObj",
                                                            *data, extended_pdf, yieldsList);

            // Get the sPlot dataset with weights
            splotData = new RooDataSet(*data);

            // Add sPlot weights to dataset
            std::cout << "  Adding sPlot weights to dataset..." << std::endl;
            for (int i = 0; i < data->numEntries(); i++)
            {
                const RooArgSet *row = data->get(i);

                // Get sPlot weights using the actual variable names
                double sigWeight = splotObj->GetSWeight(i, sig_yield->GetName());
                double bkgWeight = splotObj->GetSWeight(i, bkg_yield->GetName());

                // Add weights as variables to the dataset
                RooRealVar *sigWeightVar = new RooRealVar("sig_sWeight", "sig_sWeight", sigWeight);
                RooRealVar *bkgWeightVar = new RooRealVar("bkg_sWeight", "bkg_sWeight", bkgWeight);

                // Note: In a full implementation, you would add these weights properly to the dataset
                // For now, we'll save them to the file
            }

            // Save sPlot results to file
            std::cout << "  Saving sPlot results to file..." << std::endl;
            sFile->cd();

            // Create tree to save sPlot weights
            TTree *splotTree = new TTree(("splotTree_bin" + std::to_string(bin)).c_str(),
                                         "sPlot weights");

            double mass_val, sig_weight, bkg_weight;
            splotTree->Branch("mass", &mass_val, "mass/D");
            splotTree->Branch("sig_sWeight", &sig_weight, "sig_sWeight/D");
            splotTree->Branch("bkg_sWeight", &bkg_weight, "bkg_sWeight/D");

            // Fill the tree with sPlot weights
            for (int i = 0; i < data->numEntries(); i++)
            {
                const RooArgSet *row = data->get(i);
                RooRealVar *massVar = (RooRealVar *)row->find("tagMass");

                if (massVar)
                {
                    mass_val = massVar->getVal();
                    sig_weight = splotObj->GetSWeight(i, sig_yield->GetName());
                    bkg_weight = splotObj->GetSWeight(i, bkg_yield->GetName());

                    splotTree->Fill();
                }
            }

            splotTree->Write();

            // Flush the file to ensure data is written to disk
            sFile->Flush();

            // Force synchronization by writing the tree again if needed
            // and ensuring it's properly on disk
            sFile->cd();
            gSystem->ProcessEvents(); // Process any pending ROOT events

            // Create diagnostic plots
            std::cout << "  Creating sPlot diagnostic plots..." << std::endl;

            // Plot sPlot weights vs mass
            TCanvas *splotCanvas = new TCanvas("splotCanvas", "sPlot Weights", 800, 600);
            splotCanvas->Divide(2, 1);

            splotCanvas->cd(1);
            TH2F *sigWeightHist = new TH2F("sigWeightHist", "Signal sPlot Weights vs Mass",
                                           50, fullRange.first, fullRange.second,
                                           50, -5, 5);
            sigWeightHist->SetXTitle("Mass [GeV]");
            sigWeightHist->SetYTitle("Signal sPlot Weight");

            splotCanvas->cd(2);
            TH2F *bkgWeightHist = new TH2F("bkgWeightHist", "Background sPlot Weights vs Mass",
                                           50, fullRange.first, fullRange.second,
                                           50, -5, 5);
            bkgWeightHist->SetXTitle("Mass [GeV]");
            bkgWeightHist->SetYTitle("Background sPlot Weight");

            // Fill histograms
            for (int i = 0; i < data->numEntries(); i++)
            {
                const RooArgSet *row = data->get(i);
                RooRealVar *massVar = (RooRealVar *)row->find("tagMass");

                if (massVar)
                {
                    double mass = massVar->getVal();
                    double sigW = splotObj->GetSWeight(i, sig_yield->GetName());
                    double bkgW = splotObj->GetSWeight(i, bkg_yield->GetName());

                    sigWeightHist->Fill(mass, sigW);
                    bkgWeightHist->Fill(mass, bkgW);
                }
            }

            splotCanvas->cd(1);
            sigWeightHist->Draw("colz");
            splotCanvas->cd(2);
            bkgWeightHist->Draw("colz");

            splotCanvas->Write();

            // Flush the file again to ensure canvas is written
            sFile->Flush();

            // Cleanup sPlot objects
            delete splotObj;
            delete splotCanvas;
            delete sigWeightHist;
            delete bkgWeightHist;

            std::cout << "  sPlot analysis completed successfully" << std::endl;
        }

        // Extract fit parameters for return
        parameterArr[0] = sig_yield->getVal();
        parameterArr[1] = bkg_yield->getVal();
        parameterArr[2] = mean->getVal();
        parameterArr[3] = sigma1->getVal();
        parameterArr[4] = deltasigma->getVal();
        parameterArr[6] = n->getVal();
        parameterArr[7] = dg_frac->getVal();
        parameterArr[8] = poly0->getVal();
        parameterArr[9] = poly1->getVal();
        parameterArr[10] = NevtS_SignalRange;
        parameterArr[11] = NevtBKG_SignalRange;

        // Extract parameter errors
        if (fitTypeName != "noSig")
        {
            parameterErrArr[0] = sig_yield->getError();
        }
        parameterErrArr[1] = bkg_yield->getError();
        parameterErrArr[2] = mean->getError();
        parameterErrArr[3] = sigma1->getError();
        parameterErrArr[4] = deltasigma->getError();
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
        if (integral_fullS)
            delete integral_fullS;
        if (integral_sigS)
            delete integral_sigS;
        if (splotData)
            delete splotData;

        // Note: We're not deleting the RooFit objects created with new here
        // to avoid segmentation faults due to ROOT's object ownership model.
        // In a real implementation, proper ownership management would be needed.

        return std::make_tuple(histogram, parameterArr, parameterErrArr);
    }
    catch (std::exception &e)
    {
        std::cerr << "Error in massFit: " << e.what() << std::endl;
        return std::make_tuple(nullptr, parameterArr, parameterErrArr);
    }
}

// IP chi2 fit method that returns yield variables for sPlot
std::tuple<TH1 *, std::vector<double>, std::vector<double>, RooAbsPdf *, RooRealVar *, RooRealVar *>
Fitter::ipchi2FitWithYields(const std::string &resonance, RooDataSet *data, RooDataSet *background,
                            const std::string &figKey, int bin, const std::string &zRange,
                            bool enableSPlot, TFile *splotFile)
{
    std::cout << "\n==== Starting IP chi2 fit with yields for sPlot (bin " << bin << ") ====" << std::endl;

    // Suppress RooFit messages
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

    // Initialize return values
    std::vector<double> parameterArr(12, 0.0);
    std::vector<double> parameterErrArr(12, 0.0);
    TH1 *histogram = nullptr;
    RooAbsPdf *total_pdf = nullptr;
    // RooRealVar* promptYieldVar = nullptr;
    // RooRealVar* nonpromptYieldVar = nullptr;

    try
    {
        // Get dictionaries with parameters
        auto ipchi2_params = ipchi2Dict["Sig" + resonance];
        auto mass_params = massDict[resonance];

        // Print dataset info
        std::cout << "  Signal dataset entries: " << data->numEntries() << std::endl;
        if (background)
        {
            std::cout << "  Background dataset entries: " << background->numEntries() << std::endl;
        }

        // Create RooFit variable for log(IP Chi2)
        RooRealVar *log_ipchi2 = new RooRealVar("log_tag_ipchi2", "log(tag_ip_chi2)",
                                                ipchi2_params.logIpchi2Range.first,
                                                ipchi2_params.logIpchi2Range.second);

        // Create variables for the model
        RooRealVar *sig_yield = new RooRealVar("sig_yield", "sig_yield",
                                               mass_params.sigYield.value,
                                               mass_params.sigYield.min,
                                               mass_params.sigYield.max);

        RooRealVar *sig_yieldLim = new RooRealVar("sig_yieldLim", "sig_yieldLim",
                                                  mass_params.sigYieldLim.value,
                                                  mass_params.sigYieldLim.min,
                                                  mass_params.sigYieldLim.max);

        RooRealVar *bkg_yieldLim = new RooRealVar("bkg_yieldLim", "bkg_yieldLim",
                                                  mass_params.bkgYieldLim.value,
                                                  mass_params.bkgYieldLim.min,
                                                  mass_params.bkgYieldLim.max);

        RooRealVar *prompt_frac = new RooRealVar("prompt_frac", "prompt_frac",
                                                 ipchi2_params.promptFrac.value,
                                                 ipchi2_params.promptFrac.min,
                                                 ipchi2_params.promptFrac.max);

        // Debug: Print the prompt_frac value that was used
        std::cout << "  Created prompt_frac with value: " << prompt_frac->getVal() << std::endl;

        // Create single Bukin for prompt component
        RooRealVar *xp_prompt = new RooRealVar("xp_prompt", "xp_prompt",
                                               ipchi2_params.xpPrompt.value,
                                               ipchi2_params.xpPrompt.min,
                                               ipchi2_params.xpPrompt.max);

        RooRealVar *sigma_prompt = new RooRealVar("sigma_prompt", "sigma_prompt",
                                                  ipchi2_params.sigmaPrompt.value,
                                                  ipchi2_params.sigmaPrompt.min,
                                                  ipchi2_params.sigmaPrompt.max);

        RooRealVar *xi_prompt = new RooRealVar("xi_prompt", "xi_prompt",
                                               ipchi2_params.xiPrompt.value,
                                               ipchi2_params.xiPrompt.min,
                                               ipchi2_params.xiPrompt.max);

        RooRealVar *rho1_prompt = new RooRealVar("rho1_prompt", "rho1_prompt",
                                                 ipchi2_params.rho1Prompt.value,
                                                 ipchi2_params.rho1Prompt.min,
                                                 ipchi2_params.rho1Prompt.max);

        RooRealVar *rho2_prompt = new RooRealVar("rho2_prompt", "rho2_prompt",
                                                 ipchi2_params.rho2Prompt.value,
                                                 ipchi2_params.rho2Prompt.min,
                                                 ipchi2_params.rho2Prompt.max);

        // Create the prompt Bukin PDF
        RooBukinPdf *prompt_pdf = new RooBukinPdf("prompt_pdf", "prompt_pdf",
                                                  *log_ipchi2, *xp_prompt, *sigma_prompt, *xi_prompt,
                                                  *rho1_prompt, *rho2_prompt);

        // Create non-prompt component using a single Bukin
        RooRealVar *xp_nonprompt = new RooRealVar("xp_nonprompt", "xp_nonprompt",
                                                  ipchi2_params.xpNonprompt.value,
                                                  ipchi2_params.xpNonprompt.min,
                                                  ipchi2_params.xpNonprompt.max);

        RooRealVar *sigma_nonprompt = new RooRealVar("sigma_nonprompt", "sigma_nonprompt",
                                                     ipchi2_params.sigmaNonprompt.value,
                                                     ipchi2_params.sigmaNonprompt.min,
                                                     ipchi2_params.sigmaNonprompt.max);

        RooRealVar *xi_nonprompt = new RooRealVar("xi_nonprompt", "xi_nonprompt",
                                                  ipchi2_params.xiNonprompt.value,
                                                  ipchi2_params.xiNonprompt.min,
                                                  ipchi2_params.xiNonprompt.max);

        RooRealVar *rho1_nonprompt = new RooRealVar("rho1_nonprompt", "rho1_nonprompt",
                                                    ipchi2_params.rho1Nonprompt.value,
                                                    ipchi2_params.rho1Nonprompt.min,
                                                    ipchi2_params.rho1Nonprompt.max);

        RooRealVar *rho2_nonprompt = new RooRealVar("rho2_nonprompt", "rho2_nonprompt",
                                                    ipchi2_params.rho2Nonprompt.value,
                                                    ipchi2_params.rho2Nonprompt.min,
                                                    ipchi2_params.rho2Nonprompt.max);

        RooBukinPdf *nonprompt_pdf = new RooBukinPdf("nonprompt_pdf", "nonprompt_pdf",
                                                     *log_ipchi2, *xp_nonprompt, *sigma_nonprompt, *xi_nonprompt,
                                                     *rho1_nonprompt, *rho2_nonprompt);

        // Create yield variables for sPlot (these will be used in PDF construction)
        RooFormulaVar *prompt_yield = new RooFormulaVar("prompt_yield", "prompt_yield",
                                                        "sig_yieldLim*prompt_frac",
                                                        RooArgList(*sig_yieldLim, *prompt_frac));

        RooFormulaVar *nonprompt_yield = new RooFormulaVar("nonprompt_yield", "nonprompt_yield",
                                                           "sig_yieldLim*(1-prompt_frac)",
                                                           RooArgList(*sig_yieldLim, *prompt_frac));

        // Create extended PDFs for each component
        RooExtendPdf *prompt_pdf_ext = new RooExtendPdf("prompt_pdf_ext", "prompt_pdf_ext",
                                                        *prompt_pdf, *prompt_yield);

        RooExtendPdf *nonprompt_pdf_ext = new RooExtendPdf("nonprompt_pdf_ext", "nonprompt_pdf_ext",
                                                           *nonprompt_pdf, *nonprompt_yield);

        // Create the total PDF using the extended PDFs
        total_pdf = new RooAddPdf("ipchi2_model", "ipchi2_model",
                                  RooArgList(*prompt_pdf_ext, *nonprompt_pdf_ext));
        // Perform the fit
        std::cout << "  Performing IP chi2 fit with Bukin function..." << std::endl;

        // Check if prompt_frac is constant (fixed)
        if (prompt_frac->isConstant())
        {
            std::cout << "  WARNING: prompt_frac is constant! Setting it to variable..." << std::endl;
            prompt_frac->setConstant(false);
        }

        std::cout << "  Initial prompt_frac: " << prompt_frac->getVal() << " (range: "
                  << prompt_frac->getMin() << " - " << prompt_frac->getMax() << ")" << std::endl;
        std::cout << "  Initial sig_yieldLim: " << sig_yieldLim->getVal() << std::endl;

        RooFitResult *result = total_pdf->fitTo(*data, RooFit::Save(true),
                                                RooFit::PrintLevel(0),
                                                RooFit::SumW2Error(true),
                                                RooFit::Strategy(2),
                                                RooFit::Minos(false),
                                                RooFit::Hesse(true));

        std::cout << "  Fitted prompt_frac: " << prompt_frac->getVal() << " ± " << prompt_frac->getError() << std::endl;
        std::cout << "  Fitted sig_yieldLim: " << sig_yieldLim->getVal() << std::endl;
        std::cout << "  Fitted prompt_yield: " << prompt_yield->getVal() << std::endl;
        std::cout << "  Fitted nonprompt_yield: " << nonprompt_yield->getVal() << std::endl;

        // Create plot
        std::cout << "  Creating IP chi2 fit plot..." << std::endl;
        Plotter plotter(resonance, outfilePath, bin, false, zRange, figKey);
        histogram = plotter.ipchi2FitPlot(resonance, log_ipchi2, data, total_pdf,
                                          nonprompt_pdf, prompt_pdf, nullptr,
                                          prompt_yield, nonprompt_yield);

        // Extract fit parameters
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

        // Extract errors
        parameterErrArr[0] = sig_yield->getError();
        parameterErrArr[1] = prompt_frac->getError();

        parameterErrArr[2] = xp_prompt->getError();
        parameterErrArr[3] = sigma_prompt->getError();
        parameterErrArr[4] = xi_prompt->getError();

        parameterErrArr[7] = xp_nonprompt->getError();
        parameterErrArr[8] = sigma_nonprompt->getError();
        parameterErrArr[9] = xi_nonprompt->getError();

        std::cout << "  IP chi2 fit with yields completed successfully" << std::endl;
        std::cout << "  Prompt fraction: " << prompt_frac->getVal() << " ± " << prompt_frac->getError() << std::endl;

        // Create sPlot weights if requested
        if (enableSPlot && splotFile)
        {
            std::cout << "  Creating IP chi2 sPlot weights..." << std::endl;

            // Create a new model for SPlot using the component PDFs directly
            RooAddPdf *splot_model = new RooAddPdf("splot_model", "Model for SPlot",
                                                   RooArgList(*prompt_pdf, *nonprompt_pdf),
                                                   RooArgList(*sig_yieldLim, *bkg_yieldLim));

            // Save IP chi2 sPlot weights to ROOT file
            std::string ipSplotTreeName = "ipSplotTree_bin" + std::to_string(bin);

            std::cout << "Saving IP chi2 sPlot weights to ROOT file for bin " << bin << std::endl;


            // Check if the splot_model has parameters
            RooArgSet *params = splot_model->getParameters(*data);
            if (!params || params->getSize() == 0)
            {
                std::cerr << "Error: Model has no parameters or failed to get parameters" << std::endl;
            }

            std::cout << "  Model has " << params->getSize() << " parameters" << std::endl;

            // Get the yield variables from the splot_model coefficients
            // For a RooAddPdf, the coefficients should be the yield variables
            RooAddPdf *addPdf = dynamic_cast<RooAddPdf *>(splot_model);
            if (!addPdf)
            {
                std::cerr << "Error: Model is not a RooAddPdf" << std::endl;
                delete params;
            }

            RooArgList coefList = addPdf->coefList();

            RooRealVar *prompt_yield = dynamic_cast<RooRealVar *>(coefList.at(0));
            RooRealVar *nonprompt_yield = dynamic_cast<RooRealVar *>(coefList.at(1));

            if (!prompt_yield || !nonprompt_yield)
            {
                std::cerr << "Error: Failed to cast coefficients to RooRealVar" << std::endl;
                delete params;
            }

            std::cout << "  Using prompt yield: " << prompt_yield->getVal() << " ± " << prompt_yield->getError() << std::endl;
            std::cout << "  Using nonprompt yield: " << nonprompt_yield->getVal() << " ± " << nonprompt_yield->getError() << std::endl;

            // Create sPlot object for IP chi2 analysis
            RooArgList yieldsList;
            yieldsList.add(*prompt_yield);
            yieldsList.add(*nonprompt_yield);

            std::cout << "  Creating SPlot object..." << std::endl;
            RooStats::SPlot *splotObj = new RooStats::SPlot("ipchi2_splotObj", "ipchi2_splotObj",
                                                            *data, splot_model, yieldsList);

            std::cout << "  SPlot object created successfully" << std::endl;

            // Save to ROOT file
            splotFile->cd();
            TTree *ipSplotTree = new TTree(ipSplotTreeName.c_str(), "IP chi2 sPlot weights");

            // Declare branch variables
            double mass_val, log_ipchi2_val, prompt_sWeight, nonprompt_sWeight;
            double tagZ_val, tagY_val; // Add tagZ and tagY for event matching

            // Create branches
            ipSplotTree->Branch("mass", &mass_val, "mass/D");
            ipSplotTree->Branch("log_ipchi2", &log_ipchi2_val, "log_ipchi2/D");
            ipSplotTree->Branch("prompt_sWeight", &prompt_sWeight, "prompt_sWeight/D");
            ipSplotTree->Branch("nonprompt_sWeight", &nonprompt_sWeight, "nonprompt_sWeight/D");
            ipSplotTree->Branch("tagZ", &tagZ_val, "tagZ/D");
            ipSplotTree->Branch("tagY", &tagY_val, "tagY/D");

            std::cout << "  Filling tree with sPlot weights..." << std::endl;

            // Fill the tree with sPlot weights
            int validEntries = 0;
            for (int i = 0; i < data->numEntries(); i++)
            {
                const RooArgSet *row = data->get(i);
                if (!row)
                    continue;

                // Get event variables
                RooRealVar *massVar = (RooRealVar *)row->find("tagMass");
                RooRealVar *ipChi2Var = (RooRealVar *)row->find("log_tag_ipchi2");
                RooRealVar *tagZVar = (RooRealVar *)row->find("tagZ");
                RooRealVar *tagYVar = (RooRealVar *)row->find("tagY");

                if (massVar && ipChi2Var)
                {
                    mass_val = massVar->getVal();
                    log_ipchi2_val = ipChi2Var->getVal();
                    tagZ_val = tagZVar ? tagZVar->getVal() : -999.0;
                    tagY_val = tagYVar ? tagYVar->getVal() : -999.0;

                    // Get sPlot weights with error checking
                    try
                    {
                        prompt_sWeight = splotObj->GetSWeight(i, prompt_yield->GetName());
                        nonprompt_sWeight = splotObj->GetSWeight(i, nonprompt_yield->GetName());

                        // Check for valid weights
                        if (std::isfinite(prompt_sWeight) && std::isfinite(nonprompt_sWeight))
                        {
                            ipSplotTree->Fill();
                            validEntries++;
                        }
                    }
                    catch (const std::exception &e)
                    {
                        std::cerr << "Warning: Failed to get sPlot weight for event " << i << ": " << e.what() << std::endl;
                    }
                }
            }

            // Write tree to file
            ipSplotTree->Write();
        }

        // Clean up
        delete result;
        // Note: Do not close splotFile here as it's managed by the caller
        return std::make_tuple(histogram, parameterArr, parameterErrArr, total_pdf, sig_yieldLim, prompt_frac);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error in ipchi2FitWithYields: " << e.what() << std::endl;
        return std::make_tuple(nullptr, parameterArr, parameterErrArr, nullptr, nullptr, nullptr);
    }
}

/**
 * Creates a prompt signal tagZ distribution by applying both mass sPlot weights and IP chi2 sPlot weights.
 *
 * This method implements a two-stage statistical separation:
 * 1. Mass sPlot weights: Separate signal from combinatorial background
 * 2. IP chi2 sPlot weights: Separate prompt from nonprompt (displaced) signal
 *
 * The resulting distribution represents D0 mesons that are:
 * - True signal candidates (not combinatorial background)
 * - AND prompt signal (not displaced/nonprompt)
 *
 * Statistical interpretation:
 * - Combined weight = mass_signal_weight × prompt_weight
 * - This accounts for both background subtraction and prompt/nonprompt separation
 * - The result is the tagZ distribution for prompt signal only
 *
 * @param data: RooDataSet containing the original data
 * @param splotFileName: File containing mass sPlot weights
 * @param bin: Bin number for identification
 * @param ipChi2Weights: Map containing IP chi2 sPlot weights (prompt, nonprompt)
 * @param histName: Base name for the output histogram
 * @param nBins: Number of bins for the tagZ histogram
 * @param xMin: Minimum tagZ value
 * @param xMax: Maximum tagZ value
 * @return: Histogram of prompt signal tagZ distribution
 */
TH1D *Fitter::createPromptSignalTagZDistribution(RooDataSet *data,
                                                 const std::string &splotFileName,
                                                 int bin,
                                                 const std::map<std::pair<double, double>, std::pair<double, double>> &ipChi2Weights,
                                                 const std::string &histName,
                                                 int nBins,
                                                 double xMin,
                                                 double xMax)
{

    std::cout << "Creating prompt signal tagZ distribution for bin " << bin << std::endl;

    if (!data)
    {
        std::cerr << "Error: Null dataset provided" << std::endl;
        return nullptr;
    }

    // Create histogram
    std::string fullHistName = histName + "_bin" + std::to_string(bin);
    TH1D *promptTagZHist = new TH1D(fullHistName.c_str(),
                                    ("Prompt Signal TagZ Distribution - Bin " + std::to_string(bin)).c_str(),
                                    nBins, xMin, xMax);
    promptTagZHist->GetXaxis()->SetTitle("#it{z}_{T}");
    promptTagZHist->GetYaxis()->SetTitle("Weighted Entries");

    // Open sPlot file to read mass weights
    TFile *splotFile = TFile::Open(splotFileName.c_str(), "READ");
    if (!splotFile || splotFile->IsZombie())
    {
        std::cerr << "Error: Cannot open sPlot file: " << splotFileName << std::endl;
        delete promptTagZHist;
        return nullptr;
    }

    // Read mass sPlot tree
    TTree *massSplotTree = (TTree *)splotFile->Get(("splotTree_bin" + std::to_string(bin)).c_str());
    if (!massSplotTree)
    {
        std::cerr << "Error: Cannot find mass sPlot tree for bin " << bin << std::endl;
        splotFile->Close();
        delete splotFile;
        delete promptTagZHist;
        return nullptr;
    }

    // Variables to read from mass sPlot tree
    double mass_splot, sig_sWeight_mass, bkg_sWeight_mass;
    massSplotTree->SetBranchAddress("mass", &mass_splot);
    massSplotTree->SetBranchAddress("sig_sWeight", &sig_sWeight_mass);
    massSplotTree->SetBranchAddress("bkg_sWeight", &bkg_sWeight_mass);

    // Create map to store mass sPlot weights
    std::map<double, double> massSigWeights;

    // Read mass sPlot weights
    std::cout << "  Reading mass sPlot weights..." << std::endl;
    for (Long64_t i = 0; i < massSplotTree->GetEntries(); ++i)
    {
        massSplotTree->GetEntry(i);
        massSigWeights[mass_splot] = sig_sWeight_mass;
    }
    std::cout << "  Read " << massSigWeights.size() << " mass sPlot weight entries" << std::endl;

    splotFile->Close();
    delete splotFile;

    // Apply combined weights to create prompt signal tagZ distribution
    int processedEvents = 0;
    int weightedEvents = 0;
    int massSplotFound = 0;
    int ipChi2SplotFound = 0;
    double totalWeight = 0.0;
    const double tolerance = 1e-10; // Tolerance for floating-point comparison

    for (int i = 0; i < data->numEntries(); ++i)
    {
        const RooArgSet *row = data->get(i);
        RooRealVar *tagZVar = dynamic_cast<RooRealVar *>(row->find("tagZ"));
        RooRealVar *massVar = dynamic_cast<RooRealVar *>(row->find("tagMass"));
        RooRealVar *ipChi2Var = dynamic_cast<RooRealVar *>(row->find("log_tag_ipchi2"));

        if (tagZVar && massVar && ipChi2Var)
        {
            processedEvents++;
            double tagZ = tagZVar->getVal();
            double mass = massVar->getVal();
            double ipchi2 = ipChi2Var->getVal();

            // Find mass sPlot weight (signal weight for background subtraction)
            double massSigWeight = 0.0;
            auto massSigIt = massSigWeights.find(mass);
            if (massSigIt != massSigWeights.end())
            {
                massSigWeight = massSigIt->second;
                massSplotFound++;
            }
            else
            {
                // Try with tolerance for floating-point precision
                for (const auto &weightPair : massSigWeights)
                {
                    if (std::abs(weightPair.first - mass) < tolerance)
                    {
                        massSigWeight = weightPair.second;
                        massSplotFound++;
                        break;
                    }
                }
            }

            // Find IP chi2 sPlot weight (prompt weight for prompt/nonprompt separation)
            double promptWeight = 0.0;
            auto promptWeightIt = ipChi2Weights.find({mass, ipchi2});
            if (promptWeightIt != ipChi2Weights.end())
            {
                promptWeight = promptWeightIt->second.first;
                ipChi2SplotFound++;
            }
            else
            {
                // Try with tolerance for floating-point precision
                for (const auto &weightPair : ipChi2Weights)
                {
                    if (std::abs(weightPair.first.first - mass) < tolerance &&
                        std::abs(weightPair.first.second - ipchi2) < tolerance)
                    {
                        promptWeight = weightPair.second.first;
                        ipChi2SplotFound++;
                        break;
                    }
                }
            }

            // Combined weight: mass signal weight * prompt weight
            double combinedWeight = massSigWeight * promptWeight;

            if (combinedWeight > 0)
            {
                promptTagZHist->Fill(tagZ, combinedWeight);
                weightedEvents++;
                totalWeight += combinedWeight;
            }
        }
    }

    std::cout << "  Processed " << processedEvents << " events" << std::endl;
    std::cout << "  " << massSplotFound << " events matched mass sPlot weights" << std::endl;
    std::cout << "  " << ipChi2SplotFound << " events matched IP chi2 sPlot weights" << std::endl;
    std::cout << "  " << weightedEvents << " events had positive combined weights" << std::endl;
    std::cout << "  Total combined weight: " << totalWeight << std::endl;
    std::cout << "  Histogram entries: " << promptTagZHist->GetEntries() << std::endl;
    std::cout << "  Histogram integral: " << promptTagZHist->Integral() << std::endl;

    return promptTagZHist;
}

RooDataSet *Fitter::createWeightedDataset(RooDataSet *originalData,
                                          const std::string &splotFileName,
                                          int bin,
                                          const std::string &weightType,
                                          const std::string &datasetName)
{

    std::cout << "Creating weighted dataset from sPlot file: " << splotFileName << std::endl;

    try
    {
        // Small delay to ensure file is completely closed from previous operations
        std::this_thread::sleep_for(std::chrono::milliseconds(100));

        // Try to open the file multiple times if it's not immediately available
        TFile *splotFile = nullptr;
        int maxRetries = 5;
        int retryCount = 0;

        while (retryCount < maxRetries)
        {
            splotFile = new TFile(splotFileName.c_str(), "READ");
            if (splotFile && !splotFile->IsZombie())
            {
                break; // File opened successfully
            }

            if (splotFile)
            {
                delete splotFile;
                splotFile = nullptr;
            }

            retryCount++;
            std::cout << "  Retry " << retryCount << " to open sPlot file..." << std::endl;
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
        }

        if (!splotFile || splotFile->IsZombie())
        {
            std::cerr << "Error: Cannot open sPlot file after " << maxRetries << " retries: " << splotFileName << std::endl;
            if (splotFile)
            {
                delete splotFile;
            }
            return nullptr;
        }

        // Get the sPlot tree for this bin
        std::string treeName = "splotTree_bin" + std::to_string(bin);
        TTree *splotTree = (TTree *)splotFile->Get(treeName.c_str());
        if (!splotTree)
        {
            std::cerr << "Error: Cannot find sPlot tree: " << treeName << std::endl;
            splotFile->Close();
            delete splotFile;
            return nullptr;
        }

        // Set up branches
        double mass_val, weight_val;
        splotTree->SetBranchAddress("mass", &mass_val);
        splotTree->SetBranchAddress(weightType.c_str(), &weight_val);

        // Create a map of mass values to weights using a tolerance for floating point comparison
        std::map<double, double> massToWeight;
        Long64_t nEntries = splotTree->GetEntries();

        for (Long64_t i = 0; i < nEntries; i++)
        {
            splotTree->GetEntry(i);
            massToWeight[mass_val] = weight_val;
        }

        std::cout << "Loaded " << nEntries << " sPlot weights from tree" << std::endl;

        // Create weight variable for the dataset
        RooRealVar datasetWeight("datasetWeight", "Dataset weight", 0.0);

        // Create argument set that includes weight variable
        RooArgSet datasetVars(*originalData->get());
        datasetVars.add(datasetWeight);

        // Create new weighted dataset with weight support
        RooDataSet *weightedData = new RooDataSet(datasetName.c_str(),
                                                  ("Weighted dataset: " + datasetName).c_str(),
                                                  datasetVars,
                                                  RooFit::WeightVar(datasetWeight));

        // Apply weights to the dataset
        int matchedEntries = 0;
        double totalWeight = 0.0;

        for (int i = 0; i < originalData->numEntries(); i++)
        {
            const RooArgSet *row = originalData->get(i);
            RooRealVar *massVar = (RooRealVar *)row->find("tagMass");

            if (massVar)
            {
                double mass = massVar->getVal();

                // Find the closest mass value in the weight map (with tolerance)
                double closestMass = -1;
                double minDiff = 1e6;
                for (const auto &pair : massToWeight)
                {
                    double diff = std::abs(pair.first - mass);
                    if (diff < minDiff)
                    {
                        minDiff = diff;
                        closestMass = pair.first;
                    }
                }

                if (closestMass > 0 && minDiff < 1e-6)
                { // Tolerance for floating point comparison
                    double weight = massToWeight[closestMass];

                    // Only add entries with positive weights
                    if (weight > 0)
                    {
                        // Create a copy of the row and add the weight
                        RooArgSet weightedRow(*row);
                        datasetWeight.setVal(weight);
                        weightedRow.add(datasetWeight);
                        weightedData->add(weightedRow, weight);
                        totalWeight += weight;
                        matchedEntries++;
                    }
                }
            }
        }

        std::cout << "Created weighted dataset with " << matchedEntries << " entries" << std::endl;
        std::cout << "Total weight: " << totalWeight << std::endl;
        std::cout << "Effective entries: " << weightedData->sumEntries() << std::endl;

        // Close the file properly
        if (splotFile)
        {
            splotFile->Close();
            delete splotFile;
            splotFile = nullptr;
        }

        return weightedData;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error in createWeightedDataset: " << e.what() << std::endl;
        return nullptr;
    }
}

// Save IP chi2 sPlot weights to ROOT file
bool Fitter::saveIPChi2SPlotWeights(RooDataSet *originalData,
                                    RooAbsPdf *model,
                                    RooRealVar *sig_yieldLim,
                                    RooRealVar *prompt_frac,
                                    TFile *outputFile,
                                    int bin,
                                    const std::string &treeName)
{

    std::cout << "Saving IP chi2 sPlot weights to ROOT file for bin " << bin << std::endl;

    try
    {
        // Validate inputs
        if (!originalData || !model || !outputFile)
        {
            std::cerr << "Error: Invalid input parameters to saveIPChi2SPlotWeights" << std::endl;
            return false;
        }

        // Check if the model has parameters
        RooArgSet *params = model->getParameters(*originalData);
        if (!params || params->getSize() == 0)
        {
            std::cerr << "Error: Model has no parameters or failed to get parameters" << std::endl;
            return false;
        }

        std::cout << "  Model has " << params->getSize() << " parameters" << std::endl;

        // Get the yield variables from the model coefficients
        // For a RooAddPdf, the coefficients should be the yield variables
        RooAddPdf *addPdf = dynamic_cast<RooAddPdf *>(model);
        if (!addPdf)
        {
            std::cerr << "Error: Model is not a RooAddPdf" << std::endl;
            delete params;
            return false;
        }

        RooArgList coefList = addPdf->coefList();
        if (coefList.getSize() != 2)
        {
            std::cerr << "Error: Expected 2 coefficients (yields) in the model, got " << coefList.getSize() << std::endl;
            delete params;
            return false;
        }

        RooRealVar *prompt_yield = dynamic_cast<RooRealVar *>(coefList.at(0));
        RooRealVar *nonprompt_yield = dynamic_cast<RooRealVar *>(coefList.at(1));

        if (!prompt_yield || !nonprompt_yield)
        {
            std::cerr << "Error: Failed to cast coefficients to RooRealVar" << std::endl;
            delete params;
            return false;
        }

        std::cout << "  Using prompt yield: " << prompt_yield->getVal() << " ± " << prompt_yield->getError() << std::endl;
        std::cout << "  Using nonprompt yield: " << nonprompt_yield->getVal() << " ± " << nonprompt_yield->getError() << std::endl;

        // Create sPlot object for IP chi2 analysis
        RooArgList yieldsList;
        yieldsList.add(*prompt_yield);
        yieldsList.add(*nonprompt_yield);

        std::cout << "  Creating SPlot object..." << std::endl;
        RooStats::SPlot *splotObj = new RooStats::SPlot("ipchi2_splotObj", "ipchi2_splotObj",
                                                        *originalData, model, yieldsList);

        std::cout << "  SPlot object created successfully" << std::endl;

        // Save to ROOT file
        outputFile->cd();
        TTree *ipSplotTree = new TTree(treeName.c_str(), "IP chi2 sPlot weights");

        // Declare branch variables
        double mass_val, log_ipchi2_val, prompt_sWeight, nonprompt_sWeight;
        double tagZ_val, tagY_val; // Add tagZ and tagY for event matching

        // Create branches
        ipSplotTree->Branch("mass", &mass_val, "mass/D");
        ipSplotTree->Branch("log_ipchi2", &log_ipchi2_val, "log_ipchi2/D");
        ipSplotTree->Branch("prompt_sWeight", &prompt_sWeight, "prompt_sWeight/D");
        ipSplotTree->Branch("nonprompt_sWeight", &nonprompt_sWeight, "nonprompt_sWeight/D");
        ipSplotTree->Branch("tagZ", &tagZ_val, "tagZ/D");
        ipSplotTree->Branch("tagY", &tagY_val, "tagY/D");

        std::cout << "  Filling tree with sPlot weights..." << std::endl;

        // Fill the tree with sPlot weights
        int validEntries = 0;
        for (int i = 0; i < originalData->numEntries(); i++)
        {
            const RooArgSet *row = originalData->get(i);
            if (!row)
                continue;

            // Get event variables
            RooRealVar *massVar = (RooRealVar *)row->find("tagMass");
            RooRealVar *ipChi2Var = (RooRealVar *)row->find("log_tag_ipchi2");
            RooRealVar *tagZVar = (RooRealVar *)row->find("tagZ");
            RooRealVar *tagYVar = (RooRealVar *)row->find("tagY");

            if (massVar && ipChi2Var)
            {
                mass_val = massVar->getVal();
                log_ipchi2_val = ipChi2Var->getVal();
                tagZ_val = tagZVar ? tagZVar->getVal() : -999.0;
                tagY_val = tagYVar ? tagYVar->getVal() : -999.0;

                // Get sPlot weights with error checking
                try
                {
                    prompt_sWeight = splotObj->GetSWeight(i, prompt_yield->GetName());
                    nonprompt_sWeight = splotObj->GetSWeight(i, nonprompt_yield->GetName());

                    // Check for valid weights
                    if (std::isfinite(prompt_sWeight) && std::isfinite(nonprompt_sWeight))
                    {
                        ipSplotTree->Fill();
                        validEntries++;
                    }
                }
                catch (const std::exception &e)
                {
                    std::cerr << "Warning: Failed to get sPlot weight for event " << i << ": " << e.what() << std::endl;
                }
            }
        }

        // Write tree to file
        ipSplotTree->Write();

        std::cout << "  IP chi2 sPlot tree '" << treeName << "' saved with "
                  << validEntries << "/" << originalData->numEntries() << " valid entries" << std::endl;

        // Calculate and print summary statistics
        double totalPromptWeight = 0.0;
        double totalNonpromptWeight = 0.0;
        int positivePromptWeights = 0;
        int positiveNonpromptWeights = 0;

        std::cout << "  Calculating summary statistics..." << std::endl;
        for (int i = 0; i < originalData->numEntries(); i++)
        {
            try
            {
                double pWeight = splotObj->GetSWeight(i, prompt_yield->GetName());
                double npWeight = splotObj->GetSWeight(i, nonprompt_yield->GetName());

                if (std::isfinite(pWeight) && pWeight > 0)
                {
                    totalPromptWeight += pWeight;
                    positivePromptWeights++;
                }
                if (std::isfinite(npWeight) && npWeight > 0)
                {
                    totalNonpromptWeight += npWeight;
                    positiveNonpromptWeights++;
                }
            }
            catch (const std::exception &e)
            {
                // Skip problematic events in statistics
                continue;
            }
        }

        std::cout << "  Summary:" << std::endl;
        std::cout << "    Total prompt weight: " << totalPromptWeight
                  << " (events with positive weights: " << positivePromptWeights << ")" << std::endl;
        std::cout << "    Total nonprompt weight: " << totalNonpromptWeight
                  << " (events with positive weights: " << positiveNonpromptWeights << ")" << std::endl;

        // Clean up
        delete splotObj;
        delete params; // Clean up the parameters list

        return true;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error in saveIPChi2SPlotWeights: " << e.what() << std::endl;
        return false;
    }
}

/*!SECTION

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

        // Debug: Print the prompt_frac value that was used
        std::cout << "  Created prompt_frac with value: " << prompt_frac->getVal() << std::endl;

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
        RooAbsPdf* total_pdf = new RooAddPdf("total_pdf", "total_pdf",
                                           RooArgList(*prompt_pdf, *nonprompt_pdf),
                                           RooArgList(*prompt_yield, *nonprompt_yield));

        // Create background PDF if provided
        // if (background && background->numEntries() > 0) {
        //     RooRealVar* bkg_param1 = new RooRealVar("bkg_param1", "bkg_param1",
        //                                          ipchi2_params.bkgParam1.value,
        //                                          ipchi2_params.bkgParam1.min,
        //                                          ipchi2_params.bkgParam1.max);

        //     RooRealVar* bkg_param2 = new RooRealVar("bkg_param2", "bkg_param2",
        //                                          ipchi2_params.bkgParam2.value,
        //                                          ipchi2_params.bkgParam2.min,
        //                                          ipchi2_params.bkgParam2.max);

        //     // Use a polynomial for background
        //     background_pdf = new RooPolynomial("bkg_pdf", "bkg_pdf",
        //                                    *log_ipchi2, RooArgList(*bkg_param1, *bkg_param2));

        //     // Combined model with signal and background
        //     RooAddPdf* signal_model = new RooAddPdf("signal_model", "signal_model",
        //                                          RooArgList(*prompt_pdf, *nonprompt_pdf),
        //                                          RooArgList(*prompt_yield, *nonprompt_yield));

        //     total_pdf = new RooAddPdf("total_pdf", "total_pdf",
        //                              RooArgList(*signal_model, *background_pdf),
        //                              RooArgList(*sig_yieldLim, *bkg_yieldLim));
        // } else {
            // Signal-only model
            // total_pdf is already created above
        // }

        // Perform the fit
        std::cout << "  Performing IP chi2 fit with Bukin function..." << std::endl;
        std::cout << "  Initial prompt_frac: " << prompt_frac->getVal() << " (range: "
                  << prompt_frac->getMin() << " - " << prompt_frac->getMax() << ")" << std::endl;
        std::cout << "  Initial sig_yieldLim: " << sig_yieldLim->getVal() << std::endl;

        // Debug: Print the sig_yieldLim value that was used
        std::cout << "  Created sig_yieldLim with value: " << sig_yieldLim->getVal() << std::endl;
        std::cout << "  Dictionary mass_params.sigYieldLim.value: " << mass_params.sigYieldLim.value << std::endl;

        // Check if prompt_frac is constant (fixed)
        if (prompt_frac->isConstant()) {
            std::cout << "  WARNING: prompt_frac is constant! Setting it to variable..." << std::endl;
            prompt_frac->setConstant(false);
        }

        // Print all parameter states before fit
        std::cout << "  Parameter states before fit:" << std::endl;
        std::cout << "    prompt_frac: " << prompt_frac->getVal() << " [" << prompt_frac->getMin()
                  << ", " << prompt_frac->getMax() << "] constant=" << prompt_frac->isConstant() << std::endl;
        std::cout << "    sig_yieldLim: " << sig_yieldLim->getVal() << " [" << sig_yieldLim->getMin()
                  << ", " << sig_yieldLim->getMax() << "] constant=" << sig_yieldLim->isConstant() << std::endl;

        // Make sure sig_yieldLim is also not constant
        if (sig_yieldLim->isConstant()) {
            std::cout << "  WARNING: sig_yieldLim is constant! Setting it to variable..." << std::endl;
            sig_yieldLim->setConstant(false);
        }

        // Use improved fit strategy
        RooFitResult* result = total_pdf->fitTo(*data,
                                              RooFit::Save(true),
                                              RooFit::PrintLevel(0),
                                              RooFit::SumW2Error(true),
                                              RooFit::Strategy(2),
                                              RooFit::Minos(false),
                                              RooFit::Hesse(true));

        std::cout << "  Fitted prompt_frac: " << prompt_frac->getVal() << " ± " << prompt_frac->getError() << std::endl;
        std::cout << "  Fitted sig_yieldLim: " << sig_yieldLim->getVal() << std::endl;
        std::cout << "  Fitted prompt_yield: " << prompt_yield->getVal() << std::endl;
        std::cout << "  Fitted nonprompt_yield: " << nonprompt_yield->getVal() << std::endl;

        // Create plot
        std::cout << "  Creating IP chi2 fit plot..." << std::endl;
        Plotter plotter(resonance, outfilePath, bin, false, zRange, figKey);
        histogram = plotter.ipchi2FitPlot(resonance, log_ipchi2, data, total_pdf,
                                        nonprompt_pdf, prompt_pdf, nullptr,
                                        prompt_yield, nonprompt_yield);

        // Extract fit parameters
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

        // Extract errors
        parameterErrArr[0] = sig_yield->getError();
        parameterErrArr[1] = prompt_frac->getError();

        parameterErrArr[2] = xp_prompt->getError();
        parameterErrArr[3] = sigma_prompt->getError();
        parameterErrArr[4] = xi_prompt->getError();

        parameterErrArr[7] = xp_nonprompt->getError();
        parameterErrArr[8] = sigma_nonprompt->getError();
        parameterErrArr[9] = xi_nonprompt->getError();

        std::cout << "  IP chi2 fit with Bukin function completed successfully" << std::endl;
        std::cout << "  Prompt fraction: " << prompt_frac->getVal() << " ± " << prompt_frac->getError() << std::endl;
        std::cout << "  Prompt yield: " << prompt_yield->getVal() << std::endl;
        std::cout << "  Non-prompt yield: " << nonprompt_yield->getVal() << std::endl;

        // Clean up fit result
        delete result;

        // Create RooRealVar objects for the yields (for sPlot compatibility)
        RooRealVar* promptYieldVar = new RooRealVar("prompt_yield_var", "prompt_yield_var",
                                                   prompt_yield->getVal(), 0, 1e6);
        promptYieldVar->setError(prompt_frac->getError());

        RooRealVar* nonpromptYieldVar = new RooRealVar("nonprompt_yield_var", "nonprompt_yield_var",
                                                      nonprompt_yield->getVal(), 0, 1e6);
        nonpromptYieldVar->setError(prompt_frac->getError());

        // Return results (note: caller is responsible for cleanup of yield variables and PDF)
        return std::make_tuple(histogram, parameterArr, parameterErrArr);

    } catch (const std::exception& e) {
        std::cerr << "Error in ipchi2Fit: " << e.what() << std::endl;
        return std::make_tuple(nullptr, parameterArr, parameterErrArr);
    }
}

*/