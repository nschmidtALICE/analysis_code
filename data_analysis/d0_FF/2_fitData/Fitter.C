#include "PlotHelpers.h"
// Fitter.cpp
#include "Fitter.h"
#include "Plotter.h"
#include <chrono>
#include <thread>
#include "TSystem.h"
#include "RooStats/SPlot.h"

Fitter::Fitter(TTree *tree, const std::string &resonanceType, int numBins, bool zTObservable,
                             bool isMCData, const std::string &outputPath, bool update, const std::string &inputFile)
        : tfilePV(nullptr), TestFilename(""),
            outfilePath(outputPath), inputFileName(inputFile), isMC(isMCData), nBins(numBins), isZtObservable(zTObservable),
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
    massDict.clear();
    ipchi2Dict.clear();

    // Initialize D0 mass parameters
    MassConfig d0Config;
    d0Config.sigma1 = ParamConfig(0.008, 0.004, 0.009);
    d0Config.deltasigma = ParamConfig(1.5, 1.3, 2.2);
    d0Config.mean = ParamConfig(1.865, 1.860, 1.870);
    d0Config.alpha1 = ParamConfig(4.5, 0.5, 10.0);
    d0Config.n = ParamConfig(0.5, 0.15, 5);
    d0Config.dg_frac = ParamConfig(0.5, 0.0, 0.99999);
    d0Config.pol1 = ParamConfig(90, 40, 200);
    // d0Config.pol1 = ParamConfig(90, -1, 100);
    d0Config.pol2 = ParamConfig(58, -2e2, 5e2);
    d0Config.massRange = std::make_pair(1.815, 1.925); // 150 MeV range
    d0Config.sigYield = ParamConfig(100, 1, 15e6);
    d0Config.sigYieldLim = ParamConfig(100, 0, 15e6);
    d0Config.bkgYield = ParamConfig(1000, 0, 2e6);
    d0Config.bkgYieldLim = ParamConfig(1000, 0, 2e6);
    // d0Config.signalRegion = std::make_pair(1.845, 1.885); // 40 MeV range
    d0Config.signalRegion = std::make_pair(1.830, 1.900); // 35 MeV range = 5 sigma integral
    d0Config.sbRegion = std::make_pair(1.825, 1.910);     // Sideband region

    massDict["D0"] = d0Config;

    // Initialize D0 IP chi2 parameters
    IPChi2Config d0IPConfig;
    d0IPConfig.logIpchi2Range = std::make_pair(-3, 5);

    // Prompt component
    d0IPConfig.xpPrompt = ParamConfig(0.3, 0.23, 0.38);
    d0IPConfig.sigmaPrompt = ParamConfig(0.47, 0.43, 0.53);
    d0IPConfig.xiPrompt = ParamConfig(-0.18, -0.21, -0.15);
    d0IPConfig.rho1Prompt = ParamConfig(-0.08, -0.11, -0.05);
    d0IPConfig.rho2Prompt = ParamConfig(0.01, 0.0001, 0.98);

    // Non-prompt component
    d0IPConfig.xpNonprompt = ParamConfig(1.95, 1.87, 2.27);
    d0IPConfig.sigmaNonprompt = ParamConfig(0.45, 0.3, 0.66);
    d0IPConfig.xiNonprompt = ParamConfig(0.15, 0.05, 0.28);
    d0IPConfig.rho1Nonprompt = ParamConfig(-0.95, -0.96, -0.05);
    d0IPConfig.rho2Nonprompt = ParamConfig(0.2, 0.01, 0.98);
    // d0IPConfig.xpPrompt = ParamConfig(0.0, -0.5, 1.0);
    // d0IPConfig.sigmaPrompt = ParamConfig(0.7, 0.4, 2.0);
    // d0IPConfig.xiPrompt = ParamConfig(0.0, -0.5, 0.5);
    // d0IPConfig.rho1Prompt = ParamConfig(-0.08, -0.2, -0.01);
    // d0IPConfig.rho2Prompt = ParamConfig(0.01, 0.0001, 0.98);

    // // Non-prompt component
    // d0IPConfig.xpNonprompt = ParamConfig(1.9, 1.6, 3.0);
    // d0IPConfig.sigmaNonprompt = ParamConfig(0.4, 0.3, 0.7);
    // d0IPConfig.xiNonprompt = ParamConfig(0.1, 0.05, 0.5);
    // d0IPConfig.rho1Nonprompt = ParamConfig(-0.1, -0.95, -0.05);
    // d0IPConfig.rho2Nonprompt = ParamConfig(0.2, 0.01, 0.98);

    // Fraction - starting with a more central value to avoid boundary issues
    d0IPConfig.promptFrac = ParamConfig(0.95, 0.7, 0.999);

    // Background
    d0IPConfig.bkgParam1 = ParamConfig(0.5, 0, 1);
    d0IPConfig.bkgParam2 = ParamConfig(0.5, 0, 1);

    ipchi2Dict["SigD0"] = d0IPConfig;
}

void Fitter::resetFitDictionaries()
{
    initDictionary();
}

void Fitter::applyMassPrefitConstraints(const std::string &resonance,
                                        const std::string &fitTypeName,
                                        const std::vector<double> &fitValues,
                                        const std::vector<double> &fitErrors,
                                        double yieldScale)
{
    auto it = massDict.find(resonance);
    if (it == massDict.end())
    {
        std::cerr << "Warning: cannot apply mass pre-fit constraints for unknown resonance '"
                  << resonance << "'" << std::endl;
        return;
    }

    if (fitValues.size() < 10 || fitErrors.size() < 10)
    {
        std::cerr << "Warning: insufficient pre-fit parameters to constrain the full fit" << std::endl;
        return;
    }

    if (!(yieldScale > 0.0) || !std::isfinite(yieldScale))
    {
        yieldScale = 1.0;
    }

    auto &res = it->second;
    const bool useCrystalBall = (fitTypeName == "CBall" || fitTypeName == "DCB");

    auto clampValue = [](double value, double minValue, double maxValue) {
        return std::max(minValue, std::min(value, maxValue));
    };

    auto updateParameter = [&](ParamConfig &param, double center, double error,
                               double relativeSpan, double absoluteSpan) {
        if (!std::isfinite(center))
        {
            return;
        }

        double span = std::max(absoluteSpan, std::abs(center) * relativeSpan);
        if (std::isfinite(error) && error > 0.0)
        {
            span = std::max(span, 5.0 * error);
        }

        double newMin = std::max(param.min, center - span);
        double newMax = std::min(param.max, center + span);
        if (newMin >= newMax)
        {
            newMin = param.min;
            newMax = param.max;
        }

        param.value = clampValue(center, newMin, newMax);
        param.min = newMin;
        param.max = newMax;
    };

    auto updateYield = [&](ParamConfig &param, double value, double error) {
        if (!std::isfinite(value))
        {
            return;
        }

        double span = std::max(25.0, std::abs(value) * 0.60);
        if (std::isfinite(error) && error > 0.0)
        {
            span = std::max(span, 5.0 * error);
        }

        double newMin = std::max(param.min, value - span);
        double newMax = std::min(param.max, value + span);
        if (newMin >= newMax)
        {
            newMin = param.min;
            newMax = param.max;
        }

        param.value = clampValue(value, newMin, newMax);
        param.min = newMin;
        param.max = newMax;
    };

    updateParameter(res.mean, fitValues[2], fitErrors[2], 0.001, 0.0015);
    updateParameter(res.sigma1, fitValues[3], fitErrors[3], 0.25, 0.0005);
    updateParameter(res.deltasigma, fitValues[4], fitErrors[4], 0.30, 0.08);
    updateParameter(res.dg_frac, fitValues[7], fitErrors[7], 0.35, 0.05);
    updateParameter(res.pol1, fitValues[8], fitErrors[8], 0.30, 10.0);
    updateParameter(res.pol2, fitValues[9], fitErrors[9], 0.40, 20.0);

    if (useCrystalBall)
    {
        updateParameter(res.alpha1, fitValues[5], fitErrors[5], 0.35, 0.25);
        updateParameter(res.n, fitValues[6], fitErrors[6], 0.45, 0.20);
    }

    const double scaledSigYield = fitValues[0] * yieldScale;
    const double scaledBkgYield = fitValues[1] * yieldScale;
    const double scaledSigYieldErr = fitErrors[0] * yieldScale;
    const double scaledBkgYieldErr = fitErrors[1] * yieldScale;

    updateYield(res.sigYield, scaledSigYield, scaledSigYieldErr);
    updateYield(res.bkgYield, scaledBkgYield, scaledBkgYieldErr);

    std::cout << "Applied mass pre-fit constraints for " << resonance
              << " using yield scale " << yieldScale << std::endl;
    std::cout << "  Mean constrained to [" << res.mean.min << ", " << res.mean.max << "]"
              << " around " << res.mean.value << std::endl;
    std::cout << "  Signal yield initialized to " << res.sigYield.value
              << " with bounds [" << res.sigYield.min << ", " << res.sigYield.max << "]" << std::endl;
    std::cout << "  Background yield initialized to " << res.bkgYield.value
              << " with bounds [" << res.bkgYield.min << ", " << res.bkgYield.max << "]" << std::endl;
}

void Fitter::applyIPChi2PrefitConstraints(const std::string &resonance,
                                          const std::vector<double> &fitValues,
                                          const std::vector<double> &fitErrors)
{
    const std::string key = "Sig" + resonance;
    auto it = ipchi2Dict.find(key);
    if (it == ipchi2Dict.end())
    {
        std::cerr << "Warning: cannot apply IPChi2 pre-fit constraints for unknown resonance '"
                  << key << "'" << std::endl;
        return;
    }

    if (fitValues.size() < 12 || fitErrors.size() < 12)
    {
        std::cerr << "Warning: insufficient IPChi2 pre-fit parameters to constrain the full fit" << std::endl;
        return;
    }

    auto &cfg = it->second;

    auto clampValue = [](double value, double minValue, double maxValue) {
        return std::max(minValue, std::min(value, maxValue));
    };

    auto updateParameter = [&](ParamConfig &param, double center, double error,
                               double relativeSpan, double absoluteSpan) {
        if (!std::isfinite(center))
        {
            return;
        }

        double span = std::max(absoluteSpan, std::abs(center) * relativeSpan);
        if (std::isfinite(error) && error > 0.0)
        {
            span = std::max(span, 5.0 * error);
        }

        double newMin = std::max(param.min, center - span);
        double newMax = std::min(param.max, center + span);
        if (newMin >= newMax)
        {
            newMin = param.min;
            newMax = param.max;
        }

        param.value = clampValue(center, newMin, newMax);
        param.min = newMin;
        param.max = newMax;
    };

    updateParameter(cfg.promptFrac, fitValues[1], fitErrors[1], 0.15, 0.02);
    updateParameter(cfg.xpPrompt, fitValues[2], fitErrors[2], 0.20, 0.04);
    updateParameter(cfg.sigmaPrompt, fitValues[3], fitErrors[3], 0.20, 0.03);
    updateParameter(cfg.xiPrompt, fitValues[4], fitErrors[4], 0.35, 0.03);
    updateParameter(cfg.rho1Prompt, fitValues[5], fitErrors[5], 0.35, 0.03);
    updateParameter(cfg.rho2Prompt, fitValues[6], fitErrors[6], 0.50, 0.02);
    updateParameter(cfg.xpNonprompt, fitValues[7], fitErrors[7], 0.20, 0.05);
    updateParameter(cfg.sigmaNonprompt, fitValues[8], fitErrors[8], 0.20, 0.04);
    updateParameter(cfg.xiNonprompt, fitValues[9], fitErrors[9], 0.35, 0.04);
    updateParameter(cfg.rho1Nonprompt, fitValues[10], fitErrors[10], 0.25, 0.08);
    updateParameter(cfg.rho2Nonprompt, fitValues[11], fitErrors[11], 0.50, 0.05);

    std::cout << "Applied IPChi2 pre-fit constraints for " << key << std::endl;
    std::cout << "  prompt_frac bounds: [" << cfg.promptFrac.min << ", " << cfg.promptFrac.max
              << "] around " << cfg.promptFrac.value << std::endl;
    std::cout << "  xp_prompt bounds: [" << cfg.xpPrompt.min << ", " << cfg.xpPrompt.max
              << "] around " << cfg.xpPrompt.value << std::endl;
    std::cout << "  xp_nonprompt bounds: [" << cfg.xpNonprompt.min << ", " << cfg.xpNonprompt.max
              << "] around " << cfg.xpNonprompt.value << std::endl;
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
        keyList = {"pol0", "pol1"};
    }
    else if (fitFunc == "DGauss")
    {
        keyList = {"mean", "sigma1", "deltasigma", "dg_frac", "sig_yield",
                   "bkg_yield", "pol0", "pol1"};
    }
    else if (fitFunc == "CBall" || fitFunc == "DCB")
    {
        keyList = {"mean", "sigma1", "deltasigma", "alpha1", "n", "cb_frac",
                   "sig_yield", "bkg_yield", "pol0", "pol1"};
    }
    else if (fitFunc == "SGauss")
    {
        keyList = {"mean", "sigma1", "sig_yield", "bkg_yield", "pol0", "pol1"};
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
            else if (key == "alpha1")
                res.alpha1.value = newVal;
            else if (key == "n")
                res.n.value = newVal;
            else if (key == "dg_frac")
                res.dg_frac.value = newVal;
            else if (key == "cb_frac")
                res.dg_frac.value = newVal;
            else if (key == "sig_yield")
                res.sigYield.value = newVal;
            else if (key == "bkg_yield")
                res.bkgYield.value = newVal;
            else if (key == "pol0")
                res.pol1.value = newVal;
            else if (key == "pol1")
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
    
    // Declare variables - consolidated variable creation
    auto mass_range = res.massRange;
    
    RooRealVar *tagMass = new RooRealVar("tagMass", "tagMass", mass_range.first, mass_range.second);
    RooRealVar *tag_ipchi2 = new RooRealVar("tag_ip_chi2", "tag_ip_chi2", 0, 10000);
    RooRealVar *log_tag_ipchi2 = new RooRealVar("log_tag_ipchi2", "log_tag_ipchi2",
                                                ipchi2_params.logIpchi2Range.first,
                                                ipchi2_params.logIpchi2Range.second);

    // Standard kinematic variables
    std::vector<std::pair<std::string, std::pair<double, double>>> kinVars = {
        {"jetPt", {0, 200}}, {"tagPt", {0, 200}}, {"jetnConst", {0.0, 300.0}},
        {"tagJetdR", {0.0, 1.0}}, {"tagZ", {0.0, 1.01}}, {"tagY", {0.0, 5.01}}
    };

    RooArgSet *cutVars = new RooArgSet();
    cutVars->add(*tagMass);
    cutVars->add(*tag_ipchi2);
    cutVars->add(*log_tag_ipchi2);
    
    // Add kinematic variables
    for (const auto& var : kinVars) {
        RooRealVar *rooVar = new RooRealVar(var.first.c_str(), var.first.c_str(), 
                                           var.second.first, var.second.second);
        cutVars->add(*rooVar);
    }

    // Add distance variables
    RooRealVar *distance1 = new RooRealVar("Distance1", "Distance1", -10, 200);
    cutVars->add(*distance1);

    std::cout << "Using correction version: " << corrVer << std::endl;

    // Add efficiency weights using a loop for cleaner code
    std::vector<std::string> effVars = {
        "kaon_efficiency", "pion_efficiency", "combined_efficiency", 
        "combined_PID_efficiency", "reconstruction_efficiency", 
        "acceptance", "combined_eff_and_acceptance"
    };
    
    for (const auto& varName : effVars) {
        RooRealVar *effVar = new RooRealVar(varName.c_str(), varName.c_str(), 0, 1);
        cutVars->add(*effVar);
    }

    // Create final cut string
    std::string finalCutString = fidCutString;
    finalCutString += "&& tagPt > 2";
    finalCutString += "&& jetnConst > 1";

    if (corrVer > -2)
    {
        finalCutString += "&& Distance1 < 0.5";
    }

    // Create dataset with appropriate weight
    RooDataSet *data = nullptr;
    std::cout << "Debug: Creating dataset with corrVer = " << corrVer << std::endl;

    // Map correction versions to weight names for cleaner code
    std::map<int, std::string> correctionWeights = {
        {1, "kaon_efficiency"},
        {2, "pion_efficiency"}, 
        {3, "reconstruction_efficiency"},
        {4, "acceptance"},
        {5, "combined_PID_efficiency"},
        {6, "combined_eff_and_acceptance"}
    };

    const char *weightName = nullptr;
    auto it = correctionWeights.find(corrVer);
    if (it != correctionWeights.end()) {
        weightName = it->second.c_str();
        std::cout << "Debug: Using " << weightName << " as weight" << std::endl;
    }

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
        const bool useCrystalBall = (fitTypeName == "CBall" || fitTypeName == "DCB");
        const bool useDoubleGaussian = (fitTypeName == "DGauss");
        const bool useTwoComponentSignal = useCrystalBall || useDoubleGaussian;

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
        RooRealVar *alpha1 = new RooRealVar("alpha1", "alpha1",
                                            res.alpha1.value, res.alpha1.min, res.alpha1.max);
        RooFormulaVar *alpha2 = new RooFormulaVar("alpha2", "alpha2", "-1.0*alpha1",
                                                  RooArgList(*alpha1));
        RooRealVar *n = new RooRealVar("n", "n",
                                       res.n.value, res.n.min, res.n.max);
        const std::string fractionName = useCrystalBall ? "cb_frac" : "dg_frac";
        RooRealVar *component_frac = new RooRealVar(fractionName.c_str(), fractionName.c_str(),
                                                    res.dg_frac.value, res.dg_frac.min, res.dg_frac.max);
        RooRealVar *sig_yield = new RooRealVar("sig_yield", "sig_yield",
                                               res.sigYield.value, res.sigYield.min, res.sigYield.max);

        // Create signal PDF
        RooAbsPdf *sig_pdf = nullptr;
        if (useDoubleGaussian)
        {
            // Double Gaussian
            RooGaussian *Gauss1_pdf = new RooGaussian("Sig1_pdf", "Sig1_pdf",
                                                      *mass_tag_measured, *mean, *sigma1);
            RooGaussian *Gauss2_pdf = new RooGaussian("Sig2_pdf", "Sig2_pdf",
                                                      *mass_tag_measured, *mean, *sigma2);
            sig_pdf = new RooAddPdf("sig_pdf", "Signal",
                                    RooArgList(*Gauss2_pdf, *Gauss1_pdf), RooArgList(*component_frac));
            std::cout << "  Using Double Gaussian PDF" << std::endl;
        }
        else if (useCrystalBall)
        {
            RooCBShape *CB1_pdf = new RooCBShape("Sig1_pdf", "Sig1_pdf",
                                                 *mass_tag_measured, *mean, *sigma1, *alpha1, *n);
            RooCBShape *CB2_pdf = new RooCBShape("Sig2_pdf", "Sig2_pdf",
                                                 *mass_tag_measured, *mean, *sigma2, *alpha2, *n);
            sig_pdf = new RooAddPdf("sig_pdf", "Signal",
                                    RooArgList(*CB2_pdf, *CB1_pdf), RooArgList(*component_frac));
            std::cout << "  Using mirrored Crystal Ball PDF" << std::endl;
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

        // First: perform a binned fit to a histogram representation of the data.
        // This provides robust starting values for the shape parameters.
        std::cout << "  Performing binned fit (" << nBins << " bins) to determine shape parameters..." << std::endl;

        // Ensure variable binning is set and create a binned dataset
        // Use a larger number of bins for diagnostics if nBins is small
        int binnedN = nBins > 10 ? nBins : 50;
        mass_tag_measured->setBins(binnedN);
        RooDataHist *binnedData = new RooDataHist("binnedData", "binnedData", RooArgSet(*mass_tag_measured), *data);

        RooFitResult *binnedFitResult = nullptr;
        if (fitTypeName == "noSig")
        {
            binnedFitResult = extended_pdf->fitTo(*binnedData, RooFit::Save(),
                                                 RooFit::PrintLevel(0),
                                                 RooFit::Extended(true),
                                                 RooFit::NumCPU(8),
                                                 RooFit::Strategy(1),
                                                 RooFit::Range("SBleft,SBright"),
                                                 RooFit::DataError(RooAbsData::SumW2));
        }
        else
        {
            binnedFitResult = extended_pdf->fitTo(*binnedData, RooFit::Save(),
                                                 RooFit::PrintLevel(0),
                                                 RooFit::Extended(true),
                                                 RooFit::NumCPU(8),
                                                 RooFit::Strategy(1),
                                                 RooFit::Range("fullRange"),
                                                 RooFit::DataError(RooAbsData::SumW2));
        }

        // Create a binned TH1 for plotting (returned and/or used by plotter)
        TH1 *binnedHist = data->createHistogram("mass_hist_binned", *mass_tag_measured, RooFit::Binning(binnedN));

        // Quick diagnostic plot: binned data and fit overlay
        {
            std::string binstr = (bin >= 0) ? std::to_string(bin) : std::string("all");

            RooPlot *frame_binned = mass_tag_measured->frame(RooFit::Title("Binned mass fit"));
            binnedData->plotOn(frame_binned, RooFit::DataError(RooAbsData::SumW2));

            // Draw the fitted PDF on top
            extended_pdf->plotOn(frame_binned, RooFit::LineColor(kRed));
            if (binnedFitResult)
            {
                extended_pdf->plotOn(frame_binned, RooFit::VisualizeError(*binnedFitResult, 1), RooFit::DrawOption("F"), RooFit::FillColor(kOrange), RooFit::LineColor(kOrange));
            }

            TCanvas *c_binned = new TCanvas("c_binned", "Binned fit", 800, 600);
            frame_binned->Draw();

            std::string outdir = outfilePath;
            if (!outdir.empty() && outdir.back() != '/')
                outdir.push_back('/');
            std::string outpng = outdir + "mass_binned_fit_bin" + binstr + ".png";
            c_binned->SaveAs(outpng.c_str());

            if (sFile)
            {
                sFile->cd();
                c_binned->Write(("mass_binned_fit_bin" + binstr).c_str());
            }

            // keep the canvas in memory if user wants to view; don't delete immediately
        }

        // After binned fit, fix shape parameters and leave only yields free for the unbinned fit
        std::cout << "  Fixing shape parameters from binned fit and performing unbinned fit for yields..." << std::endl;

        // Depending on model type, set shape parameters constant
        mean->setConstant(true);
        sigma1->setConstant(true);
        deltasigma->setConstant(true);
        if (useCrystalBall)
        {
            alpha1->setConstant(true);
            n->setConstant(true);
        }
        if (useTwoComponentSignal)
        {
            component_frac->setConstant(true);
        }
        poly0->setConstant(true);
        poly1->setConstant(true);

        // Ensure yields are free to vary in the unbinned fit
        if (sig_yield)
            sig_yield->setConstant(false);
        if (bkg_yield)
            bkg_yield->setConstant(false);

        // Now perform the unbinned fit only varying yields (and any remaining free params)
        RooFitResult *fit_result = nullptr;
        if (fitTypeName == "noSig")
        {
            fit_result = extended_pdf->fitTo(*data, RooFit::Save(),
                                            RooFit::PrintLevel(0),
                                            RooFit::Extended(true),
                                            RooFit::NumCPU(8),
                                            RooFit::Strategy(1),
                                            RooFit::Range("SBleft,SBright"));
        }
        else
        {
            fit_result = extended_pdf->fitTo(*data, RooFit::Save(),
                                            RooFit::PrintLevel(0),
                                            RooFit::Extended(true),
                                            RooFit::NumCPU(8),
                                            RooFit::Strategy(1),
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

        // Create plot of the fit using the Plotter class (pass input file name for beam tagging)
        Plotter plotter(resonance, outfilePath, bin, false, zRange, inputFileName);
        histogram = plotter.individualMassFitPlotMulti(sig_yield, extended_pdf, mass_tag_measured, data, fitTypeName, isZtObservable);

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
        parameterArr[5] = useCrystalBall ? alpha1->getVal() : 0.0;
        parameterArr[6] = useCrystalBall ? n->getVal() : 0.0;
        parameterArr[7] = useTwoComponentSignal ? component_frac->getVal() : 0.0;
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
        parameterErrArr[5] = useCrystalBall ? alpha1->getError() : 0.0;
        parameterErrArr[6] = useCrystalBall ? n->getError() : 0.0;
        parameterErrArr[7] = useTwoComponentSignal ? component_frac->getError() : 0.0;
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
                            double massSigYield, double massSigYieldErr,
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
        RooRealVar *log_ipchi2 = new RooRealVar("log_tag_ipchi2", "log(tag_ipchi2)",
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

        // If a mass-fit signal yield was provided, use it to constrain the IP chi2 fit.
        // Implement constraint by allowing only ±10% variation around the mass-fit yield.
        if (massSigYield >= 0.0) {
            // Initialize to mass-fit value and set ±10% bounds
            sig_yieldLim->setVal(massSigYield);
            double lower = std::max(0.0, massSigYield * 0.90);
            double upper = massSigYield * 1.10;
            sig_yieldLim->setMin(lower);
            sig_yieldLim->setMax(upper);
            // Ensure the variable is free to vary within the ±10% window
            sig_yieldLim->setConstant(false);
            std::cout << "  Using mass-fit signal yield to initialize sig_yieldLim: " << massSigYield
                      << " (allowed range: " << lower << " - " << upper << ")" << std::endl;
        }

        // First: perform a binned fit on log_ipchi2 to determine shape parameters
        int binnedN_ip = nBins > 10 ? nBins : 50;
        std::cout << "  Performing binned IP chi2 fit (" << binnedN_ip << " bins) to determine shape parameters..." << std::endl;
        log_ipchi2->setBins(binnedN_ip);
        RooDataHist *binnedIPData = new RooDataHist("ip_binnedData", "ip_binnedData", RooArgSet(*log_ipchi2), *data);

        RooFitResult *binnedIPFitResult = nullptr;
        binnedIPFitResult = total_pdf->fitTo(*binnedIPData, RooFit::Save(),
                            RooFit::PrintLevel(0),
                            RooFit::Extended(true),
                            RooFit::NumCPU(8),
                            RooFit::Strategy(1),
                            RooFit::SumW2Error(true));

        // Diagnostic plot for binned IP chi2 fit
        {
            RooPlot *frame_ip_binned = log_ipchi2->frame(RooFit::Title("Binned IP chi2 fit"));
            binnedIPData->plotOn(frame_ip_binned, RooFit::DataError(RooAbsData::SumW2));
            total_pdf->plotOn(frame_ip_binned, RooFit::LineColor(kRed));
            if (binnedIPFitResult)
            {
                total_pdf->plotOn(frame_ip_binned, RooFit::VisualizeError(*binnedIPFitResult, 1), RooFit::DrawOption("F"), RooFit::FillColor(kOrange), RooFit::LineColor(kOrange));
            }
            TCanvas *c_ip_binned = new TCanvas("c_ip_binned", "IP binned fit", 800, 600);
            frame_ip_binned->Draw();
            std::string binstr = (bin >= 0) ? std::to_string(bin) : std::string("all");
            std::string outdir = outfilePath;
            if (!outdir.empty() && outdir.back() != '/') outdir.push_back('/');
            std::string outpng = outdir + "ipchi2_binned_fit_bin" + binstr + ".png";
            c_ip_binned->SaveAs(outpng.c_str());
            if (splotFile)
            {
                splotFile->cd();
                c_ip_binned->Write(("ipchi2_binned_fit_bin" + binstr).c_str());
            }
        }

        // Fix shape parameters from binned fit, leave yields/fractions free
        xp_prompt->setConstant(true);
        sigma_prompt->setConstant(true);
        xi_prompt->setConstant(true);
        rho1_prompt->setConstant(true);
        rho2_prompt->setConstant(true);

        xp_nonprompt->setConstant(true);
        sigma_nonprompt->setConstant(true);
        xi_nonprompt->setConstant(true);
        rho1_nonprompt->setConstant(true);
        rho2_nonprompt->setConstant(true);

        // Ensure sig_yieldLim and prompt_frac are free to vary for unbinned fit
        sig_yieldLim->setConstant(false);
        prompt_frac->setConstant(true);

        // Now perform the unbinned fit only varying yields/fraction
        RooFitResult *result = total_pdf->fitTo(*data, RooFit::Save(true),
                                                RooFit::PrintLevel(0),
                                                RooFit::SumW2Error(true),
                                                RooFit::Strategy(2),
            RooFit::Extended(true),
            RooFit::NumCPU(8),
                                                RooFit::Minos(false),
                                                RooFit::Hesse(true));

        std::cout << "  Fitted prompt_frac: " << prompt_frac->getVal() << " ± " << prompt_frac->getError() << std::endl;
        std::cout << "  Fitted sig_yieldLim: " << sig_yieldLim->getVal() << std::endl;
        std::cout << "  Fitted prompt_yield: " << prompt_yield->getVal() << std::endl;
        std::cout << "  Fitted nonprompt_yield: " << nonprompt_yield->getVal() << std::endl;

        // Create plot
        std::cout << "  Creating IP chi2 fit plot..." << std::endl;
        Plotter plotter(resonance, outfilePath, bin, false, zRange, inputFileName);
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
        parameterErrArr[5] = rho1_prompt->getError();
        parameterErrArr[6] = rho2_prompt->getError();

        parameterErrArr[7] = xp_nonprompt->getError();
        parameterErrArr[8] = sigma_nonprompt->getError();
        parameterErrArr[9] = xi_nonprompt->getError();
        parameterErrArr[10] = rho1_nonprompt->getError();
        parameterErrArr[11] = rho2_nonprompt->getError();

        std::cout << "  IP chi2 fit with yields completed successfully" << std::endl;
        std::cout << "  Results - Prompt fraction: " << prompt_frac->getVal() << " ± " << prompt_frac->getError() 
                  << ", Prompt yield: " << prompt_yield->getVal() 
                  << ", Non-prompt yield: " << nonprompt_yield->getVal() << std::endl;

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
        // Open sPlot file with retry mechanism
        TFile *splotFile = nullptr;
        for (int retryCount = 0; retryCount < 5; ++retryCount) {
            if (retryCount > 0) {
                std::this_thread::sleep_for(std::chrono::milliseconds(200));
                std::cout << "  Retry " << retryCount << " to open sPlot file..." << std::endl;
            }
            
            splotFile = new TFile(splotFileName.c_str(), "READ");
            if (splotFile && !splotFile->IsZombie()) break;
            
            delete splotFile;
            splotFile = nullptr;
        }

        if (!splotFile || splotFile->IsZombie()) {
            std::cerr << "Error: Cannot open sPlot file: " << splotFileName << std::endl;
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

        std::cout << "Creating weighted dataset: " << datasetName << std::endl;
        // Create argument set that includes weight variable
        RooArgSet datasetVars(*originalData->get());
        datasetVars.add(datasetWeight);

        // Create new weighted dataset with weight support
        RooDataSet *weightedData = new RooDataSet(datasetName.c_str(),
                                                  ("Weighted dataset: " + datasetName).c_str(),
                                                  datasetVars,
                                                  RooFit::WeightVar(datasetWeight));
        std::cout << "Initialized weighted dataset with weight variable" << std::endl;

        // Apply weights to the dataset
        int matchedEntries = 0;
        double totalWeight = 0.0;

        const double matchTolerance = 1e-6; // matching tolerance for mass values
        if (!massToWeight.empty()) {
            for (int i = 0; i < originalData->numEntries(); ++i)
            {
                const RooArgSet *row = originalData->get(i);
                RooRealVar *massVar = (RooRealVar *)row->find("tagMass");
                if (!massVar) continue;

                double mass = massVar->getVal();

                // Use lower_bound to find nearest mass key in O(log N)
                auto it = massToWeight.lower_bound(mass);
                double closestMass = 0.0;
                double minDiff = std::numeric_limits<double>::infinity();

                if (it != massToWeight.end()) {
                    double diff = std::abs(it->first - mass);
                    if (diff < minDiff) { minDiff = diff; closestMass = it->first; }
                }
                if (it != massToWeight.begin()) {
                    auto pit = std::prev(it);
                    double diff = std::abs(pit->first - mass);
                    if (diff < minDiff) { minDiff = diff; closestMass = pit->first; }
                }

                // Accept match only if within tolerance
                if (minDiff <= matchTolerance)
                {
                    double weight = massToWeight[closestMass];

                    // Only add entries with positive weights
                    if (weight > 0)
                    {
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
