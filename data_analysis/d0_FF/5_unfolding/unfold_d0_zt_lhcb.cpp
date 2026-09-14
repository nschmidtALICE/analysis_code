#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TMatrixD.h>
#include <TPad.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TTree.h>
#include <random>
#include <chrono>

#include <RooUnfoldBayes.h>
#include <RooUnfoldResponse.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
struct JetBin
{
    std::string label;
    double min = 0.0;
    double max = 0.0;
};

struct EtaBin
{
    int id = -1;
    double min = 0.0;
    double max = 0.0;
    std::string label;
};

struct D0UnfoldConfig
{
    std::string measuredFilePattern = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-04-07_pPb/TagZHistograms_%s.root";
    // std::string measuredFilePattern = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-04-01_pPb/TagZHistograms_%s.root";
    std::string responseFile = "/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pPb_MC_54plus73_EPOS_response.root";
    std::string triggeredResponseFile;
    std::string outputFile = "d0_unfolded_zt.root";
    std::string directHistogramPattern = "promptSignalTagZHist_FullyWeighted_%s";
    std::string perEtaHistogramStem = "promptSignalTagZHist_FullyWeighted_";
    std::string responseTreeName = "Response";
    std::vector<JetBin> jetBins;
    std::vector<int> etaBinIds = {0, 1, 2};
    // std::vector<int> etaBinIds = {0, 1, 2, 3, 4, 5, 6, 7};
    double etaMin = 2.5;
    double etaMax = 4.0;
    std::vector<std::string> rapidityBinLabels; // optional human-readable rapidity ranges per eta bin
    double d0Mass = 1.86484;
    double massWindow = 0.07;
    double maxJetDr = 0.4;
    double triggeredScale = 1.0;
    int minJetConstituents = 1;
    int nIterations = 6;
    bool includeFakes = false;
    bool doClosure = true;
    bool do2D = false;
    bool verbose = true;
    bool doResponseDetLevelWeighting = true;
    bool doResponseTruthLevelWeighting = true;
    bool applyResponseSmearing = true;
    double responseScaleFactor = 0.96;
    double responseSmearFactor = 0.07;
    // response bootstrap options
    bool doResponseBootstrap = true;
    int nResponseToys = 200;
    bool doPriorVariations = true;
    std::string outputFolderTag = ""; // optional tag appended to plot folder name
};

struct ResponseBundle
{
    std::unique_ptr<TH1D> recoTemplate;
    std::unique_ptr<TH1D> truthTemplate;
    std::unique_ptr<TH1D> truthTemplateUncut;
    std::unique_ptr<TH1D> missedTruth;
    std::unique_ptr<TH1D> fakeReco;
    std::unique_ptr<TH1D> measuredForUnfold;
    std::unique_ptr<TH1D> responseDetWeights;
    std::unique_ptr<TH1D> responseTruthWeights;
    std::unique_ptr<RooUnfoldResponse> response;
};

struct ResponseBundle2D
{
    std::unique_ptr<TH2D> recoTemplate;
    std::unique_ptr<TH2D> truthTemplate;
    std::unique_ptr<TH2D> truthTemplateUncut;
    std::unique_ptr<TH2D> missedTruth;
    std::unique_ptr<TH2D> fakeReco;
    std::unique_ptr<TH2D> measuredForUnfold;
    std::unique_ptr<TH2D> responseDetWeights;
    std::unique_ptr<TH2D> responseTruthWeights;
    std::unique_ptr<RooUnfoldResponse> response;
};

struct Chi2Result
{
    double chi2 = 0.0;
    int ndf = 0;
    double reduced = -1.0;
};

struct IterationScan1DResult
{
    int bestIteration = 1;
    double bestChi2 = std::numeric_limits<double>::infinity();
    std::unique_ptr<TH1D> bestUnfolded;
};

struct IterationScan2DResult
{
    int bestIteration = 1;
    double bestChi2 = std::numeric_limits<double>::infinity();
    std::unique_ptr<TH2D> bestUnfolded;
};

std::string ReplaceToken(const std::string &pattern, const std::string &replacement)
{
    const std::size_t pos = pattern.find("%s");
    if (pos == std::string::npos)
    {
        return pattern;
    }
    return pattern.substr(0, pos) + replacement + pattern.substr(pos + 2);
}

std::vector<std::string> SplitCsv(const std::string &csv)
{
    std::vector<std::string> tokens;
    std::stringstream stream(csv);
    std::string item;
    while (std::getline(stream, item, ','))
    {
        item.erase(std::remove_if(item.begin(), item.end(), ::isspace), item.end());
        if (!item.empty())
        {
            tokens.push_back(item);
        }
    }
    return tokens;
}

std::vector<JetBin> ParseJetBins(const std::string &jetBinsCsv)
{
    std::vector<JetBin> bins;
    for (const std::string &token : SplitCsv(jetBinsCsv))
    {
        const std::size_t separator = token.find('_');
        if (separator == std::string::npos)
        {
            throw std::runtime_error("Jet-bin token must be formatted as min_max: " + token);
        }
        JetBin bin;
        bin.label = token;
        bin.min = std::stod(token.substr(0, separator));
        bin.max = std::stod(token.substr(separator + 1));
        if (bin.max <= bin.min)
        {
            throw std::runtime_error("Jet-bin upper edge must exceed lower edge: " + token);
        }
        bins.push_back(bin);
    }
    return bins;
}

std::vector<double> BuildJetBinEdges(const std::vector<JetBin> &jetBins)
{
    std::vector<double> edges;
    edges.reserve(jetBins.size() + 1);
    if (jetBins.empty())
    {
        return edges;
    }
    edges.push_back(jetBins.front().min);
    for (const JetBin &bin : jetBins)
    {
        edges.push_back(bin.max);
    }
    return edges;
}

std::vector<EtaBin> BuildEtaBins(const D0UnfoldConfig &config)
{
    std::vector<EtaBin> etaBins;
    if (config.etaBinIds.empty())
    {
        return etaBins;
    }
    
    const double width = (config.etaMax - config.etaMin) / static_cast<double>(config.etaBinIds.size());
    etaBins.reserve(config.etaBinIds.size());
    for (std::size_t index = 0; index < config.etaBinIds.size(); ++index)
    {
        EtaBin etaBin;
        etaBin.id = config.etaBinIds[index];
        etaBin.min = config.etaMin + static_cast<double>(index) * width;
        etaBin.max = config.etaMin + static_cast<double>(index + 1) * width;
        if (config.rapidityBinLabels.size() == config.etaBinIds.size())
        {
            etaBin.label = config.rapidityBinLabels[index];
        }
        else
        {
            etaBin.label = "eta" + std::to_string(etaBin.id);
        }
        etaBins.push_back(etaBin);
    }
    return etaBins;
}

bool InEtaBin(double eta, const EtaBin &etaBin)
{
    return eta >= etaBin.min && eta < etaBin.max;
}

std::string EtaToken(const EtaBin &etaBin)
{
    return "eta" + std::to_string(etaBin.id);
}

std::string SanitizeLabel(std::string label)
{
    std::replace(label.begin(), label.end(), '.', 'p');
    std::replace(label.begin(), label.end(), '/', '_');
    return label;
}

/*
 * ApplyResponseRecoCorrection
 * ---------------------------
 * Purpose: optionally apply a detector-level smearing and scale to the
 * reconstructed jet pT (excluding the d0 fraction) and recompute the
 * corresponding d0 fraction so the response filling uses a reco value
 * consistent with the applied detector resolution model.
 *
 * Inputs:
 *  - config: controls `applyResponseSmearing`, `responseScaleFactor`,
 *            and `responseSmearFactor`.
 *  - jetPtDet: original reconstructed jet pT from the tree.
 *  - d0ZDet: original d0 fraction associated with the jet.
 *
 * Outputs (in-place):
 *  - jetPtDetCorr: corrected/smeared jet pT.
 *  - d0ZDetCorr: recomputed d0 fraction (kept only when jetPtDetCorr>0).
 *
 * Algorithm (when enabled):
 *  - detNoD0 = jetPtDet * (1 - d0ZDet)
 *  - mean = detNoD0 * responseScaleFactor
 *  - sigma = responseSmearFactor * responseScaleFactor * detNoD0
 *  - smeared = Gaussian(mean, sigma)  (uses ROOT `gRandom`)
 *  - jetPtDetCorr = smeared + (d0ZDet * jetPtDet)
 *  - d0ZDetCorr = (d0ZDet * jetPtDet) / jetPtDetCorr  (if jetPtDetCorr>0)
 *
 * Note: the function is stochastic when smearing is enabled; seed the
 * global RNG (`gRandom`) externally for reproducible results.
 */
void ApplyResponseRecoCorrection(const D0UnfoldConfig &config,
                                 double jetPtDet,
                                 double d0ZDet,
                                 double &jetPtDetCorr,
                                 double &d0ZDetCorr)
{
    jetPtDetCorr = jetPtDet;
    d0ZDetCorr = d0ZDet;
    if (!config.applyResponseSmearing)
    {
        return;
    }

    const double detNoD0 = jetPtDet * (1.0 - d0ZDet);
    const double mean = detNoD0 * config.responseScaleFactor;
    const double sigma = config.responseSmearFactor * config.responseScaleFactor * detNoD0;
    const double smeared = sigma > 0.0 ? gRandom->Gaus(mean, sigma) : mean;
    jetPtDetCorr = smeared + (d0ZDet * jetPtDet);
    if (jetPtDetCorr > 0.0)
    {
        d0ZDetCorr = (d0ZDet * jetPtDet) / jetPtDetCorr;
    }
}

void MapUnderOverflowToEdges(TH1 *hist)
{
    if (!hist)
    {
        return;
    }

    const int nbins = hist->GetNbinsX();
    if (hist->GetBinContent(0) != 0.0)
    {
        const double content = hist->GetBinContent(1) + hist->GetBinContent(0);
        const double error = std::hypot(hist->GetBinError(1), hist->GetBinError(0));
        hist->SetBinContent(1, content);
        hist->SetBinError(1, error);
        hist->SetBinContent(0, 0.0);
        hist->SetBinError(0, 0.0);
    }

    if (hist->GetBinContent(nbins + 1) != 0.0)
    {
        const double content = hist->GetBinContent(nbins) + hist->GetBinContent(nbins + 1);
        const double error = std::hypot(hist->GetBinError(nbins), hist->GetBinError(nbins + 1));
        hist->SetBinContent(nbins, content);
        hist->SetBinError(nbins, error);
        hist->SetBinContent(nbins + 1, 0.0);
        hist->SetBinError(nbins + 1, 0.0);
    }
}

void MapUnderOverflowToEdges(TH2 *hist)
{
    if (!hist)
    {
        return;
    }

    std::unique_ptr<TH2> copy(static_cast<TH2 *>(hist->Clone((std::string(hist->GetName()) + "_tmp").c_str())));
    copy->SetDirectory(nullptr);
    hist->Reset();

    const int nx = copy->GetNbinsX();
    const int ny = copy->GetNbinsY();
    for (int ix = 0; ix <= nx + 1; ++ix)
    {
        const int mappedX = std::min(std::max(ix, 1), nx);
        for (int iy = 0; iy <= ny + 1; ++iy)
        {
            const int mappedY = std::min(std::max(iy, 1), ny);
            const double oldContent = hist->GetBinContent(mappedX, mappedY);
            const double oldError = hist->GetBinError(mappedX, mappedY);
            const double addContent = copy->GetBinContent(ix, iy);
            const double addError = copy->GetBinError(ix, iy);
            hist->SetBinContent(mappedX, mappedY, oldContent + addContent);
            hist->SetBinError(mappedX, mappedY, std::hypot(oldError, addError));
        }
    }

    for (int ix = 0; ix <= nx + 1; ++ix)
    {
        hist->SetBinContent(ix, 0, 0.0);
        hist->SetBinContent(ix, ny + 1, 0.0);
        hist->SetBinError(ix, 0, 0.0);
        hist->SetBinError(ix, ny + 1, 0.0);
    }
    for (int iy = 0; iy <= ny + 1; ++iy)
    {
        hist->SetBinContent(0, iy, 0.0);
        hist->SetBinContent(nx + 1, iy, 0.0);
        hist->SetBinError(0, iy, 0.0);
        hist->SetBinError(nx + 1, iy, 0.0);
    }
}

Chi2Result ComputeChi2(const TH1 *data, const TH1 *model)
{
    Chi2Result result;
    if (!data || !model)
    {
        return result;
    }

    const int nbins = data->GetNbinsX();
    for (int bin = 1; bin <= nbins; ++bin)
    {
        const double error = data->GetBinError(bin);
        if (error <= 0.0)
        {
            continue;
        }
        const double delta = data->GetBinContent(bin) - model->GetBinContent(bin);
        result.chi2 += (delta * delta) / (error * error);
        ++result.ndf;
    }

    if (result.ndf > 0)
    {
        result.reduced = result.chi2 / static_cast<double>(result.ndf);
    }
    return result;
}

TH2D *BuildCovarianceHistogram(const TMatrixD &covariance, const std::string &name)
{
    const int nrows = covariance.GetNrows();
    const int ncols = covariance.GetNcols();
    TH2D *hist = new TH2D(name.c_str(), (name + ";bin;bin").c_str(), nrows, 0.5, nrows + 0.5, ncols, 0.5, ncols + 0.5);
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            hist->SetBinContent(row + 1, col + 1, covariance(row, col));
        }
    }
    return hist;
}

TH2D *BuildCorrelationHistogram(const TMatrixD &covariance, const std::string &name)
{
    const int nrows = covariance.GetNrows();
    const int ncols = covariance.GetNcols();
    TH2D *hist = new TH2D(name.c_str(), (name + ";bin;bin").c_str(), nrows, 0.5, nrows + 0.5, ncols, 0.5, ncols + 0.5);
    for (int row = 0; row < nrows; ++row)
    {
        const double sigmaRow = covariance(row, row) > 0.0 ? std::sqrt(covariance(row, row)) : 0.0;
        for (int col = 0; col < ncols; ++col)
        {
            const double sigmaCol = covariance(col, col) > 0.0 ? std::sqrt(covariance(col, col)) : 0.0;
            const double correlation = (sigmaRow > 0.0 && sigmaCol > 0.0) ? covariance(row, col) / (sigmaRow * sigmaCol) : 0.0;
            hist->SetBinContent(row + 1, col + 1, correlation);
        }
    }
    return hist;
}

// Build covariance from a vector of histograms (toys)
TMatrixD BuildCovarianceFromToys(const std::vector<std::unique_ptr<TH1D>> &toys)
{
    if (toys.empty())
    {
        return TMatrixD();
    }
    const int nbins = toys.front()->GetNbinsX();
    TMatrixD cov(nbins, nbins);
    std::vector<double> mean(nbins, 0.0);
    const int N = static_cast<int>(toys.size());
    for (const auto &h : toys)
    {
        for (int b = 1; b <= nbins; ++b)
        {
            mean[b - 1] += h->GetBinContent(b);
        }
    }
    for (int i = 0; i < nbins; ++i) mean[i] /= static_cast<double>(N);
    for (const auto &h : toys)
    {
        for (int i = 0; i < nbins; ++i)
        {
            const double vi = h->GetBinContent(i + 1) - mean[i];
            for (int j = 0; j < nbins; ++j)
            {
                const double vj = h->GetBinContent(j + 1) - mean[j];
                cov(i, j) += vi * vj;
            }
        }
    }
    if (N > 1)
    {
        cov *= 1.0 / static_cast<double>(N - 1);
    }
    return cov;
}

std::vector<double> ExtractAxisEdges(const TAxis *axis)
{
    std::vector<double> edges(axis->GetNbins() + 1, 0.0);
    for (int i = 1; i <= axis->GetNbins(); ++i)
    {
        edges[i - 1] = axis->GetBinLowEdge(i);
    }
    edges.back() = axis->GetBinUpEdge(axis->GetNbins());
    return edges;
}

void SaveCanvas(TCanvas *canvas, const std::filesystem::path &path)
{
    if (!canvas)
    {
        return;
    }
    std::filesystem::create_directories(path.parent_path());
    canvas->SaveAs(path.string().c_str());
}

TH1D *MakeRatioHistogram(const TH1 *numerator, const TH1 *denominator, const std::string &name, const std::string &yTitle)
{
    TH1D *ratio = static_cast<TH1D *>(numerator->Clone(name.c_str()));
    ratio->SetDirectory(nullptr);
    ratio->Reset();
    ratio->SetTitle((name + ";z_{T};" + yTitle).c_str());
    for (int bin = 1; bin <= numerator->GetNbinsX(); ++bin)
    {
        const double num = numerator->GetBinContent(bin);
        const double den = denominator->GetBinContent(bin);
        if (den == 0.0)
        {
            continue;
        }
        ratio->SetBinContent(bin, num / den);
        ratio->SetBinError(bin, denominator->GetBinError(bin) > 0.0 ? numerator->GetBinError(bin) / den : 0.0);
    }
    return ratio;
}

TH2D *MakeRatioHistogram2D(const TH2 *numerator, const TH2 *denominator, const std::string &name, const std::string &zTitle)
{
    TH2D *ratio = static_cast<TH2D *>(numerator->Clone(name.c_str()));
    ratio->SetDirectory(nullptr);
    ratio->Reset();
    ratio->SetTitle((name + ";z_{T};p_{T}^{jet} [GeV/c];" + zTitle).c_str());
    for (int xbin = 1; xbin <= numerator->GetNbinsX(); ++xbin)
    {
        for (int ybin = 1; ybin <= numerator->GetNbinsY(); ++ybin)
        {
            const double den = denominator->GetBinContent(xbin, ybin);
            if (den == 0.0)
            {
                continue;
            }
            ratio->SetBinContent(xbin, ybin, numerator->GetBinContent(xbin, ybin) / den);
        }
    }
    return ratio;
}

std::unique_ptr<TH1D> BuildRelativeWeightHistogram1D(const TH1 &target, const TH1 &reference, const std::string &name)
{
    std::unique_ptr<TH1D> weights(static_cast<TH1D *>(target.Clone(name.c_str())));
    weights->SetDirectory(nullptr);
    weights->Reset();

    double maxWeight = 0.0;
    for (int bin = 1; bin <= target.GetNbinsX(); ++bin)
    {
        const double den = reference.GetBinContent(bin);
        if (den <= 0.0)
        {
            continue;
        }
        const double value = target.GetBinContent(bin) / den;
        weights->SetBinContent(bin, value);
        maxWeight = std::max(maxWeight, value);
    }

    if (maxWeight > 0.0)
    {
        weights->Scale(1.0 / maxWeight);
    }
    return weights;
}

std::unique_ptr<TH2D> BuildRelativeWeightHistogram2D(const TH2 &target, const TH2 &reference, const std::string &name)
{
    std::unique_ptr<TH2D> weights(static_cast<TH2D *>(target.Clone(name.c_str())));
    weights->SetDirectory(nullptr);
    weights->Reset();

    double maxWeight = 0.0;
    for (int xbin = 1; xbin <= target.GetNbinsX(); ++xbin)
    {
        for (int ybin = 1; ybin <= target.GetNbinsY(); ++ybin)
        {
            const double den = reference.GetBinContent(xbin, ybin);
            if (den <= 0.0)
            {
                continue;
            }
            const double value = target.GetBinContent(xbin, ybin) / den;
            weights->SetBinContent(xbin, ybin, value);
            maxWeight = std::max(maxWeight, value);
        }
    }

    if (maxWeight > 0.0)
    {
        weights->Scale(1.0 / maxWeight);
    }
    return weights;
}

double LookupWeight1D(const TH1 *weights, double value)
{
    if (!weights)
    {
        return 1.0;
    }
    int bin = weights->GetXaxis()->FindBin(value);
    bin = std::max(1, std::min(bin, weights->GetNbinsX()));
    const double weight = weights->GetBinContent(bin);
    return weight > 0.0 ? weight : 1.0;
}

double LookupWeight2D(const TH2 *weights, double xValue, double yValue)
{
    if (!weights)
    {
        return 1.0;
    }
    int xbin = weights->GetXaxis()->FindBin(xValue);
    int ybin = weights->GetYaxis()->FindBin(yValue);
    xbin = std::max(1, std::min(xbin, weights->GetNbinsX()));
    ybin = std::max(1, std::min(ybin, weights->GetNbinsY()));
    const double weight = weights->GetBinContent(xbin, ybin);
    return weight > 0.0 ? weight : 1.0;
}

class D0ZTUnfolder
{
public:
    explicit D0ZTUnfolder(D0UnfoldConfig config)
        : config_(std::move(config))
    {
    }

    void Run()
    {
        PrepareOutputLayout();
        OpenFiles();

        std::vector<double> jetEdges = BuildJetBinEdges(config_.jetBins);
        const std::vector<EtaBin> etaBins = BuildEtaBins(config_);

        for (const EtaBin &etaBin : etaBins)
        {
            if (config_.verbose)
            {
                std::cout << "[unfold] Processing eta bin " << etaBin.id << " in range ["
                          << etaBin.min << ", " << etaBin.max << ")" << std::endl;
            }

            TH2D *measuredSummary = nullptr;
            TH2D *bestUnfoldedSummary = nullptr;
            std::vector<std::unique_ptr<TH1D>> kinEffHists;

            for (std::size_t jetIndex = 0; jetIndex < config_.jetBins.size(); ++jetIndex)
            {
                const JetBin &jetBin = config_.jetBins[jetIndex];
                if (config_.verbose)
                {
                    std::cout << "[unfold]   jet bin " << jetBin.label << std::endl;
                }

                std::unique_ptr<TH1D> measured = LoadMeasuredSpectrum(jetBin, etaBin);
                if (!measured)
                {
                    std::cerr << "[unfold] Skipping jet bin " << jetBin.label << " in eta bin " << etaBin.id
                              << " because no measured histogram was found." << std::endl;
                    continue;
                }

                PlotMeasuredSpectrum(jetBin, etaBin, *measured);
                MapUnderOverflowToEdges(measured.get());
                ResponseBundle bundle = BuildResponse(jetBin, etaBin, *measured);
                WriteInputs(jetBin, etaBin, *measured, bundle);

                TDirectory *etaDir = outputFile_->GetDirectory(etaBin.label.c_str());
                if (!etaDir)
                {
                    etaDir = outputFile_->mkdir(etaBin.label.c_str());
                }
                etaDir->cd();

                const std::string dirName = "jetpt_" + SanitizeLabel(jetBin.label);
                TDirectory *jetDir = etaDir->GetDirectory(dirName.c_str());
                if (!jetDir)
                {
                    jetDir = etaDir->mkdir(dirName.c_str());
                }
                jetDir->cd();

                IterationScan1DResult scanResult;
                if (config_.doResponseDetLevelWeighting)
                {
                    const TH1 *responseMeasuredNominal = bundle.response->Hmeasured();
                    std::unique_ptr<TH1D> responseMeasured1D;
                    if (responseMeasuredNominal)
                    {
                        responseMeasured1D.reset(static_cast<TH1D *>(responseMeasuredNominal->Clone(("responseMeasured_nominal_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
                        responseMeasured1D->SetDirectory(nullptr);
                        MapUnderOverflowToEdges(responseMeasured1D.get());
                    }

                    std::unique_ptr<TH1D> detWeights = responseMeasured1D ? BuildRelativeWeightHistogram1D(*bundle.measuredForUnfold, *responseMeasured1D,
                                                                                                             "responseDetWeights_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin))
                                                                           : nullptr;
                    if (detWeights)
                    {
                        outputFile_->cd();
                        detWeights->Write();
                    }

                    ResponseBundle weightedBundle = BuildWeightedResponse(jetBin, etaBin, *measured, detWeights.get(), nullptr);
                    scanResult = ScanIterations1D(jetBin, etaBin, weightedBundle);
                    bundle = std::move(weightedBundle);
                }
                else
                {
                    scanResult = ScanIterations1D(jetBin, etaBin, bundle);
                }

                if (config_.doResponseTruthLevelWeighting && scanResult.bestUnfolded)
                {
                    std::unique_ptr<TH1D> truthWeights = BuildRelativeWeightHistogram1D(*scanResult.bestUnfolded, *bundle.truthTemplate,
                                                                                        "responseTruthWeights_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin));
                    if (truthWeights)
                    {
                        outputFile_->cd();
                        truthWeights->Write();
                    }
                    ResponseBundle truthWeightedBundle = BuildWeightedResponse(jetBin, etaBin, *measured, nullptr, truthWeights.get());
                    scanResult = ScanIterations1D(jetBin, etaBin, truthWeightedBundle);
                    bundle = std::move(truthWeightedBundle);
                }

                // Compute kinematic efficiency versus zT for this jet-pt bin
                if (bundle.truthTemplate && bundle.truthTemplateUncut)
                {
                    std::string effName = "kinEff_zT_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin);
                    TH1D *effHist = MakeRatioHistogram(bundle.truthTemplate.get(), bundle.truthTemplateUncut.get(), effName, "efficiency");
                    if (effHist)
                    {
                        effHist->SetDirectory(nullptr);
                        // style and store for overlay
                        const int colors[] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 2, kOrange + 7, kCyan + 2};
                        const int ncolors = sizeof(colors) / sizeof(colors[0]);
                        effHist->SetLineColor(colors[static_cast<int>(jetIndex) % ncolors]);
                        effHist->SetMarkerColor(colors[static_cast<int>(jetIndex) % ncolors]);
                        effHist->SetLineWidth(2);
                        outputFile_->cd();
                        effHist->Write(effName.c_str());
                        kinEffHists.emplace_back(std::unique_ptr<TH1D>(effHist));
                    }
                }

                std::cout << "[unfold] Best iteration for jet bin " << jetBin.label << " in eta bin " << etaBin.id
                          << " is " << scanResult.bestIteration << " with reduced chi2 = " << scanResult.bestChi2 << std::endl;

                if (!scanResult.bestUnfolded)
                {
                    std::cerr << "[unfold] No valid unfolding result retained for jet bin " << jetBin.label
                              << " in eta bin " << etaBin.id << std::endl;
                    continue;
                }

                std::unique_ptr<TH1D> bestUnfolded(static_cast<TH1D *>(scanResult.bestUnfolded->Clone(("best_unfolded_zT_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
                bestUnfolded->SetDirectory(nullptr);
                bestUnfolded->Write();

                // Optional prior variations controlled by config flag
                if (config_.doPriorVariations)
                {
                    std::cout << "[unfold] Performing prior variations for jet bin " << jetBin.label << " in eta bin " << etaBin.id << std::endl;
                    PerformPriorVariations(jetBin, etaBin, bundle, scanResult.bestIteration, *bestUnfolded);
                }

                // Optionally compute response-bootstrap uncertainty and write covariances
                if (config_.doResponseBootstrap)
                {
                    ComputeResponseBootstrapCovariance(jetBin, etaBin, bundle, scanResult.bestIteration, *bestUnfolded);
                }

                if (config_.doClosure)
                {
                    RunClosure(jetBin, etaBin, bundle, scanResult.bestIteration);
                }

                if (!measuredSummary)
                {
                    const std::vector<double> zEdges = ExtractAxisEdges(measured->GetXaxis());
                    measuredSummary = new TH2D(("h2_measured_zT_vs_jetPt_" + EtaToken(etaBin)).c_str(),
                                               ("Measured D^{0} z_{T}, " + etaBin.label + ";z_{T};p_{T}^{jet} [GeV/c]").c_str(),
                                               static_cast<int>(zEdges.size()) - 1,
                                               zEdges.data(),
                                               static_cast<int>(jetEdges.size()) - 1,
                                               jetEdges.data());
                    bestUnfoldedSummary = new TH2D(("h2_bestUnfolded_zT_vs_jetPt_" + EtaToken(etaBin)).c_str(),
                                                   ("Best unfolded D^{0} z_{T}, " + etaBin.label + ";z_{T};p_{T}^{jet} [GeV/c]").c_str(),
                                                   static_cast<int>(zEdges.size()) - 1,
                                                   zEdges.data(),
                                                   static_cast<int>(jetEdges.size()) - 1,
                                                   jetEdges.data());
                    // per-jet 1D efficiency histograms will be collected in kinEffHists
                }

                if (measuredSummary && bestUnfoldedSummary && measured->GetNbinsX() == measuredSummary->GetNbinsX())
                {
                    for (int bin = 1; bin <= measured->GetNbinsX(); ++bin)
                    {
                        measuredSummary->SetBinContent(bin, static_cast<int>(jetIndex) + 1, measured->GetBinContent(bin));
                        measuredSummary->SetBinError(bin, static_cast<int>(jetIndex) + 1, measured->GetBinError(bin));
                        bestUnfoldedSummary->SetBinContent(bin, static_cast<int>(jetIndex) + 1, bestUnfolded->GetBinContent(bin));
                        bestUnfoldedSummary->SetBinError(bin, static_cast<int>(jetIndex) + 1, bestUnfolded->GetBinError(bin));
                    }
                }
            }

            outputFile_->cd();
            if (measuredSummary)
            {
                measuredSummary->Write();
                bestUnfoldedSummary->Write();
                // If per-jet efficiency histograms were accumulated, draw them together
                if (!kinEffHists.empty())
                {
                    // Styling colors
                    const int colors[] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 2, kOrange + 7, kCyan + 2};
                    const int ncolors = sizeof(colors) / sizeof(colors[0]);
                    TCanvas ckin(("c_kinEff_" + EtaToken(etaBin)).c_str(), "kinematic efficiency", 900, 650);
                    kinEffHists[0]->SetLineWidth(2);
                    kinEffHists[0]->SetLineColor(colors[0]);
                    kinEffHists[0]->SetMarkerColor(colors[0]);
                    kinEffHists[0]->Draw("E1");
                    TLegend leg(0.65, 0.7, 0.9, 0.9);
                    leg.SetBorderSize(0);
                    leg.AddEntry(kinEffHists[0].get(), kinEffHists[0]->GetName(), "lep");
                    for (std::size_t i = 1; i < kinEffHists.size(); ++i)
                    {
                        kinEffHists[i]->SetLineWidth(2);
                        kinEffHists[i]->SetLineColor(colors[i % ncolors]);
                        kinEffHists[i]->SetMarkerColor(colors[i % ncolors]);
                        kinEffHists[i]->Draw("E1 SAME");
                        leg.AddEntry(kinEffHists[i].get(), kinEffHists[i]->GetName(), "lep");
                    }
                    leg.Draw();
                    SaveCanvas(&ckin, plotDirectory_ / etaBin.label / ("kinematic_efficiency_overlay_" + EtaToken(etaBin) + ".png"));
                }
                PlotSummary(etaBin, *measuredSummary, *bestUnfoldedSummary);
                if (config_.do2D)
                {
                    Run2D(etaBin, *measuredSummary);
                }
            }
        }
    }

private:
    void OpenFiles()
    {
        responseInput_.reset(TFile::Open(config_.responseFile.c_str(), "READ"));
        if (!responseInput_ || responseInput_->IsZombie())
        {
            throw std::runtime_error("Failed to open response file: " + config_.responseFile);
        }

        responseTree_ = dynamic_cast<TTree *>(responseInput_->Get(config_.responseTreeName.c_str()));
        if (!responseTree_)
        {
            throw std::runtime_error("Response tree not found: " + config_.responseTreeName);
        }

        if (!config_.triggeredResponseFile.empty())
        {
            triggeredResponseInput_.reset(TFile::Open(config_.triggeredResponseFile.c_str(), "READ"));
            if (!triggeredResponseInput_ || triggeredResponseInput_->IsZombie())
            {
                throw std::runtime_error("Failed to open triggered response file: " + config_.triggeredResponseFile);
            }
            triggeredResponseTree_ = dynamic_cast<TTree *>(triggeredResponseInput_->Get(config_.responseTreeName.c_str()));
            if (!triggeredResponseTree_)
            {
                throw std::runtime_error("Triggered response tree not found: " + config_.responseTreeName);
            }
        }

        outputFile_.reset(TFile::Open(config_.outputFile.c_str(), "RECREATE"));
        if (!outputFile_ || outputFile_->IsZombie())
        {
            throw std::runtime_error("Failed to create output file: " + config_.outputFile);
        }
    }

    void FillTruthProjectionUncut(TTree *tree,
                                  const JetBin &jetBin,
                                  const EtaBin &etaBin,
                                  double scale,
                                  TH1D *truthUncut,
                                  const TH1 *detLevelWeights = nullptr,
                                  const TH1 *truthLevelWeights = nullptr) const;

    void FillTruthProjectionUncut2D(TTree *tree,
                                    const EtaBin &etaBin,
                                    double scale,
                                    TH2D *truthUncut,
                                    const TH2 *detLevelWeights = nullptr,
                                    const TH2 *truthLevelWeights = nullptr) const;

    void PrepareOutputLayout()
    {
        const std::filesystem::path outputPath(config_.outputFile);
        if (!outputPath.parent_path().empty())
        {
            std::filesystem::create_directories(outputPath.parent_path());
        }

        // Build a dated plot directory: <outputStem>_plots_<YYYY-MM-DD>
        std::time_t t = std::time(nullptr);
        std::tm tm = *std::localtime(&t);
        char dateBuf[11] = {0};
        std::strftime(dateBuf, sizeof(dateBuf), "%Y-%m-%d", &tm);
        const std::string dateStr(dateBuf);

        std::string plotStem = outputPath.stem().string() + "_plots_" + dateStr;
        if (!config_.outputFolderTag.empty()) {
            // append user tag to help distinguish systematic variations
            plotStem += std::string("_") + config_.outputFolderTag;
        }
        if (!outputPath.parent_path().empty())
        {
            plotDirectory_ = outputPath.parent_path() / plotStem;
        }
        else
        {
            plotDirectory_ = plotStem;
        }
        std::filesystem::create_directories(plotDirectory_);

        // Move the output ROOT file into the plot directory so it is saved
        // alongside the generated plot PNGs.
        const std::filesystem::path originalOutputPath(config_.outputFile);
        const std::string outputFilename = originalOutputPath.filename().string();
        config_.outputFile = (plotDirectory_ / outputFilename).string();
    }

    std::unique_ptr<TH1D> LoadMeasuredSpectrum(const JetBin &jetBin, const EtaBin &etaBin) const
    {
        const std::string measuredFileName = ReplaceToken(config_.measuredFilePattern, jetBin.label);
        std::unique_ptr<TFile> measuredFile(TFile::Open(measuredFileName.c_str(), "READ"));
        if (!measuredFile || measuredFile->IsZombie())
        {
            std::cerr << "[input] Cannot open measured file: " << measuredFileName << std::endl;
            return nullptr;
        }

        const std::string etaHistToken = jetBin.label + "_bin" + std::to_string(etaBin.id);
        const std::string directName = ReplaceToken(config_.directHistogramPattern, etaHistToken);
        if (TH1D *directHist = dynamic_cast<TH1D *>(measuredFile->Get(directName.c_str())))
        {
            std::unique_ptr<TH1D> result(static_cast<TH1D *>(directHist->Clone(("measured_zT_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
            result->SetDirectory(nullptr);
            std::cout << "[input] Loaded measured histogram for jet bin " << jetBin.label << " eta bin " << etaBin.id << std::endl;
            return result;
        }

        const std::string histName = config_.perEtaHistogramStem + etaHistToken;
        TH1D *hist = dynamic_cast<TH1D *>(measuredFile->Get(histName.c_str()));
        if (!hist)
        {
            std::cout << "[input] Cannot find histogram for eta bin " << etaBin.id << " in jet bin " << jetBin.label << std::endl;
            return nullptr;
        }

        std::unique_ptr<TH1D> result(static_cast<TH1D *>(hist->Clone(("measured_zT_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
        result->SetDirectory(nullptr);
        return result;
    }

    ResponseBundle BuildResponse(const JetBin &jetBin, const EtaBin &etaBin, const TH1D &measured)
    {
        ResponseBundle bundle;
        bundle.recoTemplate.reset(static_cast<TH1D *>(measured.Clone(("recoTemplate_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
        bundle.truthTemplate.reset(static_cast<TH1D *>(measured.Clone(("truthTemplate_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
        bundle.missedTruth.reset(static_cast<TH1D *>(measured.Clone(("missedTruth_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
        bundle.fakeReco.reset(static_cast<TH1D *>(measured.Clone(("fakeReco_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
        bundle.recoTemplate->Reset();
        bundle.truthTemplate->Reset();
        bundle.missedTruth->Reset();
        bundle.fakeReco->Reset();

        // Build an "uncut" truth projection (kinematic efficiency) and use it
        // when constructing the RooUnfoldResponse so the response carries
        // the gen-level spectrum before detector-level pT cuts (kinematic eff).
        bundle.truthTemplateUncut.reset(static_cast<TH1D *>(measured.Clone(("truthTemplateUncut_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
        bundle.truthTemplateUncut->Reset();
        // fill uncut truth (no jet pT range requirements, but keep eta/nconst filters)
        FillTruthProjectionUncut(responseTree_, jetBin, etaBin, 1.0, bundle.truthTemplateUncut.get(), nullptr, nullptr);

        if (triggeredResponseTree_)
        {
            FillTruthProjectionUncut(triggeredResponseTree_, jetBin, etaBin, config_.triggeredScale, bundle.truthTemplateUncut.get(), nullptr, nullptr);
        }

        bundle.response = std::make_unique<RooUnfoldResponse>(bundle.recoTemplate.get(), bundle.truthTemplateUncut.get());

        // Now fill the response and the (cut) truth/template diagnostics as before
        FillResponseTree(responseTree_, jetBin, etaBin, 1.0, *bundle.response, *bundle.truthTemplate, *bundle.missedTruth, *bundle.fakeReco);
        if (triggeredResponseTree_)
        {
            FillResponseTree(triggeredResponseTree_, jetBin, etaBin, config_.triggeredScale, *bundle.response, *bundle.truthTemplate, *bundle.missedTruth, *bundle.fakeReco);
        }

        if (TH2 *matrix = bundle.response->Hresponse())
        {
            MapUnderOverflowToEdges(matrix);
        }
        MapUnderOverflowToEdges(bundle.truthTemplate.get());
        MapUnderOverflowToEdges(bundle.missedTruth.get());
        MapUnderOverflowToEdges(bundle.fakeReco.get());

        bundle.measuredForUnfold.reset(static_cast<TH1D *>(measured.Clone(("measuredForUnfold_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
        bundle.measuredForUnfold->SetDirectory(nullptr);
        const TH1 *responseMeasured = bundle.response->Hmeasured();
        const double dataIntegral = bundle.measuredForUnfold->Integral();
        const double responseIntegral = responseMeasured ? responseMeasured->Integral() : 0.0;
        if (dataIntegral > 0.0 && responseIntegral > 0.0)
        {
            bundle.measuredForUnfold->Scale(responseIntegral / dataIntegral);
        }

        ReportResponseDiagnostics(jetBin, etaBin, measured, bundle, dataIntegral, responseIntegral);

        return bundle;
    }

    ResponseBundle BuildWeightedResponse(const JetBin &jetBin,
                                        const EtaBin &etaBin,
                                        const TH1D &measured,
                                        const TH1 *detLevelWeights,
                                        const TH1 *truthLevelWeights)
    {
        ResponseBundle bundle;
        bundle.recoTemplate.reset(static_cast<TH1D *>(measured.Clone(("recoTemplateWeighted_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
        bundle.truthTemplate.reset(static_cast<TH1D *>(measured.Clone(("truthTemplateWeighted_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
        bundle.missedTruth.reset(static_cast<TH1D *>(measured.Clone(("missedTruthWeighted_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
        bundle.fakeReco.reset(static_cast<TH1D *>(measured.Clone(("fakeRecoWeighted_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
        bundle.recoTemplate->Reset();
        bundle.truthTemplate->Reset();
        bundle.missedTruth->Reset();
        bundle.fakeReco->Reset();

        if (detLevelWeights)
        {
            bundle.responseDetWeights.reset(static_cast<TH1D *>(detLevelWeights->Clone(("appliedDetWeights_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
            bundle.responseDetWeights->SetDirectory(nullptr);
        }
        if (truthLevelWeights)
        {
            bundle.responseTruthWeights.reset(static_cast<TH1D *>(truthLevelWeights->Clone(("appliedTruthWeights_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
            bundle.responseTruthWeights->SetDirectory(nullptr);
        }

        // Build uncut truth projection (kinematic efficiency) taking into account any
        // det/truth-level weightings provided so the RooUnfoldResponse truth
        // reflects the gen-level counts before det-level pT selection.
        bundle.truthTemplateUncut.reset(static_cast<TH1D *>(measured.Clone(("truthTemplateUncutWeighted_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
        bundle.truthTemplateUncut->Reset();
        FillTruthProjectionUncut(responseTree_, jetBin, etaBin, 1.0, bundle.truthTemplateUncut.get(), detLevelWeights, truthLevelWeights);
        if (triggeredResponseTree_)
        {
            FillTruthProjectionUncut(triggeredResponseTree_, jetBin, etaBin, config_.triggeredScale, bundle.truthTemplateUncut.get(), detLevelWeights, truthLevelWeights);
        }

        bundle.response = std::make_unique<RooUnfoldResponse>(bundle.recoTemplate.get(), bundle.truthTemplateUncut.get());

        FillResponseTree(responseTree_, jetBin, etaBin, 1.0, *bundle.response, *bundle.truthTemplate, *bundle.missedTruth, *bundle.fakeReco,
                         detLevelWeights, truthLevelWeights);
        if (triggeredResponseTree_)
        {
            FillResponseTree(triggeredResponseTree_, jetBin, etaBin, config_.triggeredScale, *bundle.response, *bundle.truthTemplate, *bundle.missedTruth, *bundle.fakeReco,
                             detLevelWeights, truthLevelWeights);
        }

        if (TH2 *matrix = bundle.response->Hresponse())
        {
            MapUnderOverflowToEdges(matrix);
        }
        MapUnderOverflowToEdges(bundle.truthTemplate.get());
        MapUnderOverflowToEdges(bundle.missedTruth.get());
        MapUnderOverflowToEdges(bundle.fakeReco.get());

        bundle.measuredForUnfold.reset(static_cast<TH1D *>(measured.Clone(("measuredForUnfoldWeighted_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
        bundle.measuredForUnfold->SetDirectory(nullptr);
        ReportResponseDiagnostics(jetBin, etaBin, measured, bundle, bundle.measuredForUnfold->Integral(), bundle.response->Hmeasured() ? bundle.response->Hmeasured()->Integral() : 0.0);

        return bundle;
    }

    ResponseBundle2D BuildResponse2D(const EtaBin &etaBin, const TH2D &measured2D)
    {
        ResponseBundle2D bundle;
        bundle.recoTemplate.reset(static_cast<TH2D *>(measured2D.Clone(("recoTemplate2D_" + EtaToken(etaBin)).c_str())));
        bundle.truthTemplate.reset(static_cast<TH2D *>(measured2D.Clone(("truthTemplate2D_" + EtaToken(etaBin)).c_str())));
        bundle.missedTruth.reset(static_cast<TH2D *>(measured2D.Clone(("missedTruth2D_" + EtaToken(etaBin)).c_str())));
        bundle.fakeReco.reset(static_cast<TH2D *>(measured2D.Clone(("fakeReco2D_" + EtaToken(etaBin)).c_str())));
        bundle.recoTemplate->Reset();
        bundle.truthTemplate->Reset();
        bundle.missedTruth->Reset();
        bundle.fakeReco->Reset();

        // Build an uncut truth projection for kinematic-efficiency-aware response
        bundle.truthTemplateUncut.reset(static_cast<TH2D *>(measured2D.Clone(("truthTemplateUncut2D_" + EtaToken(etaBin)).c_str())));
        bundle.truthTemplateUncut->Reset();
        FillTruthProjectionUncut2D(responseTree_, etaBin, 1.0, bundle.truthTemplateUncut.get(), nullptr, nullptr);
        if (triggeredResponseTree_)
        {
            FillTruthProjectionUncut2D(triggeredResponseTree_, etaBin, config_.triggeredScale, bundle.truthTemplateUncut.get(), nullptr, nullptr);
        }

        bundle.response = std::make_unique<RooUnfoldResponse>(bundle.recoTemplate.get(), bundle.truthTemplateUncut.get());

        FillResponseTree2D(responseTree_, etaBin, 1.0, *bundle.response, *bundle.truthTemplate, *bundle.missedTruth, *bundle.fakeReco);
        if (triggeredResponseTree_)
        {
            FillResponseTree2D(triggeredResponseTree_, etaBin, config_.triggeredScale, *bundle.response, *bundle.truthTemplate, *bundle.missedTruth, *bundle.fakeReco);
        }

        if (TH2 *matrix = bundle.response->Hresponse())
        {
            MapUnderOverflowToEdges(matrix);
        }
        MapUnderOverflowToEdges(bundle.truthTemplate.get());
        MapUnderOverflowToEdges(bundle.missedTruth.get());
        MapUnderOverflowToEdges(bundle.fakeReco.get());

        bundle.measuredForUnfold.reset(static_cast<TH2D *>(measured2D.Clone(("measuredForUnfold2D_" + EtaToken(etaBin)).c_str())));
        bundle.measuredForUnfold->SetDirectory(nullptr);
        MapUnderOverflowToEdges(bundle.measuredForUnfold.get());

        const TH1 *responseMeasured = bundle.response->Hmeasured();
        const double dataIntegral = bundle.measuredForUnfold->Integral();
        const double responseIntegral = responseMeasured ? responseMeasured->Integral() : 0.0;
        // if (dataIntegral > 0.0 && responseIntegral > 0.0)
        // {
        //     bundle.measuredForUnfold->Scale(responseIntegral / dataIntegral);
        // }

        ReportResponseDiagnostics2D(etaBin, measured2D, bundle, dataIntegral, responseIntegral);

        return bundle;
    }

    ResponseBundle2D BuildWeightedResponse2D(const EtaBin &etaBin,
                                            const TH2D &measured2D,
                                            const TH2 *detLevelWeights,
                                            const TH2 *truthLevelWeights)
    {
        ResponseBundle2D bundle;
        bundle.recoTemplate.reset(static_cast<TH2D *>(measured2D.Clone(("recoTemplateWeighted2D_" + EtaToken(etaBin)).c_str())));
        bundle.truthTemplate.reset(static_cast<TH2D *>(measured2D.Clone(("truthTemplateWeighted2D_" + EtaToken(etaBin)).c_str())));
        bundle.missedTruth.reset(static_cast<TH2D *>(measured2D.Clone(("missedTruthWeighted2D_" + EtaToken(etaBin)).c_str())));
        bundle.fakeReco.reset(static_cast<TH2D *>(measured2D.Clone(("fakeRecoWeighted2D_" + EtaToken(etaBin)).c_str())));
        bundle.recoTemplate->Reset();
        bundle.truthTemplate->Reset();
        bundle.missedTruth->Reset();
        bundle.fakeReco->Reset();

        if (detLevelWeights)
        {
            bundle.responseDetWeights.reset(static_cast<TH2D *>(detLevelWeights->Clone(("appliedDetWeights2D_" + EtaToken(etaBin)).c_str())));
            bundle.responseDetWeights->SetDirectory(nullptr);
        }
        if (truthLevelWeights)
        {
            bundle.responseTruthWeights.reset(static_cast<TH2D *>(truthLevelWeights->Clone(("appliedTruthWeights2D_" + EtaToken(etaBin)).c_str())));
            bundle.responseTruthWeights->SetDirectory(nullptr);
        }

        // Build uncut weighted truth projection and use it for the RooUnfoldResponse
        bundle.truthTemplateUncut.reset(static_cast<TH2D *>(measured2D.Clone(("truthTemplateUncutWeighted2D_" + EtaToken(etaBin)).c_str())));
        bundle.truthTemplateUncut->Reset();
        FillTruthProjectionUncut2D(responseTree_, etaBin, 1.0, bundle.truthTemplateUncut.get(), detLevelWeights, truthLevelWeights);
        if (triggeredResponseTree_)
        {
            FillTruthProjectionUncut2D(triggeredResponseTree_, etaBin, config_.triggeredScale, bundle.truthTemplateUncut.get(), detLevelWeights, truthLevelWeights);
        }

        bundle.response = std::make_unique<RooUnfoldResponse>(bundle.recoTemplate.get(), bundle.truthTemplateUncut.get());

        FillResponseTree2D(responseTree_, etaBin, 1.0, *bundle.response, *bundle.truthTemplate, *bundle.missedTruth, *bundle.fakeReco,
                           detLevelWeights, truthLevelWeights);
        if (triggeredResponseTree_)
        {
            FillResponseTree2D(triggeredResponseTree_, etaBin, config_.triggeredScale, *bundle.response, *bundle.truthTemplate, *bundle.missedTruth, *bundle.fakeReco,
                               detLevelWeights, truthLevelWeights);
        }

        if (TH2 *matrix = bundle.response->Hresponse())
        {
            MapUnderOverflowToEdges(matrix);
        }
        MapUnderOverflowToEdges(bundle.truthTemplate.get());
        MapUnderOverflowToEdges(bundle.missedTruth.get());
        MapUnderOverflowToEdges(bundle.fakeReco.get());

        bundle.measuredForUnfold.reset(static_cast<TH2D *>(measured2D.Clone(("measuredForUnfoldWeighted2D_" + EtaToken(etaBin)).c_str())));
        bundle.measuredForUnfold->SetDirectory(nullptr);
        MapUnderOverflowToEdges(bundle.measuredForUnfold.get());
        ReportResponseDiagnostics2D(etaBin, measured2D, bundle, bundle.measuredForUnfold->Integral(), bundle.response->Hmeasured() ? bundle.response->Hmeasured()->Integral() : 0.0);

        return bundle;
    }

    void ReportResponseDiagnostics(const JetBin &jetBin,
                                   const EtaBin &etaBin,
                                   const TH1D &measured,
                                   const ResponseBundle &bundle,
                                   double dataIntegral,
                                   double responseIntegral) const
    {
        const TH1 *responseMeasured = bundle.response->Hmeasured();
        const TH1 *responseTruth = bundle.response->Htruth();
        const double truthIntegral = responseTruth ? responseTruth->Integral() : 0.0;
        const double missedIntegral = bundle.missedTruth ? bundle.missedTruth->Integral() : 0.0;
        const double fakeIntegral = bundle.fakeReco ? bundle.fakeReco->Integral() : 0.0;
        const double scaleFactor = dataIntegral > 0.0 ? responseIntegral / dataIntegral : 0.0;
        const double efficiency = truthIntegral > 0.0 ? responseIntegral / truthIntegral : 0.0;

        std::cout << "[diag][1D] jet=" << jetBin.label
                  << " eta=" << etaBin.id
                  << " measured_raw=" << measured.Integral()
                  << " measured_input=" << dataIntegral
                  << " response_measured=" << responseIntegral
                  << " response_truth=" << truthIntegral
                  << " missed=" << missedIntegral
                  << " fake=" << fakeIntegral
                  << " scaleFactor=" << scaleFactor
                  << " recoEff=" << efficiency
                  << std::endl;

        if (responseMeasured)
        {
            const Chi2Result beforeScaleChi2 = ComputeChi2(&measured, responseMeasured);
            std::cout << "[diag][1D] jet=" << jetBin.label
                      << " eta=" << etaBin.id
                      << " chi2(measured,responseMeasured)=" << beforeScaleChi2.chi2
                      << " ndf=" << beforeScaleChi2.ndf
                      << " reduced=" << beforeScaleChi2.reduced
                      << std::endl;
        }
    }

    void ReportResponseDiagnostics2D(const EtaBin &etaBin,
                                     const TH2D &measured2D,
                                     const ResponseBundle2D &bundle,
                                     double dataIntegral,
                                     double responseIntegral) const
    {
        const TH1 *responseMeasured = bundle.response->Hmeasured();
        const TH1 *responseTruth = bundle.response->Htruth();
        const double truthIntegral = responseTruth ? responseTruth->Integral() : 0.0;
        const double missedIntegral = bundle.missedTruth ? bundle.missedTruth->Integral() : 0.0;
        const double fakeIntegral = bundle.fakeReco ? bundle.fakeReco->Integral() : 0.0;
        const double scaleFactor = dataIntegral > 0.0 ? responseIntegral / dataIntegral : 0.0;
        const double efficiency = truthIntegral > 0.0 ? responseIntegral / truthIntegral : 0.0;

        std::cout << "[diag][2D] eta=" << etaBin.id
                  << " measured_raw=" << measured2D.Integral()
                  << " measured_input=" << dataIntegral
                  << " response_measured=" << responseIntegral
                  << " response_truth=" << truthIntegral
                  << " missed=" << missedIntegral
                  << " fake=" << fakeIntegral
                  << " scaleFactor=" << scaleFactor
                  << " recoEff=" << efficiency
                  << std::endl;

        if (responseMeasured)
        {
            std::unique_ptr<TH1D> measuredProj(measured2D.ProjectionX(("diag_measured_proj_" + EtaToken(etaBin)).c_str(), 1, measured2D.GetNbinsY()));
            std::unique_ptr<TH2D> responseMeasuredClone(static_cast<TH2D *>(responseMeasured->Clone(("diag_response_measured2D_" + EtaToken(etaBin)).c_str())));
            responseMeasuredClone->SetDirectory(nullptr);
            std::unique_ptr<TH1D> responseMeasuredProj(responseMeasuredClone->ProjectionX(("diag_response_proj_" + EtaToken(etaBin)).c_str(), 1, measured2D.GetNbinsY()));
            outputFile_->cd();
            measuredProj->Write();
            responseMeasuredClone->Write();
            responseMeasuredProj->Write();
            const Chi2Result beforeScaleChi2 = ComputeChi2(measuredProj.get(), responseMeasuredProj.get());
            std::cout << "[diag][2D] eta=" << etaBin.id
                      << " chi2(zT projection measured,responseMeasured)=" << beforeScaleChi2.chi2
                      << " ndf=" << beforeScaleChi2.ndf
                      << " reduced=" << beforeScaleChi2.reduced
                      << std::endl;
        }
    }

    // Fill response and truth using only entries in [firstEntry, lastEntry).
    // Also fills pseudoMeasured (reco level) when non-null, for building closure pseudo-data.
    void FillResponseTreeRange(TTree *tree,
                               const JetBin &jetBin,
                               const EtaBin &etaBin,
                               double scale,
                               Long64_t firstEntry,
                               Long64_t lastEntry,
                               RooUnfoldResponse *response,
                               TH1D *truthTemplate,
                               TH1D *pseudoMeasured,
                               const TH1 *detLevelWeights = nullptr,
                               const TH1 *truthLevelWeights = nullptr) const
    {
        float d0ZDet = 0.0f;
        float d0ZMc = 0.0f;
        float jetPtDet = 0.0f;
        float jetPtMc = 0.0f;
        float d0EtaDet = 0.0f;
        float d0EtaMc = 0.0f;
        float jetNconstDet = 0.0f;
        float jetNconstMc = 0.0f;
        float eventWeight = 1.0f;

        tree->SetBranchAddress("d0_z_det", &d0ZDet);
        tree->SetBranchAddress("d0_z_mc", &d0ZMc);
        tree->SetBranchAddress("jet_pt_det", &jetPtDet);
        tree->SetBranchAddress("jet_pt_mc", &jetPtMc);
        tree->SetBranchAddress("d0_eta_det", &d0EtaDet);
        tree->SetBranchAddress("d0_eta_mc", &d0EtaMc);
        tree->SetBranchAddress("jet_nconst_det", &jetNconstDet);
        tree->SetBranchAddress("jet_nconst_mc", &jetNconstMc);

        if (tree->GetBranch("event_weight"))
            tree->SetBranchAddress("event_weight", &eventWeight);
        else if (tree->GetBranch("weight"))
            tree->SetBranchAddress("weight", &eventWeight);
        else if (tree->GetBranch("evtWeight"))
            tree->SetBranchAddress("evtWeight", &eventWeight);
        else if (tree->GetBranch("eventWeight"))
            tree->SetBranchAddress("eventWeight", &eventWeight);
        else if (tree->GetBranch("totalWeight"))
            tree->SetBranchAddress("totalWeight", &eventWeight);

        const Long64_t nEntries = tree->GetEntries();
        const Long64_t end = std::min(lastEntry, nEntries);
        for (Long64_t entry = firstEntry; entry < end; ++entry)
        {
            tree->GetEntry(entry);

            double jetPtDetCorr = jetPtDet;
            double d0ZDetCorr = d0ZDet;
            ApplyResponseRecoCorrection(config_, jetPtDet, d0ZDet, jetPtDetCorr, d0ZDetCorr);

            if (jetPtDetCorr < config_.jetBins.front().min)
                continue;

            const bool passCuts = InEtaBin(d0EtaDet, etaBin) &&
                                   jetPtDetCorr >= jetBin.min && jetPtDetCorr < jetBin.max &&
                                   jetPtMc >= jetBin.min && jetPtMc < jetBin.max &&
                                   jetNconstDet > config_.minJetConstituents &&
                                   jetNconstMc > config_.minJetConstituents;
            if (!passCuts)
                continue;

            double totalWeight = scale * eventWeight;
            if (detLevelWeights)
            {
                totalWeight *= LookupWeight1D(detLevelWeights, d0ZDetCorr);
            }
            else if (truthLevelWeights)
            {
                totalWeight *= LookupWeight1D(truthLevelWeights, d0ZMc);
            }
            if (truthTemplate)
            {
                truthTemplate->Fill(d0ZMc, totalWeight);
            }
            if (response)
            {
                response->Fill(d0ZDetCorr, d0ZMc, totalWeight);
            }
            if (pseudoMeasured)
            {
                pseudoMeasured->Fill(d0ZDetCorr, totalWeight);
            }
        }
    }

    void FillResponseTreeRange2D(TTree *tree,
                                 const EtaBin &etaBin,
                                 double scale,
                                 Long64_t firstEntry,
                                 Long64_t lastEntry,
                                 RooUnfoldResponse *response,
                                 TH2D *truthTemplate,
                                 TH2D *pseudoMeasured,
                                 const TH2 *detLevelWeights = nullptr,
                                 const TH2 *truthLevelWeights = nullptr) const
    {
        float d0ZDet = 0.0f;
        float d0ZMc = 0.0f;
        float jetPtDet = 0.0f;
        float jetPtMc = 0.0f;
        float d0EtaDet = 0.0f;
        float d0EtaMc = 0.0f;
        float jetNconstDet = 0.0f;
        float jetNconstMc = 0.0f;
        float eventWeight = 1.0f;

        tree->SetBranchAddress("d0_z_det", &d0ZDet);
        tree->SetBranchAddress("d0_z_mc", &d0ZMc);
        tree->SetBranchAddress("jet_pt_det", &jetPtDet);
        tree->SetBranchAddress("jet_pt_mc", &jetPtMc);
        tree->SetBranchAddress("d0_eta_det", &d0EtaDet);
        tree->SetBranchAddress("d0_eta_mc", &d0EtaMc);
        tree->SetBranchAddress("jet_nconst_det", &jetNconstDet);
        tree->SetBranchAddress("jet_nconst_mc", &jetNconstMc);

        if (tree->GetBranch("event_weight"))
            tree->SetBranchAddress("event_weight", &eventWeight);
        else if (tree->GetBranch("weight"))
            tree->SetBranchAddress("weight", &eventWeight);
        else if (tree->GetBranch("evtWeight"))
            tree->SetBranchAddress("evtWeight", &eventWeight);
        else if (tree->GetBranch("eventWeight"))
            tree->SetBranchAddress("eventWeight", &eventWeight);
        else if (tree->GetBranch("totalWeight"))
            tree->SetBranchAddress("totalWeight", &eventWeight);

        const double jetMin = config_.jetBins.empty() ? 0.0 : config_.jetBins.front().min;
        const double jetMax = config_.jetBins.empty() ? 0.0 : config_.jetBins.back().max;
        const Long64_t nEntries = tree->GetEntries();
        const Long64_t end = std::min(lastEntry, nEntries);
        for (Long64_t entry = firstEntry; entry < end; ++entry)
        {
            tree->GetEntry(entry);

            double jetPtDetCorr = jetPtDet;
            double d0ZDetCorr = d0ZDet;
            ApplyResponseRecoCorrection(config_, jetPtDet, d0ZDet, jetPtDetCorr, d0ZDetCorr);

            if (jetPtDetCorr < jetMin)
            {
                continue;
            }

            const bool passCuts = InEtaBin(d0EtaDet, etaBin) &&
                                   jetPtDetCorr >= jetMin && jetPtDetCorr < jetMax &&
                                   jetPtMc >= jetMin && jetPtMc < jetMax &&
                                   jetNconstDet > config_.minJetConstituents &&
                                   jetNconstMc > config_.minJetConstituents;
            if (!passCuts)
            {
                continue;
            }

            double totalWeight = scale * eventWeight;
            if (detLevelWeights)
            {
                totalWeight *= LookupWeight2D(detLevelWeights, jetPtDetCorr, d0ZDetCorr);
            }
            else if (truthLevelWeights)
            {
                totalWeight *= LookupWeight2D(truthLevelWeights, jetPtMc, d0ZMc);
            }
            if (truthTemplate)
            {
                truthTemplate->Fill(d0ZMc, jetPtMc, totalWeight);
            }
            if (response)
            {
                response->Fill(d0ZDetCorr, jetPtDetCorr, d0ZMc, jetPtMc, totalWeight);
            }
            if (pseudoMeasured)
            {
                pseudoMeasured->Fill(d0ZDetCorr, jetPtDetCorr, totalWeight);
            }
        }
    }

    void FillResponseTree(TTree *tree,
                          const JetBin &jetBin,
                          const EtaBin &etaBin,
                          double scale,
                          RooUnfoldResponse &response,
                          TH1D &truthTemplate,
                          TH1D &missedTruth,
                          TH1D &fakeReco,
                          const TH1 *detLevelWeights = nullptr,
                          const TH1 *truthLevelWeights = nullptr) const
    {
        float d0ZDet = 0.0f;
        float d0ZMc = 0.0f;
        float jetPtDet = 0.0f;
        float jetPtMc = 0.0f;
        float d0EtaDet = 0.0f;
        float d0EtaMc = 0.0f;
        float jetNconstDet = 0.0f;
        float jetNconstMc = 0.0f;
        float eventWeight = 1.0f;

        tree->SetBranchAddress("d0_z_det", &d0ZDet);
        tree->SetBranchAddress("d0_z_mc", &d0ZMc);
        tree->SetBranchAddress("jet_pt_det", &jetPtDet);
        tree->SetBranchAddress("jet_pt_mc", &jetPtMc);
        tree->SetBranchAddress("d0_eta_det", &d0EtaDet);
        tree->SetBranchAddress("d0_eta_mc", &d0EtaMc);
        tree->SetBranchAddress("jet_nconst_det", &jetNconstDet);
        tree->SetBranchAddress("jet_nconst_mc", &jetNconstMc);

        if (tree->GetBranch("event_weight"))
        {
            tree->SetBranchAddress("event_weight", &eventWeight);
        }
        else if (tree->GetBranch("weight"))
        {
            tree->SetBranchAddress("weight", &eventWeight);
        }
        else if (tree->GetBranch("evtWeight"))
        {
            tree->SetBranchAddress("evtWeight", &eventWeight);
        }
        else if (tree->GetBranch("eventWeight"))
        {
            tree->SetBranchAddress("eventWeight", &eventWeight);
        }
        else if (tree->GetBranch("totalWeight"))
        {
            tree->SetBranchAddress("totalWeight", &eventWeight);
        }

        const Long64_t nEntries = tree->GetEntries();
        for (Long64_t entry = 0; entry < nEntries; ++entry)
        {
            tree->GetEntry(entry);

            double jetPtDetCorr = jetPtDet;
            double d0ZDetCorr = d0ZDet;
            ApplyResponseRecoCorrection(config_, jetPtDet, d0ZDet, jetPtDetCorr, d0ZDetCorr);

            if (jetPtDetCorr < config_.jetBins.front().min)
            {
                continue;
            }

            const bool passResponseCuts = InEtaBin(d0EtaDet, etaBin) &&
                                          jetPtDetCorr >= jetBin.min && jetPtDetCorr < jetBin.max &&
                                          jetPtMc >= jetBin.min && jetPtMc < jetBin.max &&
                                          jetNconstDet > config_.minJetConstituents &&
                                          jetNconstMc > config_.minJetConstituents;

            if (!passResponseCuts)
            {
                continue;
            }

            double totalWeight = scale * eventWeight;
            if (detLevelWeights)
            {
                totalWeight *= LookupWeight1D(detLevelWeights, d0ZDetCorr);
            }
            else if (truthLevelWeights)
            {
                totalWeight *= LookupWeight1D(truthLevelWeights, d0ZMc);
            }
            truthTemplate.Fill(d0ZMc, totalWeight);
            response.Fill(d0ZDetCorr, d0ZMc, totalWeight);
        }
    }

    void FillResponseTree2D(TTree *tree,
                            const EtaBin &etaBin,
                            double scale,
                            RooUnfoldResponse &response,
                            TH2D &truthTemplate,
                            TH2D &missedTruth,
                            TH2D &fakeReco,
                            const TH2 *detLevelWeights = nullptr,
                            const TH2 *truthLevelWeights = nullptr) const
    {
        float d0ZDet = 0.0f;
        float d0ZMc = 0.0f;
        float jetPtDet = 0.0f;
        float jetPtMc = 0.0f;
        float d0EtaDet = 0.0f;
        float d0EtaMc = 0.0f;
        float jetNconstDet = 0.0f;
        float jetNconstMc = 0.0f;
        float eventWeight = 1.0f;

        tree->SetBranchAddress("d0_z_det", &d0ZDet);
        tree->SetBranchAddress("d0_z_mc", &d0ZMc);
        tree->SetBranchAddress("jet_pt_det", &jetPtDet);
        tree->SetBranchAddress("jet_pt_mc", &jetPtMc);
        tree->SetBranchAddress("d0_eta_det", &d0EtaDet);
        tree->SetBranchAddress("d0_eta_mc", &d0EtaMc);
        tree->SetBranchAddress("jet_nconst_det", &jetNconstDet);
        tree->SetBranchAddress("jet_nconst_mc", &jetNconstMc);

        if (tree->GetBranch("event_weight"))
        {
            tree->SetBranchAddress("event_weight", &eventWeight);
        }
        else if (tree->GetBranch("weight"))
        {
            tree->SetBranchAddress("weight", &eventWeight);
        }
        else if (tree->GetBranch("evtWeight"))
        {
            tree->SetBranchAddress("evtWeight", &eventWeight);
        }
        else if (tree->GetBranch("eventWeight"))
        {
            tree->SetBranchAddress("eventWeight", &eventWeight);
        }
        else if (tree->GetBranch("totalWeight"))
        {
            tree->SetBranchAddress("totalWeight", &eventWeight);
        }

        const double jetMin = config_.jetBins.empty() ? 0.0 : config_.jetBins.front().min;
        const double jetMax = config_.jetBins.empty() ? 0.0 : config_.jetBins.back().max;
        const Long64_t nEntries = tree->GetEntries();
        for (Long64_t entry = 0; entry < nEntries; ++entry)
        {
            tree->GetEntry(entry);

            double jetPtDetCorr = jetPtDet;
            double d0ZDetCorr = d0ZDet;
            ApplyResponseRecoCorrection(config_, jetPtDet, d0ZDet, jetPtDetCorr, d0ZDetCorr);

            if (jetPtDetCorr < jetMin)
            {
                continue;
            }

            const bool passResponseCuts = InEtaBin(d0EtaDet, etaBin) &&
                                          jetPtDetCorr >= jetMin && jetPtDetCorr < jetMax &&
                                          jetPtMc >= jetMin && jetPtMc < jetMax &&
                                          jetNconstDet > config_.minJetConstituents &&
                                          jetNconstMc > config_.minJetConstituents;

            if (!passResponseCuts)
            {
                continue;
            }

            double totalWeight = scale * eventWeight;
            if (detLevelWeights)
            {
                totalWeight *= LookupWeight2D(detLevelWeights, jetPtDetCorr, d0ZDetCorr);
            }
            else if (truthLevelWeights)
            {
                totalWeight *= LookupWeight2D(truthLevelWeights, jetPtMc, d0ZMc);
            }
            truthTemplate.Fill(d0ZMc, jetPtMc, totalWeight);
            response.Fill(d0ZDetCorr, jetPtDetCorr, d0ZMc, jetPtMc, totalWeight);
        }
    }

    void PlotMeasuredSpectrum(const JetBin &jetBin, const EtaBin &etaBin, const TH1D &measured)
    {
        TCanvas measuredCanvas(("c_measured_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str(), "measured", 800, 600);
        gStyle->SetOptStat(0);
        TH1D measuredCopy(measured);
        measuredCopy.SetTitle(("Measured D^{0} z_{T}, jet p_{T} " + jetBin.label + ", " + etaBin.label + ";z_{T};counts").c_str());
        measuredCopy.Draw("E1");
        SaveCanvas(&measuredCanvas, plotDirectory_ / etaBin.label / jetBin.label / ("measured_" + jetBin.label + "_" + EtaToken(etaBin) + ".png"));
    }

    IterationScan1DResult ScanIterations1D(const JetBin &jetBin,
                                           const EtaBin &etaBin,
                                           const ResponseBundle &bundle)
    {
        IterationScan1DResult result;
        for (int iteration = 1; iteration <= config_.nIterations; ++iteration)
        {
            RooUnfoldBayes unfold(bundle.response.get(), bundle.measuredForUnfold.get(), iteration);
            TH1D *unfoldedRaw = dynamic_cast<TH1D *>(unfold.Hreco(RooUnfold::kCovToy));
            if (!unfoldedRaw)
            {
                std::cerr << "[unfold] RooUnfold failed for jet bin " << jetBin.label << " in eta bin "
                          << etaBin.id << " at iteration " << iteration << std::endl;
                continue;
            }

            std::unique_ptr<TH1D> unfolded(static_cast<TH1D *>(unfoldedRaw->Clone(("unfolded_zT_" + jetBin.label + "_" + EtaToken(etaBin) + "_iter" + std::to_string(iteration)).c_str())));
            unfolded->SetDirectory(nullptr);
            unfolded->Write();

            TH1D *refoldedRaw = dynamic_cast<TH1D *>(bundle.response->ApplyToTruth(unfolded.get()));
            std::unique_ptr<TH1D> refolded;
            if (refoldedRaw)
            {
                refolded.reset(static_cast<TH1D *>(refoldedRaw->Clone(("refolded_zT_" + jetBin.label + "_" + EtaToken(etaBin) + "_iter" + std::to_string(iteration)).c_str())));
                refolded->SetDirectory(nullptr);
                refolded->Write();
            }

            const TMatrixD covariance = unfold.Ereco(RooUnfold::kCovToy);
            std::unique_ptr<TH2D> covarianceHist(BuildCovarianceHistogram(covariance, "covariance_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin) + "_iter" + std::to_string(iteration)));
            std::unique_ptr<TH2D> correlationHist(BuildCorrelationHistogram(covariance, "correlation_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin) + "_iter" + std::to_string(iteration)));
            covarianceHist->Write();
            correlationHist->Write();

            const Chi2Result refoldChi2 = ComputeChi2(bundle.measuredForUnfold.get(), refolded.get());
            if (config_.verbose)
            {
                std::cout << "[unfold]   iter=" << iteration
                          << " measured=" << bundle.measuredForUnfold->Integral()
                          << " unfolded=" << unfolded->Integral()
                          << " refold chi2/ndf=" << refoldChi2.reduced << std::endl;
            }

            if (refoldChi2.reduced >= 0.0 && refoldChi2.reduced < result.bestChi2)
            {
                result.bestChi2 = refoldChi2.reduced;
                result.bestIteration = iteration;
                result.bestUnfolded.reset(static_cast<TH1D *>(unfolded->Clone(("best_unfolded_zT_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
                result.bestUnfolded->SetDirectory(nullptr);
            }

            PlotIteration(jetBin, etaBin, iteration, *bundle.measuredForUnfold, *unfolded, *bundle.truthTemplate, refolded.get());
        }
        return result;
    }

    IterationScan2DResult ScanIterations2D(const EtaBin &etaBin,
                                           const ResponseBundle2D &bundle)
    {
        IterationScan2DResult result;
        for (int iteration = 1; iteration <= config_.nIterations; ++iteration)
        {
            RooUnfoldBayes unfold(bundle.response.get(), bundle.measuredForUnfold.get(), iteration);
            TH2D *unfoldedRaw = dynamic_cast<TH2D *>(unfold.Hreco(RooUnfold::kCovToy));
            if (!unfoldedRaw)
            {
                std::cerr << "[unfold] 2D RooUnfold failed at iteration " << iteration << std::endl;
                continue;
            }

            std::unique_ptr<TH2D> unfolded(static_cast<TH2D *>(unfoldedRaw->Clone(("h2_unfolded_zT_vs_jetPt_" + EtaToken(etaBin) + "_iter" + std::to_string(iteration)).c_str())));
            unfolded->SetDirectory(nullptr);
            unfolded->Write();

            TH2D *refoldedRaw = dynamic_cast<TH2D *>(bundle.response->ApplyToTruth(unfolded.get()));
            std::unique_ptr<TH2D> refolded;
            if (refoldedRaw)
            {
                refolded.reset(static_cast<TH2D *>(refoldedRaw->Clone(("h2_refolded_zT_vs_jetPt_" + EtaToken(etaBin) + "_iter" + std::to_string(iteration)).c_str())));
                refolded->SetDirectory(nullptr);
                refolded->Write();
            }

            const TMatrixD covariance = unfold.Ereco(RooUnfold::kCovToy);
            std::unique_ptr<TH2D> covarianceHist(BuildCovarianceHistogram(covariance, "covariance2D_" + EtaToken(etaBin) + "_iter" + std::to_string(iteration)));
            std::unique_ptr<TH2D> correlationHist(BuildCorrelationHistogram(covariance, "correlation2D_" + EtaToken(etaBin) + "_iter" + std::to_string(iteration)));
            covarianceHist->Write();
            correlationHist->Write();

            double chi2 = std::numeric_limits<double>::infinity();
            if (refolded)
            {
                std::unique_ptr<TH1D> measuredProj(bundle.measuredForUnfold->ProjectionX(("measuredProj_iter" + std::to_string(iteration)).c_str(), 1, bundle.measuredForUnfold->GetNbinsY()));
                std::unique_ptr<TH1D> refoldedProj(refolded->ProjectionX(("refoldedProj_iter" + std::to_string(iteration)).c_str(), 1, refolded->GetNbinsY()));
                outputFile_->cd();
                measuredProj->Write();
                refoldedProj->Write();
                Chi2Result chi2Result = ComputeChi2(measuredProj.get(), refoldedProj.get());
                chi2 = chi2Result.reduced;
            }

            if (config_.verbose)
            {
                std::cout << "[unfold]   2D iter=" << iteration
                          << " measured=" << bundle.measuredForUnfold->Integral()
                          << " unfolded=" << unfolded->Integral()
                          << " refold proj chi2/ndf=" << chi2 << std::endl;
            }

            if (chi2 >= 0.0 && chi2 < result.bestChi2)
            {
                result.bestChi2 = chi2;
                result.bestIteration = iteration;
                result.bestUnfolded.reset(static_cast<TH2D *>(unfolded->Clone(("h2_bestUnfolded_zT_vs_jetPt_" + EtaToken(etaBin)).c_str())));
                result.bestUnfolded->SetDirectory(nullptr);
            }

            PlotIteration2D(etaBin, iteration, *bundle.measuredForUnfold, *unfolded, refolded.get());
        }
        return result;
    }

    void WriteInputs(const JetBin &jetBin, const EtaBin &etaBin, const TH1D &measured, const ResponseBundle &bundle)
    {
        outputFile_->cd();
        measured.Write(("measured_zT_" + jetBin.label + "_" + EtaToken(etaBin)).c_str());
        bundle.truthTemplate->Write(("truth_zT_" + jetBin.label + "_" + EtaToken(etaBin)).c_str());
        bundle.missedTruth->Write(("missedTruth_zT_" + jetBin.label + "_" + EtaToken(etaBin)).c_str());
        bundle.fakeReco->Write(("fakeReco_zT_" + jetBin.label + "_" + EtaToken(etaBin)).c_str());
        if (const TH1 *responseMeasured = bundle.response->Hmeasured())
        {
            std::unique_ptr<TH1D> responseMeasuredClone(static_cast<TH1D *>(responseMeasured->Clone(("responseMeasured_zT_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
            responseMeasuredClone->SetDirectory(nullptr);
            responseMeasuredClone->Write();
        }
        if (const TH1 *responseTruth = bundle.response->Htruth())
        {
            std::unique_ptr<TH1D> responseTruthClone(static_cast<TH1D *>(responseTruth->Clone(("responseTruth_zT_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
            responseTruthClone->SetDirectory(nullptr);
            responseTruthClone->Write();
        }
        bundle.measuredForUnfold->Write(("measuredForUnfold_zT_" + jetBin.label + "_" + EtaToken(etaBin)).c_str());
        if (TH2 *responseHist = bundle.response->Hresponse())
        {
            std::unique_ptr<TH2D> responseClone(static_cast<TH2D *>(responseHist->Clone(("responseMatrix_" + jetBin.label + "_" + EtaToken(etaBin)).c_str())));
            responseClone->SetDirectory(nullptr);
            responseClone->Write();

            TCanvas canvas(("c_response_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str(), "response", 800, 650);
            gStyle->SetOptStat(0);
            responseClone->GetXaxis()->SetTitle("reco D^{0}_{z} (det)");
            responseClone->GetYaxis()->SetTitle("truth D^{0}_{z} (mc)");
            responseClone->SetTitle("");
            responseClone->Draw("COLZ");
            TLatex* labelRespMatr = new TLatex();
            labelRespMatr->SetNDC();
            labelRespMatr->SetTextSize(0.035);
            labelRespMatr->DrawLatex(0.12, 0.92, ("Response matrix, jet p_{T}: " + jetBin.label + " GeV/c, #eta range: " + etaBin.label).c_str());
            SaveCanvas(&canvas, plotDirectory_ / etaBin.label / jetBin.label / ("response_" + jetBin.label + "_" + EtaToken(etaBin) + ".png"));
        }
    }

    void PlotIteration(const JetBin &jetBin,
                       const EtaBin &etaBin,
                       int iteration,
                       const TH1D &measured,
                       const TH1D &unfolded,
                       const TH1D &truth,
                       const TH1D *refolded)
    {
        TCanvas canvas(("c_iter_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin) + "_" + std::to_string(iteration)).c_str(), "iteration", 850, 850);
        //add dummy histogram for axis labels and range
        TPad upper("upper", "upper", 0.0, 0.30, 1.0, 1.0);
        TPad lower("lower", "lower", 0.0, 0.0, 1.0, 0.30);
        upper.SetBottomMargin(0.02);
        lower.SetTopMargin(0.03);
        lower.SetBottomMargin(0.30);
        upper.Draw();
        lower.Draw();

        upper.cd();
        upper.cd()->SetTickx(1);
        upper.cd()->SetTicky(1);
        upper.cd()->SetTopMargin(0.03);
        upper.cd()->SetRightMargin(0.01);
        upper.cd()->SetBottomMargin(0.00);
        // Draw labels for jet pT and rapidity bin
        TLatex* latex = new TLatex();
        latex->SetNDC();
        latex->SetTextSize(0.035);
        {
            std::ostringstream jl;
            jl << "jet p_{T}: " << jetBin.min << " - " << jetBin.max << " GeV/c";
            latex->DrawLatex(0.12, 0.92, jl.str().c_str());
        }
        {
            std::string rlabel = etaBin.label.empty() ? ("eta bin " + std::to_string(etaBin.id)) : etaBin.label;
            latex->DrawLatex(0.12, 0.88, rlabel.c_str());
        }
        TH1D measuredCopy(measured);
        TH1D unfoldedCopy(unfolded);
        TH1D truthCopy(truth);
        measuredCopy.SetLineColor(kBlack);
        measuredCopy.SetMarkerColor(kBlack);
        measuredCopy.SetMarkerStyle(20);
        unfoldedCopy.SetLineColor(kBlue + 1);
        unfoldedCopy.GetYaxis()->SetRangeUser(0,unfoldedCopy.GetMaximum()*1.5);
        unfoldedCopy.SetTitle("");
        unfoldedCopy.SetMarkerColor(kBlue + 1);
        unfoldedCopy.SetMarkerStyle(24);
        //scale the truthCopy to match the integral of the unfolded result for shape comparison
        const double unfoldedIntegral = unfoldedCopy.Integral();
        const double truthIntegral = truthCopy.Integral();
        if (truthIntegral > 0.0 && unfoldedIntegral > 0.0)        {
            truthCopy.Scale(unfoldedIntegral / truthIntegral);
        }
        truthCopy.SetLineColor(kGreen + 2);
        truthCopy.SetMarkerColor(kGreen + 2);
        truthCopy.SetMarkerStyle(21);
        measuredCopy.SetTitle(("D^{0} z_{T} unfolding, jet p_{T} " + jetBin.label + ", " + etaBin.label + ";z_{T};counts").c_str());
        unfoldedCopy.Draw("E1");
        measuredCopy.Draw("E1 SAME");
        truthCopy.Draw("E1 SAME");

        std::unique_ptr<TH1D> refoldedCopy;
        if (refolded)
        {
            refoldedCopy.reset(static_cast<TH1D *>(refolded->Clone(("refoldedCopy_" + jetBin.label + "_" + EtaToken(etaBin) + "_" + std::to_string(iteration)).c_str())));
            refoldedCopy->SetDirectory(nullptr);
            refoldedCopy->SetLineColor(kRed + 1);
            refoldedCopy->SetMarkerColor(kRed + 1);
            refoldedCopy->SetMarkerStyle(25);
            refoldedCopy->Draw("E1 SAME");
        }

        TLegend legend(0.13, 0.72, 0.46, 0.92);
        legend.SetBorderSize(0);
        legend.SetFillStyle(0);
        legend.AddEntry(&measuredCopy, "Measured", "lep");
        legend.AddEntry(&unfoldedCopy, ("Unfolded iter " + std::to_string(iteration)).c_str(), "lep");
        legend.AddEntry(&truthCopy, "Truth (rescaled)", "lep");
        if (refoldedCopy)
        {
            legend.AddEntry(refoldedCopy.get(), "Refolded", "lep");
        }
        legend.Draw();

        TLatex* labelJetPtAndEta = new TLatex();
        labelJetPtAndEta->SetNDC();
        labelJetPtAndEta->SetTextSize(0.04);
        labelJetPtAndEta->DrawLatex(0.52, 0.90, ("Jet p_{T}: " + jetBin.label + " GeV/c, #eta range: " + etaBin.label).c_str());

        lower.cd();
        lower.cd()->SetTopMargin(0.0);
        lower.cd()->SetRightMargin(0.01);
        lower.cd()->SetTickx(1);
        lower.cd()->SetTicky(1);

        std::unique_ptr<TH1D> ratio(refolded ? MakeRatioHistogram(refolded, &measured, "ratio_" + jetBin.label + "_" + EtaToken(etaBin) + "_iter" + std::to_string(iteration), "refolded / measured") : nullptr);
        if (ratio)
        {
            outputFile_->cd();
            ratio->Write();
            ratio->SetLineColor(kRed + 1);
            ratio->SetMarkerColor(kRed + 1);
            ratio->SetMarkerStyle(20);
            ratio->SetTitle("");
            ratio->GetYaxis()->SetRangeUser(0.81, 1.19);
            ratio->GetYaxis()->SetTitleSize(0.09);
            ratio->GetYaxis()->SetTitleOffset(0.45);
            ratio->GetYaxis()->SetLabelSize(0.09);
            ratio->GetXaxis()->SetTitleSize(0.11);
            ratio->GetXaxis()->SetLabelSize(0.10);
            ratio->Draw("E1");
            TLine* unityLine = new TLine(ratio->GetXaxis()->GetXmin(), 1.0, ratio->GetXaxis()->GetXmax(), 1.0);
            unityLine->SetLineStyle(2);
            unityLine->SetLineColor(kBlack);
            unityLine->SetLineWidth(2);
            unityLine->Draw("same");
        }

        SaveCanvas(&canvas, plotDirectory_ / etaBin.label / jetBin.label / ("unfold_iter" + std::to_string(iteration) + "_" + jetBin.label + "_" + EtaToken(etaBin) + ".png"));
    }

    void RunClosure(const JetBin &jetBin, const EtaBin &etaBin, const ResponseBundle &bundle, int bestIteration)
    {
        // Split-sample closure: first half of MC builds the response matrix, second half
        // provides the pseudo-measured (reco) input and the reference truth distribution.
        // This avoids the circularity of unfolding the same distribution used to build the response.
        std::cout << "[closure] Running split-sample closure test for jet bin " << jetBin.label
                  << " in eta bin " << etaBin.id << " with iteration " << bestIteration << std::endl;

        // Use the same binning as the nominal unfolded histogram (from bundle templates).
        if (!bundle.recoTemplate || !bundle.truthTemplate)
        {
            std::cout << "[closure] No template histograms available; cannot run closure." << std::endl;
            return;
        }

        // ---- build first-half response ----
        TH1D* closureRecoTpl = static_cast<TH1D *>(bundle.recoTemplate->Clone(("closureRecoTpl_" + jetBin.label + "_" + EtaToken(etaBin)).c_str()));
        TH1D* closureTruthTpl = static_cast<TH1D *>(bundle.truthTemplate->Clone(("closureTruthTpl_" + jetBin.label + "_" + EtaToken(etaBin)).c_str()));
        closureRecoTpl->SetDirectory(nullptr);
        closureTruthTpl->SetDirectory(nullptr);
        closureRecoTpl->Reset();
        closureTruthTpl->Reset();

        RooUnfoldResponse closureResponse(closureRecoTpl, closureTruthTpl);

        // ---- build second-half pseudo-data (reco) and reference truth ----
        TH1D* pseudoMeasured = static_cast<TH1D *>(bundle.recoTemplate->Clone(("closurePseudoData_" + jetBin.label + "_" + EtaToken(etaBin)).c_str()));
        TH1D* pseudoTruth = static_cast<TH1D *>(bundle.truthTemplate->Clone(("closurePseudoTruth_" + jetBin.label + "_" + EtaToken(etaBin)).c_str()));
        pseudoMeasured->SetDirectory(nullptr);
        pseudoTruth->SetDirectory(nullptr);
        pseudoMeasured->Reset();
        pseudoTruth->Reset();

        const TH1 *detLevelWeights = bundle.responseDetWeights.get();
        const TH1 *truthLevelWeights = bundle.responseTruthWeights.get();

        // fill from main response tree
        {
            const Long64_t n = responseTree_->GetEntries();
            const Long64_t half = n / 2;
            FillResponseTreeRange(responseTree_, jetBin, etaBin, 1.0, 0, half,
                                  &closureResponse, closureTruthTpl, nullptr,
                                  detLevelWeights, truthLevelWeights);
            FillResponseTreeRange(responseTree_, jetBin, etaBin, 1.0, half, n,
                                  nullptr, pseudoTruth, pseudoMeasured,
                                  detLevelWeights, truthLevelWeights);
        }
        // fill from triggered response tree (same split strategy, scale applied)
        if (triggeredResponseTree_)
        {
            const Long64_t n = triggeredResponseTree_->GetEntries();
            const Long64_t half = n / 2;
            FillResponseTreeRange(triggeredResponseTree_, jetBin, etaBin, config_.triggeredScale, 0, half,
                                  &closureResponse, closureTruthTpl, nullptr,
                                  detLevelWeights, truthLevelWeights);
            FillResponseTreeRange(triggeredResponseTree_, jetBin, etaBin, config_.triggeredScale, half, n,
                                  nullptr, pseudoTruth, pseudoMeasured,
                                  detLevelWeights, truthLevelWeights);
        }

        if (TH2 *matrix = closureResponse.Hresponse())
            MapUnderOverflowToEdges(matrix);
        MapUnderOverflowToEdges(closureTruthTpl);
        MapUnderOverflowToEdges(pseudoMeasured);
        MapUnderOverflowToEdges(pseudoTruth);

        const TH1 *closureMeasured = closureResponse.Hmeasured();
        const TH1 *closureTruth = closureResponse.Htruth();
        std::cout << "[closure] jet=" << jetBin.label
                  << " eta=" << etaBin.id
                  << " response_measured=" << (closureMeasured ? closureMeasured->Integral() : 0.0)
                  << " response_truth=" << (closureTruth ? closureTruth->Integral() : 0.0)
                  << " pseudo_measured=" << pseudoMeasured->Integral()
                  << " pseudo_truth=" << pseudoTruth->Integral()
                  << std::endl;
        for (int bin = 1; bin <= std::min(5, pseudoMeasured->GetNbinsX()); ++bin)
        {
            std::cout << "[closure] bin " << bin
                      << " pseudo_measured=" << pseudoMeasured->GetBinContent(bin)
                      << " pseudo_truth=" << pseudoTruth->GetBinContent(bin)
                      << std::endl;
        }

        if (pseudoMeasured->Integral() == 0.0)
        {
            std::cout << "[closure] Second-half pseudo-data is empty for jet bin " << jetBin.label
                      << " in eta bin " << etaBin.id << "; skipping closure." << std::endl;
            return;
        }
        if (pseudoTruth->Integral() == 0.0)
        {
            std::cout << "[closure] Second-half pseudo-truth is empty for jet bin " << jetBin.label
                      << " in eta bin " << etaBin.id << "; skipping closure." << std::endl;
            return;
        }

        outputFile_->cd();
        pseudoMeasured->Write();
        pseudoTruth->Write();

        RooUnfoldBayes closure(&closureResponse, pseudoMeasured, bestIteration);
        TH1D *closureRaw = dynamic_cast<TH1D *>(closure.Hreco(RooUnfold::kCovToy));
        if (!closureRaw)
        {
            std::cout << "[closure] RooUnfold failed for closure test in jet bin " << jetBin.label << " in eta bin " << etaBin.id << std::endl;
            return;
        }

        TH1D* closureUnfolded = static_cast<TH1D *>(closureRaw->Clone(("closureUnfolded_" + jetBin.label + "_" + EtaToken(etaBin)).c_str()));
        closureUnfolded->SetDirectory(nullptr);
        closureUnfolded->Write();

        std::cout << "[closure] unfolded_integral=" << closureUnfolded->Integral() << std::endl;

        // ratio is closure-unfolded / second-half truth
        std::unique_ptr<TH1D> ratio(MakeRatioHistogram(closureUnfolded, pseudoTruth, "closureRatio_" + jetBin.label + "_" + EtaToken(etaBin), "closure / truth"));
        if (!ratio)
        {
            std::cout << "[closure] Failed to build closure ratio for jet bin " << jetBin.label
                      << " in eta bin " << etaBin.id << std::endl;
            return;
        }
        ratio->Write();

        TCanvas canvas(("c_closure_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str(), "closure", 850, 850);
        TPad upper("upperClosure", "upperClosure", 0.0, 0.30, 1.0, 1.0);
        TPad lower("lowerClosure", "lowerClosure", 0.0, 0.0, 1.0, 0.30);
        upper.SetBottomMargin(0.02);
        lower.SetTopMargin(0.03);
        lower.SetBottomMargin(0.30);
        upper.Draw();
        lower.Draw();

        upper.cd();
        TH1D truthCopy(*pseudoTruth);
        TH1D closureCopy(*closureUnfolded);
        truthCopy.SetLineColor(kBlack);
        truthCopy.SetMarkerColor(kBlack);
        truthCopy.SetMarkerStyle(20);
        closureCopy.SetLineColor(kBlue + 1);
        closureCopy.SetMarkerColor(kBlue + 1);
        closureCopy.SetMarkerStyle(24);
        truthCopy.SetTitle(("Closure test (split-sample), jet p_{T} " + jetBin.label + ", " + etaBin.label + ";z_{T};counts").c_str());
        truthCopy.Draw("E1");
        closureCopy.Draw("E1 SAME");

        TLegend legend(0.58, 0.76, 0.88, 0.88);
        legend.SetBorderSize(0);
        legend.SetFillStyle(0);
        legend.AddEntry(&truthCopy, "Truth", "lep");
        legend.AddEntry(&closureCopy, ("Closure unfolded iter " + std::to_string(bestIteration)).c_str(), "lep");
        legend.Draw();

        lower.cd();
        ratio->GetYaxis()->SetRangeUser(0.7, 1.3);
        ratio->GetYaxis()->SetTitleSize(0.09);
        ratio->GetYaxis()->SetTitleOffset(0.45);
        ratio->GetYaxis()->SetLabelSize(0.09);
        ratio->GetXaxis()->SetTitleSize(0.11);
        ratio->GetXaxis()->SetLabelSize(0.10);
        ratio->Draw("E1");
        TLine unity(ratio->GetXaxis()->GetXmin(), 1.0, ratio->GetXaxis()->GetXmax(), 1.0);
        unity.SetLineStyle(2);
        unity.Draw();

        SaveCanvas(&canvas, plotDirectory_ / etaBin.label / jetBin.label / ("closure_" + jetBin.label + "_" + EtaToken(etaBin) + ".png"));
    }

    void PlotSummary(const EtaBin &etaBin, TH2D &measuredSummary, TH2D &unfoldedSummary)
    {
        TCanvas measuredCanvas(("c_summary_measured_" + EtaToken(etaBin)).c_str(), "summary measured", 900, 650);
        measuredSummary.Draw("COLZ");
        SaveCanvas(&measuredCanvas, plotDirectory_ / etaBin.label / ("summary_measured_" + EtaToken(etaBin) + ".png"));

        TCanvas unfoldedCanvas(("c_summary_unfolded_" + EtaToken(etaBin)).c_str(), "summary unfolded", 900, 650);
        unfoldedSummary.Draw("COLZ");
        SaveCanvas(&unfoldedCanvas, plotDirectory_ / etaBin.label / ("summary_unfolded_" + EtaToken(etaBin) + ".png"));
    }

    void PlotIteration2D(const EtaBin &etaBin, int iteration, const TH2D &measured, const TH2D &unfolded, const TH2D *refolded)
    {
        TCanvas canvas(("c_2d_iter_" + EtaToken(etaBin) + "_" + std::to_string(iteration)).c_str(), "2D iteration", 1500, 500);
        canvas.Divide(3, 1);

        canvas.cd(1);
        TH2D measuredCopy(measured);
        measuredCopy.SetTitle("Measured 2D;z_{T};p_{T}^{jet} [GeV/c]");
        measuredCopy.Draw("COLZ");

        canvas.cd(2);
        TH2D unfoldedCopy(unfolded);
        unfoldedCopy.SetTitle(("Unfolded 2D iter " + std::to_string(iteration) + ";z_{T};p_{T}^{jet} [GeV/c]").c_str());
        unfoldedCopy.Draw("COLZ");

        canvas.cd(3);
        if (refolded)
        {
            TH2D refoldedCopy(*refolded);
            refoldedCopy.SetTitle(("Refolded 2D iter " + std::to_string(iteration) + ";z_{T};p_{T}^{jet} [GeV/c]").c_str());
            refoldedCopy.Draw("COLZ");
        }

        SaveCanvas(&canvas, plotDirectory_ / etaBin.label / "2D" / ("unfold2D_iter" + std::to_string(iteration) + "_" + EtaToken(etaBin) + ".png"));

        if (refolded)
        {
            std::unique_ptr<TH2D> ratio(MakeRatioHistogram2D(refolded, &measured, "refoldRatio2D_" + EtaToken(etaBin) + "_iter" + std::to_string(iteration), "refolded / measured"));
            outputFile_->cd();
            ratio->Write();
            TCanvas ratioCanvas(("c_2d_ratio_" + EtaToken(etaBin) + "_" + std::to_string(iteration)).c_str(), "2D ratio", 700, 600);
            ratio->Draw("COLZ");
            SaveCanvas(&ratioCanvas, plotDirectory_ / etaBin.label / "2D" / ("refold_ratio2D_iter" + std::to_string(iteration) + "_" + EtaToken(etaBin) + ".png"));
        }
    }

    void RunClosure2D(const EtaBin &etaBin, const ResponseBundle2D &bundle, int bestIteration)
    {
        if (!bundle.recoTemplate || !bundle.truthTemplate)
        {
            return;
        }

        std::cout << "[closure] Running split-sample 2D closure test for eta bin " << etaBin.id
                  << " with iteration " << bestIteration << std::endl;

        TH2D *closureRecoTpl = static_cast<TH2D *>(bundle.recoTemplate->Clone(("closureRecoTpl2D_" + EtaToken(etaBin)).c_str()));
        TH2D *closureTruthTpl = static_cast<TH2D *>(bundle.truthTemplate->Clone(("closureTruthTpl2D_" + EtaToken(etaBin)).c_str()));
        closureRecoTpl->SetDirectory(nullptr);
        closureTruthTpl->SetDirectory(nullptr);
        closureRecoTpl->Reset();
        closureTruthTpl->Reset();

        RooUnfoldResponse closureResponse(closureRecoTpl, closureTruthTpl);

        TH2D *pseudoMeasured = static_cast<TH2D *>(bundle.recoTemplate->Clone(("closurePseudoData2D_" + EtaToken(etaBin)).c_str()));
        TH2D *pseudoTruth = static_cast<TH2D *>(bundle.truthTemplate->Clone(("closurePseudoTruth2D_" + EtaToken(etaBin)).c_str()));
        pseudoMeasured->SetDirectory(nullptr);
        pseudoTruth->SetDirectory(nullptr);
        pseudoMeasured->Reset();
        pseudoTruth->Reset();

        const TH2 *detLevelWeights = bundle.responseDetWeights.get();
        const TH2 *truthLevelWeights = bundle.responseTruthWeights.get();

        {
            const Long64_t n = responseTree_->GetEntries();
            const Long64_t half = n / 2;
            FillResponseTreeRange2D(responseTree_, etaBin, 1.0, 0, half,
                                    &closureResponse, closureTruthTpl, nullptr,
                                    detLevelWeights, truthLevelWeights);
            FillResponseTreeRange2D(responseTree_, etaBin, 1.0, half, n,
                                    nullptr, pseudoTruth, pseudoMeasured,
                                    detLevelWeights, truthLevelWeights);
        }
        if (triggeredResponseTree_)
        {
            const Long64_t n = triggeredResponseTree_->GetEntries();
            const Long64_t half = n / 2;
            FillResponseTreeRange2D(triggeredResponseTree_, etaBin, config_.triggeredScale, 0, half,
                                    &closureResponse, closureTruthTpl, nullptr,
                                    detLevelWeights, truthLevelWeights);
            FillResponseTreeRange2D(triggeredResponseTree_, etaBin, config_.triggeredScale, half, n,
                                    nullptr, pseudoTruth, pseudoMeasured,
                                    detLevelWeights, truthLevelWeights);
        }

        if (TH2 *matrix = closureResponse.Hresponse())
        {
            MapUnderOverflowToEdges(matrix);
        }
        MapUnderOverflowToEdges(closureTruthTpl);
        MapUnderOverflowToEdges(pseudoMeasured);
        MapUnderOverflowToEdges(pseudoTruth);

        if (pseudoMeasured->Integral() == 0.0 || pseudoTruth->Integral() == 0.0)
        {
            std::cout << "[closure] Empty 2D pseudo-data or pseudo-truth for eta bin " << etaBin.id
                      << "; skipping closure." << std::endl;
            return;
        }

        outputFile_->cd();
        pseudoMeasured->Write();
        pseudoTruth->Write();

        RooUnfoldBayes closure(&closureResponse, pseudoMeasured, bestIteration);
        TH2D *closureRaw = dynamic_cast<TH2D *>(closure.Hreco(RooUnfold::kCovToy));
        if (!closureRaw)
        {
            return;
        }

        std::unique_ptr<TH2D> closureUnfolded(static_cast<TH2D *>(closureRaw->Clone(("closureUnfolded2D_" + EtaToken(etaBin)).c_str())));
        closureUnfolded->SetDirectory(nullptr);
        closureUnfolded->Write();

        std::unique_ptr<TH2D> ratio(MakeRatioHistogram2D(closureUnfolded.get(), pseudoTruth, "closureRatio2D_" + EtaToken(etaBin), "closure / truth"));
        ratio->Write();

        TCanvas canvas(("c_closure2D_" + EtaToken(etaBin)).c_str(), "closure2D", 1400, 500);
        canvas.Divide(3, 1);
        canvas.cd(1);
        TH2D truthCopy(*pseudoTruth);
        truthCopy.SetTitle("Truth 2D;z_{T};p_{T}^{jet} [GeV/c]");
        truthCopy.Draw("COLZ");
        canvas.cd(2);
        TH2D closureCopy(*closureUnfolded);
        closureCopy.SetTitle("Closure unfolded 2D;z_{T};p_{T}^{jet} [GeV/c]");
        closureCopy.Draw("COLZ");
        canvas.cd(3);
        ratio->Draw("COLZ");
        SaveCanvas(&canvas, plotDirectory_ / etaBin.label / "2D" / ("closure2D_" + EtaToken(etaBin) + ".png"));
    }

    void Run2D(const EtaBin &etaBin, TH2D &measuredSummary)
    {
        if (config_.verbose)
        {
            std::cout << "[unfold] Running 2D unfolding in zT and jet pT for eta bin " << etaBin.id << std::endl;
        }

        ResponseBundle2D bundle = BuildResponse2D(etaBin, measuredSummary);
        outputFile_->cd();
        measuredSummary.Write(("h2_measured_zT_vs_jetPt_input_" + EtaToken(etaBin)).c_str());
        bundle.truthTemplate->Write(("h2_truth_zT_vs_jetPt_" + EtaToken(etaBin)).c_str());
        bundle.missedTruth->Write(("h2_missedTruth_zT_vs_jetPt_" + EtaToken(etaBin)).c_str());
        bundle.fakeReco->Write(("h2_fakeReco_zT_vs_jetPt_" + EtaToken(etaBin)).c_str());
        bundle.measuredForUnfold->Write(("h2_measuredForUnfold_zT_vs_jetPt_" + EtaToken(etaBin)).c_str());
        if (const TH1 *responseMeasured = bundle.response->Hmeasured())
        {
            std::unique_ptr<TH2D> responseMeasuredClone(static_cast<TH2D *>(responseMeasured->Clone(("h2_responseMeasured_zT_vs_jetPt_" + EtaToken(etaBin)).c_str())));
            responseMeasuredClone->SetDirectory(nullptr);
            responseMeasuredClone->Write();
        }
        if (const TH1 *responseTruth = bundle.response->Htruth())
        {
            std::unique_ptr<TH2D> responseTruthClone(static_cast<TH2D *>(responseTruth->Clone(("h2_responseTruth_zT_vs_jetPt_" + EtaToken(etaBin)).c_str())));
            responseTruthClone->SetDirectory(nullptr);
            responseTruthClone->Write();
        }
        if (TH2 *responseHist = bundle.response->Hresponse())
        {
            std::unique_ptr<TH2> responseClone(static_cast<TH2 *>(responseHist->Clone(("responseMatrix2D_flat_" + EtaToken(etaBin)).c_str())));
            responseClone->SetDirectory(nullptr);
            responseClone->Write();
        }

        // Build and save kinematic efficiency 1D histograms per jet-pt bin (from 2D response)
        if (bundle.truthTemplate && bundle.truthTemplateUncut)
        {
            std::vector<std::unique_ptr<TH1D>> kinEffHists2D;
            // derive number of jet-pt bins from measuredSummary (reference)
            const int ny = measuredSummary.GetNbinsY();
            for (int y = 1; y <= ny; ++y)
            {
                std::string projCutName = "kinEff_proj_cut_" + EtaToken(etaBin) + "_y" + std::to_string(y);
                std::string projUncutName = "kinEff_proj_uncut_" + EtaToken(etaBin) + "_y" + std::to_string(y);
                TH1D *projCut = bundle.truthTemplate->ProjectionX(projCutName.c_str(), y, y);
                TH1D *projUncut = bundle.truthTemplateUncut->ProjectionX(projUncutName.c_str(), y, y);
                if (projCut && projUncut)
                {
                    projCut->SetDirectory(nullptr);
                    projUncut->SetDirectory(nullptr);
                    std::string effName = (std::string("kinEff_zT_y") + std::to_string(y) + "_" + EtaToken(etaBin));
                    TH1D *eff = MakeRatioHistogram(projCut, projUncut, effName.c_str(), "efficiency");
                    if (eff)
                    {
                        eff->SetDirectory(nullptr);
                        // style
                        const int colors[] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 2, kOrange + 7, kCyan + 2};
                        const int ncolors = sizeof(colors) / sizeof(colors[0]);
                        eff->SetLineColor(colors[(y - 1) % ncolors]);
                        eff->SetMarkerColor(colors[(y - 1) % ncolors]);
                        eff->SetLineWidth(2);
                        outputFile_->cd();
                        eff->Write(effName.c_str());
                        kinEffHists2D.emplace_back(std::unique_ptr<TH1D>(eff));
                    }
                }
                delete projCut;
                delete projUncut;
            }
            // Overlay the per-jet efficiency histograms
            if (!kinEffHists2D.empty())
            {
                TCanvas ckin2(("c_kinEff2D_" + EtaToken(etaBin)).c_str(), "kinematic efficiency per-jetpt", 900, 650);
                kinEffHists2D[0]->Draw("E1");
                TLegend leg2(0.65, 0.7, 0.9, 0.9);
                leg2.SetBorderSize(0);
                leg2.AddEntry(kinEffHists2D[0].get(), kinEffHists2D[0]->GetName(), "lep");
                for (std::size_t i = 1; i < kinEffHists2D.size(); ++i)
                {
                    kinEffHists2D[i]->Draw("E1 SAME");
                    leg2.AddEntry(kinEffHists2D[i].get(), kinEffHists2D[i]->GetName(), "lep");
                }
                leg2.Draw();
                SaveCanvas(&ckin2, plotDirectory_ / etaBin.label / "2D" / ("kinematic_efficiency_overlay_" + EtaToken(etaBin) + ".png"));
            }
        }

        TCanvas measuredCanvas(("c_measured2D_input_" + EtaToken(etaBin)).c_str(), "measured2D", 900, 650);
        TH2D measuredCopy(measuredSummary);
        measuredCopy.Draw("COLZ");
        SaveCanvas(&measuredCanvas, plotDirectory_ / etaBin.label / "2D" / ("measured2D_" + EtaToken(etaBin) + ".png"));

        TDirectory *dir2D = outputFile_->GetDirectory(("two_dimensional_unfolding_" + EtaToken(etaBin)).c_str());
        if (!dir2D)
        {
            dir2D = outputFile_->mkdir(("two_dimensional_unfolding_" + EtaToken(etaBin)).c_str());
        }
        dir2D->cd();
        IterationScan2DResult scanResult;
        if (config_.doResponseDetLevelWeighting)
        {
            const TH1 *responseMeasuredNominal = bundle.response->Hmeasured();
            std::unique_ptr<TH2D> responseMeasured2D;
            if (responseMeasuredNominal)
            {
                responseMeasured2D.reset(static_cast<TH2D *>(responseMeasuredNominal->Clone(("responseMeasured2D_nominal_" + EtaToken(etaBin)).c_str())));
                responseMeasured2D->SetDirectory(nullptr);
                MapUnderOverflowToEdges(responseMeasured2D.get());
            }

            std::unique_ptr<TH2D> detWeights = responseMeasured2D ? BuildRelativeWeightHistogram2D(*bundle.measuredForUnfold, *responseMeasured2D,
                                                                                                     "responseDetWeights2D_" + EtaToken(etaBin))
                                                                   : nullptr;
            if (detWeights)
            {
                outputFile_->cd();
                detWeights->Write();
            }

            ResponseBundle2D weightedBundle = BuildWeightedResponse2D(etaBin, measuredSummary, detWeights.get(), nullptr);
            scanResult = ScanIterations2D(etaBin, weightedBundle);
            bundle = std::move(weightedBundle);
        }
        else
        {
            scanResult = ScanIterations2D(etaBin, bundle);
        }

        if (config_.doResponseTruthLevelWeighting && scanResult.bestUnfolded)
        {
            std::unique_ptr<TH2D> truthWeights = BuildRelativeWeightHistogram2D(*scanResult.bestUnfolded, *bundle.truthTemplate,
                                                                                "responseTruthWeights2D_" + EtaToken(etaBin));
            if (truthWeights)
            {
                outputFile_->cd();
                truthWeights->Write();
            }
            ResponseBundle2D truthWeightedBundle = BuildWeightedResponse2D(etaBin, measuredSummary, nullptr, truthWeights.get());
            scanResult = ScanIterations2D(etaBin, truthWeightedBundle);
            bundle = std::move(truthWeightedBundle);
        }

        std::cout << "[unfold] Best 2D iteration in eta bin " << etaBin.id
                  << " is " << scanResult.bestIteration
                  << " with reduced chi2 = " << scanResult.bestChi2 << std::endl;

        outputFile_->cd();
        if (scanResult.bestUnfolded)
        {
            std::unique_ptr<TH2D> bestUnfolded(static_cast<TH2D *>(scanResult.bestUnfolded->Clone(("h2_bestUnfolded_zT_vs_jetPt_" + EtaToken(etaBin)).c_str())));
            bestUnfolded->SetDirectory(nullptr);
            bestUnfolded->Write();

            TCanvas bestCanvas(("c_best2D_" + EtaToken(etaBin)).c_str(), "best2D", 900, 650);
            bestUnfolded->Draw("COLZ");
            SaveCanvas(&bestCanvas, plotDirectory_ / etaBin.label / "2D" / ("best_unfolded2D_" + EtaToken(etaBin) + ".png"));

            for (std::size_t jetIndex = 0; jetIndex < config_.jetBins.size(); ++jetIndex)
            {
                const JetBin &jetBin = config_.jetBins[jetIndex];
                const int ybin = static_cast<int>(jetIndex) + 1;
                std::unique_ptr<TH1D> measuredProj(measuredSummary.ProjectionX(("measured2D_proj_" + jetBin.label + "_" + EtaToken(etaBin)).c_str(), ybin, ybin));
                std::unique_ptr<TH1D> unfoldedProj(bestUnfolded->ProjectionX(("unfolded2D_proj_" + jetBin.label + "_" + EtaToken(etaBin)).c_str(), ybin, ybin));
                outputFile_->cd();
                measuredProj->Write();
                unfoldedProj->Write();
                TCanvas projCanvas(("c_best2D_proj_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str(), "best2Dproj", 850, 850);
                TPad upper("upper2DProj", "upper2DProj", 0.0, 0.30, 1.0, 1.0);
                TPad lower("lower2DProj", "lower2DProj", 0.0, 0.0, 1.0, 0.30);
                upper.SetBottomMargin(0.02);
                lower.SetTopMargin(0.03);
                lower.SetBottomMargin(0.30);
                upper.Draw();
                lower.Draw();
                upper.cd();
                measuredProj->SetLineColor(kBlack);
                measuredProj->SetMarkerColor(kBlack);
                measuredProj->SetMarkerStyle(20);
                unfoldedProj->SetLineColor(kBlue + 1);
                unfoldedProj->SetMarkerColor(kBlue + 1);
                unfoldedProj->SetMarkerStyle(24);
                measuredProj->SetTitle(("2D projection, jet p_{T} " + jetBin.label + ";z_{T};counts").c_str());
                measuredProj->Draw("E1");
                unfoldedProj->Draw("E1 SAME");
                TLegend legend(0.58, 0.76, 0.88, 0.88);
                legend.SetBorderSize(0);
                legend.SetFillStyle(0);
                legend.AddEntry(measuredProj.get(), "Measured", "lep");
                legend.AddEntry(unfoldedProj.get(), "Best 2D unfolded", "lep");
                legend.Draw();
                lower.cd();
                std::unique_ptr<TH1D> ratio(MakeRatioHistogram(unfoldedProj.get(), measuredProj.get(), "ratio2DProj_" + jetBin.label + "_" + EtaToken(etaBin), "unfolded / measured"));
                outputFile_->cd();
                ratio->Write();
                ratio->GetYaxis()->SetRangeUser(0.4, 1.6);
                ratio->Draw("E1");
                TLine unity(ratio->GetXaxis()->GetXmin(), 1.0, ratio->GetXaxis()->GetXmax(), 1.0);
                unity.SetLineStyle(2);
                unity.Draw();
                SaveCanvas(&projCanvas, plotDirectory_ / etaBin.label / "2D" / ("projection_" + jetBin.label + "_" + EtaToken(etaBin) + ".png"));
            }
        }

        if (config_.doClosure)
        {
            RunClosure2D(etaBin, bundle, scanResult.bestIteration);
        }
    }

    void ComputeResponseBootstrapCovariance(const JetBin &jetBin,
                                           const EtaBin &etaBin,
                                           const ResponseBundle &bundle,
                                           int bestIteration,
                                           const TH1D &bestUnfolded)
    {
        if (!bundle.response)
        {
            return;
        }
        if (config_.nResponseToys <= 1)
        {
            return;
        }

        const TH2 *responseHist = bundle.response->Hresponse();
        if (!responseHist)
        {
            return;
        }

        const int nx = responseHist->GetNbinsX();
        const int ny = responseHist->GetNbinsY();

        std::mt19937 rng(static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::vector<std::unique_ptr<TH1D>> toyUnfolds;
        toyUnfolds.reserve(config_.nResponseToys);

        for (int t = 0; t < config_.nResponseToys; ++t)
        {
            // create toy templates
            std::unique_ptr<TH1D> recoTpl(static_cast<TH1D *>(bundle.recoTemplate->Clone((std::string("recoToy_") + std::to_string(t)).c_str())));
            std::unique_ptr<TH1D> truthTpl(static_cast<TH1D *>(bundle.truthTemplate->Clone((std::string("truthToy_") + std::to_string(t)).c_str())));
            recoTpl->Reset();
            truthTpl->Reset();

            RooUnfoldResponse toyResp(recoTpl.get(), truthTpl.get());

            // For each response matrix bin, Poisson-fluctuate the original bin content
            for (int ix = 1; ix <= nx; ++ix)
            {
                const double recoCenter = responseHist->GetXaxis()->GetBinCenter(ix);
                for (int iy = 1; iy <= ny; ++iy)
                {
                    const double truthCenter = responseHist->GetYaxis()->GetBinCenter(iy);
                    const double orig = responseHist->GetBinContent(ix, iy);
                    if (orig <= 0.0)
                        continue;
                    std::poisson_distribution<int> pdist(static_cast<double>(orig));
                    const int n = pdist(rng);
                    for (int k = 0; k < n; ++k)
                    {
                        toyResp.Fill(recoCenter, truthCenter, 1.0);
                    }
                }
            }

            // Unfold measured with toy response
            RooUnfoldBayes unfoldToy(&toyResp, bundle.measuredForUnfold.get(), bestIteration);
            TH1D *toyUnfoldRaw = dynamic_cast<TH1D *>(unfoldToy.Hreco(RooUnfold::kCovToy));
            if (!toyUnfoldRaw)
            {
                continue;
            }
            std::unique_ptr<TH1D> toyUnfold(static_cast<TH1D *>(toyUnfoldRaw->Clone((std::string("toyUnfold_") + std::to_string(t)).c_str())));
            toyUnfold->SetDirectory(nullptr);
            toyUnfolds.push_back(std::move(toyUnfold));
        }

        if (toyUnfolds.empty())
        {
            return;
        }

        // Build covariance from toys
        TMatrixD cov_response = BuildCovarianceFromToys(toyUnfolds);

        // Build RooUnfold nominal covariance for bestIteration
        RooUnfoldBayes unfoldNominal(bundle.response.get(), bundle.measuredForUnfold.get(), bestIteration);
        const TMatrixD cov_meas = unfoldNominal.Ereco(RooUnfold::kCovToy);

        // Write response-only covariance
        std::unique_ptr<TH2D> covRespHist(BuildCovarianceHistogram(cov_response, "covariance_response_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)));
        outputFile_->cd();
        covRespHist->Write();

        // Add matrices (ensure same size)
        TMatrixD cov_total = cov_meas;
        if (cov_total.GetNrows() == cov_response.GetNrows())
        {
            cov_total += cov_response;
        }

        std::unique_ptr<TH2D> covTotalHist(BuildCovarianceHistogram(cov_total, "covariance_total_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)));
        outputFile_->cd();
        covTotalHist->Write();

        // Produce a clone of the best-unfolded histogram with total uncertainties
        const int nbins = bestUnfolded.GetNbinsX();
        std::unique_ptr<TH1D> bestWithTotal(static_cast<TH1D *>(bestUnfolded.Clone(("best_unfolded_with_totalErr_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str())));
        bestWithTotal->SetDirectory(nullptr);
        if (cov_total.GetNrows() == nbins)
        {
            for (int b = 1; b <= nbins; ++b)
            {
                const double err = cov_total(b - 1, b - 1) > 0.0 ? std::sqrt(cov_total(b - 1, b - 1)) : 0.0;
                bestWithTotal->SetBinError(b, err);
            }
        }
        else
        {
            // fallback: use RooUnfold nominal errors if matrix sizes mismatch
            RooUnfoldBayes unfoldNominal2(bundle.response.get(), bundle.measuredForUnfold.get(), bestIteration);
            const TMatrixD cov_meas2 = unfoldNominal2.Ereco(RooUnfold::kCovToy);
            if (cov_meas2.GetNrows() == nbins)
            {
                for (int b = 1; b <= nbins; ++b)
                {
                    const double err = cov_meas2(b - 1, b - 1) > 0.0 ? std::sqrt(cov_meas2(b - 1, b - 1)) : 0.0;
                    bestWithTotal->SetBinError(b, err);
                }
            }
        }
        outputFile_->cd();
        bestWithTotal->Write();

        // Comparison plot: nominal best_unfolded vs best_with_totalErr
        TCanvas cmpCanvas(("c_best_compare_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str(), "compare", 900, 900);
        cmpCanvas.Divide(1, 2);
        cmpCanvas.cd(1);
        gStyle->SetOptStat(0);
        std::unique_ptr<TH1D> nominalClone(static_cast<TH1D *>(bestUnfolded.Clone((std::string("best_nominal_clone_") + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str())));
        nominalClone->SetDirectory(nullptr);
        nominalClone->SetLineColor(kBlue + 1);
        nominalClone->SetMarkerColor(kBlue + 1);
        nominalClone->SetMarkerStyle(24);
        nominalClone->SetTitle((std::string("Best unfolded comparison, jet p_{T} ") + jetBin.label + ", " + etaBin.label + ";z_{T};counts").c_str());
        nominalClone->Draw("E1");
        bestWithTotal->SetLineColor(kRed + 1);
        bestWithTotal->SetMarkerColor(kRed + 1);
        bestWithTotal->SetMarkerStyle(21);
        bestWithTotal->Draw("E1 SAME");
        TLegend legComp(0.55, 0.72, 0.88, 0.88);
        legComp.SetBorderSize(0);
        legComp.SetFillStyle(0);
        legComp.AddEntry(nominalClone.get(), "Nominal best_unfolded", "lep");
        legComp.AddEntry(bestWithTotal.get(), "Best (total errors)", "lep");
        legComp.Draw();

        // Ratio panel
        cmpCanvas.cd(2);
        std::unique_ptr<TH1D> ratioComp(MakeRatioHistogram(bestWithTotal.get(), &bestUnfolded, (std::string("ratio_best_total_nominal_") + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str(), "total / nominal"));
        if (ratioComp)
        {
            ratioComp->SetLineColor(kRed + 1);
            ratioComp->SetMarkerColor(kRed + 1);
            ratioComp->SetMarkerStyle(20);
            ratioComp->GetYaxis()->SetRangeUser(0.5, 1.5);
            ratioComp->Draw("E1");
            TLine unityR(ratioComp->GetXaxis()->GetXmin(), 1.0, ratioComp->GetXaxis()->GetXmax(), 1.0);
            unityR.SetLineStyle(2);
            unityR.SetLineColor(kBlack);
            unityR.Draw("same");
        }

        SaveCanvas(&cmpCanvas, plotDirectory_ / etaBin.label / jetBin.label / ("compare_best_" + jetBin.label + "_" + EtaToken(etaBin) + ".png"));
    }

    void PerformPriorVariations(const JetBin &jetBin,
                                const EtaBin &etaBin,
                                const ResponseBundle &bundle,
                                int bestIteration,
                                const TH1D &bestUnfolded)
    {
        if (!bundle.response)
            return;
        // Read the original truth prior from the response (read-only).
        // We no longer mutate it; instead, BuildReweightedResponse
        // reweights the response matrix columns to implement a new prior.
        const TH1 *htruth = bundle.response->Htruth();
        if (!htruth)
            return;

        const int nbins = htruth->GetNbinsX();
        std::vector<double> orig(nbins + 1, 0.0);
        double origIntegral = 0.0;
        for (int b = 1; b <= nbins; ++b)
        {
            orig[b] = htruth->GetBinContent(b);
            origIntegral += orig[b];
        }

        // clone original prior for plotting and writing
        std::unique_ptr<TH1D> priorOrig(static_cast<TH1D *>(htruth->Clone(("prior_orig_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str())));
        if (priorOrig)
        {
            priorOrig->SetDirectory(nullptr);
        }

        // Helper lambda: build a RooUnfoldResponse with a reweighted
        // response matrix.  For each truth bin j the weight is
        //     w_j = newPrior_j / origPrior_j
        // Every column j of the response matrix is scaled by w_j.
        // This changes the effective prior that RooUnfoldBayes sees
        // (Htruth = column sums of R) while preserving the per-bin
        // efficiency  eps_j = sum_i R_ij / Htruth_j.
        auto BuildReweightedResponse = [&](const TH1 *newPrior) -> std::unique_ptr<RooUnfoldResponse>
        {
            TH2 *origMat = bundle.response->Hresponse();
            if (!origMat) return nullptr;
            const int nx = origMat->GetNbinsX();
            const int ny = origMat->GetNbinsY();

            // Clone the response matrix and reweight columns
            std::unique_ptr<TH2D> rwMat(static_cast<TH2D *>(origMat->Clone(
                (std::string("rwMat_") + newPrior->GetName()).c_str())));
            rwMat->SetDirectory(nullptr);

            for (int iy = 1; iy <= ny; ++iy)
            {
                const double origVal = orig[iy]; // original prior bin content
                const double newVal  = newPrior->GetBinContent(iy);
                const double weight  = (origVal > 0.0) ? (newVal / origVal) : 0.0;
                for (int ix = 1; ix <= nx; ++ix)
                {
                    rwMat->SetBinContent(ix, iy, origMat->GetBinContent(ix, iy) * weight);
                    rwMat->SetBinError(ix, iy, origMat->GetBinError(ix, iy) * weight);
                }
            }

            // Build new truth = column sums of reweighted matrix
            std::unique_ptr<TH1D> newTruth(static_cast<TH1D *>(bundle.truthTemplate->Clone(
                (std::string("rwTruth_") + newPrior->GetName()).c_str())));
            newTruth->SetDirectory(nullptr);
            newTruth->Reset();
            for (int iy = 1; iy <= ny; ++iy)
            {
                double colSum = 0.0;
                for (int ix = 1; ix <= nx; ++ix)
                    colSum += rwMat->GetBinContent(ix, iy);
                newTruth->SetBinContent(iy, colSum);
            }

            // Build new measured = row sums of reweighted matrix
            std::unique_ptr<TH1D> newMeas(static_cast<TH1D *>(bundle.recoTemplate->Clone(
                (std::string("rwMeas_") + newPrior->GetName()).c_str())));
            newMeas->SetDirectory(nullptr);
            newMeas->Reset();
            for (int ix = 1; ix <= nx; ++ix)
            {
                double rowSum = 0.0;
                for (int iy = 1; iy <= ny; ++iy)
                    rowSum += rwMat->GetBinContent(ix, iy);
                newMeas->SetBinContent(ix, rowSum);
            }

            return std::make_unique<RooUnfoldResponse>(newMeas.get(), newTruth.get(), rwMat.get());
        };

        // 2) Weighted prior: use the (rescaled) unfolded result as new prior
        std::unique_ptr<TH1D> unfW;
        std::unique_ptr<TH1D> priorWeighted;
        {
            const double unfoldedIntegral = bestUnfolded.Integral();
            if (unfoldedIntegral > 0.0)
            {
                // Build the weighted prior histogram
                priorWeighted.reset(static_cast<TH1D *>(htruth->Clone(
                    ("prior_weighted_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str())));
                priorWeighted->SetDirectory(nullptr);
                for (int b = 1; b <= nbins; ++b)
                {
                    const double v = bestUnfolded.GetBinContent(b) *
                                     (origIntegral > 0.0 ? (origIntegral / unfoldedIntegral) : 1.0);
                    priorWeighted->SetBinContent(b, v);
                }

                if (config_.verbose)
                {
                    std::cout << "[prior] weighted prior integral=" << priorWeighted->Integral()
                              << " (origIntegral=" << origIntegral << ")\n";
                }

                auto rwResp = BuildReweightedResponse(priorWeighted.get());
                if (rwResp)
                {
                    RooUnfoldBayes unfoldW(rwResp.get(), bundle.measuredForUnfold.get(), bestIteration);
                    TH1D *unfWRaw = dynamic_cast<TH1D *>(unfoldW.Hreco(RooUnfold::kCovToy));
                    if (unfWRaw)
                    {
                        unfW.reset(static_cast<TH1D *>(unfWRaw->Clone(
                            ("best_unfolded_prior_weighted_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str())));
                        unfW->SetDirectory(nullptr);
                        outputFile_->cd();
                        unfW->Write();

                        if (config_.verbose)
                        {
                            std::cout << "[prior] produced unfolded (weighted) integral=" << unfW->Integral() << "\n";
                            Chi2Result chiW = ComputeChi2(unfW.get(), &bestUnfolded);
                            std::cout << "[prior] chi2(unf_weighted, nominal) reduced=" << chiW.reduced
                                      << " ndf=" << chiW.ndf << "\n";
                        }
                    }
                }
            }
        }

        // 3) Flat prior: constant density across z_T
        std::unique_ptr<TH1D> unfFlat;
        std::unique_ptr<TH1D> priorFlat;
        {
            priorFlat.reset(static_cast<TH1D *>(htruth->Clone(
                ("prior_flat_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str())));
            priorFlat->SetDirectory(nullptr);

            double totalWidth = 0.0;
            for (int b = 1; b <= nbins; ++b)
                totalWidth += htruth->GetBinWidth(b);
            const double flatDensity = (totalWidth > 0.0 ? (origIntegral / totalWidth) : 1.0);
            for (int b = 1; b <= nbins; ++b)
                priorFlat->SetBinContent(b, flatDensity * htruth->GetBinWidth(b));

            if (config_.verbose)
            {
                std::cout << "[prior] flat prior integral=" << priorFlat->Integral()
                          << " (origIntegral=" << origIntegral << ")\n";
            }

            auto rwRespFlat = BuildReweightedResponse(priorFlat.get());
            if (rwRespFlat)
            {
                RooUnfoldBayes unfoldFlat(rwRespFlat.get(), bundle.measuredForUnfold.get(), bestIteration);
                TH1D *unfFlatRaw = dynamic_cast<TH1D *>(unfoldFlat.Hreco(RooUnfold::kCovToy));
                if (unfFlatRaw)
                {
                    unfFlat.reset(static_cast<TH1D *>(unfFlatRaw->Clone(
                        ("best_unfolded_prior_flat_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str())));
                    unfFlat->SetDirectory(nullptr);
                    outputFile_->cd();
                    unfFlat->Write();

                    if (config_.verbose)
                    {
                        std::cout << "[prior] produced unfolded (flat) integral=" << unfFlat->Integral() << "\n";
                        Chi2Result chiF = ComputeChi2(unfFlat.get(), &bestUnfolded);
                        std::cout << "[prior] chi2(unf_flat, nominal) reduced=" << chiF.reduced
                                  << " ndf=" << chiF.ndf << "\n";
                    }
                }
            }
        }

        // write prior histograms to output file
        outputFile_->cd();
        if (priorOrig) priorOrig->Write();
        if (priorFlat) priorFlat->Write();
        if (priorWeighted) priorWeighted->Write();

        // Plot the three prior shapes together: original, flat, weighted
        // Create display clones dividing by bin width so constant density plots flat.
        {
            std::unique_ptr<TH1D> dispOrig;
            std::unique_ptr<TH1D> dispFlat;
            std::unique_ptr<TH1D> dispWeighted;

            if (priorOrig)
            {
                dispOrig.reset(static_cast<TH1D *>(priorOrig->Clone((std::string("prior_orig_display_") + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str())));
                dispOrig->SetDirectory(nullptr);
                for (int b = 1; b <= dispOrig->GetNbinsX(); ++b)
                {
                    const double w = dispOrig->GetBinWidth(b);
                    if (w > 0) dispOrig->SetBinContent(b, priorOrig->GetBinContent(b) / w);
                }
            }
            if (priorFlat)
            {
                dispFlat.reset(static_cast<TH1D *>(priorFlat->Clone((std::string("prior_flat_display_") + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str())));
                dispFlat->SetDirectory(nullptr);
                for (int b = 1; b <= dispFlat->GetNbinsX(); ++b)
                {
                    const double w = dispFlat->GetBinWidth(b);
                    if (w > 0) dispFlat->SetBinContent(b, priorFlat->GetBinContent(b) / w);
                }
            }
            if (priorWeighted)
            {
                dispWeighted.reset(static_cast<TH1D *>(priorWeighted->Clone((std::string("prior_weighted_display_") + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str())));
                dispWeighted->SetDirectory(nullptr);
                for (int b = 1; b <= dispWeighted->GetNbinsX(); ++b)
                {
                    const double w = dispWeighted->GetBinWidth(b);
                    if (w > 0) dispWeighted->SetBinContent(b, priorWeighted->GetBinContent(b) / w);
                }
            }

            TH1D *base = dispOrig ? dispOrig.get() : (dispWeighted ? dispWeighted.get() : dispFlat.get());
            if (base)
            {
                TCanvas cpri(("c_prior_shapes_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str(), "Priors: orig/flat/weighted (density)", 900, 700);
                cpri.SetBottomMargin(0.08);
                cpri.SetRightMargin(0.01);
                cpri.SetTopMargin(0.01);
                cpri.SetLeftMargin(0.12);
                cpri.SetTickx(1);
                cpri.SetTicky(1);
                gStyle->SetOptStat(0);
                base->SetLineColor(kBlack);
                base->SetMarkerStyle(20);
                base->SetMarkerColor(kBlack);
                base->SetTitle("");
                base->GetYaxis()->SetTitle("density (per unit x)");
                base->GetXaxis()->SetTitle("#it{z}_{T}");
                double ymax = base->GetMaximum();
                if (dispFlat)
                {
                    ymax = std::max(ymax, dispFlat->GetMaximum());
                }
                if (dispWeighted)
                {
                    ymax = std::max(ymax, dispWeighted->GetMaximum());
                }
                base->GetYaxis()->SetRangeUser(0.0, ymax * 1.3);
                base->Draw("E1");

                if (dispFlat && dispFlat.get() != base)
                {
                    dispFlat->SetLineColor(kRed);
                    dispFlat->SetMarkerColor(kRed);
                    dispFlat->Draw("HIST same");
                }
                if (dispWeighted && dispWeighted.get() != base)
                {
                    dispWeighted->SetLineColor(kBlue);
                    dispWeighted->SetMarkerColor(kBlue);
                    dispWeighted->Draw("HIST same");
                }

                TLegend leg(0.15, 0.8, 0.52, 0.96);
                leg.SetBorderSize(0);
                leg.SetFillStyle(0);
                if (dispOrig) leg.AddEntry(dispOrig.get(), "Original prior (truth)", "lep");
                if (dispFlat) leg.AddEntry(dispFlat.get(), "Flat prior (density)", "l");
                if (dispWeighted) leg.AddEntry(dispWeighted.get(), "Weighted prior (unfold)", "l");
                leg.Draw();

                TLatex* label = new TLatex();
                label->SetNDC();
                label->SetTextSize(0.035);
                label->DrawLatex(0.5, 0.92, (std::string("Jet p_{T}: ") + jetBin.label + " GeV/c, #eta range: " + etaBin.label).c_str());

                SaveCanvas(&cpri, plotDirectory_ / etaBin.label / jetBin.label / ("priors_" + jetBin.label + "_" + EtaToken(etaBin) + ".png"));
            }
        }

        // Comparison plot: nominal vs flat and weighted prior unfolds
        {
            TCanvas cmpPrior(("c_prior_compare_" + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str(), "prior compare", 700, 900);
            cmpPrior.Divide(1, 2);
            cmpPrior.cd(1);
            cmpPrior.cd(1)->SetTickx(1);
            cmpPrior.cd(1)->SetTicky(1);
            cmpPrior.cd(1)->SetBottomMargin(0.00);
            cmpPrior.cd(1)->SetRightMargin(0.01);
            cmpPrior.cd(1)->SetTopMargin(0.02);
            gStyle->SetOptStat(0);
            // draw nominal
            std::unique_ptr<TH1D> nominalClone(static_cast<TH1D *>(bestUnfolded.Clone((std::string("nominal_clone_") + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str())));
            nominalClone->SetDirectory(nullptr);
            nominalClone->SetLineColor(kBlue + 1);
            nominalClone->SetMarkerColor(kBlue + 1);
            nominalClone->SetMarkerStyle(24);
            nominalClone->GetYaxis()->SetRangeUser(-1.0, nominalClone->GetMaximum() * 1.5);
            nominalClone->GetYaxis()->SetTitle("counts");
            // nominalClone->SetTitle((std::string("Prior variation comparison, jet p_{T} ") + jetBin.label + ", " + etaBin.label + ";z_{T};counts").c_str());
            nominalClone->SetTitle("");
            nominalClone->Draw("E1");
            TLatex* label = new TLatex();
            label->SetNDC();
            label->SetTextSize(0.04);
            label->DrawLatex(0.15, 0.92, (std::string("Jet p_{T}: ") + jetBin.label + " GeV/c, #eta range: " + etaBin.label).c_str());
            if (unfFlat)
            {
                unfFlat->SetLineColor(kRed + 1);
                unfFlat->SetMarkerColor(kRed + 1);
                unfFlat->SetMarkerStyle(21);
                unfFlat->Draw("E1 SAME");
            }
            if (unfW)
            {
                unfW->SetLineColor(kGreen + 2);
                unfW->SetMarkerColor(kGreen + 2);
                unfW->SetMarkerStyle(20);
                unfW->Draw("E1 SAME");
            }
            TLegend leg(0.55, 0.72, 0.88, 0.88);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.AddEntry(nominalClone.get(), "Nominal best unfolded", "lep");
            if (unfFlat) leg.AddEntry(unfFlat.get(), "Prior: flat (density)", "lep");
            if (unfW) leg.AddEntry(unfW.get(), "Prior: weighted (unfold)", "lep");
            leg.Draw();

            // ratio panel
            cmpPrior.cd(2);
            cmpPrior.cd(2)->SetTickx(1);
            cmpPrior.cd(2)->SetTicky(1);
            cmpPrior.cd(2)->SetTopMargin(0.0);
            cmpPrior.cd(2)->SetRightMargin(0.01);
            std::unique_ptr<TH1D> rFlat(unfFlat ? MakeRatioHistogram(unfFlat.get(), &bestUnfolded, (std::string("ratio_prior_flat_") + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str(), "flat / nominal") : nullptr);
            std::unique_ptr<TH1D> rW(unfW ? MakeRatioHistogram(unfW.get(), &bestUnfolded, (std::string("ratio_prior_weighted_") + SanitizeLabel(jetBin.label) + "_" + EtaToken(etaBin)).c_str(), "weighted / nominal") : nullptr);
            TH1D *frame = nullptr;
            if (rFlat) frame = static_cast<TH1D *>(rFlat->Clone("frame_prior_ratio"));
            else if (rW) frame = static_cast<TH1D *>(rW->Clone("frame_prior_ratio"));
            if (frame)
            {
                frame->SetDirectory(nullptr);
                frame->SetTitle(";z_{T};prior / nominal");
                frame->GetYaxis()->SetRangeUser(0.51, 1.49);
                frame->Draw("AXIS");
                if (rFlat) { rFlat->SetLineColor(kRed + 1); rFlat->SetMarkerColor(kRed + 1); rFlat->Draw("E1 SAME"); }
                if (rW) { rW->SetLineColor(kGreen + 2); rW->SetMarkerColor(kGreen + 2); rW->Draw("E1 SAME"); }
                TLine* l = new TLine(frame->GetXaxis()->GetXmin(), 1.0, frame->GetXaxis()->GetXmax(), 1.0);
                l->SetLineStyle(2);
                l->Draw();
            }

            SaveCanvas(&cmpPrior, plotDirectory_ / etaBin.label / jetBin.label / ("prior_variation_compare_" + jetBin.label + "_" + EtaToken(etaBin) + ".png"));
            outputFile_->cd();
            // also write ratio histograms
            if (rFlat) rFlat->Write();
            if (rW) rW->Write();
        }
    }

    D0UnfoldConfig config_;
    std::unique_ptr<TFile> responseInput_;
    std::unique_ptr<TFile> triggeredResponseInput_;
    std::unique_ptr<TFile> outputFile_;
    TTree *responseTree_ = nullptr;
    TTree *triggeredResponseTree_ = nullptr;
    std::filesystem::path plotDirectory_;
};
// Definitions for member helper functions moved out of BuildEtaBins
void D0ZTUnfolder::FillTruthProjectionUncut(TTree *tree,
                                           const JetBin &jetBin,
                                           const EtaBin &etaBin,
                                           double scale,
                                           TH1D *truthUncut,
                                           const TH1 *detLevelWeights,
                                           const TH1 *truthLevelWeights) const
{
    float d0ZDet = 0.0f;
    float d0ZMc = 0.0f;
    float jetPtDet = 0.0f;
    float jetPtMc = 0.0f;
    float d0EtaDet = 0.0f;
    float d0EtaMc = 0.0f;
    float jetNconstDet = 0.0f;
    float jetNconstMc = 0.0f;
    float eventWeight = 1.0f;

    tree->SetBranchAddress("d0_z_det", &d0ZDet);
    tree->SetBranchAddress("d0_z_mc", &d0ZMc);
    tree->SetBranchAddress("jet_pt_det", &jetPtDet);
    tree->SetBranchAddress("jet_pt_mc", &jetPtMc);
    tree->SetBranchAddress("d0_eta_det", &d0EtaDet);
    tree->SetBranchAddress("d0_eta_mc", &d0EtaMc);
    tree->SetBranchAddress("jet_nconst_det", &jetNconstDet);
    tree->SetBranchAddress("jet_nconst_mc", &jetNconstMc);

    if (tree->GetBranch("event_weight"))
        tree->SetBranchAddress("event_weight", &eventWeight);
    else if (tree->GetBranch("weight"))
        tree->SetBranchAddress("weight", &eventWeight);
    else if (tree->GetBranch("evtWeight"))
        tree->SetBranchAddress("evtWeight", &eventWeight);
    else if (tree->GetBranch("eventWeight"))
        tree->SetBranchAddress("eventWeight", &eventWeight);
    else if (tree->GetBranch("totalWeight"))
        tree->SetBranchAddress("totalWeight", &eventWeight);

    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry)
    {
        tree->GetEntry(entry);

        double jetPtDetCorr = jetPtDet;
        double d0ZDetCorr = d0ZDet;
        ApplyResponseRecoCorrection(config_, jetPtDet, d0ZDet, jetPtDetCorr, d0ZDetCorr);

        // Keep only basic RM-level filters (matching/eta/nconst). Require
        // the GEN-level jet pT to fall in the current jet bin so this
        // histogram represents the gen-level distribution for this jet-pt
        // bin prior to detector-level pT selection (kinematic efficiency).
        const bool baseCuts = InEtaBin(d0EtaDet, etaBin) &&
                              jetNconstDet > config_.minJetConstituents &&
                              jetNconstMc > config_.minJetConstituents;
        const bool genInBin = (jetPtMc >= jetBin.min && jetPtMc < jetBin.max);
        if (!baseCuts || !genInBin)
            continue;

        double totalWeight = scale * eventWeight;
        if (detLevelWeights)
        {
            totalWeight *= LookupWeight1D(detLevelWeights, d0ZDetCorr);
        }
        else if (truthLevelWeights)
        {
            totalWeight *= LookupWeight1D(truthLevelWeights, d0ZMc);
        }

        if (truthUncut)
        {
            truthUncut->Fill(d0ZMc, totalWeight);
        }
    }
}

void D0ZTUnfolder::FillTruthProjectionUncut2D(TTree *tree,
                                             const EtaBin &etaBin,
                                             double scale,
                                             TH2D *truthUncut,
                                             const TH2 *detLevelWeights,
                                             const TH2 *truthLevelWeights) const
{
    float d0ZDet = 0.0f;
    float d0ZMc = 0.0f;
    float jetPtDet = 0.0f;
    float jetPtMc = 0.0f;
    float d0EtaDet = 0.0f;
    float d0EtaMc = 0.0f;
    float jetNconstDet = 0.0f;
    float jetNconstMc = 0.0f;
    float eventWeight = 1.0f;

    tree->SetBranchAddress("d0_z_det", &d0ZDet);
    tree->SetBranchAddress("d0_z_mc", &d0ZMc);
    tree->SetBranchAddress("jet_pt_det", &jetPtDet);
    tree->SetBranchAddress("jet_pt_mc", &jetPtMc);
    tree->SetBranchAddress("d0_eta_det", &d0EtaDet);
    tree->SetBranchAddress("d0_eta_mc", &d0EtaMc);
    tree->SetBranchAddress("jet_nconst_det", &jetNconstDet);
    tree->SetBranchAddress("jet_nconst_mc", &jetNconstMc);

    if (tree->GetBranch("event_weight"))
        tree->SetBranchAddress("event_weight", &eventWeight);
    else if (tree->GetBranch("weight"))
        tree->SetBranchAddress("weight", &eventWeight);
    else if (tree->GetBranch("evtWeight"))
        tree->SetBranchAddress("evtWeight", &eventWeight);
    else if (tree->GetBranch("eventWeight"))
        tree->SetBranchAddress("eventWeight", &eventWeight);
    else if (tree->GetBranch("totalWeight"))
        tree->SetBranchAddress("totalWeight", &eventWeight);

    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry)
    {
        tree->GetEntry(entry);

        double jetPtDetCorr = jetPtDet;
        double d0ZDetCorr = d0ZDet;
        ApplyResponseRecoCorrection(config_, jetPtDet, d0ZDet, jetPtDetCorr, d0ZDetCorr);

        const bool baseCuts = InEtaBin(d0EtaDet, etaBin) &&
                              jetNconstDet > config_.minJetConstituents &&
                              jetNconstMc > config_.minJetConstituents;
        if (!baseCuts)
            continue;

        double totalWeight = scale * eventWeight;
        if (detLevelWeights)
        {
            totalWeight *= LookupWeight2D(detLevelWeights, jetPtDetCorr, d0ZDetCorr);
        }
        else if (truthLevelWeights)
        {
            totalWeight *= LookupWeight2D(truthLevelWeights, jetPtMc, d0ZMc);
        }

        if (truthUncut)
        {
            truthUncut->Fill(jetPtMc, d0ZMc, totalWeight);
        }
    }
}

} // namespace

void unfold_d0_zt_lhcb(
                              const std::string &measuredFilePattern = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-04-07_pPb/TagZHistograms_%s.root",
                              const std::string &responseFile = "/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pPb_MC_54plus73_EPOS_response.root",
                              const std::string &outputFile = "d0_unfolded_zt.root",
                              const std::string &outputFolderTag = "",
                            //   const std::string &jetBinsCsv = "10_15,15_20,20_30",
                              const std::string &jetBinsCsv = "10_15,15_20,20_30,30_100",
                              const std::string &rapidityLabelsCsv = "2.5-3.0,3.0-3.5,3.5-4.0",
                              const std::string &triggeredResponseFile = "",
                              double triggeredScale = 1.0,
                              int nIterations = 5,
                              bool doClosure = true,
                              bool verbose = true,
                              bool do2D = false,
                              bool doPriorVariations = true,
                              bool doResponseDetLevelWeighting = true,
                              bool doResponseTruthLevelWeighting = true,
                              bool applyResponseSmearing = true,
                              double responseScaleFactor = 0.96,
                              double responseSmearFactor = 0.07)
{
    D0UnfoldConfig config;
    config.measuredFilePattern = measuredFilePattern;
    config.responseFile = responseFile;
    config.triggeredResponseFile = triggeredResponseFile;
    config.outputFile = outputFile;
    config.triggeredScale = triggeredScale;
    config.nIterations = nIterations;
    config.doClosure = doClosure;
    config.verbose = verbose;
    config.do2D = do2D;
    config.jetBins = ParseJetBins(jetBinsCsv);
    config.rapidityBinLabels = SplitCsv(rapidityLabelsCsv);
    config.doPriorVariations = doPriorVariations;
    config.doResponseDetLevelWeighting = doResponseDetLevelWeighting;
    config.doResponseTruthLevelWeighting = doResponseTruthLevelWeighting;
    config.applyResponseSmearing = applyResponseSmearing;
    config.responseScaleFactor = responseScaleFactor;
    config.responseSmearFactor = responseSmearFactor;
    config.outputFolderTag = outputFolderTag;

    D0ZTUnfolder unfolder(std::move(config));
    unfolder.Run();
}

void unfold_d0_zt_lhcb_example()
{
    unfold_d0_zt_lhcb();
}