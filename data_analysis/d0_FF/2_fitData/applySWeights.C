#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <filesystem>
#include <iomanip>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"

// -----------------------------------------------------------------------------
// applySWeights.C
// -----------------------------------------------------------------------------
// Macro to apply sPlot weights (produced by MassFitterSimple) to zT distributions
// in the same rapidity (y) and jet pT bins.
//
// Assumptions:
//  * Original physics TTree is named "FragmNtuple" and contains branches:
//      - tagY (rapidity)
//      - jetPt (jet transverse momentum) [optional]
//      - zT   (fragmentation variable)   [REQUIRED]
//  * sWeights ROOT file contains TTrees named:
//      Mass:  sWeights_mass_y<iy>_<yLo>to<yHi>_j<ij>_<jLo>to<jHi>
//      IP:    sWeights_ipchi2_y<iy>_<yLo>to<yHi>_j<ij>_<jLo>to<jHi>
//    (matching the naming convention used in MassFitterSimple). For inclusive (no jetPt) case
//    the pattern uses j-1 internally, but we keep loops consistent with yBins and jetPtBins.
//  * Event ordering for a given bin in the sWeights TTree matches the ordering of entries
//    obtained by applying the identical selection cut on the original TTree.
//
// Usage example inside ROOT:
//   .L applySWeights.C++
//   applySWeights("/path/to/data.root", \
//                 "/path/to/SWeights_D0_DATA_pPb.root", \
//                 false,   // isMC
//                 true,    // doMassWeights
//                 true,    // doIPWeights
//                 "zT",    // zT branch name
//                 50, 0.0, 1.0); // histogram binning
//
// Output:
//   Creates an output ROOT file in working directory named:
//      WeightedZT_<DATA|MC>_weights.root
//   containing per-bin histograms:
//      hZT_sig_y<iy>_j<ij>, hZT_bkg_..., hZT_prompt_..., hZT_nonprompt_...
//   (Only those enabled by doMassWeights / doIPWeights.)
// -----------------------------------------------------------------------------

struct BinSpec {
    int iy;
    int ij;
    double yLo, yHi;
    double jLo, jHi;
    std::string suffix; // _y<iy>_<yLo>to<yHi>_j<ij>_<jLo>to<jHi>
    std::string suffix2; // _y<iy>_j<ij>
};

// Helper: format double without trailing zeros (for names) keeping original bins identical to MassFitterSimple
static std::string fmtDouble(double v) {
    std::ostringstream os; os<<std::fixed<<std::setprecision(1)<<v; return os.str();
}

// Build bin list replicating MassFitterSimple logic
static std::vector<BinSpec> buildBins(const std::vector<double>& yBins, const std::vector<double>& jetPtBins, bool haveJetPt) {
    std::vector<BinSpec> bins;
    for (size_t iy=0; iy+1<yBins.size(); ++iy) {
        double yLo=yBins[iy], yHi=yBins[iy+1];
        size_t jStart=0, jEnd = (haveJetPt ? jetPtBins.size()-1 : 1);
        for (size_t ij=jStart; ij<jEnd; ++ij) {
            double jLo = haveJetPt ? jetPtBins[ij] : jetPtBins.front();
            double jHi = haveJetPt ? jetPtBins[ij+1] : jetPtBins.back();
            std::ostringstream suf;
            suf << "_y" << iy << "_" << yLo << "to" << yHi;
            // Match MassFitterSimple suffix: always includes _j<idx>, using -1 when no jetPt
            if (haveJetPt) {
                suf << "_j" << ij << "_" << jLo << "to" << jHi;
            } else {
                suf << "_j-1"; // inclusive mode suffix in MassFitterSimple
            }
            std::ostringstream suf2;
            suf2 << "_y" << iy;
            if (haveJetPt) {
                suf2 << "_j" << ij;
            } else {
                suf2 << "_j-1";
            }
            bins.push_back({(int)iy,(int)ij,yLo,yHi,jLo,jHi,suf.str(), suf2.str()});
        }
    }
    return bins;
}

// Core function
void applySWeights(TString inputDataFile = "/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/57_FF_pPb_DATA_filtered.root",
                   TString sWeightsFile = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/20251117_D0_FF_DATA_pPb/SWeights_D0_DATA_pPb.root",
                   bool isMC=false,
                   bool doMassWeights=true,
                   bool doIPWeights=true,
                   TString zTBranchName="tagZ",
                   int zTBins=50,
                   double zTMin=0.0,
                   double zTMax=1.0,
                   bool debug=true) {

    std::cout << "[applySWeights] Starting with data: " << inputDataFile << " sWeights: " << sWeightsFile << std::endl;

    // Open input data file
    std::unique_ptr<TFile> dataFile(TFile::Open(inputDataFile, "READ"));
    if (!dataFile || dataFile->IsZombie()) { std::cerr << "[ERROR] Cannot open data file." << std::endl; return; }
    TTree* tree = dynamic_cast<TTree*>(dataFile->Get("FragmNtuple"));
    if (!tree) { std::cerr << "[ERROR] TTree 'FragmNtuple' not found in data file." << std::endl; return; }
    Long64_t nEntriesTotal = tree->GetEntries();
    if (debug) std::cout << "[DEBUG] FragmNtuple entries: " << nEntriesTotal << std::endl;

    // Check required branches
    if (!tree->GetBranch("tagY")) { std::cerr << "[ERROR] Branch tagY missing." << std::endl; return; }
    bool haveJetPt = (tree->GetBranch("jetPt")!=nullptr);
    if (!haveJetPt) std::cout << "[WARN] jetPt branch not found. Will treat as inclusive." << std::endl;
    if (!tree->GetBranch(zTBranchName)) { std::cerr << "[ERROR] zT branch '" << zTBranchName << "' not found." << std::endl; return; }

    // Quick global diagnostics: min/max for tagY and jetPt, sample a few entries
    double dbg_tagY_min= 1e9, dbg_tagY_max=-1e9;
    double dbg_jet_min= 1e9, dbg_jet_max=-1e9;
    double dbg_zT_min = 1e9, dbg_zT_max =-1e9;
    Float_t tmpY=0, tmpJet=0, tmpZT=0; // match branch types (Float_t)
    tree->SetBranchAddress("tagY", &tmpY);
    if (haveJetPt) tree->SetBranchAddress("jetPt", &tmpJet);
    tree->SetBranchAddress(zTBranchName, &tmpZT);
    for (Long64_t i=0;i<nEntriesTotal;++i) {
        tree->GetEntry(i);
        if (tmpY<dbg_tagY_min) dbg_tagY_min=tmpY;
        if (tmpY>dbg_tagY_max) dbg_tagY_max=tmpY;
        if (haveJetPt) {
            if (tmpJet<dbg_jet_min) dbg_jet_min=tmpJet;
            if (tmpJet>dbg_jet_max) dbg_jet_max=tmpJet;
        }
        if (tmpZT<dbg_zT_min) dbg_zT_min=tmpZT;
        if (tmpZT>dbg_zT_max) dbg_zT_max=tmpZT;
    }
    if (debug) {
        std::cout << std::fixed;
        std::cout << "[DEBUG] tagY range in tree: [" << dbg_tagY_min << ", " << dbg_tagY_max << "]\n";
        if (haveJetPt) std::cout << "[DEBUG] jetPt range in tree: [" << dbg_jet_min << ", " << dbg_jet_max << "]\n";
        std::cout << "[DEBUG] zT(" << zTBranchName << ") range in tree: [" << dbg_zT_min << ", " << dbg_zT_max << "]\n";
    }

    // Open sWeights file
    std::unique_ptr<TFile> swFile(TFile::Open(sWeightsFile, "READ"));
    if (!swFile || swFile->IsZombie()) { std::cerr << "[ERROR] Cannot open sWeights file." << std::endl; return; }

    // Define binning (duplicate from MassFitterSimple)
    std::vector<double> yBins = {2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0};
    std::vector<double> jetPtBins = {5,10,15,20,30,50};

    // Build bin list
    auto bins = buildBins(yBins, jetPtBins, haveJetPt);
    if (debug) {
        std::cout << "[DEBUG] Using y bins:";
        for (size_t i=0;i<yBins.size();++i) std::cout << (i?", ":" ") << yBins[i];
        std::cout << "\n";
        if (haveJetPt) {
            std::cout << "[DEBUG] Using jetPt bins:";
            for (size_t i=0;i<jetPtBins.size();++i) std::cout << (i?", ":" ") << jetPtBins[i];
            std::cout << "\n";
        } else {
            std::cout << "[DEBUG] Inclusive over jetPt (no jetPt branch)\n";
        }
    }

    // Prepare output file
    std::string outName = std::string("WeightedZT_") + (isMC?"MC":"DATA") + "_weights.root";
    std::unique_ptr<TFile> outFile(TFile::Open(outName.c_str(), "RECREATE"));
    if (!outFile || outFile->IsZombie()) { std::cerr << "[ERROR] Cannot create output file." << std::endl; return; }

    // Set up branch addresses for faster reading
    Float_t b_tagY=0.0f, b_jetPt=0.0f, b_zT=0.0f; // match branch types (Float_t)
    tree->SetBranchAddress("tagY", &b_tagY);
    if (haveJetPt) tree->SetBranchAddress("jetPt", &b_jetPt);
    tree->SetBranchAddress(zTBranchName, &b_zT);

    // Loop bins
    for (const auto& b : bins) {
        if (debug) {
            std::cout << "[DEBUG] Bin " << b.suffix << ": y in [" << b.yLo << ", " << b.yHi << ")";
            if (haveJetPt) std::cout << ", jetPt in [" << b.jLo << ", " << b.jHi << ")";
            std::cout << std::endl;
        }
        // Selection cut logic replicating MassFitterSimple
        // We evaluate selection manually rather than using TTree::Draw.
        std::vector<Long64_t> selectedIndices;
        Long64_t yOnlyCount=0, yAndJetCount=0;
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i=0;i<nEntries;++i) {
            tree->GetEntry(i);
            bool inY = (b_tagY >= b.yLo && b_tagY < b.yHi);
            if (inY) ++yOnlyCount;
            bool inJet = (!haveJetPt || (b_jetPt >= b.jLo && b_jetPt < b.jHi));
            if (inY && inJet) { ++yAndJetCount; selectedIndices.push_back(i); }
        }
        if (selectedIndices.empty()) {
            std::cout << "[ApplySWeights] Bin " << b.suffix << " empty. Skipping."
                      << " (y-only count=" << yOnlyCount << ", y+jet count=" << yAndJetCount << ")" << std::endl;
            if (debug) {
                // Additional hints if ranges are outside global ranges
                if (b.yHi <= dbg_tagY_min || b.yLo >= dbg_tagY_max) {
                    std::cout << "[DEBUG] Note: Bin y-range is outside tagY range in file." << std::endl;
                }
                if (haveJetPt && (b.jHi <= dbg_jet_min || b.jLo >= dbg_jet_max)) {
                    std::cout << "[DEBUG] Note: Bin jetPt-range is outside jetPt range in file." << std::endl;
                }
            }
            continue;
        }

        // Retrieve sWeights trees
        TTree* tMass = nullptr; TTree* tIP = nullptr;
        if (doMassWeights) {
            std::string tMassName = std::string("sWeights_mass") + b.suffix2;
            tMass = dynamic_cast<TTree*>(swFile->Get(tMassName.c_str()));
            if (!tMass) {
                std::cerr << "[WARN] Mass sWeights tree missing: " << tMassName << std::endl;
            } else if (debug) {
                std::cout << "[DEBUG] Mass sWeights entries for " << b.suffix2 << ": " << tMass->GetEntries() << std::endl;
            }
        }
        if (doIPWeights) {
            std::string tIPName = std::string("sWeights_ipchi2") + b.suffix2;
            tIP = dynamic_cast<TTree*>(swFile->Get(tIPName.c_str()));
            if (!tIP) {
                std::cerr << "[WARN] IP sWeights tree missing: " << tIPName << std::endl;
            } else if (debug) {
                std::cout << "[DEBUG] IP sWeights entries for " << b.suffix2 << ": " << tIP->GetEntries() << std::endl;
            }
        }

        // Histograms
        std::unique_ptr<TH1D> hSig, hBkg, hPrompt, hNonPrompt;
        if (tMass) {
            hSig.reset(new TH1D(("hZT_sig"+b.suffix).c_str(), ("zT signal " + b.suffix).c_str(), zTBins, zTMin, zTMax));
            hBkg.reset(new TH1D(("hZT_bkg"+b.suffix).c_str(), ("zT background " + b.suffix).c_str(), zTBins, zTMin, zTMax));
        }
        if (tIP) {
            hPrompt.reset(new TH1D(("hZT_prompt"+b.suffix).c_str(), ("zT prompt " + b.suffix).c_str(), zTBins, zTMin, zTMax));
            hNonPrompt.reset(new TH1D(("hZT_nonprompt"+b.suffix).c_str(), ("zT nonprompt " + b.suffix).c_str(), zTBins, zTMin, zTMax));
        }

        // Prepare sWeight arrays
        std::vector<double> w_sig, w_bkg, w_prompt, w_nonprompt;
        if (tMass) {
            // Only weights are needed for filling; match likely types
            Double_t wSig=0, wBkg_=0;
            tMass->SetBranchAddress("w_sig", &wSig);
            tMass->SetBranchAddress("w_bkg", &wBkg_);
            Long64_t nW = tMass->GetEntries();
            if (nW != (Long64_t)selectedIndices.size()) {
                std::cerr << "[WARN] Mass weights count ("<<nW<<") != selected events ("<<selectedIndices.size()<<") for bin "<<b.suffix<<". Will truncate to min." << std::endl;
            }
            Long64_t nUse = std::min<Long64_t>(nW, selectedIndices.size());
            w_sig.reserve(nUse); w_bkg.reserve(nUse);
            for (Long64_t i=0;i<nUse;++i) { tMass->GetEntry(i); w_sig.push_back(wSig); w_bkg.push_back(wBkg_); }
        }
        if (tIP) {
            // Only weights are needed for filling
            Double_t wP=0, wNP=0;
            tIP->SetBranchAddress("w_prompt", &wP);
            tIP->SetBranchAddress("w_nonprompt", &wNP);
            Long64_t nW = tIP->GetEntries();
            if (nW != (Long64_t)selectedIndices.size()) {
                std::cerr << "[WARN] IP weights count ("<<nW<<") != selected events ("<<selectedIndices.size()<<") for bin "<<b.suffix<<". Will truncate to min." << std::endl;
            }
            Long64_t nUse = std::min<Long64_t>(nW, selectedIndices.size());
            w_prompt.reserve(nUse); w_nonprompt.reserve(nUse);
            for (Long64_t i=0;i<nUse;++i) { tIP->GetEntry(i); w_prompt.push_back(wP); w_nonprompt.push_back(wNP); }
        }

        // Fill histograms by iterating selected indices and mapping weights by relative order
        size_t nSelUse = std::min({selectedIndices.size(), w_sig.size(), w_bkg.size()});
        size_t nSelUseIP = std::min({selectedIndices.size(), w_prompt.size(), w_nonprompt.size()});

        // Mass weights application
        if (tMass) {
            for (size_t k=0;k<nSelUse;++k) {
                tree->GetEntry(selectedIndices[k]);
                if (hSig) hSig->Fill(b_zT, w_sig[k]);
                if (hBkg) hBkg->Fill(b_zT, w_bkg[k]);
            }
        }
        // IP weights application
        if (tIP) {
            for (size_t k=0;k<nSelUseIP;++k) {
                tree->GetEntry(selectedIndices[k]);
                if (hPrompt) hPrompt->Fill(b_zT, w_prompt[k]);
                if (hNonPrompt) hNonPrompt->Fill(b_zT, w_nonprompt[k]);
            }
        }

        // Write histograms
        outFile->cd();
        if (hSig) hSig->Write();
        if (hBkg) hBkg->Write();
        if (hPrompt) hPrompt->Write();
        if (hNonPrompt) hNonPrompt->Write();

        std::cout << "[applySWeights] Bin " << b.suffix << " done. Selected=" << selectedIndices.size()
                  << " zT filled (mass/IP): " << nSelUse << "/" << nSelUseIP << std::endl;
    }

    outFile->Write();
    outFile->Close();
    std::cout << "[applySWeights] Finished. Output written to " << outName << std::endl;
}
