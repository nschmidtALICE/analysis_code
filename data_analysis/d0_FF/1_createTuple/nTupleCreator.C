#include "nTupleCreator.h"
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TParameter.h>
#include <TLegend.h>
#include <TF1.h>
#include <TLatex.h>
#include <TLine.h>
#include <TChain.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <iomanip>

int nTagPartLvlEvt = 0;
int nTagDetLvlEvt = 0;
double lastPV = 0;

// Connect branches to input tree for direct access
// Declare ALL branch pointers before the event loop
double *evt_n = nullptr;
std::vector<float> *tag_pid_vec = nullptr;
double evt_spd = 0;
std::vector<int> *mctag_pid = nullptr;
std::vector<int> *tag_l0_tos0 = nullptr;
std::vector<int> *tag_hlt1_tos0 = nullptr;
std::vector<int> *tag_idx_jet = nullptr;
std::vector<int> *tag_genAssocType = nullptr;
std::vector<int> *tag_idx_pvr = nullptr;
std::vector<int> *tag_idx_prt0 = nullptr;
std::vector<int> *tag_idx_prt1 = nullptr;
std::vector<int> *tag_idx_MCtag = nullptr;
std::vector<int> *mctag_isPrimary = nullptr;
std::vector<int> *prt_pid = nullptr;
std::vector<float> *tag_dR_jet = nullptr;
std::vector<int> *tag_n_tagsRnd = nullptr;
std::vector<float> *tag_decVtxChi2 = nullptr;
std::vector<float> *pvr_z = nullptr;
std::vector<float> *tag_z = nullptr;
std::vector<int> *jet_n_tagsComb = nullptr;
std::vector<int> *jet_n_tagsUnique = nullptr;
std::vector<int> *jet_n_neu = nullptr;
std::vector<int> *jet_n_chr = nullptr;
std::vector<float> *prt_pnn_k = nullptr;
std::vector<float> *prt_prb_ghost = nullptr;
std::vector<float> *prt_trckChi2 = nullptr;
std::vector<float> *prt_pnn_pi = nullptr;
std::vector<int> *mctag_idx_jet = nullptr;
std::vector<int> *mctag_idx_Dettag = nullptr;
std::vector<int> *mcjet_n_neu = nullptr;
std::vector<int> *mcjet_n_chr = nullptr;
std::vector<int> *mctag_idx_prt1 = nullptr;
std::vector<int> *mctag_idx_prt2 = nullptr;
std::vector<int> *mctag_idx_pvr = nullptr;

std::vector<int> *mcprt_pid = nullptr;
std::vector<float> *tag_DOCA1 = nullptr;
std::vector<float> *tag_DOCA2 = nullptr;
std::vector<float> *tag_DOCA3 = nullptr;
std::vector<float> *mcjet_px = nullptr;
std::vector<float> *mcjet_py = nullptr;
std::vector<float> *mcjet_pz = nullptr;
std::vector<float> *mcjet_e = nullptr;
std::vector<float> *mctag_px = nullptr;
std::vector<float> *mctag_py = nullptr;
std::vector<float> *mctag_pz = nullptr;

std::vector<double> *mcpvr_z = nullptr;
std::vector<double> *mctag_z = nullptr;
std::vector<float> *mctag_e = nullptr;
std::vector<float> *mcprt_px = nullptr;
std::vector<float> *mcprt_py = nullptr;
std::vector<float> *mcprt_pz = nullptr;
std::vector<float> *mcprt_e = nullptr;
std::vector<float> *jet_px = nullptr;
std::vector<float> *jet_py = nullptr;
std::vector<float> *jet_pz = nullptr;
std::vector<float> *jet_e = nullptr;
std::vector<float> *tag_px = nullptr;
std::vector<float> *tag_py = nullptr;
std::vector<float> *tag_pz = nullptr;
std::vector<float> *tag_e = nullptr;
std::vector<float> *prt_px = nullptr;
std::vector<float> *prt_py = nullptr;
std::vector<float> *prt_pz = nullptr;
std::vector<float> *prt_e = nullptr;
double pvrNumber = 0;

void debugContinue(int line, const std::string &reason)
{
    // std::cout << "DEBUG: Skipping at line " << line << ", reason: " << reason << std::endl;
}

void printLineNumber(int line)
{
    std::cout << "Line " << line << std::endl;
}

int findClosest(const std::vector<std::vector<double>> &matrix, int row, int version)
{
    if (matrix.empty() || row >= static_cast<int>(matrix.size()))
    {
        return -1;
    }

    const auto &rowValues = matrix[row];
    if (rowValues.empty())
    {
        return -1;
    }

    int closestIndex = -1;
    double minDistance = 1e9;

    for (size_t col = 0; col < rowValues.size(); ++col)
    {
        double distance = rowValues[col];
        if (distance < minDistance)
        {
            minDistance = distance;
            closestIndex = col;
        }
    }

    // Apply different matching criteria based on version
    if (version == 0)
    {
        // Return closest match regardless of distance
        return closestIndex;
    }
    else if (version == 1)
    {
        // Return match only if distance is below threshold
        return (minDistance < 0.4) ? closestIndex : -1;
    }

    return closestIndex;
}

TLorentzVector getFourVector(const std::string &type, int index)
{
    if (index < 0)
    {
        std::cerr << "This particle does not exist - please check and debug, type: " << type << std::endl;
        return TLorentzVector(0, 0, 0, 0);
    }

    double px = 0, py = 0, pz = 0, e = 0;

    if (type == "mcjet")
    {

        if (index < (int)mcjet_px->size())
        {
            px = (*mcjet_px)[index];
            py = (*mcjet_py)[index];
            pz = (*mcjet_pz)[index];
            e = (*mcjet_e)[index];
        }
    }
    else if (type == "mctag")
    {

        if (index < (int)mctag_px->size())
        {
            px = (*mctag_px)[index];
            py = (*mctag_py)[index];
            pz = (*mctag_pz)[index];
            e = (*mctag_e)[index];
        }
    }
    else if (type == "mcprt")
    {

        if (index < (int)mcprt_px->size())
        {
            px = (*mcprt_px)[index];
            py = (*mcprt_py)[index];
            pz = (*mcprt_pz)[index];
            e = (*mcprt_e)[index];
        }
    }
    else if (type == "jet")
    {
        if (index < (int)jet_px->size())
        {
            px = (*jet_px)[index];
            py = (*jet_py)[index];
            pz = (*jet_pz)[index];
            e = (*jet_e)[index];
        }
    }
    else if (type == "tag")
    {
        if (index < (int)tag_px->size())
        {
            px = (*tag_px)[index];
            py = (*tag_py)[index];
            pz = (*tag_pz)[index];
            e = (*tag_e)[index];
        }
    }
    else if (type == "prt")
    {

        if (index < (int)prt_px->size())
        {
            px = (*prt_px)[index];
            py = (*prt_py)[index];
            pz = (*prt_pz)[index];
            e = (*prt_e)[index];
        }
    }

    return TLorentzVector(px, py, pz, e);
}
double calcLifeTime(double mass, double pvrZ, double decayZ, double pz, bool debug)
{
    if (debug)
    {
        std::cout << "===== calcLifeTime Debug Info =====" << std::endl;
        std::cout << "Input parameters:" << std::endl;
        std::cout << "  mass: " << mass << " GeV/c^2" << std::endl;
        std::cout << "  pvrZ: " << pvrZ << " mm" << std::endl;
        std::cout << "  decayZ: " << decayZ << " mm" << std::endl;
        std::cout << "  pz: " << pz << " GeV/c" << std::endl;
    }

    double massD0 = 1.86483;            // PDG 2020
    double sToPicoSeconds = 1e12;       // Conversion: s->ps
    double speedOfLight = 2.99792458e8; // In m/s
    double mmToM = 1e-3;                // Conversion: mm->m

    double nominalMass = massD0;

    // Check for potential division by zero
    if (std::abs(pz) < 1e-10)
    {
        if (debug)
            std::cout << "WARNING: pz is nearly zero (" << pz << "), can't calculate lifetime!" << std::endl;
        return -999.0; // Return sentinel value
    }

    // Calculate distance in flight direction
    double flightDistance = (decayZ - pvrZ) * mmToM;
    if (debug)
        std::cout << "  Flight distance: " << flightDistance << " m" << std::endl;

    // Use relativistic time dilation formula: t = d/v * γ
    // where γ = E/m = p/(m*v)
    double lifetime = (flightDistance * nominalMass / speedOfLight) / pz;

    if (debug)
    {
        std::cout << "  Used nominal mass: " << nominalMass << " GeV/c^2" << std::endl;
        std::cout << "  Flight time (s): " << lifetime << " s" << std::endl;
    }

    // Convert to picoseconds
    lifetime *= sToPicoSeconds;

    if (debug)
    {
        std::cout << "  Final lifetime: " << lifetime << " ps" << std::endl;
        std::cout << "=============================" << std::endl;
    }

    return lifetime;
}

bool isAccepted(const TLorentzVector &particleVect, int Type)
{
    double eta = particleVect.PseudoRapidity();

    if (Type == 0)
    { // General Particles (D0, K, pi)
        if (eta < 2 || eta > 4.5)
        {
            return false;
        }
    }
    else if (Type == 1)
    { // Jets
        if (eta < 2.5 || eta > 4)
        {
            return false;
        }
    }

    // D0-specific pt cut
    if (Type == 2)
    { // D0 specific
        if (particleVect.Pt() < 2.0)
        { // Minimum pT for D0
            return false;
        }
    }

    return true; // Otherwise true
}

FilterObject::FilterObject(const std::string &fFileName, float printLvl, int inputMC)
    : isMC(inputMC == 1), printLvl(printLvl), matchingVersion(1),
      pionMapReco(nullptr), pionMapReco_Var(nullptr), pionMapSelHist(nullptr),
      pionMapSelHist_Var(nullptr), MuonMap(nullptr), MuonMap_Var(nullptr),
      TriggerSelCorrMap(nullptr), TriggerSelCorrMap_Var(nullptr),
      TriggerEffMap2D(nullptr), TriggerEffMap2D_Var(nullptr)
{

    if (isMC)
    {
        std::cout << "Start filtering for MC" << std::endl;
    }
    else
    {
        std::cout << "Start filtering for Data" << std::endl;
    }

    // Set up file name
    this->fFileName = fFileName;
    if (this->fFileName.find(".root") != std::string::npos)
    {
        this->fFileName = fFileName.substr(0, fFileName.size() - 5);
    }

    std::cout << "-> Inputfilename: " << this->fFileName << std::endl;
    this->foutName = this->fFileName + "_filterV" + std::to_string(matchingVersion) + ".root";

    // Open input file
    tfileOriginal = new TFile((this->fFileName + ".root").c_str(), "READ");
    if (!tfileOriginal || tfileOriginal->IsZombie())
    {
        std::cerr << "ERROR: Cannot open input file: " << this->fFileName + ".root" << std::endl;
        exit(1);
    }

    ttreeOriginal = dynamic_cast<TTree *>(tfileOriginal->Get("data"));
    if (!ttreeOriginal)
    {
        std::cerr << "ERROR: Cannot find 'data' tree in input file" << std::endl;
        exit(1);
    }

    // Create output file
    tfile = new TFile(foutName.c_str(), "RECREATE");
    if (!tfile || tfile->IsZombie())
    {
        std::cerr << "ERROR: Cannot create output file: " << foutName << std::endl;
        exit(1);
    }

    // Create output trees
    ttree = new TTree("FragmNtuple", "");

    if (isMC)
    {
        ttreeResponse = new TTree("Response", "");
        ttreeResponse2 = new TTree("ResponseVar", "");
        mcttree = new TTree("MCFragmNtuple", "");
    }
    else
    {
        ttreeResponse = nullptr;
        ttreeResponse2 = nullptr;
        mcttree = nullptr;
    }

    // Load efficiency maps if needed for real data
    if (!isMC)
    {
        // This would be the place to load efficiency maps as in the Python code
        // I'm omitting the map loading since it's not fully implemented in the provided code
    }
}

FilterObject::~FilterObject()
{
    // Clean up files
    if (tfile)
    {
        tfile->Close();
        delete tfile;
    }

    if (tfileOriginal)
    {
        tfileOriginal->Close();
        delete tfileOriginal;
    }

    // Clean up efficiency maps and functions
    if (pionMapReco)
        delete pionMapReco;
    if (pionMapReco_Var)
        delete pionMapReco_Var;
    if (pionMapSelHist)
        delete pionMapSelHist;
    if (pionMapSelHist_Var)
        delete pionMapSelHist_Var;
    if (MuonMap)
        delete MuonMap;
    if (MuonMap_Var)
        delete MuonMap_Var;
    if (TriggerSelCorrMap)
        delete TriggerSelCorrMap;
    if (TriggerSelCorrMap_Var)
        delete TriggerSelCorrMap_Var;
    if (TriggerEffMap2D)
        delete TriggerEffMap2D;
    if (TriggerEffMap2D_Var)
        delete TriggerEffMap2D_Var;

    for (auto &func : functionList)
    {
        delete func;
    }

    for (auto &func : functionList_Var)
    {
        delete func;
    }
}

TH2D *FilterObject::createMapVariation(TH2D *inputHisto, TEfficiency *TEffObject)
{
    TH2D *mapClone = dynamic_cast<TH2D *>(inputHisto->Clone(Form("%s_Var", inputHisto->GetName())));

    for (int binX = 0; binX < inputHisto->GetNbinsX(); binX++)
    {
        for (int binY = 0; binY < inputHisto->GetNbinsY(); binY++)
        {
            double binContent = inputHisto->GetBinContent(binX, binY);
            double binError = inputHisto->GetBinError(binX, binY);

            // This is a histogram from a TEff object. The errors are nonsensical
            // and need to be determined in a separate way
            if (TEffObject)
            {
                double errEffObj = TEffObject->GetEfficiencyErrorLow(TEffObject->GetGlobalBin(binX, binY));
                mapClone->SetBinContent(binX, binY, binContent - errEffObj);
            }
            // Normal histogram
            else
            {
                mapClone->SetBinContent(binX, binY, binContent - binError);
            }
        }
    }

    return mapClone;
}

std::vector<TF1 *> FilterObject::fitTriggerEff(TH2D *triggerHist)
{
    // Get number of eta bins == fit functions
    int etaBins = triggerHist->GetNbinsY();
    std::vector<TF1 *> fitFuncArray(etaBins, nullptr);

    std::cout << "Start fitting" << std::endl;
    triggerHist->Sumw2();

    TCanvas *canvas0 = new TCanvas("canvas0", "canvas0", 0, 0, 400, 400);
    canvas0->cd();
    triggerHist->Draw("colz");
    // canvas0->SaveAs("./Eff2D.pdf");

    TCanvas *canvas = new TCanvas("canvas", "canvas", 0, 0, 600, 400);
    canvas->Divide(2, 5);

    for (int etaBin = 0; etaBin < etaBins; etaBin++)
    {
        canvas->cd(etaBin + 1);
        fitFuncArray[etaBin] = new TF1(Form("func%d", etaBin), "[0] - ([1]*exp(-1.0*[2]*x))", 1000.0, 11000.0);
        fitFuncArray[etaBin]->SetParameter(0, 0.6);
        fitFuncArray[etaBin]->SetParLimits(0, 0.3, 0.9);
        fitFuncArray[etaBin]->SetParameter(1, 0.1);
        fitFuncArray[etaBin]->SetParLimits(1, 0.001, 5.0);
        fitFuncArray[etaBin]->SetParameter(2, 0.00001);
        fitFuncArray[etaBin]->SetParLimits(2, 0.000001, 0.01);

        // Project every etaBin to the pT axis and then fit
        TH1D *proj = triggerHist->ProjectionX(Form("effHistProjNum%d", etaBin), etaBin + 1, etaBin + 1);
        proj->Draw("E");

        // Fit - set max range to 10000
        proj->Fit(fitFuncArray[etaBin], "Q", "", 0, 10000);
        fitFuncArray[etaBin]->Draw("same");
    }

    // Save the final Canvas
    // canvas->SaveAs("./EffFits.pdf");

    delete canvas;
    delete canvas0;

    return fitFuncArray;
}

double FilterObject::getD0EffCorrWeight(const TLorentzVector &kaonVector, const TLorentzVector &pionVector, int type)
{
    // Extract kinematic properties for kaon
    double p_kaon = kaonVector.P() * 1000; // P is in GeV, Maps are in MeV
    double pt_kaon = kaonVector.Pt() * 1000;
    double eta_kaon = kaonVector.PseudoRapidity();

    // Extract kinematic properties for pion
    double p_pion = pionVector.P() * 1000;
    double pt_pion = pionVector.Pt() * 1000;
    double eta_pion = pionVector.PseudoRapidity();

    // Calculate D0 combined kinematics (for trigger)
    double pt_combined = std::sqrt(pt_kaon * pt_pion);
    double eta_d0 = (eta_kaon + eta_pion) / 2.0; // Approximate eta of D0

    // Initialize efficiency
    double eff = 1.0;

    // This function would implement all the efficiency lookup logic similar to the Python version
    // For brevity, just returning 1.0 since the original code has many dependencies on efficiency maps

    // Return the inverse efficiency as the correction weight
    if (eff > 0)
    {
        return 1.0 / eff;
    }
    else
    {
        return 0;
    }
}

double FilterObject::getEffCorrWeight(const TLorentzVector &piVector1, const TLorentzVector &piVector2,
                                      const TLorentzVector &muonVect1, const TLorentzVector &muonVect2,
                                      const TLorentzVector &triggVect, int type)
{
    // This function would implement all the efficiency lookup logic similar to the Python version
    // For brevity, just returning 1.0 since the original code has many dependencies on efficiency maps

    return 1.0;
}

bool branchExists(TTree *tree, const char *branchName)
{
    if (!tree)
        return false;
    TBranch *branch = tree->GetBranch(branchName);
    return (branch != nullptr);
}

void FilterObject::filter()
{
    // Output tree variables
    // D0 tag variables
    float v_tagdR = 0;             // Distance between D0 and jet axis
    float v_tagMass = 0;           // Reconstructed D0 mass
    float v_tagPt = 0;             // D0 transverse momentum
    float v_tagEta = 0;            // D0 pseudorapidity
    float v_tag_idx_jet = 0;       // Index of associated jet
    float v_tag_lifet = 0;         // D0 lifetime
    float v_tag_lifetWrongPV = 0;  // D0 lifetime with wrong primary vertex
    float v_tagnRnd = 0;           // Number of random tags
    float v_tag_decVtxChi2 = 0;    // D0 decay vertex chi2
    float v_tag_logdecVtxChi2 = 0; // Log of decay vertex chi2

    // Jet variables
    float v_jetnTComb = 0;   // Number of combined tags in jet
    float v_jetnTUnique = 0; // Number of unique tags in jet
    float v_jetnConst = 0;   // Number of jet constituents
    float v_jetPt = 0;       // Jet transverse momentum
    float v_jetEta = 0;      // Jet pseudorapidity
    float v_tagZ = 0;        // Fragmentation fraction (pT_D0/pT_jet)
    float v_D0Z = 0;         // Secondary fragmentation fraction

    // D0 daughter variables
    float v_mpipi = 0;  // Mass of pion pair
    float v_QValue = 0; // Q-value

    // Efficiency weight variables
    float v_EffWeight = 0;   // Combined efficiency weight
    float v_EffWeight_0 = 0; // Component 0 (pion reconstruction)
    float v_EffWeight_1 = 0; // Component 1 (pion selection)
    float v_EffWeight_2 = 0; // Component 2 (muon reconstruction)
    float v_EffWeight_3 = 0; // Component 3 (trigger line correction)
    float v_EffWeight_4 = 0; // Component 4 (trigger efficiency)
    float v_EffWeight_5 = 0; // Component 5 (data/MC ratio)

    // Variation of efficiency weights
    float v_EffWeight_0_Var = 0;
    float v_EffWeight_1_Var = 0;
    float v_EffWeight_2_Var = 0;
    float v_EffWeight_3_Var = 0;
    float v_EffWeight_4_Var = 0;
    float v_EffWeight_5_Var = 0;

    // D0 pion and kaon quality variables
    float v_piPTrckChi2 = 0;  // π+ track chi2
    float v_piPprobNNpi = 0;  // π+ PID probability
    float v_piPprobGhost = 0; // π+ ghost probability
    float v_piMTrckChi2 = 0;  // π- track chi2
    float v_piMprobNNpi = 0;  // π- PID probability
    float v_piMprobGhost = 0; // π- ghost probability
    float v_decayVtxChi2 = 0; // Decay vertex chi2
    float v_Dist1 = 0;        // Distance of closest approach 1
    float v_Dist2 = 0;        // Distance of closest approach 2
    float v_Dist3 = 0;        // Distance of closest approach 3
    float v_mD0 = 0;          // Reconstructed D0 mass
    float v_nSPD = 0;         // Number of SPD hits
    float v_isPrimary2 = 0;   // Flag for primary D0
    float v_KprobNNK = 0;     // K PID probability
    float v_KprobGhost = 0;   // K ghost probability
    float v_KTrckChi2 = 0;    // K track chi2

    // Create branches for the output tree
    ttree->Branch("tagJetdR", &v_tagdR, "tagJetdR/F");
    ttree->Branch("tagMass", &v_tagMass, "tagMass/F");
    ttree->Branch("tagPt", &v_tagPt, "tagPt/F");
    ttree->Branch("tagEta", &v_tagEta, "tagEta/F");
    ttree->Branch("tagidxjet", &v_tag_idx_jet, "tagidxjet/F");
    ttree->Branch("tag_lifet", &v_tag_lifet, "tag_lifet/F");
    ttree->Branch("tag_lifetWrongPV", &v_tag_lifetWrongPV, "tag_lifetWrongPV/F");
    ttree->Branch("tagnRnd", &v_tagnRnd, "tagnRnd/F");
    ttree->Branch("tag_ip_chi2", &v_tag_decVtxChi2, "tag_ip_chi2/F");
    ttree->Branch("log_tag_ipchi2", &v_tag_logdecVtxChi2, "log_tag_ipchi2/F");

    ttree->Branch("jetnTComb", &v_jetnTComb, "jetnTComb/F");
    ttree->Branch("jetnTUniqu", &v_jetnTUnique, "jetnTUniqu/F");
    ttree->Branch("jetnConst", &v_jetnConst, "jetnConst/F");
    ttree->Branch("jetPt", &v_jetPt, "jetPt/F");
    ttree->Branch("jetEta", &v_jetEta, "jetEta/F");
    ttree->Branch("tagZ", &v_tagZ, "tagZ/F");
    ttree->Branch("D0Z", &v_D0Z, "D0Z/F");

    ttree->Branch("mpipi", &v_mpipi, "mpipi/F");
    ttree->Branch("QValue", &v_QValue, "QValue/F");

    // Analysis cut variables for efficiency correction
    ttree->Branch("piPTrckChi2", &v_piPTrckChi2, "piPTrckChi2/F");
    ttree->Branch("piPprobNNpi", &v_piPprobNNpi, "piPprobNN/F");
    ttree->Branch("piPprobGhost", &v_piPprobGhost, "piPprobGhost/F");
    ttree->Branch("piMTrckChi2", &v_piMTrckChi2, "piMTrckChi2/F");
    ttree->Branch("piMprobNNpi", &v_piMprobNNpi, "piMprobNN/F");
    ttree->Branch("piMprobGhost", &v_piMprobGhost, "piMprobGhost/F");
    ttree->Branch("decayVtxChi2", &v_decayVtxChi2, "decayVtxChi2/F");
    ttree->Branch("Distance1", &v_Dist1, "Distance1/F");
    ttree->Branch("Distance2", &v_Dist2, "Distance2/F");
    ttree->Branch("Distance3", &v_Dist3, "Distance3/F");
    ttree->Branch("mD0", &v_mD0, "mD0/F");
    ttree->Branch("nSPD", &v_nSPD, "nSPD/F");
    ttree->Branch("EffWeight", &v_EffWeight, "EffWeight/F");
    ttree->Branch("EffWeight_0", &v_EffWeight_0, "EffWeight_0/F");
    ttree->Branch("EffWeight_1", &v_EffWeight_1, "EffWeight_1/F");
    ttree->Branch("EffWeight_2", &v_EffWeight_2, "EffWeight_2/F");
    ttree->Branch("EffWeight_3", &v_EffWeight_3, "EffWeight_3/F");
    ttree->Branch("EffWeight_4", &v_EffWeight_4, "EffWeight_4/F");
    ttree->Branch("EffWeight_5", &v_EffWeight_5, "EffWeight_5/F");
    ttree->Branch("EffWeight_0_Var", &v_EffWeight_0_Var, "EffWeight_0_Var/F");
    ttree->Branch("EffWeight_1_Var", &v_EffWeight_1_Var, "EffWeight_1_Var/F");
    ttree->Branch("EffWeight_2_Var", &v_EffWeight_2_Var, "EffWeight_2_Var/F");
    ttree->Branch("EffWeight_3_Var", &v_EffWeight_3_Var, "EffWeight_3_Var/F");
    ttree->Branch("EffWeight_4_Var", &v_EffWeight_4_Var, "EffWeight_4_Var/F");
    ttree->Branch("EffWeight_5_Var", &v_EffWeight_5_Var, "EffWeight_5_Var/F");
    ttree->Branch("isPrimary", &v_isPrimary2, "isPrimary/F");
    ttree->Branch("KprobNNK", &v_KprobNNK, "KprobNNK/F");
    ttree->Branch("KprobGhost", &v_KprobGhost, "KprobGhost/F");
    ttree->Branch("KTrckChi2", &v_KTrckChi2, "KTrckChi2/F");

    // Variables for response matrix trees (MC only)
    float v_pTDet = 0;
    float v_etaDet = 0;
    float v_phiDet = 0;
    float v_etaTagDet = 0;
    float v_phiTagDet = 0;
    float v_nConstDet = 0;
    float v_TagPtDet = 0;
    float v_TagMDet = 0;
    float v_zTDet = 0;
    float v_pTPart = 0;
    float v_etaPart = 0;
    float v_phiPart = 0;
    float v_etaTagPart = 0;
    float v_phiTagPart = 0;
    float v_nConstPart = 0;
    float v_TagPtPart = 0;
    float v_TagMPart = 0;
    float v_zTPart = 0;
    float v_dR = 0;
    float v_PartnTag = 0;
    float v_DetnTag = 0;

    // MC particle level tree variables
    float v_MCetaPart = 0;
    float v_MCphiPart = 0;
    float v_MCpTPart = 0;
    float v_MCnConstPart = 0;
    float v_MCetaTagPart = 0;
    float v_MCphiTagPart = 0;
    float v_MCTagPtPart = 0;
    float v_MCTagMPart = 0;
    float v_MCzTPart = 0;
    float v_MCtag_lifetPart = 0;
    float v_isPrimary = 0;
    float v_isDetTagRec = 0;

    if (isMC)
    {
        // Set up branches for Response tree
        ttreeResponse->Branch("jetPtDet", &v_pTDet, "jetPtDet/F");
        ttreeResponse->Branch("etaDet", &v_etaDet, "etaDet/F");
        ttreeResponse->Branch("etaTagDet", &v_etaTagDet, "etaTagDet/F");
        ttreeResponse->Branch("phiDet", &v_phiDet, "phiDet/F");
        ttreeResponse->Branch("phiTagDet", &v_phiTagDet, "phiTagDet/F");
        ttreeResponse->Branch("nConstDet", &v_nConstDet, "nConstDet/F");
        ttreeResponse->Branch("tagPtDet", &v_TagPtDet, "tagPtDet/F");
        ttreeResponse->Branch("tagMDet", &v_TagMDet, "tagMDet/F");
        ttreeResponse->Branch("zTDet", &v_zTDet, "zTDet/F");
        ttreeResponse->Branch("jetPtPart", &v_pTPart, "jetPtPart/F");
        ttreeResponse->Branch("etaPart", &v_etaPart, "etaPart/F");
        ttreeResponse->Branch("etaTagPart", &v_etaTagPart, "etaTagPart/F");
        ttreeResponse->Branch("phiPart", &v_phiPart, "phiPart/F");
        ttreeResponse->Branch("phiTagPart", &v_phiTagPart, "phiTagPart/F");
        ttreeResponse->Branch("nConstPart", &v_nConstPart, "nConstPart/F");
        ttreeResponse->Branch("tagPtPart", &v_TagPtPart, "tagPtPart/F");
        ttreeResponse->Branch("zTPart", &v_zTPart, "zTPart/F");
        ttreeResponse->Branch("tagMPart", &v_TagMPart, "tagMPart/F");
        ttreeResponse->Branch("dR", &v_dR, "dR/F");
        ttreeResponse->Branch("isPrimaryPart", &v_isPrimary2, "isPrimaryPart/F");
        ttreeResponse->Branch("PartnTag", &v_PartnTag, "PartnTag/F");
        ttreeResponse->Branch("DetnTag", &v_DetnTag, "DetnTag/F");

        // Set up branches for ResponseVar tree (same structure as Response tree)
        ttreeResponse2->Branch("jetPtDet", &v_pTDet, "jetPtDet/F");
        ttreeResponse2->Branch("etaDet", &v_etaDet, "etaDet/F");
        ttreeResponse2->Branch("etaTagDet", &v_etaTagDet, "etaTagDet/F");
        ttreeResponse2->Branch("phiDet", &v_phiDet, "phiDet/F");
        ttreeResponse2->Branch("phiTagDet", &v_phiTagDet, "phiTagDet/F");
        ttreeResponse2->Branch("nConstDet", &v_nConstDet, "nConstDet/F");
        ttreeResponse2->Branch("tagPtDet", &v_TagPtDet, "tagPtDet/F");
        ttreeResponse2->Branch("tagMDet", &v_TagMDet, "tagMDet/F");
        ttreeResponse2->Branch("zTDet", &v_zTDet, "zTDet/F");
        ttreeResponse2->Branch("jetPtPart", &v_pTPart, "jetPtPart/F");
        ttreeResponse2->Branch("etaPart", &v_etaPart, "etaPart/F");
        ttreeResponse2->Branch("etaTagPart", &v_etaTagPart, "etaTagPart/F");
        ttreeResponse2->Branch("phiPart", &v_phiPart, "phiPart/F");
        ttreeResponse2->Branch("phiTagPart", &v_phiTagPart, "phiTagPart/F");
        ttreeResponse2->Branch("nConstPart", &v_nConstPart, "nConstPart/F");
        ttreeResponse2->Branch("tagPtPart", &v_TagPtPart, "tagPtPart/F");
        ttreeResponse2->Branch("zTPart", &v_zTPart, "zTPart/F");
        ttreeResponse2->Branch("tagMPart", &v_TagMPart, "tagMPart/F");
        ttreeResponse2->Branch("dR", &v_dR, "dR/F");
        ttreeResponse2->Branch("isPrimaryPart", &v_isPrimary2, "isPrimaryPart/F");
        ttreeResponse2->Branch("PartnTag", &v_PartnTag, "PartnTag/F");
        ttreeResponse2->Branch("DetnTag", &v_DetnTag, "DetnTag/F");

        // Set up branches for MC truth tree
        mcttree->Branch("jetPtPart", &v_MCpTPart, "jetPtPart/F");
        mcttree->Branch("tagPtPart", &v_MCTagPtPart, "tagPtPart/F");
        mcttree->Branch("etaJetPart", &v_MCetaPart, "etaJetPart/F");
        mcttree->Branch("etaTagPart", &v_MCetaTagPart, "etaTagPart/F");
        mcttree->Branch("phiJetPart", &v_MCphiPart, "phiJetPart/F");
        mcttree->Branch("phiTagPart", &v_MCphiTagPart, "phiTagPart/F");
        mcttree->Branch("nConstPart", &v_MCnConstPart, "nConstPart/F");
        mcttree->Branch("zTPart", &v_MCzTPart, "zTPart/F");
        mcttree->Branch("tag_lifetPart", &v_MCtag_lifetPart, "tag_lifetPart/F");
        mcttree->Branch("isPrimary", &v_isPrimary, "isPrimary/F");
        mcttree->Branch("isDetTagRec", &v_isDetTagRec, "isDetTagRec/F");
        mcttree->Branch("tagMPart", &v_MCTagMPart, "tagMPart/F");
    }

    // Vector to store good detector tag indices
    std::vector<int> v_goodDetTags_idx;

    // To load info from input tree
    int tag_pid = 0;

    // Check if file and tree are valid
    if (!tfileOriginal || tfileOriginal->IsZombie())
    {
        std::cerr << "ERROR: Input file is invalid or could not be opened." << std::endl;
        return;
    }

    if (!ttreeOriginal)
    {
        std::cerr << "ERROR: Input tree not found." << std::endl;
        return;
    }

    ttreeOriginal->SetBranchStatus("*", 1);        // Enable all branches
    ttreeOriginal->SetCacheSize(50 * 1024 * 1024); // Set cache size to 50 MB
    bool firstEvent = true;

    // Event loop
    int totalNbOfEntries = ttreeOriginal->GetEntries();
    std::cout << "oo Start loop over ntuple to extract important quantities" << std::endl;
    std::cout << "oo Tree has " << totalNbOfEntries << " entries" << std::endl;

    if (totalNbOfEntries == 0)
        return;

    // Determine frequency of printout (should update after about 10% of events)
    int printFreq = totalNbOfEntries / printLvl;
    if (printFreq < 1)
    {
        printFreq = 1;
    }
    // Set branch addresses with proper checks
    try
    {
        // Required branches - essential for processing
        if (!branchExists(ttreeOriginal, "evt_n") ||
            !branchExists(ttreeOriginal, "tag_pid"))
        {
            std::cerr << "ERROR: Required branches missing from input file" << std::endl;
            return;
        }

        ttreeOriginal->SetBranchAddress("evt_n", &evt_n);
        ttreeOriginal->SetBranchAddress("tag_pid", &tag_pid_vec);

// Optional branches - set with error checking
// Standard branches with error handling
#define SAFE_SET_BRANCH(name, ptr)                                                              \
    if (branchExists(ttreeOriginal, name))                                                      \
    {                                                                                           \
        ttreeOriginal->SetBranchAddress(name, ptr);                                             \
    }                                                                                           \
    else                                                                                        \
    {                                                                                           \
        std::cout << "Branch '" << name << "' not found, will use default values" << std::endl; \
    }

        // Event info
        SAFE_SET_BRANCH("evt_spd", &evt_spd);
        SAFE_SET_BRANCH("pvr_n", &pvrNumber);

        // Jet properties
        SAFE_SET_BRANCH("jet_px", &jet_px);
        SAFE_SET_BRANCH("jet_py", &jet_py);
        SAFE_SET_BRANCH("jet_pz", &jet_pz);
        SAFE_SET_BRANCH("jet_e", &jet_e);
        SAFE_SET_BRANCH("jet_n_tagsComb", &jet_n_tagsComb);
        SAFE_SET_BRANCH("jet_n_tagsUnique", &jet_n_tagsUnique);
        SAFE_SET_BRANCH("jet_n_neu", &jet_n_neu);
        SAFE_SET_BRANCH("jet_n_chr", &jet_n_chr);

        // Tag properties
        SAFE_SET_BRANCH("tag_l0_tos0", &tag_l0_tos0);
        SAFE_SET_BRANCH("tag_hlt1_tos0", &tag_hlt1_tos0);
        SAFE_SET_BRANCH("tag_idx_jet", &tag_idx_jet);
        SAFE_SET_BRANCH("tag_idx_pvr", &tag_idx_pvr);
        SAFE_SET_BRANCH("tag_idx_prt0", &tag_idx_prt0);
        SAFE_SET_BRANCH("tag_idx_prt1", &tag_idx_prt1);
        SAFE_SET_BRANCH("tag_dR_jet", &tag_dR_jet);
        SAFE_SET_BRANCH("pvr_z", &pvr_z);
        SAFE_SET_BRANCH("tag_z", &tag_z);
        SAFE_SET_BRANCH("tag_pz", &tag_pz);

        // Tag kinematic properties
        SAFE_SET_BRANCH("tag_px", &tag_px);
        SAFE_SET_BRANCH("tag_py", &tag_py);
        SAFE_SET_BRANCH("tag_pz", &tag_pz);
        SAFE_SET_BRANCH("tag_e", &tag_e);

        // Particle properties
        SAFE_SET_BRANCH("prt_px", &prt_px);
        SAFE_SET_BRANCH("prt_py", &prt_py);
        SAFE_SET_BRANCH("prt_pz", &prt_pz);
        SAFE_SET_BRANCH("prt_e", &prt_e);
        SAFE_SET_BRANCH("prt_pid", &prt_pid);

        // Optional quality branches - these were causing the crash
        if (branchExists(ttreeOriginal, "prt_pnn_k"))
        {
            // Initialize these pointers only if needed
            if (!prt_pnn_k)
                prt_pnn_k = new std::vector<float>();
            ttreeOriginal->SetBranchAddress("prt_pnn_k", &prt_pnn_k);
        }

        if (branchExists(ttreeOriginal, "prt_prb_ghost"))
        {
            if (!prt_prb_ghost)
                prt_prb_ghost = new std::vector<float>();
            ttreeOriginal->SetBranchAddress("prt_prb_ghost", &prt_prb_ghost);
        }

        if (branchExists(ttreeOriginal, "prt_trckChi2"))
        {
            if (!prt_trckChi2)
                prt_trckChi2 = new std::vector<float>();
            ttreeOriginal->SetBranchAddress("prt_trckChi2", &prt_trckChi2);
        }

        if (branchExists(ttreeOriginal, "prt_pnn_pi"))
        {
            if (!prt_pnn_pi)
                prt_pnn_pi = new std::vector<float>();
            ttreeOriginal->SetBranchAddress("prt_pnn_pi", &prt_pnn_pi);
        }

        // Optional tag branches
        SAFE_SET_BRANCH("tag_n_tagsRnd", &tag_n_tagsRnd);
        SAFE_SET_BRANCH("tag_decVtxChi2", &tag_decVtxChi2);

        // MC branches if needed
        if (isMC)
        {
            // MC stuff only if needed
            SAFE_SET_BRANCH("mctag_pid", &mctag_pid);
            SAFE_SET_BRANCH("tag_genAssocType", &tag_genAssocType);
            SAFE_SET_BRANCH("tag_idx_MCtag", &tag_idx_MCtag);
            SAFE_SET_BRANCH("mctag_isPrimary", &mctag_isPrimary);
            SAFE_SET_BRANCH("mctag_idx_jet", &mctag_idx_jet);
            SAFE_SET_BRANCH("mctag_idx_Dettag", &mctag_idx_Dettag);
            SAFE_SET_BRANCH("mcjet_n_neu", &mcjet_n_neu);
            SAFE_SET_BRANCH("mcjet_n_chr", &mcjet_n_chr);
            SAFE_SET_BRANCH("mctag_idx_prt1", &mctag_idx_prt1);
            SAFE_SET_BRANCH("mctag_idx_prt2", &mctag_idx_prt2);
            SAFE_SET_BRANCH("mcprt_pid", &mcprt_pid);
            SAFE_SET_BRANCH("mctag_idx_pvr", &mctag_idx_pvr);

            // MC kinematic properties
            SAFE_SET_BRANCH("mcjet_px", &mcjet_px);
            SAFE_SET_BRANCH("mcjet_py", &mcjet_py);
            SAFE_SET_BRANCH("mcjet_pz", &mcjet_pz);
            SAFE_SET_BRANCH("mcjet_e", &mcjet_e);

            SAFE_SET_BRANCH("mctag_px", &mctag_px);
            SAFE_SET_BRANCH("mctag_py", &mctag_py);
            SAFE_SET_BRANCH("mctag_pz", &mctag_pz);
            SAFE_SET_BRANCH("mctag_z", &mctag_z);
            SAFE_SET_BRANCH("mctag_e", &mctag_e);
            SAFE_SET_BRANCH("mcpvr_z", &mcpvr_z);

            SAFE_SET_BRANCH("mcprt_px", &mcprt_px);
            SAFE_SET_BRANCH("mcprt_py", &mcprt_py);
            SAFE_SET_BRANCH("mcprt_pz", &mcprt_pz);
            SAFE_SET_BRANCH("mcprt_e", &mcprt_e);
        }

        // Try reading the first entry to make sure branches are readable
        if (ttreeOriginal->GetEntry(0) <= 0)
        {
            std::cerr << "ERROR: Failed to read first entry from tree" << std::endl;
            return;
        }
    }
    catch (std::bad_alloc &ba)
    {
        std::cerr << "ERROR: Memory allocation failed: " << ba.what() << std::endl;
        return;
    }
    catch (const std::exception &e)
    {
        std::cerr << "ERROR: Exception during branch setup: " << e.what() << std::endl;
        return;
    }
    catch (...)
    {
        std::cerr << "ERROR: Unknown exception during branch setup" << std::endl;
        return;
    }
    bool verboseOut = true;

    // Event loop
    for (int iEvt = 0; iEvt < totalNbOfEntries; iEvt++)
    {
        // std::cout << "Processing event " << iEvt + 1 << " of " << totalNbOfEntries << std::endl;
        if (firstEvent)
        {
            firstEvent = false;
        }
        // For subsequent events, we need to properly reset all vector pointers
        else
        {
            // Reset all branch pointers to avoid the segfault
            tag_pid_vec = nullptr;
            jet_px = nullptr;
            mcjet_px = nullptr;
            tag_l0_tos0 = nullptr;
            tag_hlt1_tos0 = nullptr;
            tag_idx_jet = nullptr;
            tag_genAssocType = nullptr;
            tag_idx_pvr = nullptr;
            tag_idx_prt0 = nullptr;
            tag_idx_prt1 = nullptr;
            tag_idx_MCtag = nullptr;
            mctag_isPrimary = nullptr;
            prt_pid = nullptr;
            tag_dR_jet = nullptr;
            pvr_z = nullptr;
            tag_z = nullptr;
            tag_pz = nullptr;
            jet_n_tagsComb = nullptr;
            jet_n_tagsUnique = nullptr;
            jet_n_neu = nullptr;
            jet_n_chr = nullptr;
            prt_pnn_k = nullptr;
            prt_prb_ghost = nullptr;
            prt_trckChi2 = nullptr;
            prt_pnn_pi = nullptr;
            mctag_idx_jet = nullptr;
            mctag_idx_Dettag = nullptr;
            mcjet_n_neu = nullptr;
            mcjet_n_chr = nullptr;
            mctag_idx_prt1 = nullptr;
            mctag_idx_prt2 = nullptr;
            mctag_idx_pvr = nullptr;
            mcprt_pid = nullptr;

            // Reset all branch addresses
            ttreeOriginal->ResetBranchAddresses();

            // Set up branch addresses again
            ttreeOriginal->SetBranchAddress("evt_n", &evt_n);
            ttreeOriginal->SetBranchAddress("evt_spd", &evt_spd);
            ttreeOriginal->SetBranchAddress("tag_pid", &tag_pid_vec);
            ttreeOriginal->SetBranchAddress("tag_idx_jet", &tag_idx_jet);
            ttreeOriginal->SetBranchAddress("pvr_n", &pvrNumber);
            ttreeOriginal->SetBranchAddress("tag_l0_tos0", &tag_l0_tos0);
            ttreeOriginal->SetBranchAddress("tag_hlt1_tos0", &tag_hlt1_tos0);
            ttreeOriginal->SetBranchAddress("tag_idx_pvr", &tag_idx_pvr);
            ttreeOriginal->SetBranchAddress("tag_idx_prt0", &tag_idx_prt0);
            ttreeOriginal->SetBranchAddress("tag_idx_prt1", &tag_idx_prt1);
            ttreeOriginal->SetBranchAddress("prt_pid", &prt_pid);
            ttreeOriginal->SetBranchAddress("tag_dR_jet", &tag_dR_jet);
            ttreeOriginal->SetBranchAddress("pvr_z", &pvr_z);
            ttreeOriginal->SetBranchAddress("tag_z", &tag_z);
            ttreeOriginal->SetBranchAddress("jet_n_tagsComb", &jet_n_tagsComb);
            ttreeOriginal->SetBranchAddress("jet_n_tagsUnique", &jet_n_tagsUnique);
            ttreeOriginal->SetBranchAddress("jet_n_neu", &jet_n_neu);
            ttreeOriginal->SetBranchAddress("jet_n_chr", &jet_n_chr);

            ttreeOriginal->SetBranchAddress("jet_px", &jet_px);
            ttreeOriginal->SetBranchAddress("jet_py", &jet_py);
            ttreeOriginal->SetBranchAddress("jet_pz", &jet_pz);
            ttreeOriginal->SetBranchAddress("jet_e", &jet_e);

            ttreeOriginal->SetBranchAddress("tag_px", &tag_px);
            ttreeOriginal->SetBranchAddress("tag_py", &tag_py);
            ttreeOriginal->SetBranchAddress("tag_pz", &tag_pz);
            ttreeOriginal->SetBranchAddress("tag_e", &tag_e);

            ttreeOriginal->SetBranchAddress("prt_px", &prt_px);
            ttreeOriginal->SetBranchAddress("prt_py", &prt_py);
            ttreeOriginal->SetBranchAddress("prt_pz", &prt_pz);
            ttreeOriginal->SetBranchAddress("prt_e", &prt_e);
            // Optional branches with try-catch blocks
            ttreeOriginal->SetBranchAddress("tag_n_tagsRnd", &tag_n_tagsRnd);
            ttreeOriginal->SetBranchAddress("tag_decVtxChi2", &tag_decVtxChi2);
            ttreeOriginal->SetBranchAddress("prt_pnn_k", &prt_pnn_k);
            ttreeOriginal->SetBranchAddress("prt_prb_ghost", &prt_prb_ghost);
            ttreeOriginal->SetBranchAddress("prt_trckChi2", &prt_trckChi2);
            ttreeOriginal->SetBranchAddress("prt_pnn_pi", &prt_pnn_pi);
            // Conditional MC branches
            if (isMC)
            {
                ttreeOriginal->SetBranchAddress("mctag_pid", &mctag_pid);
                ttreeOriginal->SetBranchAddress("tag_genAssocType", &tag_genAssocType);
                ttreeOriginal->SetBranchAddress("tag_idx_MCtag", &tag_idx_MCtag);
                ttreeOriginal->SetBranchAddress("mctag_isPrimary", &mctag_isPrimary);
                ttreeOriginal->SetBranchAddress("mctag_idx_jet", &mctag_idx_jet);
                ttreeOriginal->SetBranchAddress("mctag_idx_Dettag", &mctag_idx_Dettag);
                ttreeOriginal->SetBranchAddress("mcjet_n_neu", &mcjet_n_neu);
                ttreeOriginal->SetBranchAddress("mcjet_n_chr", &mcjet_n_chr);
                ttreeOriginal->SetBranchAddress("mctag_idx_prt1", &mctag_idx_prt1);
                ttreeOriginal->SetBranchAddress("mctag_idx_prt2", &mctag_idx_prt2);
                ttreeOriginal->SetBranchAddress("mctag_idx_pvr", &mctag_idx_pvr);
                ttreeOriginal->SetBranchAddress("mcprt_pid", &mcprt_pid);
                ttreeOriginal->SetBranchAddress("mcjet_px", &mcjet_px);
                ttreeOriginal->SetBranchAddress("mcjet_py", &mcjet_py);
                ttreeOriginal->SetBranchAddress("mcjet_pz", &mcjet_pz);
                ttreeOriginal->SetBranchAddress("mcjet_e", &mcjet_e);

                ttreeOriginal->SetBranchAddress("mctag_px", &mctag_px);
                ttreeOriginal->SetBranchAddress("mctag_py", &mctag_py);
                ttreeOriginal->SetBranchAddress("mctag_pz", &mctag_pz);
                ttreeOriginal->SetBranchAddress("mctag_z", &mctag_z);
                ttreeOriginal->SetBranchAddress("mctag_e", &mctag_e);
                ttreeOriginal->SetBranchAddress("mcpvr_z", &mcpvr_z);

                ttreeOriginal->SetBranchAddress("mcprt_px", &mcprt_px);
                ttreeOriginal->SetBranchAddress("mcprt_py", &mcprt_py);
                ttreeOriginal->SetBranchAddress("mcprt_pz", &mcprt_pz);
                ttreeOriginal->SetBranchAddress("mcprt_e", &mcprt_e);
            }
        }

        // std::cout << "Processing event: " << iEvt << std::endl;

        ttreeOriginal->GetEntry(iEvt);

        if (iEvt % printFreq == 0)
        {
            std::cout << " Events processed: " << iEvt << " (" << std::round(100.0 * iEvt / totalNbOfEntries) << "%)" << std::endl;
        }

        v_nSPD = evt_spd;
        int nTagDetLvl = tag_pid_vec->size();
        int nJetsDetLvl = jet_px->size();

        // Continue with the rest of the processing...
        int nTagPartLvl = 0;
        int nJetsPartLvl = 0;

        if (isMC)
        {
            nTagPartLvl = mctag_pid->size();

            ttreeOriginal->SetBranchAddress("mcjet_px", &mcjet_px);
            nJetsPartLvl = mcjet_px->size();
        }

        if (isMC && nTagPartLvl > 0)
        {
            nTagPartLvlEvt++;
        }

        if (nTagDetLvl > 0)
        {
            nTagDetLvlEvt++;
        }

        // Step 1 - Check for good events
        if (isMC && nTagPartLvl == 0)
        {
            debugContinue(__LINE__, "No particle level tags in MC");
            continue;
        }
        // Fill MC particle level tree if needed
        if (isMC)
        {
            TTree *mcTreeTemp = fillPartLvlTree(mcttree, ttreeOriginal, v_MCetaPart, v_MCphiPart, v_MCpTPart,
                                                v_MCnConstPart, v_MCetaTagPart, v_MCphiTagPart, v_MCTagPtPart,
                                                v_MCTagMPart, v_MCzTPart, v_MCtag_lifetPart, v_isPrimary, v_isDetTagRec);
        }

        if (nTagDetLvl == 0)
        {
            debugContinue(__LINE__, "No detector level tags");
            continue;
        }

        int tag_idxJetSize = tag_idx_jet->size();
        int tag_pidSize = tag_pid_vec->size();

        if (pvrNumber != 1)
        {
            debugContinue(__LINE__, "PVR number not equal to 1");
            continue;
        }
        if (tag_idxJetSize != tag_pidSize)
        {
            debugContinue(__LINE__, "Tag index jet size does not match tag PID size");
            continue;
        }

        // Step 2 - Collect all the good detector level tags in the event
        v_goodDetTags_idx.clear();

        for (int iTagDetLvl = 0; iTagDetLvl < nTagDetLvl; iTagDetLvl++)
        {
            int tag_l0_tos0_val = (*tag_l0_tos0)[iTagDetLvl];

            int tag_hlt1_tos0_val = (*tag_hlt1_tos0)[iTagDetLvl];

            int jetNumberTurbo = (*tag_idx_jet)[iTagDetLvl]; // JetID the tag is associated with

            int tag_pid_val = (*tag_pid_vec)[iTagDetLvl];

            int tag_assocType = 0;
            if (isMC)
            {
                tag_assocType = (*tag_genAssocType)[iTagDetLvl];
            }

            // Tag filtering
            if (jetNumberTurbo == -1)
            {
                debugContinue(__LINE__, "Invalid jet number (jetNumberTurbo == -1)");
                continue;
            }

            if (isMC && tag_assocType > 2)
            {
                debugContinue(__LINE__, Form("Tag assoc. type %d is larger than 2", tag_assocType));
                continue;
            }

            v_goodDetTags_idx.push_back(iTagDetLvl);
        }

        if (v_goodDetTags_idx.empty())
        {
            debugContinue(__LINE__, "No good detector tags found");
            continue;
        }

        // Step 3 - Get all info for the Tag inside Jet analysis
        for (int iTagDetLvl : v_goodDetTags_idx)
        {
            int jetNumberTurbo = (*tag_idx_jet)[iTagDetLvl]; // Jet ID the tag is associated with
            int pvrNumberTurbo = (*tag_idx_pvr)[iTagDetLvl];
            int prt0Number = (*tag_idx_prt0)[iTagDetLvl];
            int prt1Number = (*tag_idx_prt1)[iTagDetLvl];

            // Check if MC linked particle is a primary particle
            if (isMC)
            {
                int partLvlTagNumber = (*tag_idx_MCtag)[iTagDetLvl]; // Part. lvl tag ID associated with this tag
                v_isPrimary2 = (*mctag_isPrimary)[partLvlTagNumber];
            }
            else
            {
                v_isPrimary2 = -1;
            }

            std::vector<int> partList = {prt0Number, prt1Number};

            // Event observables
            int pid0 = (*prt_pid)[prt0Number];
            int pid1 = (*prt_pid)[prt1Number];

            // Tag observables
            float dRTagTurbo = (*tag_dR_jet)[iTagDetLvl];
            int nTRnd = -1;

            try
            {
                nTRnd = (*tag_n_tagsRnd)[iTagDetLvl];
            }
            catch (...)
            {
                nTRnd = -1;
            }

            TLorentzVector p4TagTurbo = getFourVector("tag", iTagDetLvl);
            if (!isAccepted(p4TagTurbo, 0))
            {
                debugContinue(__LINE__, "Tag not accepted by isAccepted()");
                continue;
            }

            v_Dist1 = -1;
            v_Dist2 = -1;
            v_Dist3 = -1;

            v_tagdR = dRTagTurbo;
            v_tagMass = p4TagTurbo.M();
            v_tagPt = p4TagTurbo.Pt();
            v_tagEta = p4TagTurbo.PseudoRapidity();
            v_tag_idx_jet = jetNumberTurbo;
            v_tagnRnd = nTRnd;

            try
            {
                v_decayVtxChi2 = (*tag_decVtxChi2)[iTagDetLvl];
                v_tag_decVtxChi2 = (*tag_decVtxChi2)[iTagDetLvl];

                // Set the logarithm of the decay vertex chi2
                v_tag_logdecVtxChi2 = std::log10(v_tag_decVtxChi2);
            }
            catch (...)
            {
                v_decayVtxChi2 = -1;
                v_tag_decVtxChi2 = -1;
                v_tag_logdecVtxChi2 = -1;
            }

            float pvrZ = (*pvr_z)[pvrNumberTurbo]; // Primary vertex of particle
            float decayZ = (*tag_z)[iTagDetLvl];   // Decay vertex of particle
            float pz = (*tag_pz)[iTagDetLvl];

            double lifetime = calcLifeTime(p4TagTurbo.M(), pvrZ, decayZ, pz);
            double lifetimeWrong = -50;

            if (lastPV != 0)
            {
                double pvrZWrong = lastPV;
                lifetimeWrong = calcLifeTime(p4TagTurbo.M(), pvrZWrong, decayZ, pz);
            }

            lastPV = pvrZ; // Set current PV as wrong PV for next loop

            v_tag_lifet = lifetime;
            v_tag_lifetWrongPV = lifetimeWrong;

            // Jet observables
            TLorentzVector p4JetTurbo; // Default constructor initializes to (0,0,0,0)

            if (jetNumberTurbo < 0 || jetNumberTurbo >= static_cast<int>(jet_px->size()))
            {
                debugContinue(__LINE__, "Invalid jet number: " + std::to_string(jetNumberTurbo));
                continue;
            }

            // Now it's safe to get the jet four-vector
            try
            {
                p4JetTurbo = getFourVector("jet", jetNumberTurbo);
                // std::cout << "Jet number: " << jetNumberTurbo << std::endl;
                // std::cout << " with pT: " << p4JetTurbo.Pt() << std::endl;

                if (!isAccepted(p4JetTurbo, 1))
                {
                    debugContinue(__LINE__, "Jet not accepted by isAccepted()");
                    continue;
                }
            }
            catch (const std::exception &e)
            {
                debugContinue(__LINE__, "Error when processing jet: " + std::string(e.what()));
                continue;
            }
            catch (...)
            {
                debugContinue(__LINE__, "Unknown error when processing jet");
                continue;
            }

            int nTComb = (*jet_n_tagsComb)[jetNumberTurbo];
            int nTUnique = (*jet_n_tagsUnique)[jetNumberTurbo];
            int nNeu = (*jet_n_neu)[jetNumberTurbo];
            int nCh = (*jet_n_chr)[jetNumberTurbo];

            // Also for jets: nConstituents in jet
            v_jetnConst = nNeu + nCh;
            v_jetnTComb = nTComb;
            v_jetnTUnique = nTUnique;
            v_jetPt = p4JetTurbo.Pt();
            v_jetEta = p4JetTurbo.PseudoRapidity();
            v_tagZ = p4TagTurbo.Pt() / p4JetTurbo.Pt();

            // Tag daughter observables
            // In case you want to know any daughter particle information
            bool foundD0 = true;
            bool foundK = false;
            bool foundPiP = false;

            for (int partNr : partList)
            {
                // Look for D0 meson (421) instead of J/Psi (443)
                if ((*prt_pid)[partNr] == 421)
                {
                    TLorentzVector p4D0Turbo = getFourVector("prt", partNr);
                    if (isAccepted(p4D0Turbo))
                    {
                        foundD0 = true;
                    }
                    v_D0Z = p4D0Turbo.Pt() / p4JetTurbo.Pt();
                }
                // Look for K- (-321) instead of mu-
                else if ((*prt_pid)[partNr] == -321)
                {
                    TLorentzVector p4KTurbo = getFourVector("prt", partNr);
                    try
                    {
                        v_KprobNNK = (*prt_pnn_k)[partNr];
                        v_KprobGhost = (*prt_prb_ghost)[partNr];
                        v_KTrckChi2 = (*prt_trckChi2)[partNr];
                    }
                    catch (...)
                    {
                        v_KprobNNK = -10;
                        v_KprobGhost = 10;
                        v_KTrckChi2 = -1;
                    }

                    if (isAccepted(p4KTurbo))
                    {
                        foundK = true;
                    }
                }
                // Keep pion (211) handling
                else if ((*prt_pid)[partNr] == 211)
                {
                    TLorentzVector p4PipTurbo = getFourVector("prt", partNr);
                    try
                    {
                        v_piPprobNNpi = (*prt_pnn_pi)[partNr];
                        v_piPprobGhost = (*prt_prb_ghost)[partNr];
                        v_piPTrckChi2 = (*prt_trckChi2)[partNr];
                    }
                    catch (...)
                    {
                        v_piPprobNNpi = -10; // Cut on > 0.9
                        v_piPprobGhost = 10; // Cut on < 0.3
                        v_piPTrckChi2 = -1;
                    }

                    if (isAccepted(p4PipTurbo))
                    {
                        foundPiP = true;
                    }
                }
            }
            // Check if all particles are found
            if (foundD0 && foundK && foundPiP)
            {
                // All particles found, fill the tree
                // std::cout << "All particles found, filling tree..." << std::endl;
                ttree->Fill();
            }
        }

        // Process MC matching for response matrix if applicable
        if (isMC)
        {
            // Step 4 - Build a matrix of Part Lvl Jet-Tag relations
            // Create detector level matrix
            std::vector<std::vector<int>> m_TagJetDetLvl(nTagDetLvl, std::vector<int>(nJetsDetLvl, 0));

            for (int DetTag_idx : v_goodDetTags_idx)
            {
                int jetNumber_idx = (*tag_idx_jet)[DetTag_idx]; // Jet ID the tag is associated with
                if (jetNumber_idx > -1)
                {
                    m_TagJetDetLvl[DetTag_idx][jetNumber_idx] = 1;
                }
            }

            // Create generator level matrix
            std::vector<std::vector<int>> m_TagJetPartLvl(nTagPartLvl, std::vector<int>(nJetsPartLvl, 0));

            for (int PartTag_idx = 0; PartTag_idx < nTagPartLvl; PartTag_idx++)
            {
                int jetNumber_idx = (*mctag_idx_jet)[PartTag_idx]; // Jet ID the tag is associated with
                if (jetNumber_idx > -1)
                {
                    m_TagJetPartLvl[PartTag_idx][jetNumber_idx] = 1;
                }
            }

            // Step 5 - Build a matrix of Jet-Jet distance relations
            if (nJetsPartLvl > 0 && nJetsDetLvl > 0)
            {
                // For matching versions 1 or 3
                if (matchingVersion == 1 || matchingVersion == 3)
                {
                    // Create distance matrix using Eigen
                    std::vector<std::vector<double>> m_JetJetDistance(nJetsPartLvl, std::vector<double>(nJetsDetLvl, 10.0));

                    for (int partJetidx = 0; partJetidx < nJetsPartLvl; partJetidx++)
                    {
                        bool containsPartLvlTag = false;

                        // Check if this jet contains a tag
                        for (int tagIdx = 0; tagIdx < nTagPartLvl; tagIdx++)
                        {
                            if (m_TagJetPartLvl[tagIdx][partJetidx] == 1)
                            {
                                containsPartLvlTag = true;
                                break;
                            }
                        }

                        if (matchingVersion == 1 && !containsPartLvlTag)
                        {
                            // debugContinue(__LINE__, "Continued due to no particle level tag in jet");
                            continue;
                        }

                        TLorentzVector p4JetPart = getFourVector("mcjet", partJetidx);

                        for (int detJetidx = 0; detJetidx < nJetsDetLvl; detJetidx++)
                        {
                            bool containsDetLvlTag = false;

                            // Check if this jet contains a tag
                            for (int tagIdx = 0; tagIdx < nTagDetLvl; tagIdx++)
                            {
                                if (m_TagJetDetLvl[tagIdx][detJetidx] == 1)
                                {
                                    containsDetLvlTag = true;
                                    break;
                                }
                            }

                            if (matchingVersion == 1 && !containsDetLvlTag)
                            {
                                // debugContinue(__LINE__, "Continued due to no detector level tag in jet");
                                continue;
                            }

                            TLorentzVector p4JetDet = getFourVector("jet", detJetidx);
                            double distance = p4JetPart.DeltaR(p4JetDet);
                            m_JetJetDistance[partJetidx][detJetidx] = distance;
                        }
                    }

                    // Step 6 - Fill Response matrix
                    for (int partJetidx = 0; partJetidx < nJetsPartLvl; partJetidx++)
                    {
                        int closestDetJetidx = -1;

                        if (matchingVersion == 1 || matchingVersion == 3)
                        {
                            closestDetJetidx = findClosest(m_JetJetDistance, partJetidx, 1);

                            if (closestDetJetidx == -1)
                            {
                                debugContinue(__LINE__, "Continued due to no closest detector jet found");
                                continue;
                            }

                            TLorentzVector p4PartJet = getFourVector("mcjet", partJetidx);
                            TLorentzVector p4DetJet = getFourVector("jet", closestDetJetidx);

                            std::vector<int> ListWithPartTags(nTagPartLvl, 0);
                            std::vector<int> ListWithDetTags(nTagDetLvl, 0);

                            // Find tags in particle level jet
                            for (int tagIdx = 0; tagIdx < nTagPartLvl; tagIdx++)
                            {
                                if (m_TagJetPartLvl[tagIdx][partJetidx] == 1)
                                {
                                    ListWithPartTags[tagIdx] = 1;
                                }
                            }

                            // Find tags in detector level jet
                            for (int tagIdx = 0; tagIdx < nTagDetLvl; tagIdx++)
                            {
                                if (m_TagJetDetLvl[tagIdx][closestDetJetidx] == 1)
                                {
                                    ListWithDetTags[tagIdx] = 1;
                                }
                            }

                            // Get particle level tag index
                            int idxTagPart = -1;
                            for (size_t i = 0; i < ListWithPartTags.size(); i++)
                            {
                                if (ListWithPartTags[i] == 1)
                                {
                                    idxTagPart = i;
                                    break;
                                }
                            }

                            if (idxTagPart == -1)
                            {
                                // debugContinue(__LINE__, "Continued due to no particle level tag found");
                                continue;
                            }

                            // Check if generator level tag fulfills fiducial requirements
                            bool particlesInAcceptance = false;
                            try
                            {
                                // Get the kaon and pion from D0
                                int prt1Number = (*mctag_idx_prt1)[idxTagPart];
                                int prt2Number = (*mctag_idx_prt2)[idxTagPart];

                                int pid1 = (*mcprt_pid)[prt1Number];
                                int pid2 = (*mcprt_pid)[prt2Number];

                                // Check if we have K(-321) and pi(211)
                                if ((pid1 == -321 && pid2 == 211) || (pid1 == 211 && pid2 == -321))
                                {
                                    TLorentzVector p4MCkaon = getFourVector("mcprt",
                                                                            (pid1 == -321) ? prt1Number : prt2Number);
                                    TLorentzVector p4MCpion = getFourVector("mcprt",
                                                                            (pid1 == 211) ? prt1Number : prt2Number);

                                    // Appropriate D0 daughter kinematic cuts
                                    if (p4MCkaon.Pt() > 0.25 && p4MCpion.Pt() > 0.25 &&
                                        p4MCkaon.P() > 2.0 && p4MCpion.P() > 2.0)
                                    {
                                        particlesInAcceptance = true;
                                    }
                                }
                            }
                            catch (...)
                            {
                                particlesInAcceptance = false;
                            }

                            // Skip if particles not in acceptance
                            if (!particlesInAcceptance)
                            {
                                debugContinue(__LINE__, "Continued due to particles not in acceptance");
                                continue;
                            }

                            // Get detector level tag index
                            int idxTagDet = -1;
                            for (size_t i = 0; i < ListWithDetTags.size(); i++)
                            {
                                if (ListWithDetTags[i] == 1)
                                {
                                    idxTagDet = i;
                                    break;
                                }
                            }

                            if (idxTagDet == -1)
                            {
                                debugContinue(__LINE__, "Continued due to no detector level tag found");
                                continue;
                            }

                            TLorentzVector p4PartTag = getFourVector("mctag", idxTagPart);
                            TLorentzVector p4DetTag = getFourVector("tag", idxTagDet);

                            // Fill variables for response tree
                            v_pTDet = p4DetJet.Pt();
                            v_etaDet = p4DetJet.Eta();
                            v_phiDet = p4DetJet.Phi();
                            v_nConstDet = (*jet_n_neu)[closestDetJetidx] + (*jet_n_chr)[closestDetJetidx];
                            v_TagPtDet = p4DetTag.Pt();
                            v_etaTagDet = p4DetTag.Eta();
                            v_phiTagDet = p4DetTag.Phi();
                            v_TagMDet = p4DetTag.M();
                            v_zTDet = p4DetTag.Pt() / p4DetJet.Pt();
                            v_pTPart = p4PartJet.Pt();
                            v_etaPart = p4PartJet.Eta();
                            v_phiPart = p4PartJet.Phi();

                            v_nConstPart = (*mcjet_n_neu)[partJetidx] + (*mcjet_n_chr)[partJetidx];
                            v_TagPtPart = p4PartTag.Pt();
                            v_etaTagPart = p4PartTag.Eta();
                            v_phiTagPart = p4PartTag.Phi();
                            v_TagMPart = p4PartTag.M();
                            v_zTPart = p4PartTag.Pt() / p4PartJet.Pt();
                            v_dR = m_JetJetDistance[partJetidx][closestDetJetidx];

                            v_isPrimary2 = (*mctag_isPrimary)[idxTagPart]; // Is the mcTag a primary particle
                            v_PartnTag = 0;
                            v_DetnTag = 0;

                            // Count number of tags in each jet
                            for (size_t i = 0; i < ListWithPartTags.size(); i++)
                            {
                                if (ListWithPartTags[i] == 1)
                                    v_PartnTag++;
                            }

                            for (size_t i = 0; i < ListWithDetTags.size(); i++)
                            {
                                if (ListWithDetTags[i] == 1)
                                    v_DetnTag++;
                            }

                            ttreeResponse->Fill();
                            std::cout << "Filled response tree for partJetidx: " << partJetidx
                                      << ", closestDetJetidx: " << closestDetJetidx << std::endl;
                        }
                    }

                    // Variation of the matching (ResponseVar tree)
                    for (int partJetidx = 0; partJetidx < nJetsPartLvl; partJetidx++)
                    {
                        // Get list with tags matched to that jet
                        std::vector<int> ListWithPartTags(nTagPartLvl, 0);

                        // Find tags in particle level jet
                        for (int tagIdx = 0; tagIdx < nTagPartLvl; tagIdx++)
                        {
                            if (m_TagJetPartLvl[tagIdx][partJetidx] == 1)
                            {
                                ListWithPartTags[tagIdx] = 1;
                            }
                        }

                        // Find first tag index
                        int idx_TagPart = -1;
                        for (size_t i = 0; i < ListWithPartTags.size(); i++)
                        {
                            if (ListWithPartTags[i] == 1)
                            {
                                idx_TagPart = i;
                                break;
                            }
                        }

                        if (idx_TagPart == -1)
                        {
                            // debugContinue(__LINE__, "Continued due to no particle level tag found in jet");
                            continue;
                        }
                        // std::cout << "idx_TagPart: " << idx_TagPart << std::endl;

                        // Get generator level properties
                        TLorentzVector p4PartJet = getFourVector("mcjet", partJetidx);
                        TLorentzVector p4PartTag = getFourVector("mctag", idx_TagPart);

                        // Find detector level tag matched to generator level tag
                        int idx_TagDet = (*mctag_idx_Dettag)[idx_TagPart]; // Det lvl tag the MC tag is associated to
                        // print all mctag_idx_Dettag
                        //  for(size_t i = 0; i < mctag_idx_Dettag->size(); i++)
                        //  {
                        //      std::cout << "mctag_idx_Dettag[" << i << "] = " << (*mctag_idx_Dettag)[i] << std::endl;
                        //  }
                        //  std::cout << "idx_TagDet: " << idx_TagDet << std::endl;
                        if (idx_TagDet < 0)
                        {
                            debugContinue(__LINE__, "Continued due to invalid detector level tag index");
                            continue;
                        }

                        int idx_jetDet = (*tag_idx_jet)[idx_TagDet]; // JetID the tag is associated to
                        // for(size_t i = 0; i < tag_idx_jet->size(); i++)
                        // {
                        //     std::cout << "tag_idx_jet[" << i << "] = " << (*tag_idx_jet)[i] << std::endl;
                        // }
                        // std::cout << "idx_jetDet: " << idx_jetDet << std::endl;
                        std::vector<int> ListWithDetTags(nTagDetLvl, 0);

                        // Find tags in detector level jet
                        for (int tagIdx = 0; tagIdx < nTagDetLvl; tagIdx++)
                        {
                            if (m_TagJetDetLvl[tagIdx][idx_jetDet] == 1)
                            {
                                ListWithDetTags[tagIdx] = 1;
                            }
                        }

                        // Get detector level properties
                        TLorentzVector p4DetJet = getFourVector("jet", idx_jetDet);
                        TLorentzVector p4DetTag = getFourVector("tag", idx_TagDet);

                        // Check if generator level tag fulfills fiducial requirements
                        bool particlesInAcceptance = false;
                        try
                        {
                            // std::cout << "idx_TagDet: " << idx_TagDet << std::endl;
                            // Get the kaon and pion from D0
                            // for(size_t i = 0; i < mctag_idx_prt1->size(); i++)
                            // {
                            //     std::cout << "mctag_idx_prt1[" << i << "] = " << (*mctag_idx_prt1)[i] << std::endl;
                            // }
                            // for(size_t i = 0; i < mctag_idx_prt2->size(); i++)
                            // {
                            //     std::cout << "mctag_idx_prt2[" << i << "] = " << (*mctag_idx_prt2)[i] << std::endl;
                            // }
                            int prt1Number = (*mctag_idx_prt1)[idx_TagPart];
                            // int prt1Number = (*mctag_idx_prt1)[idx_TagDet];
                            // std::cout << "prt1Number: " << prt1Number << std::endl;
                            int prt2Number = (*mctag_idx_prt2)[idx_TagPart];
                            // int prt2Number = (*mctag_idx_prt2)[idx_TagDet];
                            // std::cout << "prt2Number: " << prt2Number << std::endl;

                            int pid1 = (*mcprt_pid)[prt1Number];
                            int pid2 = (*mcprt_pid)[prt2Number];
                            // std::cout << "pid1: " << pid1 << ", pid2: " << pid2 << std::endl;

                            // Check if we have K(-321) and pi(211)
                            if ((pid1 == -321 && pid2 == 211) || (pid1 == 211 && pid2 == -321))
                            {
                                TLorentzVector p4MCkaon = getFourVector("mcprt",
                                                                        (pid1 == -321) ? prt1Number : prt2Number);
                                TLorentzVector p4MCpion = getFourVector("mcprt",
                                                                        (pid1 == 211) ? prt1Number : prt2Number);

                                // Appropriate D0 daughter kinematic cuts
                                if (p4MCkaon.Pt() > 0.25 && p4MCpion.Pt() > 0.25 &&
                                    p4MCkaon.P() > 2.0 && p4MCpion.P() > 2.0)
                                {
                                    particlesInAcceptance = true;
                                }
                            }
                        }
                        catch (...)
                        {
                            particlesInAcceptance = false;
                        }

                        // Skip if particles not in acceptance
                        if (!particlesInAcceptance)
                        {
                            debugContinue(__LINE__, "Continued due to particles not in acceptance");
                            continue;
                        }

                        // Fill variables for response tree
                        v_pTDet = p4DetJet.Pt();
                        v_etaDet = p4DetJet.Eta();
                        v_phiDet = p4DetJet.Phi();
                        v_nConstDet = (*jet_n_neu)[idx_jetDet] + (*jet_n_chr)[idx_jetDet];
                        v_TagPtDet = p4DetTag.Pt();
                        v_etaTagDet = p4DetTag.Eta();
                        v_phiTagDet = p4DetTag.Phi();
                        v_TagMDet = p4DetTag.M();
                        v_zTDet = p4DetTag.Pt() / p4DetJet.Pt();
                        v_pTPart = p4PartJet.Pt();
                        v_etaPart = p4PartJet.Eta();
                        v_phiPart = p4PartJet.Phi();

                        v_nConstPart = (*mcjet_n_neu)[partJetidx] + (*mcjet_n_chr)[partJetidx];
                        v_TagPtPart = p4PartTag.Pt();
                        v_etaTagPart = p4PartTag.Eta();
                        v_phiTagPart = p4PartTag.Phi();
                        v_TagMPart = p4PartTag.M();

                        if (p4PartJet.Pt() > 0)
                        {
                            v_zTPart = p4PartTag.Pt() / p4PartJet.Pt();
                        }
                        else
                        {
                            v_zTPart = 0;
                        }

                        v_dR = p4PartJet.DeltaR(p4DetJet);

                        v_isPrimary2 = (*mctag_isPrimary)[idx_TagPart];

                        // Count tags
                        v_PartnTag = 0;
                        v_DetnTag = 0;

                        for (int i = 0; i < nTagPartLvl; i++)
                        {
                            if (m_TagJetPartLvl[i][partJetidx] == 1)
                            {
                                v_PartnTag++;
                            }
                        }

                        for (int i = 0; i < nTagDetLvl; i++)
                        {
                            if (m_TagJetDetLvl[i][idx_jetDet] == 1)
                            {
                                v_DetnTag++;
                            }
                        }

                        ttreeResponse2->Fill();
                    }
                }
            }
        }
    }
    if (verboseOut)
        std::cout << "Filled PartLvl tree with " << ttree->GetEntries() << " entries." << std::endl;
    // After the event loop - save the trees
    tfile->cd();
    if (verboseOut)
        std::cout << "Writing PartLvl tree to file..." << std::endl;
    ttree->Write();
    if (verboseOut)
        std::cout << "PartLvl tree written." << std::endl;

    if (isMC)
    {
        ttreeResponse->Write();
        ttreeResponse2->Write();
        mcttree->Write();
    }

    tfile->Write();

    // Clean up - the actual deletion happens in the destructor
    std::cout << "Finished filtering file: " << fFileName << ".root" << std::endl;
    std::cout << "Saved into new file: " << foutName << std::endl;
    std::cout << "Found: " << nTagPartLvlEvt << " Evts with PartLvl tags, found: " << nTagDetLvlEvt << " Evts with DetLvl Tags" << std::endl;
}

TTree *fillPartLvlTree(TTree *mcttree, TTree *origTree,
                       float &v_MCetaPart, float &v_MCphiPart, float &v_MCpTPart,
                       float &v_MCnConstPart, float &v_MCetaTagPart, float &v_MCphiTagPart,
                       float &v_MCTagPtPart, float &v_MCTagMPart, float &v_MCzTPart,
                       float &v_MCtag_lifetPart, float &v_isPrimary, float &v_isDetTagRec)
{
    // Check for null pointers before accessing
    if (!mctag_pid || !mctag_idx_jet || !mctag_idx_Dettag ||
        !mctag_isPrimary || !mctag_idx_pvr || !mcpvr_z ||
        !mctag_z || !mctag_pz)
    {
        std::cout << "ERROR: One or more MC branch pointers are null!" << std::endl;
        return mcttree;
    }
    for (int iTagGen = 0; iTagGen < (int)mctag_pid->size(); iTagGen++)
    {
        try
        {
            int idx_jetGen = (*mctag_idx_jet)[iTagGen]; // This is the jetID the tag is associated to
            if (idx_jetGen == -1)
            {
                debugContinue(__LINE__, "Continued due to idx_jetGen == -1");
                continue;
            }

            // Check if that tag has a corresponding tag reconstructed at detector level
            int idx_tagDet = -1;
            if (iTagGen < static_cast<int>(mctag_idx_Dettag->size()))
            {
                idx_tagDet = (*mctag_idx_Dettag)[iTagGen];
            }
            else
            {
                debugContinue(__LINE__, "iTagGen out of bounds for mctag_idx_Dettag");
                continue;
            }

            int isDetRec = 0;

            if (idx_tagDet > -1 && idx_tagDet < static_cast<int>(tag_idx_MCtag->size()))
            {
                isDetRec = 1;
                int idx_MCTEST = (*tag_idx_MCtag)[idx_tagDet];
                if (iTagGen != idx_MCTEST)
                {
                    std::cout << "ERROR: these two numbers should be the same!" << std::endl;
                    std::cout << "iTagMC: " << iTagGen << std::endl;
                    std::cout << "id ass. Det Tag: " << idx_tagDet << std::endl;
                    std::cout << "id ass. MC Tag of the Det. Tag: " << idx_MCTEST << std::endl;
                }
            }

            // Safe access to mctag_isPrimary
            int isPrim = 0;
            if (iTagGen < static_cast<int>(mctag_isPrimary->size()))
            {
                isPrim = (*mctag_isPrimary)[iTagGen];
            }

            TLorentzVector p4PartTag = getFourVector("mctag", iTagGen);
            TLorentzVector p4PartJet = getFourVector("mcjet", idx_jetGen);

            // Safe access to mctag_idx_pvr
            int MCpvrNumber = -1;
            if (iTagGen < static_cast<int>(mctag_idx_pvr->size()))
            {
                MCpvrNumber = (*mctag_idx_pvr)[iTagGen];
            }
            else
            {
                debugContinue(__LINE__, "iTagGen out of bounds for mctag_idx_pvr");
                continue;
            }

            // Check if MCpvrNumber is valid for mcpvr_z
            if (MCpvrNumber < 0 || MCpvrNumber >= static_cast<int>(mcpvr_z->size()))
            {
                debugContinue(__LINE__, "MCpvrNumber out of bounds for mcpvr_z");
                continue;
            }

            double MCpvrZ = (*mcpvr_z)[MCpvrNumber]; // Primary vertex of particle
            double MCdecayZ = (*mctag_z)[iTagGen];   // Decay vertex of particle
            double MCpz = (*mctag_pz)[iTagGen];      // Z-component of momentum vector

            double lifeTime = calcLifeTime(p4PartTag.M(), MCpvrZ, MCdecayZ, MCpz);

            // Check if the particles at generator level are in the fiducial acceptance
            int pionAcc = 0;
            int kaonAcc = 0;

            try
            {
                // Get the D0 daughter particles
                int prt1Number = (*mctag_idx_prt1)[iTagGen];
                int prt2Number = (*mctag_idx_prt2)[iTagGen];

                int pid1 = (*mcprt_pid)[prt1Number];
                int pid2 = (*mcprt_pid)[prt2Number];

                // Check for pion (PID 211) and kaon (PID 321)
                TLorentzVector p4Pion;
                TLorentzVector p4Kaon;

                // Identify which particle is which
                if (std::abs(pid1) == 211 && std::abs(pid2) == 321)
                {
                    // First is pion, second is kaon
                    p4Pion = getFourVector("mcprt", prt1Number);
                    p4Kaon = getFourVector("mcprt", prt2Number);
                    pionAcc = 1;
                    kaonAcc = 1;
                }
                else if (std::abs(pid1) == 321 && std::abs(pid2) == 211)
                {
                    // First is kaon, second is pion
                    p4Kaon = getFourVector("mcprt", prt1Number);
                    p4Pion = getFourVector("mcprt", prt2Number);
                    pionAcc = 1;
                    kaonAcc = 1;
                }

                // Apply kinematic cuts to pion and kaon
                if (pionAcc == 1 && kaonAcc == 1)
                {
                    // D0 daughter acceptance cuts
                    if (!(p4Pion.Pt() > 0.25 && p4Kaon.Pt() > 0.25 &&
                          p4Pion.P() > 2.0 && p4Kaon.P() > 2.0))
                    {
                        pionAcc = 0;
                        kaonAcc = 0;
                    }
                }
            }
            catch (...)
            {
                pionAcc = -1;
                kaonAcc = -1;
            }

            // Do not accept events at generator level if the particles are not within the fiducial volume cuts
            if (pionAcc < 1 || kaonAcc < 1)
            {
                debugContinue(__LINE__, "Continued due to particles not in acceptance at generator level");
                continue;
            }

            // Set output variables for MC tree
            v_MCetaPart = p4PartJet.Eta();
            v_MCphiPart = p4PartJet.Phi();
            v_MCpTPart = p4PartJet.Pt();

            v_MCnConstPart = (*mcjet_n_neu)[iTagGen] + (*mcjet_n_chr)[iTagGen];
            v_MCetaTagPart = p4PartTag.Eta();
            v_MCphiTagPart = p4PartTag.Phi();
            v_MCTagPtPart = p4PartTag.Pt();
            v_MCTagMPart = p4PartTag.M();

            if (p4PartJet.Pt() > 0)
            {
                v_MCzTPart = p4PartTag.Pt() / p4PartJet.Pt();
            }
            else
            {
                v_MCzTPart = 0;
            }

            v_MCtag_lifetPart = lifeTime;
            v_isPrimary = isPrim;
            v_isDetTagRec = isDetRec;

            mcttree->Fill();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in fillPartLvlTree: " << e.what() << std::endl;
            continue;
        }
        catch (...)
        {
            std::cerr << "Unknown exception in fillPartLvlTree" << std::endl;
            continue;
        }
    }
    return mcttree;
}
void nTupleCreator(TString fFileName = "", int inputMC = 1)
{
    float printLvl = 10.0; // Default print level

    if (fFileName == "")
    {
        // Use default files if empty string provided
        if (inputMC)
        {
            fFileName = "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_21_MC_output_D0FF.root";
        }
        else
        {
            fFileName = "/media/niviths/SSD2/lhcb_analysis_SSD/20250501_Pbp_data/Pbp_data.root";
        }
    }

    // Convert TString to std::string for compatibility with FilterObject constructor
    std::string fileNameStr = fFileName.Data();

    // Check if this is a glob pattern with wildcards
    // if (fileNameStr.find_first_of("*?[") != std::string::npos) {
    //     // Get matching files
    //     std::vector<std::string> fileList = getFilesFromPattern(fileNameStr);
    //     if (fileList.empty()) {
    //         std::cerr << "ERROR: No files match pattern " << fileNameStr << std::endl;
    //         return;
    //     }

    //     std::cout << "Processing " << fileList.size() << " file(s) matching pattern: " << fileNameStr << std::endl;
    //     FilterObject *newFilter = new FilterObject(fileList, printLvl, inputMC);
    //     newFilter->filter();
    //     delete newFilter;
    // } else {
    // Single file - use the std::string constructor
    FilterObject *newFilter = new FilterObject(fileNameStr, printLvl, inputMC);
    newFilter->filter();
    delete newFilter;
    // }
}

// // Add an overload that directly accepts a std::string for backward compatibility
// void nTupleCreator(const std::string& fFileName = "", int inputMC = 1)
// {
//     // Just call the TString version after conversion
//     nTupleCreator(TString(fFileName.c_str()), inputMC);
// }