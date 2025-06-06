#ifndef NTUPLES_JET_ANALYSIS_RESPONSE_NO_VECT_H
#define NTUPLES_JET_ANALYSIS_RESPONSE_NO_VECT_H

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TRandom.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <array>
#include <glob.h>

// Forward declarations
class TFile;
class TTree;
class TChain;
class TLorentzVector;

/**
 * @brief Get a list of files that match a glob pattern
 *
 * @param pattern File pattern with wildcard characters
 * @return std::vector<std::string> List of matching file paths
 */
std::vector<std::string> getFilesFromPattern(const std::string& pattern);

/**
 * @brief Check if a branch exists in a tree
 * 
 * @param tree Tree to check
 * @param branchName Name of branch to look for
 * @return bool True if branch exists, false otherwise
 */
bool branchExists(TTree* tree, const char* branchName);

/**
 * @brief Function to find the closest jet in the detector level given a particle level jet
 * 
 * @param matrix Distance matrix between particle level and detector level jets
 * @param row Row index in the matrix (particle level jet index)
 * @param version Version of the matching algorithm to use
 * @return int Index of the closest detector level jet, or -1 if not found
 */
int findClosest(const std::vector<std::vector<double>>& matrix, int row, int version);

/**
 * @brief Get a four-vector for a particle in the event tree
 * 
 * @param type Type of particle ("mcjet", "mctag", "mcprt", "jet", "tag", "prt")
 * @param index Index of the particle in the corresponding collection
 * @return TLorentzVector Four-vector of the requested particle
 */
TLorentzVector getFourVector(const std::string& type, int index);

/**
 * @brief Calculate lifetime of a particle
 * 
 * @param mass Mass of the particle
 * @param pvrZ Primary vertex z-coordinate
 * @param decayZ Decay vertex z-coordinate
 * @param pz Z-component of momentum
 * @param debug Whether to print debug information
 * @return double Calculated lifetime in picoseconds
 */
double calcLifeTime(double mass, double pvrZ, double decayZ, double pz, bool debug = false);

/**
 * @brief Check if a particle is within detector acceptance
 * 
 * @param particleVect Four-vector of the particle
 * @param Type Type of the particle (0 = general, 1 = jet, 2 = D0)
 * @return bool True if particle is accepted, false otherwise
 */
bool isAccepted(const TLorentzVector& particleVect, int Type = 0);

/**
 * @brief Fill particle level tree with unbiased data
 * 
 * @param mcttree Output tree for MC truth data
 * @param origTree Input tree with original data
 * @param variables Various output variables to fill
 * @return TTree* Updated MC truth tree
 */
TTree* fillPartLvlTree(TTree* mcttree, TTree* origTree, 
                     float& MCetaPart, float& MCphiPart, float& MCpTPart, 
                     float& MCnConstPart, float& MCetaTagPart, float& MCphiTagPart,
                     float& MCTagPtPart, float& MCTagMPart, float& MCzTPart,
                     float& MCtag_lifetPart, float& isPrimary, float& isDetTagRec);

/**
 * @brief Helper function for debugging continue statements
 * 
 * @param line Source code line number
 * @param reason Reason for continuing
 */
void debugContinue(int line, const std::string& reason);

/**
 * @brief Main filter class for D0 fragmentation analysis
 */
class FilterObject {
public:
    /**
     * @brief Construct a new Filter Object from a single file
     * 
     * @param fFileName Input file name
     * @param printLvl Print level for progress reporting
     * @param inputMC Flag indicating if input is MC
     */
    FilterObject(const std::string& fFileName, float printLvl, int inputMC);
    
    /**
     * @brief Construct a new Filter Object from multiple files
     * 
     * @param fileNames Vector of input file names
     * @param printLvl Print level for progress reporting
     * @param inputMC Flag indicating if input is MC
     */
    FilterObject(const std::vector<std::string>& fileNames, float printLvl, int inputMC);
    
    /**
     * @brief Destroy the Filter Object
     */
    ~FilterObject();
    
    /**
     * @brief Main filtering function
     */
    void filter();
    
    /**
     * @brief Calculate efficiency correction weights for D0 → K- π+ decay
     * 
     * @param kaonVector Four-vector of kaon
     * @param pionVector Four-vector of pion
     * @param type Type of correction to apply
     * @return double Efficiency correction weight
     */
    double getD0EffCorrWeight(const TLorentzVector& kaonVector, const TLorentzVector& pionVector, int type = -1);
    
    /**
     * @brief Calculate efficiency correction weights for J/Psi
     * 
     * @param piVector1 Four-vector of first pion
     * @param piVector2 Four-vector of second pion
     * @param muonVect1 Four-vector of first muon
     * @param muonVect2 Four-vector of second muon
     * @param triggVect Four-vector for trigger calculations
     * @param type Type of correction to apply
     * @return double Efficiency correction weight
     */
    double getEffCorrWeight(const TLorentzVector& piVector1, const TLorentzVector& piVector2,
                          const TLorentzVector& muonVect1, const TLorentzVector& muonVect2,
                          const TLorentzVector& triggVect, int type = -1);
    
    /**
     * @brief Create map variations for uncertainty evaluation
     * 
     * @param inputHisto Input histogram
     * @param TEffObject Optional efficiency object
     * @return TH2D* Varied map
     */
    TH2D* createMapVariation(TH2D* inputHisto, TEfficiency* TEffObject = nullptr);
    
    /**
     * @brief Fit trigger efficiency histograms
     * 
     * @param triggerHist Trigger efficiency histogram
     * @return std::vector<TF1*> Fit functions for different eta bins
     */
    std::vector<TF1*> fitTriggerEff(TH2D* triggerHist);

private:
    bool isMC;
    float printLvl;
    int matchingVersion;
    
    std::string fFileName;
    std::string foutName;
    
    TFile* tfileOriginal;  // Only used for single file mode
    TTree* ttreeOriginal;  // Points to TTree or TChain depending on input
    TFile* tfile;
    TTree* ttree;
    TTree* ttreeResponse;
    TTree* ttreeResponse2;
    TTree* mcttree;
    
    // Efficiency maps
    TH2D* pionMapReco;
    TH2D* pionMapReco_Var;
    TH2D* pionMapSelHist;
    TH2D* pionMapSelHist_Var;
    TH2D* MuonMap;
    TH2D* MuonMap_Var;
    TH2D* TriggerSelCorrMap;
    TH2D* TriggerSelCorrMap_Var;
    TH2D* TriggerEffMap2D;
    TH2D* TriggerEffMap2D_Var;
    std::vector<TF1*> functionList;
    std::vector<TF1*> functionList_Var;
};

/**
 * @brief Main function for D0 fragmentation analysis with a single file
 * 
 * @param filePattern Input file name or glob pattern
 * @param inputMC Flag indicating if input is MC
 */
void nTupleCreator(int inputMC);
void nTupleCreator(TString infile, int inputMC);

/**
 * @brief Main function for D0 fragmentation analysis with multiple files
 * 
 * @param fileNames Vector of input file names
 * @param inputMC Flag indicating if input is MC
 */
void nTupleCreator(const std::vector<std::string>& fileNames, int inputMC);

// Alias for backward compatibility
void ntuplesJetAnalysis_Response_noVect(const std::string& fFileName, float printLvl, int inputMC);

#endif // NTUPLES_JET_ANALYSIS_RESPONSE_NO_VECT_H