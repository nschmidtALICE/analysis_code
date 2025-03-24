/**
 * ROOT TTree Extractor
 * 
 * This program loads TTrees named "DecayFunTuple" from different directories in a ROOT file
 * and saves each TTree to a separate output ROOT file.
 * 
 * Usage: ./extractTrees input.root
 */

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TObject.h"

// Function to check if a directory exists in the ROOT file
bool directoryExists(TFile* file, const std::string& dirName) {
    TObject* obj = file->Get(dirName.c_str());
    return (obj != nullptr && obj->InheritsFrom("TDirectory"));
}

int tree_separator(TString strinputFile = "/media/niviths/local/lhcb_analysis/20250324_LcD0AP/c2up/00280911_00000004_1.ohf_smog2_collision24c2.root") {


    // Input filename from command line
    std::string inputFileName = strinputFile.Data();
    
    // Define the directory names to process
    // Modify these names to match the directories in your ROOT file
    std::vector<std::string> directoryNames = {
        // "D0Kpi",
        "LcpKpi",
              "BmD0mu",
              "BmD0pi",
              "DpKpipi",
              "B0JPsiKpi",
              "B0JPsiKS0",
            // #   "B0Dmmup",
            // #   "B0Dmpip",
              "DsKKpi",
              "Lb0JpsiKp",
              "Lb0JpsiL0",
              "X3872JpsiPiPi"
    };

    // Try to open the input file
    std::unique_ptr<TFile> inputFile(TFile::Open(inputFileName.c_str(), "READ"));
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Cannot open input file: " << inputFileName << std::endl;
        return 1;
    }

    std::cout << "Successfully opened input file: " << inputFileName << std::endl;
    
    // Keep track of how many TTrees were successfully extracted
    int treesExtracted = 0;

    // Process each directory
    for (const auto& dirName : directoryNames) {
        // Skip directories that don't exist in the file
        if (!directoryExists(inputFile.get(), dirName)) {
            std::cout << "Directory '" << dirName << "' not found. Skipping..." << std::endl;
            continue;
        }
        
        // Navigate to the directory
        TDirectory* dir = static_cast<TDirectory*>(inputFile->Get(dirName.c_str()));
        if (!dir) {
            std::cerr << "Error accessing directory: " << dirName << std::endl;
            continue;
        }
        
        // Try to get the TTree
        TTree* inputTree = nullptr;
        dir->GetObject("DecayFunTuple", inputTree);
        
        if (!inputTree) {
            std::cerr << "Error: DecayFunTuple not found in directory: " << dirName << std::endl;
            continue;
        }
        
        // Create output filename based on directory name
        std::string outputFileName = dirName + "_DecayFunTuple.root";
        
        // Create output file
        std::unique_ptr<TFile> outputFile(TFile::Open(outputFileName.c_str(), "RECREATE"));
        if (!outputFile || outputFile->IsZombie()) {
            std::cerr << "Error: Cannot create output file: " << outputFileName << std::endl;
            continue;
        }
        
        // Clone the tree to the output file
        TTree* outputTree = inputTree->CloneTree(-1, "fast");
        if (!outputTree) {
            std::cerr << "Error: Failed to clone tree for directory: " << dirName << std::endl;
            continue;
        }
        
        // Write the tree to the output file
        outputTree->Write();
        outputFile->Close();
        
        std::cout << "Successfully extracted DecayFunTuple from '" << dirName 
                  << "' to file: " << outputFileName << std::endl;
        treesExtracted++;
    }
    
    // Summary
    std::cout << "\nExtraction complete. Processed " << treesExtracted 
              << " trees out of " << directoryNames.size() << " directories." << std::endl;
    
    return 0;
}
