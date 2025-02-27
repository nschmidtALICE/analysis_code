// this is a .root macro to plot the incident angle of the particles on the Magnet Station
// it uses the actual hit position and associated momentum vector to calculate the angle
// since the detector plane is rotated for each segment, the angle is calculated with respect to the rotated plane
// contact person: Nicolas Schmidt, schmidt_n@lanl.gov

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <iostream>
#include <vector>
#include <cmath>

void makeangleplot_SciFi_impact2(TString inputfile = "/home/niviths/Downloads/magnetStationSims/minimumBias_MS_MagDown_3200plus.root") // tilt
// void makeangleplot_SciFi_impact2(TString inputfile = "/home/niviths/Downloads/magnetStationSims/minimumBias_MS_MagDown_3120plus.root") // tilt
{
    TFile fin(inputfile, "READ");
    TTreeReader tree("ntup", &fin);
    double modSpacing = 300 * TMath::Cos(0.42);
    double minZ = 3500 + 2 * 500 + 300;
    double maxZ = 7700;

    // get ms_px, ms_py, ms_pz and ms_vx, ms_vy, ms_vz, as well as p and ms_time from the tree
    TTreeReaderArray<float> ms_px(tree, "ms_px");
    TTreeReaderArray<float> ms_py(tree, "ms_py");
    TTreeReaderArray<float> ms_pz(tree, "ms_pz");
    TTreeReaderArray<float> ms_vx(tree, "ms_vx");
    TTreeReaderArray<float> ms_vy(tree, "ms_vy");
    TTreeReaderArray<float> ms_vz(tree, "ms_vz");
    TTreeReaderArray<float> ms_energy(tree, "ms_energy");
    TTreeReaderArray<int>   ms_npart(tree, "ms_npart");
    TTreeReaderArray<int>   ms_id(tree, "ms_id");
    TTreeReaderArray<float> ft_vx(tree, "ft_vx");
    TTreeReaderArray<float> ft_vy(tree, "ft_vy");
    TTreeReaderArray<float> ft_vz(tree, "ft_vz");
    TTreeReaderArray<int>   ft_id(tree, "ft_id");
    TTreeReaderArray<int>   key(tree, "key");
    TTreeReaderArray<int>   parent_key(tree, "parent_key");
    TTreeReaderArray<int>   child1_key(tree, "child1_key");
    TTreeReaderArray<int>   child2_key(tree, "child2_key");
    TTreeReaderArray<int>   nFThits(tree, "nFThits");
    TTreeReaderArray<float> p(tree, "p");
    TTreeReaderArray<int>   pid(tree, "pid");
    TTreeReaderArray<float> pt(tree, "pt");
    TTreeReaderArray<float> eta(tree, "eta");
    TTreeReaderArray<float> vx(tree, "vx");
    TTreeReaderArray<float> vy(tree, "vy");
    TTreeReaderArray<float> vz(tree, "vz");
    TTreeReaderArray<float> ms_time(tree, "ms_time");
    TTreeReaderArray<int>   ms_bitID(tree, "ms_segment");
    const int nModules = 9;
    

    TH2F *hYZHitsinSupport = new TH2F("hYZHitsinSupport", "YZ of hits in plane", 125, minZ-1500, maxZ, 100, -1500, 1500);
    TH2F *hXYHitsinFTFromSupport = new TH2F("hXYHitsinFTFromSupport", "XY of hits in FT from support", 100, -5000, 5000, 100, -5000, 5000);
    TH2F *hXYHitsinFTFromSupportDaughter = new TH2F("hXYHitsinFTFromSupportDaughter", "XY of hits in FT from support Daughter", 100, -5000, 5000, 100, -5000, 5000);
    TH2F *hXYHitsinFTFromMS = new TH2F("hXYHitsinFTFromMS", "XY of hits in FT from MS", 100, -5000, 5000, 100, -5000, 5000);
    TH2F *hXYHitsinFTFromMSDaughter = new TH2F("hXYHitsinFTFromMSDaughter", "XY of hits in FT from MS Daughter", 100, -5000, 5000, 100, -5000, 5000);
    TH2F *hAccPartHitSuppAndFT = new TH2F("hAccPartHitSuppAndFT", "Acceptance of particles hitting MS Supports and FT tracker", 125, 0, 20, 60, 1.8, 4.5);
    TH2F *hAccPartHitSuppAndFTDaughter = new TH2F("hAccPartHitSuppAndFTDaughter", "Acceptance of particles hitting MS Supports and FT tracker from Daughter", 125, 0, 20, 60, 1.8, 4.5);
    TH2F *hAccPartHitMSAndFT = new TH2F("hAccPartHitMSAndFT", "Acceptance of particles hitting MS and FT tracker", 125, 0, 20, 60, 1.8, 4.5);
    TH2F *hAccPartHitMSAndFTDaughter = new TH2F("hAccPartHitMSAndFTDaughter", "Acceptance of particles hitting MS and FT tracker from Daughter", 125, 0, 20, 60, 1.8, 4.5);

    int curr_entry = 0;
    bool verbositysetting = 0;
    // loop over tree
    int maxbar = 0;
    bool onceMsg = 1;
    int npartwithhitinFTandMS = 0;

    int npartwithhitinAnyMS = 0;
    int npartwithhitinMS = 0;
    int npartwithhitinMSsupp = 0;

    int npartwithMotherhitinMS = 0;
    int npartwithMotherhitinMSsupp = 0;

    int npartswithhitinFT = 0;
    int npartwithouthitinMS = 0;

    int nHits_MSbars = 0;
    while (tree.Next())
    {
        // cout << "Entry: " << curr_entry << endl;
        curr_entry++;
        if (curr_entry % 100 == 0)
            cout << "Entry: " << curr_entry << endl;

        for (auto i = 0U; i < pt.GetSize(); ++i)
        {
            // if(eta[i] < 1.1 || eta[i] > 5.0){
            //     continue;
            // }
            //TODO REMOVE THIS LINE BELOW
            // if (parent_key[i] > -2000)
            //     continue;

            bool particleDaughterHitSupport = 0;
            bool particleDaughterHitMS = 0;            

            bool particleHitAnyMS = 0;
            bool particleHitMS = 0;
            bool particleMotherHitMS = 0;

            bool particleHitSupport = 0;
            bool particleMotherHitSupport = 0;
            bool particlehitinFT = 0;

            // loop over all hits in ft
            for (auto k = 0U; k < ft_id.GetSize(); ++k)
            {
                if (ft_id[k] == key[i])
                {
                    particlehitinFT = 1;
                }
            } // NOTE END LOOP OVER FT HITS

            // loop over all hits in ms
            for (auto k = 0U; k < ms_id.GetSize(); ++k)
            {
                //get isfiber and issupp
                int isfiber = (ms_bitID[k] >> 28) & 0x3; // 2 bit
                int issupp = (ms_bitID[k] >> 30) & 0x3; // 2 bit

                if (ms_id[k] == key[i])
                {
                    particleHitMS = 1;
                    if(isfiber || issupp){
                        particleHitSupport = 1;
                    }
                }
                if (parent_key[i]>-2000 && ms_id[k] == parent_key[i])// || ms_id[k] == child1_key[i] || ms_id[k] == child2_key[i])
                {
                    particleMotherHitMS = 1;
                    if(isfiber || issupp){
                        particleMotherHitSupport = 1;
                    }
                }
                if(particleHitSupport || particleMotherHitSupport || particleHitMS || particleMotherHitMS){
                    particleHitAnyMS = 1;
                }
            } // NOTE END LOOP OVER MS HITS

            if(particlehitinFT){
                npartswithhitinFT++;
            }

            if(particleHitAnyMS){
                npartwithhitinAnyMS++;
            }
            

            if(particlehitinFT && particleHitAnyMS){
                hAccPartHitMSAndFT->Fill(1 / (pt[i] / 1000), eta[i]);
                npartwithhitinFTandMS++;
            }

            if (particleHitMS)
                npartwithhitinMS++;
            if (particleMotherHitMS)
                npartwithMotherhitinMS++;
            if (particleHitSupport)
                npartwithhitinMSsupp++;
            if (particleMotherHitSupport)
                npartwithMotherhitinMSsupp++;
            if(!particleHitMS && !particleMotherHitMS && !particleHitSupport && !particleMotherHitSupport){
                npartwithouthitinMS++;
            }
        } // NOTE END LOOP OVER PARTICLES

    } // NOTE END LOOP OVER EVENTS
    // if(npartwithouthitinMS>0)cout << "Particles with hits in support: " << npartwithhitinMSsupp << " Particles without hits in MS: " << npartwithouthitinMS << " Ratio: " << (float)npartwithhitinMSsupp / (float)npartwithouthitinMS << endl;
    // if(npartwithouthitinMS>0) cout << "Particles with hits in support+MS: " << npartwithhitinMS << " Particles without hits in MS: " << npartwithouthitinMS << " Ratio: " << (float)npartwithhitinMS / (float)npartwithouthitinMS << endl;
    // if(npartswithhitinFT>0) cout << "Particles with hits in support+MS+FT: " << npartwithhitinFTandMS << " Particles with hits in FT: " << npartswithhitinFT << " Ratio: " << (float)npartwithhitinFTandMS / (float)npartswithhitinFT << endl;
    if(npartwithouthitinMS>0) cout << "Particles with hits in support+MS: " << npartwithhitinMS << " Particles without hits in MS: " << npartwithouthitinMS << " Ratio: " << (float)npartwithhitinMS / (float)npartwithouthitinMS << endl;

    float totalMShits = npartwithhitinMS + npartwithMotherhitinMS;
    float totalMShitsSupp = npartwithhitinMSsupp + npartwithMotherhitinMSsupp;

    cout << "Total MS hits: " << totalMShits << " Total MS hits in support: " << totalMShitsSupp << " Ratio: " << totalMShitsSupp / totalMShits << endl;

    //print ratio of particles that hit the FT to particles that hit the MS, Support and FT
    if(npartswithhitinFT>0) cout << "Particles with hits in support+MS+FT: " << npartwithhitinFTandMS << " Particles with hits in FT: " << npartswithhitinFT << " Ratio: " << (float)npartwithhitinFTandMS / (float)npartswithhitinFT << endl;


    //plot hAccPartHitMSAndFT
    TCanvas *c8 = new TCanvas("c8", "c8", 800, 600);
    gStyle->SetOptStat(0);
    hAccPartHitMSAndFT->GetXaxis()->SetTitle("1/pT [1/GeV]");
    hAccPartHitMSAndFT->GetYaxis()->SetTitle("#eta");
    hAccPartHitMSAndFT->Draw("colz");
    c8->SaveAs("hAccPartHitMSAndFT.pdf");


    fin.Close();

}