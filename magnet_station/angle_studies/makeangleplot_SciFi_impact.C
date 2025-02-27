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

void makeangleplot_SciFi_impact(TString inputfile = "/home/niviths/Downloads/magnetStationSims/minimumBias_MS_MagDown_3120plus.root") // tilt
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
    int npartwithhitinsupp = 0;
    int npartwithhitinFTandMS = 0;
    int npartwithhitinMS = 0;
    int npartswithhitinFT = 0;
    int npartwithouthitinMS = 0;

    int nHits_MSbars = 0;
    while (tree.Next())
    {
        // cout << "Entry: " << curr_entry << endl;
        curr_entry++;
        // if(curr_entry < 3856){
        //     continue;
        // }
        // if (curr_entry > 1500 || onceMsg)
        // {
        //     if (onceMsg)
        //     {
        //         onceMsg = 0;
        //         cout << "processing only 1500 events" << endl;
        //     } else {
        //         break;
        //     }
        // }
        if (curr_entry % 100 == 0)
            cout << "Entry: " << curr_entry << endl;

        for (auto i = 0U; i < pt.GetSize(); ++i)
        {
            // if (p[i] > 5000)
            // {
            //     // continue;
            // }
            if(eta[i] < 1.1 || eta[i] > 5.0){
                continue;
            }
            //TODO REMOVE THIS LINE BELOW
            // if (parent_key[i] > -2000)
            //     continue;

            bool particleHitSupport = 0;
            bool particleHitMS = 0;
            bool particleDaughterHitSupport = 0;
            bool particleDaughterHitMS = 0;
            // loop over size of ms_px to get all the hits per event
            for (auto j = 0U; j < ms_px.GetSize(); ++j)
            {
                if (ms_id[j] != key[i])
                    continue;
                // require hits within MS acceptance
                if (ms_vz[j] < 2500 || ms_vz[j] > maxZ+500)
                {
                    if (verbositysetting)
                        cout << "Hit outside of acceptance" << endl;
                    continue;
                }
                else
                {
                    if (verbositysetting)
                        cout << "Hit in acceptance" << endl;
                }

                int station = (ms_bitID[j] >> 8) & 0x7;  // 3 bit
                int module = (ms_bitID[j] >> 11) & 0xF;  // 4 bit
                int layer = (ms_bitID[j] >> 15) & 0x7;   // 3 bit
                int segment = (ms_bitID[j] >> 18) & 0xF; // 4 bit
                int bar = (ms_bitID[j] >> 22) & 0x3F;    // 6 bit
                // int bar = (ms_bitID[j] >> 22) & 0xFF;    // 8 bit
                int isfiber = (ms_bitID[j] >> 28) & 0x3; // 2 bit
                int issupp = (ms_bitID[j] >> 30) & 0x3; // 2 bit
                // cout << "isfiber: " << isfiber << " issupp: " << issupp << endl;
                if (bar > maxbar)
                    maxbar = bar;
                // cout << "\t\tStation: " << station << " Module: " << module << " Layer: " << layer << " Segment: " << segment << " Bar: " << bar << " Energy: " << ms_energy[j] << endl;

                if (!isfiber && !issupp){
                    nHits_MSbars++;
                }
                particleHitMS = 1;


                if (isfiber || issupp)
                {
                    // nHits_Support++;
                    hYZHitsinSupport->Fill(ms_vz[j], ms_vy[j]);
                    particleHitSupport = 1;
                }
                if(child1_key[i] > -2000 || child2_key[i] > -2000){
                    for (auto k = 0U; k < ms_px.GetSize(); ++k)
                    {
                        // cout << "ms_id[k]: " << ms_id[k] << " parent_key[i]: " << parent_key[i] << endl;
                        if (ms_id[k] != child1_key[i] && ms_id[k] != child2_key[i])
                            continue;
                        int isfiberDaughter = (ms_bitID[k] >> 28) & 0x3; // 2 bit
                        int issuppDaughter = (ms_bitID[k] >> 30) & 0x3; // 2 bit
                        // cout << "isfiberDaughter: " << isfiberDaughter << " issuppDaughter: " << issuppDaughter << endl;
                        if (isfiberDaughter || issuppDaughter){
                            particleDaughterHitSupport = 1;
                        } else {
                            particleDaughterHitMS = 1;
                        }
                    }
                }

            } // NOTE END LOOP OVER MS HITS

            // loop over ft hits and find matching hits to the current particle and fill xy histograms
            bool foundparthitinFT = 0;
            bool foundparthitinFTDaughter = 0;
            bool particlehitinFT = 0;
            // NOTE LOOP OVER FT HITS
            for (auto k = 0U; k < ft_id.GetSize(); ++k)
            {
                if (ft_id[k] == key[i])
                {
                    if(particleHitSupport)hXYHitsinFTFromSupport->Fill(ft_vx[k], ft_vy[k]);
                    
                    if(particleHitSupport || particleHitMS){
                        foundparthitinFT = 1;
                        hXYHitsinFTFromMS->Fill(ft_vx[k], ft_vy[k]);
                    }
                    particlehitinFT = 1;
                }
                //prevent code crash by check if child1_key or child2_key pointers are good
                
                try {
                    if (ft_id[k] == child1_key[i] || ft_id[k] == child2_key[i])
                    {
                        if(particleDaughterHitSupport) hXYHitsinFTFromSupportDaughter->Fill(ft_vx[k], ft_vy[k]);
                        if(particleDaughterHitSupport || particleDaughterHitMS){
                            foundparthitinFTDaughter = 1;
                            hXYHitsinFTFromMSDaughter->Fill(ft_vx[k], ft_vy[k]);
                        }

                    }
                } catch (const std::exception& e) {
                    std::cerr << e.what() << '\n';
                }
            } // NOTE END LOOP OVER FT HITS

            if(foundparthitinFT){
                hAccPartHitSuppAndFT->Fill(1 / (pt[i] / 1000), eta[i]);
            }
            if(foundparthitinFTDaughter){
                hAccPartHitSuppAndFTDaughter->Fill(1 / (pt[i] / 1000), eta[i]);
            }

            // }
            if(foundparthitinFT || foundparthitinFTDaughter){
                npartwithhitinFTandMS++;
            }
            if(particlehitinFT)
                npartswithhitinFT++;
            if (particleHitSupport || particleDaughterHitSupport)
                npartwithhitinsupp++;
            
            if (particleHitMS || particleDaughterHitMS)
                npartwithhitinMS++;
            else
                npartwithouthitinMS++;
        } // NOTE END LOOP OVER PARTICLES
    } // NOTE END LOOP OVER EVENTS
    if(npartwithouthitinMS>0)cout << "Particles with hits in support: " << npartwithhitinsupp << " Particles without hits in MS: " << npartwithouthitinMS << " Ratio: " << (float)npartwithhitinsupp / (float)npartwithouthitinMS << endl;
    if(npartwithouthitinMS>0) cout << "Particles with hits in support+MS: " << npartwithhitinMS << " Particles without hits in MS: " << npartwithouthitinMS << " Ratio: " << (float)npartwithhitinMS / (float)npartwithouthitinMS << endl;
    if(npartswithhitinFT>0) cout << "Particles with hits in support+MS+FT: " << npartwithhitinFTandMS << " Particles with hits in FT: " << npartswithhitinFT << " Ratio: " << (float)npartwithhitinFTandMS / (float)npartswithhitinFT << endl;

    //plot hYZHitsinSupport
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    gStyle->SetOptStat(0);
    hYZHitsinSupport->GetXaxis()->SetTitle("Z [mm]");
    hYZHitsinSupport->GetYaxis()->SetTitle("Y [mm]");
    hYZHitsinSupport->Draw("colz");
    c1->SaveAs("hYZHitsinSupport.pdf");

    //plot hXYHitsinFTFromSupport
    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    gStyle->SetOptStat(0);
    hXYHitsinFTFromSupport->GetXaxis()->SetTitle("X [mm]");
    hXYHitsinFTFromSupport->GetYaxis()->SetTitle("Y [mm]");
    hXYHitsinFTFromSupport->Draw("colz");
    c2->SaveAs("hXYHitsinFTFromSupport.pdf");

    //plot hXYHitsinFTFromSupportDaughter
    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    gStyle->SetOptStat(0);
    hXYHitsinFTFromSupportDaughter->GetXaxis()->SetTitle("X [mm]");
    hXYHitsinFTFromSupportDaughter->GetYaxis()->SetTitle("Y [mm]");
    hXYHitsinFTFromSupportDaughter->Draw("colz");
    c3->SaveAs("hXYHitsinFTFromSupportDaughter.pdf");

    //plot hAccPartHitSuppAndFT
    TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
    gStyle->SetOptStat(0);
    hAccPartHitSuppAndFT->GetXaxis()->SetTitle("1/pT [1/GeV]");
    hAccPartHitSuppAndFT->GetYaxis()->SetTitle("#eta");
    hAccPartHitSuppAndFT->Draw("colz");
    c4->SaveAs("hAccPartHitSuppAndFT.pdf");

    //plot hAccPartHitSuppAndFTDaughter
    TCanvas *c5 = new TCanvas("c5", "c5", 800, 600);
    gStyle->SetOptStat(0);
    hAccPartHitSuppAndFTDaughter->GetXaxis()->SetTitle("1/pT [1/GeV]");
    hAccPartHitSuppAndFTDaughter->GetYaxis()->SetTitle("#eta");
    hAccPartHitSuppAndFTDaughter->Draw("colz");
    c5->SaveAs("hAccPartHitSuppAndFTDaughter.pdf");

    //plot hXYHitsinFTFromMS
    TCanvas *c6 = new TCanvas("c6", "c6", 800, 600);
    gStyle->SetOptStat(0);
    hXYHitsinFTFromMS->GetXaxis()->SetTitle("X [mm]");
    hXYHitsinFTFromMS->GetYaxis()->SetTitle("Y [mm]");
    hXYHitsinFTFromMS->Draw("colz");
    c6->SaveAs("hXYHitsinFTFromMS.pdf");

    //plot hXYHitsinFTFromMSDaughter
    TCanvas *c7 = new TCanvas("c7", "c7", 800, 600);
    gStyle->SetOptStat(0);
    hXYHitsinFTFromMSDaughter->GetXaxis()->SetTitle("X [mm]");
    hXYHitsinFTFromMSDaughter->GetYaxis()->SetTitle("Y [mm]");
    hXYHitsinFTFromMSDaughter->Draw("colz");
    c7->SaveAs("hXYHitsinFTFromMSDaughter.pdf");

    // //plot hAccPartHitMSAndFT
    // TCanvas *c8 = new TCanvas("c8", "c8", 800, 600);
    // gStyle->SetOptStat(0);
    // hAccPartHitMSAndFT->GetXaxis()->SetTitle("1/pT [1/GeV]");
    // hAccPartHitMSAndFT->GetYaxis()->SetTitle("#eta");
    // hAccPartHitMSAndFT->Draw("colz");
    // c8->SaveAs("hAccPartHitMSAndFT.pdf");

    // //plot hAccPartHitMSAndFTDaughter
    // TCanvas *c9 = new TCanvas("c9", "c9", 800, 600);
    // gStyle->SetOptStat(0);
    // hAccPartHitMSAndFTDaughter->GetXaxis()->SetTitle("1/pT [1/GeV]");
    // hAccPartHitMSAndFTDaughter->GetYaxis()->SetTitle("#eta");
    // hAccPartHitMSAndFTDaughter->Draw("colz");
    // c9->SaveAs("hAccPartHitMSAndFTDaughter.pdf");

    cout << "Max bar: " << maxbar << endl;


    fin.Close();

}