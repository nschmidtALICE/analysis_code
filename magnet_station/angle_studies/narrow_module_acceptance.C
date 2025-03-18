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
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TString.h>
#include <TLegend.h>
#include <TMath.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TGraph.h>
#include <TStyle.h>

// void makeangleplot_occupancy(TString inputfile = "/home/niviths/Downloads/magnetStationSims/2100_wtf/minimumBias_MS_MagDown_2100plus.root") // tilt
// void makeangleplot_occupancy(TString inputfile = "/home/niviths/Downloads/magnetStationSims/minimumBias_MS_MagDown_PbPb.root") // tilt
// void makeangleplot_occupancy(TString inputfile = "/home/niviths/Downloads/magnetStationSims/20250216/minimumBias_MS_MagDown_20250216.root") // tilt
// void makeangleplot_occupancy_timeBins(TString inputfile = "/home/niviths/Downloads/magnetStationSims/20250218_8cm_inward_lastPanelsNarrow/minimumBias_MS_MagDown_20250218_8cm_inward_lastPanelsNarrow.root") // tilt
// void makeangleplot_occupancy(TString inputfile = "/home/niviths/Downloads/magnetStationSims/minimumBias_MS_MagDown_10015plus.root") // tilt
// void narrow_module_acceptance(TString inputfile = "/home/niviths/Downloads/magnetStationSims/20250227_PbPb_LayerNumbersFixed/20250227_PbPb_LayerNumbersFixed.root") // tilt
// void narrow_module_acceptance(TString inputfile = "/home/niviths/Downloads/magnetStationSims/20250311_PbPb_addUTinfo/20250311_PbPb_addUTinfo.root") // tilt
void narrow_module_acceptance(TString inputfile = "/home/niviths/Downloads/magnetStationSims/20250311_pp_newOutput/20250311_pp_newOutput.root") // tilt
{
    TFile fin(inputfile, "READ");
    cout << "Opened file " << inputfile << endl;
    TTreeReader tree("ntup", &fin);
    double modSpacing = 300 * TMath::Cos(0.42);
    double minZ = 3500 + 2 * 500 + 300;
    double maxZ = 7700;

    TString outputdir = "narrow_module_acceptance_plots/";
    // make output directory
    system("mkdir -p " + outputdir);

    // get ms_px, ms_py, ms_pz and ms_vx, ms_vy, ms_vz, as well as p and ms_time from the tree
    TTreeReaderArray<float> ms_px(tree, "ms_px");
    TTreeReaderArray<float> ms_py(tree, "ms_py");
    TTreeReaderArray<float> ms_pz(tree, "ms_pz");
    TTreeReaderArray<float> ms_vx(tree, "ms_vx");
    TTreeReaderArray<float> ms_vy(tree, "ms_vy");
    TTreeReaderArray<float> ms_vz(tree, "ms_vz");
    TTreeReaderArray<float> ms_energy(tree, "ms_energy");
    TTreeReaderArray<int> ms_npart(tree, "ms_npart");
    TTreeReaderArray<int> ms_id(tree, "ms_id");
    TTreeReaderArray<float> ft_vx(tree, "ft_vx");
    TTreeReaderArray<float> ft_vy(tree, "ft_vy");
    TTreeReaderArray<float> ft_vz(tree, "ft_vz");
    TTreeReaderArray<int> ft_id(tree, "ft_id");
    TTreeReaderArray<int> key(tree, "key");
    TTreeReaderArray<int> parent_key(tree, "parent_key");
    TTreeReaderArray<int> child1_key(tree, "child1_key");
    TTreeReaderArray<int> child2_key(tree, "child2_key");
    TTreeReaderArray<int> nFThits(tree, "nFThits");
    TTreeReaderArray<float> p(tree, "p");
    TTreeReaderArray<int> pid(tree, "pid");
    TTreeReaderArray<float> pt(tree, "pt");
    TTreeReaderArray<float> eta(tree, "eta");
    TTreeReaderArray<float> ms_time(tree, "ms_time");
    TTreeReaderArray<int> ms_bitID(tree, "ms_segment");
    TTreeReaderArray<int> ievt(tree, "ievt");
    const int nModules = 9;
    int maxbar = 56;

    // histograms with 1/pt vs pseudorapidity
    TH2F *acceptanceModuleWidth[20];
    TH2F *hHitsZXMS[20];
    TH2F *hHitsZXMSModule0[20];
    for (int i = 0; i < 20; i++)
    {
        acceptanceModuleWidth[i] = new TH2F(Form("acceptanceModuleWidth%d", i), Form("Acceptance for %d less bars", 2*i), 100, 0.0, 20, 100, 1.8, 4.5);
        hHitsZXMS[i] = new TH2F(Form("hHitsZXMS%d", i), Form("Hits in Z-X plane for %d less bars", 2*i), 100, minZ, maxZ, 600, -3000, 3000);
        hHitsZXMSModule0[i] = new TH2F(Form("hHitsZXMSModule0%d", i), Form("Hits in Z-X plane for %d less bars", 2*i), 100, 4900, 5300, 600, 1200, 2000);
    }

    int maxbarcount = 0;
    int curr_entry = 0;
    bool verbositysetting = 0;
    // loop over tree
    int currentEventNumber = 0;
    cout << "Number of entries: " << tree.GetEntries() << endl;
    int nUpdateEvents = 1000;
    if (tree.GetEntries() < 1000)
    {
        nUpdateEvents = 10;
    }
    while (tree.Next())
    {

        curr_entry++;
        if (curr_entry % nUpdateEvents == 0)
            cout << "Entry: " << curr_entry << endl;

        // if(curr_entry > 1000){
        //     break;
        // }

        for (auto i = 0U; i < pt.GetSize(); ++i)
        {
            if (eta[i] < 1.8 || eta[i] > 5.0)
            {
                continue;
            }
            bool hitMS[20] = {0};
            int nhitMS[20] = {0};
            bool hitFT = 0;
            // loop over size of ms_px to get all the hits per event
            for (auto j = 0U; j < ms_px.GetSize(); ++j)
            {
                if (ms_id[j] != key[i])
                    continue;

                int station = (ms_bitID[j] >> 8) & 0x7;  // 3 bit
                int module = (ms_bitID[j] >> 11) & 0xF;  // 4 bit
                int layer = (ms_bitID[j] >> 15) & 0x7;   // 3 bit
                int segment = (ms_bitID[j] >> 18) & 0xF; // 4 bit
                int bar = (ms_bitID[j] >> 22) & 0x3F;    // 6 bits
                int isfiber = (ms_bitID[j] >> 28) & 0x3; // 2 bits
                int issupp = (ms_bitID[j] >> 30) & 0x3;  // 2 bits

                if (isfiber || issupp)
                {
                    continue;
                }

                if (bar > maxbarcount)
                {
                    maxbarcount = bar;
                }
                for (int k = 0; k < 20; k++)
                {
                    if (ms_vy[j] > 0 && bar <= maxbar - 2*k)
                    {
                        hitMS[k] = 1;
                        nhitMS[k]++;
                        hHitsZXMS[k]->Fill(ms_vz[j], ms_vx[j]);
                        // if(module == 0 && ms_vy[j] < 0){
                        //     hHitsZXMSModule0[k]->Fill(ms_vz[j], ms_vx[j]);
                        // }
                    } else if (ms_vy[j] < 0 && bar >= 2*k){
                        // acceptanceModuleWidth[k]->Fill(1000/pt[i], eta[i]);
                        hitMS[k] = 1;
                        nhitMS[k]++;
                        hHitsZXMS[k]->Fill(ms_vz[j], ms_vx[j]);
                        if(module == 0){
                            hHitsZXMSModule0[k]->Fill(ms_vz[j], ms_vx[j]);
                        }
                    }
                }

                
            } // NOTE END LOOP OVER MS HITS
            
            for(int k = 0; k < 20; k++){
                if(hitMS[k] && nhitMS[k] > 2){
                    acceptanceModuleWidth[k]->Fill(1000/pt[i], eta[i]);
                }
            }

            // loop over ft hits as well
            for (auto k = 0U; k < ft_id.GetSize(); ++k)
            {
                if (ft_id[k] == key[i])
                {
                    hitFT = 1; 
                    break;
                    // for (int l = 0; l < 20; l++)
                    // {
                    //     acceptanceModuleWidth[l]->Fill(1000/pt[i], eta[i]);
                    // }
                }
            } // NOTE END LOOP OVER FT HITS

            if(hitFT){
                for(int k = 0; k < 20; k++){
                    acceptanceModuleWidth[k]->Fill(1000/pt[i], eta[i]);
                }
            }
        } // NOTE END LOOP OVER PARTICLES
    } // NOTE END LOOP OVER EVENTS

    std::cout << "Max bar count: " << maxbarcount << std::endl;
    gStyle->SetOptStat(0);
    // make ratio histograms to first variation
    TH2F* acceptanceModuleWidthRatio[20];
    for (int i = 0; i < 20; i++)
    {
        acceptanceModuleWidthRatio[i] = (TH2F*)acceptanceModuleWidth[i]->Clone(Form("acceptanceModuleWidthRatio%d", i));
        acceptanceModuleWidthRatio[i]->Divide(acceptanceModuleWidth[0]);
        acceptanceModuleWidthRatio[i]->SetTitle(Form("Acceptance ratio for %d less bars to default", 2*i));
    }

    // plot all histograms on a 4x5 canvas
    TCanvas *c1 = new TCanvas("c1", "c1", 1600, 1200);
    c1->Divide(4, 5);
    for (int i = 0; i < 20; i++)
    {
        c1->cd(i + 1);
        c1->SetLogz();
        acceptanceModuleWidth[i]->GetXaxis()->SetTitle("1/p_{T} [GeV^{-1}]");
        acceptanceModuleWidth[i]->GetYaxis()->SetTitle("#eta");
        acceptanceModuleWidth[i]->Draw("colz");
    }
    c1->SaveAs(outputdir + "acceptanceModuleWidth.pdf");

    TCanvas *c3 = new TCanvas("c3", "c3", 1600, 1200);
    c3->Divide(4, 5);
    for (int i = 0; i < 20; i++)
    {
        c3->cd(i + 1);
        acceptanceModuleWidthRatio[i]->GetXaxis()->SetTitle("1/p_{T} [GeV^{-1}]");
        acceptanceModuleWidthRatio[i]->GetYaxis()->SetTitle("#eta");
        acceptanceModuleWidthRatio[i]->Draw("colz");
    }
    c3->SaveAs(outputdir + "acceptanceModuleWidthRatio.pdf");

    TCanvas *c4 = new TCanvas("c4", "c4", 1600, 1200);
    c4->Divide(4, 5);
    for (int i = 0; i < 20; i++)
    {
        c4->cd(i + 1);
        hHitsZXMSModule0[i]->Draw("colz");
    }
    c4->SaveAs(outputdir + "hHitsZXMSModule0.pdf");


    TCanvas *c2 = new TCanvas("c2", "c2", 1600, 1200);
    c2->Divide(4, 5);
    for (int i = 0; i < 20; i++)
    {
        c2->cd(i + 1);
        hHitsZXMS[i]->Draw("colz");
    }
    c2->SaveAs(outputdir + "hHitsZXMS.pdf");

    fin.Close();
}
