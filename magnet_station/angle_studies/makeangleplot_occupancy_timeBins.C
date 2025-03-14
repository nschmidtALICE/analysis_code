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
void makeangleplot_occupancy_timeBins(TString inputfile = "/home/niviths/Downloads/magnetStationSims/20250311_PbPb_addUTinfo/20250311_PbPb_addUTinfo.root") // tilt
// void makeangleplot_occupancy_timeBins(TString inputfile = "/home/niviths/Downloads/magnetStationSims/20250311_pp_newOutput/20250311_pp_newOutput.root") // tilt
{
    TFile fin(inputfile, "READ");
    cout << "Opened file " << inputfile << endl;
    TTreeReader tree("ntup", &fin);
    double modSpacing = 300 * TMath::Cos(0.42);
    double minZ = 3500 + 2 * 500 + 300;
    double maxZ = 7700;

    TString outputdir = "occupancy_plots/";
    //make output directory
    system("mkdir -p " + outputdir);

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
    TTreeReaderArray<int>   ievt(tree, "ievt");
    const int nModules = 9;
    

    // TH2F *hYZHitsinMSSegments = new TH2F("hYZHitsinMSSegments", "YZ of hits in plane", 125, minZ, maxZ, 100, -1500, 1500);
    float maxNSegmentHit = 15000;
    TH1F *hNSegmentsHit = new TH1F("hNSegmentsHit", "Number of segments hit per event", 100, 0, maxNSegmentHit);
    TH1F* hNDoublySegmentsHit = new TH1F("hNDoublySegmentsHit", "Number of segments hit by multiple particles", 35, 0, 35);
    TH1F* hFractionDoublySegmentsHit = new TH1F("hFractionDoublySegmentsHit", "Fraction of segments hit by multiple particles", 80, 0, 20);

    double binningYOccupancyPlot[9] = {0, 3, 9, 18, 30, 45, 63, 92, 100};
    TH2F* hOccupancyPlotStation[4];
    TH2F* hOccupancyPlotStationTimeBins[4][10];
    for(int i = 0; i < 4; i++){
        hOccupancyPlotStation[i] = new TH2F(Form("hOccupancyPlotStation%d",i), Form("Occupancy plot for station %d",i), 57*7+2*31, 0, 57*7+2*31, 8, binningYOccupancyPlot);
        for(int j = 0; j < 10; j++){
            hOccupancyPlotStationTimeBins[i][j] = new TH2F(Form("hOccupancyPlotStation%dTimeBin%d",i,j), Form("Occupancy plot for station %d time bin %d",i,j), 57*7+2*31, 0, 57*7+2*31, 8, binningYOccupancyPlot);
        }
    }

    //Define a logarithmic binning from 0.1 to 10 GeV with 100 bins
    // const int nBins = 100;
    // double logBins[nBins + 1];
    // for (int i = 0; i <= nBins; i++)
    // {
    //     logBins[i] = pow(10, -1 + i * 0.1);
    // }

    TH2F* hMStimeVsEnergyPion = new TH2F("hMStimeVsEnergyPion", "MS time vs energy for pions", 100, 0.4, 6, 100, 0, 70);
    TH2F* hMStimeVsEnergyProton = new TH2F("hMStimeVsEnergyProton", "MS time vs energy for protons", 100, 0.4, 6, 100, 0, 70);
    TH2F* hMStimeVsEnergyElectron = new TH2F("hMStimeVsEnergyElectron", "MS time vs energy for electrons", 100, 0.4, 6, 100, 0, 70);
    TH2F* hMStimeVsEnergyKaon = new TH2F("hMStimeVsEnergyKaon", "MS time vs energy for kaons", 100, 0.4, 6, 100, 0, 70);


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

    //create vector to store the segments that are hit
    int TimeBinCounter[26] = {0};
    std::vector<std::pair<int, int>> segmentsHitTimeBins;
    std::vector<std::pair<int, int>> doublysegmentsHitTimeBins;
    std::vector<std::pair<int, int>> doublysegmentsHitKeysTimeBins;


    bool fillOnce = 1;
    int nHits_MSbars = 0;
    int currentEventNumber = 0;
    cout << "Number of entries: " << tree.GetEntries() << endl;
    int nUpdateEvents = 1000;
    if(tree.GetEntries() < 1000){
        nUpdateEvents = 10;
    }
    while (tree.Next())
    {
    
        curr_entry++;
        if (curr_entry % nUpdateEvents == 0)
            cout << "Entry: " << curr_entry << endl;
            if(fillOnce==0){
                break;
            }   
        // cout <<  pt.GetSize() << endl;
        // if( pt.GetSize()<25000) continue;
        for (auto i = 0U; i < pt.GetSize(); ++i)
        {
            if(fillOnce==0){
                break;
            }
            if(ievt[i] != currentEventNumber){
                hNSegmentsHit->Fill(segmentsHitTimeBins.size());
                // if(doublysegmentsHit.size()>0)
                hNDoublySegmentsHit->Fill(doublysegmentsHitTimeBins.size());
                if(segmentsHitTimeBins.size() > 0){
                    hFractionDoublySegmentsHit->Fill(100*(double)doublysegmentsHitTimeBins.size()/(double)segmentsHitTimeBins.size());
                        //find the time bin with the most entries from the TimeBinCounter vector
                        int maxTimeBin = 0;
                        int maxTimeBinCount = 0;
                        for(auto j = 0; j < 26; ++j){
                            int count = TimeBinCounter[j];
                            if(count > maxTimeBinCount){
                                maxTimeBinCount = count;
                                maxTimeBin = j;
                            }
                            // cout << "Time bin: " << j << " with " << count << " entries" << endl;
                        }
                        // cout << "Max time bin: " << maxTimeBin << " with " << maxTimeBinCount << " entries" << endl;

                    //loop over segmentsHit and fill the occupancy plot
                    // cout << "n segments hit: " << segmentsHit.size() << endl;
                    if(fillOnce && segmentsHitTimeBins.size() > 5000){

                        for(auto j = 0U; j < segmentsHitTimeBins.size(); ++j){
                            if(segmentsHitTimeBins[j].second == 3 || segmentsHitTimeBins[j].second == 4 || segmentsHitTimeBins[j].second == 5){
                                int station = (segmentsHitTimeBins[j].first >> 8) & 0x7;  // 3 bit
                                // if(station != 0){
                                //     continue;
                                // }
                                int module = (segmentsHitTimeBins[j].first >> 11) & 0xF;  // 4 bit
                                // int layer = (segmentsHitTimeBins[j].first >> 15) & 0x7;   // 3 bit
                                int segment = (segmentsHitTimeBins[j].first >> 18) & 0xF; // 4 bit
                                int bar = (segmentsHitTimeBins[j].first >> 22) & 0x3F;    // 6 bit
                                // int isfiber = (segmentsHitTimeBins[j].first >> 28) & 0x3; // 2 bit
                                // int issupp = (segmentsHitTimeBins[j].first >> 30) & 0x3; // 2 bit
                                int nBinXFill = 0;
                                if(module<7)
                                    nBinXFill = 57*module + bar;
                                else
                                    nBinXFill = 57*7 + 31*(module-7) + bar;
                                hOccupancyPlotStation[station]->Fill(nBinXFill, binningYOccupancyPlot[segment]);
                                // cout << "Station: " << station << " Module: " << module << " Segment: " << segment << " Bar: " << bar << endl;
                            }
                        }
                        // fillOnce = 0;
                    }
                }
                currentEventNumber = ievt[i];
                
                for(int j = 0; j < 26; ++j){
                    TimeBinCounter[j] = 0;
                }
                segmentsHitTimeBins.clear();
                doublysegmentsHitTimeBins.clear();
                doublysegmentsHitKeysTimeBins.clear();

            }
            // if (p[i] > 5000)
            // {
            //     // continue;
            // }
            if(eta[i] < 1.8 || eta[i] > 5.0){
                continue;
            }

            // loop over size of ms_px to get all the hits per event
            for (auto j = 0U; j < ms_px.GetSize(); ++j)
            {
                if (ms_id[j] != key[i])
                    continue;
                // require hits within MS acceptance
                if (ms_vz[j] < 3000 || ms_vz[j] > maxZ)
                {
                    continue;
                }

                //fill histograms for time vs energy for different particle types
                if(pid[i] == 211){
                    // cout << "pion momentum: " << p[i] << " ms_time: " << ms_time[j] << endl;
                    hMStimeVsEnergyPion->Fill(p[i]/1000, ms_time[j]);
                } else if(pid[i] == 2212){
                    hMStimeVsEnergyProton->Fill(p[i]/1000, ms_time[j]);
                } else if(pid[i] == 11){
                    hMStimeVsEnergyElectron->Fill(p[i]/1000, ms_time[j]);
                } else if(pid[i] == 321){
                    hMStimeVsEnergyKaon->Fill(p[i]/1000, ms_time[j]);
                }

                // int station = (ms_bitID[j] >> 8) & 0x7;  // 3 bit
                // int module = (ms_bitID[j] >> 11) & 0xF;  // 4 bit
                // int layer = (ms_bitID[j] >> 15) & 0x7;   // 3 bit
                // int segment = (ms_bitID[j] >> 18) & 0xF; // 4 bit
                // int bar = (ms_bitID[j] >> 22) & 0x3F;    // 6 bit
                int isfiber = (ms_bitID[j] >> 28) & 0x3; // 2 bit
                int issupp = (ms_bitID[j] >> 30) & 0x3; // 2 bit
                
                int msTimeBin =  (int)(std::fmod(ms_time[j] - 17.0,25));
                // cout << "ms time: " << ms_time[j] << " ms time bin: " << msTimeBin << endl;

                if (!isfiber && !issupp){
                    //if ms_bitID is not yet in segmentsHit for a given time bin, add it
                    if(std::find(segmentsHitTimeBins.begin(), segmentsHitTimeBins.end(), std::make_pair(ms_bitID[j], msTimeBin)) == segmentsHitTimeBins.end()){
                        segmentsHitTimeBins.push_back(std::make_pair(ms_bitID[j], msTimeBin));
                        TimeBinCounter[msTimeBin]++;
                    } else {
                        //if key is not yet in doublysegmentsHitKeys for a given time bin, add it
                        if(std::find(doublysegmentsHitKeysTimeBins.begin(), doublysegmentsHitKeysTimeBins.end(), std::make_pair(key[i], msTimeBin)) == doublysegmentsHitKeysTimeBins.end()){
                            doublysegmentsHitTimeBins.push_back(std::make_pair(ms_bitID[j], msTimeBin));
                            doublysegmentsHitKeysTimeBins.push_back(std::make_pair(key[i], msTimeBin));
                        }
                    }
                }
            } // NOTE END LOOP OVER MS HITS
        } // NOTE END LOOP OVER PARTICLES
    } // NOTE END LOOP OVER EVENTS

    //plot hNSegmentsHit
    TCanvas *c8 = new TCanvas("c8", "c8", 800, 600);
    gStyle->SetOptStat(0);
    hNSegmentsHit->Scale(100.0/(curr_entry));
    hNSegmentsHit->GetXaxis()->SetTitle("Number of segments hit");
    hNSegmentsHit->GetYaxis()->SetTitle("Fraction of events [%]");
    hNSegmentsHit->Draw();
    c8->SaveAs(Form("%shNSegmentsHit.pdf", outputdir.Data()));

    //plot hNDoublySegmentsHit
    TCanvas *c9 = new TCanvas("c9", "c9", 800, 600);
    gStyle->SetOptStat(0);
    hNDoublySegmentsHit->Scale(100.0/curr_entry);
    hNDoublySegmentsHit->GetXaxis()->SetTitle("Number of segments hit by multiple particles");
    hNDoublySegmentsHit->GetYaxis()->SetTitle("Fraction of events [%]");
    hNDoublySegmentsHit->Draw();
    c9->SaveAs(Form("%shNDoublySegmentsHit.pdf", outputdir.Data()));

    //plot hFractionDoublySegmentsHit
    TCanvas *c10 = new TCanvas("c10", "c10", 800, 600);
    gStyle->SetOptStat(0);
    hFractionDoublySegmentsHit->Scale(100.0/curr_entry);
    hFractionDoublySegmentsHit->GetXaxis()->SetTitle("Fraction of segments hit by multiple particles [%]");
    hFractionDoublySegmentsHit->GetYaxis()->SetTitle("Fraction of events [%]");
    hFractionDoublySegmentsHit->Draw();
    c10->SaveAs(Form("%shFractionDoublySegmentsHit.pdf", outputdir.Data()));

    //plot hMStimeVsEnergyPion
    TCanvas *c17 = new TCanvas("c17", "c17", 800, 600);
    gStyle->SetOptStat(0);
    hMStimeVsEnergyPion->GetXaxis()->SetTitle("p [GeV]");
    hMStimeVsEnergyPion->GetYaxis()->SetTitle("MS time [ns]");
    hMStimeVsEnergyPion->Draw("colz");
    c17->SaveAs(Form("%shMStimeVsEnergyPion.pdf", outputdir.Data()));

    //plot hMStimeVsEnergyProton
    TCanvas *c18 = new TCanvas("c18", "c18", 800, 600);
    gStyle->SetOptStat(0);
    hMStimeVsEnergyProton->GetXaxis()->SetTitle("p [GeV]");
    hMStimeVsEnergyProton->GetYaxis()->SetTitle("MS time [ns]");
    hMStimeVsEnergyProton->Draw("colz");
    c18->SaveAs(Form("%shMStimeVsEnergyProton.pdf", outputdir.Data()));

    //plot hMStimeVsEnergyElectron
    TCanvas *c19 = new TCanvas("c19", "c19", 800, 600);
    gStyle->SetOptStat(0);
    hMStimeVsEnergyElectron->GetXaxis()->SetTitle("p [GeV]");
    hMStimeVsEnergyElectron->GetYaxis()->SetTitle("MS time [ns]");
    hMStimeVsEnergyElectron->Draw("colz");
    c19->SaveAs(Form("%shMStimeVsEnergyElectron.pdf", outputdir.Data()));

    //plot hMStimeVsEnergyKaon
    TCanvas *c20 = new TCanvas("c20", "c20", 800, 600);
    gStyle->SetOptStat(0);
    hMStimeVsEnergyKaon->GetXaxis()->SetTitle("p [GeV]");
    hMStimeVsEnergyKaon->GetYaxis()->SetTitle("MS time [ns]");
    hMStimeVsEnergyKaon->Draw("colz");
    c20->SaveAs(Form("%shMStimeVsEnergyKaon.pdf", outputdir.Data()));

    //plot hOccupancyPlotStation
    for(int i = 0; i < 4; i++){
        TCanvas *c21 = new TCanvas(Form("c21_%d",i), Form("c21_%d",i), 800, 600);
        gStyle->SetOptStat(0);
        hOccupancyPlotStation[i]->GetXaxis()->SetTitle("Module and bar number");
        hOccupancyPlotStation[i]->GetYaxis()->SetTitle("Segment number");
        hOccupancyPlotStation[i]->Draw("colz");
        c21->SaveAs(Form("%shOccupancyPlotStation%d.pdf", outputdir.Data(), i));
    }

    //plot hOccupancyPlotStation on a 2x2 grid
    TCanvas *c22 = new TCanvas("c22", "c22", 1200, 800);
    c22->Divide(2, 2);
    for(int i = 0; i < 4; i++){
        c22->cd(i+1);
        gStyle->SetOptStat(0);
        hOccupancyPlotStation[i]->GetXaxis()->SetTitle("Module and bar number");
        hOccupancyPlotStation[i]->GetYaxis()->SetTitle("y-position in panel (segmented in bars)");
        hOccupancyPlotStation[i]->Draw("colz");
    }
    c22->SaveAs(Form("%shOccupancyPlotStationGrid.pdf", outputdir.Data()));

    fin.Close();

}



// if(std::find(segmentsHit.begin(), segmentsHit.end(), ms_bitID[j]) == segmentsHit.end()){
//     segmentsHit.push_back(ms_bitID[j]);
// } else {
//     //if key is not yet in doublysegmentsHitKeys, add it
//     if(std::find(doublysegmentsHitKeys.begin(), doublysegmentsHitKeys.end(), key[i]) == doublysegmentsHitKeys.end()){
//         doublysegmentsHit.push_back(ms_bitID[j]);
//         doublysegmentsHitKeys.push_back(key[i]);
//     }
// }