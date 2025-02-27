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

// void makeangleplot_occupancy(TString inputfile = "/home/niviths/Downloads/magnetStationSims/2100_wtf/minimumBias_MS_MagDown_2100plus.root") // tilt
// void makeangleplot_occupancy(TString inputfile = "/home/niviths/Downloads/magnetStationSims/minimumBias_MS_MagDown_PbPb.root") // tilt
// void makeangleplot_occupancy(TString inputfile = "/home/niviths/Downloads/magnetStationSims/20250216/minimumBias_MS_MagDown_20250216.root") // tilt
void makeangleplot_occupancy(TString inputfile = "/home/niviths/Downloads/magnetStationSims/20250218_8cm_inward_lastPanelsNarrow/minimumBias_MS_MagDown_20250218.root") // tilt
// void makeangleplot_occupancy(TString inputfile = "/home/niviths/Downloads/magnetStationSims/minimumBias_MS_MagDown_10015plus.root") // tilt
{
    TFile fin(inputfile, "READ");
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
    TH1F *hNSegmentsHitTimeCut = new TH1F("hNSegmentsHitTimeCut", "Number of segments hit per event with time cut", 100, 0, maxNSegmentHit);
    TH1F *hNSegmentsHitPrimary = new TH1F("hNSegmentsHitPrimary", "Number of segments hit per event for primary particles", 100, 0, maxNSegmentHit);
    TH1F* hNDoublySegmentsHit = new TH1F("hNDoublySegmentsHit", "Number of segments hit by multiple particles", 35, 0, 35);
    TH1F* hNDoublySegmentsHitTimeCut = new TH1F("hNDoublySegmentsHitTimeCut", "Number of segments hit by multiple particles with timing cut", 35, 0, 35);
    TH1F* hNDoublySegmentsHitPrimary = new TH1F("hNDoublySegmentsHitPrimary", "Number of segments hit by multiple primary particles", 35, 0, 35);
    TH1F* hFractionDoublySegmentsHit = new TH1F("hFractionDoublySegmentsHit", "Fraction of segments hit by multiple particles", 80, 0, 20);
    TH1F* hFractionDoublySegmentsHitTimeCut = new TH1F("hFractionDoublySegmentsHitTimeCut", "Fraction of segments hit by multiple particles with timing cut", 80, 0, 20);
    TH1F* hFractionDoublySegmentsHitPrimary = new TH1F("hFractionDoublySegmentsHitPrimary", "Fraction of segments hit by multiple primary particles", 80, 0, 20);

    double binningYOccupancyPlot[9] = {0, 3, 9, 18, 30, 45, 63, 92, 100};
    TH2F* hOccupancyPlotStation[4];
    for(int i = 0; i < 4; i++){
        hOccupancyPlotStation[i] = new TH2F(Form("hOccupancyPlotStation%d",i), Form("Occupancy plot for station %d",i), 57*7+2*31, 0, 57*7+2*31, 8, binningYOccupancyPlot);
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
    std::vector<int> segmentsHit;
    std::vector<int> doublysegmentsHit;
    std::vector<int> doublysegmentsHitKeys;
    std::vector<int> segmentsHitPrimary;
    std::vector<int> doublysegmentsHitPrimary;
    std::vector<int> doublysegmentsHitKeysPrimary;
    std::vector<int> segmentsHitTimeCut;
    std::vector<int> doublysegmentsHitTimeCut;
    std::vector<int> doublysegmentsHitKeysTimeCut;

    bool fillOnce = 1;
    int nHits_MSbars = 0;
    int currentEventNumber = 0;
    cout << "Number of entries: " << tree.GetEntries() << endl;
    while (tree.Next())
    {
    
        curr_entry++;
        if (curr_entry % 10 == 0)
            cout << "Entry: " << curr_entry << endl;
            // if(fillOnce==0){
            //     break;
            // }   
        cout <<  pt.GetSize() << endl;
        // if( pt.GetSize()<25000) continue;
        for (auto i = 0U; i < pt.GetSize(); ++i)
        {
            // if(fillOnce==0){
            //     break;
            // }
            if(ievt[i] != currentEventNumber){
                // if(segmentsHit.size()>0)
                hNSegmentsHit->Fill(segmentsHit.size());
                // if(doublysegmentsHit.size()>0)
                hNDoublySegmentsHit->Fill(doublysegmentsHit.size());
                if(segmentsHit.size() > 0){
                    hFractionDoublySegmentsHit->Fill(100*(double)doublysegmentsHit.size()/(double)segmentsHit.size());
                    //loop over segmentsHit and fill the occupancy plot
                    // cout << "n segments hit: " << segmentsHit.size() << endl;
                    if(fillOnce && segmentsHit.size() > 5000){
                        for(auto j = 0U; j < segmentsHit.size(); ++j){
                            int station = (segmentsHit[j] >> 8) & 0x7;  // 3 bit
                            // if(station != 0){
                            //     continue;
                            // }
                            int module = (segmentsHit[j] >> 11) & 0xF;  // 4 bit
                            // int layer = (segmentsHit[j] >> 15) & 0x7;   // 3 bit
                            int segment = (segmentsHit[j] >> 18) & 0xF; // 4 bit
                            int bar = (segmentsHit[j] >> 22) & 0x3F;    // 6 bit
                            // int isfiber = (segmentsHit[j] >> 28) & 0x3; // 2 bit
                            // int issupp = (segmentsHit[j] >> 30) & 0x3; // 2 bit
                            int nBinXFill = 0;
                            if(module<7)
                                nBinXFill = 57*module + bar;
                            else
                                nBinXFill = 57*7 + 31*(module-7) + bar;
                            hOccupancyPlotStation[station]->Fill(nBinXFill, binningYOccupancyPlot[segment]);
                            // cout << "Station: " << station << " Module: " << module << " Segment: " << segment << " Bar: " << bar << endl;
                        }
                    }
                    fillOnce = 0;
                }
                currentEventNumber = ievt[i];

                

                segmentsHit.clear();
                doublysegmentsHit.clear();
                doublysegmentsHitKeys.clear();

                //do the same for segments with time cut
                hNSegmentsHitTimeCut->Fill(segmentsHitTimeCut.size());
                hNDoublySegmentsHitTimeCut->Fill(doublysegmentsHitTimeCut.size());
                if(segmentsHitTimeCut.size() > 0)
                    hFractionDoublySegmentsHitTimeCut->Fill(100*(double)doublysegmentsHitTimeCut.size()/(double)segmentsHitTimeCut.size());
                segmentsHitTimeCut.clear();
                doublysegmentsHitTimeCut.clear();
                doublysegmentsHitKeysTimeCut.clear();

                //do the same for primary particles
                hNSegmentsHitPrimary->Fill(segmentsHitPrimary.size());
                hNDoublySegmentsHitPrimary->Fill(doublysegmentsHitPrimary.size());
                if(segmentsHitPrimary.size() > 0)
                    hFractionDoublySegmentsHitPrimary->Fill(100*(double)doublysegmentsHitPrimary.size()/(double)segmentsHitPrimary.size());
                segmentsHitPrimary.clear();
                doublysegmentsHitPrimary.clear();
                doublysegmentsHitKeysPrimary.clear();
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

                if (!isfiber && !issupp){
                    //if ms_bitID is not yet in segmentsHit, add it
                    if(std::find(segmentsHit.begin(), segmentsHit.end(), ms_bitID[j]) == segmentsHit.end()){
                        segmentsHit.push_back(ms_bitID[j]);
                    } else {
                        //if key is not yet in doublysegmentsHitKeys, add it
                        if(std::find(doublysegmentsHitKeys.begin(), doublysegmentsHitKeys.end(), key[i]) == doublysegmentsHitKeys.end()){
                            doublysegmentsHit.push_back(ms_bitID[j]);
                            doublysegmentsHitKeys.push_back(key[i]);
                        }
                    }

                    //do the same but require a time cut on the ms hits
                    if(ms_time[j] > 17 && ms_time[j] < 37){
                        if(std::find(segmentsHitTimeCut.begin(), segmentsHitTimeCut.end(), ms_bitID[j]) == segmentsHitTimeCut.end()){
                            segmentsHitTimeCut.push_back(ms_bitID[j]);
                        } else {
                            if(std::find(doublysegmentsHitKeysTimeCut.begin(), doublysegmentsHitKeysTimeCut.end(), key[i]) == doublysegmentsHitKeysTimeCut.end()){
                                doublysegmentsHitTimeCut.push_back(ms_bitID[j]);
                                doublysegmentsHitKeysTimeCut.push_back(key[i]);
                            }
                        }
                    }

                    //do the same but only for primary particles (particles with mother key < -1000)
                    if(parent_key[i] < -1000){
                        if(std::find(segmentsHitPrimary.begin(), segmentsHitPrimary.end(), ms_bitID[j]) == segmentsHitPrimary.end()){
                            segmentsHitPrimary.push_back(ms_bitID[j]);
                        } else {
                            if(std::find(doublysegmentsHitKeysPrimary.begin(), doublysegmentsHitKeysPrimary.end(), key[i]) == doublysegmentsHitKeysPrimary.end()){
                                doublysegmentsHitPrimary.push_back(ms_bitID[j]);
                                doublysegmentsHitKeysPrimary.push_back(key[i]);
                            }
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

    //plot hNSegmentsHitTimeCut
    TCanvas *c14 = new TCanvas("c14", "c14", 800, 600);
    gStyle->SetOptStat(0);
    hNSegmentsHitTimeCut->Scale(100.0/(curr_entry));
    hNSegmentsHitTimeCut->GetXaxis()->SetTitle("Number of segments hit with time cut");
    hNSegmentsHitTimeCut->GetYaxis()->SetTitle("Fraction of events [%]");
    hNSegmentsHitTimeCut->Draw();
    c14->SaveAs(Form("%shNSegmentsHitTimeCut.pdf", outputdir.Data()));

    //plot hNSegmentsHitPrimary
    TCanvas *c11 = new TCanvas("c11", "c11", 800, 600);
    gStyle->SetOptStat(0);
    hNSegmentsHitPrimary->Scale(100.0/(curr_entry));
    hNSegmentsHitPrimary->GetXaxis()->SetTitle("Number of segments hit by primary particles");
    hNSegmentsHitPrimary->GetYaxis()->SetTitle("Fraction of events [%]");
    hNSegmentsHitPrimary->Draw();
    c11->SaveAs(Form("%shNSegmentsHitPrimary.pdf", outputdir.Data()));

    //plot hNDoublySegmentsHit
    TCanvas *c9 = new TCanvas("c9", "c9", 800, 600);
    gStyle->SetOptStat(0);
    hNDoublySegmentsHit->Scale(100.0/curr_entry);
    hNDoublySegmentsHit->GetXaxis()->SetTitle("Number of segments hit by multiple particles");
    hNDoublySegmentsHit->GetYaxis()->SetTitle("Fraction of events [%]");
    hNDoublySegmentsHit->Draw();
    c9->SaveAs(Form("%shNDoublySegmentsHit.pdf", outputdir.Data()));

    //plot hNDoublySegmentsHitTimeCut
    TCanvas *c15 = new TCanvas("c15", "c15", 800, 600);
    gStyle->SetOptStat(0);
    hNDoublySegmentsHitTimeCut->Scale(100.0/curr_entry);
    hNDoublySegmentsHitTimeCut->GetXaxis()->SetTitle("Number of segments hit by multiple particles with timing cut");
    hNDoublySegmentsHitTimeCut->GetYaxis()->SetTitle("Fraction of events [%]");
    hNDoublySegmentsHitTimeCut->Draw();
    c15->SaveAs(Form("%shNDoublySegmentsHitTimeCut.pdf", outputdir.Data()));

    //plot hNDoublySegmentsHitPrimary
    TCanvas *c12 = new TCanvas("c12", "c12", 800, 600);
    gStyle->SetOptStat(0);
    hNDoublySegmentsHitPrimary->Scale(100.0/curr_entry);
    hNDoublySegmentsHitPrimary->GetXaxis()->SetTitle("Number of segments hit by multiple primary particles");
    hNDoublySegmentsHitPrimary->GetYaxis()->SetTitle("Fraction of events [%]");
    hNDoublySegmentsHitPrimary->Draw();
    c12->SaveAs(Form("%shNDoublySegmentsHitPrimary.pdf", outputdir.Data()));

    //plot hFractionDoublySegmentsHit
    TCanvas *c10 = new TCanvas("c10", "c10", 800, 600);
    gStyle->SetOptStat(0);
    hFractionDoublySegmentsHit->Scale(100.0/curr_entry);
    hFractionDoublySegmentsHit->GetXaxis()->SetTitle("Fraction of segments hit by multiple particles [%]");
    hFractionDoublySegmentsHit->GetYaxis()->SetTitle("Fraction of events [%]");
    hFractionDoublySegmentsHit->Draw();
    c10->SaveAs(Form("%shFractionDoublySegmentsHit.pdf", outputdir.Data()));

    //plot hFractionDoublySegmentsHitTimeCut
    TCanvas *c16 = new TCanvas("c16", "c16", 800, 600);
    gStyle->SetOptStat(0);
    hFractionDoublySegmentsHitTimeCut->Scale(100.0/curr_entry);
    hFractionDoublySegmentsHitTimeCut->GetXaxis()->SetTitle("Fraction of segments hit by multiple particles with timing cut [%]");
    hFractionDoublySegmentsHitTimeCut->GetYaxis()->SetTitle("Fraction of events [%]");
    hFractionDoublySegmentsHitTimeCut->Draw();
    c16->SaveAs(Form("%shFractionDoublySegmentsHitTimeCut.pdf", outputdir.Data()));

    //plot hFractionDoublySegmentsHitPrimary
    TCanvas *c13 = new TCanvas("c13", "c13", 800, 600);
    gStyle->SetOptStat(0);
    hFractionDoublySegmentsHitPrimary->Scale(100.0/curr_entry);
    hFractionDoublySegmentsHitPrimary->GetXaxis()->SetTitle("Fraction of segments hit by multiple primary particles [%]");
    hFractionDoublySegmentsHitPrimary->GetYaxis()->SetTitle("Fraction of events [%]");
    hFractionDoublySegmentsHitPrimary->Draw();
    c13->SaveAs(Form("%shFractionDoublySegmentsHitPrimary.pdf", outputdir.Data()));

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