#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>

std::vector<std::vector<int>> vec_segment_combinations = {
    {0, 0, 0, 0},
    {1, 1, 1, 0},
    {1, 1, 1, 1},
    {2, 2, 2, 1},
    {2, 2, 2, 2},
    {2, 3, 2, 2},
    {3, 3, 3, 2},
    {3, 3, 3, 3},
    {3, 4, 3, 3},
    {3, 4, 4, 3},
    {4, 4, 4, 3},
    {4, 4, 4, 4},
    {4, 5, 4, 4},
    {4, 5, 5, 4},
    {5, 5, 5, 4},
    {5, 5, 5, 5},
    {5, 6, 5, 5},
    {5, 6, 6, 5},
    {6, 6, 6, 5},
    {6, 6, 6, 6},
    {6, 7, 6, 6},
    {6, 7, 7, 6},
    {7, 7, 7, 6}};

class Point
{
public:
    double x, y, z, t;
    Point(double x, double y, double z, double t = 0.0) : x(x), y(y), z(z), t(t) {}
};

int track_ambiguity_finder()
{
    // bool verboseoutput = true;
    bool verboseoutput = false;
    // const char *inputfile = "/home/niviths/Downloads/magnetStationSims/20250218_8cm_inward_lastPanelsNarrow/minimumBias_MS_MagDown_20250218_8cm_inward_lastPanelsNarrow.root";
    const char *inputfile = "/home/niviths/Downloads/magnetStationSims/20250227_pp_LayerNumbersFixed/20250227_pp_LayerNumbersFixed.root";
    TFile fin(inputfile);
    TTree *ntup = (TTree *)fin.Get("ntup");
    std::cout << ntup->GetEntries() << std::endl;

    TString outputdir = "extraploation_plots/";
    // make output directory
    system("mkdir -p " + outputdir);

    TTreeReader tree(ntup);
    // TTreeReaderArray<float> ut_vx(tree, "ut_vx");
    // TTreeReaderArray<float> ut_vy(tree, "ut_vy");
    // TTreeReaderArray<float> ut_vz(tree, "ut_vz");
    // TTreeReaderArray<float> ut_tx(tree, "ut_tx");
    // TTreeReaderArray<float> ut_ty(tree, "ut_ty");
    TTreeReaderArray<float> p(tree, "p");
    TTreeReaderArray<int> pid(tree, "pid");
    TTreeReaderArray<int> key(tree, "key");
    TTreeReaderArray<float> ms_vx(tree, "ms_vx");
    TTreeReaderArray<float> ms_vy(tree, "ms_vy");
    TTreeReaderArray<float> ms_vz(tree, "ms_vz");
    TTreeReaderArray<float> ms_time(tree, "ms_time");
    TTreeReaderArray<int> ms_id(tree, "ms_id");
    TTreeReaderArray<int> ms_bitID(tree, "ms_segment");
    TTreeReaderArray<int> nUThits(tree, "nUThits");

    std::vector<int> clusterizedHits_bitID;
    std::vector<float> clusterizedHits_time;
    std::vector<int> clusterizedHits_id;

    TH2D* hRelativeBarDifferenceToFirstLayer = new TH2D("hRelativeBarDifferenceToFirstLayer", "Relative bar difference to first layer (truth)", 3, 0.5, 3.5, 21, -10.5, 10.5);
    TH1D* hRelativeBarDifferenceToFirstLayerIndiv[3];
    for (int i = 0; i < 3; i++)
    {
        hRelativeBarDifferenceToFirstLayerIndiv[i] = new TH1D(Form("hRelativeBarDifferenceToFirstLayerIndiv_%d", i+1), Form("Relative bar difference to first layer (truth) for layer %d", i+1), 21, -10.5, 10.5);
    }
    TH1D* hNGroups = new TH1D("hNGroups", "Number of groups", 10, -0.50, 9.5);
    int numevt = 0;
    while (tree.Next())
    {
        // if (numevt > 5) break;
        std::cout << "evt " << numevt << std::endl;

        if (numevt>1000) break; //NOTE break condition

        for (size_t i = 0; i < pid.GetSize(); ++i)
        {
            clusterizedHits_bitID.clear();
            clusterizedHits_time.clear();
            clusterizedHits_id.clear();
            if (p[i] > 5000)
            {
                // std::cout << "skipping track with p > 5000" << std::endl;
                continue;
            }
            if (nUThits[i] < 4)
            {
                // std::cout << "skipping track with nUThits < 5" << std::endl;
                continue;
            }
            int segment_match = -1;
            int bar_match = -1;
            Point point_match(0, 0, 0);
            float time_match = 0;
            int bitID_match = 0;
            // find the hit in MS with the lowest time and same key
            for (size_t j = 0; j < ms_vz.GetSize(); ++j)
            {
                if (ms_vz[j] < 3000 || ms_vz[j] > 7700)
                    continue;
                if (ms_bitID[j] >> 28 & 0x3 || ms_bitID[j] >> 30 & 0x3)
                    continue;

                //TODO ignoring station 1 and station 3 for now
                // if (((ms_bitID[j] >> 8) & 0x7) == 1 || ((ms_bitID[j] >> 8) & 0x7) == 3)
                // if (((ms_bitID[j] >> 8) & 0x7) == 0 || ((ms_bitID[j] >> 8) & 0x7) == 2)
                //     continue;

                Point ms_pos(ms_vx[j], ms_vy[j], ms_vz[j]);
                if (ms_id[j] == key[i])
                {
                    if (segment_match == -1 || ms_time[j] < time_match)
                    {
                        segment_match = ms_bitID[j] >> 18 & 0xF;
                        bar_match = ms_bitID[j] >> 22 & 0x3F;
                        point_match = ms_pos;
                        time_match = ms_time[j];
                        bitID_match = ms_bitID[j];
                    }
                }
                if (std::find(clusterizedHits_bitID.begin(), clusterizedHits_bitID.end(), ms_bitID[j]) == clusterizedHits_bitID.end())
                {
                    clusterizedHits_bitID.push_back(ms_bitID[j]);
                    clusterizedHits_time.push_back(ms_time[j]);
                    clusterizedHits_id.push_back(ms_id[j]);
                }
            }
            if (segment_match == -1)
            {
                if(verboseoutput)std::cout << "no matching segment found" << std::endl;
                continue;
            }
            else
            {
                if(verboseoutput)std::cout << "\tmatching segment found: " << segment_match << " layer: " << ((bitID_match >> 15) & 0x7) << " bar: " << bar_match << " module: " << ((bitID_match >> 11) & 0xF) << " station: " << ((bitID_match >> 8) & 0x7) << " time: " << time_match << std::endl;
            }

            int stationmatch = (bitID_match >> 8) & 0x7;
            int modulematch = (bitID_match >> 11) & 0xF;
            // loop again to find close by hits in MS to point_match within the same segment and similar time window, hits in at least three layers are needed
            // std::vector<Point> close_by_hits;
            std::vector<int> close_by_hits_layer;
            std::vector<int> close_by_hits_bitID;
            std::vector<float> close_by_hits_time;
            std::vector<bool> close_by_hits_true;
            // loop over clusterized hits
            for (size_t j = 0; j < clusterizedHits_bitID.size(); ++j)
            {
                if (((clusterizedHits_bitID[j] >> 8) & 0x7) != stationmatch || ((clusterizedHits_bitID[j] >> 11) & 0xF) != modulematch)
                {
                    continue;
                }
                // if (clusterizedHits_bitID[j] >> 28 & 0x3 || clusterizedHits_bitID[j] >> 30 & 0x3)
                //     continue; // NOTE not interested in support hits
                // Point ms_pos(ms_vx[j], ms_vy[j], ms_vz[j]); // in cm
                if ((clusterizedHits_bitID[j] >> 18 & 0xF) > segment_match - 1 && (clusterizedHits_bitID[j] >> 18 & 0xF) < segment_match + 2 && clusterizedHits_time[j] > time_match - 3 && clusterizedHits_time[j] < time_match + 3)
                { // && bitID_match != ms_bitID[j]) {

                    // require the bar to be within 7 bars of the matching bar
                    if (abs((clusterizedHits_bitID[j] >> 22 & 0x3F) - bar_match) < 8)
                    {
                        // close_by_hits.push_back(ms_pos);
                        close_by_hits_layer.push_back(clusterizedHits_bitID[j] >> 15 & 0x7);
                        close_by_hits_bitID.push_back(clusterizedHits_bitID[j]);
                        close_by_hits_time.push_back(clusterizedHits_time[j]);
                        if (clusterizedHits_id[j] != key[i])
                        {
                            close_by_hits_true.push_back(false);
                        }
                        else
                        {
                            close_by_hits_true.push_back(true);
                        }
                    }
                }
            }
            if (close_by_hits_layer.size() < 3)
            {
                if(verboseoutput)std::cout << "\t\tnot enough close by hits found " << close_by_hits_layer.size() << std::endl;
                continue;
            }
            else
            {
                if(verboseoutput)std::cout << "\t\tclose by hits found: " << close_by_hits_layer.size() << std::endl;
                // check how many layers are hit
                std::vector<int> layers_hit(4, 0);
                for (size_t j = 0; j < close_by_hits_layer.size(); ++j)
                {
                    layers_hit[close_by_hits_layer[j]]++;
                }
                int layers_hit_count = 0;
                for (size_t j = 0; j < layers_hit.size(); ++j)
                {
                    if (layers_hit[j] > 0)
                        layers_hit_count++;
                }
                if (layers_hit_count < 3)
                {
                    if(verboseoutput)std::cout << "\t\t\tnot enough layers hit " << layers_hit_count << std::endl;
                    continue;
                }
                else
                {
                    if(verboseoutput)std::cout << "\t\t\tlayers hit: " << layers_hit_count << std::endl;

                    // Create a vector of indices
                    std::vector<size_t> indices(close_by_hits_layer.size());
                    for (size_t j = 0; j < indices.size(); ++j)
                    {
                        indices[j] = j;
                    }

                    // Sort indices based on the values in close_by_hits_layer
                    std::sort(indices.begin(), indices.end(), [&close_by_hits_layer](size_t i1, size_t i2)
                              { return close_by_hits_layer[i1] < close_by_hits_layer[i2]; });

                    // Create copies of the original vectors
                    std::vector<int> sorted_close_by_hits_layer = close_by_hits_layer;
                    std::vector<int> sorted_close_by_hits_bitID = close_by_hits_bitID;
                    std::vector<float> sorted_close_by_hits_time = close_by_hits_time;
                    std::vector<bool> sorted_close_by_hits_true = close_by_hits_true;

                    // Apply the sorted order to the original vectors
                    for (size_t j = 0; j < indices.size(); ++j)
                    {
                        close_by_hits_layer[j] = sorted_close_by_hits_layer[indices[j]];
                        close_by_hits_bitID[j] = sorted_close_by_hits_bitID[indices[j]];
                        close_by_hits_time[j] = sorted_close_by_hits_time[indices[j]];
                        close_by_hits_true[j] = sorted_close_by_hits_true[indices[j]];
                    }

                    // print the position, segment, layer and bar of the close by hits
                    int layer0_true_bar = -1;
                    for (size_t j = 0; j < close_by_hits_layer.size(); ++j)
                    {
                        if(verboseoutput)std::cout << "\t\t\thit " << j << " layer: " << close_by_hits_layer[j] << " bar: " << (close_by_hits_bitID[j] >> 22 & 0x3F) << " segment " << (close_by_hits_bitID[j] >> 18 & 0xF) << " time: " << close_by_hits_time[j] << " true: " << close_by_hits_true[j] << std::endl;
                        // std::cout << "\t\t\thit " << j << " x: " << close_by_hits[j].x << " y: " << close_by_hits[j].y << " z: " << close_by_hits[j].z << "\tlayer: " << close_by_hits_layer[j] << " bar: " << (close_by_hits_bitID[j] >> 22 & 0x3F) << " segment " << (close_by_hits_bitID[j] >> 18 & 0xF) << " time: " << close_by_hits_time[j] << " true: " << close_by_hits_true[j] << std::endl;
                        if(close_by_hits_true[j]){
                            if(close_by_hits_layer[j] == 0){
                                layer0_true_bar = (close_by_hits_bitID[j] >> 22 & 0x3F);
                            } else {
                                hRelativeBarDifferenceToFirstLayer->Fill(close_by_hits_layer[j], (close_by_hits_bitID[j] >> 22 & 0x3F) - layer0_true_bar);
                                //if the station is 3 or 0, fill the opposite difference
                                if(((close_by_hits_bitID[j] >> 8) & 0x7) == 0 || ((close_by_hits_bitID[j] >> 8) & 0x7) == 3){
                                    hRelativeBarDifferenceToFirstLayer->Fill(close_by_hits_layer[j], layer0_true_bar - (close_by_hits_bitID[j] >> 22 & 0x3F));
                                } else {
                                    hRelativeBarDifferenceToFirstLayerIndiv[close_by_hits_layer[j]-1]->Fill((close_by_hits_bitID[j] >> 22 & 0x3F) - layer0_true_bar);
                                }
                            }
                        }
                    }
                    // find and group hits across multiple layers within a 1ns time window of each other, starting in layer0
                    std::vector<std::vector<int>> grouped_hits_layers;
                    std::vector<std::vector<int>> grouped_hits_bitIDs;
                    std::vector<std::vector<float>> grouped_hits_times;
                    std::vector<std::vector<bool>> grouped_hits_trues;
                    std::vector<std::vector<bool>> grouped_hits_removed;

                    for (size_t j = 0; j < close_by_hits_layer.size(); ++j)
                    {
                        // this loop starts from layer 0
                        if (close_by_hits_layer[j] != 0)
                        {
                            continue;
                        }
                        std::vector<int> grouped_hits_layers_temp;
                        std::vector<int> grouped_hits_bitIDs_temp;
                        std::vector<float> grouped_hits_times_temp;
                        std::vector<bool> grouped_hits_trues_temp;
                        std::vector<bool> grouped_hits_removed_temp;
                        // add layer0 hit
                        grouped_hits_layers_temp.push_back(close_by_hits_layer[j]);
                        grouped_hits_bitIDs_temp.push_back(close_by_hits_bitID[j]);
                        grouped_hits_times_temp.push_back(close_by_hits_time[j]);
                        grouped_hits_trues_temp.push_back(close_by_hits_true[j]);
                        grouped_hits_removed_temp.push_back(false);

                        for (int ilay = 1; ilay < 4; ++ilay) // loop over layers 1, 2, 3
                        {
                            for (size_t k = 0; k < close_by_hits_layer.size(); ++k)
                            {
                                if (close_by_hits_layer[k] != ilay) // select the current layer
                                    continue;
                                if (abs(close_by_hits_time[j] - close_by_hits_time[k]) < 1) // 1ns time difference
                                {
                                    //layer1 can only have a two bar difference
                                    if(ilay == 1 && abs((close_by_hits_bitID[j] >> 22 & 0x3F) - (close_by_hits_bitID[k] >> 22 & 0x3F)) > 2){
                                        continue;
                                    }
                                    //layer2 can only have a six bar difference
                                    if(ilay == 2 && abs((close_by_hits_bitID[j] >> 22 & 0x3F) - (close_by_hits_bitID[k] >> 22 & 0x3F)) > 6){
                                        continue;
                                    }
                                    //layer3 can only have a seven bar difference
                                    if(ilay == 3 && abs((close_by_hits_bitID[j] >> 22 & 0x3F) - (close_by_hits_bitID[k] >> 22 & 0x3F)) > 7){
                                        continue;
                                    }
                                    grouped_hits_layers_temp.push_back(close_by_hits_layer[k]);
                                    grouped_hits_bitIDs_temp.push_back(close_by_hits_bitID[k]);
                                    grouped_hits_times_temp.push_back(close_by_hits_time[k]);
                                    grouped_hits_trues_temp.push_back(close_by_hits_true[k]);
                                    grouped_hits_removed_temp.push_back(true);
                                }
                            }
                        }
                        grouped_hits_layers.push_back(grouped_hits_layers_temp);
                        grouped_hits_bitIDs.push_back(grouped_hits_bitIDs_temp);
                        grouped_hits_times.push_back(grouped_hits_times_temp);
                        grouped_hits_trues.push_back(grouped_hits_trues_temp);
                        grouped_hits_removed.push_back(grouped_hits_removed_temp);
                    }


                    //loop over groups and compare to vec_segment_combinations, reject any combinations that don't match
                    //the loop starts at layer0 to select the possible combinations to compare to, it then loops over the entries in grouped_hits_layers
                    
                    //create a vector to keep track of the combinations that matched with a group
                    std::vector<std::pair<int,int>> matched_combinations;

                    for (size_t j = 0; j < grouped_hits_layers.size(); ++j)
                    {
                        if (grouped_hits_layers[j].size() < 3)
                        {
                            continue;
                        }
                        int first_segment_hit = grouped_hits_bitIDs[j][0] >> 18 & 0xF;
                        for (size_t k = 0; k < vec_segment_combinations.size(); ++k)
                        {
                            if(vec_segment_combinations[k][0] != first_segment_hit){
                                continue;
                            }
                            int second_segment_combination = vec_segment_combinations[k][1];
                            int third_segment_combination = vec_segment_combinations[k][2];
                            int fourth_segment_combination = vec_segment_combinations[k][3];
                            
                            bool matchcomb[3] = {false, false, false};
                            bool hashit[3] = {false, false, false};
                            // loop over the hits in the group and check if the combination matches in the layers
                            for (size_t l = 0; l < grouped_hits_layers[j].size(); ++l)
                            {
                                if (grouped_hits_layers[j][l] == 1){
                                    hashit[0] = true;
                                    if((grouped_hits_bitIDs[j][l] >> 18 & 0xF) == second_segment_combination)
                                    {
                                        matchcomb[0] = true;
                                    }
                                }
                                if (grouped_hits_layers[j][l] == 2){
                                    hashit[1] = true;
                                    if((grouped_hits_bitIDs[j][l] >> 18 & 0xF) == third_segment_combination)
                                    {
                                        matchcomb[1] = true;
                                    }
                                }
                                if (grouped_hits_layers[j][l] == 3){
                                    hashit[2] = true;
                                    if((grouped_hits_bitIDs[j][l] >> 18 & 0xF) == fourth_segment_combination)
                                    {
                                        matchcomb[2] = true;
                                    }
                                }
                            }
                            //if there was a hit in the layer, but the combination didn't match, reject the combination
                            if((hashit[0] && !matchcomb[0]) || (hashit[1] && !matchcomb[1]) || (hashit[2] && !matchcomb[2])){
                                continue;
                            } else {
                                if(verboseoutput)std::cout << "\t\t\tgroup " << j << " matches combination " << k << " with segments " << first_segment_hit << " " << second_segment_combination << " " << third_segment_combination << " " << fourth_segment_combination << std::endl;
                                matched_combinations.push_back(std::make_pair(j,k));
                                //flag the hits that are part of the combination to not be removed
                                for (size_t l = 0; l < grouped_hits_layers[j].size(); ++l)
                                {
                                    if (grouped_hits_layers[j][l] == 1 && (grouped_hits_bitIDs[j][l] >> 18 & 0xF) == second_segment_combination){
                                        grouped_hits_removed[j][l] = false;
                                    }
                                    if (grouped_hits_layers[j][l] == 2 && (grouped_hits_bitIDs[j][l] >> 18 & 0xF) == third_segment_combination){
                                        grouped_hits_removed[j][l] = false;
                                    }
                                    if (grouped_hits_layers[j][l] == 3 && (grouped_hits_bitIDs[j][l] >> 18 & 0xF) == fourth_segment_combination){
                                        grouped_hits_removed[j][l] = false;
                                    }
                                }
                            }

                        }
                    }

                    // print groups of found hits
                    int nGroups = 0;
                    for (size_t j = 0; j < grouped_hits_layers.size(); ++j)
                    {
                        int nhitsGroup = 0;
                        if(verboseoutput)std::cout << "\t\t\tgroup " << j;
                        if(grouped_hits_layers[j].size() < 3){
                            if(verboseoutput)std::cout << " not enough hits in group";
                        }
                        if(verboseoutput)std::cout << std::endl;
                        
                        for (size_t k = 0; k < grouped_hits_layers[j].size(); ++k)
                        {
                            if(verboseoutput){
                                std::cout << "\t\t\t\thit " << k << " layer: " << grouped_hits_layers[j][k] << " bar: " << (grouped_hits_bitIDs[j][k] >> 22 & 0x3F) << " segment " << (grouped_hits_bitIDs[j][k] >> 18 & 0xF) << " time: " << grouped_hits_times[j][k] << " true: " << grouped_hits_trues[j][k];
                                if(grouped_hits_removed[j][k]){
                                    std::cout << " removed";
                                } else {
                                    nhitsGroup++;
                                }
                                std::cout << std::endl;
                            }
                            if(!grouped_hits_removed[j][k]){
                                nhitsGroup++;
                            }
                        }
                        if(nhitsGroup >= 3){
                            nGroups++;
                        }
                    }
                    hNGroups->Fill(nGroups);
                }
            }
        }

        numevt++;
    }
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    //plot hRelativeBarDifferenceToFirstLayer
    TCanvas c1("c1", "c1", 800, 600);
    c1.SetLogz();
    //set tick marks on top and right of plot
    hRelativeBarDifferenceToFirstLayer->GetXaxis()->SetTitle("Layer");
    hRelativeBarDifferenceToFirstLayer->GetYaxis()->SetTitle("Relative bar difference to first layer");
    hRelativeBarDifferenceToFirstLayer->Draw("colz");
    c1.SaveAs(Form("%shRelativeBarDifferenceToFirstLayer.pdf", outputdir.Data()));

    //make a three panel canvas and plot hRelativeBarDifferenceToFirstLayerIndiv
    TCanvas c2("c2", "c2", 1200, 800);
    c2.Divide(3, 1);
    for (int i = 0; i < 3; i++)
    {
        c2.cd(i + 1);
        hRelativeBarDifferenceToFirstLayerIndiv[i]->GetXaxis()->SetTitle("Relative bar difference to first layer");
        hRelativeBarDifferenceToFirstLayerIndiv[i]->GetYaxis()->SetTitle("Counts");
        hRelativeBarDifferenceToFirstLayerIndiv[i]->Draw();
    }
    c2.SaveAs(Form("%shRelativeBarDifferenceToFirstLayerIndiv.pdf", outputdir.Data()));

    //plot hNGroups
    TCanvas c3("c3", "c3", 800, 600);
    // c3.SetLogy();
    hNGroups->GetXaxis()->SetTitle("Number of groups");
    hNGroups->GetYaxis()->SetTitle("Counts");
    hNGroups->Draw();
    c3.SaveAs(Form("%shNGroups.pdf", outputdir.Data()));

    return 0;
}