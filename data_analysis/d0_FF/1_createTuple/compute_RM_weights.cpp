#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TF1.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <cstring>

using namespace std;

vector<string> read_list_file(const string &fname) {
    vector<string> files;
    ifstream inf(fname);
    if (!inf) return files;
    string line;
    while (getline(inf, line)) {
        // strip
        while (!line.empty() && isspace(line.back())) line.pop_back();
        size_t pos = 0;
        while (pos < line.size() && isspace(line[pos])) ++pos;
        if (pos >= line.size()) continue;
        if (line[pos] == '#') continue;
        files.push_back(line.substr(pos));
    }
    return files;
}

void add_files_to_chain(TChain* chain, const vector<string> &files) {
    for (const auto &f : files) {
        if (f.empty()) continue;
        if (gSystem->AccessPathName(f.c_str()) == 0) {
            chain->Add(f.c_str());
        } else {
            cerr << "Warning: file not found: " << f << "\n";
        }
    }
}


void compute_RM_weights(const TString &pp_input, const TString &ppb_input, const TString &output_basename) {
    gStyle->SetOptStat(0);
    // Accept two TString inputs (either a single .root path or a .txt list file)
    // Third argument: base name for outputs (folder name and .root/.json/.png names)
    string tree_name = "d0jets";
    string branch = "event_multiplicity";
    int nbins = 800;
    bool isPbp = false;
    if (output_basename.Contains("Pbp")){
        isPbp = true;
        //if output_basename contains Pbp, increase the bins to 1600 to cover the higher multiplicity tail in pPb
        nbins = 2000;
    }
    // labels used in titles/legends and filenames depending on whether input indicates 'Pbp'
    std::string ppbLabel = isPbp ? "Pbp" : "pPb"; // human-readable label
    std::string ppbFileLabel = isPbp ? "Pbp" : "pPb"; // token used in filenames
    string out_json;
    string out_root;

    // Create dated output directory: outputs/YYYY-MM-DD
    time_t now = time(nullptr);
    struct tm lt;
    #if defined(_WIN32)
    localtime_s(&lt, &now);
    #else
    localtime_r(&now, &lt);
    #endif
    char datebuf[32];
    strftime(datebuf, sizeof(datebuf), "%Y-%m-%d", &lt);
    string base_dir = "outputs";
    string date_dir = base_dir + "/" + string(datebuf);
    string output_base = string(output_basename.Data());
    if (output_base.empty()) output_base = "multiplicity_weights";
    string out_dir = date_dir + "/" + output_base;

    // create base, dated, and final output directories if necessary
    struct stat st;
    if (stat(base_dir.c_str(), &st) != 0) {
        if (mkdir(base_dir.c_str(), 0755) != 0 && errno != EEXIST) {
            cerr << "Warning: could not create directory: " << base_dir << " (" << strerror(errno) << ")\n";
        }
    }
    if (stat(date_dir.c_str(), &st) != 0) {
        if (mkdir(date_dir.c_str(), 0755) != 0 && errno != EEXIST) {
            cerr << "Warning: could not create directory: " << date_dir << " (" << strerror(errno) << ")\n";
        }
    }
    if (stat(out_dir.c_str(), &st) != 0) {
        if (mkdir(out_dir.c_str(), 0755) != 0 && errno != EEXIST) {
            cerr << "Warning: could not create directory: " << out_dir << " (" << strerror(errno) << ")\n";
        }
    }

    // place outputs in the final output directory
    out_json = out_dir + "/" + output_base + ".json";
    out_root = out_dir + "/" + output_base + ".root";
    double min_stat = 1e-12;

    TChain chain_pp(tree_name.c_str());
    TChain chain_ppb(tree_name.c_str());

    // Helper: if input ends with .txt treat as list file
    auto add_input = [&](const TString &input, TChain &chain) {
        string s = string(input.Data());
        if (s.size() >= 4 && s.substr(s.size()-4) == ".txt") {
            vector<string> tmp = read_list_file(s);
            add_files_to_chain(&chain, tmp);
        } else {
            if (gSystem->AccessPathName(s.c_str()) == 0) chain.Add(s.c_str());
            else cerr << "Warning: file not found: " << s << "\n";
        }
    };

    add_input(pp_input, chain_pp);
    add_input(ppb_input, chain_ppb);

    Long64_t entries_pp = chain_pp.GetEntries();
    Long64_t entries_ppb = chain_ppb.GetEntries();

    if (entries_pp == 0) {
        cerr << "ERROR: No entries found in pp chain. Provide files via --pp or --pp-list\n";
        return;
    }
    if (entries_ppb == 0) {
        cerr << "ERROR: No entries found in pPb chain. Provide files via --ppb or --ppb-list\n";
        return;
    }

    cout << "Filling pp histogram from " << entries_pp << " entries..." << endl;
    cout << "Filling pPb histogram from " << entries_ppb << " entries..." << endl;

    TString hppname("h_pp_raw");
    TString hppbname("h_ppb_raw");
    TH1D *h_pp_raw = new TH1D(hppname.Data(), "pp multiplicity", nbins, -0.5, nbins - 0.5);
    TH1D *h_ppb_raw = new TH1D(hppbname.Data(), "pPb multiplicity", nbins, -0.5, nbins - 0.5);
    h_pp_raw->Sumw2();
    h_ppb_raw->Sumw2();

    TString drawExpr(branch.c_str());
    TString drawCmd_pp = TString::Format("%s>>%s", drawExpr.Data(), hppname.Data());
    TString drawCmd_ppb = TString::Format("%s>>%s", drawExpr.Data(), hppbname.Data());

    chain_pp.Draw(drawCmd_pp, "", "goff");
    chain_ppb.Draw(drawCmd_ppb, "", "goff");

    double total_pp = h_pp_raw->Integral();
    double total_ppb = h_ppb_raw->Integral();
    cout << "pp total events: " << total_pp << ", pPb total events: " << total_ppb << endl;

    if (total_pp <= 0 || total_ppb <= 0) {
        cerr << "ERROR: Empty histograms, cannot compute weights" << endl;
        return;
    }

    TH1D *h_pp = (TH1D*)h_pp_raw->Clone("h_pp");
    TH1D *h_ppb = (TH1D*)h_ppb_raw->Clone("h_ppb");
    if (total_pp > 0) h_pp->Scale(1.0 / total_pp);
    if (total_ppb > 0) h_ppb->Scale(1.0 / total_ppb);

    // Optional: fit pp tail for Pbp case and use fit to extend pp fractions
    TF1 *f_pp_tail = nullptr;
    const double fitThreshold = 300.0;
    if (output_basename.Contains("Pbp")) {
        double fitMin = fitThreshold;
        double fitMax = h_pp->GetXaxis()->GetXmax();
        f_pp_tail = new TF1("f_pp_tail", "[0]*exp(-[1]*x)", fitMin, fitMax);
        int startBin = h_pp->FindBin(fitMin);
        double startVal = h_pp->GetBinContent(startBin);
        f_pp_tail->SetParameters((startVal>0?startVal:1e-6), 0.01);
        h_pp->Fit(f_pp_tail, "RQ");
    }

    // compute weights (use fit prediction for pp fraction above threshold when available)
    map<int, map<string, double>> weights;
    for (int ibin = 1; ibin <= nbins; ++ibin) {
        double c_pp = h_pp_raw->GetBinContent(ibin);
        double c_ppb = h_ppb_raw->GetBinContent(ibin);
        double n_pp = h_pp->GetBinContent(ibin);
        double n_ppb = h_ppb->GetBinContent(ibin);
        int bin_center = (int)lround(h_pp_raw->GetBinCenter(ibin));
        if (f_pp_tail && (double)bin_center > fitThreshold) {
            double pred = f_pp_tail->Eval((double)bin_center);
            if (pred > 0) n_pp = pred;
            else n_pp = 0.0;
        }
        double w = 0.0;
        if (n_pp > min_stat) w = n_ppb / n_pp;
        weights[bin_center]["weight"] = w;
        weights[bin_center]["pp_count"] = c_pp;
        weights[bin_center]["ppb_count"] = c_ppb;
        weights[bin_center]["pp_frac"] = n_pp;
        weights[bin_center]["ppb_frac"] = n_ppb;
    }

    // Create a histogram of weights for plotting
    TH1D *h_weights = new TH1D("h_weights", "Multiplicity reweighting (pPb/pp)", nbins, -0.5, nbins - 0.5);
    h_weights->Sumw2();
    for (const auto &kv : weights) {
        int mult = kv.first;
        double w = kv.second.at("weight");
        int bin = h_weights->FindBin((double)mult);
        h_weights->SetBinContent(bin, w);
        h_weights->SetBinError(bin, 0.0);
    }

    // Smooth weights by averaging over non-overlapping groups of 10 multiplicity values
    const int groupSize = 10;
    std::map<int, double> smoothed;
    int maxMult = nbins - 1;
    for (int start = 0; start <= maxMult; start += groupSize) {
        int end = std::min(start + groupSize - 1, maxMult);
        double sum = 0.0;
        int count = 0;
        for (int m = start; m <= end; ++m) {
            auto it = weights.find(m);
            if (it != weights.end()) {
                sum += it->second.at("weight");
                ++count;
            }
        }
        double avg = (count > 0) ? sum / (double)count : 0.0;
        for (int m = start; m <= end; ++m) smoothed[m] = avg;
    }

    // histogram for smoothed weights
    TH1D *h_weights_smoothed = new TH1D("h_weights_smoothed", "Smoothed multiplicity reweighting (pPb/pp)", nbins, -0.5, nbins - 0.5);
    h_weights_smoothed->Sumw2();
    for (const auto &kv : smoothed) {
        int mult = kv.first;
        double w = kv.second;
        int bin = h_weights_smoothed->FindBin((double)mult);
        h_weights_smoothed->SetBinContent(bin, w);
        h_weights_smoothed->SetBinError(bin, 0.0);
    }

    // Draw and save the weights plot to PNG in the output directory
    TCanvas *c_weights = new TCanvas("c_weights", "Multiplicity weights", 800, 600);
    c_weights->SetGrid();
    c_weights->SetLogy(1);
    h_weights->SetMarkerStyle(20);
    h_weights->Rebin(10); // Optional: rebin for smoother plot
    h_weights->GetYaxis()->SetRangeUser(0.9, h_weights->GetMaximum() * 1.2);
    {
        std::string title = std::string("Multiplicity reweighting (") + ppbLabel + "/pp);Event multiplicity;Weight (" + ppbLabel + "/pp fraction)";
        h_weights->SetTitle(title.c_str());
    }
    h_weights->SetMarkerSize(0.9);
    h_weights->Draw("E1");
    string png_path = out_dir + "/" + output_base + ".png";
    c_weights->SaveAs(png_path.c_str());

    // Overlay pp and pPb multiplicity distributions (normalized) on log-y
    TCanvas *c_mult = new TCanvas("c_mult_pp_ppb", "pp vs pPb multiplicity", 900, 700);
    c_mult->SetLogy(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    h_pp->SetLineColor(kBlue);
    h_pp->SetLineWidth(2);
    h_pp->SetMarkerStyle(20);
    h_pp->SetMarkerColor(kBlue);
    h_pp->SetTitle("Normalized multiplicity distributions;Event multiplicity;Normalized counts");
    h_pp->SetMinimum(1e-9);
    h_pp->Draw("E1");

    h_ppb->SetLineColor(kRed);
    h_ppb->SetLineWidth(2);
    h_ppb->SetMarkerStyle(21);
    h_ppb->SetMarkerColor(kRed);
    h_ppb->Draw("E1 SAME");

    TLegend legm(0.6, 0.7, 0.92, 0.88);
    legm.SetBorderSize(0);
    legm.SetFillStyle(0);
    legm.AddEntry(h_pp, "pp (normalized)", "lep");
    {
        std::string entry = ppbLabel + std::string(" (normalized)");
        legm.AddEntry(h_ppb, entry.c_str(), "lep");
    }
    if (f_pp_tail) {
        f_pp_tail->SetLineColor(kGreen+2);
        f_pp_tail->SetLineWidth(2);
        f_pp_tail->SetLineStyle(7);
        f_pp_tail->Draw("SAME");
        legm.AddEntry(f_pp_tail, "pp tail fit (exp)", "l");
    }
    legm.Draw();

    std::string png_mult = out_dir + "/multiplicity_pp_vs_" + ppbFileLabel + ".png";
    std::string pdf_mult = out_dir + "/multiplicity_pp_vs_" + ppbFileLabel + ".pdf";
    c_mult->SaveAs(png_mult.c_str());
    c_mult->SaveAs(pdf_mult.c_str());

    auto make_1d_weights = [&](const std::string &tag,
                               const std::string &plotTitle,
                               const std::string &xTitle,
                               const std::string &branchExpr,
                               int nBins,
                               double minX,
                               double maxX,
                               const std::string &rootStem) {
        TH1D *h_pp_1d = new TH1D((tag + "_pp_raw").c_str(), (plotTitle + " (pp)").c_str(), nBins, minX, maxX);
        TH1D *h_ppb_1d = new TH1D((tag + "_ppb_raw").c_str(), (plotTitle + " (" + ppbLabel + ")").c_str(), nBins, minX, maxX);
        h_pp_1d->Sumw2();
        h_ppb_1d->Sumw2();

        TString drawExpr(branchExpr.c_str());
        chain_pp.Draw(TString::Format("%s>>%s", drawExpr.Data(), h_pp_1d->GetName()), "", "goff");
        chain_ppb.Draw(TString::Format("%s>>%s", drawExpr.Data(), h_ppb_1d->GetName()), "", "goff");

        double total_pp_1d = h_pp_1d->Integral();
        double total_ppb_1d = h_ppb_1d->Integral();
        if (total_pp_1d <= 0 || total_ppb_1d <= 0) {
            cerr << "ERROR: Empty histograms for " << tag << ", cannot compute weights" << endl;
            delete h_pp_1d;
            delete h_ppb_1d;
            delete h_ppb_1d;
            return;
        }

        TH1D *h_pp_norm = (TH1D*)h_pp_1d->Clone((tag + "_pp_norm").c_str());
        TH1D *h_ppb_norm = (TH1D*)h_ppb_1d->Clone((tag + "_ppb_norm").c_str());
        h_pp_norm->Scale(1.0 / total_pp_1d);
        h_ppb_norm->Scale(1.0 / total_ppb_1d);

        TH1D *h_weights_1d = (TH1D*)h_pp_norm->Clone((tag + "_weights_histo").c_str());
        h_weights_1d->Reset("ICES");
        h_weights_1d->SetTitle((plotTitle + " reweighting;" + xTitle + ";Weight (" + ppbLabel + "/pp)").c_str());
        h_weights_1d->Sumw2();

        TTree tree_1d((tag + "_weights").c_str(), (plotTitle + " weights").c_str());
        Int_t bin_index = 0;
        Double_t bin_center = 0.0;
        Double_t weight_val = 0.0;
        Double_t pp_frac = 0.0;
        Double_t ppb_frac = 0.0;
        tree_1d.Branch("bin", &bin_index, "bin/I");
        tree_1d.Branch("center", &bin_center, "center/D");
        tree_1d.Branch("weight", &weight_val, "weight/D");
        tree_1d.Branch("pp_frac", &pp_frac, "pp_frac/D");
        tree_1d.Branch("ppb_frac", &ppb_frac, "ppb_frac/D");

        for (int ibin = 1; ibin <= nBins; ++ibin) {
            double n_pp = h_pp_norm->GetBinContent(ibin);
            double n_ppb = h_ppb_norm->GetBinContent(ibin);
            weight_val = 0.0;
            if (n_pp > min_stat) weight_val = n_ppb / n_pp;
            bin_index = ibin;
            bin_center = h_pp_norm->GetXaxis()->GetBinCenter(ibin);
            pp_frac = n_pp;
            ppb_frac = n_ppb;
            h_weights_1d->SetBinContent(ibin, weight_val);
            h_weights_1d->SetBinError(ibin, 0.0);
            tree_1d.Fill();
        }

        TCanvas *c_compare = new TCanvas((tag + "_compare_canvas").c_str(), (plotTitle + " comparison").c_str(), 900, 700);
        c_compare->SetLogy(1);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);
        h_pp_norm->SetLineColor(kBlue);
        h_pp_norm->SetLineWidth(2);
        h_pp_norm->SetMarkerStyle(20);
        h_pp_norm->SetMarkerColor(kBlue);
        h_pp_norm->GetXaxis()->SetTitle(xTitle.c_str());
        h_pp_norm->GetYaxis()->SetTitle("Normalized counts");
        h_pp_norm->SetMinimum(1e-9);
        h_pp_norm->Draw("E1");

        h_ppb_norm->SetLineColor(kRed);
        h_ppb_norm->SetLineWidth(2);
        h_ppb_norm->SetMarkerStyle(21);
        h_ppb_norm->SetMarkerColor(kRed);
        h_ppb_norm->Draw("E1 SAME");

        TLegend leg_1d(0.6, 0.7, 0.92, 0.88);
        leg_1d.SetBorderSize(0);
        leg_1d.SetFillStyle(0);
        leg_1d.AddEntry(h_pp_norm, "pp (normalized)", "lep");
        {
            std::string entry = ppbLabel + std::string(" (normalized)");
            leg_1d.AddEntry(h_ppb_norm, entry.c_str(), "lep");
        }
        leg_1d.Draw();

        TCanvas *c_weights = new TCanvas((tag + "_weights_canvas").c_str(), (plotTitle + " weights").c_str(), 900, 700);
        c_weights->SetGrid();
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);
        h_weights_1d->SetMarkerStyle(20);
        h_weights_1d->SetMarkerSize(0.9);
        h_weights_1d->SetMarkerColor(kBlue + 1);
        h_weights_1d->SetLineColor(kBlue + 1);
        h_weights_1d->GetXaxis()->SetTitle(xTitle.c_str());
        h_weights_1d->GetYaxis()->SetTitle((ppbLabel + " / pp (weight)").c_str());
        h_weights_1d->SetMinimum(0.0);
        h_weights_1d->Draw("E1");

        std::string root_path = out_dir + "/" + rootStem + ".root";
        std::string png_compare = out_dir + "/" + rootStem + "_comparison.png";
        std::string pdf_compare = out_dir + "/" + rootStem + "_comparison.pdf";
        std::string png_weights = out_dir + "/" + rootStem + ".png";
        std::string pdf_weights = out_dir + "/" + rootStem + ".pdf";

        c_compare->SaveAs(png_compare.c_str());
        c_compare->SaveAs(pdf_compare.c_str());
        c_weights->SaveAs(png_weights.c_str());
        c_weights->SaveAs(pdf_weights.c_str());

        TFile fout_1d(root_path.c_str(), "RECREATE");
        h_pp_1d->Write();
        h_ppb_1d->Write();
        h_pp_norm->Write();
        h_ppb_norm->Write();
        h_weights_1d->Write();
        tree_1d.Write();
        c_compare->Write("comparison_canvas");
        c_weights->Write("weights_canvas");
        fout_1d.Write();
        fout_1d.Close();

        delete c_compare;
        delete c_weights;
        delete h_weights_1d;
        delete h_pp_norm;
        delete h_ppb_norm;
        delete h_pp_1d;
        delete h_ppb_1d;
    };

    // Separate 1D weights and plots for D0 pT and jet pT spectra.
    make_1d_weights("d0pt",
                    "D^{0} p_{T}",
                    "p_{T} (GeV/c)",
                    "d0_pt",
                    50,
                    0.0,
                    50.0,
                    output_base + "_d0pt_weights");

    make_1d_weights("jetpt",
                    "Jet p_{T}",
                    "p_{T} (GeV/c)",
                    "jet_pt",
                    100,
                    0.0,
                    200.0,
                    output_base + "_jetpt_weights");

    // Compute jet-pT vs D0-pT weights (pPb / pp) and save to a separate ROOT file
    {
        const int nBinsJet = 100;
        const double jetMax = 200.0;
        const int nBinsD0 = 50;
        const double d0Max = 50.0;

        TH2D *h_jet_d0_pp = new TH2D("h_jet_d0_pp",
                                     "pp jet-p_{T} vs D^{0} p_{T};jet p_{T} (GeV/c);D^{0} p_{T} (GeV/c)",
                                     nBinsJet, 0.0, jetMax,
                                     nBinsD0, 0.0, d0Max);
        TH2D *h_jet_d0_ppb = new TH2D("h_jet_d0_ppb",
                                      "pPb jet-p_{T} vs D^{0} p_{T};jet p_{T} (GeV/c);D^{0} p_{T} (GeV/c)",
                                      nBinsJet, 0.0, jetMax,
                                      nBinsD0, 0.0, d0Max);

        chain_pp.Draw("d0_pt:jet_pt>>h_jet_d0_pp", "", "goff");
        chain_ppb.Draw("d0_pt:jet_pt>>h_jet_d0_ppb", "", "goff");

        TH2D *h_jet_d0_pp_norm = (TH2D*)h_jet_d0_pp->Clone("h_jet_d0_pp_norm");
        TH2D *h_jet_d0_ppb_norm = (TH2D*)h_jet_d0_ppb->Clone("h_jet_d0_ppb_norm");
        double ip = h_jet_d0_pp_norm->Integral();
        double ipb = h_jet_d0_ppb_norm->Integral();
        if (ip > 0) h_jet_d0_pp_norm->Scale(1.0 / ip);
        if (ipb > 0) h_jet_d0_ppb_norm->Scale(1.0 / ipb);

        TH2D *h_jet_d0_weights = new TH2D("h_jet_d0_weights",
                                          "Jet-p_{T} vs D^{0} p_{T} reweighting (pPb/pp);jet p_{T} (GeV/c);D^{0} p_{T} (GeV/c);Weight",
                                          nBinsJet, 0.0, jetMax,
                                          nBinsD0, 0.0, d0Max);
        h_jet_d0_weights->Sumw2();

        for (int ix = 1; ix <= nBinsJet; ++ix) {
            for (int iy = 1; iy <= nBinsD0; ++iy) {
                double n_pp = h_jet_d0_pp_norm->GetBinContent(ix, iy);
                double n_ppb = h_jet_d0_ppb_norm->GetBinContent(ix, iy);
                double w = 0.0;
                if (n_pp > min_stat) w = n_ppb / n_pp;
                h_jet_d0_weights->SetBinContent(ix, iy, w);
                h_jet_d0_weights->SetBinError(ix, iy, 0.0);
            }
        }

        TCanvas *c_jet_d0_maps = new TCanvas("c_jet_d0_maps", "Jet pT vs D0 pT maps", 1800, 600);
        c_jet_d0_maps->Divide(3, 1);

        c_jet_d0_maps->cd(1);
        gPad->SetLogz(1);
        gPad->SetRightMargin(0.14);
        h_jet_d0_pp_norm->SetStats(0);
        h_jet_d0_pp_norm->SetTitle("pp normalized jet-p_{T} vs D^{0} p_{T}");
        h_jet_d0_pp_norm->GetZaxis()->SetTitle("Normalized counts");
        h_jet_d0_pp_norm->Draw("COLZ");

        c_jet_d0_maps->cd(2);
        gPad->SetLogz(1);
        gPad->SetRightMargin(0.14);
        h_jet_d0_ppb_norm->SetStats(0);
        h_jet_d0_ppb_norm->SetTitle("pPb normalized jet-p_{T} vs D^{0} p_{T}");
        h_jet_d0_ppb_norm->GetZaxis()->SetTitle("Normalized counts");
        h_jet_d0_ppb_norm->Draw("COLZ");

        c_jet_d0_maps->cd(3);
        gPad->SetLogz(1);
        gPad->SetRightMargin(0.14);
        h_jet_d0_weights->SetStats(0);
        h_jet_d0_weights->SetTitle("pPb / pp jet-p_{T} vs D^{0} p_{T}");
        h_jet_d0_weights->GetZaxis()->SetTitle("Weight");
        h_jet_d0_weights->Draw("COLZ");

        std::string png_jet_maps = out_dir + "/" + output_base + "_jetd0pt_maps.png";
        std::string pdf_jet_maps = out_dir + "/" + output_base + "_jetd0pt_maps.pdf";
        c_jet_d0_maps->SaveAs(png_jet_maps.c_str());
        c_jet_d0_maps->SaveAs(pdf_jet_maps.c_str());

        std::string out_root_jet = out_dir + "/" + output_base + "_jetd0pt_weights.root";
        TFile foutj(out_root_jet.c_str(), "RECREATE");
        h_jet_d0_pp_norm->Write();
        h_jet_d0_ppb_norm->Write();
        h_jet_d0_weights->Write();

        TTree tj("jet_d0_pt_weights", "jet-pT vs D0-pT weights");
        Int_t jet_bin = 0;
        Int_t d0_bin = 0;
        Double_t jet_center = 0.0;
        Double_t d0_center = 0.0;
        Double_t weight_val = 0.0;
        Double_t pp_frac = 0.0;
        Double_t ppb_frac = 0.0;
        tj.Branch("jet_bin", &jet_bin, "jet_bin/I");
        tj.Branch("d0_bin", &d0_bin, "d0_bin/I");
        tj.Branch("jet_center", &jet_center, "jet_center/D");
        tj.Branch("d0_center", &d0_center, "d0_center/D");
        tj.Branch("weight", &weight_val, "weight/D");
        tj.Branch("pp_frac", &pp_frac, "pp_frac/D");
        tj.Branch("ppb_frac", &ppb_frac, "ppb_frac/D");

        for (int ix = 1; ix <= nBinsJet; ++ix) {
            for (int iy = 1; iy <= nBinsD0; ++iy) {
                jet_bin = ix;
                d0_bin = iy;
                jet_center = h_jet_d0_weights->GetXaxis()->GetBinCenter(ix);
                d0_center = h_jet_d0_weights->GetYaxis()->GetBinCenter(iy);
                weight_val = h_jet_d0_weights->GetBinContent(ix, iy);
                pp_frac = h_jet_d0_pp_norm->GetBinContent(ix, iy);
                ppb_frac = h_jet_d0_ppb_norm->GetBinContent(ix, iy);
                tj.Fill();
            }
        }
        tj.Write();
        foutj.Write();
        foutj.Close();

        TCanvas *c_jet_weights_plot = new TCanvas("c_jet_weights_plot", "Jet pT vs D0 pT Weights", 1000, 800);
        c_jet_weights_plot->SetGrid();
        h_jet_d0_weights->SetStats(0);
        h_jet_d0_weights->GetZaxis()->SetTitle("pPb / pp (weight)");
        h_jet_d0_weights->Draw("COLZ");
        std::string png_jetw = out_dir + "/" + output_base + "_jetd0pt_weights.png";
        std::string pdf_jetw = out_dir + "/" + output_base + "_jetd0pt_weights.pdf";
        c_jet_weights_plot->SaveAs(png_jetw.c_str());
        c_jet_weights_plot->SaveAs(pdf_jetw.c_str());

        delete c_jet_d0_maps;
        delete c_jet_weights_plot;
        delete h_jet_d0_weights;
        delete h_jet_d0_pp_norm;
        delete h_jet_d0_ppb_norm;
        delete h_jet_d0_pp;
        delete h_jet_d0_ppb;
    }

    // Comparison plot: original weights vs smoothed weights
    {
        // const int rebinFactor = 10; // match rebin used for original weights plotting
        TH1D *h_smoothed_reb = (TH1D*)h_weights_smoothed->Clone("h_weights_smoothed_rebinned");
        // if (rebinFactor > 1) h_smoothed_reb->Rebin(rebinFactor);

        // original h_weights was rebinned earlier for display; use it as-is
        double ymax = std::max(h_weights->GetMaximum(), h_smoothed_reb->GetMaximum());
        double ymin = 0.0;
        if (ymax <= 0) ymax = 1.0;

        TCanvas *c_weights_comp = new TCanvas("c_weights_comp", "Weights comparison", 900, 700);
        c_weights_comp->SetGrid();
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);

        h_weights->SetMarkerStyle(20);
        h_weights->SetMarkerSize(0.9);
        h_weights->SetLineColor(kBlue);
        h_weights->SetMarkerColor(kBlue);
        h_weights->GetXaxis()->SetTitle("Event multiplicity");
        {
            std::string ytitle = ppbLabel + std::string(" / pp (fraction)");
            h_weights->GetYaxis()->SetTitle(ytitle.c_str());
        }
        h_weights->GetYaxis()->SetRangeUser(ymin, ymax * 1.2);
        h_weights->Draw("E1");

        h_smoothed_reb->SetLineColor(kRed);
        h_smoothed_reb->SetLineWidth(2);
        h_smoothed_reb->SetLineStyle(2);
        h_smoothed_reb->SetMarkerStyle(24);
        h_smoothed_reb->SetMarkerColor(kRed);
        h_smoothed_reb->Draw("HIST SAME");

        TLegend legw(0.6, 0.7, 0.92, 0.88);
        legw.SetBorderSize(0);
        legw.SetFillStyle(0);
        legw.AddEntry(h_weights, "original weights (rebinned)", "lep");
        legw.AddEntry(h_smoothed_reb, "smoothed weights (rebinned)", "l");
        legw.Draw();

        std::string png_comp = out_dir + "/" + output_base + "_weights_comparison.png";
        std::string pdf_comp = out_dir + "/" + output_base + "_weights_comparison.pdf";
        c_weights_comp->SaveAs(png_comp.c_str());
        c_weights_comp->SaveAs(pdf_comp.c_str());

        delete h_smoothed_reb;
    }

    // write JSON
    ofstream jf(out_json);
    if (!jf) {
        cerr << "ERROR: cannot open " << out_json << " for writing" << endl;
        return;
    }
    jf << fixed << setprecision(6);
    jf << "{\n";
    jf << "  \"tree\": \"" << tree_name << "\",\n";
    jf << "  \"branch\": \"" << branch << "\",\n";
    jf << "  \"nbins\": " << nbins << ",\n";
    jf << "  \"total_pp\": " << total_pp << ",\n";
    jf << "  \"total_ppb\": " << total_ppb << ",\n";
    jf << "  \"weights\": {\n";
    bool first = true;
    for (const auto &kv : weights) {
        int mult = kv.first;
        const auto &info = kv.second;
        if (!first) jf << ",\n";
        first = false;
        jf << "    \"" << mult << "\": {\n";
        jf << "      \"weight\": " << info.at("weight") << ",\n";
        jf << "      \"pp_count\": " << info.at("pp_count") << ",\n";
        jf << "      \"ppb_count\": " << info.at("ppb_count") << ",\n";
        jf << "      \"pp_frac\": " << info.at("pp_frac") << ",\n";
        jf << "      \"ppb_frac\": " << info.at("ppb_frac") << "\n";
        jf << "    }";
    }
    jf << "\n  }\n}" << endl;
    jf.close();
    cout << "Wrote weights JSON to " << out_json << endl;

    // write ROOT outputs
    TFile fout(out_root.c_str(), "RECREATE");
    h_pp_raw->Write();
    h_ppb_raw->Write();
    h_pp->Write();
    h_ppb->Write();

    // also write the original weights histogram
    h_weights->Write();

    // tree of weights
    TTree tree_w("multiplicity_weights", "multiplicity weights");
    Int_t mult_i = 0;
    Double_t weight_d = 0.0, ppcount_d = 0.0, ppbcount_d = 0.0;
    tree_w.Branch("multiplicity", &mult_i, "multiplicity/I");
    tree_w.Branch("weight", &weight_d, "weight/D");
    tree_w.Branch("pp_count", &ppcount_d, "pp_count/D");
    tree_w.Branch("ppb_count", &ppbcount_d, "ppb_count/D");

    for (const auto &kv : weights) {
        mult_i = kv.first;
        weight_d = kv.second.at("weight");
        ppcount_d = kv.second.at("pp_count");
        ppbcount_d = kv.second.at("ppb_count");
        tree_w.Fill();
    }
    tree_w.Write();
    fout.Write();
    fout.Close();

    // Write second ROOT file with smoothed weights
    std::string out_root_smoothed = out_dir + "/" + output_base + "_smoothed.root";
    TFile fout2(out_root_smoothed.c_str(), "RECREATE");
    h_pp_raw->Write();
    h_ppb_raw->Write();
    h_pp->Write();
    h_ppb->Write();
    h_weights_smoothed->Write();

    // tree of smoothed weights
    TTree tree_ws("multiplicity_weights", "multiplicity weights (smoothed)");
    Int_t mult_si = 0;
    Double_t weight_sd = 0.0, ppcount_sd = 0.0, ppbcount_sd = 0.0;
    tree_ws.Branch("multiplicity", &mult_si, "multiplicity/I");
    tree_ws.Branch("weight", &weight_sd, "weight/D");
    tree_ws.Branch("pp_count", &ppcount_sd, "pp_count/D");
    tree_ws.Branch("ppb_count", &ppbcount_sd, "ppb_count/D");
    for (const auto &kv : smoothed) {
        mult_si = kv.first;
        weight_sd = kv.second;
        const auto &info = weights[mult_si];
        ppcount_sd = info.at("pp_count");
        ppbcount_sd = info.at("ppb_count");
        tree_ws.Fill();
    }
    tree_ws.Write();
    fout2.Write();
    fout2.Close();

    cout << "Wrote ROOT output with histograms and TTree to " << out_root << endl;
    cout << "Done." << endl;
}
