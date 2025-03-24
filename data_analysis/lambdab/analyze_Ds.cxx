#include "../plottingheader.h"



void analyze_Ds(bool makeTree = false) {
    // // Open the .root file
    // TFile *file = TFile::Open("/media/niviths/local/analysis_code/data_analysis/lambdab/Ds_plots/D0Kpi_DecayFunTuple_selected.root");
    TFile *file = TFile::Open("/media/niviths/local/analysis_code/data_analysis/lambdab/Ds_plots/D0Kpi_DecayFunTuple_selected.root");
    // TFile *file = TFile::Open("/media/niviths/SSD2/lhcb_analysis_SSD/20250318_fullDownloadOfAP/DsKKpi_DecayFunTuple.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }
    // Get the tree from the file
    TTree *tree;
    file->GetObject("DecayFunTuple", tree);
    if (!tree) {
        std::cerr << "Error getting tree from file" << std::endl;
        file->Close();
        return;
    }

    TString outputdir = "Ds_plots/";
    system("mkdir -p " + outputdir);

    // Set up variables to hold data from the tree
    Int_t SMOGInjectedGas;
    Int_t nPVs;
    tree->SetBranchAddress("SMOGInjectedGas", &SMOGInjectedGas);
    tree->SetBranchAddress("nPVs", &nPVs);
    TString gasNames[] = {"NONE", "HYDROGEN", "DEUTERIUM", "HELIUM", "NITROGEN", "OXYGEN", "NEON", "ARGON", "KRYPTON", "XENON"};

    UChar_t BUNCHCROSSING_TYPE;
    tree->SetBranchAddress("BUNCHCROSSING_TYPE", &BUNCHCROSSING_TYPE);

    // Float_t ALLPVZ;
    // tree->SetBranchAddress("ALLPVZ", &ALLPVZ);

    Double_t Ds_MASS;
    Float_t Ds_PT;
    Float_t Ds_BPVDIRA;
    Double_t Ds_DOCA;
    Double_t Ds_CHI2DOF;
    Float_t Ds_BPVFDCHI2;
    tree->SetBranchAddress("Ds_MASS", &Ds_MASS);
    tree->SetBranchAddress("Ds_PT", &Ds_PT);
    tree->SetBranchAddress("Ds_DOCA", &Ds_DOCA);
    tree->SetBranchAddress("Ds_BPVDIRA", &Ds_BPVDIRA);
    tree->SetBranchAddress("Ds_CHI2DOF", &Ds_CHI2DOF);
    tree->SetBranchAddress("Ds_BPVFDCHI2", &Ds_BPVFDCHI2);


    Float_t Kp_PT;
    Float_t Kp_P;
    Float_t Kp_ETA;
    Double_t Kp_CHI2DOF;
    Float_t Kp_PROBNNGHOST;
    Float_t Kp_BPVIPCHI2;
    Float_t Kp_BPVIP;
    Float_t Kp_PIDK;
    Float_t Kp_PIDP;
    tree->SetBranchAddress("Kp_PT", &Kp_PT);
    tree->SetBranchAddress("Kp_P", &Kp_P);
    tree->SetBranchAddress("Kp_ETA", &Kp_ETA);
    tree->SetBranchAddress("Kp_CHI2DOF", &Kp_CHI2DOF);
    tree->SetBranchAddress("Kp_PROBNNGHOST", &Kp_PROBNNGHOST);
    tree->SetBranchAddress("Kp_BPVIPCHI2", &Kp_BPVIPCHI2);
    tree->SetBranchAddress("Kp_BPVIP", &Kp_BPVIP);
    tree->SetBranchAddress("Kp_PIDK", &Kp_PIDK);
    tree->SetBranchAddress("Kp_PIDP", &Kp_PIDP);

    Float_t Km_PT;
    Float_t Km_P;
    Float_t Km_ETA;
    Double_t Km_CHI2DOF;
    Float_t Km_PROBNNGHOST;
    Float_t Km_BPVIPCHI2;
    Float_t Km_BPVIP;
    Float_t Km_PIDK;
    tree->SetBranchAddress("Km_PT", &Km_PT);
    tree->SetBranchAddress("Km_P", &Km_P);
    tree->SetBranchAddress("Km_ETA", &Km_ETA);
    tree->SetBranchAddress("Km_CHI2DOF", &Km_CHI2DOF);
    tree->SetBranchAddress("Km_PROBNNGHOST", &Km_PROBNNGHOST);
    tree->SetBranchAddress("Km_BPVIPCHI2", &Km_BPVIPCHI2);
    tree->SetBranchAddress("Km_BPVIP", &Km_BPVIP);
    tree->SetBranchAddress("Km_PIDK", &Km_PIDK);

    Float_t pip_PT;
    Float_t pip_P;
    Float_t pip_ETA;
    Double_t pip_CHI2DOF;
    Float_t pip_PROBNNGHOST;
    Float_t pip_BPVIPCHI2;
    Float_t pip_BPVIP;
    Float_t pip_PIDK;
    tree->SetBranchAddress("pip_PT", &pip_PT);
    tree->SetBranchAddress("pip_P", &pip_P);
    tree->SetBranchAddress("pip_ETA", &pip_ETA);
    tree->SetBranchAddress("pip_CHI2DOF", &pip_CHI2DOF);
    tree->SetBranchAddress("pip_PROBNNGHOST", &pip_PROBNNGHOST);
    tree->SetBranchAddress("pip_BPVIPCHI2", &pip_BPVIPCHI2);
    tree->SetBranchAddress("pip_BPVIP", &pip_BPVIP);
    tree->SetBranchAddress("pip_PIDK", &pip_PIDK);


    // Create output file and tree with the same structure
    TFile *fout = nullptr;
    TTree *outTree = nullptr;
    if(makeTree){
        fout = new TFile(Form("%sDsKKpi_DecayFunTuple_selected.root", outputdir.Data()), "RECREATE");
        outTree = new TTree("DecayFunTuple", "Selected D0 candidates");
        
        // Set up output tree branches with same names as input tree
        outTree->Branch("SMOGInjectedGas", &SMOGInjectedGas, "SMOGInjectedGas/I");
        outTree->Branch("nPVs", &nPVs, "nPVs/I");
        outTree->Branch("BUNCHCROSSING_TYPE", &BUNCHCROSSING_TYPE, "BUNCHCROSSING_TYPE/b");
        
        outTree->Branch("Ds_MASS", &Ds_MASS, "Ds_MASS/D");
        outTree->Branch("Ds_PT", &Ds_PT, "Ds_PT/F");
        outTree->Branch("Ds_BPVDIRA", &Ds_BPVDIRA, "Ds_BPVDIRA/F");
        outTree->Branch("Ds_DOCA", &Ds_DOCA, "Ds_DOCA/D");
        outTree->Branch("Ds_CHI2DOF", &Ds_CHI2DOF, "Ds_CHI2DOF/D");
        outTree->Branch("Ds_BPVFDCHI2", &Ds_BPVFDCHI2, "Ds_BPVFDCHI2/F");

        outTree->Branch("Kp_PT", &Kp_PT, "Kp_PT/F");
        outTree->Branch("Kp_P", &Kp_P, "Kp_P/F");
        outTree->Branch("Kp_ETA", &Kp_ETA, "Kp_ETA/F");
        outTree->Branch("Kp_CHI2DOF", &Kp_CHI2DOF, "Kp_CHI2DOF/D");
        outTree->Branch("Kp_PROBNNGHOST", &Kp_PROBNNGHOST, "Kp_PROBNNGHOST/F");
        outTree->Branch("Kp_BPVIPCHI2", &Kp_BPVIPCHI2, "Kp_BPVIPCHI2/F");
        outTree->Branch("Kp_BPVIP", &Kp_BPVIP, "Kp_BPVIP/F");
        outTree->Branch("Kp_PIDK", &Kp_PIDK, "Kp_PIDK/F");
        outTree->Branch("Kp_PIDP", &Kp_PIDP, "Kp_PIDP/F");
        
        outTree->Branch("Km_PT", &Km_PT, "Km_PT/F");
        outTree->Branch("Km_P", &Km_P, "Km_P/F");
        outTree->Branch("Km_ETA", &Km_ETA, "Km_ETA/F");
        outTree->Branch("Km_CHI2DOF", &Km_CHI2DOF, "Km_CHI2DOF/D");
        outTree->Branch("Km_PROBNNGHOST", &Km_PROBNNGHOST, "Km_PROBNNGHOST/F");
        outTree->Branch("Km_BPVIPCHI2", &Km_BPVIPCHI2, "Km_BPVIPCHI2/F");
        outTree->Branch("Km_BPVIP", &Km_BPVIP, "Km_BPVIP/F");
        outTree->Branch("Km_PIDK", &Km_PIDK, "Km_PIDK/F");
        // outTree->Branch("Km_PIDP", &Km_PIDP, "Km_PIDP/F");
        
        outTree->Branch("pip_PT", &pip_PT, "pip_PT/F");
        outTree->Branch("pip_P", &pip_P, "pip_P/F");
        outTree->Branch("pip_ETA", &pip_ETA, "pip_ETA/F");
        outTree->Branch("pip_CHI2DOF", &pip_CHI2DOF, "pip_CHI2DOF/D");
        outTree->Branch("pip_PROBNNGHOST", &pip_PROBNNGHOST, "pip_PROBNNGHOST/F");
        outTree->Branch("pip_BPVIPCHI2", &pip_BPVIPCHI2, "pip_BPVIPCHI2/F");
        outTree->Branch("pip_BPVIP", &pip_BPVIP, "pip_BPVIP/F");
        outTree->Branch("pip_PIDK", &pip_PIDK, "pip_PIDK/F");
    } else {
        fout = new TFile(Form("%sD0Kpi_histograms.root", outputdir.Data()), "RECREATE");
    }

    // Create output histograms
    // TH1F *h_Ds_MASS = new TH1F("h_Ds_MASS", "Ds_MASS", 100, 1786, 1952);
    double Ds_mass = 1968.49;
    TH1F *h_Ds_MASS_pAr = new TH1F("h_Ds_MASS_pAr", "h_Ds_MASS_pAr", 100, Ds_mass-60, Ds_mass+60);
    TH1F *h_Ds_MASS_pH = new TH1F("h_Ds_MASS_pH", "h_Ds_MASS_pH", 100, Ds_mass-60, Ds_mass+60);
    TH1F *h_Ds_MASS_pNe = new TH1F("h_Ds_MASS_pNe", "h_Ds_MASS_pNe", 100, Ds_mass-60, Ds_mass+60);
    TH1F *h_Ds_MASS_pHe = new TH1F("h_Ds_MASS_pHe", "h_Ds_MASS_pHe", 100, Ds_mass-60, Ds_mass+60);
    TH1F *h_Ds_tz = new TH1F("h_Ds_tz", "h_Ds_tz", 150, -0.8, 7);
    h_Ds_tz->GetXaxis()->SetTitle("#it{t}_{z} [ps]");
    //create cut flow histogram
    TH1F *h_cutFlow = new TH1F("h_cutFlow", "Cut Flow Histogram", 15, 0.5, 15.5);
    TH1F *h_nEvents_cutFlow = new TH1F("h_nEvents_cutFlow", "Number of events passing each cut", 15, 0.5, 15.5);

    // Create output variables
    double nEvents_SMOG[10] = {0};
    int nSelectedCandidates = 0;

    // Loop over the entries in the tree and fill the histogram
    Long64_t nEntries = tree->GetEntries();
    cout << "Number of entries: " << nEntries << endl;

    // nEntries = nEntries/50;
    // cout << "Reduced number of entries: " << nEntries << endl;
    // cout << "Reduced number of entries: " << nEntries << endl;
    // cout << "Reduced number of entries: " << nEntries << endl;
    // cout << "Reduced number of entries: " << nEntries << endl;

    Float_t Ds_pNeG_M = 1864.83;
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        //cout a progress bar in percent of all entries
        if(i>0 && nEntries>100 && i%(nEntries/(20))==0) std::cout << "//processed " << 100*(i)/nEntries << "%"  << std::endl;

        //select SMOG gas
        nEvents_SMOG[SMOGInjectedGas]++;
        // if (SMOGInjectedGas != 7) continue;
        // if(BUNCHCROSSING_TYPE != 1) continue;
        int currcut = 1;
        h_nEvents_cutFlow->Fill(currcut);
        h_cutFlow->Fill(currcut);
        currcut++; // 1
        // if(ALLPVZ > -300 || ALLPVZ < -550){ h_cutFlow->Fill(2); continue;} h_nEvents_cutFlow->Fill(2);
        

        //apply selection cuts on daughter particles
        // acceptance cuts
        if (Km_PT < 400 || pip_PT < 400 || Kp_PT < 400){ h_cutFlow->Fill(currcut); continue;} h_nEvents_cutFlow->Fill(currcut); currcut++; // 2
        if (Km_ETA < 2 || Km_ETA > 5 || pip_ETA < 2 || pip_ETA > 5 || Kp_ETA < 2 || Kp_ETA > 5){ h_cutFlow->Fill(currcut); continue;} h_nEvents_cutFlow->Fill(currcut); currcut++; // 3
        if (Km_P/1000 > 100 || Km_P/1000 < 3.2 || pip_P/1000 > 100 || pip_P/1000 < 3.2 || Kp_P/1000 > 100 || Kp_P/1000 < 3.2){ h_cutFlow->Fill(currcut); continue;} h_nEvents_cutFlow->Fill(currcut); currcut++; // 4
        // quality cuts
        if (Km_CHI2DOF > 3 || pip_CHI2DOF > 3 || Kp_CHI2DOF > 3){ h_cutFlow->Fill(currcut); continue;} h_nEvents_cutFlow->Fill(currcut); currcut++; // 5
        //ghost probability cuts
        if (Km_PROBNNGHOST > 0.3 || pip_PROBNNGHOST > 0.3 || Kp_PROBNNGHOST > 0.3){ h_cutFlow->Fill(currcut); continue;} h_nEvents_cutFlow->Fill(currcut); currcut++; // 6
        //impact parameter cuts
        if (Km_BPVIPCHI2 < 16 || pip_BPVIPCHI2 < 16 || Kp_BPVIPCHI2 < 16){ h_cutFlow->Fill(currcut); continue;} h_nEvents_cutFlow->Fill(currcut); currcut++; // 7
        // if (Km_BPVIP > 3 || pip_BPVIP > 3 || Kp_BPVIP > 3){ h_cutFlow->Fill(currcut); continue;} h_nEvents_cutFlow->Fill(currcut); currcut++; // 8
        //PID cuts
        if (Km_PIDK < 5 || pip_PIDK > 0 || Kp_PIDK < 5){ h_cutFlow->Fill(currcut); continue;} h_nEvents_cutFlow->Fill(currcut); currcut++; // 9

        // D0 cuts
        // if (Ds_DOCA > 1){ h_cutFlow->Fill(currcut); continue;} h_nEvents_cutFlow->Fill(currcut); currcut++; // 10
        if (Ds_BPVDIRA < 0.99974689 /*TMath::Cos(0.0225)*/){ h_cutFlow->Fill(currcut); continue;} h_nEvents_cutFlow->Fill(currcut); currcut++; // 11
        if (Ds_CHI2DOF > 6){ h_cutFlow->Fill(currcut); continue;} h_nEvents_cutFlow->Fill(currcut); currcut++; // 12
        if (Ds_BPVFDCHI2 < 50){ h_cutFlow->Fill(currcut); continue;} h_nEvents_cutFlow->Fill(currcut);  currcut++; // 13


        // If the candidate passes all cuts, fill output tree
        if(makeTree)outTree->Fill();
        nSelectedCandidates++;

        if(SMOGInjectedGas == 7) h_Ds_MASS_pAr->Fill(Ds_MASS);
        if(SMOGInjectedGas == 1) h_Ds_MASS_pH->Fill(Ds_MASS);
        if(SMOGInjectedGas == 6) h_Ds_MASS_pNe->Fill(Ds_MASS);
        if(SMOGInjectedGas == 3) h_Ds_MASS_pHe->Fill(Ds_MASS);

    }

    // Print the number of events in each SMOG gas
    for (int i = 0; i < 10; ++i) {
        std::cout << "SMOG gas " << gasNames[i].Data() << " events: " << nEvents_SMOG[i] << std::endl;
    }
    std::cout << "Number of candidates passing all cuts: " << nSelectedCandidates << std::endl;

    gStyle->SetOptStat(0);
    // Draw the output histograms
    TCanvas *c_Ds_MASS = new TCanvas("c_Ds_MASS", "c_Ds_MASS", 800, 600);
    c_Ds_MASS->SetTicks(1, 1);
    h_Ds_MASS_pAr->SetLineColor(kRed);
    h_Ds_MASS_pH->SetLineColor(kBlue);
    h_Ds_MASS_pNe->SetLineColor(kGreen);
    h_Ds_MASS_pHe->SetLineColor(kBlack);
    h_Ds_MASS_pAr->SetMarkerColor(kRed);
    h_Ds_MASS_pH->SetMarkerColor(kBlue);
    h_Ds_MASS_pNe->SetMarkerColor(kGreen);
    h_Ds_MASS_pHe->SetMarkerColor(kBlack);
    h_Ds_MASS_pAr->SetMarkerStyle(20);
    h_Ds_MASS_pH->SetMarkerStyle(20);
    h_Ds_MASS_pNe->SetMarkerStyle(20);
    h_Ds_MASS_pHe->SetMarkerStyle(20);
    h_Ds_MASS_pAr->SetTitle("");
    h_Ds_MASS_pAr->GetXaxis()->SetTitle("m_{D_{s}#rightarrow K^{+}K^{-}#pi^{+}} [MeV]");
    h_Ds_MASS_pAr->GetYaxis()->SetTitle("Candidates");
    h_Ds_MASS_pAr->GetYaxis()->SetRangeUser(0, 1.2*h_Ds_MASS_pAr->GetMaximum());
    h_Ds_MASS_pAr->Draw("pe");
    h_Ds_MASS_pH->Draw("same,pe");
    h_Ds_MASS_pNe->Draw("same,pe");
    h_Ds_MASS_pHe->Draw("same,pe");
    //add a legend
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(h_Ds_MASS_pAr, "#it{p}Ar", "pe");
    leg->AddEntry(h_Ds_MASS_pH, "#it{p}H", "pe");
    leg->AddEntry(h_Ds_MASS_pNe, "#it{p}Ne", "pe");
    leg->AddEntry(h_Ds_MASS_pHe, "#it{p}He", "pe");
    leg->Draw();
    //add label for center of mass energy
    TLatex *label = new TLatex();
    label->SetTextSize(0.04);
    label->DrawLatexNDC(0.15, 0.85, "D_{s} mass");

    c_Ds_MASS->SaveAs(Form("%sDs_MASS.pdf", outputdir.Data()));


    //make the mass plot for pAr in the style of LHCb. This includes a crystal ball + gaussian fit and a bottom panel with the residuals
    makeLHCbStyleMassPlot(h_Ds_MASS_pAr, outputdir, "pAr", "K^{+}K^{-}#pi^{+}", "D_{s}", "Ds");
    makeLHCbStyleMassPlot(h_Ds_MASS_pH, outputdir, "pH", "K^{+}K^{-}#pi^{+}", "D_{s}", "Ds");
    makeLHCbStyleMassPlot(h_Ds_MASS_pNe, outputdir, "pNe", "K^{+}K^{-}#pi^{+}", "D_{s}", "Ds");

    if(makeTree){
        //Draw cut flow histogram
        TCanvas *c_cutFlow = new TCanvas("c_cutFlow", "c_cutFlow", 800, 600);
        c_cutFlow->SetLogy();
        h_cutFlow->Draw();
        c_cutFlow->SaveAs(Form("%scutFlow.pdf", outputdir.Data()));
        //Draw number of events passing each cut
        TCanvas *c_nEvents_cutFlow = new TCanvas("c_nEvents_cutFlow", "c_nEvents_cutFlow", 800, 600);
        c_nEvents_cutFlow->SetLogy();
        h_nEvents_cutFlow->Draw();
        c_nEvents_cutFlow->SaveAs(Form("%snEvents_cutFlow.pdf", outputdir.Data()));
    }

    // Write output tree and close file

    h_Ds_MASS_pAr->Write();
    h_Ds_MASS_pH->Write();
    h_Ds_MASS_pNe->Write();
    if(makeTree){
        h_cutFlow->Write();
        h_nEvents_cutFlow->Write();
        outTree->Write();
    }
    fout->Close();

    std::cout << "Output file created: output_Ds_selected.root" << std::endl;
    
    // Close input file
    file->Close();
}