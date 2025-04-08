#include "../plottingheader.h"
#include "../lhcbStyle.h"

void analyze_D0(bool makeTree = false)
{
    // // Open the .root file
    TFile *file = TFile::Open("/media/niviths/local/analysis_code/data_analysis/lambdab/D0_plots/D0Kpi_DecayFunTuple_selected.root");
    // TFile *file = TFile::Open("/media/niviths/SSD2/lhcb_analysis_SSD/20250318_fullDownloadOfAP/D0Kpi_DecayFunTuple.root");
    if (!file || file->IsZombie())
    {
        std::cerr << "Error opening file" << std::endl;
        return;
    }
    // Get the tree from the file
    TTree *tree;
    file->GetObject("DecayFunTuple", tree);
    if (!tree)
    {
        std::cerr << "Error getting tree from file" << std::endl;
        file->Close();
        return;
    }

    TString outputdir = "D0_plots/";
    system("mkdir -p " + outputdir);

    setLHCbStyle();

    // Set up variables to hold data from the tree
    Int_t SMOGInjectedGas;
    Int_t nPVs;
    unsigned long long int EVENTNUMBER;
    tree->SetBranchAddress("SMOGInjectedGas", &SMOGInjectedGas);
    tree->SetBranchAddress("nPVs", &nPVs);
    tree->SetBranchAddress("EVENTNUMBER", &EVENTNUMBER);
    TString gasNames[] = {"NONE", "HYDROGEN", "DEUTERIUM", "HELIUM", "NITROGEN", "OXYGEN", "NEON", "ARGON", "KRYPTON", "XENON"};

    UChar_t BUNCHCROSSING_TYPE;
    tree->SetBranchAddress("BUNCHCROSSING_TYPE", &BUNCHCROSSING_TYPE);

    Float_t PVZ;
    tree->SetBranchAddress("PVZ", &PVZ);

    Double_t D0_Tz;
    Float_t D0_Y;
    Double_t D0_MASS;
    Float_t D0_PT;
    Float_t D0_BPVDIRA;
    Double_t D0_DOCA;
    Double_t D0_CHI2DOF;
    Float_t D0_BPVFDCHI2, D0_BPVIPCHI2;
    Float_t D0_B_PV_Z, D0_PZ, D0_END_VZ, D0_BPVLTIME;
    tree->SetBranchAddress("D0_Y", &D0_Y);
    tree->SetBranchAddress("D0_MASS", &D0_MASS);
    tree->SetBranchAddress("D0_PT", &D0_PT);
    tree->SetBranchAddress("D0_DOCA", &D0_DOCA);
    tree->SetBranchAddress("D0_BPVDIRA", &D0_BPVDIRA);
    tree->SetBranchAddress("D0_CHI2DOF", &D0_CHI2DOF);
    tree->SetBranchAddress("D0_BPVFDCHI2", &D0_BPVFDCHI2);
    tree->SetBranchAddress("D0_BPVIPCHI2", &D0_BPVIPCHI2);
    tree->SetBranchAddress("D0_B_PV_Z", &D0_B_PV_Z);
    tree->SetBranchAddress("D0_PZ", &D0_PZ);
    tree->SetBranchAddress("D0_END_VZ", &D0_END_VZ);
    tree->SetBranchAddress("D0_BPVLTIME", &D0_BPVLTIME);

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
    double tz, tzError, mass;
    if (makeTree)
    {
        tz = 0;
        tzError = 0;
        mass = 0;
    }
    if (makeTree)
    {
        fout = new TFile(Form("%sD0Kpi_DecayFunTuple_selected.root", outputdir.Data()), "RECREATE");
        outTree = new TTree("DecayFunTuple", "Selected D0 candidates");

        // Set up output tree branches with same names as input tree
        outTree->Branch("SMOGInjectedGas", &SMOGInjectedGas, "SMOGInjectedGas/I");
        outTree->Branch("nPVs", &nPVs, "nPVs/I");
        outTree->Branch("EVENTNUMBER", &EVENTNUMBER, "EVENTNUMBER/I");
        outTree->Branch("BUNCHCROSSING_TYPE", &BUNCHCROSSING_TYPE, "BUNCHCROSSING_TYPE/b");
        outTree->Branch("PVZ", &PVZ, "PVZ/F");

        outTree->Branch("D0_Tz", &D0_Tz, "D0_Tz/D");
        outTree->Branch("D0_Y", &D0_Y, "D0_Y/F");
        outTree->Branch("D0_MASS", &D0_MASS, "D0_MASS/D");
        outTree->Branch("tz", &tz, "tz/D");
        outTree->Branch("tzError", &tzError, "tzError/D");
        outTree->Branch("mass", &mass, "mass/D");
        outTree->Branch("D0_PT", &D0_PT, "D0_PT/F");
        outTree->Branch("D0_BPVDIRA", &D0_BPVDIRA, "D0_BPVDIRA/F");
        outTree->Branch("D0_DOCA", &D0_DOCA, "D0_DOCA/D");
        outTree->Branch("D0_CHI2DOF", &D0_CHI2DOF, "D0_CHI2DOF/D");
        outTree->Branch("D0_BPVFDCHI2", &D0_BPVFDCHI2, "D0_BPVFDCHI2/F");
        outTree->Branch("D0_BPVIPCHI2", &D0_BPVIPCHI2, "D0_BPVIPCHI2/F");
        outTree->Branch("D0_B_PV_Z", &D0_B_PV_Z, "D0_B_PV_Z/F");
        outTree->Branch("D0_PZ", &D0_PZ, "D0_PZ/F");
        outTree->Branch("D0_END_VZ", &D0_END_VZ, "D0_END_VZ/F");
        outTree->Branch("D0_BPVLTIME", &D0_BPVLTIME, "D0_BPVLTIME/F");

        outTree->Branch("Km_PT", &Km_PT, "Km_PT/F");
        outTree->Branch("Km_P", &Km_P, "Km_P/F");
        outTree->Branch("Km_ETA", &Km_ETA, "Km_ETA/F");
        outTree->Branch("Km_CHI2DOF", &Km_CHI2DOF, "Km_CHI2DOF/D");
        outTree->Branch("Km_PROBNNGHOST", &Km_PROBNNGHOST, "Km_PROBNNGHOST/F");
        outTree->Branch("Km_BPVIPCHI2", &Km_BPVIPCHI2, "Km_BPVIPCHI2/F");
        outTree->Branch("Km_BPVIP", &Km_BPVIP, "Km_BPVIP/F");
        outTree->Branch("Km_PIDK", &Km_PIDK, "Km_PIDK/F");

        outTree->Branch("pip_PT", &pip_PT, "pip_PT/F");
        outTree->Branch("pip_P", &pip_P, "pip_P/F");
        outTree->Branch("pip_ETA", &pip_ETA, "pip_ETA/F");
        outTree->Branch("pip_CHI2DOF", &pip_CHI2DOF, "pip_CHI2DOF/D");
        outTree->Branch("pip_PROBNNGHOST", &pip_PROBNNGHOST, "pip_PROBNNGHOST/F");
        outTree->Branch("pip_BPVIPCHI2", &pip_BPVIPCHI2, "pip_BPVIPCHI2/F");
        outTree->Branch("pip_BPVIP", &pip_BPVIP, "pip_BPVIP/F");
        outTree->Branch("pip_PIDK", &pip_PIDK, "pip_PIDK/F");
    }
    else
    {
        fout = new TFile(Form("%sD0Kpi_histograms.root", outputdir.Data()), "RECREATE");
    }

    // Create output histograms
    // TH1F *h_D0_MASS = new TH1F("h_D0_MASS", "D0_MASS", 100, 1786, 1952);
    TH1F *h_D0_MASS_pAr = new TH1F("h_D0_MASS_pAr", "h_D0_MASS_pAr", 100, 1865 - 50, 1865 + 50);
    TH1F *h_D0_MASS_pH = new TH1F("h_D0_MASS_pH", "h_D0_MASS_pH", 100, 1865 - 50, 1865 + 50);
    TH1F *h_D0_MASS_pNe = new TH1F("h_D0_MASS_pNe", "h_D0_MASS_pNe", 100, 1865 - 50, 1865 + 50);
    // make additional mass histograms in 15 rapidity bins from 2.0 to 5.0
    TH1F *h_D0_MASS_pAr_rapidity[15];
    TH1F *h_D0_MASS_pH_rapidity[15];
    TH1F *h_D0_MASS_pNe_rapidity[15];
    for (int i = 0; i < 15; i++)
    {
        h_D0_MASS_pAr_rapidity[i] = new TH1F(Form("h_D0_MASS_pAr_rapidity_%d", i), Form("h_D0_MASS_pAr_rapidity_%d", i), 100, 1865 - 50, 1865 + 50);
        h_D0_MASS_pH_rapidity[i] = new TH1F(Form("h_D0_MASS_pH_rapidity_%d", i), Form("h_D0_MASS_pH_rapidity_%d", i), 100, 1865 - 50, 1865 + 50);
        h_D0_MASS_pNe_rapidity[i] = new TH1F(Form("h_D0_MASS_pNe_rapidity_%d", i), Form("h_D0_MASS_pNe_rapidity_%d", i), 100, 1865 - 50, 1865 + 50);
    }

    TH1F *h_D0_tz = new TH1F("h_D0_tz", "h_D0_tz", 150, -0.8, 7);
    h_D0_tz->GetXaxis()->SetTitle("#it{t}_{z} [ps]");
    // create cut flow histogram
    TH1F *h_cutFlow = new TH1F("h_cutFlow", "Cut Flow Histogram", 15, 0.5, 15.5);
    TH1F *h_nEvents_cutFlow = new TH1F("h_nEvents_cutFlow", "Number of events passing each cut", 15, 0.5, 15.5);
    TString binlabels[15] = {"All (1)", "#it{p}_{T} (2)", "#eta (3)", "#it{p} (4)", "#chi^{2} (5)", "Ghost prob (6)", "BPVIPCHI2 (7)", "BPVIP (8)", "PID (9)", "D0 DOCA (10)", "D0 BPVDIRA (11)", "D0 CHI2DOF (12)", "D0 cuts (13)", "D0 cuts (14)", "D0 cuts (15)"};
    for (int i = 1; i <= 15; i++)
        h_cutFlow->GetXaxis()->SetBinLabel(i, Form("%s", binlabels[i - 1].Data()));
    for (int i = 1; i <= 15; i++)
        h_nEvents_cutFlow->GetXaxis()->SetBinLabel(i, Form("%s", binlabels[i - 1].Data()));
    // make a proper time histogram
    TH1F *h_D0_tz_pAr = new TH1F("h_D0_tz_pAr", "h_D0_tz_pAr", 150, 0, 1);
    TH1F *h_D0_tz_pH = new TH1F("h_D0_tz_pH", "h_D0_tz_pH", 150, 0, 1);
    TH1F *h_D0_tz_pNe = new TH1F("h_D0_tz_pNe", "h_D0_tz_pNe", 150, 0, 1);
    h_D0_tz_pAr->GetXaxis()->SetTitle("#it{t}_{z} [ps]");
    h_D0_tz_pH->GetXaxis()->SetTitle("#it{t}_{z} [ps]");
    h_D0_tz_pNe->GetXaxis()->SetTitle("#it{t}_{z} [ps]");
    // make log(D0_BPVIPCHI2) histogram
    TH1F *h_D0_BPVIPCHI2_pAr = new TH1F("h_D0_BPVIPCHI2_pAr", "h_D0_BPVIPCHI2_pAr", 200, -5, 4);
    TH1F *h_D0_BPVIPCHI2_pH = new TH1F("h_D0_BPVIPCHI2_pH", "h_D0_BPVIPCHI2_pH", 200, -5, 4);
    TH1F *h_D0_BPVIPCHI2_pNe = new TH1F("h_D0_BPVIPCHI2_pNe", "h_D0_BPVIPCHI2_pNe", 200, -5, 4);

    // Create output variables
    double nEvents_SMOG[10] = {0};
    int nSelectedCandidates = 0;

    // Loop over the entries in the tree and fill the histogram
    Long64_t nEntries = tree->GetEntries();
    // nEntries = nEntries/50;
    int scaledownfactor = 1;
    cout << "Number of entries: " << nEntries << endl;
    if (scaledownfactor != 1)
        cout << "Inspecting every " << scaledownfactor << "th entry" << endl;

    // cout << "Reduced number of entries: " << nEntries << endl;
    // cout << "Reduced number of entries: " << nEntries << endl;
    // cout << "Reduced number of entries: " << nEntries << endl;
    // cout << "Reduced number of entries: " << nEntries << endl;

    Float_t D0_pNeG_M = 1864.83;
    unsigned long long int lastEventNumber = -1;
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        if (i > 0 && nEntries > 100 && i % (nEntries / (20)) == 0)
            std::cout << "//processed " << 100 * (i) / nEntries << "%" << std::endl;
        if (scaledownfactor != 1 && i % scaledownfactor != 0)
            continue;
        tree->GetEntry(i);
        // cout a progress bar in percent of all entries

        // Check if the event number is the same as the previous one
        // if (EVENTNUMBER == lastEventNumber) std::cout << "Event number " << EVENTNUMBER << " is duplicated" << std::endl;
        // lastEventNumber = EVENTNUMBER; //std::cout << "Event number " << EVENTNUMBER << std::endl;

        // select SMOG gas
        nEvents_SMOG[SMOGInjectedGas]++;
        // if (SMOGInjectedGas != 7) continue;
        // if(BUNCHCROSSING_TYPE != 1) continue;
        h_nEvents_cutFlow->Fill(1);
        h_cutFlow->Fill(1);
        // if(ALLPVZ > -300 || ALLPVZ < -550){ h_cutFlow->Fill(2); continue;} h_nEvents_cutFlow->Fill(2);

        // apply selection cuts on daughter particles
        //  acceptance cuts
        if (Km_PT < 400 || pip_PT < 400)
        {
            h_cutFlow->Fill(2);
            continue;
        }
        h_nEvents_cutFlow->Fill(2);
        if (Km_ETA < 2 || Km_ETA > 5 || pip_ETA < 2 || pip_ETA > 5)
        {
            h_cutFlow->Fill(3);
            continue;
        }
        h_nEvents_cutFlow->Fill(3);
        if (Km_P / 1000 > 100 || Km_P / 1000 < 3.2 || pip_P / 1000 > 100 || pip_P / 1000 < 3.2)
        {
            h_cutFlow->Fill(4);
            continue;
        }
        h_nEvents_cutFlow->Fill(4);
        // quality cuts
        if (Km_CHI2DOF > 3 || pip_CHI2DOF > 3)
        {
            h_cutFlow->Fill(5);
            continue;
        }
        h_nEvents_cutFlow->Fill(5);
        // ghost probability cuts
        if (Km_PROBNNGHOST > 0.3 || pip_PROBNNGHOST > 0.3)
        {
            h_cutFlow->Fill(6);
            continue;
        }
        h_nEvents_cutFlow->Fill(6);
        // impact parameter cuts
        if (Km_BPVIPCHI2 < 16 || pip_BPVIPCHI2 < 16)
        {
            h_cutFlow->Fill(7);
            continue;
        }
        h_nEvents_cutFlow->Fill(7);
        if (Km_BPVIP > 3 || pip_BPVIP > 3)
        {
            h_cutFlow->Fill(8);
            continue;
        }
        h_nEvents_cutFlow->Fill(8);
        // PID cuts
        if (Km_PIDK < 5 || pip_PIDK > 0)
        {
            h_cutFlow->Fill(9);
            continue;
        }
        h_nEvents_cutFlow->Fill(9);

        // D0 cuts
        if (D0_DOCA > 1)
        {
            h_cutFlow->Fill(10);
            continue;
        }
        h_nEvents_cutFlow->Fill(10);
        if (D0_BPVDIRA < 0.99985036 /*TMath::Cos(0.0173)*/)
        {
            h_cutFlow->Fill(11);
            continue;
        }
        h_nEvents_cutFlow->Fill(11);
        if (D0_CHI2DOF > 10)
        {
            h_cutFlow->Fill(12);
            continue;
        }
        h_nEvents_cutFlow->Fill(12);
        if (D0_BPVFDCHI2 < 49)
        {
            h_cutFlow->Fill(13);
            continue;
        }
        h_nEvents_cutFlow->Fill(13);

        // If the candidate passes all cuts, fill output tree
        if (makeTree)
            outTree->Fill();
        nSelectedCandidates++;

        // calculate t_z
        //  double t_z = (D0_MASS - D0_pNeG_M) * 1.115683 / D0_MASS;
        double t_z = (D0_END_VZ - D0_B_PV_Z) * D0_MASS / D0_PZ;
        D0_Tz = t_z;
        tz = t_z;
        tzError = t_z / 100;
        mass = D0_MASS;

        // Fill the mass histogram for each rapidity bin
        int rapidityBin = (int)((D0_Y - 2.0) / 0.2);
        if (rapidityBin > 14)
            rapidityBin = 14; // make sure the bin is within range
        // std::cout << "rapidityBin: " << rapidityBin << " D0_Y: " << D0_Y << std::endl;

        if (SMOGInjectedGas == 7)
        {
            h_D0_MASS_pAr->Fill(D0_MASS);
            h_D0_MASS_pAr_rapidity[rapidityBin]->Fill(D0_MASS);
            h_D0_tz_pAr->Fill(t_z);
            h_D0_BPVIPCHI2_pAr->Fill(log10(D0_BPVIPCHI2));
        }
        if (SMOGInjectedGas == 1)
        {
            h_D0_MASS_pH->Fill(D0_MASS);
            h_D0_MASS_pH_rapidity[rapidityBin]->Fill(D0_MASS);
            h_D0_tz_pH->Fill(t_z);
            h_D0_BPVIPCHI2_pH->Fill(log10(D0_BPVIPCHI2));
        }
        if (SMOGInjectedGas == 6)
        {
            h_D0_MASS_pNe->Fill(D0_MASS);
            h_D0_MASS_pNe_rapidity[rapidityBin]->Fill(D0_MASS);
            h_D0_tz_pNe->Fill(t_z);
            h_D0_BPVIPCHI2_pNe->Fill(log10(D0_BPVIPCHI2));
        }
    }

    // Print the number of events in each SMOG gas
    for (int i = 0; i < 10; ++i)
    {
        std::cout << "SMOG gas " << gasNames[i].Data() << " events: " << nEvents_SMOG[i] << std::endl;
    }
    std::cout << "Number of candidates passing all cuts: " << nSelectedCandidates << std::endl;

    gStyle->SetOptStat(0);
    // Draw the output histograms
    TCanvas *c_D0_MASS = new TCanvas("c_D0_MASS", "c_D0_MASS", 800, 600);
    h_D0_MASS_pAr->SetLineColor(kRed);
    h_D0_MASS_pH->SetLineColor(kBlue);
    h_D0_MASS_pNe->SetLineColor(kGreen);
    h_D0_MASS_pAr->Draw();
    h_D0_MASS_pH->Draw("same");
    h_D0_MASS_pNe->Draw("same");
    // add a legend
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(h_D0_MASS_pAr, "#it{p}Ar", "l");
    leg->AddEntry(h_D0_MASS_pH, "#it{p}H", "l");
    leg->AddEntry(h_D0_MASS_pNe, "#it{p}Ne", "l");
    leg->Draw();
    // add label for center of mass energy
    TLatex *label = new TLatex();
    label->SetTextSize(0.04);
    label->DrawLatexNDC(0.15, 0.85, "D^{0} mass");

    c_D0_MASS->SaveAs(Form("%sD0_MASS.pdf", outputdir.Data()));

    // //make the mass plot for pAr in the style of LHCb. This includes a crystal ball + gaussian fit and a bottom panel with the residuals
    // // combinedMassAndTzFit(h_D0_MASS_pAr, h_D0_tz_pAr, outputdir, "pAr", "K^{-}#pi^{+}", "D^{0}", "D0");

    // std::vector<double> pAr_results = makeLHCbStyleMassPlot(h_D0_MASS_pAr, outputdir, "pAr", "K^{-}#pi^{+}", "D^{0}", "D0");
    // std::vector<double> pH_results = makeLHCbStyleMassPlot(h_D0_MASS_pH, outputdir, "pH", "K^{-}#pi^{+}", "D^{0}", "D0");
    std::vector<double> pNe_results = makeLHCbStyleMassPlot(h_D0_MASS_pNe, outputdir, "pNe", "K^{-}#pi^{+}", "D^{0}", "D0");
    TGraphErrors *g_pAr = new TGraphErrors();
    TGraphErrors *g_pH = new TGraphErrors();
    TGraphErrors *g_pNe = new TGraphErrors();
    // for(int irap = 0; irap < 15; irap++){
    //     std::vector<double> pAr_rapidity_results = makeLHCbStyleMassPlot(h_D0_MASS_pAr_rapidity[irap], outputdir, Form("pAr_rapidity_%d", irap), "K^{-}#pi^{+}", "D^{0}", "D0");
    //     g_pAr->SetPoint(irap, 2.0 + irap*0.2, pAr_rapidity_results[0]);
    //     g_pAr->SetPointError(irap, 0.1, pAr_rapidity_results[1]);
    //     std::vector<double> pH_rapidity_results = makeLHCbStyleMassPlot(h_D0_MASS_pH_rapidity[irap], outputdir, Form("pH_rapidity_%d", irap), "K^{-}#pi^{+}", "D^{0}", "D0");
    //     g_pH->SetPoint(irap, 2.0 + irap*0.2, pH_rapidity_results[0]);
    //     g_pH->SetPointError(irap, 0.1, pH_rapidity_results[1]);
    //     std::vector<double> pNe_rapidity_results = makeLHCbStyleMassPlot(h_D0_MASS_pNe_rapidity[irap], outputdir, Form("pNe_rapidity_%d", irap), "K^{-}#pi^{+}", "D^{0}", "D0");
    //     g_pNe->SetPoint(irap, 2.0 + irap*0.2, pNe_rapidity_results[0]);
    //     g_pNe->SetPointError(irap, 0.1, pNe_rapidity_results[1]);
    // }

    // double nsig_pAr = pAr_results[0];
    // double nsig_err_pAr = pAr_results[1];
    // double nbkg_pAr = pAr_results[2];
    // double nbkg_err_pAr = pAr_results[3];

    // fitD0Lifetime(h_D0_tz_pAr, nsig_pAr, nsig_err_pAr, nbkg_pAr, nbkg_err_pAr, outputdir, "pAr");
    // makeLHCbStyleMassPlot(h_D0_MASS_pAr, outputdir, "pAr", "K^{-}#pi^{+}", "D^{0}", "D0");
    // makeLHCbStyleMassPlot(h_D0_MASS_pH, outputdir, "pH", "K^{-}#pi^{+}", "D^{0}", "D0");
    // makeLHCbStyleMassPlot(h_D0_MASS_pNe, outputdir, "pNe", "K^{-}#pi^{+}", "D^{0}", "D0");

    // d0MassAndLifetimeFit(file, "D0_fits", "pAr", "K^{-}#pi^{+}", "D^{0}", "D0");
    // d0MassAndLifetimeFit(h_D0_MASS_pAr, h_D0_tz_pAr, "D0_fits", "pAr", "K^{-}#pi^{+}", "D^{0}", "D0");

    // plot tz histograms
    TCanvas *c_D0_tz = new TCanvas("c_D0_tz", "c_D0_tz", 800, 600);
    h_D0_tz_pAr->SetLineColor(kRed);
    h_D0_tz_pH->SetLineColor(kBlue);
    h_D0_tz_pNe->SetLineColor(kGreen);
    h_D0_tz_pAr->Draw();
    h_D0_tz_pH->Draw("same");
    h_D0_tz_pNe->Draw("same");
    // add a legend
    TLegend *leg_tz = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg_tz->AddEntry(h_D0_tz_pAr, "#it{p}Ar", "l");
    leg_tz->AddEntry(h_D0_tz_pH, "#it{p}H", "l");
    leg_tz->AddEntry(h_D0_tz_pNe, "#it{p}Ne", "l");
    leg_tz->Draw();
    // add label for center of mass energy
    TLatex *label_tz = new TLatex();
    label_tz->SetTextSize(0.04);
    label_tz->DrawLatexNDC(0.15, 0.85, "#it{t}_{z}");

    c_D0_tz->SaveAs(Form("%sD0_tz.pdf", outputdir.Data()));

    // make a multi panel plot with all the rapidity bins
    TCanvas *c_D0_MASS_rapidity = new TCanvas("c_D0_MASS_rapidity", "c_D0_MASS_rapidity", 2400, 1800);
    c_D0_MASS_rapidity->Divide(3, 5);
    for (int i = 0; i < 15; i++)
    {
        c_D0_MASS_rapidity->cd(i + 1);
        h_D0_MASS_pAr_rapidity[i]->SetLineColor(kRed);
        h_D0_MASS_pH_rapidity[i]->SetLineColor(kBlue);
        h_D0_MASS_pNe_rapidity[i]->SetLineColor(kGreen);
        h_D0_MASS_pAr_rapidity[i]->Draw();
        h_D0_MASS_pH_rapidity[i]->Draw("same");
        h_D0_MASS_pNe_rapidity[i]->Draw("same");
        // add a legend
        //  int rapidityBin = (int)((D0_Y - 2.0) / 0.2);
        TLegend *leg_rapidity = new TLegend(0.7, 0.65, 0.9, 0.9);
        leg_rapidity->AddEntry((TObject *)0, Form("y = %.1f", 2.0 + i * 0.2), "");
        leg_rapidity->AddEntry(h_D0_MASS_pAr_rapidity[i], "#it{p}Ar", "l");
        leg_rapidity->AddEntry(h_D0_MASS_pH_rapidity[i], "#it{p}H", "l");
        leg_rapidity->AddEntry(h_D0_MASS_pNe_rapidity[i], "#it{p}Ne", "l");
        leg_rapidity->Draw();
    }
    // add label for center of mass energy
    TLatex *label_rapidity = new TLatex();
    label_rapidity->SetTextSize(0.04);
    label_rapidity->DrawLatexNDC(0.15, 0.85, "D^{0} mass");
    c_D0_MASS_rapidity->SaveAs(Form("%sD0_MASS_rapidity.pdf", outputdir.Data()));

    // plot g_pAr, g_pH, g_pNe
    TCanvas *c_g_pAr = new TCanvas("c_g_pAr", "c_g_pAr", 800, 600);
    g_pAr->SetMarkerStyle(20);
    g_pAr->SetMarkerSize(1);
    g_pAr->SetMarkerColor(kRed);
    g_pAr->SetLineColor(kRed);
    g_pNe->SetMarkerStyle(20);
    g_pNe->SetMarkerSize(1);
    g_pNe->SetMarkerColor(kGreen);
    g_pNe->SetLineColor(kGreen);
    g_pH->SetMarkerStyle(20);
    g_pH->SetMarkerSize(1);
    g_pH->SetMarkerColor(kBlue);
    g_pH->SetLineColor(kBlue);

    g_pAr->SetTitle("D^{0} mass vs rapidity");
    g_pAr->GetXaxis()->SetTitle("y");
    g_pAr->GetYaxis()->SetTitle("D^{0} signal counts");
    g_pAr->GetYaxis()->SetTitleOffset(1.5);
    g_pAr->GetXaxis()->SetLimits(2.0, 5.0);
    // g_pAr->GetYaxis()->SetRangeUser(1865-50, 1865+50);
    g_pAr->Draw("AP");
    g_pH->Draw("P,same");
    g_pNe->Draw("P,same");
    // add a legend
    TLegend *leg_g_pAr = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg_g_pAr->SetTextSize(0.04);
    leg_g_pAr->AddEntry(g_pAr, "#it{p}Ar", "p");
    leg_g_pAr->AddEntry(g_pH, "#it{p}H", "p");
    leg_g_pAr->AddEntry(g_pNe, "#it{p}Ne", "p");
    leg_g_pAr->Draw();
    // add label for center of mass energy
    //  TLatex *label_g_pAr = new TLatex();
    //  label_g_pAr->SetTextSize(0.04);
    //  label_g_pAr->DrawLatexNDC(0.15, 0.85, "D^{0} mass vs rapidity");
    c_g_pAr->SaveAs(Form("%sD0_vs_rap.pdf", outputdir.Data()));

    // plot log(D0_BPVIPCHI2) histograms
    TCanvas *c_D0_BPVIPCHI2 = new TCanvas("c_D0_BPVIPCHI2", "c_D0_BPVIPCHI2", 800, 600);
    c_D0_BPVIPCHI2->SetLogy();
    // normalize histograms
    //  h_D0_BPVIPCHI2_pAr->Scale(1./h_D0_BPVIPCHI2_pAr->Integral());
    //  h_D0_BPVIPCHI2_pH->Scale(1./h_D0_BPVIPCHI2_pH->Integral());
    //  h_D0_BPVIPCHI2_pNe->Scale(1./h_D0_BPVIPCHI2_pNe->Integral());
    h_D0_BPVIPCHI2_pAr->SetLineColor(kRed);
    h_D0_BPVIPCHI2_pH->SetLineColor(kBlue);
    h_D0_BPVIPCHI2_pNe->SetLineColor(kGreen);
    h_D0_BPVIPCHI2_pAr->SetMarkerColor(kRed);
    h_D0_BPVIPCHI2_pH->SetMarkerColor(kBlue);
    h_D0_BPVIPCHI2_pNe->SetMarkerColor(kGreen);
    // make all histograms have open circle markers
    h_D0_BPVIPCHI2_pAr->SetMarkerStyle(24);
    h_D0_BPVIPCHI2_pH->SetMarkerStyle(24);
    h_D0_BPVIPCHI2_pNe->SetMarkerStyle(24);
    h_D0_BPVIPCHI2_pAr->GetXaxis()->SetTitle("log(#chi^{2}_{IP}(D^{0}))");
    h_D0_BPVIPCHI2_pH->GetXaxis()->SetTitle("log(#chi^{2}_{IP}(D^{0}))");
    h_D0_BPVIPCHI2_pNe->GetXaxis()->SetTitle("log(#chi^{2}_{IP}(D^{0}))");
    h_D0_BPVIPCHI2_pAr->GetYaxis()->SetTitle("counts");
    h_D0_BPVIPCHI2_pH->GetYaxis()->SetTitle("counts");
    h_D0_BPVIPCHI2_pNe->GetYaxis()->SetTitle("counts");
    // h_D0_BPVIPCHI2_pAr->GetYaxis()->SetTitle("Normalized counts");
    // h_D0_BPVIPCHI2_pH->GetYaxis()->SetTitle("Normalized counts");
    // h_D0_BPVIPCHI2_pNe->GetYaxis()->SetTitle("Normalized counts");
    h_D0_BPVIPCHI2_pAr->Draw();
    h_D0_BPVIPCHI2_pH->Draw("same");
    h_D0_BPVIPCHI2_pNe->Draw("same");
    // add a legend
    TLegend *leg_D0_BPVIPCHI2 = new TLegend(0.2, 0.7, 0.4, 0.9);
    leg_D0_BPVIPCHI2->AddEntry(h_D0_BPVIPCHI2_pAr, "#it{p}Ar", "l");
    leg_D0_BPVIPCHI2->AddEntry(h_D0_BPVIPCHI2_pH, "#it{p}H", "l");
    leg_D0_BPVIPCHI2->AddEntry(h_D0_BPVIPCHI2_pNe, "#it{p}Ne", "l");
    leg_D0_BPVIPCHI2->Draw();
    // add label for center of mass energy
    //  TLatex *label_D0_BPVIPCHI2 = new TLatex();
    //  label_D0_BPVIPCHI2->SetTextSize(0.04);
    //  label_D0_BPVIPCHI2->DrawLatexNDC(0.15, 0.85, "log(D^{0} BPVIPCHI2)");
    c_D0_BPVIPCHI2->SaveAs(Form("%sD0_BPVIPCHI2.pdf", outputdir.Data()));

    // Fit the log(BPVIPCHI2) distributions
    // std::vector<double> pAr_ipchi2_results = fitLogIPChi2Distribution(h_D0_BPVIPCHI2_pAr, outputdir, "pAr", "D^{0}", "D0");
    // std::vector<double> pH_ipchi2_results = fitLogIPChi2Distribution(h_D0_BPVIPCHI2_pH, outputdir, "pH", "D^{0}", "D0");
    // std::vector<double> pNe_ipchi2_results = fitLogIPChi2Distribution(h_D0_BPVIPCHI2_pNe, outputdir, "pNe", "D^{0}", "D0");

    // // Compare fit parameters across collision systems
    // std::cout << "CB Mean comparison:" << std::endl;
    // std::cout << "  pAr: " << pAr_ipchi2_results[0] << " ± " << pAr_ipchi2_results[1] << std::endl;
    // std::cout << "  pH:  " << pH_ipchi2_results[0] << " ± " << pH_ipchi2_results[1] << std::endl;
    // std::cout << "  pNe: " << pNe_ipchi2_results[0] << " ± " << pNe_ipchi2_results[1] << std::endl;

    if (makeTree)
    {
        // Draw cut flow histogram
        TCanvas *c_cutFlow = new TCanvas("c_cutFlow", "c_cutFlow", 800, 600);
        c_cutFlow->SetLogy();
        h_cutFlow->Draw();
        c_cutFlow->SaveAs(Form("%scutFlow.pdf", outputdir.Data()));
        // Draw number of events passing each cut
        TCanvas *c_nEvents_cutFlow = new TCanvas("c_nEvents_cutFlow", "c_nEvents_cutFlow", 800, 600);
        c_nEvents_cutFlow->SetLogy();
        h_nEvents_cutFlow->Draw();
        c_nEvents_cutFlow->SaveAs(Form("%snEvents_cutFlow.pdf", outputdir.Data()));
    }

    // Write output tree and close file

    h_D0_MASS_pAr->Write();
    h_D0_MASS_pH->Write();
    h_D0_MASS_pNe->Write();
    if (makeTree)
    {
        h_cutFlow->Write();
        h_nEvents_cutFlow->Write();
        outTree->Write();
    }
    h_D0_BPVIPCHI2_pAr->Write();
    h_D0_BPVIPCHI2_pH->Write();
    h_D0_BPVIPCHI2_pNe->Write();
    fout->Close();

    std::cout << "Output file created: output_D0_selected.root" << std::endl;

    // Close input file
    file->Close();
}
