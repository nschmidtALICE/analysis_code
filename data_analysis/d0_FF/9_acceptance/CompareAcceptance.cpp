void CompareAcceptance() {
    // Implementation of acceptance comparison
TString inputFile1 = "/media/niviths/local/analysis_code/data_analysis/d0_FF/9_acceptance/D0AcceptanceMap_pPb_53to56_2026-01-15/D0AcceptanceMap_pPb_53to56.root";
TString inputFile2 = "/media/niviths/local/analysis_code/data_analysis/d0_FF/9_acceptance/D0AcceptanceMap_Pbp_2025-10-14/D0AcceptanceMap_Pbp.root";
// TString inputFile2 = "/media/niviths/local/analysis_code/data_analysis/d0_FF/9_acceptance/D0AcceptanceMap_pPb_54_2026-01-06/D0AcceptanceMap_pPb_54.root";

    TFile* inFile1 = TFile::Open(inputFile1, "READ");
    if (!inFile1 || inFile1->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile1 << std::endl;
        return;
    }

    TH2D* hAcceptance = dynamic_cast<TH2D*>(inFile1->Get("hAcceptance"));
    if (!hAcceptance) {
        std::cerr << "Error: Cannot find histogram 'hAcceptance' in file " << inputFile1 << std::endl;
        inFile1->Close();
        return;
    }

    TFile* inFile2 = TFile::Open(inputFile2, "READ");
    if (!inFile2 || inFile2->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile2 << std::endl;
        inFile1->Close();
        return;
    }

    TH2D* hAcceptance2 = dynamic_cast<TH2D*>(inFile2->Get("hAcceptance"));
    if (!hAcceptance2) {
        std::cerr << "Error: Cannot find histogram 'hAcceptance' in file " << inputFile2 << std::endl;
        inFile2->Close();
        inFile1->Close();
        return;
    }

    TH2D* hAcceptanceRatio = dynamic_cast<TH2D*>(hAcceptance->Clone("hAcceptanceRatio"));
    hAcceptanceRatio->SetTitle("#eta;p_{T} [GeV];Acceptance Ratio");
    hAcceptanceRatio->Divide(hAcceptance2);

    TCanvas* cRatio = new TCanvas("cRatio", "Acceptance Ratio", 800, 600);
    hAcceptanceRatio->SetStats(0);
    hAcceptanceRatio->GetZaxis()->SetTitle("Acceptance Ratio");
    hAcceptanceRatio->Draw("COLZ");
    cRatio->SetRightMargin(0.15);
    cRatio->SaveAs("AcceptanceRatio.png");
    cRatio->SaveAs("AcceptanceRatio.pdf");

    // create slices in eta and plot them all on one canvas with different colors
    int nEtaBins = hAcceptance->GetYaxis()->GetNbins();
    TCanvas* cSlices = new TCanvas("cSlices", "Acceptance Ratio Slices", 1000, 700);
    cSlices->SetRightMargin(0.15);
    TLegend* leg = new TLegend(0.75, 0.25, 0.95, 0.92);
    leg->SetBorderSize(0);

    // color palette to cycle through
    int colors[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kCyan, kOrange+7, kViolet, kAzure, kGray+2};
    const int nColors = sizeof(colors)/sizeof(int);

    bool first = true;
    TList* owned = new TList();
    owned->SetOwner(kTRUE);

    for (int i = 1; i <= nEtaBins; ++i) {
        TH1D* hSlice1 = hAcceptance->ProjectionX(Form("hSlice1_etaBin%d", i), i, i);
        //rebin for better visibility
        hSlice1->Rebin(2);
        TH1D* hSlice2 = hAcceptance2->ProjectionX(Form("hSlice2_etaBin%d", i), i, i);
        hSlice2->Rebin(2);
        TH1D* hSliceRatio = dynamic_cast<TH1D*>(hSlice1->Clone(Form("hSliceRatio_etaBin%d", i)));
        hSliceRatio->Divide(hSlice2);
        int col = colors[(i-1) % nColors];
        hSliceRatio->SetLineColor(col);
        hSliceRatio->SetMarkerColor(col);
        hSliceRatio->SetMarkerStyle(20);
        hSliceRatio->SetStats(0);
        hSliceRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
        hSliceRatio->GetXaxis()->SetRangeUser(0, 15);

        if (first) {
            hSliceRatio->SetTitle("Acceptance Ratio;p_{T} [GeV];Acceptance Ratio");
            hSliceRatio->Draw("hist");
            first = false;
        } else {
            hSliceRatio->Draw("hist same");
        }

        leg->AddEntry(hSliceRatio, Form("#eta bin %d", i), "lp");

        // keep ownership so we can clean up later
        owned->Add(hSlice1);
        owned->Add(hSlice2);
        owned->Add(hSliceRatio);
    }

    leg->Draw();
    cSlices->SaveAs("AcceptanceRatio_slices.png");
    cSlices->SaveAs("AcceptanceRatio_slices.pdf");

    delete cSlices;
    delete leg;
    delete owned;

    // TString outputFile = "output.root";

}