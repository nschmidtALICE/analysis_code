// Overlay zT distributions for multiple jet pT intervals.
// For each y-bin, produce a plot with all jet pT intervals:
//  - pPb histograms in graduated blues
//  - Pbp histograms in graduated reds
// Assumes histogram naming: unfolded_zT_<jetLow_jetHigh>_<yLow>_<yHigh>_iter<iter>
// Adjust input file paths & iteration tag as needed.

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TString.h>
#include <ctime>
#include <TColor.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <map>
#include <algorithm>

void applyCorrectionsAndPlotFinalResults() {
    //------------------------------------------------------------------
    // Configuration
    //------------------------------------------------------------------
    const TString inputFile_pPb = "/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_zT_pPb_2025-10-14/unfolded_output.root";
    // const TString inputFile_pPb = "/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_zT_2025-10-06-pPb/unfolded_output.root";
    const TString inputFile_Pbp = "/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_zT_Pbp_2025-10-14/unfolded_output.root";
    // const TString inputFile_Pbp = "/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_zT_2025-10-06-Pbp/unfolded_output.root";
    const TString iterUsed      = "6";   // iteration label in histogram names
    const bool    normalize     = true;  // if true, area-normalize each histo (after efficiency correction)
    const bool    useLogY       = false; // enable logY per plot
    const TString outDirBase    = "./final_zT_overlays";
    // compute today's date string
    time_t _t = time(nullptr);
    struct tm _tm = *localtime(&_t);
    char _datebuf[32];
    strftime(_datebuf, sizeof(_datebuf), "%Y-%m-%d", &_tm);
    TString outDir = Form("%s_%s", outDirBase.Data(), _datebuf);
    // Jet-finding efficiency correction using pol2 TF1 fits stored in efficiency file
    const bool    enablePol2EfficiencyCorrection = true; // set false to disable
    const TString effFilePath = "/media/niviths/local/analysis_code/data_analysis/d0_FF/7_efficiency/jetFindingEffMinimal_fullMC_outputs_2025-10-07/jetFindingEffMinimal_fullMC.root";
    // Expected TF1 names: f_pol2_pt5_10, f_pol2_pt10_15, etc. Matching jetMomentumBins below.

    gSystem->Exec( Form( "[ -d %s ] || mkdir -p %s", outDir.Data(), outDir.Data() ) );
    gStyle->SetOptStat(0);
    //------------------------------------------------------------------
    // Binnings (keep in sync with unfolding macro)
    //------------------------------------------------------------------
    const TString jetMomentumBins[] = {"5_10", "10_15", "15_20", "20_30", "30_50"};
    const TString jetMomentumBinsLabels[] = {"5 < #it{p}_{T}^{jet} < 10", "10 < #it{p}_{T}^{jet} < 15", "15 < #it{p}_{T}^{jet} < 20", "20 < #it{p}_{T}^{jet} < 30", "30 < #it{p}_{T}^{jet} < 50"};
    const int     nJet = sizeof( jetMomentumBins ) / sizeof( jetMomentumBins[0] );
    const std::vector<double> yBinBorders = {2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0};
    const int nY = static_cast<int>( yBinBorders.size() ) - 1;

    //------------------------------------------------------------------
    // Styles
    //------------------------------------------------------------------
    int markerStyle_pPb[nJet] = {20, 21, 22, 47, 33};
    int markerStyle_Pbp[nJet] = {24, 25, 26, 46, 27};
    int colorBlue[nJet]       = {kBlue - 9, kBlue - 7, kBlue - 5, kBlue - 3, kBlue - 1};
    int colorRed[nJet]        = {kRed - 9, kRed - 7, kRed - 5, kRed - 3, kRed - 1};

    //------------------------------------------------------------------
    // Open input files
    //------------------------------------------------------------------
    std::unique_ptr<TFile> f_pPb( TFile::Open( inputFile_pPb, "READ" ) );
    std::unique_ptr<TFile> f_Pbp( TFile::Open( inputFile_Pbp, "READ" ) );
    if ( !f_pPb || f_pPb->IsZombie() ) {
        std::cerr << "[ERROR] Cannot open pPb file: " << inputFile_pPb << std::endl;
        return;
    }
    if ( !f_Pbp || f_Pbp->IsZombie() ) {
        std::cerr << "[ERROR] Cannot open Pbp file: " << inputFile_Pbp << std::endl;
        return;
    }
    // Efficiency TF1 file
    std::unique_ptr<TFile> f_eff;
    std::map<std::string, TF1*> pol2Map;
    if ( enablePol2EfficiencyCorrection ) {
        f_eff.reset( TFile::Open( effFilePath, "READ" ) );
        if ( !f_eff || f_eff->IsZombie() ) {
            std::cerr << "[WARN] Cannot open efficiency file: " << effFilePath << " -- disabling pol2 efficiency correction." << std::endl;
        } else {
            // Build mapping of jet bin labels to TF1 objects
            for ( const TString &jb : jetMomentumBins ) {
                int low=0, high=0;
                if ( sscanf( jb.Data(), "%d_%d", &low, &high ) != 2 ) {
                    std::cerr << "[WARN] Could not parse jet bin label '" << jb << "' for efficiency mapping." << std::endl;
                    continue;
                }
                TString fname = Form( "f_pol2_pt%d_%d", low, high );
                TF1 *f = dynamic_cast<TF1*>( f_eff->Get( fname ) );
                if ( !f ) {
                    std::cerr << "[WARN] Missing TF1 '" << fname << "' in efficiency file." << std::endl;
                } else {
                    // Clone to detach so we can close file later
                    TF1 *fClone = (TF1*)f->Clone( Form("%s_clone", fname.Data()) );
                    pol2Map[jb.Data()] = fClone;
                    std::cout << "[INFO] Loaded efficiency fit " << fname << std::endl;
                }
            }
        }
    }
    
    //------------------------------------------------------------------
    // Prepare containers to accumulate mean zT per jet-pT bin as a function of rapidity
    // We'll loop over y-bins below and append the mean from each y-slice into these
    //------------------------------------------------------------------
    std::vector<std::vector<double>> yCentersPerJet(nJet);
    std::vector<std::vector<double>> mean_pPb_perJet(nJet);
    std::vector<std::vector<double>> meanErr_pPb_perJet(nJet);
    std::vector<std::vector<double>> mean_Pbp_perJet(nJet);
    std::vector<std::vector<double>> meanErr_Pbp_perJet(nJet);

    //------------------------------------------------------------------
    // Loop over y-bins; one figure per y-slice
    //------------------------------------------------------------------
    for ( int iy = 0; iy < nY; ++iy ) {
        double yLow  = yBinBorders[iy];
        double yHigh = yBinBorders[iy + 1];

        // Containers for this y-bin
        std::vector<TH1D*> h_pPb, h_Pbp;              // corrected (if enabled)
        std::vector<TH1D*> h_pPb_raw, h_Pbp_raw;      // raw (uncorrected)
    std::vector<std::string> jetLabelsLoaded;     // keep track of which jet bins were successfully loaded
        h_pPb.reserve( nJet );
        h_Pbp.reserve( nJet );
        h_pPb_raw.reserve( nJet );
        h_Pbp_raw.reserve( nJet );

        // Load histograms per jet pT bin
        for ( int ij = 0; ij < nJet; ++ij ) {
            TString hnamepPb = Form( "unfolded_zT_pPb_%s_%.1f_%.1f_iter%s", jetMomentumBins[ij].Data(), yLow, yHigh, iterUsed.Data() );
            TString hnamePbp = Form( "unfolded_zT_Pbp_%s_%.1f_%.1f_iter%s", jetMomentumBins[ij].Data(), yLow, yHigh, iterUsed.Data() );
            TH1D* hp = dynamic_cast<TH1D*>( f_pPb->Get( hnamepPb ) );
            TH1D* hr = dynamic_cast<TH1D*>( f_Pbp->Get( hnamePbp ) );
            if ( !hp ) { std::cerr << "[WARN] Missing pPb histo: " << hnamepPb << std::endl; continue; }
            if ( !hr ) { std::cerr << "[WARN] Missing Pbp histo: " << hnamePbp << std::endl; continue; }

            TH1D* hpC = static_cast<TH1D*>( hp->Clone( Form( "%s_pPb_clone", hp->GetName() ) ) );
            TH1D* hrC = static_cast<TH1D*>( hr->Clone( Form( "%s_Pbp_clone", hr->GetName() ) ) );
            hpC->SetDirectory( nullptr );
            hrC->SetDirectory( nullptr );

            // Raw clones BEFORE efficiency correction or normalization
            TH1D* hpRaw = static_cast<TH1D*>( hp->Clone( Form( "%s_pPb_raw", hp->GetName() ) ) );
            TH1D* hrRaw = static_cast<TH1D*>( hr->Clone( Form( "%s_Pbp_raw", hr->GetName() ) ) );
            hpRaw->SetDirectory( nullptr );
            hrRaw->SetDirectory( nullptr );

            // Apply pol2 efficiency correction (divide by efficiency) BEFORE normalization
            if ( enablePol2EfficiencyCorrection && !pol2Map.empty() ) {
                auto itF = pol2Map.find( jetMomentumBins[ij].Data() );
                if ( itF != pol2Map.end() && itF->second ) {
                    TF1 *fEff = itF->second;
                    int nBins = hpC->GetNbinsX();
                    for ( int b=1; b<=nBins; ++b ) {
                        double zCenter = hpC->GetXaxis()->GetBinCenter( b );
                        double eff = fEff->Eval( zCenter );
                        if ( eff > 0 && eff < 2.0 ) { // sanity guard
                            double c1 = hpC->GetBinContent( b );
                            double e1 = hpC->GetBinError( b );
                            double c2 = hrC->GetBinContent( b );
                            double e2 = hrC->GetBinError( b );
                            // Propagate only statistical bin error (ignore fit parameter uncertainty): scale errors by 1/eff
                            hpC->SetBinContent( b, c1 / eff );
                            hpC->SetBinError( b, e1 / eff );
                            hrC->SetBinContent( b, c2 / eff );
                            hrC->SetBinError( b, e2 / eff );
                        } else if ( eff <= 0 ) {
                            hpC->SetBinContent( b, 0 ); hpC->SetBinError( b, 0 );
                            hrC->SetBinContent( b, 0 ); hrC->SetBinError( b, 0 );
                        }
                    }
                } else {
                    std::cerr << "[WARN] No efficiency TF1 for jet bin '" << jetMomentumBins[ij] << "'; skipping efficiency correction for this bin." << std::endl;
                }
            }

            if ( normalize ) {
                double ip = hpC->Integral( "width" );
                double ir = hrC->Integral( "width" );
                if ( ip > 0 ) hpC->Scale( 1.0 / ip );
                if ( ir > 0 ) hrC->Scale( 1.0 / ir );
                // Normalize raw for shape comparison
                double ipRaw = hpRaw->Integral( "width" );
                double irRaw = hrRaw->Integral( "width" );
                if ( ipRaw > 0 ) hpRaw->Scale( 1.0 / ipRaw );
                if ( irRaw > 0 ) hrRaw->Scale( 1.0 / irRaw );
            }

            hpC->SetMarkerStyle( markerStyle_pPb[ij] );
            hpC->SetMarkerColor( colorBlue[ij] );
            hpC->SetLineColor( colorBlue[ij] );
            hpC->SetLineWidth( 2 );
            hrC->SetMarkerStyle( markerStyle_Pbp[ij] );
            hrC->SetMarkerColor( colorRed[ij] );
            hrC->SetLineColor( colorRed[ij] );
            hrC->SetLineWidth( 2 );
            hpC->SetMarkerSize( 1.5 );
            hrC->SetMarkerSize( 1.5 );

            h_pPb.push_back( hpC );
            h_Pbp.push_back( hrC );
            jetLabelsLoaded.push_back( std::string( jetMomentumBins[ij].Data() ) );
            // Style raw variants (dashed)
            hpRaw->SetLineColor( colorBlue[ij] );
            hpRaw->SetMarkerColor( colorBlue[ij] );
            hpRaw->SetLineStyle( 2 );
            hpRaw->SetLineWidth( 2 );
            hpRaw->SetMarkerStyle( markerStyle_pPb[ij] );
            hpRaw->SetMarkerSize( 1.2 );
            hrRaw->SetLineColor( colorRed[ij] );
            hrRaw->SetMarkerColor( colorRed[ij] );
            hrRaw->SetLineStyle( 2 );
            hrRaw->SetLineWidth( 2 );
            hrRaw->SetMarkerStyle( markerStyle_Pbp[ij] );
            hrRaw->SetMarkerSize( 1.2 );
            h_pPb_raw.push_back( hpRaw );
            h_Pbp_raw.push_back( hrRaw );
        }

        if ( h_pPb.empty() || h_Pbp.empty() ) {
            std::cerr << "[INFO] Skipping y-bin (" << yLow << "," << yHigh << ") due to missing histograms." << std::endl;
            for ( auto* h : h_pPb ) delete h;
            for ( auto* h : h_Pbp ) delete h;
            for ( auto* h : h_pPb_raw ) delete h;
            for ( auto* h : h_Pbp_raw ) delete h;
            continue;
        }

        //------------------------------------------------------------------
        // Accumulate mean zT per jet-pT for this y-slice so we can later plot mean vs rapidity
        //------------------------------------------------------------------
        if ( !h_pPb.empty() ) {
            double yCenter = 0.5 * ( yLow + yHigh );
            for ( size_t ip=0; ip < h_pPb.size(); ++ip ) {
                // determine which jet bin this corresponds to
                std::string label = jetLabelsLoaded[ip];
                int jetIdx = -1;
                for ( int k=0; k<nJet; ++k ) {
                    if ( label == std::string( jetMomentumBins[k].Data() ) ) { jetIdx = k; break; }
                }
                if ( jetIdx < 0 ) continue;
                double mean_pPb = h_pPb[ip]->GetMean();
                double meanErr_pPb = h_pPb[ip]->GetMeanError();
                double mean_Pbp = h_Pbp[ip]->GetMean();
                double meanErr_Pbp = h_Pbp[ip]->GetMeanError();
                yCentersPerJet[jetIdx].push_back( yCenter );
                mean_pPb_perJet[jetIdx].push_back( mean_pPb );
                meanErr_pPb_perJet[jetIdx].push_back( meanErr_pPb );
                mean_Pbp_perJet[jetIdx].push_back( mean_Pbp );
                meanErr_Pbp_perJet[jetIdx].push_back( meanErr_Pbp );
            }
        }

        //------------------------------------------------------------------
        // Determine y-range
        //------------------------------------------------------------------
        double ymax = 0.0;
        for ( auto* h : h_pPb ) ymax = std::max( ymax, h->GetMaximum() );
        for ( auto* h : h_Pbp ) ymax = std::max( ymax, h->GetMaximum() );
        if ( ymax <= 0 ) ymax = 1.0;
        double ymin = 0.0;

        //------------------------------------------------------------------
        // Canvas & frame
        //------------------------------------------------------------------
        TString cname = Form( "c_zT_ybin_%d", iy );
        std::unique_ptr<TCanvas> c( new TCanvas( cname, cname, 900, 800 ) );
        if ( useLogY ) c->SetLogy();

        // Create an empty frame using the first histogram binning
        TH1D* frame = static_cast<TH1D*>( h_pPb.front()->Clone( "frame" ) );
        frame->Reset();
        frame->SetTitle( Form( "z_{T} distributions;z_{T};%s", normalize ? "1/N dN/dz_{T}" : "Entries" ) );
        frame->GetYaxis()->SetRangeUser( useLogY ? ymax * 1e-4 : ymin, ymax * ( useLogY ? 10 : 1.25 ) );
        frame->Draw();

        //------------------------------------------------------------------
        // Draw histograms (pPb first, then Pbp) to avoid reds hiding blues
        //------------------------------------------------------------------
        for ( size_t ij = 0; ij < h_pPb.size(); ++ij ) h_pPb[ij]->Draw( "E1,hist SAME" );
        for ( size_t ij = 0; ij < h_Pbp.size(); ++ij ) h_Pbp[ij]->Draw( "E1,hist SAME" );

        //------------------------------------------------------------------
        // Legend
        //------------------------------------------------------------------
        std::unique_ptr<TLegend> leg( new TLegend( 0.58, 0.60, 0.89, 0.89 ) );
        leg->SetBorderSize( 0 );
        leg->SetFillStyle( 0 );
        leg->SetTextSize( 0.032 );
        for ( int ij = 0; ij < nJet; ++ij ) {
            // Find corresponding loaded index (may have skipped if missing)
            if ( ij >= static_cast<int>( h_pPb.size() ) || ij >= static_cast<int>( h_Pbp.size() ) ) break;
            TString jetLabel = jetMomentumBins[ij];
            jetLabel.ReplaceAll( "_", "-" );
            leg->AddEntry( h_pPb[ij], Form( "pPb %s GeV", jetLabel.Data() ), "lep" );
            leg->AddEntry( h_Pbp[ij], Form( "Pbp %s GeV", jetLabel.Data() ), "lep" );
        }
        leg->Draw();

        //------------------------------------------------------------------
        // Annotation
        //------------------------------------------------------------------
        TLatex lat;
        lat.SetNDC();
        lat.SetTextSize( 0.038 );
        lat.DrawLatex( 0.16, 0.88, Form( "#bf{Preliminary}" ) );
        lat.DrawLatex( 0.16, 0.83, Form( "%.1f < y < %.1f", yLow, yHigh ) );
        lat.DrawLatex( 0.16, 0.78, Form( "Iter %s" , iterUsed.Data() ) );
        if ( enablePol2EfficiencyCorrection && !pol2Map.empty() ) {
            lat.DrawLatex( 0.16, 0.73, "Jet-finding eff. (pol2) corrected" );
        }

        //------------------------------------------------------------------
        // Save outputs
        //------------------------------------------------------------------
        TString foutBase = Form( "%s/zT_overlay_y%.1f_%.1f", outDir.Data(), yLow, yHigh );
        c->SaveAs( foutBase + ".png" );

        // Also save species-separated overlays: pPb only and Pbp only
        // pPb-only
        {
            TString cname_pPb = Form( "c_zT_ybin_%d_pPb", iy );
            std::unique_ptr<TCanvas> c_pPb( new TCanvas( cname_pPb, cname_pPb, 900, 800 ) );
            if ( useLogY ) c_pPb->SetLogy();
            TH1D* frame_pPb = static_cast<TH1D*>( h_pPb.front()->Clone( "frame_pPb" ) );
            frame_pPb->Reset();
            frame_pPb->SetTitle( Form( "z_{T} distributions (pPb only);z_{T};%s", normalize ? "1/N dN/dz_{T}" : "Entries" ) );
            frame_pPb->GetYaxis()->SetRangeUser( useLogY ? ymax * 1e-4 : ymin, ymax * ( useLogY ? 10 : 1.25 ) );
            frame_pPb->Draw();
            for ( size_t ij = 0; ij < h_pPb.size(); ++ij ) h_pPb[ij]->Draw( "E1,hist SAME" );
            std::unique_ptr<TLegend> leg_pPb( new TLegend( 0.58, 0.60, 0.89, 0.89 ) );
            leg_pPb->SetBorderSize( 0 ); leg_pPb->SetFillStyle( 0 ); leg_pPb->SetTextSize( 0.032 );
            for ( int ij = 0; ij < nJet; ++ij ) {
                if ( ij >= static_cast<int>( h_pPb.size() ) ) break;
                TString jetLabel = jetMomentumBins[ij]; jetLabel.ReplaceAll( "_", "-" );
                leg_pPb->AddEntry( h_pPb[ij], Form( "pPb %s GeV", jetLabel.Data() ), "lep" );
            }
            leg_pPb->Draw();
            TLatex lat_p; lat_p.SetNDC(); lat_p.SetTextSize( 0.038 );
            lat_p.DrawLatex( 0.16, 0.88, Form( "#bf{Preliminary}" ) );
            lat_p.DrawLatex( 0.16, 0.83, Form( "%.1f < y < %.1f", yLow, yHigh ) );
            TString fout_pPb = foutBase + "_pPb";
            c_pPb->SaveAs( fout_pPb + ".png" );
            delete frame_pPb;
        }

        // Pbp-only
        {
            TString cname_Pbp = Form( "c_zT_ybin_%d_Pbp", iy );
            std::unique_ptr<TCanvas> c_Pbp( new TCanvas( cname_Pbp, cname_Pbp, 900, 800 ) );
            if ( useLogY ) c_Pbp->SetLogy();
            TH1D* frame_Pbp = static_cast<TH1D*>( h_Pbp.front()->Clone( "frame_Pbp" ) );
            frame_Pbp->Reset();
            frame_Pbp->SetTitle( Form( "z_{T} distributions (Pbp only);z_{T};%s", normalize ? "1/N dN/dz_{T}" : "Entries" ) );
            frame_Pbp->GetYaxis()->SetRangeUser( useLogY ? ymax * 1e-4 : ymin, ymax * ( useLogY ? 10 : 1.25 ) );
            frame_Pbp->Draw();
            for ( size_t ij = 0; ij < h_Pbp.size(); ++ij ) h_Pbp[ij]->Draw( "E1,hist SAME" );
            std::unique_ptr<TLegend> leg_Pbp( new TLegend( 0.58, 0.60, 0.89, 0.89 ) );
            leg_Pbp->SetBorderSize( 0 ); leg_Pbp->SetFillStyle( 0 ); leg_Pbp->SetTextSize( 0.032 );
            for ( int ij = 0; ij < nJet; ++ij ) {
                if ( ij >= static_cast<int>( h_Pbp.size() ) ) break;
                TString jetLabel = jetMomentumBins[ij]; jetLabel.ReplaceAll( "_", "-" );
                leg_Pbp->AddEntry( h_Pbp[ij], Form( "Pbp %s GeV", jetLabel.Data() ), "lep" );
            }
            leg_Pbp->Draw();
            TLatex lat_r; lat_r.SetNDC(); lat_r.SetTextSize( 0.038 );
            lat_r.DrawLatex( 0.16, 0.88, Form( "#bf{Preliminary}" ) );
            lat_r.DrawLatex( 0.16, 0.83, Form( "%.1f < y < %.1f", yLow, yHigh ) );
            TString fout_Pbp = foutBase + "_Pbp";
            c_Pbp->SaveAs( fout_Pbp + ".png" );
            delete frame_Pbp;
        }
        // c->SaveAs( foutBase + ".pdf" );

        // Comparison canvas (raw vs corrected) if efficiency correction active
        if ( enablePol2EfficiencyCorrection && !pol2Map.empty() ) {
            double ymaxCmp = 0.0;
            for ( auto* h : h_pPb )      ymaxCmp = std::max( ymaxCmp, h->GetMaximum() );
            for ( auto* h : h_Pbp )      ymaxCmp = std::max( ymaxCmp, h->GetMaximum() );
            for ( auto* h : h_pPb_raw )  ymaxCmp = std::max( ymaxCmp, h->GetMaximum() );
            for ( auto* h : h_Pbp_raw )  ymaxCmp = std::max( ymaxCmp, h->GetMaximum() );
            if ( ymaxCmp <= 0 ) ymaxCmp = 1.0;
            TString cnameCmp = Form( "c_zT_ybin_%d_rawVsCorr", iy );
            std::unique_ptr<TCanvas> cCmp( new TCanvas( cnameCmp, cnameCmp, 900, 800 ) );
            if ( useLogY ) cCmp->SetLogy();
            TH1D* frameCmp = static_cast<TH1D*>( h_pPb.front()->Clone( "frameCmp" ) );
            frameCmp->Reset();
            frameCmp->SetTitle( Form( "z_{T} distributions: raw vs corrected;z_{T};%s", normalize ? "1/N dN/dz_{T}" : "Entries" ) );
            frameCmp->GetYaxis()->SetRangeUser( useLogY ? ymaxCmp * 1e-4 : 0.0, ymaxCmp * ( useLogY ? 10 : 1.25 ) );
            frameCmp->Draw();
            // Raw first
            for ( size_t ij=0; ij<h_pPb_raw.size(); ++ij ) h_pPb_raw[ij]->Draw( "E1,hist SAME" );
            for ( size_t ij=0; ij<h_Pbp_raw.size(); ++ij ) h_Pbp_raw[ij]->Draw( "E1,hist SAME" );
            // Corrected
            for ( size_t ij=0; ij<h_pPb.size(); ++ij ) h_pPb[ij]->Draw( "E1,hist SAME" );
            for ( size_t ij=0; ij<h_Pbp.size(); ++ij ) h_Pbp[ij]->Draw( "E1,hist SAME" );
            std::unique_ptr<TLegend> legCmp( new TLegend( 0.50, 0.52, 0.89, 0.89 ) );
            legCmp->SetBorderSize(0); legCmp->SetFillStyle(0); legCmp->SetTextSize(0.028);
            for ( size_t ij=0; ij<h_pPb.size(); ++ij ) {
                if ( ij >= h_pPb_raw.size() || ij >= h_Pbp_raw.size() ) break;
                TString jetLabel = jetMomentumBins[ij]; jetLabel.ReplaceAll("_","-");
                legCmp->AddEntry( h_pPb[ij],     Form("pPb %s GeV (corr)", jetLabel.Data()), "l" );
                legCmp->AddEntry( h_pPb_raw[ij], Form("pPb %s GeV (raw)",  jetLabel.Data()), "l" );
                legCmp->AddEntry( h_Pbp[ij],     Form("Pbp %s GeV (corr)", jetLabel.Data()), "l" );
                legCmp->AddEntry( h_Pbp_raw[ij], Form("Pbp %s GeV (raw)",  jetLabel.Data()), "l" );
            }
            legCmp->Draw();
            TLatex latCmp; latCmp.SetNDC(); latCmp.SetTextSize(0.038);
            latCmp.DrawLatex(0.16,0.88,"#bf{Preliminary}");
            latCmp.DrawLatex(0.16,0.83,Form("%.1f < y < %.1f", yLow, yHigh));
            latCmp.DrawLatex(0.16,0.78,Form("Iter %s", iterUsed.Data()));
            latCmp.DrawLatex(0.16,0.73,"Pol2 eff. applied");
            TString foutBaseCmp = Form( "%s/zT_overlay_y%.1f_%.1f_rawVsCorr", outDir.Data(), yLow, yHigh );
            cCmp->SaveAs( foutBaseCmp + ".png" );
            // cCmp->SaveAs( foutBaseCmp + ".pdf" );
            delete frameCmp;
        }

        // Cleanup per y-bin
        delete frame;
        for ( auto* h : h_pPb ) delete h;
        for ( auto* h : h_Pbp ) delete h;
        for ( auto* h : h_pPb_raw ) delete h;
        for ( auto* h : h_Pbp_raw ) delete h;
    }

    // Cleanup TF1 clones
    for ( auto &kv : pol2Map ) delete kv.second;
    std::cout << "[INFO] Finished producing zT overlay plots." << std::endl;

    //------------------------------------------------------------------
    // After looping over rapidity slices: produce mean zT vs rapidity per jet bin
    //------------------------------------------------------------------
    // create a vector to hold all the graphs, indexed by jet bin
    std::vector<TGraphErrors*> vecGraphsPpb(nJet, nullptr), vecGraphsPbp(nJet, nullptr);
    for ( int ij = 0; ij < nJet; ++ij ) {
        if ( yCentersPerJet[ij].empty() ) continue;
        int npts = static_cast<int>( yCentersPerJet[ij].size() );
        TGraphErrors *g_pPb_y = new TGraphErrors();
        TGraphErrors *g_Pbp_y = new TGraphErrors();
        g_pPb_y->SetName( Form("g_pPb_mean_vs_y_%s", jetMomentumBins[ij].Data() ) );
        g_Pbp_y->SetName( Form("g_Pbp_mean_vs_y_%s", jetMomentumBins[ij].Data() ) );
        for ( int ip=0; ip<npts; ++ip ) {
            double y = yCentersPerJet[ij][ip];
            double my = mean_pPb_perJet[ij][ip];
            double mye = meanErr_pPb_perJet[ij][ip];
            double sy = mean_Pbp_perJet[ij][ip];
            double sye = meanErr_Pbp_perJet[ij][ip];
            // find matching y-bin to get bin width for x-error
            double xerr = 0.0;
            for ( size_t ib=0; ib < yBinBorders.size()-1; ++ib ) {
                if ( y >= yBinBorders[ib] && y <= yBinBorders[ib+1] ) {
                    xerr = 0.5 * ( yBinBorders[ib+1] - yBinBorders[ib] );
                    break;
                }
            }
            g_pPb_y->SetPoint( ip, y, my );
            g_pPb_y->SetPointError( ip, xerr, mye );
            g_Pbp_y->SetPoint( ip, y, sy );
            g_Pbp_y->SetPointError( ip, xerr, sye );
        }
        // Apply rapidity shift to Pbp points
        const double rapidityShift = 0.465;
        for ( int ipt = 0; ipt < g_Pbp_y->GetN(); ++ipt ) {
            double xpt, ypt; g_Pbp_y->GetPoint(ipt, xpt, ypt);
            double ex = g_Pbp_y->GetErrorX(ipt);
            double ey = g_Pbp_y->GetErrorY(ipt);
            g_Pbp_y->SetPoint(ipt, xpt + rapidityShift, ypt);
            g_Pbp_y->SetPointError(ipt, ex, ey);
        }

    // Style and draw; store graphs for combined plotting
    g_pPb_y->SetMarkerStyle(20); g_pPb_y->SetMarkerColor(kBlue); g_pPb_y->SetLineColor(kBlue);
    g_Pbp_y->SetMarkerStyle(21); g_Pbp_y->SetMarkerColor(kRed);  g_Pbp_y->SetLineColor(kRed);
    // store in indexed vectors so we can overlay later (keep ownership)
    vecGraphsPpb[ij] = g_pPb_y;
    vecGraphsPbp[ij] = g_Pbp_y;
    TCanvas *cY = new TCanvas( Form("c_mean_vs_y_%s", jetMomentumBins[ij].Data()), Form("Mean zT vs y %s", jetMomentumBins[ij].Data()), 800, 600 );
        // set up frame
        double xmin = 2.4, xmax = 4.8;
        double ymin = 0.11, ymax = 0.79;
        TH2D *frameY = new TH2D( Form("frameY_%s", jetMomentumBins[ij].Data()), "Mean z_{T} vs rapidity; y; mean z_{T}", 100, xmin, xmax, 100, ymin, ymax );
        frameY->Draw();
        g_pPb_y->Draw( "P same" );
        g_Pbp_y->Draw( "P same" );
        TLegend *lgy = new TLegend( 0.62, 0.75, 0.9, 0.9 ); lgy->SetFillStyle(0); lgy->SetBorderSize(0);
        lgy->AddEntry( g_pPb_y, "pPb", "p" ); lgy->AddEntry( g_Pbp_y, "Pbp (shifted)", "p" ); lgy->Draw();
        TLatex t; t.SetNDC(); t.SetTextSize(0.035); t.DrawLatex( 0.16, 0.92, Form( "%s GeV", jetMomentumBins[ij].Data() ) );
        TString foutY = Form( "%s/mean_zT_vs_y_%s", outDir.Data(), jetMomentumBins[ij].Data() );
        cY->SaveAs( foutY + ".png" );
        delete frameY; delete lgy; delete cY; // keep g_pPb_y/g_Pbp_y alive in vecGraphs
    }

    //------------------------------------------------------------------
    // Combined overlay of all mean-vs-y graphs (all jet pT bins)
    //------------------------------------------------------------------
    {
        const double rapidityShift = 0.465;
        double xmin = yBinBorders.front();
        double xmax = yBinBorders.back() + rapidityShift;

        

        TCanvas *cAll = new TCanvas("c_all_mean_vs_y", "Mean zT vs y (all jet bins)", 800, 700);
        cAll->SetTickx(); cAll->SetTicky();
        cAll->SetLeftMargin(0.09); cAll->SetRightMargin(0.01); cAll->SetTopMargin(0.01); cAll->SetBottomMargin(0.08);
        TH2D *frameAll = new TH2D("frameAll", "; #it{y}; <#it{z}_{T}>", 100, xmin-0.2, xmax+0.2, 100, 0.17, 0.89);
        frameAll->Draw();
        TLegend *legAll = new TLegend( 0.115, 0.67, 0.56, 0.84 ); legAll->SetFillStyle(0); legAll->SetBorderSize(0); legAll->SetTextSize(0.03); legAll->SetNColumns(3);
        for ( int ij=0; ij<nJet-1; ++ij ) {
            legAll->AddEntry( (TObject*)nullptr, Form("%s GeV/#it{c}: ", jetMomentumBinsLabels[ij].Data()), "" );
            if ( vecGraphsPpb[ij] ) {
                vecGraphsPpb[ij]->SetMarkerStyle( markerStyle_pPb[ij] ); vecGraphsPpb[ij]->SetMarkerColor( colorBlue[ij] ); vecGraphsPpb[ij]->SetLineColor( colorBlue[ij] ); vecGraphsPpb[ij]->SetMarkerSize(1.5);
                vecGraphsPpb[ij]->Draw("P same");
                legAll->AddEntry( vecGraphsPpb[ij], TString::Format("pPb"), "p" );
            }
            if ( vecGraphsPbp[ij] ) {
                vecGraphsPbp[ij]->SetMarkerStyle( markerStyle_Pbp[ij] ); vecGraphsPbp[ij]->SetMarkerColor( colorRed[ij] ); vecGraphsPbp[ij]->SetLineColor( colorRed[ij] ); vecGraphsPbp[ij]->SetMarkerSize(1.5);
                vecGraphsPbp[ij]->Draw("P same");
                legAll->AddEntry( vecGraphsPbp[ij], TString::Format("Pbp"), "p" );
            }
        }
        legAll->Draw();
        TLatex tax; 
        tax.SetNDC(); 
        tax.SetTextSize(0.035); 
        tax.DrawLatex(0.15, 0.92, "#font[22]{LHCb} #bf{in progress}");
        tax.DrawLatex(0.15, 0.87, "#it{p}Pb, Pb#it{p} #sqrt{s_{NN}} = 8.16 TeV, D^{0}-tagged jets");
        TString foutAll = Form( "%s/mean_zT_vs_y_allJetBins", outDir.Data() );
        cAll->SaveAs( foutAll + ".png" );
        cAll->SaveAs( foutAll + ".pdf" );
        delete frameAll; delete legAll; delete cAll;
    }

}