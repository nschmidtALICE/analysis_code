#!/usr/bin/env python3
# filepath: /media/niviths/local/analysis_code/data_analysis/d0_FF/d0-code/signalFit/MassFitterScript_TestReduce.py
"""
D0 Meson Analysis Framework

This script performs mass and lifetime fits for D0 mesons in different kinematic regions.
It can be used for both data and Monte Carlo analysis.
"""

import os
import sys
import ROOT
import Fitter as fit
import Plotter as plt
import array
import numpy
import argparse
import math

ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

class FitSpectraObject:
    """
    Main class for fitting mass spectra of D0 mesons in different kinematic regions.
    It handles the creation of datasets, fitting, and visualization of the results.
    """
    
    def __init__(self, resonance, pt_range, is_binned, is_mc, z_bins, is_zt_observable):
        """
        Initialize the fit object with configuration parameters.
        
        Args:
            resonance: String identifying the particle to analyze ("D0")
            pt_range: Tuple of (min_pt, max_pt) defining the jet pT range
            is_binned: Boolean flag for binned vs unbinned fitting
            is_mc: Boolean flag for MC vs data analysis
            z_bins: Array of bin edges for z or dR binning
            is_zt_observable: Boolean flag for zT vs dR binning
        """
        # Core configuration parameters
        self.binned = is_binned
        self.isMC = is_mc
        self.jetPt = pt_range
        self.zBins = z_bins
        self.nzTBins = len(self.zBins) - 1
        self.zTObservable = is_zt_observable
        
        # Configure the resonance specific parameters
        self._configure_resonance(resonance)
        
        # Configure file paths
        self._configure_file_paths()
        
    def _configure_resonance(self, resonance):
        """Set up parameters specific to the resonance type"""
        if "D0" in resonance:
            self.dictKey = "D0"
            self.sideBandLimits = (1.840, 1.890)  # D0 mass sideband region
        else:
            # Default to D0 if not specified
            self.dictKey = "D0"
            self.sideBandLimits = (1.840, 1.890)
            
    def _configure_file_paths(self):
        """Set up input and output file paths"""
        # Add MC/data and observable tags for file naming
        mc_tag = "MC" if self.isMC else ""
        obs_tag = "zT" if self.zTObservable else "dR"
        
        # Base path for data files
        base_path = "/media/niviths/local/analysis_code/data_analysis/d0_FF/test/"
        
        # Select appropriate input file based on particle and data/MC
        if self.isMC:
            self.ntupleFile = "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_MC_output_D0FF_filterV1.root"
        else:
            # self.ntupleFile = "/media/niviths/local/download_ganga/2025_05_Pbp/1103152_filterV1.root"
            self.ntupleFile = "/media/niviths/SSD2/lhcb_analysis_SSD/20250501_Pbp_data/Pbp_data_filterV1.root"
            # self.ntupleFile = "/media/niviths/SSD2/lhcb_analysis_SSD/20250430_ganga_output/ntuple_filterV1.root"
            # self.ntupleFile = "/media/niviths/local/download_ganga/2025_05_Pbp/202505Pbp_filterV1.root"
            # self.ntupleFile = "/media/niviths/local/analysis_code/data_analysis/d0_FF/Output01May2025_MeasuredData_D0NewPF_III_filterV1.root"
            # self.ntupleFile = "/media/niviths/SSD2/lhcb_analysis_SSD/20250430_ganga_output/ntuple_filterV1.root"
            # self.ntupleFile = "/media/niviths/local/analysis_code/data_analysis/d0_FF/outputs/3/ntuple_JetFrag-1099031_filterV1.root"
            # self.ntupleFile = "/media/niviths/local/analysis_code/data_analysis/d0_FF/Output29Apr2025_MeasuredData_D0NewPF_III_filterV1.root"
        
        # Set up output directories
        #set a different directory for data and MC
        if self.isMC:
            output_dir = f"/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/{self.dictKey}_MC"
        else:
            output_dir = f"/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/{self.dictKey}_FF"
        self.outfilePath = f"{output_dir}/{self.jetPt[0]}_{self.jetPt[1]}/"
        self.fOutDataName = f"{output_dir}/FitParametersUnBinned{self.dictKey}{obs_tag}_{self.jetPt[0]}_{self.jetPt[1]}"
        self.fOutDataNameB = f"{output_dir}/BinnedSpectra{self.dictKey}{obs_tag}_{self.jetPt[0]}_{self.jetPt[1]}"
        
        # Print configuration information
        print(f"Create fit object for {'binned' if self.binned else 'unbinned'} fit")
        print(f"Output directory: {output_dir}")
        print(f"Output data file: {self.fOutDataName}")
                
    def startFitting(self):
        """Main method to perform mass and lifetime fits"""
        print("\n===== STARTING FITTING PROCESS =====")
        
        # Initialize arrays for results
        print("Initializing result arrays...")
        self._initialize_result_arrays()
        
        # Calculate bin centers and widths
        print("Calculating bin properties...")
        zBinsCent, zBinsWidth = self._calculate_bin_properties()
        print(f"Bin centers: {zBinsCent}")
        print(f"Bin widths: {zBinsWidth}")
        
        # Create fitter object and prepare datasets
        print("\nCreating fitter object...")
        fitter = self._create_fitter()
        
        print("\nPreparing master dataset...")
        try:
            data_master = self._prepare_master_dataset(fitter)
            print(f"Master dataset created with {data_master.numEntries()} entries")
        except Exception as e:
            print(f"ERROR creating master dataset: {str(e)}")
            import traceback
            traceback.print_exc()
            return
        
        # Process correction factors if requested
        print("\nProcessing correction factors...")
        self._process_correction_factors(fitter, data_master)
        
        # If running just for correction factors, return early
        if hasattr(self, 'correction_only') and self.correction_only:
            print("Correction-only mode, returning early")
            return
        
        print("\n===== STARTING BIN-BY-BIN FITTING =====")
        # Process fit for each bin
        self._process_fits_by_bin(fitter, data_master, zBinsCent, zBinsWidth)
        
    def _initialize_result_arrays(self):
        """Initialize arrays to store fit results"""
        # Number of parameters to store for each fit value (value, error, min, max, etc.)
        self.num_fit_items = 5
        
        # Create arrays for mass fit results
        self.FitMRes_SYield = numpy.zeros((self.nzTBins, self.num_fit_items))
        self.FitMRes_BYield = numpy.zeros((self.nzTBins, self.num_fit_items))
        self.FitMRes_Mean = numpy.zeros((self.nzTBins, self.num_fit_items))
        self.FitMRes_Sig1 = numpy.zeros((self.nzTBins, self.num_fit_items))
        self.FitMRes_deltaSig = numpy.zeros((self.nzTBins, self.num_fit_items))
        self.FitMRes_Sig2 = numpy.zeros((self.nzTBins, self.num_fit_items))
        self.FitMRes_alpha = numpy.zeros((self.nzTBins, self.num_fit_items))
        self.FitMRes_n = numpy.zeros((self.nzTBins, self.num_fit_items))
        self.FitMRes_CBFrac = numpy.zeros((self.nzTBins, self.num_fit_items))
        self.FitMRes_pol0 = numpy.zeros((self.nzTBins, self.num_fit_items))
        self.FitMRes_pol1 = numpy.zeros((self.nzTBins, self.num_fit_items))
        self.FitMRes_SYieldLim = numpy.zeros((self.nzTBins, self.num_fit_items))
        self.FitMRes_BYieldLim = numpy.zeros((self.nzTBins, self.num_fit_items))
        self.FitMRes_SYieldSG = numpy.zeros((self.nzTBins, self.num_fit_items))
        self.FitMRes_SYieldDCB = numpy.zeros((self.nzTBins, self.num_fit_items))
        
        # Create arrays for IP chi2 fit results - updated for double Gaussian parameters
        self.FitIPRes_SYield = numpy.zeros((self.nzTBins, self.num_fit_items))
        # Replace MeanP with MeanP1 and add MeanP2 (double Gaussian means)
        self.FitIPRes_MeanP1 = numpy.zeros((self.nzTBins, self.num_fit_items))  # First Gaussian mean
        self.FitIPRes_MeanP2 = numpy.zeros((self.nzTBins, self.num_fit_items))  # Second Gaussian mean
        self.FitIPRes_MeanNP = numpy.zeros((self.nzTBins, self.num_fit_items))  # Non-prompt mean
        # Replace SigmaP with SigmaP1 and add SigmaP2 (double Gaussian widths)
        self.FitIPRes_SigmaP1 = numpy.zeros((self.nzTBins, self.num_fit_items))  # First Gaussian sigma
        self.FitIPRes_SigmaP2 = numpy.zeros((self.nzTBins, self.num_fit_items))  # Second Gaussian sigma
        self.FitIPRes_GaussFrac = numpy.zeros((self.nzTBins, self.num_fit_items))  # Mixing fraction
        self.FitIPRes_SigmaNP = numpy.zeros((self.nzTBins, self.num_fit_items))  # Non-prompt sigma
        self.FitIPRes_PromptFrac = numpy.zeros((self.nzTBins, self.num_fit_items))  # Prompt fraction
        self.FitIPRes_PromptYield = numpy.zeros((self.nzTBins, self.num_fit_items))  # Prompt yield
        self.FitIPRes_NonPromptYield = numpy.zeros((self.nzTBins, self.num_fit_items))  # Non-prompt yield
        
        # Create arrays for histograms
        self.ipchi2HistoArray = [None] * self.nzTBins
        self.massHistoArray = [None] * self.nzTBins
        
    def _calculate_bin_properties(self):
        """Calculate bin centers and widths"""
        zBinsCent = [0.0] * self.nzTBins
        zBinsWidth = [0.0] * self.nzTBins
        
        for bin in range(self.nzTBins):
            zBinsCent[bin] = 0.5 * (self.zBins[bin] + self.zBins[bin + 1])
            zBinsWidth[bin] = 0.5 * (self.zBins[bin + 1] - self.zBins[bin])
            
        return zBinsCent, zBinsWidth
        
    def _create_fitter(self):
        print("Create and initialize the Fitter object")
        fitter = fit.Fitter(
            rootFileName=self.ntupleFile,
            resonance=self.dictKey,
            nBins=self.nzTBins,
            isBinned=self.binned,
            isMC=self.isMC,
            outfilePath=self.outfilePath,
            update=True
        )
        return fitter
        
    def _prepare_master_dataset(self, fitter):
        """Prepare the master dataset with appropriate cuts"""
        # Get mass range from the fitter's dictionary
        res_sig = fitter.massDict[self.dictKey]
        mass_range = list(res_sig["mass_range"])
        mass_cut = (mass_range[0], mass_range[1])
        
        # Create fiducial cut string and master dataset
        fid_cut = eval(f'fitter.fiducialCutString(jetPt=self.jetPt, tagMass=mass_cut)')
        data_master = fitter.createDataSet(
            self.dictKey,
            "Inclusive data filtered",
            fidCutString=fid_cut,
            isMass=True,
            corrVer=-1
        )
        
        return data_master
        
    def _create_roofit_variables(self):
        """Define RooRealVars for the analysis"""
        # Define variables with appropriate ranges for D0
        tagMass = ROOT.RooRealVar("tagMass", "tagMass", 1.800, 1.950)
        
        # Replace lifetime with IP chi2 variables
        tag_ipchi2 = ROOT.RooRealVar("tag_ip_chi2", "tag_ip_chi2", 0, 10000)
        log_tag_ipchi2 = ROOT.RooRealVar("log_tag_ipchi2", "log(tag_ip_chi2)", -3, 10)
        
        # Other variables remain unchanged
        dRValue = ROOT.RooRealVar("tagJetdR", "tagJetdR", 0.0, 1.0)
        tagZ = ROOT.RooRealVar("tagZ", "tagZ", 0.0, 1.01)
        
        # D0 mass window for signal region
        mass_tag_measured = ROOT.RooRealVar("tagMass", "tagMass", 1.845, 1.885)
        
        # Efficiency weight
        effWeight_3 = ROOT.RooRealVar("EffWeight_3", "EffWeight_3", -1, 2)
        
        return {
            'tagMass': tagMass,
            'tag_ipchi2': tag_ipchi2,
            'log_tag_ipchi2': log_tag_ipchi2,
            'dRValue': dRValue,
            'tagZ': tagZ,
            'mass_tag_measured': mass_tag_measured,
            'effWeight_3': effWeight_3
        }
        
    def _process_correction_factors(self, fitter, data_master):
        """Process correction factors for efficiency studies"""
        # Set correction_only flag to indicate we're only processing correction factors
        self.correction_only = False
        
        # Create RooFit variables
        variables = self._create_roofit_variables()
        
        # Number of correction factors to process
        num_corr_factors = 8
        corr_factor_list = [None] * num_corr_factors
        
        # Create arrays for canvas objects
        canvas = [None] * num_corr_factors
        canvas_norm = [None] * num_corr_factors
        canvas_div = [None] * num_corr_factors
        
        # Define grid layout
        rows = 2
        columns = int(self.nzTBins / rows)
        
        # Arrays for storing yields and errors
        yield_array = numpy.zeros((num_corr_factors, self.nzTBins))
        yield_array_err = numpy.zeros((num_corr_factors, self.nzTBins))
        yield_array_rnd = numpy.zeros((num_corr_factors, self.nzTBins))
        yield_array_err_rnd = numpy.zeros((num_corr_factors, self.nzTBins))
        
        # Set up legends for plots
        legend_list = [None] * num_corr_factors
        for i in range(num_corr_factors):
            legend_list[i] = [None] * self.nzTBins
        
        # Set up canvases and legends
        self._setup_canvases(canvas, canvas_norm, canvas_div, legend_list, columns, rows)
        
        # Create different correction datasets
        fid_cut = eval(f'fitter.fiducialCutString(jetPt=self.jetPt, tagMass=list(fitter.massDict[self.dictKey]["mass_range"]))')
        
        # DOCA cuts dataset
        data_corr_doca = fitter.createDataSet(
            self.dictKey, "DOCA_Cuts", 
            fidCutString=fid_cut, isMass=True, corrVer=-2
        )
        
        # Random selection dataset
        data_corr_rnd = fitter.createDataSet(
            self.dictKey, "EffCorr_Rnd",
            fidCutString=fid_cut, isMass=True, corrVer=0
        )
        
        # Create datasets for all other correction types
        data_corr_1 = fitter.createDataSet(
            self.dictKey, "EffCorr_totalEff",
            fidCutString=fid_cut, isMass=True, corrVer=1  # Total efficiency and acceptance
        )
        
        data_corr_2 = fitter.createDataSet(
            self.dictKey, "EffCorr_pionReco",
            fidCutString=fid_cut, isMass=True, corrVer=2  # Pion reconstruction efficiency
        )
        
        data_corr_3 = fitter.createDataSet(
            self.dictKey, "EffCorr_pionSel",
            fidCutString=fid_cut, isMass=True, corrVer=3  # Pion selection efficiency
        )
        
        data_corr_4 = fitter.createDataSet(
            self.dictKey, "EffCorr_muonReco",
            fidCutString=fid_cut, isMass=True, corrVer=4  # Muon reconstruction efficiency
        )
        
        data_corr_5 = fitter.createDataSet(
            self.dictKey, "EffCorr_stripLine",
            fidCutString=fid_cut, isMass=True, corrVer=5  # Stripping line efficiency
        )
        
        data_corr_6 = fitter.createDataSet(
            self.dictKey, "EffCorr_trigg",
            fidCutString=fid_cut, isMass=True, corrVer=6  # Trigger efficiency
        )
        
        # Extract random efficiency correction
        plotter = plt.Plotter(
            resonance=self.dictKey, bin=0, range="z_{T}=0-1",
            basepath=self.outfilePath, binned=self.binned
        )
        
        # Process each zT bin
        for bin_idx in range(self.nzTBins):
            print(f"Processing correction factors for bin {bin_idx}/{self.nzTBins-1}")
            bin_range = f"z_{{T}}={self.zBins[bin_idx]:.2f}-{self.zBins[bin_idx+1]:.2f}"
            
            # Get data in this bin
            if self.zTObservable:
                bin_cut = f"tagZ>{self.zBins[bin_idx]} && tagZ<{self.zBins[bin_idx+1]}"
            else:
                bin_cut = f"tagJetdR>{self.zBins[bin_idx]} && tagJetdR<{self.zBins[bin_idx+1]}"
            
            # Reduce datasets for this bin
            bin_data_master = data_master.reduce(bin_cut)
            bin_data_corr_1 = data_corr_1.reduce(bin_cut)
            bin_data_corr_2 = data_corr_2.reduce(bin_cut)
            bin_data_corr_3 = data_corr_3.reduce(bin_cut)
            bin_data_corr_4 = data_corr_4.reduce(bin_cut)
            bin_data_corr_5 = data_corr_5.reduce(bin_cut)
            bin_data_corr_6 = data_corr_6.reduce(bin_cut)
            bin_data_corr_doca = data_corr_doca.reduce(bin_cut)
            
            # Extract correction factors - random efficiency
            yield_val_rnd, error_val_rnd = plotter.extractCorParamRnd(
                "Corr_Rnd", variables['mass_tag_measured'], 
                legend_list[0], data_master, data_corr_rnd
            )
            
            # Extract other correction factors with proper weight selections
            # Total efficiency (Type1)
            data_sel = bin_data_master.reduce("EffWeight_0 > 0 && EffWeight_1 > 0 && EffWeight_2 > 0 && EffWeight_3 > 0 && EffWeight_4 > 0")
            plotter.extractCorParam("Corr_totalEff", bin_idx, 1, canvas, canvas_norm, canvas_div, 
                                  yield_array, yield_array_err, legend_list[1], 
                                  variables['mass_tag_measured'], data_sel, bin_data_corr_1)
            
            # Pion reco efficiency (Type2)
            data_sel = bin_data_master.reduce("EffWeight_0 > 0")
            plotter.extractCorParam("Corr_pionReco", bin_idx, 2, canvas, canvas_norm, canvas_div, 
                                  yield_array, yield_array_err, legend_list[2], 
                                  variables['mass_tag_measured'], data_sel, bin_data_corr_2)
            
            # Pion selection efficiency (Type3)
            data_sel = bin_data_master.reduce("EffWeight_1 > 0")
            plotter.extractCorParam("Corr_pionSel", bin_idx, 3, canvas, canvas_norm, canvas_div, 
                                  yield_array, yield_array_err, legend_list[3], 
                                  variables['mass_tag_measured'], data_sel, bin_data_corr_3)
            
            # Muon reco efficiency (Type4)
            data_sel = bin_data_master.reduce("EffWeight_2 > 0")
            plotter.extractCorParam("Corr_muonReco", bin_idx, 4, canvas, canvas_norm, canvas_div, 
                                  yield_array, yield_array_err, legend_list[4], 
                                  variables['mass_tag_measured'], data_sel, bin_data_corr_4)
            
            # Stripping line efficiency (Type5)
            data_sel = bin_data_master.reduce("EffWeight_3 > 0")
            plotter.extractCorParam("Corr_stripLine", bin_idx, 5, canvas, canvas_norm, canvas_div, 
                                  yield_array, yield_array_err, legend_list[5], 
                                  variables['mass_tag_measured'], data_sel, bin_data_corr_5)
            
            # Trigger efficiency (Type6)
            data_sel = bin_data_master.reduce("EffWeight_4 > 0")
            plotter.extractCorParam("Corr_trigg", bin_idx, 6, canvas, canvas_norm, canvas_div, 
                                  yield_array, yield_array_err, legend_list[6], 
                                  variables['mass_tag_measured'], data_sel, bin_data_corr_6)
            
            # DOCA cuts (Type7)
            plotter.extractCorParam("Corr_DOCA", bin_idx, 7, canvas, canvas_norm, canvas_div, 
                              yield_array, yield_array_err, legend_list[7], 
                              variables['mass_tag_measured'], bin_data_corr_doca, bin_data_master)
        
        # Calculate bin centers and widths
        bin_centers = [0.5 * (self.zBins[i] + self.zBins[i+1]) for i in range(self.nzTBins)]
        bin_widths = [0.5 * (self.zBins[i+1] - self.zBins[i]) for i in range(self.nzTBins)]
        
        # Create correction factor graphs for each type
        for corr_type in range(num_corr_factors):
            if corr_type == 0:
                # For Type0 (random selection), use the global yield_val_rnd value for all bins
                for i in range(self.nzTBins):
                    yield_array_rnd[corr_type][i] = yield_val_rnd
                    yield_array_err_rnd[corr_type][i] = error_val_rnd
                    
                corr_factor_list[corr_type] = ROOT.TGraphErrors(
                    self.nzTBins,
                    numpy.array(bin_centers),
                    yield_array_rnd[corr_type],
                    numpy.array(bin_widths),
                    yield_array_err_rnd[corr_type]
                )
                corr_factor_list[corr_type].SetName("CorrHist_TypeRnd")
            else:
                # Create graphs for Types 1-7
                corr_factor_list[corr_type] = ROOT.TGraphErrors(
                    self.nzTBins,
                    numpy.array(bin_centers),
                    yield_array[corr_type],
                    numpy.array(bin_widths),
                    yield_array_err[corr_type]
                )
                corr_factor_list[corr_type].SetName(f"CorrHist_Type{corr_type}")
        
        # Write all correction factors to ROOT file
        print(f"Writing correction factors to {self.fOutDataName}_CorrFact.root")
        f_out = ROOT.TFile(f"{self.fOutDataName}_CorrFact.root", "RECREATE")
        for item in corr_factor_list:
            if item:
                print(f"Writing correction histogram: {item.GetName()}")
                item.Write()
        f_out.Close()
        
    def _setup_canvases(self, canvas, canvas_norm, canvas_div, legend_list, columns, rows):
        """Set up ROOT canvases and legends for visualization"""
        for i in range(len(canvas)):
            # Create main canvas
            canvas[i] = ROOT.TCanvas(
                f"cMulti{self.jetPt[0]}_{self.jetPt[1]}_{i}", 
                f"cMultiCorr{self.jetPt[0]}_{self.jetPt[1]}", 800*4, 400*4
            )
            canvas[i].Divide(columns, rows)
            
            # Create normalized canvas
            canvas_norm[i] = ROOT.TCanvas(
                f"cMultiNorm{self.jetPt[0]}_{self.jetPt[1]}_{i}", 
                f"cMultiCorrNorm{self.jetPt[0]}_{self.jetPt[1]}", 800*4, 400*4
            )
            canvas_norm[i].Divide(columns, rows)
            
            # Create ratio canvas
            canvas_div[i] = ROOT.TCanvas(
                f"cMultiDiv{self.jetPt[0]}_{self.jetPt[1]}_{i}", 
                f"cMultiCorrDiv{self.jetPt[0]}_{self.jetPt[1]}", 800*4, 400*4
            )
            canvas_div[i].Divide(columns, rows)
            
            # Create legends for each bin
            for bin in range(self.nzTBins):
                legend_list[i][bin] = ROOT.TLegend(0.15, 0.65, 0.25, 0.88)
                legend_list[i][bin].SetTextFont(42)
                legend_list[i][bin].SetBorderSize(0)
                legend_list[i][bin].SetFillStyle(0)
                legend_list[i][bin].SetFillColor(0)
                legend_list[i][bin].SetMargin(0.25)
                legend_list[i][bin].SetTextSize(0.03)
                
    def _process_fits_by_bin(self, fitter, data_master, bin_centers, bin_widths):
        """Process mass and IP chi2 fits for each bin"""
        print("\nProcessing fits bin by bin...")
        
        # Define RooFit variables
        print("Creating RooFit variables...")
        variables = self._create_roofit_variables()
        
        # Process each bin
        for bin_idx in range(self.nzTBins):
            print(f"\n----- Processing bin {bin_idx}/{self.nzTBins-1} -----")
            # Reset dictionary for fresh fit parameters
            fitter.initDictionary()
            
            # Get bin range and label
            z_bin = (self.zBins[bin_idx], self.zBins[bin_idx + 1])
            bin_label = self._get_bin_label(bin_idx)
            print(f"Bin range: {z_bin}, Label: {bin_label}")
            
            # Process the bin data
            print(f"Preparing data for bin {bin_idx}...")
            try:
                bin_data = self._prepare_bin_data(data_master, bin_idx, variables)
                print(f"  Bin data created with {bin_data.numEntries()} entries")
            except Exception as e:
                print(f"ERROR preparing bin data: {str(e)}")
                import traceback
                traceback.print_exc()
                continue
            
            # Mass fit
            print(f"Performing mass fit for bin {bin_idx}...")
            try:
                self._perform_mass_fit(fitter, bin_data, bin_idx, bin_label)
                print(f"  Mass fit completed for bin {bin_idx}")
            except Exception as e:
                print(f"ERROR in mass fit: {str(e)}")
                import traceback
                traceback.print_exc()
                continue
            
            # IP chi2 fit 
            print(f"Performing IP chi2 fit for bin {bin_idx}...")
            try:
                self._perform_ipchi2_fit(fitter, bin_data, bin_idx, bin_label)
                print(f"  IP chi2 fit completed for bin {bin_idx}")
            except Exception as e:
                print(f"ERROR in IP chi2 fit: {str(e)}")
                import traceback
                traceback.print_exc()
                continue
        
        # Create and save plots with fit results
        self._create_result_plots(bin_centers, bin_widths)
        
    def _get_bin_label(self, bin_idx):
        """Get label for the current bin"""
        if self.zTObservable:
            return f"z_{{T}}={self.zBins[bin_idx]:.1f}-{self.zBins[bin_idx+1]:.1f}"
        else:
            return f"#Delta R={self.zBins[bin_idx]:.2f}-{self.zBins[bin_idx+1]:.2f}"
            
    def _prepare_bin_data(self, data_master, bin_idx, variables):
        """Prepare data for the current bin"""
        # Create a RooArgSet with all necessary variables including tagMass
        arg_set = ROOT.RooArgSet()
        
        # Always include mass variable needed for plotting
        arg_set.add(variables['tagMass'])
        
        # Define cut string based on observable type
        if self.zTObservable:
            cut_string = f"tagZ>{self.zBins[bin_idx]} && tagZ<{self.zBins[bin_idx+1]}"
            # Add tagZ variable used in the cut
            arg_set.add(variables['tagZ'])
        else:
            cut_string = f"tagJetdR>{self.zBins[bin_idx]} && tagJetdR<{self.zBins[bin_idx+1]}"
            # Add dRValue variable used in the cut
            arg_set.add(variables['dRValue'])
        
        # Add IP chi2 variables that are needed for IP chi2 fits
        if 'tag_ipchi2' in variables:
            arg_set.add(variables['tag_ipchi2'])
        if 'log_tag_ipchi2' in variables:
            arg_set.add(variables['log_tag_ipchi2'])
            
        # Log the operation for debugging
        print(f"  Reducing dataset with cut: {cut_string}")
        # print(f"  Variables included: {[arg_set.at(i).GetName() for i in range(arg_set.getSize())]}")
        
        # Reduce the dataset with the complete set of variables
        reduced_data = data_master.reduce(arg_set, cut_string)
        
        # Verify the reduced dataset has the required variables
        if reduced_data.numEntries() > 0:
            if not reduced_data.get(0).find('tagMass'):
                print("WARNING: Reduced dataset is missing tagMass variable!")
                
        return reduced_data
        
    def _perform_mass_fit(self, fitter, bin_data, bin_idx, bin_label):
        """Perform mass fit for the current bin"""
        print(f"  Starting mass fit for bin {bin_idx}")

        # Check output directories exist
        mass_fit_dir = f"{self.outfilePath}MassFits_{self.zTObservable and 'zT' or 'dR'}/"
        if not os.path.exists(mass_fit_dir):
            print(f"  Creating output directory: {mass_fit_dir}")
            os.makedirs(mass_fit_dir, exist_ok=True)
        
        # Single Gaussian fit
        print("  Performing Single Gaussian fit...")
        try:
            histo_sg, param_sg, param_err_sg = fitter.massFit(
                self.dictKey, bin_data, fitTypeName="SGauss", 
                bin=bin_idx, zRange=bin_label
            )
            print("    SG fit completed successfully")
        except Exception as e:
            print(f"    ERROR in SG fit: {str(e)}")
            import traceback
            traceback.print_exc()
            return
        
        # Double Gaussian fit
        print("  Performing Double Gaussian fit...")
        try:
            fitter.massFit(
                self.dictKey, bin_data, fitTypeName="DGauss", 
                bin=bin_idx, zRange=bin_label
            )
            print("    DG fit completed successfully")
        except Exception as e:
            print(f"    ERROR in DG fit: {str(e)}")
            import traceback
            traceback.print_exc()
            return
        
        # Double Crystal Ball fit
        print("  Performing Double Crystal Ball fit...")
        try:
            histo_dcb, param_dcb, param_err_dcb = fitter.massFit(
                self.dictKey, bin_data, fitTypeName="DCB", 
                bin=bin_idx, zRange=bin_label
            )
            print("    DCB fit completed successfully")
        except Exception as e:
            print(f"    ERROR in DCB fit: {str(e)}")
            import traceback
            traceback.print_exc()
            return
        
        # Fix alpha and n parameters from MC if needed
        if self.nzTBins > 1 and not self.isMC:
            fitter.fixNAlphaValue(self.dictKey, self.jetPt[0])
            
        # Final DCB fit with fixed parameters
        print("  Performing final DCB fit with fixed parameters...")
        try:
            histo, param, param_err = fitter.massFit(
                self.dictKey, bin_data, fitTypeName="DCBFixed", 
                bin=bin_idx, zRange=bin_label
            )
            print("    Final DCB fit completed successfully")
        except Exception as e:
            print(f"    ERROR in final DCB fit: {str(e)}")
            import traceback
            traceback.print_exc()
            return
        
        # Store results from the appropriate fit
        final_param = param if not self.isMC else param_dcb
        final_param_err = param_err if not self.isMC else param_err_dcb
        
        # Store mass histogram
        self.massHistoArray[bin_idx] = histo
        
        # Store fit parameters
        self._store_mass_fit_parameters(bin_idx, final_param, final_param_err)
        
        # Store SG and DCB yields for comparison
        self.FitMRes_SYieldSG[bin_idx, 0] = param_sg[0]
        self.FitMRes_SYieldSG[bin_idx, 1] = param_err_sg[0]
        self.FitMRes_SYieldDCB[bin_idx, 0] = param_dcb[0]
        self.FitMRes_SYieldDCB[bin_idx, 1] = param_err_dcb[0]
        
    def _load_mass_parameters(self, fitter, bin_idx):
        """Load mass fit parameters from file"""
        z_bin = (self.zBins[bin_idx], self.zBins[bin_idx + 1])
        param, param_err = fitter.loadExternalParameters(
            self.dictKey, self.fOutDataName, bin_idx, z_bin
        )
        
        # Store loaded parameters
        self._store_mass_fit_parameters(bin_idx, param, param_err)
        
    def _store_mass_fit_parameters(self, bin_idx, param, param_err):
        """Store mass fit parameters in result arrays"""
        # Basic fit results
        self.FitMRes_SYield[bin_idx, 0] = param[0]
        self.FitMRes_BYield[bin_idx, 0] = param[1]
        self.FitMRes_Mean[bin_idx, 0] = param[2]
        self.FitMRes_Sig1[bin_idx, 0] = param[3]
        self.FitMRes_deltaSig[bin_idx, 0] = param[4]
        self.FitMRes_Sig2[bin_idx, 0] = self.FitMRes_Sig1[bin_idx, 0] * param[4]
        self.FitMRes_alpha[bin_idx, 0] = param[5]
        self.FitMRes_n[bin_idx, 0] = param[6]
        self.FitMRes_CBFrac[bin_idx, 0] = param[7]
        self.FitMRes_pol0[bin_idx, 0] = param[8]
        self.FitMRes_pol1[bin_idx, 0] = param[9]
        self.FitMRes_SYieldLim[bin_idx, 0] = param[10]
        self.FitMRes_BYieldLim[bin_idx, 0] = param[11]
        
        # Error values
        self.FitMRes_SYield[bin_idx, 1] = param_err[0]
        self.FitMRes_BYield[bin_idx, 1] = param_err[1]
        self.FitMRes_Mean[bin_idx, 1] = param_err[2]
        self.FitMRes_Sig1[bin_idx, 1] = param_err[3]
        self.FitMRes_deltaSig[bin_idx, 1] = param_err[4]
        # Propagate error for Sig2
        self.FitMRes_Sig2[bin_idx, 1] = self._propagate_error(
            self.FitMRes_Sig1[bin_idx, 0], param[4], 
            param_err[3], param_err[4], 1
        )
        self.FitMRes_alpha[bin_idx, 1] = param_err[5]
        self.FitMRes_n[bin_idx, 1] = param_err[6]
        self.FitMRes_CBFrac[bin_idx, 1] = param_err[7]
        self.FitMRes_pol0[bin_idx, 1] = param_err[8]
        self.FitMRes_pol1[bin_idx, 1] = param_err[9]
        
    def _perform_ipchi2_fit(self, fitter, bin_data, bin_idx, bin_label):
        """Perform IP chi2 fit for the current bin"""
        # Get signal region from fitter
        res_sig = fitter.massDict[self.dictKey]
        sig_region = res_sig["signal_region"]
        
        # Update prompt fraction based on previous fit or initialization
        if hasattr(fitter, 'ipchi2Dict') and hasattr(fitter.ipchi2Dict[f"Sig{self.dictKey}"], "prompt_frac"):
            prompt_frac = fitter.ipchi2Dict[f"Sig{self.dictKey}"]["prompt_frac"]
            print(f"  Using prompt fraction from dictionary: {prompt_frac[0]} (range: {prompt_frac[1]}-{prompt_frac[2]})")
        
        # Update signal yield limits based on mass fit
        fitter.updateSigYield(self.dictKey)
        
        # Prepare datasets for signal and background regions
        sig_data, bkg_data = self._prepare_ipchi2_datasets(bin_data, sig_region)
        
        # Fit the IP chi2 distribution
        histo, param, param_err = fitter.ipchi2Fit(
            self.dictKey, data=sig_data, background=bkg_data, 
            figKey="BKGincluded", bin=bin_idx, zRange=bin_label
        )
        
        # Store histogram
        self.ipchi2HistoArray[bin_idx] = histo
        
        # Store IP chi2 fit parameters
        self._store_ipchi2_fit_parameters(bin_idx, param, param_err, fitter)
        
        # Print main results
        print(f"  IP chi2 fit results:")
        print(f"    Prompt fraction: {param[5]:.3f} ± {param_err[5]:.3f}")
        print(f"    Prompt yield: {param[6]:.1f}")
        print(f"    Non-prompt yield: {param[7]:.1f}")
        
    def _prepare_ipchi2_datasets(self, bin_data, sig_region):
        """Prepare datasets for IP chi2 fit with signal and background regions"""
        # Signal region
        sig_cut = f"tagMass>{sig_region[0]} && tagMass<{sig_region[1]}"
        sig_data = bin_data.reduce(sig_cut)
        sig_data.SetName("Sig")
        
        # Left side background
        left_cut = f"tagMass<{self.sideBandLimits[0]}"
        bkg_left = bin_data.reduce(left_cut)
        bkg_left.SetName("BKG1")
        
        # Right side background
        right_cut = f"tagMass>{self.sideBandLimits[1]}"
        bkg_right = bin_data.reduce(right_cut)
        bkg_right.SetName("BKG2")
        
        # Combine backgrounds
        bkg_data = bkg_left.Clone("BKG")
        bkg_data.append(bkg_right)
        
        # Check if dataset has IP chi2 variable
        ip_var = sig_data.get().find("tag_ip_chi2")
        log_ip_var = sig_data.get().find("log_tag_ipchi2")
        
        if not ip_var and not log_ip_var:
            print("WARNING: IP chi2 variable not found in dataset, adding calculated column")
            # Create custom column with log(IP chi2) in both datasets
            sig_data = self._add_log_ipchi2_column(sig_data)
            bkg_data = self._add_log_ipchi2_column(bkg_data)
        
        return sig_data, bkg_data
        
    def _create_result_plots(self, bin_centers, bin_widths):
        """Create and save plots with fit results"""
        # Calculate yields for prompt and non-prompt components from IP chi2 fits
        prompt_yield = [self.FitIPRes_PromptYield[i, 0] for i in range(self.nzTBins)]
        nonprompt_yield = [self.FitIPRes_NonPromptYield[i, 0] for i in range(self.nzTBins)]
        
        # Create fragmentation function graphs
        incl_frag_func = ROOT.TGraphErrors(
            self.nzTBins, numpy.array(bin_centers), numpy.array(self.FitMRes_SYield[:, 0]),
            numpy.array(bin_widths), numpy.array([0.0] * self.nzTBins)
        )
        prompt_frag_func = ROOT.TGraphErrors(
            self.nzTBins, numpy.array(bin_centers), numpy.array(prompt_yield),
            numpy.array(bin_widths), numpy.array([0.0] * self.nzTBins)
        )
        nonprompt_frag_func = ROOT.TGraphErrors(
            self.nzTBins, numpy.array(bin_centers), numpy.array(nonprompt_yield),
            numpy.array(bin_widths), numpy.array([0.0] * self.nzTBins)
        )
        
        # Create parameter graphs
        mass_graphs = self._create_parameter_graphs(bin_centers, bin_widths)
        ipchi2_graphs = self._create_ipchi2_parameter_graphs(bin_centers, bin_widths)
        
        # Save results to ROOT file
        self._save_results_to_file(
            incl_frag_func, nonprompt_frag_func, prompt_frag_func,
            mass_graphs, ipchi2_graphs
        )
        
    def _create_parameter_graphs(self, bin_centers, bin_widths):
        """Create graphs for mass fit parameters"""
        # Create graphs for each parameter
        graphs = {
            'SYield': create_graphs("FitMSYield", self.nzTBins, bin_centers, self.FitMRes_SYield, bin_widths),
            'BYield': create_graphs("FitMBYield", self.nzTBins, bin_centers, self.FitMRes_BYield, bin_widths),
            'Mean': create_graphs("FitMMean", self.nzTBins, bin_centers, self.FitMRes_Mean, bin_widths),
            'Sig1': create_graphs("FitMSig1", self.nzTBins, bin_centers, self.FitMRes_Sig1, bin_widths),
            'DeltaSig': create_graphs("FitMDeltaSig", self.nzTBins, bin_centers, self.FitMRes_deltaSig, bin_widths),
            'Sig2': create_graphs("FitMSig2", self.nzTBins, bin_centers, self.FitMRes_Sig2, bin_widths),
            'Alpha': create_graphs("FitMAlpha", self.nzTBins, bin_centers, self.FitMRes_alpha, bin_widths),
            'N': create_graphs("FitMN", self.nzTBins, bin_centers, self.FitMRes_n, bin_widths),
            'CBFrac': create_graphs("FitMCBFrac", self.nzTBins, bin_centers, self.FitMRes_CBFrac, bin_widths),
            'Pol0': create_graphs("FitMPol0", self.nzTBins, bin_centers, self.FitMRes_pol0, bin_widths),
            'Pol1': create_graphs("FitMPol1", self.nzTBins, bin_centers, self.FitMRes_pol1, bin_widths),
            'SYieldLim': create_graphs("FitMSYieldLim", self.nzTBins, bin_centers, self.FitMRes_SYieldLim, bin_widths),
            'BYieldLim': create_graphs("FitMBYieldLim", self.nzTBins, bin_centers, self.FitMRes_BYieldLim, bin_widths),
            'SYieldSG': create_graphs("FitMSYieldSG", self.nzTBins, bin_centers, self.FitMRes_SYieldSG, bin_widths),
            'SYieldDCB': create_graphs("FitMSYieldDCB", self.nzTBins, bin_centers, self.FitMRes_SYieldDCB, bin_widths),
        }
        
        # Create range graphs
        range_graphs = {
            'MeanR': create_graphs_asymm_err("FitMMeanRange", self.nzTBins, bin_centers, self.FitMRes_Mean, bin_widths),
            'Sig1R': create_graphs_asymm_err("FitMSig1Range", self.nzTBins, bin_centers, self.FitMRes_Sig1, bin_widths),
            'DeltaSigR': create_graphs_asymm_err("FitMDeltaSigRange", self.nzTBins, bin_centers, self.FitMRes_deltaSig, bin_widths),
            'Sig2R': create_graphs_asymm_err("FitMSig2Range", self.nzTBins, bin_centers, self.FitMRes_Sig2, bin_widths),
            'AlphaR': create_graphs_asymm_err("FitMAlphaRange", self.nzTBins, bin_centers, self.FitMRes_alpha, bin_widths),
            'NR': create_graphs_asymm_err("FitMNRange", self.nzTBins, bin_centers, self.FitMRes_n, bin_widths),
            'CBFracR': create_graphs_asymm_err("FitMCBFracRange", self.nzTBins, bin_centers, self.FitMRes_CBFrac, bin_widths)
        }
        
        return {'graphs': graphs, 'range_graphs': range_graphs}
        
    def _create_ipchi2_parameter_graphs(self, bin_centers, bin_widths):
        """Create graphs for IP chi2 fit parameters with single Bukin model"""
        # Create graphs for each parameter
        graphs = {
            'SYield': create_graphs("FitIPSYield", self.nzTBins, bin_centers, self.FitIPRes_SYield, bin_widths),
            # Prompt Bukin parameters
            'XpPrompt': create_graphs("FitIPXpPrompt", self.nzTBins, bin_centers, self.FitIPRes_XpPrompt, bin_widths),
            'SigmaPrompt': create_graphs("FitIPSigmaPrompt", self.nzTBins, bin_centers, self.FitIPRes_SigmaPrompt, bin_widths),
            'XiPrompt': create_graphs("FitIPXiPrompt", self.nzTBins, bin_centers, self.FitIPRes_XiPrompt, bin_widths),
            # Non-prompt Bukin parameters
            'XpNonprompt': create_graphs("FitIPXpNonprompt", self.nzTBins, bin_centers, self.FitIPRes_XpNonprompt, bin_widths),
            'SigmaNonprompt': create_graphs("FitIPSigmaNonprompt", self.nzTBins, bin_centers, self.FitIPRes_SigmaNonprompt, bin_widths),
            'XiNonprompt': create_graphs("FitIPXiNonprompt", self.nzTBins, bin_centers, self.FitIPRes_XiNonprompt, bin_widths),
            # Fractions and yields
            'PromptFrac': create_graphs("FitIPPromptFrac", self.nzTBins, bin_centers, self.FitIPRes_PromptFrac, bin_widths),
            'PromptYield': create_graphs("FitIPPromptYield", self.nzTBins, bin_centers, self.FitIPRes_PromptYield, bin_widths),
            'NonPromptYield': create_graphs("FitIPNonPromptYield", self.nzTBins, bin_centers, self.FitIPRes_NonPromptYield, bin_widths),
        }
        
        # Create range graphs
        range_graphs = {
            'XpPromptR': create_graphs_asymm_err("FitIPXpPromptRange", self.nzTBins, bin_centers, self.FitIPRes_XpPrompt, bin_widths),
            'SigmaPromptR': create_graphs_asymm_err("FitIPSigmaPromptRange", self.nzTBins, bin_centers, self.FitIPRes_SigmaPrompt, bin_widths),
            'XiPromptR': create_graphs_asymm_err("FitIPXiPromptRange", self.nzTBins, bin_centers, self.FitIPRes_XiPrompt, bin_widths),
            'XpNonpromptR': create_graphs_asymm_err("FitIPXpNonpromptRange", self.nzTBins, bin_centers, self.FitIPRes_XpNonprompt, bin_widths),
            'SigmaNonpromptR': create_graphs_asymm_err("FitIPSigmaNonpromptRange", self.nzTBins, bin_centers, self.FitIPRes_SigmaNonprompt, bin_widths),
            'XiNonpromptR': create_graphs_asymm_err("FitIPXiNonpromptRange", self.nzTBins, bin_centers, self.FitIPRes_XiNonprompt, bin_widths),
            'PromptFracR': create_graphs_asymm_err("FitIPPromptFracRange", self.nzTBins, bin_centers, self.FitIPRes_PromptFrac, bin_widths)
        }
        
        return {'graphs': graphs, 'range_graphs': range_graphs}
        
    def _save_results_to_file(self, incl_func, nonprompt_func, prompt_func, mass_graphs, ipchi2_graphs):
        """Save fit results to ROOT files"""
        # Create output file for parameters
        f_out_data = ROOT.TFile(f"{self.fOutDataName}.root", "RECREATE")
        
        # Save fragmentation functions
        incl_func.SetName("ginclFragFunc")
        incl_func.Write()
        nonprompt_func.SetName("nonpromptFragFunc")
        nonprompt_func.Write()
        prompt_func.SetName("promptFragFunc")
        prompt_func.Write()
        
        # Save mass parameter graphs
        for graph_set in mass_graphs['graphs'].values():
            for i in range(len(graph_set)):
                if graph_set[i]:
                    graph_set[i].Write()
                    
        # Save IP chi2 parameter graphs
        for graph_set in ipchi2_graphs['graphs'].values():
            for i in range(len(graph_set)):
                if graph_set[i]:
                    graph_set[i].Write()
                    
        # Save range graphs
        for graph in mass_graphs['range_graphs'].values():
            graph.Write()
            
        for graph in ipchi2_graphs['range_graphs'].values():
            graph.Write()
            
        f_out_data.Close()
        
        # Create output file for histograms
        f_out_histo = ROOT.TFile(f"{self.fOutDataNameB}.root", "RECREATE")
        
        # Save fragmentation functions
        incl_func.SetName("ginclFragFunc")
        incl_func.Write()
        nonprompt_func.SetName("nonpromptFragFunc")
        nonprompt_func.Write()
        prompt_func.SetName("promptFragFunc")
        prompt_func.Write()
        
        # Save histograms if not loading parameters
        for i in range(self.nzTBins):
            # Replace timeHistoArray with ipchi2HistoArray
            if hasattr(self, 'ipchi2HistoArray') and self.ipchi2HistoArray[i]:
                self.ipchi2HistoArray[i].SetName(f"hIPChi2_{i}")
                self.ipchi2HistoArray[i].Write()
            if hasattr(self, 'massHistoArray') and self.massHistoArray[i]:
                self.massHistoArray[i].SetName(f"hMassSpectr_{i}")
                self.massHistoArray[i].Write()
                    
        f_out_histo.Close()
        
    @staticmethod
    def _propagate_error(factor_a, factor_b, error_a, error_b, error_type=0):
        """Propagate error for derived quantities"""
        # A/B error propagation
        if error_type == 0:
            if factor_b > 0:
                error = math.pow(factor_a/factor_b, 2) * (
                    math.pow(error_a/factor_a, 2) + math.pow(error_b/factor_b, 2)
                )
            else:
                error = 0
        # A*B error propagation
        elif error_type == 1:
            if factor_a > 0 and factor_b > 0:
                error = math.pow(factor_a*factor_b, 2) * (
                    math.pow(error_a/factor_a, 2) + math.pow(error_b/factor_b, 2)
                )
            else:
                error = 0
        # A*B (B has no error)
        elif error_type == 2:
            error = math.pow(error_a * factor_b, 2)
        else:
            print("Error propagation type not implemented!")
            return 0
            
        return math.sqrt(error)
        
    def _add_log_ipchi2_column(self, dataset):
        """Add log(IP chi2) column to dataset if it doesn't exist"""
        if not dataset:
            return None
        
        # Check if we need to add the column
        if dataset.get().find("log_tag_ipchi2"):
            return dataset  # Already has the column
        
        # Find IP chi2 variable or use proxy
        ip_var = dataset.get().find("tag_ip_chi2")
        if not ip_var:
            print("ERROR: Cannot find tag_ip_chi2 variable in dataset")
            return dataset
        
        # Create new dataset with log transform
        log_range = (-3, 10)  # Default range for log(IP chi2)
        
        # Create new RooRealVar for log(IP chi2)
        log_ipchi2 = ROOT.RooRealVar("log_tag_ipchi2", "log(tag_ip_chi2)", log_range[0], log_range[1])
        
        # Create variables list for new dataset
        vars_list = ROOT.RooArgSet(dataset.get())
        vars_list.add(log_ipchi2)
        
        # Create new dataset
        new_dataset = ROOT.RooDataSet(f"{dataset.GetName()}_with_log", 
                                     f"{dataset.GetTitle()} with log(IP chi2)",
                                     vars_list)
        
        # Fill new dataset with transformed values
        for i in range(dataset.numEntries()):
            row = dataset.get(i)
            ip_chi2_val = row.getRealValue("tag_ip_chi2")
            
            # Calculate log(IP chi2) with safety check
            if ip_chi2_val <= 0:
                log_val = log_range[0]  # Use minimum value for invalid inputs
            else:
                log_val = math.log(ip_chi2_val)
                # Clip to valid range
                log_val = max(log_range[0], min(log_range[1], log_val))
            
            # Set values in the row
            log_ipchi2.setVal(log_val)
            
            # Add row to new dataset
            new_dataset.add(vars_list)
        
        return new_dataset
        
    def _store_ipchi2_fit_parameters(self, bin_idx, param, param_err, fitter=None):
        """Store IP chi2 fit parameters in result arrays for simplified Bukin function fits"""
        # Basic fit results - handling simplified Bukin function parameters
        self.FitIPRes_SYield[bin_idx, 0] = param[0]  # Total signal yield
        self.FitIPRes_PromptFrac[bin_idx, 0] = param[1]  # Prompt fraction
        
        # Prompt Bukin parameters
        self.FitIPRes_XpPrompt = getattr(self, 'FitIPRes_XpPrompt', 
                                        numpy.zeros((self.nzTBins, self.num_fit_items)))
        self.FitIPRes_SigmaPrompt = getattr(self, 'FitIPRes_SigmaPrompt', 
                                           numpy.zeros((self.nzTBins, self.num_fit_items)))
        self.FitIPRes_XiPrompt = getattr(self, 'FitIPRes_XiPrompt', 
                                        numpy.zeros((self.nzTBins, self.num_fit_items)))
        self.FitIPRes_Rho1Prompt = getattr(self, 'FitIPRes_Rho1Prompt', 
                                          numpy.zeros((self.nzTBins, self.num_fit_items)))
        self.FitIPRes_Rho2Prompt = getattr(self, 'FitIPRes_Rho2Prompt', 
                                          numpy.zeros((self.nzTBins, self.num_fit_items)))
        
        self.FitIPRes_XpPrompt[bin_idx, 0] = param[2]  # Peak position
        self.FitIPRes_SigmaPrompt[bin_idx, 0] = param[3]  # Width
        self.FitIPRes_XiPrompt[bin_idx, 0] = param[4]  # Asymmetry
        self.FitIPRes_Rho1Prompt[bin_idx, 0] = param[5]  # Left tail
        self.FitIPRes_Rho2Prompt[bin_idx, 0] = param[6]  # Right tail
        
        # Non-prompt Bukin parameters
        self.FitIPRes_XpNonprompt = getattr(self, 'FitIPRes_XpNonprompt', 
                                          numpy.zeros((self.nzTBins, self.num_fit_items)))
        self.FitIPRes_SigmaNonprompt = getattr(self, 'FitIPRes_SigmaNonprompt', 
                                             numpy.zeros((self.nzTBins, self.num_fit_items)))
        self.FitIPRes_XiNonprompt = getattr(self, 'FitIPRes_XiNonprompt', 
                                          numpy.zeros((self.nzTBins, self.num_fit_items)))
        self.FitIPRes_Rho1Nonprompt = getattr(self, 'FitIPRes_Rho1Nonprompt', 
                                            numpy.zeros((self.nzTBins, self.num_fit_items)))
        self.FitIPRes_Rho2Nonprompt = getattr(self, 'FitIPRes_Rho2Nonprompt', 
                                            numpy.zeros((self.nzTBins, self.num_fit_items)))
        
        self.FitIPRes_XpNonprompt[bin_idx, 0] = param[7]  # Peak position
        self.FitIPRes_SigmaNonprompt[bin_idx, 0] = param[8]  # Width
        self.FitIPRes_XiNonprompt[bin_idx, 0] = param[9]  # Asymmetry
        self.FitIPRes_Rho1Nonprompt[bin_idx, 0] = param[10]  # Left tail
        self.FitIPRes_Rho2Nonprompt[bin_idx, 0] = param[11]  # Right tail
        
        # Derived quantities (yields)
        self.FitIPRes_PromptYield[bin_idx, 0] = param[0] * param[1]  # Calculate prompt yield
        self.FitIPRes_NonPromptYield[bin_idx, 0] = param[0] * (1 - param[1])  # Calculate non-prompt yield
        
        # Store errors for key parameters
        self.FitIPRes_SYield[bin_idx, 1] = param_err[0]
        self.FitIPRes_PromptFrac[bin_idx, 1] = param_err[1]
        
        # Prompt Bukin parameter errors
        if len(param_err) > 2:
            self.FitIPRes_XpPrompt[bin_idx, 1] = param_err[2]
        if len(param_err) > 3:
            self.FitIPRes_SigmaPrompt[bin_idx, 1] = param_err[3]
        if len(param_err) > 4:
            self.FitIPRes_XiPrompt[bin_idx, 1] = param_err[4]
        
        # Non-prompt Bukin parameter errors
        if len(param_err) > 7:
            self.FitIPRes_XpNonprompt[bin_idx, 1] = param_err[7]
        if len(param_err) > 8:
            self.FitIPRes_SigmaNonprompt[bin_idx, 1] = param_err[8]
        if len(param_err) > 9:
            self.FitIPRes_XiNonprompt[bin_idx, 1] = param_err[9]
        
        # Store parameter ranges from fitter dictionary (if available)
        if fitter and hasattr(fitter, 'ipchi2Dict'):
            ipchi2_params = fitter.ipchi2Dict[f"Sig{self.dictKey}"]
            
            # Set parameter ranges from dictionary
            params = [
                # Main parameters
                (self.FitIPRes_PromptFrac, ipchi2_params.get("prompt_frac", (0.95, 0.8, 0.999))),
                
                # Prompt Bukin
                (self.FitIPRes_XpPrompt, ipchi2_params.get("xp_prompt", (0.0, -0.5, 1.0))),
                (self.FitIPRes_SigmaPrompt, ipchi2_params.get("sigma_prompt", (0.7, 0.4, 2.0))),
                (self.FitIPRes_XiPrompt, ipchi2_params.get("xi_prompt", (0.0, -0.5, 0.5))),
                (self.FitIPRes_Rho1Prompt, ipchi2_params.get("rho1_prompt", (0.0, -1.0, 1.0))),
                (self.FitIPRes_Rho2Prompt, ipchi2_params.get("rho2_prompt", (0.0, -1.0, 1.0))),
                
                # Non-prompt Bukin
                (self.FitIPRes_XpNonprompt, ipchi2_params.get("xp_nonprompt", (2.5, 1.0, 4.0))),
                (self.FitIPRes_SigmaNonprompt, ipchi2_params.get("sigma_nonprompt", (2.0, 0.4, 4.0))),
                (self.FitIPRes_XiNonprompt, ipchi2_params.get("xi_nonprompt", (0.3, -0.5, 0.5))),
                (self.FitIPRes_Rho1Nonprompt, ipchi2_params.get("rho1_nonprompt", (0.2, -1.0, 1.0))),
                (self.FitIPRes_Rho2Nonprompt, ipchi2_params.get("rho2_nonprompt", (0.2, -1.0, 1.0)))
            ]
            
            # Store parameter ranges (start, min, max)
            for arr, val_list in params:
                arr[bin_idx, 2] = val_list[0]  # Start value
                arr[bin_idx, 3] = val_list[1]  # Min value
                arr[bin_idx, 4] = val_list[2]  # Max value
                
        # Print summary of key Bukin parameters
        print(f"  IP Chi2 Bukin fit parameters (bin {bin_idx}):")
        print(f"    Prompt fraction: {self.FitIPRes_PromptFrac[bin_idx, 0]:.3f}")
        print(f"    Prompt Bukin: peak={self.FitIPRes_XpPrompt[bin_idx, 0]:.3f}, "
              f"width={self.FitIPRes_SigmaPrompt[bin_idx, 0]:.3f}, "
              f"asymm={self.FitIPRes_XiPrompt[bin_idx, 0]:.3f}")
        print(f"    Non-prompt Bukin: peak={self.FitIPRes_XpNonprompt[bin_idx, 0]:.3f}, "
              f"width={self.FitIPRes_SigmaNonprompt[bin_idx, 0]:.3f}, "
              f"asymm={self.FitIPRes_XiNonprompt[bin_idx, 0]:.3f}")


def create_graphs(key, n_bins, x_pos, y_pos_arr, x_width):
    """Create TGraphErrors objects from arrays of values"""
    graphs = [None] * 4
    endings = ["F", "Start", "LL", "HL"]
    
    for i in range(5):
        if i == 1:  # Skip the error array
            continue
            
        temp_array = numpy.array(y_pos_arr[:, i])
        temp_array_err = numpy.array(y_pos_arr[:, 1])
        
        if i == 0:
            # Main graph with values and errors
            graphs[0] = ROOT.TGraphErrors(
                n_bins,
                numpy.array(x_pos),
                temp_array,
                numpy.array(x_width),
                temp_array_err
            )
            if len(temp_array) > 0:  # Avoid errors with empty arrays
                graphs[0].SetMinimum(min(temp_array))
                graphs[0].SetMaximum(max(temp_array))
            graphs[0].SetName(f"{key}{endings[0]}")
        else:
            # Graphs for range limits
            idx = i-1
            if idx < len(graphs) and idx >= 0:  # Check if index is valid
                graphs[idx] = ROOT.TGraphErrors(
                    n_bins,
                    numpy.array(x_pos),
                    temp_array,
                    numpy.array(x_width)
                )
                if len(temp_array) > 0:  # Avoid errors with empty arrays
                    graphs[idx].SetMinimum(min(temp_array))
                    graphs[idx].SetMaximum(max(temp_array))
                graphs[idx].SetName(f"{key}{endings[idx]}")
    
    return graphs


def create_graphs_asymm_err(key, n_bins, x_pos, y_pos_arr, x_width):
    """Create TGraphAsymmErrors with asymmetric error bars"""
    # Get arrays for start values, lower and upper limits
    start_array = numpy.array(y_pos_arr[:, 2])
    ll_array = numpy.array(y_pos_arr[:, 3])
    hl_array = numpy.array(y_pos_arr[:, 4])
    
    # Calculate lower and upper errors
    err_low = [start_array[i] - ll_array[i] for i in range(len(start_array))]
    err_high = [hl_array[i] - start_array[i] for i in range(len(start_array))]
    
    # Create graph
    graph = ROOT.TGraphAsymmErrors(
        n_bins,
        numpy.array(x_pos),
        start_array,
        numpy.array(x_width),
        numpy.array(x_width),
        numpy.array(err_low),
        numpy.array(err_high)
    )
    
    graph.SetMinimum(min(ll_array))
    graph.SetMaximum(max(hl_array))
    graph.SetName(key)
    
    return graph


def mass_fitter_script(resonance, is_mc):
    """Main function to run the mass fitter"""
    # Convert string to boolean if needed
    if isinstance(is_mc, str):
        is_mc = (is_mc != "False")
    
    # Configuration
    is_fit_single_bin = False  # Single bin fitting mode
    is_zt_observable = True  # Use zT bins instead of dR bins
    is_binned = False  # Use binned data
    
    if is_fit_single_bin:
        # Single bin configuration
        jet_pt = (5, 60)
        z_bins = [0.2,0.5,0.65,0.75,0.85,0.95,1]  # D0 zT bins
        # z_bins = [0, 0.3, 0.5, 0.7, 0.85, 1.0]  # D0 zT bins
        
        # Create and run fit object
        fitter = FitSpectraObject(
            resonance, jet_pt, is_binned, is_mc, z_bins, 
            is_zt_observable
        )
        fitter.startFitting()
    else:
        # Multi-bin configuration
        z_bins = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]  # D0 zT bins
        # z_bins = [0.2, 0.5, 0.65, 0.75, 0.85, 1.0]  # D0 zT bins
        # z_bins = [0, 0.2, 0.5, 0.65, 0.75, 0.85, 1.0]  # D0 zT bins
        r_bins = [0, 0.015, 0.03, 0.06, 0.1, 0.2, 0.5]  # D0 R bins
        
        # pT binning for jets
        start_pt = [5, 10, 15, 20, 30]
        end_pt = [10, 15, 20, 30, 50]
        
        # Process each pT bin
        for jet_bin in range(len(start_pt)):
            jet_pt = (start_pt[jet_bin], end_pt[jet_bin])
            
            # Choose bin array based on observable type
            bin_array = z_bins if is_zt_observable else r_bins
            
            # Create and run fit object
            fitter = FitSpectraObject(
                resonance, jet_pt, is_binned, is_mc, bin_array,
                is_zt_observable
            )
            fitter.startFitting()


if __name__ == '__main__':
    print("===== D0 Mass Fitter Script Starting =====")
    print(f"Python version: {sys.version}")
    print(f"ROOT version: {ROOT.gROOT.GetVersion()}")
    print(f"Current directory: {os.getcwd()}")
    
    parser = argparse.ArgumentParser(description="D0 Mass Fitter Script")
    parser.add_argument(
        "-r", "--resonance", 
        action="store", type=str,
        default="D0", choices=["D0", "X3872", "Psi2S"],
        help="Resonance to analyze (D0, X3872, or Psi2S)"
    )
    parser.add_argument(
        "-m", "--isMC", 
        action="store", default=True,
        help="Boolean flag for MC vs data analysis"
    )
    
    args = parser.parse_args()
    print(f"Running with resonance: {args.resonance}, isMC: {args.isMC}")
    mass_fitter_script(args.resonance, args.isMC)
    print("===== D0 Mass Fitter Script Complete =====")