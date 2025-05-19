import os;
#import yaml
#import parameters as pars
import sys, ROOT
import array
import math
#from ordered_set import OrderedSet
import numpy
import argparse

ROOT.gROOT.SetBatch(True)

#>> This plots graphs for the pptx presentation summarizing the fits and the fit parameters!
#
# How to run: python plot_Graphs.py -r Psi2S -b False -p 5_10 -z True
# How to run: python plot_Graphs.py -r X3872 -p 20_50 -z True
#
# How to run: python plot_Graphs.py -r X3872 -z True


class PlotGraphsObject:
    def __init__(self, resonance, ptRange, iszT, basedir):

        print("- - - - - - - - - - - - - - - - -")
        print("Start Plots for Resonance: {} in pT range: {}".format(resonance,ptRange))
        print("- - - - - - - - - - - - - - - - -")
        self.resonance = resonance
        
        self.ptRange = ptRange.replace("_", "-")
        self.binTag       ="UnBinned"
        print("Bin Tag: {}".format(self.binTag))
        
        if iszT=="True" or iszT=="TRUE":
          self.obsTag ="zT"
        else:
          self.obsTag ="dR"
       
        self.minPlotRange   = 0
        if "D0" in self.resonance:
          self.minPlotRange = 0.2#??
          self.minFitRange  = self.minPlotRange
          self.maxFitRange  = 1
        # elif "Psi" in self.resonance:
        #   self.minPlotRange = 0 #0.4
        #   self.minFitRange  = self.minPlotRange
        #   self.maxFitRange  = 1

        #self.baseDir     = "/Users/eliane/LHCb/x3872-code/signalFit/{}_Cuts".format(resonance)
        self.baseDir     = basedir
        self.OutfilePath = "{}/{}/".format(self.baseDir,ptRange)

        #-if only the time fit is performed by loading previous mass fit parameters
        rootFileName = "{}/FitParametersUnBinned{}{}_{}.root".format(self.baseDir, resonance,self.obsTag,ptRange)
        #rootFileName = "{}/FitParametersUnBinned{}{}_{}IIgenFun.root".format(self.baseDir, resonance,self.obsTag,ptRange)
        #rootFileName = "{}/FitParametersUnBinned{}{}_{}IIgenHis.root".format(self.baseDir, resonance,self.obsTag,ptRange)

        fInFileHisto = ROOT.TFile(rootFileName)
        if fInFileHisto:
          print("o Open data file: {}".format(fInFileHisto.GetName()))
          print(fInFileHisto)
        else:
          print("o File: {} does not exist, fix before continuing".format(rootFileName))
          return

        keyList = ["F","Start","LL","HL"]
        
        self.hMYield     = [None]*4
        self.hMYieldSG   = [None]*4
        self.hMYieldDCB  = [None]*4
        self.hMYieldFix  = [None]*4
        self.hMMean      = [None]*4
        self.hMsigma1    = [None]*4
        self.hMsigma2    = [None]*4
        self.hMDeltaSig  = [None]*4
        self.hMCBFraction= [None]*4
        self.hMalpha     = [None]*4
        self.hMN         = [None]*4
        self.hMPol0      = [None]*4
        self.hMPol1      = [None]*4
        
        self.hTYield          = [None]*4
        self.hTMeanP          = [None]*4
        self.hTMeanNP         = [None]*4
        self.hTsigma1         = [None]*4
        self.hTsigma2         = [None]*4
        self.hTgFraction      = [None]*4
        self.hTnpDecayConst   = [None]*4
        self.hTbDecayFraction = [None]*4
        
        for i in range(0,4):
          print("round: {}".format(i))
          self.hMYield[i]       = fInFileHisto.Get("FitMSYield{}".format(keyList[i]))    #Final result
          self.hMYieldSG[i]     = fInFileHisto.Get("FitMSYieldSG{}".format(keyList[i]))  #Single gaussian fit
          self.hMYieldDCB[i]    = fInFileHisto.Get("FitMSYieldDCB{}".format(keyList[i])) #DCB fit
          self.hMYieldFix[i]    = fInFileHisto.Get("FitMSYieldFix{}".format(keyList[i])) #DCB fit with fixed alpha and n
          self.hMMean[i]        = fInFileHisto.Get("FitMMean{}".format(keyList[i]))
          self.hMsigma1[i]      = fInFileHisto.Get("FitMSig1{}".format(keyList[i]))
          self.hMsigma2[i]      = fInFileHisto.Get("FitMSig2{}".format(keyList[i]))
          self.hMDeltaSig[i]    = fInFileHisto.Get("FitMDeltaSig{}".format(keyList[i]))
          self.hMCBFraction[i]  = fInFileHisto.Get("FitMCBFrac{}".format(keyList[i]))
          self.hMalpha[i]       = fInFileHisto.Get("FitMalpha{}".format(keyList[i]))
          self.hMN[i]           = fInFileHisto.Get("FitMN{}".format(keyList[i]))
          self.hMPol0[i]        = fInFileHisto.Get("FitMPol0{}".format(keyList[i]))
          self.hMPol1[i]        = fInFileHisto.Get("FitMPol1{}".format(keyList[i]))

          self.hTYield[i]          = fInFileHisto.Get("FitTSYield{}".format(keyList[i]))
          self.hTMeanP[i]          = fInFileHisto.Get("FitTMeanP{}".format(keyList[i]))
          self.hTMeanNP[i]         = fInFileHisto.Get("FitTMeanNP{}".format(keyList[i]))
          self.hTsigma1[i]         = fInFileHisto.Get("FitTSig1{}".format(keyList[i]))
          self.hTsigma2[i]         = fInFileHisto.Get("FitTSig2{}".format(keyList[i]))
          self.hTgFraction[i]      = fInFileHisto.Get("FitTGFrac{}".format(keyList[i]))
          self.hTnpDecayConst[i]   = fInFileHisto.Get("FitTRes_Dec{}".format(keyList[i]))
          self.hTbDecayFraction[i] = fInFileHisto.Get("FitTRes_BDec{}".format(keyList[i]))

          # Now look for IP Chi2 parameters in the file
          self.hIPYield = [None]*4
          self.hIPPromptFrac = [None]*4
          self.hIPXpPrompt = [None]*4
          self.hIPSigmaPrompt = [None]*4
          self.hIPXiPrompt = [None]*4
          self.hIPXpNonprompt = [None]*4
          self.hIPSigmaNonprompt = [None]*4
          self.hIPXiNonprompt = [None]*4
          
          # Try to load the IP Chi2 parameters using different possible naming conventions
          for i in range(0, 4):
              print(f"Loading IP Chi2 parameters (round {i})...")
              
              # First try FitIP naming convention
              self.hIPYield[i] = fInFileHisto.Get(f"FitIPSYield{keyList[i]}")
              self.hIPPromptFrac[i] = fInFileHisto.Get(f"FitIPPromptFrac{keyList[i]}")
              self.hIPXpPrompt[i] = fInFileHisto.Get(f"FitIPXpPrompt{keyList[i]}")
              self.hIPSigmaPrompt[i] = fInFileHisto.Get(f"FitIPSigmaPrompt{keyList[i]}")
              self.hIPXiPrompt[i] = fInFileHisto.Get(f"FitIPXiPrompt{keyList[i]}")
              self.hIPXpNonprompt[i] = fInFileHisto.Get(f"FitIPXpNonprompt{keyList[i]}")
              self.hIPSigmaNonprompt[i] = fInFileHisto.Get(f"FitIPSigmaNonprompt{keyList[i]}")
              self.hIPXiNonprompt[i] = fInFileHisto.Get(f"FitIPXiNonprompt{keyList[i]}")
              
              # If not found, try alternative naming
              if not self.hIPPromptFrac[i]:
                  self.hIPPromptFrac[i] = fInFileHisto.Get(f"FitIPRes_PromptFrac{keyList[i]}")
          
          # Load range graphs for IP Chi2 parameters
          self.hIPYieldR = fInFileHisto.Get("FitIPSYieldRange")
          self.hIPPromptFracR = fInFileHisto.Get("FitIPPromptFracRange")
          self.hIPXpPromptR = fInFileHisto.Get("FitIPXpPromptRange")
          self.hIPSigmaPromptR = fInFileHisto.Get("FitIPSigmaPromptRange")
          self.hIPXiPromptR = fInFileHisto.Get("FitIPXiPromptRange")
          self.hIPXpNonpromptR = fInFileHisto.Get("FitIPXpNonpromptRange")
          self.hIPSigmaNonpromptR = fInFileHisto.Get("FitIPSigmaNonpromptRange")
          self.hIPXiNonpromptR = fInFileHisto.Get("FitIPXiNonpromptRange")
          
          # Try alternative naming if not found
          if not self.hIPPromptFracR:
              self.hIPPromptFracR = fInFileHisto.Get("FitIPRes_PromptFracRange")
        
        self.hMYieldR       = fInFileHisto.Get("FitMSYieldRange")
        self.hMMeanR        = fInFileHisto.Get("FitMMeanRange")
        self.hMsigma1R      = fInFileHisto.Get("FitMSig1Range")
        self.hMsigma2R      = fInFileHisto.Get("FitMSig2Range")
        self.hMDeltaSigR    = fInFileHisto.Get("FitMDeltaSigRange")
        self.hMCBFractionR  = fInFileHisto.Get("FitMCBFracRange")
        self.hMalphaR       = fInFileHisto.Get("FitMalphaRange")
        self.hMNR           = fInFileHisto.Get("FitMNRange")
        self.hMPol0R        = fInFileHisto.Get("FitMPol0Range")
        self.hMPol1R        = fInFileHisto.Get("FitMPol1Range")
        
        self.hTYieldR          = fInFileHisto.Get("FitTSYieldRange")
        self.hTMeanPR          = fInFileHisto.Get("FitTMeanPRange")
        self.hTMeanNPR         = fInFileHisto.Get("FitTMeanNPRange")
        self.hTsigma1R         = fInFileHisto.Get("FitTSig1Range")
        self.hTsigma2R         = fInFileHisto.Get("FitTSig2Range")
        self.hTgFractionR      = fInFileHisto.Get("FitTGFracRange")
        self.hTnpDecayConstR   = fInFileHisto.Get("FitTRes_DecRange")
        self.hTbDecayFractionR = fInFileHisto.Get("FitTRes_BDecRange")

        self.gInclFragFunc  = fInFileHisto.Get("ginclFragFunc")
        self.gDecayFragFunc = fInFileHisto.Get("decayFragFunc")
        self.gPromptFragFunc= fInFileHisto.Get("promptFragFunc")

        
        
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def plotAll(self):
        """Plot all fit parameters and yields, checking for existence of required graphs first."""
        
        # Helper function to check if both graphs needed for plotting exist and are valid graph types
        def can_plot(graph1, graph2, name):
            if graph1 is None or graph2 is None:
                print(f"Warning: Skipping {name} plot - one or both required graphs do not exist")
                return False
            
            # Check if graphs have the necessary methods
            required_methods = ['GetMinimum', 'GetMaximum', 'Draw', 'GetN']
            try:
                valid1 = all(hasattr(graph1, method) for method in required_methods)
                valid2 = all(hasattr(graph2, method) for method in required_methods)
                
                if not valid1 or not valid2:
                    print(f"Warning: Skipping {name} plot - one or both graphs are not valid graph types")
                    return False
                    
                # Check if graphs have any points
                if graph1.GetN() == 0 or graph2.GetN() == 0:
                    print(f"Warning: Skipping {name} plot - one or both graphs have no points")
                    return False
                    
                return True
            except Exception as e:
                print(f"Warning: Error checking graph for {name}: {e}")
                return False
        
        print("Plotting all graphs")
        
        # Mass parameter plots
        if can_plot(self.hMYield[0], self.hMYieldR, "Mass Yield"):
            self.plotHistParamR(self.hMYield[0], self.hMYieldR, 
                              f"{self.OutfilePath}MParam_YieldR{self.binTag}.png", "Z P")
        
        if can_plot(self.hMMean[0], self.hMMeanR, "Mass Mean"):
            self.plotHistParamR(self.hMMean[0], self.hMMeanR, 
                              f"{self.OutfilePath}MParam_MeanR{self.binTag}.png", "Z P")
        
        if can_plot(self.hMsigma1[0], self.hMsigma1R, "Mass Sigma1"):
            self.plotHistParamR(self.hMsigma1[0], self.hMsigma1R, 
                              f"{self.OutfilePath}MParam_Sig1R{self.binTag}.png", "Z P")
        
        if can_plot(self.hMsigma2[0], self.hMsigma2R, "Mass Sigma2"):
            self.plotHistParamR(self.hMsigma2[0], self.hMsigma2R, 
                              f"{self.OutfilePath}MParam_Sig2R{self.binTag}.png", "Z P")
        
        if can_plot(self.hMDeltaSig[0], self.hMDeltaSigR, "Mass DeltaSig"):
            self.plotHistParamR(self.hMDeltaSig[0], self.hMDeltaSigR, 
                              f"{self.OutfilePath}MParam_DeltaSigR{self.binTag}.png", "Z P")
        
        if can_plot(self.hMCBFraction[0], self.hMCBFractionR, "Mass CBFraction"):
            self.plotHistParamR(self.hMCBFraction[0], self.hMCBFractionR, 
                              f"{self.OutfilePath}MParam_CBFracR{self.binTag}.png", "Z P")
        
        if can_plot(self.hMalpha[0], self.hMalphaR, "Mass Alpha"):
            self.plotHistParamR(self.hMalpha[0], self.hMalphaR, 
                              f"{self.OutfilePath}MParam_alphaR{self.binTag}.png", "Z P")
        
        if can_plot(self.hMN[0], self.hMNR, "Mass N"):
            self.plotHistParamR(self.hMN[0], self.hMNR, 
                              f"{self.OutfilePath}MParam_NR{self.binTag}.png", "Z P")
        
        # Only attempt to plot polynomial parameters if they exist
        if hasattr(self, 'hMPol0') and hasattr(self, 'hMPol0R') and can_plot(self.hMPol0[0], self.hMPol0R, "Mass Pol0"):
            self.plotHistParamR(self.hMPol0[0], self.hMPol0R, 
                              f"{self.OutfilePath}MParam_pol0R{self.binTag}.png", "Z P")
        
        if hasattr(self, 'hMPol1') and hasattr(self, 'hMPol1R') and can_plot(self.hMPol1[0], self.hMPol1R, "Mass Pol1"):
            self.plotHistParamR(self.hMPol1[0], self.hMPol1R, 
                              f"{self.OutfilePath}MParam_pol1R{self.binTag}.png", "Z P")

        # Plot IP Chi2 parameters instead of lifetime parameters
        ip_params = [
            # Format: (attribute_name, display_name, output_filename)
            ('hIPYield', 'IP Chi2 Yield', 'IPParam_YieldR'),
            ('hIPPromptFrac', 'Prompt Fraction', 'IPParam_PromptFracR'),
            ('hIPXpPrompt', 'Prompt Peak Position', 'IPParam_XpPromptR'),
            ('hIPSigmaPrompt', 'Prompt Width', 'IPParam_SigmaPromptR'),
            ('hIPXiPrompt', 'Prompt Asymmetry', 'IPParam_XiPromptR'),
            ('hIPXpNonprompt', 'Non-Prompt Peak Position', 'IPParam_XpNonpromptR'),
            ('hIPSigmaNonprompt', 'Non-Prompt Width', 'IPParam_SigmaNonpromptR'),
            ('hIPXiNonprompt', 'Non-Prompt Asymmetry', 'IPParam_XiNonpromptR')
        ]
        
        # Look for IP Chi2 parameters using various naming conventions
        for base_attr, display_name, output_file in ip_params:
            # Try different attribute name variations
            attr_variations = [
                base_attr,                      # Base name
                f"FitIP{base_attr[3:]}",        # With 'FitIP' prefix
                f"hIP{base_attr[3:]}",          # With 'hIP' prefix
                base_attr.replace('hIP', 'h')   # Without 'IP' in name
            ]
            
            # Also look for range versions
            range_variations = [f"{attr}R" for attr in attr_variations]
            
            # Try to find matching attributes
            for attr in attr_variations:
                if hasattr(self, attr) and hasattr(self, f"{attr}R"):
                    try:
                        graph1 = getattr(self, attr)[0]
                        graph2 = getattr(self, f"{attr}R")
                        if can_plot(graph1, graph2, display_name):
                            self.plotHistParamR(graph1, graph2, 
                                              f"{self.OutfilePath}{output_file}{self.binTag}.png", "Z P")
                            break
                    except Exception as e:
                        print(f"Error accessing IP Chi2 parameter {attr}: {e}")
        
        # Use updated method to create yield results with proper error checking
        try:
            # Check for required histograms needed for yield calculations
            required_yield_hists = [self.hMYield[0], self.hMYieldSG[0], self.hMYieldDCB[0]]
            
            # Add IP Chi2 prompt fraction parameter - try different attribute names
            ip_prompt_frac = None
            for attr_name in ['hIPPromptFraction', 'hIPPromptFrac', 'FitIPRes_PromptFrac']:
                if hasattr(self, attr_name) and getattr(self, attr_name)[0] is not None:
                    ip_prompt_frac = getattr(self, attr_name)[0]
                    break
            
            if all(hist is not None for hist in required_yield_hists) and ip_prompt_frac is not None:
                # Create yield results using IP prompt fraction instead of B decay fraction
                print("Creating yield results using IP Chi2 prompt fraction...")
                self.createYieldResultWithIPChi2(True, ip_prompt_frac)
                tot, prompt, nonPrompt = self.createYieldResultWithIPChi2(False, ip_prompt_frac)
                
                # Save plots in file if the graphs were created successfully
                if all(g is not None and len(g) > 0 and g[0] is not None for g in [tot, prompt, nonPrompt]):
                    print(f"Saving yield graphs to {self.OutfilePath}YieldPerzTGraphs{self.binTag}.root")
                    
                    outputFilename = f"{self.OutfilePath}YieldPerzTGraphs{self.binTag}.root"
                    fOutData = ROOT.TFile(outputFilename, "RECREATE")
                    
                    tot[0].SetName(f"totalYieldPer_{self.obsTag}")
                    tot[0].Write()
                    prompt[0].SetName(f"promptYieldPer_{self.obsTag}")
                    prompt[0].Write()
                    nonPrompt[0].SetName(f"nonpromptYieldPer_{self.obsTag}")
                    nonPrompt[0].Write()
                  
                    fOutData.Close()
                else:
                    print("Warning: Skipping yield graph output - one or more yield graphs were not created successfully")
            else:
                missing = []
                if any(hist is None for hist in required_yield_hists):
                    missing.append("mass histograms")
                if ip_prompt_frac is None:
                    missing.append("IP Chi2 prompt fraction")
                print(f"Warning: Skipping yield result creation - missing {', '.join(missing)}")
        except Exception as e:
            print(f"Error in yield result creation: {e}")
            import traceback
            traceback.print_exc()

    ###########################################################################
    # Creat the prompt and non prompt components     ##########################
    ###########################################################################
    def createYieldResult(self, Absolute):
      
      self.setOptions()
      
      hmassyield     = self.hMYield[0]   #Final result
      hmassyieldSG   = self.hMYieldSG[0]
      #self.hMYieldDCB[0]
      #hTmassyieldFix = self.hMYieldFix[0]
      hTmassyieldFix = self.hMYieldDCB[0]
      hBdecFrac      = self.hTbDecayFraction[0]
      gpromt, gNpromt, gTotYield = self.createPNPFractions(hmassyield,hmassyieldSG,hTmassyieldFix,hBdecFrac, Absolute)
      
      self.plotYieldResult(gpromt, gNpromt, gTotYield, Absolute,0) #Nominal result
      self.plotYieldResult(gpromt, gNpromt, gTotYield, Absolute,1) #Version1
      self.plotYieldResult(gpromt, gNpromt, gTotYield, Absolute,-1)#All together
      
      return gTotYield, gpromt, gNpromt

    ###########################################################################
    def createYieldResultWithIPChi2(self, Absolute, ip_prompt_frac):
        """
        Create yield results using IP Chi2 prompt fraction instead of lifetime B-decay fraction.
        
        Args:
            Absolute: Whether to use absolute or normalized yields
            ip_prompt_frac: Graph with prompt fraction from IP Chi2 fit
        """
        self.setOptions()
        
        # Get mass yield histograms
        hmassyield = self.hMYield[0]   # Final result
        hmassyieldSG = self.hMYieldSG[0]
        hmassyieldDCB = self.hMYieldDCB[0]
        
        # Use the IP Chi2 prompt fraction instead of B-decay fraction
        gpromt, gNpromt, gTotYield = self.createPNPFractionsWithIPChi2(
            hmassyield, hmassyieldSG, hmassyieldDCB, ip_prompt_frac, Absolute)
        
        # Plot the yield results
        if gpromt and gNpromt and gTotYield:
            self.plotYieldResult(gpromt, gNpromt, gTotYield, Absolute, 0)  # Nominal result
            self.plotYieldResult(gpromt, gNpromt, gTotYield, Absolute, 1)  # Version1
            self.plotYieldResult(gpromt, gNpromt, gTotYield, Absolute, -1) # All together
        
        return gTotYield, gpromt, gNpromt

    ###########################################################################
    def createPNPFractionsWithIPChi2(self, hyield, hyieldVar1, hyieldVar2, hPromptFrac, Absolute):
        """
        Create prompt and non-prompt fraction graphs using IP Chi2 prompt fraction.
        
        Args:
            hyield: Total yield graph
            hyieldVar1, hyieldVar2: Yield graphs from different fit methods 
            hPromptFrac: Prompt fraction from IP Chi2 fit
            Absolute: Whether to use absolute or normalized yields
        """
        # Check if input objects are valid graphs
        if not all(x is not None for x in [hyield, hyieldVar1, hyieldVar2, hPromptFrac]):
            print("Error: One or more input graphs are None")
            # Return empty graphs as fallback
            empty_graph = ROOT.TGraphErrors(0)
            return [empty_graph], [empty_graph], [empty_graph]
        
        # Check if objects have the GetN method (are they proper graphs?)
        for obj, name in zip([hyield, hyieldVar1, hyieldVar2, hPromptFrac], 
                            ["hyield", "hyieldVar1", "hyieldVar2", "hPromptFrac"]):
            if not hasattr(obj, 'GetN'):
                print(f"Error: {name} is not a valid graph object (no GetN method)")
                # Return empty graphs as fallback
                empty_graph = ROOT.TGraphErrors(0)
                return [empty_graph], [empty_graph], [empty_graph]
        
        try:
            pointsA = hyield.GetN()
            pointsB = hPromptFrac.GetN()
            
            if pointsA != pointsB:
                print(f"Warning: Different number of points - yield has {pointsA}, prompt fraction has {pointsB}")
                # Use the minimum number of points to avoid out-of-bounds errors
                points = min(pointsA, pointsB)
            else:
                points = pointsA
                
            if points == 0:
                print("Error: One or both graphs have zero points")
                empty_graph = ROOT.TGraphErrors(0)
                return [empty_graph], [empty_graph], [empty_graph]
            
            NbOfVariations = 3
            graphNP = [None]*NbOfVariations
            graphP = [None]*NbOfVariations
            graphyield = [None]*NbOfVariations
            
            npYield = numpy.zeros((points, NbOfVariations))
            pYield = numpy.zeros((points, NbOfVariations))
            npYieldE = numpy.zeros((points, NbOfVariations))
            pYieldE = numpy.zeros((points, NbOfVariations))
            
            x = hyield.GetX()
            xErr = hyield.GetEX()
            yieldVal = hyield.GetY()
            yErr = hyield.GetEY()
            yieldValVar1 = hyieldVar1.GetY()
            yErrVar1 = hyieldVar1.GetEY()
            yieldValVar2 = hyieldVar2.GetY()
            yErrVar2 = hyieldVar2.GetEY()
            
            # Get prompt fractions from IP Chi2 fit
            frac_x = hPromptFrac.GetX()
            fracVal = hPromptFrac.GetY()
            fracValE = hPromptFrac.GetEY()
            
            for i in range(points):
                if Absolute:
                    scale = 1.
                else:
                    # Scaling for the Radial area of the dR slice
                    if "dR" in self.obsTag:
                        Area = math.pi*(math.pow((x[i]+xErr[i]), 2)-math.pow((x[i]-xErr[i]), 2))
                        scale = Area
                    else:
                        scale = xErr[i]*2
                
                # For prompt fractions from IP Chi2, we use the value directly
                # Note: B-decay fraction would be 1-prompt_fraction
                prompt_fraction = fracVal[i]
                prompt_fraction_err = fracValE[i]
                
                # Calculate yields with proper scaling
                yieldVal[i] *= 1./scale
                yErr[i] *= 1./scale
                
                # For prompt and non-prompt yields using prompt fraction
                npYield[i, 0] = yieldVal[i] * (1.0 - prompt_fraction)  # Non-prompt is 1-prompt_fraction
                pYield[i, 0] = yieldVal[i] * prompt_fraction
                
                # Error propagation
                npYieldE[i, 0] = math.sqrt(
                    pow(yieldVal[i] * prompt_fraction_err, 2) + 
                    pow((1.0 - prompt_fraction) * yErr[i], 2)
                )
                pYieldE[i, 0] = math.sqrt(
                    pow(yieldVal[i] * prompt_fraction_err, 2) + 
                    pow(prompt_fraction * yErr[i], 2)
                )
                
                # Variation 1
                yieldValVar1[i] *= 1./scale
                yErrVar1[i] *= 1./scale
                npYield[i, 1] = yieldValVar1[i] * (1.0 - prompt_fraction)
                pYield[i, 1] = yieldValVar1[i] * prompt_fraction
                npYieldE[i, 1] = math.sqrt(
                    pow(yieldValVar1[i] * prompt_fraction_err, 2) + 
                    pow((1.0 - prompt_fraction) * yErrVar1[i], 2)
                )
                pYieldE[i, 1] = math.sqrt(
                    pow(yieldValVar1[i] * prompt_fraction_err, 2) + 
                    pow(prompt_fraction * yErrVar1[i], 2)
                )
                
                if NbOfVariations > 2:
                    # Variation 2
                    yieldValVar2[i] *= 1./scale
                    yErrVar2[i] *= 1./scale
                    
                    npYield[i, 2] = yieldValVar2[i] * (1.0 - prompt_fraction)
                    pYield[i, 2] = yieldValVar2[i] * prompt_fraction
                    npYieldE[i, 2] = math.sqrt(
                        pow(yieldValVar2[i] * prompt_fraction_err, 2) + 
                        pow((1.0 - prompt_fraction) * yErrVar2[i], 2)
                    )
                    pYieldE[i, 2] = math.sqrt(
                        pow(yieldValVar2[i] * prompt_fraction_err, 2) + 
                        pow(prompt_fraction * yErrVar2[i], 2)
                    )
            
            # Find maximum for scaling 
            max_vals = [max(yieldVal) if len(yieldVal) > 0 else 0,
                       max(yieldValVar1) if len(yieldValVar1) > 0 else 0]
            
            if len(max_vals) > 0:
                maxTot = max(max_vals)
                if maxTot > 0:
                    hyield.SetMaximum(maxTot)
                    hyieldVar1.SetMaximum(maxTot)
            else:
                print("Warning: Could not determine maximum yield value")
            
            graphyield[0] = hyield
            graphyield[1] = hyieldVar1
            
            if NbOfVariations > 2:
                if hyieldVar2 is not None:
                    hyieldVar2.SetMaximum(maxTot)
                    graphyield[2] = hyieldVar2
            
            # Create prompt and non-prompt graphs
            try:
                for i in range(NbOfVariations):
                    # Check for NaN or infinite values
                    valid_points = [(j, npYield[j, i], pYield[j, i]) 
                                    for j in range(points) 
                                    if not (math.isnan(npYield[j, i]) or 
                                           math.isnan(pYield[j, i]) or
                                           math.isinf(npYield[j, i]) or
                                           math.isinf(pYield[j, i]))]
                    
                    if not valid_points:
                        print(f"Warning: No valid points for variation {i}")
                        graphNP[i] = ROOT.TGraphErrors(0)
                        graphP[i] = ROOT.TGraphErrors(0)
                        continue
                    
                    # Extract arrays of valid points
                    valid_indices = [p[0] for p in valid_points]
                    x_vals = numpy.array([x[j] for j in valid_indices])
                    np_vals = numpy.array([npYield[j, i] for j in valid_indices])
                    p_vals = numpy.array([pYield[j, i] for j in valid_indices])
                    x_errs = numpy.array([xErr[j] for j in valid_indices])
                    np_errs = numpy.array([npYieldE[j, i] for j in valid_indices])
                    p_errs = numpy.array([pYieldE[j, i] for j in valid_indices])
                    
                    # Create graphs with valid points
                    graphNP[i] = ROOT.TGraphErrors(len(valid_indices), x_vals, np_vals, x_errs, np_errs)
                    graphP[i] = ROOT.TGraphErrors(len(valid_indices), x_vals, p_vals, x_errs, p_errs)
            except Exception as e:
                print(f"Error creating prompt/non-prompt graphs: {e}")
                import traceback
                traceback.print_exc()
            
            return graphP, graphNP, graphyield
        
        except Exception as e:
            print(f"Error in createPNPFractionsWithIPChi2: {e}")
            import traceback
            traceback.print_exc()
            empty_graph = ROOT.TGraphErrors(0)
            return [empty_graph], [empty_graph], [empty_graph]

    ###########################################################################
    # Draw the Absolute yield, different versions and NP and P components of it
    ###########################################################################
    def plotYieldResult(self, gpromt, gNpromt, gTotYield, Absolute, Version):
        # Existing code...
        
        # Safe access of data for plotting
        if Version != -1:
            if Version < len(gpromt) and Version < len(gNpromt) and Version < len(gTotYield):
                gpromt_use = gpromt[Version]
                gNpromt_use = gNpromt[Version]
                hmassyield = gTotYield[Version]
                
                if not gpromt_use or not gNpromt_use or not hmassyield:
                    print(f"Warning: Missing yield graphs for version {Version}, skipping plot")
                    return
                
                # File naming
                if Absolute:
                    outputFilename = f"{self.OutfilePath}YieldAbs{self.binTag}_V{Version}.png"
                else:
                    outputFilename = f"{self.OutfilePath}YieldPer_{self.obsTag}{self.binTag}_V{Version}.png"
            else:
                print(f"Warning: Version {Version} out of range, skipping plot")
                return
        else:
            # All versions together
            if len(gTotYield) == 0 or not gTotYield[0]:
                print("Warning: Missing total yield graph, skipping plot")
                return
            
            hmassyield = gTotYield[0]
            
            if Absolute:
                outputFilename = f"{self.OutfilePath}YieldAbs{self.binTag}_All.png"
            else:
                outputFilename = f"{self.OutfilePath}YieldPer_{self.obsTag}{self.binTag}_All.png"
        
        # The rest of the function remains the same...
      
        if Version!=-1:
            gpromt    = gpromt[Version]
            gNpromt   = gNpromt[Version]
            hmassyield= gTotYield[Version]
        
            if Absolute:
                outputFilename = "{}YieldAbs{}_V{}.png".format(self.OutfilePath,self.binTag,Version)
            else:
                outputFilename = "{}YieldPer_{}{}_V{}.png".format(self.OutfilePath,self.obsTag,self.binTag,Version)
        else:
            hmassyield= gTotYield[0]

            if Absolute:
                outputFilename = "{}YieldAbs{}_All.png".format(self.OutfilePath,self.binTag)
            else:
                outputFilename = "{}YieldPer_{}{}_All.png".format(self.OutfilePath,self.obsTag,self.binTag)
        
        if "dR" in self.obsTag:
            myLegend1 = ROOT.TLegend(0.3,0.73,0.4,0.9)
        else:
            myLegend1 = ROOT.TLegend(0.2,0.73,0.4,0.9)
        #myLegend1 = ROOT.TLegend(0.6,0.75,0.8,0.92)
        myLegend1.SetTextFont(42)
        myLegend1.SetBorderSize(0)
        myLegend1.SetFillStyle(0)
        myLegend1.SetFillColor(0)
        myLegend1.SetMargin(0.25)
        myLegend1.SetTextSize(0.04)
        #myLegend1.SetEntrySeparation(0.5)
        
        c = ROOT.TCanvas("c","c: hist",600,450)
        c.cd()
        #if Absolute and "dR" in self.obsTag:
    
        ROOT.TGaxis.SetMaxDigits(3)
        # Set pad and histo arrangement
        myPad2 = ROOT.TPad("myPad", "The pad",0,0,1,1)
        #myPad2.SetLeftMargin(0.21)
        myPad2.SetLeftMargin(0.15)
        myPad2.SetTopMargin(0.06)
        myPad2.SetRightMargin(0.04)
        myPad2.SetBottomMargin(0.15)
        myPad2.Draw()
        if "dR" in self.obsTag:
            myPad2.cd().SetLogy()
        else:
            myPad2.cd()

        max = hmassyield.GetMaximum()
        print("This is the maximum: {}".format(max))
        myBlankHisto2 = ROOT.TH1F("myBlankHisto2","Blank Histogram",20,0,1)
        myBlankHisto2.SetNdivisions(505)
        if "dR" in self.obsTag:
            myBlankHisto2.SetXTitle("#DeltaR")
        else:
            # if "X3872" in self.resonance:
                myBlankHisto2.SetXTitle("p_{T}^{D^{0}}/p_{T}^{jet}")
            # else:
                # myBlankHisto2.SetXTitle("p_{T}^{#psi(2S)}/p_{T}^{jet}")

        myBlankHisto2.GetYaxis().SetTitleOffset(1.8)
        myBlankHisto2.GetXaxis().SetRangeUser(0,1)
        
        #Absolute yield
        if Absolute:
            myBlankHisto2.SetYTitle("Yield in {} range".format(self.obsTag))
            if "zT" in self.obsTag:
                if "D0" in self.resonance and max>15000:
                    myBlankHisto2.GetYaxis().SetRangeUser(0,50000)
            #   if "Psi2S" in self.resonance and max>13000:
            #     myBlankHisto2.GetYaxis().SetRangeUser(0,10000)
            # if "Psi2S" in self.resonance and max>600000:
            #   myBlankHisto2.GetYaxis().SetRangeUser(0,250000)

        #Scaled yield
        else:
            if "dR" in self.obsTag:
                myBlankHisto2.SetYTitle("dN/dA")
            else:
                myBlankHisto2.SetYTitle("dN/dz_{T}")
        
        if "dR" in self.obsTag:
            myBlankHisto2.GetYaxis().SetRangeUser(100,max*4)
            myBlankHisto2.GetXaxis().SetRangeUser(0,0.5)
            if Absolute:
                myBlankHisto2.GetYaxis().SetRangeUser(10,max*4)
        else:
            myBlankHisto2.GetYaxis().SetRangeUser(0,max*1.3)
            myBlankHisto2.GetXaxis().SetRangeUser(self.minPlotRange,1)
        myBlankHisto2.SetLineColor(0)
        myBlankHisto2.Draw("E")
        
        #Final Value
        hmassyield.SetMarkerSize(1.3)
        hmassyield.SetMarkerStyle(20)
        hmassyield.SetMarkerColor(ROOT.kGreen+2)
        hmassyield.SetLineStyle(1)
        hmassyield.SetLineWidth(1)
        hmassyield.SetLineColor(ROOT.kGreen+2)
        hmassyield.Draw("same EP")

        if Version!=-1:
            #Start Value
            gpromt.SetMarkerStyle(4)
            gpromt.SetMarkerColor(ROOT.kBlue)
            gpromt.SetLineColor(ROOT.kBlue)
            gpromt.SetMarkerSize(1.2)
            gpromt.Draw("same EP")

            #Lower Limit
            gNpromt.SetMarkerStyle(4)
            gNpromt.SetMarkerColor(ROOT.kBlue-9)
            gNpromt.SetLineColor(ROOT.kBlue-9)
            gNpromt.SetMarkerSize(1.2)
            gNpromt.Draw("same EP")
        
            myLegend1.AddEntry(hmassyield," total yield","LP")
            myLegend1.AddEntry(gpromt," prompt yield","LP")
            myLegend1.AddEntry(gNpromt," non-prompt yield","LP")
        else:
            colorList =[ROOT.kGreen+2,ROOT.kBlue,ROOT.kRed]
            #LegendList=["Double CrystalBall","Single Gauss","Fixed #alpha"]
            LegendList=["Double CB, #alpha and n from MC","Single Gauss","Double CB, #alpha free"]
            myLegend1.AddEntry(hmassyield,LegendList[0],"LP")

            for i in range(1,len(gTotYield)):
                gTotYield[i].SetLineColor(colorList[i])
                gTotYield[i].SetMarkerColor(colorList[i])
                gTotYield[i].SetMarkerSize(1.3)
                gTotYield[i].SetMarkerStyle(20)
                gTotYield[i].Draw("same EP")
                myLegend1.AddEntry(gTotYield[i],LegendList[i],"LP")
            hmassyield.Draw("same EP")

        myLegend1.Draw()
     
     
        myLegend0 = ROOT.TLegend(0.7,0.85,0.9,0.9)
        #myLegend1 = ROOT.TLegend(0.6,0.75,0.8,0.92)
        myLegend0.SetTextFont(42)
        myLegend0.SetBorderSize(0)
        myLegend0.SetFillStyle(0)
        myLegend0.SetFillColor(0)
        myLegend0.SetMargin(0.25)
        myLegend0.SetTextSize(0.04)
        myLegend0.AddEntry(myBlankHisto2,"p_{T}^{Jet}=%s GeV"%(self.ptRange),"")
        myLegend0.Draw()

        c.SaveAs(outputFilename)
        c.Close()
    
    ###########################################################################
    def createPNPFractions(self, hyield, hyieldVar1, hyieldVar2, hBdecFrac, Absolute):
        """Create prompt and non-prompt fraction graphs with proper error handling"""
        
        # Check if input objects are valid graphs
        if not all(x is not None for x in [hyield, hyieldVar1, hyieldVar2, hBdecFrac]):
            print("Error: One or more input graphs are None")
            # Return empty graphs as fallback
            empty_graph = ROOT.TGraphErrors(0)
            return [empty_graph], [empty_graph], [empty_graph]
        
        # Check if objects have the GetN method (are they proper graphs?)
        for obj, name in zip([hyield, hyieldVar1, hyieldVar2, hBdecFrac], 
                             ["hyield", "hyieldVar1", "hyieldVar2", "hBdecFrac"]):
            if not hasattr(obj, 'GetN'):
                print(f"Error: {name} is not a valid graph object (no GetN method)")
                # Return empty graphs as fallback
                empty_graph = ROOT.TGraphErrors(0)
                return [empty_graph], [empty_graph], [empty_graph]
        
        pointsA = hyield.GetN()
        pointsB = hBdecFrac.GetN()
        
        if pointsA != pointsB:
            print(f"Error: Incompatible graphs - hyield has {pointsA} points but hBdecFrac has {pointsB} points")
            # Return empty graphs as fallback
            empty_graph = ROOT.TGraphErrors(0)
            return [empty_graph], [empty_graph], [empty_graph]
        
        NbOfVariations = 3
        graphNP = [None]*NbOfVariations
        graphP = [None]*NbOfVariations
        graphyield = [None]*NbOfVariations
        
        try:
            npYield = numpy.zeros((pointsA, NbOfVariations))
            pYield = numpy.zeros((pointsA, NbOfVariations))
            npYieldE = numpy.zeros((pointsA, NbOfVariations))
            pYieldE = numpy.zeros((pointsA, NbOfVariations))
            
            x = hyield.GetX()
            xErr = hyield.GetEX()
            yieldVal = hyield.GetY()
            yErr = hyield.GetEY()
            yieldValVar1 = hyieldVar1.GetY()
            yErrVar1 = hyieldVar1.GetEY()
            yieldValVar2 = hyieldVar2.GetY()
            yErrVar2 = hyieldVar2.GetEY()
            
            fracVal = hBdecFrac.GetY()
            fracValE = hBdecFrac.GetEY()
            
            for i in range(0, pointsA):
                if Absolute:
                    scale = 1.
                else:
                    # Scaling for the Radial area of the dR slice
                    if "dR" in self.obsTag:
                        Area = math.pi*(math.pow((x[i]+xErr[i]), 2)-math.pow((x[i]-xErr[i]), 2))
                        scale = Area
                    else:
                        scale = xErr[i]*2
                    
                yieldVal[i] *= 1./scale
                yErr[i] *= 1./scale
                npYield[i, 0] = yieldVal[i]*fracVal[i]
                pYield[i, 0] = yieldVal[i]*(1-fracVal[i])
                npYieldE[i, 0] = yieldVal[i]*fracValE[i]+fracVal[i]*yErr[i]
                pYieldE[i, 0] = -yieldVal[i]*fracValE[i]+(1-fracVal[i])*yErr[i]
                
                # Variation 1
                yieldValVar1[i] *= 1./scale
                yErrVar1[i] *= 1./scale
                npYield[i, 1] = yieldValVar1[i]*fracVal[i]
                pYield[i, 1] = yieldValVar1[i]*(1-fracVal[i])
                npYieldE[i, 1] = yieldValVar1[i]*fracValE[i]+fracVal[i]*yErrVar1[i]
                pYieldE[i, 1] = -yieldValVar1[i]*fracValE[i]+(1-fracVal[i])*yErrVar1[i]
                
                if NbOfVariations > 2:
                    # Variation 2
                    yieldValVar2[i] *= 1./scale
                    yErrVar2[i] *= 1./scale
                    
                    npYield[i, 2] = yieldValVar2[i]*fracVal[i]
                    pYield[i, 2] = yieldValVar2[i]*(1-fracVal[i])
                    npYieldE[i, 2] = yieldValVar2[i]*fracValE[i]+fracVal[i]*yErrVar2[i]
                    pYieldE[i, 2] = -yieldValVar2[i]*fracValE[i]+(1-fracVal[i])*yErrVar2[i]
            
            # Find maximum for scaling
            max1 = max(yieldVal)
            max2 = max(yieldValVar1)
            maxTot = max(max1, max2)
            
            hyield.SetMaximum(maxTot)
            hyieldVar1.SetMaximum(maxTot)
            
            graphyield[0] = hyield
            graphyield[1] = hyieldVar1
            
            if NbOfVariations > 2:
                hyieldVar2.SetMaximum(maxTot)
                graphyield[2] = hyieldVar2
            
            # Create two new graphs
            for i in range(0, NbOfVariations):
                graphNP[i] = ROOT.TGraphErrors(pointsA, numpy.array(x), numpy.array(npYield[:, i]), 
                                              numpy.array(xErr), numpy.array(npYieldE[:, i]))
                graphP[i] = ROOT.TGraphErrors(pointsA, numpy.array(x), numpy.array(pYield[:, i]), 
                                             numpy.array(xErr), numpy.array(pYieldE[:, i]))
            
            return graphP, graphNP, graphyield
            
        except Exception as e:
            print(f"Error in createPNPFractions: {e}")
            # Return empty graphs as fallback
            empty_graph = ROOT.TGraphErrors(0)
            return [empty_graph], [empty_graph], [empty_graph]

    ###########################################################################
    # Draw the fit parameter ranges     #######################################
    ###########################################################################
    def plotHistParam(self,graphArray, outputFilename="T.pdf", drawOptions = "", setLogy = False, setLogz = False):
      
      self.setOptions()
      
      c = ROOT.TCanvas("c","c: hist",600,450)
      c.cd()
      if setLogy:
        c.SetLogy()
      if setLogz:
        c.SetLogz()
        
      # Set pad and histo arrangement
      myPad2 = ROOT.TPad("myPad", "The pad",0,0,1,1)
      #myPad2.SetLeftMargin(0.21)
      myPad2.SetLeftMargin(0.15)
      myPad2.SetTopMargin(0.04)
      myPad2.SetRightMargin(0.04)
      myPad2.SetBottomMargin(0.15)
      myPad2.Draw()
      myPad2.cd()

      min = graphArray[0].GetMinimum()
      if graphArray[2].GetMinimum()<min:
        min = graphArray[2].GetMinimum()

      #max = 10
      max = graphArray[0].GetMaximum()
      if graphArray[3].GetMaximum()>max:
        max = graphArray[3].GetMaximum()

      #diff = max-min
      print("min: {}, max: {}".format(min,max))
      
      myBlankHisto2 = ROOT.TH1F("myBlankHisto2","Blank Histogram",20,0,1)
      myBlankHisto2.SetNdivisions(505)
      myBlankHisto2.SetXTitle("z_{T}")
      myBlankHisto2.GetYaxis().SetTitleOffset(1.8)
      myBlankHisto2.SetYTitle("")
      if min<0:
        myBlankHisto2.GetYaxis().SetRangeUser(min*1.1,max*1.1)
      else:
        myBlankHisto2.GetYaxis().SetRangeUser(min*0.9,max*1.1)
      if "dR" in self.obsTag:
        myBlankHisto2.GetXaxis().SetRangeUser(0,0.5)
      else:
        myBlankHisto2.GetXaxis().SetRangeUser(self.minPlotRange,1)
      myBlankHisto2.SetLineColor(0)
      myBlankHisto2.Draw("E")
        
      #Final Value
      graphArray[0].SetMarkerSize(1.8)
      graphArray[0].SetMarkerStyle(28)
      graphArray[0].SetMarkerColor(ROOT.kGreen+2)
      graphArray[0].SetLineStyle(1)
      graphArray[0].SetLineWidth(2)
      graphArray[0].SetLineColor(1)
      graphArray[0].Draw("same {} ".format(drawOptions))

      #Start Value
      graphArray[1].SetMarkerStyle(4)
      graphArray[1].SetMarkerColor(ROOT.kGreen+2)
      graphArray[1].SetMarkerSize(1.2)
      graphArray[1].Draw("same {}".format(drawOptions))

      #Lower Limit
      graphArray[2].SetMarkerStyle(26)
      graphArray[2].SetMarkerColor(ROOT.kGray+1)
      graphArray[2].SetMarkerSize(1.4)
      graphArray[2].Draw("same {}".format(drawOptions))
      
      #Upper Limit
      graphArray[3].SetMarkerStyle(32)
      graphArray[3].SetMarkerColor(ROOT.kGray+1)
      graphArray[3].SetMarkerSize(1.4)
      graphArray[3].Draw("same {}".format(drawOptions))
  
      #g.Draw("E2 ")
      #gAlice0204_275.Draw("PE X0 same") #Plot without x-Error
      #gAlice0204_275.Draw("Z P same")
      
      text = graphArray[1].GetName()
      if text:
        textFit = ROOT.TLatex()
        textFit.SetTextSize(0.05)
        textFit.SetNDC()
        textFit.DrawLatex(0.3,0.8,text)
      
      myLegend0 = ROOT.TLegend(0.7,0.85,0.9,0.9)
      #myLegend1 = ROOT.TLegend(0.6,0.75,0.8,0.92)
      myLegend0.SetTextFont(42)
      myLegend0.SetBorderSize(0)
      myLegend0.SetFillStyle(0)
      myLegend0.SetFillColor(0)
      myLegend0.SetMargin(0.25)
      myLegend0.SetTextSize(0.04)
      myLegend0.AddEntry(myBlankHisto2,"p_{T}^{Jet}=%s GeV"%(self.ptRange),"")
      myLegend0.Draw()
  
      c.SaveAs(outputFilename)
      c.Close()

    ###########################################################################
    # Draw the fit parameter ranges     #######################################
    ###########################################################################
    def plotHistParamR(self, graphFin, tgraphR, outputFilename="T.pdf", drawOptions="", setLogy=False, setLogz=False):
        """
        Draw the fit parameter ranges with proper type checking
        """
        self.setOptions()
        
        # Check if tgraphR is a valid graph type with required methods
        if not tgraphR or not hasattr(tgraphR, 'GetMinimum') or not hasattr(tgraphR, 'GetMaximum'):
            print(f"ERROR: Invalid graph object for {outputFilename}. Cannot get min/max values.")
            return
        
        c = ROOT.TCanvas("c", "c: hist", 600, 450)
        c.cd()
        if setLogy:
            c.SetLogy()
        if setLogz:
            c.SetLogz()
        
        # Set pad and histo arrangement
        myPad2 = ROOT.TPad("myPad", "The pad", 0, 0, 1, 1)
        myPad2.SetLeftMargin(0.15)
        myPad2.SetTopMargin(0.05)
        myPad2.SetRightMargin(0.04)
        myPad2.SetBottomMargin(0.13)
        myPad2.Draw()
        myPad2.cd()

        # Safely get min and max values
        try:
            min_val = tgraphR.GetMinimum()
            max_val = tgraphR.GetMaximum()
            diff = (max_val - min_val)
            mean = min_val + 0.5 * diff
            
            plotMin = mean - 1.1 * 0.5 * diff
            plotMax = mean + 1.1 * 0.5 * diff
            print(f"min: {min_val}, max: {max_val}")
        except Exception as e:
            print(f"Error getting min/max values from graph: {e}")
            # Use default values
            plotMin = -1
            plotMax = 10
            print("Using default plot range due to error")
        
        myBlankHisto2 = ROOT.TH1F("myBlankHisto2", "Blank Histogram", 20, 0, 1)
        myBlankHisto2.SetNdivisions(505)
        myBlankHisto2.SetXTitle(self.obsTag)
        
        # Check if graphFin has a name, otherwise use a default
        if hasattr(graphFin, 'GetName') and graphFin.GetName():
            myBlankHisto2.SetYTitle(graphFin.GetName())
        else:
            myBlankHisto2.SetYTitle("Parameter Value")
        
        myBlankHisto2.GetYaxis().SetTitleOffset(1.5)
        myBlankHisto2.GetYaxis().SetRangeUser(plotMin, plotMax)
        
        if "dR" in self.obsTag:
            myBlankHisto2.GetXaxis().SetRangeUser(0, 0.5)
        else:
            myBlankHisto2.GetXaxis().SetRangeUser(self.minPlotRange, 1)
        
        myBlankHisto2.SetLineColor(0)
        myBlankHisto2.Draw("E")
        
        # Check if final graph exists
        if not graphFin or not hasattr(graphFin, 'Draw'):
            print(f"ERROR: Invalid final graph object for {outputFilename}")
            c.Close()
            return
        
        # Final Value - with safety checks
        try:
            graphFin.SetMarkerSize(1.8)
            graphFin.SetMarkerStyle(28)
            graphFin.SetMarkerColor(ROOT.kGreen+2)
            graphFin.SetLineStyle(1)
            graphFin.SetLineWidth(2)
            graphFin.SetLineColor(1)
            graphFin.Draw("same EP")
        except Exception as e:
            print(f"Error drawing graphFin: {e}")
        
        # Start Value - with safety checks
        try:
            tgraphR.SetMarkerStyle(4)
            tgraphR.SetMarkerColor(ROOT.kGreen+2)
            tgraphR.SetMarkerSize(1.2)
            tgraphR.Draw("same {}".format("PE X0 same"))  # Plot only the point
            
            tgraphR.SetLineColor(0)
            tgraphR.SetFillColor(ROOT.kBlue-6)
            tgraphR.SetFillColorAlpha(ROOT.kBlue-6, 0.2)
            tgraphR.SetFillStyle(1001)
            tgraphR.Draw("same E2")  # Plot only the point
        except Exception as e:
            print(f"Error drawing tgraphR: {e}")
        
        # Plot the final values after the fit (again to overlay)
        try:
            graphFin.Draw("same PE")
        except Exception as e:
            print(f"Error re-drawing graphFin: {e}")
        
        if "alpha" in outputFilename or "N" in outputFilename:
            try:
                # Fit the final values with a line
                self.SetMinMaxFitRange()
                
                polyFunc0 = ROOT.TF1("lineFit", "[0]", self.minFitRange, self.maxFitRange)
                polyFunc1 = ROOT.TF1("lineFit", "[0]+[1]*x", self.minFitRange, self.maxFitRange)
                graphFin.Fit(polyFunc0, "QNRE")
                polyFunc0.SetLineColor(ROOT.kBlue-8)
                polyFunc0.SetLineStyle(3)
                polyFunc0.Draw("same")
                graphFin.Fit(polyFunc1, "QNRE")
                polyFunc1.SetLineColor(ROOT.kBlue)
                
                myLegend1 = ROOT.TLegend(0.2, 0.8, 0.4, 0.9)
                myLegend1.SetTextFont(42)
                myLegend1.SetBorderSize(0)
                myLegend1.SetFillStyle(0)
                myLegend1.SetFillColor(0)
                myLegend1.SetMargin(0.25)
                myLegend1.SetTextSize(0.04)
                myLegend1.AddEntry(polyFunc0, "%1.3f #pm %1.3f" % (polyFunc0.GetParameter(0), polyFunc0.GetParError(0)), "L")
                myLegend1.Draw()
            except Exception as e:
                print(f"Error in fitting section: {e}")
        
        # Add the pt range legend
        myLegend0 = ROOT.TLegend(0.7, 0.85, 0.9, 0.9)
        myLegend0.SetTextFont(42)
        myLegend0.SetBorderSize(0)
        myLegend0.SetFillStyle(0)
        myLegend0.SetFillColor(0)
        myLegend0.SetMargin(0.25)
        myLegend0.SetTextSize(0.04)
        myLegend0.AddEntry(myBlankHisto2, "p_{T}^{Jet}=%s GeV" % (self.ptRange), "")
        myLegend0.Draw()
        
        # Save the canvas
        try:
            c.SaveAs(outputFilename)
        except Exception as e:
            print(f"Error saving canvas to {outputFilename}: {e}")
        
        c.Close()
    
    ###########################################################################
    # Extract min and max fit range for alpha parameter     ###################
    ###########################################################################
    def SetMinMaxFitRange(self):
      nBins     = self.hMYield[0].GetN()
      x         = self.hMYield[0].GetX()
      xErr      = self.hMYield[0].GetEX()
      yieldVal  = self.hMYield[0].GetY()
      yErr      = self.hMYield[0].GetEY()
      minYield  = 40
            
      meanYield = numpy.mean(yieldVal)
      for bin in range(0,nBins):
        print("Value: {}".format(yieldVal[bin]))
        #Ony accdept bins that have at least 10% of mean yield
        if yieldVal[bin]>meanYield*0.1:
          print("reached limit")
          self.minFitRange=x[bin]-xErr[bin]
          break
 
      for bin in range(nBins-1,0,-1):
        print("Max bin counts: {}".format(yieldVal[bin]))
        #Ony accdept bins that have at least 10% of mean yield
        if yieldVal[bin]>meanYield*0.1:
          print("reached limit")
          self.maxFitRange=x[bin]+xErr[bin]
          break
      print("Found min and max for fitting: {}, {}".format(self.minFitRange,self.maxFitRange))
      
    ###########################################################################
    # Configure style options     #############################################
    ###########################################################################
    def setOptions(self):

      font = 42
      
      ROOT.gStyle.SetFrameBorderMode(0)
      ROOT.gStyle.SetFrameFillColor(0)
      ROOT.gStyle.SetCanvasBorderMode(0)
      ROOT.gStyle.SetPadBorderMode(0)
      ROOT.gStyle.SetPadColor(10)
      ROOT.gStyle.SetCanvasColor(10)
      ROOT.gStyle.SetTitleFillColor(10)
      ROOT.gStyle.SetTitleBorderSize(1)
      ROOT.gStyle.SetStatColor(10)
      ROOT.gStyle.SetStatBorderSize(1)
      ROOT.gStyle.SetLegendBorderSize(1)
      
      ROOT.gStyle.SetDrawBorder(0)
      ROOT.gStyle.SetTextFont(font)
      ROOT.gStyle.SetStatFont(font)
      ROOT.gStyle.SetStatFontSize(0.05)
      ROOT.gStyle.SetStatX(0.97)
      ROOT.gStyle.SetStatY(0.98)
      ROOT.gStyle.SetStatH(0.03)
      ROOT.gStyle.SetStatW(0.3)
      ROOT.gStyle.SetTickLength(0.02,"y")
      ROOT.gStyle.SetEndErrorSize(3)
      ROOT.gStyle.SetLabelSize(0.05,"xyz")
      ROOT.gStyle.SetLabelFont(font,"xyz")
      ROOT.gStyle.SetLabelOffset(0.01,"xyz")
      ROOT.gStyle.SetTitleFont(font,"xyz")
      ROOT.gStyle.SetTitleOffset(1.2,"xyz")
      ROOT.gStyle.SetTitleSize(0.045,"xyz")
      ROOT.gStyle.SetMarkerSize(1)
      ROOT.gStyle.SetPalette(1)
      
      ROOT.gStyle.SetOptTitle(0)
      ROOT.gStyle.SetOptStat(0)
      ROOT.gStyle.SetOptFit(0)
      

#________________________________________________________________________
#________________________________________________________________________
def PlotGraphs(resonance, ptRange,isZt,baseDir):

  print("is binned var")
  #-Call the object
  if ptRange!="":
    graphs = PlotGraphsObject(resonance, ptRange,isZt,baseDir)
    graphs.plotAll()
  else:
    # if resonance=="X3872":
    startpT = [ 5,10,15,20]
    endpT   = [10,15,20,30]
    # else:
    #   startpT = [ 5,10,15,20,30, 40]
    #   endpT   = [10,15,20,30,40,60]
    for jetBin in range(0,len(startpT)):
        graphs = PlotGraphsObject(resonance, "{}_{}".format(startpT[jetBin],endpT[jetBin]),isZt, baseDir)
        graphs.plotAll()
  
#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  # Define arguments def plotPowerpoint(resKey, ptRange):

  parser = argparse.ArgumentParser(description="Create a run-by-run QA presentation")
  parser.add_argument("-r", "--resonance", action="store",
                      type=str, metavar="resonance",
                      default="D0",
                      help="List of runs to iterate over")
  parser.add_argument("-p", "--ptRange", action="store",
                      type=str, metavar="ptRange",
                      default="",
                      help="pt range for the jet")
  parser.add_argument("-z", "--iszT", action="store",
                      type=str, metavar="iszT",
                      #default=True,
                      #metavar="isBinned",
                      help="is it for the zT observable")
  parser.add_argument("-b", "--baseDir", action="store",
                      type=str, metavar="baseDir",
                      #default=True,
                      #metavar="isBinned",
                      help="base dir for the figures")

  # Parse the arguments
  args = parser.parse_args()
  
  PlotGraphs(args.resonance, args.ptRange, args.iszT, args.baseDir)
