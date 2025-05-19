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

#
# How to run: python plot_RawSignalYields.py -r D0 -p 5_10 -z True
# How to run: python plot_RawSignalYields.py -r D0 -z True
# How to run: python plot_RawSignalYields.py -r D0 -p 15_20 -z True
#
    
def is_valid_graph(graph):
    """Check if an object is a valid graph with required methods"""
    return (graph is not None and 
            hasattr(graph, 'GetN') and 
            hasattr(graph, 'GetX') and 
            hasattr(graph, 'GetY'))  
class PlotGraphsObject:
    def __init__(self, resonance, ptRange, iszT):

        print("- - - - - - - - - - - - - - - - -")
        print("Start Plots for Resonance: {} in pT range: {}".format(resonance,ptRange))
        print("- - - - - - - - - - - - - - - - -")
        self.resonance = resonance
        self.ptString  = ptRange
        self.ptRange   = ptRange.replace("_", "-")
  
        if iszT=="True" or iszT=="TRUE":
          self.obsTag ="zT"
          self.xTitle="p_{T}^{D^{0}}/p_{T}^{jet}"
        else:
          self.obsTag ="dR"
          self.xTitle="#DeltaR"

        self.minPlotRange     =0
        #self.minPlotRange = 0.2#??
        self.minPlotRange = 0 # for consitency among plots
        self.minFitRange  = self.minPlotRange
        self.maxFitRange  = 1

        #- - - - - - - - - - - - - -
        #rootFileName = "/Users/Eliane/LHCb/PostProcessing/SignalFit/{}/FitParameters{}{}{}_{}IIgenFun.root".format(resonance,self.binTag,resonance,self.obsTag,ptRange)#One version extracts prompr and non prompt fit to the t-spectrum with a histogram the other one with a function
        
        #-Prepare input
        basepath     = "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF"

        
        #basepath     = "/Users/eliane/LHCb/x3872-code/signalFit/Psi2S_Cuts"   #These results should be equivalent to the previous Fitting results
        #basepath    = "/Users/eliane/LHCb/x3872-code/signalFit/Psi2S_noCuts"
        rootFileName = "{}/FitParametersUnBinned{}{}_{}.root".format(basepath,resonance,self.obsTag,self.ptString)#One version extracts prompr and non prompt fit to the t-spectrum with a histogram the other one with a function
        fInFileHisto  = ROOT.TFile(rootFileName) 

        #-Create output
        self.OutfilePath  ="{}/RawSignalYields_{}/".format(basepath,resonance)
        if not os.path.exists(self.OutfilePath):
          os.makedirs(self.OutfilePath)
          print("make new directory: {}".format(self.OutfilePath))
          
        
        if fInFileHisto:
          print("o Open data file: {}".format(fInFileHisto.GetName()))
        else:
          print("o File: {} does not exist, fix before continuing".format(rootFileName))
          return
      
        #- - - - - - - - - - - - - -
        #-Get the relevant histograms
        self.hMYield          = fInFileHisto.Get("FitMSYieldF")
        self.hMYieldSG        = fInFileHisto.Get("FitMSYieldSGF")
        self.hMYieldFix       = fInFileHisto.Get("FitMSYieldFixF")
        if not self.hMYieldFix:
          self.hMYieldFix       = fInFileHisto.Get("FitMSYieldDCBF")
        else:
          print("->->->Still used the FitMSYieldFixF histogram, please update to  FitMSYieldDCBF, aka run the fitting code again!")
        self.hTbDecayFraction = fInFileHisto.Get("FitTRes_BDecF")
        if not is_valid_graph(self.hTbDecayFraction):
            print("Warning: FitTRes_BDecF is not a valid graph, creating a fallback")
            # Try alternative name patterns that might exist in the file
            alt_names = ["FitIPPromptFrac", "hIPPromptFrac", "FitIPRes_PromptFrac", "BdecayFrac"]
            for name in alt_names:
                test_graph = fInFileHisto.Get(name)
                if is_valid_graph(test_graph):
                    print(f"Found alternative B-decay fraction graph: {name}")
                    self.hTbDecayFraction = test_graph
                    break
            
            # If still not valid, create a default graph with 50% B-decay fraction
            if not is_valid_graph(self.hTbDecayFraction):
                print("Creating default B-decay fraction graph with 50% fraction")
                self.hTbDecayFraction = self.create_constant_graph(0.5, "FitTRes_BDecF_default")

        self.gInclFragFunc    = fInFileHisto.Get("ginclFragFunc")
        self.gDecayFragFunc   = fInFileHisto.Get("decayFragFunc")
        self.gPromptFragFunc  = fInFileHisto.Get("promptFragFunc")
        
        #-Get the histogram for the acceptance correction with error checking
        try:
            self.gAccCorr0 = fInFileHisto.Get("CorrHist_TypeRnd") #Rnd selection correction
            if not is_valid_graph(self.gAccCorr0):
                print("Warning: CorrHist_TypeRnd is not a valid graph, creating empty graph")
                self.gAccCorr0 = self.create_empty_graph()
                
            self.gAccCorr1 = fInFileHisto.Get("CorrHist_Type1") #total sel & reco & trigger efficiency
            if not is_valid_graph(self.gAccCorr1):
                print("Warning: CorrHist_Type1 is not a valid graph, creating empty graph")
                self.gAccCorr1 = self.create_empty_graph()
                
            self.gAccCorr2 = fInFileHisto.Get("CorrHist_Type2") #pion reco
            if not is_valid_graph(self.gAccCorr2):
                print("Warning: CorrHist_Type2 is not a valid graph, creating empty graph")
                self.gAccCorr2 = self.create_empty_graph()
                
            self.gAccCorr3 = fInFileHisto.Get("CorrHist_Type3") #pion selection
            if not is_valid_graph(self.gAccCorr3):
                print("Warning: CorrHist_Type3 is not a valid graph, creating empty graph")
                self.gAccCorr3 = self.create_empty_graph()
                
            self.gAccCorr4 = fInFileHisto.Get("CorrHist_Type4") #muon reco
            if not is_valid_graph(self.gAccCorr4):
                print("Warning: CorrHist_Type4 is not a valid graph, creating empty graph")
                self.gAccCorr4 = self.create_empty_graph()
                
            self.gAccCorr5 = fInFileHisto.Get("CorrHist_Type5") #stipping line correction
            if not is_valid_graph(self.gAccCorr5):
                print("Warning: CorrHist_Type5 is not a valid graph, creating empty graph")
                self.gAccCorr5 = self.create_empty_graph()
                
            self.gAccCorr6 = fInFileHisto.Get("CorrHist_Type6") #trigger efficiency
            if not is_valid_graph(self.gAccCorr6):
                print("Warning: CorrHist_Type6 is not a valid graph, creating empty graph")
                self.gAccCorr6 = self.create_empty_graph()
        except Exception as e:
            print(f"Error loading correction histograms: {e}")
            # Create empty graphs as fallbacks
            self.gAccCorr0 = self.create_empty_graph()
            self.gAccCorr1 = self.create_empty_graph()
            self.gAccCorr2 = self.create_empty_graph()
            self.gAccCorr3 = self.create_empty_graph()
            self.gAccCorr4 = self.create_empty_graph()
            self.gAccCorr5 = self.create_empty_graph()
            self.gAccCorr6 = self.create_empty_graph()
        
        #-Total combined factor
        self.gAccCorr = self.multiplyGraphs(self.gAccCorr0,self.gAccCorr1)
        #self.gAccCorr = self.gAccCorr1
        #-Convert the factor to an efficiency
        self.gAccCorr  = self.invertGraphs(self.gAccCorr)
        self.gAccCorr0 = self.invertGraphs(self.gAccCorr0)
        self.gAccCorr1 = self.invertGraphs(self.gAccCorr1)
        self.gAccCorr2 = self.invertGraphs(self.gAccCorr2)
        self.gAccCorr3 = self.invertGraphs(self.gAccCorr3)
        self.gAccCorr4 = self.invertGraphs(self.gAccCorr4)
        self.gAccCorr5 = self.invertGraphs(self.gAccCorr5)
        self.gAccCorr6 = self.invertGraphs(self.gAccCorr6)
        
        self.gAccCorr5_Ext = None
        if ptRange=="40_60" and resonance=="Psi2S":
            print("extrapolate missing points")
            self.gAccCorr5_Ext = self.extrapolateIfNecessary(self.gAccCorr5)
        
        #-Create
        #-create total acc corr graph
        self.gAccCorrTotal = self.multiplyGraphs(self.gAccCorr0,self.gAccCorr2)
        self.gAccCorrTotal = self.multiplyGraphs(self.gAccCorrTotal,self.gAccCorr3)
        self.gAccCorrTotal = self.multiplyGraphs(self.gAccCorrTotal,self.gAccCorr4)
        if self.gAccCorr5_Ext:
            self.gAccCorrTotal = self.multiplyGraphs(self.gAccCorrTotal,self.gAccCorr5_Ext)
        else:
             self.gAccCorrTotal = self.multiplyGraphs(self.gAccCorrTotal,self.gAccCorr5)
        self.gAccCorrTotal = self.multiplyGraphs(self.gAccCorrTotal,self.gAccCorr6)
                
        #-Separate the prompt and non-promt fractions
        self.createPNPFractions(self.hMYield,self.hMYieldSG, self.hMYieldFix, self.hTbDecayFraction,False)
  
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def plotPt(self,isCorrected=False):
  
        #- - - - - - - - - - -
        #- -Plot extracted values for specific pT range
        #- - - - - - - - - - -
        self.plotYieldResult(isCorrected)
        self.plotCorrFacAcceptance()
        gBdecay = self.plotBdecayFraction()
  
        return gBdecay
    ###########################################################################
    # Extrapolate factor, if necessary
    ###########################################################################
    def extrapolateIfNecessary(self, graph):
        
        extrapolatedGraph = ROOT.TGraphErrors(graph)
        nBins   = graph.GetN()
    
        x            = graph.GetX()
        xErr         = graph.GetEX()
        yieldVal     = graph.GetY()
        yErr         = graph.GetEY()
     
        lastYield=0
        lastError=0
        for bin in range(0,nBins):
            graphYield  = yieldVal[bin]
            graphYieldE = yErr[bin]
            if graphYield==0:
                extrapolatedGraph.SetPoint(bin,x[bin],lastYield)
                extrapolatedGraph.SetPointError(bin,xErr[bin], lastError)
            else:
                #extrapolatedGraph.SetPoint(bin,x[bin],0)
                lastYield = graphYield
                lastError = graphYieldE
        
        return extrapolatedGraph
    
    ###########################################################################
    # Plot the acceptance correction factor per bin
    ###########################################################################
    def plotCorrFacAcceptance(self):
       
        self.setOptions()
        outputFilename = "{}FinFig_AccCorrFactor_{}_{}.png".format(self.OutfilePath,self.obsTag,self.ptString)
        #outputFilename = "{}FinFig_AccCorrFactor_{}_{}.C".format(self.OutfilePath,self.obsTag,self.ptString)
        
        c = ROOT.TCanvas("c","c: hist",500*2,450*2)
        c.cd()
        ROOT.TGaxis.SetMaxDigits(3)
        
        # Set pad and histo arrangement
        myPad2 = ROOT.TPad("myPad", "The pad",0,0,1,1)
        #myPad2.SetLeftMargin(0.21)
        myPad2.SetLeftMargin(0.15)
        myPad2.SetTopMargin(0.06)
        myPad2.SetRightMargin(0.04)
        myPad2.SetBottomMargin(0.15)
        myPad2.SetTicks()
        myPad2.Draw()
        myPad2.cd()

        myBlankHisto2 = ROOT.TH1F("myBlankHisto2","Blank Histogram",20,-0.3,1)
        myBlankHisto2.GetXaxis().SetNdivisions(505)
        myBlankHisto2.SetXTitle(self.xTitle)
        myBlankHisto2.GetXaxis().SetTitleSize(0.05)
        myBlankHisto2.GetXaxis().SetRangeUser(-0.3,1)
        if "dR" in self.obsTag:
          myBlankHisto2.GetXaxis().SetRangeUser(0,0.5)
        myBlankHisto2.GetXaxis().SetNdivisions(405)
        myBlankHisto2.GetYaxis().SetTitleOffset(1.35)
        myBlankHisto2.GetYaxis().SetTitleSize(0.055)
        myBlankHisto2.SetLineColor(0)

        #myBlankHisto2.SetYTitle("#epsilon_{acc}x#epsilon_{trig}")
        myBlankHisto2.SetYTitle("corr Facor")
        #myBlankHisto2.GetYaxis().SetRangeUser(0.5,6)
        myBlankHisto2.GetYaxis().SetRangeUser(0,1)
        #myBlankHisto2.GetYaxis().SetRangeUser(0,7)
        myBlankHisto2.Draw("E")
          
        MarkerScale=1.6

        self.setHisto(self.gAccCorr0,MarkerScale,ROOT.kGreen-8) #random sel
        self.gAccCorr0.Draw("same EP")
        #self.setHisto(self.gAccCorr1,MarkerScale,ROOT.kBlack)
        #self.gAccCorr1.Draw("same EP")
        self.setHisto(self.gAccCorr2,MarkerScale,ROOT.kRed-10)
        self.gAccCorr2.Draw("same EP")
        self.setHisto(self.gAccCorr3,MarkerScale,ROOT.kRed-7)
        self.gAccCorr3.Draw("same EP")
        self.setHisto(self.gAccCorr4,MarkerScale,ROOT.kBlue-9)
        self.gAccCorr4.Draw("same EP")
        if self.gAccCorr5_Ext:
            self.setHisto(self.gAccCorr5_Ext,MarkerScale,ROOT.kBlue-4)
            self.gAccCorr5_Ext.SetMarkerStyle(24)
            self.gAccCorr5_Ext.Draw("same EP")
        self.setHisto(self.gAccCorr5,MarkerScale,ROOT.kBlue-4)
        self.gAccCorr5.Draw("same EP")
        self.setHisto(self.gAccCorr6,MarkerScale,ROOT.kBlue+2)
        self.gAccCorr6.Draw("same EP")
        
        #self.gAccCorrTotal
        self.setHisto(self.gAccCorrTotal,MarkerScale,ROOT.kGreen+2)
        self.gAccCorrTotal.Draw("same EP")
        #self.setHisto(self.gAccCorr,MarkerScale,ROOT.kGreen+2)
        #self.gAccCorr.Draw("same EP")
        #. . . . . . . . . . . . . . . . . . . . . .
        #. . . . . . . . . . . . . . . . . . . . . .
        #-Legend
        #-Jet and tag details
        if "dR" in self.obsTag:
          #myLegend0 = ROOT.TLegend(0.15,0.78,0.4,0.9)#2rows
          myLegend0 = ROOT.TLegend(0.5,0.62,0.7,0.8)
        else:
            myLegend0 = ROOT.TLegend(0.15,0.76,0.4,0.9)
        myLegend0.SetTextFont(42)
        myLegend0.SetBorderSize(0)
        myLegend0.SetFillStyle(0)
        myLegend0.SetFillColor(0)
        myLegend0.SetMargin(0.25)
        myLegend0.SetTextSize(0.04)

        #myLegend0.AddEntry(myBlankHisto2,"Anti-#it{k}_{T} #it{R} = 0.5, #it{#eta}_{jet}= 2.5-4","")
        myLegend0.AddEntry(myBlankHisto2,"#it{p}_{T}^{jet}=%s (GeV/#it{c})"%(self.ptRange),"")
        myLegend0.AddEntry(myBlankHisto2,"#it{p}_{T}^{D^{0}}>2 (GeV/#it{c})","")

        #collaboration = ROOT.TLatex(0.69,0.88,"#bf{LHCb}")
        #collaboration.SetNDC()
        #collaboration.SetTextSize(0.044)
        #system = ROOT.TLatex(0.69,0.83,"pp  #sqrt{#it{s}} = 13 TeV")
        #system.SetNDC()
        #system.SetTextSize(0.044)

        myLegend0.Draw()

        
        #myLegend1 = ROOT.TLegend(0.6,0.5,0.8,0.9)
        myLegend1 = ROOT.TLegend(0.17,0.4,0.45,0.72)
        myLegend1.SetTextFont(42)
        myLegend1.SetBorderSize(0)
        myLegend1.SetFillStyle(0)
        myLegend1.SetFillColor(0)
        myLegend1.SetMargin(0.25)
        myLegend1.SetTextSize(0.04)
      
        #myLegend1.AddEntry(self.gAccCorr,"total corr factor","pl")
        myLegend1.AddEntry(self.gAccCorrTotal,"total corr factor","pl")
        myLegend1.AddEntry(self.gAccCorr0,"Rnd sel. corr","pl")
        myLegend1.AddEntry(self.gAccCorr2,"#pi reco","pl")
        myLegend1.AddEntry(self.gAccCorr3,"#pi sel.","pl")
        myLegend1.AddEntry(self.gAccCorr4,"#mu reco","pl")
        myLegend1.AddEntry(self.gAccCorr5,"stripping line corr","pl")
        myLegend1.AddEntry(self.gAccCorr6,"trigger eff","pl")
        myLegend1.AddEntry(self.gAccCorr6,"?? tag selection","")
        myLegend1.Draw()


        c.SaveAs(outputFilename)
        c.Close()
        
        ##self.fOutData.cd()
        ##self.gAccCorr.SetName("{}_{}".format(self.gAccCorr.GetName(),self.ptString))
        ##self.gAccCorr.Write()
        #return self.hTbDecayFraction
        
    ###########################################################################
    def setHisto(self, histo, MarkerScale, color):
        
        histo.SetMarkerSize(1.3*MarkerScale)
        histo.SetMarkerStyle(20)
        histo.SetMarkerColor(color)
        histo.SetLineStyle(1)
        histo.SetLineWidth(2)
        histo.SetLineColor(color)
       
    ###########################################################################
    # Draw the inclusive yield for all pT ranges
    ###########################################################################
    def plotBdecayFraction(self,isCorrected=False):
       
        self.setOptions()
        corrTag=""
        if isCorrected:
          corrTag="_Corr"
        
        outputFilename = "{}FinFig_BdecayFrac_{}_{}{}.png".format(self.OutfilePath,self.obsTag,self.ptString,corrTag)
        
        c = ROOT.TCanvas("c","c: hist",500*2,450*2)
        c.cd()
        ROOT.TGaxis.SetMaxDigits(3)
        
        # Set pad and histo arrangement
        myPad2 = ROOT.TPad("myPad", "The pad",0,0,1,1)
        #myPad2.SetLeftMargin(0.21)
        myPad2.SetLeftMargin(0.15)
        myPad2.SetTopMargin(0.06)
        myPad2.SetRightMargin(0.04)
        myPad2.SetBottomMargin(0.15)
        myPad2.SetTicks()
        myPad2.Draw()
        myPad2.cd()

        myBlankHisto2 = ROOT.TH1F("myBlankHisto2","Blank Histogram",20,0,1)
        myBlankHisto2.GetXaxis().SetNdivisions(505)
        myBlankHisto2.SetXTitle(self.xTitle)
        myBlankHisto2.GetXaxis().SetTitleSize(0.05)
        myBlankHisto2.GetXaxis().SetRangeUser(self.minPlotRange,1)
        myBlankHisto2.GetXaxis().SetNdivisions(405)
        myBlankHisto2.GetYaxis().SetTitleOffset(1.35)
        myBlankHisto2.GetYaxis().SetTitleSize(0.055)
        myBlankHisto2.SetLineColor(0)

        myBlankHisto2.SetYTitle("Non-Prompt Fraction")
        myBlankHisto2.GetYaxis().SetRangeUser(0,1)
        if "dR" in self.obsTag:
          myBlankHisto2.GetXaxis().SetRangeUser(0,0.5)
        myBlankHisto2.Draw("E")
          
        MarkerScale=1.6
        #Final Value
        self.hTbDecayFraction.SetMarkerSize(1.3*MarkerScale)
        self.hTbDecayFraction.SetMarkerStyle(20)
        self.hTbDecayFraction.SetMarkerColor(ROOT.kGreen+2)
        self.hTbDecayFraction.SetLineStyle(1)
        self.hTbDecayFraction.SetLineWidth(2)
        self.hTbDecayFraction.SetLineColor(ROOT.kGreen+2)
        self.hTbDecayFraction.Draw("same EP")
 
        #. . . . . . . . . . . . . . . . . . . . . .
        #. . . . . . . . . . . . . . . . . . . . . .
        #-Legend
        #-Jet and tag details
        if "dR" in self.obsTag:
          #myLegend0 = ROOT.TLegend(0.15,0.78,0.4,0.9)#2rows
          myLegend0 = ROOT.TLegend(0.5,0.62,0.7,0.8)
        myLegend0 = ROOT.TLegend(0.15,0.72,0.4,0.9)
        myLegend0.SetTextFont(42)
        myLegend0.SetBorderSize(0)
        myLegend0.SetFillStyle(0)
        myLegend0.SetFillColor(0)
        myLegend0.SetMargin(0.25)
        myLegend0.SetTextSize(0.04)

        myLegend0.AddEntry(myBlankHisto2,"Anti-#it{k}_{T} #it{R} = 0.5, #it{#eta}_{jet}= 2.5-4","")
        myLegend0.AddEntry(myBlankHisto2,"#it{p}_{T}^{jet}=%s (GeV/#it{c})"%(self.ptRange),"")
        myLegend0.AddEntry(myBlankHisto2,"#it{p}_{T}^{D^{0}}>2 (GeV/#it{c})","")

        collaboration = ROOT.TLatex(0.67,0.88,"#bf{LHCb} in progress")
        collaboration.SetNDC()
        collaboration.SetTextSize(0.044)
        system = ROOT.TLatex(0.67,0.83,"pp  #sqrt{#it{s}} = 13 TeV")
        system.SetNDC()
        system.SetTextSize(0.044)

        myLegend0.Draw()
        collaboration.Draw()
        system.Draw()
        
        
        c.SaveAs(outputFilename)
        c.Close()
        
        return self.hTbDecayFraction
    ###########################################################################
    # Draw the inclusive yield for all pT ranges
    ###########################################################################
    def plotYieldSummary(self, yieldArray, pTArray, normType=0, Seltag="incl", isCorrected=False,corrArray=None):
  
        ColorArray = [ROOT.kAzure, ROOT.kAzure-4, ROOT.kCyan-6,  ROOT.kGreen-3, ROOT.kTeal-6,ROOT.kGreen+4,1,1]
        if len(ColorArray)<len(yieldArray) or len(pTArray)<len(yieldArray):
          print("Adapt size of color array!")

        #-apply acceptance correction
        if isCorrected:
         corrTag="_Corr"
        else:
          corrTag=""
          
        titleAddon = ""
        self.setOptions()
        if normType==0:
          outputFilename = "{}FinFig_YieldSummary_{}{}{}.png".format(self.OutfilePath,self.obsTag,Seltag,corrTag)
        elif normType==1:
          outputFilename = "{}FinFig_YieldSummary_{}{}Norm{}.png".format(self.OutfilePath,self.obsTag,Seltag,corrTag)
          titleAddon = "/d#sigma"
        elif normType==2:
          outputFilename = "{}FinFig_BdecaySummary_{}{}{}.png".format(self.OutfilePath,self.obsTag,Seltag,corrTag)
        
        c = ROOT.TCanvas("c","c: hist",500*2,450*2)
        c.cd()
        ROOT.TGaxis.SetMaxDigits(3)
        
        # Set pad and histo arrangement
        myPad2 = ROOT.TPad("myPad", "The pad",0,0,1,1)
        #myPad2.SetLeftMargin(0.21)
        myPad2.SetLeftMargin(0.15)
        myPad2.SetTopMargin(0.06)
        myPad2.SetRightMargin(0.04)
        myPad2.SetBottomMargin(0.15)
        myPad2.SetTicks()
        myPad2.Draw()
        if "dR" in self.obsTag and normType<2:
          myPad2.cd().SetLogy()
        else:
          myPad2.cd()

        
        max  = yieldArray[0].GetMaximum()
        max1 = yieldArray[1].GetMaximum()
        max2 = yieldArray[2].GetMaximum()
        if max1>max:
            max=max1
        if max2>max:
            max=max2
        print("This is the maximum: {}".format(max))
        if max<0:
            max=10000
        myBlankHisto2 = ROOT.TH1F("myBlankHisto2","Blank Histogram",20,0,1)
        myBlankHisto2.GetXaxis().SetNdivisions(505)
        myBlankHisto2.SetXTitle(self.xTitle)
        myBlankHisto2.GetXaxis().SetTitleSize(0.05)
        myBlankHisto2.GetXaxis().SetRangeUser(self.minPlotRange,1)
        myBlankHisto2.GetXaxis().SetNdivisions(405)
        myBlankHisto2.GetYaxis().SetTitleOffset(1.35)
        myBlankHisto2.GetYaxis().SetTitleSize(0.055)
        myBlankHisto2.SetLineColor(0)

        #Absolute yield
        if "dR" in self.obsTag:
          myBlankHisto2.SetYTitle("dN/dA")
          myBlankHisto2.GetYaxis().SetRangeUser(10,max*4)
          myBlankHisto2.GetXaxis().SetRangeUser(0,0.5)
          if normType==1:
            myBlankHisto2.GetYaxis().SetRangeUser(1,100000)

        else:
          myBlankHisto2.SetYTitle("dN/dz_{T}%s"%(titleAddon))
          myBlankHisto2.GetYaxis().SetRangeUser(0,max*1.2)
          if normType==1:
            myBlankHisto2.GetYaxis().SetRangeUser(0,0.45)
        #-B-decay frac
        if normType==2:
          myBlankHisto2.SetYTitle("Non-Prompt Fraction")
          myBlankHisto2.GetYaxis().SetRangeUser(0,1)
        myBlankHisto2.Draw("E")
          
        MarkerScale=1.6
        #Final Value
        for i in range(0,len(yieldArray)):

          print("Histo Name: {}".format(yieldArray[i].GetName()))
          yieldArray[i].SetMarkerSize(1.3*MarkerScale)
          yieldArray[i].SetMarkerStyle(20)
          yieldArray[i].SetMarkerColor(ColorArray[i])
          yieldArray[i].SetLineStyle(1)
          yieldArray[i].SetLineWidth(2)
          yieldArray[i].SetLineColor(ColorArray[i])
          if normType==1:
            self.normalize(yieldArray[i])
          yieldArray[i].Draw("same EP")
 

        #. . . . . . . . . . . . . . . . . . . . . .
        #. . . . . . . . . . . . . . . . . . . . . .
        #-Legend
        #-about different contributions
        if "dR" in self.obsTag:
          #myLegend1 = ROOT.TLegend(0.55,0.43,0.6,0.59)
          myLegend1 = ROOT.TLegend(0.4,0.7,0.6,0.9)
        else:
          myLegend1 = ROOT.TLegend(0.2,0.7,0.4,0.9)
        if len(yieldArray)==4:
          myLegend1 = ROOT.TLegend(0.2,0.7,0.4,0.9)
        if len(yieldArray)==6:
          myLegend1 = ROOT.TLegend(0.2,0.6,0.4,0.9)

        myLegend1.SetTextFont(42)
        myLegend1.SetBorderSize(0)
        myLegend1.SetFillStyle(0)
        myLegend1.SetFillColor(0)
        myLegend1.SetMargin(0.25)
        myLegend1.SetTextSize(0.04)

        myLegend1.AddEntry(yieldArray[0]," #it{p}_{T}^{jet}:","")
        for i in range(0,len(yieldArray)):
          myLegend1.AddEntry(yieldArray[i],"  %s (GeV/%s)" % (pTArray[i].replace("_", "-"),"#it{c}"),"LP")

        collaboration = ROOT.TLatex(0.67,0.88,"#bf{LHCb} in progress")
        collaboration.SetNDC()
        collaboration.SetTextSize(0.044)
        system = ROOT.TLatex(0.67,0.83,"pp  #sqrt{#it{s}} = 13 TeV")
        system.SetNDC()
        system.SetTextSize(0.044)
        if normType!=2:
            component = ROOT.TLatex(0.67,0.73,"{} component".format(Seltag))
            component.SetNDC()
            component.SetTextSize(0.042)
            component.Draw()

        myLegend1.Draw()
        collaboration.Draw()
        system.Draw()
        
        c.SaveAs(outputFilename)
        c.Close()
    ###########################################################################
    # Draw the Absolute yield, different versions and NP and P components of it
    ###########################################################################
    def normalize(self, graph):
    
      sum=0
      for bin in range(0,graph.GetN()):
        sum+=graph.GetY()[bin]
      
      if "dR" in self.obsTag:
        sum*=0.00001
    
      if sum==0: return
      #print("sum1: {}".format(sum))
      for bin in range(0,graph.GetN()):
        graph.GetY()[bin] *= 1/sum
        graph.GetEY()[bin] *= 1/sum

    ###########################################################################
    def applyAccCorrection(self, graph, gAccCorrInput=None):
      
      print("Divide GRAPHS")
      if not gAccCorrInput:
        print("Try to get the histogram from the constructor")
        gAccCorrInput=self.gAccCorrTotal
        
      nBinsAcc   = gAccCorrInput.GetN()
      nBinsgraph = graph.GetN()
      if nBinsAcc!=nBinsgraph:
        print("!!!This can not work fix this error!!!")
        return None
      
      x            = graph.GetX()
      xErr         = graph.GetEX()
      yieldVal     = graph.GetY()
      yErr         = graph.GetEY()
      corrFact     = gAccCorrInput.GetY()
      corrFactErr  = gAccCorrInput.GetEY()

      CorrectedYield    = numpy.zeros(nBinsAcc)
      CorrectedYieldErr = numpy.zeros(nBinsAcc)
      
      for bin in range(0,graph.GetN()):
        graphYield  = yieldVal[bin]
        graphYieldE = yErr[bin]
        corrFactVal = corrFact[bin]
        corrFactE   = corrFactErr[bin]
        if corrFactVal>0:
          CorrectedYield[bin]   = graphYield/corrFactVal
          CorrectedYieldErr[bin]  = self.propagateError(graphYield, 1./corrFactVal, graphYieldE, corrFactE, Type=2)
        else:
          CorrectedYield[bin]    = 0
          CorrectedYieldErr[bin] = 0
        #-proper error propagation
        #CorrectedYieldErr[bin]  = self.propagateError(graphYield, corrFactVal, graphYieldE, corrFactE, Type=0) #This propagates the stat. error of the corr factor

      maximum = max(CorrectedYield)
      ratioGraph = ROOT.TGraphErrors(nBinsgraph, numpy.array(x), numpy.array(CorrectedYield), numpy.array(xErr), numpy.array(CorrectedYieldErr))
      ratioGraph.SetTitle("{}AccCorr{}".format(graph.GetTitle(),self.ptString))
      ratioGraph.SetName("{}AccCorr{}".format(graph.GetName(),self.ptString))
      ratioGraph.SetMaximum(maximum)
 
      return ratioGraph
      
    ###########################################################################
    def multiplyGraphs(self, graph1, graph2):
      
      print("Multiply GRAPHS")
       
      # Check if both graphs are valid
      if not is_valid_graph(graph1):
        print(f"Warning: First graph is not a valid graph object")
        return self.create_empty_graph()
        
      if not is_valid_graph(graph2):
        print(f"Warning: Second graph is not a valid graph object")
        return self.create_empty_graph()
    
      nBins1 = graph1.GetN()
      nBins2 = graph2.GetN()
    
      if nBins1 != nBins2:
        print(f"Warning: Graphs have different numbers of bins ({nBins1} vs {nBins2})")
        return self.create_empty_graph()
    
      if nBins1 == 0:
        print("Warning: Graphs have zero points")
        return self.create_empty_graph()
    
      x = graph1.GetX()
      xErr = graph1.GetEX()
      yieldVal = graph1.GetY()
      yErr = graph1.GetEY()
      yieldVal2 = graph2.GetY()
      yErr2 = graph2.GetEY()

      CorrectedYield = numpy.zeros(nBins1)
      CorrectedYieldErr = numpy.zeros(nBins1)
      
      for bin in range(0, nBins1):
        graphYield = yieldVal[bin]
        graphYieldE = yErr[bin]
        graphYield2 = yieldVal2[bin]
        graphYieldE2 = yErr2[bin]
        CorrectedYield[bin] = graphYield * graphYield2
        # Proper error propagation
        CorrectedYieldErr[bin] = self.propagateError(graphYield, graphYield2, graphYieldE, graphYieldE2, Type=1)

      maximum = max(CorrectedYield) if len(CorrectedYield) > 0 else 0
      ratioGraph = ROOT.TGraphErrors(nBins1, numpy.array(x), numpy.array(CorrectedYield), 
                                    numpy.array(xErr), numpy.array(CorrectedYieldErr))
      ratioGraph.SetTitle(f"{graph1.GetTitle()}Multipl{self.ptString}")
      ratioGraph.SetName(f"{graph1.GetName()}Multipl{self.ptString}")
      ratioGraph.SetMaximum(maximum)

      return ratioGraph
    ###########################################################################s
    def create_empty_graph(self):
        """Create an empty graph as a fallback"""
        try:
            # Create an empty graph with proper initialization
            empty_graph = ROOT.TGraphErrors()
            empty_graph.SetName(f"empty_graph_{self.ptString}")
            # Set point at 0,0 to ensure it's not completely empty
            empty_graph.SetPoint(0, 0.0, 0.0)
            empty_graph.SetPointError(0, 0.0, 0.0)
            return empty_graph
        except Exception as e:
            print(f"Error creating empty graph: {e}")
            # Last resort fallback
            return None
    ###########################################################################
    def invertGraphs(self, graph1):
      """Invert a graph with proper error handling"""
      
      print("Invert GRAPHS")
      
      # Check if graph is valid
      if not is_valid_graph(graph1):
          print(f"Warning: Graph is not a valid graph object")
          return self.create_empty_graph()
      
      nBins1 = graph1.GetN()
      if nBins1 == 0:
          print("Warning: Graph has zero points")
          return self.create_empty_graph()
      
      x = graph1.GetX()
      xErr = graph1.GetEX()
      yieldVal = graph1.GetY()
      yErr = graph1.GetEY()

      CorrectedYield = numpy.zeros(nBins1)
      CorrectedYieldErr = numpy.zeros(nBins1)
      
      for bin in range(0, nBins1):
          graphYield = yieldVal[bin]
          graphYieldE = yErr[bin]
          
          if graphYield != 0:
              CorrectedYield[bin] = 1.0 / graphYield
              # Proper error propagation
              CorrectedYieldErr[bin] = self.propagateError(1.0, graphYield, 0, graphYieldE, Type=0)
          else:
              CorrectedYield[bin] = 0
              CorrectedYieldErr[bin] = 0

      maximum = max(CorrectedYield) if len(CorrectedYield) > 0 else 0
      ratioGraph = ROOT.TGraphErrors(nBins1, numpy.array(x), numpy.array(CorrectedYield), 
                                    numpy.array(xErr), numpy.array(CorrectedYieldErr))
      ratioGraph.SetTitle(f"{graph1.GetTitle()}Multipl{self.ptString}")
      ratioGraph.SetName(f"{graph1.GetName()}Multipl{self.ptString}")
      ratioGraph.SetMaximum(maximum)

      return ratioGraph

    ###########################################################################
    def propagateError(self, factorA, factorB, factorAErr, factorBErr, Type=0):
      #Calculate error of f(A,B) = A/B
      #                sig(f(A,B))^2 =
      
      #-In this case the scaling applied also has an error
      #A/B
      if Type==0:
        if factorB>0:
          error = math.pow(factorA/factorB,2)*(pow((factorAErr/factorA),2)+pow((factorBErr/factorB),2)) #No correlation assumed
        else:
          error = 0
      #A*B
      elif Type==1:
        if factorA>0 and factorB>0:
          error = math.pow(factorA*factorB,2)*(pow((factorAErr/factorA),2)+pow((factorBErr/factorB),2))
        else:
          error = 0
      #A*B (B has no error)
      elif Type==2:
        error = math.pow(factorAErr*factorB,2)
      else:
        print("Not implemented!")
        return 0
        
      return math.sqrt(error)
    ###########################################################################
    # Draw the Absolute yield, different versions and NP and P components of it
    ###########################################################################
    def plotYieldResult(self,isCorrected=False):
      
      self.setOptions()
      
      # Make sure attributes exist
      if not hasattr(self, 'graphPRaw') or not hasattr(self, 'graphNPRaw') or not hasattr(self, 'graphInclRaw'):
          print("Warning: Required yield graphs not found - they may have failed to be created")
          # Create empty fallback graphs
          if not hasattr(self, 'graphPRaw'):
              self.graphPRaw = self.create_empty_graph()
              self.graphPRaw.SetName(f"PromptRaw_{self.ptString}")
          if not hasattr(self, 'graphNPRaw'):
              self.graphNPRaw = self.create_empty_graph()
              self.graphNPRaw.SetName(f"NonPromptRaw_{self.ptString}")
          if not hasattr(self, 'graphInclRaw'):
              self.graphInclRaw = self.create_empty_graph()
              self.graphInclRaw.SetName(f"inclRaw_{self.ptString}")
    
      gpromtOrig    = self.graphPRaw
      gNpromtOrig   = self.graphNPRaw
      hmassyieldOrig= self.graphInclRaw
      
      corrTag=""
      if isCorrected:
        corrTag="_Corr"
        gpromt     = self.applyAccCorrection(gpromtOrig)    #This is applying the acceptance correction
        gNpromt    = self.applyAccCorrection(gNpromtOrig)   #This is applying the acceptance correction
        hmassyield = self.applyAccCorrection(hmassyieldOrig)#This is applying the acceptance correction
        self.graphPCorr    = gpromt
        self.graphPCorr.SetName("PromptCorr_{}".format(self.ptString))
        self.graphNPCorr   = gNpromt
        self.graphNPCorr.SetName("NonPromptCorr_{}".format(self.ptString))
        self.graphInclCorr = hmassyield
        self.graphInclCorr.SetName("inclCorr_{}".format(self.ptString))
      else:
        gpromt     = gpromtOrig
        gNpromt    = gNpromtOrig
        hmassyield = hmassyieldOrig
        
      #outputFilename = "{}FinFig_Yield_{}_{}{}.pdf".format(self.OutfilePath,self.obsTag,self.ptString,corrTag)
      outputFilename = "{}FinFig_Yield_{}_{}{}.png".format(self.OutfilePath,self.obsTag,self.ptString,corrTag)
      

      #-Jet and tag details
      if "dR" in self.obsTag:
        #myLegend0 = ROOT.TLegend(0.15,0.78,0.4,0.9)#2rows
        myLegend0 = ROOT.TLegend(0.5,0.62,0.7,0.8)
      else:
        myLegend0 = ROOT.TLegend(0.15,0.72,0.4,0.9)
      myLegend0.SetTextFont(42)
      myLegend0.SetBorderSize(0)
      myLegend0.SetFillStyle(0)
      myLegend0.SetFillColor(0)
      myLegend0.SetMargin(0.25)
      myLegend0.SetTextSize(0.04)

      #-Legend
      #-about different contributions
      if "dR" in self.obsTag:
        myLegend1 = ROOT.TLegend(0.55,0.47,0.6,0.59)
      else:
        myLegend1 = ROOT.TLegend(0.2,0.59,0.4,0.7)
      #myLegend1 = ROOT.TLegend(0.6,0.75,0.8,0.92)
      myLegend1.SetTextFont(42)
      myLegend1.SetBorderSize(0)
      myLegend1.SetFillStyle(0)
      myLegend1.SetFillColor(0)
      myLegend1.SetMargin(0.25)
      myLegend1.SetTextSize(0.04)
      #myLegend1.SetEntrySeparation(0.5)
        
      scale=1
      MarkerScale=1.6
      #labelOffset1 = 0.16
      labelOffset2 = 0.23+0.2-0.3
      collaboration = ROOT.TLatex(0.67,0.88,"#bf{LHCb} in progress")
      collaboration.SetNDC()
      collaboration.SetTextSize(0.044*scale)
      system = ROOT.TLatex(0.67,0.83,"pp  #sqrt{#it{s}} = 13 TeV")
      system.SetNDC()
      system.SetTextSize(0.044*scale)
      
      c = ROOT.TCanvas("c","c: hist",500*2,450*2)
      c.cd()
      ROOT.TGaxis.SetMaxDigits(3)
      
      # Set pad and histo arrangement
      myPad2 = ROOT.TPad("myPad", "The pad",0,0,1,1)
      #myPad2.SetLeftMargin(0.21)
      myPad2.SetLeftMargin(0.15)
      myPad2.SetTopMargin(0.06)
      myPad2.SetRightMargin(0.04)
      myPad2.SetBottomMargin(0.15)
      myPad2.SetTicks()
      myPad2.Draw()
      if "dR" in self.obsTag:
        myPad2.cd().SetLogy()
      else:
        myPad2.cd()

      max = hmassyield.GetMaximum()
      #print("This is the maximum: {}".format(max))
      myBlankHisto2 = ROOT.TH1F("myBlankHisto2","Blank Histogram",20,0,1)
      myBlankHisto2.GetXaxis().SetNdivisions(505)
      myBlankHisto2.SetXTitle(self.xTitle)
      myBlankHisto2.GetXaxis().SetTitleSize(0.05)
      myBlankHisto2.GetXaxis().SetRangeUser(self.minPlotRange,1)
      myBlankHisto2.GetXaxis().SetNdivisions(405)
      myBlankHisto2.GetYaxis().SetTitleOffset(1.35)
      myBlankHisto2.GetYaxis().SetTitleSize(0.055)
      myBlankHisto2.SetLineColor(0)

      #Absolute yield
      if "dR" in self.obsTag:
        myBlankHisto2.SetYTitle("dN/dA")
        #myBlankHisto2.GetYaxis().SetRangeUser(100,max*4)
        myBlankHisto2.GetYaxis().SetRangeUser(10,max*4)
        myBlankHisto2.GetXaxis().SetRangeUser(0,0.5)
      else:
        myBlankHisto2.SetYTitle("dN/dz_{T}")
        myBlankHisto2.GetYaxis().SetRangeUser(0,max*1.3)
        #myBlankHisto2.GetYaxis().SetRangeUser(0,max*0.8)
      myBlankHisto2.Draw("E")
        
      #-Total yield
      hmassyield.SetMarkerSize(1.3*MarkerScale)
      hmassyield.SetMarkerStyle(20)
      hmassyield.SetMarkerColor(ROOT.kGreen+2)
      hmassyield.SetLineStyle(1)
      hmassyield.SetLineWidth(2)
      hmassyield.SetLineColor(ROOT.kGreen+2)
      hmassyield.Draw("same EP")

      
      #-Prompt fraction
      gpromt.SetMarkerSize(1.3*MarkerScale)
      #gpromt.SetMarkerStyle(4)
      gpromt.SetMarkerStyle(4)
      gpromt.SetLineWidth(2)
      gpromt.SetMarkerColor(ROOT.kBlue)
      gpromt.SetLineColor(ROOT.kBlue)
      gpromt.Draw("same EP")

      #-Non prompt fraction
      gNpromt.SetMarkerSize(1.3*MarkerScale)
      gNpromt.SetMarkerStyle(4)
      gNpromt.SetLineWidth(2)
      gNpromt.SetMarkerColor(ROOT.kBlue-9)
      gNpromt.SetLineColor(ROOT.kBlue-9)
      gNpromt.Draw("same EP")
      
      myLegend1.AddEntry(hmassyield," inclusive yield","LP")
      myLegend1.AddEntry(gpromt," prompt yield","LP")
      myLegend1.AddEntry(gNpromt," non-prompt yield","LP")


      myLegend0.AddEntry(myBlankHisto2,"Anti-#it{k}_{T} #it{R} = 0.5, #it{#eta}_{jet}= 2.5-4","")
      myLegend0.AddEntry(myBlankHisto2,"#it{p}_{T}^{jet}=%s (GeV/#it{c})"%(self.ptRange),"")
      myLegend0.AddEntry(myBlankHisto2,"#it{p}_{T}^{D^{0}}>2 (GeV/#it{c})","")


      myLegend1.Draw()
      collaboration.Draw()
      system.Draw()

      myLegend0.Draw()

      c.SaveAs(outputFilename)
      c.Close()

      
    ###########################################################################
    def create_constant_graph(self, value, name, n_points=10):
        """Create a graph with constant y value as fallback"""
        try:
            # Create bins from 0 to 1 for zT
            x_vals = numpy.linspace(0, 1, n_points)
            x_errs = numpy.array([0.05] * n_points)  # Fixed bin width
            y_vals = numpy.array([value] * n_points)
            y_errs = numpy.array([0.1] * n_points)  # Reasonable uncertainty
            
            # Create the graph
            graph = ROOT.TGraphErrors(n_points, x_vals, y_vals, x_errs, y_errs)
            graph.SetName(name)
            return graph
        except Exception as e:
            print(f"Error creating constant graph: {e}")
            return self.create_empty_graph()

    ###########################################################################
    def createPNPFractions(self, hyield, hyieldVar1, hyieldVar2, hBdecFrac, Absolute):
      """Create prompt and non-prompt fraction graphs with proper error handling"""
      
      # Check if input objects are valid graphs
      if not all(is_valid_graph(x) for x in [hyield, hyieldVar1, hyieldVar2, hBdecFrac]):
          print("Error: One or more input graphs are not valid")
          missing = []
          if not is_valid_graph(hyield): missing.append("hyield")
          if not is_valid_graph(hyieldVar1): missing.append("hyieldVar1")
          if not is_valid_graph(hyieldVar2): missing.append("hyieldVar2") 
          if not is_valid_graph(hBdecFrac): missing.append("hBdecFrac")
          print(f"Missing valid graphs: {', '.join(missing)}")
          
          # Create empty graphs as fallback
          self.graphPRaw = self.create_empty_graph()
          self.graphPRaw.SetName("PromptRaw_empty")
          self.graphNPRaw = self.create_empty_graph()
          self.graphNPRaw.SetName("NonPromptRaw_empty")
          self.graphInclRaw = self.create_empty_graph()
          self.graphInclRaw.SetName("inclRaw_empty")
          return
    
      pointsA = hyield.GetN()
      pointsB = hBdecFrac.GetN()
      
      if pointsA != pointsB:
          print(f"Warning: Different number of points between graphs: {pointsA} vs {pointsB}")
          print("Attempting to adapt by using common bins...")
          
          # Use the minimum number of points to avoid out-of-bounds
          common_points = min(pointsA, pointsB)
          if common_points == 0:
              print("Error: No common bins to process")
              # Create empty fallback graphs
              self.graphPRaw = self.create_empty_graph()
              self.graphPRaw.SetName("PromptRaw_empty")
              self.graphNPRaw = self.create_empty_graph() 
              self.graphNPRaw.SetName("NonPromptRaw_empty")
              self.graphInclRaw = self.create_empty_graph()
              self.graphInclRaw.SetName("inclRaw_empty")
              return
      else:
          common_points = pointsA
    
      # Proceed with the calculation using common_points instead of pointsA
      NbOfVariations = 3
      graphNP    = [None]*NbOfVariations
      graphP     = [None]*NbOfVariations
      graphyield = [None]*NbOfVariations

      npYield   = numpy.zeros((common_points, NbOfVariations))
      pYield    = numpy.zeros((common_points, NbOfVariations))
      npYieldE  = numpy.zeros((common_points, NbOfVariations))
      pYieldE   = numpy.zeros((common_points, NbOfVariations))
      
      # Extract arrays with safety bounds checking
      x_arr = []
      xErr_arr = []
      yieldVal_arr = []
      yErr_arr = []
      yieldValVar1_arr = []
      yErrVar1_arr = []
      yieldValVar2_arr = []
      yErrVar2_arr = []
      fracVal_arr = []
      fracValE_arr = []
      
      for i in range(common_points):
          x_arr.append(hyield.GetX()[i])
          xErr_arr.append(hyield.GetEX()[i])
          yieldVal_arr.append(hyield.GetY()[i])
          yErr_arr.append(hyield.GetEY()[i])
          yieldValVar1_arr.append(hyieldVar1.GetY()[i])
          yErrVar1_arr.append(hyieldVar1.GetEY()[i])
          yieldValVar2_arr.append(hyieldVar2.GetY()[i])
          yErrVar2_arr.append(hyieldVar2.GetEY()[i])
          fracVal_arr.append(hBdecFrac.GetY()[i])
          fracValE_arr.append(hBdecFrac.GetEY()[i])
      
      # Now process the data
      npYield = [[0 for _ in range(NbOfVariations)] for _ in range(common_points)]
      pYield = [[0 for _ in range(NbOfVariations)] for _ in range(common_points)]
      npYieldE = [[0 for _ in range(NbOfVariations)] for _ in range(common_points)]
      pYieldE = [[0 for _ in range(NbOfVariations)] for _ in range(common_points)]
      
      for i in range(common_points):
          if Absolute:
              scale = 1.0
          else:
              # Scaling for the Radial area of the dR slice
              if "dR" in self.obsTag:
                  Area = math.pi*(math.pow((x_arr[i]+xErr_arr[i]), 2)-math.pow((x_arr[i]-xErr_arr[i]), 2))
                  scale = Area
              else:
                  # Scale by the bin width
                  scale = xErr_arr[i]*2
        
          # Apply scaling
          yieldVal_arr[i] *= 1./scale
          yErr_arr[i] *= 1./scale
          npYield[i][0] = yieldVal_arr[i]*fracVal_arr[i]
          pYield[i][0] = yieldVal_arr[i]*(1-fracVal_arr[i])
          npYieldE[i][0] = self.propagateError(yieldVal_arr[i], fracVal_arr[i], yErr_arr[i], fracValE_arr[i], 1)
          pYieldE[i][0] = self.propagateError(yieldVal_arr[i], (1-fracVal_arr[i]), yErr_arr[i], fracValE_arr[i], 1)
          
          # Variation 1
          yieldValVar1_arr[i] *= 1./scale
          yErrVar1_arr[i] *= 1./scale
          npYield[i][1] = yieldValVar1_arr[i]*fracVal_arr[i]
          pYield[i][1] = yieldValVar1_arr[i]*(1-fracVal_arr[i])
          npYieldE[i][1] = self.propagateError(yieldValVar1_arr[i], fracVal_arr[i], yErrVar1_arr[i], fracValE_arr[i], 1)
          pYieldE[i][1] = self.propagateError(yieldValVar1_arr[i], (1-fracVal_arr[i]), yErrVar1_arr[i], fracValE_arr[i], 1)
          
          if NbOfVariations > 2:
              # Variation 2
              yieldValVar2_arr[i] *= 1./scale
              yErrVar2_arr[i] *= 1./scale
              npYield[i][2] = yieldValVar2_arr[i]*fracVal_arr[i]
              pYield[i][2] = yieldValVar2_arr[i]*(1-fracVal_arr[i])
              npYieldE[i][2] = self.propagateError(yieldValVar2_arr[i], fracVal_arr[i], yErrVar2_arr[i], fracValE_arr[i], 1)
              pYieldE[i][2] = self.propagateError(yieldValVar2_arr[i], (1-fracVal_arr[i]), yErrVar2_arr[i], fracValE_arr[i], 1)
    
      # Find max values safely
      max1 = max(yieldVal_arr) if len(yieldVal_arr) > 0 else 0
      max2 = max(yieldValVar1_arr) if len(yieldValVar1_arr) > 0 else 0
      maxTot = max(max1, max2)
    
      # Set maximums
      try:
          hyield.SetMaximum(maxTot)
          hyieldVar1.SetMaximum(maxTot)
          if NbOfVariations > 2:
              hyieldVar2.SetMaximum(maxTot)
      except Exception as e:
          print(f"Warning: Error setting maximums: {e}")
    
      graphyield[0] = hyield
      graphyield[1] = hyieldVar1
      if NbOfVariations > 2:
          graphyield[2] = hyieldVar2
    
      # Create the new graphs using a different constructor approach
      try:
          # Create graphs point by point to avoid numpy array conversion issues
          for i in range(NbOfVariations):
              graphNP[i] = ROOT.TGraphErrors()
              graphP[i] = ROOT.TGraphErrors()
              
              for j in range(common_points):
                  graphNP[i].SetPoint(j, x_arr[j], npYield[j][i])
                  graphNP[i].SetPointError(j, xErr_arr[j], npYieldE[j][i])
                  
                  graphP[i].SetPoint(j, x_arr[j], pYield[j][i])
                  graphP[i].SetPointError(j, xErr_arr[j], pYieldE[j][i])
        
          # Store the graphs as class attributes
          self.graphPRaw = graphP[0]
          self.graphPRaw.SetName(f"PromptRaw_{self.ptString}")
          self.graphNPRaw = graphNP[0]
          self.graphNPRaw.SetName(f"NonPromptRaw_{self.ptString}")
          self.graphInclRaw = graphyield[0]
          self.graphInclRaw.SetName(f"inclRaw_{self.ptString}")
      except Exception as e:
          print(f"Error creating output graphs: {e}")
          import traceback
          traceback.print_exc()
          # Create fallback empty graphs
          self.graphPRaw = self.create_empty_graph()
          self.graphPRaw.SetName(f"PromptRaw_{self.ptString}")
          self.graphNPRaw = self.create_empty_graph()
          self.graphNPRaw.SetName(f"NonPromptRaw_{self.ptString}")
          self.graphInclRaw = self.create_empty_graph()
          self.graphInclRaw.SetName(f"inclRaw_{self.ptString}")
        
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
def plotRawSignalYields(resonance, ptRange,isZt):

  print("is binned var")
  
  if resonance=="Psi2S":
    #pTRangeArray = ["5_10","10_15","15_20","20_30","30_40","40_100"]#pre March 22
    #pTRangeArray = ["5_10","10_15","15_20","20_30","30_40","40_60"]
    pTRangeArray = ["2_5","5_10","10_15","15_20","20_30","30_40"]
  else:
    #pTRangeArray = ["5_10","10_15","15_20","20_50"]#pre march 22
    #pTRangeArray = ["5_10","10_15","15_20","20_30"]
    pTRangeArray = ["2_5","5_10","10_15","15_20","20_30"]

  yieldArray   = []
  yieldArrayP  = []
  yieldArrayNP = []
  CyieldArray   = []
  CyieldArrayP  = []
  CyieldArrayNP = []
  CyieldArrayUnNorm = []
  accArray     = []
  bFracArray   = []
  
  gObj = [None]*len(pTRangeArray)
  #-Call the object
  if ptRange=="":
    #
    #plotAccCorr = False  #do not apply acceptance correction factors
    plotAccCorr = True  #Apply acceptance correction factors

    for i in range(0,len(pTRangeArray)):
        gObj[i] = PlotGraphsObject(resonance,pTRangeArray[i] ,isZt)
        gObj[i].plotPt(plotAccCorr)
        yieldArray.append(gObj[i].graphInclRaw)
        yieldArrayP.append(gObj[i].graphPRaw)
        yieldArrayNP.append(gObj[i].graphNPRaw)
        #-Corrected distributions (folded with acc,trigg,eff corrections)
        #CyieldArrayUnNorm.append(gObj[i].graphInclNoNormCorr) #don'T remember what that is for
        if plotAccCorr:
            CyieldArray.append(gObj[i].graphInclCorr)
            CyieldArrayP.append(gObj[i].graphPCorr)
            CyieldArrayNP.append(gObj[i].graphNPCorr)
        #accArray.append(hacc1)
        bFracArray.append(gObj[i].hTbDecayFraction)
        
    
    #- - - - - - - - - - - - - -
    #-Create an output file to save all final histograms in
    if plotAccCorr:
        fOutData  = ROOT.TFile("{}CorrectedFinalHistograms_{}.root".format(gObj[0].OutfilePath,resonance), "RECREATE")
        for i in range(0,len(yieldArray)):
            #CyieldArrayUnNorm[i].Write() #don'T remember what that is for
            CyieldArray[i].Write()
            CyieldArrayP[i].Write()
            CyieldArrayNP[i].Write()


    #-Uncorrected plots
    #graphs1.plotYieldSummary(yieldArray,0)
    #graphs1.plotYieldSummary(yieldArray,1) #Line A1
    #-Corrected plots
    #-When you want to use this you need to comment line A1
    if plotAccCorr:
        gObj[0].plotYieldSummary(CyieldArray, pTRangeArray, 0,"incl",plotAccCorr)
        gObj[0].plotYieldSummary(CyieldArray, pTRangeArray, 1,"incl",plotAccCorr)
        gObj[0].plotYieldSummary(CyieldArrayP, pTRangeArray, 0,"promt",plotAccCorr)
        gObj[0].plotYieldSummary(CyieldArrayP, pTRangeArray, 1,"promt",plotAccCorr)
        gObj[0].plotYieldSummary(CyieldArrayNP, pTRangeArray, 0,"Npromt",plotAccCorr)
        gObj[0].plotYieldSummary(CyieldArrayNP, pTRangeArray, 1,"Npromt",plotAccCorr)
    else:
        gObj[0].plotYieldSummary(yieldArray, pTRangeArray, 0,"incl",plotAccCorr)
        gObj[0].plotYieldSummary(yieldArray, pTRangeArray, 1,"incl",plotAccCorr)
        gObj[0].plotYieldSummary(yieldArrayP, pTRangeArray, 0,"promt",plotAccCorr)
        gObj[0].plotYieldSummary(yieldArrayP, pTRangeArray, 1,"promt",plotAccCorr)
        gObj[0].plotYieldSummary(yieldArrayNP, pTRangeArray, 0,"Npromt",plotAccCorr)
        gObj[0].plotYieldSummary(yieldArrayNP, pTRangeArray, 1,"Npromt",plotAccCorr)
      
    gObj[0].plotYieldSummary(bFracArray, pTRangeArray, 2)

  
  #-Plot just one single pT bin
  else:
    graphs = PlotGraphsObject(resonance, ptRange, isZt)
    graphs.plotPt()


 

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
                      help="is it for the zT observable")

  # Parse the arguments
  args = parser.parse_args()
  
  plotRawSignalYields(args.resonance, args.ptRange, args.iszT)
#
# How to run:
# python plot_RawSignalYields.py -r D0 -p 5_10 -z True
# python plot_RawSignalYields.py -r D0 -z True
# python plot_RawSignalYields.py -r D0 -p 15_20 -z True
#
