#! /usr/bin/env python

# Macro to unfold a jet spectrum.
# To run:
#     python unfoldzTDistributionLHCb.py
#  on lxplus
#     lb-run DaVinci/v46r0 python unfoldzTDistributionLHCb.py -v 0

# General
import os
import sys
import argparse
import itertools
import math
from array import *
import numpy
# ROOT
import ROOT
from ROOT import TH1
from ROOT import TH2
# RooUnfold
#ROOT.gSystem.Load("/afs/cern.ch/user/e/efranzis/public/RooUnfold/libRooUnfold")
# ROOT.gSystem.Load("/Users/eliane/LHCb/RooUnfold/libRooUnfold")
from ROOT import RooUnfoldResponse
from ROOT import RooUnfold
from ROOT import RooUnfoldBayes
from ROOT import RooUnfoldSvd

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)
    
ROOT.gErrorIgnoreLevel = ROOT.kPrint
ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 1; " ) #does not work somehow

#should we normalize the RM?
#Todo
#* done - Change prior swap from RM swap to ACTUAL Prior swap -> include method to provide external prior
#include kinematic efficiency
#prev. used pT=7+3-2 as lower input and reported starting from 20
#check selection cuts in PrepareResponseMatrix2D and PrepareResponseMatrix3D
#is the prior/truth with or without matching eg under/overlow bins
class UnfoldSpectraClass :
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__(self, promptFlag, inFileNameRM, inFileNameData, resonance, ptRangeArray):
        #general variables
        
        self.useTagPt = False  #-Bin in tag pT instead of jet pT

        #LOCALLY
        #outputDir="/Users/Eliane/LHCb/PostProcessing/UnfoldJetFragmentation/UnfoldedDistributions"
        self.inPathRM   = "/Users/eliane/LHCb/OutFiles/X3872_MC/"
        
        
        
        #- - - - - - - - - - - - - - - - -
        #-Standard result
        #inPathData = "/Users/eliane/LHCb/x3872-code/signalFit/{}_ProbNN07/RawSignalYields_{}/".format(resonance,resonance)
        #outpath    = "/Users/eliane/LHCb/x3872-code/unfolding/OutputApril_"

        #- - - - - - - - - - - - - - - - -
        #-Special attempt to plot >20GeV yield only unfold to only two jet pT bins.
        #-You only need to change:  self.pTBinsArrayTruth  = array('i', [5,20,30])  or  array('i', [5,20,60])   -- THATS ALL :-)
        #inPathData = "/Users/eliane/LHCb/x3872-code/signalFit/{}_ProbNN07/RawSignalYields_{}/".format(resonance,resonance)
        #outpath    = "/Users/eliane/LHCb/x3872-code/unfolding/OutputAprilSingular_"
        
        #- - - - - - - - - - - - - - - - -
        #-tag pt result
        self.useTagPt = True
        inPathData = "/Users/eliane/LHCb/x3872-code/signalFit/{}_ProbNN07_TagpT/RawSignalYields_{}/".format(resonance,resonance)
        outpath    = "/Users/eliane/LHCb/x3872-code/unfolding/OutputApril_TagpT_"
        
        
        
        self.dictKey    = resonance
        if self.dictKey=="Psi2S":
            self.figureTag    = "PII_{}".format(promptFlag)
        else:
            self.figureTag    = "X_{}".format(promptFlag)

        outpath   +=resonance
        self.outpathBase = outpath
        if promptFlag=="P":
            self.isPrompt     = True
            outpath    += "/Prompt_0/"
        elif promptFlag=="NP":
            self.isPrompt     = False
            outpath    += "/NonPrompt_0/"
        
        self.applyRMCut   = True   #pose certain selections on the RM
        self.inFileNRM    = "{}{}.root".format(self.inPathRM,inFileNameRM)
        self.inFileNData  = "{}{}.root".format(inPathData,inFileNameData)

        #-Set binning for input and output spectrum
        self.zBinsArrayDet    = []   #This is the binning of the measured input sprectum. Is set in getRawSpectra()
        if self.dictKey=="Psi2S":
            #self.zBinsArrayTruth   = ([0,0.4,0.6,0.84,0.94,1])#Half of Psi2S bins
            self.zBinsArrayTruth   = ([0,0.3,0.4,0.5,0.6,0.75,0.84,0.9,0.94,0.97,1])#Same as det lvl of Psi2S bins
        else:
            self.zBinsArrayTruth   = ([0.2,0.5,0.65,0.75,0.85,0.95,1])#Same as det lvl of X3872 bins
        #
        #self.pTBinsArrayTruth  = array('i', [5,20,30]) #Special plot
        #self.pTBinsArrayTruth  = array('i', [5,20,60]) #Special plot
        self.pTBinsArrayTruth  = array('i',ptRangeArray) #Same as input, for the time being
        self.pTBinsArrayDet    = array('i',ptRangeArray) #This is the array with the pT ranges of the spectra
        #-Lists of unfolded output for different number of iterations
        self.unfoldedArr2D     = []  #Is set in unfold2D()
        self.unfoldedArrPerBin = []  #Is set in unfold1D()
        self.hOriginalPrior    = None #Is set in unfold2D()
        self.hCurrentPrior     = None #Is set in unfold2D()
        self.externalPrior     = None #Is set in provideExtPrior()
        self.errorType         = RooUnfold.kCovToy  # RooUnfold.kCovToy , ROOT.RooUnfold.kCovariance
        
        for i in range(0,len(ptRangeArray)-1):
            self.unfoldedArrPerBin.append([]) # store an empty list for each pT bin (this will later contain the unfolding results for several iterations)
        self.unfoldLabel = ""
        #self.RM          = None         #Is set in PrepareResponseMatrix()

        #- - - - - - - - - - - - - - - - - - - - - -
        #-Prepare Output Dir
        self.outPath     = outpath
        if not self.outPath.endswith("/"):
          self.outPath = self.outPath + "/"
        if not os.path.exists(self.outPath):
          os.makedirs(self.outPath)
          
        #- - - - - - - - - - - - - - - - - - - - - -
        #-GetInput Files
        self.fResponse  = ROOT.TFile(self.inFileNRM)
        self.fData      = ROOT.TFile(self.inFileNData)

        print("o Open input File with Response Matrix: {}".format(self.fResponse.GetName()))
        print("o Open input File with Data: {}".format(self.fData.GetName()))
        self.measuredSpectraArray = None
        self.measuredSpectraArray = self.getRawSpectra()
        self.measuredSpectra2D    = self.getRawSpectra(True)

          
          
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def unfold1D(self, ptBin, regParam, iteration, powerLawOffset=0):
        
        self.unfoldLabel = "1DBayes{}".format(regParam)
        self.unfoldLabel+= "_Round{}".format(iteration)
 
        #-Create specific output directory
        outputDirUnfolding = os.path.join(self.outPath, self.unfoldLabel)
        print(os.path.isdir(outputDirUnfolding))
        if not os.path.isdir(outputDirUnfolding): # Create current label dir if necessary
            print("directory does not exist. Make directory: {}".format(outputDirUnfolding))
            os.makedirs(outputDirUnfolding)
     
        #-Prepare the 2D response matrix and plot QA
        responseMatrix  = self.PrepareResponseMatrix2D()  #Returns a RooUnfoldResponse object
        responseMatrix1 = self.PrepareResponseMatrix2D(1) #Returns a RooUnfoldResponse object
        responseMatrix2 = self.PrepareResponseMatrix2D(2) #Returns a RooUnfoldResponse object

        #-Call RooUnfold
        for i in range(1,regParam+3):
            
            #-Set up the Bayesian unfolding object
            unfoldBayes = ROOT.RooUnfoldBayes(responseMatrix, self.measuredSpectraArray[ptBin], i)
            #unfoldBayes.SetNToys(1000)
            
            #-Perform the unfolding
            hSpectrumUnfoldedPerBin = unfoldBayes.Hreco(self.errorType)
            hSpectrumUnfoldedPerBin.SetName("UnfoldedSprectraPerBin_nIter{}".format(i))
            self.plotHist(hSpectrumUnfoldedPerBin,"UnfoldedSprectraBin{}_nIter{}".format(ptBin,i),"E")
            self.unfoldedArrPerBin[ptBin].append(hSpectrumUnfoldedPerBin)
            # Produces the truth distribution, with errors, PerBin (will scale by bin width below)
            # 1 -- (default) sqrt of cov matrix diagonals (for SVD, uses toy MCs to account for response matrix errors)
            # 3 -- sqrt of cov matrix from toy MC tests
            
            #-Plot Pearson correlation coeffs for each nIter, to get a measure of the correlation between the bins
            covarianceMatrix = unfoldBayes.Ereco(errorType) # Get the covariance matrix
            self.plotCorrelationCoefficients(covarianceMatrix, i, "CorMatr_Bin{}".format(ptBin))
        
            #. . . . . . . . . . . . . . . . . . . .
            #-1) Refolding Test -- unfold measured spectrum with response1, then apply response2 to unfolded result, and compare to the measured spectrum.
            self.RefoldingTest(ptBin, i, self.measuredSpectraArray[ptBin], responseMatrix1,responseMatrix2)
            
            #. . . . . . . . . . . . . . . . . . . .
            #-2) Unfolding test -- unfold the smeared det-level result with response, and compare to truth-level MC.
            # That just tests the RM itself, no measured spectrum involved
            self.UnfoldingTest(responseMatrix, i)
            
        
        #. . . . . . . . . . . . . . . . . . . .
        #-3) Stability of the number of iterations
        self.StabilityTest(ptBin, regParam)
        self.StabilityTest(ptBin, regParam+1)
        self.StabilityTest(ptBin, regParam, True)

        #. . . . . . . . . . . . . . . . . . . .
        #-4) Test of the rubustness of the RM in terms of statistic
        self.StatTestRM(responseMatrix,ptBin)
        
        #. . . . . . . . . . . . . . . . . . . .
        #-5) Plot the result before and after unfolding
        self.plotUnfoldingEffect(ptBin, regParam)
        #self.plotUnfoldingEffect(ptBin, regParam+1)

    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def unfold2D(self, regParam, iteration, tag=""):

        self.unfoldLabel = "2DBayes{}_Round{}".format(regParam,iteration)

        #-update output path
        if self.isPrompt:
            self.outPath =  self.outpathBase+"/Prompt{}_{}/".format(tag,iteration)
        else:
            self.outPath =  self.outpathBase+"/NonPrompt{}_{}/".format(tag,iteration)
        if not os.path.exists(self.outPath):
            os.makedirs(self.outPath)
         
        #-Prepare RM
        hweights=None
        if iteration>0:
            hweights = self.prepareRMWeights(regParam,iteration)
     
        responseMatrix3D          = self.PrepareResponseMatrix3D(0,hweights) #Returns a RooUnfoldResponse object
        if iteration==0:
            self.hOriginalPrior   = responseMatrix3D.Htruth() # set the current prior
        else:
            self.hCurrentPrior    = responseMatrix3D.Htruth() # set the current prior
        responseMatrix3D_Pt1  = self.PrepareResponseMatrix3D(1,hweights)  #Returns a RooUnfoldResponse object
        responseMatrix3D_Pt2  = self.PrepareResponseMatrix3D(2,hweights)  #Returns a RooUnfoldResponse object
        
        if iteration==0 and self.externalPrior:
            #prepare all the RM again with the new  external prior weight
            #need to get the original response first to be able to weight properly in the first round
            hweights = self.prepareRMWeights(regParam,iteration)
            responseMatrix3D      = self.PrepareResponseMatrix3D(0,hweights)  #Returns a RooUnfoldResponse
            responseMatrix3D_Pt1  = self.PrepareResponseMatrix3D(1,hweights)  #Returns a RooUnfoldResponse
            responseMatrix3D_Pt2  = self.PrepareResponseMatrix3D(2,hweights)  #Returns a RooUnfoldResponse

        
        
        #-Call RooUnfold
        self.plotHist(self.measuredSpectra2D,"hMeasuredData2D","colz","#it{p}_{T,jet} (GeV/#it{c})","#it{z}_{T}")
        
        self.unfoldedArr2D.clear()
        print("length of array: {}".format(len(self.unfoldedArr2D)))
        for i in range(1,regParam+3):
        #for i in range(1,regParam+5):
        #for i in range(regParam,regParam+1):
            
            #-Set up the Bayesian unfolding object
            unfoldBayes2D = ROOT.RooUnfoldBayes(responseMatrix3D, self.measuredSpectra2D, i)
            #unfoldBayes2D.SetVerbose(1000)#print out all details of the unfolding

            #-Perform the unfolding
            #-Produces the truth distribution, with errors, PerBin (will scale by bin width below)
            h2DUnfoldedPerBin = unfoldBayes2D.Hreco(self.errorType)
            h2DUnfoldedPerBin.SetName("hUnfoldedData2DPerBin_nIter{}".format(i))
            self.plotHist(h2DUnfoldedPerBin,"hUnfoldedData2D_nIter{}".format(i),"colz","#it{p}_{T,jet} (GeV/#it{c})","#it{z}_{T}")
            self.unfoldedArr2D.append(h2DUnfoldedPerBin)
            
            #. . . . . . . . . . . . . . . . . . . .
            #-0) Plot Pearson correlation coeffs for each nIter, to get a measure of the correlation between the bins
            covarianceMatrix = unfoldBayes2D.Ereco(errorType) # Get the covariance matrix
            self.plotCorrelationCoefficients(covarianceMatrix, i, "CorMatr2D")

            #. . . . . . . . . . . . . . . . . . . .
            #-1) Refolding Test -- unfold measured spectrum with response1, then apply response2 to unfolded result, and compare to the measured spectrum.
            self.RefoldingTest2D(i, responseMatrix3D_Pt1, responseMatrix3D_Pt2)
        
        #. . . . . . . . . . . . . . . . . . . .
        #-2) Unfolding test -- unfold the det-level MC spectrum with the response, and compare to truth-level MC.
        # That just tests the RM itself, no measured spectrum involved
        # Also called closure test
        self.UnfoldingTest2D(responseMatrix3D, regParam)
        
        
        #. . . . . . . . . . . . . . . . . . . .
        #-3) Stability of the number of iterations
        self.StabilityTest2D(regParam)
        self.StabilityTest2D(regParam+1)
        #self.StabilityTest2D(regParam, True)
        
        #. . . . . . . . . . . . . . . . . . . .
        #-4) Evaluate regularization parameter
        #``The optimal regularisation parameter can be selected by finding the largest
        #value up to which the errors remain reasonable´´
        self.TestRegParam2D()
        
        #. . . . . . . . . . . . . . . . . . . .
        #-5) Test of the robustness of the RM in terms of statistic
        if iteration==0:
            self.StatTestRM2D(responseMatrix3D)
        
        #. . . . . . . . . . . . . . . . . . . .
        #-6) Plot the prior as compared to the unfolded result
        self.plotPrior2D(regParam, responseMatrix3D)

        #. . . . . . . . . . . . . . . . . . . .
        #-7) Plot the result before and after unfolding
        self.plotUnfoldingEffect2D(regParam)
         
        #. . . . . . . . . . . . . . . . . . . .
        #-8) save the output spectra in a root file for paper figure plotting and sys err. evaluation
        self.saveResult(regParam-1, iteration, tag)
        self.saveResult(regParam, iteration, tag)
        self.saveResult(regParam+1, iteration, tag)
        
        #. . . . . . . . . . . . . . . . . . . .
        #-9) plot the kinematic efficiency. To so how much percentage of the projected region is covered
        #-by measured staatistic
        self.getKinEfficiency(bin)
 
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def provideExtPrior(self, tag, FileName=None, specialType=""):
        
        if FileName=="" or FileName==None: return 0
        inFileNameExt_RM  = "{}{}.root".format(self.inPathRM,FileName)
        inFileExt_RM      = ROOT.TFile(inFileNameExt_RM)
      
        #-Provide prior from file
        self.externalPrior = self.getResponseMatrix(0,"jetPtPart","zTPart", False,0, False, inFileExt_RM)
        #-Provide a flat prior
        if specialType=="Flat":
            #-Set all bin entries to 1
            for binx in range(self.externalPrior.GetNbinsX()+1):
                for biny in range(self.externalPrior.GetNbinsY()+1):
                    #self.externalPrior.SetBinContent(binx,biny,1) #Flat per bin
                    xWidth = self.externalPrior.GetXaxis().GetBinWidth(binx)
                    yWidth = self.externalPrior.GetYaxis().GetBinWidth(biny)
                    self.externalPrior.SetBinContent(binx,biny,xWidth*yWidth) #Flat per Area
                    #self.externalPrior.SetBinContent(binx,biny,yWidth) #Flat in zT
        
        
        #-update output path
        #-Necessary to update the output path so that the prior is saved in the correct folder
        if self.isPrompt:
            self.outPath =  self.outpathBase+"/Prompt{}_0/".format(tag)
        else:
            self.outPath =  self.outpathBase+"/NonPrompt{}_0/".format(tag)
        if not os.path.exists(self.outPath):
            os.makedirs(self.outPath)
         
        self.plotHist(self.externalPrior,"h_extPrior","colz","#it{p}_{T,jet} (GeV/#it{c})","#it{z}_{T}")

    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def prepareRMWeights(self, regParam, iteration):
        #This function finds the weights that need to be used while filling the RM
      
        hweights=None
        if self.externalPrior and iteration==0:
            #-Build the 2D ratio between original prior and externally provided prior
            hweights = self.externalPrior.Clone("hweightsForPrior_Ext")
            hweights.Divide(self.hOriginalPrior)
            #-scale by maximum bin
            #hweights.Scale(1.0/hweights.Integral())
            oldMax = self.hOriginalPrior.GetMaximum()
            newMax = hweights.GetMaximum()
            hweights.Scale(oldMax/newMax)
            self.plotHist(hweights,"hweights_extPrior","colz","#it{p}_{T,jet} (GeV/#it{c})","#it{z}_{T}")
        else:
            #-This version will adapt the prior to the mesured spectrum and repeat the unfolding
            #-with the improved RM
            #-Build the 2D ratio between prior and reconstructed data after unfolding
            hweights = self.unfoldedArr2D[regParam].Clone("hweightsForPrior_{}".format(iteration))
            hweights.Divide(self.hOriginalPrior)
            #-scale by maximum bin
            #hweights.Scale(1.0/hweights.Integral())
            oldMax = self.hOriginalPrior.GetMaximum()
            newMax = hweights.GetMaximum()
            hweights.Scale(oldMax/newMax)
            self.plotHist(hweights,"hweights_Iter{}".format(iteration),"colz","#it{p}_{T,jet} (GeV/#it{c})","#it{z}_{T}")
            
        return hweights
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getEventWeight(self, hweights, zTValue, pTValue):
        weight=1
        
        xBin = hweights.GetXaxis().FindBin(pTValue)
        yBin = hweights.GetYaxis().FindBin(zTValue)
        weight = hweights.GetBinContent(xBin,yBin)
        

        #print("zt: {}, pT: {}, weight: {}".format(zTValue, pTValue,weight))
        return weight

    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def RefoldingTest(self, ptBin, nIter, histo, RM1, RM2):
        
        #-unfold the measured spectrum with half the RM statistic and refold it
        #-back with the other half of the RM statistic.
        unfoldBayes1     = ROOT.RooUnfoldBayes(RM1, histo, nIter)
        unfoldedSpectrum = unfoldBayes1.Hreco(self.errorType)
        refoldedSpectrum = RM2.ApplyToTruth(unfoldedSpectrum)
        
        #-Draw the result
        c = ROOT.TCanvas("c","c: pT",800,850)
        c.cd()
        c.cd().SetLeftMargin(0.15)
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0)
        pad1.SetLeftMargin(0.15)
        pad1.SetRightMargin(0.05)
        pad1.SetTopMargin(0.05)
        #pad1.SetLogy()
        pad1.Draw()
        pad1.cd()
     
        histo.GetYaxis().SetTitle("counts per bin")
        histo.GetYaxis().SetTitleSize(0.06)
        histo.GetYaxis().SetLabelFont(43)
        histo.GetYaxis().SetLabelSize(20)
        histo.SetNdivisions(505)
        histo.SetLineColor(4)
        histo.SetLineWidth(3)
        histo.SetLineStyle(1)
        histo.Draw("hist E")
        refoldedSpectrum.SetLineColor(1)
        refoldedSpectrum.SetLineWidth(2)
        refoldedSpectrum.SetLineStyle(1)
        refoldedSpectrum.DrawCopy("hist E same")
     
        
        leg2 = ROOT.TLegend(0.2,0.8,0.5,0.93,"")
        leg2.SetFillColor(10)
        leg2.SetBorderSize(0)
        leg2.SetFillStyle(0)
        leg2.SetTextSize(0.04)
        leg2.AddEntry(refoldedSpectrum,   "measured spectrum unfolded (RM1)+refolded (RM2)", "l")
        leg2.AddEntry(histo,  "measured spectrum, detLvlv", "l")
        leg2.Draw("same")
  
        c.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.35)
        pad2.SetLeftMargin(0.15)
        pad2.SetRightMargin(0.05)
        pad2.Draw()
        pad2.cd()

        refoldedSpectrumC = refoldedSpectrum.Clone("{}Clone".format(refoldedSpectrum.GetName()))
        refoldedSpectrumC.Divide(histo)
        refoldedSpectrumC.SetMarkerStyle(21)
        refoldedSpectrumC.SetMarkerColor(2)

        refoldedSpectrumC.GetXaxis().SetTitleSize(30)
        refoldedSpectrumC.GetXaxis().SetTitleFont(43)
        refoldedSpectrumC.GetXaxis().SetTitleOffset(3.) #was 4
        refoldedSpectrumC.GetXaxis().SetLabelFont(43)
        refoldedSpectrumC.GetXaxis().SetLabelSize(20)
        refoldedSpectrumC.GetYaxis().SetRangeUser(0,2)
        refoldedSpectrumC.GetYaxis().SetTitle("Ratio to orig, nIter={}".format(nIter))
        refoldedSpectrumC.GetYaxis().SetTitleSize(20)
        refoldedSpectrumC.GetYaxis().SetTitleFont(43)
        refoldedSpectrumC.GetYaxis().SetTitleOffset(2.2)
        refoldedSpectrumC.GetYaxis().SetLabelFont(43)
        refoldedSpectrumC.GetYaxis().SetLabelSize(20)
        refoldedSpectrumC.GetYaxis().SetNdivisions(505)
        refoldedSpectrumC.GetXaxis().SetTitle("z_{T}") #ELI make this flexible
        refoldedSpectrumC.DrawCopy("hist E")

        lines = self.drawLines()
        for line in lines:
            line.Draw("same")
 
        c.cd()
  
        FileName="{}RefoldingTestBin{}_nIter{}_{}.png".format(self.outPath,ptBin,nIter,self.figureTag)
        #print("Safe figure under: {}".format(FileName))
        c.SaveAs(FileName)
        c.Close()
        
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def RefoldingTest2D(self, nIter, RM1, RM2):
        
        #-unfold the measured spectrum with half the RM statistic and refold it
        #-back with the other half of the RM statistic.
        unfoldBayes1     = ROOT.RooUnfoldBayes(RM1, self.measuredSpectra2D, nIter)
        unfoldedSpectrum = unfoldBayes1.Hreco(self.errorType)
        refoldedSpectrum = RM2.ApplyToTruth(unfoldedSpectrum)
        
        #-Draw the result
        #The result is a 2D histogram with x-Axis==Jet pT, y-Axis==zT
        #make a canvas pad for each pT bin that was unfolded.
        npTBins = len(self.pTBinsArrayTruth)-1
        c = ROOT.TCanvas("c","c: pT",800*npTBins,850)
        c.Divide(npTBins,2,0,0)
        for ptBin in range(1,npTBins+1):
            c.cd(ptBin)
            c.cd(ptBin).SetLeftMargin(0.15)
            c.cd(ptBin).SetRightMargin(0.05)
            c.cd(ptBin).SetTopMargin(0.05)
            c.cd(ptBin).SetBottomMargin(0)
            OrigPTBinProjection      = self.measuredSpectra2D.ProjectionY("_pyOrig",ptBin,ptBin)
            RefoldedPTBinProjection  = refoldedSpectrum.ProjectionY("_pyRefolded",ptBin,ptBin)
            OrigPTBinProjection.GetYaxis().SetTitle("counts per bin")
            OrigPTBinProjection.GetYaxis().SetTitleSize(0.06)
            OrigPTBinProjection.GetYaxis().SetLabelFont(43)
            OrigPTBinProjection.GetYaxis().SetLabelSize(20)
            OrigPTBinProjection.SetNdivisions(505)
            OrigPTBinProjection.SetLineColor(4)
            OrigPTBinProjection.SetLineWidth(3)
            OrigPTBinProjection.SetLineStyle(1)
            OrigPTBinProjection.DrawCopy("hist E")
            RefoldedPTBinProjection.SetLineColor(1)
            RefoldedPTBinProjection.SetLineWidth(2)
            RefoldedPTBinProjection.SetLineStyle(1)
            RefoldedPTBinProjection.DrawCopy("hist E same")
     
        
            leg2 = ROOT.TLegend(0.2,0.8,0.5,0.93,"")
            leg2.SetFillColor(10)
            leg2.SetBorderSize(0)
            leg2.SetFillStyle(0)
            leg2.SetTextSize(0.04)
            leg2.AddEntry(RefoldedPTBinProjection,   "measured spectrum unfolded (RM1)+refolded (RM2)", "l")
            leg2.AddEntry(OrigPTBinProjection,  "measured spectrum, detLvlv", "l")
            leg2.Draw("same")
  
    
            c.cd(npTBins+ptBin)
            c.cd(npTBins+ptBin).SetLeftMargin(0.15)
            c.cd(npTBins+ptBin).SetRightMargin(0.05)
            c.cd(npTBins+ptBin).SetTopMargin(0)
            c.cd(npTBins+ptBin).SetBottomMargin(0.35)
     

            refoldedSpectrumC = RefoldedPTBinProjection.Clone("{}1DClone".format(refoldedSpectrum.GetName()))
            refoldedSpectrumC.Divide(OrigPTBinProjection)
            refoldedSpectrumC.SetMarkerStyle(21)
            refoldedSpectrumC.SetMarkerColor(2)

            refoldedSpectrumC.GetXaxis().SetTitleSize(30)
            refoldedSpectrumC.GetXaxis().SetTitleFont(43)
            refoldedSpectrumC.GetXaxis().SetTitleOffset(3.) #was 4
            refoldedSpectrumC.GetXaxis().SetLabelFont(43)
            refoldedSpectrumC.GetXaxis().SetLabelSize(20)
            refoldedSpectrumC.GetYaxis().SetRangeUser(0,2)
            refoldedSpectrumC.GetYaxis().SetTitle("Ratio to orig, nIter={}".format(nIter))
            refoldedSpectrumC.GetYaxis().SetTitleSize(20)
            refoldedSpectrumC.GetYaxis().SetTitleFont(43)
            refoldedSpectrumC.GetYaxis().SetTitleOffset(2.2)
            refoldedSpectrumC.GetYaxis().SetLabelFont(43)
            refoldedSpectrumC.GetYaxis().SetLabelSize(20)
            refoldedSpectrumC.GetYaxis().SetNdivisions(505)
            refoldedSpectrumC.GetXaxis().SetTitle("z_{T}") #ELI make this flexible
            refoldedSpectrumC.DrawCopy("hist E")

            lines = self.drawLines()
            for line in lines:
                cLine = line.Clone("c")
                cLine.Draw("same")
 

 
        c.cd()
        
        
        FileName="{}RefoldingTest2D_nIter{}{}.png".format(self.outPath,nIter,self.figureTag)
        #print("Safe figure under: {}".format(FileName))
        c.SaveAs(FileName)
        c.Close()

    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def UnfoldingTest2D(self, response, nIter):
    
        TH2DRM = response.Hresponse()

        # Get the truth-level jet spectrum (matched) from response matrix
        hSpectrumTruePerBin  = response.Htruth()
        hSpectrumTruePerBin.SetName("hSpectrumTruePerBin")
    
        hSpectrumMCDetPerBin = response.Hmeasured()
        hSpectrumMCDetPerBin.SetName("hSpectrumMCDetPerBin")

        #-This is a closure test
        #-Take MC det lvl spectrum to see if one can reliably unfold to the gen lvl spectrum
        bayesUnfold = RooUnfoldBayes(response, hSpectrumMCDetPerBin, nIter)
        hSpectrumMCDetPerBinUnfolded = bayesUnfold.Hreco(self.errorType) # Produces the truth distribution, with errors, PerBin (will scale by bin width below)
        
        #-Draw the result
        #The result is a 2D histogram with x-Axis==Jet pT, y-Axis==zT
        #make a canvas pad for each pT bin that was unfolded.
        npTBins = len(self.pTBinsArrayTruth)-1
        
        leg2 = ROOT.TLegend(0.2,0.8,0.5,0.93,"")
        leg2.SetFillColor(10)
        leg2.SetBorderSize(0)
        leg2.SetFillStyle(0)
        leg2.SetTextSize(0.04)
            
        c = ROOT.TCanvas("c","c: pT",800*npTBins,850)
        c.Divide(npTBins,2,0,0)
        for ptBin in range(1,npTBins+1):
            c.cd(ptBin)
            c.cd(ptBin).SetLeftMargin(0.15)
            c.cd(ptBin).SetRightMargin(0.05)
            c.cd(ptBin).SetTopMargin(0.05)
            c.cd(ptBin).SetBottomMargin(0)
            MCTruthPTBinProjection             = hSpectrumTruePerBin.ProjectionY("_pyOrig",ptBin,ptBin)
            MCMeasuredUnfoldedPTBinProjection  = hSpectrumMCDetPerBinUnfolded.ProjectionY("_pyRefolded",ptBin,ptBin)

     
            MCMeasuredUnfoldedPTBinProjection.GetYaxis().SetTitle("counts per bin")
            MCMeasuredUnfoldedPTBinProjection.GetYaxis().SetTitleSize(0.06)
            MCMeasuredUnfoldedPTBinProjection.GetYaxis().SetLabelFont(43)
            MCMeasuredUnfoldedPTBinProjection.GetYaxis().SetLabelSize(20)
            hSpectrumMCDetPerBinUnfolded.SetNdivisions(505)
            MCMeasuredUnfoldedPTBinProjection.SetLineColor(4)
            MCMeasuredUnfoldedPTBinProjection.SetLineWidth(3)
            MCMeasuredUnfoldedPTBinProjection.SetLineStyle(1)
            MCMeasuredUnfoldedPTBinProjection.DrawCopy("hist E")
            MCTruthPTBinProjection.SetLineColor(1)
            MCTruthPTBinProjection.SetLineWidth(2)
            MCTruthPTBinProjection.SetLineStyle(1)
            MCTruthPTBinProjection.DrawCopy("hist E same")
         
            if ptBin==1:
                leg2.AddEntry(MCTruthPTBinProjection,   "MC Gen. lvl spectrum", "l")
                leg2.AddEntry(hSpectrumMCDetPerBinUnfolded,  "MC Det. lvl sprectrum unfolded", "l")
            leg2.Draw("same")
  
            c.cd(npTBins+ptBin)
            c.cd(npTBins+ptBin).SetLeftMargin(0.15)
            c.cd(npTBins+ptBin).SetRightMargin(0.05)
            c.cd(npTBins+ptBin).SetTopMargin(0)
            c.cd(npTBins+ptBin).SetBottomMargin(0.35)
  

            MCMeasuredUnfoldedPTBinProjectionC = MCMeasuredUnfoldedPTBinProjection.Clone("{}Clone".format(MCMeasuredUnfoldedPTBinProjection.GetName()))
            MCMeasuredUnfoldedPTBinProjectionC.Divide(MCTruthPTBinProjection)
            MCMeasuredUnfoldedPTBinProjectionC.SetMarkerStyle(21)
            MCMeasuredUnfoldedPTBinProjectionC.SetMarkerColor(2)

            MCMeasuredUnfoldedPTBinProjectionC.GetXaxis().SetTitleSize(30)
            MCMeasuredUnfoldedPTBinProjectionC.GetXaxis().SetTitleFont(43)
            MCMeasuredUnfoldedPTBinProjectionC.GetXaxis().SetTitleOffset(3.) #was 4
            MCMeasuredUnfoldedPTBinProjectionC.GetXaxis().SetLabelFont(43)
            MCMeasuredUnfoldedPTBinProjectionC.GetXaxis().SetLabelSize(20)
            #MCMeasuredUnfoldedPTBinProjectionC.GetYaxis().SetRangeUser(0,2)
            MCMeasuredUnfoldedPTBinProjectionC.GetYaxis().SetRangeUser(0.5,1.5)
            MCMeasuredUnfoldedPTBinProjectionC.GetYaxis().SetTitle("Ratio Unfolded to MC Gen lvl, nIter={}".format(nIter))
            MCMeasuredUnfoldedPTBinProjectionC.GetYaxis().SetTitleSize(20)
            MCMeasuredUnfoldedPTBinProjectionC.GetYaxis().SetTitleFont(43)
            MCMeasuredUnfoldedPTBinProjectionC.GetYaxis().SetTitleOffset(2.2)
            MCMeasuredUnfoldedPTBinProjectionC.GetYaxis().SetLabelFont(43)
            MCMeasuredUnfoldedPTBinProjectionC.GetYaxis().SetLabelSize(20)
            MCMeasuredUnfoldedPTBinProjectionC.GetYaxis().SetNdivisions(505)
            MCMeasuredUnfoldedPTBinProjectionC.GetXaxis().SetTitle("z_{T}") #ELI make this flexible
            MCMeasuredUnfoldedPTBinProjectionC.DrawCopy("hist E")

            lines = self.drawLines()
            for line in lines:
                cLine = line.Clone("c")
                cLine.Draw("same")
        c.cd()
  
        FileName="{}UnfoldingClosureTest2D_nIter{}{}.png".format(self.outPath,nIter,self.figureTag)
        #print("Safe figure under: {}".format(FileName))
        c.SaveAs(FileName)
        c.Close()
    
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def UnfoldingTest(self, response, nIter):
    
        TH2DRM = response.Hresponse()
        #self.plotHist(TH2DRM,"TTT{}".format(TH2DRM.GetName()),"colz")

        # Get the truth-level jet spectrum (matched) from response matrix (already re-binned)
        hSpectrumTruePerBin  = TH2DRM.ProjectionY("_py",1,TH2DRM.GetNbinsX()) # Do exclude under and overflow bins
        hSpectrumTruePerBin.SetName("hSpectrumTruePerBin")
    
        hSpectrumMCDetPerBin = TH2DRM.ProjectionX("_px",1,TH2DRM.GetNbinsY()) # Do exclude under and overflow bins
        hSpectrumMCDetPerBin.SetName("hSpectrumMCDetPerBin")

        #-This is a closure test
        #-Take MC det lvl spectrum to see if one can reliably unfold to the gen lvl spectrum
        bayesUnfold = RooUnfoldBayes(response, hSpectrumMCDetPerBin, nIter)
        #bayesUnfold = ROOT.RooUnfoldSvd(response, hSpectrumMCDetPerBin, nIter)
        hSpectrumMCDetPerBinUnfolded = bayesUnfold.Hreco(self.errorType) # Produces the truth distribution, with errors, PerBin (will scale by bin width below)
        
 
        #-Draw the result
        c = ROOT.TCanvas("c","c: pT",800,850)
        c.cd()
        c.cd().SetLeftMargin(0.15)
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0)
        pad1.SetLeftMargin(0.15)
        pad1.SetRightMargin(0.05)
        pad1.SetTopMargin(0.05)
        #pad1.SetLogy()
        pad1.Draw()
        pad1.cd()
     
        hSpectrumMCDetPerBinUnfolded.GetYaxis().SetTitle("counts per bin")
        hSpectrumMCDetPerBinUnfolded.GetYaxis().SetTitleSize(0.06)
        hSpectrumMCDetPerBinUnfolded.GetYaxis().SetLabelFont(43)
        hSpectrumMCDetPerBinUnfolded.GetYaxis().SetLabelSize(20)
        hSpectrumMCDetPerBinUnfolded.SetNdivisions(505)
        hSpectrumMCDetPerBinUnfolded.SetLineColor(4)
        hSpectrumMCDetPerBinUnfolded.SetLineWidth(3)
        hSpectrumMCDetPerBinUnfolded.SetLineStyle(1)
        hSpectrumMCDetPerBinUnfolded.Draw("hist E")
        hSpectrumTruePerBin.SetLineColor(1)
        hSpectrumTruePerBin.SetLineWidth(2)
        hSpectrumTruePerBin.SetLineStyle(1)
        hSpectrumTruePerBin.DrawCopy("hist E same")
     
        leg2 = ROOT.TLegend(0.2,0.8,0.5,0.93,"")
        leg2.SetFillColor(10)
        leg2.SetBorderSize(0)
        leg2.SetFillStyle(0)
        leg2.SetTextSize(0.04)
        leg2.AddEntry(hSpectrumTruePerBin,   "MC Gen. lvl spectrum", "l")
        leg2.AddEntry(hSpectrumMCDetPerBinUnfolded,  "MC Det. lvl sprectrum unfolded", "l")
        leg2.Draw("same")
  
        c.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.35)
        pad2.SetLeftMargin(0.15)
        pad2.SetRightMargin(0.05)
        pad2.Draw()
        pad2.cd()

        hSpectrumMCDetPerBinUnfoldedC = hSpectrumMCDetPerBinUnfolded.Clone("{}Clone".format(hSpectrumMCDetPerBinUnfolded.GetName()))
        hSpectrumMCDetPerBinUnfoldedC.Divide(hSpectrumTruePerBin)
        hSpectrumMCDetPerBinUnfoldedC.SetMarkerStyle(21)
        hSpectrumMCDetPerBinUnfoldedC.SetMarkerColor(2)

        hSpectrumMCDetPerBinUnfoldedC.GetXaxis().SetTitleSize(30)
        hSpectrumMCDetPerBinUnfoldedC.GetXaxis().SetTitleFont(43)
        hSpectrumMCDetPerBinUnfoldedC.GetXaxis().SetTitleOffset(3.) #was 4
        hSpectrumMCDetPerBinUnfoldedC.GetXaxis().SetLabelFont(43)
        hSpectrumMCDetPerBinUnfoldedC.GetXaxis().SetLabelSize(20)
        hSpectrumMCDetPerBinUnfoldedC.GetYaxis().SetRangeUser(0,2)
        hSpectrumMCDetPerBinUnfoldedC.GetYaxis().SetTitle("Ratio Unfolded to MC Gen lvl, nIter={}".format(nIter))
        hSpectrumMCDetPerBinUnfoldedC.GetYaxis().SetTitleSize(20)
        hSpectrumMCDetPerBinUnfoldedC.GetYaxis().SetTitleFont(43)
        hSpectrumMCDetPerBinUnfoldedC.GetYaxis().SetTitleOffset(2.2)
        hSpectrumMCDetPerBinUnfoldedC.GetYaxis().SetLabelFont(43)
        hSpectrumMCDetPerBinUnfoldedC.GetYaxis().SetLabelSize(20)
        hSpectrumMCDetPerBinUnfoldedC.GetYaxis().SetNdivisions(505)
        hSpectrumMCDetPerBinUnfoldedC.GetXaxis().SetTitle("z_{T}") #ELI make this flexible
        hSpectrumMCDetPerBinUnfoldedC.DrawCopy("hist E")

        lines = self.drawLines()
        for line in lines:
            line.Draw("same")
 
        c.cd()
  
        FileName="{}UnfoldingClosureTest_nIter{}{}.png".format(self.outPath,nIter,self.figureTag)
        #print("Safe figure under: {}".format(FileName))
        c.SaveAs(FileName)
        c.Close()
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def plotPrior2D(self, nIter, RooUnfoldRM):
   
        ROOT.gStyle.SetOptTitle(0)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(0)
 
        if self.dictKey=="Psi2S":
            xTitle="p_{T}^{#psi(2S)}/p_{T}^{jet}"
        else:
            xTitle="p_{T}^{#chi_{c1}(3872)}/p_{T}^{jet}"
  
        mainUnfolded   = self.unfoldedArr2D[nIter-1].Clone("TUnfold") #Array position 0 is iteration number 1,
        original       = self.measuredSpectra2D.Clone("TOrig")
        prior          = RooUnfoldRM.Htruth().Clone("TTruth")
        
        colorArray          = [ROOT.kRed-9,ROOT.kGreen-9,ROOT.kGreen-8,ROOT.kBlue-9, ROOT.kGreen, ROOT.kSpring+8,ROOT.kSpring]
        colorArrayFinalFigs = [ROOT.kAzure, ROOT.kAzure-4, ROOT.kCyan-6,  ROOT.kGreen-3, ROOT.kTeal-6,ROOT.kGreen+4,1,1]
        maxArrayOrig   = []
        maxArrayUnfold = []
        hArrayOrig   = []
        hArrayUnfold = []
        leg1 = ROOT.TLegend(0.2,0.65,0.5,0.85,"")
        leg1.SetFillColor(10)
        leg1.SetBorderSize(0)
        leg1.SetFillStyle(0)
        leg1.SetTextSize(0.06)
        #-Draw the result
        #The result is a 2D histogram with x-Axis==Jet pT, y-Axis==zT
        #make a canvas pad for each pT bin that was unfolded.
        npTBins = len(self.pTBinsArrayTruth)-1
        if self.dictKey=="Psi2S":
            c1 = ROOT.TCanvas("c1","c1: pT",800*3,800*2)
            c1.Divide(3,2,0,0)
        else:
            c1 = ROOT.TCanvas("c1","c1: pT",800*npTBins,800)
            c1.Divide(npTBins,1,0,0)
        for ptBin in range(1,npTBins+1):
            c1.cd(ptBin)
            c1.cd(ptBin).SetLeftMargin(0.15)
            c1.cd(ptBin).SetRightMargin(0.05)
            c1.cd(ptBin).SetTopMargin(0.05)
            c1.cd(ptBin).SetBottomMargin(0.15)

            mainUnfoldedProj  = mainUnfolded.ProjectionY("_pyMainUnfolded_Bin{}".format(ptBin),ptBin,ptBin)
            originalProj      = original.ProjectionY("_pyOrig_Bin{}".format(ptBin),ptBin,ptBin)
            priorProj         = prior.ProjectionY("_pyPrior_Bin{}".format(ptBin),ptBin,ptBin)
            mainUnfoldedS     = self.scaleHistogram(mainUnfoldedProj) # divides entries by bin width
            originalS         = self.scaleHistogram(originalProj)     # divides entries by bin width
            priorS            = self.scaleHistogram(priorProj)        # divides entries by bin width
            if priorS.Integral()>0:
                priorS.Scale(mainUnfoldedS.Integral()/priorS.Integral())
  
            originalS.SetNdivisions(505)
            maxList = [originalS.GetMaximum(),priorS.GetMaximum(),mainUnfoldedS.GetMaximum()]

 
            originalS.GetYaxis().SetTitle("dN/dz_{T}")
            originalS.GetYaxis().SetTitleSize(0.06)
            originalS.GetYaxis().SetLabelFont(43)
            originalS.GetYaxis().SetLabelSize(20)
            originalS.GetXaxis().SetTitle(xTitle)
            originalS.GetXaxis().SetTitleSize(0.05)
            originalS.GetXaxis().SetRangeUser(0,1)
            originalS.GetXaxis().SetNdivisions(405)
  
            originalS.GetYaxis().SetRangeUser(0,max(maxList)*1.3)
            originalS.SetLineColor(4)
            originalS.SetMarkerColor(4)
            originalS.SetLineWidth(3)
            originalS.SetMarkerStyle(20)
            originalS.SetLineStyle(1)
            originalS.DrawCopy("hist E")
            mainUnfoldedS.SetLineColor(2)#1
            mainUnfoldedS.SetMarkerColor(2)
            mainUnfoldedS.SetLineWidth(2)
            mainUnfoldedS.SetMarkerStyle(20)
            mainUnfoldedS.SetLineStyle(1)
            mainUnfoldedS.DrawCopy("same hist E")
            
            priorS.SetLineColor(ROOT.kSpring+3)#1
            priorS.SetMarkerColor(ROOT.kSpring+3)#1
            priorS.SetLineWidth(2)
            priorS.SetLineStyle(1)
            priorS.SetMarkerStyle(20)
            priorS.DrawCopy("hist E same")

            maxArrayOrig.append(originalS.GetMaximum())
            maxArrayUnfold.append(mainUnfoldedS.GetMaximum())
            hArrayOrig.append(originalS)
            hArrayUnfold.append(mainUnfoldedS)
            
            if ptBin==1:
                leg1.AddEntry(originalS, "measured", "l")
                leg1.AddEntry(mainUnfoldedS, "unfolded nIter={}".format(nIter), "l")
                leg1.AddEntry(priorS, "prior", "l")
            leg1.Draw("same")

            textFit = ROOT.TLatex()
            textFit.SetTextSize(0.06)
            textFit.SetNDC()
            textFit.DrawLatex(0.2,0.87," #it{p}_{T}^{jet}: %d-%d GeV"%(self.pTBinsArrayTruth[ptBin-1],self.pTBinsArrayTruth[ptBin]))

        c1.cd()
        FileName="{}PriorComparision2D_Yaxis_nIter{}{}.png".format(self.outPath, nIter,self.figureTag)
        c1.SaveAs(FileName)
        c1.Close()
 
        leg2 = ROOT.TLegend(0.5,0.7,0.8,0.93,"")
        leg2.SetFillColor(10)
        leg2.SetBorderSize(0)
        leg2.SetFillStyle(0)
        leg2.SetTextSize(0.06)
        
        nzTBins = len(self.zBinsArrayTruth)-1
        
        if self.dictKey=="Psi2S":
            c2 = ROOT.TCanvas("c1","c1: pT",800*4,850*3)
            c2.Divide(4,3)
        else:
            c2 = ROOT.TCanvas("c1","c1: pT",800*3,850*2)
            c2.Divide(3,2)

        for zTBin in range(1,nzTBins+1):
    
            c2.cd(zTBin)
            c2.cd(zTBin).SetLeftMargin(0.15)
            c2.cd(zTBin).SetRightMargin(0.05)
            c2.cd(zTBin).SetTopMargin(0.05)
            c2.cd(zTBin).SetBottomMargin(0.15)

            mainUnfoldedProj  = mainUnfolded.ProjectionX("_pxMainUnfolded_Bin{}".format(zTBin),zTBin,zTBin)
            originalProj      = original.ProjectionX("_pxOrig_Bin{}".format(zTBin),zTBin,zTBin)
            priorProj         = prior.ProjectionX("_pxPrior_Bin{}".format(zTBin),zTBin,zTBin)
            mainUnfoldedS     = self.scaleHistogram(mainUnfoldedProj) # divides entries by bin width
            originalS         = self.scaleHistogram(originalProj)     # divides entries by bin width
            priorS            = self.scaleHistogram(priorProj)        # divides entries by bin width

            if priorS.Integral()>0:
                priorS.Scale(mainUnfoldedS.Integral()/priorS.Integral())
 
            maxList = [originalS.GetMaximum(),priorS.GetMaximum(),mainUnfoldedS.GetMaximum()]
            mainUnfoldedS.GetYaxis().SetRangeUser(0,max(maxList)*1.3)

            mainUnfoldedS.GetYaxis().SetTitle("dN/dp_{T}")
            mainUnfoldedS.GetYaxis().SetTitleSize(0.06)
            mainUnfoldedS.GetYaxis().SetLabelFont(43)
            mainUnfoldedS.GetYaxis().SetLabelSize(20)
            mainUnfoldedS.GetXaxis().SetTitle("jet p_{T}")
            mainUnfoldedS.GetXaxis().SetTitleSize(0.05)
            mainUnfoldedS.GetXaxis().SetRangeUser(0,1)
            mainUnfoldedS.GetXaxis().SetNdivisions(405)

            mainUnfoldedS.SetNdivisions(505)
            mainUnfoldedS.SetLineColor(2)#1
            mainUnfoldedS.SetMarkerColor(2)
            mainUnfoldedS.SetLineWidth(2)
            mainUnfoldedS.SetMarkerStyle(20)
            mainUnfoldedS.SetLineStyle(1)
            mainUnfoldedS.DrawCopy("hist E")
            originalS.SetLineColor(4)
            originalS.SetMarkerColor(4)
            originalS.SetLineWidth(3)
            originalS.SetMarkerStyle(20)
            originalS.SetLineStyle(1)
            originalS.DrawCopy("hist E same")
            
            priorS.SetLineColor(ROOT.kSpring+3)#1
            priorS.SetMarkerColor(ROOT.kSpring+3)#1
            priorS.SetLineWidth(2)
            priorS.SetLineStyle(1)
            priorS.SetMarkerStyle(20)
            priorS.DrawCopy("hist E same")
            
            if ptBin==1:
                leg2.AddEntry(originalS, "measured", "l")
                leg2.AddEntry(mainUnfoldedS, "unfolded nIter={}".format(nIter), "l")
                leg2.AddEntry(priorS, "prior", "l")
            leg2.Draw("same")

            textFit = ROOT.TLatex()
            textFit.SetTextSize(0.06)
            textFit.SetNDC()
            textFit.DrawLatex(0.4,0.8," #it{z}_{T}: %d-%d GeV"%(self.zBinsArrayTruth[zTBin-1],self.zBinsArrayTruth[zTBin]))

        c2.cd()
        FileName="{}PriorComparision2D_Xaxis_nIter{}{}.png".format(self.outPath, nIter,self.figureTag)
        c2.SaveAs(FileName)
        c2.Close()
        
        '''
        #- . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . -
        #-Plot summary after unfolding
        c3 = ROOT.TCanvas("c2","c: hist",500*2,450*2)
        c3.cd()
        ROOT.TGaxis.SetMaxDigits(3)
        
        # Set pad and histo arrangement
        myPad3 = ROOT.TPad("myPad", "The pad",0,0,1,1)
        #myPad2.SetLeftMargin(0.21)
        myPad3.SetLeftMargin(0.15)
        myPad3.SetTopMargin(0.06)
        myPad3.SetRightMargin(0.04)
        myPad3.SetBottomMargin(0.15)
        myPad3.SetTicks()
        myPad3.Draw()
        myPad3.cd()

        maximum = max(maxArrayUnfold)
        #print("This is the maximum: {}".format(max))
        myBlankHisto3 = ROOT.TH1F("myBlankHisto3","Blank Histogram",20,0,1)
        myBlankHisto3.GetXaxis().SetNdivisions(505)
        myBlankHisto3.SetXTitle(xTitle)
        myBlankHisto3.GetXaxis().SetTitleSize(0.05)
        myBlankHisto3.GetXaxis().SetRangeUser(0,1)
        myBlankHisto3.GetXaxis().SetNdivisions(405)
        myBlankHisto3.GetYaxis().SetTitleOffset(1.35)
        myBlankHisto3.GetYaxis().SetTitleSize(0.055)
        myBlankHisto3.SetLineColor(0)

        #Absolute yield
        myBlankHisto3.SetYTitle("dN/dz_{T} unfolded")
        if len(self.zBinsArrayTruth)==len(self.zBinsArrayDet):
            myBlankHisto3.GetYaxis().SetRangeUser(0,0.45) #When drawn normaluzed
        else:
            myBlankHisto3.GetYaxis().SetRangeUser(0,0.8) #When drawn normaluzed
        myBlankHisto3.Draw("E")

        #-Legend
        #-about different contributions
        myLegend3 = ROOT.TLegend(0.2,0.7,0.4,0.9)
        if len(hArrayUnfold)==4:
          myLegend3 = ROOT.TLegend(0.2,0.7,0.4,0.9)
        if len(hArrayUnfold)==6:
          myLegend3 = ROOT.TLegend(0.2,0.6,0.4,0.9)

        myLegend3.SetTextFont(42)
        myLegend3.SetBorderSize(0)
        myLegend3.SetFillStyle(0)
        myLegend3.SetFillColor(0)
        myLegend3.SetMargin(0.25)
        myLegend3.SetTextSize(0.04)
        myLegend3.AddEntry(hArrayUnfold[0]," #it{p}_{T}^{jet}:","")
  
        MarkerScale=1.6
        for i in range(0,len(hArrayUnfold)):
            hArrayUnfold[i].SetMarkerSize(1.3*MarkerScale)
            hArrayUnfold[i].SetMarkerStyle(20)
            hArrayUnfold[i].SetMarkerColor(colorArrayFinalFigs[i])
            hArrayUnfold[i].SetLineStyle(1)
            hArrayUnfold[i].SetLineWidth(2)
            hArrayUnfold[i].SetLineColor(colorArrayFinalFigs[i])
            #if normType==1:
            #  self.normalize(yieldArray[i])
            #hArrayOrig[i].Draw("same EP")
            hArrayUnfold[i].DrawNormalized("same EP")
            
            myLegend3.AddEntry(hArrayUnfold[i],"  %s-%s (GeV/%s)" % (self.pTBinsArrayTruth[i],self.pTBinsArrayTruth[i+1],"#it{c}"),"LP")

        myLegend3.Draw()
         
        c3.cd()
        FileName="{}UnfoldingEffect2DFinal.png".format(self.outPath)
        c3.SaveAs(FileName)
        c3.Close()
        '''
      
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def plotUnfoldingEffect2D(self, nIter):
   
        ROOT.gStyle.SetOptTitle(0)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(0)
 
        if self.dictKey=="Psi2S":
            xTitle="p_{T}^{#psi(2S)}/p_{T}^{jet}"
        else:
            xTitle="p_{T}^{#chi_{c1}(3872)}/p_{T}^{jet}"
  
        mainUnfolded   = self.unfoldedArr2D[nIter-1].Clone("TUnfold") #Array position 0 is iteration number 1,
        original       = self.measuredSpectra2D.Clone("TOrig")
        
        colorArray          = [ROOT.kRed-9,ROOT.kGreen-9,ROOT.kGreen-8,ROOT.kBlue-9, ROOT.kGreen, ROOT.kSpring+8,ROOT.kSpring]
        colorArrayFinalFigs = [ROOT.kAzure, ROOT.kAzure-4, ROOT.kCyan-6,  ROOT.kGreen-3, ROOT.kTeal-6,ROOT.kGreen+4,1,1]
        maxArrayOrig   = []
        maxArrayUnfold = []
        hArrayOrig   = []
        hArrayUnfold = []
        leg1 = ROOT.TLegend(0.2,0.7,0.5,0.85,"")
        leg1.SetFillColor(10)
        leg1.SetBorderSize(0)
        leg1.SetFillStyle(0)
        leg1.SetTextSize(0.06)
        #-Draw the result
        #The result is a 2D histogram with x-Axis==Jet pT, y-Axis==zT
        #make a canvas pad for each pT bin that was unfolded.
        npTBins = len(self.pTBinsArrayTruth)-1
        if self.dictKey=="Psi2S":
            c1 = ROOT.TCanvas("c1","c1: pT",800*3,800*2)
            c1.Divide(3,2)
        else:
            c1 = ROOT.TCanvas("c1","c1: pT",800*npTBins,800)
            c1.Divide(npTBins,1)
        for ptBin in range(1,npTBins+1):
            c1.cd(ptBin)
            c1.cd(ptBin).SetLeftMargin(0.15)
            c1.cd(ptBin).SetRightMargin(0.05)
            c1.cd(ptBin).SetTopMargin(0.05)
            c1.cd(ptBin).SetBottomMargin(0.15)
     
            mainUnfoldedProj  = mainUnfolded.ProjectionY("_pyMainUnfolded_Bin{}".format(ptBin),ptBin,ptBin)
            originalProj      = original.ProjectionY("_pyOrig_Bin{}".format(ptBin),ptBin,ptBin)
            mainUnfoldedS     = self.scaleHistogram(mainUnfoldedProj) #scale by bin width
            originalS         = self.scaleHistogram(originalProj)     #scale by bin width

            maxList = [originalS.GetMaximum(),mainUnfoldedS.GetMaximum()]
            originalS.GetYaxis().SetRangeUser(0,max(maxList)*1.3)

            originalS.GetYaxis().SetTitle("dN/dz_{T}")
            originalS.GetYaxis().SetTitleSize(0.06)
            originalS.GetYaxis().SetLabelFont(43)
            originalS.GetYaxis().SetLabelSize(20)
            originalS.GetXaxis().SetTitle(xTitle)
            originalS.GetXaxis().SetTitleSize(0.05)
            originalS.GetXaxis().SetRangeUser(0,1)
            originalS.GetXaxis().SetNdivisions(405)
            originalS.SetNdivisions(505)
            originalS.SetLineColor(4)
            originalS.SetLineWidth(3)
            originalS.SetLineStyle(1)
            originalS.DrawCopy("hist E")
            mainUnfoldedS.SetLineColor(2)#1
            mainUnfoldedS.SetLineWidth(2)
            mainUnfoldedS.SetLineStyle(1)
            mainUnfoldedS.DrawCopy("hist E same")

            maxArrayOrig.append(originalS.GetMaximum())
            maxArrayUnfold.append(mainUnfoldedS.GetMaximum())
            hArrayOrig.append(originalS)
            hArrayUnfold.append(mainUnfoldedS)
            
            if ptBin==1:
                leg1.AddEntry(originalS, "measured", "l")
                leg1.AddEntry(mainUnfoldedS, "unfolded nIter={}".format(nIter), "l")
            leg1.Draw("same")

            textFit = ROOT.TLatex()
            textFit.SetTextSize(0.06)
            textFit.SetNDC()
            textFit.DrawLatex(0.2,0.87," #it{p}_{T}^{jet}: %d-%d GeV"%(self.pTBinsArrayTruth[ptBin-1],self.pTBinsArrayTruth[ptBin]))
    
        c1.cd()
        FileName="{}UnfoldingEffect2D_nIter{}{}.png".format(self.outPath, nIter,self.figureTag)
        c1.SaveAs(FileName)
        c1.Close()
        
        #- . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . -
        #-Plot summary before unfolding for testing equivalent to: PlotFinalFigures.py
        c2 = ROOT.TCanvas("c2","c: hist",500*2,450*2)
        c2.cd()
        ROOT.TGaxis.SetMaxDigits(3)
        
        # Set pad and histo arrangement
        myPad2 = ROOT.TPad("myPad2", "The pad",0,0,1,1)
        #myPad2.SetLeftMargin(0.21)
        myPad2.SetLeftMargin(0.15)
        myPad2.SetTopMargin(0.06)
        myPad2.SetRightMargin(0.04)
        myPad2.SetBottomMargin(0.15)
        myPad2.SetTicks()
        myPad2.Draw()
        myPad2.cd()

        maximum = max(maxArrayOrig)
        #print("This is the maximum: {}".format(max))
        myBlankHisto2 = ROOT.TH1F("myBlankHisto2","Blank Histogram",20,0,1)
        myBlankHisto2.GetXaxis().SetNdivisions(505)

        myBlankHisto2.SetXTitle(xTitle)
        myBlankHisto2.GetXaxis().SetTitleSize(0.05)
        myBlankHisto2.GetXaxis().SetRangeUser(0,1)
        myBlankHisto2.GetXaxis().SetNdivisions(405)
        myBlankHisto2.GetYaxis().SetTitleOffset(1.35)
        myBlankHisto2.GetYaxis().SetTitleSize(0.055)
        myBlankHisto2.SetLineColor(0)

        #Absolute yield
        myBlankHisto2.SetYTitle("dN/dz_{T} measured")
        #myBlankHisto2.GetYaxis().SetRangeUser(0,maximum*1.2)
        myBlankHisto2.GetYaxis().SetRangeUser(0,0.45) #When drawn normaluzed
        myBlankHisto2.Draw("E")

        #-Legend
        #-about different contributions
        myLegend2 = ROOT.TLegend(0.2,0.7,0.4,0.9)
        if len(hArrayOrig)==4:
          myLegend2 = ROOT.TLegend(0.2,0.7,0.4,0.9)
        if len(hArrayOrig)==6:
          myLegend2 = ROOT.TLegend(0.2,0.6,0.4,0.9)

        myLegend2.SetTextFont(42)
        myLegend2.SetBorderSize(0)
        myLegend2.SetFillStyle(0)
        myLegend2.SetFillColor(0)
        myLegend2.SetMargin(0.25)
        myLegend2.SetTextSize(0.04)
        myLegend2.AddEntry(hArrayOrig[0]," #it{p}_{T}^{jet}:","")
  
        MarkerScale=1.6
        for i in range(0,len(hArrayOrig)):
            hArrayOrig[i].SetMarkerSize(1.3*MarkerScale)
            hArrayOrig[i].SetMarkerStyle(20)
            hArrayOrig[i].SetMarkerColor(colorArrayFinalFigs[i])
            hArrayOrig[i].SetLineStyle(1)
            hArrayOrig[i].SetLineWidth(2)
            hArrayOrig[i].SetLineColor(colorArrayFinalFigs[i])
            #hArrayOrig[i].Draw("same EP")
            hArrayOrig[i].DrawNormalized("same EP")
            
            myLegend2.AddEntry(hArrayOrig[i],"  %s-%s (GeV/%s)" % (self.pTBinsArrayDet[i],self.pTBinsArrayDet[i+1],"#it{c}"),"LP")

        myLegend2.Draw()
         
        c2.cd()
        FileName="{}UnfoldingEffect2DOrig_{}.png".format(self.outPath,self.figureTag)
        c2.SaveAs(FileName)
        c2.Close()
   
        #- . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . -
        #-Plot summary after unfolding
        c3 = ROOT.TCanvas("c2","c: hist",500*2,450*2)
        c3.cd()
        ROOT.TGaxis.SetMaxDigits(3)
        
        # Set pad and histo arrangement
        myPad3 = ROOT.TPad("myPad", "The pad",0,0,1,1)
        #myPad2.SetLeftMargin(0.21)
        myPad3.SetLeftMargin(0.15)
        myPad3.SetTopMargin(0.06)
        myPad3.SetRightMargin(0.04)
        myPad3.SetBottomMargin(0.15)
        myPad3.SetTicks()
        myPad3.Draw()
        myPad3.cd()

        maximum = max(maxArrayUnfold)
        #print("This is the maximum: {}".format(max))
        myBlankHisto3 = ROOT.TH1F("myBlankHisto3","Blank Histogram",20,0,1)
        myBlankHisto3.GetXaxis().SetNdivisions(505)
        myBlankHisto3.SetXTitle(xTitle)
        myBlankHisto3.GetXaxis().SetTitleSize(0.05)
        myBlankHisto3.GetXaxis().SetRangeUser(0,1)
        myBlankHisto3.GetXaxis().SetNdivisions(405)
        myBlankHisto3.GetYaxis().SetTitleOffset(1.35)
        myBlankHisto3.GetYaxis().SetTitleSize(0.055)
        myBlankHisto3.SetLineColor(0)

        #Absolute yield
        myBlankHisto3.SetYTitle("dN/dz_{T} unfolded")
        if len(self.zBinsArrayTruth)==len(self.zBinsArrayDet):
            myBlankHisto3.GetYaxis().SetRangeUser(0,0.45) #When drawn normaluzed
        else:
            myBlankHisto3.GetYaxis().SetRangeUser(0,0.8) #When drawn normaluzed
        myBlankHisto3.Draw("E")

        #-Legend
        #-about different contributions
        myLegend3 = ROOT.TLegend(0.2,0.7,0.4,0.9)
        if len(hArrayUnfold)==4:
          myLegend3 = ROOT.TLegend(0.2,0.7,0.4,0.9)
        if len(hArrayUnfold)==6:
          myLegend3 = ROOT.TLegend(0.2,0.6,0.4,0.9)

        myLegend3.SetTextFont(42)
        myLegend3.SetBorderSize(0)
        myLegend3.SetFillStyle(0)
        myLegend3.SetFillColor(0)
        myLegend3.SetMargin(0.25)
        myLegend3.SetTextSize(0.04)
        myLegend3.AddEntry(hArrayUnfold[0]," #it{p}_{T}^{jet}:","")
  
        MarkerScale=1.6
        for i in range(0,len(hArrayUnfold)):
            hArrayUnfold[i].SetMarkerSize(1.3*MarkerScale)
            hArrayUnfold[i].SetMarkerStyle(20)
            hArrayUnfold[i].SetMarkerColor(colorArrayFinalFigs[i])
            hArrayUnfold[i].SetLineStyle(1)
            hArrayUnfold[i].SetLineWidth(2)
            hArrayUnfold[i].SetLineColor(colorArrayFinalFigs[i])
            #if normType==1:
            #  self.normalize(yieldArray[i])
            #hArrayOrig[i].Draw("same EP")
            hArrayUnfold[i].DrawNormalized("same EP")
            
            myLegend3.AddEntry(hArrayUnfold[i],"  %s-%s (GeV/%s)" % (self.pTBinsArrayTruth[i],self.pTBinsArrayTruth[i+1],"#it{c}"),"LP")

        myLegend3.Draw()
         
        c3.cd()
        FileName="{}UnfoldingEffect2DFinal_{}.png".format(self.outPath,self.figureTag)
        c3.SaveAs(FileName)
        c3.Close()
      
        
        #- . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . -
        #-Plot summary after unfolding
        firstBin  = self.pTBinsArrayTruth.index(20)+1
        lastpTBin = len(self.pTBinsArrayTruth)-1
        
        print("project from bin {} to bin {}".format(firstBin,lastpTBin))
        mainUnfolded20Proj  = mainUnfolded.ProjectionY("_pyMainUnfolded_Bin20",firstBin,lastpTBin)
        mainUnfolded20S     = self.scaleHistogram(mainUnfolded20Proj)
  
        
        
        c4 = ROOT.TCanvas("c4","c: hist",500*2,450*2)
        c4.cd()
        ROOT.TGaxis.SetMaxDigits(3)
        
        # Set pad and histo arrangement
        myPad4 = ROOT.TPad("myPad4", "The pad",0,0,1,1)
        #myPad2.SetLeftMargin(0.21)
        myPad4.SetLeftMargin(0.15)
        myPad4.SetTopMargin(0.06)
        myPad4.SetRightMargin(0.04)
        myPad4.SetBottomMargin(0.15)
        myPad4.SetTicks()
        myPad4.Draw()
        myPad4.cd()

        maximum = max(maxArrayUnfold)
        #print("This is the maximum: {}".format(max))
        myBlankHisto4 = ROOT.TH1F("myBlankHisto4","Blank Histogram",20,0,1)
        myBlankHisto4.GetXaxis().SetNdivisions(505)
        myBlankHisto4.SetXTitle(xTitle)
        myBlankHisto4.GetXaxis().SetTitleSize(0.05)
        myBlankHisto4.GetXaxis().SetRangeUser(0,1)
        myBlankHisto4.GetXaxis().SetNdivisions(405)
        myBlankHisto4.GetYaxis().SetTitleOffset(1.35)
        myBlankHisto4.GetYaxis().SetTitleSize(0.055)
        myBlankHisto4.SetLineColor(0)

        #Absolute yield
        myBlankHisto4.SetYTitle("dN/dz_{T} unfolded")
        if self.dictKey=="Psi2S":
            myBlankHisto4.GetYaxis().SetRangeUser(0,0.5)
        else:
            myBlankHisto4.GetYaxis().SetRangeUser(0,0.7)
        myBlankHisto4.Draw("E")

        #-Legend
        #-about different contributions
        myLegend4 = ROOT.TLegend(0.2,0.7,0.4,0.9)
        myLegend4.SetTextFont(42)
        myLegend4.SetBorderSize(0)
        myLegend4.SetFillStyle(0)
        myLegend4.SetFillColor(0)
        myLegend4.SetMargin(0.25)
        myLegend4.SetTextSize(0.04)
        myLegend4.AddEntry(hArrayUnfold[0]," #it{p}_{T}^{jet}:","")
  
        MarkerScale=1.6
        #Add all unfolded spectra for jet_pT>20 GeV
        mainUnfolded20S.SetMarkerSize(1.3*MarkerScale)
        mainUnfolded20S.SetMarkerStyle(20)
        mainUnfolded20S.SetMarkerColor(colorArrayFinalFigs[0])
        mainUnfolded20S.SetLineStyle(1)
        mainUnfolded20S.SetLineWidth(2)
        mainUnfolded20S.SetLineColor(colorArrayFinalFigs[0])
        mainUnfolded20S.DrawNormalized("same EP")
        
        myLegend4.AddEntry(mainUnfolded20S,"  %s-%s (GeV/%s)" % (self.pTBinsArrayTruth[firstBin-1],self.pTBinsArrayTruth[lastpTBin],"#it{c}"),"LP")
        myLegend4.Draw()
         
        c4.cd()
        FileName="{}UnfoldingEffect_SpecialJPsiComp_{}.png".format(self.outPath,self.figureTag)
        c4.SaveAs(FileName)
        c4.Close()
      
           
        
     #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def saveResult(self, regParam, nIter, tag):
        
        #-Create an outputRoot file to save the histograms
        #-Will later be used in a plotting Macro
        #prompt/nonpromt, regparam, nIter, tag
        FileName="{}Spectra_reg{}.root".format(self.outPath,regParam)
        fOutData  = ROOT.TFile(FileName, "RECREATE")
        #CyieldArray[i].Write()
       
        if self.isPrompt:
            flag="P"
        else:
            flag="NP"
        mainUnfolded   = self.unfoldedArr2D[regParam-1].Clone("TUnfold_Save") #Array position 0 is iteration number 1,
        original       = self.measuredSpectra2D.Clone("TOrig_Save")
      
        for ptBin in range(1,len(self.pTBinsArrayTruth)):
        
            #-Get histos and save them
            mainUnfoldedProj  = mainUnfolded.ProjectionY("_pyMainUnfolded_Bin{}".format(ptBin),ptBin,ptBin)
            originalProj      = original.ProjectionY("_pyOrig_Bin{}".format(ptBin),ptBin,ptBin)
            #-Make sure they are normalized by bin width
            mainUnfoldedS     = self.scaleHistogram(mainUnfoldedProj)
            originalS         = self.scaleHistogram(originalProj)
        
            #-1) Prompt/NonPrompt Raw Results
            originalS.SetName("hOrig{}_{}_{}".format(flag, self.pTBinsArrayTruth[ptBin-1],self.pTBinsArrayTruth[ptBin]))
            originalS.Write()
        
            #-2) Prompt/NonPrompt Reco Results
            mainUnfoldedS.SetName("hReco{}_{}_{}".format(flag, self.pTBinsArrayTruth[ptBin-1],self.pTBinsArrayTruth[ptBin]))
            mainUnfoldedS.Write()
        
        for ztBin in range(1,len(self.zBinsArrayTruth)):
        
            #-Get histos and save them
            mainUnfoldedProj  = mainUnfolded.ProjectionX("_pxMainUnfolded_Bin{}".format(ztBin),ztBin,ztBin)
            originalProj      = original.ProjectionX("_pxOrig_Bin{}".format(ztBin),ztBin,ztBin)
            #-Make sure they are normalized by bin width
            mainUnfoldedSII     = self.scaleHistogram(mainUnfoldedProj)
            originalSII         = self.scaleHistogram(originalProj)
            #-1) Inclusive Raw Results
        
            #-2) Prompt/NonPrompt Raw Results
            originalSII.SetName("hpTOrig{}_{}_{}".format(flag, self.zBinsArrayTruth[ztBin-1],self.zBinsArrayTruth[ztBin]))
            originalSII.Write()
        
            #-5) Prompt/NonPrompt Reco Results
            mainUnfoldedSII.SetName("hpTReco{}_{}_{}".format(flag, self.zBinsArrayTruth[ztBin-1],self.zBinsArrayTruth[ztBin]))
            mainUnfoldedSII.Write()
        
        
       
        
        
        
        fOutData.Close()


    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def StabilityTest2D(self, nIter, plotAll=False):
        
        ROOT.gStyle.SetOptTitle(0)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(0)
 
        main   = self.unfoldedArr2D[nIter-1] #Array position 0 is iteration number 1,
        lower  = self.unfoldedArr2D[nIter-2]
        higher = self.unfoldedArr2D[nIter]
       
        colorArray = [ROOT.kRed-9,ROOT.kGreen-9,ROOT.kGreen-8,ROOT.kBlue-9, ROOT.kGreen, ROOT.kSpring+8,ROOT.kSpring]
        
        leg2 = ROOT.TLegend(0.2,0.7,0.5,0.93,"")
        leg2.SetFillColor(10)
        leg2.SetBorderSize(0)
        leg2.SetFillStyle(0)
        leg2.SetTextSize(0.04)
        #-Draw the result
        #The result is a 2D histogram with x-Axis==Jet pT, y-Axis==zT
        #make a canvas pad for each pT bin that was unfolded.
        npTBins = len(self.pTBinsArrayTruth)-1
        c = ROOT.TCanvas("c","c: pT",800*npTBins,850)
        c.Divide(npTBins,2,0,0)
        for ptBin in range(1,npTBins+1):
            c.cd(ptBin)
            c.cd(ptBin).SetLeftMargin(0.15)
            c.cd(ptBin).SetRightMargin(0.05)
            c.cd(ptBin).SetTopMargin(0.05)
            c.cd(ptBin).SetBottomMargin(0)
            hMainPTBinProjection    = main.ProjectionY("_pymain",ptBin,ptBin)
            hlowerPTBinProjection   = lower.ProjectionY("_pylower",ptBin,ptBin)
            hhigherPTBinProjection  = higher.ProjectionY("_pyhigher",ptBin,ptBin)
        

            hMainPTBinProjection.GetYaxis().SetTitle("counts per bin")
            hMainPTBinProjection.GetYaxis().SetTitleSize(0.06)
            #hMainPTBinProjection.GetXaxis().SetRangeUser(xRangeMin, xRangeMax)
            #hMainPTBinProjection.GetYaxis().SetRangeUser(2e-10,2e-3)
            hMainPTBinProjection.GetYaxis().SetLabelFont(43)
            hMainPTBinProjection.GetYaxis().SetLabelSize(20)
            hMainPTBinProjection.SetNdivisions(505)
            hMainPTBinProjection.SetLineColor(4)
            hMainPTBinProjection.SetLineWidth(3)
            hMainPTBinProjection.SetLineStyle(1)
            hMainPTBinProjection.DrawCopy("hist E")
            hlowerPTBinProjection.SetLineColor(1)
            hlowerPTBinProjection.SetLineWidth(2)
            hlowerPTBinProjection.SetLineStyle(1)
            hlowerPTBinProjection.DrawCopy("hist E same")
            hhigherPTBinProjection.SetLineColor(2)
            hhigherPTBinProjection.SetLineWidth(2)
            hhigherPTBinProjection.SetLineStyle(1)
            hhigherPTBinProjection.DrawCopy("hist E same")
            
            if ptBin==1:
                leg2.AddEntry(hlowerPTBinProjection,  "nIter={}".format(nIter-1), "l")
                leg2.AddEntry(hMainPTBinProjection,   "final nIter={}".format(nIter), "l")
                leg2.AddEntry(hhigherPTBinProjection, "nIter={}".format(nIter+1), "l")
     
            '''
            if plotAll:
                for i in range(1, len(self.unfoldedArrPerBin[ptBin])+1):
                    if i==nIter or i==(nIter-1) or i==(nIter+1): continue
                    print("iteration: {}".format(i))
                    self.unfoldedArr2D[i-1].SetLineColor(colorArray[i-1])
                    self.unfoldedArr2D[i-1].SetMarkerColor(colorArray[i-1])
                    self.unfoldedArr2D[i-1].DrawCopy("hist same")
                    #leg2.AddEntry(self.unfoldedArrPerBin[ptBin][i-1], "nIter={}".format(i), "l")
            '''
            leg2.Draw("same")
  
            c.cd(npTBins+ptBin)
            c.cd(npTBins+ptBin).SetLeftMargin(0.15)
            c.cd(npTBins+ptBin).SetRightMargin(0.05)
            c.cd(npTBins+ptBin).SetTopMargin(0)
            c.cd(npTBins+ptBin).SetBottomMargin(0.35)
     

            hhigherPTBinProjectionC = hhigherPTBinProjection.Clone("{}Clone".format(hhigherPTBinProjection.GetName()))
            hhigherPTBinProjectionC.Divide(hMainPTBinProjection)
            hhigherPTBinProjectionC.SetMarkerStyle(21)
            hhigherPTBinProjectionC.SetMarkerColor(2)

            hlowerPTBinProjectionC = hlowerPTBinProjection.Clone("{}Clone".format(hlowerPTBinProjection.GetName()))
            hlowerPTBinProjectionC.Divide(hMainPTBinProjection)
            hlowerPTBinProjectionC.SetMarkerStyle(21)
            hlowerPTBinProjectionC.SetMarkerColor(1)

            hlowerPTBinProjectionC.GetXaxis().SetTitleSize(30)
            hlowerPTBinProjectionC.GetXaxis().SetTitleFont(43)
            hlowerPTBinProjectionC.GetXaxis().SetTitleOffset(3.) #was 4
            hlowerPTBinProjectionC.GetXaxis().SetLabelFont(43)
            hlowerPTBinProjectionC.GetXaxis().SetLabelSize(20)
            hlowerPTBinProjectionC.GetYaxis().SetRangeUser(0,2)
            hlowerPTBinProjectionC.GetYaxis().SetTitle("Ratio to nIter={}".format(nIter))
            hlowerPTBinProjectionC.GetYaxis().SetTitleSize(20)
            hlowerPTBinProjectionC.GetYaxis().SetTitleFont(43)
            hlowerPTBinProjectionC.GetYaxis().SetTitleOffset(2.2)
            hlowerPTBinProjectionC.GetYaxis().SetLabelFont(43)
            hlowerPTBinProjectionC.GetYaxis().SetLabelSize(20)
            hlowerPTBinProjectionC.GetYaxis().SetNdivisions(505)
            hlowerPTBinProjectionC.GetXaxis().SetTitle("z_{T}") #ELI make this flexible
      
            hlowerPTBinProjectionC.DrawCopy("hist E")
            hhigherPTBinProjectionC.DrawCopy("hist E same")
            
            '''
            if plotAll:
                for i in range(1, len(self.unfoldedArrPerBin[ptBin])+1):
                    if i==nIter or i==(nIter-1) or i==(nIter+1): continue
                    cloneT = self.unfoldedArrPerBin[ptBin][i-1].Clone("{}Clone".format(self.unfoldedArrPerBin[ptBin][i-1].GetName()))
                    cloneT.SetMarkerStyle(21)
                    cloneT.Divide(main)
                    cloneT.DrawCopy("hist E same")
            '''
            lines = self.drawLines()
            for line in lines:
                cLine = line.Clone("c")
                cLine.Draw("same")
     
        c.cd()
        if plotAll:
            FileName="{}IterationTest2D_{}.png".format(self.outPath,self.figureTag)
        else:
            FileName="{}StabilityTest2D_nIter{}_{}.png".format(self.outPath, nIter,self.figureTag)
        #print("Safe figure under: {}".format(FileName))
        c.SaveAs(FileName)
        c.Close()
        
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def StabilityTest(self, ptBin, nIter, plotAll=False):
        
        ROOT.gStyle.SetOptTitle(0)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(0)
 
        main   = self.unfoldedArrPerBin[ptBin][nIter-1] #Array position 0 is iteration number 1,
        lower  = self.unfoldedArrPerBin[ptBin][nIter-2]
        higher = self.unfoldedArrPerBin[ptBin][nIter]
       
        colorArray = [ROOT.kRed-9,ROOT.kGreen-9,ROOT.kGreen-8,ROOT.kBlue-9, ROOT.kGreen, ROOT.kSpring+8,ROOT.kSpring]
        c = ROOT.TCanvas("c","c: pT",800,850)
        c.cd()
        c.cd().SetLeftMargin(0.15)
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0)
        pad1.SetLeftMargin(0.15)
        pad1.SetRightMargin(0.05)
        pad1.SetTopMargin(0.05)
        #pad1.SetLogy()
        pad1.Draw()
        pad1.cd()
     
        main.GetYaxis().SetTitle("counts per bin")
        main.GetYaxis().SetTitleSize(0.06)
        #main.GetXaxis().SetRangeUser(xRangeMin, xRangeMax)
        #main.GetYaxis().SetRangeUser(2e-10,2e-3)
        main.GetYaxis().SetLabelFont(43)
        main.GetYaxis().SetLabelSize(20)
        main.SetNdivisions(505)
        main.SetLineColor(4)
        main.SetLineWidth(3)
        main.SetLineStyle(1)
        main.Draw("hist E")
        lower.SetLineColor(1)
        lower.SetLineWidth(2)
        lower.SetLineStyle(1)
        lower.DrawCopy("hist E same")
        higher.SetLineColor(2)
        higher.SetLineWidth(2)
        higher.SetLineStyle(1)
        higher.DrawCopy("hist E same")
        
        
        leg2 = ROOT.TLegend(0.2,0.7,0.5,0.93,"")
        leg2.SetFillColor(10)
        leg2.SetBorderSize(0)
        leg2.SetFillStyle(0)
        leg2.SetTextSize(0.04)
        leg2.AddEntry(lower,  "nIter={}".format(nIter-1), "l")
        leg2.AddEntry(main,   "final nIter={}".format(nIter), "l")
        leg2.AddEntry(higher, "nIter={}".format(nIter+1), "l")
 
        if plotAll:
            for i in range(1, len(self.unfoldedArrPerBin[ptBin])+1):
                if i==nIter or i==(nIter-1) or i==(nIter+1): continue
                print("iteration: {}".format(i))
                self.unfoldedArrPerBin[ptBin][i-1].SetLineColor(colorArray[i-1])
                self.unfoldedArrPerBin[ptBin][i-1].SetMarkerColor(colorArray[i-1])
                self.unfoldedArrPerBin[ptBin][i-1].DrawCopy("hist same")
                leg2.AddEntry(self.unfoldedArrPerBin[ptBin][i-1], "nIter={}".format(i), "l")

        leg2.Draw("same")
  
        c.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.35)
        pad2.SetLeftMargin(0.15)
        pad2.SetRightMargin(0.05)
        pad2.Draw()
        pad2.cd()

        higherC = higher.Clone("{}Clone".format(higher.GetName()))
        higherC.Divide(main)
        higherC.SetMarkerStyle(21)
        higherC.SetMarkerColor(2)

        lowerC = lower.Clone("{}Clone".format(lower.GetName()))
        lowerC.Divide(main)
        lowerC.SetMarkerStyle(21)
        lowerC.SetMarkerColor(1)

        lowerC.GetXaxis().SetTitleSize(30)
        lowerC.GetXaxis().SetTitleFont(43)
        lowerC.GetXaxis().SetTitleOffset(3.) #was 4
        lowerC.GetXaxis().SetLabelFont(43)
        lowerC.GetXaxis().SetLabelSize(20)
        lowerC.GetYaxis().SetRangeUser(0,2)
        lowerC.GetYaxis().SetTitle("Ratio to nIter={}".format(nIter))
        lowerC.GetYaxis().SetTitleSize(20)
        lowerC.GetYaxis().SetTitleFont(43)
        lowerC.GetYaxis().SetTitleOffset(2.2)
        lowerC.GetYaxis().SetLabelFont(43)
        lowerC.GetYaxis().SetLabelSize(20)
        lowerC.GetYaxis().SetNdivisions(505)
        lowerC.GetXaxis().SetTitle("z_{T}") #ELI make this flexible
  
        lowerC.DrawCopy("hist E")
        higherC.DrawCopy("hist E same")

        if plotAll:
            for i in range(1, len(self.unfoldedArrPerBin[ptBin])+1):
                if i==nIter or i==(nIter-1) or i==(nIter+1): continue
                cloneT = self.unfoldedArrPerBin[ptBin][i-1].Clone("{}Clone".format(self.unfoldedArrPerBin[ptBin][i-1].GetName()))
                cloneT.SetMarkerStyle(21)
                cloneT.Divide(main)
                cloneT.DrawCopy("hist E same")
 
        lines = self.drawLines()
        for line in lines:
            line.Draw("same")
 

        c.cd()
        if plotAll:
            FileName="{}IterationTestBin{}_{}.png".format(self.outPath,ptBin,self.figureTag)
        else:
            FileName="{}StabilityTestBin{}_nIter{}_{}.png".format(self.outPath,ptBin, nIter,self.figureTag)
        #print("Safe figure under: {}".format(FileName))
        c.SaveAs(FileName)
        c.Close()
    
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def TestRegParam2D(self):
    
        colorArray = [ROOT.kRed-9,ROOT.kGreen-9,ROOT.kGreen-8,ROOT.kBlue-9, ROOT.kGreen, ROOT.kSpring+8,ROOT.kSpring,ROOT.kRed,ROOT.kRed+2,ROOT.kRed+5]

        #-Declare empty histograms to store relative error as a function of zT for all iterations
        ##hArrpTBins = []*len(self.pTBinsArrayTruth) #Fill list of iterations per pT bin
        hArrpTBins = [[] for i in range(len(self.pTBinsArrayTruth))]
        dummyHisto = self.unfoldedArr2D[0].ProjectionY("_pymain", 1, 1)
        dummyHisto.Reset("ICESM")
        npTBins = len(self.pTBinsArrayTruth)-1
        for iterStep in range (0,len(self.unfoldedArr2D)):
            for pTBin in range(0,npTBins+1):
                hArrpTBins[pTBin].append(dummyHisto.Clone("hpTBin{}_ier{}".format(pTBin,iterStep)))
  
        #-Evaluate the error change for different steps of iterations
        #-Build histograms with the relative growth of errors for each iter step as function of zT
        for iterStep in range (1,len(self.unfoldedArr2D)+1):
            for pTBin in range (1,len(self.pTBinsArrayTruth)+1):
                proj = self.unfoldedArr2D[iterStep-1].ProjectionY("_pymain", pTBin, pTBin)
                for zTBin in range(1,proj.GetNbinsX()+1):
                    cont = proj.GetBinContent(zTBin)
                    err  = proj.GetBinError(zTBin)
                    if cont!=0:
                        relErr = err/cont
                        #print("pt {}, zT {}, rel {}".format(pTBin, zTBin, relErr))
                        hArrpTBins[pTBin][iterStep-1].SetBinContent(zTBin,relErr)
        
        leg2 = ROOT.TLegend(0.2,0.7,0.5,0.93,"")
        leg2.SetFillColor(10)
        leg2.SetBorderSize(0)
        leg2.SetFillStyle(0)
        leg2.SetTextSize(0.04)
        leg2.AddEntry(hArrpTBins[1][0],  "nIter=1", "l")

        #-Draw the result
        #The result is a 2D histogram with x-Axis==Jet pT, y-Axis==zT
        #make a canvas pad for each pT bin that was unfolded.
        c = ROOT.TCanvas("c","c: pT",800*npTBins,850)
        c = ROOT.TCanvas("c","c: pT",800*npTBins,800)
        #c.Divide(npTBins,2,0,0)
        c.Divide(npTBins)
        for pTBin in range(1,npTBins+1):
            c.cd(pTBin)
            c.cd(pTBin).SetLeftMargin(0.15)
            c.cd(pTBin).SetRightMargin(0.05)
            c.cd(pTBin).SetTopMargin(0.05)
            #c.cd(pTBin).SetBottomMargin(0)
            c.cd(pTBin).SetBottomMargin(0.15)

   
            hArrpTBins[pTBin][0].GetYaxis().SetTitle("relative error of iteration")
            hArrpTBins[pTBin][0].GetYaxis().SetTitleSize(0.06)
            hArrpTBins[pTBin][0].GetYaxis().SetLabelFont(43)
            hArrpTBins[pTBin][0].GetYaxis().SetLabelSize(20)
            hArrpTBins[pTBin][0].SetNdivisions(505)
            hArrpTBins[pTBin][0].SetLineColor(4)
            hArrpTBins[pTBin][0].SetLineWidth(3)
            hArrpTBins[pTBin][0].SetLineStyle(1)
            hArrpTBins[pTBin][0].SetMarkerStyle(21)
            hArrpTBins[pTBin][0].Draw("hist")
      
            for iterStep in range (2,len(self.unfoldedArr2D)+1):
                hArrpTBins[pTBin][iterStep-1].SetLineColor(colorArray[iterStep])
                hArrpTBins[pTBin][iterStep-1].Draw("same hist")
                if pTBin==1:
                    leg2.AddEntry(hArrpTBins[pTBin][iterStep-1],  "nIter={}".format(iterStep), "l")
 
            leg2.Draw("same")

            '''
            c.cd(npTBins+pTBin)
            c.cd(npTBins+pTBin).SetLeftMargin(0.15)
            c.cd(npTBins+pTBin).SetRightMargin(0.05)
            c.cd(npTBins+pTBin).SetTopMargin(0)
            c.cd(npTBins+pTBin).SetBottomMargin(0.35)
 
            #hvariationError.GetYaxis().SetRangeUser(0.8,1.2)
            hvariationError[pTBin].GetYaxis().SetRangeUser(0.5,1.8)
            hvariationError[pTBin].SetLineColor(4)
            hvariationError[pTBin].SetLineWidth(3)
            hvariationError[pTBin].SetLineStyle(1)
            hvariationError[pTBin].SetMarkerStyle(21)
            hvariationError[pTBin].DrawCopy("hist E")
        
            lines = self.drawLines()
            for line in lines:
                cLine = line.Clone("c")
                cLine.Draw("same")
            '''
        c.cd()
        FileName="{}RegParamTest2D_{}.png".format(self.outPath,self.figureTag)
        #print("Safe figure under: {}".format(FileName))
        c.SaveAs(FileName)
        c.Close()
        
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def StatTestRM2D(self,response):
        
        histo  = self.measuredSpectra2D
        #-This function varies the RM within the statistical uncertainty and then unfolds the measured spectra with this variation
        #-This checks how much the statistical errors in the RM contribute
        #-to uncertainties in the unfolded result
        #-Ergo if you should stuff more events into your RM
        TH2DRM    = response.Hresponse()
        TH2DTruth = response.Htruth()
        TH2DMeas  = response.Hmeasured()

        npTBins         = len(self.pTBinsArrayTruth)-1
        RegulForStatTest= 3
        errorType       = RooUnfold.kCovariance #only for this test
        fRandom         = ROOT.TRandom3(0)
        unfoldTesthistograms = []
        factorStatistic = 1 #Test how the errors would look like if we had X-times the statitic in the RM
        for i in range(1, 1000):  #500
            
            #-Smear the RM and build a new RooUnfoldResponseObject
            RMVariation    = self.SmearPoints(TH2DRM,fRandom,factorStatistic)
            
            #Test to see if this works :-0
            #setup empty object
            RooUnfoldRMVar = ROOT.RooUnfoldResponse("hResponseMatrix3D{}".format(i), "hResponseMatrix3D{}".format(i))
            RooUnfoldRMVar.Setup(TH2DMeas,TH2DTruth,RMVariation)
            #RooUnfoldRMVar.Setup(TH2DMeas,TH2DTruth,TH2DRM)  #Original RM - sanity check

    
            #-Set up the unfolding object
            unfoldVar = RooUnfoldBayes(RooUnfoldRMVar, histo, RegulForStatTest)#
            #-Perform the unfolding
            hunf = unfoldVar.Hreco(errorType) # Produces the truth distribution, with errors, PerBin
            unfoldTesthistograms.append(hunf)
        
    
        #-Create empy histograms for reach pT bin to store mean and std of all histograms
        hvariationResult = []
        hvariationError  = []
        hvariationResultT = []
        hvariationErrorT  = []
        histo = ROOT.TH1D(unfoldTesthistograms[0].ProjectionY("_pymain", 1, 1))
        histo.Reset("ICESM")
        for ptBin in range(0,npTBins+1):
            hvariationResult.append(histo.Clone("pTBinCont_{}".format(ptBin)))
            hvariationError.append(histo.Clone("pTBinErr_{}".format(ptBin)))
            hvariationResultT.append(histo.Clone("pTBinContT_{}".format(ptBin)))
            hvariationErrorT.append(histo.Clone("pTBinErrT_{}".format(ptBin)))
        
        nzTBins     = unfoldTesthistograms[0].ProjectionY("_pymain", 1, 1).GetNbinsX()
        nVariations = len(unfoldTesthistograms)
        
        '''
        #TEST
        c1 = ROOT.TCanvas("c1","c1: pT",800*npTBins,850)
        c1.Divide(npTBins,2,0,0)
        for pTBin in range(1,npTBins+1):
            c1.cd(pTBin)
            c1.cd(pTBin).SetLeftMargin(0.15)
            c1.cd(pTBin).SetRightMargin(0.05)
            c1.cd(pTBin).SetTopMargin(0.05)
            c1.cd(pTBin).SetBottomMargin(0)
            #print("pTBin: {}".format(pTBin))

            for variation in range (0,len(unfoldTesthistograms)):
                
                proj = unfoldTesthistograms[variation].ProjectionY("_pymain", pTBin, pTBin)
                #print("content: {}".format(proj.GetBinContent(3)))
                if variation==0:
                    proj.DrawCopy("hist E")
                else:
                    proj.DrawCopy("same hist E"
        c1.cd()
        FileName="{}RMStatTest2DTTT_nIter{}.png".format(self.outPath, RegulForStatTest)
        #print("Safe figure under: {}".format(FileName))
        c1.SaveAs(FileName)
        c1.Close()
        #TEST
        '''
        
        Version=2   #determine mean and std in each bin from the plain values
        #Version=2   #determine the width of a gaussian distribution of the values
        #if Version==1:
        #-loop over all unfolded histograms and project to the desired pT bin
        for zTBin in range(1,nzTBins+1):
            for pTBin in range(1,npTBins+1):
                contentList=[]
                for variation in range (0,len(unfoldTesthistograms)):
                    proj = unfoldTesthistograms[variation].ProjectionY("_pymain", pTBin, pTBin)
                    contentList.append(proj.GetBinContent(zTBin))
                    #print("content: {}".format(proj.GetBinContent(1)))
                #-when all variations are collected
                contentArr = numpy.array(contentList)
                std  = numpy.std(contentArr, axis = None)
                mean = numpy.median(contentArr, axis = None) #a bit better than mean
                #zTBin=1
                #print("pTBin: {}, zTBin: {}, mean: {}, error: {}".format(pTBin,zTBin,mean,std))
                #print("currently in zTbin: {}".format(zTBin))
                hvariationResult[pTBin-1].SetBinContent(zTBin,mean)
                hvariationResult[pTBin-1].SetBinError(zTBin,std)#1sigma or 3sigma
                hvariationError[pTBin-1].SetBinContent(zTBin,1)
                if mean>0:
                    hvariationError[pTBin-1].SetBinError(zTBin,std/mean)

        #-Use the information from the first round to creat appropriate
        #-histogram ranges
        hvariationResult, hvariationError  = self.calculateDispersion(unfoldTesthistograms, hvariationResult, hvariationError, nzTBins, npTBins)
                
        #-Draw the result
        #The result is a 2D histogram with x-Axis==Jet pT, y-Axis==zT
        #make a canvas pad for each pT bin that was unfolded.
        c = ROOT.TCanvas("c","c: pT",800*npTBins,850)
        c.Divide(npTBins,2,0,0)
        for pTBin in range(1,npTBins+1):
            c.cd(pTBin)
            c.cd(pTBin).SetLeftMargin(0.15)
            c.cd(pTBin).SetRightMargin(0.05)
            c.cd(pTBin).SetTopMargin(0.05)
            c.cd(pTBin).SetBottomMargin(0)
                
            hvariationResult[pTBin-1].GetYaxis().SetTitle("counts per bin")
            hvariationResult[pTBin-1].GetYaxis().SetTitleSize(0.06)
            hvariationResult[pTBin-1].GetYaxis().SetLabelFont(43)
            hvariationResult[pTBin-1].GetYaxis().SetLabelSize(20)
            hvariationResult[pTBin-1].SetNdivisions(505)
            hvariationResult[pTBin-1].SetLineColor(4)
            hvariationResult[pTBin-1].SetLineWidth(3)
            hvariationResult[pTBin-1].SetLineStyle(1)
            hvariationResult[pTBin-1].SetMarkerStyle(21)
            hvariationResult[pTBin-1].Draw("hist E")
      
            c.cd(npTBins+pTBin)
            c.cd(npTBins+pTBin).SetLeftMargin(0.15)
            c.cd(npTBins+pTBin).SetRightMargin(0.05)
            c.cd(npTBins+pTBin).SetTopMargin(0)
            c.cd(npTBins+pTBin).SetBottomMargin(0.35)
 
            #hvariationError.GetYaxis().SetRangeUser(0.8,1.2)
            hvariationError[pTBin-1].GetYaxis().SetRangeUser(0.5,1.8)
            hvariationError[pTBin-1].SetLineColor(4)
            hvariationError[pTBin-1].SetLineWidth(3)
            hvariationError[pTBin-1].SetLineStyle(1)
            hvariationError[pTBin-1].SetMarkerStyle(21)
            hvariationError[pTBin-1].DrawCopy("hist E")
        
            lines = self.drawLines()
            for line in lines:
                cLine = line.Clone("c")
                cLine.Draw("same")
                
        c.cd()
        FileName="{}RMStatTest2D_nIter{}_{}.png".format(self.outPath, RegulForStatTest,self.figureTag)
        #print("Safe figure under: {}".format(FileName))
        c.SaveAs(FileName)
        c.Close()
  
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    def calculateDispersion(self, histogramArray2D, hvariationResult, hvariationError, nzTBins, npTBins):
        
        
        for pTBin in range(1,npTBins+1):
            #-for the first pT bin create several histograms filled with the results of
            #-the unfolding
            hArray1D_v1 = [ROOT.TH1D]*nzTBins
            hArray1D_v2 = [ROOT.TH1D]*nzTBins
         
            for i_zTBin in range(1,nzTBins+1):
                #-Iteration 1
                for i_var in range (0,len(histogramArray2D)):
                    proj = histogramArray2D[i_var].ProjectionY("_pymain", pTBin, pTBin)
                    entry = proj.GetBinContent(i_zTBin)
                    if i_var==0:
                        mean  = hvariationResult[pTBin-1].GetBinContent(i_zTBin)
                        error = hvariationResult[pTBin-1].GetBinError(i_zTBin)
                        #print("ROUND 1: zTBin = {}, mean: {}, error: {}".format(i_zTBin,mean,error))
                        hArray1D_v1[i_zTBin-1] = ROOT.TH1D("Var_{}_zTBin{}".format(i_var,i_zTBin),"",200,mean-error,mean+error)
                    hArray1D_v1[i_zTBin-1].Fill(entry)
                #-Iteration 2
                for i_var in range (0,len(histogramArray2D)):
                    proj = histogramArray2D[i_var].ProjectionY("_pymain", pTBin, pTBin)
                    entry = proj.GetBinContent(i_zTBin)
                    if i_var==0:
                        meanBin  = hArray1D_v1[i_zTBin-1].GetMaximumBin()
                        mean2    = hArray1D_v1[i_zTBin-1].GetBinCenter(meanBin)
                        error    = hArray1D_v1[i_zTBin-1].GetRMS()
                        #print("ROUND 2: zTBin = {}, mean: {}, error: {}".format(i_zTBin,mean2,error))
                        hArray1D_v2[i_zTBin-1] = ROOT.TH1D("Var_{}_zTBin{}".format(i_var,i_zTBin),"",100,mean2-4*error,mean2+4*error)
                    hArray1D_v2[i_zTBin-1].Fill(entry)
                #-Fitting step. Define gauss function for the histogram
                gaussFunc = ROOT.TF1("fb","gaus(0)",hArray1D_v2[i_zTBin-1].GetBinCenter(1),hArray1D_v2[i_zTBin-1].GetBinCenter(hArray1D_v2[i_zTBin-1].GetNbinsX()));
                hArray1D_v2[i_zTBin-1].Fit(gaussFunc)
                #print("ROUND 3: zTBin = {}, mean: {}, error: {}".format(i_zTBin,gaussFunc.GetParameter(1),gaussFunc.GetParameter(2)))
                #-Set the new Values
                hvariationResult[pTBin-1].SetBinContent(i_zTBin, gaussFunc.GetParameter(1))
                hvariationResult[pTBin-1].SetBinError(i_zTBin, gaussFunc.GetParameter(2))
                hvariationError[pTBin-1].SetBinContent(i_zTBin, 1)
                if gaussFunc.GetParameter(1)>0:
                    hvariationError[pTBin-1].SetBinError(i_zTBin, gaussFunc.GetParameter(2)/gaussFunc.GetParameter(1))
            #canvasArray = ROOT.TCanvas("c".format(),"c: pT",800,1600)
            cA = ROOT.TCanvas("c_{}".format(pTBin),"c: pT",1600,1600)
            #canvasArray.Divide(int(0.5*nzTBins),2) #Divide the canvas into pads. One pad per zT bin
            cA.Divide(3,4) #Divide the canvas into pads. One pad per zT bin
            for i_zTBin in range(1,nzTBins+1):
                cA.cd(i_zTBin)
                hArray1D_v2[i_zTBin-1].Draw()
            
            #-
            cA.cd()
            FileName="{}RMStatTest2D_TestAttempt_pT{}_{}.png".format(self.outPath,pTBin,self.figureTag)
            #print("Safe figure under: {}".format(FileName))
            cA.SaveAs(FileName)
            cA.Close()
            
            
        return hvariationResult, hvariationError

   
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def StatTestRM(self,response, ptBin):
        
        histo  = self.measuredSpectraArray[ptBin]
        #-This function varies the RM within the statistical uncertainty and unfold
        #-This checks how much the statistical errors in the RM contribute
        #-to uncertainties in the unfolded result
        #-Ergo if you should stuff more events into your RM
        TH2DRM = response.Hresponse()

        RegulForStatTest= 3
        errorType       = RooUnfold.kCovariance #only for this test
        fRandom         = ROOT.TRandom3(0)
        unfoldTesthistograms = []
        for i in range(1, 1000):  #500
            
            #-Smear the RM and build a new RooUnfoldResponseObject
            RMVariation    = self.SmearPoints(TH2DRM,fRandom)
            hGenLevelMatchedUncutPerBinVar = RMVariation.ProjectionY()
            hGenLevelMatchedUncutPerBinVar.SetName("hGenLevelMatchedUncutPerBin_Var{}".format(i))
            RooUnfoldRMVar = ROOT.RooUnfoldResponse(0, hGenLevelMatchedUncutPerBinVar, RMVariation, "hResponseMatrixMain{}_V{}".format(self.unfoldLabel,i), "hResponseMatrixMain{}_V{}".format(self.unfoldLabel,i))
            RooUnfoldRMVar.UseOverflow(False)
               
            #-Set up the unfolding object
            unfoldVar = RooUnfoldBayes(RooUnfoldRMVar, histo, RegulForStatTest)#
            #-Perform the unfolding
            hunf = unfoldVar.Hreco(errorType) # Produces the truth distribution, with errors, PerBin
            unfoldTesthistograms.append(hunf)
        
        
        Version=1   #determine mean and std in each bin from the plain values
        #Version=2   #determine the width of a gaussian distribution of the values
        if Version==1:
            #-Create empy histogram to store mean and std of all histograms
            hvariationResult = ROOT.TH1D(unfoldTesthistograms[0])
            hvariationError  = ROOT.TH1D(unfoldTesthistograms[0])
            #-Calculate the standard deviation
            for bin in range(1,unfoldTesthistograms[0].GetNbinsX()+1):
                contentList=[]
                for j in range (0,len(unfoldTesthistograms)):
                    contentList.append(unfoldTesthistograms[j].GetBinContent(bin))
                #-In this version just calculate the mean and std from the datapoints
                contentArr = numpy.array(contentList)
                std  = numpy.std(contentArr, axis = None)
                mean = numpy.mean(contentArr, axis = None)
                hvariationResult.SetBinContent(bin,mean)
                #hvariationResult.SetBinError(bin,std)
                hvariationResult.SetBinError(bin,std)#3sigma
                hvariationError.SetBinContent(bin,1)
                hvariationError.SetBinError(bin,std/mean)

        
        #-Draw the solutions
        c = ROOT.TCanvas("c","c: pT",800,850)
        c.cd()
        c.cd().SetLeftMargin(0.15)
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0)
        pad1.SetLeftMargin(0.15)
        pad1.SetRightMargin(0.05)
        pad1.SetTopMargin(0.05)
        #pad1.SetLogy()
        pad1.Draw()
        pad1.cd()
     
        hvariationResult.GetYaxis().SetTitle("counts per bin")
        hvariationResult.GetYaxis().SetTitleSize(0.06)
        hvariationResult.GetYaxis().SetLabelFont(43)
        hvariationResult.GetYaxis().SetLabelSize(20)
        hvariationResult.SetNdivisions(505)
        hvariationResult.SetLineColor(4)
        hvariationResult.SetLineWidth(3)
        hvariationResult.SetLineStyle(1)
        hvariationResult.SetMarkerStyle(21)
        hvariationResult.Draw("hist E")
  
        c.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.35)
        pad2.SetLeftMargin(0.15)
        pad2.SetRightMargin(0.05)
        pad2.Draw()
        pad2.cd()

        #hvariationError.GetYaxis().SetRangeUser(0.8,1.2)
        #hvariationError.GetYaxis().SetRangeUser(0.7,1.3)
        hvariationError.SetLineColor(4)
        hvariationError.SetLineWidth(3)
        hvariationError.SetLineStyle(1)
        hvariationError.SetMarkerStyle(21)
        hvariationError.DrawCopy("hist E")
        c.cd()
        FileName="{}RMStatTestBin{}_nIter{}_{}.png".format(self.outPath, ptBin, RegulForStatTest,self.figureTag)
        #print("Safe figure under: {}".format(FileName))
        c.SaveAs(FileName)
        c.Close()
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def SmearPoints(self, hist,fRandom, factor=1):

        hnew = hist.Clone("hnew")
        
        #-The factor is an experimental factor.
        #-What would happen if we had factor times more statistic
        #-The new error with factor would be = srt(factor)*oldError
        #-The bin counts would be factor*counts
        #-So if we scale down to the old counts to have a comparativley RM
        #-We have new bin content = content*factor/factor
        #-        new error       = error sqrt(factor)/factor

        for binx in range(1, hnew.GetNbinsX()):
            for biny in range(1, hnew.GetNbinsY()):
                cont    = hist.GetBinContent(binx,biny)
                err     = hist.GetBinError(binx,biny)*math.sqrt(factor)/factor
                NewContent = fRandom.Gaus(cont, err)  #1Sigma variation
                #NewContent = fRandom.Poisson(cont)   #in the MC RM the errors should be purely poisson like, but doesn't make a big difference
                if NewContent>0:
                    hnew.SetBinContent(binx,biny,NewContent)
                    hnew.SetBinError(binx,biny,math.sqrt(NewContent))
                else:
                    hnew.SetBinContent(binx,biny,0)
                    hnew.SetBinError(binx,biny,0)

        return hnew
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def PrepareResponseMatrix2D(self,part=0):
    
        #-Get the RM from the input file - binned as /0.05
        RM = self.getResponseMatrix(part,"zTDet","zTPart")
        self.plotHist(RM,"{}".format(RM.GetName()),"colz")

        #print(self.RM)
        #print(type(self.RM))

        #-Get truth-level spectrum (matched) from response matrix projection, before cutting ranges, and rebinning
        #-Needed for kinematic efficiency correction
        hGenLevelMatchedUncutPerBin = RM.ProjectionY() #ELI should this include under/overflow bins??? according to old script alice this is a spectrum with! u/o bins and it has no limits in detlvl jet pt for the projection
        hGenLevelMatchedUncutPerBin.SetName("hGenLevelMatchedUncutPerBin_part{}".format(part))
        hDetLevelMatchedUncutPerBin = RM.ProjectionX()
        hDetLevelMatchedUncutPerBin.SetName("hDetLevelMatchedUncutPerBin{}".format(part))
        self.plotHist(hGenLevelMatchedUncutPerBin,"hGenLevelMatchedUncutPerBin{}".format(part),"E")
        self.plotHist(hDetLevelMatchedUncutPerBin,"hDetLevelMatchedUncutPerBin{}".format(part),"E")

        #-Finally Plot the extracted response matrix
        #-Plot a couple of QA plots
        #-*JER, *JES shift, *Jet finding Efficiency, *kinematic Efficiency
        self.plotJetResponse()
        
        
        #-Prepare a RooUnfold RM object
        #-The first entry is empty since we do not account for fakes
        #-For the second entry, we pass in the projection before the RM was cut to a specific  range, in order to account for the kinematic efficiency
        RooUnfoldRM = ROOT.RooUnfoldResponse(0, hGenLevelMatchedUncutPerBin.Clone("T1"), RM.Clone("T"), "hResponseMatrixMain{}".format(self.unfoldLabel), "hResponseMatrixMain{}".format(self.unfoldLabel))
        RooUnfoldRM.UseOverflow(False)
        

        return RooUnfoldRM
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def PrepareResponseMatrix3D(self, part, weightHistogram):
    
        #-Empty constructor
        RooUnfoldRM = ROOT.RooUnfoldResponse("hResponseMatrix3DMain{}_part{}".format(self.unfoldLabel,part), "hResponseMatrix3DMain{}".format(self.unfoldLabel))

        RM2DTrue = self.getResponseMatrix(part,"jetPtPart","zTPart")
        RM2DDet  = self.getResponseMatrix(part,"jetPtDet","zTDet")
        #-Provide 2DMeasured and 2DTrue, only used for binning
        RooUnfoldRM.Setup(RM2DDet,RM2DTrue) #Is overwritten in Fill function
 
        #self.plotHist(RM2DTrue,"h2DRMtrue","colz")
        #self.plotHist(RM2DDet,"h2DRMdet","colz")

        #-Fill the RM by iterating though the DataFrame
        if "X3872" in self.dictKey:
            MassCut    = (3.840,3.920)
        elif "Psi" in self.dictKey:
            MassCut    = (3.640,3.720)
        else:
            MassCut    = (3.5,4)
            
        MassCut    = (3.5,4)

        jetPtMin=5

        R = ROOT.TRandom()
        #-Get the Ttree from the file
        tTree = self.fResponse.Get("Response")
        for iEvt in range(0, tTree.GetEntries()+1):
            weight  =1
            b_weight=1
            tTree.GetEntry(iEvt)

            rndm = R.Uniform(1)
            zTDet      = tTree.zTDet
            zTPart     = tTree.zTPart
            jetPtDet   = tTree.jetPtDet
            jetPtPart  = tTree.jetPtPart
            tagPtDet   = tTree.tagPtDet
            tagPtPart  = tTree.tagPtPart
            nConstDet  = tTree.nConstDet
            nConstPart = tTree.nConstPart
            etaDet     = tTree.etaDet
            if self.isPrompt==False:
                #-use the event weight from the b-decay merging so that the
                #-ratio of B+ to B0 and hyperon to non-hyperon events is correct!
                b_weight = tTree.weight
                #print("b-weight: {}".format(b_weight))
            if jetPtDet<5: continue    #EEEEELLLIII do we need this cut?
            #.....Add filter here
            #if (self.applyRMCut==False) or (self.applyRMCut==True and nConstDet>1 and etaDet>2.5 and  etaDet<4): #just for testing
            if (self.applyRMCut==False) or (self.applyRMCut==True and nConstDet>1 and nConstPart>1 and etaDet>2.5 and  etaDet<4):
                
                if weightHistogram:
                    evt_weight = self.getEventWeight(weightHistogram, zTPart, jetPtPart)
                    weight     = b_weight*evt_weight
                else:
                    weight = b_weight

                #- - - - - - - - - - -
                if part==0:
                    if self.useTagPt:
                        RooUnfoldRM.Fill(tagPtDet, zTDet, tagPtPart, zTPart, weight)
                    else:
                        RooUnfoldRM.Fill(jetPtDet, zTDet, jetPtPart, zTPart, weight)
                #-Fill RM with half of the statistic (for testing purposes)
                elif part==1 and rndm>0.5:
                    if self.useTagPt:
                        RooUnfoldRM.Fill(tagPtDet, zTDet, tagPtPart, zTPart, weight)
                    else:
                        RooUnfoldRM.Fill(jetPtDet, zTDet, jetPtPart, zTPart, weight)
                #-Fill RM with other half of the statistic (for testing purposes)
                elif part==2 and rndm<0.5:
                    if self.useTagPt:
                        RooUnfoldRM.Fill(tagPtDet, zTDet, tagPtPart, zTPart, weight)
                    else:
                        RooUnfoldRM.Fill(jetPtDet, zTDet, jetPtPart, zTPart, weight)


        RooUnfoldRM.UseOverflow(False)
        
        #for Cross check get truth, and measured and RM from RooUnfoldResponse object
        self.plotHist(RooUnfoldRM.Htruth(),"h2DRMtrue_RM_{}".format(part),"colz","#it{p}_{T,jet} (GeV/#it{c})","#it{z}_{T}")
        self.plotHist(RooUnfoldRM.Hmeasured(),"h2DRMdet_RM_{}".format(part),"colz","#it{p}_{T,jet} (GeV/#it{c})","#it{z}_{T}")
        self.plotHist(RooUnfoldRM.Hresponse(),"h2DRM_RM_{}".format(part),"colz", " "," ")
        
        return RooUnfoldRM
     
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getKinEfficiency(self, bin):
        
        colorArrayFinalFigs = [ROOT.kAzure, ROOT.kAzure-4, ROOT.kCyan-6,  ROOT.kGreen-3, ROOT.kTeal-6,ROOT.kGreen+4,1,1]
        kinEffArr=[]
        for bin in range(0,len(self.pTBinsArrayTruth)-1):
            Full = self.getResponseMatrix(0,"zTDet","zTPart", True, bin, False)  #no cuts
            Cut  = self.getResponseMatrix(0,"zTDet","zTPart", True, bin, True)   #cuts
 
            pojFull = Full.ProjectionY("_pyFull_{}".format(bin),1,Full.GetNbinsX())
            pojCut  = Cut.ProjectionY("_pyCut_{}".format(bin),1,Cut.GetNbinsX())
            pojCut.SetName("ProjectionBin_{}".format(bin))
            pojCut.Divide(pojFull)
            kinEffArr.append(pojCut)
    
        if self.dictKey=="Psi2S":
            xTitle="p_{T}^{#psi(2S)}/p_{T}^{jet}"
        else:
            xTitle="p_{T}^{#chi_{c1}(3872)}/p_{T}^{jet}"
  
   
        #- . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . -
        #-Plot summary after unfolding
        c3 = ROOT.TCanvas("c2","c: hist",500*2,450*2)
        c3.cd()
        ROOT.TGaxis.SetMaxDigits(3)
        
        # Set pad and histo arrangement
        myPad3 = ROOT.TPad("myPad", "The pad",0,0,1,1)
        #myPad2.SetLeftMargin(0.21)
        myPad3.SetLeftMargin(0.15)
        myPad3.SetTopMargin(0.06)
        myPad3.SetRightMargin(0.04)
        myPad3.SetBottomMargin(0.15)
        myPad3.SetTicks()
        myPad3.Draw()
        myPad3.cd()

        #print("This is the maximum: {}".format(max))
        myBlankHisto3 = ROOT.TH1F("myBlankHisto3","Blank Histogram",20,0,1)
        myBlankHisto3.GetXaxis().SetNdivisions(505)
        myBlankHisto3.SetXTitle(xTitle)
        myBlankHisto3.GetXaxis().SetTitleSize(0.05)
        myBlankHisto3.GetXaxis().SetRangeUser(0,1)
        myBlankHisto3.GetXaxis().SetNdivisions(405)
        myBlankHisto3.GetYaxis().SetTitleOffset(1.35)
        myBlankHisto3.GetYaxis().SetTitleSize(0.055)
        myBlankHisto3.SetLineColor(0)
        myBlankHisto3.SetYTitle("kinematic efficiency")
        myBlankHisto3.GetYaxis().SetRangeUser(0.7,1.3)
        myBlankHisto3.Draw("P")

        #-Legend
        myLegend3 = ROOT.TLegend(0.2,0.6,0.4,0.9)
        myLegend3.SetTextFont(42)
        myLegend3.SetBorderSize(0)
        myLegend3.SetFillStyle(0)
        myLegend3.SetFillColor(0)
        myLegend3.SetMargin(0.25)
        myLegend3.SetTextSize(0.04)
        myLegend3.AddEntry(kinEffArr[0]," #it{p}_{T}^{jet}:","")
  
        MarkerScale=1.6
        for i in range(0,len(kinEffArr)):
            kinEffArr[i].SetMarkerSize(1.3*MarkerScale)
            kinEffArr[i].SetMarkerStyle(20)
            kinEffArr[i].SetMarkerColor(colorArrayFinalFigs[i])
            kinEffArr[i].SetLineStyle(1)
            kinEffArr[i].SetLineWidth(2)
            kinEffArr[i].SetLineColor(colorArrayFinalFigs[i])
            kinEffArr[i].DrawCopy("hist P same")
            
            myLegend3.AddEntry(kinEffArr[i],"  %s-%s (GeV/%s)" % (self.pTBinsArrayTruth[i],self.pTBinsArrayTruth[i+1],"#it{c}"),"LP")

        myLegend3.Draw()
         
        c3.cd()
        FileName="{}kinematicEff_{}.png".format(self.outPath,self.figureTag)
        c3.SaveAs(FileName)
        c3.Close()
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getResponseMatrix(self, part, xAxisVar, yAxisVar, fineBinning=False,bin=0,isCut=False, externalFile=None):
    
        if "X3872" in self.dictKey:
            MassCut    = (3.840,3.920)
        elif "Psi" in self.dictKey:
            MassCut    = (3.640,3.720)
        else:
            MassCut    = (3.5,4)
            
        MassCut    = (3.5,4)

        #-Cuts for the Response creation
        resonanceCutString  = "tagMPart > {} && tagMPart < {} ".format(MassCut[0],MassCut[1])
        jetCuts             = "nConstPart>1 && etaDet>2.5 && etaDet<4"
        #jetCuts             = "etaDet>2.5 && etaDet<4"
 
        #-Get the DataFrame from the file
        if externalFile:
            DF = ROOT.RDataFrame("Response",externalFile)
        else:
            DF = ROOT.RDataFrame("Response",self.fResponse)
        
        #-Filter the RM for certain conditions
        #filterString = "{} && {} && {}".format(association,resonanceCutString,jetbinCutString)
        #filterString = "{} && {}".format(resonanceCutString,jetbinCutString)
        filterString = "{} && {}".format(resonanceCutString, jetCuts)
        if not self.isPrompt and externalFile:
            print("Ext filenname:{}".format(externalFile.GetName()))
            #for the external file in case self.isPrompt=False the prompt external file is provided and
            #then is Primary needs to be selected
            filterString = "{} && {} && isPrimaryPart==1".format(resonanceCutString, jetCuts)

        DFfiltered = DF.Filter(filterString)
            
        if fineBinning:
            filterString = "jetPtPart > {} && jetPtPart < {} ".format(self.pTBinsArrayTruth[bin], self.pTBinsArrayTruth[bin+1])
            if isCut:
                filterString+="&& jetPtDet>{} &&  jetPtDet<{}".format(self.pTBinsArrayDet[0],self.pTBinsArrayDet[-1])
                filterString+="&& zTDet>{}".format(self.zBinsArrayTruth[0])
            #- - -
            #-This is the NP montecarlo
            if (self.isPrompt and externalFile) or (not self.isPrompt and externalFile==None):
                RM = DFfiltered.Filter(filterString).Histo2D(("FineBinned_Bin{}".format(bin),"FineBinned",20,0,1,20,0,1), xAxisVar,yAxisVar,"weight")
            else:
                RM = DFfiltered.Filter(filterString).Histo2D(("FineBinned_Bin{}".format(bin),"FineBinned",20,0,1,20,0,1), xAxisVar,yAxisVar)
        else:
            nBinsDet = len(self.zBinsArrayDet) - 1
            if self.zBinsArrayTruth==None:
                self.zBinArrayTruth=self.zBinsArrayDet
            yBinArr  = array('d',self.zBinsArrayTruth)
            nBinsGen = len(self.zBinsArrayTruth)-1
            
            #-Binning for part lvl 2D RM matrix
            #-X-axis
            if xAxisVar=="jetPtPart":
                xBinArr    = array('d',self.pTBinsArrayTruth)
                nBinsXaxis = len(self.pTBinsArrayTruth)-1
            if xAxisVar=="jetPtDet":
                xBinArr    = array('d',self.pTBinsArrayDet)
                nBinsXaxis = len(self.pTBinsArrayDet)-1
            elif xAxisVar=="zTDet":
                xBinArr  = array('d',self.zBinsArrayDet)
                nBinsXaxis = len(self.zBinsArrayDet)-1
            #-Y-axis
            if yAxisVar=="zTPart":
                yBinArr  = array('d',self.zBinsArrayTruth)
                nBinsYaxis = len(self.zBinsArrayTruth)-1
            elif yAxisVar=="zTDet":
                yBinArr  = array('d',self.zBinsArrayDet)
                nBinsYaxis = len(self.zBinsArrayDet)-1

            #-Build the response from the DataFrame
            RMName = "RM_{}_{}".format(xAxisVar,yAxisVar)
            if part==0:
                #-Full statistic for Response Matrix
                #-This is the NP montecarlo
                if (self.isPrompt and externalFile) or (not self.isPrompt and externalFile==None):
                    RM = DFfiltered.Histo2D(("{}Full".format(RMName),"RMFull",nBinsXaxis,xBinArr,nBinsYaxis,yBinArr),xAxisVar,yAxisVar,"weight")
                else:
                    RM = DFfiltered.Histo2D(("{}Full".format(RMName),"RMFull",nBinsXaxis,xBinArr,nBinsYaxis,yBinArr),xAxisVar,yAxisVar)
            elif part==1:
                #-First 1/2 statistic for Response Matrix
                #-Filter all events with uneven rdfentry
                #-This is the NP montecarlo
                if (self.isPrompt and externalFile) or (not self.isPrompt and externalFile==None):
                    RM = DFfiltered.Filter("rdfentry_ % 2 != 0").Histo2D(("{}Half1".format(RMName),"RMHalf1",nBinsXaxis,xBinArr,nBinsYaxis,yBinArr),xAxisVar,yAxisVar,"weight")
                else:
                    RM = DFfiltered.Filter("rdfentry_ % 2 != 0").Histo2D(("{}Half1".format(RMName),"RMHalf1",nBinsXaxis,xBinArr,nBinsYaxis,yBinArr),xAxisVar,yAxisVar)
            elif part==2:
                #-Second 1/2 statistic for Response Matrix
                #-Filter all events with even rdfentry
                #-This is the NP montecarlo
                if (self.isPrompt and externalFile) or (not self.isPrompt and externalFile==None):
                    RM = DFfiltered.Filter("rdfentry_ % 2 == 0").Histo2D(("{}Half2".format(RMName),"RMHalf2",nBinsXaxis,xBinArr,nBinsYaxis,yBinArr),xAxisVar,yAxisVar,"weight")
                else:
                    RM = DFfiltered.Filter("rdfentry_ % 2 == 0").Histo2D(("{}Half2".format(RMName),"RMHalf2",nBinsXaxis,xBinArr,nBinsYaxis,yBinArr),xAxisVar,yAxisVar)

        #-Set sumW2 for correct error propagation
        RM.Sumw2()

        return RM.Clone("C{}".format(RM.GetName())) #To get a proper type and not only the pointer
     
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getRawSpectra(self, mergeTo2D=False):
        
        histoArray = []
        if self.measuredSpectraArray==None:
            #-Raw spectra should be presented as counts per bin (not normalized per bin width)
            #-Open the data file and get the spectra
            for i in range(0,len(self.pTBinsArrayDet)-1):
                if self.isPrompt:
                    histname = "PromptCorr_{}_{}".format(self.pTBinsArrayDet[i],self.pTBinsArrayDet[i+1])#ELI fix this
                else:
                    histname = "NonPromptCorr_{}_{}".format(self.pTBinsArrayDet[i],self.pTBinsArrayDet[i+1])
                #histname = "inclCorr_{}_{}".format(self.pTBinsArrayDet[i],self.pTBinsArrayDet[i+1])
                histObj = self.fData.Get(histname)
                if type(histObj) is not ROOT.TH1:
                    self.plotHist(histObj,"HmeasuredGraph{}".format(i),"E")
                    #-first undo the scaling by bin width so that the counts are per bin
                    self.deScaleGraph(histObj)
                    #-create array of low-edges for each bin
                    arrErrX  = histObj.GetEX()
                    arrX     = histObj.GetX()
                    arrErrY  = histObj.GetEY()
                    arrY     = histObj.GetY()

                    xArrLowEdge = []
                    for j in range(0,len(arrX)):
                          #print("Add {}".format(float(arrX[j])-float(arrErrX[j])))
                          xArrLowEdge.append(arrX[j]-arrErrX[j])
                    #print("Add {}".format(arrX[-1]+arrErrX[-1]))
                    xArrLowEdge.append(arrX[-1]+arrErrX[-1])
                    
                    #-Create a new histogram with the lower bin edges
                    xArrLowEdgeNP = array('d',xArrLowEdge)
                    histogramTransform = ROOT.TH1D(histObj.GetName(),histObj.GetName(),histObj.GetN(),xArrLowEdgeNP)
                    self.zBinsArrayDet = xArrLowEdgeNP #ELI be careful numbers are not limited to 2digits but have rounding errors in 8th digit or so
                    for j in range(0,len(arrX)):
                        histogramTransform.SetBinContent(j+1,arrY[j])
                        histogramTransform.SetBinError(j+1,arrErrY[j])
                        #print("Add content: {}, bin {}".format(arrY[j],j+1))
         
                    #-Set sumW2 for correct errorpropagation
                    #print("1histogramTransform  sumw2")
                    #histogramTransform.Sumw2()
                    #print("2histogramTransform  sumw2")
                    #histogramTransform.Sumw2()
                    self.plotHist(histogramTransform,"HmeasuredHisto{}".format(i),"E")
                    histoArray.append(histogramTransform)
                else:
                    self.plotHist(histObj,"HmeasuredHisto{}".format(i),"E")
                    histogram.Sumw2()
                    histoArray.append(histogram)
                    print("need to still define zT for this case!")
                    self.zBinsArrayDet   = None
        
    
        if mergeTo2D==False:
            return histoArray
        else:
            #-Create a 2D histogram from all the 1D histograms,  x-Axis=jet pT bins, y-Axis=zT bins
            histo2DMeasurement = ROOT.TH2D("mesured2D", "mesured2D", len(self.pTBinsArrayDet)-1, array('d',self.pTBinsArrayDet), len(self.zBinsArrayDet)-1, self.zBinsArrayDet)
            #histo2DMeasurement = ROOT.TH2D("mesured2D", "mesured2D", len(self.zBinsArrayDet)-1, self.zBinsArrayDet, len(self.zBinsArrayDet)-1, self.zBinsArrayDet) self.pTBinsArrayDet    = array('i',ptRangeArray)
            for binx in range(1,histo2DMeasurement.GetNbinsX()+1):     #-pT bins
                for biny in range(1,histo2DMeasurement.GetNbinsY()+1): #-zT bins
                    content = self.measuredSpectraArray[binx-1].GetBinContent(biny)
                    error   = self.measuredSpectraArray[binx-1].GetBinError(biny)
                    histo2DMeasurement.SetBinContent(binx,biny,content)
                    histo2DMeasurement.SetBinError(binx,biny,error)
            #histo2DMeasurement.Sumw2()
            #print("1histo2DMeasurement  sumw2")
            #histo2DMeasurement.Sumw2()
            #print("2histo2DMeasurement  sumw2")
            return histo2DMeasurement

    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def deScaleGraph(self, graph):
        #This function takes a graph that is normalized by the bin width and undoes the normalization
        arrX     = graph.GetX()
        arrErrX  = graph.GetEX()
        arrY     = graph.GetY()
        arrErrY  = graph.GetEY()
        for i in range(0,len(arrErrX)):
            scale=arrErrX[i]*2
            arrY[i]*=scale     #Scale bin content
            arrErrY[i]*=scale  #Scale bin error
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #-This function takes a histogram and scales the bin content by bin width
    def scaleHistogram(self, histo):
       
        histoScale = histo.Clone("{}Scaled".format(histo.GetName()))
        histoScale.Reset("ICESM")
        for i in range (1, histo.GetNbinsX()+1):
            num = histo.GetBinContent(i)
            err = histo.GetBinError(i)
            den = histo.GetBinWidth(i)
            value = 0
            if den!=0:
                nValue = num/den;
                nError = err/den;
                histoScale.SetBinContent(i,nValue)
                histoScale.SetBinError(i,nError)
        #print("scale histogram")
        #histoScale.Sumw2()
        return histoScale

    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def plotCorrelationCoefficients(self, covarianceMatrix, i, Name):
  
        #outputDir = self.outPath + "CorrelationCoefficients/"
        #if not os.path.exists(outputDir):
        #    os.makedirs(outputDir)

        nBinsX = covarianceMatrix.GetNrows()
        nBinsY = covarianceMatrix.GetNcols()
          
        correlationCoefficientMatrix = ROOT.TH2D("correlationCoefficientMatrix", "correlationCoefficientMatrix", nBinsX, 0, nBinsX, nBinsY, 0, nBinsY)

        for xbin in range(0, nBinsX):
            varianceX = covarianceMatrix(xbin, xbin)
            sigmaX = math.sqrt(varianceX)

            for ybin in range(0, nBinsY):
              varianceY = covarianceMatrix(ybin, ybin)
              sigmaY = math.sqrt(varianceY)

              covXY = covarianceMatrix(xbin, ybin)
              if sigmaX > 0 and sigmaY > 0:
                Cxy = covXY / (sigmaX * sigmaY)
                correlationCoefficientMatrix.SetBinContent(xbin+1, ybin+1, Cxy)

                #print "sigma x: {}, sigmay: {}".format(sigmaX, sigmaY)
                #print "cov (x,y) = {}".format(covXY)
                #print "Cxy = {}".format(Cxy)

        self.plotHist(correlationCoefficientMatrix, "hCorrelationCoefficientMatrix{}".format(i), "colz")

    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def plotJetResponse(self):
        #-Create subdirectory for the jet response matrix QA
        subdirpath = self.makeOutDir("Response")


    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def drawLines(self):
        
        line = ROOT.TLine(0,1,1,1)
        line.SetLineColor(1)
        line.SetLineStyle(1)
        
        line2 = ROOT.TLine(0,1.2,1,1.2)
        line2.SetLineColor(1)
        line2.SetLineStyle(2)
        
        line3 = ROOT.TLine(0,0.8,1,0.8)
        line3.SetLineColor(1)
        line3.SetLineStyle(2)
        
        lineList = [line, line2, line3]
        return lineList
        
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def makeOutDir(self, subDirName):
        
        if not self.outPath.endswith("/"):
            subDirPath = self.outPath + "/{}".format(subDirName)
        else:
            subDirPath = self.outPath + subDirName

        if not os.path.exists(subDirPath): # Create current label dir if necessary
          os.makedirs(subDirPath)
        
        return subDirPath
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def plotHist(self, hist, outputFilename, drawOptions = "", xTitle = "", yTitle = ""):
  
        ROOT.gStyle.SetOptTitle(0)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(0)
        
        c = ROOT.TCanvas("c","c: hist",600,450)
        c.cd()
        c.cd().SetLeftMargin(0.15)

        #if setLogy:
        #  c.SetLogy()
        #if setLogz:
        #  c.SetLogz()
        #obj.InheritsFrom(ROOT.TH1.Class()): ELI
        if type(hist) is ROOT.TH1 or type(hist) is ROOT.TH1D:
            #print("is 1D histo")
            hist.SetNdivisions(505)
            #hist.GetYaxis().SetRangeUser(0,hist.GetMaximum()*1.2)
            #hist.GetYaxis().SetRangeUser(hist.GetMinimum(),hist.GetMaximum()*1.2)
            hist.Draw(drawOptions)
        elif type(hist) is ROOT.TH2 or type(hist) is ROOT.TH2D:
            #print("is 2D histo")
            hist.SetNdivisions(505)
            hist.SetXTitle(xTitle)
            hist.SetYTitle(yTitle)
            #hist.Scale(1.0/hist.Integral())###ELI comment out later
            hist.Draw(drawOptions)
            
            #textFit = ROOT.TLatex()
            #textFit.SetTextSize(0.04)
            #textFit.SetNDC()
            #textFit.DrawLatex(0.6,0.8,"counts: {}".format(hist.Integral()))
            
        else:
            #print("is not a 1D or 2D histo")
            myBlankHisto = hist.GetHistogram()
            #myBlankHisto = ROOT.TH1F("myBlankHisto{}".format(hist.GetName()),"Blank Histogram{}".format(hist.GetName()),10,0,1)
            myBlankHisto.SetNdivisions(505)
            #myBlankHisto.SetXTitle("#it{p}_{T,jet} (GeV/#it{c})")
            #myBlankHisto.GetYaxis().SetTitleOffset(1.8)
            #myBlankHisto.SetYTitle("#it{R}_{AA}")
            #myBlankHisto.GetYaxis().SetRangeUser(0,hist.GetMaximum())
            myBlankHisto.GetXaxis().SetRangeUser(0,1)
            myBlankHisto.GetYaxis().SetRangeUser(0,hist.GetMaximum()*1.2)
            myBlankHisto.Draw("E")
            hist.Draw("same E")

        #if setLogy:
        #    FileName="{}{}.C".format(self.outPath,outputFilename)
        #else:
        FileName="{}{}.png".format(self.outPath,outputFilename)
        #FileName="{}{}.C".format(self.outPath,outputFilename)
        #print("Safe figure under: {}".format(FileName))
        c.SaveAs(FileName)
        c.Close()

###############################################################################
def unfoldzTDistributionLHCb(variation):
        
    maxIter=4
    #maxIter=1

    extFilename_P  = ""
    extFilename_NP = ""
    priorType      = ""
    
 
    tag=""
    #FileRM_P ="AllPromptPsi2S_X3872"
    #FileRM_P ="job425_426_filteredV1"
    FileRM_P ="All_PromptX3872"
    #----
    FileRM_NP="16_17_18AllBdecay_ScaledNom_V2" #All MC  2016+2017+2018

    #- - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - -
    #Variations
    #1. Swap prior
    if variation==1:
        tag="SWP_prior"
        extFilename_P  = FileRM_NP
        extFilename_NP = FileRM_P
        priorType      = ""
        maxIter=1
    #2. Flat prior
    if variation==2:
        tag="Flat_prior"
        extFilename_P  = FileRM_NP
        extFilename_NP = FileRM_P
        priorType      = "Flat"
        maxIter=1
    #3. hyperon fraction in B-decay +20%
    if variation==3:
        tag="moreHyperons"
        FileRM_NP="16_17_18AllBdecay_Scaled_p20pc_V2" #All MC  2016+2017+2018
        maxIter=1
    #4. hyperon fraction in B-decay -20%
    if variation==4:
        tag="lessHyperons"
        FileRM_NP="16_17_18AllBdecay_Scaled_m20pc_V2" #All MC  2016+2017+2018
        maxIter=1
    '''
    #5. variation of time fit Standard is timeFfit (functional fit vs. graph fit)
    if var==5:
        tag="timeHfit"
        DataHistograms=""
    #6. acc Corr variation
    if var==6:
        tag="accCorrVar"
        DataHistograms=""
    #7. Variation of the bDecay??? Isn't this included in the erros that propagate?
    if var==7:
        tag="BfracVar"
        DataHistograms=""
    '''
 
    '''
    #- - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - -
    #Psi(2S)
    #pTBinArray = [5,10,15,20,30,40,100]
    #pTBinArray = [5,10,15,20,30,40,60]
    pTBinArray = [2,5,10,15,20,30,40]  #for Tag pT binning
    
    
    regularizationParam=3
    unfoldObjectPrompt = UnfoldSpectraClass("P",FileRM_P,"CorrectedFinalHistograms_Psi2S","Psi2S",pTBinArray)
    unfoldObjectPrompt.provideExtPrior(tag, extFilename_P,priorType)
    for iter in range(0,maxIter):
        unfoldObjectPrompt.unfold2D(regularizationParam,iter, tag)
    
    
    regularizationParam=3
    unfoldObjectNonPrompt = UnfoldSpectraClass("NP",FileRM_NP,"CorrectedFinalHistograms_Psi2S","Psi2S",pTBinArray)
    unfoldObjectNonPrompt.provideExtPrior(tag, extFilename_NP,priorType)
    for iter in range(0,maxIter):
        unfoldObjectNonPrompt.unfold2D(regularizationParam,iter, tag)

    
    '''
    #- - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - -
    #X3872
    #pTBinArray = [5,10,15,20,50]
    #pTBinArray = [5,10,15,20,30]
    pTBinArray = [2,5,10,15,20,30]
    
    
    regularizationParam=3
    unfoldObjectPrompt = UnfoldSpectraClass("P",FileRM_P,"CorrectedFinalHistograms_X3872","X3872",pTBinArray)
    unfoldObjectPrompt.provideExtPrior(tag, extFilename_P,priorType)
    for iter in range(0,maxIter):
        unfoldObjectPrompt.unfold2D(regularizationParam,iter, tag)
    

    regularizationParam=3
    unfoldObjectNonPrompt = UnfoldSpectraClass("NP",FileRM_NP,"CorrectedFinalHistograms_X3872","X3872",pTBinArray)
    unfoldObjectNonPrompt.provideExtPrior(tag, extFilename_NP,priorType)
    for iter in range(0,maxIter):
        unfoldObjectNonPrompt.unfold2D(regularizationParam,iter, tag)
    

###############################################################################
def unfoldzTAllDistributionLHCb(variation):
 
    unfoldzTDistributionLHCb(0) #0. MAIN
    unfoldzTDistributionLHCb(1) #1. Swap prior
    unfoldzTDistributionLHCb(2) #2. Flat prior
    unfoldzTDistributionLHCb(3) #3. hyperon fraction in B-decay +20%
    unfoldzTDistributionLHCb(4) #4. hyperon fraction in B-decay -20%
 
#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description="- -")
 
  parser.add_argument("-v", "--variationArg", action="store",
                      type=int, metavar="hyperon Fraction",
                      default=0,
                      help="variation=0 ->> main result, variation>0 ->> systematic variations")

  # Parse the arguments
  args = parser.parse_args()

  #unfoldzTDistributionLHCb(variation = args.variationArg)
  unfoldzTAllDistributionLHCb(variation = args.variationArg)
#  on lxplus
#     lb-run DaVinci/v46r0 python unfoldzTDistributionLHCb.py -v 0
#     lb-run DaVinci/v46r0 python unfoldzTDistributionLHCb.py -v 0

