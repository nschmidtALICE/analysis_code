import os
import time
import math
import sys, ROOT

import numpy as np
import array
from collections import OrderedDict
import random
import argparse
from ROOT import gDirectory
import pickle
import ratio

sys.stdout.flush()
ROOT.gROOT.SetBatch(True)

#- - - - - -  - - - - - - -
# how to run: python ntuplesJetAnalysis_Response_noVect.py -f filename -p 10 -m 0
# On lxplus: lb-run DaVinci/v45r6 python ntuplesJetAnalysis_Response_noVect.py -f /eos/lhcb/user/e/efranzis/Output/246/merge246 -p 10 -m 0
#This creates a fitlered dataset and the response matrix for MC
#The data structure is a simple ntuple
#- - - - - - - - - - - - -

#Todo:
#-** Add Matching Criteria to the RM ntuple (1 or 4)
#-** created ResponseVar - Add a type of RM not crated by physical distance but by finding the MC match (only 1)
#-** created separate Tree - Consider the caveat that only MC events that have a tag at Det and Part lvl are accepted
#-     That means the part. lvl spectrum is biased!!!

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
def findClosest(matrix,row,version):
  
  #rows:     PartLvl Jet idx
  #collumns: DetLvl Jet idx
  #content:  distance between jets
  position   = -1
  
  minPartLvl = np.amin(matrix[row])  #Find the value    of smallest element (minimum distance) for this particle level jet
  collumnMin = np.argmin(matrix[row])#Find the position of the minimum in that row (=collumn==DetLvlv-jet)

  if minPartLvl==10: return -1 #default value of matrix elements
  #- - - - - - - - - - - - - -
  #-just find the smallest distance for this PartLvl jet
  #-This means one detLvl jet can be linked to multiple PartLvl jets
  if version==0:
    position = collumnMin
  #- - - - - - - - - - - - - -
  #-Only use the smallest distance, if this PartLvl jet is also
  #-the closest jet for the Det Lvl jet
  #-Uniquley closest distance
  if version==1:
    minDetLvl = np.amin(matrix[:,collumnMin])
    if minDetLvl==10: return -1 #default value of matrix elements

    if minDetLvl==minPartLvl:
      position = collumnMin
   
  return int(position) #collumn=DetLvl jet idx
#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
def getFourVector(tree,type,index):
    
    if index<0:
        print("This particle does not exist - please check and debug, type: {}".format(type))
    
    if type=="mcjet":
    #FindBranch (const char *name)
    #GetBranch (const char *name)
      #px = tree.GetBranch('{}_px'.format(type))[index]
      px = tree.mcjet_px[index]
      py = tree.mcjet_py[index]
      pz = tree.mcjet_pz[index]
      e  = tree.mcjet_e[index]
    if type=="mctag":
      px = tree.mctag_px[index]
      py = tree.mctag_py[index]
      pz = tree.mctag_pz[index]
      e  = tree.mctag_e[index]
    if type=="mcprt":
      px = tree.mcprt_px[index]
      py = tree.mcprt_py[index]
      pz = tree.mcprt_pz[index]
      e  = tree.mcprt_e[index]
    if type=="jet":
      px = tree.jet_px[index]
      py = tree.jet_py[index]
      pz = tree.jet_pz[index]
      e  = tree.jet_e[index]
    if type=="tag":
      px = tree.tag_px[index]
      py = tree.tag_py[index]
      pz = tree.tag_pz[index]
      e  = tree.tag_e[index]
    if type=="prt":
      px = tree.prt_px[index]
      py = tree.prt_py[index]
      pz = tree.prt_pz[index]
      e  = tree.prt_e[index]
      
    fourVect = ROOT.TLorentzVector(px,py,pz,e)
    return fourVect

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
def calcLifeTime(mass,pvrZ,decayZ,pz):
  
  massX      = 3.87160   #PDG 2020
  massPsi2s  = 3.686097  #PDG 2020

  sToPicoSeconds = 1e12           #conversion: s->ps
  speedOfLight   = 2.99792458e8   #in m/s
  mmTom          = 1e-3           #conversion: mm->m
        
  if mass>3.779:
    NominalMass=massX
  else:
    NominalMass=massPsi2s
    
  lifetime = ((decayZ - pvrZ)*mmTom*(NominalMass)/(speedOfLight))/pz
  lifetime*=sToPicoSeconds #Change the liftime from s to ps to avoid these large exponents, better for fitting

  return lifetime

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
def isAccepted(particleVect, Type=0):
    eta = particleVect.PseudoRapidity()
    if Type == 0:  # General Particles (D0, K, pi)
        if eta < 2 or eta > 4.5:
            return False
    elif Type == 1:  # Jets
        if eta < 2.5 or eta > 4:
            return False
    
    # D0-specific pt cut
    if Type == 2:  # D0 specific
        if particleVect.Pt() < 2.0:  # Minimum pT for D0
            return False
            
    # Otherwise True
    return True


########################################################################################################
# Main I                                  ##############################################################
########################################################################################################
class filterObject():
 
    def __init__(self,fFileName,printLvl,inputMC): #
        self.isMC=False
        if inputMC==1:
            self.isMC=True

        #print(isMC)
        #print(type(isMC))
        if self.isMC:
            print("Start filtering for MC")
        else:
            print("Start filtering for Data")
            
        self.printLvl = printLvl
        #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        #-Select the way you want to match your jets
        self.matchingVersion=-1

        #o Version 1
        #-match all particle and detector level jets to each other
        #-that contain a tag
        #-(exclusivley closest pair wins)
        self.matchingVersion=1

        #o Version 2
        #-find the associated Generator Tag to the Detector level tag
        #-and match the two jets which contain them
        #self.matchingVersion=2

        #o Version 3
        #- match all part and det level jets with each other and then afterwards
        #- ask if the matched pair both contain also a tag
        #self.matchingVersion= 3


        #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        #-Load the efficiency maps
        # if not self.isMC:
        #     #. . . . . . . . . . . . . . . . . . . . . .
        #     #- Pions
        #     base="/eos/lhcb/user/e/efranzis/Efficiencies"  #on lxplus
        #     #base="/Users/eliane/LHCb/x3872-code/filterGangaOutput"  #-local
        #     pathPion="{}/NewEffCorrection_pion".format(base)
        #     #-
        #     fInFilePionReco  = ROOT.TFile("{}/eff2DMap_recocombinedPtPion.root".format(pathPion),"READ")
        #     self.pionMapReco = fInFilePionReco.Get("eff_histo")
        #     self.pionMapReco.SetDirectory(0) #Verified

        #     fInFilePionReco_Var  = ROOT.TFile("{}/effHist_recoHalfBins.root".format(pathPion),"READ")
        #     self.pionMapReco_Var = fInFilePionReco_Var.Get("eff_histo")
        #     self.pionMapReco_Var.SetName("pionRecoHalfBin")
        #     self.pionMapReco_Var.SetDirectory(0) #Verified



        #     #fInFilePionSel  = ROOT.TFile("{}/effboth9_pionProbNN05.root".format(pathPion),"READ")
        #     fInFilePionSel  = ROOT.TFile("{}/effboth10_pionProbNN07.root".format(pathPion),"READ")
        #     #fInFilePionSel  = ROOT.TFile("{}/effboth11_pionProbNN09.root".format(pathPion),"READ")
        #     pionMapSel      = fInFilePionSel.Get("total_clone")  #is a TEfficiency object
        #     pionMapSel.SetDirectory(0)
        #     self.pionMapSelHist  = pionMapSel.CreateHistogram()
        #     self.pionMapSelHist.SetName("effHistoPionSel")
        #     self.pionMapSelHist.SetDirectory(0) #Verified

        #     self.pionMapSelHist_Var = self.createMapVariation(self.pionMapSelHist,pionMapSel) #-subtract one sigma of stat error from the map value
        #     self.pionMapSelHist_Var.SetDirectory(0)

        #     #. . . . . . . . . . . . . . . . . . . . . .
        #     #- Muons
        #     pathMuon="{}/NewEffCorrection_muonJPsi".format(base)
        #     fInFileMuon = ROOT.TFile("{}/effHist_reco.root".format(pathMuon),"READ")
        #     self.MuonMap = fInFileMuon.Get("eff_histo")
        #     self.MuonMap.SetDirectory(0)   #Verified

        #     fInFileMuon_Var  = ROOT.TFile("{}/effHist_recoHalfBins.root".format(pathMuon),"READ")
        #     self.MuonMap_Var = fInFileMuon_Var.Get("eff_histo")
        #     self.MuonMap_Var.SetDirectory(0)   #Verified


        #     #. . . . . . . . . . . . . . . . . . . . . .
        #     #- J/Psi / Trigger
        #     fInFileTriggerSel = ROOT.TFile("{}/pidmu0_ismuon_magup.root".format(pathMuon),"READ")
        #     self.TriggerSelCorrMap = fInFileTriggerSel.Get("efficiency_avg")
        #     self.TriggerSelCorrMap.SetDirectory(0) #Verified

        #     self.TriggerSelCorrMap_Var = self.createMapVariation(self.TriggerSelCorrMap) #-subtract one sigma of stat error from the map value
        #     self.TriggerSelCorrMap_Var.SetDirectory(0)
        #     #-
        #     fInFileTriggerEff = ROOT.TFile("{}/effHist_trigger_Data.root".format(pathMuon),"READ")
        #     self.TriggerEffMap2D = fInFileTriggerEff.Get("eff_histo")
        #     self.TriggerEffMap2D.SetDirectory(0)
        #     self.functionList = self.fitTriggerEff(self.TriggerEffMap2D)

        #     fInFileTriggerEff_Var = ROOT.TFile("{}/effHist_triggerHalfBins.root".format(pathMuon),"READ")
        #     self.TriggerEffMap2D_Var = fInFileTriggerEff_Var.Get("eff_histo")
        #     self.TriggerEffMap2D_Var.SetDirectory(0)
        #     self.TriggerEffMap2D_Var.SetName("effHalfBin")
        #     self.functionList_Var = self.fitTriggerEff(self.TriggerEffMap2D_Var)

        #     #. . . . . . . . . . . . . . . . . . . . . .
        #     #- R_Data/MC
        #     ratioPath="{}/R_Data_MC/".format(base)
        #     self.Ratio_D_MC = ratio.Ratio(ratioPath)

        #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        self.fFileName = fFileName
        if self.fFileName.endswith(".root"):
            self.fFileName = fFileName[:-5]
        print("-> Inputfilename: {}".format(self.fFileName))
        self.foutName='{}_filterV{}.root'.format(self.fFileName,self.matchingVersion)

        self.tfileOriginal = ROOT.TFile('{}.root'.format(self.fFileName),'READ')
        self.ttreeOriginal = self.tfileOriginal.data

        self.tfile = ROOT.TFile(self.foutName, 'RECREATE')
        self.ttree = ROOT.TTree('FragmNtuple', '')



    ########################################################################################################
    def fitTriggerEff(self, triggerHist):
        
        #-get number of eta bins == fit functions
        etaBins = triggerHist.GetNbinsY()
        FitFuncArray  = [None]*etaBins
        
        print("start fitting")
        triggerHist.Sumw2()
        canvas0 = ROOT.TCanvas("canvas0","canvas0",0,0,400,400)
        canvas0.cd()
        triggerHist.Draw("colz")
        #canvas0.SaveAs("./Eff2D.pdf")
             
        canvas = ROOT.TCanvas("canvas","canvas",0,0,600,400)
        canvas.Divide(2,5);

        for etaBin in range(0,etaBins):
            canvas.cd(etaBin+1)
            FitFuncArray[etaBin] = ROOT.TF1("func{0}".format(etaBin),"[0] - ([1]*exp(-1.0*[2]*x))",1000.0,11000.0);
            FitFuncArray[etaBin].SetParameter(0,0.6)
            FitFuncArray[etaBin].SetParLimits(0,0.3,0.9)
            FitFuncArray[etaBin].SetParameter(1,0.1)
            FitFuncArray[etaBin].SetParLimits(1,0.001,5.0)
            FitFuncArray[etaBin].SetParameter(2,0.00001)
            FitFuncArray[etaBin].SetParLimits(2,0.000001,0.01);
     
            #-Project every etaBin to the pT axis and then fit
            proj = triggerHist.ProjectionX("effHistProjNum{0}".format(etaBin),etaBin+1,etaBin+1)
            proj.Draw("E")
            #-Fit
            #set max range to 10000
            proj.Fit(FitFuncArray[etaBin],"Q","",0,10000)
            FitFuncArray[etaBin].Draw("same")
        
        #-Save the final Canvas
        #canvas.SaveAs("./EffFits.pdf")
        
        return FitFuncArray
    
    ########################################################################################################
    def createMapVariation(self, inputHisto, TEffObject=None):
        
        mapClone = inputHisto.Clone("{}_Var".format(inputHisto.GetName()))
        
        for binx in range (0,inputHisto.GetNbinsX()):
            for biny in range (0,inputHisto.GetNbinsY()):
                binContent = inputHisto.GetBinContent(binx,biny)
                binError   = inputHisto.GetBinError(binx,biny)
                #-This is a histogram from a TEff object. The errors are nonsensical
                #-and need to be determined in a separate way
                if TEffObject:
                    errEffObj = TEffObject.GetEfficiencyErrorLow(TEffObject.GetGlobalBin(binx,biny))
                    mapClone.SetBinContent(binx,biny,binContent-errEffObj)
                #-Normal histogram
                else:
                    mapClone.SetBinContent(binx,biny,binContent-binError)
        return mapClone

    ########################################################################################################
    def filter(self):
        # D0 tag variables
        v_tagdR = array.array('f', [0])        # Distance between D0 and jet axis
        v_tagMass = array.array('f', [0])      # Reconstructed D0 mass
        v_tagPt = array.array('f', [0])        # D0 transverse momentum
        v_tagEta = array.array('f', [0])       # D0 pseudorapidity
        v_tag_idx_jet = array.array('f', [0])  # Index of associated jet
        v_tag_lifet = array.array('f', [0])    # D0 lifetime
        v_tag_lifetWrongPV = array.array('f', [0])  # D0 lifetime with wrong primary vertex
        v_tagnRnd = array.array('f', [0])      # Number of random tags
        v_tag_decVtxChi2 = array.array('f', [0])  # D0 decay vertex chi2
        v_tag_logdecVtxChi2 = array.array('f', [0])  # Log of decay vertex chi2



        # Jet variables
        v_jetnTComb = array.array('f', [0])    # Number of combined tags in jet
        v_jetnTUnique = array.array('f', [0])  # Number of unique tags in jet
        v_jetnConst = array.array('f', [0])    # Number of jet constituents
        v_jetPt = array.array('f', [0])        # Jet transverse momentum
        v_jetEta = array.array('f', [0])       # Jet pseudorapidity
        v_tagZ = array.array('f', [0])         # Fragmentation fraction (pT_D0/pT_jet)
        v_D0Z = array.array('f', [0])          # Secondary fragmentation fraction
        
        # D0 daughter variables
        v_mpipi = array.array('f', [0])        # Mass of pion pair
        v_QValue = array.array('f', [0])       # Q-value
    
        # Efficiency weight variables
        v_EffWeight = array.array('f', [0])     # Combined efficiency weight
        v_EffWeight_0 = array.array('f', [0])   # Component 0 (pion reconstruction)
        v_EffWeight_1 = array.array('f', [0])   # Component 1 (pion selection)
        v_EffWeight_2 = array.array('f', [0])   # Component 2 (muon reconstruction)
        v_EffWeight_3 = array.array('f', [0])   # Component 3 (trigger line correction)
        v_EffWeight_4 = array.array('f', [0])   # Component 4 (trigger efficiency)
        v_EffWeight_5 = array.array('f', [0])   # Component 5 (data/MC ratio)
        
        # Variation of efficiency weights
        v_EffWeight_0_Var = array.array('f', [0])
        v_EffWeight_1_Var = array.array('f', [0])
        v_EffWeight_2_Var = array.array('f', [0])
        v_EffWeight_3_Var = array.array('f', [0])
        v_EffWeight_4_Var = array.array('f', [0])
        v_EffWeight_5_Var = array.array('f', [0])

        # D0 pion and kaon quality variables
        v_piPTrckChi2 = array.array('f', [0])   # π+ track chi2
        v_piPprobNNpi = array.array('f', [0])   # π+ PID probability
        v_piPprobGhost = array.array('f', [0])  # π+ ghost probability
        v_piMTrckChi2 = array.array('f', [0])   # π- track chi2
        v_piMprobNNpi = array.array('f', [0])   # π- PID probability
        v_piMprobGhost = array.array('f', [0])  # π- ghost probability
        v_decayVtxChi2 = array.array('f', [0])  # Decay vertex chi2
        v_Dist1 = array.array('f', [0])         # Distance of closest approach 1
        v_Dist2 = array.array('f', [0])         # Distance of closest approach 2
        v_Dist3 = array.array('f', [0])         # Distance of closest approach 3
        v_mD0 = array.array('f', [0])           # Reconstructed D0 mass
        v_nSPD = array.array('f', [0])          # Number of SPD hits
        v_isPrimary2 = array.array('f', [0])    # Flag for primary D0
        v_KprobNNK = array.array('f', [0])      # K PID probability
        v_KprobGhost = array.array('f', [0])    # K ghost probability
        v_KTrckChi2 = array.array('f', [0])     # K track chi2


        self.ttree.Branch('tagJetdR', v_tagdR, "tagJetdR/F")
        self.ttree.Branch('tagMass', v_tagMass, "tagMass/F")
        self.ttree.Branch('tagPt', v_tagPt, "tagPt/F")
        self.ttree.Branch('tagEta', v_tagEta, "tagEta/F")
        self.ttree.Branch('tagidxjet', v_tag_idx_jet, "tagidxjet/F")
        self.ttree.Branch('tag_lifet', v_tag_lifet, "tag_lifet/F")
        self.ttree.Branch('tag_lifetWrongPV', v_tag_lifetWrongPV, "tag_lifetWrongPV/F")
        self.ttree.Branch('tagnRnd', v_tagnRnd, "tagnRnd/F")   #n_tagsRnd
        self.ttree.Branch('tag_ip_chi2', v_tag_decVtxChi2, "tag_ip_chi2/F")
        self.ttree.Branch('log_tag_ipchi2', v_tag_logdecVtxChi2, "log_tag_ipchi2/F")

        self.ttree.Branch('jetnTComb', v_jetnTComb, "jetnTComb/F")
        self.ttree.Branch('jetnTUniqu', v_jetnTUnique, "jetnTUniqu/F")
        self.ttree.Branch('jetnConst', v_jetnConst, "jetnConst/F")
        self.ttree.Branch('jetPt', v_jetPt, "jetPt/F")
        self.ttree.Branch('jetEta', v_jetEta, "jetEta/F")
        self.ttree.Branch('tagZ', v_tagZ, "tagZ/F")
        self.ttree.Branch('D0Z', v_D0Z, "D0Z/F")  # Change from JPsiZ

        self.ttree.Branch('mpipi', v_mpipi, "mpipi/F")
        self.ttree.Branch('QValue', v_QValue, "QValue/F")

        #-Analysis cut variables for efficiency correction
        self.ttree.Branch('piPTrckChi2', v_piPTrckChi2, "piPTrckChi2/F")
        self.ttree.Branch('piPprobNNpi', v_piPprobNNpi, "piPprobNN/F")     #--
        self.ttree.Branch('piPprobGhost', v_piPprobGhost, "piPprobGhost/F") #--
        self.ttree.Branch('piMTrckChi2', v_piMTrckChi2, "piMTrckChi2/F")
        self.ttree.Branch('piMprobNNpi', v_piMprobNNpi, "piMprobNN/F")     #--
        self.ttree.Branch('piMprobGhost', v_piMprobGhost, "piMprobGhost/F") #--
        self.ttree.Branch('decayVtxChi2', v_decayVtxChi2 , "decayVtxChi2/F")
        self.ttree.Branch('Distance1', v_Dist1, "Distance1/F")  #--
        self.ttree.Branch('Distance2', v_Dist2, "Distance2/F")  #--
        self.ttree.Branch('Distance3', v_Dist3, "Distance3/F")  #--
        self.ttree.Branch('mD0', v_mD0, "mD0/F")  # Change from mJPsi
        self.ttree.Branch('nSPD', v_nSPD, "nSPD/F")
        self.ttree.Branch('EffWeight', v_EffWeight, "EffWeight/F")
        self.ttree.Branch('EffWeight_0', v_EffWeight_0, "EffWeight_0/F")
        self.ttree.Branch('EffWeight_1', v_EffWeight_1, "EffWeight_1/F")
        self.ttree.Branch('EffWeight_2', v_EffWeight_2, "EffWeight_2/F")
        self.ttree.Branch('EffWeight_3', v_EffWeight_3, "EffWeight_3/F")
        self.ttree.Branch('EffWeight_4', v_EffWeight_4, "EffWeight_4/F")
        self.ttree.Branch('EffWeight_5', v_EffWeight_5, "EffWeight_5/F")  #To account for Data/MC difference of track correction
        self.ttree.Branch('EffWeight_0_Var', v_EffWeight_0_Var, "EffWeight_0_Var/F")
        self.ttree.Branch('EffWeight_1_Var', v_EffWeight_1_Var, "EffWeight_1_Var/F")
        self.ttree.Branch('EffWeight_2_Var', v_EffWeight_2_Var, "EffWeight_2_Var/F")
        self.ttree.Branch('EffWeight_3_Var', v_EffWeight_3_Var, "EffWeight_3_Var/F")
        self.ttree.Branch('EffWeight_4_Var', v_EffWeight_4_Var, "EffWeight_4_Var/F")
        self.ttree.Branch('EffWeight_5_Var', v_EffWeight_5_Var, "EffWeight_5_Var/F")
        self.ttree.Branch('isPrimary', v_isPrimary2, "isPrimary/F")
        self.ttree.Branch('KprobNNK', v_KprobNNK, "KprobNNK/F")
        self.ttree.Branch('KprobGhost', v_KprobGhost, "KprobGhost/F") 
        self.ttree.Branch('KTrckChi2', v_KTrckChi2, "KTrckChi2/F")

        if self.isMC==True:

            ################
            #Prepare components of RM tree
            #v_isPionAcc   = array.array('f',[0])
            #v_isMuonAcc   = array.array('f',[0])
            v_pTDet       = array.array('f',[0])
            v_etaDet      = array.array('f',[0])
            v_phiDet      = array.array('f',[0])
            v_etaTagDet   = array.array('f',[0])
            v_phiTagDet   = array.array('f',[0])
            v_nConstDet   = array.array('f',[0])
            v_TagPtDet    = array.array('f',[0])
            v_TagMDet     = array.array('f',[0])
            v_zTDet       = array.array('f',[0])
            v_pTPart      = array.array('f',[0])
            v_etaPart     = array.array('f',[0])
            v_phiPart     = array.array('f',[0])
            v_etaTagPart  = array.array('f',[0])
            v_phiTagPart  = array.array('f',[0])
            v_nConstPart  = array.array('f',[0])
            v_TagPtPart   = array.array('f',[0])
            v_TagMPart    = array.array('f',[0])
            v_zTPart      = array.array('f',[0])
            v_dR          = array.array('f',[0])
            v_PartnTag    = array.array('f',[0])
            v_DetnTag     = array.array('f',[0])

            #Response Matrix
            self.ttreeResponse = ROOT.TTree('Response', '')
            #self.ttreeResponse.Branch('isPionAcc',v_isPionAcc, "isPionAcc/F" )
            #self.ttreeResponse.Branch('isMuonAcc',v_isMuonAcc, "isMuonAcc/F" )
            self.ttreeResponse.Branch('jetPtDet', v_pTDet,    "jetPtDet/F")
            self.ttreeResponse.Branch('etaDet',   v_etaDet,   "etaDet/F")
            self.ttreeResponse.Branch('etaTagDet',   v_etaTagDet,   "etaTagDet/F")
            self.ttreeResponse.Branch('phiDet',   v_phiDet,   "phiDet/F")
            self.ttreeResponse.Branch('phiTagDet',   v_phiTagDet,   "phiTagDet/F")
            self.ttreeResponse.Branch('nConstDet',v_nConstDet,"nConstDet/F")
            self.ttreeResponse.Branch('tagPtDet', v_TagPtDet, "tagPtDet/F")
            self.ttreeResponse.Branch('tagMDet',  v_TagMDet,  "tagMDet/F")
            self.ttreeResponse.Branch('zTDet',    v_zTDet,    "zTDet/F")
            self.ttreeResponse.Branch('jetPtPart',v_pTPart,   "jetPtPart/F")
            self.ttreeResponse.Branch('etaPart',  v_etaPart,  "etaPart/F")
            self.ttreeResponse.Branch('etaTagPart',  v_etaTagPart,  "etaTagPart/F")
            self.ttreeResponse.Branch('phiPart',  v_phiPart,  "phiPart/F")
            self.ttreeResponse.Branch('phiTagPart',  v_phiTagPart,  "phiTagPart/F")
            self.ttreeResponse.Branch('nConstPart',v_nConstPart,"nConstPart/F")
            self.ttreeResponse.Branch('tagPtPart',v_TagPtPart,"tagPtPart/F")
            self.ttreeResponse.Branch('zTPart',   v_zTPart,   "zTPart/F")
            self.ttreeResponse.Branch('tagMPart', v_TagMPart, "tagMPart/F")
            self.ttreeResponse.Branch('dR', v_dR, "dR/F")
            self.ttreeResponse.Branch('isPrimaryPart', v_isPrimary2,   "isPrimaryPart/F")
            self.ttreeResponse.Branch('PartnTag',     v_PartnTag,     "PartnTag/F")
            self.ttreeResponse.Branch('DetnTag',      v_DetnTag,      "DetnTag/F")

            #Response Matrix
            self.ttreeResponse2 = ROOT.TTree('ResponseVar', '')
            #self.ttreeResponse2.Branch('isPionAcc',v_isPionAcc, "isPionAcc/F" )
            #self.ttreeResponse2.Branch('isMuonAcc',v_isMuonAcc, "isMuonAcc/F" )
            self.ttreeResponse2.Branch('jetPtDet', v_pTDet,    "jetPtDet/F")
            self.ttreeResponse2.Branch('etaDet',   v_etaDet,   "etaDet/F")
            self.ttreeResponse2.Branch('etaTagDet',   v_etaTagDet,   "etaTagDet/F")
            self.ttreeResponse2.Branch('phiDet',   v_phiDet,   "phiDet/F")
            self.ttreeResponse2.Branch('phiTagDet',   v_phiTagDet,   "phiTagDet/F")
            self.ttreeResponse2.Branch('nConstDet',v_nConstDet,"nConstDet/F")
            self.ttreeResponse2.Branch('tagPtDet', v_TagPtDet, "tagPtDet/F")
            self.ttreeResponse2.Branch('tagMDet',  v_TagMDet,  "tagMDet/F")
            self.ttreeResponse2.Branch('zTDet',    v_zTDet,    "zTDet/F")
            self.ttreeResponse2.Branch('jetPtPart',v_pTPart,   "jetPtPart/F")
            self.ttreeResponse2.Branch('etaPart',  v_etaPart,  "etaPart/F")
            self.ttreeResponse2.Branch('etaTagPart',  v_etaTagPart,  "etaTagPart/F")
            self.ttreeResponse2.Branch('phiPart',  v_phiPart,  "phiPart/F")
            self.ttreeResponse2.Branch('phiTagPart',  v_phiTagPart,  "phiTagPart/F")
            self.ttreeResponse2.Branch('nConstPart',v_nConstPart,"nConstPart/F")
            self.ttreeResponse2.Branch('tagPtPart',v_TagPtPart,"tagPtPart/F")
            self.ttreeResponse2.Branch('zTPart',   v_zTPart,   "zTPart/F")
            self.ttreeResponse2.Branch('tagMPart', v_TagMPart, "tagMPart/F")
            self.ttreeResponse2.Branch('dR', v_dR, "dR/F")
            self.ttreeResponse2.Branch('isPrimaryPart',v_isPrimary2,   "isPrimaryPart/F")
            self.ttreeResponse2.Branch('PartnTag',     v_PartnTag,     "PartnTag/F")
            self.ttreeResponse2.Branch('DetnTag',      v_DetnTag,      "DetnTag/F")

            ################
            #Particle level distribution
            #v_MCisPionAcc   = array.array('f',[0])
            #v_MCisMuonAcc   = array.array('f',[0])
            v_MCpTPart      = array.array('f',[0])
            v_MCetaPart     = array.array('f',[0])
            v_MCphiPart     = array.array('f',[0])
            v_MCetaTagPart  = array.array('f',[0])
            v_MCphiTagPart  = array.array('f',[0])
            v_MCnConstPart  = array.array('f',[0])
            v_MCTagPtPart   = array.array('f',[0])
            v_MCTagMPart    = array.array('f',[0])
            v_MCzTPart      = array.array('f',[0])
            v_MCtag_lifetPart = array.array('f',[0])
            v_isPrimary       = array.array('f',[0])
            v_isDetTagRec     = array.array('f',[0])

            self.mcttree = ROOT.TTree('MCFragmNtuple', '')
            #self.mcttree.Branch('isPionAcc',  v_MCisPionAcc,   "isPionAcc/F")
            #self.mcttree.Branch('isMuonAcc',  v_MCisMuonAcc,   "isMuonAcc/F")
            self.mcttree.Branch('jetPtPart',  v_MCpTPart,   "jetPtPart/F")
            self.mcttree.Branch('tagPtPart',  v_MCTagPtPart,"tagPtPart/F")
            self.mcttree.Branch('etaJetPart', v_MCetaPart,  "etaJetPart/F")
            self.mcttree.Branch('etaTagPart', v_MCetaTagPart, "etaTagPart/F")
            self.mcttree.Branch('phiJetPart', v_MCphiPart,  "phiJetPart/F")
            self.mcttree.Branch('phiTagPart', v_MCphiTagPart,  "phiTagPart/F")
            self.mcttree.Branch('nConstPart', v_MCnConstPart,"nConstPart/F")
            self.mcttree.Branch('zTPart',     v_MCzTPart,   "zTPart/F")
            self.mcttree.Branch('tag_lifetPart', v_MCtag_lifetPart, "tag_lifetPart/F")
            self.mcttree.Branch('isPrimary',  v_isPrimary,   "isPrimary/F")
            self.mcttree.Branch('isDetTagRec',v_isDetTagRec, "isDetTagRec/F")
            self.mcttree.Branch('tagMPart',   v_MCTagMPart,  "tagMPart/F")



        ################
        v_goodDetTags_idx     = ROOT.std.vector('int')()

        ################
        #To load the info from the input tree
        evt_n         = array.array('f',[0])
        tag_pid       = array.array('f',[0])
        tag_idx_jet   = array.array('f',[0])

        #========================================================================================
        # Event loop
        #========================================================================================
        totalNbOfEntries = self.ttreeOriginal.GetEntries()
        #totalNbOfEntries = 10000 #TODO this can be changed for lower statistics
        print("oo Start loop over ntuple to extract important quantities")
        print("oo Tree has {} entries".format(totalNbOfEntries))

        if totalNbOfEntries==0: return
        #-Determine frequency of printout (printout should update after about 10% of events)
        printFreq = int(totalNbOfEntries/self.printLvl)
        if printFreq<1:
            printFreq=1

        nTagPartLvlEvt=0
        nTagDetLvlEvt =0
        lastPV = None
 
 
        for iEvt in range(0, totalNbOfEntries+1):
      
            self.ttreeOriginal.GetEntry(iEvt)
            if iEvt % printFreq==0:
                print(" Events processed: {} ({}%)".format(iEvt,round(100*iEvt/totalNbOfEntries,0)))

            evt_n[0]     = self.ttreeOriginal.evt_n #do I need this for something?
            v_nSPD[0]    = int(self.ttreeOriginal.evt_spd)
            nTagDetLvl   = int(self.ttreeOriginal.tag_pid.size())
            nJetsDetLvl  = int(self.ttreeOriginal.jet_px.size())
            if self.isMC==True:
                nTagPartLvl  = int(self.ttreeOriginal.mctag_pid.size())
                nJetsPartLvl = int(self.ttreeOriginal.mcjet_px.size())

            if self.isMC and nTagPartLvl:
                nTagPartLvlEvt+=1
            if nTagDetLvl:
                nTagDetLvlEvt+=1

            # print("\n nJetsDetLvlv: {}, nTags DetLvl: {}".format(nJetsDetLvl,nTagDetLvl))
            # print("n Jets PartLvl: {}, nTagsPartLvl: {}, nJetsDetLvlv: {}, nTags DetLvl: {}".format(nJetsPartLvl,nTagPartLvl,nJetsDetLvl,nTagDetLvl))

            #- - - - - - - - - - - - - - - - - - - - - - - - - - -
            #- Step 1 - Check for good events
            #- - - - - - - - - - - - - - - - - - - - - - - - - - -
            if self.isMC and nTagPartLvl==0: continue

            #Fill unbiased PartLvl spectrum - for Eff determination
            #mcttree = fillPartLvlTree(mcttree,ttreeOriginal)
            if self.isMC:
                self.mcttree = fillPartLvlTree(self.mcttree,self.ttreeOriginal,v_MCetaPart,v_MCphiPart,v_MCpTPart,v_MCnConstPart,v_MCetaTagPart, v_MCphiTagPart,v_MCTagPtPart,v_MCTagMPart,v_MCzTPart, v_MCtag_lifetPart, v_isPrimary, v_isDetTagRec)
 
      
            if nTagDetLvl==0:
                # print("No tags in event")
                continue

            tag_idxJetSize  = int(self.ttreeOriginal.tag_idx_jet.size())
            tag_pidSize     = int(self.ttreeOriginal.tag_pid.size())
            pvrNumber       = int(self.ttreeOriginal.pvr_n)
      
            if pvrNumber!= 1: continue
            if (tag_idxJetSize!=tag_pidSize):
                # print("length tag_idx_jet={}, length tag_pid={}".format(tag_idxJetSize,tag_pidSize))
                continue
            # print("length tag_idx_jet_after={}, length tag_pid_after={}".format(tag_idxJetSize,tag_pidSize))
         
            #- - - - - - - - - - - - - - - - - - - - - - - - - - -
            #- Step 2 - Collect all the good detector level tags in the event
            #- - - - - - - - - - - - - - - - - - - - - - - - - - -
            v_goodDetTags_idx.clear()
            for iTagDetLvl in range(0, nTagDetLvl):

                tag_l0_tos0    = int(self.ttreeOriginal.tag_l0_tos0[iTagDetLvl])
                tag_hlt1_tos0  = int(self.ttreeOriginal.tag_hlt1_tos0[iTagDetLvl])
                jetNumberTurbo = int(self.ttreeOriginal.tag_idx_jet[iTagDetLvl]) #This is the jetID the tag is associated to
                ##tagNumberGen   = int(ttreeOriginal.tag_idx_MCtag[iTagDetLvl]) #This is the generator tag ID the Detector level tag is associated to
                tag_pid        = int(self.ttreeOriginal.tag_pid[iTagDetLvl])
                if self.isMC:
                    tag_assocType  = int(self.ttreeOriginal.tag_genAssocType[iTagDetLvl])

                #. . . . . . . . . . . . . . .
                # Tag filtering
                #. . . . . . . . . . . . . . .
                #if nprtPID != nprtpnn_mu: continue #Q what is this for?
                if jetNumberTurbo == -1: 
                    # print("jetNumberTurbo == -1")
                    continue
                # else:
                #     print("jetNumberTurbo != -1")

                if self.isMC and tag_assocType    > 2: 
                    # print("tag_assocType > 2")
                    continue #the way how the MC det lvl particle was matched to the gen. lvl. particle
                # if self.isMC==False: #TODO add this line back with proper trigger
                    #print("tag_l0_tos: {}, tag_hlt1_tos:{}, jetNumberTurbo:{}, tag_assocType:{}".format(tag_l0_tos0,tag_hlt1_tos0,jetNumberTurbo,tag_assocType))
                    # if tag_l0_tos0     != 1: continue #TODO add this line back with proper trigger
                    # if tag_hlt1_tos0   != 1: continue #TODO add this line back with proper trigger
              
                v_goodDetTags_idx.push_back(int(iTagDetLvl))
                #-only use the event if it contains a good det Lvl tag.
            if v_goodDetTags_idx.size()==0:
                # print("No good tags in event")
                continue

            #- - - - - - - - - - - - - - - - - - - - - - - - - - -
            #- Step 3 - Get all info for the Tag inside Jet analysis
            #- - - - - - - - - - - - - - - - - - - - - - - - - - -
            for iTagDetLvl in v_goodDetTags_idx:

                jetNumberTurbo   = int(self.ttreeOriginal.tag_idx_jet[iTagDetLvl]) #This is the jet ID the tag is associated to
                pvrNumberTurbo = int(self.ttreeOriginal.tag_idx_pvr[iTagDetLvl])
                prt0Number     = int(self.ttreeOriginal.tag_idx_prt0[iTagDetLvl])
                prt1Number     = int(self.ttreeOriginal.tag_idx_prt1[iTagDetLvl])
                
                #Check if the MC linked particle is a primary particle
                if self.isMC:
                    partLvlTagNumber = int(self.ttreeOriginal.tag_idx_MCtag[iTagDetLvl]) #This is the partLvl tag ID the tag is associated to
                    v_isPrimary2[0]  = self.ttreeOriginal.mctag_isPrimary[partLvlTagNumber]
                else:
                    v_isPrimary2[0] = -1

                partList = [prt0Number, prt1Number]
                
                #. . . . . . . . . . . . . . .
                # Event observables
                #. . . . . . . . . . . . . . .
                pid0 = self.ttreeOriginal.prt_pid[prt0Number]
                pid1 = self.ttreeOriginal.prt_pid[prt1Number]
                
                #. . . . . . . . . . . . . . .
                # Tag observables
                #. . . . . . . . . . . . . . .
                dRTagTurbo = self.ttreeOriginal.tag_dR_jet[iTagDetLvl]
                try:
                    nTRnd = self.ttreeOriginal.tag_n_tagsRnd[iTagDetLvl]
                except:
                    nTRnd = -1
                    
                p4TagTurbo = getFourVector(self.ttreeOriginal, "tag", iTagDetLvl)
                if not isAccepted(p4TagTurbo, 0):
                    # print("Tag not accepted")
                    continue
               
                try: #This only works for data and newer MC Davinci files
                    v_Dist1[0] = self.ttreeOriginal.tag_DOCA1[iTagDetLvl]
                    v_Dist2[0] = self.ttreeOriginal.tag_DOCA2[iTagDetLvl]
                    v_Dist3[0] = self.ttreeOriginal.tag_DOCA3[iTagDetLvl]
                except:
                    v_Dist1[0] = -1
                    v_Dist2[0] = -1
                    v_Dist3[0] = -1
                    
                v_tagdR[0]          =dRTagTurbo
                v_tagMass[0]        =p4TagTurbo.M()
                v_tagPt[0]          =p4TagTurbo.Pt()
                v_tagEta[0]         =p4TagTurbo.PseudoRapidity()
                v_tag_idx_jet[0]    =jetNumberTurbo
                v_tagnRnd[0]        =nTRnd
                try:
                    v_decayVtxChi2[0]   =self.ttreeOriginal.tag_decVtxChi2[iTagDetLvl]
                    v_tag_decVtxChi2[0] =self.ttreeOriginal.tag_decVtxChi2[iTagDetLvl]
                    #set the logarithm of the decay vertex chi2 via calculation
                    v_tag_logdecVtxChi2[0] = math.log10(v_tag_decVtxChi2[0])
                except:
                    v_decayVtxChi2[0]   =-1
                    v_tag_decVtxChi2[0] =-1
                    v_tag_logdecVtxChi2[0] =-1

          
                pvrZ       = self.ttreeOriginal.pvr_z[pvrNumberTurbo] #-primary vertex of particle
                decayZ     = self.ttreeOriginal.tag_z[iTagDetLvl]      #-decay vertex of particle
                pz         = self.ttreeOriginal.tag_pz[iTagDetLvl]

                lifetime        = calcLifeTime(p4TagTurbo.M(),pvrZ,decayZ,pz)
                if lastPV!=None:
                    pvrZWrong     = lastPV
                    lifetimeWrong = calcLifeTime(p4TagTurbo.M(),pvrZWrong,decayZ,pz)
                else:
                    lifetimeWrong = -50
                lastPV     = pvrZ #Set the current primary vertex as the wrong PV to be used in next loop

                v_tag_lifet[0]        =lifetime
                v_tag_lifetWrongPV[0] =lifetimeWrong
              
                #. . . . . . . . . . . . . . .
                # Jet observables
                #. . . . . . . . . . . . . . .
                p4JetTurbo = getFourVector(self.ttreeOriginal,"jet",jetNumberTurbo)
                if not isAccepted(p4JetTurbo,1):
                    # print("Jet not accepted")
                    continue

                nTComb     = self.ttreeOriginal.jet_n_tagsComb[jetNumberTurbo]
                nTUnique   = self.ttreeOriginal.jet_n_tagsUnique[jetNumberTurbo]
                nNeu       = self.ttreeOriginal.jet_n_neu[jetNumberTurbo]
                nCh        = self.ttreeOriginal.jet_n_chr[jetNumberTurbo]

                #print("jetNumberTurbo: {}, pTJet:{}".format(jetNumberTurbo,p4JetTurbo.Pt()))
                #print("    pT jet in DetLvl jet position {}".format(v_SelDetJet.at(jetNumberTurbo).Pt()))

                #also for jets: nConstituents in jet
                v_jetnConst[0]   =nNeu+nCh
                v_jetnTComb[0]   =nTComb
                v_jetnTUnique[0] =nTUnique
                v_jetPt[0]       =p4JetTurbo.Pt()
                #v_jetY[0]        =p4JetTurbo.Rapidity()
                v_jetEta[0]      =p4JetTurbo.PseudoRapidity()
                v_tagZ[0]        =p4TagTurbo.Pt()/p4JetTurbo.Pt()

                #. . . . . . . . . . . . . . .
                # Tag daughter observables
                #. . . . . . . . . . . . . . .
                #In case you want to know any daughter particle information
                foundD0=True
                foundK=False
                foundPiP=False
      
                for partNr in partList:
                    # Look for D0 meson (421) instead of J/Psi (443)
                    if self.ttreeOriginal.prt_pid[partNr] == 421:
                        p4D0Turbo = getFourVector(self.ttreeOriginal, "prt", partNr)
                        if isAccepted(p4D0Turbo):
                            foundD0 = True
                        v_D0Z[0] = p4D0Turbo.Pt()/p4JetTurbo.Pt()
                        
                    # Look for K- (-321) instead of mu-
                    elif self.ttreeOriginal.prt_pid[partNr] == -321:
                        p4KTurbo = getFourVector(self.ttreeOriginal, "prt", partNr)
                        try:
                            v_KprobNNK[0] = self.ttreeOriginal.prt_pnn_k[partNr]
                            v_KprobGhost[0] = self.ttreeOriginal.prt_prb_ghost[partNr]
                            v_KTrckChi2[0] = self.ttreeOriginal.prt_trckChi2[partNr]
                        except:
                            v_KprobNNK[0] = -10
                            v_KprobGhost[0] = 10
                            v_KTrckChi2[0] = -1
                        if isAccepted(p4KTurbo):
                            foundK = True
                            
                    # Keep pion (211) handling
                    elif self.ttreeOriginal.prt_pid[partNr] == 211:
                        p4PipTurbo = getFourVector(self.ttreeOriginal, "prt", partNr)
                        try:
                            v_piPprobNNpi[0] = self.ttreeOriginal.prt_pnn_pi[partNr]
                            v_piPprobGhost[0] = self.ttreeOriginal.prt_prb_ghost[partNr]
                            v_piPTrckChi2[0] = self.ttreeOriginal.prt_trckChi2[partNr]
                        except:
                            v_piPprobNNpi[0] = -10 #cut on > 0.9
                            v_piPprobGhost[0] = 10 #cut on < 0.3
                            v_piPTrckChi2[0] = -1
                        if isAccepted(p4PipTurbo):
                            foundPiP = True

                # 3. Updated condition for filling event data - check for D0, K and pi
                if foundD0 and foundK and foundPiP:
                    pass #TODO this needs to be fixed
                    # print("Found all particles!")
                    # Fill D0 information
                    # v_mD0[0] = p4D0Turbo.M() #TODO uncomment this!!
                    # do nothing
                    # Remove all the J/Psi and muon specific code
                    # No need for v_mpipi or v_QValue
                    # Update the efficiency calculation to handle D0 → K pi instead of J/Psi→μμ decays
                    # if not self.isMC:
                    #     v_EffWeight[0] = self.getD0EffCorrWeight(p4KTurbo, p4PipTurbo, -1) #TODO uncomment this!!
                    #     v_EffWeight_0[0] = self.getD0EffCorrWeight(p4KTurbo, p4PipTurbo, 0) #TODO uncomment this!!
                        # etc. for the other efficiency weights
                else:
                    # print("Did not find all particles!")
                    # print("D0 {}, K- {}, pi+ {}".format(foundD0, foundK, foundPiP))
                    continue


                self.ttreeOriginal.GetEntry(iEvt)#Q: what is that??
                self.ttree.Fill()

            if self.isMC==True:
                #- - - - - - - - - - - - - - - - - - - - - - - - - - -
                #- Step 4 - Build a matrix of Part Lvl Jet-Tag relations
                #- - - - - - - - - - - - - - - - - - - - - - - - - - -
                #           One jet (although rare) can have multiple tags
                #           One Tag can only relate to one jet
                #           Thus: You can uniquley relate the tags to the jet but jet->tag relation
                #                 is not always clear because there COULD be multiple per jet
                #
                #           Example:
                #                  Jet1    Jet2    Jet 3   Jet 4
                #           Tag1             1
                #           Tag2    1
                #           Tag3             1
                #           Tag4                             1
                #- - - - - - - - - - - - - - - - - - - - - - - - - - -

                #-Detector level
                m_TagJetDetLvl = np.zeros([nTagDetLvl,nJetsDetLvl])
                for DetTag_idx in v_goodDetTags_idx:
                    jetNumber_idx = int(self.ttreeOriginal.tag_idx_jet[DetTag_idx]) #This is the jetID the tag is associated to
                    if jetNumber_idx>-1:
                      m_TagJetDetLvl[DetTag_idx,jetNumber_idx] = 1

                #-Generator level
                if self.isMC==True:
                    m_TagJetPartLvl = np.zeros([nTagPartLvl,nJetsPartLvl])
                    for PartTag_idx in range(0, nTagPartLvl):
                        jetNumber_idx = int(self.ttreeOriginal.mctag_idx_jet[PartTag_idx]) #This is the jetID the tag is associated to
                        if jetNumber_idx>-1:
                            m_TagJetPartLvl[PartTag_idx,jetNumber_idx] = 1
                
                    
                #- - - - - - - - - - - - - - - - - - - - - - - - - - -
                #- Step 5 - Build a matrix Jet-Jet distance relations
                #- regard speacial conditions for matching
                #- - - - - - - - - - - - - - - - - - - - - - - - - - -
                if nJetsPartLvl>0 and nJetsDetLvl>0:
                  
                  if self.matchingVersion==1 or  self.matchingVersion==3:
                      #-Declare a matrix of size nSelPartLvlJets and size nSelDetLvlJets
                      #-and fill it with all distances between the two sets
                      m_JetJetDistance = np.ones([nJetsPartLvl,nJetsDetLvl])
                      m_JetJetDistance*=10


                      for partJetidx in range(0, nJetsPartLvl):
                          containsPartLvlTag=False
                          #check if this jet contains a tag:
                          ListWithPartTags = m_TagJetPartLvl[:,partJetidx]
                          if 1 in ListWithPartTags: containsPartLvlTag=True
                          if self.matchingVersion==1 and containsPartLvlTag==False: continue
                          
                          p4JetPart = getFourVector(self.ttreeOriginal,"mcjet",partJetidx)
                          for detJetidx in range(0, nJetsDetLvl):
                              containsDetLvlTag=False
                              #check if this jet contains a tag:
                              ListWithDetTags = m_TagJetDetLvl[:,detJetidx]
                              if 1 in ListWithDetTags:  containsDetLvlTag=True
                              if self.matchingVersion==1 and containsDetLvlTag==False: continue

                              p4JetDet  = getFourVector(self.ttreeOriginal,"jet",detJetidx)
                              distance  = p4JetPart.DeltaR(p4JetDet)
                              m_JetJetDistance[partJetidx,detJetidx] = distance
                  if self.matchingVersion==2:
                    print("not implemented yet")
                
                #- - - - - - - - - - - - - - - - - - - - - - - - - - -
                #- Step 5 - Fill Response matrix
                #- - - - - - - - - - - - - - - - - - - - - - - - - - -
                #- - - - - - - - - - - - - - - - - - - - - - - -
                #-Retrieve the:
                #               (in case 1 and 3) closest jet-jet pair
                #               (in case 2) 2 jets that contain the matched tag at part and detector level
                for partJetidx in range(0, nJetsPartLvl):

                    closestDetJetidx=-1
                    if self.matchingVersion==1 or self.matchingVersion==3:
                        closestDetJetidx = findClosest(m_JetJetDistance,row=partJetidx,version=1) #check if they are uniquely the closest to each other
 
                        if closestDetJetidx==-1: continue #No closest jet found
                        p4PartJet = getFourVector(self.ttreeOriginal,"mcjet",partJetidx)
                        p4DetJet  = getFourVector(self.ttreeOriginal,"jet",closestDetJetidx)

                        ListWithPartTags = m_TagJetPartLvl[:,partJetidx]
                        ListWithDetTags  = m_TagJetDetLvl[:,closestDetJetidx]
                        #-get Part Lvl Tag index
                        arrayWithInicesP = np.where(ListWithPartTags==1) #-find first index with value (1)
                        idxTagPart       = int(arrayWithInicesP[0][0])

                        #-Check if the generator level tag fulfills all the fiducial requirements
                        #-pt and p for the individual particles
                        particlesInAcceptance = False
                        try:
                            # The kaon and pion from D0
                            prt1Number = int(self.ttreeOriginal.mctag_idx_prt1[idxTagPart])
                            prt2Number = int(self.ttreeOriginal.mctag_idx_prt2[idxTagPart])
                            pid1 = self.ttreeOriginal.mcprt_pid[prt1Number]
                            pid2 = self.ttreeOriginal.mcprt_pid[prt2Number]
                            
                            # Check if we have K(-321) and pi(211)
                            if (pid1 == -321 and pid2 == 211) or (pid1 == 211 and pid2 == -321):
                                p4MCkaon = getFourVector(self.ttreeOriginal, "mcprt", 
                                                       prt1Number if pid1 == -321 else prt2Number)
                                p4MCpion = getFourVector(self.ttreeOriginal, "mcprt", 
                                                       prt1Number if pid1 == 211 else prt2Number)
                                
                                # Appropriate D0 daughter kinematic cuts
                                if p4MCkaon.Pt() > 0.25 and p4MCpion.Pt() > 0.25 and \
                                   p4MCkaon.P() > 2.0 and p4MCpion.P() > 2.0:
                                    particlesInAcceptance = True
                        except:
                            particlesInAcceptance = False

                        #-Do not accept events at generator level if the particles are not within the fiducial volume cuts
                        if not particlesInAcceptance: continue

                        #v_isPionAcc[0]     = pionsAcc
                        #v_isMuonAcc[0]     = muonsAcc

                        #-get Det Lvl Tag index
                        arrayWithInicesD  = np.where(ListWithDetTags==1)  #-find first index with value (1)
                        idxTagDet = int(arrayWithInicesD[0][0])
                        #for tag in ListWithDetTags:
                        #  print("Tag idx: {}".format(tag))
                        #print("index det: {}, index part: {}".format(idxTagDet,idxTagPart))
                        p4PartTag = getFourVector(self.ttreeOriginal,"mctag",idxTagPart)
                        p4DetTag  = getFourVector(self.ttreeOriginal,"tag",idxTagDet)

                        v_pTDet[0]       = p4DetJet.Pt()
                        v_etaDet[0]      = p4DetJet.Eta()
                        v_phiDet[0]      = p4DetJet.Phi()
                        v_nConstDet[0]   = self.ttreeOriginal.jet_n_neu[closestDetJetidx]+self.ttreeOriginal.jet_n_chr[closestDetJetidx]
                        v_TagPtDet[0]    = p4DetTag.Pt()
                        v_etaTagDet[0]   = p4DetTag.Eta()
                        v_phiTagDet[0]   = p4DetTag.Phi()
                        v_TagMDet[0]     = p4DetTag.M()
                        v_zTDet[0]       = p4DetTag.Pt()/p4DetJet.Pt()
                        v_pTPart[0]      = p4PartJet.Pt()
                        v_etaPart[0]     = p4PartJet.Eta()
                        v_phiPart[0]     = p4PartJet.Phi()
                        v_nConstPart[0]  = self.ttreeOriginal.mcjet_n_neu[partJetidx]+self.ttreeOriginal.mcjet_n_chr[partJetidx]
                        v_TagPtPart[0]   = p4PartTag.Pt()
                        v_etaTagPart[0]  = p4PartTag.Eta()
                        v_phiTagPart[0]  = p4PartTag.Phi()
                        v_TagMPart[0]    = p4PartTag.M()
                        v_zTPart[0]      = p4PartTag.Pt()/p4PartJet.Pt()
                        v_dR[0]          = m_JetJetDistance[partJetidx,closestDetJetidx]
                        v_isPrimary2[0]  = self.ttreeOriginal.mctag_isPrimary[idxTagPart] #-is the mcTag a primary particle
                        v_PartnTag[0]    = len(ListWithPartTags)
                        v_DetnTag[0]     = len(ListWithDetTags)
                        #print("distance: {}, pT MC: {}, pT Data: {}".format(distance,mcJet.Pt(),dataJet.Pt()))

                        self.ttreeResponse.Fill()

                #- - - - - - - - - - - - - - - - - - - - - - - -
                #- - - - - - - - - - - - - - - - - - - - - - - -
                #VARIATION of the matching
                for partJetidx in range(0, nJetsPartLvl):
  
                    #Get list with tags matched to that jet
                    ListWithPartTags = m_TagJetPartLvl[:,partJetidx]
                    arrayWithInicesP = np.where(ListWithPartTags==1) #-find first index with value (1)
                    if len(arrayWithInicesP[0])==0: continue
                    idx_TagPart       = int(arrayWithInicesP[0][0])
                    #-Get Gen lvl properties
                    p4PartJet = getFourVector(self.ttreeOriginal,"mcjet",partJetidx)
                    p4PartTag = getFourVector(self.ttreeOriginal,"mctag",idx_TagPart)

                    #-Find the det lvl. tag that is matched to the gen lvl. tag
                    idx_TagDet = int(self.ttreeOriginal.mctag_idx_Dettag[idx_TagPart]) #This is the Det lvl tag the MC tag is associated to
                    if idx_TagDet<0: continue  # There is no match for that Gen lvl. tag
                    idx_jetDet = int(self.ttreeOriginal.tag_idx_jet[idx_TagDet])       #This is the jetID the tag is associated to
                    ListWithDetTags  = m_TagJetDetLvl[:,idx_jetDet]

                    #-Get Det lvl properties
                    p4DetJet   = getFourVector(self.ttreeOriginal,"jet",idx_jetDet)
                    p4DetTag   = getFourVector(self.ttreeOriginal,"tag",idx_TagDet)


                    #-Check if the generator level tag fulfills all the fiducial requirements
                    #- pt and p for the individual particles
                    particlesInAcceptance = False
                    try:
                        # The kaon and pion from D0
                        prt1Number = int(self.ttreeOriginal.mctag_idx_prt1[idx_TagDet])
                        prt2Number = int(self.ttreeOriginal.mctag_idx_prt2[idx_TagDet])
                        pid1 = self.ttreeOriginal.mcprt_pid[prt1Number]
                        pid2 = self.ttreeOriginal.mcprt_pid[prt2Number]
                        
                        # Check if we have K(-321) and pi(211)
                        if (pid1 == -321 and pid2 == 211) or (pid1 == 211 and pid2 == -321):
                            p4MCkaon = getFourVector(self.ttreeOriginal, "mcprt", 
                                                   prt1Number if pid1 == -321 else prt2Number)
                            p4MCpion = getFourVector(self.ttreeOriginal, "mcprt", 
                                                   prt1Number if pid1 == 211 else prt2Number)
                            
                            # Appropriate D0 daughter kinematic cuts
                            if p4MCkaon.Pt() > 0.25 and p4MCpion.Pt() > 0.25 and \
                               p4MCkaon.P() > 2.0 and p4MCpion.P() > 2.0:
                                particlesInAcceptance = True
                    except:
                        particlesInAcceptance = False

                    #-Do not accept events at generator level if the particles are not within the fiducial volume cuts
                    if not particlesInAcceptance: continue

                    #v_isPionAcc[0]   = pionsAcc
                    #v_isMuonAcc[0]   = muonsAcc

                    v_pTDet[0]       = p4DetJet.Pt()
                    v_etaDet[0]      = p4DetJet.Eta()
                    v_phiDet[0]      = p4DetJet.Phi()
                    v_nConstDet[0]   = self.ttreeOriginal.jet_n_neu[idx_jetDet]+self.ttreeOriginal.jet_n_chr[idx_jetDet]
                    v_TagPtDet[0]    = p4DetTag.Pt()
                    v_etaTagDet[0]   = p4DetTag.Eta()
                    v_phiTagDet[0]   = p4DetTag.Phi()
                    v_TagMDet[0]     = p4DetTag.M()
                    v_zTDet[0]       = p4DetTag.Pt()/p4DetJet.Pt()
                    #
                    v_pTPart[0]      = p4PartJet.Pt()
                    v_etaPart[0]     = p4PartJet.Eta()
                    v_phiPart[0]     = p4PartJet.Phi()
                    v_nConstPart[0]  = self.ttreeOriginal.mcjet_n_neu[partJetidx]+self.ttreeOriginal.mcjet_n_chr[partJetidx]
                    v_TagPtPart[0]   = p4PartTag.Pt()
                    v_etaTagPart[0]  = p4PartTag.Eta()
                    v_phiTagPart[0]  = p4PartTag.Phi()
                    v_TagMPart[0]    = p4PartTag.M()
                    if p4PartJet.Pt()>0:
                        v_zTPart[0]      = p4PartTag.Pt()/p4PartJet.Pt()
                    else:
                        v_zTPart[0]      = 0
                    #
                    v_dR[0]          = p4PartJet.DeltaR(p4DetJet)
                    v_isPrimary2[0]  = self.ttreeOriginal.mctag_isPrimary[idx_TagPart]
                    v_PartnTag[0]    = len(ListWithPartTags)
                    v_DetnTag[0]     = len(ListWithDetTags)
                    #print("distance: {}, pT MC: {}, pT Data: {}".format(distance,mcJet.Pt(),dataJet.Pt()))

                    self.ttreeResponse2.Fill()


        
        
        #After the event loop - save the trees
        self.tfile.cd()
        self.ttree.Write()
        if self.isMC==True:
            self.ttreeResponse.Write()
            self.ttreeResponse2.Write()
            self.mcttree.Write()
        self.tfile.Write()

        self.tfileOriginal.Close()
        self.tfile.Purge()
        self.tfile.Close()

        print("Finished filtering file: {}.root".format(self.fFileName))
        print("Saved into new file: {}".format(self.foutName))
        print("Found: {} Evts with PartLvl tags, found: {} Evts with DetLvl Tags".format(nTagPartLvlEvt,nTagDetLvlEvt))
    
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
    def getEffCorrWeight(self, piVector1,piVector2,muonVect1,muonVect2,triggVect, type=-1):


        #. . . . . . . . . . . . . . . . . . . . . .
        #-get info from pion particles
        p1     = piVector1.P()*1000#P is in GeV ,Maps are in MeV
        pt1    = piVector1.Pt()*1000#Pt is in GeV ,Maps are in MeV
        eta1   = piVector1.PseudoRapidity()
        
        p2     = piVector2.P()*1000#P is in GeV ,Maps are in MeV
        pt2    = piVector2.Pt()*1000#Pt is in GeV ,Maps are in MeV
        eta2   = piVector2.PseudoRapidity()

        #. . . . . . . . . . . . . . . . . . . . . .
        #-get info from muon particles
        pt3    = muonVect1.Pt()*1000#Pt is in GeV ,Maps are in MeV
        p3     = muonVect1.P()*1000#Pt is in GeV ,Maps are in MeV
        eta3   = muonVect1.PseudoRapidity()
        pt4    = muonVect2.Pt()*1000#Pt is in GeV ,Maps are in MeV
        p4     = muonVect2.P()*1000#Pt is in GeV ,Maps are in MeV
        eta4   = muonVect2.PseudoRapidity()
        
        #. . . . . . . . . . . . . . . . . . . . . .
        #-trigger efficiency
        #-get info from the trigger (J/Psi)
        ptComb   = math.sqrt(pt3*pt4)
        etaTrig  = triggVect.PseudoRapidity()
        

        #. . . . . . . . . . . . . . . . . . . . . .
        #-Return efficiencies
        #-Type =-1 -> all combined
        #-Type =0 -> pion reconstr efficiencies
        #-Type =1 -> pion selection efficiencies
        #-Type =2 -> muon reconstr efficiencies
        #-Type =3 -> trig line selection correction
        #-Type =4 -> Trigger efficiency  (trigger eff)
        #-Type =5 -> Correction R_Data/MC applied to each particle
        eff = 1
        if type==-1:
            #eff = (piRecoEff1*piRecoEff2*piSelEff1*piSelEff2*muEff1*muEff2*effStrippLineCorr1*effStrippLineCorr2*triggEff)  #*triggerEff * triggerSelEff
            eff = 1 #-not used
        if type==0:
            piRecoEff1 = self.pionMapReco.GetBinContent(self.pionMapReco.FindBin(pt1,eta1))
            piRecoEff2 = self.pionMapReco.GetBinContent(self.pionMapReco.FindBin(pt2,eta2))
            eff = (piRecoEff1*piRecoEff2)
        #-This is the variation of type=0
        if type==10:
            piRecoEff1 = self.pionMapReco_Var.GetBinContent(self.pionMapReco_Var.FindBin(pt1,eta1))
            piRecoEff2 = self.pionMapReco_Var.GetBinContent(self.pionMapReco_Var.FindBin(pt2,eta2))
            eff = (piRecoEff1*piRecoEff2)
        if type==1:
            piSelEff1  = self.pionMapSelHist.GetBinContent(self.pionMapSelHist.FindBin(p1,eta1))
            piSelEff2  = self.pionMapSelHist.GetBinContent(self.pionMapSelHist.FindBin(p2,eta2))
            eff = (piSelEff1*piSelEff2)
        #-This is the variation of type=1
        if type==11:
            piSelEff1  = self.pionMapSelHist_Var.GetBinContent(self.pionMapSelHist_Var.FindBin(p1,eta1))
            piSelEff2  = self.pionMapSelHist_Var.GetBinContent(self.pionMapSelHist_Var.FindBin(p2,eta2))
            eff = (piSelEff1*piSelEff2)
        if type==2:
            muEff1 = self.MuonMap.GetBinContent(self.MuonMap.FindBin(pt3,eta3))
            muEff2 = self.MuonMap.GetBinContent(self.MuonMap.FindBin(pt4,eta4))
            eff = (muEff1*muEff2)
        #-This is the variation of type=2
        if type==12:
            muEff1 = self.MuonMap_Var.GetBinContent(self.MuonMap_Var.FindBin(pt3,eta3))
            muEff2 = self.MuonMap_Var.GetBinContent(self.MuonMap_Var.FindBin(pt4,eta4))
            eff = (muEff1*muEff2)
        if type==3:
            effStrippLineCorr1 = self.TriggerSelCorrMap.GetBinContent(self.TriggerSelCorrMap.FindBin(eta3,pt3))
            effStrippLineCorr2 = self.TriggerSelCorrMap.GetBinContent(self.TriggerSelCorrMap.FindBin(eta4,pt4))
            if effStrippLineCorr1>1: effStrippLineCorr1=0
            if effStrippLineCorr2>1: effStrippLineCorr2=0
            eff = (effStrippLineCorr1*effStrippLineCorr2)
        #-This is the variation of type=3
        if type==13:
            effStrippLineCorr1 = self.TriggerSelCorrMap_Var.GetBinContent(self.TriggerSelCorrMap_Var.FindBin(eta3,pt3))
            effStrippLineCorr2 = self.TriggerSelCorrMap_Var.GetBinContent(self.TriggerSelCorrMap_Var.FindBin(eta4,pt4))
            if effStrippLineCorr1>1: effStrippLineCorr1=0
            if effStrippLineCorr2>1: effStrippLineCorr2=0
            eff = (effStrippLineCorr1*effStrippLineCorr2)
        if type==4:
            #Get the eta bin number
            etaBin   = self.TriggerEffMap2D.GetYaxis().FindBin(etaTrig)
            triggEff = self.functionList[etaBin-1].Eval(ptComb)
            eff = triggEff
        #-This is the variation of type=4
        if type==14:
            #Get the eta bin number
            etaBin   = self.TriggerEffMap2D_Var.GetYaxis().FindBin(etaTrig)
            triggEff = self.functionList_Var[etaBin-1].Eval(ptComb)
            eff = triggEff
        if type==5:
            #-R_Data/MC correction
            pion1DataMCcorr = self.Ratio_D_MC.getCorr(eta = eta1, pt = pt1, p = p1)
            pion2DataMCcorr = self.Ratio_D_MC.getCorr(eta = eta2, pt = pt2, p = p2)
            muon1DataMCcorr = self.Ratio_D_MC.getCorr(eta = eta3, pt = pt3, p = p3)
            muon2DataMCcorr = self.Ratio_D_MC.getCorr(eta = eta4, pt = pt4, p = p4)

            eff = pion1DataMCcorr*pion2DataMCcorr*muon1DataMCcorr*muon2DataMCcorr
        #-This is the variation of type=5
        if type==15:
            #-R_Data/MC correction
            pion1DataMCcorr = self.Ratio_D_MC.getCorr_Var(eta = eta1, pt = pt1, p = p1)
            pion2DataMCcorr = self.Ratio_D_MC.getCorr_Var(eta = eta2, pt = pt2, p = p2)
            muon1DataMCcorr = self.Ratio_D_MC.getCorr_Var(eta = eta3, pt = pt3, p = p3)
            muon2DataMCcorr = self.Ratio_D_MC.getCorr_Var(eta = eta4, pt = pt4, p = p4)

            eff = pion1DataMCcorr*pion2DataMCcorr*muon1DataMCcorr*muon2DataMCcorr
        
        #- - - - - - -
        if eff>0:
            effCorrection = 1./eff
        else:
            effCorrection = 0
        return effCorrection

    def getD0EffCorrWeight(self, kaonVector, pionVector, type=-1):
        """
        Calculate efficiency correction weights for D0 → K- π+ decay
        
        Args:
            kaonVector: TLorentzVector of the kaon
            pionVector: TLorentzVector of the pion
            type: Type of correction to apply:
                -1: Combined correction (returns 1 as placeholder)
                0: Kaon/pion reconstruction efficiency
                1: Kaon/pion selection efficiency
                2: Kaon/pion reconstruction tracking efficiency
                3: Stripping line correction
                4: Trigger efficiency
                5: Data/MC ratio correction
                10-15: Variations of the above corrections
        
        Returns:
            float: Efficiency correction weight (1/eff)
        """
        # Extract kinematic properties for kaon
        p_kaon = kaonVector.P() * 1000  # P is in GeV, Maps are in MeV
        pt_kaon = kaonVector.Pt() * 1000
        eta_kaon = kaonVector.PseudoRapidity()
        
        # Extract kinematic properties for pion
        p_pion = pionVector.P() * 1000
        pt_pion = pionVector.Pt() * 1000
        eta_pion = pionVector.PseudoRapidity()
        
        # Calculate D0 combined kinematics (for trigger)
        pt_combined = math.sqrt(pt_kaon * pt_pion)
        eta_d0 = (eta_kaon + eta_pion) / 2.0  # Approximate eta of D0
        
        # Initialize efficiency
        eff = 1.0
        
        # Apply appropriate correction based on type
        if type == -1:
            # Combined correction - not used, return placeholder
            eff = 1.0
            
        elif type == 0 or type == 10:
            # Kaon/pion reconstruction efficiency
            map_to_use = self.pionMapReco if type == 0 else self.pionMapReco_Var
            
            # Apply reconstruction efficiency for both kaon and pion
            # Using same map for K as for π is an approximation
            k_reco_eff = map_to_use.GetBinContent(map_to_use.FindBin(pt_kaon, eta_kaon))
            pi_reco_eff = map_to_use.GetBinContent(map_to_use.FindBin(pt_pion, eta_pion))
            eff = k_reco_eff * pi_reco_eff
            
        elif type == 1 or type == 11:
            # Kaon/pion selection efficiency
            map_to_use = self.pionMapSelHist if type == 1 else self.pionMapSelHist_Var
            
            # Apply selection efficiency for both kaon and pion
            # Using pion map for K is an approximation
            k_sel_eff = map_to_use.GetBinContent(map_to_use.FindBin(p_kaon, eta_kaon))
            pi_sel_eff = map_to_use.GetBinContent(map_to_use.FindBin(p_pion, eta_pion))
            eff = k_sel_eff * pi_sel_eff
            
        elif type == 2 or type == 12:
            # Tracking efficiency - using muon map as approximation
            map_to_use = self.MuonMap if type == 2 else self.MuonMap_Var
            
            # Apply tracking efficiency for both K and π
            k_track_eff = map_to_use.GetBinContent(map_to_use.FindBin(pt_kaon, eta_kaon))
            pi_track_eff = map_to_use.GetBinContent(map_to_use.FindBin(pt_pion, eta_pion))
            eff = k_track_eff * pi_track_eff
            
        elif type == 3 or type == 13:
            # Stripping line correction
            map_to_use = self.TriggerSelCorrMap if type == 3 else self.TriggerSelCorrMap_Var
            
            # Apply stripping line correction for both K and π
            k_strip_eff = map_to_use.GetBinContent(map_to_use.FindBin(eta_kaon, pt_kaon))
            pi_strip_eff = map_to_use.GetBinContent(map_to_use.FindBin(eta_pion, pt_pion))
            
            # Handle invalid efficiencies (>1)
            if k_strip_eff > 1: k_strip_eff = 0
            if pi_strip_eff > 1: pi_strip_eff = 0
            
            eff = k_strip_eff * pi_strip_eff
            
        elif type == 4 or type == 14:
            # Trigger efficiency
            func_list = self.functionList if type == 4 else self.functionList_Var
            
            # Get appropriate eta bin number
            if type == 4:
                eta_bin = self.TriggerEffMap2D.GetYaxis().FindBin(eta_d0)
            else:
                eta_bin = self.TriggerEffMap2D_Var.GetYaxis().FindBin(eta_d0)
                
            # Evaluate trigger efficiency
            if eta_bin > 0 and eta_bin <= len(func_list):
                trig_eff = func_list[eta_bin-1].Eval(pt_combined)
                eff = trig_eff
            else:
                # Handle out-of-range eta
                eff = 0
                
        elif type == 5 or type == 15:
            # Data/MC ratio correction
            corr_func = self.Ratio_D_MC.getCorr if type == 5 else self.Ratio_D_MC.getCorr_Var
            
            # Apply Data/MC correction for both K and π
            k_dm_corr = corr_func(eta=eta_kaon, pt=pt_kaon, p=p_kaon)
            pi_dm_corr = corr_func(eta=eta_pion, pt=pt_pion, p=p_pion)
            
            eff = k_dm_corr * pi_dm_corr
        
        # Return the inverse efficiency as the correction weight
        if eff > 0:
            eff_correction = 1.0 / eff
        else:
            eff_correction = 0
            
        return eff_correction
  
#- - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - -
#Fill unbiased part lvl tree
#def fillPartLvlTree(mcttree,origTree):
def fillPartLvlTree(mcttree,origTree, v_MCetaPart,v_MCphiPart,v_MCpTPart, v_MCnConstPart,v_MCetaTagPart,v_MCphiTagPart,v_MCTagPtPart,v_MCTagMPart,v_MCzTPart, v_MCtag_lifetPart, v_isPrimary, v_isDetTagRec):
    
    #for iTag in range(0, int(origTree.tag_pid.size())):
    for iTagGen in range(0, int(origTree.mctag_pid.size())):

        #idx_tagGen   = int(origTree.tag_idx_MCtag[iTagGen]) #This is the generator tag ID the Detector level tag is associated to   -> idx_tagGen = iTagGen
                
        idx_jetGen = int(origTree.mctag_idx_jet[iTagGen])    #This is the jetID the tag is associated to
        if idx_jetGen==-1: continue # no jet was reconstructed
        
        #-check if that tag has a corresponding tag reconstructed at detector level, regardless of jet situation
        idx_tagDet = int(origTree.mctag_idx_Dettag[iTagGen]) #This is the Det lvl tag the MC tag is associated to
        isDetRec=0
        #print("idx_tagDet: {}".format(idx_tagDet))
        if idx_tagDet>-1:
            isDetRec=1
            idx_MCTEST   = int(origTree.tag_idx_MCtag[idx_tagDet])
            if(iTagGen!=idx_MCTEST):
                print("ERROR: these two numbers should be the same!")
                print("iTagMC: {}".format(iTagGen))
                print("id ass. Det Tag: {}".format(idx_tagDet))
                print("id ass. MC Tag of the Det. Tag: {}".format(idx_MCTEST))
            
            
        isPrim         = origTree.mctag_isPrimary[iTagGen]
        p4PartTag      = getFourVector(origTree,"mctag",iTagGen)
        p4PartJet      = getFourVector(origTree,"mcjet",idx_jetGen)
        MCpvrNumber    = int(origTree.mctag_idx_pvr[iTagGen])
          
        MCpvrZ     = origTree.mcpvr_z[MCpvrNumber]  #-primary vertex of particle
        MCdecayZ   = origTree.mctag_z[iTagGen]      #-decay vertex of particle
        MCpz       = origTree.mctag_pz[iTagGen]     #-z-component of momentum vector

        lifeTime   = calcLifeTime(p4PartTag.M(),MCpvrZ,MCdecayZ,MCpz)

        #-Check if the particles at generator level are in the fiducial acceptance
        pionsAcc = 0
        muonsAcc = 0
        try:
            #-The two pions
            prt1Number     = int(origTree.mctag_idx_prt1[iTagGen])
            prt2Number     = int(origTree.mctag_idx_prt2[iTagGen])
            pid1 =origTree.mcprt_pid[prt1Number]
            pid2 =origTree.mcprt_pid[prt2Number]
            if math.fabs(pid1)==211 and math.fabs(pid2)==211:
                p4MCpion1 = getFourVector(origTree,"mcprt",prt1Number)
                p4MCpion2 = getFourVector(origTree,"mcprt",prt2Number)
                if p4MCpion1.Pt()>0.5 and p4MCpion2.Pt()>0.5 and p4MCpion1.P()>3 and p4MCpion2.P()>3:
                    pionsAcc = 1
        except:
            pionsAcc = -1

        try:
            #-The two muons
            prt3Number     = int(origTree.mctag_idx_prt3[iTagGen])
            prt4Number     = int(origTree.mctag_idx_prt4[iTagGen])

            pid3 =origTree.mcprt_pid[prt3Number]
            pid4 =origTree.mcprt_pid[prt4Number]
            #if pid0==443 and fabs(pid1)==211 and fabs(pid2)==211:
            if math.fabs(pid3)==13 and math.fabs(pid4)==13:
                p4MCmuon1 = getFourVector(origTree,"mcprt",prt3Number)
                p4MCmuon2 = getFourVector(origTree,"mcprt",prt4Number)
                if p4MCmuon1.Pt()>0.5 and p4MCmuon2.Pt()>0.5 and p4MCmuon1.P()>6 and p4MCmuon2.P()>6:
                    muonsAcc = 1
        except:
            muonsAcc = -1

        #-Do not accept events at generator level if the particles are not within the fiducial volume cuts
        if pionsAcc<1 or muonsAcc<1: continue


        #v_MCisPionAcc[0]  = pionsAcc
        #v_MCisMuonAcc[0]  = muonsAcc
        v_MCetaPart[0]    = p4PartJet.Eta()
        v_MCphiPart[0]    = p4PartJet.Phi()
        v_MCpTPart[0]     = p4PartJet.Pt()
        v_MCnConstPart[0] = origTree.mcjet_n_neu[iTagGen]+origTree.mcjet_n_chr[iTagGen]
        v_MCetaTagPart[0] = p4PartTag.Eta()
        v_MCphiTagPart[0] = p4PartTag.Phi()
        v_MCTagPtPart[0]  = p4PartTag.Pt()
        v_MCTagMPart[0]   = p4PartTag.M()
        if p4PartJet.Pt()>0:
            v_MCzTPart[0] = p4PartTag.Pt()/p4PartJet.Pt()
        else:
           v_MCzTPart[0]  = 0
        v_MCtag_lifetPart[0] = lifeTime
        v_isPrimary[0]    = isPrim
        v_isDetTagRec[0]  = isDetRec


        mcttree.Fill()
    return mcttree
    
    
    
#________________________________________________________________________
#________________________________________________________________________
def ntuplesJetAnalysis_Response_noVect(fFileName,printLvl,inputMC):
  
  newFilter = filterObject(fFileName,printLvl,inputMC)
  newFilter.filter()
  
#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description="- -")
 
  parser.add_argument("-f", "--fFileName", action="store",
                      type=str, metavar="fFileName",
                      default="AnalysisResults.root (-m)",
                      help="Path to root file containing the gh NTuple")
  parser.add_argument("-p", "--printLvl", action="store",
                      type=float, metavar="printLvl",
                      default="10",
                      help="Rate at which progress is printed default 10->every 10% of the final goal")
  parser.add_argument("-m", "--isMC", action="store",
                      type=int, metavar="isMC",
                      default="0",
                      help="is this for mc, yes=1, no=0")

  # Parse the arguments
  args = parser.parse_args()
  
  
  ntuplesJetAnalysis_Response_noVect(fFileName = args.fFileName, printLvl = args.printLvl, inputMC=args.isMC)

