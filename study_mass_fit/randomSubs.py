import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("dimusel"   , help = "Define if keep or remove dimuon resonances. You can choose: keepPsiP, keepJpsi, rejectPsi, keepPsi")
parser.add_argument("year"      , help = "choose among:2016,2017,2018", default = '2018')
parser.add_argument('--ntoys', type=int, default = 1)
args = parser.parse_args()

import os, sys, inspect
from os import path
sys.path.insert(0, os.environ['HOME'] + '/.local/lib/python2.7/site-packages')
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import ROOT
from ROOT import gSystem, gStyle
ROOT.gROOT.SetBatch(True)

gStyle.SetOptFit(True) 
gSystem.Load('libRooFit')
from ROOT import RooFit, RooRealVar, RooDataSet, RooArgList, RooArgSet, RooAddPdf, RooFormulaVar, TFile
import sys, math, pdb, random
import numpy as np

sys.path.append("../utils")
from utils.utils import *
nSigma_psiRej = 3.

tData = ROOT.TChain('ntuple')
# tData.Add('/gwpool/users/fiorendi/p5prime/miniAOD/CMSSW_10_2_14/src/miniB0KstarMuMu/miniKstarMuMu/bdt/final_ntuples/2017MC_LMNR_100k.root')
tData.Add('/gwteray/users/fiorendi/final_ntuples_p5prime_allyears/%sMC_LMNR.root'%args.year)

## numbers for 2018, will be read from common function in the next iteration
n_bin = n_data[args.year]
q2binning = [
                1,
                2, 
                4.3,
                6,
                8.68,
                10.09,
                12.86,
                14.18,
                16,
]

cut_base      = applyB0PsiCut(args.dimusel, nSigma_psiRej)



def findFrt(fulldata, ibin, h_truef):  

    cut  = cut_base + '&& (mumuMass*mumuMass > %s && mumuMass*mumuMass < %s)'%(q2binning[ibin], q2binning[ibin+1])
    fulldata_v2 = fulldata.reduce(RooArgSet(tagged_mass,mumuMass,mumuMassE), cut)

    nDataEntries = fulldata_v2.sumEntries()
    nDesired = n_bin[ibin]/nDataEntries

    for it in range(args.ntoys):

        print it
        htag = ROOT.TH1F('htag', 'htag', 4,0,4)
        htag.Reset()
        tData.Draw("tagB0+genSignal >> htag", "rndm() < %s"%nDesired, "goff")##; // draw random sample of about 20% of entries

        n_rt = htag.GetBinContent(3)
        n_wt = htag.GetBinContent(2) + htag.GetBinContent(4)

        h_truef.Fill(float(n_rt)/(n_rt + n_wt))
        
    

tagged_mass     = RooRealVar("tagged_mass" , "#mu^{+}#mu^{-}K#pi mass", 5., 5.6, "GeV")
mumuMass        = RooRealVar("mumuMass"    , "mumuMass" , 0, 6);
mumuMassE       = RooRealVar("mumuMassE"   , "mumuMassE", 0, 10000);

tagged_mass.setRange("full",   5.0,5.6) ;
thevars = RooArgSet()
thevars.add(tagged_mass)
thevars.add(mumuMass)
thevars.add(mumuMassE)

fulldata   = RooDataSet('fulldata', 'fulldataset', tData,  RooArgSet(thevars))
## add to the input tree the combination of the variables, to be used for the cuts on the dimuon mass
deltaB0Mfunc = RooFormulaVar("deltaB0M", "deltaB0M", "@0 - @1", RooArgList(tagged_mass,B0Mass) )
deltaJMfunc  = RooFormulaVar("deltaJpsiM" , "deltaJpsiM" , "@0 - @1", RooArgList(mumuMass,JPsiMass) )
deltaPMfunc  = RooFormulaVar("deltaPsiPM" , "deltaPsiPM" , "@0 - @1", RooArgList(mumuMass,PsiPMass) )
deltaB0M     = fulldata.addColumn(deltaB0Mfunc) ;
deltaJpsiM   = fulldata.addColumn(deltaJMfunc) ;
deltaPsiPM   = fulldata.addColumn(deltaPMfunc) ;
thevars.add(deltaB0M)
thevars.add(deltaJpsiM)
thevars.add(deltaPsiPM)


out_f = TFile ("checkfrt%s.root"%args.year,"RECREATE") 
out_w = ROOT.RooWorkspace("toy_w")

for ibin in range(len(q2binning)-1):

    h_truef      = ROOT.TH1F('bin%s'%ibin,    'bin%s;frt_{real}; '%ibin, 250, 0.82,0.92)
    if args.dimusel == 'rejectPsi' and \
       (q2binning[ibin] == 8.68 or q2binning[ibin] == 12.86): 
           continue

    print ibin
    findFrt(fulldata, ibin, h_truef)

    c2 = ROOT.TCanvas('c2','c2',400,400)
    h_truef.Draw()
    h_truef.Fit('gaus')
    c2.SaveAs('frt_toy_bin%s_%s_%s.pdf' %(ibin,args.ntoys,args.year))
    out_f.cd()
    h_truef.Write()

out_f.Close()
