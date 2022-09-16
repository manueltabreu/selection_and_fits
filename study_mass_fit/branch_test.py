import argparse
from ast import Yield
from calendar import c
from tkinter import N
parser = argparse.ArgumentParser(description="")
#parser.add_argument("inputfile" , help = "Path to the input ROOT file")
parser.add_argument("dimusel"   , help = "Define if keep or remove dimuon resonances. You can choose: keepPsiP, keepJpsi, rejectPsi, keepPsi")
parser.add_argument("year"      , help = "choose among:2016,2017,%s", default = '2018')
args = parser.parse_args()


'''
code to get the branching fraction:
'''

import os, sys, inspect
from os import path
sys.path.insert(0, os.environ['HOME'] + '/.local/lib/python2.7/site-packages')

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import ROOT
from ROOT import gSystem
ROOT.gROOT.SetBatch(True)

gSystem.Load('libRooFit')
gSystem.Load('../utils/func_roofit/libRooDoubleCBFast')
gSystem.Load('../utils/func_roofit/libRooGaussDoubleSidedExp')
from ROOT import RooFit, RooRealVar, RooAbsReal, RooDataSet, RooArgList, RooTreeData, RooArgSet, RooAddPdf, RooFormulaVar
from ROOT import RooGaussian, RooExponential, RooChebychev, RooProdPdf, RooCBShape, TFile, TGraphErrors, RooPolynomial, RooExtendPdf
import sys, math, pdb
#from uncertainties import ufloat
import array
import random
import numpy as np
import time

ROOT.RooMsgService.instance().setGlobalKillBelow(4)
ROOT.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(50000)
T = time.time()


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
#                 19,
]

num_real_events = [
                        60600671,
                        113787898,
                        86443130,
                        157793256,
                        776202791,
                        185031164,
                        85358493,
                        110280298,
]

MC_num_entries = [104824.0, 239995.0, 201743.0, 360605.0, 0, 631800.0, 0, 390447.0]
MC_RT_num_entries = [0, 0, 0, 0, 0, 0, 0, 0]
MC_WT_num_entries = [0, 0, 0, 0, 0, 0, 0, 0]
MC_num_entries_error = [0, 0, 0, 0, 0, 0, 0, 0]
eff = [0, 0, 0, 0, 0, 0, 0, 0]
eff_error = [0, 0, 0, 0, 0, 0, 0, 0]
tgraph_y = [266.24961202947884, 457.4291188016444, 360.8221403914036, 792.7427824641064, 0, 1200.496289495756, 0, 1.3179568547627696e-05]
tgraph_x = [0, 0, 0, 0, 0, 0, 0, 0]
tgraph_ex = [0, 0, 0, 0, 0, 0, 0, 0]
tgraph_ey = [19, 26, 23, 33.5, 0, 42, 0, 0.5]


#create Tgraph
for ibin in range(len(q2binning)-1):

#   yields
    middle_value = (q2binning[ibin+1] +  q2binning[ibin]) / 2

#   efficiencies
#    MC_num_entries[ibin] = MC_RT_num_entries[ibin] + MC_WT_num_entries[ibin]
    eff[ibin] = MC_num_entries[ibin] / num_real_events[ibin]
    MC_num_entries_error[ibin] = math.sqrt(MC_num_entries[ibin])
    eff_error[ibin] = MC_num_entries_error[ibin] / num_real_events[ibin]

#   q2bins
    tgraph_x[ibin] = middle_value
    tgraph_ex[ibin] = middle_value - q2binning[ibin]


arr_tgraph_y = array.array('d',tgraph_y)
arr_tgraph_x = array.array('d',tgraph_x)
arr_tgraph_ex = array.array('d',tgraph_ex)
arr_tgraph_y = array.array('d',tgraph_y)
arr_tgraph_ey = array.array('d',tgraph_ey)
arr_eff = array.array('d',eff)
arr_eff_error = array.array('d',eff_error)

tgraph_yield = TGraphErrors(8, np.array(arr_tgraph_x), np.array(arr_tgraph_y), np.array(arr_tgraph_ex), np.array(arr_tgraph_ey))
#tgraph_yield = TGraphErrors(8, tgraph_x, tgraph_y, tgraph_ex, tgraph_ey)

tgraph_efficiency = TGraphErrors(8, np.array(arr_tgraph_x), np.array(arr_eff), np.array(arr_tgraph_ex), np.array(arr_eff_error))

print 'tgraph_y'
print(tgraph_y)

print 'eff'
print(eff)

print 'eff_error'
print(eff_error)

c2 = ROOT.TCanvas() 
tgraph_yield.SetName("Yields")
tgraph_yield.SetTitle("Yields")
tgraph_yield.Draw("AP")
c2.SaveAs("tgraph_yields1.pdf")

print 'tgraph_yield created'
print (tgraph_yield)

c3 = ROOT.TCanvas() 
tgraph_efficiency.SetName("Efficiencies_num")
tgraph_efficiency.SetTitle("Efficiencies_num")
tgraph_efficiency.Draw("AP")
c3.SaveAs("tgraph_efficiencies.pdf")

print 'tgraph_efficiency created'
print (tgraph_efficiency)

out_gname = "fit_results_mass_checkOnMC/newbdt_puw/yields_efficiencies_%s%s.root"%(args.year, '_Jpsi'*(args.dimusel=='keepJpsi') + '_Psi'*(args.dimusel=='keepPsiP'))
out_g = TFile (out_gname, "RECREATE") 
tgraph_yield.Write("tgraph_yields_%s%s"%(args.year, '_Jpsi'*(args.dimusel=='keepJpsi') + '_Psi'*(args.dimusel=='keepPsiP')))
#out_g.WriteObject(tgraph_yield, "Yields")
tgraph_efficiency.Write("tgraph_efficiencies_%s%s"%(args.year, '_Jpsi'*(args.dimusel=='keepJpsi') + '_Psi'*(args.dimusel=='keepPsiP')))
#out_g.WriteObject(tgraph_efficiency, "Efficiencies_num")

print 'It took', time.time() - T, 'seconds'
print (time.time()- T)/60, 'minutes'
print (time.time() - T)/3600, 'hours'

out_g.Close()
