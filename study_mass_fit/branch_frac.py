import argparse
from ast import Yield
from calendar import c
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

yields = [0, 0, 0, 0, 0, 0, 0, 0]
yields_error = [0, 0, 0, 0, 0, 0, 0, 0]

efficiency = [0, 0, 0, 0, 0, 0, 0, 0]
efficiency_error = [0, 0, 0, 0, 0, 0, 0, 0]

num_real_events = [158292, 328218, 300020, 409302, 837837, 725094, 7273466, 453982]

tgraph_x = [0, 0, 0, 0, 0, 0, 0, 0]
tgraph_ex = [0, 0, 0, 0, 0, 0, 0, 0]

branch_f = [0, 0, 0, 0, 1.46, 0, 0, 0]
branch_f_error = [0, 0, 0, 0, 0, 0, 0, 0]

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

#get the efficiencies and the yields
#t_MC_Data = ROOT.TChain('ntuple')
#t_MC_Data.Add("fit_results_mass_checkOnMC/newbdt_puw/yields_%s%s.root"%(args.year, '_Jpsi'*(args.dimusel=='keepJpsi') + '_Psi'*(args.dimusel=='keepPsiP')))
#data = RooDataSet("data", "yields_and_efficiencies",t_MC_Data,RooArgSet(Yields))


#print 'mc file name:',  t_MC_Data.GetFile().GetName()
fname = "fit_results_mass_checkOnMC/newbdt_puw/yields_efficiencies_%s%s.root"%(args.year, '_Jpsi'*(args.dimusel=='keepJpsi') + '_Psi'*(args.dimusel=='keepPsiP'))
f = TFile(fname, "READ")
#f = TFile.Open("fit_results_mass_checkOnMC/newbdt_puw/yields_2018_Psi.root") 

bf_file_name = "fit_results_mass_checkOnMC/newbdt_puw/branching_fraction_%s%s.root"%(args.year, '_Jpsi'*(args.dimusel=='keepJpsi') + '_Psi'*(args.dimusel=='keepPsiP'))
bf_file = TFile(bf_file_name, "RECREATE")

print 'the objects are'
f.Print()
print 'size'
#print(f.GetSize())
#print (f.Sizeof())
#graph_Y = TGraphErrors()
#f.GetObject("Yields", graph_Y)
#graph_Y = f.Yields
graph_Y = f.Get("tgraph_yields_%s%s"%(args.year, '_Jpsi'*(args.dimusel=='keepJpsi') + '_Psi'*(args.dimusel=='keepPsiP')))
graph_e = f.Get("tgraph_efficiencies_%s%s"%(args.year, '_Jpsi'*(args.dimusel=='keepJpsi') + '_Psi'*(args.dimusel=='keepPsiP')))

#graph_Y = TGraphErrors("fit_results_mass_checkOnMC/newbdt_puw/yields_%s%s.root"%(args.year, '_Jpsi'*(args.dimusel=='keepJpsi') + '_Psi'*(args.dimusel=='keepPsiP')), "%lg %lg %lg %lg", "")

print(type(graph_Y))
print(graph_Y)
print 'here'
#print(graph_Y.GetEY()[3])
print(graph_Y.GetPointY(3))
print 'errors'
print(graph_Y.GetErrorX(2))
print(graph_Y.GetErrorX(3))
print 'yield'
yields = graph_Y.GetEY()

print 'FINALLY'
for ibin in range(len(q2binning)-1): 
    yields[ibin] = graph_Y.GetPointY(ibin)
    efficiency[ibin] = graph_e.GetPointY(ibin)
    efficiency[ibin] = efficiency[ibin] / num_real_events[ibin]
    tgraph_x[ibin] = graph_Y.GetPointX(ibin)
    tgraph_ex = graph_Y.GetErrorX(ibin)

efficiency[4] = 0.9
efficiency[6] = 0.9
yields[4] = 883249
yields[6] = 772182
print 'yields'
print(yields)
print 'what'
print(yields[3])
#for ibin in range(len(q2binning)-1):
#    print(yields[ibin])
print 'efficiency'
print(efficiency)
#efficiency

#get the branching fraction
for ibin in range(len(q2binning)-1):
   delta_qi2 = (q2binning[ibin + 1] - q2binning[ibin]) * (q2binning[ibin + 1] - q2binning[ibin])
   branch_f[ibin] = yields[ibin] * efficiency[4] * branch_f[4] / ( yields[4] * efficiency[ibin] * delta_qi2)
   branch_f[ibin] = branch_f[ibin] * 1000

print 'branch_f'
print (branch_f)

arr_branch_f = array.array('d',branch_f)
arr_branch_f_error = array.array('d',branch_f_error)
arr_tgraph_x = array.array('d',tgraph_x)
arr_tgraph_ex = array.array('d',tgraph_ex)

tgraph_bf = TGraphErrors(8, np.array(arr_tgraph_x), np.array(arr_branch_f), np.array(arr_tgraph_ex), np.array(arr_branch_f_error))
c1 = ROOT.TCanvas()
tgraph_bf.Draw("AP")
c1.SaveAs("branching_fraction.pdf")

tgraph_bf.Write("tgraph_branching_fraction%s%s"%(args.year, '_Jpsi'*(args.dimusel=='keepJpsi') + '_Psi'*(args.dimusel=='keepPsiP')))

print 'It took', time.time() - T, 'seconds'
print (time.time()- T)/60, 'minutes'
print (time.time() - T)/3600, 'hours'

f.Close()
bf_file.Close()