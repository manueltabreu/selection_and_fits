import argparse
from ast import Yield
from calendar import c
parser = argparse.ArgumentParser(description="")
#parser.add_argument("inputfile" , help = "Path to the input ROOT file")
#parser.add_argument("dimusel"   , help = "Define if keep or remove dimuon resonances. You can choose: keepPsiP, keepJpsi, rejectPsi, keepPsi")
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
from ROOT import RooGaussian, RooExponential, RooChebychev, RooProdPdf, RooCBShape, TFile, TGraphErrors, TGraphMultiErrors, RooPolynomial, RooExtendPdf
import sys, math, pdb
#from uncertainties import ufloat
import array
import random
import numpy as np
import time

ROOT.RooMsgService.instance().setGlobalKillBelow(4)
ROOT.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(50000)
T = time.time()

dimusel = 'rejectPsi'

yields = [0, 0, 0, 0, 0, 0, 0, 0]
yields_error = [0, 0, 0, 0, 0, 0, 0, 0]

efficiency = [0, 0, 0, 0, 0, 0, 0, 0]
efficiency_error = [0, 0, 0, 0, 0, 0, 0, 0]

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

tgraph_x = [0, 0, 0, 0, 0, 0, 0, 0]
tgraph_ex = [0, 0, 0, 0, 0, 0, 0, 0]

d_branch_f = [0, 0, 0, 0, 0, 0, 0, 0]
d_branch_f_stat_error = [0, 0, 0, 0, 0, 0, 0, 0]
d_branch_f_syst_error = [0, 0, 0, 0, 0, 0, 0, 0]

branch_norm = 0.0593
branch_norm_error = 0.0006

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

#get the efficiencies and the yields from the 3 files

#names of the files, it must change with the option "dimusel"
f_NPsi_name = "fit_results_mass_checkOnMC/newbdt_puw/yields_efficiencies_%s%s.root"%(args.year, '_Jpsi'*(dimusel=='keepJpsi') + '_Psi'*(dimusel=='keepPsiP'))
graphic_yields_name = "tgraph_yields_%s%s"%(args.year, '_Jpsi'*(dimusel=='keepJpsi') + '_Psi'*(dimusel=='keepPsiP'))
graphic_efficiencies_name = "tgraph_efficiencies_%s%s"%(args.year, '_Jpsi'*(dimusel=='keepJpsi') + '_Psi'*(dimusel=='keepPsiP'))

dimusel = 'keepJpsi'
f_JPsi_name = "fit_results_mass_checkOnMC/newbdt_puw/yields_efficiencies_%s%s.root"%(args.year, '_Jpsi'*(dimusel=='keepJpsi') + '_Psi'*(dimusel=='keepPsiP'))
graphic_yields_Jpsi_name = "tgraph_yields_%s%s"%(args.year, '_Jpsi'*(dimusel=='keepJpsi') + '_Psi'*(dimusel=='keepPsiP'))
graphic_efficiencies_Jpsi_name = "tgraph_efficiencies_%s%s"%(args.year, '_Jpsi'*(dimusel=='keepJpsi') + '_Psi'*(dimusel=='keepPsiP'))

dimusel = 'keepPsiP'
f_PsiP_name = "fit_results_mass_checkOnMC/newbdt_puw/yields_efficiencies_%s%s.root"%(args.year, '_Jpsi'*(dimusel=='keepJpsi') + '_Psi'*(dimusel=='keepPsiP'))
graphic_yields_PsiP_name = "tgraph_yields_%s%s"%(args.year, '_Jpsi'*(dimusel=='keepJpsi') + '_Psi'*(dimusel=='keepPsiP'))
graphic_efficiencies_PsiP_name = "tgraph_efficiencies_%s%s"%(args.year, '_Jpsi'*(dimusel=='keepJpsi') + '_Psi'*(dimusel=='keepPsiP'))

#name of the final file
bf_file_name = "fit_results_mass_checkOnMC/newbdt_puw/branching_fraction_%s.root"%(args.year)
bf_file = TFile(bf_file_name, "RECREATE")

#For the non-resonant regions
dimusel = 'rejectPsi'

f_NPsi = TFile(f_NPsi_name, "READ")

graph_Y = f_NPsi.Get(graphic_yields_name)
graph_e = f_NPsi.Get(graphic_efficiencies_name)

for ibin in range(len(q2binning)-1): 
    yields[ibin] = graph_Y.GetPointY(ibin)
    yields_error[ibin] = graph_Y.GetErrorY(ibin)

    efficiency[ibin] = graph_e.GetPointY(ibin)
    efficiency_error[ibin] = graph_e_JPsi.GetErrorY(ibin)
    efficiency[ibin] = efficiency[ibin] / num_real_events[ibin]

    tgraph_x[ibin] = graph_Y.GetPointX(ibin)
    tgraph_ex[ibin] = graph_Y.GetErrorX(ibin)

#For the Jpsi
dimusel = 'keepJpsi'
print (f_JPsi_name)
f_JPsi = TFile(f_JPsi_name, "READ")

graph_Y_JPsi = f_JPsi.Get(graphic_yields_Jpsi_name)
graph_e_JPsi = f_JPsi.Get(graphic_efficiencies_Jpsi_name)

# in case we still didn't run for the jpsi, 1,92E+06 for the numerator

yields[4] = graph_Y_JPsi.GetPointY(4)
yields_error[4] = graph_Y_JPsi.GetErrorY(4)

efficiency[4] = graph_e_JPsi.GetPointY(4)
efficiency_error[4] = graph_e_JPsi.GetErrorY(4)
#efficiency[4] = 1920000
efficiency[4] = efficiency[4] / num_real_events[4]

tgraph_x[4] = graph_Y_JPsi.GetPointX(0)
tgraph_ex[4] = graph_Y_JPsi.GetErrorX(0)

#For the Psi Prime
dimusel = 'keepPsiP'
f_PsiP = TFile(f_PsiP_name, "READ")

graph_Y_PsiP = f_PsiP.Get(graphic_yields_PsiP_name)
graph_e_PsiP = f_PsiP.Get(graphic_efficiencies_PsiP_name)

# in case we still didn't run for the PsiP, 226602 for the numerator
yields[6] = graph_Y_PsiP.GetPointY(6)
yields_error[6] = graph_Y_PsiP.GetErrorY(6)

efficiency[6] = graph_e_PsiP.GetPointY(6)
efficiency_error[6] = graph_e_PsiP.GetErrorY(6)
#efficiency[6] = 226602
efficiency[6] = efficiency[6] / num_real_events[6]

tgraph_x[6] = graph_Y_PsiP.GetPointX(6)
tgraph_ex[6] = graph_Y_PsiP.GetErrorX(6)

#protecting from errors (eg.: divided by zero)
print (yields[4])
if yields[4] == 0.:
    print 'no value obtained for normalized yield'
    print 'assuming value of 8801283'
    yields[4] = 3000

if yields[6] == 0.:
    print 'no value obtained for the yield of the Psi Prime'
    print 'assuming value of 972182'
    yields[6] = 2000 

if efficiency[4] == 0.:
    print 'no value calculated for normalized efficiency'
    print 'assuming value of 0.00247'
    efficiency[4] = 0.00247

if efficiency[6] == 0.:
    print 'no value calculated for efficiency of the Psi Prime'
    print 'assuming value of 0.00265'
    efficiency[6] = 0.00265   

print 'yields'
print(yields)
for ibin in range(len(q2binning)-1):
    print(yields[ibin])
print 'efficiency'
print(efficiency)

#get the branching fraction
for ibin in range(len(q2binning)-1):
    delta_qi2 = (q2binning[ibin + 1] - q2binning[ibin]) * (q2binning[ibin + 1] - q2binning[ibin])
    print 'delta_qi2 = ', delta_qi2
    print 'yield[',ibin,'] = ', yields[ibin]
    print 'efficiency[',ibin,'] = ', efficiency[ibin]
    print 'yield[',4,'] = ', yields[4]
    print 'efficiency[',4,'] = ', efficiency[4]
    print 'd_branch_f[',4,'] = ',branch_norm

    d_branch_f[ibin] = 100 * yields[ibin] * efficiency[4] *branch_norm / ( yields[4] * efficiency[ibin] * delta_qi2)
    print 'd_branch_f[',ibin,'] = ', d_branch_f[ibin]
    print '----------------'

    fst_stat_term = yields_error[ibin] / yields[4]
    snd_stat_term = yields_error[4] * yields[ibin] / ( yields[4] * yields[4] )
    squared_sum = fst_stat_term * fst_stat_term + snd_stat_term * snd_stat_term
    d_branch_f_stat_error[ibin] = math.sqrt(squared_sum) 

    fst_syst_term = efficiency_error[ibin] / efficiency[4]
    snd_syst_term = efficiency_error[4] * efficiency[ibin] / ( efficiency[4] * efficiency[4] )
    squared_sum = fst_syst_term * fst_syst_term + snd_syst_term * snd_syst_term
    fst_syst_rel_error = math.sqrt(squared_sum) * efficiency[ibin] / ( branch_norm * efficiency[4] )
    snd_syst_rel_error = branch_norm_error * efficiency[ibin] / ( branch_norm * efficiency[4] )
    squared_sum = fst_syst_rel_error * fst_syst_rel_error + snd_syst_rel_error * snd_syst_rel_error
    d_branch_f_syst_rel_error = math.sqrt(squared_sum) 
    d_branch_f_syst_error[ibin] = d_branch_f_syst_rel_error * branch_norm * efficiency[4] / efficiency[ibin]
# d_branch_f[ibin] = d_branch_f[ibin] * 1000

print 'd_branch_f'
print (d_branch_f)
#d_branch_f[4] =branch_norm/1000
#d_branch_f[5] = d_branch_f[5]/1000
#d_branch_f[6] = d_branch_f[6]/1000000
#print (d_branch_f)

arr_d_branch_f = array.array('d',d_branch_f)
arr_d_branch_f_stat_error = array.array('d',d_branch_f_stat_error)
arr_d_branch_f_syst_error = array.array('d',d_branch_f_syst_error)
arr_tgraph_x = array.array('d',tgraph_x)
arr_tgraph_ex = array.array('d',tgraph_ex)

c1 = ROOT.TCanvas()
tgraph_bf = TGraphMultiErrors( "tgraph_bf", "Branching Fraction", 8, np.array(arr_tgraph_x), np.array(arr_d_branch_f), np.array(arr_tgraph_ex), np.array(arr_tgraph_ex), np.array(arr_d_branch_f_stat_error), np.array(arr_d_branch_f_stat_error))
tgraph_bf.AddYError(8, np.array(arr_d_branch_f_syst_error), np.array(arr_d_branch_f_syst_error))
tgraph_bf.SetMarkerStyle(20)
tgraph_bf.GetAttLine(0).SetLineColor(kRed)
tgraph_bf.GetAttLine(1).SetLineColor(kBlue)


tgraph_bf.Draw("AP")
c1.SaveAs("branching_fraction.pdf")

tgraph_bf.Write("tgraph_branching_fraction%s"%(args.year))

print 'It took', time.time() - T, 'seconds'
print (time.time()- T)/60, 'minutes'
print (time.time() - T)/3600, 'hours'

f_NPsi.Close()
f_JPsi.Close()
f_PsiP.Close()
bf_file.Close()