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
from ROOT import gStyle
from ROOT import gROOT
from ROOT import TStyle

gSystem.Load('libRooFit')
gSystem.Load('../utils/func_roofit/libRooDoubleCBFast')
gSystem.Load('../utils/func_roofit/libRooGaussDoubleSidedExp')
from ROOT import RooFit, RooRealVar, RooAbsReal, RooDataSet, RooArgList, RooTreeData, RooArgSet, RooAddPdf, RooFormulaVar
from ROOT import RooGaussian, RooExponential, RooChebychev, RooProdPdf, RooCBShape, TFile, TGraphPainter, TGraphErrors, TGraphMultiErrors, RooPolynomial, RooExtendPdf
from ROOT import TLegend
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
corrected_yields = [0, 0, 0, 0, 0, 0, 0]

num_events = [0, 0, 0, 0, 1, 0, 1, 0]
efficiency = [0, 0, 0, 0, 0, 0, 0, 0]
efficiency_error = [0., 0., 0., 0., 0., 0., 0., 0.]
efficiency_error_squared = [0., 0., 0., 0., 0., 0., 0., 0.]
eff_uncert_term = [0, 0, 0, 0, 0, 0, 0, 0]

num_real_events = [
                        60600671,
                        113787898,
                        86443130,
                        157793256,
                        #776202791,
                        517468528,
                        185031164,
                        #85358493,
                        56905662,
                        110280298,
]

tgraph_x = [0, 0, 0, 0, 0, 0, 0, 0]
tgraph_ex = [0, 0, 0, 0, 0, 0, 0, 0]
tgraph_temp = [0, 0, 0, 0, 0, 0, 0, 0]
tgraph_corrected_x = [0, 0, 0, 0, 0, 0, 0]

d_branch_f = [0, 0, 0, 0, 0, 0, 0, 0]
yields_uncert_term = [0, 0, 0, 0, 0, 0, 0, 0]
d_branch_f_stat_error = [0, 0, 0, 0, 0, 0, 0, 0]

bf_uncert_term = 0.

d_branch_f_tot_rel_error = [0, 0, 0, 0, 0, 0, 0, 0]
d_branch_f_tot_error = [0, 0, 0, 0, 0, 0, 0, 0]

branch_norm1 = 0.0593
branch_norm2 = 0.00127
branch_norm = 0.0593 * 0.00127
branch_norm_error = 0
branch_norm_1st_error = 0.0006
branch_norm_2nd_error = 0.00005

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

f_NPsi = TFile(f_NPsi_name, "READ")

graph_Y = f_NPsi.Get(graphic_yields_name)
graph_e = f_NPsi.Get(graphic_efficiencies_name)

for ibin in range(len(q2binning)-1): 
    yields[ibin] = graph_Y.GetPointY(ibin)
    yields_error[ibin] = graph_Y.GetErrorY(ibin)

    num_events[ibin] = graph_e.GetPointY(ibin)
#    efficiency_error_squared[ibin] = graph_e.GetErrorY(ibin)
#   to avoid division by zero 
    num_events[4] = 1
    num_events[6] = 1
    efficiency[ibin] = num_events[ibin] / num_real_events[ibin]
    efficiency_error_squared[ibin] = 1 / num_events[ibin]

    tgraph_x[ibin] = graph_Y.GetPointX(ibin)
    tgraph_ex[ibin] = graph_Y.GetErrorX(ibin)

#For the Jpsi

print (f_JPsi_name)
f_JPsi = TFile(f_JPsi_name, "READ")

graph_Y_JPsi = f_JPsi.Get(graphic_yields_Jpsi_name)
graph_e_JPsi = f_JPsi.Get(graphic_efficiencies_Jpsi_name)

# in case we still didn't run for the jpsi, 1,92E+06 for the numerator

yields[4] = graph_Y_JPsi.GetPointY(4)
yields_error[4] = graph_Y_JPsi.GetErrorY(4)

num_events[4] = graph_e_JPsi.GetPointY(4)
#efficiency_error_squared[4] = graph_e_JPsi.GetErrorY(4)
print 'eff'
print (num_events[4])
print (graph_e_JPsi.GetPointY(4))
#efficiency[4] = 1920000
efficiency[4] = num_events[4] / num_real_events[4]
efficiency_error_squared[4] = 1 / num_events[4]

tgraph_x[4] = graph_Y_JPsi.GetPointX(4)
tgraph_ex[4] = graph_Y_JPsi.GetErrorX(4)

#For the Psi Prime
dimusel = 'keepPsiP'
f_PsiP = TFile(f_PsiP_name, "READ")

graph_Y_PsiP = f_PsiP.Get(graphic_yields_PsiP_name)
graph_e_PsiP = f_PsiP.Get(graphic_efficiencies_PsiP_name)

# in case we still didn't run for the PsiP, 226602 for the numerator
yields[6] = graph_Y_PsiP.GetPointY(6)
yields_error[6] = graph_Y_PsiP.GetErrorY(6)

num_events[6] = graph_e_PsiP.GetPointY(6)
#efficiency[6] = 226602
efficiency[6] = num_events[6] / num_real_events[6]
efficiency_error_squared[6] = 1 / num_events[6]

tgraph_x[6] = graph_Y_PsiP.GetPointX(6)
tgraph_ex[6] = graph_Y_PsiP.GetErrorX(6)

#protecting from errors (eg.: divided by zero)
print (yields[4])
if yields[4] == 0.:
    print 'no value obtained for normalized yield'
    print 'assuming value of 6264'
    yields[4] = 6264
    yields_error[4] = 87

if yields[6] == 0.:
    print 'no value obtained for the yield of the Psi Prime'
    print 'assuming value of 429'
    yields[6] = 429 
    yields_error[6] = 23

if efficiency[4] == 0.:
    print 'no value calculated for normalized efficiency'
    print 'assuming value of 0.00247'
    efficiency[4] = 0.00247
    efficiency_error_squared[4] = 1 / num_events[4]

if efficiency[6] == 0.:
    print 'no value calculated for efficiency of the Psi Prime'
    print 'assuming value of 0.00265'
    efficiency[6] = 0.00265   
    efficiency_error_squared[6] = math.sqrt(226602) / num_real_events[6]

print 'yields'
print(yields)
for ibin in range(len(q2binning)-1):
    print(yields[ibin])
print 'efficiency'
print(efficiency)

#get the branching fraction
for ibin in range(len(q2binning)-1):
    delta_qi2 = (q2binning[ibin + 1] - q2binning[ibin]) 
    print 'delta_qi2 = ', delta_qi2
    print 'yield[',ibin,'] = ', yields[ibin]
    print 'efficiency[',ibin,'] = ', efficiency[ibin]
    print 'yield[',4,'] = ', yields[4]
    print 'efficiency[',4,'] = ', efficiency[4]
    print 'd_branch_f[',4,'] = ',branch_norm

    d_branch_f[ibin] = yields[ibin] * efficiency[4] * branch_norm / ( yields[4] * efficiency[ibin] * delta_qi2)
    print 'd_branch_f[',ibin,'] = ', d_branch_f[ibin]
    print '----------------'

    efficiency_error[ibin] = math.sqrt(efficiency_error_squared[ibin]) * efficiency[ibin]
    print 'eff_error'
    print efficiency_error[ibin]

    print 'yields_error[ibin]'
    print(yields_error[ibin])
    print 'yields_error[4]'
    print(yields_error[4])
    fst_stat_term = yields_error[ibin] / yields[ibin]
    snd_stat_term = yields_error[4] / yields[4] 
    yields_uncert_term[ibin] = fst_stat_term * fst_stat_term + snd_stat_term * snd_stat_term
    d_branch_f_stat_error[ibin] = math.sqrt( yields_uncert_term[ibin] ) * d_branch_f[ibin]

    eff_uncert_term[ibin] = efficiency_error_squared[ibin] + efficiency_error_squared[4]

    print '1st_norm_err'
    print branch_norm_1st_error
    print (branch_norm_1st_error * branch_norm_1st_error)
    print 'branch_norm'
    print branch_norm
    print (branch_norm * branch_norm)
    branch_norm_rel_1st_error_squared = (branch_norm_1st_error * branch_norm_1st_error) / (branch_norm1 * branch_norm1)
    branch_norm_rel_2nd_error_squared = (branch_norm_2nd_error * branch_norm_2nd_error) / (branch_norm2 * branch_norm2)

    bf_uncert_term = branch_norm_rel_1st_error_squared + branch_norm_rel_2nd_error_squared

    branch_norm_error = math.sqrt(bf_uncert_term) * branch_norm
    print '1st err squared'
    print branch_norm_rel_1st_error_squared
    print '2nd err squared'
    print branch_norm_rel_2nd_error_squared
    print 'bf_uncert_term'
    print bf_uncert_term
    d_branch_f_tot_rel_error[ibin] = yields_uncert_term[ibin] + eff_uncert_term[ibin] + bf_uncert_term
    d_branch_f_tot_error[ibin] = math.sqrt( d_branch_f_tot_rel_error[ibin] ) * d_branch_f[ibin]
    #    fst_syst_rel_error = math.sqrt(squared_sum) 
    #    snd_syst_rel_error = branch_norm_error * efficiency[ibin] / ( branch_norm * efficiency[4] )
    #   print 'fst_syst_rel_error'
    #   print(fst_syst_rel_error)
    #    squared_sum = fst_syst_rel_error * fst_syst_rel_error + snd_syst_rel_error * snd_syst_rel_error
    #    d_branch_f_syst_rel_error = math.sqrt(squared_sum) 
    #    d_branch_f_syst_error[ibin] = d_branch_f_syst_rel_error * branch_norm * efficiency[4] / efficiency[ibin]
    # d_branch_f[ibin] = d_branch_f[ibin] * 1000

#temp4_ex_err = tgraph_ex[4]
#tgraph_ex[4] = 0
#temp6_ex_err = tgraph_ex[6]
#tgraph_ex[6] = 0
d_branch_f[4] = -(1E-9)
d_branch_f[6] = -(1E-9)
d_branch_f[7] = -(1E-9)
#d_branch_f[7] = 0
#d_branch_f_stat_error[4] = 0
#d_branch_f_tot_error[4] = 0
#d_branch_f_stat_error[6] = 0
#d_branch_f_tot_error[6] = 0

print (tgraph_x)
print 'x error'
print (tgraph_ex)
print 'd_branch_f'
print (d_branch_f)
#d_branch_f[4] =branch_norm/1000
#d_branch_f[5] = d_branch_f[5]/1000
#d_branch_f[6] = d_branch_f[6]/1000000
#print (d_branch_f)

print 'stat_errors'
print (d_branch_f_stat_error)
print 'tot_errors'
print (d_branch_f_tot_error)

print 'num_events'
print num_events

for i in range(len(tgraph_corrected_x) - 1):
    if i < 4:
        corrected_yields[i] = yields[i]
        tgraph_corrected_x[i] = tgraph_x[i]
    if i > 4:
        corrected_yields[i] = yields[i+1]
        tgraph_corrected_x[i] = tgraph_x[i+1]

print 'NORM_BRANCH'
print branch_norm
print 'NORM_BRANCH_ERROR'
print branch_norm_error

#temp7_x = tgraph_x[7]
#temp6_x = tgraph_x[6]
#temp4_x = tgraph_x[4]
#tgraph_x[7] = tgraph_x[0]
#tgraph_x[6] = tgraph_x[0]
#tgraph_x[4] = tgraph_x[0]
#temp7_ex_err = tgraph_ex[7]
#tgraph_ex[7] = 0
arr_d_branch_f = array.array('d',d_branch_f)
arr_d_branch_f_stat_error = array.array('d',d_branch_f_stat_error)
arr_d_branch_f_tot_error = array.array('d',d_branch_f_tot_error)
arr_tgraph_x = array.array('d',tgraph_x)
arr_tgraph_ex = array.array('d',tgraph_ex)

c1 = ROOT.TCanvas()
tgraph_bf = TGraphMultiErrors( "tgraph_bf", "", 8, np.array(arr_tgraph_x), np.array(arr_d_branch_f), np.array(arr_tgraph_ex), np.array(arr_tgraph_ex), np.array(arr_d_branch_f_stat_error), np.array(arr_d_branch_f_stat_error))
tgraph_bf.AddYError(8, np.array(arr_d_branch_f_tot_error), np.array(arr_d_branch_f_tot_error))
tgraph_bf.SetMinimum(0)
tgraph_bf.SetMaximum(9E-8)
tgraph_bf.SetMarkerStyle(20)
tgraph_bf.SetLineColor(2)
# 632 = kRed
tgraph_bf.GetAttLine(0).SetLineColor(632)
tgraph_bf.GetAttLine(0).SetLineWidth(2)
# 600 = kBlue
tgraph_bf.GetAttLine(1).SetLineColor(600)
tgraph_bf.GetAttFill(1).SetFillStyle(0)
tgraph_bf.GetYaxis().SetTitle("Differential Branching Fraction [ 1 / GeV^{2}]")
tgraph_bf.GetXaxis().SetTitle("q^{2} bins [GeV^{2}]")
tgraph_bf.Draw("APS")
legend = TLegend(0.6,0.7,0.9,0.9)
arr_tgraph_temp = array.array('d',tgraph_temp)
tgraph_legend1 = TGraphErrors(8, np.array(arr_tgraph_x), np.array(arr_d_branch_f), np.array(arr_tgraph_temp), np.array(arr_d_branch_f_stat_error))
#tgraph_bf.SetMarkerStyle(0)
tgraph_legend1.SetLineColor(632)
tgraph_legend1.SetLineWidth(2)
legend.AddEntry(tgraph_legend1,"Statistical Error","lef")
#tgraph_bf = TGraphMultiErrors( "tgraph_bf", "", 8, np.array(arr_tgraph_x), np.array(arr_d_branch_f), np.array(arr_tgraph_temp), np.array(arr_tgraph_temp), np.array(arr_d_branch_f_tot_error), np.array(arr_d_branch_f_tot_error))
#tgraph_bf.GetAttLine(0).SetLineColor(600)
tgraph_legend2 = TGraphErrors(8, np.array(arr_tgraph_x), np.array(arr_d_branch_f), np.array(arr_tgraph_temp), np.array(arr_d_branch_f_tot_error))
tgraph_legend2.SetLineColor(600)
tgraph_legend2.SetLineWidth(1)
legend.AddEntry(tgraph_legend2,"Total Error","lef")
#print legend.GetEntry()
legend.Draw()
c1.SaveAs("branching_fraction.pdf")
#tgraph_bf.Write("tgraph_branching_fraction%s"%(args.year))
#yields[4] = 0
#yields[6] = 0
#yields_error[4] = temp4_ex_err
#yields_error[6] = temp6_ex_err
#tgraph_ex[4] = temp4_ex_err
#tgraph_ex[6] = temp6_ex_err
#yields_error[7] = 0
#tgraph_x[4] = temp4_x
#tgraph_x[6] = temp6_x
#yields_error
print 'yields'
print yields
print 'yields_error'
print yields_error
print 'x_err'
print tgraph_ex
#arr_tgraph_ex = array.array('d',tgraph_ex)
arr_yields = array.array('d', yields)
arr_yields_error = array.array('d',yields_error)
#arr_tgraph_corrected_x = array.array('d', tgraph_corrected_x)

c2 = ROOT.TCanvas()
c2.SetLogy()
#tgraph_yield = TGraphErrors(7, np.array(arr_tgraph_corrected_x), np.array(arr_yields), np.array(arr_tgraph_ex), np.array(arr_yields_error))
tgraph_yield = TGraphErrors(8, np.array(arr_tgraph_x), np.array(arr_yields), np.array(arr_tgraph_ex), np.array(arr_yields_error))
tgraph_yield.SetMinimum(1)
tgraph_yield.SetMaximum(1E7)
tgraph_yield.SetTitle()
tgraph_yield.SetMarkerStyle(20)
tgraph_yield.GetYaxis().SetTitle("Yield")
tgraph_yield.GetXaxis().SetTitle("q^{2} bins [GeV^{2}]")
tgraph_yield.Draw("AP")
c2.SaveAs("yields.pdf")


#tgraph_ex[4] = temp4_ex_err
#tgraph_ex[6] = temp6_ex_err
#tgraph_ex[7] = temp7_ex_err
arr_tgraph_ex = array.array('d',tgraph_ex)
arr_eff = array.array('d', efficiency)
arr_eff_error = array.array('d',efficiency_error)
print 'efficiency'
print efficiency
print 'eff_error'
print efficiency_error
print 'x_err'
print tgraph_ex
c3 = ROOT.TCanvas()
ef_title = TStyle("ef_title","ef_title")
tgraph_efficiency = TGraphErrors(8, np.array(arr_tgraph_x), np.array(arr_eff), np.array(arr_tgraph_ex), np.array(arr_eff_error))
tgraph_efficiency.SetTitle()
tgraph_efficiency.SetMarkerStyle(20)
ef_title.SetTitleOffset(0.01,"y")
tgraph_efficiency.GetYaxis().SetTitle("Efficiency")
tgraph_efficiency.GetXaxis().SetTitle("q^{2} bins [GeV^{2}]")
tgraph_efficiency.Draw("AP")
c3.SaveAs("efficiencies.pdf")

print 'It took', time.time() - T, 'seconds'
print (time.time()- T)/60, 'minutes'
print (time.time() - T)/3600, 'hours'

f_NPsi.Close()
f_JPsi.Close()
f_PsiP.Close()
bf_file.Close()