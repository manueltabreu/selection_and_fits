import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("year"      , help = "choose among:2016,2017,2018", default = '2018')
parser.add_argument("--tag"      , help = "", default = '_punzi_noTkMu')
args = parser.parse_args()
year = args.year

import ROOT
import numpy, math, array, sys
from os import path
from ROOT import TGraph, TFile, TMultiGraph

NORM      = True
n_samples = 11

norm_factor = 1
sys.path.append( path.dirname( path.dirname( path.abspath('../../utils/utils.py') ) ) )
from utils.eras_allYears import *
data_lumi   = lumi_eras[args.year]
mc_lumi     = lumi_mc['LMNR' + args.year]
norm_factor = norm_fact[args.year]
scale = data_lumi/(mc_lumi * norm_factor)

print 'data lumi: ', data_lumi
print 'mc lumi: ', mc_lumi
print 'norm_factor: ', norm_factor
norm_factor = 2.5
print 'scale: ', scale
tag = args.tag # '_sign_yesTkMu'

infile = TFile('optNov2020/outcome_wp_finding_%s_lmnrPlusCharm%s.root'%(args.year,tag),'read')

graphs_vy = []
for i in range(n_samples):
    graphs_vy_sig_tmp = infile.Get('graph_signal_sample%s'%i).GetY()
    graphs_vy_bkg_tmp = infile.Get('graph_background_sample%s'%i).GetY()

    graphs_vy_sig_tmp.SetSize(infile.Get('graph_signal_sample%s'%i).GetN())
    graphs_vy_bkg_tmp.SetSize(infile.Get('graph_background_sample%s'%i).GetN())

    ## now convert to usual python objects
    vy = []
    vy_sig = array.array('f',graphs_vy_sig_tmp)
    vy_bkg = array.array('f',graphs_vy_bkg_tmp)
    for j in range(len(vy_sig)):
        vy.append(scale*vy_sig[j] / math.sqrt( scale*vy_sig[j] + vy_bkg[j]))
    
    graphs_vy_tmp2 = numpy.asarray(vy )
    graphs_vy.append(graphs_vy_tmp2)

    if i==0: 
        graphs_vx_tmp = infile.Get('graph_signal_sample%s'%i).GetX()
        graphs_vx_tmp.SetSize(infile.Get('graph_signal_sample%s'%i).GetN())
        vx = array.array('f',graphs_vx_tmp)
        graphs_vx = numpy.asarray(vx )


## find point with max average s/s+B
ave = []
for ibdt in range(len(graphs_vx)):
  theave = 0
  for i in range(n_samples):  
    theave = theave+graphs_vy[i][ibdt]
  ave.append(theave/n_samples)  
ind = numpy.argmax(ave)
bestavebdt = graphs_vx[ind]
print 'WP with highest ave BDT: ', bestavebdt

best = []
## now prepare plots
norm_graphs = []
for i in range(n_samples):
    graphs_vy_array = numpy.asarray(graphs_vy[i] )
    ## find highest BDT per subsample 
    ind = numpy.argmax(graphs_vy_array)
    bestbdt = graphs_vx[ind]
    best.append(bestbdt)
    
    norm_graphs.append( ROOT.TGraph(len(graphs_vx),graphs_vx, graphs_vy_array )) # convert lists to arrays
    norm_graphs[i].GetXaxis().SetTitle('bdt cut')
    norm_graphs[i].GetYaxis().SetTitle('S/#sqrt{S+B} a.u.')
    norm_graphs[i].SetTitle('')
    norm_graphs[i].SetMarkerStyle(8)
    norm_graphs[i].SetName('graph_sample%s'%i)
    norm_graphs[i].SetMarkerColor(ROOT.kViolet+i)
    norm_graphs[i].SetMarkerStyle(20+i)
#     norm_graphs[i].GetYaxis().SetRangeUser(8,15)

    print 'after plotds WP with highest ave BDT: ', bestbdt

## add graph with ave BDT
ave_array = numpy.asarray(ave )
graph_ave = ROOT.TGraph(len(graphs_vx),graphs_vx, ave_array ) 
graph_ave.GetXaxis().SetTitle('bdt cut')
graph_ave.GetYaxis().SetTitle('S/#sqrt{S+B} a.u.')
graph_ave.SetTitle('')
graph_ave.SetMarkerStyle(8)
graph_ave.SetName('average')
graph_ave.SetMarkerColor(ROOT.kAzure+7)
graph_ave.SetMarkerStyle(47)
# graph_ave.GetYaxis().SetRange(8,15)


mg = TMultiGraph()
mg.SetTitle('')
# mg.SetTitle('%s'%tag.replace('_', ' - ').replace('yes', 'keep ').replace('- p', ' p').replace('- s', ' s').replace('no', 'remove '))

for i in range(n_samples):
  mg.Add(norm_graphs[i])

mg.Add(graph_ave)

canv = ROOT.TCanvas()
mg.Draw('AP')
mg.GetYaxis().SetRangeUser(5,11)
mg.GetXaxis().SetTitle('BDT score')
mg.GetYaxis().SetTitle('S/#sqrt{S+B}')

leg = ROOT.TLegend(0.35,0.14,0.7,0.46)
for i in range(n_samples):
    leg.AddEntry(norm_graphs[i], 'BDT %s,  highest sign.WP: %.3f'%(i, best[i]), 'p')
leg.AddEntry(graph_ave, 'ave sign.,   highest sign.WP: %.3f'%bestavebdt, 'p')
leg.SetBorderSize(0)
leg.Draw()
canv.SaveAs('optNov2020/bdt_wps_%s_%.4f_%s.pdf'%(args.year,scale,tag))

