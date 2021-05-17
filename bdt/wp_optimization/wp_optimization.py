import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("year"      , help = "choose among:2016,2017,2018", default = '2018')
args = parser.parse_args()
year = args.year

import ROOT, math, sys
import numpy as np
# import pandas, root_numpy

from collections import OrderedDict
from ROOT import gSystem
gSystem.Load('libRooFit')
from ROOT import RooFit, RooRealVar, RooDataSet, RooArgList, RooTreeData, RooArgSet, RooAddPdf, RooFormulaVar, RooExtendPdf
from ROOT import RooGaussian, RooExponential, RooChebychev

ROOT.gROOT.SetBatch(True)
# ROOT.gErrorIgnoreLevel = ROOT.kFatal
# ROOT.RooMsgService.instance().setSilentMode(True)

from os import path
sys.path.append( path.dirname( path.dirname( path.abspath('../../utils/utils.py') ) ) )
from utils.utils import *

ROOT.RooMsgService.instance().setGlobalKillBelow(4)

nSigma_psiRej = 3.

outfile   = ROOT.TFile('optNov2020/outcome_wp_finding_%s_lmnrPlusCharm_punzi_noTkMu_addHLT.root'%args.year, 'recreate')

n_bdt_points = 20
sample_range = 11


folder = '/gwpool/users/fiorendi/p5prime/miniAOD/CMSSW_10_2_14/src/miniB0KstarMuMu/miniKstarMuMu/selection_and_fits/bdt/'


def fitNSig(ibdt, fulldata,isample):
    mean        = RooRealVar ("mass"         , "mean"          ,  B0Mass_,   3,    7, "GeV")
    sigma       = RooRealVar ("#sigma_{1}"   , "sigma"         ,  0.028,     0,   10, "GeV")
    signalGauss = RooGaussian("signalGauss"  , "signal gauss"  ,  theBMass,  mean,sigma)
    
    sigma2       = RooRealVar ("#sigma_{2}"       , "sigma2"         ,  0.048,     0,   0.09, "GeV")
    signalGauss2 = RooGaussian("signalGauss2"  , "signal gauss2"  ,  theBMass,  mean,sigma2)
    f1           = RooRealVar ("f1"            , "f1"             ,  0.8  ,     0.,   1.)
    gaus         = RooAddPdf  ("gaus"          , "gaus1+gaus2"    , RooArgList(signalGauss,signalGauss2), RooArgList(f1))
    
    pol_c1      = RooRealVar ("p1"           , "coeff x^0 term",    0.5,   -10, 10);
    pol_c2      = RooRealVar ("p2"           , "coeff x^1 term",    0.5,   -10, 10);
    # pol_c3      = RooRealVar ("p3"           , "coeff x^2 term",    0.5,   -10, 10);
    # slope       = RooRealVar ("slope"        , "slope"         ,    0.5,   -10, 10);
    # bkg_exp     = RooExponential("bkg_exp"   , "exponential"   ,  slope,   theBMass  );
    bkg_pol     = RooChebychev("bkg_pol"     , "2nd order pol" ,  theBMass, RooArgList(pol_c1));
    
    nsig        = RooRealVar("Yield"         , "signal frac"   ,   40000,     0,   1000000);
    nbkg        = RooRealVar("nbkg"          , "bkg fraction"  ,    1000,     0,   550000);
    
    cut = cut_base + '&& bdt_prob > %s'%(ibdt)
    
    data       = fulldata.reduce(RooArgSet(theBMass,mumuMass,mumuMassE), cut)
    fitFunction = RooAddPdf ("fitfunction" , "fit function"  ,  RooArgList(gaus, bkg_pol), RooArgList(nsig, nbkg))
    r = fitFunction.fitTo(data, RooFit.Extended(True), RooFit.Save(), RooFit.Range(4.9,5.6), RooFit.PrintLevel(-1))
    
    frame = theBMass.frame()
    data.plotOn(frame, RooFit.Binning(70), RooFit.MarkerSize(.7))
    fitFunction.plotOn(frame, );
    fitFunction.plotOn(frame, RooFit.Components("bkg_pol")     , RooFit.LineStyle(ROOT.kDashed));
    fitFunction.plotOn(frame, RooFit.Components("signalGauss") , RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+1));
    fitFunction.plotOn(frame, RooFit.Components("signalGauss2"), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kOrange+1));
    
    parList = RooArgSet (nsig,sigma,sigma2, mean)
    ###### fitFunction.plotOn(frame, RooFit.Components("signalGauss2"), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+2));
    
    fitFunction.paramOn(frame, RooFit.Parameters(parList), RooFit.Layout(0.62,0.86,0.88))
    canv  = ROOT.TCanvas()
    frame.Draw()
#     canv.SaveAs('sig_fit_bdt%f_sample%i.pdf'%(ibdt,isample))
    
    dict_s_v1[ibdt]  = [nsig.getVal(), nsig.getError()]
    dict_sigma[ibdt] = math.sqrt(f1.getVal()*(sigma.getVal()**2) + (1-f1.getVal())*(sigma2.getVal()**2))


def fitNBkg(ibdt, fullbkg,isample):

    slope       = RooRealVar    ("slope"     , "slope"         ,    0.5,   -10, 10);
    bkg_exp     = RooExponential("bkg_exp"   , "exponential"   ,  slope,   theBMass  );
    
    cut = cut_base + '&& bdt_prob > %s'%(ibdt)
    
    theBMass.setRange('sigRangeMC', B0Mass_ - 3*dict_sigma[ibdt], B0Mass_ + 3*dict_sigma[ibdt])
    
    databkg = fullbkg.reduce(RooArgSet(theBMass,mumuMass,mumuMassE), cut)
    r       = bkg_exp.fitTo(databkg, RooFit.Save(), ROOT.RooFit.Range('left,right'),  RooFit.PrintLevel(-1))
    
    frame = theBMass.frame()
    databkg.plotOn(frame, RooFit.Binning(70), RooFit.MarkerSize(.7))
    bkg_exp.plotOn(frame, );
    canv  = ROOT.TCanvas()
    frame.Draw()
    
    nbkg = RooRealVar  ('nbkg', 'bkg n'  ,  1000,       0,   550000)
    ebkg = RooExtendPdf('ebkg','ebkg'    ,  bkg_exp, nbkg, 'sigRangeMC')  ## here imposing the range to calculate bkg yield
    ebkg.fitTo(databkg, ROOT.RooFit.Range('left,right'),  RooFit.PrintLevel(-1))
    ebkg.plotOn(frame, RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+1), RooFit.Range(4.9,5.6));
    frame.Draw()
#     canv.SaveAs('bkg_fit_bdt%f_sample%i.pdf'%(ibdt,isample))
    
    dict_b_v1[ibdt] =  [nbkg.getVal(), nbkg.getError()]


for isample in range(sample_range):
    print 'now working on subsample', isample
    ifileMC   = folder + '/sub_samples/sample_%s_MC_LMNR_%s_newphi_addBDT.root'%(year,isample)
#     ifileMC   = folder + '/sub_samples/sample_%s_MC_LMNR_%s_newphi_addBDT_addHLT.root'%(year,isample)

    dict_s_v1  = OrderedDict()
    dict_b_v1  = OrderedDict()
    dict_sigma = OrderedDict()
    
    ### retrieve S from fitting the MC sample
    bMass     = RooRealVar("bMass"             , "#mu^{+}#mu^{-}K* mass", 2, 20, "GeV")
    bBarMass  = RooRealVar("bBarMass"          , "#mu^{+}#mu^{-}K* mass", 2, 20, "GeV")
    mumuMass  = RooRealVar("mumuMass"          , "mumuMass"             , 0, 1000)
    mumuMassE = RooRealVar("mumuMassE"         , "mumuMassE"            , 0, 10000)
    tagB0     = RooRealVar("tagB0"             , "tagB0"                , -1, 2)
    bdt_prob  = RooRealVar("bdt_prob"          , "bdt_prob"             ,  0.7, 2) ## already cut here on BDT !!!
    pass_pre  = RooRealVar("pass_preselection" , "pass_preselection"    ,  1  , 2) ## cut on preselection 
    
    thevars = RooArgSet()
    thevars.add(bMass)
    thevars.add(bBarMass)
    thevars.add(mumuMass)
    thevars.add(mumuMassE)
    thevars.add(tagB0)
    thevars.add(bdt_prob)



    tree = ROOT.TChain('ntuple')
    tree.AddFile(ifileMC)
    fulldata   = RooDataSet('fulldata', 'fulldataset', tree,  RooArgSet(thevars))
    
    ## add to the input tree the combination of the variables for the B0 arb. mass
    theBMassfunc = RooFormulaVar("theBMass", "#mu^{+}#mu^{-}K^{#pm}#pi^{#mp} mass [GeV]", "@0*@1 + (1-@0)*@2", RooArgList(tagB0,bMass,bBarMass) )
    ## add to the input tree the combination of the variables, to be used for the cuts on the dimuon mass
    theBMass     = fulldata.addColumn(theBMassfunc) ;
    theBMass.setRange(4.9,5.6);

    cut_base = '( mumuMass < 2.702 || mumuMass > 4.)'
#     cut_base = '(mumuMass*mumuMass > 4.3 && mumuMass*mumuMass < 6)'
    
    
    for i in range(n_bdt_points):
#         ibdt = i*0.002+0.97
#         ibdt = i*0.005+0.85
        ibdt = i*0.005+0.9
#         ibdt = i*0.001+0.985
        print 'sig, ibdt ', ibdt
        fitNSig(ibdt, fulldata,isample)
    
    
    ### retrieve B from fitting the data sample
    treeBkg = ROOT.TChain('ntuple')
    treeBkg.Add(folder + 'sub_samples/sample_%s_data_LMNR_%s_newphi_addBDT.root'%(args.year,isample))
    treeBkg.Add(folder + 'sub_samples/sample_%s_data_Charmonium_%s_newphi_addBDT.root'%(args.year,isample))

    fullbkg     = RooDataSet('fullbkg', 'fullbkg', treeBkg,  RooArgSet(thevars))
    theBMass    = fullbkg.addColumn(theBMassfunc) 
    
    theBMass.setRange(4.9,5.6);
    theBMass.setRange('left'    , 4.9, 5.13)
    theBMass.setRange('right'   , 5.4, 5.6 )
    
    
    for i in range(n_bdt_points):
#         ibdt = i*0.002+0.97
#         ibdt = i*0.005+0.92
        ibdt = i*0.005+0.9
#         ibdt = i*0.001+0.985
        print 'bkg, ibdt ', ibdt
        fitNBkg(ibdt, fullbkg, isample)
        
    
    sign_vec_v1 = []
    bkg_vec_v1  = []
    bdt_vec     = []  
    # 
    for k,v in dict_s_v1.iteritems():
        bdt_vec.append(k)
        sign_vec_v1.append(v[0])
    for k,v in dict_b_v1.iteritems():
        bkg_vec_v1.append(v[0])
    
    
    sign_array_v1 = np.asarray(sign_vec_v1)
    bkg_array_v1  = np.asarray(bkg_vec_v1)
    bdt_array     = np.asarray(bdt_vec )
    

    graph_v1 = ROOT.TGraph(len(sign_array_v1),bdt_array, sign_array_v1 ) # convert lists to arrays
    graph_v1.GetXaxis().SetTitle('bdt cut')
    graph_v1.GetYaxis().SetTitle('S')
    graph_v1.SetTitle('')
    graph_v1.SetMarkerStyle(8)
    graph_v1.SetName('graph_signal_sample%s'%isample)

    canv  = ROOT.TCanvas()
    graph_v1.Draw('AP')

#     canv.SaveAs('bdt_optimize_wp_sample%s_newSB.pdf'%isample)

    graph_v2 = ROOT.TGraph(len(bkg_array_v1),bdt_array, bkg_array_v1 ) # convert lists to arrays
    graph_v2.GetXaxis().SetTitle('bdt cut')
    graph_v2.GetYaxis().SetTitle('B')
    graph_v2.SetTitle('')
    graph_v2.SetMarkerStyle(8)
    graph_v2.SetName('graph_background_sample%s'%isample)


    outfile.cd()
    graph_v1.Write()
    graph_v2.Write()


outfile.Close()