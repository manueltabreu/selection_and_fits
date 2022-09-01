import argparse
from ast import Yield
from calendar import c
parser = argparse.ArgumentParser(description="")
# parser.add_argument("inputfile" , help = "Path to the input ROOT file")
parser.add_argument("dimusel"   , help = "Define if keep or remove dimuon resonances. You can choose: keepPsiP, keepJpsi, rejectPsi, keepPsi")
parser.add_argument("year"      , help = "choose among:2016,2017,%s", default = '2018')
args = parser.parse_args()

'''
code to fit the B0 mass distribution:
- unbinned fit
- possibility to apply cuts on the dimuon mass [B0&Psi cut in RunI analysis] (e.g. to exclude the Jpsi mass region, or the psi) via the parameter dimusel
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

def _getFittedVar(varName, w=None):
#    if w is not None:
#        return ufloat (w.var(varName).getVal() , w.var(varName).getError())
#    else :
#        return ufloat (varName.getVal()        , varName.getError())
    if w is not None:
	val = w.var(varName).getVal()
	err = w.var(varName).getError()
        return 
    else :
	val = varName.getVal()
        err = varName.getError()
        return 

def _goodFit(r):
    return (r.status()==0 and r.covQual() == 3)

def _accFit(r):
    return (r.status()==4 and r.covQual() == 3)

def _writeFitStatus(r):
    str_status = "GOOD" if r.status()==0 else "NOT CONV"
    txt = ROOT.TLatex(.16,.7, "fit status: " + str_status + ", covQ = %s" %r.covQual() )
    txt . SetNDC() ;
    txt . SetTextSize(0.033) ;
    txt . SetTextFont(42)
    return txt

def _writeChi2(chi2):
    txt = ROOT.TLatex(.16,.8, "fit #chi^{2}: %.1f "%chi2 )
    txt . SetNDC() ;
    txt . SetTextSize(0.033) ;
    txt . SetTextFont(42)
    return txt
    
#def _constrainVar(var, nsigma):
def _constrainVar(var, nsigma):
  
    if w is not None:
        val = w.var(var.GetName()).getVal()
        err = w.var(var.GetName()).getError()
    else :
        val = (var.GetName()).getVal()
        err = (var.GetName()).getError()

#    constr = _getFittedVar(var.GetName(), w)
#    gauss_constr = RooGaussian(  "c_%s" %var.GetName() , 
#                                 "c_%s" %var.GetName() , 
#                                var         ,  
#                                ROOT.RooFit.RooConst( constr.n ), 
#                                ROOT.RooFit.RooConst( nsigma*constr.s )
#                                ) 
#    print 'constraining var',   var.GetName(), ': ',     constr.n , ' with uncertainty ' , nsigma*constr.s                          
    constr_val = val
    constr_err = err
    gauss_constr = RooGaussian(  "c_%s" %var.GetName() ,
                                 "c_%s" %var.GetName() ,
                                var         ,
                                ROOT.RooFit.RooConst(constr_val),
                                ROOT.RooFit.RooConst( nsigma*constr_err )
                                )
    print 'constraining var',   var.GetName(), ': ',     constr_val , ' with uncertainty ' , nsigma*constr_err

    return gauss_constr                        


sys.path.append("../utils")
from utils.utils import *
from utils.fit_functions import *

nbins = 80
nSigma_psiRej = 3.
# cut_base      = applyB0PsiCut(args.dimusel, nSigma_psiRej)
cut_base      = 'passB0Psi_lmnr == 1 '
## entries per bin in data
n_bin = n_data[args.year]
## uncertainty on fRT
fm_sigma = fM_sigmas[args.year]

MC_num_entries = [0, 0, 0, 0, 0, 0, 0, 0]
MC_RT_num_entries = [0, 0, 0, 0, 0, 0, 0, 0]
MC_WT_num_entries = [0, 0, 0, 0, 0, 0, 0, 0]
MC_num_entries_error = [0, 0, 0, 0, 0, 0, 0, 0]

tgraph_y = [0, 0, 0, 0, 0, 0, 0, 0]
tgraph_x = [0, 0, 0, 0, 0, 0, 0, 0]
tgraph_ex = [0, 0, 0, 0, 0, 0, 0, 0]
tgraph_ey = [0, 0, 0, 0, 0, 0, 0, 0]

q2binningone = [
    10.09,
    12.86
]

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

from collections import OrderedDict
fitStats = OrderedDict()
covStats = OrderedDict()
chi2s    = OrderedDict()
print("start function declaration")
  
def fitMC(fulldata, correctTag, ibin):

    print 'now fitting: ', ibin, ' for ', correctTag*'correctTag ', (1-correctTag)*'wrongTag'  
    cut = cut_base + '&& (mumuMass*mumuMass > %s && mumuMass*mumuMass < %s)'%(q2binning[ibin], q2binning[ibin+1])
    fulldata_v2        = fulldata.reduce(RooArgSet(thevarsMC), cut)

    cut = 'rand < 100' ## valid for 2018
    data = fulldata_v2.reduce(RooArgSet(tagged_mass,mumuMass,mumuMassE), cut)


    pol_c1      = RooRealVar ("p1"           , "coeff x^0 term" ,  -0.5,   -10, 10);
    pol_c2      = RooRealVar ("p2"           , "coeff x^0 term" ,  0.5,   -10, 10);
    bkg_pol     = RooChebychev("bkg_pol"     , "bkg_pol" ,  tagged_mass, RooArgList(pol_c1));
    slope       = RooRealVar    ("slope"   , "slope"           ,   -4,   -10, 10);
    bkg_exp     = RooExponential("bkg_exp" , "exponential"     ,  slope,   tagged_mass  );
    fsig        = RooRealVar("fsig"         , "fsig"   ,   0.9, 0, 1)
#     signalFunction = bkg_pol ### just a placeholder

    nsig        = RooRealVar("Yield"         , "nsig"   ,   500000,     0,    1E7)
    nbkg        = RooRealVar("nbkg"          , "nbkg"   ,     1000,     0,    1E6 )
    
    doextended = False
    fitrange   = "mcrange"
    tag        = 'RT' if correctTag else 'WT'
    if correctTag:
        fitrange   = "datarange"
        mean        = RooRealVar ("mean_{RT}^{%s}"%ibin,        "massRT"         , B0Mass_,     5.2,    5.4, "GeV")
#        mean        = RooRealVar ("mean_{RT}^{%s}"%ibin,        "massRT"         , B0Mass_,     5.2,    5.8, "GeV")
       
        MC_RT_num_entries[ibin] = fulldata_v2.sumEntries()
        if ibin < 4:  #4  
            ### double CB, 1 sigma
            sigmaCB     = RooRealVar ("#sigma_{RT1}^{%s}"%ibin, "sigmaCB"        ,  0.03  ,     0,   1  )#1
            alpha1      = RooRealVar ("#alpha_{RT1}^{%s}"%ibin,  "#alpha_{1}"     ,  0.5   ,    0,  10  )#10
            alpha2      = RooRealVar ("#alpha_{RT2}^{%s}"%ibin,  "#alpha_{2}"     ,  2.5   ,    0,  10  )#10
            n1          = RooRealVar ("n_{RT1}^{%s}"%ibin,       "n_1"            ,  1     ,    0.01,  20  )
            n2          = RooRealVar ("n_{RT2}^{%s}"%ibin,       "n_2"            ,  1     ,    0.01,  100  )
            doublecb = ROOT.RooDoubleCBFast("doublecb_RT%s"%ibin, "doublecb", tagged_mass, mean, sigmaCB, alpha1, n1, alpha2, n2)	
            signalFunction = doublecb

        elif ibin == 7:  #4  
            ### CB + gaus
            crystalBall( mean  , initial_sigma1,  1.5,  1, tagged_mass, w, 'RT1', ibin, [0, 10])
            singleG(mean, initial_sigma2, tagged_mass, w, 'RT2' , ibin)
            gausCB( w.pdf("cbshape_RT1_%s"%ibin) , w.pdf("gaus_RT2_%s"%ibin), 0.5, tagged_mass, w, 'RT', ibin )
            signalFunction = w.pdf("gauscb_RT_%s"%ibin)   

        else:
            ### double CB, 2 sigmas
            crystalBall( mean  , initial_sigma1,  1.5,  1, tagged_mass, w, 'RT1', ibin, [0, 10])
            crystalBall( mean  , initial_sigma2,  -2 ,  1, tagged_mass, w, 'RT2', ibin, [-10, 0])
            doubleCB     ( w.pdf("cbshape_RT1_%s"%(ibin)), w.pdf("cbshape_RT2_%s"%ibin), 0.8  , tagged_mass, w, "RT%s"%ibin)
            signalFunction = w.pdf("doublecb_RT%s"%ibin)   
          
        fitFunction = signalFunction  ## sara 04.11 to uncomment

    else:
        MC_WT_num_entries[ibin] = fulldata_v2.sumEntries()
        ### double CB, 1 sigma
        mean        = RooRealVar ("mean_{WT}^{%s}"%ibin,     "massWT"         , B0Mass_,     5,    6, "GeV")
        sigmaCB     = RooRealVar ("#sigma_{WT1}^{%s}"%ibin,  "sigmaCB"        ,  0.03  ,    0,   1  )
        alpha1      = RooRealVar ("#alpha_{WT1}^{%s}"%ibin,  "#alpha_{1}"     ,  0.5   ,    0,  10  )
        alpha2      = RooRealVar ("#alpha_{WT2}^{%s}"%ibin,  "#alpha_{2}"     ,  2.5   ,    0,  10  )
        n1          = RooRealVar ("n_{WT1}^{%s}"%ibin,       "n_1"            ,  1     ,    0.01,  100  )
        n2          = RooRealVar ("n_{WT2}^{%s}"%ibin,       "n_2"            ,  1     ,    0.01,  100  )
        doublecb = ROOT.RooDoubleCBFast("doublecb_%s"%ibin, "doublecb", tagged_mass, mean, sigmaCB, alpha1, n1, alpha2, n2)	
        signalFunction = doublecb
        fitFunction    = doublecb  

#     fitFunction = RooAddPdf("fitfunction%s"%ibin, "fitfunction%s"%ibin , RooArgList(signalFunction, bkg_exp), RooArgList(fsig))

    getattr(w,"import")(fitFunction)
    print 'fitting the mc...'
    r = fitFunction.fitTo(data, RooFit.Extended(doextended), RooFit.Save(), RooFit.Range(fitrange))
    print 'fit status: ', r.status(), r.covQual() 
#    total_integral = signalFunction.createIntegral(RooArgSet(tagged_mass), RooFit.NormSet(RooArgSet(tagged_mass)), RooFit.Range('fitrange')).getVal()
    center = mean.getVal()
    initial_sigma = 0.03
    for dsigma in np.linspace(0,0.3,300):
        tmp_sigma = initial_sigma + dsigma
        tagged_mass.setRange('tmp_range', center-tmp_sigma, center+tmp_sigma)
#        tmp_int = signalFunction.createIntegral(RooArgSet(tagged_mass), RooFit.NormSet(RooArgSet(tagged_mass)), RooFit.Range('tmp_range')).getVal()
#        if tmp_int/total_integral > 0.68:
#          print tag, ' reached sigma:', tmp_sigma, ' for which we have ', tmp_int/total_integral 
#          break
    
    ## draw everything 
    params = fitFunction.getParameters(RooArgSet(tagged_mass)) 
    w.saveSnapshot("reference_fit_%s_%s"%(tag, ibin),params,ROOT.kTRUE) 
    frame = tagged_mass.frame(RooFit.Range(fitrange))#, RooFit.Title('correctly'*correctTag + 'wrongly'*(1-correctTag) + ' tagged events'))
    data.plotOn(frame, RooFit.Binning(nbins), RooFit.MarkerSize(.7))
    
    drawPdfComponents(fitFunction, frame, ROOT.kGreen if correctTag else ROOT.kViolet, RooFit.NormRange(fitrange), RooFit.Range(fitrange), isData=False)
    fitFunction.plotOn(frame, RooFit.NormRange(fitrange), RooFit.Range(fitrange) )
    fitFunction.paramOn(frame,  RooFit.Layout(0.62,0.86,0.88))

    frame.Draw()
    niceFrame(frame, 'correctly'*correctTag + 'wrongly'*(1-correctTag) + ' tagged events')
    frame. addObject(_writeFitStatus(r))
    
    ## evaluate sort of chi2 and save number of RT/WT events
    observables = RooArgSet(tagged_mass)
    flparams    = fitFunction.getParameters(observables)
    nparam      = int(flparams.selectByAttrib("Constant",ROOT.kFALSE).getSize())
    pdfstring = ''
#     pdb.set_trace()
    if correctTag:
        pdfstring = "doublecb_RT%s_Norm[tagged_mass]_Comp[doublecb_RT%s]_Range[datarange]_NormRange[datarange]"%(ibin,ibin)
	nRT_fromMC = data.sumEntries()
#         pdfstring = "doublegaus_RT%s_Norm[tagged_mass]_Comp[doublegaus_RT%s]_Range[mcrange]_NormRange[mcrange]"%(ibin,ibin)
#         pdfstring = "gauscb_RT%s_Norm[tagged_mass]_Comp[gauscb_RT%s]_Range[mcrange]_NormRange[mcrange]"%(ibin,ibin)
#         pdfstring = "expGaussExp_RT%s_Norm[tagged_mass]_Comp[expGaussExp_RT%s]_Range[mcrange]_NormRange[mcrange]"%(ibin,ibin)

        if ibin == 7:   
            pdfstring = "gauscb_RT_%s_Norm[tagged_mass]_Comp[gauscb_RT_%s]_Range[datarange]_NormRange[datarange]"%(ibin,ibin)

        if doextended:
            dict_s_rt[ibin]   = _getFittedVar(nsig)
#        else:
#            dict_s_rt[ibin]    = ufloat(data.sumEntries(), math.sqrt(data.sumEntries()))
        nRT_fromMC = data.sumEntries()
        nRT_fromMC_err = math.sqrt(data.sumEntries())
        nRT = RooRealVar ("nRT_%s"%ibin, "yield of RT signal",0,1.E9)
#        nRT.setVal(  dict_s_rt[ibin].n)
#        nRT.setError(dict_s_rt[ibin].s)
        nRT.setVal(nRT_fromMC)
        nRT.setError(nRT_fromMC_err)
#        print 'setting nRT to ', dict_s_rt[ibin].n
	print 'setting nRT to ', nRT_fromMC
        getattr(w,"import")(nRT)

    else:
        pdfstring = "doublecb_%s_Norm[tagged_mass]_Comp[doublecb_%s]_Range[mcrange]_NormRange[mcrange]"%(ibin,ibin)     
#        dict_s_wt[ibin]    = ufloat(data.sumEntries(), math.sqrt(data.sumEntries()))
        nWT_fromMC = data.sumEntries()
        nWT_fromMC_err = math.sqrt(data.sumEntries())
        nWT = RooRealVar ("nWT_%s"%ibin, "yield of WT signal",0,1.E7)
#        nWT.setVal(  dict_s_wt[ibin].n)
#        nWT.setError(dict_s_wt[ibin].s)
#        print 'setting nWT to ', dict_s_wt[ibin].n
	nWT.setVal(nWT_fromMC)
        nWT.setError(nWT_fromMC_err)
        print 'setting nWT to ', nWT_fromMC
        getattr(w,"import")(nWT)
    
#     pdfstring = "fitfunction%s_Norm[tagged_mass]_Range[datarange]_NormRange[datarange]"%(ibin)  ## sara 04.11 to comment

    ## eval and save goodness of fit indicators
    chi2s['%s%s'%(tag,ibin)] = frame.chiSquare(pdfstring, "h_fullmc",  nparam)
    frame. addObject(_writeChi2( chi2s['%s%s'%(tag,ibin)] ))
    fitStats['%s%s'%(tag,ibin)] = r.status()
    covStats['%s%s'%(tag,ibin)] = r.covQual()

    c1 = ROOT.TCanvas() 
    upperPad = ROOT.TPad('upperPad' , 'upperPad' , 0., 0.35 , 1.,  1.    )  
    lowerPad = ROOT.TPad('lowerPad' , 'lowerPad' , 0., 0.0  , 1.,  0.345 )  
    upperPad.SetBottomMargin(0.012)
    lowerPad.SetTopMargin(0)
    lowerPad.SetBottomMargin(0.2)

    upperPad.Draw()
    lowerPad.Draw()

    upperPad.cd()
    if not args.year=='test':  writeCMS(frame, args.year, [ q2binning[ibin], q2binning[ibin+1] ], -1)
    frame.Draw()
    
    ## add plot of pulls
    lowerPad.cd()
    hpull  = frame.pullHist("h_fullmc", pdfstring)
    frame2 = tagged_mass.frame(RooFit.Range(fitrange), RooFit.Title(''))
    frame2.addPlotable(hpull,"P") 
    niceFrameLowerPad(frame2, 'pull')
    frame2.Draw()
    line = ROOT.TLine(4.9,0,5.7,0)
    if correctTag:
        line = ROOT.TLine(5.0,0,5.6,0)
    line.SetLineColor(ROOT.kGreen+3)
    line.Draw()
 
   
    ## save to pdf and root files
    for ilog in [False]:
        upperPad.SetLogy(ilog)
        c1.SaveAs('fit_results_mass_checkOnMC/newbdt_puw/save_fit_mc_%s_%s_%s_newSigmaFRT%s%s_asInANv8.pdf'%(ibin, args.year, tag, '_logScale'*ilog, '_Jpsi'*(args.dimusel=='keepJpsi')))
    out_f.cd()
    r.Write('results_%s_%s'%(tag, ibin))
    
#    return dict_s_rt[ibin].n if correctTag else dict_s_wt[ibin].n
    return nRT_fromMC if correctTag else nWT_fromMC
   

## outdated, not used   
def fitData(fulldata, ibin, nRT_fromMC, nWT_fromMC):

    cut  = cut_base + '&& (mumuMass*mumuMass > %s && mumuMass*mumuMass < %s)'%(q2binning[ibin], q2binning[ibin+1])
    fulldata_v2 = fulldata.reduce(RooArgSet(tagged_mass,mumuMass,mumuMassE, randVar), cut)

    ## reduce to data-like statistics
    nDataEntries = fulldata_v2.sumEntries()
    nDesired = n_bin[ibin]/nDataEntries
    cut = 'rand < %f'%nDesired
    data = fulldata_v2.reduce(RooArgSet(tagged_mass,mumuMass,mumuMassE), cut)

#    fraction = dict_s_wt[ibin] / (dict_s_rt[ibin] + dict_s_wt[ibin])
    nRT_fromMC_err = math.sqrt(nRT_fromMC)
    nWT_fromMC_err = math.sqrt(nWT_fromMC)
    fraction_val = nWT_fromMC / (nRT_fromMC + nWT_fromMC)
    fraction_err = nWT_fromMC_err / (nRT_fromMC_err + nWT_fromMC_err)
#    print 'mistag fraction on MC for bin ', ibin , ' : ' , fraction.n , '+/-', fraction.s 
    print 'mistag fraction on MC for bin ', ibin , ' : ' , fraction_val , '+/-', fraction_err
    
    ### creating RT component
    w.loadSnapshot("reference_fit_RT_%s"%ibin)
    mean_rt      = w.var("mean_{RT}^{%s}"%ibin)
    sigma_rt1    = w.var("#sigma_{RT1}^{%s}"%ibin)

    sigma_rt2 = RooRealVar()
    alpha_rt1 = RooRealVar()
    alpha_rt2 = RooRealVar()
    n_rt1     = RooRealVar()
    n_rt2     = RooRealVar()
    f1rt     = RooRealVar()

    ## double cb fast
    if ibin < 4:
        alpha_rt1    = w.var("#alpha_{RT1}^{%s}"%ibin)
        alpha_rt2    = w.var("#alpha_{RT2}^{%s}"%ibin)
        n_rt1        = w.var("n_{RT1}^{%s}"%ibin)
        n_rt2        = w.var("n_{RT2}^{%s}"%ibin)
    ## for gauscb
    elif ibin == 7:
        sigma_rt2    = w.var("#sigma_{RT2}^{%s}"%ibin)
        alpha_rt1    = w.var("#alpha_{RT1}^{%s}"%ibin)
        n_rt1        = w.var("n_{RT1}^{%s}"%ibin)
        f1rt         = w.var("f^{RT%s}"%ibin)
        print (Yield)
    ## double cb old
    else:
        sigma_rt2    = w.var("#sigma_{RT2}^{%s}"%ibin)
        alpha_rt1    = w.var("#alpha_{RT1}^{%s}"%ibin)
        alpha_rt2    = w.var("#alpha_{RT2}^{%s}"%ibin)
        n_rt1        = w.var("n_{RT1}^{%s}"%ibin)
        n_rt2        = w.var("n_{RT2}^{%s}"%ibin)
        f1rt         = w.var("f^{RT%s}"%ibin)        

    theRTgauss  = w.pdf("doublecb_RT%s"%ibin)   

    ### creating WT component
    w.loadSnapshot("reference_fit_WT_%s"%ibin)
    mean_wt      = w.var("mean_{WT}^{%s}"%ibin)
    sigma_wt     = w.var("#sigma_{WT1}^{%s}"%ibin)
    alpha_wt1    = w.var("#alpha_{WT1}^{%s}"%ibin)
    alpha_wt2    = w.var("#alpha_{WT2}^{%s}"%ibin)
    n_wt1        = w.var("n_{WT1}^{%s}"%ibin)
    n_wt2        = w.var("n_{WT2}^{%s}"%ibin)
    theWTgauss   = w.pdf("doublecb_%s"%ibin)   

    ### creating constraints
    c_vars = RooArgSet()
    c_pdfs_rt = RooArgSet()
    c_pdfs_wt = RooArgSet()

    c_sigma_rt1   = _constrainVar(sigma_rt1, 1)
    c_alpha_rt1   = _constrainVar(alpha_rt1, 1)
    c_n_rt1       = _constrainVar(n_rt1, 1)
    c_sigma_wt    = _constrainVar(sigma_wt,  1)
    c_alpha_wt1   = _constrainVar(alpha_wt1, 1)
    c_alpha_wt2   = _constrainVar(alpha_wt2, 1)
    c_n_wt1       = _constrainVar(n_wt1, 1)
    c_n_wt2       = _constrainVar(n_wt2, 1)

    if ibin < 4:
        c_alpha_rt2   = _constrainVar(alpha_rt2, 1)
        c_n_rt2       = _constrainVar(n_rt2, 1)

        c_pdfs_rt = RooArgSet(c_sigma_rt1, c_alpha_rt1, c_alpha_rt2, c_n_rt1, c_n_rt2)
        c_vars = RooArgSet(sigma_rt1,     alpha_rt1,   alpha_rt2,   n_rt1,   n_rt2)

    elif ibin == 7:
        c_sigma_rt2   = _constrainVar(sigma_rt2, 1)
        c_f1rt        = _constrainVar(f1rt, 1)

        c_pdfs_rt = RooArgSet(c_sigma_rt1, c_sigma_rt2, c_alpha_rt1, c_n_rt1, c_f1rt)
        c_vars = RooArgSet(  sigma_rt1,   sigma_rt2,   alpha_rt1,  n_rt1,   f1rt)
        
    else:
        c_alpha_rt2   = _constrainVar(alpha_rt2, 1)
        c_n_rt2       = _constrainVar(n_rt2, 1)
        c_sigma_rt2   = _constrainVar(sigma_rt2, 1)
        c_f1rt        = _constrainVar(f1rt, 1)

        c_pdfs_rt = RooArgSet(c_sigma_rt1, c_sigma_rt2, c_alpha_rt1, c_alpha_rt2, c_n_rt1, c_n_rt2, c_f1rt)
        c_vars = RooArgSet(  sigma_rt1,   sigma_rt2,   alpha_rt1,   alpha_rt2,   n_rt1,   n_rt2,   f1rt)
    
    c_pdfs_wt = RooArgSet(c_sigma_wt, c_alpha_wt1, c_alpha_wt2, c_n_wt1, c_n_wt2)
    c_vars.add(sigma_wt)
    c_vars.add(alpha_wt1)
    c_vars.add(alpha_wt2)
    c_vars.add(n_wt1)
    c_vars.add(n_wt2)

    constr_list_rt = RooArgList(c_pdfs_rt)
    constr_list_rt.add(theRTgauss)
    c_theRTgauss = RooProdPdf ("c_theRTgauss%s"%ibin, "c_theRTgauss", constr_list_rt)     
    constr_list_wt = RooArgList(c_pdfs_wt)
    constr_list_wt.add(theRTgauss)
    c_theWTgauss = RooProdPdf ("c_theWTgauss%s"%ibin, "c_theWTgauss", constr_list_wt)     


    ### creating variable for the difference between the two peaks
    deltaPeaks = RooFormulaVar("deltaPeaks%s"%ibin, "@0 - @1", RooArgList(mean_rt, mean_wt))  
#    fm              = RooRealVar ("f_{M}^{%s}"%ibin , "fm"             , fraction.n , 0, 1)
    fm              = RooRealVar ("f_{M}^{%s}"%ibin , "fm"             , fraction_val , 0, 1) 
    signalFunction  = RooAddPdf  ("sumgaus%s"%ibin  , "rt+wt"          , RooArgList(c_theWTgauss,c_theRTgauss), RooArgList(fm))

    ### now create background parametrization
    slope         = RooRealVar    ("slope"      , "slope"           ,    0.5,   -10, 10);
    bkg_exp       = RooExponential("bkg_exp"    , "exponential"     ,  slope,   tagged_mass  );
    pol_c1        = RooRealVar    ("p1"         , "coeff x^0 term"  ,    0.5,   -10, 10);
    pol_c2        = RooRealVar    ("p2"         , "coeff x^1 term"  ,    0.5,   -10, 10);
    bkg_pol       = RooChebychev  ("bkg_pol"    , "2nd order pol"   ,  tagged_mass, RooArgList(pol_c1));
   
    nsig          = RooRealVar("Yield"         , "signal frac"    ,     nRT_fromMC + nWT_fromMC,     0,   1000000);
    nbkg          = RooRealVar("nbkg"          , "bkg fraction"   ,     1000,     0,   550000);

    print 'nsig is'
    print nsig.getVal()

    #trying to fit with the 2 gaussians
#    frame4 = tagged_mass.frame()
#    meanG = RooRealVar("meanG", "meanG", 13, 12, 14)
#    stdG1 = RooRealVar("stdG1", "stdG1", 0.1*2,0.,0.5*2)
#    stdG2 = RooRealVar("stdG2", "stdG2", 1, 0, 3)
#    G1 = RooGaussian("G1", "G1", tagged_mass, meanG, stdG1)
#    G2 = RooGaussian("G2", "G2", tagged_mass, meanG, stdG2) 
#    modelG = RooAddPdf("modelG", "modelG", RooArgList(G1,G2))

    ### creating constraints for the difference between the two peaks
    c_deltaPeaks = RooGaussian("c_deltaPeaks%s"%ibin , "c_deltaPeaks", deltaPeaks, ROOT.RooFit.RooConst( deltaPeaks.getVal() ), 
                                ROOT.RooFit.RooConst( 0.0005 )  ## value to be checked
                                ) 
    
# signal pdf
    c_pdfs = RooArgSet()
#     c_pdfs.add(c_deltaPeaks)
#     c_vars.add(deltaPeaks)

#    c_fm            = RooGaussian("c_fm%s"%ibin    , "c_fm"           , fm,  ROOT.RooFit.RooConst(fraction.n) , ROOT.RooFit.RooConst(fm_sigma[ibin]) )
    c_fm            = RooGaussian("c_fm%s"%ibin    , "c_fm"           , fm,  ROOT.RooFit.RooConst(fraction_val) , ROOT.RooFit.RooConst(fm_sigma[ibin]) ) 
    c_pdfs.add(c_fm)
    c_vars.add(fm)

    constr_list = RooArgList(c_pdfs)
#    constr_list.add(signalFunction)
    sig_bkg_function = RooAddPdf  ("totpdf%s"%ibin  , "wt+rt+bkg"          , RooArgList(signalFunction, bkg_exp), RooArgList(nsig,nbkg))
    constr_list.add(sig_bkg_function)
    totalFunction = RooProdPdf ("totalFunction%s"%ibin, "totalFunction", constr_list)   

#background pdf
#defined above

#total pdf (data)
    fitFunction = totalFunction   
    #fitFunction = RooAddPdf  ("totpdf%s"%ibin  , "wt+rt+bkg"          , RooArgList(signalFunction, bkg_exp), RooArgList(nsig,nbkg))
    #fitFunction = RooAddPdf  ("totpdf%s"%ibin  , "wt+rt+bkg"          , RooArgList(totalFunction, bkg_exp), RooArgList(nsig,nbkg))
    
    print 'here'
    
    print 'total pdf created'
    r = fitFunction.fitTo(data, 
                          RooFit.Range("datarange"), 
                          ROOT.RooFit.Constrain(c_vars),
                          ROOT.RooFit.Minimizer("Minuit2","migrad"),
                          ROOT.RooFit.Hesse(True),
                          ROOT.RooFit.Strategy(2),
                          ROOT.RooFit.Minos(False),
                         )
    print 'fit with Hesse strategy 2 done, now Minos'    
    r = fitFunction.fitTo(data, 
                          RooFit.Save(), 
                          RooFit.Range("datarange"), 
                          RooFit.Verbose(False),
                          ROOT.RooFit.Constrain(c_vars),
                          ROOT.RooFit.Minimizer("Minuit2","migrad"),
                          ROOT.RooFit.Hesse(True),
                          ROOT.RooFit.Strategy(2),
                          ROOT.RooFit.Minos(True),
                         )
                        
    r.Print()
    r.correlationMatrix().Print()
    fitStats['data%s'%(ibin)] = r.status()
    covStats['data%s'%(ibin)] = r.covQual()
    frame = tagged_mass.frame( RooFit.Range("datarange") )
    data.plotOn(frame, RooFit.Binning(nbins), RooFit.MarkerSize(.7))
    fitFunction.plotOn(frame, RooFit.NormRange("datarange"), RooFit.Range("datarange"))

    ## evaluate sort of chi2 and save number of RT/WT events
    observables = RooArgSet(tagged_mass)
    flparams = fitFunction.getParameters(observables)
    nparam = int(flparams.selectByAttrib("Constant",ROOT.kFALSE).getSize())
    pdfstring = "totalFunction%s_Norm[tagged_mass]_Range[datarange]_NormRange[datarange]"%ibin
#     pdfstring = "totalFunction_Norm[tagged_mass]_Range[full]_NormRange[full]"
    chi2s['data%s'%ibin] = frame.chiSquare(pdfstring, "h_fulldata",  nparam)
    frame. addObject(_writeChi2( chi2s['data%s'%ibin] ))

#   get the Yields
    print 'ibin and Yield'
    print (ibin)
    print (nsig)
    tgraph_y[ibin] = nsig.getVal()
    tgraph_ey[ibin] = nsig.getError()

    drawPdfComponents(fitFunction, frame, ROOT.kAzure, RooFit.NormRange("datarange"), RooFit.Range("datarange"), isData = True)

    parList1 = RooArgSet (nsig, alpha_wt1, alpha_wt2,  alpha_rt1, alpha_rt2, sigma_wt, sigma_rt1, sigma_rt2)
    parList2 = RooArgSet (mean_rt, n_wt1, n_wt2, n_rt1, n_rt2, f1rt, fm, nbkg)

#    parList = RooArgSet (nsig, mean_rt, mean_wt, sigma_rt1, sigma_rt2, alpha_rt1, alpha_rt2, mean_wt, sigma_wt)
#    parList.add(alpha_wt1)
#    parList.add(alpha_wt2)
#    parList.add(n_wt1)
#    parList.add(n_wt2)
#    parList.add(fm)
#    fitFunction.paramOn(frame, RooFit.Parameters(parList), RooFit.Layout(0.56,0.86,0.89))

#    fitFunction.paramOn(frame, RooFit.Parameters(parList2), RooFit.Layout(0.72,0.86,0.89)) 
  
    fitFunction.paramOn(frame, RooFit.Parameters(parList1), RooFit.Layout(0.56,0.70,0.89))
    frame.getAttText().SetTextSize(0.023)
    fitFunction.paramOn(frame, RooFit.Parameters(parList2), RooFit.Layout(0.72,0.80,0.89))
    frame.getAttText().SetTextSize(0.023)

    frame.Draw()
    niceFrame(frame, '')
    frame. addObject(_writeFitStatus(r))
    
    c1 = ROOT.TCanvas() 
    upperPad = ROOT.TPad('upperPad' , 'upperPad' , 0., 0.35 , 1.,  1.    )  
    lowerPad = ROOT.TPad('lowerPad' , 'lowerPad' , 0., 0.0  , 1.,  0.345 )  
    upperPad.SetBottomMargin(0.012)
    lowerPad.SetTopMargin(0)
    lowerPad.SetBottomMargin(0.2)

    upperPad.Draw(); lowerPad.Draw()
    upperPad.cd()
    frame.Draw()
    if not args.year=='test':  writeCMS(frame, args.year, [ q2binning[ibin], q2binning[ibin+1] ], 1)
    frame.Draw()

    ## add plot of pulls
    lowerPad.cd()
    hpull = frame.pullHist("h_fulldata", pdfstring)
    frame2 =  tagged_mass.frame(RooFit.Range("datarange"), RooFit.Title(''))
    frame2.addPlotable(hpull,"P") 
    niceFrameLowerPad(frame2, 'pull')
    frame2.Draw()
    line = ROOT.TLine(5.0,0,5.6,0)
    line.SetLineColor(ROOT.kGreen+3)
    line.Draw()

    for ilog in [False]:
        upperPad.SetLogy(ilog)
        c1.SaveAs('fit_results_mass_checkOnMC/newbdt_puw/save_fit_data_%s_%s_LMNR_newSigmaFRT_%s_asInANv8.pdf'%(ibin, args.year, '_logScale'*ilog))

    out_f.cd()
    r.Write('results_data_%s'%(ibin))

    params = fitFunction.getParameters(RooArgSet(tagged_mass)) 
    w.saveSnapshot("reference_fit_data_%s"%(ibin),params,ROOT.kTRUE) 
    getattr(w, 'import')(fitFunction)




tData = ROOT.TChain('ntuple')
tMC = ROOT.TChain('ntuple')

if args.year == 'test':
    tData.Add('/gwteray/users/fiorendi/final_ntuples_p5prime_allyears/2016Data_100k.root')
    tMC.Add('/gwteray/users/fiorendi/final_ntuples_p5prime_allyears/2016MC_LMNR_100k.root')
else:    
    if args.dimusel == 'rejectPsi':
        tData.Add('/lstore/cms/boletti/Run2-BdToKstarMuMu/%sdata_newphi_punzi_removeTkMu_fixBkg_B0Psicut_addxcutvariable.root'%args.year)
        tMC.Add('/lstore/cms/boletti/Run2-BdToKstarMuMu/%sMC_LMNR_newphi_punzi_removeTkMu_fixBkg_B0Psicut_addxcutvariable.root'%args.year)
#         tMC.Add('/eos/cms/store/user/fiorendi/p5prime/2016/skims/newphi/2016MC_LMNR.root')
    elif args.dimusel == 'keepJpsi':
        tData.Add('/lstore/cms/boletti/Run2-BdToKstarMuMu/%sdata_newphi_punzi_removeTkMu_fixBkg_B0Psicut_addxcutvariable.root'%args.year)
        tMC.Add('/lstore/cms/boletti/Run2-BdToKstarMuMu/%sMC_JPSI_newphi_punzi_removeTkMu_fixBkg_B0Psicut_addxcutvariable.root'%args.year)
    elif args.dimusel == 'keepPsiP':
        tData.Add('/lstore/cms/boletti/Run2-BdToKstarMuMu/%sdata_newphi_punzi_removeTkMu_fixBkg_B0Psicut_addxcutvariable.root'%args.year)
        tMC.Add('/lstore/cms/boletti/Run2-BdToKstarMuMu/%sMC_PSI_newphi_punzi_removeTkMu_fixBkg_B0Psicut_addxcutvariable.root'%args.year)

print 'mc file name:',  tMC.GetFile().GetName()
print 'data file name:', tData.GetFile().GetName()

tagged_mass     = RooRealVar("tagged_mass" , "#mu^{+}#mu^{-}K#pi mass", 4.9, 5.7, "GeV")
mumuMass        = RooRealVar("mumuMass"    , "mumuMass" , 0, 6);
mumuMassE       = RooRealVar("mumuMassE"   , "mumuMassE", 0, 10000);
tagB0           = RooRealVar("tagB0"       , "tagB0"    , 0, 2);
passB0Psi_lmnr  = RooRealVar("passB0Psi_lmnr" , "passB0Psi_lmnr", -200, 2);
passB0Psi_jpsi  = RooRealVar("passB0Psi_jpsi" , "passB0Psi_jpsi", -200, 2);
passB0Psi_psip  = RooRealVar("passB0Psi_psip" , "passB0Psi_psip", -200, 2);
runN            = RooRealVar("runN"   ,   "runN", 0, 10000000)
weight          = RooRealVar("weight" , "weight", 0,10);

tagged_mass.setRange("datarange", 5.0,5.6) ;
# tagged_mass.setRange("datarange"  , 4.9,5.7) ;
tagged_mass.setRange("mcrange"  , 4.9,5.7) ;
thevars = RooArgSet()
thevars.add(tagged_mass)
thevars.add(mumuMass)
thevars.add(mumuMassE)
thevars.add(tagB0)
thevars.add(passB0Psi_lmnr)
thevars.add(passB0Psi_jpsi)
thevars.add(passB0Psi_psip)
thevars.add(runN)
thevars.add(weight)

print 'reading data...'
#fulldata   = RooDataSet('fulldata', 'fulldataset', tData,  RooArgSet(thevars), "weight")
#fulldataall   = RooDataSet('fulldata', 'fulldataset', tData,  RooArgSet(thevars), "weight")
# = fulldataall.reduce("runN>316000 && runN<317000")
#fulldata.Print()

# fastest run to Jpsi (~10 minutes) "runN > 316050 && runN < 316060"
# to the PsiP "runN > 316000 && runN <316100"
# fastest run to PsiP (~10 minutes) "runN > 316050 && runN < 316060"
#fulldata   = RooDataSet('fulldata', 'fulldataset', tData,  RooArgSet(thevars))

if args.dimusel == 'rejectPsi': 
    fulldata   = RooDataSet('fulldata', 'fulldataset', tData,  RooArgSet(thevars))       
if args.dimusel == 'keepJpsi': 
    fulldata   = RooDataSet('fulldata', 'fulldataset', tData,  RooArgSet(thevars), "runN > 316050 && runN < 316060")   
if args.dimusel == 'keepPsiP':  
    fulldata   = RooDataSet('fulldata', 'fulldataset', tData,  RooArgSet(thevars), "runN > 316050 && runN <316060")    

print 'it worked :)'
#fulldata.printValue()

## add to the data dataset a random variable, in order to scale it to desired stat
nDataEntries = fulldata.sumEntries()
randVar  = RooRealVar("rand","rand",0,1) 
p0       = RooPolynomial("px","px",randVar) ;
rDataset = p0.generate(RooArgSet(randVar),int(nDataEntries))
fulldata.merge(rDataset) 

## add to the input tree the combination of the variables, to be used for the cuts on the dimuon mass
# deltaB0Mfunc = RooFormulaVar("deltaB0M", "deltaB0M", "@0 - @1", RooArgList(tagged_mass,B0Mass) )
# deltaJMfunc  = RooFormulaVar("deltaJpsiM" , "deltaJpsiM" , "@0 - @1", RooArgList(mumuMass,JPsiMass) )
# deltaPMfunc  = RooFormulaVar("deltaPsiPM" , "deltaPsiPM" , "@0 - @1", RooArgList(mumuMass,PsiPMass) )
# deltaB0M     = fulldata.addColumn(deltaB0Mfunc) ;
# deltaJpsiM   = fulldata.addColumn(deltaJMfunc) ;
# deltaPsiPM   = fulldata.addColumn(deltaPMfunc) ;

print 'reading monte carlo...'
genSignal       = RooRealVar("genSignal"      , "genSignal"      , 0, 10);
thevarsMC   = thevars; 
thevarsMC.add(genSignal)
fullmc      = RooDataSet('fullmc', 'fullmc', tMC,  RooArgSet(thevarsMC))

nMCEntries = fullmc.sumEntries()
rDataset = p0.generate(RooArgSet(randVar),int(nMCEntries))
fullmc.merge(rDataset) 
thevarsMC.add(randVar)
# deltaB0M    = fullmc.addColumn(deltaB0Mfunc) 
# deltaJpsiM  = fullmc.addColumn(deltaJMfunc)  
# deltaPsiPM  = fullmc.addColumn(deltaPMfunc)  
# 
# thevars.add(deltaB0M)
# thevars.add(deltaJpsiM)
# thevars.add(deltaPsiPM)
# thevarsMC.add(deltaB0M)
# thevarsMC.add(deltaJpsiM)
# thevarsMC.add(deltaPsiPM)

### define correct and wrong tag samples
rt_mc       = fullmc.reduce(RooArgSet(thevarsMC), '((tagB0==1 && genSignal==1) || (tagB0==0 && genSignal==2))')
wt_mc       = fullmc.reduce(RooArgSet(thevarsMC), '((tagB0==0 && genSignal==1) || (tagB0==1 && genSignal==2))')

dict_s_rt  = {}
dict_s_wt  = {}

out_fname = "fit_results_mass_checkOnMC/newbdt_puw/results_fits_%s_fM%s_asInANv8.root"%(args.year, '_Jpsi'*(args.dimusel=='keepJpsi') + '_Psi'*(args.dimusel=='keepPsiP'))
out_f = TFile (out_fname, "RECREATE") 
w = ROOT.RooWorkspace("w")
initial_n_1 =  3.
initial_n_2 =  1.
initial_a_1 =  1.
initial_a_2 = -1.
initial_sigma1 = 0.028
initial_sigma2 = 0.048
initial_sigmaCB = 0.048


for ibin in range(len(q2binning)-1):
    print 'dimuon selection: ', args.dimusel
    if args.dimusel == 'rejectPsi' and \
       (q2binning[ibin] == 8.68 or q2binning[ibin] == 12.86): 
           continue
    if args.dimusel == 'keepJpsi' and \
       (q2binning[ibin] < 8.68 or q2binning[ibin] > 10): 
           continue
    if args.dimusel == 'keepPsiP' and \
       (q2binning[ibin] < 12 or q2binning[ibin] > 14): 
           continue
    if   ibin ==4:  cut_base = 'passB0Psi_jpsi== 1 '
    elif ibin ==6:  cut_base = 'passB0Psi_psip== 1 '

    nRT_fromMC = fitMC(rt_mc, True, ibin)
    nWT_fromMC = fitMC(wt_mc, False, ibin)
#    fitMC(rt_mc, True, ibin)
#    fitMC(wt_mc, False, ibin)
    fitData(fulldata, ibin, nRT_fromMC, nWT_fromMC)
    print ' --------------------------------------------------------------------------------------------------- '

#create Tgraph
for ibin in range(len(q2binning)-1):

#   yields
    middle_value = (q2binning[ibin+1] +  q2binning[ibin]) / 2

#   efficiencies
    MC_num_entries[ibin] = MC_RT_num_entries[ibin] + MC_WT_num_entries[ibin]
    MC_num_entries_error[ibin] = math.sqrt(MC_num_entries[ibin])

#   q2bins
    tgraph_x[ibin] = middle_value
    tgraph_ex[ibin] = middle_value - q2binning[ibin]


#arr_tgraph_y = array.array('d',tgraph_y)
#arr_tgraph_x = array.array('d',tgraph_x)
#arr_tgraph_ex = array.array('d',tgraph_ex)
#arr_tgraph_y = array.array('d',tgraph_y)

tgraph_yield = TGraphErrors(8, np.array(tgraph_x), np.array(tgraph_y), np.array(tgraph_ex), np.array(tgraph_ey))
#tgraph_yield = TGraphErrors(8, tgraph_x, tgraph_y, tgraph_ex, tgraph_ey)

tgraph_efficiency = TGraphErrors(8, np.array(tgraph_x), np.array(MC_num_entries), np.array(tgraph_ex), np.array(MC_num_entries_error))

print 'tgraph_y'
print(tgraph_y)

print 'MC_num_entries'
print(MC_num_entries)

print 'MC_num_entries_error'
print(MC_num_entries_error)

c2 = ROOT.TCanvas() 
tgraph_yield.SetName("Yields")
tgraph_yield.Draw("AP")
c2.SaveAs("tgraph_yields1.pdf")

print 'tgraph_yield created'
print (tgraph_yield)

c3 = ROOT.TCanvas() 
tgraph_efficiency.SetName("Efficiencies_num")
tgraph_efficiency.Draw("AP")
c3.SaveAs("tgraph_efficiencies.pdf")

print 'tgraph_efficiency created'
print (tgraph_efficiency)

out_gname = "fit_results_mass_checkOnMC/newbdt_puw/yields_efficiencies_%s%s.root"%(args.year, '_Jpsi'*(args.dimusel=='keepJpsi') + '_Psi'*(args.dimusel=='keepPsiP'))
out_g = TFile (out_gname, "RECREATE") 
tgraph_yield.Write("tgraph_yields_%s%s"%(args.year, '_Jpsi'*(args.dimusel=='keepJpsi') + '_Psi'*(args.dimusel=='keepPsiP')))
tgraph_efficiency.Write("tgraph_efficiencies_%s%s"%(args.year, '_Jpsi'*(args.dimusel=='keepJpsi') + '_Psi'*(args.dimusel=='keepPsiP')))

w.Print()
print '--------------------------------------------------------------------------------------------------- '
print 'bin\t\t fit status \t cov. matrix \t\t chi2'
# for i,k in enumerate(fitStats.keys()):    
#     if i%3==0:  print '------------------------------------------------------'
#     print k , '\t\t', fitStats[k], '\t\t', covStats[k], '\t\t', chi2s[k]
# 
# for i,k in enumerate(dict_s_rt.keys()):    
#     print 'bin ' , i , '\t\t n rt events: ', dict_s_rt[k], '\t\t n wt events: ', dict_s_wt[k], '\t\t n tot events: ', dict_s_wt[k]+dict_s_rt[k]

print 'It took', time.time() - T, 'seconds'
print (time.time()- T)/60, 'minutes'
print (time.time() - T)/3600, 'hours'

out_g.Close()
out_f.Close()
w.writeToFile(out_f.GetName(), False)
