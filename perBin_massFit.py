import argparse
parser = argparse.ArgumentParser(description="")
# parser.add_argument("inputfile" , help = "Path to the input ROOT file")
parser.add_argument("dimusel"   , help = "Define if keep or remove dimuon resonances. You can choose: keepPsiP, keepJpsi, rejectPsi, keepPsi")
parser.add_argument("year"      , help = "choose among:2016,2017,2018", default = '2018')
args = parser.parse_args()


'''
code to fit the B0 mass distribution:
- unbinned fit
- possibility to apply cuts on the dimuon mass [B0&Psi cut in RunI analysis] (e.g. to exclude the Jpsi mass region, or the psi) via the parameter dimusel
'''

import os, sys
from os import path
sys.path.insert(0, os.environ['HOME'] + '/.local/lib/python2.7/site-packages')

import ROOT
from ROOT import gSystem
ROOT.gROOT.SetBatch(True)

gSystem.Load('libRooFit')
from ROOT import RooFit, RooRealVar, RooDataSet, RooArgList, RooTreeData, RooArgSet, RooAddPdf, RooFormulaVar
from ROOT import RooGaussian, RooExponential, RooChebychev, RooProdPdf, RooCBShape, TFile, RooPolynomial
import sys, math
from uncertainties import ufloat
import random

ROOT.RooMsgService.instance().setGlobalKillBelow(4)
ROOT.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(50000)


def _getFittedVar(varName, w=None):
    if w is not None:
        return ufloat (w.var(varName).getVal() , w.var(varName).getError())
    else :
        return ufloat (varName.getVal()        , varName.getError())

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
    
def _constrainVar(var):
    
    constr = _getFittedVar(var.GetName(), w)
    gauss_constr = RooGaussian(  "c_%s" %var.GetName() , 
                                 "c_%s" %var.GetName() , 
                                var         ,  
                                ROOT.RooFit.RooConst( constr.n ), 
                                ROOT.RooFit.RooConst( constr.s )
                                ) 
    print 'constraining var',   var.GetName(), ': ',     constr.n , ' with uncertainty ' , constr.s                          
    return gauss_constr                        


from utils.utils import *
from utils.fit_functions import *


nSigma_psiRej = 3.
cut_base      = applyB0PsiCut(args.dimusel, nSigma_psiRej)

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
                19,
]







def fitMC(fulldata, correctTag, ibin):

    print 'now fitting: ', ibin, ' for ', correctTag*'correctTag ', (1-correctTag)*'wrongTag'  
    cut = cut_base + '&& (mumuMass*mumuMass > %s && mumuMass*mumuMass < %s)'%(q2binning[ibin], q2binning[ibin+1])
    data        = fulldata.reduce(RooArgSet(thevarsMC), cut)

    pol_c1      = RooRealVar ("p1"           , "coeff x^0 term" ,  -0.5,   -10, 10);
    bkg_pol     = RooChebychev("bkg_pol"     , "bkg_pol" ,  tagged_mass, RooArgList(pol_c1));
    signalFunction = bkg_pol ### just a placeholder

    nsig        = RooRealVar("Yield"         , "nsig"   ,   10000,     0,    1000000)
    nbkg        = RooRealVar("nbkg"          , "nbkg"   ,      10,     0,    100000 )
    
    doextended = False
    fitrange   = "mcrange"

    if correctTag:
        doubleG( B0Mass_ , initial_sigma1 , initial_sigma2,  0.8, tagged_mass, w, "RT%s"%ibin)    ## (mean_   , sigma1_, sigma2_, f1_)
        signalFunction = w.pdf("doublegaus_RT%s"%ibin)   
        fitFunction    = RooAddPdf ("fitfunction" , "fit function"  ,  RooArgList(signalFunction, bkg_pol), RooArgList(nsig, nbkg))
        doextended = True
        fitrange   = "full"

    else:
        mean         = RooRealVar ("mean^{WT%s}"%ibin, "massWT",   B0Mass_    ,     5,    6, "GeV")
        crystalBall(  mean, initial_sigmaCB     , initial_a_1   ,  initial_n_1 , tagged_mass, w, "1", "WT%s"%ibin, [0, 10])    ## (mean    , sigma_ , alpha_ ,  n_)
        crystalBall(  mean, initial_sigmaCB     , initial_a_2   ,  initial_n_2 , tagged_mass, w, "2", "WT%s"%ibin, [-10,0])    
        doubleCB (    w.pdf("cbshape_1_WT%s"%ibin) , w.pdf("cbshape_2_WT%s"%ibin)  , 0.8  ,  tagged_mass, w, "%s"%ibin    )
        signalFunction = w.pdf("doublecb_%s"%ibin)   
        fitFunction    = signalFunction

    r = fitFunction.fitTo(data, RooFit.Extended(doextended), RooFit.Save(), RooFit.Range(fitrange))
    print 'fit status: ', r.status(), r.covQual() 
    
    if not _goodFit(r):
        r = fitFunction.fitTo(data, RooFit.Extended(doextended), RooFit.Save(), RooFit.Range(fitrange))
        print 'fit status (redo): ', r.status(), r.covQual() 

    if not _goodFit(r) and correctTag:
        r = fitFunction.fitTo(data, RooFit.Extended(doextended), RooFit.Save(), RooFit.Range(fitrange))
        print 'fit status (redo2): ', r.status(), r.covQual() 

    if not correctTag and (not _goodFit(r) or \
                                w.var("n_{1}^{WT%s}"%ibin).getVal() > 12 or \
                                w.var("n_{2}^{WT%s}"%ibin).getVal() > 12 or \
                                w.var("n_{1}^{WT%s}"%ibin).getVal() < 0.001 or \
                                w.var("n_{2}^{WT%s}"%ibin).getVal() < 0.001\
                                ):
        mean = RooRealVar ("mean^{WT%s}"%ibin, "massWT",   B0Mass_, 5, 6, "GeV")
        singleG( mean , 0.028 , tagged_mass, w, "WT%s"%ibin)   
        gausCB ( w.pdf("cbshape_1_WT%s"%ibin) , w.pdf("gaus_WT%s"%ibin) , 0.3  ,  tagged_mass, w, "%s"%ibin    ) 
        signalFunction = w.pdf("gauscb_%s"%ibin)   
        fitFunction    = signalFunction
        r = fitFunction.fitTo(data, RooFit.Extended(doextended), RooFit.Save(), RooFit.Range(fitrange))
        print 'fit status (gaus + cb): ', r.status(), r.covQual() 

        
            
    params = signalFunction.getParameters(RooArgSet(tagged_mass)) 
    w.saveSnapshot("reference_fit_%s_%s"%('RT'*correctTag + 'WT'*(1-correctTag), ibin),params,ROOT.kTRUE) 
    
    frame = tagged_mass.frame(RooFit.Range(fitrange))
    data.plotOn(frame, RooFit.Binning(60), RooFit.MarkerSize(.7))
    
    drawPdfComponents(fitFunction, frame, ROOT.kGreen if correctTag else ROOT.kViolet, RooFit.NormRange(fitrange), RooFit.Range(fitrange), isData=False)
    fitFunction.plotOn(frame, RooFit.NormRange(fitrange), RooFit.Range(fitrange) )
    
##     parList = RooArgSet (nsig,sigma,sigma2, sigmaCB, mean)
##     fitFunction.paramOn(frame, RooFit.Parameters(parList), RooFit.Layout(0.62,0.86,0.88))
    fitFunction.paramOn(frame,  RooFit.Layout(0.62,0.86,0.88))

    frame.Draw()
    niceFrame(frame, '')
    frame. addObject(_writeFitStatus(r))
    frame.Draw()

    if correctTag:
        dict_s_rt[ibin]   = _getFittedVar(nsig)
    else:
        dict_s_wt[ibin]    = ufloat(data.sumEntries(), math.sqrt(data.sumEntries()))
    

    frame.SetTitle('correctly'*correctTag + 'wrongly'*(1-correctTag) + ' tagged events')
    c1.SaveAs('fit_results_mass/save_fit_mc_%s_%s_%sT.pdf'%(ibin, args.year, "R"*correctTag + "W"*(1-correctTag)))
    out_f.cd()
    r.Write('results_%s_%s'%(correctTag*'RT' + (1-correctTag)*'WT', ibin))

   
   
def fitData(fulldata, ibin):

    cut  = cut_base + '&& (mumuMass*mumuMass > %s && mumuMass*mumuMass < %s)'%(q2binning[ibin], q2binning[ibin+1])
    data = fulldata.reduce(RooArgSet(tagged_mass,mumuMass,mumuMassE), cut)

    fraction = dict_s_rt[ibin] / (dict_s_rt[ibin] + dict_s_wt[ibin])
    print 'mistag fraction on MC for bin ', ibin , ' : ' , fraction.n , '+/-', fraction.s 
    
    ### creating RT component
    w.loadSnapshot("reference_fit_RT_%s"%ibin)
    sigmart1    = w.var("#sigma_{1}^{RT%s}"%ibin  )
    sigmart2    = w.var("#sigma_{2}^{RT%s}"%ibin  )
    massrt      = w.var("mean^{RT%s}"%ibin  )
    f1rt        = w.var("f^{RT%s}"%ibin)

    theRTgauss  = w.pdf("doublegaus_RT%s"%ibin)   
    c_sigma_rt1 = _constrainVar(sigmart1)
    c_sigma_rt2 = _constrainVar(sigmart2)
    c_mean_rt   = _constrainVar(massrt)
    c_f1rt      = _constrainVar(f1rt)

    ### creating WT component
    w.loadSnapshot("reference_fit_WT_%s"%ibin)

    sigmawt1    = w.var("#sigma_{CB1}^{WT%s}"%ibin)
    alphawt1    = w.var("#alpha_{1}^{WT%s}"%ibin)
    nwt1        = w.var("n_{1}^{WT%s}"%ibin)
    sigmawt2    = w.var("#sigma^{WT%s}"%ibin)
    f3wt        = w.var("f^{WT%s}"%ibin)
    meanwt      = w.var("mean^{WT%s}"%ibin)

    theWTgauss  = w.pdf("gauscb_%s"%ibin)   
    if theWTgauss==None:
        theWTgauss = w.pdf("doublecb_%s"%ibin)   
        sigmawt2   = w.var("#sigma_{CB2}^{WT%s}"%ibin)
        alphawt2   = w.var("#alpha_{2}^{WT%s}"%ibin)
        nwt2       = w.var("n_{2}^{WT%s}"%ibin)
        
        c_alpha_wt2   = _constrainVar(alphawt2)
        c_n_wt2       = _constrainVar(nwt2)

    c_sigma_wt1   = _constrainVar(sigmawt1)
    c_alpha_wt1   = _constrainVar(alphawt1)
    c_n_wt1       = _constrainVar(nwt1)
    c_sigma_wt2   = _constrainVar(sigmawt2)
    c_f3          = _constrainVar(f3wt)
    c_mean_wt     = _constrainVar(meanwt)
    

    ### creating constraints for the RT component
    c_RTgauss   = RooProdPdf  ("c_RTgauss" , "c_RTgauss" , RooArgList(theRTgauss, c_sigma_rt1, c_sigma_rt2, c_mean_rt, c_f1rt  ) )     

    c_vars = RooArgSet(c_sigma_rt1, c_sigma_rt2, c_f1rt, c_mean_rt)
    c_vars.add(c_sigma_wt1)
    c_vars.add(c_sigma_wt2)
    c_vars.add(c_mean_wt)
    c_vars.add(c_alpha_wt1)
    c_vars.add(c_n_wt1)
    c_vars.add(c_f3)

    ### creating constraints for the WT component
    if theWTgauss.GetName() == "gauscb_%s"%ibin:
        c_WTgauss  = RooProdPdf  ("c_WTgauss" , "c_WTgauss" , RooArgList(theWTgauss, c_sigma_wt1, c_alpha_wt1, c_n_wt1, c_sigma_wt2, c_f3, c_mean_wt  ) )     
    else:
        c_WTgauss  = RooProdPdf  ("c_WTgauss" , "c_WTgauss" , RooArgList(theWTgauss, c_sigma_wt1, c_alpha_wt1, c_n_wt1, c_sigma_wt2, c_f3, c_mean_wt, c_alpha_wt2, c_n_wt2  ) )     
        c_vars.add(c_alpha_wt2)
        c_vars.add(c_n_wt2)


    frt              = RooRealVar ("F_{RT}"          , "frt"             , fraction.n , 0, 1)
    signalFunction   = RooAddPdf  ("sumgaus"         , "rt+wt"           , RooArgList(c_RTgauss,c_WTgauss), RooArgList(frt))
    c_frt            = RooGaussian("c_frt"           , "c_frt"           , frt,  ROOT.RooFit.RooConst(fraction.n) , ROOT.RooFit.RooConst(fraction.s) )
    c_signalFunction = RooProdPdf ("c_signalFunction", "c_signalFunction", RooArgList(signalFunction, c_frt))     
    c_vars.add(frt)

    ### now create background parametrization
    slope         = RooRealVar    ("slope"      , "slope"           ,    0.5,   -10, 10);
    bkg_exp       = RooExponential("bkg_exp"    , "exponential"     ,  slope,   tagged_mass  );
    pol_c1        = RooRealVar    ("p1"         , "coeff x^0 term"  ,    0.5,   -10, 10);
    pol_c2        = RooRealVar    ("p2"         , "coeff x^1 term"  ,    0.5,   -10, 10);
    bkg_pol       = RooChebychev  ("bkg_pol"    , "2nd order pol"   ,  tagged_mass, RooArgList(pol_c1,pol_c2));
   
    nsig          = RooRealVar("Yield"         , "signal frac"    ,    4000,     0,   1000000);
    nbkg          = RooRealVar("nbkg"          , "bkg fraction"   ,    1000,     0,   550000);
    

#     fitFunction = RooAddPdf ("fitfunction" , "fit function"  ,  RooArgList(signalFunction, bkg_pol), RooArgList(nsig, nbkg))
    fitFunction = RooAddPdf ("fitfunction" , "fit function"  ,  RooArgList(c_signalFunction, bkg_exp), RooArgList(nsig, nbkg))

    r = fitFunction.fitTo(data, 
                          RooFit.Extended(True), 
                          RooFit.Save(), 
                          RooFit.Range("full"), 
                          RooFit.Verbose(False),
                          ROOT.RooFit.Constrain(c_vars)
                         )

    frame = tagged_mass.frame( RooFit.Range("full") )
    data.plotOn(frame, RooFit.Binning(35), RooFit.MarkerSize(.7))
    fitFunction.plotOn(frame);
    drawPdfComponents(fitFunction, frame, ROOT.kAzure, RooFit.NormRange("full"), RooFit.Range("full"), isData = True)

    fitFunction.paramOn(frame,  RooFit.Layout(0.62,0.86,0.88))
    frame.Draw()
    niceFrame(frame, '')
    frame. addObject(_writeFitStatus(r))

    if not args.year=='test':  writeCMS(frame, args.year, [ q2binning[ibin], q2binning[ibin+1] ])
    frame.Draw()
    c1.SaveAs('fit_results_mass/save_fit_data_%s_%s_LMNR.pdf'%(ibin, args.year))








tData = ROOT.TChain('ntuple')
tMC = ROOT.TChain('ntuple')

if args.year == 'test':
    tData.Add('/gwteray/users/fiorendi/final_ntuples_p5prime_allyears/2016Data_100k.root')
    tMC.Add('/gwteray/users/fiorendi/final_ntuples_p5prime_allyears/2016MC_LMNR_100k.root')
else:    
    tData.Add('/gwteray/users/fiorendi/final_ntuples_p5prime_allyears/%sData_All_finalSelection.root'%args.year)
    tMC.Add('/gwteray/users/fiorendi/final_ntuples_p5prime_allyears/%sMC_LMNR.root'%args.year)


tagged_mass     = RooRealVar("tagged_mass" , "#mu^{+}#mu^{-}K#pi mass", 4.9, 5.6, "GeV")
mumuMass        = RooRealVar("mumuMass"    , "mumuMass" , 0, 6);
mumuMassE       = RooRealVar("mumuMassE"   , "mumuMassE", 0, 10000);
tagB0           = RooRealVar("tagB0"       , "tagB0"    , 0, 2);

tagged_mass.setRange("full",   5. ,5.6) ;
tagged_mass.setRange("mcrange",4.9,5.6) ;
thevars = RooArgSet()
thevars.add(tagged_mass)
thevars.add(mumuMass)
thevars.add(mumuMassE)
thevars.add(tagB0)

fulldata   = RooDataSet('fulldata', 'fulldataset', tData,  RooArgSet(thevars))



## add to the input tree the combination of the variables, to be used for the cuts on the dimuon mass
deltaB0Mfunc = RooFormulaVar("deltaB0M", "deltaB0M", "@0 - @1", RooArgList(tagged_mass,B0Mass) )
deltaB0M     = fulldata.addColumn(deltaB0Mfunc) ;
deltaJMfunc  = RooFormulaVar("deltaJpsiM" , "deltaJpsiM" , "@0 - @1", RooArgList(mumuMass,JPsiMass) )
deltaJpsiM   = fulldata.addColumn(deltaJMfunc) ;
deltaPMfunc  = RooFormulaVar("deltaPsiPM" , "deltaPsiPM" , "@0 - @1", RooArgList(mumuMass,PsiPMass) )
deltaPsiPM   = fulldata.addColumn(deltaPMfunc) ;

genSignal       = RooRealVar("genSignal"      , "genSignal"      , 0, 10);
thevarsMC   = thevars; 
thevarsMC.add(genSignal)
fullmc      = RooDataSet('fullmc', 'fullmc', tMC,  RooArgSet(thevarsMC))

deltaB0M    = fullmc.addColumn(deltaB0Mfunc) 
deltaJpsiM  = fullmc.addColumn(deltaJMfunc)  
deltaPsiPM  = fullmc.addColumn(deltaPMfunc)  

thevars.add(deltaB0M)
thevars.add(deltaJpsiM)
thevars.add(deltaPsiPM)

thevarsMC.add(deltaB0M)
thevarsMC.add(deltaJpsiM)
thevarsMC.add(deltaPsiPM)

### define correct and wrong tag samples
rt_mc       = fullmc.reduce(RooArgSet(thevarsMC), '((tagB0==1 && genSignal==1) || (tagB0==0 && genSignal==2))')
wt_mc       = fullmc.reduce(RooArgSet(thevarsMC), '((tagB0==0 && genSignal==1) || (tagB0==1 && genSignal==2))')


c1 = ROOT.TCanvas() 
dict_s_rt  = {}
dict_s_wt  = {}

out_f = TFile ("results_fits_%s.root"%args.year,"RECREATE") 

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
           
    fitMC(rt_mc, True, ibin)
    fitMC(wt_mc, False, ibin)
    fitData(fulldata, ibin)

