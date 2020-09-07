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
gSystem.Load('utils/func_roofit/libRooDoubleCBFast')
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

def _writeChi2(chi2):
    txt = ROOT.TLatex(.16,.6, "fit #chi^{2}: %.1f "%chi2 )
    txt . SetNDC() ;
    txt . SetTextSize(0.033) ;
    txt . SetTextFont(42)
    return txt
    
def _constrainVar(var, nsigma):
    
    constr = _getFittedVar(var.GetName(), w)
    gauss_constr = RooGaussian(  "c_%s" %var.GetName() , 
                                 "c_%s" %var.GetName() , 
                                var         ,  
                                ROOT.RooFit.RooConst( constr.n ), 
                                ROOT.RooFit.RooConst( nsigma*constr.s )
                                ) 
    print 'constraining var',   var.GetName(), ': ',     constr.n , ' with uncertainty ' , nsigma*constr.s                          
    return gauss_constr                        

from utils.utils import *
from utils.fit_functions import *

nbins = 60
nSigma_psiRej = 3.
cut_base      = applyB0PsiCut(args.dimusel, nSigma_psiRej)
## entries per bin in data
n_bin = [520, 1000, 850, 1900, 0, 3166, 0, 1860]

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




def fitMC(fulldata, correctTag, ibin):

    print 'now fitting: ', ibin, ' for ', correctTag*'correctTag ', (1-correctTag)*'wrongTag'  
    cut = cut_base + '&& (mumuMass*mumuMass > %s && mumuMass*mumuMass < %s)'%(q2binning[ibin], q2binning[ibin+1])
    data        = fulldata.reduce(RooArgSet(thevarsMC), cut)

    pol_c1      = RooRealVar ("p1"           , "coeff x^0 term" ,  -0.5,   -10, 10);
    pol_c2      = RooRealVar ("p2"           , "coeff x^0 term" ,  0.5,   -10, 10);
    bkg_pol     = RooChebychev("bkg_pol"     , "bkg_pol" ,  tagged_mass, RooArgList(pol_c1, pol_c2));
    signalFunction = bkg_pol ### just a placeholder

    nsig        = RooRealVar("Yield"         , "nsig"   ,   500000,     0,    1E7)
    nbkg        = RooRealVar("nbkg"          , "nbkg"   ,     1000,     0,    1E6 )
    
    doextended = False
    fitrange   = "mcrange"
    nbins = 70

    if correctTag:
        doubleG( B0Mass_ , initial_sigma1 , initial_sigma2,  0.8, tagged_mass, w, "RT%s"%ibin)    ## (mean_   , sigma1_, sigma2_, f1_)
        signalFunction = w.pdf("doublegaus_RT%s"%ibin)   
        fitFunction    = RooAddPdf ("fitfunction" , "fit function"  ,  RooArgList(signalFunction, bkg_pol), RooArgList(nsig, nbkg))
        doextended = True
        fitrange   = "full"
        nbins = 60

    else:
        mean        = RooRealVar ("mean^{WT%s}"%ibin,        "massWT"         , B0Mass_,     5,    6, "GeV")
        sigmaCB     = RooRealVar ("#sigma_{CB}^{WT%s}"%ibin, "sigmaCB"        ,  0.03  ,    0,   1  )
        alpha1      = RooRealVar ("#alpha_{1}^{WT%s}"%ibin,  "#alpha_{1}"     ,  0.5   ,    0,  10  )
        alpha2      = RooRealVar ("#alpha_{2}^{WT%s}"%ibin,  "#alpha_{2}"     ,  2.5   ,    0,  10  )
        n1          = RooRealVar ("n_{1}^{WT%s}"%ibin,       "n_1"            ,  1     ,    0,  20  )
        n2          = RooRealVar ("n_{2}^{WT%s}"%ibin,       "n_2"            ,  1     ,    0,  20  )
        doublecb = ROOT.RooDoubleCBFast("doublecb_%s"%ibin, "doublecb", tagged_mass, mean, sigmaCB, alpha1, n1, alpha2, n2)	
#         doublecb      = RooCBShape ("doublecb_%s"%(bin)       , "cbshape_%s"%(bin) ,  tagged_mass, mean, sigmaCB, alpha1, n1)
#         getattr(w, 'import')(doublecb)
        signalFunction = doublecb
        fitFunction    = doublecb
        
    getattr(w,"import")(signalFunction)

    r = fitFunction.fitTo(data, RooFit.Extended(doextended), RooFit.Save(), RooFit.Range(fitrange))
    print 'fit status: ', r.status(), r.covQual() 
    
    
    if not _goodFit(r):
        r = fitFunction.fitTo(data, RooFit.Extended(doextended), RooFit.Save(), RooFit.Range(fitrange))
        print 'fit status (redo): ', r.status(), r.covQual() 

    if not _goodFit(r) and correctTag:
        r = fitFunction.fitTo(data, RooFit.Extended(doextended), RooFit.Save(), RooFit.Range(fitrange))
        print 'fit status (redo2): ', r.status(), r.covQual() 

         
    params = signalFunction.getParameters(RooArgSet(tagged_mass)) 
    w.saveSnapshot("reference_fit_%s_%s"%('RT'*correctTag + 'WT'*(1-correctTag), ibin),params,ROOT.kTRUE) 
    frame = tagged_mass.frame(RooFit.Range(fitrange))
    print 'n bins now is ', nbins, '(',correctTag, ')' 
    data.plotOn(frame, RooFit.Binning(nbins), RooFit.MarkerSize(.7))
    
    drawPdfComponents(fitFunction, frame, ROOT.kGreen if correctTag else ROOT.kViolet, RooFit.NormRange(fitrange), RooFit.Range(fitrange), isData=False)
    fitFunction.plotOn(frame, RooFit.NormRange(fitrange), RooFit.Range(fitrange) )
    fitFunction.paramOn(frame,  RooFit.Layout(0.62,0.86,0.88))

    frame.Draw()
    niceFrame(frame, '')
    frame. addObject(_writeFitStatus(r))
    
    ## evaluate sort of chi2 and save number of RT/WT events
    observables = RooArgSet(tagged_mass)
    flparams = fitFunction.getParameters(observables)
    nparam = int(flparams.selectByAttrib("Constant",ROOT.kFALSE).getSize())
    if correctTag:
        frame. addObject(_writeChi2( frame.chiSquare("fitfunction_Norm[tagged_mass]_Range[full]_NormRange[full]", "h_fullmc",  nparam) ))
        dict_s_rt[ibin]   = _getFittedVar(nsig)
        nRT = RooRealVar ("nRT_%s"%ibin, "yield of RT signal",0,1.E6)
        nRT.setVal(  dict_s_rt[ibin].n)
        nRT.setError(dict_s_rt[ibin].s)
        print 'setting nRT to ', dict_s_rt[ibin].n
        nRT_forLater = dict_s_rt[ibin].n
        getattr(w,"import")(nRT)
    else:
        frame. addObject(_writeChi2( frame.chiSquare("doublecb_%s_Norm[tagged_mass]_Comp[doublecb_%s]_Range[mcrange]_NormRange[mcrange]"%(ibin,ibin), "h_fullmc",  nparam) ))
        dict_s_wt[ibin]    = ufloat(data.sumEntries(), math.sqrt(data.sumEntries()))
        nWT = RooRealVar ("nWT_%s"%ibin, "yield of WT signal",0,1.E6)
        nWT.setVal(  dict_s_wt[ibin].n)
        nWT.setError(dict_s_wt[ibin].s)
        nWT_forLater = dict_s_wt[ibin].n
        print 'setting nWT to ', dict_s_wt[ibin].n
        print 'nWT events in data: ', data.sumEntries()
        getattr(w,"import")(nWT)

#         chi2 = frame.chiSquare("doublecb_%s_Norm[tagged_mass]_Comp[doublecb_%s]_Range[mcrange]_NormRange[mcrange]"%(ibin,ibin), "h_fullmc",  nparam)
#         if chi2 == -1:
#             chi2 = frame.chiSquare("gauscb_%s_Norm[tagged_mass]_Comp[gauscb_%s]_Range[mcrange]_NormRange[mcrange]"%(ibin,ibin), "h_fullmc",  nparam) 
#         frame. addObject(_writeChi2( chi2 ) )
            
    frame.Draw()
    frame.SetTitle('correctly'*correctTag + 'wrongly'*(1-correctTag) + ' tagged events')
    c1.SaveAs('fit_results_mass_checkOnMC/save_fit_mc_%s_%s_%sT_pol2.pdf'%(ibin, args.year, "R"*correctTag + "W"*(1-correctTag)))
    c1.SetLogy()
    c1.SaveAs('fit_results_mass_checkOnMC/save_fit_mc_%s_%s_%sT_pol2_log.pdf'%(ibin, args.year, "R"*correctTag + "W"*(1-correctTag)))
    c1.SetLogy(False)
    out_f.cd()
    r.Write('results_%s_%s'%(correctTag*'RT' + (1-correctTag)*'WT', ibin))
    
    if correctTag:
      return nRT_forLater
    else:
      return nWT_forLater
   
   
def fitData(fulldata, ibin, nRT_forLater, nWT_forLater):

    cut  = cut_base + '&& (mumuMass*mumuMass > %s && mumuMass*mumuMass < %s)'%(q2binning[ibin], q2binning[ibin+1])
    fulldata_v2 = fulldata.reduce(RooArgSet(tagged_mass,mumuMass,mumuMassE, randVar), cut)

#     import pdb; pdb.set_trace()
    nDataEntries = fulldata_v2.sumEntries()
    nDesired = n_bin[ibin]/nDataEntries
    cut = 'rand < %f'%nDesired
    data = fulldata_v2.reduce(RooArgSet(tagged_mass,mumuMass,mumuMassE), cut)


    fraction = dict_s_rt[ibin] / (dict_s_rt[ibin] + dict_s_wt[ibin])
    print 'mistag fraction on MC for bin ', ibin , ' : ' , fraction.n , '+/-', fraction.s 
    
    ### creating RT component
    w.loadSnapshot("reference_fit_RT_%s"%ibin)
    sigmart1    = w.var("#sigma_{1}^{RT%s}"%ibin  )
    sigmart2    = w.var("#sigma_{2}^{RT%s}"%ibin  )
    massrt      = w.var("mean^{RT%s}"%ibin  )
    f1rt        = w.var("f^{RT%s}"%ibin)

    theRTgauss  = w.pdf("doublegaus_RT%s"%ibin)   
    c_sigma_rt1 = _constrainVar(sigmart1, 5)
    c_sigma_rt2 = _constrainVar(sigmart2, 5)
    c_mean_rt   = _constrainVar(massrt, 3)
    c_f1rt      = _constrainVar(f1rt, 3)

    ### creating WT component
    w.loadSnapshot("reference_fit_WT_%s"%ibin)
    meanwt      = w.var("mean^{WT%s}"%ibin)
    sigmawt     = w.var("#sigma_{CB}^{WT%s}"%ibin)
    alphawt1    = w.var("#alpha_{1}^{WT%s}"%ibin)
    alphawt2    = w.var("#alpha_{2}^{WT%s}"%ibin)
    nwt1        = w.var("n_{1}^{WT%s}"%ibin)
    nwt2        = w.var("n_{2}^{WT%s}"%ibin)

    theWTgauss  = w.pdf("doublecb_%s"%ibin)   
    c_mean_wt     = _constrainVar(meanwt, 3)
    c_sigma_wt    = _constrainVar(sigmawt, 3)
    c_alpha_wt1   = _constrainVar(alphawt1, 3)
    c_alpha_wt2   = _constrainVar(alphawt2, 3)
    c_n_wt1       = _constrainVar(nwt1, 3)
    c_n_wt2       = _constrainVar(nwt2, 3)


    ### creating constraints for the RT component
#     c_RTgauss   = RooProdPdf  ("c_RTgauss" , "c_RTgauss" , RooArgList(theRTgauss, c_mean_rt  ) )     
    c_RTgauss   = RooProdPdf  ("c_RTgauss" , "c_RTgauss" , RooArgList(theRTgauss, c_sigma_rt1, c_sigma_rt2, c_mean_rt, c_f1rt  ) )     

#     c_vars = RooArgSet(c_mean_rt)
    c_vars = RooArgSet(c_sigma_rt1, c_sigma_rt2, c_f1rt, c_mean_rt)
    c_vars.add(c_sigma_wt)
    c_vars.add(c_mean_wt)
    c_vars.add(c_alpha_wt1)
    c_vars.add(c_alpha_wt2)
    c_vars.add(c_n_wt1)
    c_vars.add(c_n_wt2)

    ### creating constraints for the WT component
#     c_WTgauss  = RooProdPdf  ("c_WTgauss" , "c_WTgauss" , RooArgList(theWTgauss, c_mean_wt  ) )     
#     import pdb; pdb.set_trace()
    c_WTgauss  = RooProdPdf  ("c_WTgauss" , "c_WTgauss" , RooArgList(theWTgauss, c_alpha_wt1, c_n_wt1, c_sigma_wt, c_mean_wt, c_alpha_wt2, c_n_wt2  ) )     

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
    bkg_pol       = RooChebychev  ("bkg_pol"    , "2nd order pol"   ,  tagged_mass, RooArgList(pol_c1));
   
    nsig          = RooRealVar("Yield"         , "signal frac"    ,     nRT_forLater + nWT_forLater,     0,   1000000);
    nbkg          = RooRealVar("nbkg"          , "bkg fraction"   ,     1000,     0,   550000);
#     nsig.setConstant(True)

    print nsig.getVal()
    fitFunction = c_signalFunction
#     fitFunction = RooAddPdf ("fitfunction" , "fit function"  ,  RooArgList(c_signalFunction, bkg_pol), RooArgList(nsig, nbkg))
#     fitFunction = RooAddPdf ("fitfunction" , "fit function"  ,  RooArgList(c_signalFunction, bkg_exp), RooArgList(nsig, nbkg))
    
#     meanwt  .setConstant(True)
#     sigmawt .setConstant(True)
#     alphawt1.setConstant(True)
#     alphawt2.setConstant(True)
#     nwt1    .setConstant(True)
#     nwt2    .setConstant(True)
#     fitFunction.fitTo(data, 
#                               RooFit.Extended(True), 
#                               RooFit.Range("full"), 
#                               RooFit.Verbose(False),
#                               ROOT.RooFit.Constrain(c_vars)
#                          )     
#     meanwt  .setConstant(False)
#     sigmawt .setConstant(False)
#     alphawt1.setConstant(False)
#     alphawt2.setConstant(False)
#     nwt1    .setConstant(False)
#     nwt2    .setConstant(False)
#     sigmart1.setConstant(True)
#     sigmart2.setConstant(True) 
#     massrt  .setConstant(True)    
#     f1rt    .setConstant(True)   
#     fitFunction.fitTo(data, 
#                               RooFit.Extended(True), 
#                               RooFit.Range("full"), 
#                               RooFit.Verbose(False),
#                               ROOT.RooFit.Constrain(c_vars)
#                          )     
#     sigmart1.setConstant(False)
#     sigmart2.setConstant(False) 
#     massrt  .setConstant(False)    
#     f1rt    .setConstant(False)   

    r = fitFunction.fitTo(data, 
#                           RooFit.Extended(True), 
                          RooFit.Range("full"), 
                          RooFit.Verbose(False),
                          ROOT.RooFit.Constrain(c_vars),
                          ROOT.RooFit.Minimizer("Minuit2","migrad"),
                          ROOT.RooFit.Hesse(True),
                          ROOT.RooFit.Strategy(0),
                          ROOT.RooFit.Minos(False),
#                           ROOT.RooFit.Offset(True)
                         )
    print 'fit with strategy 0 done, now strategy 2'    
    r = fitFunction.fitTo(data, 
#                           RooFit.Extended(True), 
                          RooFit.Range("full"), 
                          RooFit.Verbose(False),
                          ROOT.RooFit.Constrain(c_vars),
                          ROOT.RooFit.Minimizer("Minuit2","migrad"),
                          ROOT.RooFit.Hesse(True),
                          ROOT.RooFit.Strategy(2),
                          ROOT.RooFit.Minos(False),
#                           ROOT.RooFit.Offset(True)
                         )
    print 'fit with Hesse strategy 2 done, now Minos'    
    r = fitFunction.fitTo(data, 
#                           RooFit.Extended(True), 
                          RooFit.Save(), 
                          RooFit.Range("full"), 
#                           RooFit.Verbose(False),
                          ROOT.RooFit.Constrain(c_vars),
                          ROOT.RooFit.Minimizer("Minuit2","migrad"),
#                           ROOT.RooFit.Hesse(True),
                          ROOT.RooFit.Strategy(2),
                          ROOT.RooFit.Minos(True),
#                           ROOT.RooFit.Offset(True)
                         )
                        
    r.Print()
    r.covarianceMatrix().Print()
    r.correlationMatrix().Print()
    frame = tagged_mass.frame( RooFit.Range("full") )
    data.plotOn(frame, RooFit.Binning(nbins), RooFit.MarkerSize(.7))
    fitFunction.plotOn(frame, RooFit.NormRange("full"), RooFit.Range("full"))
    drawPdfComponents(fitFunction, frame, ROOT.kAzure, RooFit.NormRange("full"), RooFit.Range("full"), isData = True)

    parList = RooArgSet (nsig, massrt, sigmart1, sigmart2, f1rt, meanwt, sigmawt, alphawt1)
    parList.add(alphawt2)
    parList.add(nwt1)
    parList.add(nwt2)
    parList.add(frt)
    fitFunction.paramOn(frame, RooFit.Parameters(parList), RooFit.Layout(0.62,0.86,0.89))
    frame.Draw()
    niceFrame(frame, '')
    frame. addObject(_writeFitStatus(r))

    if not args.year=='test':  writeCMS(frame, args.year, [ q2binning[ibin], q2binning[ibin+1] ], 1)
    frame.Draw()
    c1.SaveAs('fit_results_mass_checkOnMC/save_fit_data_%s_%s_LMNR_pol2.pdf'%(ibin, args.year))
    c1.SetLogy()
    c1.SaveAs('fit_results_mass_checkOnMC/save_fit_data_%s_%s_LMNR_pol2_logScale.pdf'%(ibin, args.year))
    c1.SetLogy(False)

#     c1.SaveAs('fit_results_mass_checkOnMC/save_fit_data_%s_%s_LMNR_newBDT_inclusiveQ.pdf'%(ibin, args.year))








tData = ROOT.TChain('ntuple')
tMC = ROOT.TChain('ntuple')

if args.year == 'test':
    tData.Add('/gwteray/users/fiorendi/final_ntuples_p5prime_allyears/2016Data_100k.root')
    tMC.Add('/gwteray/users/fiorendi/final_ntuples_p5prime_allyears/2016MC_LMNR_100k.root')
else:    
    tData.Add('/gwpool/users/fiorendi/p5prime/miniAOD/CMSSW_10_2_14/src/miniB0KstarMuMu/miniKstarMuMu/bdt/final_ntuples/2018MC_LMNR_100k.root')
#     tMC.Add('/gwteray/users/fiorendi/final_ntuples_p5prime_allyears/%sMC_JPSI.root'%args.year)
#     tMC.Add('/gwteray/users/fiorendi/final_ntuples_p5prime_allyears/%sMC_LMNR.root'%args.year)
    tMC.Add('/gwpool/users/fiorendi/p5prime/miniAOD/CMSSW_10_2_14/src/miniB0KstarMuMu/miniKstarMuMu/bdt/final_ntuples/2018MC_LMNR.root')
#     tMC.Add('/gwpool/users/fiorendi/p5prime/CMSSW_8_0_24/src/B0KstarMM/B0KstMuMu/bdt/feb5_ntuples_fixPUW/final_ntuples/%sMC_LMNR_NoL1Selection.root'%args.year)


tagged_mass     = RooRealVar("tagged_mass" , "#mu^{+}#mu^{-}K#pi mass", 4.9, 5.6, "GeV")
mumuMass        = RooRealVar("mumuMass"    , "mumuMass" , 0, 6);
mumuMassE       = RooRealVar("mumuMassE"   , "mumuMassE", 0, 10000);
tagB0           = RooRealVar("tagB0"       , "tagB0"    , 0, 2);

tagged_mass.setRange("full",   5.0,5.6) ;
tagged_mass.setRange("mcrange",5.0,5.6) ;
thevars = RooArgSet()
thevars.add(tagged_mass)
thevars.add(mumuMass)
thevars.add(mumuMassE)
thevars.add(tagB0)

fulldata   = RooDataSet('fulldata', 'fulldataset', tData,  RooArgSet(thevars))

# import pdb; pdb.set_trace()
## add to the data dataset a random variable, in order to scale it to desired stat
nDataEntries = fulldata.sumEntries()
randVar = RooRealVar("rand","rand",0,1) 
p0      = RooPolynomial("px","px",randVar) ;
rDataset = p0.generate(RooArgSet(randVar),int(nDataEntries))
fulldata.merge(rDataset) 

## add to the input tree the combination of the variables, to be used for the cuts on the dimuon mass
deltaB0Mfunc = RooFormulaVar("deltaB0M", "deltaB0M", "@0 - @1", RooArgList(tagged_mass,B0Mass) )
deltaJMfunc  = RooFormulaVar("deltaJpsiM" , "deltaJpsiM" , "@0 - @1", RooArgList(mumuMass,JPsiMass) )
deltaPMfunc  = RooFormulaVar("deltaPsiPM" , "deltaPsiPM" , "@0 - @1", RooArgList(mumuMass,PsiPMass) )
deltaB0M     = fulldata.addColumn(deltaB0Mfunc) ;
deltaJpsiM   = fulldata.addColumn(deltaJMfunc) ;
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

out_f = TFile ("fit_results_mass_checkOnMC/results_fits_%s_newCB_newBDT.root"%args.year,"RECREATE") 
# out_f = TFile ("fit_results_mass_checkOnMC/results_fits_%s_newCB_newBDT_inclusiveQ.root"%args.year,"RECREATE") 

w = ROOT.RooWorkspace("w")
initial_n_1 =  3.
initial_n_2 =  1.
initial_a_1 =  1.
initial_a_2 = -1.
initial_sigma1 = 0.028
initial_sigma2 = 0.048
initial_sigmaCB = 0.048
nRT_forLater = 10
nWT_forLater = 10


for ibin in range(len(q2binning)-1):

    print 'dimuon selection: ', args.dimusel
    if args.dimusel == 'rejectPsi' and \
       (q2binning[ibin] == 8.68 or q2binning[ibin] == 12.86): 
           continue
           
    nRT_forLater = fitMC(rt_mc, True, ibin)
    nWT_forLater = fitMC(wt_mc, False, ibin)

    fitData(fulldata, ibin, nRT_forLater, nWT_forLater)
    print ' --------------------------------------------------------------------------------------------------- '


out_f.Close()
w.writeToFile(out_f.GetName(), False)
