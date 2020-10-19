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
from ROOT import RooFit, RooRealVar, RooDataSet, RooArgList, RooTreeData, RooArgSet, RooAddPdf, RooFormulaVar
from ROOT import RooGaussian, RooExponential, RooChebychev, RooProdPdf, RooCBShape, TFile, RooPolynomial
import sys, math, pdb
from uncertainties import ufloat

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
    txt = ROOT.TLatex(.16,.8, "fit #chi^{2}: %.1f "%chi2 )
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


sys.path.append("../utils")
from utils.utils import *
from utils.fit_functions import *

nbins = 60
nSigma_psiRej = 3.
cut_base      = applyB0PsiCut(args.dimusel, nSigma_psiRej)

## uncertainty on fRT
frt_sigma = frt_sigmas[args.year]

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




   
   
   
def fitData(fulldata, ibin, w):

    cut  = cut_base + '&& (mumuMass*mumuMass > %s && mumuMass*mumuMass < %s)'%(q2binning[ibin], q2binning[ibin+1])
    data = fulldata.reduce(RooArgSet(tagged_mass,mumuMass,mumuMassE), cut)

    nrt_mc = _getFittedVar("nRT_%s"%ibin, w)
    nwt_mc = _getFittedVar("nWT_%s"%ibin, w)
    fraction = nrt_mc / (nrt_mc + nwt_mc)
    print 'mistag fraction on MC for bin ', ibin , ' : ' , fraction.n , '+/-', fraction.s 
    
    ### creating RT component
    w.loadSnapshot("reference_fit_RT_%s"%ibin)
    mean_rt      = w.var("mean^{RT%s}"%ibin)
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

    ## double cb old
    else:
        sigma_rt2    = w.var("#sigma_{RT2}^{%s}"%ibin)
        alpha_rt1    = w.var("#alpha_{RT1}^{%s}"%ibin)
        alpha_rt2    = w.var("#alpha_{RT2}^{%s}"%ibin)
        n_rt1        = w.var("n_{RT1}^{%s}"%ibin)
        n_rt2        = w.var("n_{RT2}^{%s}"%ibin)
        f1rt        = w.var("f^{RT%s}"%ibin)

    theRTgauss  = w.pdf("doublecb_RT%s"%ibin)   

    ### creating WT component
    w.loadSnapshot("reference_fit_WT_%s"%ibin)
    mean_wt      = w.var("mean_{WT}^{%s}"%ibin)
    sigma_wt     = w.var("#sigma_{WT1}^{%s}"%ibin)
    alpha_wt1    = w.var("#alpha_{WT1}^{%s}"%ibin)
    alpha_wt2    = w.var("#alpha_{WT2}^{%s}"%ibin)
    n_wt1        = w.var("n_{WT1}^{%s}"%ibin)
    n_wt2        = w.var("n_{WT2}^{%s}"%ibin)
    theWTgauss  = w.pdf("doublecb_%s"%ibin)   

    ### creating variable for the difference between the two peaks
    deltaPeaks = RooFormulaVar("deltaPeaks%s"%ibin, "@0 - @1", RooArgList(mean_rt, mean_wt))  
    frt              = RooRealVar ("F_{RT}%s"%ibin   , "frt"             , fraction.n , 0, 1)
    signalFunction   = RooAddPdf  ("sumgaus%s"%ibin  , "rt+wt"           , RooArgList(theRTgauss,theWTgauss), RooArgList(frt))

    ### now create background parametrization
    slope         = RooRealVar    ("slope_%s"%ibin   , "slope"           ,    0.5,   -10, 10);
    bkg_exp       = RooExponential("bkg_exp%s"%ibin , "exponential"     ,  slope,   tagged_mass  );
    pol_c1        = RooRealVar    ("p1_%s"%ibin      , "coeff x^0 term"  ,    0.5,   -10, 10);
    pol_c2        = RooRealVar    ("p2_%s"%ibin      , "coeff x^1 term"  ,    0.5,   -10, 10);
    bkg_pol       = RooPolynomial ("bkg_pol%s"%ibin  , "2nd order pol"   ,  tagged_mass, RooArgList(pol_c1, pol_c2));

    fsig          = RooRealVar("fsig%s"%ibin   , "fsig"    ,     0.9,     0,   1);
#     nsig          = RooRealVar("Yield%s"%ibin  , "signal frac"    ,     1000,     0,   10000);
#     nbkg          = RooRealVar("nbkg%s"%ibin   , "bkg fraction"   ,     1000,     0,   500000);
    nsig        = RooRealVar("Yield"         , "signal frac"   ,   600000,     0,   5000000);
    nbkg        = RooRealVar("nbkg"          , "bkg fraction"  ,   100000,     0,   2000000);

#     if ibin==4:  
#         nsig.setRange(500000,1500000)
#         nsig.setVal(900000)
#         nbkg.setRange(80000,1000000)
#         nbkg.setVal(100000)


    ### creating constraints
    c_vars = RooArgSet()
    c_pdfs = RooArgSet()

    c_sigma_rt1   = _constrainVar(sigma_rt1, 1)
    c_alpha_rt1   = _constrainVar(alpha_rt1, 1)
    c_alpha_rt2   = _constrainVar(alpha_rt2, 1)
    c_n_rt1       = _constrainVar(n_rt1, 1)
    c_n_rt2       = _constrainVar(n_rt2, 1)

    c_sigma_wt    = _constrainVar(sigma_wt,  1)
    c_alpha_wt1   = _constrainVar(alpha_wt1, 1)
    c_alpha_wt2   = _constrainVar(alpha_wt2, 1)
    c_n_wt1       = _constrainVar(n_wt1, 1)
    c_n_wt2       = _constrainVar(n_wt2, 1)


    if ibin < 4:
        c_pdfs = RooArgSet(c_sigma_rt1, c_alpha_rt1, c_alpha_rt2, c_n_rt1, c_n_rt2)
        c_vars = RooArgSet(sigma_rt1,     alpha_rt1,   alpha_rt2,   n_rt1,   n_rt2)
    else:
        c_sigma_rt2   = _constrainVar(sigma_rt2, 1)
        c_pdfs = RooArgSet(c_sigma_rt1, c_sigma_rt2, c_alpha_rt1, c_alpha_rt2, c_n_rt1, c_n_rt2)
        c_vars = RooArgSet(  sigma_rt1,   sigma_rt2,   alpha_rt1,   alpha_rt2,   n_rt1,   n_rt2)
    
    c_pdfs.add(c_sigma_wt);  c_vars.add(sigma_wt)
    c_pdfs.add(c_alpha_wt1); c_vars.add(alpha_wt1)
    c_pdfs.add(c_alpha_wt2); c_vars.add(alpha_wt2)
    c_pdfs.add(c_n_wt1);     c_vars.add(n_wt1)
    c_pdfs.add(c_n_wt2);     c_vars.add(n_wt2)

    c_deltaPeaks = RooGaussian("c_deltaPeaks%s"%ibin , "c_deltaPeaks", deltaPeaks, ROOT.RooFit.RooConst( deltaPeaks.getVal() ), 
                                ROOT.RooFit.RooConst( 0.0005 )  ## value to be checked
                                ) 
    c_pdfs.add(c_deltaPeaks)
    c_vars.add(deltaPeaks)
    c_frt            = RooGaussian("c_frt%s"%ibin    , "c_frt"           , frt,  ROOT.RooFit.RooConst(fraction.n) , ROOT.RooFit.RooConst(frt_sigma[ibin]) )
    c_pdfs.add(c_frt)
    c_vars.add(frt)
    
    constr_list = RooArgList(c_pdfs)
    constr_list.add(signalFunction)
    c_signalFunction = RooProdPdf ("c_signalFunction", "c_signalFunction", constr_list)     

#     mean        = RooRealVar ("mass"         , "mean"          ,  B0Mass_,   3,    7, "GeV")
#     sigma       = RooRealVar ("#sigma_{1}"   , "sigma"         ,  0.028,     0,   10, "GeV")
#     signalGauss = RooGaussian("signalGauss"  , "signal gauss"  ,  tagged_mass,  mean,sigma)
#     
#     sigma2       = RooRealVar ("#sigma_{2}"    , "sigma2"         ,  0.048,     0,   0.07, "GeV")
#     signalGauss2 = RooGaussian("signalGauss2"  , "signal gauss2"  ,  tagged_mass,  mean,sigma2)
#     f1           = RooRealVar ("f1"            , "f1"             ,  0.8  ,     0.,   1.)
#     gaus         = RooAddPdf  ("gaus"          , "gaus1+gaus2"    , RooArgList(signalGauss,signalGauss2), RooArgList(f1))
    
#     pol_c1      = RooRealVar ("p1"           , "coeff x^0 term",    0.5,   -10, 10);
#     pol_c2      = RooRealVar ("p2"           , "coeff x^1 term",    0.5,   -10, 10);
#     pol_c3      = RooRealVar ("p3"           , "coeff x^2 term",    0.5,   -10, 10);
#     slope       = RooRealVar ("slope"        , "slope"         ,    0.5,   -10, 10);
#     bkg_exp     = RooExponential("bkg_exp"   , "exponential"   ,  slope,   tagged_mass  );
#     bkg_pol     = RooChebychev("bkg_pol"     , "2nd order pol" ,  tagged_mass, RooArgList(pol_c1,pol_c2));

    fitFunction = RooAddPdf ("fitfunction%s"%ibin , "fit function"  ,  RooArgList(c_signalFunction, bkg_exp), RooArgList(nsig, nbkg))
    
#     r = fitFunction.fitTo(data, 
# #                           RooFit.Extended(True), 
#                           RooFit.Range("full"), 
#                           ROOT.RooFit.Constrain(c_vars),
#                           ROOT.RooFit.Minimizer("Minuit2","migrad"),
#                           ROOT.RooFit.Hesse(True),
#                           ROOT.RooFit.Strategy(2),
#                           ROOT.RooFit.Minos(False),
#                          )
    print 'fit with Hesse strategy 2 done, now Minos'    
    r = fitFunction.fitTo(data, 
                          RooFit.Extended(True), 
                          RooFit.Save(), 
                          RooFit.Range("full"), 
                          RooFit.Verbose(False),
                          ROOT.RooFit.Constrain(c_vars),
#                           ROOT.RooFit.Minimizer("Minuit2","migrad"),
#                           ROOT.RooFit.Hesse(True),
#                           ROOT.RooFit.Strategy(2),
#                           ROOT.RooFit.Minos(False),
                         )

    r.Print()
    r.correlationMatrix().Print()
    fitStats['data%s'%(ibin)] = r.status()
    covStats['data%s'%(ibin)] = r.covQual()
    frame = tagged_mass.frame( RooFit.Range("full") )
    data.plotOn(frame, RooFit.Binning(nbins), RooFit.MarkerSize(.7))
    fitFunction.plotOn(frame, RooFit.NormRange("full"), RooFit.Range("full"))

    ## evaluate sort of chi2 and save number of RT/WT events
    observables = RooArgSet(tagged_mass)
    flparams = fitFunction.getParameters(observables)
    nparam = int(flparams.selectByAttrib("Constant",ROOT.kFALSE).getSize())
    pdfstring = "fitfunction%s_Norm[tagged_mass]_Range[full]_NormRange[full]"%ibin
    chi2s['data%s'%ibin] = frame.chiSquare(pdfstring, "h_fulldata",  nparam)
    frame. addObject(_writeChi2( chi2s['data%s'%ibin] ))

    drawPdfComponents(fitFunction, frame, ROOT.kAzure, RooFit.NormRange("full"), RooFit.Range("full"), isData = True)
#     fitFunction.paramOn(frame, RooFit.Layout(0.62,0.86,0.89))

    parList = RooArgSet (nsig, mean_rt, sigma_rt, alpha_rt1, alpha_rt2, n_rt1, n_rt2, mean_wt, sigma_wt)
#     parList.add(alphawt1)
#     parList.add(alphawt2)
#     parList.add(nwt1)
#     parList.add(nwt2)
    parList.add(frt)
    fitFunction.paramOn(frame, RooFit.Parameters(parList), RooFit.Layout(0.62,0.86,0.89))

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
    if not args.year=='test':  writeCMS(frame, args.year, [ q2binning[ibin], q2binning[ibin+1] ], 0)
    frame.Draw()

    ## add plot of pulls
    lowerPad.cd()
    hpull = frame.pullHist("h_fulldata", pdfstring)
    frame2 =  tagged_mass.frame(RooFit.Range("full"), RooFit.Title(''))
    frame2.addPlotable(hpull,"P") 
    niceFrameLowerPad(frame2, 'pull')
    frame2.Draw()
    line = ROOT.TLine(5.0,1,5.6,1)
    line.SetLineColor(ROOT.kGreen+3)
    line.Draw()

    for ilog in [True,False]:
        upperPad.SetLogy(ilog)
        c1.SaveAs('fit_results_mass/save_fit_data_%s_%s_LMNR_Update%s_newSigmaFRT_pars_Jpsi.pdf'%(ibin, args.year, '_logScale'*ilog))


    out_f.cd()
    r.Write('results_data_%s'%(ibin))

    params = fitFunction.getParameters(RooArgSet(tagged_mass)) 
    out_w.saveSnapshot("reference_fit_data_%s"%(ibin),params,ROOT.kTRUE) 
    getattr(out_w, 'import')(fitFunction)
#     getattr(out_w, 'import')(signalFunction)
#     getattr(out_w, 'import')(bkg_pol)






tData = ROOT.TChain('ntuple')

if args.year == 'test':
    tData.Add('/gwteray/users/fiorendi/final_ntuples_p5prime_allyears/2016Data_100k.root')
else:    
    tData.Add('/gwteray/users/fiorendi/final_ntuples_p5prime_allyears/%sData_All_finalSelection.root'%args.year)


tagged_mass     = RooRealVar("tagged_mass" , "#mu^{+}#mu^{-}K#pi mass", 5., 5.6, "GeV")
mumuMass        = RooRealVar("mumuMass"    , "mumuMass" , 0, 6);
mumuMassE       = RooRealVar("mumuMassE"   , "mumuMassE", 0, 10000);
tagB0           = RooRealVar("tagB0"       , "tagB0"    , 0, 2);

tagged_mass.setRange("full",   5.0,5.6) ;
thevars = RooArgSet()
thevars.add(tagged_mass)
thevars.add(mumuMass)
thevars.add(mumuMassE)
thevars.add(tagB0)

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

if args.year == 'test':
    fname_mcresults = 'fit_results_mass_checkOnMC/results_fits_2016_newSigmaFRT_Jpsi.root'
else:
    fname_mcresults = 'fit_results_mass_checkOnMC/results_fits_%s_newSigmaFRT_Jpsi.root'%args.year
fo = ROOT.TFile()
try:
  fo = ROOT.TFile(fname_mcresults,'open')
except:
  print ('file %s or not found'%(fo))
w = fo.Get('w')


out_f = TFile ("results_data_fits_%s_newSigmaFRT_Jpsi.root"%args.year,"RECREATE") 
out_w = ROOT.RooWorkspace("data_w")

for ibin in range(len(q2binning)-1):

    print 'dimuon selection: ', args.dimusel
    if args.dimusel == 'rejectPsi' and \
       (q2binning[ibin] == 8.68 or q2binning[ibin] == 12.86): 
           continue
    if q2binning[ibin] < 8:  continue       
    if q2binning[ibin] > 10:  continue       
    fitData(fulldata, ibin, w)
    print ' --------------------------------------------------------------------------------------------------- '


print '--------------------------------------------------------------------------------------------------- '
print 'bin\t\t fit status \t cov. matrix \t\t chi2'
for i,k in enumerate(fitStats.keys()):    
    if i%3==0:  print '------------------------------------------------------'
    print k , '\t\t', fitStats[k], '\t\t', covStats[k], '\t\t', chi2s[k]

out_f.Close()
out_w.writeToFile(out_f.GetName(), False)
