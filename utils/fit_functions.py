import ROOT
from ROOT import gSystem

gSystem.Load('libRooFit')
from ROOT import RooFit, RooRealVar, RooDataSet, RooArgList, RooTreeData, RooArgSet, RooAddPdf, RooFormulaVar
from ROOT import RooGaussian, RooExponential, RooChebychev, RooProdPdf, RooCBShape, TFile, RooPolynomial, RooVoigtian, RooBreitWigner, RooFFTConvPdf
#from uncertainties import ufloat, umath

from .utils import *


def _import(wsp, obj):
    getattr(wsp, 'import')(obj)


def singleG(mean, sigma_, tagged_mass, w, fn, bin):

#     mean         = RooRealVar ("mean^{%s}"%fn   , "massSG"        ,  mean_     ,      5,    6, "GeV")
    sigma        = RooRealVar ("#sigma_{%s}^{%s}"%(fn,bin) , "sigmaSG"       ,  sigma_    ,     0,    1, "GeV")
    singlegaus   = RooGaussian("gaus_%s_%s"%(fn,bin)        , "singlegaus"     , tagged_mass,  mean, sigma)
    _import(w,singlegaus)


def doubleG(mean_, sigma1_, sigma2_, f1_, tagged_mass, w, fn):

    mean         = RooRealVar ("mean^{%s}"%fn          , "massDG"         ,  mean_      ,      5,    6, "GeV")
    sigma1       = RooRealVar ("#sigma_{1}^{%s}"%fn    , "sigmaDG1"       ,  sigma1_    ,      0,    1, "GeV")
    signalGauss1 = RooGaussian("dg_firstGauss_%s"%fn   , "firstGauss"     ,  tagged_mass,   mean, sigma1)

    sigma2       = RooRealVar ("#sigma_{2}^{%s}"%fn    , "sigmaDG2"       ,  sigma2_    ,      0,   0.12, "GeV")
    signalGauss2 = RooGaussian("dg_secondGauss_%s"%fn  , "secondGauss"    ,  tagged_mass,   mean, sigma2)

    f1           = RooRealVar ("f^{%s}"%fn             , "f1"             ,  f1_        ,     0.,    1. )
    doublegaus   = RooAddPdf  ("doublegaus_%s"%fn      , "gaus1+gaus2"    ,  RooArgList(signalGauss1,signalGauss2), RooArgList(f1))
    _import(w,doublegaus)


def tripleG(doublegaus, mean, sigma3_, f2_, tagged_mass, w):

    sigma3       = RooRealVar ("#sigma_{TG3}"  , "sigmaTG3"        ,  sigma3_    ,      0,   0.2, "GeV")
    signalGauss3 = RooGaussian("thirdGauss"    , "thirdGauss"      ,  tagged_mass,   mean, sigma3)
    f2           = RooRealVar ("f2"            , "f2"              ,  f2_        ,     0.,    1. )
    triplegaus   = RooAddPdf  ("triplegaus"    , "doublegaus+gaus3",  RooArgList(doublegaus,signalGauss3), RooArgList(f2))
    _import(w,triplegaus)


def crystalBall(mean, sigma_, alpha_, n_, tagged_mass, w, fn, bin, rangeAlpha):

    sigmaCB      = RooRealVar ("#sigma_{%s}^{%s}"%(fn, bin)   , "sigmaCB_%s"%fn        ,  sigma_  ,     0,   1  )
    alpha        = RooRealVar ("#alpha_{%s}^{%s}"%(fn, bin)   , "#alpha_{%s}^{%s}"%(fn, bin) ,  alpha_  ,    rangeAlpha[0],  rangeAlpha[1] ) # was 0 - 5
    n            = RooRealVar ("n_{%s}^{%s}"%(fn, bin)        , "n_%s"%fn              ,  n_      ,      0.001,   200	 )
    cbshape      = RooCBShape ("cbshape_%s_%s"%(fn,bin)       , "cbshape_%s_%s"%(fn, bin)        ,  tagged_mass, mean, sigmaCB, alpha, n)
    _import(w,cbshape)


def doubleCB(cbshape1, cbshape2, f3_, tagged_mass, w, fn):

    f3           = RooRealVar ("f^{%s}"%fn      , "f3"       ,  f3_  ,     0.,   1.)
    doublecb     = RooAddPdf  ("doublecb_%s"%fn, "doublecb"  ,  RooArgList(cbshape1,cbshape2), RooArgList(f3))
    _import(w,doublecb)


def gausCB(cbshape, gaus, f3_, tagged_mass, w, fn, bin):

    f3           = RooRealVar ("f^{%s%s}"%(fn, bin)  , "f3"      ,  f3_  ,     0.,   1.)
    gauscb       = RooAddPdf  ("gauscb_%s_%s"%(fn,bin)  , "gauscb"  ,  RooArgList(gaus,cbshape), RooArgList(f3))
    _import(w,gauscb)


def doubleGausCB(cbshape, doublegaus, f3_, tagged_mass, w):
    f4           = RooRealVar ("f4"            , "f4"            ,  f4_  ,     0.,   1.)
    doublegauscb = RooAddPdf  ("doublegauscb"  , "doublegauscb"  ,  RooArgList(doublegaus,cbshape), RooArgList(f4))
    _import(w,doublegauscb)


def voigtian(mean, width_, sigma_, tagged_mass, w):

    sigmaV      = RooRealVar ("#sigma_{V}"   , "sigmaV"        ,  sigma_  ,     0,   10 )
    widthV      = RooRealVar ("widthV"       , "widthV"        ,  width_  ,     0,    5 )
    vgshape     = RooVoigtian ("vgshape"     , "vgshape"       ,  tagged_mass, mean, widthV, sigmaV)
    _import(w,vgshape)


def bwcb(mean_, width_, sigma_, alpha_, n_, fn, tagged_mass, w):

## Breit-Wigner
    meanBW       = RooRealVar ("massBW_%s"%fn        , "massBW_%s"%fn   ,  mean_   ,      3,    7, "GeV")
    widthBW      = RooRealVar ("widthBW_%s"%fn       , "widthBW_%s"%fn  ,  width_  ,     0,   10 )
    bwshape      = RooBreitWigner ("bwshape_%s"%fn   , "bwshape_%s"%fn  ,  tagged_mass, meanBW, widthBW)

    meanCB       = RooRealVar ("massBW_%s"%fn          , "massBW_%s"%fn       ,  0.)
    sigmabwCB    = RooRealVar ("#sigma_{bwCB}_%s"%fn   , "sigmabwCB_%s"%fn    ,  sigma_  ,     0,   1  )
    alphabw      = RooRealVar ("#alphabw_%s"%fn        , "alphabw_%s"%fn      ,  alpha_  ,     0,    10 ) # was 0 - 5
    nbw          = RooRealVar ("nbw_%s"%fn             , "nbw_%s"%fn          ,  n_      ,     0,   25 )
    cbshape      = RooCBShape ("cbshapebw_%s"%fn       , "cbshapebw_%s"%fn        ,  tagged_mass, meanCB, sigmabwCB, alphabw, nbw)
    
    cbbw = RooFFTConvPdf( "cbbw_%s"%fn, "cbbw_%s"%fn, tagged_mass, bwshape, cbshape);
    _import(w,cbbw)



# def calculateTotSigma(sigma_vec, f_vec):
#     
#     if len(sigma_vec) != len(f_vec) + 1:
#         print 'calculateTotSigma: Warning! wrong vector lenghts'
#         return 0
#     totSigma = ufloat(0., 0.)
#     sum_f    = ufloat(0., 0.)
# 
#     s1 = sigma_vec[0]
#     s2 = sigma_vec[1]
#     f1 = f_vec[0]
#     totSigma = pow(s1, 2) * f1 + (1-f1)* pow(s2,2)
#     for i in range(len(f_vec)-1):
#         s1 = totSigma
#         s2 = sigma_vec[i+2]
#         f1 = f_vec[i+1]
#         totSigma = pow(s1, 2) * f1 + (1-f1)* pow(s2,2)
#     
# #     print 'calculateTotSigma: totSigma = ', totSigma, ' sum_f = ',  sum_f
# #     totSigma += (1-sum_f) *  pow(sigma_vec[-1],2)  
#     print 'calculateTotSigma: ', umath.sqrt(totSigma)
#     return umath.sqrt(totSigma) 
    



from itertools import product
def drawPdfComponents(fitFunction, frame, base_color, normrange, range, isData = False):

    from ROOT import gStyle
#     gStyle.SetPalette(109);
    colorlist = []
    gStyle.SetPalette(55);
    colorlist.append(ROOT.TColor.GetColorPalette(20))
#     colorlist.append(ROOT.TColor.GetColorPalette(50))
    colorlist.append(ROOT.TColor.GetColorPalette(80))
#     colorlist.append(ROOT.TColor.GetColorPalette(100))
    colorlist.append(ROOT.TColor.GetColorPalette(130))
    colorlist.append(ROOT.TColor.GetColorPalette(150))
    colorlist.append(ROOT.TColor.GetColorPalette(170))
    colorlist.append(ROOT.TColor.GetColorPalette(190))
    colorlist.append(ROOT.TColor.GetColorPalette(210))
    colorlist.append(ROOT.TColor.GetColorPalette(250))
    colorlist.append(ROOT.TColor.GetColorPalette(2))
    pdf_components = fitFunction.getComponents()
    iter = pdf_components.createIterator()
    var = iter.Next();  color = 0
    list_to_plot      = ['fitfunction', 'c_theRTgauss', 'c_theWTgauss', 'bkg_exp', 'bkg_pol', 'cbshape_bs',
                         'c_signalFunction', 
                         #'doublecb_', 'doublecb_RT', 'gauscb_RT_', 
                         #'cbshape_RT1_', 'cbshape_RT2_', 'gaus_RT2_', 
                         'bs_shape_kk', 'bs_shape_phi', 'bs_shape_kst', 'erf_pdf' ]
    list_to_plot_bins = ['%s%s'%(i,ibin) for i,ibin in product(list_to_plot,list(xrange(8)))]
#     second_list_to_plot = ['cbshape_bs4_kk', 'n_cb_bs_kk']
#     import pdb; pdb.set_trace()
    while var :
        ### https://root-forum.cern.ch/t/roofit-normalization/23644/5
        if isData and var.GetName() not in list_to_plot and var.GetName() not in list_to_plot_bins:  
            var = iter.Next()
            continue
        print var.GetName(), ' -----> plotting  '
        fitFunction.plotOn(frame, RooFit.Components(var.GetName()), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(colorlist[color]), normrange, range)
#         fitFunction.plotOn(frame, RooFit.Components(var.GetName()), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(base_color+color), normrange, range)
        var = iter.Next()
        color += 1




### maybe this can be used somehow
#     params = doublegaus.getParameters(tagged_mass) ;
#     w.defineSet("doublegaus_parameters",*params) ;
### save parameters (maybe can be used)
#     w.saveSnapshot("reference_fit",*params,kTRUE) ;
