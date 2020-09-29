import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("bin", help = "choose q2 bin range", default = -1, type = float)
parser.add_argument("ns_min", help = "choose nsigma close to signal", default = 3, type = float)
parser.add_argument("ns_max", help = "choose nsigma far from signal", default = 3, type = float)
# # parser.add_argument("year", help = "choose year [format:2016, 20162017]", default = '2016')
args = parser.parse_args()

import ROOT
from ROOT import RooFit, gSystem
from math import sqrt, sin, cos, pow
import math, itertools
import pdb
from pdb import set_trace
from uncertainties import ufloat
from uncertainties.umath import sqrt as usqrt
from uncertainties.umath import sin as usin
from uncertainties.umath import cos as ucos
from uncertainties.umath import pow as upow

from ROOT import RooFit, RooRealVar, RooDataSet, RooArgList, RooTreeData, RooArgSet, RooAddPdf, RooFormulaVar

from collections import OrderedDict
ROOT.gROOT.SetBatch(True)

ROOT.RooMsgService.instance().setGlobalKillBelow(4)

gSystem.Load('libRooFit')
gSystem.Load('utils/func_roofit/libRooDoubleCBFast')

ns_min = args.ns_min
ns_max = args.ns_max

class fit_pars(object):
    '''
    '''
    def __init__(self):
        self.Reset()

    def Reset(self):
        self.mean            = ufloat(-10,0)
        self.sigmaRT1        = ufloat(-10,0)
        self.sigmaRT2        = ufloat(-10,0)
        self.fRT1            = ufloat(-10,0)
        self.sigmaWT         = ufloat(-10,0)
        self.fRT             = ufloat(-10,0)
        self.RTSigma         = ufloat(-10,0)
        self.totSigma        = ufloat(-10,0)

    def __str__(self):
        toWrite = ' \
mean     = %f +/- %f \n \
sigmaRT1 = %f +/- %f \n \
sigmaRT2 = %f +/- %f \n \
fRT1     = %f +/- %f \n \
sigmaWT  = %f +/- %f \n \
fRT      = %f +/- %f \n \
\n \
RTSigma  = %f +/- %f \n \
totSigma = %f +/- %f \n '%(self.mean.n     , self.mean.s, 
                           self.sigmaRT1.n , self.sigmaRT1.s,
                           self.sigmaRT2.n , self.sigmaRT2.s,
                           self.fRT1.n     , self.fRT1.s,
                           self.sigmaWT.n  , self.sigmaWT.s,
                           self.fRT.n      , self.fRT.s,
                           self.RTSigma.n  , self.RTSigma.s,
                           self.totSigma.n , self.totSigma.s
                           )
        return toWrite 

    def setSigma(self):
        if self.sigmaRT2.n > 0:
            self.RTSigma = usqrt(self.fRT1    * upow(self.sigmaRT1,2) + \
                                (1-self.fRT1) * upow(self.sigmaRT2,2) )
            self.totSigma = usqrt(self.fRT    * upow(self.RTSigma,2) + \
                                 (1-self.fRT) * upow(self.sigmaWT,2) )
        else:
            self.RTSigma  = self.sigmaRT1
            self.totSigma = usqrt(self.fRT    * upow(self.sigmaRT1,2) + \
                                  (1-self.fRT)* upow(self.sigmaWT,2) )


bin_range = [0,1,2,3,5,7]
if args.bin != -1:
    bin_range = [int(args.bin)]  

range_sh1 = [1]
range_sh2 = [1]
range_c   = [1]

fname_oldr = 'results_data_fits_2018_expbkg.root'

try:
  fo = ROOT.TFile(fname_oldr,'open')
except:
  print ('file %s or %s not found'%(fo))
in_ws = fo.Get('data_w') 


for ii, ibin in enumerate(bin_range):
    
    par_list = {
      "mean^{RT%s}"%ibin        : 'mean'    ,
      "#sigma_{CB}^{WT%s}"%ibin : 'sigmaWT' ,
      "F_{RT}%s"%ibin           : 'fRT'
    }
    if ibin < 5:
      par_list["#sigma_{CB}^{RT%s}"%ibin] = 'sigmaRT1'
    else:
      par_list['#sigma_{CBRT0}^{%s}'%ibin] = 'sigmaRT1'
      par_list['#sigma_{CBRT1}^{%s}'%ibin] = 'sigmaRT2'
      par_list['f^{RT%s}'%ibin] = 'fRT1'
        
    par_list = OrderedDict(sorted(par_list.items(), key=lambda t: t[0]))
    
    ## retrieve parameter values and calculate sigma
    keys = fo.GetListOfKeys()
    pars_o = fit_pars()
    for i in keys:
        if 'results_data_%s'%ibin in i.GetName():   
            fitres = fo.Get(i.GetName()) 
            for ipar,ii in par_list.iteritems():
                setattr(pars_o, ii, ufloat (fitres.floatParsFinal().find( ipar).getVal(), 
                                            fitres.floatParsFinal().find( ipar).getError()))
            pars_o.setSigma()
            print 'bin: ', ibin, '\n',  pars_o

#     pdb.set_trace()        
    
    ## retrieve signal and total pdf and calculate integral in the sidebands

    in_ws.loadSnapshot("reference_fit_data_%s"%ibin)
#     fit_pdf = in_ws.pdf('fitfunction%s'%ibin)
    sig_pdf = in_ws.pdf('c_signalFunction%s'%ibin)
    bkg_pdf = in_ws.pdf('bkg_exp_%s'%ibin)
    x = in_ws.var('tagged_mass')
    x.setRange('lSB', (pars_o.mean - ns_max*pars_o.totSigma).n, (pars_o.mean - ns_min*pars_o.totSigma).n )
    x.setRange('rSB', (pars_o.mean + ns_min*pars_o.totSigma).n, (pars_o.mean + ns_max*pars_o.totSigma).n )
    x.setRange('sig', (pars_o.mean - 3*pars_o.totSigma).n,      (pars_o.mean + 3*pars_o.totSigma).n )
    x.setRange('all', 5.0, 5.6)
    xset = RooArgSet(x)
    normSet = ROOT.RooFit.NormSet(xset)
    lInt = sig_pdf.createIntegral(RooArgSet(x), normSet, RooFit.Range('lSB')) 
    rInt = sig_pdf.createIntegral(RooArgSet(x), normSet, RooFit.Range('rSB')) 
    aInt = sig_pdf.createIntegral(RooArgSet(x), normSet, RooFit.Range('all')) 
    sInt = sig_pdf.createIntegral(RooArgSet(x), normSet, RooFit.Range('sig')) 
    print 'remaining signal events in the left SB: %.3f' %lInt.getVal(), \
          '\nremaining of signal events in the right SB: %.3f' %rInt.getVal() 

    n_sig_fit = in_ws.var('Yield%s'%ibin).getVal()
    n_bkg_fit = in_ws.var('nbkg%s'%ibin).getVal()

    lIntb = bkg_pdf.createIntegral(RooArgSet(x), normSet, RooFit.Range('lSB')) 
    rIntb = bkg_pdf.createIntegral(RooArgSet(x), normSet, RooFit.Range('rSB')) 
    print 'remaining bkg events in the left SB: %.3f' %lIntb.getVal(), \
          '\nremaining bkg events in the right SB: %.3f' %rIntb.getVal() 

    
    fit_pdf = in_ws.pdf('fitfunction%s'%ibin)
    lallInt = fit_pdf.createIntegral(RooArgSet(x), normSet, RooFit.Range('lSB')) 
    rallInt = fit_pdf.createIntegral(RooArgSet(x), normSet, RooFit.Range('rSB')) 
    print 'percentage of signal events in the left SB: %.3f'    %(lInt.getVal()*n_sig_fit/(lallInt.getVal()*(n_bkg_fit+n_sig_fit))), \
          '\npercentage of signal events in the right SB: %.3f' %(rInt.getVal()*n_sig_fit/(rallInt.getVal()*(n_bkg_fit+n_sig_fit)))

    print 'L:', (pars_o.mean - ns_max*pars_o.totSigma).n, (pars_o.mean - ns_min*pars_o.totSigma).n
    print 'R:', (pars_o.mean + ns_min*pars_o.totSigma).n, (pars_o.mean + ns_max*pars_o.totSigma).n

#     print 'signal region:', sInt.getVal() 
    print ' --------'

# (fit_pdf.createIntegral(xset, RooFit.Range('sig')).getVal() / fit_pdf.createIntegral(xset).getVal())


#     pdb.set_trace()        

fo.Close()















