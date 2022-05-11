import ROOT
from ROOT import gSystem

gSystem.Load('libRooFit')
from ROOT import RooFit, RooRealVar


from eras_allYears import *
# from utils.eras_allYears import *

B0Mass_   = 5.27958
JPsiMass_ = 3.096916
PsiPMass_ = 3.686109
KStMass_  = 0.896

B0Mass     = RooRealVar("B0Mass"    , "B0Mass"  , B0Mass_  )
JPsiMass   = RooRealVar("JPsiMass"  , "JPsiMass", JPsiMass_)
PsiPMass   = RooRealVar("PsiPMass"  , "PsiPMass", PsiPMass_)
KStMass    = RooRealVar("KStMass"   , "KStMass" , KStMass_ )

n_data = {}
n_data['2016'] = [205, 454, 391,  689,  627365, 1174,  37827,  704]
n_data['2017'] = [307, 581, 495, 1013,  787149, 1524,  48372,  835]
n_data['2018'] = [500, 981, 821, 1608, 1694952, 3018, 102424, 1836]
n_data['test'] = [520, 1000, 850, 1900,  100000, 3166, 10000, 1860]

# std::map<int,std::vector<float>> nbkg_years = {
#   {2016, {162, 535, 462,  810, 0.005, 1342, 0.006, 467}},
#   {2017, {185, 496, 441,  711, 0.004, 1363, 0.005, 379}},
#   {2018, {288, 842, 734, 1270, 0.002, 2954, 0.003, 779}},
# };
# 
# std::map<int,std::vector<float>> nsig_years = {
#   {2016, {205, 454, 391,  689, 0.005, 1174, 0.006,  704}},
#   {2017, {307, 581, 495, 1013, 0.004, 1524, 0.005,  835}},
#   {2018, {500, 981, 821, 1608, 0.002, 3018, 0.003, 1836}},
# };




frt_sigmas = {}
frt_sigmas['2016'] = [0.021, 0.015, 0.016, 0.011, 0.009, 0.009, 0.008, 0.011]
frt_sigmas['2017'] = [0.018, 0.013, 0.014, 0.010, 0.008, 0.008, 0.007, 0.011]
frt_sigmas['2018'] = [0.013, 0.010, 0.011, 0.007, 0.006, 0.006, 0.006, 0.007]
frt_sigmas['test'] = [0.013, 0.010, 0.011, 0.007, 0.006, 0.006, 0.006, 0.007]

fM_sigmas = {}
fM_sigmas['2016'] = [0.023, 0.015, 0.017, 0.013, 0.0005 , 0.010, 0.0018, 0.013]
fM_sigmas['2017'] = [0.018, 0.014, 0.015, 0.010, 0.0004 , 0.008, 0.0016, 0.011]
fM_sigmas['2018'] = [0.015, 0.010, 0.011, 0.008, 0.00027, 0.006, 0.0011, 0.008]
fM_sigmas['test'] = [0.013, 0.010, 0.011, 0.007, 0.006  , 0.006, 0.006 , 0.007]
## check sigmas for Jpsi and Psi


# mc_scale = {}
# mc_scale['2016'] = 2.2## 2.5
# mc_scale['2017'] = 2.1
# mc_scale['2018'] = 1.9


# q2binning_base = [
#                 1,
#                 2, 
#                 4.3,
#                 6,
#                 8.68,
#                 10.09,
#                 12.86,
#                 14.18,
#                 16,
#                 19,
# ]


def applyB0PsiCut(dimusel, nSigma_psiRej, triangular=False):

    cut_base = 'mumuMass > 0'
    if dimusel == 'keepJpsi':
      cut_base = '(abs(mumuMass - {JPSIM}) < {CUT}*mumuMassE)'.format( JPSIM=JPsiMass_, CUT=nSigma_psiRej)
    elif dimusel == 'keepPsiP':
      cut_base = '(abs(mumuMass - {PSIM}) < {CUT}*mumuMassE)'.format( PSIM=PsiPMass_, CUT=nSigma_psiRej)
    elif dimusel == 'rejectPsi' and not triangular:
      cut_base = '( abs(mumuMass - {JPSIM}) > {CUT}*mumuMassE && abs(mumuMass - {PSIM}) > {CUT}*mumuMassE &&  \
               (( mumuMass < {JPSIM} && !( abs(deltaB0M - deltaJpsiM) < 0.18 || abs(deltaB0M - deltaPsiPM) < 0.0) ) || \
                ( mumuMass > {PSIM}  && !( abs(deltaB0M - deltaJpsiM) < 0.0  || abs(deltaB0M - deltaPsiPM) < 0.08) ) || \
                ( mumuMass > {JPSIM} && mumuMass < {PSIM} && !( abs(deltaB0M - deltaJpsiM) < 0.08 || abs(deltaB0M - deltaPsiPM) < 0.09 ))))'.format(JPSIM=JPsiMass_, PSIM=PsiPMass_,  CUT=nSigma_psiRej)  
    elif dimusel == 'rejectPsi' and triangular:
      cut_base = '( abs(mumuMass - {JPSIM}) > {CUT}*mumuMassE && abs(mumuMass - {PSIM}) > {CUT}*mumuMassE &&  \
               (( mumuMass < {JPSIM} && !( (deltaB0M - deltaJpsiM) < 0.18 || abs(deltaB0M - deltaPsiPM) < 0.0) ) || \
                ( mumuMass > {PSIM}  && !( abs(deltaB0M - deltaJpsiM) < 0.0  || (deltaB0M - deltaPsiPM) > -0.08) ) || \
                ( mumuMass > {JPSIM} && mumuMass < {PSIM} && !( (deltaB0M - deltaJpsiM) >- 0.08 || (deltaB0M - deltaPsiPM) < 0.09 ))))'.format(JPSIM=JPsiMass_, PSIM=PsiPMass_,  CUT=nSigma_psiRej)  
    elif dimusel == 'keepPsi':
      cut_base = '(abs(mumuMass - {JPSIM}) < {CUT}*mumuMassE || abs(mumuMass - {PSIM}) < {CUT}*mumuMassE)'.format( JPSIM=JPsiMass_, PSIM=PsiPMass_, CUT=nSigma_psiRej)
    elif dimusel == 'nocut':
      cut_base = cut_base
    else:
      print ('\nYou should define which dimuon mass to consider. Please choose between following options: \nkeepPsiP, keepJpsi, rejectPsi, keepPsi')
 
    return cut_base






def writeCMS(frame, year, ibin = [-1,-1], toy=0):

    txt = ROOT.TLatex(.11,.91,"CMS") 
    txt. SetNDC() 
    txt. SetTextSize(0.045) 
    frame.addObject(txt) 
    
    lumiYear = lumi_eras[str(year)]
    if toy==0:     txt2 = ROOT.TLatex(.75,.91,"%.1f fb^{-1}, 13 TeV"%lumiYear) 
    elif toy==1:   txt2 = ROOT.TLatex(.55,.91,"simulation, equivalent of %.1f fb^{-1}, 13 TeV"%lumiYear) 
    elif toy==-1:  txt2 = ROOT.TLatex(.55,.93,"simulation, 13 TeV") 
    else:          txt2 = ROOT.TLatex(.55,.91,"simulation, equivalent of %.1f fb^{-1}, 13 TeV"%toy) 
    txt2 . SetNDC() ;
    txt2 . SetTextSize(0.03) ;
    txt2 . SetTextFont(42) ;
    frame. addObject(txt2) ;
    
    if ibin[0] > -1:
        txtq = ROOT.TLatex(.16,.6, "%s GeV^{2} < q^{2} < %s GeV^{2}" %(ibin[0],ibin[1])) 
        txtq . SetNDC() ;
        txtq . SetTextSize(0.036) ;
        txtq . SetTextFont(42) 
        frame. addObject(txtq) 
       
    
def niceFrame(frame, title):
    frame.GetYaxis().SetTitleOffset(1.35)
    frame.SetTitle(title)
    try:
        frame.getAttText().SetTextSize(0.024) 
        frame.getAttText().SetTextFont(42) 
        frame.getAttLine().SetLineColor(0) 
    except:
        pass
            
def niceFrameLowerPad(frame, titleY):
    frame.GetYaxis().SetTitle(titleY)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetYaxis().SetTitleSize(0.06)
    frame.GetYaxis().SetLabelSize(0.06)
    frame.GetXaxis().SetLabelSize(0.06)
    frame.GetXaxis().SetTitleSize(0.07)
    frame.SetTitle('')

    