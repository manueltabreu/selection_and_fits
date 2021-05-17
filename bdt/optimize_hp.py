import os, sys
# sys.path.insert(0, os.environ['HOME'] + '/.local/lib/python2.7/site-packages')

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
# mpl.use('TkAgg')
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt
import pandas, root_numpy
from copy import deepcopy
import argparse
import pickle
import sklearn
from   sklearn.externals import joblib
from   sklearn.ensemble  import GradientBoostingClassifier
from   sklearn.metrics   import roc_curve
from   sklearn.model_selection import train_test_split
from   pdb import set_trace
from   array       import array

from   xgboost import XGBClassifier, plot_importance

from ROOT  import TFile, TTree, TH1F, gROOT
from math  import sqrt 

from bayes_opt import BayesianOptimization
from bayes_opt.observer import JSONLogger
from bayes_opt.event import Events

# mpl.interactive(True)

mc_sigma = 0.0400 ## was 49
mc_mass  = 5.27783 

a_punzi = 2
b_punzi = 5

lumi_mc   = 8256 * (10./11) * 1.95
lumi_data = 61.1 * (10./11)

gROOT.SetBatch(True)


def xgboost_fom(max_depth,
                learning_rate,
                n_estimators,
                gamma,
                min_child_weight,
                max_delta_step,
                subsample,
                colsample_bytree,
              ):

    clf = XGBClassifier(
        max_depth        = int(max_depth),
        learning_rate    = learning_rate,
        n_estimators     = int(n_estimators),
        subsample        = subsample,
        colsample_bytree = colsample_bytree,
        min_child_weight = min_child_weight,
        gamma            = gamma,     ## optimized
        max_delta_step   = max_delta_step,
#         reg_alpha        = 0,
#         reg_lambda       = 1,
        seed             = 1986,
        silent           = True,
        )
    
    clf.fit(
        train[features], 
        train.target,
        eval_set              = [(train[features], train.target), (test[features], test.target)],
        early_stopping_rounds = 50,
        eval_metric           = 'auc',
        verbose               = False,
        sample_weight         = train['normfactor'],
    )

    true_bkg = test[(test.target==0) & (test.pass_preselection==1)]
    true_sig = test[(test.target==1) & (test.pass_preselection==1)]
    
    true_bkg['bdt_score']  = clf.predict_proba(true_bkg[features])[:, 1] ## should be probability to be signal (bdt_score) for each data entry
    true_sig['bdt_score']  = clf.predict_proba(true_sig[features])[:, 1] ## should be probability to be signal (bdt_score) for each data entry
    

    punzis = {}
    print ('score \t S (on test only) \t B (on test s.) \t punzi')
    for wp in np.arange(0.8,1,0.01):
        n_sig = float( len(true_sig[true_sig.bdt_score > wp])) 
        n_bkg = float( len(true_bkg[true_bkg.bdt_score > wp])) 
        den = pow(a_punzi,2)/8 + 9*pow(b_punzi,2)/13 + a_punzi*sqrt(n_bkg) + 0.5*b_punzi*sqrt(pow(b_punzi,2) + 4*a_punzi*sqrt(n_bkg) + 4*n_bkg)
        if den > 0:
            punzis[wp] =  n_sig/ den
        else:
            punzis[wp] = 0.0
        
        print ('%.2f\t'%wp, n_sig, '\t', n_bkg, '\t', punzis[wp])
    return max(punzis.values())

#     significances = {}
#     print 'score \t S (on test only) \t B (on test s.) \t significance'
#     for wp in np.arange(0.5,1,0.02):
#         n_sig = float( len(true_sig[true_sig.bdt_score > wp])) * lumi_data/lumi_mc 
#         n_bkg = float( len(true_bkg[true_bkg.bdt_score > wp])) * 5./8 ## to account that we use 8 sigma region for bkg and 5 sigma regin for signal
#         if n_sig > 0:
#             significances[wp] =  n_sig/ sqrt(n_sig + n_bkg) 
#         else:
#             significances[wp] = 0.0
#         
#         print wp, '\t', n_sig, '\t', n_bkg, '\t', significances[wp]
#             
#     return max(significances.values())

##########################################################################################
#####   SIGNAL AND BACKGROUND SELECTION
##########################################################################################
sig_selection_cutbased = '((tagB0==1 && (bMass    > {M} - 2.5*{S} && bMass    < {M} + 2.5*{S})) || \
                           (tagB0==0 && (bBarMass > {M} - 2.5*{S} && bBarMass < {M} + 2.5*{S}))) &&\
                            truthMatchMum==1 && truthMatchMup ==1 && truthMatchTrkm==1 && truthMatchTrkp ==1 && \
                           ((tagB0==1 && genSignal==1) || (tagB0==0 && genSignal==2)) && \
                           ( kstTrkmTrackerMuon == 0 ) && ( kstTrkpTrackerMuon == 0 ) && \
                          trig == 0 && \
                           (mumuMass < 2.702) '.format( M=mc_mass,S=mc_sigma)

bkg_selection_cutbased = '((tagB0==1 && ((bMass    > {M}-7.*{S} && bMass    < {M}-3.*{S}) || (bMass    > {M}+3.*{S} && bMass    < {M}+7*{S}))) || \
                           (tagB0==0 && ((bBarMass > {M}-7.*{S} && bBarMass < {M}-3.*{S}) || (bBarMass > {M}+3.*{S} && bBarMass < {M}+7*{S}))) && \
                           ( kstTrkmTrackerMuon == 0 ) && ( kstTrkpTrackerMuon == 0 ) && \
                           (mumuMass < 2.702)) '.format(M=mc_mass,S=mc_sigma)


sig_list = [
#           'sub_samples/sample_MC_LMNR_1.root', ### need to include PU profile
          'sub_samples/sample_2018_MC_LMNR_0_newphi.root',  ### need to include PU profile
          'sub_samples/sample_2018_MC_LMNR_2_newphi.root',  ### need to include PU profile
          'sub_samples/sample_2018_MC_LMNR_3_newphi.root',  ### need to include PU profile
          'sub_samples/sample_2018_MC_LMNR_4_newphi.root',  ### need to include PU profile
          'sub_samples/sample_2018_MC_LMNR_5_newphi.root',  ### need to include PU profile
          'sub_samples/sample_2018_MC_LMNR_6_newphi.root',  ### need to include PU profile
          'sub_samples/sample_2018_MC_LMNR_7_newphi.root',  ### need to include PU profile
          'sub_samples/sample_2018_MC_LMNR_8_newphi.root',  ### need to include PU profile
          'sub_samples/sample_2018_MC_LMNR_9_newphi.root',  ### need to include PU profile
          'sub_samples/sample_2018_MC_LMNR_10_newphi.root', ### need to include PU profile
]

bkg_list = [
#           'sub_samples/sample_data_LMNR_1.root',
          'sub_samples/sample_2018_data_LMNR_0_newphi.root',
          'sub_samples/sample_2018_data_LMNR_2_newphi.root',
          'sub_samples/sample_2018_data_LMNR_3_newphi.root',
          'sub_samples/sample_2018_data_LMNR_4_newphi.root',
          'sub_samples/sample_2018_data_LMNR_5_newphi.root',
          'sub_samples/sample_2018_data_LMNR_6_newphi.root',
          'sub_samples/sample_2018_data_LMNR_7_newphi.root',
          'sub_samples/sample_2018_data_LMNR_8_newphi.root',
          'sub_samples/sample_2018_data_LMNR_9_newphi.root',
          'sub_samples/sample_2018_data_LMNR_10_newphi.root',
]


tag = '_optimize_splitIso'

##########################################################################################
#####   FEATURES AND BRANCHES
##########################################################################################
features = [
    'bCosAlphaBS',
    'bLBS/bLBSE',
    'kstTrkmDCABS/kstTrkmDCABSE',
    'kstTrkpDCABS/kstTrkpDCABSE',
    'bVtxCL',
    'bDCABS/bDCABSE',
    'kstTrkmMinIP2D', 
    'kstTrkpMinIP2D', 
]

branches = features + [
    'bMass',
    'bBarMass',
    'tagB0',
    'mumTMOneStationTight',
    'mupTMOneStationTight',
    'mumuVtxCL',
    'kstVtxCL',
    'kstTrkpPt',
    'kstTrkmPt',
    'kstTrkmDCABSE',
    'kstTrkmDCABS',
    'kstTrkpDCABS',
    'kstTrkpDCABSE',
    'kstMass',
    'kstBarMass',
    'kkMass',
    'bLBS',
    'bLBSE',    
    'kkMass',  ## could be moved to preselection
    'kstTrkmGlobalMuon',
    'kstTrkmNTrkLayers',
    'kstTrkmNPixHits',
    'kstTrkpGlobalMuon',
    'kstTrkpNTrkLayers',
    'kstTrkpNPixHits',
    'mumuMass',
    'charge_trig_matched',
    'mumIsoPt_dr04',
    'mupIsoPt_dr04',
    'kstTrkmIsoPt_dr04',
    'kstTrkpIsoPt_dr04',
    'mupPt',
    'mumPt',
    'kstTrkpPt',
    'kstTrkmPt',
    'kstTrkpTrackerMuon',
    'kstTrkmTrackerMuon',
]

branches = list(set(branches))

sig = pandas.DataFrame(
    root_numpy.root2array(
        sig_list, 
        'ntuple',
        branches  = branches + ['weight'],
        selection = sig_selection_cutbased,
    )
)

bkg = pandas.DataFrame(
    root_numpy.root2array(
        bkg_list, 
        'ntuple',
        branches  = branches,
        selection = bkg_selection_cutbased,
    )
)

##########################################################################################
#####   DEFINE THE TARGETS
##########################################################################################

sig['target'] = np.ones (sig.shape[0]).astype(np.int)
bkg['target'] = np.zeros(bkg.shape[0]).astype(np.int)

## use only one combination for K* mass
sig['kstarmass']  = sig.tagB0*sig.kstMass +(1- sig.tagB0)*sig.kstBarMass
bkg['kstarmass']  = bkg.tagB0*bkg.kstMass +(1- bkg.tagB0)*bkg.kstBarMass
features.append('kstarmass')

## define jsonlation: # tracks with pt in a cone
sig['isopt_mum_04']  = sig.mumIsoPt_dr04/sig.mumPt
sig['isopt_mup_04']  = sig.mupIsoPt_dr04/sig.mupPt
sig['isopt_trkm_04'] = sig.kstTrkmIsoPt_dr04/sig.kstTrkmPt
sig['isopt_trkp_04'] = sig.kstTrkpIsoPt_dr04/sig.kstTrkpPt
sig['sum_isopt_04']  = sig.isopt_mum_04 + sig.isopt_mup_04 + sig.isopt_trkm_04 + sig.isopt_trkp_04
bkg['isopt_mum_04']  = bkg.mumIsoPt_dr04/bkg.mumPt
bkg['isopt_mup_04']  = bkg.mupIsoPt_dr04/bkg.mupPt
bkg['isopt_trkm_04'] = bkg.kstTrkmIsoPt_dr04/bkg.kstTrkmPt
bkg['isopt_trkp_04'] = bkg.kstTrkpIsoPt_dr04/bkg.kstTrkpPt
bkg['sum_isopt_04']  = bkg.isopt_mum_04 + bkg.isopt_mup_04 + bkg.isopt_trkm_04 + bkg.isopt_trkp_04
features.append('sum_isopt_04')

## add b mass calculation for plotting only
bkg['themass']  = bkg.bMass*bkg.tagB0 + bkg.bBarMass*(1-bkg.tagB0)
sig['themass']  = sig.bMass*sig.tagB0 + sig.bBarMass*(1-sig.tagB0)
branches.append('themass')

## add b mass calculation for plotting only
bkg['trkpDCASign']  = abs(bkg.kstTrkpDCABS/bkg.kstTrkpDCABSE)
sig['trkpDCASign']  = abs(sig.kstTrkpDCABS/sig.kstTrkpDCABSE)
bkg['trkmDCASign']  = abs(bkg.kstTrkmDCABS/bkg.kstTrkmDCABSE)
sig['trkmDCASign']  = abs(sig.kstTrkmDCABS/sig.kstTrkmDCABSE)
branches.append('trkmDCASign')
branches.append('trkpDCASign')

sig['normfactor'] = sig.weight
bkg['normfactor'] = 1. 


##########################################################################################
#####   SPLIT TRAIN AND TEST SAMPLES
##########################################################################################
data_all = pandas.concat([sig, bkg])

## add column for pass-preselection
data_all['pass_preselection'] = ( data_all.mumTMOneStationTight == 1 ) & ( data_all.mupTMOneStationTight == 1 ) & \
                                ( data_all.kkMass > 1.035 ) & \
                                (~((data_all.kstTrkmGlobalMuon == 1) & ( data_all.kstTrkmNTrkLayers > 5 ) & ( data_all.kstTrkmNPixHits > 0))) & \
                                (~((data_all.kstTrkpGlobalMuon == 1) & ( data_all.kstTrkpNTrkLayers > 5 ) & ( data_all.kstTrkpNPixHits > 0))) & \
                                ( data_all.kstTrkmTrackerMuon == 0 ) & ( data_all.kstTrkpTrackerMuon == 0 ) & \
                                (( (data_all.charge_trig_matched ==  1) & (data_all.kstTrkpPt > 1.2) & (data_all.trkpDCASign > 2) ) | \
                                 ( (data_all.charge_trig_matched == -1) & (data_all.kstTrkmPt > 1.2) & (data_all.trkmDCASign > 2) ) )


data = data_all[data_all.pass_preselection==1]
train, test = train_test_split(data, test_size=0.3, random_state = 17)

xgboostBO = BayesianOptimization(xgboost_fom,
                                 {
                                  'max_depth': (2, 10),
                                  'learning_rate': (0.01, 0.4),
                                  'n_estimators': (200, 1000),
                                  'subsample': (0.4, 0.8),
                                  'colsample_bytree' :(0.5, 0.99),
                                  'min_child_weight': (2, 10),
                                  'gamma': (0., 1.),
                                  'max_delta_step': (0, 0.2),
                                  }
                                )
## save optimization steps
logger = JSONLogger(path="logs_punzi_noTrkMuTrks.json")
xgboostBO.subscribe(Events.OPTMIZATION_STEP, logger)

xgboostBO.maximize(
    n_iter = 30,
    init_points = 20,
)

print('-'*53)

print('Final Results')
print(xgboostBO.max)

print '\n ---- all iterations ------ \n'
for i, res in enumerate(xgboostBO.res):
    print("Iteration {}: \n\t{}".format(i, res))