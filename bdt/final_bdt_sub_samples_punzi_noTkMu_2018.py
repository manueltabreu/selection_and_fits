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

from ROOT      import TFile, TTree, TH1F, gROOT, TChain

# mpl.interactive(True)

mc_sigma = 0.0400 ## was 49
mc_mass  = 5.27783 

gROOT.SetBatch(True)

## weight input background mass
min_m  = 4.8
max_m  = 5.8
n_bins = 500


##########################################################################################
#####   SIGNAL AND BACKGROUND SELECTION
##########################################################################################
sig_mass_str = '( (tagB0==1 && (bMass > {M} - 2.5*{S} && bMass< {M} + 2.5*{S})) || \
                  (tagB0==0 && (bBarMass > {M} - 2.5*{S} && bBarMass < {M} + 2.5*{S})))'.format( M=mc_mass,S=mc_sigma)     
truth_match_str = '(truthMatchMum==1 && truthMatchMup ==1 && truthMatchTrkm==1 && truthMatchTrkp ==1 )'
ct_str = '( (tagB0==1 && genSignal==1) || (tagB0==0 && genSignal==2))'

bkg_mass_str = '( (tagB0==1 && ((bMass> {M}-7.*{S} && bMass < {M}-3.*{S}) || (bMass > {M}+3.*{S} && bMass < {M}+7*{S}))) || \
                  (tagB0==0 && ((bBarMass > {M}-7.*{S} && bBarMass < {M}-3.*{S}) || (bBarMass > {M}+3.*{S} && bBarMass < {M}+7*{S})))) ' .format( M=mc_mass,S=mc_sigma) 


print '----- consider mass preselection which is different in 2016 and 2018 --------'
sig_selection_cutbased = sig_mass_str + ' && ' + \
                         truth_match_str + ' && ' +\
                         ct_str + ' && ' +\
                         'trig == 0 && (mumuMass < 2.702)'

bkg_selection_cutbased = bkg_mass_str + ' && ' + \
                         '(mumuMass < 2.702) '



# sig_selection_cutbased = '(((tagB0==1 && (bMass    > {M} - 2.5*{S} && bMass    < {M} + 2.5*{S})) || \
#                             (tagB0==0 && (bBarMass > {M} - 2.5*{S} && bBarMass < {M} + 2.5*{S}))) &&\
#                              truthMatchMum==1 && truthMatchMup ==1 && truthMatchTrkm==1 && truthMatchTrkp ==1 && \
#                             ((tagB0==1 && genSignal==1) | (tagB0==0 && genSignal==2)) && \
#                             trig == 0 && \
#                             ((mumuMass < 2.702)  ))'.format( M=mc_mass,S=mc_sigma)     
# # 
# bkg_selection_cutbased = '((tagB0==1 && ((bMass    > {M}-7.*{S} && bMass    < {M}-3.*{S}) || (bMass    > {M}+3.*{S} && bMass    < {M}+7*{S}))) || \
#                            (tagB0==0 && ((bBarMass > {M}-7.*{S} && bBarMass < {M}-3.*{S}) || (bBarMass > {M}+3.*{S} && bBarMass < {M}+7*{S}))) && \
#                            ((mumuMass < 2.702)  )) '.format(M=mc_mass,S=mc_sigma)
# 

sig_list = [
          'sub_samples/sample_2018_MC_LMNR_0_newphi.root', 
          'sub_samples/sample_2018_MC_LMNR_1_newphi.root',
          'sub_samples/sample_2018_MC_LMNR_2_newphi.root', 
          'sub_samples/sample_2018_MC_LMNR_3_newphi.root', 
          'sub_samples/sample_2018_MC_LMNR_4_newphi.root', 
          'sub_samples/sample_2018_MC_LMNR_5_newphi.root', 
          'sub_samples/sample_2018_MC_LMNR_6_newphi.root', 
          'sub_samples/sample_2018_MC_LMNR_7_newphi.root', 
          'sub_samples/sample_2018_MC_LMNR_8_newphi.root', 
          'sub_samples/sample_2018_MC_LMNR_9_newphi.root', 
          'sub_samples/sample_2018_MC_LMNR_10_newphi.root',
]

bkg_list = [
          'sub_samples/sample_2018_data_LMNR_0_newphi.root',
          'sub_samples/sample_2018_data_LMNR_1_newphi.root',
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

for isample in range(11):

    tag = '_punzi_removeTkMu_fixBkg_2018_%s'%isample

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
        'kstTrkmMinIP2D', ## min impact parameter of the track from any PV (in 2D)
        'kstTrkpMinIP2D', ## min impact parameter of the track from any PV (in 2D)
    ]
    
    branches = features + [
        'bMass',
        'bBarMass',
        'tagB0',
        'mumNTrkLayers',
        'mupNTrkLayers',
        'mumNPixLayers',
        'mupNPixLayers',
        'mupdxyBS',
        'mumdxyBS',
        'mumdzBS',
        'mupdzBS',
        'mupHighPurity',
        'mumHighPurity',
        'mumTMOneStationTight',
        'mupTMOneStationTight',
        'mupPt',
        'mumPt',
        'mumuVtxCL',
        'kstVtxCL',
        'kstTrkpTrackerMuonArbitrated',
        'kstTrkmTrackerMuonArbitrated',
        'kstTrkpTrackerMuon',
        'kstTrkmTrackerMuon',
        'kstTrkpHighPurity',
        'kstTrkmHighPurity',
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
        'mumIsoPt_dr04',
        'mupIsoPt_dr04',
        'kstTrkmIsoPt_dr04',
        'kstTrkpIsoPt_dr04',
        'charge_trig_matched'
    ]

    branches = list(set(branches))
    
    sig = pandas.DataFrame(
        root_numpy.root2array(
            sig_list[:isample] + sig_list[(isample + 1):], 
            'ntuple',
            branches  = branches + ['weight'],
            selection = sig_selection_cutbased,
        )
    )
    
    bkg = pandas.DataFrame(
        root_numpy.root2array(
            bkg_list[:isample] + bkg_list[(isample + 1):], 
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
    
    bkg['trkpDCASign']  = abs(bkg.kstTrkpDCABS/bkg.kstTrkpDCABSE)
    sig['trkpDCASign']  = abs(sig.kstTrkpDCABS/sig.kstTrkpDCABSE)
    bkg['trkmDCASign']  = abs(bkg.kstTrkmDCABS/bkg.kstTrkmDCABSE)
    sig['trkmDCASign']  = abs(sig.kstTrkmDCABS/sig.kstTrkmDCABSE)
    branches.append('trkmDCASign')
    branches.append('trkpDCASign')
    
    ## define isolation: # tracks with pt in a cone
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
    
    
    ## use only one combination for K* mass
    sig['kstarmass']  = sig.tagB0*sig.kstMass +(1- sig.tagB0)*sig.kstBarMass
    bkg['kstarmass']  = bkg.tagB0*bkg.kstMass +(1- bkg.tagB0)*bkg.kstBarMass
    features.append('kstarmass')
    
    ## add b mass calculation for plotting only
    bkg['themass']  = bkg.bMass*bkg.tagB0 + bkg.bBarMass*(1-bkg.tagB0)
    sig['themass']  = sig.bMass*sig.tagB0 + sig.bBarMass*(1-sig.tagB0)
    branches.append('themass')
    
    bkg.hist(column='themass',bins=400)
    plt.savefig('mass_bkg.pdf')
    sig.hist(column='themass',bins=400)   
    plt.savefig('mass_sig.pdf')
    
#     bkg['m_weight']  = add_m_weights_bkg(bkg.themass)
#     sig['m_weight']  = add_m_weights_sig(sig.themass)
    ## add normfactor here
    
    sig['normfactor'] = sig.weight#*sig.m_weight #/len(sig) *10000
    bkg['normfactor'] = 1. ##bkg.m_weight# /len(bkg) #* 10000
    
#     sig['normfactor'] = sig.weight*sig.m_weight /len(sig)*1000
#     bkg['normfactor'] = bkg.m_weight /len(bkg)*1000
    
    ##########################################################################################
    #####   SPLIT TRAIN AND TEST SAMPLES
    ##########################################################################################
    data_all = pandas.concat([sig, bkg])
    
    ## add column for pass-preselection
    data_all['pass_preselection'] = ( data_all.mumTMOneStationTight == 1 ) & ( data_all.mupTMOneStationTight == 1 ) & \
                                    ( data_all.kkMass > 1.035 ) & \
                                    ( data_all.kstTrkmTrackerMuon == 0 ) & ( data_all.kstTrkpTrackerMuon == 0 ) & \
                                    (~((data_all.kstTrkmGlobalMuon == 1) & ( data_all.kstTrkmNTrkLayers > 5 ) & ( data_all.kstTrkmNPixHits > 0))) & \
                                    (~((data_all.kstTrkpGlobalMuon == 1) & ( data_all.kstTrkpNTrkLayers > 5 ) & ( data_all.kstTrkpNPixHits > 0))) & \
                                    (( (data_all.charge_trig_matched ==  1) & (data_all.kstTrkpPt > 1.2) & (data_all.trkpDCASign > 2) ) | \
                                     ( (data_all.charge_trig_matched == -1) & (data_all.kstTrkmPt > 1.2) & (data_all.trkmDCASign > 2) ) )

    
    data = data_all[data_all.pass_preselection==1]
    train, test = train_test_split(data, test_size=0.30, random_state = 17)
    
#     import pdxb; pdb.set_trace()
# {'params': {'colsample_bytree': 0.5819271135450363, 'learning_rate': 0.2157230517930242, 
# 'max_delta_step': 0.061506541557647364, 'min_child_weight': 7.033350795238926, 
# 'n_estimators': 577.721545092955, 'subsample': 0.5919674072752289, 'max_depth': 6.375707865120728, 
# 'gamma': 0.776500014626179}, 'target': 570.7153757775731}
    clf = XGBClassifier(
        colsample_bytree = 0.580,
        learning_rate    = 0.21,
        max_delta_step   = 0.06150,
        min_child_weight = 7.03,
        n_estimators     = 578,
        max_depth        = 6,
        subsample        = 0.59,
        gamma            = 0.78,     
        seed             = 1986,
        silent           = False,
        )

    clf.fit(
        train[features], 
        train.target,
        eval_set              = [(train[features], train.target), (test[features], test.target)],
        early_stopping_rounds = 50,
        eval_metric           = 'auc',
        verbose               = True,
        sample_weight         = train['normfactor'],
    )
    
    joblib.dump(clf, 'results/classifier_%s.pkl' %(tag), compress=True)
    
    ##########################################################################################
    # #####   PREDICT ON THE TEST UNBIAS SAMPLE
    # ##########################################################################################
    pred = clf.predict_proba(test[features])[:, 1] #### !!! changed from test_unbias
    
    train_sig = clf.predict_proba(train[features][train.target>0.5])[:,1]
    train_bkg = clf.predict_proba(train[features][train.target<0.5])[:,1]
    
    test_sig = clf.predict_proba(test[features][test.target>0.5])[:,1]
    test_bkg = clf.predict_proba(test[features][test.target<0.5])[:,1]
    
    ##########################################################################################
    #####   ROC CURVE
    ##########################################################################################
    plt.clf()
    
    import itertools
    xy = [i*j for i,j in itertools.product([10.**i for i in range(-8, 0)], [1,2,4,8])]+[1]
    plt.plot(xy, xy, color='grey', linestyle='--')
    plt.xlim([10**-3, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    
    ## draw baseline point
    #   plt.plot([fpr], [tpr], label='baseline', markerfacecolor='red', marker='o', markersize=10)
    # 
    # ## draw ROC         
    fpr, tpr, threshold = roc_curve(test.target, (pred * test.pass_preselection) ) ### !! changed from test unbias
    plt.plot(fpr, tpr, color='r', label='ROC test')
    
    pred2 = clf.predict_proba(train[features])[:, 1]
    fpr2, tpr2, sara2 = roc_curve(train.target, pred2)
    plt.plot(fpr2, tpr2, color='b', label='ROC train')

    plt.xscale('log')
    plt.grid()
    
    plt.legend(loc='best')
    plt.grid()
    plt.title('ROC')
    plt.tight_layout()
    plt.savefig('results/roc_train_test_%s.pdf' %(tag))
    plt.clf()

    roc_file = open('roc_%s.pck' %(tag), 'w+')
    pickle.dump((tpr, fpr), roc_file)
    roc_file.close()
    
    ### add plot of BDT vs tpr      ############################
#     plt.plot(threshold, tpr, color='b')
#     
#     plt.legend(loc='best')
#     plt.grid()
#     plt.title('turn on')
#     plt.tight_layout()
#     plt.xlabel('BDT output')
#     plt.ylabel('True Positive Rate')
#     plt.savefig('results/turn_on_%s.pdf' %(tag))
#     plt.clf()
#     
#     roc_file = open('results/turnon_%s.pck' %(tag), 'w+')
#     pickle.dump((threshold, tpr), roc_file)
#     roc_file.close()
    
    ### add plot of BDT vs fpr      ############################
#     plt.plot(threshold, fpr, color='b')
#     plt.legend(loc='best')
#     plt.grid()
#     plt.title('rate')
#     plt.xlabel('BDT output')
#     plt.ylabel('False Positive Rate')
#     plt.yscale('log')
#     
#     plt.tight_layout()
#     plt.savefig('results/rate_%s.pdf' %(tag))
#     plt.clf()
#     
#     roc_file = open('results/rate_%s.pck' %(tag), 'w+')
#     pickle.dump((threshold, fpr), roc_file)
#     roc_file.close()
    ##########################################################################################
    #####   OVERTRAINING TEST
    ##########################################################################################
    
    low  = 0
    high = 1
    low_high = (low,high)
    bins = 50
    
    #################################################
    plt.hist(
        train_sig,
        color='r', 
        alpha=0.5, 
        range=low_high, 
        bins=bins,
        histtype='stepfilled', 
        normed=True,
        log=True,
        label='S (train)'
    )
    
    #################################################
    plt.hist(
        train_bkg,
        color='b', 
        alpha=0.5, 
        range=low_high, 
        bins=bins,
        histtype='stepfilled', 
        normed=True,
        log=True,
        label='B (train)'
    )
    
    #################################################
    hist, bins = np.histogram(
        test_sig,
        bins=bins, 
        range=low_high, 
        normed=True,
    )
    
    width  = (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    scale  = len(test_sig) / sum(hist)
    err    = np.sqrt(hist * scale) / scale
    
    plt.errorbar(
        center, 
        hist, 
        yerr=err, 
        fmt='o', 
        c='r', 
        label='S (test)'
    )
    
    #################################################
    hist, bins = np.histogram(
        test_bkg,
        bins=bins, 
        range=low_high, 
        normed=True
    )
    
    width  = (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    scale  = len(test_bkg) / sum(hist)
    err    = np.sqrt(hist * scale) / scale
    
    plt.errorbar(
        center, 
        hist, 
        yerr=err, 
        fmt='o', 
        c='b', 
        label='B (test)'
    )
    
    #################################################
    plt.xlabel('BDT output')
    plt.ylabel('Arbitrary units')
    plt.legend(loc='best')
    ks_sig = ks_2samp(train_sig, test_sig)
    ks_bkg = ks_2samp(train_bkg, test_bkg)
    plt.suptitle('KS p-value: sig = %.3f%s - bkg = %.2f%s' %(ks_sig.pvalue * 100., '%', ks_bkg.pvalue * 100., '%'))
    
    # plt.tight_layout()
    plt.savefig('results/overtrain_%s.pdf' %(tag))
    
    ##########################################################################################
    #####   FEATURE IMPORTANCE
    ##########################################################################################
    plot_importance(clf)
    plt.tight_layout()
    plt.savefig('results/feat_importance_%s.pdf' %(tag))
    
    plt.close()        
    
    ##########################################################################################
    #####   ROC TRAIN
    ##########################################################################################
#     plt.clf()
#     
#     xy = [i*j for i,j in itertools.product([10.**i for i in range(-8, 0)], [1,2,4,8])]+[1]
#     plt.plot(xy, xy, color='grey', linestyle='--')
#     plt.xlim([10**-3, 1.0])
#     plt.ylim([0.0, 1.05])
#     plt.xlabel('False Positive Rate')
#     plt.ylabel('True Positive Rate')
#     
#     plt.xscale('log')
#     
#     pred = clf.predict_proba(train[features])[:, 1]
#     # fpr, tpr, _ = roc_curve(train.target, pred)
#     fpr, tpr, sara = roc_curve(train.target, pred)
#     plt.plot(fpr, tpr, label='BDT', color='b')
#     
#     plt.legend(loc='best')
#     plt.grid()
#     plt.title('ROC')
#     plt.tight_layout()
#     plt.savefig('results/roc_train_%s.pdf' %(tag))
    
    ##########################################################################################
    #####   OVERTRAINING SCORE
    ##########################################################################################
    plt.clf()
    
    auc_train = clf.evals_result()['validation_0']['auc']
    auc_test  = clf.evals_result()['validation_1']['auc']
    
    n_estimators = np.arange(len(auc_train))
    
    plt.plot(n_estimators, auc_train, color='r', label='AUC train')
    plt.plot(n_estimators, auc_test , color='b', label='AUC test' )
    
    plt.xlabel('# tree')
    plt.ylabel('Area Under ROC')
    
    plt.xscale('log')
    plt.grid()
    
    # plt.xlim([1, 1000])
    # plt.ylim([0.985, 1.0])
    
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('results/auc_score_%s.pdf' %(tag))
    
    
    
    ######################################
    # Correlation Matrix Plot
    ######################################
    sig  ['bdt'] = clf.predict_proba(sig[features])[:, 1] 
    bkg  ['bdt'] = clf.predict_proba(bkg[features])[:, 1] 

    sig_corr = sig[features+ ['bdt', 'themass']]
    bkg_corr = bkg[features+ ['bdt', 'themass']]
    
    sig_correlations = sig_corr.corr()
    bkg_correlations = bkg_corr.corr()
    
    fig = plt.figure()
    ticks = np.arange(0,len(features)+2,1)
    
    ax  = fig.add_subplot(121)
    cax = ax.matshow(sig_correlations, vmax=1., vmin=-1, cmap=plt.cm.get_cmap('Blues', 9))


    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(features+ ['bdt', 'themass'],fontsize=5)
    ax.set_yticklabels(features+ ['bdt', 'themass'],fontsize=5)
    plt.xticks(rotation=90)
    
    # plot correlation matrix
    ax2  = fig.add_subplot(122)
    cax2 = ax2.matshow(bkg_correlations, vmax=1., vmin=-1, cmap=plt.cm.get_cmap('Blues', 11))
    ax2.set_xticks(ticks)
    ax2.set_yticks(ticks)
    ax2.set_xticklabels(features+ ['bdt', 'themass'],fontsize=5)
    ax2.set_yticklabels(features+ ['bdt', 'themass'],fontsize=0)
    plt.xticks(rotation=90)
    
    cbar = fig.colorbar(cax, orientation="horizontal", pad=0.045)
    cbar.ax.tick_params(labelsize=5)
    # cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8]) 
    # cb = plt.colorbar(cax, cax = cbaxes)  
    # fig.colorbar(cax)
    # fig.colorbar(cax2)
    
    plt.show()
    plt.savefig('results/bkg_correlation_%s.pdf' %(tag))
    
    plt.close()
    