#!/bin/tcsh
source  /cvmfs/cms.cern.ch/cmsset_default.csh
eval `scram runtime -csh`

set year = 2016
set tag_string='_newphi_punzi_removeTkMu_fixBkg_B0Psicut_fixPres'
foreach i ( 'MC_JPSI' 'MC_LMNR' 'MC_PSI' 'data' )
# foreach i ( 'MC_JPSI' 'MC_LMNR' 'MC_PSI' 'MC_BS' 'MC_BSJPSIPHI' 'MC_BSJPSIKST' 'data' )
# foreach i ( 'MC_JPSI' 'MC_LMNR' 'MC_PSI' 'MC_HBJPSIX' 'data' )
#     set nsub = `echo $nsub_all:q | sed 's/ / /g'`
#     set nroot="`ls ntuples/charmonium_fixl144_${year}${i}/*.root | wc -l`"
#     if ($nsub[2] != $nroot)  then
#         echo 'probably a root file is missing'
#     endif

    echo "hadd ../final_ntuples/${year}${i}${tag_string}.root ../final_ntuples/${year}${i}${tag_string}_part*.root"
    hadd ../final_ntuples/${year}${i}${tag_string}.root ../final_ntuples/${year}${i}${tag_string}_part*.root
#     hadd /gwteray/users/fiorendi/p5prime/data20${year}/flat_ntuples/Oct2020/ntuple_CharmoniumR_20${year}${i}_fixl144.root ntuples/charmonium_fixl144_${year}${i}/*.root
end
