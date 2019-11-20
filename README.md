# selection_and_fits

## Train BDT and add bdt score to the ntuples  
Instructions to come

## Fit mass spectrum
The [perBin_massFit.py](https://github.com/CMSKStarMuMu/selection_and_fits/blob/master/perBin_massFit.py) script performs the fit:  
- correctly tagged MC events  
- wrongly tagged MC events  
- data, constraining the parameters (Gaussian constrain) to the results from the MC fits above.  

NB: in order to run this script, you first need to compile the RooDoubleCBFast pdf (once for all).
```
cd utils/func_roofit/  
make
cd ../../
python perBin_massFit.py [mumuSelection] [year]
```
