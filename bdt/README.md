Scripts to perform signal/background discrimination via bdt, single candidate per-event selection, and application of the B0Psi cut.  
The following steps should be performed:

#### BDT training/testing/WP optimisation
create_subsample.C  
          |   
optimize_hp.py  (only needed once to optimize the bdt hyperparameters)  
          |  
final_bdt_sub_samples.py     (only needed once to train the bdt)  
          |  
add_bdt_subsamples.py       
          |    
wp_optimization/wp_optimization.py  (only needed when searching for the optimal bdt working point)  
          |  
wp_optimization/draw_wp_optimization.py  (only needed when searching for the optimal bdt working point)   
          |  
selectSingleCand.py            
          |     
applyB0PsiCut.py   
          |      
addVarsData.py    

----------------------------------------
More specifically, the workflow to add the bdt score to your samples and to select your B0 candidate is the following :  
```
mkdir sub_samples
root -l 'create_subsample(year, type, mc).C'
```
where year can be 2016, 2017 or 2018,  
type is defined [here](https://github.com/CMSKStarMuMu/selection_and_fits/blob/master/bdt/create_subsample.C#L26-L30),  
and mc is 0 or 1 for data or MC respectively.  

```
python add_bdt_subsamples.py year
python selectSingleCand.py year
python applyB0PsiCut.py year
python addVarsData.py year
```

Please note that in each script there is a list defining which samples you would like to process, e.g. [here](https://github.com/CMSKStarMuMu/selection_and_fits/blob/master/bdt/add_bdt_subsamples.py#L23-L31), so you may wish to modify the list/add new samples there.
