Scripts to perform signal/background discrimination via bdt, single candidate per-event selection, and application of the B0Psi cut.  
The following steps should be performed:

#### BDT training/testing/WP optimisation
create_subsample.C  
          |  
          |   
optimize_hp.py  
          |  
          |  
final_bdt_sub_samples.py     
          |  
          |  
add_bdt_subsamples.py       
          |    
          |  
wp_optimization/wp_optimization.py  
          |  
          |  
wp_optimization/draw_wp_optimization.py  
     

#### Apply BDT selection/select one candidate per event/apply B0PsiCut
selectSingleCand.py  
          |    
          |    
applyB0PsiCut.py  

  
