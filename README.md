# AP-scheme-for-LFP-equation
This project provides code for paper (). 
# Precompute fractional laplacian matrix
One can use following code to compute the computatioanl matrix for different parameters; we also provide all the precomputed matrix we used in the paper.
## s=0.5 (alpha=1)
```
M_pre_pl_alp1.m
```
## s!=0.5 and 0<s<1
```
M_pre_pl_new.m
```
# Compute fractional laplacian
Fo example
## Fig 1,2
```
compute_fl(128,3,0.4,300,2)
```
## Fig 3
```
FL_Error_vs_N(0.4,3,300,2)
```

# Spatially homogeneuous case
## For Fig~4,5,6, following example is about s=0.4
```
LF_hom(128,0.8,1,3,300,1)
```
# Spatially inhomogeneuous case
## When \epsilon =1, compare AP scheme with IMEX scheme
## For Fig~7,8, following example is about s=0.8, gaussian initial condition
```
LFP_AP_NHE_ssn(100,128,pi,3,0.8,1,300,1,0.05,0.5,2,1)
```

## Energy stability
For Fig~10, following example is for s=0.4, \epsilon=1
```
LF_hom(128,0.8,1,3,300,1)
LFP_AP_NHE_ssn(100,128,5,3,0.4,1,300,1,0.01,0.1,1,2)
LFP_Energy(100,128,0.8,0.01,0.1,1,1)
```

## AP property
For Fig~11,12, following example is for s=0.8, gaussian initial condition
```
LFP_AP_check_new(100,128,5,3,1,0.8,300,0.01,0.1,1)
```

