# AP-scheme-for-LFP-equation
This project provides code for paper (). 
# Precompute fractional laplacian matrix
## s=0.5 (alpha=1)
```
M_pre_pl_alp1.m
```
## s!=0.5 and 0<s<1
```
M_pre_pl_new.m
```
# Compute fractional laplacian
## Fig 1,2
```
compute_fl
```
## Fig 3
```
FL_Error_vs_N(s,L,l_lim,IC)
```

# Spatially inhomogeneuous case
## generate Fig~4,5,6, following example is about s=0.4
```
LF_hom(128,0.8,1,3,300,1)
```

