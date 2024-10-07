# Deuteron Peak Fit for NMR Signal Analysis

Simple fitting program for spin-1 NMR signals with quadrupolar splitting, based on C code by C. Dulya. [[1]](#1)

Uses Non-Linear Least-Squares Minimization and Curve-Fitting for Python (LMFIT), https://lmfit.github.io/lmfit-py/


### Use
May be called as an object, passing a list of frequency points, a list signal magnitudes as each point, and a parameter dictionary.
```
fit_results = DFits(freqs, sweep, params)

r = fit_results.result.params['r'].value
pol = (r * r - 1) / (r * r + r + 1)
```

The success of the fit is highly dependent on the initial parameters passed.

## References

<a id="1">[1]</a> 
Dulya, C. et. al.  
"A line-shape analysis for spin-1 NMR signals"
NIM A, 398, 109-125 (1997). (https://doi.org/10.1016/S0168-9002(97)00317-3)
