# PowerSwitchingTrial
==========
* [Introduction](#introduction)
* [Installation](#installation)
* [Citation](#citation)
* [Example](#example)

<a name="introduction"/>

# Introduction

PowerSwitchingTrial is a R package to estimate power and sample sizes in clinical trials with treatment switching.

<a name="installation"/>

# Installation

```R

library(devtools)
install_github("darwin-hub/PowerSwitchingTrial")
```


<a name="citation"/>

# Citation

 Lejun Deng, Hsu Chih-Yuan, Shyr Yu. https://doi.org/10.1101/2022.09.07.22279686
 

<a name="example"/>

# Example

To find the power of a clinical trial where the extimated median survival time of the experimental treatment group is 1.5x the control group. It is assumed that all patients are recruitted at the start of the clinical trial, and the clinical trial doesn't end. The proportion of patients who switch is 0.6, and the switching time is also the estimated median survival time of the control group. The power required is 0.8

```R
LogRankTestMix2PowerMedian(n=500, m1=1, m2=1.5, r=1, Ta=0, Te=999, proportion=0.6, alpha=0.05, 
                                       s.dist="beta", pt=0.3, rho=0.5,
                                       censor.rate=c("AC.only"),
                                       reps=5000)
```
```R
LogRankTestMix2NMedian(power=0.8, lower=100, upper=300, m1=1, m2=1.5, r=1, Ta=0, Te=4, proportion=0.3, alpha=0.05, s.dist="beta", pt=0.3, rho=0.5,
                       censor.rate=c("AC.only"),
                       reps=5000)
```
