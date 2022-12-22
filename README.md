# PowerSwitchingTrial
==========
* [Introduction](#introduction)
* [Installation](#installation)
* [Citation](#citation)
* [Usage](#Usage)

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

 Lejun Deng, Chih-Yuan Hsu, Yu Shyr . https://doi.org/10.1101/2022.09.07.22279686
 

<a name="Usage"/>

# Usage

To estimate the power of a clinical trial with treatment switching at a significance level of 0.05 (alpha=0.05), where the trial has 500 patients in each group (n=500), the allocation ratio r of 1, the median survival time of the control group and the experimental group is 1 and 1.5 (m1=1, m2=1.5), respectively, the switching probability of 0.3 (proportion=0.3), the option of switching time to be "beta" (s.dist="beta"), the ratio of average switching time and the event time to be 0.3, the correlation between the switching time and the event time of 0.5. 

```R
library(PowerSwitchingTrial)
LogRankTestMix2PowerMedian(n=500, m1=1, m2=1.5, r=1, Ta=0, Te=4, proportion=0.3, alpha=0.05, 
                                       s.dist="beta", pt=0.3, rho=0.5,
                                       censor.rate=c("AC.only"),
                                       reps=5000)
```

To determine the sample size required to achieve 80% power at a significance level of 0.05 and other settings described as above. 
```R
LogRankTestMix2NMedian(power=0.8, lower=100, upper=300, m1=1, m2=1.5, r=1, Ta=0, Te=4, 
                       proportion=0.3, alpha=0.05, s.dist="beta", pt=0.3, rho=0.5,
                       censor.rate=c("AC.only"),
                       reps=5000)
```
