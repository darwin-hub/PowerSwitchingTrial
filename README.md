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
LogRankTestMedian(m1=1,m2=1.5,Ta=0,Te=999,reps=5000,alpha=0.05,r=1,lower=50,upper=300,power=0.8)

 N  K
98 11
In the above, the initial number of trial patients using treatment 1 that is required is 98. The function performed binary search K times
LogRankTestMedian(m1=1,m2=1.5,Ta=0,Te=999,reps=5000,alpha=0.05,r=1,lower=50,upper=70,power=0.8)
pwLower   0.5190
pwUpper   0.6556
lower    50.0000
upper    70.0000
In the above, the range of the lower and upper numbers did not include the power value we specified (0.8) Looking at the results, we can see that they were both smaller than 0.8. If pwLower is greater than the specified power, decrease the lower argument. If pwUpper is smaller than the specified power, increase the upper argument.
