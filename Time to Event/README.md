# Time to Event Model with Smooth Dose-Response Curve

This [R walkthrough](https://berryconsultants.github.io/code-libraries/Time%20to%20Event/TTE_smooth.html) guides through data generation, prior calibration and posterior simulation for a piecewise exponential model with a NDLM for the underlying dose-response curve.

Highlights from this walkthrough include: 

* Review of survival analysis
* NDLM for smooth curves
* How to generate data from a piecewise exponential model
* Prior calibration for a piecewise exponential model
* Efficient model fitting

***Note***: the piecewise exponential model can be fitted as a Poisson regression with a reparametrization trick. 
However, in practice I have observed that the model fitting procedure proposed here (using the likelihood explicitly) is not only faster but is also mixing better when run with the same number of iterations. 
