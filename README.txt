bs_orthoreg_lgcen.pro is an IDL program that determines the confidence interval of the slope, intercept, and intrinsic scatter of a line fitted to logarithmic values (e.g., luminosities) using an orthogonal regression methodology in the presence of censored data. 

Orthogonal regression minimizes the perpendicular distance between the fitted line and the data. It is best suited to problems where either X or Y could be the dependent or independent variable. 

A bootstrap analysis is undertaken to estimate the confidence interval.

MPFIT is used as the algorithm to minimize the inverse of the likelihood functions.

The likelihood function is calculated using the function orthreg_lgcens.pro, which implements the likelihood functions from Pihajoki 2017 arxiv:1704.05466v2, specifically equations 33-35 and B5-B6, depending on the type of censoring found in each data point. While likelihood functions are available for censoring in either X or Y, we do not have one for a data point that is simultaneously censored in both directions and therefore effectively exclude it by assigning it a high likelihood so that it contributes very little to the inverse likelihood.

The inputs for bs_orthoreg_lgcen are: 
- dfile: an ASCII data file containing up to 5 pairs of columns containing value 
		and flag (0=value, 1 = upper limit)
- nboots: an integer indicating the number of bootstrap trials to run
- nme_out: a string identifier for this dataset

It creates as output an IDL savefile containing 
  - fit_bs: the results of each mpfit run
  - fit_stats: a distillation of the best fit and confidence intervals
  - the necessary inputs to reproduce the fit if necessary

Fit_stats is also output to the screen and, optionally, the fits can be plotted.


Acknowledgements: 

If you use this software for your work, please cite the following two papers

- Lanz, L., Hickox, R.C., Balokovic, M. et al. 2018, Astrophysical Journal, 
	"A Joint Study of X-ray and IR Reprocessing in Obscured Swift/BAT AGN"

- Pihajoki, P. 2017 MNRAS, 472, 3407, "A geometric approach to non-linear correlations with intrinsic scatter"
