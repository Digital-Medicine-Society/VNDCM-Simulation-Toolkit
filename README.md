# Welcome to the Simulation Toolkit for Digital Clinical Measure Analytical Validation

This toolkit contains a selection of tools that complement and support the Framework For Novel Digital Clinical Measures, and allow a user to:

1) Investigate the performance and characteristics of a method for assessing agreement between measures when they do not have directly comparable units, through a simulated analytical validation study.
2) Understand how this method compares to established methods of assessing agreement in a simulated analytical validation study.
3) Examine the influence of common study design factors on the estimates generated by this approach.
4) Experiment with model factors to simulate your own novel digital clinical measure scenarios.

## Contents of the Toolkit

There are four components to the toolkit:

1) R code to run a simulation study for analytical validation of a novel digital clinical measure, where a digital clinical measure that captures step count data is validated against a combination of COA-based reference measures that do not have directly comparable units to the digital clinical measure.
2) An app that visualizes the results of the simulation study code using the default simulation parameters.
3) The full dataset produced by running the simulation study code using the default simulation parameters. This is the data that was analyzed to produce the manuscript linked below.

The methodology behind this toolkit is described in [this manuscript](https://todo).

More information on this project that generated this toolkit can be accessed [on the website of the Digital Medicine Society](https://datacc.dimesociety.org/validating-novel-digital-clinical-measures/).

## Code for simulation study

The toolkit contains the following R scripts:

1) `Final Simulation for Manuscript.R` - This builds and runs the simulation study. The design condition matrix is created, the data generation mechanisms are created and the data is generated, the data is then analyzed using three methodologies (Pearson Correlation, Confirmatory Factor Analysis, Simple/Multiple Linear Regression), the performance measures are calculated (including empirical bias and empirical standard error), and the estimates and performance measures are collated, and finally saved to CSV files for convenience of data manipulation. 
> NOTE: results from each individual simulation condition are saved to a local location on your device, depending on your R environment working path.
2) `Functions for Simulation Data Gen.R` - This is a collection of support functions for the data generation mechanisms used in the main simulation script. 
> IMPORTANT: This script must be run first to load key functions required for the main simulation script into your R environment.
3) `Streamlined Plotting Collection.R` - This produces a number of ggplots that illustrate key features of the simulation study results using the default simulation parameters.

## How to run the simulation study

1) Run `Functions for Simulation Data Gen.R`
2) Run `Final Simulation for Manuscript.R`
3) Run `Streamlined Plotting Collection.R`

Note that with the default settings, the simulation will complete only 10 replications per simulation condition. We start with 10 replication to allow you to run an initial pass of the simulation and code in a reasonable timeframe, to test the code in your enviroment, and familiarze yourself with the outputs. We recommend 500 replications in order for the Monte Carlo Standard Error of the simulation to be acceptable, which is likely to take several hours to complete.

## Parameters of the simulation, and modifying the simulation to fit your analytical validation scenario

The simulation study code contains several parameters that can be adjusted to modify the simulation and allow you to investigate additional scenarios. 

The default scenario is a digital clinical measure capturing daily summary step count data, validated against three COA-based reference measures that do not have directly comparable units to the digital clinical measure.

First, consider the main simulation parameters that are varied and tested in the full factorial simulation (`parameter name`: parameter description). 

1) `N`: sample size
2) `meas_error_mag`: magnitude of the measurement error in the digital clinical measure
3) `n_assess`: number of repeated assessments included in the analysis
4) `missing_method`: the digital clinical measure data missingness method
5) `missing_rate`: the proportion of data missingness in the digital clinical measure data

To modify the way these parameters are varied, edit the design condition matrix named `Design` in `Functions for Simulation Data Gen.R`.

> NOTE: When modifying the data missingness method, the design matrix only lists strings that describe the missingness method - to enact your missingness method, you will also need to add code to the if-else block that deals with data missingness. For further details on how these parameters control the data generation mechanism, please see the [manuscript](https://todo), and the comments in main simulation script.

Second, there several additional parameters which are fixed across all simulation conditions and each repetition in the simulation. These parameters control various aspects of the digital clinical measure and reference measure data generation. Adjusting these parameters will allow you to tailor the simulation study to your specific analyical validation scenario.

1) `fluct_sd`: Controls the daily fluctuation in an individual's latent trait. Larger values increase the daily fluctuations.
2) `meth_filt_sd`: Controls the method filter, the ability of the digital clinical measure to observe the latent trait. Larger values decrease the digital clinical measure's observational ability.
3) `base_rate`: The hypothesized mean of the digital clinical measure for an individual from the population with mean physical ability. The default setting is 10000 to model step count.
4) `latent_effect`: The proportional effect of an individual's latent physical ability on their expected digital clinical measure count. Increasing this value increases the effect of being "non-average" on the digital clinical measure data.
5) `per_filt_sd`: Controls the perception filter for a given reference measure, i.e., the imperfect ability of an individual to perceive their mean latent trait value. Larger values decrease the reference measure's observational ability.
6) `b0/b2`: Controls the difficulty threshold parameters for the simulated weekly recall COA data. This data is simulated using Item Response Theory, in a Graded Response model.
7) `reliability`: Controls the reliability of the weekly recall COA instruments. To control additional parameters of the graded response models, adjust the functions `Generate_4_12_IRT_parameters` and `Generate_5_7_ClinRO_IRT_parameters` in `Functions for Simulation Data Gen.R`.
8) `imic/thrds`: Controls the distribution of the transition thresholds between response categories for the daily recall PRO reference measure.
     
For full details on how these parameterss control the simulation, please see the [manuscript](https://todo) and the comments in main simulation script.
    
## How to use the visualization app

1) Download and run `Visualization app.R`.
2) Select sample size; note that CFA functions poorly with small sample sizes, so we advise setting N>=35.
3) Select magnitude of measurement error (as a fraction of the latent trait's effect - enter values between 0.5 and 2.0).
4) Select data missingness rate from the drop down menu.
5) Select number of repetitions for each simulation condition; note: large values may take a long time to run.
6) Click `Run simulation`.
7) Summary statistics and plots of the mean empirical bias for each statistical method are displayed. To change how the plots are grouped, use the `Group By:` drop down menu.



