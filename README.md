# Welcome to the Simulation Toolkit for Digital Clinical Measure Validation.

This toolkit contains a selection of tools that complement and support the Framework For Novel Digital Clinical Measures, and allow a user to:

1) Investigate the performance and characteristics of a novel method for assessing agreement between measures, through a simulated analytical validation study where the measures do nt have directly comparable units.
2) Understand how the novel method compares to established methods of assessing agreement in a simulated analytical validation study.
3) Examine the influence of common study design factors on the estimates gathered from these methods.
4) Experiment with factors of a realistic data generation mechanism to simulate your own novel digital measure scenarios.

# Contents of the Toolkit

There are five aspects to the toolkit:

1) R code for a simulation study of methods for analytical validation of a novel digital measure, where a digital measure that captures step count data is validated against a combination of COA-based reference measures that do not have directly comparable units to the digital measure.
2) A Shiny app that visualizes the results of running the simulation study code using the default simulation setup.
3) A manuscript concerning the results of the simulation study that results from running the simluation study code was coded using using the default setup. This manuscript details the methodology choices and underlying ideas that the simulation study was built on, in addition to the results of the study itself.
4) The full dataset produced by running the simulation study code using the default simulation setup. This is the data that was analyzed to produce the above manuscript.
5) A note on statistical methodology considerations for analytical validation studies. This note makes suggestions and recommendations for scenarios where measures do or do not have directly comparable units, builds on ideas explored in the simulation study manuscript, brings in additional validation concepts that may be familiar to investigators, and explores ideas that may bring a fresh perspective to analytical validation of digital measures.

# Code for simulation study

There are a number of R scripts that constitute the simulation study code. These are:

1) Final Simulation for Manuscript.R - This builds and runs the simulation study. The design condition matrix is created, the data generation mechanisms are created and the data is generated, the data is then analyzed using three methodologies (Pearson Correlation, Confirmatory Factor Analysis, Simple/Multiple Linear Regression), the performance measures are calculated (including empirical bias and empirical standard error), and the estimates and performance measures are collated, thensaved to CSV files for convenience of data manipulation. Note: results from each individual simulation condition are saved to a local location on your device, depending on your R environment working path.
2) Functions for Simulation Data Gen.R - This is a collection of support functions for the data generation mechanisms used in the main simulation script. NOTE: THIS SCRIPT MUST BE RUN FIRST, TO LOAD KEY FUNCTIONS REQUIRED FOR THE MAIN SIMULATION SCRIPT INTO YOUR R ENVIRONMENT.
3) Streamlined Plotting Collection.R - This produces a number of ggplots that illustrate key features of the simulation study results using the default simulation setup.

# How to run the simulation study

1) Run "Functions for Simulation Data Gen.R"
2) Run "Final Simulation for Manuscript.R"
3) Run "Streamlined Plotting Collection.R"

Note that with the default settings, the simulation will complete only 10 replications per simulation condition. We start with 10 replication to allow you to run an initial pass of the simulation and code in a reasonable timeframe, to test the code in your enviroment and get a sense of the outputs. We recommend 500 replications in order for the Monte Carlo Standard Error of the simulation to be acceptable, but note that running the full simulation with 500 replications is likely to take several hours to complete.

# Levers of the simulation, and modifying the simulation to fit your analytical validation scenario

The simulation study code contains a suite of parameters that can be adjusted, to modify the simulation and allow you to investigate broader scenarios than the default scenario presented here (a digital measure capturing daily summary step count data, validated against three COA-based reference measures that do not have directly comparable units to the digital measure).

Firstly, we have the parameters of the simulation: the parameters that are varied and tested in the full factorial simulation. They are listed below, in the format of variable name : description of variable. 

1) N: sample size
2) meas_error_mag: magnitude of the measurement error in the digital measure
3) n_assess: number of repeated assessments being included in the analysis,
4) missing_method: the digital measure data missingness method
5) missing_rate: the proportion of data missingness in the digital measure data.

To modify the way these parameters are varied, edit the design condition matrix named "Design" in "Functions for Simulation Data Gen.R". (Note: when modifying the data missingness method, the design matrix only lists strings that describe the missingness method - to enact your missingness method, you will also need to add code to the if-else block that deals with data missingness.). For further details on how these parameters control the data generation mechanism, please see the manuscript in this Toolkit, or the comments in main simulation script in this Toolkit.

Further to the parameters that make up the simulation design condition matrix, there are a number of other parameters which are fixed across all simulation conditions and each repetition in the simulation. These parameters control various aspects of the digital measure and reference measure data generation. Adjusting these parameters will allow you to tailor the simulation study to your specific analyical validation scenario. A list of these parameters can be found below, in the format of variable name : description of variable. 

So-called "fixed" values that can also be modified by a user for additional flexibility:
    1. fluct_sd : Controls the daily fluctuation in an individual's latent trait. Larger values increase the daily fluctuations.
    2. meth_filt_sd : Controls the method filter, the ability of the digital measure to observe the latent trait. Larger values decrease the digital measure's observational ability.
    3) base_rate : the hypothesized mean of the digital measure for an individual from the population with mean physical ability. Default setting is 10000 to model step count.
    3) latent_effect : the proportional effect of an individual's latent physical ability on their expected digital measure count. Increasing this value increases the effect of being "non-average" on the digital measure data.
    4) per_filt_sd : Controls the perception filter for a given reference measure, i.e. the imperfect ability of an individual to perceive their mean latent trait value. Larger values decrease the reference measure's observational ability.
    5) b0/b2 : Controls the difficulty threshold parameters for the simulated weekly recall COA data. This data is simulated using Item Response Theory, in a Graded Response model.
    6) reliability: Controls the reliability of the weekly recall COA instruments (Note: to control further parameters of the graded response models, adjust the functions "Generate_4_12_IRT_parameters" and "Generate_5_7_ClinRO_IRT_parameters" in "Functions for Simulation Data Gen.R"
    7) imic/thrds : Controls the distribution of the transition thresholds between response categories for the daily recall PRO reference measure.
     
For full details on how these parameterss control the simulation, please see the manuscript in this Toolkit, or the comments in main simulation script in this Toolkit.
    
# How to use the Shiny visualization app

1) Downloead and run "Shiny Visualization app.R"
2) Select sample size (Note that CFA functions poorly with small sample sizes, so we advise setting N>=35
3) Select magnitude of measurement error (as a fraction of the latent trait's effect - enter values between 0.5 and 2.0)
4) Select data missingness rate from the drop down menu
5) Select number of repetitions for each simulation condition (Note: large values may take a long time to run)
6) Click "Run simulation"
7) Summary statistics and plots of the mean empirical bias for each statistical method are displayed. To change how the plots are grouped, use the "Group By:" drop down menu



