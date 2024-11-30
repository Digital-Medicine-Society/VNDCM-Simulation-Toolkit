# Welcome to the Simulation Toolkit for Digital Clinical Measure Validation.

This toolkit contains a selection of tools that complement and support the Framework For Novel Digital Clinical Measures, and allow a user to:

1) Investigate the performance and characteristics of a novel method for assessing agreement between measures, through a simulated analytical validation study where the measures do nt have directly comparable units.
2) Understand how the novel method compares to established methods of assessing agreement in a simulated analytical validation study.
3) Examine the influence of common study design factors on the estimates gathered from these methods.
4) Experiment with factors of a realistic data generation mechanism to simulate your own novel digital measure scenarios.

# Contents of the Toolkit

There are four aspects to the toolkit:

1) R code for a simulation study of methods for analytical validation of a novel digital measure, where a digital measure that captures step count data is validated against a combination of COA-based reference measures that do not have directly comparable units to the digital measure.
2) A Shiny app that visualizes the results of the simulation study in the code from point #1, using the default simulation setup.
3) A manuscript concerning the results of the simulation study that was coded using point #1. This manuscript details the methodology choices and underlying ideas that the simulation study was built on, in addition to the results of the study itself. 
4) A note on statistical methodology considerations for analytical validation studies. This note makes suggestions and recommendations for scenarios where measures do or do not have directly comparable units, builds on ideas explored in the simulation study manuscript, brings in additional validation concepts that may be familiar to investigators, and explores ideas that may bring a fresh perspective to analytical validation of digital measures.

# Code for simulation study

There are a number of R scripts that constitute the simulation study code. These are:

1) Final Simulation for Manuscript.R - This builds and runs the simulation study. The design condition matrix is created, the data generation mechanisms are created and the data is generated, the data is then analyzed using three methodologies (Pearson Correlation, Confirmatory Factor Analysis, Simple/Multiple Linear Regression), the performance measures are calculated (including empirical bias and empirical standard error), and the estimates and performance measures are collated and saved to CSV files.
2) Functions for Simulation Data Gen.R - This is a collection of support functions for the data generation mechanisms used in the main simulation script. NOTE: THIS SCRIPT MUST BE RUN FIRST, TO LOAD KEY FUNCTIONS REQUIRED FOR THE MAIN SIMULATION SCRIPT INTO YOUR R ENVIRONMENT.
3) Streamlined Plotting Collection.R - This produces a number of ggplots that illustrate key features of the simulation study results using the default simulation setup.

# How to run the simulation study

1) Run "Functions for Simulation Data Gen.R"
2) Run "Final Simulation for Manuscript.R"
3) Run "Streamlined Plotting Collection.R"

Note that with the default settings, the simulation will complete only 10 replications per simulation condition. We start with 10 replication to allow you to run an intial pass of the simulation and code in a reasonable timeframe. We recommend 500 replications in order for the Monte Carlo Standard Error of the simulation to be acceptable, but note that running the full simulation with 500 replications is likely to take several hours to complete.



Explain the files: Functions... contains a collection of helper functions that need to be loaded into your environment prior to running your sim.
                   Final Sim... runs the simulation
                   Collating the estimates... does what it says - it loads the final results (which are saved to a local location by dim condition) and does two things: 1) collates the estimates from each individual sim repetition into one data table, and 2) creates a summary table of e.g. mean bias by sim condition. This is useful for data manipulation.

Expain the levers of the simulation, and what those variables are in the code.
    Sample size
    Missing data - mechanism and rate
    Number of repeated assessments
    Magnitude of measurement error (as a proportion of the Poisson rate)

So-called "fixed" values that can also be modified by a user for additional flexibility:
    Daily latent fluctuation factor
    Method Filter
    Perception filter
