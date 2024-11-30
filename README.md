# VNDCM-Simulation-Toolkit

Welcome to the Simulation Toolkit for Digital Clinical Measure Validation.

This toolkit contains a selection of tools that complement and support the Framework For Novel Digital Clinical Measures, and allow a user to:

1) Investigate the performance and characteristics of a novel method for assessing agreement between measures, through a simulated analytical validation study where the measures do nt have directly comparable units.
2) Understand how the novel method compares to established methods of assessing agreement in a simulated analytical validation study.
3) Examine the influence of common study design factors on the estimates gathered from these methods.
4) Experiment with factors of a realistic data generation mechanism to simulate your own novel digital measure scenarios.

Contents of the Toolkit
=======================

There are four aspects to the toolkit:

1) R code for a simulation study of methods for analytical validation of a novel digital measure, where a digital measure that captures step count data is validated against a combination of COA-based reference measures that do not have directly comparable units to the digital measure.
2) A Shiny app that visualizes the results of the simulation study in the code from point #1, using the default simulation settings.
3) A manuscript concerning the results of the simulation study that was coded using point #1. This manuscript details the methodology choices and underlying ideas that the simulation study was built on, in addition to the results of the study itself. 
4) A note on statistical methodology considerations for analytical validation studies. This note makes suggestions and recommendations for scenarios where measures do or do not have directly comparable units, builds on ideas explored in the simulation study manuscript, brings in additional validation concepts that may be familiar to investigators, and explores ideas that may bring a fresh perspective to analytical validation of digital measures.


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
