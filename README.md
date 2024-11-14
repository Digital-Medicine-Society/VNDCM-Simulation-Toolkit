# VNDCM-Simulation-Toolkit


This repository allows a user to run a simulation for an analytical validation study with a digital measure that captures step count data, and is validated against a combination of COA-based reference measures. 

Explain the files: Functions... contains a collection of helper functions that need to be loaded into your environment prior to running your sim.
                   Final Sim... runs the simulation
                   Collating the estimates... does what it says - it loads the final results (which are saved to a local location by dim condition) and does two things: 1) collates the estimates from each individual sim repetition into one data table, and 2) creates a summary table of e.g. mean bias by sim condition. This is useful for data manipulation.

Expain the levers of the simulation, and what those variables are in the code.
    Sample size
    Missing data - mechanism and rate
    Number of repeated assessments
    Magnitude of measurement error (as a proportion of the Poisson rate)

So-called "fixed" values that can also be modified by a user for additional flexibility
    Daily latent fluctuation factor
    Method Filter
    Perception filter
