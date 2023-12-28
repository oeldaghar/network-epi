# network-epi


This is the accompanying repo for: **Multi-scale Local Network Structure Critically Impacts
Epidemic Spread and Interventions**

It contains functions for performing SIR and SEIR simulations with a neighbor quarantining policy. We model stochastic simulations using a Discrete Event Simulation that only tracks updates in node states. 

The other notable computation is for computing the Network Community Profile (NCP) for both seeded PageRank and epidemic simulations. 

## Usage
See juypter notebooks in the **examples/** directory for examples of code usage. 

## Citations 
Scripts have citations where the code is based on the work of others. If you believe some portion requires a citation, please bring this to our attention. Some notable portions of code that deserve mentioning are:
1. Computing the NCP from seeded PageRank solution vectors
2. Code for producing HexBins

The first is mostly code that I modifed from David Gleich (https://github.com/dgleich). 
The second was modified from https://github.com/RalphAS/HexBinPlots.jl for the missed sets computation.

## Data Sources
For sources of data, please see the accompanying paper for appropriate citations.
