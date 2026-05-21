# SRWAggregationFormation

Scripts and data accompanying manuscript submitted to PNAS,  investigating the spatiotemporal dynamics of SRWs' coastal habitat use and site fidelity during five decades of population recovery 


# README

Dataset accompanying research article 'Post-collapse habitat selection: density-dependent behavioral shifts governed the persistence of a migratory mammal population’ submitted to PNAS. Please familiarize yourself with the contents of the article and its contexts before using attached datasets and R code. Any use of the data or code should refer to the original research article:

A. Grundlehner, [ADD CITATION]

 * Common abbreviations *
CC = Cow-calf pair
UA = Unaccompanied individual


# R Scripts

1) '...R' For calculating and plotting Ripley's K coefficient for spatial clustering
2) '...R' To simulate/calculate the maximum population growth rate based on fixed assumptions.
3) '...R' Contains code to run cluster analyses in order to identify aggregation areas
4) '...R' To model temporal trends in the annual abundance of CC, UA and Totals, and calculate growth rates.


# Data files

1)	'SRWAbundance_byArea.csv'

This file contains annual SRW abundance as used for statistical analysis in the manuscript. It contains the following columns:

AggregationArea
Description: Name of the aggregation area, corresponding to the sites distinguished in the manuscript (see Figure 1).

Year
Description: Year of observation

Sum_FemaleCalfPairs 
Description: The sum of female-calf pairs (number of CC pairs) counted in the aggregation area in the focal year. Note that this is the sum of female-calf pairs, not the number of individuals (which can be derived by multiplying the number of pairs by two).

Sum_UnaccompaniedAnimals
Description: The sum of unaccompanied animals (number of individuals) counted in the aggregation area in the focal year. Unaccompanied animals are all SRW individuals that are not a calf or female with a calf. This includes males, juveniles, and females without a calf.


