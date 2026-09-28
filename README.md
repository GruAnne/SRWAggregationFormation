# SRWAggregationFormation

Scripts and data accompanying manuscript accepted for publication in PNAS (2026),  investigating the spatiotemporal dynamics of SRWs' coastal habitat use and site fidelity during five decades of population recovery 


# README

Dataset accompanying research article 'Post-collapse habitat selection: density-dependent behavioral shifts governed the persistence of a migratory mammal population’ submitted to PNAS. Please familiarize yourself with the contents of the article and its contexts before using attached datasets and R code. Any use of the data or code should refer to the original research article:

[[CITATION TBA]]


Common abbreviations

CC = Cow-calf pair

UA = Unaccompanied individual


# R Scripts

1) 'RipleyK.R' For calculating and plotting Ripley's K coefficient for spatial clustering
2) 'ClusterAnalysis.R' Contains code to run cluster analyses in order to identify aggregation areas
3) 'GAMs_GrowthAbundanceDensDep.R' To model temporal trends in the annual abundance of CC, UA and Totals, and calculate growth rates; and to test density dependent affects

#' Description of the scripts and required datasets

1) 'RipleyK.R'
Description: Calculate Ripley's K 

Required datasets:
"fg_grid0.1_AggrAreaAttributed.xlsx"

2) 'ClusterAnalysis.R'
Description: Run 2-step cluster analysis (see Main text and Supporting information of research article) to identify aggregation areas

Required datasets:
"fg_grid0.1_AggrAreaAttributed.xlsx"

3) 'GAMs_GrowthAbundanceDensDep.R'
Description: This script contains code to fit and validate the GAM models for the temporal trends of individual aggregation areas; it contains code for plotting area-specific trends in abundances and growth rates; And to run statistical tests to check whether annual abundances (regional and local; tested for different groups) affect the occupation of the areas (density-dependent effects).

Required datasets:
"clusteranalysis_polygons_aggregation_areas.shp"; "fg_area.csv"


# Data files

1)	'SRWAbundance_byArea.csv'

This file contains annual SRW abundance as used for statistical analysis in the manuscript. It contains the following columns:

Column: AggregationArea

Description: Name of the aggregation area, corresponding to the sites distinguished in the manuscript (see Figure 1).


Column: Year

Description: Year of observation


Column: Sum_FemaleCalfPairs 

Description: The sum of female-calf pairs (number of CC pairs) counted in the aggregation area in the focal year. Note that this is the sum of female-calf pairs, not the number of individuals (which can be derived by multiplying the number of pairs by two).


Column: Sum_UnaccompaniedAnimals

Description: The sum of unaccompanied animals (number of individuals) counted in the aggregation area in the focal year. Unaccompanied animals are all SRW individuals that are not a calf or female with a calf. This includes males, juveniles, and females without a calf.


2) 'fg_grid0.1_AggrAreaAttributed.xlsx'

This file contains gridded SRW abundance data. This data is already filtered to only contain annual maximum values per grid cell per year.


Column: Year

Description: Year of observation


Column: AggregationArea

Description: Name of the aggregation area, corresponding to the sites distinguished in the manuscript (see Figure 1).


Column: LonBin

Description: Longitude of grid cell


Column: LatBin

Description: Latitude of grid cell


Column: calf	

Description: Total number of calves in the grid cell


Column: total	

Description: Total number of SRWs in the grid cell


Column: unacc

Description: Total number of unaccompanied individuals in the grid cell


4) 'clusteranalysis_polygons_aggregation_areas.shp'

This is a shapefile containing the outline of the identified aggregation areas; can also be made using ClusterAnalysis.R

3) 'fg_area_GH.csv'

This file contains annual maximum abundances for each individual aggregation area (always observed in Aug/Sept, see Methods section in the manuscript). This data is already filtered to only contain annual maximum values per area per year.


Column: Year

Description: Year of observation


Column: AggregationArea

Description: Name of the aggregation area, corresponding to the sites distinguished in the manuscript (see Figure 1).


Column: calf	

Description: Total number of calves in the grid cell


Column: total	

Description: Total number of SRWs in the grid cell


Column: unacc

Description: Total number of unaccompanied individuals in the grid cell

