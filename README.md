# Genetic-Variation-and-Mutation-Patterns-in-Giant-Kelp
Scripts for article "Genetic Variation and Mutation Patterns in Giant Kelp (Macrocystis pyrifera) Insights from a Comprehensive Analysis"

This project contains the main scripts to analyze the VCF file result from SnpEff after the file has been decomposed. 

Mutation_table_from_SnpEff_output.py:
Creates an easy-to-handle mutation table from the decomposed VCF file output from SnpEff. 
This table is used to generate several other files that become the input for "Figures_and_Analysis.R"

Statistics_pop_gen_per_gene:
This is a bash/R script that is used to calculate all the population genetics statistics used in the paper. 
It starts by defining the limits of a gene, generating a VCF file for that particular gene using BCFtools, then it utilizes R to calculate statistics.

associate_gene_function_to_mutations.py:
This file helps analyze mutations in different KOG annotated pathways. 
Utilizing the Mutation table, separates different KOG pathways and their genes, then calculates the amount of every type of mutation in each individual.

Figures_and_Analysis.R:
R file to make every figure and graph on the paper. 
It describes the input files required for each figure, all of them made form the Mutation Table.


