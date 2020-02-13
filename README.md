# The ancestral *C. elegans* cuticle suppresses *rol-1*

![rol-1 tracks](assets/BE8_longest_circgt0.9_pretty.png?raw=true)

This repository contains `R` code and data to reproduce analyses and figures. Here's the [preprint](https://www.biorxiv.org/content/10.1101/2020.02.07.938696v1).

S1: This file list. 

data/S2_UncRolRecombinantGenotypes.txt: N2/CB4856 genotype frequencies, described in Methods, `Observation of segregation distortion`.

videos/S3-S9: sample videos of mutant and suppressed worms.

data/S10_col182_orthologs.csv: coding sequence coordinates based on [*Caenorhabditis* Genomes Project](http://caenorhabditis.org/) data and denovo assembly for *elegans* wild isolates. See Methods, `Gene model`.

data/S11_col182_CDS.mfa: multiple alignment of coding sequences, for the species in File S10. See Methods, `Gene model`.

data/S12_merged_roller_MWT.rda: `R` data archive containing raw, aggregated Multi-Worm Tracker data. See Methods, `Quantitative locomotion analysis`, and references therein. 

data/S13_merged_roller_MWT_trackStat.rda: `R` data archive containing track-level summary statistics derived from the raw Multi-Worm Tracker data. See Methods, `Quantitative locomotion analysis`, and references therein. 

code/S14_col182_mwt.R: `R` code to analyse Multi-Worm Tracker data. See Methods, `Quantitative locomotion analysis`, and references therein. 

data/S15_CeMEE_phenoGeno.rda: `R` data archive containing CeMEE genotypes, phenotypes and utility functions used for testing genetic interactions with the N2 *col-182* allele. See Methods, `Epistasis analysis in the CeMEE`, and references therein.

code/S16_col182_CeMEE_interactions.R: `R` code containing functions used for testing genetic interactions with the N2 *col-182* allele. See Methods, `Epistasis analysis in the CeMEE`, and references therein.



