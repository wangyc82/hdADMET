# hdADMET

**H**eterogeneous network and **D**eep neural network - based **ADMET** prediction. 

## Status

Active development

## Introduction

Systematical prediction of ADMET is a critical step in drug discovery. Among various obstacles hindering clinical translation, lacking effective multi-task methods for integration inner-relationships among ADMET features has become a bottleneck. hdADMET incorporated compound-drug relationships, the experimental validated ADMET features of drugs, and the interrelationships among ADMET properties to transform the embedding representations into a kernel-based DNN for boosting the learning ability of prediction model.

## Usage

1. Installation

   Prerequisites of hdADMET includes the following: 

   - R is properly installed; 

   - Rscript is available in your system path ($PATH);

   - git (2.21.1)

    Installation of hdADMET includes the following steps:

    - step 1: git clone https://github.com/wangyc82/hdADMET;

    - step 2: download the example data (example-data.RData) from hdADMET repository, and put it in the hdADMET folder.

    Dependencies of hdADMET includes the following: 

    - word2vec package and its dependencies.
    - h2o package and its dependencies.

    Testing of successful installation by running the following commands in R:
     
       > library(word2vec)
       > library(h2o)


2. Preparation of the input files
The folliwng three matrix need be to prepared for running hdADMET

Acd (adjacent matrix representing the interactions between new compounds and known clinical drugs)

Adp (adjacent matrix representing the validated ADMET properties for clinical drugs)

App (adjacent matrix representing ADMET inner-relationships)

take example data for instance

     > load('~/hdADMET/example-data.RData')

example data includes the MACCs fingerprints for new compounds (MACCs.fp.mat) and anticancer drugs (MACCs.fp.mat.CA)
adjacent matrix representing drug-ADMET interaction,adj_ADMET_CA1
adjacent matrix representing ADMET inter-relationships, ADMET_cor that were obtained by the ADMET profiles of drugs

Following the instruction of hdADMET.R to prepare Acd, Adp, and App, and make sure compound and ADMET features has no any delimeter

set len1, len2, len3, len4, len5 accoding the dimension of Acd, Adp, App

For this example data, set len1=10000;len2=1e5;len3=1e5;len4=1000;len5=5*1e+5.

3. Running hdADMET

The main function of hdADMET is hdADMET.R. Get your input files prepared, and run it like this:

Usage example:

    > source('~/hdADMET/hdADMET.R')
    > predictions<-CIPHEN(Acd,Adp,App,len1,len2,len3,len4,len5)
    

## Contact

For technical issues please send an email to wangyongcui@mail.kib.ac.cn.
