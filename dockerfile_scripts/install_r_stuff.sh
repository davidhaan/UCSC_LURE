#!/bin/bash

#chrisw 20190709
#install various R packages for running LURE

R -e 'install.packages(c("cowplot","gridExtra","survival","survminer","ggpubr","magrittr","igraph","RcppGreedySetCover"))'
R -e 'install.packages(c("optparse","ggplot2","PRROC","dplyr","doMC","iterators","data.table","matrixStats"))'
R -e 'install.packages(c("plyr","glmnet","foreach","Matrix"))'

#R -e 'source("https://bioconductor.org/biocLite.R")'
#R -e 'BiocInstaller::biocLite(c("ComplexHeatmap","FDb.InfiniumMethylation.hg19"))'

