# LURE_positive_controls.R
# This script will run LURE on the positive controls for the LURE manuscript



run_timestamp<-strftime(Sys.time(),"%Y_%m_%d_%H_%M")

print("Running LURE Positive Controls")

INPUT<-"./input/"
TEMP<-"./temp/"
SCRIPTS<-"./scripts/"
OUTPUT<-"./output/"
source(paste(sep="",SCRIPTS,"LURE_functions.R"))
registerDoMC(detectCores()-1)

# Run the SF3B1 in UVM Positive Control

# Load feature data for UVM (gene expression data)
feature_data<-data.frame(fread(paste(INPUT,"pancan_RNAexp_UVM",sep=""),stringsAsFactors = FALSE),row.names = 1)
tissue<-"UVM"
# Run LURE using the default settings
LURE(bait_gene="SF3B1-SET1_missense",
     gmt_file="positive_control_SF3B1_missense.gmt",
     feature_data=feature_data,
     num_permutations=5,
     max_num_events=5,
     percent_overlap=0,
     LURE_pvalue_threshold=.05,
     min_gene_set_size=5,
     gsea_pvalue_threshold=.05,
     gsea_fdr_threshold=.25,
     max_tree_length = 5,
     folds=10,
     enrichment_analysis_only=FALSE,
     output_file_prefix="Pos_Ctrl")

# Run the IDH1 Positive Control
# load the feature data for LGG (gene expression data)
feature_data<-data.frame(fread(paste(INPUT,"pancan_RNAexp_LGG",sep=""),stringsAsFactors = FALSE),row.names = 1)
tissue<-"LGG"
LURE(bait_gene="IDH1-SET1_missense",
       gmt_file="positive_control_IDH1_missense.gmt",
       feature_data=feature_data,
       num_permutations=5,
       max_num_events=5,
       percent_overlap=0,
       LURE_pvalue_threshold=.05,
       min_gene_set_size=5,
       gsea_pvalue_threshold=.05,
       gsea_fdr_threshold=.25,
       max_tree_length = 5,
       folds=10,
       enrichment_analysis_only=FALSE,
       output_file_prefix="Pos_Ctrl")



# For creating the gmt files. These are already provided so this code is not necessary to run
# gene<-"SF3B1_missense"
# orig_gmt_file<-"mutation_specific_cytoband_fusion_5_15_2019.gmt"
# # divide up the SF3B1 UVM mutants
# setwd(DD_HOME)
# system(paste("grep SF3B1",orig_gmt_file,">","new_temp.gmt"))
#
# # change these variables to "mutant_samples"
# SF3B1_mutants<-as.character(t(read.table(paste0(DD_HOME,"SF3B1_mutants"))[,-c(1,2)]))
# SF3B1_mutants<-as.character(SF3B1_mutants[which(SF3B1_mutants %in% rownames(feature_data))])
# SF3B1_mutant_set1<-sample(SF3B1_mutants,8)
# SF3B1_mutant_set2<-sample(SF3B1_mutants[-(which(SF3B1_mutants %in% SF3B1_mutant_set1))],5)
# SF3B1_mutant_set3<-sample(SF3B1_mutants[-(which(SF3B1_mutants %in% c(SF3B1_mutant_set1,SF3B1_mutant_set2)))],5)
#
# # create gmt file
# output_file<-data.frame(matrix(ncol=length(SF3B1_mutant_set1)+2,nrow=0))
# output_file<-rbind(output_file,data.frame(matrix(c("SF3B1-SET1_missense","N/A",SF3B1_mutant_set1),nrow=1)))
# output_file<-rbind.fill(output_file,data.frame(matrix(c("SF3B1-SET2_missense","N/A",SF3B1_mutant_set2),nrow=1)))
# output_file<-rbind.fill(output_file,data.frame(matrix(c("SF3B1-SET3_missense","N/A",SF3B1_mutant_set3),nrow=1)))
# write.table(output_file,file=paste0(DD_HOME,"mutants",gene,".gmt"),sep="\t",quote=FALSE,na="",row.names=FALSE,col.names = FALSE)
#
# gmt_file<-"positive_control_SF3B1.gmt"
#
#
# # create mutant set file for IDH1
# gene<-"IDH1_missense"
# orig_gmt_file<-"mutation_specific_cytoband_fusion_5_15_2019.gmt"
# # divide up the IDH1 LGG mutants
# system(paste("grep IDH1_missense",orig_gmt_file,">","new_temp.gmt"))
#
#
#
# feature_data<-data.frame(fread(paste(PANCAN_DATA,"pancan_RNAexp_LGG",sep=""),stringsAsFactors = FALSE),row.names = 1)
#
#
# mutant_samples<-as.character(t(read.table(paste0(DD_HOME,"new_temp.gmt"))[,-c(1,2)]))
# mutant_samples<-as.character(mutant_samples[which(mutant_samples %in% rownames(feature_data))])
# length(unique(mutant_samples))
# mutant_samples<-unique(mutant_samples)
# mutant_sample_set1<-sample(mutant_samples,317)
# mutant_sample_set2<-sample(mutant_samples[-(which(mutant_samples %in% mutant_sample_set1))],20)
# mutant_sample_set3<-sample(mutant_samples[-(which(mutant_samples %in% c(mutant_sample_set1,mutant_sample_set2)))],20)
# mutant_sample_set4<-sample(mutant_samples[-(which(mutant_samples %in% c(mutant_sample_set1,mutant_sample_set2,mutant_sample_set3 )))],20)
#
# # create gmt file
# output_file<-data.frame(matrix(ncol=length(mutant_sample_set1)+2,nrow=0))
# output_file<-rbind(output_file,data.frame(matrix(c("IDH1-SET1_missense","N/A",mutant_sample_set1),nrow=1)))
# output_file<-rbind.fill(output_file,data.frame(matrix(c("IDH1-SET2_missense","N/A",mutant_sample_set2),nrow=1)))
# output_file<-rbind.fill(output_file,data.frame(matrix(c("IDH1-SET3_missense","N/A",mutant_sample_set3),nrow=1)))
# output_file<-rbind.fill(output_file,data.frame(matrix(c("IDH1-SET4_missense","N/A",mutant_sample_set4),nrow=1)))
#
# write.table(output_file,file=paste0(DD_HOME,"mutants",gene,".gmt"),sep="\t",quote=FALSE,na="",row.names=FALSE,col.names = FALSE)
#
# setwd(DD_HOME)
# system(paste0("cat mutants",gene,".gmt ",orig_gmt_file," > positive_control_",gene,".gmt"))
#
#
