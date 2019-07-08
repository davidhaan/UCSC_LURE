
run_timestamp<-strftime(Sys.time(),"%Y_%m_%d_%H_%M")

if ((Sys.info()["nodename"] == "Davids-MacBook-Pro.local") | 
    (grepl("eduroam",Sys.info()["nodename"])==1) |
    (grepl("dhcp",Sys.info()["nodename"])==1)) {
  print("Running Locally...")
  PCAWG_DATA<<-"~/TCGA_data/pancan/pcawg_data/"
  PANCAN_DATA<<-"~/TCGA_data/pancan/pancan_data/"
  OUTPUT_DATA<<-"~/TCGA_data/pancan/output/"
  DD_HOME<<-"~/TCGA_data/"
  SCRIPTS<<-"~/GoogleDriveUCSC/DARKMATTER/"
  # LOCAL SETTINGS
  max_folds<-20
  # set to 20 for production run
  num_permutations<-5
  min_mutants<-5
  min_score<-.75
  max_num_events<-15
  folds<-10
  # the flag is either set to expression or a splicing phenotype to analyze
  tissue<-"UVM"
  #percent_overlap<-.33
  percent_overlap<-.5
  
  dumpsterdiverflag<-"N"
  source(paste(sep="",SCRIPTS,"dd_functions.R"))
  registerDoMC(detectCores()/4)
  
  
} else {
  print("Running on Server")
  PCAWG_DATA<<-Sys.getenv("PCAWG_DATA")
  PANCAN_DATA<<-Sys.getenv("PANCAN_DATA")
  OUTPUT_DATA<<-Sys.getenv("OUTPUT_DATA")
  DD_HOME<<-Sys.getenv("DD_HOME")
  SCRIPTS<<-Sys.getenv("SCRIPTS")
  max_folds<-20
  folds<-10
  # set to 50 for production run
  num_permutations<-5
  # used to find the initial baits, only genes with 10 or more mutated samples are considered
  min_mutants<-5
  # used to find the initial baits, only mutations with this score or better are considered
  min_initial_pr_auc_score<-.5
  min_initial_recall_scores<-.7
  # used to limit the number of "catches" found by GSEA and considered for FDR analysis.  The events are sorted by GSEA
  # score so the top events will be chosen. 
  max_num_events<-5
  # if  <percent_overlap> of the new event samples are in the existing positive set already then we skip it.  
  # A smaller number is more restrictive, for telomere analysis I do .5, for the larger TCGA analysis I do .3
  percent_overlap<-.3
  # parameter for GSEA, only events with n or more mutated samples are considered
  min_gene_set_size<-5
  args = commandArgs(trailingOnly=TRUE)
  tissue<<-args[1]
  dumpsterdiverflag<<-args[2]
  source(paste(sep="",SCRIPTS,"dd_functions.R"))
  registerDoMC(detectCores()/2)
  
  
  
}







# for SF3B1 positive control use UVM
# for IDH1 use LGG

tissue<-"UVM"

if (is.na(tissue)) {
  print("Must supply a tissue")
  q()
}


print(paste0("Loading Data and Analyzing...",tissue))

feature_data<-data.frame(fread(paste(PANCAN_DATA,"pancan_RNAexp_",tissue,sep=""),stringsAsFactors = FALSE),row.names = 1)


gene<-"SF3B1_missense"
orig_gmt_file<-"mutation_specific_cytoband_fusion_5_15_2019.gmt"
# divide up the SF3B1 UVM mutants
setwd(DD_HOME)
system(paste("grep SF3B1",orig_gmt_file,">","new_temp.gmt"))

# change these variables to "mutant_samples"
SF3B1_mutants<-as.character(t(read.table(paste0(DD_HOME,"SF3B1_mutants"))[,-c(1,2)]))
SF3B1_mutants<-as.character(SF3B1_mutants[which(SF3B1_mutants %in% rownames(feature_data))])
SF3B1_mutant_set1<-sample(SF3B1_mutants,8)
SF3B1_mutant_set2<-sample(SF3B1_mutants[-(which(SF3B1_mutants %in% SF3B1_mutant_set1))],5)
SF3B1_mutant_set3<-sample(SF3B1_mutants[-(which(SF3B1_mutants %in% c(SF3B1_mutant_set1,SF3B1_mutant_set2)))],5)

# create gmt file
output_file<-data.frame(matrix(ncol=length(SF3B1_mutant_set1)+2,nrow=0)) 
output_file<-rbind(output_file,data.frame(matrix(c("SF3B1-SET1_missense","N/A",SF3B1_mutant_set1),nrow=1)))
output_file<-rbind.fill(output_file,data.frame(matrix(c("SF3B1-SET2_missense","N/A",SF3B1_mutant_set2),nrow=1)))
output_file<-rbind.fill(output_file,data.frame(matrix(c("SF3B1-SET3_missense","N/A",SF3B1_mutant_set3),nrow=1)))
write.table(output_file,file=paste0(DD_HOME,"mutants",gene,".gmt"),sep="\t",quote=FALSE,na="",row.names=FALSE,col.names = FALSE)


# create gmt file
output_file<-data.frame(matrix(ncol=length(SF3B1_mutants)+3)) 
output_file<-rbind(output_file,data.frame(matrix(c(paste0(gene,"_set1"),"N/A",mutant_sample_set1),nrow=1)))
output_file<-rbind.fill(output_file,data.frame(matrix(c(paste0(gene,"_set2"),"N/A",mutant_sample_set2),nrow=1)))
output_file<-rbind.fill(output_file,data.frame(matrix(c(paste0(gene,"_set3"),"N/A",mutant_sample_set3),nrow=1)))
output_file<-rbind.fill(output_file,data.frame(matrix(c(paste0(gene,"_set4"),"N/A",mutant_sample_set4),nrow=1)))


gmt_file<-"positive_control_SF3B1.gmt"



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












# create mutant set file for IDH1
gene<-"IDH1_missense"
orig_gmt_file<-"mutation_specific_cytoband_fusion_5_15_2019.gmt"
# divide up the IDH1 LGG mutants
system(paste("grep IDH1_missense",orig_gmt_file,">","new_temp.gmt"))



feature_data<-data.frame(fread(paste(PANCAN_DATA,"pancan_RNAexp_LGG",sep=""),stringsAsFactors = FALSE),row.names = 1)


mutant_samples<-as.character(t(read.table(paste0(DD_HOME,"new_temp.gmt"))[,-c(1,2)]))
mutant_samples<-as.character(mutant_samples[which(mutant_samples %in% rownames(feature_data))])
length(unique(mutant_samples))
mutant_samples<-unique(mutant_samples)
mutant_sample_set1<-sample(mutant_samples,317)
mutant_sample_set2<-sample(mutant_samples[-(which(mutant_samples %in% mutant_sample_set1))],20)
mutant_sample_set3<-sample(mutant_samples[-(which(mutant_samples %in% c(mutant_sample_set1,mutant_sample_set2)))],20)
mutant_sample_set4<-sample(mutant_samples[-(which(mutant_samples %in% c(mutant_sample_set1,mutant_sample_set2,mutant_sample_set3 )))],20)

# create gmt file
output_file<-data.frame(matrix(ncol=length(mutant_sample_set1)+2,nrow=0)) 
output_file<-rbind(output_file,data.frame(matrix(c("IDH1-SET1_missense","N/A",mutant_sample_set1),nrow=1)))
output_file<-rbind.fill(output_file,data.frame(matrix(c("IDH1-SET2_missense","N/A",mutant_sample_set2),nrow=1)))
output_file<-rbind.fill(output_file,data.frame(matrix(c("IDH1-SET3_missense","N/A",mutant_sample_set3),nrow=1)))
output_file<-rbind.fill(output_file,data.frame(matrix(c("IDH1-SET4_missense","N/A",mutant_sample_set4),nrow=1)))

write.table(output_file,file=paste0(DD_HOME,"mutants",gene,".gmt"),sep="\t",quote=FALSE,na="",row.names=FALSE,col.names = FALSE)

setwd(DD_HOME)
system(paste0("cat mutants",gene,".gmt ",orig_gmt_file," > positive_control_",gene,".gmt"))




gene<-"IDH1-SET1_missense"

  print(paste0("Initial Processing of ",gene))
  

  #gmt_file<-"positive_control_SF3B1.gmt"

  print(paste0("gmt_file:",gmt_file))
  

  
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
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
