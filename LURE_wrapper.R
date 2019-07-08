# LURE WRAPPER
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
  source(paste(sep="",SCRIPTS,"dd_functions.R"))
  registerDoMC(detectCores()/4)
  
  
} else {
  print("Running on Server")
  PCAWG_DATA<<-Sys.getenv("PCAWG_DATA")
  PANCAN_DATA<<-Sys.getenv("PANCAN_DATA")
  OUTPUT_DATA<<-Sys.getenv("OUTPUT_DATA")
  DD_HOME<<-Sys.getenv("DD_HOME")
  SCRIPTS<<-Sys.getenv("SCRIPTS")
  source(paste(sep="",SCRIPTS,"dd_functions.R"))
  registerDoMC(detectCores())
  }

option_list = list(
  make_option(c("--folds"), type="numeric", default=10, 
              help="Number of Cross Validation Folds", metavar="character"),
  
  make_option(c("--num_permutations"), type="numeric", default=25, 
              help="Number of Permutations/Iterations", metavar="character"),
  
  make_option(c("--min_gene_set_size"), type="numeric", default=5, 
              help="prey event minimum size: parameter for GSEA, only events with n or more mutated samples are considered", metavar="character"),
  
  make_option(c("--percent_overlap"), type="numeric", default=.3, 
              help="If <percent_overlap> of the new event samples are in the existing positive set already then we skip it. A smaller number is more restrictive, for telomere analysis I do .5, for the larger TCGA analysis I do .3",
              metavar="character"),
  
  make_option(c("--max_tree_length"), type="numeric", default=3, 
              help="Max Tree Length: Here we set the max length of the Event Discovery Tree", metavar="character"),

  make_option(c("--bait_gene"), type="character", default="custom_to_test.tsv", 
              help=".tsv file containing bait genes to test", metavar="character"),
   
  make_option(c("--gmt_file"), type="character", default="mutation_specific_cytoband_fusion_5_15_2019.gmt", 
              help=".gmt file containing seed data", metavar="character"),
  
  make_option(c("--gsea_fdr_threshold"), type="numeric", default=.25, 
              help="FDR value threshold for GSEA step, .25 for PANCAN, .5 for telomere", metavar="character"),
  
  make_option(c("--gsea_pvalue_threshold"), type="numeric", default=.05, 
              help="P value threshold for GSEA step", metavar="character"),
  
  make_option(c("--LURE_pvalue_threshold"), type="numeric", default=.05, 
              help="P value threshold for LURE AUC score step", metavar="character"),
  
  make_option(c("--max_num_events"), type="numeric", default=3, 
              help="Used to limit the number of prey events found by GSEA and considered for LURE's classifier AUC score step.  The events are sorted by GSEA NES score so the top events will be chosen. The larger this parameter the longer the runtime.",
              metavar="character"),
  make_option(c("--feature_data_file"), type="character", default="both", 
              help="Feature Data file", metavar="character"),
  
  make_option(c("--target_gmt_file"), type="character", default="LUAD_functional_coding_non_coding_amp_del_fusion.gmt", 
              help="This argument only pertains when LURE is run with enrichment only.  It is the gmt file for the test/target dataset.  The original gmt_file argument is for the bait", metavar="character"),
  
  make_option(c("--target_feature_file"), type="character", default="mutation_specific_cytoband_fusion_5_15_2019.gmt", 
              help="This argument only pertains when LURE is run with enrichment only.  It is the gmt file for the test/target dataset.  The original gmt_file argument is for the bait", metavar="character"),
  
  make_option(c("--output_file_prefix"), type="character", default="V50",
              help="This is the file prefix assigned to all the output files.  For multiple runs it helps keep track of each run.")
  
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
folds<<-opt$folds
print(opt)


feature_data<-data.frame(fread(paste(PANCAN_DATA,opt$feature_data_file,sep=""),stringsAsFactors = FALSE),row.names = 1)

tissue<-"NA"


LURE(bait_gene=opt$bait_gene,
     gmt_file=opt$gmt_file,
     feature_data=feature_data,
     num_permutations=opt$num_permutations,
     max_num_events=opt$max_num_events,
     percent_overlap=opt$percent_overlap,
     LURE_pvalue_threshold=opt$LURE_pvalue_threshold,
     min_gene_set_size=opt$min_gene_set_size,
     gsea_pvalue_threshold=opt$gsea_pvalue_threshold,
     gsea_fdr_threshold=opt$gsea_fdr_threshold,
     max_tree_length =opt$max_tree_length,
     folds=opt$folds,
     enrichment_analysis_only=FALSE,
     output_file_prefix=opt$output_file_prefix)



