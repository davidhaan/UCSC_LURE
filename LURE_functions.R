## LURE_functions.R
## David Haan 2015-2019
## Functions associated with LURE




list.of.packages <- c("glmnet", "plyr","matrixStats","data.table","methods","doMC","tools","dplyr","PRROC","ggplot2","optparse","RcppGreedySetCover","igraph","survminer","survival")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  print("Packages not found:")
  print(new.packages)
  print("Installing now")
  install.packages(new.packages,repos = "http://cran.us.r-project.org")
}

lapply(list.of.packages, require, character.only = TRUE)


#### --- LURE FUNCTIONS
LURE<-function(bait_gene,
               gmt_file,
               feature_data,
               num_permutations,
               max_num_events,
               percent_overlap,
               LURE_pvalue_threshold,
               min_gene_set_size,
               gsea_pvalue_threshold,
               gsea_fdr_threshold,
               max_tree_length,
               folds,
               enrichment_analysis_only=TRUE,
               target_gmt_file,
               target_feature_data,
               output_file_prefix) {
  #set defaults if none given
  if(is.na(bait_gene)) {
    print("Must provide a bait gene")
    return(1)
  }
  if(is.na(gmt_file)) {
    print("Must provide a gmt file location")
    return(1)
  }
  if(!(exists("feature_data"))) {
    print("Must provide feature data in data.frame")
    return(1)
  }
  if(enrichment_analysis_only & !(exists("target_feature_data"))) {
    print("Enrichment Analysis is set to TRUE and no TARGET feature file is provided!")
    return(1)
  }
  if(enrichment_analysis_only & !(exists("target_gmt_file"))) {
    print("Enrichment Analysis is set to TRUE and no TARGET feature file is provided!")
    return(1)
  }
  
  
  if(is.na(num_permutations))
    num_permutations<-25
  if(is.na(max_num_events))
    max_num_events<-3
  if(is.na(percent_overlap))
    percent_overlap<-.25
  if(is.na(LURE_pvalue_threshold))
    LURE_pvalue_threshold<-.05
  if(is.na(min_gene_set_size))
    min_gene_set_size<-4
  if(is.na(gsea_pvalue_threshold))
    gsea_pvalue_threshold<-.05
  if(is.na(gsea_fdr_threshold))
    gsea_fdr_threshold<-.25
  if(is.na(max_tree_length))
    max_tree_length<-3
  if(is.na(folds))
    folds<-10
  if(is.na(enrichment_analysis_only))
    enrichment_analysis_only<-FALSE
  
  # load gmt file
  no_col <- max(count.fields(paste0(DD_HOME,gmt_file), sep = "\t"))
  gmt_file_data <- read.table(paste0(DD_HOME,gmt_file),sep="\t",fill=TRUE,header = F,col.names=1:no_col,stringsAsFactors = FALSE)
  
  print(paste0("Initial Processing of ",bait_gene))
  
  print(paste0("gmt_file:",gmt_file))
  
  
  # Here we make our feature set X and label set Y
  # multiple bait genes are allowed, they are separated by a semicolon ";"
  # the function create_XY only works with one alteration at a time, so we must loop through
  positive_set<-c()
  num_samples<-c()
  bait_gene_list<-strsplit(bait_gene,";")[[1]]
  for (gene in bait_gene_list) {
    print(gene)
    gene_prefix<-unlist(strsplit(gene,"_"))[c(TRUE,FALSE)]
    gene_mut_type<-unlist(strsplit(gene,"_"))[c(FALSE,TRUE)]
    print(paste0("Selecting Mutation Type: ",gene_mut_type))
    temp_set_XY<-createXY_gmt(gene, gmt_file_data, feature_data, NA)
    if (length(rownames(temp_set_XY$Y)[which(temp_set_XY$Y==1)]) == 0)
      print("WARNING: Bait Gene Alterations not found in provided gmt file")
    num_samples<-c(num_samples,length(rownames(temp_set_XY$Y)[which(temp_set_XY$Y==1)]))
    positive_set<-c(positive_set,rownames(temp_set_XY$Y)[which(temp_set_XY$Y==1)])
  }
  # set full_set_XY to our final X and Y to use for our classification model
  full_set_XY<-temp_set_XY
  full_set_XY$Y$sample<-0
  full_set_XY$Y$sample[which(rownames(full_set_XY$Y) %in% positive_set)]<-1
  
  # check to see if there are at least 5 positive samples
  if (length(full_set_XY$Y[full_set_XY$Y==1]) < 5) {
    print("At least 5 mutant samples are needed to establish a 'bait'")
    return()
  } 
  
  # Run logistic regression classifier
  initial_bait_model<-runClassifier_V2(full_set_XY$X,full_set_XY$Y,folds,num_permutations,alpha=1,weighted = TRUE,stratified = TRUE)
  
  print(paste("ORIGINAL PR AUC:",mean(initial_bait_model$pr_auc)))
  print(paste("ORIGINAL PRECISION:",mean(initial_bait_model$precision)))
  print(paste("ORIGINAL RECALL:",mean(initial_bait_model$recall)))
  print(paste("ORIGINAL AUC LAMBDA:",log(mean(initial_bait_model$min_lambda))))
  print(paste("ORIGINAL PRECISION LAMBDA:",log(mean(initial_bait_model$min_precision_lambda ))))
  print(paste("ORIGINAL RECALL LAMBDA:",log(mean(initial_bait_model$min_recall_lambda ))))
  

  # make initial predictions on expected positive test data
  if (enrichment_analysis_only==TRUE) {
    pred_resp<-((predict(initial_bait_model$model, s=mean(initial_bait_model$min_lambda), newx=data.matrix(target_feature_data), type="response")))[,1]
  } else {
    pred_resp<-((predict(initial_bait_model$model, s=mean(initial_bait_model$min_lambda), newx=data.matrix(full_set_XY$test), type="response")))[,1]
  }
  pred_resp<-sort(pred_resp,decreasing = TRUE)
  
  # make a pdf plot of the initial predictions
  scores_labels<-data.frame("Sample_Names"=names(pred_resp), "Classifier_Score"=pred_resp, "Color"="black",stringsAsFactors = FALSE)
  scores_labels$Color[(scores_labels$Classifier_Score > .5)]<-"red"

  scores_labels$Sample_Names<-factor(scores_labels$Sample_Names, levels=scores_labels$Sample_Names[order(scores_labels$Classifier_Score, decreasing=TRUE)])
  positive_samples<-length(scores_labels$Classifier_Score[scores_labels$Classifier_Score>.5])
  pdf(paste0(PANCAN_DATA,output_file_prefix,"_Initial_Classifier_Scores_",positive_samples,"_positive_",round(mean(initial_bait_model$pr_auc),2),"_AUC_",tissue,"_",bait_gene,".pdf"),onefile = FALSE,width=8,height = 5)
  print(ggplot(data=scores_labels, aes(x=Sample_Names,y=Classifier_Score,fill=Color)) +
          geom_bar(stat="identity", width=0.5) +
          scale_fill_manual("", values = c("black" = "black", "red" = "red")) +
          xlab("Samples") +
          ylab("Predicted Mutation Status") +
          theme_bw() +
          geom_hline(yintercept = .5) +
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
          theme(legend.position = "none") +
          theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=.3)) +
          labs(title = "Initial LURE Sample Scores", subtitle = paste0("Bait Gene: ",bait_gene,"; CV Test PR AUC: ",round(mean(initial_bait_model$pr_auc),2))))
  dev.off()
  
  #create the global classifier tree variable and start the LURE steps by recursively run through the tree paths
  tree_path<<-bait_gene
  tree_path_list<<-c()
  STG_output<<-data.frame(matrix(nrow=0,ncol=5))

  # If enrichment analysis only is TRUE, then only one iteration of LURE is executed.  Essentially the initial_bait_classifier is applied to 
  # the user provided test set
  if (enrichment_analysis_only==TRUE) {
    results<-positive_recursive_iterate(bait=paste0("Initial_",bait_gene), 
                                        pred_resp, 
                                        tree_path, 
                                        target_gmt_file, 
                                        oncogene_neg_expression_data, 
                                        "origY"=data.frame(matrix(0,ncol=1,nrow=nrow(oncogene_neg_expression_data))), 
                                        orig_score_vector=initial_bait_model$pr_auc, 
                                        num_permutations, 
                                        max_num_events,
                                        percent_overlap, 
                                        pvalue_threshold=LURE_pvalue_threshold, 
                                        min_gene_set_size=min_gene_set_size, 
                                        gsea_pvalue_threshold=gsea_pvalue_threshold, 
                                        gsea_fdr_threshold=gsea_fdr_threshold, 
                                        tree_length = 3, 
                                        max_tree_length = max_tree_length, 
                                        folds=folds, 
                                        enrichment_analysis_only=TRUE)
  } else { 
    results<-positive_recursive_iterate(bait=paste0("Initial_",bait_gene),
                                        pred_resp, 
                                        tree_path,
                                        gmt_file,
                                        full_set_XY$X, 
                                        "origY"=full_set_XY$Y, 
                                        orig_score_vector=initial_bait_model$pr_auc, 
                                        num_permutations, 
                                        max_num_events,
                                        percent_overlap, 
                                        pvalue_threshold=LURE_pvalue_threshold, 
                                        min_gene_set_size=min_gene_set_size,
                                        gsea_pvalue_threshold=gsea_pvalue_threshold, 
                                        gsea_fdr_threshold=gsea_fdr_threshold, 
                                        tree_length = 1, 
                                        max_tree_length = max_tree_length, 
                                        folds=folds, 
                                        enrichment_analysis_only=FALSE)
  }
  # LURE iterative analysis is done, now processing the results
  
  # IF ENRICHMENT ONLY ANALYSIS
  if (enrichment_analysis_only==TRUE) {
    filename_gene<-gsub("\\(","_",bait_gene)
    filename_gene<-gsub("\\)","",filename_gene)
    gmt_prefix<-strsplit(gmt_file,"\\.")[[1]][1]
    filename1<-paste0(PANCAN_DATA,output_file_prefix,"_TCGA_",tissue,"_Enrich_Only_Analysis_Output_",run_timestamp,"_",gmt_prefix,"_",filename_gene,".tsv")
    filename4<-paste0(PANCAN_DATA,output_file_prefix,"_TCGA_",tissue,"_Enrich_Only_classifier_scores_",run_timestamp,"_",gmt_prefix,"_",filename_gene,".tsv")
    filename5<-paste0(PANCAN_DATA,output_file_prefix,"_TCGA_",tissue,"_Enrich_Only_mutation_matrix_",run_timestamp,"_",gmt_prefix,"_",filename_gene,".tsv")
    filename6<-paste0(PANCAN_DATA,output_file_prefix,"_TCGA_",tissue,"_Enrich_Only_Oncoprint_Plot_",run_timestamp,"_",gmt_prefix,"_",filename_gene,".pdf")
    filename7<-paste0(PANCAN_DATA,output_file_prefix,"_TCGA_",tissue,"_Enrich_Only_Original_PR_AUC_scores",run_timestamp,"_",gmt_prefix,"_",filename_gene,".tsv")
    filename8<-paste0(PANCAN_DATA,output_file_prefix,"_TCGA_",tissue,"_Initial_Classifier_Coefficients",run_timestamp,"_",gmt_prefix,"_",filename_gene,".tsv")
    
    print(paste0("Routine Complete for ",bait_gene,"_","____",gmt_prefix,"...writing classifier score file: "))
    write.csv(pred_resp,file=filename4)
    write.csv(initial_bait_model$pr_auc,file=filename7)
    write.csv(coef(initial_bait_model, s = mean(initial_bait_model$min_lambda)),file=filename8)
    
    # Create output files for Yanni's oncoprint plot Python function
    #load gmt_file
    if (nrow(results)>0) {
      write.table(results,file=filename1,quote=FALSE,sep="\t",col.names = FALSE,row.names = FALSE)
      
      no_col <- max(count.fields(paste0(DD_HOME,gmt_file), sep = "\t"))
      gmt_file_data <- read.table(paste0(DD_HOME,gmt_file),sep="\t",fill=TRUE,header = F,col.names=1:no_col,stringsAsFactors = FALSE)
      
      gene_list<-c(bait_gene,rownames(results))
      oncoprint_mutation_data<-data.frame(matrix(NA,ncol=nrow(oncogene_neg_expression_data),nrow=length(gene_list)))
      colnames(oncoprint_mutation_data)<-names(pred_resp)
      output_gene_list<-c()
      for (raw_gene in gene_list) {
        print(raw_gene)
        mutant_samples<-gmt_file_data[which(gmt_file_data$X1==raw_gene),]
        mutant_samples<-mutant_samples[,-c(1,2)]
        print(length(mutant_samples))
        split_list<-strsplit(raw_gene,"_")
        gene_name<-split_list[[1]][1]
        mutation<-paste(split_list[[1]][2:length(split_list[[1]])],collapse="_")
        output_gene_list<-c(output_gene_list,gene_name)
        oncoprint_mutation_data[which(gene_list==raw_gene),which(colnames(oncoprint_mutation_data) %in% as.character(mutant_samples))]<-mutation
      }
      rownames(oncoprint_mutation_data)<-make.names(output_gene_list,unique = TRUE)
      write.csv(oncoprint_mutation_data,file=filename5)
      # print oncoprint using yanni's oncoprint plot:
      system(paste0("python3 ",SCRIPTS,"oncoprint.py -memo -c ",filename4," -m ",filename5," -o ",filename6))
      
    }
    print("No results found in the samples which recieved high scores")
    
    
  } # end of ENRICHMENT_ONLY LOOP
  
  #  
  # ENRICHMENT ONLY IS FALSE, Running full LURE and set coverage 
  #
  else {
    
    filename_gene<-gsub("\\(","_",bait_gene)
    filename_gene<-gsub("\\)","",filename_gene)
    gmt_prefix<-strsplit(gmt_file,"\\.")[[1]][1]
    filename1<-paste0(PANCAN_DATA,output_file_prefix,"_TCGA_",tissue,"_",filename_gene,"_tree_path_list_",run_timestamp,"_",gmt_prefix,".tsv")
    filename2<-paste0(PANCAN_DATA,output_file_prefix,"_TCGA_",tissue,"_",filename_gene,"_STG_data_",run_timestamp,"_",gmt_prefix,".tsv")
    filename3<-paste0(PANCAN_DATA,output_file_prefix,"_TCGA_",tissue,"_",filename_gene,"_set_cover_solution_",run_timestamp,"_",gmt_prefix,".tsv")
    filename4<-paste0(PANCAN_DATA,output_file_prefix,"_TCGA_",tissue,"_",filename_gene,"_set_cover_classifier_scores_",run_timestamp,"_",gmt_prefix,".tsv")
    filename5<-paste0(PANCAN_DATA,output_file_prefix,"_TCGA_",tissue,"_",filename_gene,"_Set_Cover_mutation_matrix_",run_timestamp,"_",gmt_prefix,".tsv")
    filename6<-paste0(PANCAN_DATA,output_file_prefix,"_TCGA_",tissue,"_",filename_gene,"_Oncoprint_Plot_",run_timestamp,"_",gmt_prefix,".pdf")
    filename7<-paste0(PANCAN_DATA,output_file_prefix,"_TCGA_",tissue,"_",filename_gene,"_Original_PR_AUC_scores",run_timestamp,"_",gmt_prefix,".tsv")
    filename8<-paste0(PANCAN_DATA,output_file_prefix,"_TCGA_",tissue,"_",filename_gene,"_set_cover_PR_AUC_scores_",run_timestamp,"_",gmt_prefix,".tsv")
    
    
    # Run the Set Coverage Algorithm  
    gene_list<-unique(c(bait_gene,as.character(STG_output$target)))
    if (length(gene_list) < 2) {
      print("No Catch Events Found...")
      tree_path_list<-ldply(tree_path_list, rbind)
      write.table(tree_path_list,file=filename1,quote=FALSE,sep="\t",col.names = FALSE,row.names = FALSE)
      return()
    }

    # create classifier using all of the new catches and original bait
    positive_set<-c()
    gene<-gene_list[1]
    for (gene in gene_list) {
      #print(gene)
      gene_prefix<-unlist(strsplit(gene,"_"))[c(TRUE,FALSE)]
      gene_mut_type<-unlist(strsplit(gene,"_"))[c(FALSE,TRUE)]
      print(paste0("Selecting Mutation Type: ",gene_mut_type))
      #full_set_XY<-createXY_pos(gene_prefix,mutation_data,expression_data,NA,mutation_type=gene_mut_type, cnv_data=cnv_data, fusion_data = fusion_data, pcawg_mutation_data=pcawg_mutation_data)
      full_set_XY<-createXY_gmt(gene, gmt_file_data, feature_data, NA)
      
      positive_set<-c(positive_set,rownames(full_set_XY$Y)[which(full_set_XY$Y==1)])
    }
    final_set_XY<-full_set_XY
    final_set_XY$Y$sample<-0
    final_set_XY$Y$sample[which(rownames(final_set_XY$Y) %in% positive_set)]<-1
    
    # run logistic regression classifier
    initial_model<-runClassifier_V2(final_set_XY$X,final_set_XY$Y,nfolds_input=folds,iterations=num_permutations,alpha=1,weighted = TRUE,stratified = TRUE)
    
    print(paste("NEW CATCH PR AUC:",mean(initial_model$pr_auc)))
    print(paste("NEW CATCH PRECISION:",mean(initial_model$precision)))
    print(paste("NEW CATCH RECALL:",mean(initial_model$recall)))
    
    # Score the samples using new classifier
    classifier_scores<-((predict(initial_model$model, s=min((initial_model$min_lambda)[which(initial_model$pr_auc==max(initial_model$pr_auc))]), newx=data.matrix(final_set_XY$X), type="response")))[,1]
    classifier_scores<-classifier_scores[order(classifier_scores,decreasing = TRUE)]

    # shorten bait gene size, just for filename output    
    bait_gene_for_filename<-substr(bait_gene,1,200)

    # load GMT data file
    no_col <- max(count.fields(paste0(DD_HOME,gmt_file), sep = "\t"))
    gmt_file_data <- read.table(paste0(DD_HOME,gmt_file),sep="\t",fill=TRUE,header = F,col.names=1:no_col,stringsAsFactors = FALSE)
    
    # create mutation matrix
    oncoprint_mutation_data<-data.frame(matrix(NA,ncol=length(classifier_scores),nrow=length(gene_list)))
    colnames(oncoprint_mutation_data)<-names(classifier_scores)
    gene_list_for_df<-c()
    row<-1
    for (row in 1:length(gene_list)) {
      print(row)
      mutant_samples<-gmt_file_data[which(toupper(gmt_file_data$X1)==toupper(gene_list[row])),]
      mutant_samples<-mutant_samples[,-c(1,2)]
      split_list<-strsplit(gene_list[row],"_")
      print(toupper(gene_list[row]))
      gene<-split_list[[1]][1]
      mutation<-toupper(paste(split_list[[1]][2:length(split_list[[1]])],collapse="_"))
      gene_list_for_df<-c(gene_list_for_df,toupper(gene_list[row]))
      oncoprint_mutation_data[row,which(colnames(oncoprint_mutation_data) %in% as.character(mutant_samples))]<-mutation
      
    }
    rownames(oncoprint_mutation_data)<-gene_list_for_df
    
    
    # print bipartite graph
    # Create a adjacency matrix for the bipartite graph and set coverage algorithm
    # Remove samples with classifier scores lower than .5
    adj_matrix<-data.frame(oncoprint_mutation_data[,which(classifier_scores > .5)],stringsAsFactors = FALSE)
    adj_matrix[adj_matrix!=""]<-1
    adj_matrix[adj_matrix==""]<-0
    adj_matrix[is.na(adj_matrix)]<-0
    adj_matrix<-data.matrix(adj_matrix)
    
    # remove the samples which do not have any mutation in the genes in question
    adj_matrix<-adj_matrix[,(which(colSums(adj_matrix)>0))]
    
    # create a bipartite data structure
    bipartite_data<-(matrix(ncol=2,nrow=0))
    for (col in 1:ncol(adj_matrix)) {
      for (row in 1:nrow(adj_matrix)) {
        if (adj_matrix[row,col]==1)
          bipartite_data<-rbind(bipartite_data,c(rownames(adj_matrix)[row],colnames(adj_matrix)[col]))
        
      }
    }
    bipartite_data<-data.frame(bipartite_data)
    
    # run a greedy(quick) Set Coverage Algorithm on the remaining samples
    df_output<-greedySetCover(bipartite_data, data.table = FALSE)
    set_coverage_solution<-toupper(gene_list[(toupper(gene_list) %in% unique(as.character(df_output$X1)))])
    not_in_set_coverage_solution<-unique(bipartite_data$X1[!(bipartite_data$X1 %in% set_coverage_solution)])
    print(paste0("Set Coverage Solution:",set_coverage_solution))
    
    # check to see if bait gene is in the solution
    if (!(toupper(bait_gene) %in% set_coverage_solution)) {
      print("The bait gene is not in the set coverage solution, skipping this result...")
      return()
    }
    # make bipartite graph with shaded samples 
    g<-graph_from_data_frame(data.frame(bipartite_data), directed=FALSE)
    
    V(g)$type<-bipartite_mapping(g)$type
    V(g)$shape<-ifelse(V(g)$type, "circle", "square")
    V(g)$label.dist<-2
    V(g)$label.degree<-ifelse(V(g)$type, pi/2, 3*pi/2)
    # node/vertex color
    V(g)$color <-ifelse(V(g)$type, "#E69F00","#E69F00")
    # set the vertex and edge color to red for those not in the covered set 
    V(g)$color[which(vertex_attr(g)$name %in% not_in_set_coverage_solution)] <- "red"
    E(g)$color<-"grey"
    E(g)$color[which(!(bipartite_data$X1 %in% set_coverage_solution))]<-"red"
    # remove (change color to white) for the sample names
    V(g)$label.color<-ifelse(V(g)$type, "#FFFFFF","black")
    
    # space out the mutation vertices so the label names don't overlap
    coords<-layout_as_bipartite(g)
    num_mutations<-length(which(coords[,2]==1))
    max_length<-(max(coords[,1]))
    max_length/num_mutations
    
    rank_index<-rank(coords[which(coords[,2]==1),1])
    increment<-max_length/(num_mutations)
    new_coords<-seq(0+increment/2, max_length-increment/2, max_length/(num_mutations))
    coords[which(coords[,2]==1),1]<-new_coords[rank_index]
    
    # abbreviate vertex names for plot
    vertex_attr(g)$name<-gsub("TRUNCATING","TRUNC",vertex_attr(g)$name)
    vertex_attr(g)$name<-gsub("MISSENSE","MISS",vertex_attr(g)$name)
    
    # output bipartite graph for viewing
    pdf(paste0(PANCAN_DATA,tissue,"_",gsub("\\|","_",bait_gene_for_filename),"_bipartite.pdf"))
    #print(plot(g, layout=layout.bipartite, vertex.size=7, vertex.label.cex=0.5, vertex.label.color = "black"))
    print(plot(g, layout=coords, vertex.size=10, vertex.label.cex=0.5, asp = 0.75) +
            text(x=0,y=-1.2,"Samples") +
            text(x=0, y=1.3, "Mutations"))
    dev.off()
    write.csv(bipartite_data,file=paste0(PANCAN_DATA,tissue,"_",gsub("\\|","_",bait_gene_for_filename),"_bipartite.csv"),row.names = FALSE)
    
    # generate classifier scores again but with the smaller gene set after set coverage algorithm
    positive_set<-c()
    num_samples<-c()
    gene<-gene_list[1]
    for (gene in set_coverage_solution) {
      print(gene)
      gene_prefix<-unlist(strsplit(gene,"_"))[c(TRUE,FALSE)]
      gene_mut_type<-unlist(strsplit(gene,"_"))[c(FALSE,TRUE)]
      if (gene_mut_type=="SNV") {
        print("Combining all SNV mutations")
        temp_set_XY<-createXY_gmt(gene, gmt_file_data, feature_data, NA)
      } else {
        print(paste0("Selecting Mutation Type: ",gene_mut_type))
        temp_set_XY<-createXY_gmt(gene, gmt_file_data, feature_data, NA)
        
        #temp_set_XY<-createXY_pos(gene_prefix,mutation_data,expression_data,NA,mutation_type=gene_mut_type, cnv_data=cnv_data, fusion_data = fusion_data, pcawg_mutation_data=pcawg_mutation_data)
      }
      num_samples<-c(num_samples,length(rownames(temp_set_XY$Y)[which(temp_set_XY$Y==1)]))
      positive_set<-c(positive_set,rownames(temp_set_XY$Y)[which(temp_set_XY$Y==1)])
      
    }
    set_coverage_final_solution<-temp_set_XY
    set_coverage_final_solution$Y$sample<-0
    set_coverage_final_solution$Y$sample[which(rownames(set_coverage_final_solution$Y) %in% positive_set)]<-1
    
    # creating classifier scores for new set coverage solution
    new_model<-runClassifier_V2(set_coverage_final_solution$X,set_coverage_final_solution$Y,folds,num_permutations,alpha=1,weighted = TRUE,stratified = TRUE)
    print(paste("New Set Coverage Solution PR AUC:",mean(new_model$pr_auc)))
    print(paste("New Set Coverage Solution PRECISION:",mean(new_model$precision)))
    print(paste("New Set Coverage Solution RECALL:",mean(new_model$recall)))
    
    set_cover_classifier_scores<-((predict(new_model$model, s=min((new_model$min_lambda)[which(new_model$pr_auc==max(new_model$pr_auc))]), newx=data.matrix(set_coverage_final_solution$X), type="response")))[,1]
    set_cover_classifier_scores<-set_cover_classifier_scores[order(set_cover_classifier_scores,decreasing = TRUE)]
    
    set_cover_oncoprint_mutation_data<-data.frame(oncoprint_mutation_data[which(rownames(oncoprint_mutation_data) %in% toupper(set_coverage_solution)),],stringsAsFactors = FALSE)
    set_cover_oncoprint_mutation_data[is.na(set_cover_oncoprint_mutation_data)]<-""
    
    print(paste0("Routine Complete for ",gene,"_","____",gmt_file,"...writing files: "))
    # convert list of lists to dataframe for writing
    write.table(set_coverage_solution,file=filename3,quote=FALSE,sep="\t",row.names = FALSE)
    write.csv(set_cover_classifier_scores,file=filename4)
    write.csv(set_cover_oncoprint_mutation_data,file=filename5)
    write.csv(initial_bait_model$pr_auc,file=filename7)
    write.csv(new_model$pr_auc,file=filename8)
    
    # print oncoprint using yanni's oncoprint plot:
    system(paste0("python3 ",SCRIPTS,"oncoprint.py -memo -c ",filename4," -m ",filename5," -o ",filename6))
    
  } # END OF SET COVERAGE LOOP
  
} # END OF LURE FUNCTION


# This is the core function in the LURE Method
# It is a recursive function, the bait gene is essentially the base case
# At each node/iteration the KS test (GSEA preranked) is ran and this function is called on every event
# Most variables are just passed through from LURE function, please read LURE function for details
# tree_path is an internal variable that keeps track of each individual path.  
#   : Consider it a stack variable that is pushed new leaves to go down the tree and pops off leaves to back up the tree
#   : Set it initially to the starting oncogene
# tree_path_list is an external variable which is a list of the tree paths.  At the bottom of each tree, the current tree_path

positive_recursive_iterate<-function(bait, # name of gene
                                     pred_resp, # sample classifier scores
                                     tree_path, # current tree path
                                     gmt_file, # gmt file (mutation data)
                                     X, # feature file
                                     origY, # original labels
                                     orig_score_vector, # original CV test scores
                                     num_permutations, # number of permutations/iterations to run classifier
                                     num_events, # number of maximum events
                                     percent_overlap, # max percent of overlap allowed
                                     pvalue_threshold, 
                                     min_gene_set_size, 
                                     gsea_pvalue_threshold, 
                                     gsea_fdr_threshold,
                                     tree_length,
                                     max_tree_length, 
                                     folds, # CV folds
                                     enrichment_analysis_only) {
  print(paste0("Tree Length:",tree_length))
  # Run the GSEA Preranked Analysis to identify enriched events in high scored samples
  results<-run_gsea_V2(bait, gmt_file, X, origY, pred_resp, orig_score_vector=orig_score_vector, num_permutations, num_events=num_events, pvalue_threshold=pvalue_threshold, min_gene_set_size=min_gene_set_size,gsea_pvalue_threshold=gsea_pvalue_threshold, gsea_fdr_threshold=gsea_fdr_threshold, percent_overlap, folds = folds, enrichment_analysis_only=enrichment_analysis_only)
  
  # if enrichment analysis only is set to TRUE, then we do NOT run the LURE recursive process.  
  # We are returning the initial set of results
  if (enrichment_analysis_only==TRUE) {
    print("Enrichment Analysis Only...Returning results now")
    return(results)
  }
  # Check for results and if there are results, call the add event to positive class and re-run function
  if ((nrow(results$df_results) > 0) && (tree_length < max_tree_length)) {
    print(paste("Number of Events to analyze:",nrow(results$df_results)))
    i<-1  
    df_results<-results$df_results
    results$models
    for (i in 1:nrow(df_results)) {
      gene<-df_results$event[i]
      
      print(paste("Adding samples containing event",gene,"to the positive set to recreate a classifier and rerun GSEA and look for new events:"))
      # Didn't figure out how to pass a model variable so I just re-run the classifier here
      # add event samples to 1 in Y (add event samples to the positive set)
      # must create a new Y, otherwise it gets changed
      newY<-origY
      event_sample_names_inner_iteration<-strsplit(df_results$sample_names[i]," ")[[1]]
      pr_auc_scores<-as.numeric(strsplit(df_results$pr_auc_score[i]," ")[[1]])
      newY[which(rownames(newY) %in% event_sample_names_inner_iteration) , ]<-1
      
      print(paste("Length of Old Positive Sample Names:",length(origY[origY==1])))
      print(paste("Length of New Positive Sample Names:",length(event_sample_names_inner_iteration)))
      print(paste("Length of New and Old Positive Sample Names:",length(newY[newY==1])))
      pos_set<-rownames(newY)[which(newY==1)]
      # create a new test not containing the event samples
      new_XY<-createXY_pos("test","nada",X,pos_set=pos_set)
      
      big_model<-runClassifier_V2(X,newY,folds,num_permutations,alpha=1,weighted = TRUE,stratified = TRUE)
      # creating new response vector with all the samples without known event samples
      # old method was to choose the best lambda from all the iterations
      #new_pred_resp<-((predict(big_model$model, s=min((big_model$min_lambda)[which(big_model$pr_auc==max(big_model$pr_auc))]), newx=data.matrix(new_XY$test), type="response")))[,1]
      # new method averages the best lambdas from all the iterations.  It should remove some variability in the results
      new_pred_resp<-((predict(big_model$model, s=mean(big_model$min_lambda), newx=data.matrix(new_XY$test), type="response")))[,1]
      
      
      print(paste("Length of test samples:",length(new_pred_resp)))
      # new score vector is the previous(orig) classifier scores
      new_score_vector<-as.numeric(strsplit(df_results$pr_auc_score[i]," ")[[1]])
      
      # add gene to tree path and increment tree path length
      tree_path<-c(tree_path,gene)
      tree_length<-tree_length+1
      
      # output to source-target spreadsheet
      df_output<-data.frame("source"=df_results$bait[i],
                            "interaction"="positive",
                            "target"=gene,
                            "pvalue_over_orig"=df_results$pvalue_over_orig[i],
                            "pr_auc"=mean(pr_auc_scores))
      assign("STG_output", rbind(STG_output,df_output), envir = .GlobalEnv)
      
      # call recursive function
      positive_recursive_iterate(bait=gene,
                                 pred_resp=new_pred_resp,
                                 tree_path=tree_path, 
                                 gmt_file=gmt_file, 
                                 X=X, 
                                 origY=newY,
                                 orig_score_vector=new_score_vector,
                                 num_permutations=num_permutations,
                                 num_events=num_events,
                                 percent_overlap = percent_overlap,
                                 pvalue_threshold=pvalue_threshold, 
                                 min_gene_set_size=min_gene_set_size,
                                 gsea_pvalue_threshold=gsea_pvalue_threshold, 
                                 gsea_fdr_threshold=gsea_fdr_threshold,
                                 tree_length=tree_length,
                                 max_tree_length=max_tree_length,
                                 folds=folds,
                                 enrichment_analysis_only=enrichment_analysis_only
      )
      tree_path<-head(tree_path,-1)
      tree_length<-tree_length-1
    } # end of event for loop
  } else {
    # base case
    print("Bottom of Tree or Max Tree Length achieved")
    print(tree_path)
    assign("tree_path_list", c(tree_path_list, list(tree_path)), envir = .GlobalEnv)
  }
  return(tree_path_list)
}


# createXY_gmt
# creates X and Y variables, the features(X) and labels(Y)
# function takes a gene name, matches it to the gene name in the gmt file and gets the sample names
createXY_gmt<- function(gene_name_mutation_type, gmt_file_data, expression_file, pos_set) {
  NA_pos_set<-""
  if (is.na(pos_set[1])) {
    NA_pos_set<-NA
    pos_set<-gmt_file_data[which(toupper(gmt_file_data$X1)==toupper(gene_name_mutation_type)),]
    pos_set<-pos_set[,-c(1,2)]
  }

  X_pos<-expression_file[(row.names(expression_file) %in% pos_set) , ,drop=FALSE]
  X_pos<-X_pos[!(is.na(rownames(X_pos))) , ,drop=FALSE]
  old_X_pos<-X_pos
  X_pos_num<-length(X_pos[,1])
  
  print(paste("Number of Positive Samples: ", length(X_pos[,1])))
  # this was originally wrote so that it would run faster while initially creating classifier scores across Pancan
  # if (!(is.na(min_mutants))) {
  #   if (length(X_pos[,1]) < min_mutants) {
  #     print("Too few mutants")
  #     return(list("X" = NA, "Y" = NA,"test" = NA, "X_total" = NA, "X_pos_num" = NA))
  #   }
  # }

  # load negative training expression data for A
  neg_set<-row.names(expression_file)[!(row.names(expression_file) %in% pos_set)]
  X_neg<-expression_file[(row.names(expression_file) %in% neg_set) , ,drop=FALSE]
  X_neg_num<-length(X_neg[,1])
  print(paste("Number of Negative Samples: ", X_neg_num))
  
  X_total<-X_pos_num+X_neg_num
  Y_pos<-data.frame(matrix(1,length(X_pos[,1]),1))
  colnames(Y_pos)<-"sample"
  rownames(Y_pos)<-(rownames(X_pos))
  Y_neg<-data.frame(matrix(0,length(X_neg[,1]),1))
  colnames(Y_neg)<-"sample"
  rownames(Y_neg)<-(rownames(X_neg))
  # combine pos and neg
  X<-rbind(X_pos, X_neg)
  Y<-rbind(Y_pos, Y_neg)
  # here we set test to whatever is left
  if (is.na(NA_pos_set)) {
    #test<-expression_file[(row.names(expression_file) %in% unique(mutation_file$Tumor_Sample_Barcode)) , ]
    test<- expression_file[ ((row.names(expression_file) %in% neg_set)) , ]
  }
  else
    test<- expression_file[ (row.names(expression_file) %in% neg_set) , ]
  
  print(paste("Number of Remaining(test) Samples: ", length(test[,1])))
  # return variables
  # replace the period in rownames Y and replace with -
  return(list("X" = X, "Y" = Y,"test" = test, "X_total" = X_total, "X_pos_num" = X_pos_num))
}



# Runs a logistic regression classifier on X and Y, features(X) and labels(Y)

runClassifier_V2<- function(X, # features
                            Y, # labels
                            nfolds_input, # the number of cross validation folds 
                            iterations, # the number of times to run the model. (better to run more to get average lambda)
                            alpha, # elastic-net penalty, see below
                            weighted, # TRUE/FALSE weight the labels
                            stratified) { # TRUE/FALSE, stratify the folds
  # The elastic-net penalty is controlled by alpha, and bridges the gap between lasso (α=1, the default) and ridge (α=0).
  # lasso results in sparse coefficients
  # ridge regression results in more coefficients
  
  # calculate number of folds
  # to guarantee that at least one positive/negative label is in each fold, otherwise at least 10 samples per fold
  # first check if 10 observations exist per fold if 10 folds is chosen, if not set to min number of folds to keep at least 10 observations in each fold
  if (nrow(Y)/nfolds_input < 10) {
    nfolds_input<-floor(nrow(Y)/10)
  }
  
  nfolds<-min(max(3,nfolds_input),sum(Y==1),sum(Y==0))
  print(paste("Performing cross-validation using",nfolds,"folds"))
  print(paste("Alpha set to:",alpha,"(1=LASSO; fast; less coefficients, 0=Ridge Regression; slow; all coefficients)"))
  # create weights, weights are 1 minus the fraction of the pos or neg over the total number of labels.  
  # Basically smaller classes are given proportionally bigger weights
  fraction_0<-rep(1-sum(Y==0)/nrow(Y),sum(Y==0))
  fraction_1<-rep(1-sum(Y==1)/nrow(Y),sum(Y==1))
  # assign 1 - that value to a "weights" vector
  weights<-numeric(nrow(Y))
  if (weighted==TRUE) {
    weights[Y==0]<-fraction_0
    weights[Y==1]<-fraction_1
  } else {
    weights<-rep(1,nrow(Y))
  }
  
  # create an initial model and get a lambda, this is somewhat cheating as I am using glmnet's hotstart, but then do my own CV
  lambda_model<-glmnet(as.matrix(X), as.factor(Y[,1]), family = "binomial", weights = weights, nlambda = 100)
  lambda_start<-lambda_model$lambda
  
  
  
  # begin iterations
  pr_score_list<-as.numeric(list())
  lambda_list<-as.numeric(list())
  precision_score_list<-as.numeric(list())
  recall_score_list<-as.numeric(list())
  min_precision_lambda<-c()
  min_recall_lambda<-c()
  overall_recall_list<-c()
  overall_precision_list<-c()
  min_overall_recall_lambda<-c()
  min_overall_precision_lambda<-c()
  pr<-c()
  i<-1
  for (i in 1:iterations) {
    print(paste("Iteration:",i,"of",iterations))
    # here we set the seed, but allow it to change according to the iteration.
    # we will get consistent results, but increasing the number of iterations will still be beneficial as it will produce different
    # cross validation folds
    set.seed(69/i)
    if (stratified==TRUE) {
      # assign folds evenly using the mod operator
      fold0 <- sample.int(sum(Y==0)) %% nfolds
      fold1 <- sample.int(sum(Y==1)) %% nfolds
      foldid <- numeric(nrow(Y))
      foldid[Y==0] <- fold0
      foldid[Y==1] <- fold1
      foldid <- foldid + 1
    } else {
      foldid <- numeric(nrow(Y))
      foldid <- sample.int(nrow(Y)) %% nfolds
      foldid <- foldid+1
    }
    
    
    pr_auc_matrix<-matrix(ncol=100,nrow=nfolds)
    fold<-1
    # start parallel cross validation here
    cv_pr_score_list<-foreach(fold=1:(nfolds)) %dopar% {
      
      # train on all data this is not in the fold
      model<-glmnet(as.matrix(X[which(foldid!=fold) , ]), lambda=lambda_start, as.factor(Y[which(foldid!=fold),1]), family = "binomial", weights = weights[which(foldid!=fold)])
      # get the response of the data in the fold
      resp<-data.frame(predict(model, s=lambda_start, data.matrix(as.matrix(X[which(foldid==fold) , ])), type="response"))
      
      pr_auc_list<-rep(0,ncol(resp))
      precision_list<-rep(0,ncol(resp))
      recall_list<-rep(0,ncol(resp))
      col<-22
      for (col in 1:ncol(resp)) {
        resp_column<-resp[,col]
        if (abs(max(as.numeric(resp_column)) - min(as.numeric(resp_column))) == 0) {
          pr$auc.integral<-0
          pr$auc.davis.goadrich<-0
          pr_auc<-0
        }
        else {
          # calculate PR AUC 
          pr <- pr.curve( scores.class0=as.numeric(resp_column),weights.class0=as.numeric(Y[which(foldid==fold),1]),curve=FALSE)
          pr_auc<-pr$auc.integral
        }
        pr_auc_list[col]<-pr_auc
        
        # Calculate Precision and Recall
        # THIS CREATES A CUTOFF, better to use the PR AUC
        resp_column[resp_column>=.5]<-1
        resp_column[resp_column<.5]<-0
        
        TP<-length(which(resp_column==1 & Y[which(foldid==fold),1]==1))
        FP<-length(which(resp_column==1 & Y[which(foldid==fold),1]==0))
        FN<-length(which(resp_column==0 & Y[which(foldid==fold),1]==1))
        if ((TP+FP) > 0) { precision<-TP/(TP+FP) } else { precision<-0 }
        if ((TP+FN) > 0) { recall<-TP/(TP+FN) } else { recall <- 0 }
        precision_list[col]<-precision
        recall_list[col]<-recall
      }
      
      
      list(pr_auc_list,precision_list,recall_list)
      
    } # end of parallel loop
    # add up PR AUC scores
    list_of_auc_scores<-lapply(cv_pr_score_list, `[[`, 1)
    
    pr_auc_matrix<-matrix(unlist(list_of_auc_scores), ncol = nfolds, byrow = FALSE)
    # replace Nan's with zeros
    pr_auc_matrix[is.nan(pr_auc_matrix)]<-0
    average_auc_per_fold<-rowMeans(pr_auc_matrix)
    #plot(log(lambda_start), average_auc_per_fold)
    # calculate minimum lambda for the maximum PR
    min_lambda<-min(lambda_start[which(average_auc_per_fold==max(average_auc_per_fold))])
    lambda_list<-c(lambda_list,min_lambda)
    pr_score_list<-c(pr_score_list,max(average_auc_per_fold))
    
    list_of_precision<-lapply(cv_pr_score_list, `[[`, 2)
    precision_matrix<-matrix(unlist(list_of_precision), ncol = nfolds, byrow = FALSE)
    average_prec_per_fold<-rowMeans(precision_matrix)
    precision_score_list<-c(precision_score_list,max(average_prec_per_fold[which(average_auc_per_fold==max(average_auc_per_fold))]))
    
    list_of_recall<-lapply(cv_pr_score_list, `[[`, 3)
    recall_matrix<-matrix(unlist(list_of_recall), ncol = nfolds, byrow = FALSE)
    average_recall_per_fold<-rowMeans(recall_matrix)
    recall_score_list<-c(recall_score_list,max(average_recall_per_fold[which(average_auc_per_fold==max(average_auc_per_fold))]))
    
    min_recall_lambda<-c(min_recall_lambda,min(lambda_start[which(average_recall_per_fold==max(average_recall_per_fold))]))
    min_precision_lambda<-c(min_precision_lambda,min(lambda_start[which(average_prec_per_fold==max(average_prec_per_fold))]))
    
    overall_recall_list<-c(overall_recall_list,max(average_recall_per_fold))
    overall_precision_list<-c(overall_precision_list,max(average_prec_per_fold))
    
    min_overall_recall_lambda<-c(min_overall_recall_lambda,min(lambda_start[which(overall_recall_list==max(overall_recall_list))]))
    min_overall_precision_lambda<-c(min_overall_precision_lambda,min(lambda_start[which(overall_precision_list==max(overall_precision_list))]))
    
    
  } # end of iterations
  
  return(list("model" = lambda_model,
              "pr_auc" = pr_score_list,
              "precision_pr_auc"=precision_score_list,
              "recall_pr_auc"=recall_score_list,
              "min_lambda"=lambda_list,
              "min_recall_lambda"=min_recall_lambda,
              "min_precision_lambda"=min_precision_lambda,
              "overall_precision"=overall_precision_list,
              "overall_recall"=overall_recall_list,
              "overall_precision_lambda"=min_overall_precision_lambda,
              "overall_recall_lambda"=min_overall_recall_lambda))
  
  
  
} # end of RunClassifier_V2


# This function runs the GSEA preranked tool for LURE and performs a permutation test by adding in random samples to see if it does better
run_gsea_V2<-function(bait, 
                      gmt_file, 
                      X, # features
                      Y, # labels
                      pred_resp, # predicted response vector...passed to GSEA
                      orig_score_vector, # original set of classifier scores
                      num_permutations, # number of times to run classifier
                      num_events, # max number of events
                      pvalue_threshold, 
                      min_gene_set_size,
                      gsea_pvalue_threshold, 
                      gsea_fdr_threshold, 
                      percent_overlap, 
                      folds, 
                      enrichment_analysis_only) {
  no_col<-max(count.fields(paste0(DD_HOME,gmt_file),sep = "\t"))
  gmt_mutation_data<-(read.csv(paste0(DD_HOME,gmt_file),sep="\t",col.names=1:no_col,header = FALSE,stringsAsFactors = FALSE))
  
  print(paste("Bait:",bait))
  print(paste("LengthX:",length(X[,1])))
  print(paste("LengthY:",length(Y[,1])))

  print(paste("Percent Overlap",percent_overlap))
  print(paste("gmt file:",gmt_file))
  print(paste("num_permutations:",num_permutations))
  print(paste("num_events:",num_events))
  print(paste("pvalue_threshold:",pvalue_threshold))
  print(paste("min_gene_set_size:",min_gene_set_size))
  set.seed(NULL)
  rand<-sample(1:100000000,1)
  print(paste("Using Rand Num:",rand))
  print(paste("Original Scores:",orig_score_vector))
  
  # initialize output variable
  result_output<-data.frame(matrix(ncol=9,nrow=0))
  models_output<-c()

  # running GSEA on the pred_resp, note this is the "test" dataset, everything in Y that is not 1
  pred_resp<-pred_resp[order(pred_resp, decreasing=TRUE)]
  
  # check to see if there is at least the minumum positively classed sample
  if (sum(pred_resp>.5) < min_gene_set_size) {
    print("Not enough positively classified samples, not running GSEA...")
    return(list("df_results"=result_output,"models"=models_output))
  }
  
  # write data file for GSEA java app to read
  write.table(pred_resp, file=paste(DD_HOME,"rankedfile_",rand,".rnk",sep=""), quote=FALSE, sep="\t", col.names=FALSE)
  print("Running Pre-ranked GSEA for enrichment test...")
  gsea_cmd<-paste("java -cp ",DD_HOME,"gsea2-2.2.2.jar -Xmx15000m xtools.gsea.GseaPreranked -gmx ",DD_HOME,gmt_file," -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ",DD_HOME,"rankedfile_",rand,".rnk -scoring_scheme weighted -rpt_label my_analysis -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed 1234 -set_max 500 -set_min ",min_gene_set_size," -zip_report false -out ",DD_HOME,"GSEA_reports/",rand," -gui false > ",OUTPUT_DATA,"/GSEA_logfile_",rand,".log", sep="")

  # check return code from GSEA 
  return_code<-system(paste0(gsea_cmd, " 2> ",DD_HOME,"logfile",rand,".txt"))
  print(paste0("GSEA Return code:",return_code))
  if (return_code > 0) {
    print("Java failure!!!!")
    print(paste0("logfile located at ",DD_HOME,"logfile",rand,".txt"))
    logfile<-system(paste0("cat ",DD_HOME,"logfile",rand,".txt"), intern = TRUE)
    print(logfile)
    if (!(grepl("none of the gene sets passed size thresholds",logfile[1]))) {
      print("JAVA Failure unknown, stopping process")
      return(list("df_results"=result_output,"models"=models_output))
    }
    print("No Positive Found!")
    return(list("df_results"=result_output,"models"=models_output))
    
  }
  # process GSEA results, read output files
  setwd(paste(DD_HOME,"GSEA_reports/",rand,sep=""))
  dir_results<-list.files()
  dir_results[length(dir_results)]
  setwd(paste(DD_HOME,"GSEA_reports/",rand,"/",dir_results[length(dir_results)],sep=""))
  print("Loading GSEA results...")
  GSEA_output<-read.csv(list.files(pattern=".*gsea_report_for_na_pos_.*xls.*")   , sep="\t", row.names=1, stringsAsFactors = FALSE)
  GSEA_output<-GSEA_output[ (GSEA_output$NOM.p.val < gsea_pvalue_threshold) , ]
  GSEA_output<-GSEA_output[ (GSEA_output$FDR.q.val < gsea_fdr_threshold) , ]
  pos_events<-data.frame()
  # check to see if enrichment_analysis equals true, if so return just the GSEA analysis
  if (enrichment_analysis_only==TRUE)
    return(GSEA_output)
  
  # running normal LURE process 
  if (length(GSEA_output[,1]) == 0) {
    print("NO POSITIVE EVENTS FOUND!")
    return(list("df_results"=result_output,"models"=models_output))
  }
  # events found, proceeding with LURE
  print(paste(length(GSEA_output[,1]), "positive event(s) found!"))
  
  # load gmt data
  no_col<-max(count.fields(paste0(DD_HOME,gmt_file),sep = "\t"))
  gmt_mutation_data<-(read.csv(paste0(DD_HOME,gmt_file),sep="\t",col.names=1:no_col,header = FALSE,stringsAsFactors = FALSE))
  
  
  # take top 50 events events, make sure the GSEA setting to build the xls files is at 50 or more...
  GSEA_output<-GSEA_output[1:min(nrow(GSEA_output), 50) , ]
  print(row.names(GSEA_output))

  # Check for overlap between new event samples and pos_set. percent_overlap is passed in as an argument. 

  event_sample_names<-data.frame(matrix(nrow=0,ncol=2),stringsAsFactors = FALSE)
  colnames(event_sample_names)<-c("event","samples")
  row<-2
  
  for (row in 1:length(rownames(GSEA_output))) {
    print(paste("Processing:",toupper(rownames(GSEA_output)[row])))
    orig_pos_set<-rownames(Y)[Y==1]
    
    all_mutant_samples<-gmt_mutation_data[which(toupper(gmt_mutation_data$X1) == toupper(rownames(GSEA_output)[row])),]
    all_mutant_samples<-c(unlist(all_mutant_samples[,-c(1,2)]))
    all_mutant_samples<-unique(all_mutant_samples)
    num_mutants_in_orig_pos_set<-length(which(all_mutant_samples %in% orig_pos_set))
    percent_of_orig_pos_set_overlapping<-num_mutants_in_orig_pos_set/length(orig_pos_set)
    percent_of_mutants_in_orig_pos_set<-num_mutants_in_orig_pos_set/(length(all_mutant_samples)+num_mutants_in_orig_pos_set)
    
    all_mutant_samples<-all_mutant_samples[all_mutant_samples %in% rownames(Y)[Y==0]]
    print(paste0("New Mutant Samples:",length(all_mutant_samples)))

    print(paste0("Overlapping (co-mutated) Samples:",num_mutants_in_orig_pos_set))
    if (percent_of_orig_pos_set_overlapping > percent_overlap) {
      print("(Overlapping New-Orig Mutants)/(Original Mutants) ratio larger than set_overlap, skipping event...")
      next()
    }
    if (percent_of_mutants_in_orig_pos_set > percent_overlap) {
      print("(Overlapping New-Orig Mutants)/(New Mutant) ratio larger than set_overlap, skipping event...")
      next()
    }
    
    sample_names_for_output<-paste((all_mutant_samples),collapse = " ")
    # move our new set of events into our event_sample dataframe and limit the event size to num_events
    event_sample_names<-rbind(event_sample_names,data.frame("event"=rownames(GSEA_output)[row], "samples"=sample_names_for_output,stringsAsFactors = FALSE))
  }
  # limit the event size to num_events
  event_sample_names<-event_sample_names[1:min(nrow(event_sample_names),num_events) , ]
  print(paste("Now considering",nrow(event_sample_names),"event(s)..."))
  
  row<-1
  
  for (row in 1:nrow(event_sample_names)) {
    if (is.na(event_sample_names$event[row])) {
      print("No events left...")
      next()
    }
    print(paste("Processing:",event_sample_names$event[row]))
    #print(paste("old_sample_names",rownames((Y[Y==1]))))
    sample_names<-strsplit(event_sample_names$samples[row]," ")[[1]]
    print(paste("Length of Old Positive Sample Names:",length(Y[Y==1])))
    print(paste("Length of New Positive Sample Names:",length(sample_names)))
    pos_set<-c(rownames(Y)[which(Y==1)],sample_names)
    
    # create new X and Y 
    new_XY<-createXY_pos("test","nada",X,pos_set=pos_set)
    
    # run classifier for new positive set
    print("Running Classifier for Actual Event Samples")
    actual_model<-runClassifier_V2(new_XY$X,new_XY$Y,folds,num_permutations,alpha=1,weighted = TRUE,stratified = TRUE)
    actual_scores<-actual_model$pr_auc
    print(paste("ACTUAL NEW EVENT PR AUC:",mean(actual_scores)))
    print(paste("ACTUAL PR AUC:",actual_scores))
    random_scores<-c()
    for (i in 1:num_permutations) {
      print(paste("Permutation:",i))
      
      set.seed(69/i)
      
      random_set<-c(rownames(Y)[which(Y==1)],sample(rownames(Y)[which(Y==0)],length(sample_names)))
      random_XY<-createXY_pos("test","nada",X,pos_set=random_set)
      
      random_model<-runClassifier_V2(random_XY$X,random_XY$Y,folds,1,alpha=1, weighted=TRUE,stratified = TRUE)
      print(paste("RANDOM PR AUC:",mean(random_model$pr_auc)))
      random_scores<-c(random_scores,mean(random_model$pr_auc))
    }

    # Perform t test to determine if actual scores did better than random
    actual_vs_random<-t.test(actual_scores,random_scores,alternative="greater")
    print(paste("Actual Greater than Random Significance Value:",actual_vs_random$p.value))
    
    # Perform t test to determine if actual scores did not decrease from the original score
    
    actual_vs_orig<-t.test(actual_scores,orig_score_vector,alternative="greater")
    print(paste("Actual Greater than Orig Score Significance Value:",actual_vs_orig$p.value))
    print(paste("t statistic:(must be positive)",actual_vs_orig$statistic))
    
    
    if (actual_vs_random$p.value > (pvalue_threshold/num_permutations)) {
      print("Event did not do better than random...SKIPPING EVENT")
      next()
    }
    if (actual_vs_orig$statistic < 0) {
      print("Event lowered classifier score...SKIPPING EVENT")
      next()
    }
    actual_vs_random_pvalue<-actual_vs_random$p.value
    actual_vs_orig_pvalue<-actual_vs_orig$p.value
    
    
    # store data in dataframe
    print("Event passed thresholds...")
    print("Saving classifier scores and events")
    # create a response variable for new event
    #((predict(big_model$model, s=mean(big_model$min_recall_lambda), data.matrix(full_set_XY$test), type="response")))[,1]
    resp_output<-predict(actual_model$model, s=mean(actual_model$min_recall_lambda), data.matrix(new_XY$X), type="response")[,1]
    

    result_output<-rbind(result_output,data.frame("bait"=bait,
                                                  "event"=event_sample_names$event[row],
                                                  "sample_names"=event_sample_names$samples[row],
                                                  "actual_vs_random_pvalue"=actual_vs_random_pvalue,
                                                  "actual_vs_orig_pvalue"=actual_vs_orig_pvalue,
                                                  "pr_auc_score"=paste(actual_model$pr_auc,collapse=" "),
                                                  "mean_lambda"=mean(actual_model$min_recall_lambda),
                                                  "pvalue_over_random"=actual_vs_random$p.value,
                                                  "pvalue_over_orig"=actual_vs_orig$p.value,
                                                  stringsAsFactors = FALSE)) 
    models_output<-c(models_output, actual_model)
  } # end of event loop
  return(list("df_results"=result_output,"models"=models_output))
} # end of run_gseaV2 function

# Takes two matrices and the intersects the feature names (gene names) and outputs two matrices using only the intersected features
create_common_gene_matrix<-function(input_matrix1, input_matrix2) {
  # determine common gene names between 1 and 2
  common_genes<-intersect(colnames(input_matrix1),colnames(input_matrix2))
  # create new exp matrices 
  common_input_matrix1<-input_matrix1[ , colnames(input_matrix1) %in% common_genes]
  ordered_input_matrix1<-common_input_matrix1[ , order(colnames(common_input_matrix1)) ]
  
  common_input_matrix2<-input_matrix2[ , colnames(input_matrix2) %in% common_genes]
  ordered_input_matrix2<-common_input_matrix2[ , order(colnames(common_input_matrix2)) ]
  
  return(list("matrix1" = common_input_matrix1, "matrix2" = common_input_matrix2))
}



# t_test function, this is will prevent the error thrown the data is basically the same
run_t_test<-function(A, B) {
  myt_test <- try(t.test(A,B))
  if (inherits(myt_test, "try-error"))
  {
    cat(myt_test)
    myt_test$p.value <- 1
  }
  return(myt_test)
}


# function for Bonferroni FDR calculation
fdr <-
  function(NES,NESobs,NESnull) { 
    
    nrot <- length(NESnull)
    
    Nnp <- numeric(nrot)
    Nnl <- numeric(nrot)
    Nnn <- numeric(nrot)
    Nns <- numeric(nrot)
    
    #For each permutation, counting the number of positive/negative NES and the number of NES more 
    #extreme than the given NES in the null distribution
    
    if(NES >= 0) {
      Nnp <- sum(NESnull >= 0)
      Nnl <- sum(NESnull >= NES)
    } else {
      Nnn[i] <- sum(NESnull < 0) + 1
      Nns[i] <- sum(NESnull <= NES) + 1
    }
    #print(Nnl)
    #print(Nnp)
    #Counting the number of positive/negative NES and NES more extreme than the given NES,
    #among the NES of all gene sets to be tested. Finally, calculating the q-value.
    if(NES >= 0) {
      
      Np <- sum(NESobs >= 0)
      Nl <- sum(NESobs >= NES)
      Nnpl <- Nnl/Nnp 
      #print(Nl)
      #print(Nnpl)
      q <- (Nnpl)/(Nl/Np)
      
    } else {
      
      Nn <- sum(NESobs < 0)
      Ns <- sum(NESobs <= NES)
      Nnns <- Nns/Nnn
      q <- mean(Nnns)/(Ns/Nn)
      
    }
    if( q > 1) 
      q <- 1
    q
    #print("q=")
    #print(q)
  }





# This function performs a memo sort, the sorting algorithm used in oncoprint plots
memoSort <- function(M) {
  geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
  scoreCol <- function(x) {
    score <- 0;
    for(i in 1:length(x)) {
      if(x[i]) {
        score <- score + 2^(length(x)-i);
      }
    }
    return(score);
  }
  scores <- apply(M[geneOrder, ], 2, scoreCol);
  sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
  return(list("geneOrder" = geneOrder, "sampleOrder" = sampleOrder));
}

# This is an orphaned function, replaced by create_XY_gmt
createXY_pos<- function(gene_name, mutation_file, expression_file, pos_set, mutation_type, cnv_data, fusion_data, pcawg_mutation_data) {
  NA_pos_set<-""
  # CREATE X AND Y 
  if (is.na(pos_set[1])) {
    NA_pos_set<-NA
    if (is.na(mutation_type)) {
      pos_set<-unique(mutation_file$Tumor_Sample_Barcode[which(mutation_file$Hugo_Symbol %in% gene_name & !(mutation_file$Variant_Classification %in% c("Intron","Silent","3'UTR","5'UTR","3'Flank","5'Flank","IGR")))])
    } else {
      if (mutation_type=="SNV")
        pos_set<-unique(mutation_file$Tumor_Sample_Barcode[which(mutation_file$Hugo_Symbol %in% gene_name & (mutation_file$Variant_Classification %in% c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","Nonstop_Mutation")))])
      else if (toupper(mutation_type)=="MISSENSE")
        pos_set<-unique(mutation_file$Tumor_Sample_Barcode[which(mutation_file$Hugo_Symbol %in% gene_name & (mutation_file$Variant_Classification %in% c("Missense_Mutation")))])
      else if (toupper(mutation_type)=="TRUNCATING")
        pos_set<-unique(mutation_file$Tumor_Sample_Barcode[which(mutation_file$Hugo_Symbol %in% gene_name & (mutation_file$Variant_Classification %in% c("Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","Nonstop_Mutation")))])
      else if (toupper(mutation_type)=="SPLICESITE")
        pos_set<-unique(mutation_file$Tumor_Sample_Barcode[which(mutation_file$Hugo_Symbol %in% gene_name & (mutation_file$Variant_Classification %in% c("Splice_Site")))])
      else if (toupper(mutation_type)=="SPLICE")
        pos_set<-unique(mutation_file$Tumor_Sample_Barcode[which(mutation_file$Hugo_Symbol %in% gene_name & (mutation_file$Variant_Classification %in% c("Splice_Site")))])
      else if (toupper(mutation_type)=="AMP") {
        #sample_names<-colnames(cnv_data)[2:length(colnames(cnv_data))] 
        #colnames(cnv_data)[2:length(colnames(cnv_data))]<-gsub("-",".",strtrim(sample_names,12))
        gene_data<-cnv_data[which(toupper(cnv_data$`Gene Symbol`) == toupper(gene_name)), ]
        pos_set<-unique(colnames(gene_data)[which(gene_data %in% c(2))])
      }
      else if (toupper(mutation_type)=="DEL") {
        print("Deletion Event")
        #sample_names<-colnames(cnv_data)[2:length(colnames(cnv_data))] 
        #colnames(cnv_data)[2:length(colnames(cnv_data))]<-gsub("-",".",strtrim(sample_names,12))
        gene_data<-cnv_data[which(toupper(cnv_data$`Gene Symbol`) == toupper(gene_name)), ]
        pos_set<-unique(colnames(gene_data)[which(gene_data %in% c(-2))])
      }
      else if (toupper(mutation_type)=="FUSION") {
        print("Fusion Event")
        sample_names<-as.character(fusion_data[which(toupper(fusion_data$X1) == toupper(gene)),3:ncol(fusion_data)])
        pos_set<-sample_names[!(sample_names=="")]
      }
      else if (toupper(mutation_type)=="3UTR") {
        full_gene_name<-paste0(gene_name,"(",mutation_type,")")
        pos_set<-as.character(pcawg_mutation_data$V2[which(toupper(pcawg_mutation_data$gene)==full_gene_name)])
      }
      else if (toupper(mutation_type)=="5UTR") {
        full_gene_name<-paste0(gene_name,"(",mutation_type,")")
        pos_set<-as.character(pcawg_mutation_data$V2[which(toupper(pcawg_mutation_data$gene)==full_gene_name)])
      }
      else if (toupper(mutation_type)=="CDS") {
        full_gene_name<-paste0(gene_name,"(",mutation_type,")")
        pos_set<-as.character(pcawg_mutation_data$V2[which(toupper(pcawg_mutation_data$gene)==full_gene_name)])
      }
      else if (toupper(mutation_type)=="PROMCORE") {
        full_gene_name<-paste0(gene_name,"(",mutation_type,")")
        pos_set<-as.character(pcawg_mutation_data$V2[which(toupper(pcawg_mutation_data$gene)==full_gene_name)])
      }
      else if (toupper(mutation_type)=="ENH") {
        full_gene_name<-paste0(gene_name,"(",mutation_type,")")
        pos_set<-as.character(pcawg_mutation_data$V2[which(toupper(pcawg_mutation_data$gene)==full_gene_name)])
      }
      else if (toupper(mutation_type)=="PROMDOMAIN") {
        full_gene_name<-paste0(gene_name,"(",mutation_type,")")
        pos_set<-as.character(pcawg_mutation_data$V2[which(toupper(pcawg_mutation_data$gene)==full_gene_name)])
      }
      else {
        print("Mutation Type not recognized")
        return(list("X" = NA, "Y" = NA,"test" = NA, "X_total" = NA, "X_pos_num" = NA))
        
      }
    }
  }
  
  # new calculation for X_pos
  X_pos<-expression_file[(row.names(expression_file) %in% pos_set) , ,drop=FALSE]
  X_pos<-X_pos[!(is.na(rownames(X_pos))) , ,drop=FALSE]
  old_X_pos<-X_pos
  X_pos_num<-length(X_pos[,1])
  
  print(paste("Number of Positive Samples: ", length(X_pos[,1])))
  # load negative training expression data for A
  neg_set<-row.names(expression_file)[!(row.names(expression_file) %in% pos_set)]
  X_neg<-expression_file[(row.names(expression_file) %in% neg_set) , ,drop=FALSE]
  X_neg_num<-length(X_neg[,1])
  print(paste("Number of Negative Samples: ", X_neg_num))
  
  X_total<-X_pos_num+X_neg_num
  Y_pos<-data.frame(matrix(1,length(X_pos[,1]),1))
  colnames(Y_pos)<-"sample"
  rownames(Y_pos)<-(rownames(X_pos))
  Y_neg<-data.frame(matrix(0,length(X_neg[,1]),1))
  colnames(Y_neg)<-"sample"
  rownames(Y_neg)<-(rownames(X_neg))
  # combine pos and neg
  X<-rbind(X_pos, X_neg)
  Y<-rbind(Y_pos, Y_neg)
  # here we set test to whatever is left
  if (is.na(NA_pos_set)) {
    #test<-expression_file[(row.names(expression_file) %in% unique(mutation_file$Tumor_Sample_Barcode)) , ]
    test<- expression_file[ ((row.names(expression_file) %in% neg_set)) , ]
  }
  else
    test<- expression_file[ (row.names(expression_file) %in% neg_set) , ]
  
  print(paste("Number of Remaining(test) Samples: ", length(test[,1])))
  # return variables
  # replace the period in rownames Y and replace with -
  return(list("X" = X, "Y" = Y,"test" = test, "X_total" = X_total, "X_pos_num" = X_pos_num))
}
#pos_set<-unique(juri_file$V2[which(juri_file$gene==gene)])

#pos_set, feature_data, pcawg_clinical, .05)
#percent_neg_set<-.05
#feature_data<-junc_subset
downsample_pos_neg_set <- function(pos_set, junc_subset, pcawg_clinical, percent_neg_set) {
  
  pos_set<-sapply(strsplit(pos_set,","), "[[", 1)
  pos_set<-pos_set[which(pos_set %in% rownames(junc_subset))]
  
  tissues<-unique(pcawg_clinical$dcc_project_code[match(pos_set,pcawg_clinical$tumor_wgs_aliquot_id)])
  samples<-pcawg_clinical$tumor_wgs_aliquot_id[which(pcawg_clinical$dcc_project_code %in% tissues)]
  pos_set<-pcawg_clinical$tumor_wgs_aliquot_id[which((pcawg_clinical$dcc_project_code %in% tissues) & (pcawg_clinical$tumor_wgs_aliquot_id %in% pos_set))]
  
  
  set<-c()
  neg_set<-c()
  running_neg_set<-c()
  running_rand_pos_set<-c()
  running_rand_neg_set<-c()
  tissue<-tissues[1]
  for (tissue in tissues) {
    #which((pcawg_clinical$dcc_project_code == tissue))
    #which((pcawg_clinical$tumor_wgs_aliquot_id %in% pos_set))
    
    
    num_pos_set<-length(pcawg_clinical$tumor_wgs_aliquot_id[which((pcawg_clinical$dcc_project_code == tissue) & (pcawg_clinical$tumor_wgs_aliquot_id %in% pos_set))])
    all_set_tissue<-(pcawg_clinical$tumor_wgs_aliquot_id[which((pcawg_clinical$dcc_project_code == tissue) & (pcawg_clinical$tumor_wgs_aliquot_id %in% rownames(junc_subset)))])
    neg_set_tissue<-(pcawg_clinical$tumor_wgs_aliquot_id[which((pcawg_clinical$dcc_project_code == tissue) & (pcawg_clinical$tumor_wgs_aliquot_id %in% rownames(junc_subset)) & !(pcawg_clinical$tumor_wgs_aliquot_id %in% pos_set))])
    num_neg_set<-length(neg_set_tissue)
    num_all_set<-num_pos_set/percent_neg_set
    
    if (num_all_set<length(all_set_tissue)) {
      neg_set<-sample(neg_set_tissue,num_neg_set,replace = FALSE)
      random_pos_set<-sample(all_set_tissue,num_pos_set,replace=FALSE)
      random_neg_set<-sample(all_set_tissue[!(all_set_tissue %in% random_pos_set)], num_neg_set)
    } else {
      neg_set<-neg_set_tissue
      random_pos_set<-sample(all_set_tissue,num_pos_set,replace=FALSE)
      random_neg_set<-all_set_tissue[!(all_set_tissue %in% random_pos_set)]
      
      
      
    }
    #print(paste("ACTUAL:",tissue,num_pos_set,length(neg_set),length(all_set_tissue)))
    #print(paste("RANDOM:",tissue,length(random_pos_set),length(random_neg_set),length(all_set_tissue)))
    running_neg_set<-c(running_neg_set,neg_set)
    running_rand_pos_set<-c(running_rand_pos_set,random_pos_set)
    running_rand_neg_set<-c(running_rand_neg_set,random_neg_set)
    
  }
  set<-c(pos_set,running_neg_set)
  random_set<-c(running_rand_pos_set,running_rand_neg_set)
  random_pos_set<-sapply(strsplit(running_rand_pos_set,","), "[[", 1)
  random_set<-sapply(strsplit(random_set,","), "[[", 1)
  set<-sapply(strsplit(set,","), "[[", 1)
  junc_tissue_subset<-junc_subset[which(rownames(junc_subset) %in% set),]
  return(list("pos_set"=pos_set, "random_pos_set"=random_pos_set, "junc_tissue_subset"=junc_tissue_subset, "tissues"=tissues))
} # end of downsampling function

nPr <- function(n,r) {
  return(factorial(n)/factorial(n-r))
}



# run_mutual_exclusive_test<-function(seed, input_tree_path_file, all_sample_names, gmt_file) {
#   file<-input_tree_path_file
#   no_col <- max(count.fields(file))
#   if (no_col < 2) {
#     print(paste0("Seed: ",seed))
#     print("No results")
#     return(data.frame("Seed"=seed,"Num_Tree_Paths"=0,"Best_Tree_Path"=NA,"Score"=NA))
#   }
#   
#   tree_path_file<-read.csv(file,header=FALSE, sep="\t",col.names=1:no_col)
#   tree_path_file<-cbind("mutual_exclusive_ratio"=NA, tree_path_file)
#   tree_path_file<-cbind("comet_score"=NA, tree_path_file)
#   
#   #best_score<-.01
#   #tree<-7
#   #tree<-tree+1
#   
#   
#   #tree<-1
#   rand_dir<-sample(10^10,1)
#   dir.create(paste0("~/wext/",rand_dir))
#   for (tree in 1:nrow(tree_path_file)) {
#     #for (tree in 1:(nrow(tree_path_file)-4)) {
#     
#     print(tree)
#     line<-tree_path_file[tree,]
#     genes<-line[!(is.na(line))]
#     print(genes)
#     
#     # create mutation file for WEXT
#     output_file<-matrix(NA,ncol=length(genes)+1,nrow=length(all_sample_names))
#     rownames(output_file)<-all_sample_names
#     output_file[,1]<-all_sample_names
#     colnames(output_file)<-c("sample",genes)
#     row<-2
#     for (row in 1:length(genes)) {
#       return_code<-system(paste0("grep -i ",genes[row]," ",gmt_file, " > ~/capture_file"))
#       mutated_samples<-as.vector(t(read.table("~/capture_file"))[,1])[-c(1:2)]
#       tissue_mutated_samples<-mutated_samples[mutated_samples %in% all_sample_names]
#       sample<-tissue_mutated_samples[1]
#       for (sample in tissue_mutated_samples) {
#         output_file[which(rownames(output_file) == sample),which(colnames(output_file)==genes[row])]<-genes[row]
#       }
#     }
#     system("rm ~/wext/out.txt")
#     apply(output_file,1,function(x) cat(paste(x[!is.na(x)], collapse="\t"), sep="\n", file="~/wext/out.txt", append=TRUE))
#     
#     # process mutation data
#     system(paste0(DD_HOME,"venv/bin/python ",DD_HOME,"wext-master/process_mutations.py -m ~/wext/out.txt -o ~/wext/data.json -ct NA"))
#     
#     # create mutation probabilities
#     system(paste0(DD_HOME,"venv/bin/python ",DD_HOME,"wext-master/compute_mutation_probabilities.py     -mf ~/wext/data.json     -np 1000     -nc 4     -wf ~/wext/weights.npy     -s  12345     -v  1"))
#     
#     # find exclusive sets
#     rand_file2<-sample(10^10,1)
#     
#     system(paste0(DD_HOME,"venv/bin/python ",DD_HOME,"wext-master/find_exclusive_sets.py -mf ~/wext/data.json -o ~/wext/",rand_dir,"/output_",tree,"_",rand_file2," -ks ",length(genes)," -s Enumerate WRE -m Saddlepoint -wf ~/wext/weights.npy"))
#     
#   }
#   
#   
#   gene_set_output<-data.frame(matrix(ncol=2,nrow=0),stringsAsFactors = FALSE)
#   
#   setwd(paste0("~/wext/",rand_dir))
#   files<-list.files(paste0("~/wext/",rand_dir))
#   for (file in files) {
#     data_read<-read.table(file,sep="\t",stringsAsFactors = FALSE)
#     gene_set_output<-rbind(gene_set_output,data.frame(data_read$V1,as.numeric(data_read$V2)))
#   }
#   
#   gene_set_output<-gene_set_output[order(gene_set_output$as.numeric.data_read.V2.),]
#   best_tree_path<-as.character(gene_set_output[1,1])
#   best_tree_path_score<-as.numeric(gene_set_output[1,2])
#   gene_set_output<-cbind("Seed"=seed, gene_set_output)
#   #df_best_tree_path<-data.frame("Seed"=seed,"Num_Tree_Paths"=nrow(tree_path_file),"Best_Tree_Path"=best_tree_path,"Score"=best_tree_path_score)
# 
# return(gene_set_output)
# } # end of run_mutual_exclusive_test







### pcawg.colour.palette.R #########################################################################
#
# Authors: Jennifer Aguiar & Constance Li (constance.li@oicr.on.ca)
#

### DESCRIPTION ####################################################################################
#
# Return standard PCAWG colour palettes. Case insensitive.
# 
# To see all available schemes, set:			scheme='all', return.scheme = FALSE
# To return all full schemes, set: 			scheme='all', return.scheme = TRUE
# To return specific full schemes, set: 	scheme=<wanted scheme>, return.scheme = TRUE
#	 Note: x will be ignored when scheme='all' OR return.scheme = TRUE 
# To return colours for specific values, 

### ARGUMENTS ######################################################################################
# x 				Chracter vector with terms to be mapped to colours. Ignored if scheme='all' or 
#						return.scheme=TRUE
# scheme 			String specifying desired colour scheme. To see all available schemes, use 
#						scheme='all', returns.scheme=FALSE
# fill.colour 		Unrecognized output will be filled with this colour. Default to 'slategrey'
# return.scheme 	TRUE/FALSE. Set to true to return full specified scheme. Set to false to map
#						x to colours 
# 

### MAIN ###########################################################################################

pcawg.colour.palette <- function(
  x = NULL,
  scheme = NULL,
  fill.colour = 'slategrey',
  return.scheme = FALSE
) {
  
  # Define all colours 
  # Coding SNV mutation subtypes & consequences 
  nonsynonymous <- '#698B69'
  synonymous <- '#FFD700'
  stop.gain <- '#8B4789'
  stop.loss <- '#DA70D6'
  indel.frameshift <- '#FF8C00'
  indel.nonframeshift <- '#003366'
  splicing <- '#00CED1'
  # Non-coding SNV mutation subtypes, consequences & gene types
  non.coding <- '#A80015'
  promoter <- '#4C191E'
  enhancer <- '#7F000F'
  operator <- '#A84955'
  silencer <- '#E78A96'
  insulator <- '#FFC1C9'
  lncRNA <- '#331900'
  sncRNA <- '#594027'
  tRNA <- '#A87849'
  rRNA <- '#E7B98A'
  miRNA <- '#FFE0C1'
  utr5.utr3 <- '#1A1A1A'
  intronic <- '#4D4D4D'
  intergenic <- '#7F7F7F'
  telomeres <- '#B3B3B3'
  # Structral variant mutation subtypes
  cna.gain <- '#FF0000'
  cna.loss <- '#0000FF'
  inversion <- '#FFA500'
  transposition <- '#7300E7'
  translocation <- '#458B00'
  # Chromosomes
  chr1 <- '#DE47AB'
  chr2 <- '#72BE97'
  chr3 <- '#F7F797'
  chr4 <- '#7C749B'
  chr5 <- '#E85726'
  chr6 <- '#B395F8'
  chr7 <- '#DC8747'
  chr8 <- '#96D53D'
  chr9 <- '#DC85EE'
  chr10 <- '#7D32B3'
  chr11 <- '#88DB68'
  chr12 <- '#78AAF1'
  chr13 <- '#D9C6CA'
  chr14 <- '#336C80'
  chr15 <- '#F7CA44'
  chr16 <- '#32C7C7'
  chr17 <- '#D4C5F2'
  chr18 <- '#995493'
  chr19 <- '#F88B78'
  chr20 <- '#475ECC'
  chr21 <- '#E0BD8C'
  chr22 <- '#9E2800'
  chrX <- '#F2BBD2'
  chrY <- '#B6EBEA'
  # Sex 
  male <- '#B6EBEA'
  female <- '#F2BBD2'
  # Tumour stage 
  st.one <- '#FFFFFF'
  st.two <- '#FFFF00'
  st.three <- '#FFA500'
  st.four <- '#FF0000'
  st.one.two <- '#FFE4B5'
  st.one.three <- '#EEE8AA'
  st.two.one <- '#FFD700'
  st.two.three <- '#F4A460'
  # TNM Cat
  tnm.zero <- '#FFFFFF'
  tn.one <- '#FFD399'
  tn.two <- '#FFAE45'
  tn.three <- '#B87217'
  tn.four <- '#774607'
  m.one <- '#000000'
  tnm.x <- '#708090'
  # Grade 
  gr.one <- '#FFFFFF'
  gr.two <- '#9CF0FC'
  gr.three <- '#335FE5'
  gr.four <- '#003366'
  gr.well <- '#FFFFFF'
  gr.mod <- '#9CE750'
  gr.poor <- '#00CC00'
  gr.un <- '#005900'
  pr.three.three <- '#FFFFFF'
  pr.three.four <- '#FFFF00'
  pr.three.five <- '#CD2990'
  pr.four.three <- '#FFA500'
  pr.four.four <- '#FF0000'
  pr.four.five <- '#A52A2A'
  pr.five.three <- '#8B008B'
  pr.five.four <- '#0000CD'
  pr.five.five <- '#000000'
  # Primary or Met 
  primary <- '#FFFFFF'
  metastatic <- '#7217A5'
  # Generic 
  other <- '#E5E5E5'
  unknown <- '#708090'
  # Tumour Subtype 
  biliary.adenoca <- '#00CD66'
  bladder.tcc <- '#EEAD0E' 
  bone.osteosarc <- '#FFD700'
  softtissue.leiomyo <- '#FFEC8B'
  softtissue.liposarc <- '#CDCB50'
  bone.epith <- '#ADAC44'
  breast.adenoca <- '#CD6090' 
  cervix.scc <- '#79CDCD' 
  cns.medullo <- '#D8BFD8'
  cns.piloastro <- '#B0B0B0' 
  cns.gbm <- '#3D3D3D' 
  cns.gbm.alt <- '#4A4A4A' 
  cns.oligo <- '#787878'
  colorect.adenoca <- '#191970' 
  eso.adenoca <- '#1E90FF'
  head.scc <- '#8B2323'
  kidney.rcc <- '#FF4500' 
  kidney.chrcc <- '#B32F0B' 
  liver.hcc <- '#006400' 
  lung.scc <- '#FDF5E6'
  lung.adenoca <- '#FFFFFF' 
  lymph.bnhl <- '#698B22'
  lymph.cll <- '#698B22' 
  myeloid.mpn <- '#FFC100' 
  myeloid.aml <- '#CD6600' 
  ovary.adenoca <- '#008B8B' 
  panc.adenoca <- '#7A378B'
  panc.endocrine <- '#E066FF' 
  prost.adenoca <- '#87CEFA'
  skin.melanoma <- '#000000' 
  stomach.adenoca <- '#BFEFFF' 
  thy.adenoca <- '#9370DB'
  uterus.adenoca <- '#FF8C69'
  bone.cart <- '#DDCDCD' 
  bone.cart.alt <- '#F0EE60'
  breast.lobularca <- '#DDCDCD' 
  breast.lobularca.alt <- '#F095BD'
  breast.dcis <- '#DDCDCD'
  lymph.nos <- '#DDCDCD'
  lymph.nos.alt <- '#698B22'
  myeloid.mds <- '#DDCDCD'
  cervix.adenoca <- '#DDCDCD'
  
  #-----------------------------------------------------------------------------------------------
  # Some input checking & processing 
  if (class(x) == 'factor') {
    stop('x cannot be a factor: please coerce to character before passing')
  }
  # Some parameters override provided input x. 
  if ((scheme == 'all' || return.scheme) && length(x) != 0) {
    warning('Input x ignored when scheme = \'all\' OR return.scheme = TRUE. Returning all schemes')
  }
  scheme <- tolower(scheme)
  x.input <- tolower(x)
  x.input <- gsub('-', '.', x.input)
  if (return.scheme || scheme == 'all') {
    x.input <- NULL
  }
  
  colour.schemes <- list(
    coding.snv = list(
      levels = c(
        'nonsynonymous',
        'synonymous',
        'stopgain',
        'stoploss',
        'indel.frameshift',
        'indel.nonframeshift',
        'splicing'
      ),
      colours = c(
        nonsynonymous,
        synonymous,
        stop.gain,
        stop.loss,
        indel.frameshift,
        indel.nonframeshift,
        splicing
      )
    ),
    noncoding.snv = list(
      levels = c(
        'noncoding',
        'promoter',
        'enhancer',
        'operator',
        'silencer',
        'insulator',
        'lncrna',
        'sncrna',
        'trna',
        'rrna',
        'mirna',
        'utr5.utr3',
        'intronic',
        'intergenic',
        'telomeres'
      ),
      colours = c(
        non.coding,
        promoter,
        enhancer,
        operator,
        silencer,
        insulator,
        lncRNA,
        sncRNA,
        tRNA,
        rRNA,
        miRNA,
        utr5.utr3,
        intronic,
        intergenic,
        telomeres
      )
    ),
    # ADD TFBS
    structural.variants = list(
      levels = c(
        'cna.gain',
        'cna.loss',
        'inversion',
        'transposition',
        'translocation'
      ),
      colours = c(
        cna.gain,
        cna.loss,
        inversion,
        transposition,
        translocation
      )
    ),
    chromosomes = list(
      levels	= c(
        '1',
        '2',
        '3',
        '4',
        '5',
        '6',
        '7',
        '8',
        '9',
        '10',
        '11',
        '12',
        '13',
        '14',
        '15',
        '16',
        '17',
        '18',
        '19',
        '20',
        '21',
        '22',
        'x',
        'y'
      ),
      colours = c(
        chr1,
        chr2,
        chr3,
        chr4,
        chr5,
        chr6,
        chr7,
        chr8,
        chr9,
        chr10,
        chr11,
        chr12,
        chr13,
        chr14,
        chr15,
        chr16,
        chr17,
        chr18,
        chr19,
        chr20,
        chr21,
        chr22,
        chrX,
        chrY
      )
    ),
    sex = list(
      levels = c(
        'male',
        'female'
      ),
      colours = c(
        male,
        female
      )
    ),
    stage.arabic = list(
      levels = c(
        '1',
        '2',
        '3',
        '4'
      ),
      colours = c(
        st.one,
        st.two,
        st.three,
        st.four
      )
    ),
    stage.roman = list(
      levels = c(
        'i',
        'i.ii',
        'i.iii',
        'ii',
        'ii.i',
        'ii.iii',
        'iii',
        'iv'
      ),
      colours = c(
        st.one,
        st.one.two,
        st.one.three,
        st.two,
        st.two.one,
        st.two.three,
        st.three,
        st.four
      )
    ),
    t.category = list(
      levels	= c(
        '0',
        '1',
        '2',
        '3',
        '4',
        'x'
      ),
      colours = c(
        tnm.zero,
        tn.one,
        tn.two,
        tn.three,
        tn.four,
        tnm.x
      )
    ),
    n.category = list(
      levels	= c(
        '0',
        '1',
        '2',
        '3',
        '4',
        'x'
      ),
      colours = c(
        tnm.zero,
        tn.one,
        tn.two,
        tn.three,
        tn.four,
        tnm.x
      )
    ),
    m.category = list(
      levels = c(
        '0',
        '1',
        'x'
      ),
      colours = c(
        tnm.zero,
        m.one,
        tnm.x
      )
    ),
    grade = list(
      levels = c(
        'G1',
        'G2',
        'G3',
        'G4'
      ),
      colours = c(
        gr.one,
        gr.two,
        gr.three,
        gr.four
      )
    ),
    grade.word = list(
      levels = c(
        'well.differentiated',
        'moderately.differentiated',
        'poorly.differentiated',
        'undifferentiated'
      ),
      colours = c(
        gr.well,
        gr.mod,
        gr.poor,
        gr.un
      )
    ),
    prostate.grade = list(
      levels	= c(
        '3+3',
        '3+4',
        '3+5',
        '4+3',
        '4+4',
        '4+5',
        '5+3',
        '5+4',
        '5+5'
      ),
      colours = c(
        pr.three.three,
        pr.three.four,
        pr.three.five,
        pr.four.three,
        pr.four.four,
        pr.four.five,
        pr.five.three,
        pr.five.four,
        pr.five.five
      )
    ),
    primary.met = list(
      levels = c(
        'primary',
        'metastatic'
      ),
      colours = c(
        primary,
        metastatic
      )
    ),
    tumour.subtype = list(
      levels	= c(
        'biliary.adenoca',
        'bladder.tcc',
        'bone.osteosarc',
        'softtissue.leiomyo', 
        'softtissue.liposarc',
        'bone.epith',
        'breast.adenoca',
        'cervix.scc',
        'cns.medullo',
        'cns.piloastro',
        'cns.gbm',
        'cns.gbm.alt',
        'cns.oligo',
        'colorect.adenoca',
        'eso.adenoca',
        'head.scc',
        'kidney.rcc',
        'kidney.chrcc',
        'liver.hcc',
        'lung.scc',
        'lung.adenoca',
        'lymph.bnhl',
        'lymph.cll',
        'myeloid.mpn',
        'myeloid.aml',
        'ovary.adenoca',
        'panc.adenoca',
        'panc.endocrine',
        'prost.adenoca',
        'skin.melanoma',
        'stomach.adenoca',
        'thy.adenoca',
        'uterus.adenoca',
        'bone.cart',
        'bone.cart.alt',
        'breast.lobularca',
        'breast.lobularca.alt',
        'breast.dcis',
        'lymph.nos',
        'lymph.nos.alt',
        'myeloid.mds',
        'cervix.adenoca'
      ),
      colours = c(
        biliary.adenoca,
        bladder.tcc,
        bone.osteosarc,
        softtissue.leiomyo,
        softtissue.liposarc,
        bone.epith,
        breast.adenoca,
        cervix.scc,
        cns.medullo,
        cns.piloastro,
        cns.gbm,
        cns.gbm.alt,
        cns.oligo,
        colorect.adenoca,
        eso.adenoca,
        head.scc,
        kidney.rcc,
        kidney.chrcc,
        liver.hcc,
        lung.scc,
        lung.adenoca,
        lymph.bnhl,
        lymph.cll,
        myeloid.mpn,
        myeloid.aml,
        ovary.adenoca,
        panc.adenoca,
        panc.endocrine,
        prost.adenoca,
        skin.melanoma,
        stomach.adenoca,
        thy.adenoca,
        uterus.adenoca,
        bone.cart,
        bone.cart.alt,
        breast.lobularca,
        breast.lobularca.alt,
        breast.dcis,
        lymph.nos,
        lymph.nos.alt,
        myeloid.mds,
        cervix.adenoca
      )
    ),
    organ.system = list(
      levels = c(
        'biliary',
        'bladder',
        'bone.softtissue',
        'breast',
        'cervix',
        'cns',
        'colon.rectum',
        'esophagus',
        'head.neck',
        'kidney',
        'liver',
        'lung',
        'lymphoid',
        'myeloid',
        'ovary',
        'pancreas',
        'prostate',
        'skin',
        'stomach',
        'thyroid',
        'uterus'
      ),
      colours = c(
        biliary.adenoca,
        bladder.tcc,
        softtissue.leiomyo,
        breast.adenoca,
        cervix.scc,
        cns.oligo,
        colorect.adenoca,
        eso.adenoca,
        head.scc,
        kidney.rcc,
        liver.hcc,
        lung.scc,
        lymph.bnhl,
        myeloid.aml,
        ovary.adenoca,
        panc.adenoca,
        prost.adenoca,
        skin.melanoma,
        stomach.adenoca,
        thy.adenoca,
        uterus.adenoca
      )
    )
  )
  
  # Error if wanted scheme doesn't match existing schemes
  if (is.null(colour.schemes[[scheme]]) && scheme != 'all'){
    stop('Scheme not found!')
  }
  # Return full specified schemes if return.scheme is TRUE
  if (return.scheme & 'all' == scheme) {
    return(colour.schemes)
  } else if (return.scheme & 'all' != scheme) {
    return(colour.schemes[[scheme]])
  } else if (!return.scheme & 'all' == scheme) {
    return(names(colour.schemes))
  }
  
  # Form output colours 
  matched <- match(x.input, colour.schemes[[scheme]]$levels);
  x.colours <- colour.schemes[[scheme]]$colours[matched]
  names(x.colours) <- colour.schemes[[scheme]]$levels[matched]
  
  # Deal with unrecognized input by setting to fill colour (slategray by default)
  if (any(is.na(x.colours))) {
    warning('Unrecognized input value for x. Default to fill.colour.')
  }
  x.colours[which(is.na(x.colours))] <- fill.colour;
  
  return(x.colours)
}


create_output_matrix<-function(gene_list, sample_names, mutation_all, non_coding_mutations, cnv_data) {
  # creates output matrix for oncoprint for the telomere marker paper
  #cluster should = "Tumor Cluster1" or "Tumor Cluster2", etc
  
  
  #sample_names<-as.character(tumor_telomere_data2$Tumor_Barcode[which(tumor_telomere_data2$new_cluster ==cluster)])
  output_matrix<-matrix(nrow=0,ncol=length(sample_names))
  colnames(output_matrix)<-sample_names
  i<-4
  for (i in 1:length(gene_list)) {
    print(gene_list[i])
    new_row<-rep("",length(colnames(output_matrix)))
    
    
    # I'm removing the mutation_all analysis....it is giving different results than juri's file
    # SNVs: trunc
    all_gene_mutations<-mutation_all[(which(mutation_all$Hugo_Symbol==gene_list[i])) , ]
    trunc_mutations<-all_gene_mutations[ c(which(all_gene_mutations$Variant_Classification=="Frame_Shift_Ins") , which(all_gene_mutations$Variant_Classification=="Frame_Shift_Del"), which(all_gene_mutations$Variant_Classification=="De_novo_Start_OutOfFrame"), which(all_gene_mutations$Variant_Classification=="Nonsense_Mutation")) , ]
    new_row[which(colnames(output_matrix)%in%trunc_mutations$Tumor_Sample_Barcode)]<-paste0(new_row[which(colnames(output_matrix)%in%trunc_mutations$Tumor_Sample_Barcode)],"TRUNC;")
    
    # SNVs: splice site
    splice_mut<-all_gene_mutations[ c(which(all_gene_mutations$Variant_Classification=="Splice_Site")) , ]
    new_row[which(colnames(output_matrix)%in%splice_mut$Tumor_Sample_Barcode)]<-paste0(new_row[which(colnames(output_matrix)%in%splice_mut$Tumor_Sample_Barcode)],"SPLICE;")
    
    # structural variants
    if (gene_list[i]=="ATRX") {
      new_row[which(colnames(output_matrix)%in%ATRX_sv_samples)]<-paste0(new_row[which(colnames(output_matrix)%in%ATRX_sv_samples)],"SV;")
    }
    if (gene_list[i]=="RB1") {
      new_row[which(colnames(output_matrix)%in%RB1_sv_samples)]<-paste0(new_row[which(colnames(output_matrix)%in%RB1_sv_samples)],"SV;")
    }
    if (gene_list[i]=="DAXX") {
      new_row[which(colnames(output_matrix)%in%DAXX_sv_samples)]<-paste0(new_row[which(colnames(output_matrix)%in%DAXX_sv_samples)],"SV;")
    }
    
    
    # SNVs: missense
    misssense_mut<-all_gene_mutations[ c(which(all_gene_mutations$Variant_Classification=="Missense_Mutation")) , ]
    new_row[which(colnames(output_matrix)%in%misssense_mut$Tumor_Sample_Barcode)]<-paste0(new_row[which(colnames(output_matrix)%in%misssense_mut$Tumor_Sample_Barcode)],"SNV;")
    # copy number amp or del
    cnv_gene_data<-cnv_data[which(cnv_data$`Gene Symbol` == gene_list[i]) , ]
    cnv_gene_loss<-colnames(cnv_gene_data)[which(cnv_gene_data %in% c(-2))]
    new_row[which(colnames(output_matrix)%in%cnv_gene_loss)]<-paste0(new_row[which(colnames(output_matrix)%in%cnv_gene_loss)],"DEEPDEL;")
    if (gene_list[i]=="RB1") {
      cnv_gene_loss<-colnames(cnv_gene_data)[which(cnv_gene_data %in% c(-1))]
      new_row[which(colnames(output_matrix)%in%cnv_gene_loss)]<-paste0(new_row[which(colnames(output_matrix)%in%cnv_gene_loss)],"SHALLDEL;")
    }
    
    # promoter mutation
    if (!(is.na(non_coding_mutations))) {
      unique(non_coding_mutations$region_type)
      prom_core_mutations<-non_coding_mutations[which(non_coding_mutations$region_type=="gc19_pc.promCore.bed")]
      prom_core_mutations<-rbind(prom_core_mutations,non_coding_mutations[which(non_coding_mutations$region_type=="gc19_pc.promDomain.bed")])
      
      prom_core_gene_mutations_samples<-prom_core_mutations$V2[(grep(paste0("::",gene_list[i],"::"),prom_core_mutations$reg_id))]
      new_row[which(colnames(output_matrix)%in%prom_core_gene_mutations_samples)]<-paste0(new_row[which(colnames(output_matrix)%in%prom_core_gene_mutations_samples)],"PROMSNV;")
      
      # 3' UTR mutations
      unique(non_coding_mutations$region_type)
      prom_core_mutations<-non_coding_mutations[which(non_coding_mutations$region_type=="gc19_pc.3utr.bed")]
      prom_core_gene_mutations_samples<-prom_core_mutations$V2[(grep(paste0("::",gene_list[i],"::"),prom_core_mutations$reg_id))]
      new_row[which(colnames(output_matrix)%in%prom_core_gene_mutations_samples)]<-paste0(new_row[which(colnames(output_matrix)%in%prom_core_gene_mutations_samples)],"UTR3SNV;")
      
      # 5' UTR mutations
      unique(non_coding_mutations$region_type)
      prom_core_mutations<-non_coding_mutations[which(non_coding_mutations$region_type=="gc19_pc.5utr.bed")]
      prom_core_gene_mutations_samples<-prom_core_mutations$V2[(grep(paste0("::",gene_list[i],"::"),prom_core_mutations$reg_id))]
      new_row[which(colnames(output_matrix)%in%prom_core_gene_mutations_samples)]<-paste0(new_row[which(colnames(output_matrix)%in%prom_core_gene_mutations_samples)],"UTR5SNV;")
    }
    
    # end of gene loop
    output_matrix<-rbind(output_matrix, new_row)
    rownames(output_matrix)[i]<-gene_list[i]
  }  
  return("output_matrix"=output_matrix)
}



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# bait_gene="SF3B1_missense_set1"
# gmt_file="positive_control_SF3B1_missense.gmt"
# feature_data=feature_data
# num_permutations=5
# max_num_events=5
# percent_overlap=0
# LURE_pvalue_threshold=.05
# min_gene_set_size=5
# gsea_pvalue_threshold=.05
# gsea_fdr_threshold=.25
# max_tree_length = 5
# folds=10
# enrichment_analysis_only=FALSE
# output_file_prefix="Pos_Ctrl"


# make tumor type legend
make_pancan_tumor_type_legend<- function() {
  tissues <- c("BRCA",
               "PRAD","TGCT","KICH","KIRP","KIRC","BLCA",
               "OV","UCS","CESC","UCEC",
               "THCA","PCPG","ACC",
               "SKCM",
               "UVM","HNSC",
               "SARC",
               "ESCA","STAD","COAD","READ",
               "CHOL","PAAD","LIHC",
               "MESO","LUSC","LUAD",
               "GBM","LGG",
               "DLBC","LAML","THYM")
  tissue_colors = c("#ED2891","#7E1918","#BE1E2D","#ED1C24","#EA7075","#F8AFB3","#FAD2D9",
                    "#D97D25","#F89420","#F6B667","#FBE3C7",
                    "#F9ED32","#E8C51D","#C1A72F",
                    "#BBD642",
                    "#009444","#97D1A9",
                    "#00A99D",
                    "#007EB5","#00AEEF","#9EDDF9","#DAF1FC",
                    "#104A7F","#6E7BA2","#CACCDB",
                    "#542C88","#A084BD","#D3C3E0",
                    "#B2509E","#D49DC7",
                    "#3953A4","#754C29","#CEAC8F"
  )
  names(tissue_colors) <- tissues
  pdf("~/pancan_legend.pdf", height=12)
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("topleft", legend =names(tissue_colors), pch=15, pt.cex=3, cex=1.5, bty='n',
         col = tissue_colors)
  mtext("Tumor Type", at=0.14, cex=2)
  dev.off()
}



