#load necessary packages
library(glmnet)
library(plyr)
library(matrixStats)
library(data.table)
library(methods)
library(doMC)
library(tools)
library(dplyr)
library(PRROC)
library(ggplot2)
library(preprocessCore)
library(optparse)
library(RcppGreedySetCover)
library(igraph)



# dumpster diver functions

# load PCAWG clinical
#clinical<-read.csv(paste("~/Downloads/release_may2016.v1.4.tsv"), sep="\t", row.names=1, stringsAsFactors = FALSE)


#mutation_file<-mutation_pancan
#expression_file<-t(expression)

#### ---- FUNCTIONS
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


#mutation_file<-mutation_data
#expression_file<-expression_data
#XY<-createXY_pos(gene_name, TERT_mutant_file, pcawg_exp)
#pos_set<-NA
#gene_name<-"TP53"
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


createXY_pos_neg_list<- function(gene_name, pos_set, neg_set, expression_file) {
  
  # CREATE X AND Y 
  X_pos<-expression_file[(row.names(expression_file) %in% pos_set) , ,drop=FALSE]
  X_pos<-X_pos[!(is.na(rownames(X_pos))) , ,drop=FALSE]
  X_pos_num<-length(X_pos[,1])
  
  print(paste("Number of Positive Samples: ", length(X_pos[,1])))
  # load negative training expression data for A
  X_neg<-expression_file[(row.names(expression_file) %in% neg_set) , ,drop=FALSE]
  X_neg_num<-length(X_neg[,1])
  print(paste("Number of Negative Samples: ", X_neg_num))
  
  # check and adjust positive and negative sets to decrease chances of overfitting... we want 30% / 70% minimum
  X_total<-X_pos_num+X_neg_num
  #if ((X_neg_num/X_total < .75) && (X_neg_num/X_total > .25)) {
  #  print("Ratio of pos/neg is ok...")
  #} else if (X_neg_num/X_total > .75) {
  #  new_X_neg<-(.75/.25*X_pos_num)
  #  X_neg<-data.frame(X_neg[!(is.na(rownames(X_neg))) , ])
  #  X_neg<-data.frame(X_neg[1:new_X_neg,])
  #  print (paste("Modified Number of Negatives:", length(X_neg[,1])))
  #} else if (X_neg_num/X_total < .25) {
  #  new_X_pos<-(.75/.25*X_neg_num)
  #  X_pos<-X_pos[!(is.na(rownames(X_pos))) , ]
  #  X_pos<-X_pos[1:new_X_pos,]
  #  print (paste("Modified Number of Positives:", length(X_pos[,1])))
  #}
  Y_pos<-data.frame(matrix(1,length(X_pos[,1]),1))
  colnames(Y_pos)<-"sample"
  rownames(Y_pos)<-make.names(rownames(X_pos), unique=TRUE)
  
  
  Y_neg<-data.frame(matrix(0,length(X_neg[,1]),1))
  colnames(Y_neg)<-"sample"
  rownames(Y_neg)<-make.names(rownames(X_neg), unique=TRUE)
  # combine pos and neg
  X<-rbind(X_pos, X_neg)
  Y<-rbind(Y_pos, Y_neg)
  # here we set test to whatever is left
  #test<-expression_file[(row.names(expression_file) %in% unique(mutation_file$Tumor_Sample_Barcode)) , ]
  #test<- test[ !((row.names(test) %in% row.names(X_pos))) , ]
  
  
  # return variables
  return(list("X" = X, "Y" = Y, "X_total" = X_total, "X_pos_num" = X_pos_num))
}







createXY_pos_neg<- function(pos_gene_name, neg_gene_name, mutation_file, expression_file) {
  
  # CREATE X AND Y 
  pos_set<-unique(mutation_file$sample[which(mutation_file$gene == pos_gene_name & mutation_file$effect != "Silent")])
  # old calculation for X_pos
  X_pos<-(expression_file[match(pos_set, row.names(expression_file)) ,])
  # new calculation for X_pos
  X_pos<-expression_file[(row.names(expression_file) %in% pos_set) , ]
  X_pos<-X_pos[!(is.na(rownames(X_pos))) , ]
  Y_pos<-data.frame(matrix(1,length(X_pos[,1]),1))
  colnames(Y_pos)<-"sample"
  rownames(Y_pos)<-rownames(X_pos)
  print(paste("Number of Positive Samples: ", length(X_pos[,1])))
  # load negative training expression data for A
  #neg_set<-unique(mutation_file$sample[!(mutation_file$sample %in% pos_set)])
  neg_set<-unique(mutation_file$sample[which(mutation_file$gene == neg_gene_name & mutation_file$effect != "Silent")])
  
  # old
  #X_neg<-(expression_file[match(neg_set, row.names(expression_file)) ,])
  # new
  X_neg<-expression_file[(row.names(expression_file) %in% neg_set) , ]
  print(paste("Number of Negative Samples: ", length(X_neg[,1])))
  
  X_neg<-X_neg[!(is.na(rownames(X_neg))) , ]
  Y_neg<-data.frame(matrix(0,length(X_neg[,1]),1))
  colnames(Y_neg)<-"sample"
  rownames(Y_neg)<-rownames(X_neg)
  # combine pos and neg
  X<-rbind(X_pos, X_neg)
  Y<-rbind(Y_pos, Y_neg)
  # here we determine
  test<-expression_file[(row.names(expression_file) %in% unique(mutation_file$sample)) , ]
  test<- test[ -(which(row.names(test) %in% row.names(X))) , ]
  print(paste("Number of Test(remaining) Samples: ", length(test[,1])))
  
  # here we set test to all the samples...
  #test<-
  # return variables
  return(list("X" = X, "Y" = Y,"test" = test))
}

#X<-output_pos$matrix1
#Y<-output$Y
#nfolds_input=20


# # runClassifier(rna_XY$X,rna_XY$Y,folds,num_iterations,alpha=1)
# #X<-full_set_XY$X
# #Y<-full_set_XY$Y
# #nfolds_input<-10
# #iterations<-1
# #alpha=1
# runClassifier<- function(X, Y, nfolds_input, iterations, alpha ,weighted, stratified) {
#   # The elastic-net penalty is controlled by alpha, and bridges the gap between lasso (α=1, the default) and ridge (α=0).
#   # lasso results in sparse coefficients
#   # ridge regression results in more coefficients
#   
#   # calculate number of folds
#   # to guarantee that at least one positive/negative label is in each fold, otherwise at least 10 samples per fold
#   # first check if 10 observations exist per fold if 10 folds is chosen, if not set to min number of folds to keep at least 10 observations in each fold
#   if (nrow(Y)/nfolds_input < 10) {
#     nfolds_input<-floor(nrow(Y)/10)
#   }
#   
#   nfolds<-min(max(3,nfolds_input),sum(Y==1),sum(Y==0))
#   print(paste("Performing cross-validation using",nfolds,"folds"))
#   print(paste("Alpha set to:",alpha,"(1=LASSO; fast; less coefficients, 0=Ridge Regression; slow; all coefficients)"))
#   # create weights, weights are a fraction of the pos/neg over the total number of labels
#   fraction_0<-rep(1-sum(Y==0)/nrow(Y),sum(Y==0))
#   fraction_1<-rep(1-sum(Y==1)/nrow(Y),sum(Y==1))
#   # assign 1 - that value to a "weights" vector
#   weights<-numeric(nrow(Y))
#   if (weighted==TRUE) {
#     weights[Y==0]<-fraction_0
#     weights[Y==1]<-fraction_1
#   } else {
#     weights<-rep(1,nrow(Y))
#   }
#   # begin iterations
#   pr_score_list<-as.numeric(list())
#   lambda_list<-as.numeric(list())
#   score_list<-as.numeric(list())
#   precision_list<-as.numeric(list())
#   recall_list<-as.numeric(list())
#   pr<-c()
#   i<-1
#   for (i in 1:iterations) {
#     print(paste("Iteration:",i,"of",iterations))
#     
#     if (stratified==TRUE) {
#       # assign folds evenly using the mod operator
#       fold0 <- sample.int(sum(Y==0)) %% nfolds
#       fold1 <- sample.int(sum(Y==1)) %% nfolds
#       foldid <- numeric(nrow(Y))
#       foldid[Y==0] <- fold0
#       foldid[Y==1] <- fold1
#       foldid <- foldid + 1
#     } else {
#       foldid <- numeric(nrow(Y))
#       foldid <- sample.int(nrow(Y)) %% nfolds
#       foldid <- foldid+1
#     }
#     
#     cv<-cv.glmnet(as.matrix(X), as.factor(Y[,1]), alpha=alpha, foldid = foldid, family = "binomial", type.measure='auc', parallel=TRUE, weights = weights)
#     
#     
#     # INSTEAD I AM OVERWRITTING THAT WITH A LAMBDA MIN FROM ROC MAXIMIZATION
#     min_lambda<-cv$lambda[which(cv$cvm==max(cv$cvm))]
#     
#     # Precision/Recall and PR AUC calculation
#     resp<-data.frame(predict(cv, s=min_lambda, data.matrix(X), type="response"))
#     resp$binary[resp$X1>.5]<-1
#     resp$binary[resp$X1<.5]<-0
#     
#     TP<-length(which(resp$binary==1 & Y[,1]==1))
#     FP<-length(which(resp$binary==1 & Y[,1]==0))
#     FN<-length(which(resp$binary==0 & Y[,1]==1))
#     if ((TP+FP) > 0) { precision<-TP/(TP+FP) } else { precision<-0 }
#     if ((TP+FN) > 0) { recall<-TP/(TP+FN) } else { recall <- 0 }
#     #print(paste("Precision:",precision))
#     #print(paste("Recall:",recall))
#     
#     #pos_resp
#     if (abs(max(as.numeric(resp$X1)) - min(as.numeric(resp$X1))) == 0) {
#       pr$auc.integral<-0
#       pr$auc.davis.goadrich<-0
#       pr_auc<-0
#     }
#     else {
#       pr <- pr.curve( scores.class0=as.numeric(resp$binary),weights.class0=as.numeric(Y[,1]),curve=FALSE)
#       pr_auc<-pr$auc.integral
#       #print(paste("PR AUC:",pr_auc))
#     }
#     
#     
#     precision_list<-c(precision_list,precision)
#     recall_list<-c(recall_list,recall)
#     
#     pr_score_list<-c(pr_score_list, pr_auc)
#     score_list<-c(score_list, max(cv$cvm))
#     lambda_list<-c(lambda_list, min_lambda)
#   }
#   return(list("cv" = cv, 
#               "lambda.min" = mean(lambda_list),
#               "score_list" = (score_list),
#               "score.var" = var(score_list),
#               "pr_auc" = (pr_score_list),
#               "recall_list" = (recall_list),
#               "precision_list" = (precision_list)))
# }
# 
# 
# #X<-full_set_XY$X
# #Y<-full_set_XY$Y
# #nfolds_input<-10
# #iterations<-5
# #alpha=1
# #weighted=TRUE
# #stratified=TRUE



runClassifier_V2<- function(X, Y, nfolds_input, iterations, alpha, weighted, stratified) {
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
  # create weights, weights are a fraction of the pos/neg over the total number of labels
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
  
  # create an initial model and get a lambda
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
    #set.seed(123)
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
      
      #loop contents here
      
      
      #for (fold in 0:(nfolds-1)) {
      model<-glmnet(as.matrix(X[which(foldid!=fold) , ]), lambda=lambda_start, as.factor(Y[which(foldid!=fold),1]), family = "binomial", weights = weights[which(foldid!=fold)])
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
          pr <- pr.curve( scores.class0=as.numeric(resp_column),weights.class0=as.numeric(Y[which(foldid==fold),1]),curve=FALSE)
          pr_auc<-pr$auc.integral
          #print(paste("PR AUC:",pr_auc))
        }
        pr_auc_list[col]<-pr_auc
        
        # Calculate Precision and Recall
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
      
      #plot(log(model$lambda), pr_auc_list)
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
  
  
  
} # end of Run_Classifier_PR




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



# output<-run_quick_gsea("All_Prom_UTR_enh_cds.tsv", "LGG", resp_matrix,10, 4)
#file<-"All_Prom_UTR_enh_cds.tsv"
#subtype_cde<-"LGG"
#num_events<-10
#min_events_per_sample<-4
run_quick_gsea<-function(file, subtype_cde, resp_matrix, num_events, min_events_per_sample, p_value_threshold, q_value_threshold) {
  # initialize output variables
  gene_name<-"ALL"
  neg_events<-data.frame()
  pos_events<-data.frame()
  pos_event_sample_names<-data.frame()
  neg_event_sample_names<-data.frame()
  
  
  # copy mutation element list (prey events)
  system(paste(sep="","cp ",PANCAN_DATA,file," ",DD_HOME,gene_name,subtype_cde,"geneset.gmt"))
  
  
  resp_ordered<-resp_matrix[order(resp_matrix, decreasing=TRUE)]
  write.table(resp_ordered, file=paste(DD_HOME,"rankedfile_",subtype_cde,"_",gene_name,".rnk",sep=""), quote=FALSE, sep="\t", col.names=FALSE)
  
  print("Running GSEA...")
  # generate random number for output file directory
  rand<-sample(1:100000000,1)
  print(paste("Using Random Directory",rand))  
  
  # options
  # -rnd_type equalize_and_balance 
  gsea_cmd<-paste("java -cp ",DD_HOME,"gsea2-2.2.2.jar -Xmx15000m xtools.gsea.GseaPreranked -gmx ",DD_HOME,gene_name,subtype_cde,"geneset.gmt -rnd_type equalize_and_balance -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ",DD_HOME,"rankedfile_",subtype_cde,"_",gene_name,".rnk -scoring_scheme weighted_p2 -rpt_label my_analysis -include_only_symbols true -make_sets true -plot_top_x ",num_events+25," -rnd_seed 149 -set_max 500 -set_min ",min_events_per_sample, "-zip_report false -out ",DD_HOME,"GSEA_reports/",rand," -gui false > ",OUTPUT_DATA,"/GSEA_logfile_",gene_name,"_",subtype_cde,".log", sep="")
  
  
  
  return_code<-system(gsea_cmd, intern=FALSE)
  print(return_code)
  if (return_code > 0) {
    print("Java failure!!!!")
    return(list("pos_events" = pos_events, "neg_events" = neg_events))
  }
  setwd(paste(DD_HOME,"GSEA_reports/",rand,sep=""))
  dir_results<-list.files()
  dir_results[length(dir_results)]
  setwd(paste(DD_HOME,"GSEA_reports/",rand,"/",dir_results[length(dir_results)],sep=""))
  
  # positive evrand
  GSEA_output<-read.csv(list.files(pattern=".*gsea_report_for_na_pos_.*xls.*")   , sep="\t", row.names=1, stringsAsFactors = FALSE)
  GSEA_output<-GSEA_output[ (GSEA_output$NOM.p.val < p_value_threshold) , ]
  # I removed the FDR calculation from GSEA as a lot of my events did had high FDR rates
  GSEA_output<-GSEA_output[ (GSEA_output$FDR.q.val < q_value_threshold) , ]
  pos_events<-data.frame()
  if (length(GSEA_output[,1]) == 0) {
    print("NO POS EVENTS FOUND!")
  } else {
    pos_events<-GSEA_output   
    print(paste(length(pos_events[,1]), " positive events found!"))
    
    # take top num_events (passed from wrapper script) events, make sure the GSEA setting to build the xls files is at 10 or more...
    pos_events<-pos_events[1:min(length(pos_events[,1]), num_events) , ]
    print(pos_events[,1])
    # get sample names for each event
    event<-data.frame()
    events<-data.frame()
    for (i in 1:length(rownames(pos_events))) {
      filename<-paste(rownames(pos_events)[i],".xls", sep="")
      event_output<-read.csv(filename, sep="\t", row.names=1, stringsAsFactors = FALSE)
      events<-rbind.fill(events, data.frame(cbind(rownames(pos_events)[i], t(event_output$PROBE))))
    }
    rownames(events)<-events[,1]
    events<-events[,-1]
    pos_event_sample_names<-t(events)
    
    
    
    
    
    
    
  }
  
  # NEGATIVE EVENTS
  GSEA_output_neg<-read.csv(list.files(pattern=".*gsea_report_for_na_neg_.*xls.*")   , sep="\t", row.names=1, stringsAsFactors = FALSE)
  GSEA_output_neg<-GSEA_output_neg[ (GSEA_output_neg$NOM.p.val < p_value_threshold) , ]
  GSEA_output_neg<-GSEA_output_neg[ (GSEA_output_neg$FDR.q.val < q_value_threshold) , ]
  neg_events<-data.frame()
  neg_final_hello<-matrix(ncol=6,nrow=0)
  neg_final_hello_data<-matrix(ncol=5,nrow=0)
  if (length(GSEA_output_neg[,1]) == 0) {
    print("NO NEG EVENTS FOUND!")
  } else {
    neg_events<-GSEA_output_neg   
    print(paste(length(neg_events[,1]), "negative events found!"))
    
    # take top 10 events, make sure the GSEA setting to build the xls file is at 10 or more...
    neg_events<-neg_events[1:min(length(neg_events[,1]), num_events) , ]
    print(neg_events[,1])
    # get sample names for each event
    event<-data.frame()
    events<-data.frame()
    for (i in 1:length(rownames(neg_events))) {
      filename<-paste(rownames(neg_events)[i],".xls", sep="")
      event_output<-read.csv(filename, sep="\t", row.names=1, stringsAsFactors = FALSE)
      events<-rbind.fill(events, data.frame(cbind(rownames(neg_events)[i], t(event_output$PROBE))))
    }
    rownames(events)<-events[,1]
    events<-events[,-1]
    neg_event_sample_names<-t(events)
    
  }
  
  return(list("pos_events" = pos_events, "neg_events" = neg_events, "pos_event_samples" = pos_event_sample_names, "neg_event_samples" = neg_event_sample_names))
}


#gmt_file<-"splicing_enhancers.gmt"
#num_permutations<-5
#num_events<-5
##pvalue_threshold<-.05
#min_gene_set_size<-5

# positive_recursive_iterate(bait=paste0("Initial_",gene), pred_resp, tree_path, gmt_file, full_set_XY$X, "origY"=full_set_XY$Y, orig_score_vector=big_model$pr_auc, num_permutations, num_events)
# bait=paste0("Initial_",gene)
# X<-full_set_XY$X
# Y<-full_set_XY$Y
# orig_score_vector=big_model$pr_auc
# gsea_pvalue_threshold<-.05
# gsea_fdr_threshold<-.05
# #orig_score_vector<-c(.5,.5,.5,.5,.5)
# pvalue_threshold<-.05
# min_gene_set_size<-5
# # test run
# #results<-run_gsea_V2(gmt_file, X, Y, pred_resp, orig_score_vector, num_permutations, num_events, pvalue_threshold, min_gene_set_size,gsea_pvalue_threshold, gsea_fdr_threshold)

run_gsea_V2<-function(bait, gmt_file, X, Y, pred_resp, orig_score_vector, num_permutations, num_events, pvalue_threshold, min_gene_set_size,
                      gsea_pvalue_threshold, gsea_fdr_threshold, percent_overlap, folds, enrichment_analysis_only) {
  print(paste("Bait:",bait))
  print(paste("LengthX:",length(X[,1])))
  print(paste("LengthY:",length(Y[,1])))
  #print(paste("Y SAMPLES:",rownames(Y)))
  
  print(paste("Percent Overlap",percent_overlap))
  print(paste("gmt file:",gmt_file))
  print(paste("num_permutations:",num_permutations))
  print(paste("num_events:",num_events))
  print(paste("pvalue_threshold:",pvalue_threshold))
  print(paste("min_gene_set_size:",min_gene_set_size))
  rand<-sample(1:100000000,1)
  print(paste("Using Rand Num:",rand))
  print(paste("Original Scores:",orig_score_vector))
  
  # initialize output variable
  result_output<-data.frame(matrix(ncol=9,nrow=0))
  models_output<-c()
  #result_output<-list("df_results"=results,"models"=models_output)
  
  # running GSEA on the pred_resp, note this is the "test" dataset, everything in Y that is not 1
  pred_resp<-pred_resp[order(pred_resp, decreasing=TRUE)]
  write.table(pred_resp, file=paste(DD_HOME,"rankedfile_",rand,".rnk",sep=""), quote=FALSE, sep="\t", col.names=FALSE)
  print("Running GSEA...")
  gsea_cmd<-paste("java -cp ",DD_HOME,"gsea2-2.2.2.jar -Xmx15000m xtools.gsea.GseaPreranked -gmx ",DD_HOME,gmt_file," -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ",DD_HOME,"rankedfile_",rand,".rnk -scoring_scheme weighted -rpt_label my_analysis -include_only_symbols true -make_sets true -plot_top_x 200 -rnd_seed timestamp -set_max 500 -set_min ",min_gene_set_size," -zip_report false -out ",DD_HOME,"GSEA_reports/",rand," -gui false > ",OUTPUT_DATA,"/GSEA_logfile_",rand,".log", sep="")
  
  return_code<-system(paste0(gsea_cmd, " 2> ",DD_HOME,"logfile",rand,".txt"))
  print(return_code)
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
  
  setwd(paste(DD_HOME,"GSEA_reports/",rand,sep=""))
  dir_results<-list.files()
  dir_results[length(dir_results)]
  setwd(paste(DD_HOME,"GSEA_reports/",rand,"/",dir_results[length(dir_results)],sep=""))
  
  GSEA_output<-read.csv(list.files(pattern=".*gsea_report_for_na_pos_.*xls.*")   , sep="\t", row.names=1, stringsAsFactors = FALSE)
  GSEA_output<-GSEA_output[ (GSEA_output$NOM.p.val < gsea_pvalue_threshold) , ]
  GSEA_output<-GSEA_output[ (GSEA_output$FDR.q.val < gsea_fdr_threshold) , ]
  pos_events<-data.frame()
  # check to see if enrichment_analysis equals true, if so return just the GSEA analysis
  if (enrichment_analysis_only==TRUE)
    return(GSEA_output)
  
  # running regular LURE
  if (length(GSEA_output[,1]) == 0) {
    print("NO POSITIVE EVENTS FOUND!")
    return(list("df_results"=result_output,"models"=models_output))
  }
  # events found, proceeding with dumpster diver
  print(paste(length(GSEA_output[,1]), "positive event(s) found!"))
  
  
  # take top 50 events events, make sure the GSEA setting to build the xls files is at 50 or more...
  GSEA_output<-GSEA_output[1:min(nrow(GSEA_output), 50) , ]
  print(row.names(GSEA_output))
  # get sample names for each event
  event_sample_names<-data.frame(matrix(nrow=0,ncol=2),stringsAsFactors = FALSE)
  i<-1
  for (i in 1:length(rownames(GSEA_output))) {
    filename<-paste(rownames(GSEA_output)[i],".xls", sep="")
    file_data<-read.csv(filename, sep="\t", row.names=1, stringsAsFactors = FALSE)
    sample_names<-paste(sort(file_data$PROBE),collapse = " ")
    event_sample_names<-rbind(event_sample_names,data.frame("event"=rownames(GSEA_output)[i], "sampes"=sample_names,stringsAsFactors = FALSE))
  }
  
  # need to add something here to remove duplicate events (really only need it in CNV data)
  #duplicated(event_sample_names$sampes)
  
  # Check for overlap between new event samples and pos_set. percent_overlap is passed in as an argument. 
  # If the percent of samples of the new event in the old pos_set is greater than the percent_overlap, we do not consider that event
  # to start I am putting percent overlap at only 33%.  
  # So if 33% of the mutant samples are in the original pos_set have the mutation under consideration, we skip it.
  no_col<-max(count.fields(paste0(DD_HOME,gmt_file)))
  gmt_mutation_data<-(read.csv(paste0(DD_HOME,gmt_file),sep="\t",col.names=1:no_col,header = FALSE,stringsAsFactors = FALSE))
  
  new_event_sample_names<-data.frame(matrix(nrow=0,ncol=2),stringsAsFactors = FALSE)
  colnames(new_event_sample_names)<-colnames(event_sample_names)
  #row<-2
  
  for (row in 1:nrow(event_sample_names)) {
    print(paste("Processing:",event_sample_names$event[row]))
    orig_pos_set<-rownames(Y)[Y==1]
    sample_names<-strsplit(event_sample_names$sampes[row]," ")[[1]]
    all_mutant_samples<-gmt_mutation_data[which(toupper(gmt_mutation_data$X1) == toupper(event_sample_names$event[row])),]
    all_mutant_samples<-c(unlist(all_mutant_samples[,-c(1,2)]))
    
    num_mutants_in_orig_pos_set<-length(which(all_mutant_samples %in% orig_pos_set))
    percent_of_mutants_in_orig_pos_set<-num_mutants_in_orig_pos_set/(length(sample_names)+num_mutants_in_orig_pos_set)
    if (percent_of_mutants_in_orig_pos_set > percent_overlap) {
      print("Too many mutants in the orig positive set...too much overlap...skipping event")
      next()
    }
    new_event_sample_names<-rbind(new_event_sample_names,event_sample_names[row,])
  }
  # move our new set of events back into the original array and limit the event size to num_events
  event_sample_names<-new_event_sample_names[1:min(nrow(new_event_sample_names),num_events) , ]
  print(paste("Now considering",nrow(event_sample_names),"events..."))
  
  #row<-1
  
  for (row in 1:nrow(event_sample_names)) {
    if (is.na(event_sample_names$event[row])) {
      print("No events left...")
      next()
    }
    print(paste("Processing:",event_sample_names$event[row]))
    print(paste("old_sample_names",rownames((Y[Y==1]))))
    sample_names<-strsplit(event_sample_names$sampes[row]," ")[[1]]
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
      
      random_set<-c(rownames(Y)[which(Y==1)],sample(rownames(Y)[which(Y==0)],length(sample_names)))
      random_XY<-createXY_pos("test","nada",X,pos_set=random_set)
      
      random_model<-runClassifier_V2(random_XY$X,random_XY$Y,folds,1,alpha=1, weighted=TRUE,stratified = TRUE)
      print(paste("RANDOM PR AUC:",mean(random_model$pr_auc)))
      random_scores<-c(random_scores,mean(random_model$pr_auc))
    }
    # check scores
    print("Variance Actual:",var(actual_scores))
    print("Variance Orig:",var(orig_score_vector))
    print("Variance Random:",var(random_scores))
    
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
    
    
    #saveRDS(new_XY$Y, paste0(PANCAN_DATA,"Labels_",bait,"_",event_sample_names$event[row],".rds"))
    #saveRDS(resp_output, paste0(PANCAN_DATA,"Scores_",bait,"_",event_sample_names$event[row],".rds"))
    
    
    
    result_output<-rbind(result_output,data.frame("bait"=bait,
                                                  "event"=event_sample_names$event[row],
                                                  "sample_names"=event_sample_names$sampes[row],
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




run_t_test<-function(A, B) {
  myt_test <- try(t.test(A,B))
  if (inherits(myt_test, "try-error"))
  {
    cat(myt_test)
    myt_test$p.value <- 1
  }
  return(myt_test)
}



# Bonferroni FDR calculation
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


# Bonferroni FDR calculation
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


# tree_path is an internal variable that keeps track of each individual path.  
#   : Consider it a stack variable that is pushed new leaves to go down the tree and pops off leaves to back up the tree
#   : Set it initially to the starting oncogene
# tree_path_list is an external variable which is a list of the tree paths.  At the bottom of each tree, the current tree_path
# TEST RUN
#tree_path_list<-c()
#positive_recursive_iterate(pred_resp, tree_path, gmt_file, X, "origY"=Y, orig_score_vector, num_permutations)

#origY<-full_set_XY$Y
#X<-full_set_XY$X
#orig_score_vector=rep(0,num_permutations)
#permutations<-5
#bait="ATRX_truncating"
positive_recursive_iterate<-function(bait, pred_resp, tree_path, gmt_file, X, origY, orig_score_vector, num_permutations, num_events, percent_overlap, pvalue_threshold, min_gene_set_size, gsea_pvalue_threshold, gsea_fdr_threshold,tree_length,max_tree_length, folds, enrichment_analysis_only) {
  print(paste0("Tree Length:",tree_length))
  results<-run_gsea_V2(bait, gmt_file, X, origY, pred_resp, orig_score_vector=orig_score_vector, num_permutations, num_events=num_events, pvalue_threshold=pvalue_threshold, min_gene_set_size=min_gene_set_size,gsea_pvalue_threshold=gsea_pvalue_threshold, gsea_fdr_threshold=gsea_fdr_threshold, percent_overlap, folds = folds, enrichment_analysis_only=enrichment_analysis_only)
  
  if (enrichment_analysis_only==TRUE) {
    print("Enrichment Analysis Only...Returning results now")
    return(results)
  }
  #iteration_events<-run_gsea(gmt_file, "pos",X,origY,test, "LGG", pos_resp, .8, 10, num_iterations, 25, .05)  
  # OUTPUT FROM RUN GSEA:
  #result_output<-rbind(result_output,data.frame("bait"=bait
  #                                            "event"=event_sample_names$event[row],
  #                                              "sample_names"=event_sample_names$sampes[row],
  #                                              "actual_vs_random_pvalue"=actual_vs_random_pvalue,
  #                                              "actual_vs_orig_pvalue"=actual_vs_orig_pvalue,
  #                                              "pr_auc_score"=paste(actual_model$pr_auc,collapse=" "),
  #                                              "mean_lambda"=actual_model$lambda.min,
  #                                              "pvalue_over_random"=actual_vs_random$p.value,
  #                                              "pvalue_over_orig"=actual_vs_orig$p.value,
  
  
  
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
      
      #positive_recursive_iterate<-function(pred_resp, tree_path, gmt_file, X, origY, orig_score_vector, num_permutations) {
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
    } # end of for loop
  } else {
    # base case
    print("Bottom of Tree or Max Tree Length achieved")
    print(tree_path)
    assign("tree_path_list", c(tree_path_list, list(tree_path)), envir = .GlobalEnv)
  }
  return(tree_path_list)
}


# best_classifier_iterate<-function(bait, pred_resp, tree_path, gmt_file, X, origY, orig_score_vector, num_permutations, num_events, percent_overlap, pvalue_threshold, min_gene_set_size, gsea_pvalue_threshold, gsea_fdr_threshold) {
#   event_list<-data.frame(matrix(nrow=0,ncol=2))
#   event_list<-rbind(event_list,data.frame("event"=bait,"pr_auc"=mean(orig_score_vector)))
#   
#   #results<-run_gsea_V2(gmt_file, X, Y, pred_resp, orig_score_vector, num_permutations, num_events, pvalue_threshold, min_gene_set_size,gsea_pvalue_threshold, gsea_fdr_threshold)
#   results<-run_gsea_V2(bait, gmt_file, X, origY, pred_resp, orig_score_vector=orig_score_vector, num_permutations, num_events=num_events, pvalue_threshold=pvalue_threshold, min_gene_set_size=min_gene_set_size,gsea_pvalue_threshold=gsea_pvalue_threshold, gsea_fdr_threshold=gsea_fdr_threshold, percent_overlap)
#   #iteration_events<-run_gsea(gmt_file, "pos",X,origY,test, "LGG", pos_resp, .8, 10, num_iterations, 25, .05)  
#   # OUTPUT FROM RUN GSEA:
#   #result_output<-rbind(result_output,data.frame("bait"=bait
#   #                                            "event"=event_sample_names$event[row],
#   #                                              "sample_names"=event_sample_names$sampes[row],
#   #                                              "actual_vs_random_pvalue"=actual_vs_random_pvalue,
#   #                                              "actual_vs_orig_pvalue"=actual_vs_orig_pvalue,
#   #                                              "pr_auc_score"=paste(actual_model$pr_auc,collapse=" "),
#   #                                              "mean_lambda"=actual_model$lambda.min,
#   #                                              "pvalue_over_random"=actual_vs_random$p.value,
#   #                                              "pvalue_over_orig"=actual_vs_orig$p.value,
#   if (nrow(results$df_results) == 0)
#     return(event_list)
#   tree_length<-1
#   while (nrow(results$df_results) > 0) {
#     print(paste0("Tree Length: ",tree_length))
#     tree_length<-tree_length+1
#     df_results<-results$df_results
#     best_event<-which(min(results$df_results$actual_vs_orig_pvalue) == results$df_results$actual_vs_orig_pvalue)
#     gene<-results$df_results$event[best_event]
#     print(paste("Adding samples containing event",gene,"to the positive set to recreate a classifier and rerun GSEA and look for new events:"))
#     # add event samples to 1 in Y (add event samples to the positive set)
#     # must create a new Y, otherwise it gets changed
#     newY<-origY
#     event_sample_names_inner_iteration<-strsplit(df_results$sample_names[best_event]," ")[[1]]
#     pr_auc_scores<-as.numeric(strsplit(df_results$pr_auc_score[best_event]," ")[[1]])
#     newY[which(rownames(newY) %in% event_sample_names_inner_iteration) , ]<-1
#     
#     print(paste("Length of Old Positive Sample Names:",length(origY[origY==1])))
#     print(paste("Length of New Positive Sample Names:",length(event_sample_names_inner_iteration)))
#     print(paste("Length of New and Old Positive Sample Names:",length(newY[newY==1])))
#     pos_set<-rownames(newY)[which(newY==1)]
#     # create a new test not containing the event samples
#     new_XY<-createXY_pos("test","nada",X,pos_set=pos_set)
#     
#     big_model<-runClassifier_V2(X,newY,folds,num_permutations,alpha=1,weighted = TRUE,stratified = TRUE)
#     # creating new response vector with all the samples without known event samples
#     new_pred_resp<-((predict(big_model$model, s=min((big_model$min_lambda)[which(big_model$pr_auc==max(big_model$pr_auc))]), newx=data.matrix(new_XY$test), type="response")))[,1]
#     
#     
#     print(paste("Length of test samples:",length(new_pred_resp)))
#     # new score vector is the previous(orig) classifier scores
#     new_score_vector<-as.numeric(strsplit(df_results$pr_auc_score[best_event]," ")[[1]])
#     
#     # add to event list
#     event_list<-rbind(event_list,data.frame("event"=gene,"pr_auc"=mean(new_score_vector)))
#     
# 
#     # output to source-target spreadsheet
#     df_output<-data.frame("source"=df_results$bait[best_event],
#                           "interaction"="positive",
#                           "target"=gene,
#                           "pvalue_over_orig"=df_results$pvalue_over_orig[best_event],
#                           "pr_auc"=mean(pr_auc_scores))
#     assign("STG_output", rbind(STG_output,df_output), envir = .GlobalEnv)
#     tempY<-newY
#     results<-run_gsea_V2(bait=gene, gmt_file, X, newY, pred_resp=new_pred_resp, orig_score_vector=new_score_vector, num_permutations, num_events=num_events, pvalue_threshold=pvalue_threshold, min_gene_set_size=min_gene_set_size,gsea_pvalue_threshold=gsea_pvalue_threshold, gsea_fdr_threshold=gsea_fdr_threshold, percent_overlap)
#     origY<-tempY
#   }
#     
# 
#     
#     
# return(event_list)
# }




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






