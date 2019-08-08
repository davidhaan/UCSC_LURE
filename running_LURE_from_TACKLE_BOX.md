### What is LURE?

LURE is a semi-supervised machine learning tool, Learning UnRealized Events (LURE), that finds new drivers sharing molecular signatures with known drivers. The LURE paper is available on [BioRxiv](https://www.biorxiv.org/content/10.1101/727891v1). The citation is:


    Using Transcriptional Signatures to Find Cancer Drivers with LURE
    David Haan, Ruikang Tao, Verena Friedl, Ioannis Nikolaos Anastopoulos, Christopher K Wong, Alana S Weinstein, Joshua M Stuart
    bioRxiv 727891; doi: https://doi.org/10.1101/727891


The LURE website is [here](https://sysbiowiki.soe.ucsc.edu/lure).

### LURE Scripts

LURE_positive_controls.R: 

Runs the positive controls for the LURE manuscript. There are two positive controls: SF3B1 in Uveal Melanoma (UVM) and IDH1 in Lower Grade Gliomas (LGG). The positive controls are further described in the manuscript.

LURE_functions.R: 

This script contains R functions associated with the LURE method. LURE can be ran as a function in R or using the LURE_wrapper.R script.

LURE_wrapper.R: 

A wrapper command line script which calls LURE using command line arguments. The command line arguments are descibed in the "LURE Parameters" section below.


### What is Tackle Box?

Tackle Box is the Docker image for running LURE. While it is possible to run the LURE code from directly within a local installation of R, we recommend running LURE from within Tackle Box.

---
### What is the directory structure for running LURE with Tackle Box?

To run LURE from Tackle Box, you'll need to have a specific directory structure in your working directory. The structure of the `/example_jobs` directory is shown below.


    /example_jobs
    ├── input
    │   ├── pancan_RNAexp_LGG
    │   ├── pancan_RNAexp_UVM
    │   ├── positive_control_IDH1_missense.gmt
    │   └── positive_control_SF3B1_missense.gmt
    ├── output
    └── temp
        └── GSEA_reports


- LURE reads the input files from `input`. The `/example_jobs/input` directory has the data required for two jobs. Two files are required for each LURE job.
  1. expression matrix
  2. gmt file where the set name is a feature and the set members are samples that are positive for the named feature

- LURE output will be saved to the `output` directory.

- LURE writes temporary files and GSEA reports in the `temp` directory.

---
### How can I run the example LURE job using Tackle Box?

1. Get an interactive session in the Docker container with:
    ```
    docker run -ti --entrypoint /bin/bash -v `pwd`:/data stuartlab/tackle_box
    ```


2. Go to `/example_jobs` in the Docker container. This is the working directory for the example LURE jobs.


3. Start R.


4. In R, do `source("/lure_scripts/LURE_positive_controls.R")`. This will load some R libraries, load some input files, and then run the LURE analysis pipeline for two jobs, UVM positive controls and LGG positive controls. The UVM job should complete in a few minutes. The LGG job is bigger, so it will take more time.

5. Exit R with `q()`.

5. Move the results to `/data` directory with `cp -r /example_jobs /data/.`. This will make the files available for review outside of Tackle Box.

6. Exit the Docker container with `exit`.

---
### What are all of those output files in the `example_jobs` directory?

This is the output from running the UVM positive controls example job:

    output
    ├── Pos_Ctrl_Initial_Classifier_Scores_10_positive_0.71_AUC_UVM_SF3B1-SET1_missense.pdf
    ├── Pos_Ctrl_TCGA_UVM_SF3B1-SET1_missense_Oncoprint_Plot_2019_08_08_03_46_positive_control_SF3B1_missense.pdf
    ├── Pos_Ctrl_TCGA_UVM_SF3B1-SET1_missense_Original_PR_AUC_scores2019_08_08_03_46_positive_control_SF3B1_missense.tsv
    ├── Pos_Ctrl_TCGA_UVM_SF3B1-SET1_missense_Set_Cover_mutation_matrix_2019_08_08_03_46_positive_control_SF3B1_missense.tsv
    ├── Pos_Ctrl_TCGA_UVM_SF3B1-SET1_missense_set_cover_PR_AUC_scores_2019_08_08_03_46_positive_control_SF3B1_missense.tsv
    ├── Pos_Ctrl_TCGA_UVM_SF3B1-SET1_missense_set_cover_classifier_scores_2019_08_08_03_46_positive_control_SF3B1_missense.tsv
    ├── Pos_Ctrl_TCGA_UVM_SF3B1-SET1_missense_set_cover_solution_2019_08_08_03_46_positive_control_SF3B1_missense.tsv
    ├── UVM_SF3B1-SET1_missense_bipartite.csv
    └── UVM_SF3B1-SET1_missense_bipartite.pdf

- `Pos_Ctrl_TCGA_UVM_SF3B1-SET1_missense_Oncoprint_Plot_2019_08_08_03_46_positive_control_SF3B1_missense.pdf` is an oncoprint figure that shows which samples are positive for the bait and catch events. The samples are ordered by the classifier score, the top row in the figure.

- `UVM_SF3B1-SET1_missense_bipartite.pdf` shows a bipartite graph representing the result of the set coverage step in LURE. The edges of the graph connect events with samples that have the event.

- `Pos_Ctrl_Initial_Classifier_Scores_10_positive_0.71_AUC_UVM_SF3B1-SET1_missense.pdf` is a waterfall plot of each sample's LURE bait score.

### LURE Parameters

--folds<br />
Default: 10<br />
Description: Number of Cross Validation Folds<br />

--num_permutations<br />
Default: 5<br />
Description: Number of Permutations/Iterations<br />

--min_gene_set_size<br />
Default: 4<br />
Description: prey event minimum size: parameter for GSEA, only events with n or more mutated samples are considered<br />

--percent_overlap<br />
Default: 0.5<br />
Description: If <percent_overlap> of the new event samples are in the existing positive set already then we skip it. A smaller number is more restrictive<br />

--max_tree_length<br />
Default: 5<br />
Description: Max Tree Length: Here we set the max length of the Event Discovery Tree<br />

--bait_gene<br />
Default: "custom_to_test.tsv"<br />
Description: .tsv file containing bait genes to test<br />

--gmt_file<br />
Default: "mutation_specific_cytoband_fusion_5_15_2019.gmt"<br />
Description: .gmt file containing catch data<br />

--gsea_fdr_threshold<br />
Default: .25<br />
Description: FDR value threshold for GSEA step<br />

--gsea_pvalue_threshold<br />
Default: .05<br />
Description: P value threshold for GSEA step<br />

--LURE_pvalue_threshold<br />
Default: .05<br />
Description: P value threshold for LURE AUC score step<br />

--max_num_events<br />
Default: 5<br />
Description: Used to limit the number of prey events found by GSEA and considered for LURE's classifier AUC score step.  The events are sorted by GSEA NES score so the top events will be chosen. The larger this parameter the longer the runtime<br />

--feature_data_file<br />
Default: "../input/pancan_RNAexp_UVM"<br />
Description: Feature Data File (Must be located in input directory)<br />

--target_gmt_file<br />
Default: ""<br />
Description: This argument only pertains when LURE is run with enrichment only.  It is the gmt file for the test/target dataset.  The original gmt_file argument is for the bait<br />

--target_feature_file<br />
Default: ""<br />
Description: This argument only pertains when LURE is run with enrichment only.  It is the gmt file for the test/target dataset.  The original gmt_file argument is for the bait<br />

--output_file_prefix<br />
Default="V10"<br />
Description: This is the file prefix assigned to all the output files.  For multiple runs it helps keep track of each run.<br />

--tissue<br />
Default=""<br />
Description: Tissue or Tumor Type, used for additional filename prefix for large pancan runs<br />
