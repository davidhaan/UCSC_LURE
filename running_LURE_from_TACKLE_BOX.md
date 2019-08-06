### What is Tackle Box?

Tackle Box is the Docker image for running LURE, a machine learning method developed by David Haan. The LURE paper is available on BioRxiv at [https://www.biorxiv.org/](https://www.biorxiv.org/).


### What is the directory structure for running LURE with Tackle Box?

To run LURE from Tackle Box, you'll need to have the following directory structure in your working directory:

```
.
├── GSEA_reports
├── LURE_functions.R
├── LURE_positive_controls.R
├── LURE_wrapper.R
├── gsea2-2.2.2.jar
├── input
│   ├── pancan_RNAexp_UVM
│   ├── positive_control_IDH1_missense.gmt
│   ├── positive_control_SF3B1_missense.gmt
│   └── readme
├── oncoprint.py
```

- The LURE scripts (`LURE_*.R`, `oncoprint.py`, and `gsea2-2.2.2.jar`) are in the outermost directory, `./`.
- The input data files are in the `input` directory.
- GSEA, `gsea2-2.2.2.jar`, will write output to `GSEA_reports`.


### How can I run the example LURE job using Tackle Box?

1. Get an interactive session in the Docker container with:
```
docker run -ti --entrypoint /bin/bash -v `pwd`:/data stuartlab/tackle_box
```


2. Go to `/data` in the Docker container. This is the working directory for the LURE job.


3. Start R.


4. In R, do `source("LURE_positive_controls.R")`. This will load some R libraries, load some input files, and then run the LURE analysis pipeline.


5. (Check the results.)
