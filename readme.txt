Python/R code files

- readme.txt

- Source_functions_other: contains eight source functions for FPCA-gY, FPCA-gX, FGLasso, and PSKL:
  A.prellim.func.R, ADMM.new.R, B.prellim.func.R, FPCA.score.R, ProxAlg_FGM.R, auc.R, base.func.R, and prec.rec.R.
  All functions are downloaded from: https://github.com/PercyZhai/FGM_Neighborhood

- Source_functions_npFGM: contains two source functions for the npFGM method:
  main1.R and main2.R

- Simulation:
  Folders A, B, C, and D contain simulation codes for Models A, B, C, and D, respectively.
  Within each folder, there are three subfolders corresponding to different dimensions.
  For example:
  - p30 → p = 30
  - p50 → p = 50
  - p100 → p = 100

- Within sub-folder of Model A (JCGS_Code → Simulation → A → p30), there are four files:
  classifier.sh, classifier.py, npFGM.R, and data_generator.R.

- Within each subfolder of Models B, C, and D (e.g., JCGS_Code → Simulation → B → p30), there are five files:
  classifier.sh, classifier.py, npFGM.R, other.R, and data_generator.R.

- Real_data: contains five folders including real data (ADHD), source files, R/Python scripts, and analysis results.
  - data: time.series.ADHD.Rdata is the raw ADHD dataset.
  - npFGM: R script run_ADHD.R is used for data analysis of npFGM method.
    The Results sub-folder contains:
    - ADHD1_npFGM.Rdata for ADHD group
    - ADHD2_npFGM.Rdata for control group
  - fgDNN: contains two Python scripts (fgDNN_ADHD.sh and fgDNN_ADHD.py) and two sub-folders (get_score and results).
    Outputs are saved in the results folder.
  - Other_alternatives: Results sub-folder contains six .Rdata files (results of FPCA-gY, FPCA-gX, FGLasso, and PSKL methods).
    Results Analysis sub-folder contains three R scripts for analyzing these results.
    All files are downloaded from: https://github.com/PercyZhai/FGM_Neighborhood
  - plot: Run plot.R for final plots and analysis results. AAL.info.R contains coordinate information for the atlas.

---

How to run the simulation code:
1. Save all source files from Source_functions_npFGM and Source_functions_other into your working directory.
2. For each data-generating setting, go to the corresponding sub-folder in Simulation. Example: JCGS_Code → Simulation → A → p30 and run data_generator.R for Model A with p = 30.
3. For fgDNN method: run classifier.py (Python script). The bash script classifier.sh can be used to submit jobs on HPC.
4. For npFGM method: go to sub-folder npFGM and run npFGM.R.
5. For Model A, code for alternative methods (FPCA-gY, FPCA-gX, FGLasso, and PSKL) is omitted; results can be taken from: https://github.com/PercyZhai/FGM_Neighborhood.
6. For Models B, C, and D: run other.R for FPCA-gY, FPCA-gX, FGLasso, and PSKL methods.

---

How to run the real data analysis code:
1. Go to folder Real_data.
2. Save time.series.ADHD.Rdata (from data) and all source files from Source_functions_npFGM and Source_functions_other into your working directory.
3. For npFGM method: go to npFGM and run run_ADHD.R.
4. For FPCA-gY, FPCA-gX, FGLasso, and PSKL methods: results are in six .Rdata files in Other_alternatives → Results. Analysis scripts are in Other_alternatives → Results Analysis.
5. For fgDNN method:
   a. Run ADHD_score.R in get_score to extract scores (saved as ADHD_1.csv and ADHD_2.csv).
   b. Run fgDNN_ADHD.sh to call fgDNN_ADHD.py to obtain results.
6. For final plots: go to plot and run plot.R. AAL.info.R contains atlas coordinate information.
