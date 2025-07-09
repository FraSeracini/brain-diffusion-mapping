# COMP0118 – Coursework 1: Parametric Models

## Overview

This project was completed as part of the COMP0118 module *Computational Modelling for Biomedical Imaging* at UCL. The coursework focuses on implementing and analyzing parametric models such as the linear diffusion tensor and the ball-and-stick model for processing and interpreting diffusion MRI data.

## Repository Structure

The repository includes the following components:

/
├── COMP0118_CW1_Report.pdf # Main report (max 3 pages)
├── COMP0118_Coursework1_Figures/ # Directory containing figures and tables
│ └── fig1.png, fig2.png, ... # Numbered and captioned figures referenced in the report
├── COMP0118_Coursework1_Code/ # MATLAB code directory
│ ├── Coursework1_Code.m # Main script covering all core tasks
│ ├── BallStickSSD_constraints.m # Ball-and-stick model with parameter constraints
│ ├── BallStickSSD_constraints_signal.m # Signal prediction function
│ ├── find_optimal_param_newdataset.m # Global minimum search with multiple initializations
│ ├── DT_starting_point.m # Linear tensor model used to generate starting points
│ └── ... # Additional helper scripts and function


## Requirements

- MATLAB (R2021a or later recommended)
- Optimization Toolbox (for `fminunc`, `fmincon`)
- Parallel Computing Toolbox (for `parfor`, if used)

## Running the Code

1. **Load the data**: Use the provided instructions to load `dwis`, `bvecs`, and generate `bvals`.
2. **Execute main script**: Run `Coursework1_Code.m` to reproduce the results discussed in the report.
3. **Check output**: Output parameter maps and plots will be displayed as figures.

## Submission Components

As required by the coursework specification, the submission consists of:

1. **Written report**: `COMP0118_CW1_Report.pdf`
2. **Figures document**: All figures referenced in the report are located in the `COMP0118_Coursework1_Figures/` directory.
3. **Code listing**: The full MATLAB code is contained within the `COMP0118_Coursework1_Code/` directory.

## Author

Francesco Seracini  
MSc Machine Learning  
University College London  

