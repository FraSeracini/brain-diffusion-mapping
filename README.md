# COMP0118 – Coursework 1: Parametric Models

## Overview

This project was completed as part of the COMP0118 module *Computational Modelling for Biomedical Imaging* at UCL. The coursework focuses on implementing and analyzing parametric models such as the linear diffusion tensor and the ball-and-stick model for processing and interpreting diffusion MRI data.

## Repository Structure

The repository contains the following main components:

- **Coursework1.pdf**  
  This is the coursework specification document provided by the instructor. It outlines the core and advanced tasks, required deliverables, and submission instructions.

- **COMP0118_Coursework1.pdf**  
  This is the final written report (maximum 3 pages) that describes the methodology, analysis, and results for all the required tasks. It references both the code and figures.

- **COMP0118_Coursework1_Figures**  
  This document contains all figures and tables referenced in the report. Each figure is numbered and includes a short caption, as required by the coursework instructions.

- **COMP0118_Coursework1_Code/**  
  This folder includes all MATLAB scripts and functions used to perform the parameter estimation, model fitting, and mapping tasks described in the report. It contains:

## Requirements

- MATLAB (R2021a or later recommended)
- Optimization Toolbox (for `fminunc`, `fmincon`)
- Parallel Computing Toolbox (for `parfor`, if used)

## Running the Code

1. **Load the data**: Use the provided instructions to load `dwis`, `bvecs`, and generate `bvals`.
2. **Execute scripts for each question**: Run `Q.1.1`, `Q.1.2`, ...,  to reproduce the results discussed in the report.
3. **Check output**: Output parameter maps and plots will be displayed as figures.

## Submission Components

As required by the coursework specification, the submission consists of:

1. **Written report**: `COMP0118_Coursework1.pdf`
2. **Figures document**: All figures referenced in the report are inside the `COMP0118_Coursework1_Figures` file.
3. **Code listing**: The full MATLAB code is contained within the `COMP0118_Coursework1_Code/` directory.

## Author

Francesco Seracini  
MSc Computer Science - Artificial Intelligence  
Politecnico di Milano - University College London  

