# Zymomonas_Modeling
This repository contains code and documentation for reconstructing and refining a genome-scale metabolic model of *Zymomonas mobilis* using MATLAB and the COBRA Toolbox.

## Project Overview

The goal of this project is to simulate the metabolism of *Z. mobilis* under anaerobic conditions using flux balance analysis (FBA). The model is based on the iZM478 template and includes:

- Manual curation and reaction edits
- Gap-filling for biomass production
- Oxygen closure for anaerobic simulation
- FBA to predict growth and metabolite yields

## Folder Structure

- `/sdsu_model/`: MAT files of draft and refined models  
- `/scripts/`: MATLAB scripts for editing and running the model  
- `/results/`: FBA outputs, figures, and logs  
- `/docs/`: Background information and methods

## How to Run

1. Install MATLAB with COBRA Toolbox  
2. Load your model file (e.g., `SDSU.mat`)  
3. Run simulation using `optimizeCbModel`  
4. Adjust exchange reactions (e.g., block O₂) as needed  

## Author
Nhi Nguyen  
CSU ARI NEXTGen Fellow, San Diego State University
