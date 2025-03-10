# Project Overview

This repository contains the code related to Chai et al., 2025, "A survival-critical role for Drosophila giant interneurons during predation" The scripts are known to run with MATLAB R2022a.

## Directory Structure

```
Chai_et_al_2025_ProcB/
│
├── Analysis/                       # Contains time and tracking data
│   ├── Time of Contact                       
│       └── Time 
│           └── Time_data.xlsx      # Time stamp for either contract or escape in each event 
│   ├── Time_fps              
│       └── Time 
│           └── Time_data.xlsx      # Time matrix corresponding to tracked points
│   ├── dl(wild) escape             # Damselfly and fruit fly tracked points for predation attempts with dl(wild) flies that escaped
│   ├── dl(wild) no escape          # Damselfly and fruit fly tracked points for predation attempts with dl(wild) flies that did not escape
│   ├── kir2.1 escape               # Damselfly and fruit fly tracked points for predation attempts with GF silenced flies that escaped
│   ├── kir2.1 no escape            # Damselfly and fruit fly tracked points for predation attempts with GF silenced flies that did not escape
│
├── Box Analysis/                   # Initial fruit fly tracked head and tail points for each predation attempt
│   ├── DL(Wild) Escape
│   ├── DL(Wild) No Escape
│   ├── kir2.1 escape
│   ├── kir2.1 no escape                   
│
├── FlyPEZ Data/                    # Contains data files from FlyPEZ escape timing experiments                                                      
│
├── Damselfly_analysis.m            # Main analysis script that generates all plots in the paper figures                     
│
├── align_axislabel.m               # Support function that aligns axes labels to axes in 3D plots              
│
├── axislabel_translation.m         # Support function that controls distance between axes labels and axes               
│
├── axislabel_translation_slider.m  # Support function that enables dynamic manipulation of axislabel_translation.m              
│
├── circ_mean.m                     # Support function that computes the mean direction for circular data    
│
├── get_error_bars.m                # Support function that computes the Wilson Score Interval for binary data
│
├── golayDifferentiate.m            # Support function that filters data and calculates differentials         
│
├── padcat.m                        # Support function that concatenate vectors with different lengths by padding with NaN
│
└── shadedErrorBar.m                # Support function that generates error bars around a line plot

```

## Data Files

Describe types of data files, variables in each column, and units

## Scripts Overview

### 1. `Damselfly_analysis`
The main analysis script that generates all plots in the paper figures 

- **Features**:


- **Inputs and Outputs**:

- **Dependencies**:

---


## Contact Information

