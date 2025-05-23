# Project Overview

This repository contains the code related to Chai et al., 2025, "Shorter-duration escapes driven by Drosophila giant interneurons promote survival during predation". The scripts are known to run with MATLAB R2022a.

# Project Abstract

Large axon-diameter descending neurons are metabolically costly but transmit information rapidly from sensory neurons in the brain to motor neurons in the nerve cord. They have thus endured as a common feature of escape circuits in many animal species where speed is paramount. Though often considered isolated command neurons triggering fast-reaction-time, all-or-none escape responses, giant neurons are one of multiple parallel pathways enabling selection between behavioral alternatives. Such degeneracy among escape circuits makes it unclear if and how giant neurons benefit prey fitness. Here we competed Drosophila melanogaster flies with genetically silenced Giant Fibers (GFs) against flies with functional GFs in an arena with wild-caught damselfly predators and find that GF silencing decreases prey survival. Kinematic analysis of damselfly attack trajectories shows that decreased prey survival results from predator capture of GF-silenced flies during some attack speeds and approach distances that would normally elicit successful escapes. In previous studies with a virtual looming stimulus, we proposed a model in which GFs enforce selection of a short-duration takeoff sequence as opposed to reducing reaction time. Our findings here demonstrate that, during real predation scenarios, the GFs indeed promote prey survival by influencing action selection as a means to increase escape probability. 


# Directory Structure

```
Chai_et_al_2025_ProcB/
│
├── Analysis/                       # Contains time and tracking data
│   ├── Time of Contact                       
│       └── Time 
│           └── Time_data.xlsx      # Time stamp for either contact or escape in each event 
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
├── Figure Data Files/              # Contains csv files with the data from each figure exported from Damselfly_analysis.m                                                     
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
├── padcat.m                        # Support function that concatenates vectors with different lengths by padding with NaN
│
└── shadedErrorBar.m                # Support function that generates error bars around a line plot

```
# Analysis

## Time of Contact/Time/Time_data.xlsx
Time from start of predation event until fly escape or capture by damselfly for all predation attempts analyzed (n=31). Column A is the time until fly escape or capture in seconds (ms), Column C is a constant (1000) used to convert between milliseconds and seconds, and Column E is the time until fly escape or capture in milliseconds (s).

## Time_fps/Time/Time_data.xlsx
Time matrix used for video analysis. Column A is time in seconds (s). Column E is the time interval between frames (0.001 s from sampling rate of 1000fps).

## *Predation_data_xyzpts.csv
Four subfolders in the Analysis folder contain .csv files of damselfly and fruit fly tracked points for each predation attempt generated using DLTdv software: The two DL(Wild) (escape or no escape) subfolders contain tracking results for GF1>+ GF-wildtype flies and the two kir2.1 (escape or no escape) subfolders contain tracking results for GF1>Kir2.1 GF-silenced flies. "*Predation_data_xyzpts.csv" below is a representative file description that also applies to all other files in the four subfolders.

Each column heading denotes the coordinates (X, Y or Z) of the three tracked points in the video ('pt1' = damselfly head centroid, 'Tail' = damselfly tail centroid, 'Fly' = fly centroid). Each row of tracked position values correspond to individual video frames throughout the recorded predation event.
*Different letters in the filename represent individual predation attempts.

---

# Box Analysis
This folder contains .csv files with initial fruit fly tracked head and tail points for each predation attempt generated using DLTdv software. The two DL(Wild) (escape or no escape) subfolders contain results for GF1>+ GF-wdildtype flies and the two kir2.1 (escape or no escape) subfolders contain results for GF1>Kir2.1 GF-silenced flies. The following is a representative file description that also applies to all other files in the "Box Analysis" folder.

## *FlyAxis_data_xyzpts.csv
Each column heading denotes the coordinates (X, Y or Z) of the two tracked points on the fly's body ('pt1' = fly head and 'pt2' = fly tail). Only the coordinates for the first video frame of the predation event was tracked. 
*Different letters in the filename represent individual predation attempts.

---

# Figure Data Files
This folder contains .csv files with the values of all data points plotted in the paper's figures.

## Damselfly_Fig1B.csv
Percentages (%) of flies assayed that took off upon stimulus presentation in the FlyPEZ for different looming stimulus speeds for each GF-specific split-GAL4 driver pairing. Column headings denote GF-specific split-GAL4 driver (GF1> or GF2>) and stimulus looming speed (r/v=40 or r/v=10) combination. Row headings denote UAS effector lines (DL wild-type +, Kir2.1;Gal80ts, or Kir2.1) that were crossed to split-GAL4 driver lines denoted in column headings. 

## Damselfly_Fig1C.csv
Escape sequence durations in milliseconds (ms) for flies with functional (GF>+) or silenced GFs (GF>Kir2.1). Column headings denote fly genotypes (GF>+ or GF>Kir2.1). Data were pooled from GF1>+ and GF2>+ fly genotypes presented with fast (r/v=10) and slow (r/v=40) looming speeds into GF>+ column. Data were pooled from GF1>Kir2.1 and GF2>Kir2.1 fly genotypes presented with fast (r/v=10) and slow (r/v=40) looming speeds into GF>Kir2.1 column.

## Damselfly_Fig1D.csv
Percentages (%) of escape takeoffs upon stimulus presentation in the FlyPEZ that were short mode takeoffs for different looming stimulus speeds for each GF-specific split-GAL4 driver pairing. Column headings denote GF-specific split-GAL4 driver (GF1> or GF2>) and stimulus looming speed (r/v=40 or r/v=10) combination. Row headings denote UAS effector lines (DL wild-type +, Kir2.1;Gal80ts, or Kir2.1) that were crossed to split-GAL4 driver lines denoted in column headings. 

## Damselfly_Fig1F.csv
Number of flies per genotype that were eaten for a subset of Wing-Clipping Bias (WCB) and Prey Consumption Index (PCI) survival competition assays. Column headings denote side of wing clipping (Right or Left) and fly genotype (GF1>+ or GF1>Kir2.1). Columns A to E are for Wing-Clipping Bias (WCB) plots and Columns G to K are for Prey Consumption Index (PCI) plots.

## Damselfly_Fig1G.csv
Columns A to B: Wing-Clipping Bias (WCB) index scores for different fly genotypes (row headings). Index values are unitless. 
Columns E to F: Prey Consumption Index (PCI) scores for different fly genotype matchups (row headings). Index values are unitless. 

## Damselfly_Fig2B.csv
Percentages (%) of damselfly predation events that resulted in fly capture or escape for flies with functional (GF1>+) or silenced GFs (GF1>Kir2.1). Column headings denote fly genotype (GF1>+ or GF1>Kir2.1) and fly escape or capture outcome combination.

## Damselfly_Fig2E.csv
Coordinates of damselflies in the X-Z plane for all predation events throughout attack time course. Columns A to J are for GF1>+ flies that were captured (n=5 flies, X (mm) and Z (mm) coordinates for each fly are reported in adjacent columns). Columns L to AK are for GF1>Kir2.1 flies that were captured (n=13 flies, X (mm) and Z (mm) coordinates for each fly are reported in adjacent columns). Columns AM to AZ are for GF1>+ flies that escaped (n=7 flies, X (mm) and Z (mm) coordinates for each fly are reported in adjacent columns). Columns BB to BM are for GF1>Kir2.1 flies that escaped (n=6 flies, X (mm) and Z (mm) coordinates for each fly are reported in adjacent columns). Time interval between rows is 1 millisecond (ms). 

## Damselfly_Fig2F.csv
Coordinates of damselflies in the X-Y plane for all predation events throughout attack time course. Columns A to J are for GF1>+ flies that were captured (n=5 flies, X (mm) and Y (mm) coordinates for each fly are reported in adjacent columns). Columns L to AK are for GF1>Kir2.1 flies that were captured (n=13 flies, X (mm) and Y (mm) coordinates for each fly are reported in adjacent columns). Columns AM to AZ are for GF1>+ flies that escaped (n=7 flies, X (mm) and Y (mm) coordinates for each fly are reported in adjacent columns). Columns BB to BM are for GF1>Kir2.1 flies that escaped (n=6 flies, X (mm) and Y (mm) coordinates for each fly are reported in adjacent columns). Time interval between rows is 1 millisecond (ms). 

## Damselfly_Fig2G.csv
Average damselfly elevation angles (°) relative to fly prey for predation events that resulted in fly capture (Column A to B) or escape (Column C to D) for flies with functional (GF1>+) or silenced GFs (GF1>Kir2.1). Column headings denote fly genotype (GF1>+ or GF1>Kir2.1) and predation event outcome (fly capture or fly escape) combination.

## Damselfly_Fig2H.csv
Average damselfly azimuth angles (°) relative to fly prey for predation events that resulted in fly capture (Column A to B) or escape (Column C to D) for flies with functional (GF1>+) or silenced GFs (GF1>Kir2.1). Column headings denote fly genotype (GF1>+ or GF1>Kir2.1) and predation event outcome (fly capture or fly escape) combination.

## Damselfly_Fig3A-Dist.csv
Damselfly distance from fly prey (mm) throughout attack time course for predation events that resulted in fly capture or escape for flies with functional (GF1>+) or silenced GFs (GF1>Kir2.1). Columns A to J are for predation events resulting in GF1>+ flies being captured (n=5 damselflies, Time (ms) and Distance from fly (mm) values for each damselfly are reported in adjacent columns). Columns L to AK are for predation events resulting in GF1>Kir2.1 flies being captured (n=13 damselflies, Time (ms) and Distance from fly (mm) values for each damselfly are reported in adjacent columns). Columns AM to AZ are for predation events resulting in GF1>+ flies escaping (n=7 damselflies, Time (ms) and Distance from fly (mm) values for each damselfly are reported in adjacent columns). Columns BB to BM are for predation events resulting in GF1>Kir2.1 flies escaping (n=6 damselflies, Time (ms) and Distance from fly (mm) values for each damselfly are reported in adjacent columns). Time= 0 ms is the time of either fly prey capture or escape.

## Damselfly_Fig3A-Speed.csv
Damselfly speed (m/s) throughout attack time course for predation events that resulted in fly capture or escape for flies with functional (GF1>+) or silenced GFs (GF1>Kir2.1). Columns A to J are for predation events resulting in GF1>+ flies being captured (n=5 damselflies, Time (ms) and Speed (m/s) values for each damselfly are reported in adjacent columns). Columns L to AK are for predation events resulting in GF1>Kir2.1 flies being captured (n=13 damselflies, Time (ms) and Speed (m/s) values for each damselfly are reported in adjacent columns). Columns AM to AZ are for predation events resulting in GF1>+ flies escaping (n=7 damselflies, Time (ms) and Speed (m/s) values for each damselfly are reported in adjacent columns). Columns BB to BM are for predation events resulting in GF1>Kir2.1 flies escaping (n=6 damselflies, Time (ms) and Speed (m/s) values for each damselfly are reported in adjacent columns). Time= 0 ms is the time of either fly prey capture or escape.

## Damselfly_Fig3A-Acceleration.csv 
Damselfly acceleration (m/s^2) throughout attack time course for predation events that resulted in fly capture or escape for flies with functional (GF1>+) or silenced GFs (GF1>Kir2.1). Columns A to J are for predation events resulting in GF1>+ flies being captured (n=5 damselflies, Time (ms) and Acceleration (m/s^2) values for each damselfly are reported in adjacent columns). Columns L to AK are for predation events resulting in GF1>Kir2.1 flies being captured (n=13 damselflies, Time (ms) and Acceleration (m/s^2) values for each damselfly are reported in adjacent columns). Columns AM to AZ are for predation events resulting in GF1>+ flies escaping (n=7 damselflies, Time (ms) and Acceleration (m/s^2) values for each damselfly are reported in adjacent columns). Columns BB to BM are for predation events resulting in GF1>Kir2.1 flies escaping (n=6 damselflies, Time (ms) and Acceleration (m/s^2) values for each damselfly are reported in adjacent columns). Time= 0 ms is the time of either fly prey capture or escape.

## Damselfly_Fig3A-n.csv
Number of predation events contributing to each time point for predation events that resulted in fly capture or escape for flies with functional (GF1>+) or silenced GFs (GF1>Kir2.1). Columns A to B are for predation events resulting in GF1>+ flies being captured (Time (ms) and Number of predation events at each timepoint n are reported in adjacent columns, total number of damselflies=5). Columns D to E are for predation events resulting in GF1>Kir2.1 flies being captured (Time (ms) and Number of predation events at each timepoint n are reported in adjacent columns, total number of damselflies=13). Columns G to H are for predation events resulting in GF1>+ flies escaping (Time (ms) and Number of predation events at each timepoint n are reported in adjacent columns, total number of damselflies=7). Columns J to K are for predation events resulting in GF1>Kir2.1 flies escaping (Time (ms) and Number of predation events at each timepoint n are reported in adjacent columns, total number of damselflies=6). Time= 0 ms is the time of either fly prey capture or escape.

## Damselfly_Fig3B.csv
Damselfly peak speeds (m/s) during predation events that resulted in fly capture (Column A to B) or escape (Column C to D) for flies with functional (GF1>+) or silenced GFs (GF1>Kir2.1).

## Damselfly_Fig3C.csv
Damselfly peak accelerations (m/s^2) during predation events that resulted in fly capture (Column A to B) or escape (Column C to D) for flies with functional (GF1>+) or silenced GFs (GF1>Kir2.1). 

## Damselfly_Fig3D.csv
Percentages of flies (%) with functional (GF1>+) or silenced GFs (GF1>Kir2.1) that escaped damselfly attacks with peak speeds slower or faster than 0.167 m/s.

## Damselfly_Fig3E.csv
Peak damselfly speed (m/s) was binned with the center bin value reported in Column A. Column B and C are the number of damselflies in each peak speed bin for the GF>+ and GF1>Kir2.1 cases respectively. 

## Damselfly_Fig3F.csv
Angular sizes of damselflies on fly’s retina (°) at time of capture or escape for predation events that resulted in fly capture (Column A to B) or escape (Column C to D) for flies with functional (GF1>+) or silenced GFs (GF1>Kir2.1). Column headings denote fly genotype (GF1>+ or GF1>Kir2.1) and predation event outcome (fly capture or fly escape) combination.

## Damselfly_Fig4A.csv
Loom angular sizes (°) throughout stimulus expansion time course in milliseconds (ms) in the FlyPEZ for different stimulus looming speeds (r/v=40: Columns A to B and r/v=10: Columns D to E).

## Damselfly_Fig4B.csv
Time (ms) after start of looming stimulus expansion at which fly executed wing lift and takeoff behaviors for different looming stimulus speeds in the FlyPEZ for flies with functional (GF1>+) or silenced GFs (GF1>Kir2.1). Columns A and B are time of wing lift (ms) and takeoff (ms) for GF1>+ flies challenged with looming stimulus speed r/v=40. Columns D and E are time of wing lift (ms) and takeoff (ms) for GF1>+ flies challenged with looming stimulus speed r/v=10. Columns G and H are time of wing lift (ms) and takeoff (ms) for GF1>Kir2.1 flies challenged with looming stimulus speed r/v=40. Columns J and K are time of wing lift (ms) and takeoff (ms) for GF1>Kir2.1 flies challenged with looming stimulus speed r/v=10. Each row represents behavior time points for an individual fly.

## Damselfly_Fig4C.csv
Fly reaction times (ms) for different looming stimulus speeds (r/v=40: Columns A to B and r/v=10: Columns C to D) in the FlyPEZ for flies with functional (GF1>+) or silenced GFs (GF1>Kir2.1). Column headings denote looming stimulus speed (r/v=40 or r/v=10) and fly genotype (GF1>+ or GF1>Kir2.1) combination. Each row represents the reaction time for an individual fly.

## Damselfly_Fig4D.csv
Fly latencies (ms) to the end of the takeoff sequence for different looming stimulus speeds (r/v=40: Columns A to B and r/v=10: Columns C to D) in the FlyPEZ for flies with functional (GF1>+) or silenced GFs (GF1>Kir2.1). Column headings denote looming stimulus speed (r/v=40 or r/v=10) and fly genotype (GF1>+ or GF1>Kir2.1) combination. Each row represents the latency to the end of the takeoff sequence for an individual fly.

## Damselfly_FigS1C.csv
Wing-Clipping Bias (WCB) index scores for different fly genotypes (GF1>+, GF1>Kir2.1, GF1>Kir2.1;Gal80ts, GF2>+, GF2>Kir2.1). Index values are unitless. 

## Damselfly_FigS2A.csv
Elevation angles (°) of damselflies relative to fly prey throughout attack time course for predation events that resulted in fly capture or escape for flies with functional (GF1>+) or silenced GFs (GF1>Kir2.1). Columns A to J are for predation events resulting in GF1>+ flies being captured (n=5 damselflies, Time (ms) and Elevation angle (°) values for each damselfly are reported in adjacent columns). Columns L to AK are for predation events resulting in GF1>Kir2.1 flies being captured (n=13 damselflies, Time (ms) and Elevation angle (°) values for each damselfly are reported in adjacent columns). Columns AM to AZ are for predation events resulting in GF1>+ flies escaping (n=7 damselflies, Time (ms) and Elevation angle (°) values for each damselfly are reported in adjacent columns). Columns BB to BM are for predation events resulting in GF1>Kir2.1 flies escaping (n=6 damselflies, Time (ms) and Elevation angle (°) values for each damselfly are reported in adjacent columns). Time= 0 ms is the time of either fly prey capture or escape.

## Damselfly_FigS2B.csv
Azimuth angles (°) of damselflies relative to fly prey throughout attack time course for predation events that resulted in fly capture or escape for flies with functional (GF1>+) or silenced GFs (GF1>Kir2.1). Columns A to J are for predation events resulting in GF1>+ flies being captured (n=5 damselflies, Time (ms) and Azimuth angle (°) values for each damselfly are reported in adjacent columns). Columns L to AK are for predation events resulting in GF1>Kir2.1 flies being captured (n=13 damselflies, Time (ms) and Azimuth angle (°) values for each damselfly are reported in adjacent columns). Columns AM to AZ are for predation events resulting in GF1>+ flies escaping (n=7 damselflies, Time (ms) and Azimuth angle (°) values for each damselfly are reported in adjacent columns). Columns BB to BM are for predation events resulting in GF1>Kir2.1 flies escaping (n=6 damselflies, Time (ms) and Azimuth angle (°) values for each damselfly are reported in adjacent columns). Time= 0 ms is the time of either fly prey capture or escape.

## Damselfly_FigS2C.csv
Number of predation events contributing to each time point throughout attack time course for predation events that resulted in fly capture or escape for flies with functional (GF1>+) or silenced GFs (GF1>Kir2.1). Columns A to B are for predation events resulting in GF1>+ flies being captured (Time (ms) and Number of predation events at each timepoint n are reported in adjacent columns, total number of damselflies=5). Columns D to E are for predation events resulting in GF1>Kir2.1 flies being captured (Time (ms) and Number of predation events at each timepoint n are reported in adjacent columns, total number of damselflies=13). Columns G to H are for predation events resulting in GF1>+ flies escaping (Time (ms) and Number of predation events at each timepoint n are reported in adjacent columns, total number of damselflies=7). Columns J to K are for predation events resulting in GF1>Kir2.1 flies escaping (Time (ms) and Number of predation events at each timepoint n are reported in adjacent columns, total number of damselflies=6). Time= 0 ms is the time of either fly prey capture or escape.

## Damselfly_FigS3A.csv
Damselfly approach speeds (m/s) throughout attack time course for predation events ending in GF1>Kir2.1 fly capture. Columns A to N are for predation events with peak speeds that are <0.167 m/s (n=7 damselflies, Time (ms) and Speed (m/s) values for each damselfly are reported in adjacent columns). Columns P to AA are for predation events with peak speeds that are >0.167 m/s (n=6 damselflies, Time (ms) and Speed (m/s) values for each damselfly are reported in adjacent columns). Time= 0 ms is the time of fly prey capture.

## Damselfly_FigS3B.csv
Distances between damselflies and flies (mm) at time of capture or escape for predation events that resulted in fly capture (Column A to B) or escape (Column C to D) for flies with functional (GF1>+) or silenced GFs (GF1>Kir2.1). Column headings denote fly genotype (GF1>+ or GF1>Kir2.1) and predation event outcome (fly capture or fly escape) combination.

## Damselfly_FigS3C.csv
Damselfly head widths in millimeters (mm).

---

# FlyPEZ Data
This folder contains files from many flyPEZ experiments, but there are two general file types, `*_dataForVisualization.mat` and `*_manualAnnotations.mat` which are prepended with experiment ID codes specific to the flyPEZ data organization structure. 

## *_dataForVisualization.mat
MATLAB structure that results from a flyPEZ data processing pipeline. It contains variables used in common flyPEZ analyses. Relevant variables for this project are: manualJumpTest, a manual annotation of whether a fly took off (1), or did not take off (0), and stimStart, the frame number the looming stimulus began. 

## *_manualAnnotations.mat
MATLAB structure that contains manual annotation data from flyPEZ experiments related to takeoff timing. Relevant variables for this project are: frame_of_wing_movement, and frame_of_leg_push, frame_of_take_off. These capture the frame number at which the fly first moved its wings in preparation of takeoff, first started pushing off with its legs, and lifted its legs from the platform, respectively. These variables were used in this project in combination with stimStart from `*_dataForVisualization.mat` to calculate takeoff timing parameters.

---

# Scripts Overview

## 1. `Damselfly_analysis`
The main analysis script that generates all plots in the paper figures 

- **Inputs and Outputs**:
    - Inputs: Fruit fly and damselfly tracked positions, time series information, flyPEZ data
    - Outputs: Paper figure plots

- **Dependencies**:
    - `align_axislabel.m`: Aligns axes labels to axes in 3D plots              
    - `axislabel_translation.m`: Controls distance between axes labels and axes               
    - `axislabel_translation_slider.m`: Enables dynamic manipulation of `axislabel_translation.m`              
    - `circ_mean.m`: Computes the mean direction for circular data    
    - `get_error_bars.m`: Computes the Wilson Score Interval for binary data
    - `golayDifferentiate.m`: Filters data and calculates differentials         
    - `padcat.m`: Concatenates vectors with different lengths by padding with NaN
    - `shadedErrorBar.m`: Generates error bars around a line plot

---


# Author and Contact Information
AUTHORS: Cynthia M. Chai [1,2], Carmen M. Morrow [1], Dhyey D. Parikh [1], Catherine R. von Reyn [1,3], Anthony Leonardo [1,4], Gwyneth M. Card [1,5,6,*]

AFFILIATIONS:
1 Janelia Research Campus, Howard Hughes Medical Institute, Ashburn, VA, USA
2 Current address: Department of Biological Sciences, Columbia University, New York, NY, USA
3 Current address: School of Biomedical Engineering, Science & Health Systems, Drexel University, Philadelphia, PA, USA 
4 Current address: Wispr, 444 Townsend St, San Francisco, CA
5 Current address: Mortimer B. Zuckerman Mind Brain Behavior Institute, Columbia University, New York, NY, USA
6 Current address: Howard Hughes Medical Institute, Chevy Chase, MD, USA 
* Corresponding author: gc3017@columbia.edu (G.M.C.)

