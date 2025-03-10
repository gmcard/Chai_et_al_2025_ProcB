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

## Overview

These scripts automate the following workflow:  
1. **Real-Time Monitoring**: Tracks fish position using cameras and detectors.  
2. **Behavior-Triggered Rewards**: Dispenses rewards (eggs) based on fish behavior and randomized feeder selection.  
3. **Data Logging**: Records egg delivery timing, feeder selection, and video footage for analysis.  

## Scripts Overview

### 1. `Danionella_reward_experiment_2_feeders_feederROI.m` (previous version: `Danionella_reward_experiment_2_feeders.m`)
The main script for running reward-based experiments. Integrates hardware controls and manages experimental workflows.

- **Features**:
  - Detects fish position using pre-defined ROIs.
  - Controls solenoid valves for feeder operation.
  - Manages syringe pumps to dispense rewards.
  - Logs metadata and generates output files for analysis.
  - Saves video footage for each experimental session.

- **Inputs and Outputs**:
  - Inputs: Configuration parameters for feeders, detectors, and hardware.
  - Outputs: Experimental metadata and video recordings.

- **Dependencies**:
  - `Danionella_camera_initialization.m`: Initializes cameras.
  - `camera_control.m`: Manages BIAS cameras.
  - `view_detector_cameras.m`: Displays live feeds from egg detection cameras.

---

### 2. `syringe_pump_egg_prime.m`
Used to prime feeders and control syringe pumps prior to experiments.

- **Features**:
  - GUI for controlling syringe pump loading, solenoid valve opening and closing, egg reservoir pressure, and egg priming.  
  - Live feedback on system status.
  - Integrates with `view_detector_cameras.m` for live monitoring of egg detection cameras.

- **Usage**:
  - Load syringes with Hatch Brine eggs.
  - Prime feeders before experiments.

---

### 3. `view_detector_cameras.m`
Displays live feeds from the two egg detection cameras.  

- **Features**:
  - Real-time monitoring of the egg detection process.  
  - GUI displaying both detector feeds.  
  - Configurable serial ports for camera connections. 
  - Video display at 30 Hz.  

- **Usage**:
  - Called by `syringe_pump_egg_prime.m` and `Danionella_reward_experiment_2_feeders.m` for live camera visualization.  

---

### 4. `MovieConversion.m`
Processes and renames video files recorded during experiments. Converts videos from MJPG to MP4 format for analysis.

- **Key Functions**:
  - Renames files to match experiment metadata.
  - Validates folder naming conventions (`yyyy-mm-dd_runXXX_camX_movie`).
  - Uses `JAABA/filehandling/ReadIndexedMJPGHeader.m` for index checks.

---

### 5. `camera_control.m`
Manages BIAS camera operations, such as starting/stopping capture and saving configurations.

- **Key Commands**:
  - `startCapture`: Begins video recording.
  - `stopCapture`: Stops recording.
  - `saveConfiguration`: Saves camera settings to a file.

---

### 6. `Danionella_camera_initialization.m`
Sets up cameras for real-time fish tracking and video recording.

- **Features**:
  - Configures frame rate, ROI, and trigger mode for cameras.
  - Communicates with the BIAS system via `bias_gui.bat`.
  - Supports multiple cameras for simultaneous monitoring.

- **Referenced Files**:
  - **Configuration Files**:
    - `bias_config0.json`, `bias_config1.json`: Contain camera-specific settings.
  - **Background Image**:
    - `bg_image.png`: Used for background subtraction during fish detection.

- **Dependencies**:
  - **BiasControl Class**:
    - Facilitates communication with the BIAS system.
    - Provides functions like connecting/disconnecting cameras, adjusting settings, and managing ROIs.

---

### 7. `get_bias_params.m`
Provides pre-configured URLs and commands for interacting with the BIAS camera system.

- **Purpose**:
  - Interfaces with the BIAS software to manage fish detection and system state.

- **Key Commands**:
  - `statusurl`: Retrieves fish detection status.
  - `reseturl`: Resets the fish detection trigger.

---

### 8. `bias_gui.bat`
Launches BIAS camera GUI.

- **Referenced Files**:  
  - Contains the file path of the Bias .exe file.  
  
- **Usage**:
  - Called by `Danionella_camera_initialization.m`

---

### 9. `prfsfiletransfer.bat`
Automates file transfer from the local system to the PRFS server.

- **Usage**:
  - Schedule this script to run daily using Task Scheduler or a similar tool.
  - Ensures experiment data is backed up on the server.

---


### 10. `Aux Board Scripts`
Contains auxiliary scripts for hardware control:
- **`find_eggs_detector1_H7_R1.py`** and **`find_eggs_detector2_H7_R1.py`**:  
  - Control egg detectors and update ROIs.  
  - Used to adjust the detection threshold for experimental precision.  

- **`BubblePumpV2_CM.ino`**:
  - Arduino script for live reservoir control.
  - Adjusts pressure and flow rates during experiments.

- **`fishFeeding_miniBCS_2feeders.ino`**:
  - Arduino script for miniBCS control
  - Controls camera acquisition through TTL pulses sent to cameras.
  - Controls feeder valve closure after egg detection signal is received.
  - Sends signal to speaker.
---

## Libraries and Dependencies

### `encoder/`
Contains utilities for encoding and decoding data formats used in video processing.

### `jsonlab/`
Provides JSON file handling capabilities in MATLAB.

### `JAABA/filehandling/ReadIndexedMJPGHeader.m`
Used to validate and process MJPG file headers during video conversion.

---

## Usage Notes

1. **Preparing the System**:
   - Prime feeders with `syringe_pump_egg_prime.m`.

2. **Running the Experiment**:
   - Execute `reward_experiment.m` to begin the session.
   - Monitor fish behavior and ensure hardware is functioning properly.

3. **Post-Experiment Processing**:
   - Use `MovieConversion.m` to process video files.
   - Run `prfsfiletransfer.bat` to back up data to the server.

4. **Hardware Adjustments**:
   - Use scripts in `Aux Board Scripts` to fine-tune egg detectors and reservoir settings.
