# Cervical tSCS EMG Artifact Rejection Pipeline

## Overview
This repository provides a MATLAB-based signal processing pipeline designed to extract physiological electromyography (EMG) signals during spinal cord stimulation (SCS). 

Because SCS artifacts can mask the underlying EMG signals. This toolkit can enable manual or dynamical detection of these stimulation artifacts and utilizes targeted notch filtering suppress them, preserving the integrity and time-domain morphology of the EMG data.

## Citation and Data Availability
The sample data (`SingleCh_StimOFF.mat` and `SingleCh_StimON.mat`) provided in the `Data` directory is derived from published clinical trials and the pipeline have been explained in the extended figure of the following publication:

* **Authors**: Nikhil Verma, Jeonghoon Oh, Ernesto Bedoy, Nikole Chetty, Alexander G. Steele, Seo Jeong Park, Jaime R. Guerrero, Amir H. Faraji, Douglas Weber, and Dimitry G. Sayenko
* **Title**: Transcutaneous stimulation of the cervical spinal cord facilitates motoneuron firing and improves hand-motor function after spinal cord injury
* **Journal**: Journal of Neurophysiology, Volume 134, Issue 1
* **DOI**: [https://doi.org/10.1152/jn.00422.2024](https://doi.org/10.1152/jn.00422.2024)

## Repository Structure & Core Methodologies
This repository contains two distinct pipeline architectures for artifact rejection, depending on the complexity of your noise floor.

### 1. Manual Fixed-Threshold Pipeline
**Script**: `Stim_artifact_removal_pipeline_Manual_fixed_threshold.m`
* **How it works**: A semi-automated approach. It uses static notch filters to remove known ambient line noise (e.g., 60 Hz and harmonics) from the Stim-OFF baseline. For the Stim-ON data, it prompts the user with an interactive plot to manually select a peak height threshold and minimum peak distance. **Crucially, if the stimulation artifact peaks are smaller than the high-power, low-frequency physiological EMG band, the user can set a lower frequency cutoff (`lower_freq_cutoff`). This restricts the detection and notching algorithm exclusively to frequencies above this cutoff, perfectly preserving the broadband EMG while still eliminating higher-frequency stimulation harmonics.**

### 2. Adaptive Savitzky-Golay Pipeline
**Script**: `Artifact_Rejection_Pipeline_Adaptive_threshold_SavGol_filter_on_baseline_.m`
* **How it works**: This is an advanced, automated pipeline. It first computes the frequency spectrum of a baseline, Stim-OFF recording. It applies a Savitzky-Golay smoothing filter to model the $1/f$ physiological noise floor. It then calculates the local Median Absolute Deviation (MAD) to create a dynamic, frequency-dependent threshold. When applied to the Stim-ON data, it isolates and notches only the non-physiological peaks that pierce this dynamic threshold.


## Repository Structure 
* `/Data/`: Contains sample `.mat` files for Stim-ON and Stim-OFF conditions.
* `Artifact_Rejection_Pipeline_Adaptive_threshold_SavGol_filter_on_baseline_.m`: The primary, fully automated pipeline using Savitzky-Golay baseline modeling.
* `Stim_artifact_removal_pipeline_Manual_fixed_threshold.m`: A semi-automated pipeline for manual peak detection.
* `find_notch_peaks_for_stim_artifact.m`: Helper function for the manual pipeline to interactively select notch targets.
* `make_pretty.m`: Helper function for publication-ready plot formatting.

## Installation & Requirements
* **MATLAB**: R2021a or newer is recommended.
* **Toolboxes**: The **Signal Processing Toolbox** is required (relies on `filtfilt`, `butter`, `findpeaks`).
