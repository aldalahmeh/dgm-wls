# DGM-WLS: Distributed Gaussian Mixture - Weighted Least Squares

This repository contains the MATLAB source code for the proposed Distributed Gaussian Mixture Weighted Least Squares (DGM-WLS) algorithm. This localized, edge-computing architecture is designed for resilient power system state estimation under high-density impulsive noise conditions.

## Overview
Standard global estimators (like Ordinary WLS) fail systematically under severe impulsive noise, while computationally heavy robust estimators (like IRWLS) suffer from the "swamping effect" under extreme noise densities. 
The DGM-WLS framework leverages a localized DBSCAN-GMM filtering stage at the sensor level to isolate and remove bad data *prior* to the central estimation step. This ensures the central WLS matrix is populated exclusively with pre-filtered data, yielding high scalability, lower computational iterations, and superior error bounding.

## Prerequisites
To run this code, you will need:
* **MATLAB** (Tested on version 2025a)
* **MATPOWER** (Open-source power system simulation package)

## Repository Structure
* `PowerGrid.m`: Class handling the IEEE 30-bus grid instantiation, power flow ground truth generation, and sensor strategy mapping.
* `MonteCarloPerformance.m`: Evaluation class for computing Mean Absolute Error (MAE), convergence statistics, F1 scores, and generating publication-ready IEEE-compliant figures.
* `[Insert_Your_Main_Script_Name].m`: The top-level batch script to run the varying noise probability scenarios ($p=0.001$ to $p=0.1$).

## Usage
1. Ensure MATPOWER is installed and added to your MATLAB path.
2. Clone or download this repository.
3. Open `[Insert_Your_Main_Script_Name].m` in MATLAB and run the script. 
4. The script will automatically generate the results and save the corresponding `.eps` figures to the `\Figures` directory.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation
If you use this code in your research, please cite the corresponding paper:
> [Author Names], "[Your Paper Title]," in *IET Smart Cities*, 2026. (DOI pending)
