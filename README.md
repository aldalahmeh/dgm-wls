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

### Core Classes
* `PowerGrid.m`: Class handling the IEEE test case instantiation, power flow ground truth generation, and sensor strategy mapping.
* `MonteCarloPerformance.m`: Evaluation class for computing Mean Absolute Error (MAE), convergence statistics, F1 scores, and generating publication-ready figures.

### Top-Level Scripts
* `genSimPowerDataIEEEtestCases.m`: Generates the baseline power flow data and ground truth states for the test grid.
* `simDbscanGmmTuning.m`: Executes hyperparameter tuning for the DBSCAN-GMM filtering stage to optimize the F1 detection score.
* `simGridWSNSystemAlg.m`: The primary simulation script that runs the WLS, DGM-WLS, and IRWLS state estimation algorithms across the grid network.
* `computeMonteCarloPerformance.m`: Evaluates the results of a single Monte Carlo simulation run.
* `computeMonteCarloPerformanceBatch.m`: Processes and evaluates multiple simulation results in batches across varying noise probability scenarios.
* `computeMonteCarloPerformanceAve.m`: Computes the grid-wide average performance metrics and generates the final trend figures.

## Usage Workflow
1. Ensure MATPOWER is installed and added to your MATLAB path.
2. Clone or download this repository.
3. **Generate Data:** Run `genSimPowerDataIEEEtestCases.m` to establish the ground truth.
4. **Tune Filters (Optional):** Run `simDbscanGmmTuning.m` if you wish to replicate the hyperparameter selection process.
5. **Run Estimation:** Execute `simGridWSNSystemAlg.m` to perform the state estimation under impulsive noise.
6. **Evaluate:** Use the `computeMonteCarloPerformance` scripts to process the saved data and generate the `.eps` performance figures.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation
If you use this code in your research, please cite the corresponding paper:
> [Author Names], "[Your Paper Title]," in *IET Smart Cities*, 2026. (DOI pending)
