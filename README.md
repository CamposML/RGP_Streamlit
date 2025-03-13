# README for Monte Carlo Simulation for PTA and CFR Analysis

## Overview
This repository contains Python code to perform Monte Carlo simulations for analyzing pharmacokinetics/pharmacodynamics (PK/PD) data. Specifically, it calculates:

- **Probability of Target Attainment (PTA):** The likelihood that a given dosing regimen achieves a specific pharmacodynamic target.
- **Cumulative Fraction Response (CFR):** A measure of overall effectiveness against a distribution of Minimum Inhibitory Concentrations (MICs).

The analysis supports drug development and dosing optimization by simulating diverse patient populations and drug regimens.

### Key Features
- Lognormal distributions for pharmacokinetic parameters.
- Configurable dosing regimens and MIC values.
- PTA results saved as CSV matrix.
- CFR calculations for MIC distributions.
- Visualization of PTA vs. MIC and CFR.

## Getting Started

### Prerequisites

- Python 3.8+
- Required Python libraries:
  - `math`
  - `numpy`
  - `matplotlib`
  - `csv`

Install required libraries using:
```bash
pip install numpy matplotlib
```

### Clone Repository
Clone the repository from GitHub:
```bash
git clone https://github.com/CamposML/RGP.git
cd RGP
```

### Usage

The script is designed to be executed as a standalone file. To run the simulation:

```bash
python <filename>.py
```

#### Input Parameters
The following parameters can be modified in the script under the `if __name__ == "__main__":` block:

- **Dosing Regimens**: `doses_intervals` specifies dose (mg) and interval (h):
  ```python
  doses_intervals = [(2000, 24), (2000, 12), (4000, 24), (6000, 12)]
  ```

- **PK Parameters**:
  - Fraction unbound: `fu_mean`, `fu_std`
  - Volume of distribution (L): `vd_mean`, `vd_std`
  - Clearance (L/h): `clt_mean`, `clt_std`

- **MIC Values**: `mics` contains a list of MIC values (mg/L):
  ```python
  mics = [0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
  ```

- **MIC Distribution**: Distribution fractions for MIC values:
  ```python
  mic_distribution = {0.125: 0.0, 0.25: 0.0, 0.5: 0.0, 1.0: 0.011, ...}
  ```

- **Number of Patients**: `num_patients` defines the number of simulated patients:
  ```python
  num_patients = 5000
  ```

- **PK/PD Target**: `target` specifies the %fT > MIC required to achieve the target:
  ```python
  target = 55
  ```

#### Outputs
1. **PTA Results (CSV):**
   PTA results are saved as `pta_results_matrix.csv`, presenting MIC values vs. dosing regimens.

2. **CFR Results (CSV):**
   CFR results are saved as `cfr_results.csv` with corresponding CFR percentages for each dosing regimen.

3. **Visualizations:**
   - **PTA vs. MIC**: Plotted for all dosing regimens.
   - **CFR Bar Chart**: Summarizing CFR percentages.

#### Example Run
Example configurations and results are already set in the script for ease of use. Simply execute the file to generate CSV outputs and visualizations.

## Code Structure

1. **Functions:**
   - `calculate_ft_mic`: Calculates %fT > MIC for a given regimen.
   - `calculate_lognormal_params`: Converts mean and std to lognormal parameters.
   - `monte_carlo_simulation`: Performs Monte Carlo simulations for PTA calculations.
   - `calculate_cfr`: Computes CFR based on MIC distributions.
   - `save_pta_to_csv`: Exports PTA results to CSV.
   - `save_cfr_to_csv`: Exports CFR results to CSV.
   - `plot_pta_vs_mic`: Visualizes PTA vs. MIC results.
   - `plot_cfr_bar_chart`: Creates a CFR bar chart.

2. **Main Script:**
   Configurable input parameters for direct execution.

## Licensing
This project is currently under a **Proprietary License**. Redistribution and use are prohibited without explicit permission.

## Contributing
Contributions are currently not open for this repository. For inquiries, please contact the repository owner through GitHub.

## Support
For issues or questions, submit an inquiry on the GitHub repository: [RGP Issues](https://github.com/CamposML/RGP/issues).

---

Happy simulating!

