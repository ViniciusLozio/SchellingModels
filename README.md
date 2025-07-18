# SchellingModels
Python project that implements a mean‑field approach to derive differential equations that simulate extended Schelling models and generate iterative time‑series data for the most important parameters.

This project was carried out over the course of my undergraduate studies in physics as a scientific initiation project under the supervision of Professor Dr. André de Pinho Vieira. It was completed in 2024 and provided several novel insights into the study of Schelling models, such as delineating the limits of a mean‑field approach and demonstrating the existence of phase transitions when specific parameters are varied.

## The Schelling Models




## Extendind Schelling Models



## The Mean Field Approach



## Generating Differential Equations




## Project Structure

.Directories:
  - Python: Python project
  - Misc: Reports and data analysis


.bi: Scripts with non‑vectorized functions for integration using the built‑in RK4 method. (Does not use NumPy or SciPy; compatible with PyPy)

 - MAIN_bi.py: Executes the dynamics using the built‑in RK4 integrator.

 - functions_bi.py: Utility functions used by the non‑vectorized scripts.

.Folder sp: Scripts with vectorized functions leveraging NumPy and SciPy for integration via SciPy’s solve_ivp.
