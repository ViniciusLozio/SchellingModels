# SchellingModels
Python project that implements a mean‑field approach to derive differential equations that simulate extended Schelling models and generate iterative time‑series data for the most important parameters.

This project was carried out over the course of my undergraduate studies in physics as a scientific initiation project under the supervision of Professor Dr. André de Pinho Vieira. It was completed in 2024 and provided several novel insights into the study of Schelling models, such as delineating the limits of a mean‑field approach and demonstrating the existence of phase transitions when specific parameters are varied.

## The Schelling Models

In 1969, the American economist and 2005 Nobel Prize laureate in Economics, Thomas Schelling (1921–2016), published the article “Models of Segregation” [1], in which—aiming to study racial segregation in American cities—he introduced one of the first mathematical models of social agents. The model consisted of two types of cells occupying predefined sites on a lattice; these cells could swap positions to increase their “satisfaction”, quantified by the number of like‑type cells in their eight‑cell Moore neighborhood.

<img width="414" height="195" alt="schelling" src="https://github.com/user-attachments/assets/d559e7cc-6529-4eef-a231-e091b9a848f2" />


## Extendind Schelling Models

Since Schelling's article, studies have been conducted on Schelling models, exploring different dynamics, neighbourhoods and cells fractions. 

The model in this project works as follows:

- There are two types of agents (type 0 and type 1), whose fractions of the total are denoted by ρ₀ and ρ₁ (with ρ₀ + ρ₁ = 1).
- Each agent has a neighborhood of 4 agents (a von Neumann neighborhood).
- Satisfaction rules are encoded as a binary string s₀ s₁ s₂ s₃ s₄, where each sₖ is a Boolean indicating satisfaction (s = 1) or dissatisfaction (s = 0) when there are k same-type neighbors.
  For example, the rule 10001 means agents are satisfied with 0 or 4 same-type neighbors (s₀ = s₄ = 1) and dissatisfied with 1, 2, or 3 same-type neighbors (s₁ = s₂ = s₃ = 0). In this study the same satisfaction rule was applied to both agent types.
- Dynamics proceed in discrete steps: at each step, a randomly chosen dissatisfied type-0 agent swaps positions with a randomly chosen dissatisfied type-1 agent, even if both remain dissatisfied after the swap.

 <img width="222" height="402" alt="model" src="https://github.com/user-attachments/assets/df4cbf7c-e325-401f-bb6a-6a9bb0781dfa" />


## The Mean Field Approach

One way to simplify the calculations of a Schelling model is by using the mean‑field approach, widely employed in Statistical Mechanics. This allows one to describe the system’s overall average behavior by treating each agent’s interactions as a single average interaction. In a Schelling model, agents seek satisfaction based on the number of same‑type neighbors they have. Under the mean‑field approximation, an agent’s satisfaction depends only on the average fraction of similar neighbors, thereby greatly simplifying the computations.

## Generating Differential Equations

The core element of the code is the computation of five 4×4 matrices. The calculation of these matrices is directly tied to the very definition of the system’s dynamics and the satisfaction rule being used. 

From these matrices, a set of differential equations is formulated, allowing for the numerical computation of the rate of change over time of ρₖ, the fraction of agents of a given type with k neighbors of type 1 in their neighborhood. Through numerical integration, like Runge-Kutta method, one can then analyze the evolution of these fractions as well as the total number of dissatisfied agents of each type.


## Project Structure

### Directories:
  - Python: Python project
  - Misc: Reports and data analysis


.bi*: Scripts with non‑vectorized functions for integration using the built‑in RK4 method. (Does not use NumPy or SciPy; compatible with PyPy)

 - MAIN_bi.py: Executes the dynamics using the built‑in RK4 integrator.

 - functions_bi.py: Utility functions used by the non‑vectorized scripts.

.sp*: Scripts with vectorized functions leveraging NumPy and SciPy for integration via SciPy’s solve_ivp.
