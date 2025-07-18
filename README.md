# A Mean-Field Approach for Extended Schelling Models
Python project that implements a mean‑field approach to derive differential equations that simulate extended Schelling models and generate iterative time‑series data for the most important parameters.

This project was carried out over the course of my undergraduate studies in physics as a scientific initiation project under the supervision of Professor Dr. André de Pinho Vieira. It was completed in 2024 and provided several novel insights into the study of Schelling models, such as delineating the limits of a mean‑field approach and demonstrating the existence of phase transitions when specific parameters are varied.

## The Schelling Models

In 1969, the American economist and 2005 Nobel Prize laureate in Economics, Thomas Schelling (1921–2016), published the article [Models of Segregation [1]](https://www.jstor.org/stable/1823701), in which, aiming to study racial segregation in American cities, he introduced one of the first mathematical models of social agents. The model consisted of two types of cells occupying predefined sites on a lattice; these cells could swap positions to increase their “satisfaction”, quantified by the number of like‑type cells in their eight‑cell Moore neighborhood.

![Schelling Diagram](https://github.com/user-attachments/assets/d559e7cc-6529-4eef-a231-e091b9a848f2)

*Figure 1. Schelling, T.C, 1974.*



## Extendind Schelling Models

Since Schelling's article, several studies have been conducted on Schelling models, exploring different dynamics, neighbourhoods and cells fractions. 

The model in this project works as follows:

- There are two types of agents (type 0 and type 1), whose fractions of the total are denoted by ρ₀ and ρ₁ (with ρ₀ + ρ₁ = 1).
- Each agent has a neighborhood of 4 agents (a von Neumann neighborhood).
- Satisfaction rules are encoded as a binary string s₀ s₁ s₂ s₃ s₄, where each sₖ is a Boolean indicating satisfaction (s = 1) or dissatisfaction (s = 0) when there are k same-type neighbors.
  For example, the rule 10001 means agents are satisfied with 0 or 4 same-type neighbors (s₀ = s₄ = 1) and dissatisfied with 1, 2, or 3 same-type neighbors (s₁ = s₂ = s₃ = 0). In this study the same satisfaction rule was applied to both agent types.
- Dynamics proceed in discrete steps: at each step, a randomly chosen dissatisfied type-0 agent swaps positions with a randomly chosen dissatisfied type-1 agent, even if both remain dissatisfied after the swap.

![Neighbourhood](https://github.com/user-attachments/assets/df4cbf7c-e325-401f-bb6a-6a9bb0781dfa) 

*Figure 2: Possible configurations for a Von Neumann neighbourhood with rules of the class s₀ s₁10 s₄.*


## The Mean Field Approach

One way to simplify the calculations of a Schelling model is by using the mean‑field approach, widely employed in Statistical Mechanics. This allows one to describe the system’s overall average behavior by treating each agent’s interactions as a single average interaction. In a Schelling model, agents seek satisfaction based on the number of same‑type neighbors they have. Under the mean‑field approximation, an agent’s satisfaction depends only on the average fraction of similar neighbors, thereby greatly simplifying the computations.

## Generating Differential Equations

The core element of the code is the computation of five 4×4 matrices. The calculation of these matrices is directly tied to the very definition of the system’s dynamics and the satisfaction rule being used. 

From these matrices, a set of differential equations is formulated, allowing for the numerical computation of the rate of change over time of ρₖ, the fraction of agents of a given type with k neighbors of type 1 in their neighborhood. Through numerical integration, like Runge-Kutta method, one can then analyze the evolution of these fractions as well as the total number of dissatisfied agents of each type.


## Project Files

- main.py: simulates the dynamic once for a specific input rule and fixed initial cell fractions;
- function.py: contains most of the main math based functions used;
- pcs.py: simulates the dynamic several times, changing the initial cell fractions to find critical points of phase transitions;
- deviations.py: simulates the dynamic several times, changing the initial cell fractions to investigate the dynamic when initial cell fractions are close to critical points.


## Bibliography:

[1] SCHELLING, Thomas. Models of Segregation. The American Economic
Review, Vol. 59, No. 2, 1969, pp. 488-493;

[2] André P Vieira et al. Dynamics of extended Schelling models J. Stat. Mech.
(2020) 013212

[3] André P Vieira et al. Phase transitions in a conservative game of life. PHYSICAL
REVIEW, E 103, (2021) 012132





