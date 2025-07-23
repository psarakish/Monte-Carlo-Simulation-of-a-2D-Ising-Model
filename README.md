# Monte Carlo Simulation of the 2D Ising Model

This project simulates a 2D Ising model using the Metropolis Monte Carlo algorithm. The goal is to study how temperature affects the magnetization and energy of a spin lattice system, as well as calculate physical quantities like heat capacity and magnetic susceptibility.

## 🧪 Overview

- Simulates a square \( N \times N \) lattice of spin-1 particles (\( s = \pm1 \))
- Each spin interacts with its four nearest neighbors
- Uses the Metropolis algorithm to evolve the system toward thermal equilibrium
- Calculates temperature-dependent properties:
  - Total energy
  - Magnetization
  - Heat capacity
  - Magnetic susceptibility

## ⚙️ Parameters

- Lattice sizes: `N = 10, 20, 30`
- Temperature range: `T ∈ (0, 4)` with step `ΔT = 0.1`
- Number of Monte Carlo steps: `10,000`
- Magnetic field: `0` (i.e., zero-field case)
- Spin value: `±1`

## 📊 Results

The simulation demonstrates a phase transition around the critical temperature. As temperature increases, spontaneous magnetization decreases, and the system transitions from an ordered to a disordered state. The larger the grid, the sharper the transition.

## 💻 Implementation

- Language: C++
- Random number generation for thermal noise
- Uses periodic boundary conditions
- Observables are averaged over Monte Carlo steps for each temperature
