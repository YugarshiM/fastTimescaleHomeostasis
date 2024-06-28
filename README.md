# Burster Assembly and Analysis

This code is designed to generate bursters, measure their activity, and plot the half-(in)activations and maximal conductances as the model assembles a burster that satisfies the calcium targets. Follow the steps below to run the simulations and generate the plots:

## Steps to Run the Code

1. **Run `main.jl` in Julia**
   - Ensure the required packages have been installed.
   - This script produces two simulation files in the `./simulation` directory.
   - The first file shows the assembly of the burster from specified initial conditions.
   - At the end of this simulation, the homeostatic mechanism is turned off, and the result is saved in the file with the suffix `_RESTART`.

2. **Run `plotOutput.m` in MATLAB**
   - Once the files from the previous step are generated, run this script in MATLAB.
   - This script measures various aspects of the activity, including Period, Interburst Interval, Burst Duration, Max Hyperpolarization, Spike Height, and Slow Wave Amplitude, along with their Coefficients of Variation (CVs).
   - The measurements are output to the terminal.
   - A plot is generated showing how the model changes the half-(in)activations and maximal conductances in time to create a burster that satisfies the calcium targets

## Output

- **Simulation Files**: Generated in the `./simulation` directory.
- **Activity Measurements**: Displayed in the terminal.
- **Plots**: Showing time evolution of half-(in)activations and maximal conductances.

Make sure to have both Julia and MATLAB installed, with the necessary packages, to run the respective scripts.
