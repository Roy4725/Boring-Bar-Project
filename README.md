# Anti-Vibration-Boring-Bar-Project
Vibration Damping in a Two-Mass-Forced System – Simulation & Optimization
This project simulates a mechanical system composed of two connected masses with springs and dampers, analyzing their dynamic response to an external force. It also includes parameter optimization to minimize undesired vibrations.

🔧 Key Files and Code Snippets
1. optimize_parametrs.py – Geometry & Mechanical Parameter Optimization
Performs optimization of internal geometry (e.g., inner diameter, auxiliary mass dimensions) to minimize the vibration amplitude of the main body.

Simulation is done using numerical solution of coupled differential equations (Euler method). The goal is to find parameter combinations that minimize maximum displacement.

Key Code Sections:
-  [Inertia, mass, stiffness, damping calculation](https://github.com/troy-piz/vibration-damping/blob/main/optimize_parametrs.py#L27-L34)
-  [Numerical simulation loop](https://github.com/troy-piz/vibration-damping/blob/main/optimize_parametrs.py#L71-L78)
-  [Optimal configuration logic](https://github.com/troy-piz/vibration-damping/blob/main/optimize_parametrs.py#L95-L97)


2. damper.py – Damper Stiffness Calculation Based on Geometry
Uses theoretical formulas to calculate damper stiffness (K) based on length, inner and outer radii, and viscosity coefficient.
Includes two approaches to find the damper length that achieves the target K value:

Method 1: Optimization by minimizing error.
Method 2: Root-finding to solve K(L) = Target K.

Key Code Sections:
-  [Theoretical calculation functions for S, β, K](https://github.com/troy-piz/vibration-damping/blob/main/damper.py#L12-L28)
-  [Optimization to find L such that K ≈ target_K](https://github.com/troy-piz/vibration-damping/blob/main/damper.py#L35-L37)
-  [Plotting K vs L graph with annotations](https://github.com/troy-piz/vibration-damping/blob/main/damper.py#L83-L104)


3. graph.py – System Response Simulation With and Without Damper
This script simulates the time response of the two-mass system with and without the auxiliary damping mass. It compares displacement to show how the damper reduces vibration amplitude.

Key Code Sections:
-  [Define sinusoidal external force input](https://github.com/troy-piz/vibration-damping/blob/main/graph.py#L27)
-  [Numerical simulation of bar + damper system using Euler’s method](https://github.com/troy-piz/vibration-damping/blob/main/graph.py#L44-L51)
-  [Plotting displacement with and without auxiliary mass](https://github.com/troy-piz/vibration-damping/blob/main/graph.py#L74-L86)


🧠 Technologies Used
Python: numpy, matplotlib, scipy, pandas
Numerical ODE solving (Euler method)
Optimization and root-finding techniques
Visualization and displacement comparison
