# Pilot-Wave-Dynamics
This repository is dedicated to the computational investigation of pilot wave theory, also known as de Broglie–Bohm theory. Through numerical simulations and visualizations, we explore the underlying mechanics of quantum phenomena as described by this interpretation. 

schrodingerequation_numerical.m

Explanation of the MATLAB Code for Solving the 2D Time-Dependent Schrödinger Equation
This MATLAB code simulates the evolution of a 2D Gaussian wave packet according to the time-dependent Schrödinger equation using a finite difference method. The Schrödinger equation is a fundamental equation in quantum mechanics that describes how the quantum state of a physical system changes over time.
Key Sections of the Code:
1. Initialization and Setup:
•	Clearing Workspace: The code starts by clearing the workspace using clear; and clc; to remove any previous data and output, ensuring a clean environment.
•	Grid and Spatial Parameters:
o	NP: Number of grid points along each axis (x and y).
o	xmax, xmin: Define the boundaries of the simulation area.
o	dx: Spatial step size.
o	x, y: Spatial grids in x and y directions.
o	X, Y: 2D meshgrid combining the x and y grids for easier matrix operations.
•	Physical Constants: The constants hbar (reduced Planck's constant) and m (particle mass) are set to 1 in atomic units for simplicity.
2. Time Step Calculation:
•	CFL Condition: The maximum allowable time step dt_max is calculated using the Courant–Friedrichs–Lewy (CFL) condition, critical for ensuring the stability of the numerical method.
•	Adjusted Time Step: The actual time step dt is set to 50% of this maximum value to further ensure numerical stability, a common practice in simulations.
3. Initial Wave Packet:
•	Wave Packet Parameters:
o	x0, y0: Center of the Gaussian in x and y.
o	rho_x, rho_y: Widths of the Gaussian in x and y.
o	kx, ky: Wavenumbers (related to the momentum) in x and y.
•	Construction: The initial wave packet is constructed as a 2D Gaussian modulated by plane waves in both x and y directions.
4. Potential Setup:
•	Free Particle Potential: A simple potential V is defined, set to zero everywhere, representing a free particle with no external forces acting on it.
5. Time Evolution Loop:
•	Crank-Nicolson Method: The wave function psi is evolved over time using the Crank-Nicolson method, a common implicit scheme used to solve time-dependent differential equations.
•	Kinetic Energy Operator: The kinetic energy operators T_x and T_y are calculated using finite differences, representing the second derivatives in the x and y directions.
•	Hamiltonian Operator: The Hamiltonian operator H_psi combines the kinetic energy terms with the potential.
•	Wave Function Update: The wave function psi is updated at each time step based on the Hamiltonian and the current state of psi.
•	Live Animation: The code updates the plot of the probability density |\psi(x, y, t)|^2 in real-time, providing a live visualization of the wave packet's evolution.
•	Video Saving: Simultaneously, the animation is saved as an MP4 video, allowing for later review and analysis. This is handled by the VideoWriter object, which captures and saves each frame.
6. Final Visualization:
•	Comparison: After completing the time evolution, the final wave packet is plotted and compared with the initial wave packet to visualize how the wave packet has evolved over time.
Tailored Explanation for Students:
This code provides a practical demonstration of how quantum mechanics can be simulated on a computer. By using a grid to represent space and evolving a wave packet over time, you can observe how quantum states change. The Crank-Nicolson method used here is stable and preserves the probability density, which is crucial in quantum simulations. The real-time animation and video saving features help visualize this evolution, making it easier to build an intuition for quantum particle behavior.
Tailored Explanation for Experts:
This MATLAB implementation demonstrates a finite difference approach to solving the 2D time-dependent Schrödinger equation with implicit periodic boundary conditions applied via the circshift function. The choice of a Gaussian wave packet allows for clear visualization of dispersion effects, and the Crank-Nicolson scheme ensures unitarity, maintaining the norm of the wave function over time. The time step is carefully chosen according to the CFL condition to ensure the stability of the simulation. The code is modular, allowing for easy modification to include different potentials or initial conditions, making it a robust tool for exploring quantum dynamics in two dimensions. The added functionality of live animation and video saving enhances its utility for both educational purposes and detailed analysis.
Suggestions for Next Steps:
•	Custom Potentials: Introduce different potentials (e.g., harmonic oscillator, double-well) to explore various quantum phenomena.
•	Adaptive Time-Stepping: Implement adaptive time-stepping to optimize computational efficiency and accuracy further.
•	Boundary Conditions: Experiment with different boundary conditions to study their effects on the wave packet's evolution.

https://github.com/user-attachments/assets/884a8566-455a-4342-961d-252074ff3aae

![wave_packet_evolution](https://github.com/user-attachments/assets/3c26034f-3120-4874-ad64-240cd0237aeb)


potentialbarriers.m

Explanation of the MATLAB Code for Visualizing Square Well, Elliptical, and Cylindrical Potential Barriers
This MATLAB code is designed to visualize three different types of potential barriers—square well, elliptical, and cylindrical (circular)—in a 2D grid. Potential barriers are important in quantum mechanics as they can represent regions where particles may be confined or reflected.
Key Sections of the Code:
1.	Initialization and Setup:
o	The number of grid points NP is defined, determining the resolution of the grid on which the potentials will be visualized.
o	The height of the potential barriers, potHeight, is set to a large value (10000), simulating an impenetrable barrier for a particle in quantum mechanical terms.
2.	Grid Definition:
o	A 1D grid x and y is created using linspace from 0 to 1 with NP points.
o	The meshgrid function generates a 2D grid X, Y, which allows for the creation of potential barriers in a 2D space.
3.	Square Well Potential Barrier:
o	The potential matrix V1 is initialized as a zero matrix, representing no potential.
o	A square well is created by setting a rectangular region within the grid (from X > 0.5 & X < 0.7 and Y > 0.2 & Y < 0.8) to the value of potHeight, creating a high potential barrier within that region.
4.	Elliptical Potential Barrier:
o	Similarly, the potential matrix V2 is initialized as a zero matrix.
o	An elliptical barrier is created by defining an ellipse centered at (0.5, 0.5) with a semi-major axis along the x-axis and a semi-minor axis along the y-axis. The ellipse is defined by the condition sqrt((X - 0.5).^2 + (0.5 * (Y - 0.5)).^2) < 0.1.
5.	Cylindrical (Circular) Potential Barrier:
o	The potential matrix V3 is again initialized as a zero matrix.
o	A cylindrical (circular) barrier is created by defining a circle centered at (0.5, 0.5) with a radius of 0.1. The circle is defined by the condition sqrt((X - 0.5).^2 + (Y - 0.5).^2) < 0.1.
6.	Visualization:
o	The imagesc function is used to display the potential barriers for each case:
	V1 (Square Well), V2 (Elliptical), and V3 (Cylindrical) are visualized in separate subplots within a single figure.
o	The axis xy command is used to set the y-axis direction so that the origin is at the bottom-left corner, matching typical Cartesian coordinates.
o	The colorbar function is added to each subplot to show the potential scale.
o	Titles and axis labels are added for clarity.
o	The overall figure is adjusted in size and layout using sgtitle to add a common title and set to modify the figure's dimensions.
Tailored Explanation for Students:
This code helps you understand how different types of potential barriers can be visualized in 2D space. Potential barriers are regions where particles encounter a "wall" that they cannot easily pass through. In this example, you see three different shapes of these barriers: a square, an ellipse, and a circle. The code uses simple mathematical conditions to create these barriers on a grid and then shows you how they look using color maps.
Each type of barrier can have different effects on particles in quantum mechanics, and this visualization gives you an intuitive sense of how these barriers might confine or restrict particles in different ways.
Tailored Explanation for Experts:
This MATLAB script efficiently generates and visualizes three canonical potential barriers commonly used in quantum mechanics and related fields. The potential barriers are represented on a uniform 2D grid using logical conditions to define the regions of high potential (potHeight). The square well and elliptical potentials are created using straightforward rectangular and elliptical inequalities, while the circular potential is defined using a radial distance condition.
The visualization approach leverages imagesc to render the potential landscapes, allowing for easy comparison between the different geometries. The script is structured to be easily extendable for more complex potential geometries or higher-dimensional studies. The overall layout is optimized for comparative analysis, making it a useful tool for both educational and research purposes.
Suggestions for Next Steps: a. Modify the code to simulate the dynamics of a particle interacting with these potential barriers.
b. Implement a user-defined function to create more complex or custom potential shapes.

![comparisonofbarriers](https://github.com/user-attachments/assets/31ef520f-5acf-49cb-9585-cdae5048f6d9)

gaussianwaveinteractingpotentialbarrier.m

Explanation of the MATLAB Code for Simulating 2D TDSE with Gaussian Wave Packet and Cylindrical Potential
This MATLAB code simulates the evolution of a 2D Gaussian wave packet interacting with a cylindrical potential barrier by solving the time-dependent Schrödinger equation (TDSE) using a finite difference method. The code also visualizes the wave packet's interaction with the potential and creates an animation of the time evolution.

Key Sections of the Code:
Initialization and Setup:

Grid and Wave Packet Parameters:

NP: Number of grid points, defining the resolution of the simulation.
rho: Initial width of the Gaussian wave packet.
x0, y0: Initial position of the wave packet on the grid.
kx, ky: Wavenumbers in the x and y directions, which define the initial momentum of the wave packet.
xmax, xmin, ymax, ymin: Define the spatial domain of the simulation.
dx: Spatial step size, determined by the number of grid points and the size of the domain.
dt: Time step size, calculated based on dx to maintain numerical stability.
lambda: CFL-like parameter, ensuring stability in the numerical solution.
Potential Parameters:

ellip: Ellipticity parameter; set to 1 to create a cylindrical potential barrier.
potheight: Height of the potential barrier, representing an impenetrable region.
Grid Initialization:

The spatial grid x, y is created using linspace, and the 2D meshgrid X, Y is generated using meshgrid to facilitate matrix operations.
Wave Function Initialization:

The initial wave function Inwave is defined as a Gaussian wave packet, centered at (x0, y0) with widths defined by rho, and modulated by plane waves with wavenumbers kx and ky.
Potential Barrier Initialization:

The potential matrix V1 is initialized to represent a cylindrical barrier centered at (0.55, 0.55) with a radius of 0.05. The arrayfun function applies the cylindrical potential condition across the grid.
Time Evolution:

The main time evolution loop iteratively updates the wave function according to the TDSE:
Alpha and Gamma Coefficients: These are intermediate variables used to propagate the wave function using a tridiagonal matrix approach. The coefficients are calculated during two passes over the grid, one for each spatial dimension.
Forward Passes: The wave function is updated in two forward passes, first to propagate phi and then psi. The potential barrier is applied to the wave function between these passes, introducing the interaction between the wave packet and the cylindrical barrier.
Wave Function Update: The updated wave function phi is stored at each time step for later visualization.
Visualization:

The code visualizes the wave packet's probability density |\psi|^2 at four distinct time steps. The plots include the potential barrier, allowing for a clear visualization of the wave packet's interaction with the barrier.
Animation Creation:

An optional section of the code creates an animated GIF of the wave packet's time evolution. The imagesc function is used to generate frames, and imwrite writes these frames to a GIF file.
Tailored Explanation for Students:
This code simulates how a quantum particle, represented as a Gaussian wave packet, behaves when it encounters a cylindrical barrier in a 2D space. The wave packet starts with a certain position and momentum, and the code shows how it evolves over time when it hits a barrier. You can think of this as simulating how a water wave would interact with an obstacle, but in a quantum context.

The code also visualizes the simulation results, allowing you to see how the wave packet spreads out and interacts with the barrier. Additionally, it creates an animation to help you better understand the time evolution of the system.

Tailored Explanation for Experts:
This MATLAB implementation solves the 2D TDSE for a Gaussian wave packet interacting with a cylindrical potential using a finite difference method with an explicit time-stepping scheme. The simulation employs a tridiagonal matrix approach to efficiently handle the wave function propagation, with stability maintained via a CFL-like parameter lambda. The potential is applied between directional sweeps, ensuring an accurate representation of the wave packet's interaction with the barrier.

The code is structured to facilitate both visualization and further analysis, providing detailed snapshots of the wave packet at multiple time steps and offering an option to generate a full animation. The use of array operations and meshgrid ensures scalability and efficiency, making the script suitable for more complex simulations with different potential configurations.

Suggestions for Next Steps:
a. Modify the initial conditions (e.g., wave packet width or momentum) to explore different scattering behaviors.
b. Implement absorbing boundary conditions to reduce reflections from the grid edges, simulating an open system.

![2DCylindricalPotential3D](https://github.com/user-attachments/assets/93431e58-58f5-46fd-be62-6585786e1ab4)

![wave_function_interaction_barrier](https://github.com/user-attachments/assets/61578d33-774f-40c4-a17f-f944ae9ec4e5)


