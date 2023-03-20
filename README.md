
# <div align="center">Seminar "Simulation of electrical machines" <br> SoSe 23 <br> v1.2</div>


## Model problem

We want to solve the nonlinear magnetostatics problem 

$$\begin{aligned}
\operatorname{curl}h &= j, \qquad  \text{in } \Omega, \\
\operatorname{div}b &= 0, \qquad \text{in } \Omega, \\
n \cdot b &= 0, \qquad \text{on } \partial\Omega, 
\end{aligned}$$

To complete the model, a general material law is
provided in the form $$h = h(b)$$

The precise form of this relation for the different
materials will be specified below.

## Model specifications

<div align="center"><img src="https://raw.githubusercontent.com/radu-bogdan/SeminarEM/main/geometry.png" width="500" height="500">

<b>Figure [1] - Geometry of the motor</b>
</div>


Let us describe the exact geometry, see Figure [1], which is comprised of several materials:

-   16 permanent magnets $\Omega_{M_i}$, depicted by the blue areas $$h = \nu_0 b - m_i,$$ where $\nu_0 = 10^7/(4\pi)$ is the vacuum reluctivity. The magnetization of the individual magnets is chosen as $|m_i|=1$ and directed perpendicular (close to axial direction) to their elongation direction; the two magnets of each pair point into the same direction (in/outwards), and the direction changes between neighbouring pairs.

-   48 slots (windings) $\Omega_{C_i}$, depicted by the green areas, filled by copper, with $$h = \nu_0 b$$ The slots carry the windings with impressed current densities $j_3$, which are different constants for each slot; specified in files available at git.

-   The iron stator $\Omega_{\text{S}}$, iron rotor $\Omega_{\text{R}}$ (red areas), $$h = \nu(|b|)b.$$ We use is the Brauer model where $\nu(|b|) = k_1e^{k_2|b|^2}+k_3$ and the parameters $$k_1 = 49.4,\quad k_2 = 1.46,\quad k_3 = 520.6,$$ from [@Brauer1975], which represent cast iron.

-   The iron shaft $\Omega_{\text{SH}}$ (orange) is supposed to be made up of non-conductive material, and hence $$h = \nu_0 b.$$

-   48 air gap regions $\Omega_{AC_i}$ located below the coils;

    32 air gap regions $\Omega_{AM_i}$; 2 around each of the 16 magnets depicted by the purple areas; 
    
    3 air gap regions between stator and rotor $\Omega_{AG_i}$. 
    
    For all air-gaps, we use the linear relation $$h = \nu_0b$$
    
Also check Peter's thesis [@Gangl2017 Chapter 2] for a detailed description; a good resource for the material laws is the book by Stratton, see [@Stratton1941 Sec. 1.6].

## The repository

The geometry information and the exact values of the right-hand sides are given in this repository and are available for either $\texttt{Python}$, $\texttt{Matlab}$ or $\texttt{NGSolve}$. Go to <https://www.radubogdan.de/fem/motor.html> if you wish to look at the geometry (and the mesh) in more detail.

## Solution of the vector potential formulation

My solution to the vector potential formulation can be viewed here: https://www.radubogdan.de/fem/solution_magnetostatics.html. You can compare point evaluations to check if you obtained the same answer.

## Data format

We use matlab notation to explain the main data provided:

         p          ... 2 x np array with vertex coordinates
         e          ... 3 x ne array with edge points plus subboundary indices 
         t          ... 4 x nt array with triangle points plus subdomain indices
         regions_2d ... cell array with subdomain descriptinos
         regions_1d ... cell array with boundary desciptions
         j3         ... 1 x ncoils array with current densities
         m          ... 2 x nmagnets array with (per magnet) magnetization values
         
## Questions? Found any mistakes?
Come to my office, you know where to find me!

## List of changes
v1.2 : added geometry files for $\texttt{Matlab}$ (Thank you Josef) and an IGS format
       added solution to the magnetostatic problem for comparison
       
v1.1 : updated magnetization values

## Bibliography

[@Brauer1975] J. Brauer. Simple equations for the magnetization and reluctivity curves
of steel. , 11(1):81--81, 1975.

[@Gangl2017] Peter Gangl. Topology and Shape Optimization with Application to Electrical
  Machines, volume 43 of *Schriftenreihe Advances in Mechatronics*.
Trauner Verlag, Austria, 2017.

[@Stratton1941] Julius Adams Stratton. . McGraw-Hill Companies (New York), 1 edition,
1941. The URL and misc. info are for a re-issue published in 2007 by
Wiley-IEEE Press.
