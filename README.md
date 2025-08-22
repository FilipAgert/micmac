# Macroscopic-Microscopic Nuclear Model

A **work-in-progress** program for calculating nuclear ground state masses, deformation, fission barrier heights and deformation using a macroscopic-microscopic method.
There is likely an error in the implementation of the shell model as the spectra does not look right.

Microscopic corrections are from the Strutinsky method on single particle levels from a Wood-Saxon + Coulomb + spin-orbit Hamiltonian, diagonalised with an axially deformed harmonic oscillator basis.
The microscopic method attempts to replicate [1], while the implementation uses calculations from [2]. Quadrature nodes and weights are computed using [this](https://github.com/FilipAgert/fquad/) program.

> [1] **P. Jachimowicz, M. Kowal, and J. Skalski.**
>     *Properties of heaviest nuclei with 98≤Z≤126 and 134≤N≤192,*
>     At. Data Nucl. Data Tables 138, 101393 (2021).
  
> [2] **H. C. Pauli.**
> *On the shell model and its application to the deformation energy of heavy nuclei,*
> Physics Reports 7, 35 (1973)



---
## Requirements
This program uses the libraries LAPACK and BLAS for diagonalisation and solving linear systems. These need to be located on your computer. 
In the makefile, change `DLIB` to the directory where these libraries are contained.


## Build Instructions
### Compile shell model program
```bash
make sh
```
which builds the executable at:
```
app/shell
```
### Compile Strutinsky shell correction program
```bash
make st
```
which builds the executable at:
```
app/strut
```

---
## Running the program

### Run shell model calculation
The shell model program is run with the following syntax:
```bash
app/shell <Z> <A> <beta2> <beta4>
```
where
- `<Z>`: Proton number of nucleus
- `<A>`: Mass number of nucleus
- `<beta2>`: Quadrupole deformation of nucleus (optional, defaults to 0)
- `<beta4>`: Hexadecapole deformation of nucleus (optional, defaults to 0)

This generates two files containing the proton and neutron single-particle levels
```
data/out/<Z>_<A>/levels_n.dat
data/out/<Z>_<A>/levels_p.dat
```
which contains single particle levels for protons and neutrons. 
Each file contains three columns:
| Column | Description        |
|--------|--------------------|
| 1      | m_j |
| 2      | parity             |
| 3      | Energy (MeV)       |

m_j is in half integer format such that value 3 in column one means m_j = 3/2.
Parity is 0 in column 2 if parity is not a good quantum number with the current settings. This is only if `beta3`, `beta5` or `beta7` is nonzero.

The format for reading the files is: `(i3,1x,i3,1x,f10.5)`.

Note: Each level has degeneracy 2.

### Run Strutinsky shell correction
The strutinsky shell correction program is run with the following syntax:
```bash
app/strut <Z> <A> <beta2> <beta4>
```
where
- `<Z>`: Proton number of nucleus
- `<A>`: Mass number of nucleus
- `<beta2>`: Quadrupole deformation of nucleus (optional, defaults to 0)
- `<beta4>`: Hexadecapole deformation of nucleus (optional, defaults to 0)

This outputs the shell correction to console for the provided nucleus.
## Settings and parameters
In the file src/settings.f90 the settings for the shell model can be found. These are the physical parameters such as potential depth, or the technical parameters, such as `N_max` for the maximum harmonic oscillator shell to include in the states. From all the states up to `N_max` in the axially deformed harmonic oscillator, we pick the lowest energy states to include in the diagonalisation. The number of states to include is decided by the parameters `num_p_states` and `num_n_states`. 
## TODO
### Microscopic
- Verify single particle levels to verify it's implemented correctly.
- Verify strutinsky shell correction and pairing correction is implemented correctly.
- Implement higher order deformations `beta6`, `beta8` and reflection symmetry breaking `beta3`, `beta5`, `beta7`.

### Macroscopic
- Method for finding fission barrier in multidimensional deformation space.
- Incorporate Strutinsky shell correction program into microscopic correction instead of current simplified shell correction calculation.
- Implement higher order deformations.
- Change the Coulomb and surface macroscopic deformation energy from and expansion to performing the integral directly. 
