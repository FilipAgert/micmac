# Macroscopic-Microscopic Nuclear Model

A **work-in-progress** program for calculating nuclear ground state masses, deformation, fission barrier heights and deformation using a macroscopic-microscopic method.

Microscopic corrections are from the Strutinsky method on single particle levels from a Wood-Saxon + Coulomb + spin-orbit Hamiltonian, diagonalised with an axially deformed harmonic oscillator basis.
The microscopic method attempts to replicate [1], while the implementation uses calculations from [2]. Quadrature nodes and weights are computed using [this](https://github.com/FilipAgert/fquad/) program.

> [1] **P. Jachimowicz, M. Kowal, and J. Skalski**
>     *Properties of heaviest nuclei with 98≤Z≤126 and 134≤N≤192*
>     At. Data Nucl. Data Tables 138, 101393 (2021).
  
> [2] **H. C. Pauli**
> *On the shell model and its application to the deformation energy of heavy nuclei*
> Physics Reports 7, 35 (1973)



---
## Requirements
This program uses the libraries LAPACK and BLAS for diagonalisation and solving linear systems. These need to be located on your computer. 
In the makefile, change `DLIB` to the directory where these libraries are contained.


## Build Instructions

Use `make` to compile different parts of the code:
### Compile macroscopic + microscopic fitting program
```bash
make
```
### Compile shell model program
```bash
make sh
```
### Compile Strutinsky shell correction program
```bash
make st
```
---
## Running the program
### Run macroscopic-microscopic calculation
```bash
make run
```

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

This generates a file containing the proton and neutron single-particle levels
```
data/out/levels_<Z>_<A>.dat
```
which contains single particle levels for protons in column 1 and neutrons in proton 2.

The format for reading the levels is: `f10.5,3x,f10.5`

Note: Each level has degeneracy 2

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
