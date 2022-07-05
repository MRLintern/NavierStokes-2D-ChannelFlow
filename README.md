## 2D Navier-Stokes Equations Solver for 2D Chanel Flow
This code solves the 2D Navier-Stokes equations.
The equations are discretized using the finite volume method.
The System of Equations is solved using the Successive Over-Relaxation (SOR) Method

## Operating System
Operating system used is Ubuntu 20.04.

## Make
This project uses Make to generate the executable: [Get Make from the GNU Project](https://www.gnu.org/software/make/).
You may need to install this: 
1. Updates: `$ sudo apt-get update`
2. Install: `$ sudo apt-get install make`


## Compiler
The compiler used for this software is gfortran: [Get the gcc compiler from GNU](https://gcc.gnu.org/fortran/).
To install it: 
1. Updates: `$ sudo apt-get update`
2. Install: `$ sudo apt-get install gfortran-9`

## Instructions
1. Clone the repository: `https://github.com/MRLintern/NavierStokes-2D-ChanelFlow.git`
2. `$ make`
3. `$./NavierStokes-2D-ChanelFlow`

## Results
Once the program has been run, a file called NS_RESULTS.csv will be generated
with the solution. For visualization, [ParaView](https://www.paraview.org/) is a good choice.
If you need help go to [ParaView Tutorial](https://www.paraview.org/Wiki/images/b/bc/ParaViewTutorial56.pdf).
If you don't want to run the program, NS_RESULTS.csv is included.
