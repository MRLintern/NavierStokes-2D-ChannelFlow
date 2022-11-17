
## 2D Navier-Stokes Equations Solver for 2D Chanel Flow
This code solves the Navier-Stokes Equations for 2D Open Chanel Flow.
The equations are discretized using the finite volume method.
The Resulting System of Algebraic Equations is solved using the Successive Over-Relaxation (SOR) Method

## Operating System
`Ubuntu 20.04`.

## Make
This project uses `make` to generate the executable: [Get Make from the GNU Project](https://www.gnu.org/software/make/).
You may need to install this: 
1. Updates: `$ sudo apt-get update`
2. Install: `$ sudo apt-get install make`


## Compiler
The compiler used for this software is `gfortran`: [Get gfortran from GNU](https://gcc.gnu.org/fortran/).
To install it: 
1. Updates: `$ sudo apt-get update`
2. Install: `$ sudo apt-get install gfortran-9`

## Running the Application
1. Clone the repository: `$ git clone https://github.com/MRLintern/NavierStokes-2D-ChanelFlow.git`
2. `$ make`
3. `$./NavierStokes-2D-ChanelFlow`

## Results
Once the program has been run, a file called `RESULTS.dat` will be generated
with the solution. For visualization, [ParaView](https://www.paraview.org/) is a good choice.
If you need help go to [ParaView Tutorial](https://www.paraview.org/Wiki/images/b/bc/ParaViewTutorial56.pdf).
If you don't want to run the program, `RESULTS.dat` is included.
The original source code generated the solution to a .csv file, called `NS_RESULTS.csv` which is still available.
If you want to run the program and generated a .csv file, replace the file extension from .dat to .csv.
