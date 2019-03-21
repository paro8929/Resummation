##########################################
Numerical Solvers of scalar phi^4 in 1+1 dimensions
R2 and R3 approximations

by Paul Romatschke, January 2019

If you are using (parts of) this code, it would be great
if you can reference the corresponding paper(s):

https://arxiv.org/pdf/1901.05483.pdf

File list:

* 2d-R2.cpp [main part of the R2-level solver]
* 2d-R3.cpp [main part of the R3-level solver]
* 2d-R3.h   [header file for the R3-level solver (vertex class)]
* compile.sh [bash shell script for compiling executable]
* linear_int.cpp [linear interpolation routines]
* linear_int.h   [linear interpolation header file]
* mf-2d-R2.cpp [measurement routines for R2-level solver]
* mf-2d-R2.cpp [measurement routines for R3-level solver]
* rootfinding.cpp [rootfinding routines]
* structs.h [header file for data structures]
* R2.exe [executable R2-level]
* R3.exe [executable R3-level]

Directories:
data [output for fields for each iteration (for debugging, mostly)]
LegendreData [stencils and weights for Gauss-Legendre quadrature]
store [main output directory]
DataFiles [Data text files matching published manuscript]

---------------------------------------

Setting parameters: modify the files R2.cpp (R3.cpp)

- Value of interaction: change the parameter "interaction" (which is g=lambda/mR^2) to an appropriate value

- Resolution of integration: change the parameter "N0"  to an appropriate value; note that corresponding to your choice of N0, an appropriate file "Leg-Quad-NXX.dat" must be present in the "LegendreData" subdirectory where xx=N0-1
  (R3 only: change the parameter "N1" to an appropriate value for adjusting the resolution of the angular interpolation in the vertex function)

- Maximum number of iterations: change the parameter "mIT" to an appropriate value

- Maximum value of K^2 in self-energies: change the parameter "cut" to an appropriate value

--------------------------------------

Compiling the code: use ./compile

Output: R2.exe (R3.exe)

Note: shell script assumes that gsl and openmp libraries are present and installed at your machine, and your compiler is c++11 standard compliant.


Single thread compiling: It is possible to remove dependence on openmp by uncommenting all "#pragma" pre-processor commands in the code if you insist on running on single thread


-------------------------------------

Running the code and inspecting results

Use ./R2.exe (./R3.exe) to run the code

Code runs multiple iterations (governed by "mIT" parameter, see above) to converge to correct results for self-energies and vertices. Code outputs progress to std and results to "history" (and "final") files in the "store" subdirectory.

---------------------------------------

Output files: "history-xxx.dat"

Columns correspond to iteration number, delta m^2, delta Pi(K=0), Lambda/mR^2, M/mR, c2 and error in c2. If successfully converged, the content of rows will be identical after a sufficient number of iterations
