# SCUC
SCUC is a  MIP(Mixed Integer Programming) solver for security constrained unit commitment problem. It aims to develop a customized approach that can include  specialized techniques such as dynamic cut and column generation, specialized branching, or heuristics to help find feasible solutions more quickly. Benefit from knowledge about the problem being solved, the solver is expected to solve the large scale day-ahead security constrained unit commitment performance more efficiently than a general MIP solver such as CPLEX, Gurobi.

SCUC is written in C++ and is released as open source under the [Eclipse Public License 2.0](https://opensource.org/licenses/EPL-2.0).

## Dependencies 
* [CBC](https://github.com/coin-or/Cbc) a general mixed integer programming solver 
* [Clp](https://github.com/coin-or/Clp) the default solver for LP relaxations 
* [Cgl](https://github.com/coin-or/Cgl) for cut generation
* [CoinUtils](https://github.com/coin-or/CoinUtils) for reading input files and various utilitiesBranching

## Usage
   ### Binaries
   The   binaries of dependencies for most platforms are available as part of  SCUC.
   * Windows MSVisualStudio: The static libraries for x64 and x86 platforms used in Microsoft Visual studio for debug mode and release mode are in the win_lib  folder.
   
   * Windows Msys2: The static libraries used in Msys2 shell are in the  mingw_lib folder.
   
   * Linux: On Debian/Ubuntu, the dynamic libraries  are in folder.
   
   ### Build
   
   * With Microsoft Visual Studio: 
   
   The easiest way to build SCUC on Windows is through MSVisualStudio.
   For Microsoft Visual C++ users, there are project files for version 10 and version 19 available in the MSVisualStudio directory. First, obtain the source code using either a Windows git client or download a snapshot. In MSVC++, open the solution file (this should be converted to whatever version of MSVC+ you are using) and build the SCUC project. The code should build out of the box with default settings.
   
   Assumptions:

   - A VS solution with all necessary libs (libCbc, libClp, libCbcSolver, libCgl, libCoinUtils, libOsi, libOsiCbc, libOsiClp). The libraries files can be found inside the win_lib folders.
   
   Steps (based on VS 2019):

1.  add `..\..\..\include\coin;` under Properties -> Configuration Properties -> C/C++ -> General -> Additional Include Directories 
   
2.  add `..\..\..\win_lib\Vxx_VSxxxx\x??-vxx-Debug(Release);` under Properties -> Configuration Properties -> Linker-> General -> Additional Library Directories 
   
3.  add   `libCbc.lib;libCbcSolver.lib;libCgl.lib;libClp.lib;libCoinUtils.lib;libOsi.lib;`    `libOsiCbc.lib;libOsiClp.lib;`  under Properties -> Configuration Properties -> Linker ->
   Input ->  Additional Dependencies 
   
   * With Msys2: 
   
   A Makefile is in src folder. Users can modify it in terms of the files which will be built.
   
   * For Linux: 
   
   A Makefile_Linux is in src folder. Users can modify it in terms of the files which will be built.
   
 
