# Julia Workshop @ YoMos 2025 

This repository contains material for the workshop. 

## Installing Julia 
You can install an up-to-date version of Julia by following the instructions on
the official webpage <https://julialang.org/install/>.

If you use VSCode or VSCodium, you can install an extension that brings a lot
of nice features and language tooling to the editor:
<https://www.julia-vscode.org/>. 
However, Julia is supported by other editors as well (see the 'Editors and
IDEs' section at <https://julialang.org/>).

## Running Code 
In order to install the packages needed for this workshop, pull this repository
and navigate to the directory. 
Start a Julia session there (or navigate to the directory using the function `cd("<path-to-dir>")`). 
Then, type `]` in the REPL to activate the package mode. 
Since this directory contains a project and manifest file, you can install the
necessary packages by simply typing 
```julia REPL
(@v1.11) pkg> activate .
```
to activate the project environment, followed by 
```julia REPL
(JuliaWorkshop_YoMos2025) pkg> instantiate 
```
The latter command obtains and installs the packages that belong to the
environment defined by the project and manifest file.

## Workshop Content 
The practical part of the workshop consists of three parts:
1. The fundamentals of Julia (with focus on multiple dispatch and its features)
2. Numerical Simulations of ODEs and plotting; an example of how easy it is to integrate different packages in Julia
3. Building an interactive GUI that displays numerical solutions of an ODE
   system (with sliders for the parameters and initial condition)

You find the source code in the file `yomos.jl`. 
The static html file `yomos.html` is generated with `Weave.jl` from the Julia
markdown file `yomos.jmd`.  
This captures the output of the Julia code chunks, so have a look at this file
if you don't want to install Julia.
