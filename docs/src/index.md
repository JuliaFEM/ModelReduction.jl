# ModelReduction.jl documentation

```@contents
```

ModelReduction.jl is a Julia package to perform model reduction methods for i.e. multibody dynamics problems. The packcage includes model order reduction methods such as the Guyan reduction and the Craig-Bampton method.

Reducing the sizes of stiffness and mass matrices of the model will greatly decrease the computation resources needed when performing dynamic analyses.

ModelReduction.jl is a part of JuliaFEM. All codes are MIT licensed.

## Installing and testing the package

Download and unzip the packge files into your Julia folder and then install the package the same way other Julia packages are installed.

```julia
julia> Pkg.add("ModelReduction")

```

Test the package with ```Pkg.test``` etc.

```julia
julia> Pkg.test("ModelReduction")

```


## Usage example

This example demonstrates how to use the Craig-Bampton method function.

Problem setup:

```julia
K =
[2 -1 0 0;
-1 2 -1 0;
 0 -1 2 -1;
 0 0 -1 1]

M =
[2 0 0 0;
 0 2 0 0;
 0 0 2 0;
 0 0 0 1]

r = [4]
l = [1, 2, 3]
n = 1

```
K = original stiffness matrix, M = original mass matrix, r = retained DOF:s, l = internal DOF:s, n = the number of the internal modes to keep.
Calculate the reduced mass and stiffness matrices Mred and Kred.

```julia
craig_bampton(K, M, r, l, n)

# output

Mred, Kred = ([2.75 -1.20711; -1.20711 1.0], [0.25 0.0; 0.0 0.292893])

```


## Contributing



Have a new great idea and want to share it with the open source community? 
From [here](https://github.com/JuliaLang/julia/blob/master/CONTRIBUTING.md)
and [here](https://juliadocs.github.io/Documenter.jl/stable/man/contributing/)
you can look for coding styles. [Here](https://docs.julialang.org/en/stable/manual/packages/#Making-changes-to-an-existing-package-1) it is explained how to contribute to 
open source project, in general.

