# ModelReduction.jl
[![DOI](https://zenodo.org/badge/97968807.svg)](https://zenodo.org/badge/latestdoi/97968807)[![Build Status](https://travis-ci.org/JuliaFEM/ModelReduction.jl.svg?branch=master)](https://travis-ci.org/JuliaFEM/ModelReduction.jl)[![Coverage Status](https://coveralls.io/repos/github/JuliaFEM/ModelReduction.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaFEM/ModelReduction.jl?branch=master)[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliafem.github.io/ModelReduction.jl/stable)[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliafem.github.io/ModelReduction.jl/latest)[![Issues](https://img.shields.io/github/issues/JuliaFEM/ModelReduction.jl.svg)](https://github.com/JuliaFEM/ModelReduction.jl/issues)

ModelReduction.jl is a Julia package to perform model reduction methods for i.e. multibody dynamics problems. The packcage includes model order reduction methods such as the Guyan reduction and the Craig-Bampton method.

Reducing the sizes of stiffness and mass matrices of the model will greatly decrease the computation resources needed when performing dynamic analyses.

ModelReduction.jl is a part of JuliaFEM. All codes are MIT licensed.

## Installing and testing the package

Install the package the same way other Julia packages are installed.

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
julia> K = [2 -1  0  0;
           -1  2 -1  0;
            0 -1  2 -1;
            0  0 -1  1]

julia> M = [2 0 0 0;
            0 2 0 0;
            0 0 2 0;
            0 0 0 1]

julia> r = [4]
julia> l = [1, 2, 3]
julia> n = 1

```
K = original stiffness matrix, M = original mass matrix, r = retained DOF:s, l = internal DOF:s, n = the number of the internal modes to keep.
Calculate the reduced mass and stiffness matrices Mred and Kred.

```julia
julia> using ModelReduction

julia> Mred, Kred = ModelReduction.craig_bampton(K, M, r, l, n)
([2.75 -1.20711; -1.20711 1.0], [0.25 0.0; 0.0 0.292893])

```

## Citing

If you like using our package, please consider citing our [article](https://rakenteidenmekaniikka.journal.fi/article/view/69026/35912):
```
@article{rapo2018implementing,
 title={Implementing model reduction to the JuliaFEM platform},
 volume={51},
 url={https://rakenteidenmekaniikka.journal.fi/article/view/69026},
 doi={10.23998/rm.69026},
 number={1},
 journal={Rakenteiden Mekaniikka},
 author={Rapo, Marja and Aho, Jukka and Koivurova, Hannu and Frondelius, Tero},
 year={2018},
 pages={36-54}
}
```


## Contributing

Have a new great idea and want to share it with the open source community? 
From [here](https://github.com/JuliaLang/julia/blob/master/CONTRIBUTING.md)
and [here](https://juliadocs.github.io/Documenter.jl/stable/man/contributing/)
you can look for coding styles. [Here](https://docs.julialang.org/en/stable/manual/packages/#Making-changes-to-an-existing-package-1) it is explained how to contribute to 
open source project, in general.


