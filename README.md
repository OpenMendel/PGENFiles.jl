# PGEN.jl
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://OpenMendel.github.io/PGEN.jl/dev)
<!--[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://OpenMendel.github.io/PGEN.jl/stable)-->
[![codecov](https://codecov.io/gh/OpenMendel/PGEN.jl/branch/main/graph/badge.svg?token=W28QPREGC7)](https://codecov.io/gh/OpenMendel/PGEN.jl)
[![build Actions Status](https://github.com/OpenMendel/PGEN.jl/workflows/CI/badge.svg)](https://github.com/OpenMendel/PGEN.jl/actions)

Routines for reading compressed storage of genotyped or imputed markers

[*Genome-wide association studies (GWAS)*](https://en.wikipedia.org/wiki/Genome-wide_association_study) data with imputed markers are often saved in the [**PGEN format**](https://www.cog-genomics.org/plink/2.0/input#pgen) in `.pgen` file.
It can store both hard calls and imputed data, unphased genotypes and phased haplotypes. Each variant is compressed separately. This is the central data format for [PLINK 2](https://www.cog-genomics.org/plink/2.0/). 

## Installation

This package requires Julia v1.6 or later, which can be obtained from
https://julialang.org/downloads/ or by building Julia from the sources in the
https://github.com/JuliaLang/julia repository.


This package is registered in the default Julia package registry, and can be installed through standard package installation procedure: e.g., running the following code in Julia REPL.
```julia
using Pkg
pkg"add https://github.com/OpenMendel/PGEN.jl"
```

## Citation

If you use [OpenMendel](https://openmendel.github.io) analysis packages in your research, please cite the following reference in the resulting publications:

*OPENMENDEL: a cooperative programming project for statistical genetics. Zhou H, Sinsheimer JS, Bates DM, Chu BB, German CA, Ji SS, Keys KL, Kim J, Ko S, Mosher GD, Papp JC, Sobel EM, Zhai J, Zhou JJ, Lange K. Hum Genet. 2019 Mar 26. doi: 10.1007/s00439-019-02001-z. [Epub ahead of print] PMID: 30915546*

## Acknowledgments

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
