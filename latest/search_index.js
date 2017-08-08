var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#ModelReduction.jl-documentation-1",
    "page": "Introduction",
    "title": "ModelReduction.jl documentation",
    "category": "section",
    "text": "ModelReduction.jl is (write here)"
},

{
    "location": "index.html#Installing-and-testing-package-1",
    "page": "Introduction",
    "title": "Installing and testing package",
    "category": "section",
    "text": "( write here how to install and test package )"
},

{
    "location": "index.html#Usage-example-1",
    "page": "Introduction",
    "title": "Usage example",
    "category": "section",
    "text": "( demonstrate here how to use package, simple example )"
},

{
    "location": "index.html#Contributing-1",
    "page": "Introduction",
    "title": "Contributing",
    "category": "section",
    "text": "( write here how to contribute to package )"
},

{
    "location": "theory.html#",
    "page": "Theory",
    "title": "Theory",
    "category": "page",
    "text": ""
},

{
    "location": "theory.html#Theory-1",
    "page": "Theory",
    "title": "Theory",
    "category": "section",
    "text": "( write here what is the theory behing package )"
},

{
    "location": "theory.html#References-1",
    "page": "Theory",
    "title": "References",
    "category": "section",
    "text": "( add here the list of references, if any, I leave the list below as example )Wikipedia contributors. \"Mortar methods.\" Wikipedia, The Free Encyclopedia. Wikipedia, The Free Encyclopedia.\nMaday, Yvon, Cathy Mavriplis, and Anthony Patera. \"Nonconforming mortar element methods: Application to spectral discretizations.\" (1988).\nYang, Bin, Tod A. Laursen, and Xiaonong Meng. \"Two dimensional mortar contact methods for large deformation frictional sliding.\" International journal for numerical methods in engineering 62.9 (2005): 1183-1225.\nYang, Bin, and Tod A. Laursen. \"A contact searching algorithm including bounding volume trees applied to finite sliding mortar formulations.\" Computational Mechanics 41.2 (2008): 189-205.\nWohlmuth, Barbara I. \"A mortar finite element method using dual spaces for the Lagrange multiplier.\" SIAM journal on numerical analysis 38.3 (2000): 989-1012."
},

{
    "location": "api.html#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api.html#API-documentation-1",
    "page": "API",
    "title": "API documentation",
    "category": "section",
    "text": ""
},

{
    "location": "api.html#ModelReduction.craig_bampton-NTuple{5,Any}",
    "page": "API",
    "title": "ModelReduction.craig_bampton",
    "category": "Method",
    "text": "craig_bampton(K, M, r, l, n)\n\nReduce the stiffness and mass matrices with the Craig-Bampton method. K = original stiffness matrix, M = original mass matrix, r = retained DOF:s, l = internal DOF:s, n = the number of modes to keep.\n\n\n\n"
},

{
    "location": "api.html#ModelReduction.guyan_reduction-Tuple{Any,Any,Any}",
    "page": "API",
    "title": "ModelReduction.guyan_reduction",
    "category": "Method",
    "text": "guyan_reduction(K, m, s)\n\nReduce the stiffness matrix by Guyan Reduction. K = original stiffness matrix, m = master nodes, s= slave nodes.\n\n\n\n"
},

{
    "location": "api.html#Index-1",
    "page": "API",
    "title": "Index",
    "category": "section",
    "text": "DocTestSetup = quote\n    using ModelReduction\nendModules = [ModelReduction]"
},

]}
