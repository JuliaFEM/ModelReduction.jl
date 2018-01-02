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
    "text": "ModelReduction.jl is a Julia package to perform model reduction methods for i.e. multibody dynamics problems. The packcage includes model order reduction methods such as the Guyan reduction and the Craig-Bampton method.Reducing the sizes of stiffness and mass matrices of the model will greatly decrease the computation resources needed when performing dynamic analyses.ModelReduction.jl is a part of JuliaFEM. All codes are MIT licensed."
},

{
    "location": "index.html#Installing-and-testing-the-package-1",
    "page": "Introduction",
    "title": "Installing and testing the package",
    "category": "section",
    "text": "Download and unzip the packge files into your Julia folder and then install the package the same way other Julia packages are installed.julia> Pkg.add(\"ModelReduction\")\nTest the package with Pkg.test etc.julia> Pkg.test(\"ModelReduction\")\n"
},

{
    "location": "index.html#Usage-example-1",
    "page": "Introduction",
    "title": "Usage example",
    "category": "section",
    "text": "This example demonstrates how to use the Craig-Bampton method function.Problem setup:K =\n[2 -1 0 0;\n-1 2 -1 0;\n 0 -1 2 -1;\n 0 0 -1 1]\n\nM =\n[2 0 0 0;\n 0 2 0 0;\n 0 0 2 0;\n 0 0 0 1]\n\nr = [4]\nl = [1, 2, 3]\nn = 1\nK = original stiffness matrix, M = original mass matrix, r = retained DOF:s, l = internal DOF:s, n = the number of the internal modes to keep. Calculate the reduced mass and stiffness matrices Mred and Kred.craig_bampton(K, M, r, l, n)\n\n# output\n\nMred, Kred = ([2.75 -1.20711; -1.20711 1.0], [0.25 0.0; 0.0 0.292893])\n"
},

{
    "location": "index.html#Contributing-1",
    "page": "Introduction",
    "title": "Contributing",
    "category": "section",
    "text": "Have a new great idea and want to share it with the open source community?  From here and here you can look for coding styles. Here it is explained how to contribute to  open source project, in general."
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
    "text": "The Craig-Bampton method is a dynamic reduction technique that reduces the mass and stiffness matrices of the model by expressing the boundary modes in physical coordinates and the elastic modes in modal coordinates.The equation of motion is:beginequation\nboldsymbolMddotboldsymbolu+boldsymbolKboldsymbolu=boldsymbolf\nendequationThe matrices are partitioned into boundary nodes R and the independent elastic nodes L:beginequation\nboldsymbolu=beginbmatrixboldsymbolu_mathrmR\nboldsymbolu_mathrmL\nendbmatrix\nendequationEquation (1) becomes:beginequation\nbeginbmatrixboldsymbolM_mathrmRR  boldsymbolM_mathrmmathrmRL\nboldsymbolM_mathrmRR  boldsymbolM_mathrmLL\nendbmatrixbeginbmatrixddotboldsymbolu_mathrmR\nddotboldsymbolu_mathrmL\nendbmatrix+beginbmatrixboldsymbolK_mathrmRR  boldsymbolK_mathrmRL\nboldsymbolK_mathrmLR  boldsymbolK_mathrmLL\nendbmatrixbeginbmatrixboldsymbolu_mathrmR\nboldsymbolu_mathrmL\nendbmatrix=boldsymbolf\nendequationThe degrees of freedom are are transformed to hybrid coordinatesbeginequation\nbeginbmatrixboldsymbolu_mathrmR\nboldsymbolu_mathrmL\nendbmatrix=beginbmatrixboldsymbolI  boldsymbol0\nboldsymbolX_mathrmR  boldsymbolX_mathrmL\nendbmatrixbeginbmatrixboldsymbolu_mathrmR\nboldsymbolq_mathrmm\nendbmatrix\nendequationEquation (1) can be rewritten asbeginequation\nbeginbmatrixboldsymbolM_mathrmRR  boldsymbolM_mathrmmathrmRL\nboldsymbolM_mathrmLR  boldsymbolM_mathrmLL\nendbmatrixbeginbmatrixboldsymbolI  boldsymbol0\nboldsymbolX_mathrmR  boldsymbolX_mathrmL\nendbmatrixbeginbmatrixddotboldsymbolu_mathrmR\nddotboldsymbolq_mathrmm\nendbmatrix+beginbmatrixboldsymbolK_mathrmRR  boldsymbolK_mathrmRL\nboldsymbolK_mathrmLR  boldsymbolK_mathrmLL\nendbmatrixbeginbmatrixboldsymbolI  boldsymbol0\nboldsymbolX_mathrmR  boldsymbolX_mathrmL\nendbmatrixbeginbmatrixboldsymbolu_mathrmR\nboldsymbolq_mathrmm\nendbmatrix=beginbmatrixboldsymbolf_mathrmR\nboldsymbol0\nendbmatrix\nendequationEquation (1) reduces tobeginequation\nboldsymbolK_mathrmLRboldsymbolK_mathrmLRboldsymbolu_mathrmR+boldsymbolK_mathrmLLboldsymbolu_mathrmL\nendequationThe internal degrees of freedom can be expressed asbeginequation\nboldsymbolu_mathrmL=-boldsymbolK_mathrmLL^-1boldsymbolK_mathrmLRboldsymbolu_mathrmR=boldsymbolX_mathrmRboldsymbolu_mathrmR\nendequationwherebeginequation\nboldsymbolX_mathrmR=-boldsymbolK_mathrmLL^-1boldsymbolK_mathrmLR\nendequationTo determine mathitboldsymbolX_mathrmL the retained degrees of freedom are fixed. The equation of motion reduces tobeginequation\nboldsymbolM_mathrmLLddotboldsymbolu_mathrmL+boldsymbolK_mathrmLLboldsymbolu_mathrmL=0\nendequationBy assuming harmonic response and substituting the coordinate transformation (4)beginequation\n(-omega^2boldsymbolM_mathrmLL+boldsymbolK_mathrmLL)boldsymbolX_mathrmLboldsymbolq_mathrmme^iomega t=0\nendequationThe eigenvectors can be normalized:beginequation\nboldsymbolX_mathrmL^mathrmTboldsymbolM_mathrmLLboldsymbolX_mathrmL=boldsymbolI\nendequationbeginequation\nboldsymbolX_mathrmL^mathrmTboldsymbolK_mathrmLLboldsymbolX_mathrmL=boldsymbolLambda\nendequationSince boldsymbolX_mathrmR in (9) contains boldsymbolK_mathrmLL^-1, an inverse of boldsymbolK_mathrmLL, determining it will require lots of computing resources. This can be avoided by determining the boldsymbolK_mathrmLL inverse as follows.beginequation\n-boldsymbolK_mathrmLL^-1=boldsymbolX_mathrmLboldsymbolLambda^-1boldsymbolX_mathrmL^mathrmT\nendequationIn order to get the dynamic equations of the system, equation (6) is multiplied with the coordination transformation matrix.beginmultline\nbeginbmatrixboldsymbolI  boldsymbolX_mathrmR^mathrmT\nboldsymbol0  boldsymbolX_mathrmL^mathrmT\nendbmatrixbeginbmatrixboldsymbolM_mathrmRR  boldsymbolM_mathrmmathrmRL\nboldsymbolM_mathrmLR  boldsymbolM_mathrmLL\nendbmatrixbeginbmatrixboldsymbolI  boldsymbol0\nboldsymbolX_mathrmR  boldsymbolX_mathrmL\nendbmatrixbeginbmatrixddotboldsymbolu_mathrmR\nddotboldsymbolq_mathrmm\nendbmatrix\n+beginbmatrixboldsymbolI  boldsymbolX_mathrmR^mathrmT\nboldsymbol0  boldsymbolX_mathrmL^mathrmT\nendbmatrixbeginbmatrixboldsymbolK_mathrmRR  boldsymbolK_mathrmRL\nboldsymbolK_mathrmLR  boldsymbolK_mathrmLL\nendbmatrixbeginbmatrixboldsymbolI  boldsymbol0\nboldsymbolX_mathrmR  boldsymbolX_mathrmL\nendbmatrixbeginbmatrixmathbfmathitboldsymbolu_mathrmR\nboldsymbolq_mathrmm\nendbmatrix=beginbmatrixboldsymbolI  boldsymbolX_mathrmR^mathrmT\nboldsymbol0  boldsymbolX_mathrmL^mathrmT\nendbmatrixbeginbmatrixboldsymbolf_mathrmR\nboldsymbol0\nendbmatrix\nendmultlineBy simplifying the equation of motion (1) becomesbeginmultline\nbeginbmatrixboldsymbolM_mathrmRR+boldsymbolX_mathrmR^mathrmTboldsymbolM_mathrmLR+boldsymbolX_mathrmR^mathrmTboldsymbolM_mathrmLLboldsymbolX_mathrmR  boldsymbolM_mathrmmathrmRLboldsymbolX_mathrmL+boldsymbolX_mathrmR^mathrmTboldsymbolM_mathrmLLboldsymbolX_mathrmL\nboldsymbolX_mathrmL^mathrmTboldsymbolM_mathrmLR+boldsymbolX_mathrmL^mathrmTboldsymbolM_mathrmLLboldsymbolX_mathrmR  boldsymbolI\nendbmatrixbeginbmatrixddotboldsymbolu_mathrmR\nddotboldsymbolq_mathrmm\nendbmatrix\n+beginbmatrixboldsymbolK_mathrmRR+boldsymbolK_mathrmRLboldsymbolX_mathrmR  boldsymbol0\nboldsymbol0  boldsymbolLambda\nendbmatrixbeginbmatrixboldsymbolu_mathrmR\nboldsymbolq_mathrmm\nendbmatrix=beginbmatrixboldsymbolf_mathrmR\nboldsymbol0\nendbmatrix\nendmultline"
},

{
    "location": "theory.html#References-1",
    "page": "Theory",
    "title": "References",
    "category": "section",
    "text": "Qu, Zu-Qing. Model Order Reduction Techniques (2004). p. 322 - 329.\nYoung, John T. Prime on the Craig-Bampton Method (2000). p. 5 - 17."
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
    "location": "api.html#ModelReduction.global_mass-Tuple{Any,Any}",
    "page": "API",
    "title": "ModelReduction.global_mass",
    "category": "Method",
    "text": "global_mass(m, N)\n\nCalculate global mass matrix of size N.\n\n\n\n"
},

{
    "location": "api.html#ModelReduction.global_stiffness-Tuple{Any,Any}",
    "page": "API",
    "title": "ModelReduction.global_stiffness",
    "category": "Method",
    "text": "global_stiffness(k, N)\n\nCalculate global stiffness matrix of size N.\n\n\n\n"
},

{
    "location": "api.html#ModelReduction.guyan_reduction-Tuple{Any,Any,Any}",
    "page": "API",
    "title": "ModelReduction.guyan_reduction",
    "category": "Method",
    "text": "guyan_reduction(K, m, s)\n\nReduce the stiffness matrix by Guyan Reduction. K = original stiffness matrix, m = master nodes, s= slave nodes.\n\n\n\n"
},

{
    "location": "api.html#ModelReduction.sort_nodes-Tuple{Any,Any}",
    "page": "API",
    "title": "ModelReduction.sort_nodes",
    "category": "Method",
    "text": "sort_nodes(nodes, node_sets)\n\nCreate the r and l arrays for the retained and internal nodes. nodes = all nodes of the model node_sets = node sets containing the nodes that are to be retained\n\n\n\n"
},

{
    "location": "api.html#Index-1",
    "page": "API",
    "title": "Index",
    "category": "section",
    "text": "DocTestSetup = quote\n    using ModelReduction\nendModules = [ModelReduction]"
},

]}
