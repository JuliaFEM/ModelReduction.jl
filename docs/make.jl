# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

using Documenter, ModelReduction

makedocs(modules=[ModelReduction],
         format = :html,
         sitename = "ModelReduction",
         pages = [
                  "Introduction" => "index.md",
                  "Theory" => "theory.md",
                  "API" => "api.md"
                 ])
