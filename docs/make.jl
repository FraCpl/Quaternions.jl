using Quaternions
using Documenter

DocMeta.setdocmeta!(Quaternions, :DocTestSetup, :(using Quaternions); recursive=true)

makedocs(;
    modules=[Quaternions],
    authors="F. Capolupo",
    repo="https://github.com/FraCpl/Quaternions.jl/blob/{commit}{path}#{line}",
    sitename="Quaternions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true", canonical="https://FraCpl.github.io/Quaternions.jl", edit_link="master", assets=String[]
    ),
    pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/FraCpl/Quaternions.jl", devbranch="master")
