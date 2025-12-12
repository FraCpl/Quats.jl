using Quats
using Documenter

DocMeta.setdocmeta!(Quats, :DocTestSetup, :(using Quats); recursive=true)

makedocs(;
    modules=[Quats],
    authors="F. Capolupo",
    sitename="Quats.jl",
    format=Documenter.HTML(;
        canonical="https://FraCpl.github.io/Quats.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/FraCpl/Quats.jl",
    devbranch="master",
)
