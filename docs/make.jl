push!(LOAD_PATH, "../src/")
using Documenter, BracedExcavation

format = Documenter.HTML(edit_link="master", prettyurls=get(ENV, "CI", nothing) == "true")

About = "Introduction" => "index.md"

Tutorials = "Tutorials" => ["tutorials/01_triaxial.md", "tutorials/02_constitutive.md", "tutorials/03_excavation.md"]

Notes = "Developer Notes" => ["notes/01_elemgeo.md", "notes/02_equivforce.md", "notes/03_excav_nointer.md", "notes/04_interface.md"]

API = "API Reference" => "api.md"

PAGES = [About, Tutorials, Notes, API]

makedocs(sitename="Braced Excavation",
    authors="Mao Ouyang",
    format=format,
    checkdocs=:exports,
    pages=PAGES,
)

deploydocs(
    repo="github.com/einoo/BracedExcavation.jl.git",
)
