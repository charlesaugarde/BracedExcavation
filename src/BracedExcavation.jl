module BracedExcavation

include("common.jl")
include("shapefunc.jl")
include("constitutive.jl")
export interpretElement, engtocauchy, cauchytoeng, PARAMSDT, MESH, GP

# make the geometry
include("geomaker.jl")
export unithex, unittet, unitwed, unittetwed, excavation

# input parameters
include("parameter.jl")

# read the mesh data
include("readgmsh.jl")
export readgmsh_triaxial, readgmsh_excavation, readgmsh_excavnointer

# initialisation
include("gpmaterial.jl")
include("geostatic.jl")
export geostatic_triaxial, geostatic_excavation, geostatic_test, gpmaterial

# boundary conditions
include("boundary.jl")
export boundary_consolidation, boundary_compression, boundary_excavation

# loading steps
include("loading_triaxial.jl")
include("solve.jl")
include("loading_excavinter.jl")
export loading_triaxial, solve, loading_excav

# postprocess
include("plotting.jl")
include("lineplot.jl")
export postpro, lineplot_triaxial, lineplot_excavation

# excavation
include("firstexcav.jl")
include("secondexcav.jl")
export firstexcavation, secondexcavation

end
