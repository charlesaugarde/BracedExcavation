using AppleAccelerate, Ipopt
push!(LOAD_PATH, "./src/")
using BracedExcavation

project_name = "BracedExcavation"
outdir = project_name * "_out"
if isdir(outdir)
    rm(outdir, recursive=true)
end
mkpath(outdir)

# excavation(project_name, true)

# For Drucker-Prager, set THETA=π/6, TOL=1e-4, (geomaker.jl N=Nu=4) for demostration of plastic behaviours of first excavation
# For Drucker-Prager, set THETA = π/3, TOL=1e-9, (geomaker.jl, N=Nu==4) for demonstration of elastic behaviours of first excavation
PARAMS = PARAMSDT(CONSTM="DPSSO", TOL=1e-4, THETA=π/3, E0=1e8, E1=1e10)  
# PARAMS = PARAMSDT(CONSTM="DPconst", TOL=1e-4, THETA=π/6)  
# PARAMS = PARAMSDT(CONSTM="Hooke3d")  
# PARAMS = PARAMSDT(CONSTM="LSSO")  

mesh = MESH()

mesh0, mesh1, mesh2 = readgmsh_excavnointer(project_name * ".msh");
Snnodes = 0  # no additional nodes added to interface elements

mesh.coord = mesh0[1]
mesh.mat = getindex.(mesh0[2], 1)
mesh.etpl = getindex.(mesh0[2], 2)
mesh.ngp = getindex.(mesh0[2], 3)

gps = gpmaterial(mesh, PARAMS)
println("    Initialise the materials.")

fexc1, fexc2, kindex = geostatic_excavation(mesh, gps, mesh1, mesh2, PARAMS.DENS, PARAMS.GRAV, PARAMS.K0, PARAMS.H)
postpro(mesh.uvw, mesh.lstp, mesh.coord, mesh.etpl, gps, outdir)
# lineplot_excavation(mesh.lstp, mesh.uvw, mesh.coord, Snnodes, outdir)
println("    Geostatic state.")

gps, kindex, remesh = firstexcavation(mesh, kindex, fexc1, mesh0, mesh1, gps, mesh2)
println("    Conduct the first excavation.")

mesh.bc = boundary_excavation(mesh.coord, PARAMS.W, PARAMS.L)

fexv = loading_excav(mesh, gps, PARAMS.CONSTM, PARAMS.A, PARAMS.GAM07, PARAMS.THETA, PARAMS.C, -1, PARAMS.THETAI, PARAMS.CI, PARAMS.DENS, PARAMS.GRAV, remesh, PARAMS.TOL)
postpro(100 * mesh.uvw, mesh.lstp, mesh.coord, mesh.etpl, gps, outdir)
# lineplot_excavation(mesh.lstp, mesh.uvw, mesh.coord, Snnodes, outdir)
println("    Deformation due to first excavation.")

gps, kindex = secondexcavation(mesh, kindex, fexv, mesh0, mesh1, mesh2, gps, fexc2)
println("    Conduct the second excavation.")

mesh.bc = boundary_excavation(mesh.coord, PARAMS.W, PARAMS.L)

fexv = loading_excav(mesh, gps, PARAMS.CONSTM, PARAMS.A, PARAMS.GAM07, PARAMS.THETA, PARAMS.C, -1, PARAMS.THETAI, PARAMS.CI, PARAMS.DENS, PARAMS.GRAV, remesh, PARAMS.TOL)
postpro(1 * mesh.uvw, mesh.lstp, mesh.coord, mesh.etpl, gps, outdir)
lineplot_excavation(mesh.lstp, mesh.uvw, mesh.coord, Snnodes, outdir)
println("    Deformation due to second excavation.")
