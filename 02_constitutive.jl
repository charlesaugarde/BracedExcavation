using AppleAccelerate, Ipopt
push!(LOAD_PATH, "./src/")
using BracedExcavation

project_name = "ElementCode"
outdir = project_name * "_out"
if isdir(outdir)
    rm(outdir, recursive=true)
end
mkpath(outdir)

# unithex(project_name, true)

# PARAMS = PARAMSDT(CONSTM="LSSO")  # Linear small strain overlay SSO model
# PARAMS = PARAMSDT(CONSTM="DPconst")  # Elasto-perfect plastic Drucker-Prager (DP) model
PARAMS = PARAMSDT(CONSTM="DPSSO", THETA=π / 14)  # combination of SSO and DP model (DPSSO)

mesh = MESH()

mesh.coord, etpl, benode = readgmsh_triaxial(project_name * ".msh")
mesh.mat = getindex.(etpl, 1)
mesh.etpl = getindex.(etpl, 2)
mesh.ngp = getindex.(etpl, 3)

gps = gpmaterial(mesh, PARAMS)
results = zeros(Float64, 5, PARAMS.LSTPC + PARAMS.LSTPX + 1)
apex = 7

geostatic_triaxial(mesh, gps, PARAMS.DENS, PARAMS.GRAV, PARAMS.K0, PARAMS.H)
postpro(mesh.uvw, mesh.lstp, mesh.coord, mesh.etpl, gps, outdir)

bc0, fext0 = boundary_consolidation(mesh.coord, benode[1], benode[2], benode[3], PARAMS.P)
bc1, fext1 = boundary_compression(mesh.coord, benode[2], benode[3], PARAMS.P, PARAMS.U, PARAMS.LSTPX)

for lstp = 1:PARAMS.LSTPC
    mesh.lstp = lstp
    mesh.bc = bc0
    mesh.fext = lstp / PARAMS.LSTPC * fext0
    loading_triaxial(mesh, gps, PARAMS.CONSTM, PARAMS.A, PARAMS.GAM07, PARAMS.THETA, PARAMS.C, 1, PARAMS.TOL)
    postpro(mesh.uvw, mesh.lstp, mesh.coord, mesh.etpl, gps, outdir)
    results[1:3, lstp+1] = mesh.uvw[apex*3-2:apex*3]
    results[4, lstp+1] = -sum(gps[1][1].sigP[1:3]) / 3.0 / 1e3
    results[5, lstp+1] = (gps[1][1].sigP[1] - gps[1][1].sigP[3]) / 1e3
end

for lstp = PARAMS.LSTPC+1:PARAMS.LSTPX+PARAMS.LSTPC
    mesh.lstp = lstp
    mesh.bc = bc1
    mesh.fext = fext1
    loading_triaxial(mesh, gps, PARAMS.CONSTM, PARAMS.A, PARAMS.GAM07, PARAMS.THETA, PARAMS.C, 1, PARAMS.TOL)
    postpro(mesh.uvw, mesh.lstp, mesh.coord, mesh.etpl, gps, outdir)
    results[1:3, lstp+1] = mesh.uvw[apex*3-2:apex*3]
    results[4, lstp+1] = -sum(gps[1][1].sigP[1:3]) / 3.0 / 1e3
    results[5, lstp+1] = (gps[1][1].sigP[1] - gps[1][1].sigP[3]) / 1e3
end

lineplot_triaxial(PARAMS.LSTPC + PARAMS.LSTPX, PARAMS.CONSTM, results, bc1, PARAMS.LSTPC, PARAMS.LSTPX, PARAMS.THETA, PARAMS.C, 1, outdir)
