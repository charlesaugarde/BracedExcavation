using Plots

"""
    lineplot_triaxial(lstps, constm, results, bc1, lstpc, lstpx, outdir)

Plot the results of triaxial tests: load step and displacement in x, y, z-directions, stress strain relationship, and effective stress paths.

# Arguments:
- `lstps`: total load steps.
- `constm`: constitutive models.
- `results`: results data.
- `bc1`: boundary condition for the compression test.
- `lstpc`: load steps for the isotropic consolidation.
- `lstpx`: load steps for the drained compression.
- `theta`: friction angle for the Drucker-Prager yield surface.
- `c`: cohesion for the Drucker-Prager yield surface.
- `direct`: loading direction: 1 for monotonic loading; -1 for unloading/reloading.
- `outdir`: directory for saving the figures.
"""
function lineplot_triaxial(lstps, constm, results, bc1, lstpc, lstpx, theta, c, direct, outdir)
    tstep = collect(0:lstps)
    disx = results[1, :]
    disy = results[2, :]
    disz = results[3, :]
    sigp = results[4, :]
    sigq = results[5, :]

    # displacement
    fig = plot!(tstep, disx, label="x-displacement", markers=:circle, markercolor=:white, ms=5, linewidth=2, dpi=200)
    fig = plot!(tstep, disy, label="y-displacement", markers=:circle, markercolor=:white, ms=5, linewidth=2, dpi=200)
    fig = plot!(tstep, disz, label="z-displacement", markers=:circle, markercolor=:white, ms=5, linewidth=2, dpi=200)
    # xlims!(0, lstps)
    # ylims!(0, 100)
    xlabel!("Load step (-)")
    ylabel!("Displacement (-)")
    savefig(fig, "./" * outdir * "/" * constm * "_displacement.png")

    # effective stress path
    fig = plot(sigp, sigq, label="Effective stress path", markers=:circle, color=:black, ms=5, linewidth=2, markercolor=:white, dpi=200)
    if occursin("DP", constm)
        # yield surface
        a = collect(0:500)
        eta = 6 * sin(theta) / (sqrt(3) * (3 - sin(theta)))
        xi = 6 * cos(theta) / (sqrt(3) * (3 - sin(theta)))
        b = -(xi * c .- eta * a) * sqrt(3)
        plot!(a, b, label="Drucker-Prager yield surface")
    end
    # xlims!(100, 340)
    # ylims!(0, 100)
    xlabel!("Hydrostatic stress (kPa)")
    ylabel!("Deviatoric stress (kPa)")
    savefig(fig, "./" * outdir * "/" * constm * "_effective_stress_path.png")

    # stress-strain relationship
    if direct == 1
        da = -minimum(bc1[:, 2]) * collect(0:lstpx)  # for compression
    elseif direct == -1
        da = -maximum(bc1[:, 2]) * collect(0:lstpx)   # for extension
    end
    fig = plot(abs.(da), sigq[lstpc+1:1:end], label="Stress strain relationship", markers=:circle, color=:black, ms=5, linewidth=2, markercolor=:white, dpi=200)
    # xlims!(0, 0.1)
    # ylims!(0, 100)
    xlabel!("Axial strain (-)")
    ylabel!("Deviatoric stress (kPa)")
    savefig(fig, "./" * outdir * "/" * constm * "_stress_strain.png")

    return nothing
end

"""
    lineplot_excavation(lstp, uvw, coord, maxnode, outdir)

Plot the wall deflection due to excavation at different excavation stages.

# Arguments:
- `lstp`: load step.
- `uvw`: displacement.
- `coord`: coordinates of the model nodes. 
- `maxnode`: maximum tag of original nodes.
- `outdir`: directory for saving the figures. 
"""
function lineplot_excavation(lstp, uvw, coord, maxnode::Int64, outdir)
    # 1. wall deflection
    tol = 1e-5
    z = []  # the depth of the nodes 
    dispxc = []  # wall x-deflection in the corner
    dispyc = []  # wall y-deflection in the corner
    dispxx = []  # wall x-deflection in the x-direction
    dispyx = []  # wall y-deflection in the x-direction
    dispxy = []  # wall x-deflection in the y-direction
    dispyy = []  # wall y-deflection in the y-direction
    nodes, nD = size(coord)
    if maxnode > 0
        snode = nodes - maxnode
    else
        snode = 1
    end
    for node = snode:nodes
        # for the wall in the corner 
        if (abs(coord[node, 1] - 0.48) <= tol && abs(coord[node, 2] - 0.48) <= tol)
            push!(z, 1.0 - coord[node, 3])
            push!(dispxc, uvw[node*nD-2])
            push!(dispyc, uvw[node*nD-1])
        end
        if (abs(coord[node, 1] - 0.48) <= tol && abs(coord[node, 2] - 0.24) <= tol)
            push!(dispxx, uvw[node*nD-2])
            push!(dispyx, uvw[node*nD-1])
        end
        if (abs(coord[node, 1] - 0.24) <= tol && abs(coord[node, 2] - 0.48) <= tol)
            push!(dispxy, uvw[node*nD-2])
            push!(dispyy, uvw[node*nD-1])
        end
    end
    # sort by the depth 
    perm = sortperm(z)

    fig = plot!(dispxc[perm], z[perm], markers=:circle, ms=5, linewidth=2, markercolor=:white, label="Wall x-deflection @ center", dpi=300)
    fig = plot!(dispyc[perm], z[perm], markers=:circle, ms=5, linewidth=2, markercolor=:white, label="Wall y-deflection @ center", dpi=300)
    fig = plot!(dispxx[perm], z[perm], markers=:circle, ms=5, linewidth=2, markercolor=:white, label="Wall x-deflection @ x-axis", dpi=300)
    # fig = plot!(dispyx[perm], z[perm], markers=:circle, ms=5, linewidth=2, markercolor=:white, label="Wall y-deflection @ x-axis", dpi=300)
    # fig = plot!(dispxy[perm], z[perm], markers=:circle, ms=5, linewidth=2, markercolor=:white, label="Wall x-deflection @ y-axis", dpi=300)
    fig = plot!(dispyy[perm], z[perm], markers=:circle, ms=5, linewidth=2, markercolor=:white, label="Wall y-deflection @ y-axis", dpi=300)
    # ylims!(0, 50)
    # xlims!(-0.002, 0.002)
    yflip!(true)
    xlabel!("Wall deflection (m)")
    ylabel!("Depth of the diaphragm wall (m)")
    savefig(fig, outdir * "/walldefelction_$lstp.png")
    # savefig(fig, outdir * "/walldefelction_$lstp.svg")

    return nothing
end
