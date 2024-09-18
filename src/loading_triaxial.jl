"""
    loading_triaxial(mesh, gps, constm, G0, a, Gam07, direct)

Loading step for the FE analysis. 
All the MESH and GP variables are updated as the global properties. 

# Arguments:
- `mesh`: MESH struct datatype from [`BracedExcavation.MESH`](@ref).
- `gps`: GP struct datatype from [`BracedExcavation.GP`](@ref).
- `constm`: constitutive model.
- `a`: parameter for the small strain overlay model.
- `Gam07`: strain at 70% initial shear stiffness.
- `theta`: friction angle for the Drucker-Prager yield surface.
- `c`: cohesion for the Drucker-Prager yield surface.
- `direct`: loading direction, affecting the order the eigen values. 
- `NRtol`: tolerance for the Newton_Ralphson iteration. 
"""
function loading_triaxial(mesh, gps, constm::String, a, Gam07, theta, c, direct, NRtol)
    # NRtol = 1e-9
    NRitmax = 40
    # data need to be copied in the calculation
    uvw = copy(mesh.uvw)
    uvwold = copy(mesh.uvw)
    fext, react, fint, fb = copy(mesh.fext), copy(mesh.react), copy(mesh.fint), copy(mesh.fb)
    kval = copy(mesh.kval)

    coord = copy(mesh.coord)
    nels = length(mesh.etpl)
    nodes, nD = size(coord)
    kindex = zeros(Int64, nels + 1)
    nDoF = nodes * nD

    NRit = 0
    oobf = fext + react - fint + fb
    oobfnorm = 2 * NRtol
    while ((NRit < NRitmax) && (oobfnorm > NRtol))
        NRit += 1
        # d, R = solve(mesh.kval, mesh.krow, mesh.kcol, oobf, mesh.bc, NRit, nDoF, mesh.lstp)
        d, R = solve(kval, mesh.krow, mesh.kcol, oobf, mesh.bc, nDoF, NRit)
        uvw .+= d
        react .+= R
        duvw = uvw - uvwold
        fint .*= 0  # zero the internal force, be careful on the gravitational force
        for nel = 1:nels
            nen = length(mesh.etpl[nel])
            neDoF = (nen * nD)^2
            kindex[nel+1] = kindex[nel] + neDoF
            ed = ones(Int, nD, 1) * reshape(mesh.etpl[nel] * nD, 1, nen) - collect(nD-1:-1:0) * ones(Int, 1, nen)
            ed = reshape(ed, 1, nen * nD)
            ke = zeros(Float64, nen * nD, nen * nD)  # zero the element stiffness matrix
            felem = zeros(Float64, nen * nD)
            if nen == 12
                ngp, N, dNr, wp = shapefunc(2, nen)
                if occursin(r"matzy[p,n]", mesh.mat[nel])
                    JT = dNr[:, 1:6] * coord[mesh.etpl[nel][1:6], [1, 3]]
                    T = [1 0 0; 0 0 1; 0 1 0]
                elseif occursin(r"matzx[p,n]", mesh.mat[nel])
                    JT = dNr[:, 1:6] * coord[mesh.etpl[nel][1:6], [2, 3]]
                    T = [0 1 0; 0 0 1; 1 0 0]
                end
                for gp in 1:ngp
                    B = zeros(Float64, 3, nen * nD)
                    indx = 2 * gp .- collect(1:-1:0)
                    detJ = abs(det(JT[indx, :]))
                    for (i, j) in Base.product(1:3, 1:6)
                        B[i, 3*j-3+i] = N[j, gp]
                    end
                    B[:, 19:end] = -copy(B[:, 1:18])
                    # coordinate transformation (use a simplified way here)
                    B = T * B
                    deps = B * duvw[ed]'  # strain increment
                    # mechanical behaviours
                    Ks, Kn = gps[nel][gp].Ks, gps[nel][gp].Kn  # no separation/slip
                    D = [Ks 0 0; 0 Ks 0; 0 0 Kn]
                    dsig = D * deps
                    gps[nel][gp].sig[1:3] = vec(gps[nel][gp].sigG[1:3] + gps[nel][gp].sigP[1:3] + dsig)
                    gps[nel][gp].eps[1:3] = vec(gps[nel][gp].epsP[1:3] + deps)
                    ke .+= (B' * D * B * detJ * wp[gp])
                    felem .+= (B' * gps[nel][gp].sig[1:3] * detJ * wp[gp])
                end
            else
                ngp, _, dNr, wp = shapefunc(nD, nen)
                JT = dNr * coord[mesh.etpl[nel], :]
                for gp = 1:mesh.ngp[nel]
                    B = zeros(Float64, 6, nen * nD)
                    indx = nD * gp .- collect(nD-1:-1:0)
                    detJ = abs(det(JT[indx, :]))
                    dNx = JT[indx, :] \ dNr[indx, :]
                    B[[1 4 6], 1:nD:end] = dNx
                    B[[4 2 5], 2:nD:end] = dNx
                    B[[6 5 3], 3:nD:end] = dNx
                    deps = B * duvw[ed]'
                    if occursin("Hooke3d", constm)
                        D, dsig = Hooke3d(deps, gps[nel][gp].E, gps[nel][gp].v)
                    elseif occursin("LSSO", constm)
                        D, dsig, gps[nel][gp].H, _ = smallstrain(deps, gps[nel][gp].H0, gps[nel][gp].G, gps[nel][gp].v, a, Gam07, direct)
                    elseif occursin("DPconst", constm)
                        D, dsig, deps, _ = DPconst(deps, gps[nel][gp].E, gps[nel][gp].v, theta, c, gps[nel][gp].sigP, gps[nel][gp].epsP)
                    elseif occursin("DPSSO", constm)
                        if gps[nel][gp].Y == false   # elastic region
                            _, _, gps[nel][gp].H, gps[nel][gp].E = smallstrain(deps, gps[nel][gp].H0, gps[nel][gp].G, gps[nel][gp].v, a, Gam07, direct)
                            D, dsig, deps, gps[nel][gp].Y = DPconst(deps, gps[nel][gp].E, gps[nel][gp].v, theta, c, gps[nel][gp].sigP, gps[nel][gp].epsP)
                        else  # plastic region, not update the shear stiffness
                            D, dsig, deps, gps[nel][gp].Y = DPconst(deps, gps[nel][gp].E, gps[nel][gp].v, theta, c, gps[nel][gp].sigP, gps[nel][gp].epsP)
                        end
                    else
                        println("Please select the constitutive model from: Hooke3d, LSSO, DPconst, DPSSO.")
                        exit()
                    end
                    # D, dsig = Hooke3d(deps, gps[nel][gp].E, gps[nel][gp].v)
                    gps[nel][gp].sig = vec(gps[nel][gp].sigG + gps[nel][gp].sigP + dsig)
                    gps[nel][gp].eps = vec(gps[nel][gp].epsP + deps)
                    ke .+= (B' * D * B * detJ * wp[gp])
                    felem .+= (B' * gps[nel][gp].sig * detJ * wp[gp])
                end
            end
            fint[ed] .+= felem'
            kval[kindex[nel]+1:kindex[nel+1]] = reshape(ke, neDoF, 1)
        end
        oobf = fext + react - fint + fb
        if norm(fext + react + fint + fb) < NRtol
            oobfnorm = norm(oobf)
        else
            oobfnorm = norm(oobf) / norm(fext + react + fint + fb)
        end
    end
    println("Time step: $(mesh.lstp); Number of iterations: $NRit; Residual error: $oobfnorm")

    # Update the physical states
    mesh.kval = kval
    mesh.uvw = uvw
    mesh.fint = fint
    mesh.react = react
    mesh.fb = fb
    for nel in 1:nels
        for gp in 1:mesh.ngp[nel]
            gps[nel][gp].sigP = copy(gps[nel][gp].sig - gps[nel][gp].sigG)
            gps[nel][gp].epsP = copy(gps[nel][gp].eps)
            # update the H (strain history), and Y (yield criterion)
            gps[nel][gp].H0 = copy(gps[nel][gp].H)
        end
    end

    return nothing
end
