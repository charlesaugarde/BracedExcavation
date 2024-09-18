"""
    loading_excav(mesh, gps, constm::String, a, Gam07, theta, c, direct, phii, cohi, dens, grav, remesh, NRtol)

Loading steps for the braced excavation. 
The main difference with the [`BracedExcavation.loading_triaxial`](@ref) is that the *fint*, *fb*, and *react* induced by the further excavated soils need to be subtracted for the subsequent calculations. 

# Arguments:
- `mesh`: MESH struct datatype from [`BracedExcavation.MESH`](@ref).
- `gps`: GP struct datatype from [`BracedExcavation.GP`](@ref).
- `constm`: constitutive model.
- `a`: parameter for the small strain overlay model.
- `Gam07`: strain at 70% initial shear stiffness.
- `theta`: friction angle for the Drucker-Prager yield surface.
- `c`: cohesion for the Drucker-Prager yield surface.
- `direct`: loading direction, affecting the order the eigen values. 
- `phii`: friction angle for the interface element.
- `cohi`: cohesion for the interface element. 
- `dens`: density of materials.
- `grav`: gravitational forces.
- `remesh`: *elements* and *boundary nodes* to be excavated at the next stage. 
- `NRtol`: tolerance for the Newton_Ralphson iteration. 
"""
function loading_excav(mesh, gps, constm::String, a, Gam07, theta, c, direct, phii, cohi, dens, grav, remesh, NRtol)
    etple = remesh[1]
    benode = remesh[2]
    # NRtol = 1e-5
    NRitmax = 40
    # data need to be copied in the calculation
    uvw = copy(mesh.uvw)
    uvwold = copy(mesh.uvw)
    fext, react, fint, fb = copy(mesh.fext), copy(mesh.react), copy(mesh.fint), copy(mesh.fb)
    kval, krow, kcol = copy(mesh.kval), copy(mesh.krow), copy(mesh.kcol)

    coord = copy(mesh.coord)
    nels = length(mesh.etpl)
    nodes, nD = size(coord)
    kindex = zeros(Int64, nels + 1)
    nDoF = nodes * nD
    finte = zeros(Float64, nDoF)
    fbe = zeros(Float64, nDoF)

    NRit = 0
    oobf = fext + react - fint + fb
    oobfnorm = 2 * NRtol
    while ((NRit < NRitmax) && (oobfnorm > NRtol))
        NRit += 1
        # d, R = solve(mesh.kval, mesh.krow, mesh.kcol, oobf, mesh.bc, NRit, nDoF, mesh.lstp)
        d, R = solve(kval, krow, kcol, oobf, mesh.bc, nDoF, NRit)
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
                    Ks, Kn = gps[nel][gp].Ks, gps[nel][gp].Kn
                    # allow separation
                    if occursin(r"matz[x,y]p", mesh.mat[nel]) && deps[3] > 0.0
                        Ks = Kn = 0
                    end
                    if occursin(r"matz[x,y]n", mesh.mat[nel]) && deps[3] < 0.0
                        Ks = Kn = 0
                    end
                    D = [Ks 0 0; 0 Ks 0; 0 0 Kn]
                    dsig = D * deps
                    # allow slip (to close at first for accelerating the computation)
                    sigmax = cohi + abs(tan(phii) * (gps[nel][gp].sigP[3] + dsig[3]))  # maximum interface shear stress
                    if abs(gps[nel][gp].sigP[1] + dsig[1]) > sigmax
                        dsig[1] = sign(gps[nel][gp].sigP[1] + dsig[1]) * (sigmax - gps[nel][gp].sigP[1])
                        Ks = 0
                        D = [Ks 0 0; 0 Ks 0; 0 0 Kn]
                    end
                    if abs(gps[nel][gp].sigP[2] + dsig[2]) > sigmax
                        dsig[2] = sign(gps[nel][gp].sigP[2] + dsig[2]) * (sigmax - gps[nel][gp].sigP[2])
                        Ks = 0
                        D = [Ks 0 0; 0 Ks 0; 0 0 Kn]
                    end
                    gps[nel][gp].sig[1:3] = vec(gps[nel][gp].sigG[1:3] + gps[nel][gp].sigP[1:3] + dsig)
                    gps[nel][gp].eps[1:3] = vec(gps[nel][gp].epsP[1:3] + deps)
                    ke .+= (B' * D * B * detJ * wp[gp])
                    felem .+= (B' * gps[nel][gp].sig[1:3] * detJ * wp[gp])
                end
            elseif nen == 2
                B = zeros(Float64, 6, nen * nD)
                len = norm(mesh.coord[mesh.etpl[nel][1], :] - mesh.coord[mesh.etpl[nel][2], :])  # length
                dx = mesh.coord[mesh.etpl[nel][2], 1] - mesh.coord[mesh.etpl[nel][1], 1]
                dy = mesh.coord[mesh.etpl[nel][2], 2] - mesh.coord[mesh.etpl[nel][1], 2]
                B[1, 1] = -dx / len^2
                B[2, 2] = -dy / len^2
                B[4, 4] = dx / len^2
                B[5, 5] = dy / len^2
                detJ = len / 2.0
                deps = B * duvw[ed]'
                gp = 1  # one Gauss point
                D = elasticD(gps[nel][gp].E, gps[nel][gp].v)
                dsig = D * deps
                sig = gps[nel][gp].sigG + gps[nel][gp].sigP + dsig  # update stress
                gps[nel][gp].sig = vec(sig)
                gps[nel][gp].eps = vec(gps[nel][gp].epsP + deps)
                ke .+= B' * D * B * detJ * 2
                felem .+= B' * sig * detJ * 2
            else
                ngp, N, dNr, wp = shapefunc(nD, nen)
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
                    if mesh.mat[nel] == "LondonClay"  # soils are modelled as nonlinear materials
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
                    else   # structures are modelled as linear elastic behaviours
                        D, dsig = Hooke3d(deps, gps[nel][gp].E, gps[nel][gp].v)
                    end
                    gps[nel][gp].sig = vec(gps[nel][gp].sigG + gps[nel][gp].sigP + dsig)
                    gps[nel][gp].eps = vec(gps[nel][gp].epsP + deps)
                    ke .+= (B' * D * B * detJ * wp[gp])
                    felem .+= (B' * gps[nel][gp].sig * detJ * wp[gp])
                    if mesh.etpl[nel] in etple  # obtain the equivalent nodal force after the first excavation
                        for node = 1:nen
                            if sum(mesh.etpl[nel][node] .∈ benode) > 0
                                fbe[mesh.etpl[nel][node]*nD] += N[node, gp] * dens * grav * wp[gp] * detJ
                            end
                        end
                    end
                end
            end
            fint[ed] .+= felem'
            kval[kindex[nel]+1:kindex[nel+1]] = reshape(ke, neDoF, 1)
            if mesh.etpl[nel] in etple
                for node = 1:nen
                    if sum(mesh.etpl[nel][node] .∈ benode) > 0
                        if nen == 12
                            finte[mesh.etpl[nel][node]*nD-2:mesh.etpl[nel][node]*nD-1] .-= felem[node*nD-2:node*nD-1]
                            finte[mesh.etpl[nel][node]*nD] += felem[node*nD]
                            # finte[mesh.etpl[nel][node]*nD-2:mesh.etpl[nel][node]*nD] .-= felem[node*nD-2:node*nD]
                        else
                            finte[mesh.etpl[nel][node]*nD-2:mesh.etpl[nel][node]*nD] .+= felem[node*nD-2:node*nD]
                        end
                    end
                end
            end
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
    mesh.krow = krow
    mesh.kcol = kcol
    mesh.uvw = uvw
    mesh.fint = fint
    mesh.react = react
    mesh.fb = fb
    for nel in 1:nels
        for gp in 1:mesh.ngp[nel]
            gps[nel][gp].sigP = copy(gps[nel][gp].sig - gps[nel][gp].sigG)
            gps[nel][gp].epsP = copy(gps[nel][gp].eps)
        end
    end

    return [finte, fbe]
end
