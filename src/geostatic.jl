"""
    geostatic_triaxial(mesh, gps, dens, grav, K0, H)

Set the geostatic conditions for the triaxial tests using the simplified K0 approach. 

# Arguments:
- `mesh`: MESH structure defined by [`BracedExcavation.MESH`](@ref).
- `gps`: GP structure defined by [`BracedExcavation.GP`](@ref).
- `dens`: density of materials.
- `grav`: gravitational forces.
- `K0`: K0 value.
- `H`: height of the model.
"""
function geostatic_triaxial(mesh, gps, dens, grav, K0, H)
    mesh.lstp = 0
    coord = mesh.coord  # coordinate
    etpl = mesh.etpl  # topography
    mat = mesh.mat  # materials

    nodes, nD = size(coord)
    nels = length(etpl)
    nDoF = nodes * nD
    fb = zeros(Float64, nDoF)
    fint = zeros(Float64, nDoF)
    kindex = zeros(Int64, nels + 1)
    neDoF = 0
    for e in etpl
        neDoF += (length(e) * nD)^2
    end
    kval = zeros(Float64, neDoF) # values for the global stiffness matrix
    krow = zeros(Int64, neDoF)  # [sparse] row location for the stiffness matrix values
    kcol = zeros(Int64, neDoF)  # [sparse] column location for the stiffness matrix values

    for nel in 1:nels
        nen = length(etpl[nel])
        neDoF = (nen * nD)^2
        kindex[nel+1] = kindex[nel] + neDoF
        ed = ones(Int, nD, 1) * reshape(etpl[nel] * nD, 1, nen) - collect(nD-1:-1:0) * ones(Int, 1, nen)
        ed = reshape(ed, 1, nen * nD)   # element degrees of freedom
        ke = zeros(Float64, nen * nD, nen * nD)
        felem = zeros(Float64, nen * nD)   # zero the internal force matrix
        if nen == 12
            ngp, N, dNr, wp = shapefunc(2, nen)
            if occursin(r"matzy[p,n]", mat[nel])
                JT = dNr[:, 1:6] * coord[etpl[nel][1:6], [1, 3]]
                T = [1 0 0; 0 0 1; 0 1 0]
            elseif occursin(r"matzx[p,n]", mat[nel])
                JT = dNr[:, 1:6] * coord[etpl[nel][1:6], [2, 3]]
                T = [0 1 0; 0 0 1; 1 0 0]
            end
            for gp in 1:ngp
                B = zeros(Float64, 3, nen * nD)
                Ks, Kn = gps[nel][gp].Ks, gps[nel][gp].Kn
                D = [Ks 0 0; 0 Ks 0; 0 0 Kn]
                indx = 2 * gp .- collect(1:-1:0)
                detJ = abs(det(JT[indx, :]))
                for (i, j) in Base.product(1:3, 1:6)
                    B[i, 3*j-3+i] = N[j, gp]
                end
                B[:, 19:end] = -copy(B[:, 1:18])
                # coordinate transformation (use a simplified way here)
                B = T * B
                z = N[1:6, gp]' * coord[etpl[nel][1:6], 3] - H  # depth of the Gauss point
                gps[nel][gp].sigG[3] = -z * dens * grav
                gps[nel][gp].sigG[1] = -K0 * z * dens * grav
                gps[nel][gp].sigG[2] = -K0 * z * dens * grav
                gps[nel][gp].sig = copy(gps[nel][gp].sigG)
                felem .+= (B' * gps[nel][gp].sig[1:3] * detJ * wp[gp])
                ke .+= (B' * D * B * detJ * wp[gp])
            end
        else
            ngp, N, dNr, wp = shapefunc(nD, nen)
            JT = dNr * coord[etpl[nel], :]
            for gp in 1:ngp                      # loop over the gauss points
                B = zeros(Float64, 6, nen * nD)
                D = elasticD(gps[nel][gp].E, gps[nel][gp].v)
                indx = nD * gp .- collect(nD-1:-1:0)  # index for the Gauss point entries
                detJ = det(JT[indx, :])               # determinant of the Jacobian matrix
                dNx = JT[indx, :] \ dNr[indx, :]
                # Assemble the strain-displacement matrix 
                B[[1 4 6], 1:nD:end] = dNx
                B[[4 2 5], 2:nD:end] = dNx
                B[[6 5 3], 3:nD:end] = dNx
                @debug "Strain-displacement matrix" B
                z = N[:, gp]' * coord[etpl[nel], 3] - H
                @debug "The global height of the Gauss point is" z
                gps[nel][gp].sigG[3] = -z * dens * grav
                gps[nel][gp].sigG[1] = -K0 * z * dens * grav
                gps[nel][gp].sigG[2] = -K0 * z * dens * grav
                gps[nel][gp].sig = copy(gps[nel][gp].sigG)
                felem .+= (B' * gps[nel][gp].sigG * detJ * wp[gp])
                ke .+= (B' * D * B * detJ * wp[gp])
                fb[ed[nD:nD:end]] += N[:, gp] * dens * grav * wp[gp] * detJ  # body force
            end
        end
        fint[ed] .+= felem'  # internal force for each element
        krow[kindex[nel]+1:kindex[nel+1]] = reshape(repeat(ed, 1, nen * nD), neDoF, 1)
        kcol[kindex[nel]+1:kindex[nel+1]] = reshape(repeat(ed, nen * nD, 1), neDoF, 1)
        kval[kindex[nel]+1:kindex[nel+1]] = reshape(ke, neDoF, 1)
    end
    uvw = zeros(Float64, nDoF)  # zero displacement at the geostatic state
    react = fint - fb

    mesh.uvw = uvw
    mesh.fint, mesh.fb, mesh.react = fint, fb, react
    mesh.kval, mesh.krow, mesh.kcol = kval, krow, kcol

    # initialise the stresses
    for nel in 1:nels
        for gp in 1:mesh.ngp[nel]
            # initialisation
            gps[nel][gp].H0 = zeros(Float64, 6)  # consistent with engnieering tensor, Bu
            gps[nel][gp].H = zeros(Float64, 6)
            gps[nel][gp].Y = false
            gps[nel][gp].sigP = zeros(Float64, 6)
            gps[nel][gp].eps = zeros(Float64, 6)
            gps[nel][gp].epsP = zeros(Float64, 6)
        end
    end

    return nothing
end

"""
    geostatic_excavation(mesh, gps, mesh1, mesh2, dens, grav, K0, H)

Set the geostatic conditions for the braced excavation using the simplified K0 approach. 
Return the equivalent nodal forces for each stage of excavation. 

# Arguments:
- `mesh`: MESH structure defined by [`BracedExcavation.MESH`](@ref).
- `gps`: GP structure defined by [`BracedExcavation.GP`](@ref).
- `mesh1`: mesh information for the first excavation.
- `mesh2`: mesh information for the second excavation.
- `dens`: density of materials.
- `grav`: gravitational forces.
- `K0`: K0 value.
- `H`: height of the model.
"""
function geostatic_excavation(mesh, gps, mesh1, mesh2, dens, grav, K0, h)
    mesh.lstp = 0
    coord = mesh.coord  # coordinate
    etpl = mesh.etpl  # topography
    mat = mesh.mat  # materials

    nodes, nD = size(coord)
    nels = length(etpl)
    nDoF = nodes * nD
    fb = zeros(Float64, nDoF)
    fint = zeros(Float64, nDoF)
    kindex = zeros(Int64, nels + 1)
    neDoF = 0
    for e in etpl
        neDoF += (length(e) * nD)^2
    end
    kval = zeros(Float64, neDoF) # values for the global stiffness matrix
    krow = zeros(Int64, neDoF)  # [sparse] row location for the stiffness matrix values
    kcol = zeros(Int64, neDoF)  # [sparse] column location for the stiffness matrix values

    # to obtain the equivalent nodal forces for each excavation stage
    etple1 = mesh1[1]  # elements to be excavated in stage 1
    benode1 = mesh1[2]  # excavation boundary of stage 1
    etple2 = mesh2[1]  # elements to be excavated in stage 2
    benode2 = mesh2[2]  # excavation boundary of stage 2
    fbe1 = zeros(Float64, nDoF)
    finte1 = zeros(Float64, nDoF)
    fbe2 = zeros(Float64, nDoF)
    finte2 = zeros(Float64, nDoF)

    for nel in 1:nels
        nen = length(etpl[nel])
        neDoF = (nen * nD)^2
        kindex[nel+1] = kindex[nel] + neDoF
        ed = ones(Int, nD, 1) * reshape(etpl[nel] * nD, 1, nen) - collect(nD-1:-1:0) * ones(Int, 1, nen)
        ed = reshape(ed, 1, nen * nD)   # element degrees of freedom
        ke = zeros(Float64, nen * nD, nen * nD)
        felem = zeros(Float64, nen * nD)   # zero the internal force matrix
        if nen == 12
            ngp, N, dNr, wp = shapefunc(2, nen)
            if occursin(r"matzy[p,n]", mat[nel])
                JT = dNr[:, 1:6] * coord[etpl[nel][1:6], [1, 3]]
                T = [1 0 0; 0 0 1; 0 1 0]
            elseif occursin(r"matzx[p,n]", mat[nel])
                JT = dNr[:, 1:6] * coord[etpl[nel][1:6], [2, 3]]
                T = [0 1 0; 0 0 1; 1 0 0]
            end
            for gp in 1:ngp
                B = zeros(Float64, 3, nen * nD)
                Ks, Kn = gps[nel][gp].Ks, gps[nel][gp].Kn
                D = [Ks 0 0; 0 Ks 0; 0 0 Kn]
                indx = 2 * gp .- collect(1:-1:0)
                detJ = abs(det(JT[indx, :]))
                for (i, j) in Base.product(1:3, 1:6)
                    B[i, 3*j-3+i] = N[j, gp]
                end
                B[:, 19:end] = -copy(B[:, 1:18])
                # coordinate transformation (use a simplified way here)
                B = T * B
                z = N[1:6, gp]' * coord[etpl[nel][1:6], 3] - h  # depth of the Gauss point
                gps[nel][gp].sigG[3] = -z * dens * grav
                gps[nel][gp].sigG[1] = -K0 * z * dens * grav
                gps[nel][gp].sigG[2] = -K0 * z * dens * grav
                gps[nel][gp].sig = copy(gps[nel][gp].sigG)
                felem .+= (B' * gps[nel][gp].sig[1:3] * detJ * wp[gp])
                ke .+= (B' * D * B * detJ * wp[gp])
            end
        else
            ngp, N, dNr, wp = shapefunc(nD, nen)
            JT = dNr * coord[etpl[nel], :]
            for gp in 1:ngp                      # loop over the gauss points
                B = zeros(Float64, 6, nen * nD)
                D = elasticD(gps[nel][gp].E, gps[nel][gp].v)
                indx = nD * gp .- collect(nD-1:-1:0)  # index for the Gauss point entries
                detJ = det(JT[indx, :])               # determinant of the Jacobian matrix
                dNx = JT[indx, :] \ dNr[indx, :]
                # Assemble the strain-displacement matrix 
                B[[1 4 6], 1:nD:end] = dNx
                B[[4 2 5], 2:nD:end] = dNx
                B[[6 5 3], 3:nD:end] = dNx
                @debug "Strain-displacement matrix" B
                z = N[:, gp]' * coord[etpl[nel], 3] - h
                @debug "The global height of the Gauss point is" z
                gps[nel][gp].sigG[3] = -z * dens * grav
                gps[nel][gp].sigG[1] = -K0 * z * dens * grav
                gps[nel][gp].sigG[2] = -K0 * z * dens * grav
                gps[nel][gp].sig = copy(gps[nel][gp].sigG)
                felem .+= (B' * gps[nel][gp].sigG * detJ * wp[gp])
                ke .+= (B' * D * B * detJ * wp[gp])
                fb[ed[nD:nD:end]] += N[:, gp] * dens * grav * wp[gp] * detJ  # body force
                if etpl[nel] in etple1  # obtain the equivalent nodal force after the first excavation
                    for node = 1:nen
                        if sum(etpl[nel][node] .∈ benode1) > 0
                            fbe1[etpl[nel][node]*nD] += N[node, gp] * dens * grav * wp[gp] * detJ
                        end
                    end
                end
                if etpl[nel] in etple2  # obtain the equivalent nodal force after the second excavation
                    for node = 1:nen
                        if sum(etpl[nel][node] .∈ benode2) > 0
                            fbe2[etpl[nel][node]*nD] += N[node, gp] * dens * grav * wp[gp] * detJ
                        end
                    end
                end
            end
        end
        fint[ed] .+= felem'  # internal force for each element
        if etpl[nel] in etple1
            for node = 1:nen
                if sum(etpl[nel][node] .∈ benode1) > 0
                    if nen == 12
                        finte1[etpl[nel][node]*nD-2:etpl[nel][node]*nD-1] .-= felem[node*nD-2:node*nD-1]
                        finte1[etpl[nel][node]*nD] += felem[node*nD]
                        # finte1[etpl[nel][node]*nD-2:etpl[nel][node]*nD] .-= felem[node*nD-2:node*nD]
                    else
                        finte1[etpl[nel][node]*nD-2:etpl[nel][node]*nD] .+= felem[node*nD-2:node*nD]
                    end
                end
            end
        end
        if etpl[nel] in etple2
            for node = 1:nen
                if sum(etpl[nel][node] .∈ benode2) > 0
                    if nen == 12
                        finte2[etpl[nel][node]*nD-2:etpl[nel][node]*nD-1] .-= felem[node*nD-2:node*nD-1]
                        finte2[etpl[nel][node]*nD] += felem[node*nD]
                        # finte2[etpl[nel][node]*nD-2:etpl[nel][node]*nD] .-= felem[node*nD-2:node*nD]
                    else
                        finte2[etpl[nel][node]*nD-2:etpl[nel][node]*nD] .+= felem[node*nD-2:node*nD]
                    end
                end
            end
        end
        krow[kindex[nel]+1:kindex[nel+1]] = reshape(repeat(ed, 1, nen * nD), neDoF, 1)
        kcol[kindex[nel]+1:kindex[nel+1]] = reshape(repeat(ed, nen * nD, 1), neDoF, 1)
        kval[kindex[nel]+1:kindex[nel+1]] = reshape(ke, neDoF, 1)
    end
    uvw = zeros(Float64, nDoF)  # zero displacement at the geostatic state
    react = fint - fb

    mesh.uvw = uvw
    mesh.fint, mesh.fb, mesh.react = fint, fb, react
    mesh.kval, mesh.krow, mesh.kcol = kval, krow, kcol

    # initialise the stresses
    for nel in 1:nels
        for gp in 1:mesh.ngp[nel]
            # initialisation
            gps[nel][gp].H0 = zeros(Float64, 6)  # consistent with engnieering tensor, Bu
            gps[nel][gp].H = zeros(Float64, 6)
            gps[nel][gp].Y = false
            gps[nel][gp].sigP = zeros(Float64, 6)
            gps[nel][gp].eps = zeros(Float64, 6)
            gps[nel][gp].epsP = zeros(Float64, 6)
        end
    end

    # For checking the equivalent nodal force
    # fv = finte1 - fbe1
    # @show "Equivalent nodal force of the first excavation"
    # @show sum(fv[1:3:end])
    # @show sum(fv[2:3:end])
    # @show sum(fv[3:3:end])
    # fv = finte2 - fbe2
    # @show "Equivalent nodal force of the second excavation"
    # @show sum(fv[1:3:end])
    # @show sum(fv[2:3:end])
    # @show sum(fv[3:3:end])

    return [finte1, fbe1], [finte2, fbe2], kindex
end

