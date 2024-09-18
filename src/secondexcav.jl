function secondexcavation(mesh, kindex0, fexc, mesh0, mesh1, mesh2, gps, fexcgeo)
    mesh.lstp = 2  # set the load step for the excavation

    # update COORD related parameters
    coord20 = mesh2[4]  # remaining nodes after the excavation
    _, nD = size(mesh.coord)  # dimension (3) to be used for the f values
    nnodes = length(coord20)  # total number of nodes after excavation

    fint0 = fexc[1]  # internal force of the remaining model
    fb0 = fexc[2]  # body force of the remaining model

    fext0 = fexcgeo[1] - fexcgeo[2] # net equivalent nodal force due to excavation

    react0 = mesh.react  # reaction force
    uvw0 = mesh.uvw  # initial displacement

    # initialise the coordinates, f values, and displacement (uvw)
    coord = zeros(Float64, nnodes, nD)
    fext = zeros(Float64, nnodes * nD)
    fint = zeros(Float64, nnodes * nD)
    fb = zeros(Float64, nnodes * nD)
    react = zeros(Float64, nnodes * nD)
    uvw = zeros(Float64, nnodes * nD)

    # transform the coord-related parameters
    for (i, node) in enumerate(coord20)
        # mesh1[4] remaining nodes after the first excavation
        indx = findfirst(item -> item == node, mesh1[4])
        coord[i, :] = mesh.coord[indx, :]  # the coordinates are decided from previous load step
        @assert mesh0[1][node, :] == mesh.coord[indx, :] "coordinate transformation is wrong, $node"
        fext[i*nD-nD+1:i*nD] = fext0[node*nD-2:node*nD]  # transform from initial order to current order
        fint[i*nD-nD+1:i*nD] = fint0[indx*nD-2:indx*nD]
        fb[i*nD-nD+1:i*nD] = fb0[indx*nD-2:indx*nD]
        react[i*nD-nD+1:i*nD] = react0[indx*nD-2:indx*nD]
        uvw[i*nD-nD+1:i*nD] = uvw0[indx*nD-2:indx*nD]
    end

    # redefine the MESH struct data for the FE calculation
    mesh.coord = coord
    mesh.fext = fext
    mesh.fint = fint
    mesh.react = react
    mesh.fb = fb
    mesh.uvw = uvw

    # update ELEMENT related parameters
    etpl10map = mesh1[end]  # mapping the 3D elements of the previous step
    etpl20map = mesh2[end]  # mapping the 3D elements of the remaining model
    etpl20 = mesh2[3]  # remaining elements after excavation
    nels = length(etpl20)  # total number of elements after the excavation

    # calculation the element DoF for the stiffness matrix
    neDoF = 0
    for e in etpl20
        neDoF += (length(e) * nD)^2
    end

    # initialise the element stiffness matrices
    kval = zeros(Float64, neDoF) # values for the global stiffness matrix
    krow = zeros(Int64, neDoF)  # [sparse] row location for the stiffness matrix values
    kcol = zeros(Int64, neDoF)  # [sparse] column location for the stiffness matrix values

    kindex = zeros(Int64, nels + 1)  # location for each element
    indx = Int64[]  # remaining elements after the first excavation

    # transform the element-related parameters
    for nel = 1:nels
        nen = length(etpl20[nel])
        neDoF = (nen * nD)^2
        kindex[nel+1] = kindex[nel] + neDoF  # location of the stiffness matrix for the nel element
        for node = 1:nen
            # node number is sorted based on "coord20"
            etpl20[nel][node] = findfirst(item -> item == etpl20[nel][node], coord20)
        end
        ed = ones(Int, nD, 1) * reshape(etpl20[nel] * nD, 1, nen) - collect(nD-1:-1:0) * ones(Int, 1, nen)
        ed = reshape(ed, 1, nen * nD)   # element degrees of freedom
        krow[kindex[nel]+1:kindex[nel+1]] = reshape(repeat(ed, 1, nen * nD), neDoF, 1)
        kcol[kindex[nel]+1:kindex[nel+1]] = reshape(repeat(ed, nen * nD, 1), neDoF, 1)
        # convert the element orders
        localindx = findfirst(item -> item == etpl20map[nel], etpl10map)
        if nen == 10
            @assert mesh.ngp[Int(localindx)] == 4 "MESH npg is wrong, $nel"
            @assert length(gps[Int(localindx)]) == 4 "GP ngp is wrong, $nel"
        elseif nen == 15
            @assert mesh.ngp[Int(localindx)] == 6 "MESH npg is wrong, $nel"
            @assert length(gps[Int(localindx)]) == 6 "GP ngp is wrong, $nel"
        end
        push!(indx, localindx)
        kval[kindex[nel]+1:kindex[nel+1]] = mesh.kval[Int64(kindex0[localindx] + 1):Int64(kindex0[localindx+1])]
    end

    mats = copy(mesh.mat)
    ngps = copy(mesh.ngp)
    mat = mats[indx]
    ngp = ngps[indx]

    # redefine the MESH struct data
    mesh.etpl = etpl20
    mesh.mat = mat
    mesh.ngp = ngp
    mesh.krow = krow
    mesh.kcol = kcol
    mesh.kval = kval

    # redefine the element-related stress parameters
    gps1 = gps[indx]

    return gps1, kindex

end
