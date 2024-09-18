function firstexcavation(mesh, kindex0, fexc1, mesh0, mesh1, gps, mesh2)
    mesh.lstp = 1  # set the load step for the excavation

    # update COORD related parameters
    coord10 = mesh1[4]  # remaining nodes after the excavation
    _, nD = size(mesh.coord)  # dimension (3) to be used for the f values
    nnodes = length(coord10)  # total number of nodes after excavation

    finte0 = fexc1[1]  # equivalent internal force due to excavation
    fbe0 = fexc1[2]  # equivalent body force due to excavation
    fext0 = finte0 - fbe0  # net equivalent nodal force due to excavation

    fint0 = mesh.fint - finte0  # internal force minus the excavation force
    fb0 = mesh.fb - fbe0  # body force minus the excavation force
    uvw0 = mesh.uvw  # initial displacement

    # initialise the coordinates, f values and the displacement (uvw)
    coord = zeros(Float64, nnodes, nD)
    fext = zeros(Float64, nnodes * nD)
    fint = zeros(Float64, nnodes * nD)
    fb = zeros(Float64, nnodes * nD)
    uvw = zeros(Float64, nnodes * nD)

    # transform the coord-related parameters
    for (i, node) in enumerate(coord10)
        coord[i, :] = mesh.coord[node, :]
        fext[i*nD-nD+1:i*nD] = fext0[node*nD-2:node*nD]
        fint[i*nD-nD+1:i*nD] = fint0[node*nD-2:node*nD]
        fb[i*nD-nD+1:i*nD] = fb0[node*nD-2:node*nD]
        uvw[i*nD-nD+1:i*nD] = uvw0[node*nD-2:node*nD]
    end

    # redefine the MESH struct data for the FE calculation
    mesh.coord = coord
    mesh.fext = fext
    mesh.fint = fint
    mesh.fb = fb
    mesh.react = fint - fb
    mesh.uvw = uvw

    # update ELEMENT related parameters
    etpl0map = mesh0[end]  # mapping the 3D elements of the whole model
    etpl10map = mesh1[end]  # mapping the 3D elements of the remaining model
    etpl10 = mesh1[3]  # remaining elements after excavation
    nels = length(etpl10)  # total number of elements after the excavation

    # calculate the element DoF for the stiffness matrix
    neDoF = 0
    for e in etpl10
        neDoF += (length(e) * nD)^2
    end

    # initialise the element stiffness matrices
    kval = zeros(Float64, neDoF) # values for the global stiffness matrix
    krow = zeros(Int64, neDoF)  # [sparse] row location for the stiffness matrix values
    kcol = zeros(Int64, neDoF)  # [sparse] column location for the stiffness matrix values

    kindex = zeros(Int64, nels + 1)  # location for each element
    indx = Int64[]  # remaining elements after the excavation

    # transform the element-related parameters
    for nel = 1:nels
        nen = length(etpl10[nel])
        neDoF = (nen * nD)^2
        kindex[nel+1] = kindex[nel] + neDoF  # location of the stiffness matrix for the nel element
        for node = 1:nen
            # node number is sorted based on "coord10"
            etpl10[nel][node] = findfirst(item -> item == etpl10[nel][node], coord10)
        end
        ed = ones(Int, nD, 1) * reshape(etpl10[nel] * nD, 1, nen) - collect(nD-1:-1:0) * ones(Int, 1, nen)
        ed = reshape(ed, 1, nen * nD)   # element degrees of freedom
        krow[kindex[nel]+1:kindex[nel+1]] = reshape(repeat(ed, 1, nen * nD), neDoF, 1)
        kcol[kindex[nel]+1:kindex[nel+1]] = reshape(repeat(ed, nen * nD, 1), neDoF, 1)
        # convert the element orders
        localindx = findfirst(item -> item == etpl10map[nel], etpl0map) # remaining element number after excavation
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

    # redefine the MESH struct data
    mesh.etpl = etpl10
    mesh.mat = mesh.mat[indx]
    mesh.ngp = mesh.ngp[indx]
    mesh.krow = krow
    mesh.kcol = kcol
    mesh.kval = kval

    # redefine the element-related stress parameters
    gps1 = gps[indx]

    # change the process to record the fint and fb of remaining model
    # This part is necessary for conducting the subsequent excavation to get the corrent fb and fint on the excavation interface
    # pay attention to the "deepcopy" to avoid change the values in the mesh2
    benode1 = deepcopy(mesh2[4])  # remaning nodes of the second excavation
    benode = zeros(Int64, length(benode1))  # initialise the nodes number
    for (i, node) in enumerate(benode1)
        benode[i] = findfirst(item -> item == node, coord10) # remaining node number after excavation
    end
    etpl1 = deepcopy(mesh2[3])  # remaining elements of the second excavation  
    nels = length(etpl1)
    for nel = 1:nels
        nen = length(etpl1[nel])
        for node = 1:nen
            etpl1[nel][node] = findfirst(item -> item == etpl1[nel][node], coord10)  # remaining nodes after excavation
        end
    end

    return gps1, kindex, [etpl1, benode]
end
