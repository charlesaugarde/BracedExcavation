"""
    boundary_excavation(coordinate, W, L)

Set the Dirichlet (displacement) boundary conditions for the excavation model.

# Arguments:
- `coordinate`: coordinates of model nodes.
- `W`: width of the model.
- `L`: length of the model. 
"""
function boundary_excavation(coordinate, W, L)
    # roller boundary condition for the symmetric condition
    nnodes, nD = size(coordinate)
    nDoF = Int(nnodes * nD)
    bc = zeros(Int64, nDoF, 2)  # boundary condition
    for node = 1:nnodes
        if coordinate[node, 1] <= 1e-5 || abs(coordinate[node, 1] - W) <= 1e-5
            bc[node*nD-2, :] = [node * nD - 2, 0]
        end
        if coordinate[node, 2] <= 1e-5 || abs(coordinate[node, 2] - L) <= 1e-5
            bc[node*nD-1, :] = [node * nD - 1, 0]
        end
        if coordinate[node, 3] <= 1e-5
            bc[node*nD-0, :] = [node * nD - 0, 0]
            bc[node*nD-1, :] = [node * nD - 1, 0]
            bc[node*nD-2, :] = [node * nD - 2, 0]
        end
    end

    # reduce the storage
    bc = bc[vec(mapslices(col -> any(col .!= 0), bc, dims=2)), :]

    return bc
end

"""
    boundary_compression(coord, be2, be3, p, u, lstpx)

Set the boundary conditions for the triaxial tests: the bottom and two sides are roller conditions. 
The boundaries with names *sidex*, and *sidey* are kept constnat with pressure **p**. 
During the compression stage, the *top* boundary is compressed to **u** within **lstpx** steps. 

# Arguments:
- `coord`: coordinates of model nodes.
- `be2`: boundary nodes of *sidex*.
- `be3`: boundary nodes of *sidey*.
- `p`: confining pressure **p**. 
- `u`: total displacement of *top* boundary. 
- `lstpx`: load steps for the compression. 
"""
function boundary_compression(coord, be2, be3, p, u, lstpx)
    nnodes, nD = size(coord)
    nDoF = Int(nnodes * nD)
    bc = zeros(Float64, nDoF, 2)
    for node = 1:nnodes
        if coord[node, 1] == 0
            bc[node*nD-2, :] = [node * nD - 2, 0]
        end
        if coord[node, 2] == 0
            bc[node*nD-1, :] = [node * nD - 1, 0]
        end
        if coord[node, 3] == 0
            bc[node*nD-0, :] = [node * nD - 0, 0]
        end
        if coord[node, 3] == 1
            bc[node*nD-0, :] = [node * nD - 0, u / lstpx]
        end
    end
    bc = bc[vec(mapslices(col -> any(col .!= 0), bc, dims=2)), :]

    # external force
    fext = zeros(Float64, nDoF)
    for elem in be2
        nen = length(elem)
        ngp, N, dNr, wp = shapefunc(2, nen)  # 2D shape function
        JT = dNr * coord[elem, 2:3]  # side pressure
        for node in 1:nen
            nodf = 0.0
            for gp = 1:ngp
                indx = 2 * gp .- collect(1:-1:0)
                detJ = abs(det(JT[indx, :]))
                nodf += N[node, gp] * p * detJ * wp[gp]
            end
            fext[elem[node]*nD-2] += nodf
        end
    end
    for elem in be3
        nen = length(elem)
        ngp, N, dNr, wp = shapefunc(2, nen)  # 2D shape function
        JT = dNr * coord[elem, [1, 3]]  # side pressure
        for node in 1:nen
            nodf = 0.0
            for gp = 1:ngp
                indx = 2 * gp .- collect(1:-1:0)
                detJ = abs(det(JT[indx, :]))
                nodf += N[node, gp] * p * detJ * wp[gp]
            end
            fext[elem[node]*nD-1] += nodf
        end
    end

    return bc, fext
end

"""
    boundary_consolidation(coord, be1, be2, be3, p)

Set the boundary conditions for the triaxial tests: the bottom and two sides are roller conditions. 
The boundaries with names *top*, *sidex*, and *sidey* are consolidated to the prescribed pressure **p**. 

# Arguments:
- `coord`: coordinates of model nodes.
- `be1`: boundary nodes of *top*. 
- `be2`: boundary nodes of *sidex*.
- `be3`: boundary nodes of *sidey*.
- `p`: confining pressure **p**. 
"""
function boundary_consolidation(coord, be1, be2, be3, p)
    nnodes, nD = size(coord)
    nDoF = Int(nnodes * nD)
    bc = zeros(Int64, nDoF, 2)
    for node = 1:nnodes
        if coord[node, 1] == 0
            bc[node*nD-2, :] = [node * nD - 2, 0]
        end
        if coord[node, 2] == 0
            bc[node*nD-1, :] = [node * nD - 1, 0]
        end
        if coord[node, 3] == 0
            bc[node*nD-0, :] = [node * nD - 0, 0]
        end
    end
    bc = bc[vec(mapslices(col -> any(col .!= 0), bc, dims=2)), :]

    # external force
    fext = zeros(Float64, nDoF)
    for elem in be1
        nen = length(elem)
        ngp, N, dNr, wp = shapefunc(2, nen)  # 2D shape function
        JT = dNr * coord[elem, 1:2]  # top pressure
        for node in 1:nen
            nodf = 0.0
            for gp = 1:ngp
                indx = 2 * gp .- collect(1:-1:0)
                detJ = abs(det(JT[indx, :]))
                nodf += N[node, gp] * p * detJ * wp[gp]
            end
            fext[elem[node]*nD] += nodf
        end
    end
    for elem in be2
        nen = length(elem)
        ngp, N, dNr, wp = shapefunc(2, nen)  # 2D shape function
        JT = dNr * coord[elem, 2:3]  # side pressure
        for node in 1:nen
            nodf = 0.0
            for gp = 1:ngp
                indx = 2 * gp .- collect(1:-1:0)
                detJ = abs(det(JT[indx, :]))
                nodf += N[node, gp] * p * detJ * wp[gp]
            end
            fext[elem[node]*nD-2] += nodf
        end
    end
    for elem in be3
        nen = length(elem)
        ngp, N, dNr, wp = shapefunc(2, nen)  # 2D shape function
        JT = dNr * coord[elem, [1, 3]]  # side pressure
        for node in 1:nen
            nodf = 0.0
            for gp = 1:ngp
                indx = 2 * gp .- collect(1:-1:0)
                detJ = abs(det(JT[indx, :]))
                nodf += N[node, gp] * p * detJ * wp[gp]
            end
            fext[elem[node]*nD-1] += nodf
        end
    end

    return bc, fext
end
