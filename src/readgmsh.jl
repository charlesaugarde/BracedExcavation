using Gmsh: gmsh

"""
    readgmsh_triaxial(argv)

Read the mesh data for the triaxial test, needs *top*, *sidex*, *sidey* boundaries. 

# Arguments:
- `argv`: the input mesh file with suffix, e.g., "ElementCode.msh"
"""
function readgmsh_triaxial(argv)
    gmsh.initialize()
    gmsh.open(argv)

    # read all the nodes
    nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
    Tnnodes = length(nodeTags)
    coord0 = zeros(Float64, Tnnodes, 3)
    for i = 1:Tnnodes
        coord0[nodeTags[i], :] = nodeCoords[(i-1)*3+1:i*3]
    end

    # read all the elements with material properties
    etpl0 = Vector[]   # total 3D elements for the initial body force
    entities = gmsh.model.getEntities(3)
    for e in entities
        physicalTags = gmsh.model.getPhysicalGroupsForEntity(e[1], e[2])
        for p in physicalTags
            name = gmsh.model.getPhysicalName(e[1], p)
            if occursin("tet", name)
                _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                numTelem, element = interpretElement(elemTags, elemNodeTags)
                element = view(element, :, [1, 2, 3, 4, 5, 6, 7, 8, 10, 9])
                for mat in ["LondonClay"]
                    if occursin(mat, name)
                        for i = 1:numTelem
                            push!(etpl0, [mat, element[i, :], 4])
                        end
                    end
                end
            elseif occursin("wed", name)
                _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                numTelem, element = interpretElement(elemTags, elemNodeTags)
                element = view(element, :, [1, 2, 3, 4, 5, 6, 7, 10, 8, 13, 15, 14, 9, 11, 12])
                for mat in ["LondonClay"]
                    if occursin(mat, name)
                        for i = 1:numTelem
                            push!(etpl0, [mat, element[i, :], 6])
                        end
                    end
                end
            elseif occursin("hex", name)
                _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                numTelem, element = interpretElement(elemTags, elemNodeTags)
                element = view(element, :, [1, 2, 3, 4, 5, 6, 7, 8])
                for mat in ["LondonClay"]
                    if occursin(mat, name)
                        for i = 1:numTelem
                            push!(etpl0, [mat, element[i, :], 8])
                        end
                    end
                end
            end
        end
    end

    # return the 2D boundary conditions
    be1 = Vector{Int}[]
    be2 = Vector{Int}[]
    be3 = Vector{Int}[]
    entities = gmsh.model.getEntities(2)
    for e in entities
        # get the physical groups
        physicalTags = gmsh.model.getPhysicalGroupsForEntity(e[1], e[2])
        if length(physicalTags) > 0
            for p in physicalTags
                name = gmsh.model.getPhysicalName(2, p)
                if occursin("top", name)
                    _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                    numTelem, element = interpretElement(elemTags, elemNodeTags)
                    for i = 1:numTelem
                        push!(be1, element[i, :])
                    end
                elseif occursin("sidex", name)
                    _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                    numTelem, element = interpretElement(elemTags, elemNodeTags)
                    for i = 1:numTelem
                        push!(be2, element[i, :])
                    end
                elseif occursin("sidey", name)
                    _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                    numTelem, element = interpretElement(elemTags, elemNodeTags)
                    for i = 1:numTelem
                        push!(be3, element[i, :])
                    end
                end
            end
        end
    end

    return coord0, etpl0, [be1, be2, be3]
end

"""
    readgmsh_excavation(argv)

Read the mesh data for the excavation test, need the material physical groups and interface elements.

# Arguments:
- `argv`: the input mesh file with suffix, e.g., "BracedExcavation.msh"
"""
function readgmsh_excavation(argv)
    gmsh.initialize()
    gmsh.open(argv)

    # read the nodes based on the interface physical groups
    coordi = Int[]  # interface nodes
    nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(2, 10)  # interface
    Snnodes = length(nodeTags)
    for i = 1:Snnodes
        push!(coordi, reinterpret(Int, nodeTags[i]))
    end

    # get the maximum node tags
    maxNodeTag = reinterpret(Int, gmsh.model.mesh.getMaxNodeTag())
    @debug "The max tag of node in original model is $maxNodeTag"

    # duplicate the interface nodes from the original model
    nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
    Tnnodes = length(nodeTags)
    @debug "The total nodes in the model: $Tnnodes"
    @debug "The nodes in the interface is: $Snnodes"
    coord0 = zeros(Float64, Tnnodes + Snnodes, 3)
    for i = 1:Tnnodes
        coord0[nodeTags[i], :] = nodeCoords[(i-1)*3+1:i*3]
    end
    for i = 1:Snnodes
        coord0[maxNodeTag+i, :] = coord0[coordi[i], :]
    end
    @debug "The nodes adding the interface elements are" coord0
    @debug "The additional nodes are" coord0[end-Snnodes:end, :]

    # get the 2D elements in the interface
    inter0 = Vector{Int}[]  # the total interface elements
    # the interface is read based on the geomaker
    inter_01_xp = Vector{Int}[]  # excavated interface at the first stage in x-p
    inter_01_yp = Vector{Int}[]  # excavated interface at the first stage in y-p
    inter_02_xp = Vector{Int}[]  # excavated interface at the second stage in x-p
    inter_02_yp = Vector{Int}[]  # excavated interface at the second stage in y-p
    inter_xp = Vector{Int}[]  # remaining interface in x-p direction of excavation side
    inter_yp = Vector{Int}[]  # remaining interface in y-p direction of excavation side
    inter_xn = Vector{Int}[]  # remaining interface in x-n direction of unexcavation side
    inter_yn = Vector{Int}[]  # remaining interface in y-n direction of unexcavation side
    entities = gmsh.model.getEntities(2)
    for e in entities
        physicalTags = gmsh.model.getPhysicalGroupsForEntity(e[1], e[2])
        if length(physicalTags) > 0
            for p in physicalTags
                name = gmsh.model.getPhysicalName(2, p)
                if occursin("inter0", name)
                    @debug "Boundary conditions - Physical groups: $name"
                    _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                    numTelem, element = interpretElement(elemTags, elemNodeTags)
                    @debug "reshaped elements" element
                    for i = 1:numTelem
                        push!(inter0, element[i, :])
                    end
                elseif occursin("inter_01_xp", name)
                    _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                    numTelem, element = interpretElement(elemTags, elemNodeTags)
                    for i = 1:numTelem
                        push!(inter_01_xp, element[i, :])
                    end
                elseif occursin("inter_01_yp", name)
                    _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                    numTelem, element = interpretElement(elemTags, elemNodeTags)
                    for i = 1:numTelem
                        push!(inter_01_xp, element[i, :])
                    end
                elseif occursin("inter_02_xp", name)
                    _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                    numTelem, element = interpretElement(elemTags, elemNodeTags)
                    for i = 1:numTelem
                        push!(inter_02_xp, element[i, :])
                    end
                elseif occursin("inter_02_yp", name)
                    _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                    numTelem, element = interpretElement(elemTags, elemNodeTags)
                    for i = 1:numTelem
                        push!(inter_02_yp, element[i, :])
                    end
                elseif occursin("inter_xp", name)
                    _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                    numTelem, element = interpretElement(elemTags, elemNodeTags)
                    for i = 1:numTelem
                        push!(inter_xp, element[i, :])
                    end
                elseif occursin("inter_yp", name)
                    _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                    numTelem, element = interpretElement(elemTags, elemNodeTags)
                    for i = 1:numTelem
                        push!(inter_yp, element[i, :])
                    end
                elseif occursin("inter_xn", name)
                    _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                    numTelem, element = interpretElement(elemTags, elemNodeTags)
                    for i = 1:numTelem
                        push!(inter_xn, element[i, :])
                    end
                elseif occursin("inter_yn", name)
                    _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                    numTelem, element = interpretElement(elemTags, elemNodeTags)
                    for i = 1:numTelem
                        push!(inter_yn, element[i, :])
                    end
                end
            end
        end
    end
    @debug "The 2D element in the interface is" inter0
    @debug "The order of the newly generated nodes is" coordi

    # 2) connect the original interface elements with the newly generated elements
    numISelem = length(inter0)
    inter1 = Int[]
    for i in inter0
        for j in i
            indx = findfirst(item -> item == j, coordi)
            push!(inter1, indx + maxNodeTag)
        end
    end
    inter1 = reshape(inter1, Int(length(inter1) / length(inter0)), length(inter0))
    inter1 = inter1'
    @debug "The newly generated 2D mesh" inter1

    # 3) read the 3D elements with material properties
    etpl0 = Vector[]   # total 3D elements for the initial body force
    etpl0map = Int[]   # element tags of the total model

    etpl1 = Vector[]   # elements of the first excavation 
    etpl10 = Vector[]  # remaining elements after the first excavation
    etpl1map = Int[]  # remaining element tags after the first excavation

    etpl2 = Vector[]   # elements of the second excavation 
    etpl20 = Vector[]  # remaining elements after the second excavation
    etpl2map = Int[]  # remaining element tags after the second excavation

    entities = gmsh.model.getEntities(3)
    for e in entities
        physicalTags = gmsh.model.getPhysicalGroupsForEntity(e[1], e[2])
        for p in physicalTags
            name = gmsh.model.getPhysicalName(e[1], p)
            if occursin("tet", name)
                _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                numTelem, element = interpretElement(elemTags, elemNodeTags)
                element = view(element, :, [1, 2, 3, 4, 5, 6, 7, 8, 10, 9])
                for mat in ["LondonClay"]
                    if occursin(mat, name)
                        for i = 1:numTelem
                            push!(etpl0, [mat, element[i, :], 4])
                            push!(etpl0map, reinterpret(Int, elemTags[1][i]))
                        end
                    end
                end
                if occursin("excav01", name)  # excavated soils are modelled as tetrahedron
                    for i = 1:numTelem
                        push!(etpl1, element[i, :])
                    end
                end
                if occursin("rem01", name)
                    for i = 1:numTelem
                        push!(etpl10, element[i, :])
                        push!(etpl1map, reinterpret(Int, elemTags[1][i]))
                    end
                end
                if occursin("excav02", name)  # excavated soils are modelled as tetrahedron
                    for i = 1:numTelem
                        push!(etpl2, element[i, :])
                    end
                end
                if occursin("rem02", name)
                    for i = 1:numTelem
                        push!(etpl20, element[i, :])
                        push!(etpl2map, reinterpret(Int, elemTags[1][i]))
                    end
                end
            elseif occursin("wed", name)
                _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                numTelem, element = interpretElement(elemTags, elemNodeTags)
                element = view(element, :, [1, 2, 3, 4, 5, 6, 7, 10, 8, 13, 15, 14, 9, 11, 12])
                # revise the elemental nodes for the interface
                for i = 1:numTelem
                    for k in 1:15
                        indx = findfirst(item -> item == element[i, k], coordi)
                        if indx !== nothing
                            element[i, k] = maxNodeTag + indx
                        end
                    end
                end
                for mat in ["DiaphragmWall"]
                    if occursin(mat, name)
                        for i = 1:numTelem
                            push!(etpl0, [mat, element[i, :], 6])
                            push!(etpl0map, reinterpret(Int, elemTags[1][i]))
                        end
                    end
                end
                if occursin("rem01", name)
                    for i = 1:numTelem
                        push!(etpl10, element[i, :])
                        push!(etpl1map, reinterpret(Int, elemTags[1][i]))
                    end
                end
                if occursin("rem02", name)
                    for i = 1:numTelem
                        push!(etpl20, element[i, :])
                        push!(etpl2map, reinterpret(Int, elemTags[1][i]))
                    end
                end
            end
        end
    end

    # 4) add the interface elements to the model
    # add the mappings to the interface elements
    maxElemTag = reinterpret(Int, gmsh.model.mesh.getMaxElementTag())
    @debug "The max tag of node in original model is $maxElemTag"
    for i in 1:numISelem
        push!(etpl0map, maxElemTag + i)
        push!(etpl1map, maxElemTag + i)
        push!(etpl2map, maxElemTag + i)
        if inter0[i] in inter_01_xp  # to be excavated at the first stage in x-p direction
            push!(etpl0, ["matzxp", reduce(vcat, (inter0[i], inter1[i, :])), 3])
            push!(etpl1, reduce(vcat, (inter0[i], inter1[i, :])))
            push!(etpl2, reduce(vcat, (inter0[i], inter1[i, :])))
        elseif inter0[i] in inter_01_yp  # to be excavated at the first stage in y-p direction
            push!(etpl0, ["matzyp", reduce(vcat, (inter0[i], inter1[i, :])), 3])
            push!(etpl1, reduce(vcat, (inter0[i], inter1[i, :])))
            push!(etpl2, reduce(vcat, (inter0[i], inter1[i, :])))
        elseif inter0[i] in inter_02_xp  # to be excavated at the second stage in x-p direction
            push!(etpl0, ["matzxp", reduce(vcat, (inter0[i], inter1[i, :])), 3])
            push!(etpl2, reduce(vcat, (inter0[i], inter1[i, :])))
            push!(etpl10, reduce(vcat, (inter0[i], inter1[i, :])))
        elseif inter0[i] in inter_02_yp  # to be excavated at the second stage in y-p direction
            push!(etpl0, ["matzyp", reduce(vcat, (inter0[i], inter1[i, :])), 3])
            push!(etpl2, reduce(vcat, (inter0[i], inter1[i, :])))
            push!(etpl10, reduce(vcat, (inter0[i], inter1[i, :])))
        elseif inter0[i] in inter_xp  # remain at all stages in x-p direction
            push!(etpl0, ["matzxp", reduce(vcat, (inter0[i], inter1[i, :])), 3])
            push!(etpl10, reduce(vcat, (inter0[i], inter1[i, :])))
            push!(etpl20, reduce(vcat, (inter0[i], inter1[i, :])))
        elseif inter0[i] in inter_yp  # remain at all stages in y-p direction
            push!(etpl0, ["matzyp", reduce(vcat, (inter0[i], inter1[i, :])), 3])
            push!(etpl10, reduce(vcat, (inter0[i], inter1[i, :])))
            push!(etpl20, reduce(vcat, (inter0[i], inter1[i, :])))
        elseif inter0[i] in inter_xn  # remain at all stages in x-n direction
            push!(etpl0, ["matzxn", reduce(vcat, (inter0[i], inter1[i, :])), 3])
            push!(etpl10, reduce(vcat, (inter0[i], inter1[i, :])))
            push!(etpl20, reduce(vcat, (inter0[i], inter1[i, :])))
        elseif inter0[i] in inter_yn  # remain at all stages in y-n direction
            push!(etpl0, ["matzyn", reduce(vcat, (inter0[i], inter1[i, :])), 3])
            push!(etpl10, reduce(vcat, (inter0[i], inter1[i, :])))
            push!(etpl20, reduce(vcat, (inter0[i], inter1[i, :])))
        else
            println("The interface generation is not correct, exit.")
            exit()
        end
    end

    # return the nodes in the boundaries of first excavation with interface
    coord1into = Int[]
    nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(2, 17)  # interface nodes of first excavation in y-p direction
    nnodes = length(nodeTags)
    for i = 1:nnodes
        push!(coord1into, reinterpret(Int, nodeTags[i]))
    end
    nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(2, 18)  # interface nodes of first excavation in x-p direction
    nnodes = length(nodeTags)
    for i = 1:nnodes
        push!(coord1into, reinterpret(Int, nodeTags[i]))
    end
    coord1 = Int[]
    nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(2, 102)  # bottom soil nodes
    nnodes = length(nodeTags)
    for i = 1:nnodes
        push!(coord1, reinterpret(Int, nodeTags[i]))
    end
    for node in coord1into
        indx = findfirst(item -> item == node, coordi)
        push!(coord1, indx + maxNodeTag)  # transform original nodes to interface nodes on wedge
    end
    coord1 = unique(coord1, dims=1)
    sort!(coord1)

    # return the nodes after the first excavation
    coord10 = Int[]
    nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(3, 103)  # tet elements
    nnodes = length(nodeTags)
    for i = 1:nnodes
        push!(coord10, reinterpret(Int, nodeTags[i]))
    end
    nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(3, 104)  # wedge elements
    nnodes = length(nodeTags)
    for i = 1:nnodes
        node = reinterpret(Int, nodeTags[i])
        indx = findfirst(item -> item == node, coordi)
        if isnothing(indx)
            push!(coord10, node)
        else
            push!(coord10, indx + maxNodeTag)
        end
    end
    coord10 = unique(coord10, dims=1)
    sort!(coord10)

    # return the nodes in the boundaries of second excavation with interface
    coord2into = Int[]
    nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(2, 15)  # interface nodes of second excavation in y-p direction
    nnodes = length(nodeTags)
    for i = 1:nnodes
        push!(coord2into, reinterpret(Int, nodeTags[i]))
    end
    nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(2, 16)  # interface nodes of second excavation in x-p direction
    nnodes = length(nodeTags)
    for i = 1:nnodes
        push!(coord2into, reinterpret(Int, nodeTags[i]))
    end
    coord2 = Int[]  # boundaries of second excavation
    nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(2, 202)  # bottom soil nodes
    nnodes = length(nodeTags)
    for i = 1:nnodes
        push!(coord2, reinterpret(Int, nodeTags[i]))
    end
    for node in coord2into
        indx = findfirst(item -> item == node, coordi)
        push!(coord2, indx + maxNodeTag)  # transform original nodes to interface nodes on wedge
    end
    coord2 = unique(coord2, dims=1)
    sort!(coord2)

    # return the nodes after the second excavation
    coord20 = Int[]
    nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(3, 203)  # tet elements
    nnodes = length(nodeTags)
    for i = 1:nnodes
        push!(coord20, reinterpret(Int, nodeTags[i]))
    end
    nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(3, 204)  # wedge elements
    nnodes = length(nodeTags)
    for i = 1:nnodes
        node = reinterpret(Int, nodeTags[i])
        indx = findfirst(item -> item == node, coordi)
        if isnothing(indx)
            push!(coord20, node)
        else
            push!(coord20, indx + maxNodeTag)
        end
    end
    coord20 = unique(coord20, dims=1)
    sort!(coord20)

    gmsh.finalize()

    @show Snnodes, numISelem, maxNodeTag, maxElemTag

    return [coord0, etpl0, numISelem, etpl0map], [etpl1, coord1, etpl10, coord10, etpl1map], [etpl2, coord2, etpl20, coord20, etpl2map], Snnodes
end

"""
    readgmsh_excavnointer(argv)

Read the mesh data for the excavation test without interface to demonstrate the feasiblity of the code.

# Arguments:
- `argv`: the input mesh file with suffix, e.g., "BracedExcavation.msh"
"""
function readgmsh_excavnointer(argv)
    gmsh.initialize()
    gmsh.open(argv)

    # read all the nodes
    nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
    Tnnodes = length(nodeTags)
    coord0 = zeros(Float64, Tnnodes, 3)
    for i = 1:Tnnodes
        coord0[nodeTags[i], :] = nodeCoords[(i-1)*3+1:i*3]
    end

    # read all the elements with material properties
    etpl0 = Vector[]   # total 3D elements for the initial body force
    etpl0map = Int[]   # element tags of the total model

    etpl1 = Vector[]   # elements of the first excavation 
    etpl10 = Vector[]  # remaining elements after the first excavation
    etpl1map = Int[]  # remaining element tags after the first excavation

    etpl2 = Vector[]   # elements of the second excavation 
    etpl20 = Vector[]  # remaining elements after the second excavation
    etpl2map = Int[]  # remaining element tags after the second excavation

    entities = gmsh.model.getEntities(3)
    for e in entities
        physicalTags = gmsh.model.getPhysicalGroupsForEntity(e[1], e[2])
        for p in physicalTags
            name = gmsh.model.getPhysicalName(e[1], p)
            if occursin("tet", name)
                _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                numTelem, element = interpretElement(elemTags, elemNodeTags)
                element = view(element, :, [1, 2, 3, 4, 5, 6, 7, 8, 10, 9])
                for mat in ["LondonClay"]
                    if occursin(mat, name)
                        for i = 1:numTelem
                            push!(etpl0, [mat, element[i, :], 4])
                            push!(etpl0map, reinterpret(Int, elemTags[1][i]))
                        end
                    end
                end
                if occursin("excav01", name)  # excavated soils are modelled as tetrahedron
                    for i = 1:numTelem
                        push!(etpl1, element[i, :])
                    end
                end
                if occursin("rem01", name)
                    for i = 1:numTelem
                        push!(etpl10, element[i, :])
                        push!(etpl1map, reinterpret(Int, elemTags[1][i]))
                    end
                end
                if occursin("excav02", name)  # excavated soils are modelled as tetrahedron
                    for i = 1:numTelem
                        push!(etpl2, element[i, :])
                    end
                end
                if occursin("rem02", name)
                    for i = 1:numTelem
                        push!(etpl20, element[i, :])
                        push!(etpl2map, reinterpret(Int, elemTags[1][i]))
                    end
                end
            elseif occursin("wed", name)
                _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[1], e[2])
                numTelem, element = interpretElement(elemTags, elemNodeTags)
                element = view(element, :, [1, 2, 3, 4, 5, 6, 7, 10, 8, 13, 15, 14, 9, 11, 12])
                for mat in ["DiaphragmWall"]
                    if occursin(mat, name)
                        for i = 1:numTelem
                            push!(etpl0, [mat, element[i, :], 6])
                            push!(etpl0map, reinterpret(Int, elemTags[1][i]))
                        end
                    end
                end
                if occursin("rem01", name)
                    for i = 1:numTelem
                        push!(etpl10, element[i, :])
                        push!(etpl1map, reinterpret(Int, elemTags[1][i]))
                    end
                end
                if occursin("rem02", name)
                    for i = 1:numTelem
                        push!(etpl20, element[i, :])
                        push!(etpl2map, reinterpret(Int, elemTags[1][i]))
                    end
                end
            end
        end
    end

    # Nodes for the first excavation
    coord1 = Int[]  # boundary of first excavation
    nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(2, 105)  # read the coordinates
    nnodes = length(nodeTags)
    # println("Total nodes in the first excavation boundary is: $nnodes")
    for i = 1:nnodes
        push!(coord1, reinterpret(Int, nodeTags[i]))
    end

    # return the nodes after the first excavation
    coord10 = Int[]
    nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(3, 103)  # read the coordinates
    nnodes = length(nodeTags)
    for i = 1:nnodes
        push!(coord10, reinterpret(Int, nodeTags[i]))
    end
    nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(3, 104)  # read the coordinates
    nnodes = length(nodeTags)
    for i = 1:nnodes
        push!(coord10, reinterpret(Int, nodeTags[i]))
    end
    coord10 = unique(coord10, dims=1)
    sort!(coord10)

    # Nodes for the second excavation
    coord2 = Int[]  # boundary of second excavation
    nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(2, 205)  # read the coordinates
    nnodes = length(nodeTags)
    # println("Total nodes in the second excavation boundary is: $nnodes")
    for i = 1:nnodes
        push!(coord2, reinterpret(Int, nodeTags[i]))
    end

    # return the nodes after the second excavation
    coord20 = Int[]
    nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(3, 203)  # tet
    nnodes = length(nodeTags)
    for i = 1:nnodes
        push!(coord20, reinterpret(Int, nodeTags[i]))
    end
    nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(3, 204)  # wed 
    nnodes = length(nodeTags)
    for i = 1:nnodes
        push!(coord20, reinterpret(Int, nodeTags[i]))
    end
    coord20 = unique(coord20, dims=1)
    sort!(coord20)

    return [coord0, etpl0, etpl0map], [etpl1, coord1, etpl10, coord10, etpl1map], [etpl2, coord2, etpl20, coord20, etpl2map]
end
