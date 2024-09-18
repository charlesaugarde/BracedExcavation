using DelimitedFiles

"""
    postpro(uvw, lstp, coord, etpl, gps, outdir)

Generate the mesh vtk files, deformed mesh vtk files, and the properties in each Gauss point.

# Arguments:
- `uvw`: displacement field in the MESH struct.
- `lstp`: load step in the MESH struct.
- `coord`: coordinates of nodes in the MESH struct.
- `etpl`: element topography in the MESH struct.
- `gps`: GP properties of elements.
- `outdir`: output directory to save the *.vtk* files.
"""
function postpro(uvw, lstp, coord, etpl, gps, outdir)
    meshName = outdir * "/mesh_$lstp.vtk"
    makeVtk(coord, etpl, uvw, meshName)
    nodes, nD = size(coord)
    coordDef = coord + reshape(uvw, nD, nodes)'
    meshName = outdir * "/meshDef_$lstp.vtk"
    makeVtk(coordDef, etpl, uvw, meshName)

    gpC = Vector{Float64}[]
    gpDat = Vector{Float64}[]
    nels = length(etpl)
    for nel in 1:nels
        nen = length(etpl[nel])
        if nen == 12
            ngp, N, _, _ = shapefunc(2, nen)
        elseif nen == 2
            ngp, N, _, _ = shapefunc(1, nen)
        else
            ngp, N, _, _ = shapefunc(nD, nen)
        end
        for gp in 1:ngp
            push!(gpDat, [gps[nel][gp].sig; gps[nel][gp].eps])
            if nen == 12
                push!(gpC, round.((N[1:6, gp]' * coordDef[etpl[nel][1:6], :])', digits=6))
            else
                push!(gpC, round.((N[:, gp]' * coordDef[etpl[nel], :])', digits=6))
            end
        end
    end
    gpDatFields = ["sig_xx", "sig_yy", "sig_zz", "sig_xy", "sig_yz", "sig_zx",
        "epsE_xx", "epsE_yy", "epsE_zz", "epsE_xy", "epsE_yz",
        "epsE_zx"]                                                  # Gauss point data field description

    gpDataName = outdir * "/gpData_$lstp.vtk"
    makeVtkGP(gpDataName, gpC, gpDat, gpDatFields)

    return nothing
end

"""
    makeVtkGP(mpFileName, gpC, gpDat, datNames)

Generate the GP field data in each Gauss point.

# Arguments:
- `mpFileName`: output filename.
- `gpC`: total number of Gauss points.
- `gpDat`: stress and strain properties in each Gauss point.
- `datNames`: names of the GP properties.
"""
function makeVtkGP(mpFileName, gpC, gpDat, datNames)
    ngp = size(gpC, 1)
    nD = length(gpC[1])
    open(mpFileName, "w") do f
        write(f, "# vtk DataFile Version 2.0\n")
        write(f, "Julia generated vtk file\n")
        write(f, "ASCII\n")
        write(f, "DATASET UNSTRUCTURED_GRID\n")
        write(f, "POINTS $ngp double\n")

        # position output 
        if nD < 3
            gpC = [gpC zeros(ngp, 3 - nD)]
        end
        writedlm(f, gpC, " ")
        write(f, "\n")
        write(f, "POINT_DATA $ngp \n")
        nDat = length(gpDat[1])
        for i = 1:nDat
            write(f, "SCALARS ", string(datNames[i]), " FLOAT 1\n")
            write(f, "LOOKUP_TABLE default\n")
            writedlm(f, getindex.(gpDat, i), " ")
            write(f, "\n")
        end
    end

    return nothing
end

"""
    makeVtk(coord, etpl, uvw, meshName)

Generate the mesh data *.vtk* files.

# Arguments:
- `coord`: coordinates of nodes in MESH struct.
- `etpl`: element topography in MESH struct.
- `uvw`: displacement field in the MESH struct.
- `meshName`: output file name. 
"""
function makeVtk(coord, etpl, uvw, meshName)
    nodes, nD = size(coord)
    nels = length(etpl)
    coord = round.(coord, digits=6)
    uvw = round.(uvw, digits=6)

    # Generation of vtk file
    open(meshName, "w") do f
        write(f, "# vtk DataFile Version 2.0\n")
        write(f, "Julia generated vtk file\n")
        write(f, "ASCII\n")
        write(f, "DATASET UNSTRUCTURED_GRID\n")
        write(f, "POINTS $nodes double\n")

        # nodal coordinates
        if nD < 3
            coord = [coord zeros(nodes, 3 - nD)]
        end
        writedlm(f, coord, " ")
        write(f, "\n")

        # element topology
        nn = sum(length, etpl) + nels  # total number of elements in the etpl
        @debug "The total No. of elements in the etpl is" nn
        write(f, "CELLS $nels $nn \n")
        etplOutput = []
        eltype = []
        for nel in etpl
            nen = length(nel)
            tvtk, elemId = etpl_vtk(nD, nen)
            push!(eltype, elemId)
            elem = view(nel .- 1, tvtk)
            push!(etplOutput, [nen elem])
        end
        writedlm(f, etplOutput, " ")
        write(f, "\n")

        # element types
        write(f, "CELL_TYPES $nels \n")
        writedlm(f, eltype, " ")
        write(f, "\n")

        # displacement output
        write(f, "POINT_DATA $nodes\n")
        if nD == 3
            write(f, "SCALARS u_x FLOAT 1\n")
            write(f, "LOOKUP_TABLE default\n")
            writedlm(f, uvw[1:nD:end], " ")
            write(f, "\n")

            write(f, "SCALARS u_y FLOAT 1\n")
            write(f, "LOOKUP_TABLE default\n")
            writedlm(f, uvw[2:nD:end], " ")
            write(f, "\n")

            write(f, "SCALARS u_z FLOAT 1\n")
            write(f, "LOOKUP_TABLE default\n")
            writedlm(f, uvw[3:nD:end], " ")
            write(f, "\n")
        elseif nD == 2
            write(f, "SCALARS u_x FLOAT 1\n")
            write(f, "LOOKUP_TABLE default\n")
            writedlm(f, uvw[1:nD:end], " ")
            write(f, "\n")

            write(f, "SCALARS u_y FLOAT 1\n")
            write(f, "LOOKUP_TABLE default\n")
            writedlm(f, uvw[2:nD:end], " ")
            write(f, "\n")
        elseif nD == 1
            write(f, "SCALARS u_x FLOAT 1\n")
            write(f, "LOOKUP_TABLE default\n")
            writedlm(f, uvw, " ")
            write(f, "\n")
        end
    end

    return nothing
end

"""
    etpl_vtk(nD, nen)

Transform the nodal numbering in the FE etpl to paraview vtk format. 

# Arguments:
- `nD`: dimension.
- `nen`: number of nodes in the element. 
"""
function etpl_vtk(nD, nen)
    if nD == 3 && nen == 8
        tvtk = [1 2 3 4 5 6 7 8]
        elemId = 12
    elseif nD == 3 && nen == 10
        tvtk = [1 2 3 4 5 6 7 8 9 10]
        elemId = 24
    elseif nD == 3 && nen == 15
        tvtk = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]
        elemId = 26
    elseif nen == 12
        tvtk = [1 2 3 7 8 9 4 5 6 10 11 12]
        elemId = 31
    elseif nen == 2
        tvtk = [1 2]
        elemId = 3
    end

    return tvtk, elemId
end
