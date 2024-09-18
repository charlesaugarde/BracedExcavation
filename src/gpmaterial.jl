"""
    gpmaterial(mesh, params)

Assign the initial mechanical properties to the Gauss points of elements. 

# Arguments:
- `mesh`: MESH structure defined in the [`BracedExcavation.MESH`](@ref).
- `params`: Input parameters defined by the [`BracedExcavation.PARAMSDT`](@ref). 
"""
function gpmaterial(mesh, params)
    gps = Vector[]
    nels = length(mesh.etpl)
    for nel in 1:nels
        elems = []
        for _ in 1:mesh.ngp[nel]
            elemd = GP()
            if mesh.mat[nel] == "LondonClay"
                elemd.E = params.E0
                elemd.v = params.v0
                elemd.G = params.E0 / 2 / (1 + params.v0)
                elemd.Kn = 0
                elemd.Ks = 0
            elseif mesh.mat[nel] == "DiaphragmWall"
                elemd.E = params.E1
                elemd.v = params.v1
                elemd.G = 0  # not consider the SSO model for diaphragm walls
                elemd.Kn = 0
                elemd.Ks = 0
            elseif occursin("matz", mesh.mat[nel])
                elemd.E = 0
                elemd.v = 0
                elemd.G = 0
                elemd.Kn = params.KN
                elemd.Ks = params.KS
            end
            elemd.sigP = zeros(Float64, 6)
            elemd.sigG = zeros(Float64, 6)
            elemd.sig = zeros(Float64, 6)
            elemd.eps = zeros(Float64, 6)
            elemd.epsP = zeros(Float64, 6)
            push!(elems, elemd)
        end
        push!(gps, elems)
    end

    return gps
end
