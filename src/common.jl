"""
    interpretElement(elemTags, elemNodeTags)

Reinterpret the Gmsh element data to the Int64 type. 
Return (a) the total number of elements in one physical group and (b) the nodes array of each element. 

# Arguments
- `elemTags`: Tags of element in the Gmsh format.
- `elemNodeTags`: Tags of node in the element in the Gmsh format.
"""
function interpretElement(elemTags, elemNodeTags)
    numTelem = 0
    numElem = sum(length(i) for i in elemTags)
    numTelem += numElem
    element = reinterpret(Int64, elemNodeTags[1])
    element = reshape(element, Int64(length(elemNodeTags[1]) / numTelem), numTelem)
    element = element'
    return numTelem, element
end


"""
    engtocauchy(strain::Array)

Transform the engineering strain tensor to Cauchy strain tensor.

# Arguments
- `strain`: input strain in 1D array (1 x 6).
"""
function engtocauchy(strain)
    @assert length(strain) == 6 "The size of strain is not correct."
    return [strain[1] strain[4]/2 strain[6]/2; strain[4]/2 strain[2] strain[5]/2; strain[6]/2 strain[5]/2 strain[3]]
end

"""
    cauchytoeng(strain::Array)

Transform the Cauchy strain tensor to the engineering strain tensor.

# Arguments
- `strain`: input strian in the tensor form (3 x 3).
"""
function cauchytoeng(strain)
    @assert size(strain) == (3, 3) "The size of strain tensor is not correct"
    return [strain[1, 1] strain[2, 2] strain[3, 3] strain[1, 2] + strain[2, 1] strain[1, 3] + strain[3, 1] strain[2, 3] + strain[3, 2]]
end


"""
    MESH()

Geometry properties of the mesh (**Global Properties**).

# Arguments:
- `lstp`: load step.
- `coord`: nodal coordinates.
- `etpl`: element topography.
- `mat`: material names.
- `ngp`: number of Gauss points.
- `bc`: boundary conditions.
- `krow`: global stiffness row indices.
- `kcol`: global stiffness column indices.
- `kval`: global stiffness values.
- `fint`: internal forces.
- `fb`: body forces.
- `fext`: external forces.
- `react`: reaction forces.
- `uvw`: mesh deformation.
"""
mutable struct MESH
    lstp::Int64                 # load step
    coord::Array{Float64,2}     # nodal coordinates
    etpl::Array{Vector{Int64}}  # element topography
    mat::Array{String}          # material names
    ngp::Array{Int8}            # number of Gauss points
    bc::Array{Float64,2}        # boundary conditions
    krow::Vector{Int64}         # global stiffness row indices
    kcol::Vector{Int64}         # global stiffness column indices
    kval::Vector{Float64}       # global stiffness values
    fint::Vector{Float64}       # internal forces
    fb::Vector{Float64}         # body forces
    fext::Vector{Float64}       # external forces
    react::Vector{Float64}      # reaction forces
    uvw::Vector{Float64}        # mesh deformation
    MESH() = new()
end

"""
    GP()

Mechanical behaviours at each Gauss point in an element (**Global Properties**).

# Arguments:
- `E`: Young's modulus.
- `G`: Initial shear modulus.
- `v`: Poisson's ratio.
- `Ks`: shear stiffness of interface element.
- `Kn`: normal stiffness of interface element.
- `H0`: deviatoric strain history of previous load step.
- `H`: deviatoric strian history of current load step.
- `Y`: yield flag.
- `sig`: stress of current load step.
- `sigP`: stress of previous load step.
- `sigG`: stress of geostatic state.
- `eps`: elastic strain of current load step.
- `epsP`: elastic strain of previous load step.
"""
mutable struct GP
    E::Float64                  # Young's modulus
    G::Float64                  # Initial shear modulus (constant, for SSO model)
    v::Float64                  # Poisson's ratio
    Ks::Float64                 # shear stiffness of interface element
    Kn::Float64                 # normal stiffness of interface element
    H0::Array{Float64}          # deviatoric strain history of previous load step
    H::Array{Float64}           # deviatoric strain history of current load step
    Y::Bool                     # yield flag
    sig::Vector{Float64}        # stress of current load step 
    sigP::Vector{Float64}       # stress of previous load step
    sigG::Vector{Float64}       # geostatic stress
    eps::Vector{Float64}        # elastic strain
    epsP::Vector{Float64}       # strain of previous load step
    GP() = new()
end
