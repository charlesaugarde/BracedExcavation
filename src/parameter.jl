"""
    PARAMSDT()

Declare the type of input parameters, and the default values.

# Arguments
- `E0`: Young's modulus for the London clay (*default: 1.0e6 N/m2*).
- `v0`: Poisson's ratio for the London clay (*default: 0.2*).
- `E1`: Young's modulus for the retaining wall (*default: 1.0e9 N/m2*).
- `v1`: Poisson's ratio for the retaining wall (*default: 0.3*).
- `E2`: Young's modulus for the prop (*default: 2.0e9 N/m2*).
- `v2`: Poisson's ratio for the prop (*default: 0.3*).
- `DENS`: density (*default: 2000.0 kg/m2*).
- `GRAV`: gravitational froce (*default: -9.81 N/kg*).
- `K0`: K0 value (*default: 1.0*).
- `W`: width of the model (*default: 1.0 m*).
- `L`: length of the model (*default: 1.0 m*).
- `H`: height of the model (*default: 1.0 m*).
- `P`: isotropic consolidation pressure (*default: -200e3 N/m2*).
- `U`: axial compressure displacement (*default: -0.4 m*).
- `LSPTC`: load steps for the consolidation (*default: 20*).
- `LSTPX`: load steps for the compression (*default: 20*).
- `CONSTM`: constitutive model (*default: "Hooked3d"*).
- `A`: parameters for the small-strain overlay model (*default: 0.385*).
- `GAM07`: strain at 70% initial shear stiffness (*default: 1e-1*).
- `THETA`: friction angle for the Drucker-Prager yield surface (*default:  ``\\pi/9``*).
- `C`: cohesion for the Drucker-Prager yield surface (*default: 0.0 N/m2*).
- `KN`: normal stiffness for the interface element (*default: 1.0e8 N/m2*).
- `KS`: shear stiffness for the interface element (*default: 1.0e6 N/m2*).
- `THETAI`: friction angle for the interface element (*default:  ``\\pi/9``*).
- `CI`: cohesion for the interface element (*default: 0.0 N/m2*).
- `TOL`: tolerance for the Newton-Ralphson iteration (*default: 1e-9*).
"""
@kwdef struct PARAMSDT
    E0::Float64 = 1.0e6          # Young's modulus for the London clay
    v0::Float64 = 0.2            # Poisson's ratio for the London clay
    E1::Float64 = 1.0e9          # Young's modulus for the retaining wall
    v1::Float64 = 0.3            # Poisson's ratio for the retaining wall
    E2::Float64 = 2.0e9          # Young's modulus for the prop
    v2::Float64 = 0.3            # Poisson's ratio for the prop
    DENS::Float64 = 2000.0       # density
    GRAV::Float64 = -9.81        # gravitational force
    K0::Float64 = 1.0            # K0 value
    W::Float64 = 1.0             # width of the model
    L::Float64 = 1.0             # length of the model
    H::Float64 = 1.0             # height of the model
    P::Float64 = -200e3          # isotropic consolidation pressure
    U::Float64 = -0.4            # axial compression displacement 
    LSTPC::Int8 = 20             # load steps for the consolidation
    LSTPX::Int8 = 20             # load steps for the compression
    CONSTM::String = "Hooke3d"   # constitutive model 
    A::Float64 = 0.385           # parameters for small-strain overlay model
    GAM07::Float64 = 1e-1        # strain at 70% initial shear stiffenss
    THETA::Float64 = π / 9       # friction angle for the Drucker-Prager yield surface
    C::Float64 = 0.0             # cohesion for the Drucker-Prager yield surface
    KN::Float64 = 1.0e8          # normal stiffness for the interface element
    KS::Float64 = 1.0e6          # shear stiffness for the interface element
    THETAI::Float64 = π / 9      # friction angle for the interface element
    CI::Float64 = 0.0            # cohesion for the interface element
    TOL::Float64 = 1e-9          # tolerance for the Newton-Ralphson iteration
end
