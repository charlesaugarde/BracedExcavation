using LinearAlgebra, Infiltrator
"""
    elasticD(E, v)

Return the elastic stiffness matrix ``\\mathbf{D}``:

``
\\mathbf{D} = \\frac{E}{(1+\\nu)(1-2\\nu)} \\left[  \\begin{array}{cccccc} 1-\\nu & \\nu & \\nu & 0 & 0 & 0 \\\\ \\nu & 1-\\nu & \\nu & 0 & 0 & 0 \\\\ \\nu & \\nu & 1-\\nu &  0 & 0 & 0 \\\\  0 & 0 & 0 &  \\frac{1-2\\nu}{2} & 0 & 0 \\\\ 0 & 0 & 0 &  0 & \\frac{1-2\\nu}{2}  & 0 \\\\ 0 & 0 & 0 &  0 & 0 & \\frac{1-2\\nu}{2}   \\end{array} \\right].
``

# Arguments:
- `E`: Young's modulus.
- `v`: Poisson's ratio.
"""
function elasticD(E, v)
    D = zeros(Float64, 6, 6)
    D[1:3, 1:3] = E / (1 + v) / (1 - 2 * v) * (v * ones(Float64, 3, 3) + (1 - 2 * v) * diagm(ones(Float64, 3)))
    D[4:end, 4:end] = E / (1 + v) / 2 * diagm(ones(Float64, 3))

    return D
end

"""
    Hooke3d(deps, E, v)

Calculate the stress induced by the elastic strain. 

# Arguments:
- `deps`: elastic strain rate, or trial elastic strain rate for the return method.
- `E`: Young's modulus.
- `v`: Poisson's ratio.
"""
function Hooke3d(deps, E, v)
    D = elasticD(E, v)
    dsig = D * deps

    return D, dsig
end

"""
    smallstrain(deps, H0, G0, v, a, Gam07, direct)

Small strain overlay model ([Benz et al. 2009](https://onlinelibrary.wiley.com/doi/abs/10.1002/nag.701)). 
The shear modulus will be affected by the deviatoric strain history through the Hardin-Drnevich relationships. 

# Arguments:
- `deps`: elastic strain rate, or trial elastic strain rate for the return method.
- `H0`: deviatoric strain history.
- `G0`: initial shear modulus.
- `v`: Poisson's ratio.
- `a`: parameters for the small-strain overlay model.
- `Gam07`: strain at 70% initial shear stiffness.
- `direct`: direction of the loading: 1 for the monotonic loading, -1 for the unloading/reloading. 
"""
function smallstrain(deps, H0, G0, v, a, Gam07, direct)
    # assert the loading condition
    @assert direct in [1, -1] "Please assign the load direction to be 1 (compression) or -1 (extention)"
    H0 = engtocauchy(H0)  # transform the strain history vector to tensor
    te = engtocauchy(deps)  # transform the totoal strain vector to tensor
    de = te - tr(te) / 3 * I  # deviatoric strain
    if direct == 1
        eigvals, eigvecs = eigen(de, sortby=-)  # get the princinple direction of deviatoric strain (compression)
    elseif direct == -1
        eigvals, eigvecs = eigen(de, sortby=+)  # get the princinple direction of deviatoric strain (extension)
    end
    H1 = eigvecs' * H0 * eigvecs  # map the strain history to the princinple direction
    gam1 = sqrt(3) * norm(H1 * eigvals) / norm(eigvals)  # calculate the monotonic shear strain at previous load step
    T = diagm(ones(Float64, 3))  # transformation matrix for the loading direction
    for direct = 1:3
        if H1[direct, direct] * eigvals[direct] < 0.0
            T[direct, direct] = 1.0 / sqrt(H1[direct, direct] + 1.0)
        end
    end
    H2 = T * (H1 + I) * T' + diagm(eigvals) - I  # update the strain history tensor
    uH = cauchytoeng(H2)  # transform the updated strain history tensor to the vector for storage
    gam2 = sqrt(3) * norm(H2 * eigvals) / norm(eigvals)  # calculate the monotonic shear strain at current load step
    if abs(gam2 - gam1) < 1e-5  # Almost no change in the strain history (isotropic deformation)
        G = G0
    else
        G = G0 / (gam2 - gam1) * (gam2 / (1 + a * gam2 / Gam07) - gam1 / (1 + a * gam1 / Gam07))  # incremental form of stiffenss
    end
    E = 2 * G * (1 + v)
    D, dsig = Hooke3d(deps, E, v)  # update the stress strain relation based on linear elastic law

    return D, dsig, uH, E
end

"""
    yieldfunction(stress, frictionangle, cohesion)

Yield surface for the Drucker-Prage model. 

# Arguments:
- `stress`: stress at current load step.
- `frictionangle`: friction angle. 
- `cohesion`: cohesion.
"""
function yieldfunction(stress, frictionangle, cohesion)
    stress = -vec(stress)
    # Drucker-Prager yield function
    I1 = sum(stress[1:3])
    p = I1 / 3.0
    s = stress - p * [1, 1, 1, 0, 0, 0]
    J2 = 1 / 2 * s' * s
    eta = 6 * sin(frictionangle) / (sqrt(3) * (3 - sin(frictionangle)))
    xi = 6 * cos(frictionangle) / (sqrt(3) * (3 - sin(frictionangle)))
    r = sqrt(J2) - eta * abs(p) + xi * cohesion

    return r, s, J2, eta
end

"""
    DPconst(deps, E, v, frictionangle, cohesion, sigma, strain)

Elasto-perfect plastic Drucker-Prager model. 

# Arguments:
- `deps`: elastic strain rate, or trial elastic strain rate for the return method.
- `E`: Young's modulus.
- `v`: Poisson's ratio.
- `frictionangle`: friction angle controlling the yield surface.
- `cohesion`: cohession controlling the yield surface.
- `sigma`: stress of previous load step. 
- `strain`: strain of previous load step.
"""
function DPconst(deps, E, v, frictionangle, cohesion, sigma, strain)
    # parameters for the Drucker-Prager model
    # frictionangle = pi / 20  # For DPSSO model
    # frictionangle = pi / 9  # For demo of DPconst model
    # dilationangle = pi / 36
    dilationangle = frictionangle
    # cohesion = 0
    tol = 1e-9  # tolerance
    # elastic trial
    De = elasticD(E, v)
    dsigma = De * deps
    sigE = sigma + dsigma
    ftrial, strial, J2trial, eta = yieldfunction(sigE, frictionangle, cohesion)
    if ftrial > tol # plastic behaviour
        # @show "comes to here"
        # @show sigE
        # @show sigma
        # @infiltrate
        flag = 1
        etabar = 6 * sin(dilationangle) / (sqrt(3) * (3 - sin(dilationangle)))
        dpsi = strial / (2 * sqrt(J2trial)) + etabar / 3 * [1, 1, 1, 0, 0, 0]
        dphi = strial / (2 * sqrt(J2trial)) + eta / 3 * [1, 1, 1, 0, 0, 0]
        FACT1 = De' * dpsi * dphi' * De
        FACT2 = dphi' * De * dpsi
        Dalg = De - FACT1 / FACT2
        dsigma = Dalg * deps
        sigE = sigma + dsigma
        λ = dphi' * De * deps / FACT2
        dεp = λ .* dpsi
        epsE = strain + deps - dεp
    else # elastic behaviour
        flag = 0
        epsE = strain + deps
        Dalg = De
    end

    return Dalg, dsigma, epsE, flag
end
