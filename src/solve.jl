using SparseArrays

"""
    solve(kval, krow, kcol, f, bc, nd::Int, NRit::Int)

Solve the algebra equations for the global FE analysis. 

# Arguments:
- `kval`: global stiffness values.
- `krow`: global stiffness row indices.
- `kcol`: global stiffness column indices.
- `f`: right-hand side for the equation ``Ku=f``.
- `bc`: boundary conditions.
- `nd`: model dimensions.
- `NRit`: iteration number.
"""
function solve(kval, krow, kcol, f, bc, nd::Int, NRit::Int)
    d = zeros(Float64, nd)
    fdof = zeros(Int64, nd)
    pdof = zeros(Int64, nd)
    K = sparse(krow, kcol, kval, nd, nd)
    fdof = collect(1:nd)
    pdof = Int.(bc[:, 1])
    deleteat!(fdof, pdof)
    np = size(bc, 1)
    dp = zeros(Float64, np)
    if NRit == 1
        dp = bc[:, 2]
    end
    s = K[fdof, fdof] \ (f[fdof] - K[fdof, pdof] * dp)
    d[pdof] = copy(dp)
    d[fdof] = copy(s)
    R = K * d - f

    return d, R
end
