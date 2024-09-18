"""
    shapefunc(nD, nen)

Return (a) the total number of Guass points in the elements; (b) the shape function at each Gauss point; (c) the derivation of shape function at each Gauss point; and (d) the weights of Gauss points. 

The shape functions of current code contain: 

### 1D element

- 2-node line element with 1 Gauss point

### 2D element

- 4-node rectangle element with 1 Gauss point
- 6-node triangle element with 3 Gauss points
- 8-node rectangle element with 4 Gauss points
- 12-node interface element with 6 Gauss points

### 3D element

- 8-node hexahedron element with 8 Gauss points
- 10-node tetrahedron element with 4 Gauss points
- 15-node wedge element with 6 Gauss points

# Arguments:
- `nD`: dimension of the element
- `nen`: number of nodes in each element
"""
function shapefunc(nD, nen)
    if nD == 1
        if nen == 2  # 2-node line
            ngp = 1
            N = zeros(Float64, nen, ngp)
            dNr = zeros(Float64, ngp * nD, nen)

            # coordinate of Gauss points
            xsi = [0]

            # weight 
            wp = [2]

            # shape function
            N[1, :] = 1 / 2 * (1 .- xsi)
            N[2, :] = 1 / 2 * (1 .+ xsi)

            # derivation of shape function
            r2 = ngp * nD
            dNr[1:2:r2, 1] .= -1 / 2
            dNr[2:2:r2, 1] .= 1 / 2
        end
    elseif nD == 2
        if nen == 4  # 4-node rectangle
            ngp = 1
            N = zeros(Float64, nen, ngp)
            dNr = zeros(Float64, ngp * nD, nen)

            # coordinate of Gauss points
            xsi = [0]
            eta = [0]

            # weight 
            wp = [4]

            # shape function
            N[1, :] = 1 / 4 * (1 .- xsi) .* (1 .- eta)
            N[2, :] = 1 / 4 * (1 .- xsi) .* (1 .+ eta)
            N[3, :] = 1 / 4 * (1 .+ xsi) .* (1 .+ eta)
            N[4, :] = 1 / 4 * (1 .+ xsi) .* (1 .- eta)

            # derivation of shape function
            r2 = ngp * nD
            dNr[1:2:r2, 1] = 1 / 4 * (eta .- 1)
            dNr[1:2:r2, 2] = -1 / 4 * (eta .+ 1)
            dNr[1:2:r2, 3] = 1 / 4 * (eta .+ 1)
            dNr[1:2:r2, 4] = -1 / 4 * (eta .- 1)

            dNr[2:2:r2, 1] = 1 / 4 * (xsi .- 1)
            dNr[2:2:r2, 2] = -1 / 4 * (xsi .- 1)
            dNr[2:2:r2, 3] = 1 / 4 * (xsi .+ 1)
            dNr[2:2:r2, 4] = -1 / 4 * (xsi .+ 1)
        elseif nen == 6  # 6-node triangle
            ngp = 3
            N = zeros(Float64, nen, ngp)
            dNr = zeros(Float64, ngp * nD, nen)

            # coordinate of Gauss points
            g2 = 1 / 6
            xsi = [g2, 4 * g2, g2]
            eta = [g2, g2, 4 * g2]

            # weight 
            wp = ones(ngp, 1) / 6

            # shape function
            N[1, :] = (1 .- 2 * xsi - 2 * eta) .* (1 .- xsi - eta)
            N[2, :] = (2 * xsi .- 1) .* xsi
            N[3, :] = (2 * eta .- 1) .* eta
            N[4, :] = 4 * (1 .- xsi - eta) .* xsi
            N[5, :] = 4 * xsi .* eta
            N[6, :] = 4 * (1 .- xsi - eta) .* eta

            # derivation of shape function
            r2 = ngp * nD
            dNr[1:2:r2, 1] = -4 * (1 .- xsi - eta) .+ 1
            dNr[1:2:r2, 2] = 4 .* xsi .- 1
            # dNr[1:2:r2, 3] = 0
            dNr[1:2:r2, 4] = 4 * (1 .- 2 * xsi - eta)
            dNr[1:2:r2, 5] = 4 * eta
            dNr[1:2:r2, 6] = -4 * eta

            dNr[2:2:r2, 1] = -4 * (1 .- xsi - eta) .+ 1
            # dNr[2:2:r2, 2] = 0
            dNr[2:2:r2, 3] = 4 * eta .- 1
            dNr[2:2:r2, 4] = -4 * xsi
            dNr[2:2:r2, 5] = 4 * xsi
            dNr[2:2:r2, 6] = 4 * (1 .- xsi - 2 * eta)
        elseif nen == 8  # 8-node rectangle
            ngp = 4
            N = zeros(Float64, nen, ngp)
            dNr = zeros(Float64, ngp * nD, nen)

            # coordinate of Gauss points
            g2 = 1 / sqrt(3)
            xsi = [-1 -1 1 1]' * g2
            eta = [-1 1 -1 1]' * g2

            # weight 
            wp = ones(ngp, 1)

            # shape function
            N[1, :] = 1 / 4 * (1 .- xsi) .* (1 .- eta) .* (-xsi - eta .- 1)
            N[2, :] = 1 / 4 * (1 .+ xsi) .* (1 .- eta) .* (xsi - eta .- 1)
            N[3, :] = 1 / 4 * (1 .+ xsi) .* (1 .+ eta) .* (xsi + eta .- 1)
            N[4, :] = 1 / 4 * (1 .- xsi) .* (1 .+ eta) .* (-xsi + eta .- 1)
            N[5, :] = 1 / 2 * (1 .- xsi .^ 2) .* (1 .- eta)
            N[6, :] = 1 / 2 * (1 .+ xsi) .* (1 .- eta .^ 2)
            N[7, :] = 1 / 2 * (1 .- xsi .^ 2) .* (1 .+ eta)
            N[8, :] = 1 / 2 * (1 .- xsi) .* (1 .- eta .^ 2)

            # derivation of shape function
            r2 = ngp * nD
            dNr[1:2:r2, 1] = 1 / 4 * (eta .- 1) .* (-2 * xsi - eta)
            dNr[1:2:r2, 2] = 1 / 4 * (eta .- 1) .* (-2 * xsi + eta)
            dNr[1:2:r2, 3] = 1 / 4 * (eta .+ 1) .* (2 * xsi + eta)
            dNr[1:2:r2, 4] = 1 / 4 * (eta .+ 1) .* (2 * xsi - eta)
            dNr[1:2:r2, 5] = (eta .- 1) .* xsi
            dNr[1:2:r2, 6] = 1 / 2 * (1 .- eta .^ 2)
            dNr[1:2:r2, 7] = -(eta .+ 1) .* xsi
            dNr[1:2:r2, 8] = -1 / 2 * (1 .- eta .^ 2)

            dNr[2:2:r2, 1] = 1 / 4 * (xsi .- 1) .* (-2 * eta - xsi)
            dNr[2:2:r2, 2] = 1 / 4 * (xsi .+ 1) .* (2 * eta - xsi)
            dNr[2:2:r2, 3] = 1 / 4 * (xsi .+ 1) .* (2 * eta + xsi)
            dNr[2:2:r2, 4] = 1 / 4 * (xsi .- 1) .* (-2 * eta + xsi)
            dNr[2:2:r2, 5] = -1 / 2 * (1 .- xsi .^ 2)
            dNr[2:2:r2, 6] = -(xsi .+ 1) .* eta
            dNr[2:2:r2, 7] = 1 / 2 * (1 .- xsi .^ 2)
            dNr[2:2:r2, 8] = (xsi .- 1) .* eta
        elseif nen == 12
            ngp = 3
            N = zeros(Float64, nen, ngp)
            dNr = zeros(Float64, ngp * nD, nen)

            # coordinate of Gauss points
            g2 = 1 / 6
            xsi = [g2, 4 * g2, g2]
            eta = [g2, g2, 4 * g2]

            # weight 
            wp = ones(ngp, 1) / 6

            # shape function
            N[1, :] = (1 .- 2 * xsi - 2 * eta) .* (1 .- xsi - eta)
            N[2, :] = (2 * xsi .- 1) .* xsi
            N[3, :] = (2 * eta .- 1) .* eta
            N[4, :] = 4 * (1 .- xsi - eta) .* xsi
            N[5, :] = 4 * xsi .* eta
            N[6, :] = 4 * (1 .- xsi - eta) .* eta
            N[7:end, :] = copy(N[1:6, :])  # assign the same shape functions

            # derivation of shape function
            r2 = ngp * nD
            dNr[1:2:r2, 1] = -4 * (1 .- xsi - eta) .+ 1
            dNr[1:2:r2, 2] = 4 .* xsi .- 1
            # dNr[1:2:r2, 3] = 0
            dNr[1:2:r2, 4] = 4 * (1 .- 2 * xsi - eta)
            dNr[1:2:r2, 5] = 4 * eta
            dNr[1:2:r2, 6] = -4 * eta

            dNr[2:2:r2, 1] = -4 * (1 .- xsi - eta) .+ 1
            # dNr[2:2:r2, 2] = 0
            dNr[2:2:r2, 3] = 4 * eta .- 1
            dNr[2:2:r2, 4] = -4 * xsi
            dNr[2:2:r2, 5] = 4 * xsi
            dNr[2:2:r2, 6] = 4 * (1 .- xsi - 2 * eta)
            dNr[:, 7:end] = copy(dNr[:, 1:6])
        end
    elseif nD == 3
        if nen == 8  # 8-node hex
            ngp = 8  # No. of Gauss points
            N = zeros(Float64, nen, ngp)
            dNr = zeros(Float64, ngp * nD, nen)

            # coordinate of Gauss points
            xsi = zeros(Float64, ngp)
            eta = zeros(Float64, ngp)
            zet = zeros(Float64, ngp)
            g2 = 1 / sqrt(3)
            xsi = [-1, -1, 1, 1, -1, -1, 1, 1] * g2
            eta = [-1, -1, -1, -1, 1, 1, 1, 1] * g2
            zet = [-1, 1, 1, -1, -1, 1, 1, -1] * g2

            # weight
            wp = ones(ngp, 1)

            # shape function
            N[1, :] = (1 .- xsi) .* (1 .- eta) .* (1 .- zet) / 8
            N[2, :] = (1 .- xsi) .* (1 .- eta) .* (1 .+ zet) / 8
            N[3, :] = (1 .+ xsi) .* (1 .- eta) .* (1 .+ zet) / 8
            N[4, :] = (1 .+ xsi) .* (1 .- eta) .* (1 .- zet) / 8
            N[5, :] = (1 .- xsi) .* (1 .+ eta) .* (1 .- zet) / 8
            N[6, :] = (1 .- xsi) .* (1 .+ eta) .* (1 .+ zet) / 8
            N[7, :] = (1 .+ xsi) .* (1 .+ eta) .* (1 .+ zet) / 8
            N[8, :] = (1 .+ xsi) .* (1 .+ eta) .* (1 .- zet) / 8

            # derivation of shape function
            r2 = ngp * nD
            dNr[1:3:r2, 1] = -1 / 8 * (1.0 .- eta) .* (1.0 .- zet)
            dNr[1:3:r2, 2] = -1 / 8 * (1.0 .- eta) .* (1.0 .+ zet)
            dNr[1:3:r2, 3] = 1 / 8 * (1.0 .- eta) .* (1.0 .+ zet)
            dNr[1:3:r2, 4] = 1 / 8 * (1.0 .- eta) .* (1.0 .- zet)
            dNr[1:3:r2, 5] = -1 / 8 * (1.0 .+ eta) .* (1.0 .- zet)
            dNr[1:3:r2, 6] = -1 / 8 * (1.0 .+ eta) .* (1.0 .+ zet)
            dNr[1:3:r2, 7] = 1 / 8 * (1.0 .+ eta) .* (1.0 .+ zet)
            dNr[1:3:r2, 8] = 1 / 8 * (1.0 .+ eta) .* (1.0 .- zet)

            dNr[2:3:r2, 1] = -1 / 8 * (1.0 .- xsi) .* (1.0 .- zet)
            dNr[2:3:r2, 2] = -1 / 8 * (1.0 .- xsi) .* (1.0 .+ zet)
            dNr[2:3:r2, 3] = -1 / 8 * (1.0 .+ xsi) .* (1.0 .+ zet)
            dNr[2:3:r2, 4] = -1 / 8 * (1.0 .+ xsi) .* (1.0 .- zet)
            dNr[2:3:r2, 5] = 1 / 8 * (1.0 .- xsi) .* (1.0 .- zet)
            dNr[2:3:r2, 6] = 1 / 8 * (1.0 .- xsi) .* (1.0 .+ zet)
            dNr[2:3:r2, 7] = 1 / 8 * (1.0 .+ xsi) .* (1.0 .+ zet)
            dNr[2:3:r2, 8] = 1 / 8 * (1.0 .+ xsi) .* (1.0 .- zet)

            dNr[3:3:r2, 1] = -1 / 8 * (1.0 .- xsi) .* (1.0 .- eta)
            dNr[3:3:r2, 2] = 1 / 8 * (1.0 .- xsi) .* (1.0 .- eta)
            dNr[3:3:r2, 3] = 1 / 8 * (1.0 .+ xsi) .* (1.0 .- eta)
            dNr[3:3:r2, 4] = -1 / 8 * (1.0 .+ xsi) .* (1.0 .- eta)
            dNr[3:3:r2, 5] = -1 / 8 * (1.0 .- xsi) .* (1.0 .+ eta)
            dNr[3:3:r2, 6] = 1 / 8 * (1.0 .- xsi) .* (1.0 .+ eta)
            dNr[3:3:r2, 7] = 1 / 8 * (1.0 .+ xsi) .* (1.0 .+ eta)
            dNr[3:3:r2, 8] = -1 / 8 * (1.0 .+ xsi) .* (1.0 .+ eta)
        elseif nen == 10
            ngp = 4
            N = zeros(Float64, nen, ngp)
            dNr = zeros(Float64, ngp * nD, nen)

            # coordinate of Gauss points
            xsi = zeros(Float64, ngp)
            eta = zeros(Float64, ngp)
            zet = zeros(Float64, ngp)
            g1 = (5 - sqrt(5)) / 20
            g2 = (5 + 3 * sqrt(5)) / 20
            xsi = [g1 g2 g1 g1]'
            eta = [g1 g1 g2 g1]'
            zet = [g1 g1 g1 g2]'

            # weight
            wp = ones(ngp, 1) / 24

            # shape function
            N[1, :] = (2 * (1 .- xsi - eta - zet) .- 1) .* (1 .- xsi - eta - zet)
            N[2, :] = (2 * xsi .- 1) .* xsi
            N[3, :] = (2 * eta .- 1) .* eta
            N[4, :] = (2 * zet .- 1) .* zet
            N[5, :] = 4 * (1 .- xsi - eta - zet) .* xsi
            N[6, :] = 4 * xsi .* eta
            N[7, :] = 4 * (1 .- xsi - eta - zet) .* eta
            N[8, :] = 4 * (1 .- xsi - eta - zet) .* zet
            N[9, :] = 4 * xsi .* zet
            N[10, :] = 4 * eta .* zet

            # derivation of shape function
            r2 = ngp * nD
            dNr[1:3:r2, 1] = -4 * (1 .- xsi - eta - zet) .+ 1
            dNr[1:3:r2, 2] = 4 * xsi .- 1
            #dNr[1:3:r2, 3]= 0
            #dNr[1:3:r2, 4]= 0
            dNr[1:3:r2, 5] = 4 * (1 .- 2 * xsi - eta - zet)
            dNr[1:3:r2, 6] = 4 * eta
            dNr[1:3:r2, 7] = -4 * eta
            dNr[1:3:r2, 8] = -4 * zet
            dNr[1:3:r2, 9] = 4 * zet
            #dNr[1:3:r2,10]= 0

            dNr[2:3:r2+1, 1] = -4 * (1 .- xsi - eta .- zet) .+ 1
            #dNr[2:3:r2+1, 2]= 0
            dNr[2:3:r2+1, 3] = 4 * eta .- 1
            #dNr[2:3:r2+1, 4]= 0
            dNr[2:3:r2+1, 5] = -4 * xsi
            dNr[2:3:r2+1, 6] = 4 * xsi
            dNr[2:3:r2+1, 7] = 4 * (1 .- xsi - 2 * eta - zet)
            dNr[2:3:r2+1, 8] = -4 * zet
            #dNr[2:3:r2+1, 9]= 0
            dNr[2:3:r2+1, 10] = 4 * zet

            dNr[3:3:r2+2, 1] = -4 * (1 .- xsi - eta .- zet) .+ 1
            #dNr[3:3:r2+2, 2]= 0 
            #dNr[3:3:r2+2, 3]= 0
            dNr[3:3:r2+2, 4] = 4 * zet .- 1
            dNr[3:3:r2+2, 5] = -4 * xsi
            #dNr[3:3:r2+2, 6]= 0
            dNr[3:3:r2+2, 7] = -4 * eta
            dNr[3:3:r2+2, 8] = 4 * (1 .- xsi - eta - 2 .* zet)
            dNr[3:3:r2+2, 9] = 4 * xsi
            dNr[3:3:r2+2, 10] = 4 * eta
        elseif nen == 15
            ngp = 6
            N = zeros(Float64, nen, ngp)
            dNr = zeros(Float64, ngp * nD, nen)

            # coordinate of Gauss points
            xsi = zeros(Float64, ngp)
            eta = zeros(Float64, ngp)
            zet = zeros(Float64, ngp)
            g1 = 1 / 6
            g2 = 1 / sqrt(3)
            xsi = [g1; 4 * g1; g1; g1; 4 * g1; g1]
            eta = [g1; g1; 4 * g1; g1; g1; 4 * g1]
            zet = [-g2; -g2; -g2; g2; g2; g2]

            # weight
            wp = ones(ngp, 1) / 6

            # shape function
            f1 = 1 .- xsi - eta
            f2 = (1 .- zet .* zet)
            N[1, :] = 0.5 * (f1 .* (2 .* f1 .- 1) .* (1 .- zet) .- f1 .* f2)
            N[2, :] = 0.5 * (xsi .* (2 * xsi .- 1) .* (1 .- zet) .- xsi .* f2)
            N[3, :] = 0.5 * (eta .* (2 * eta .- 1) .* (1 .- zet) .- eta .* f2)
            N[4, :] = 0.5 * (f1 .* (2 * f1 .- 1) .* (1 .+ zet) .- f1 .* f2)
            N[5, :] = 0.5 * (xsi .* (2 * xsi .- 1) .* (1 .+ zet) .- xsi .* f2)
            N[6, :] = 0.5 * (eta .* (2 * eta .- 1) .* (1 .+ zet) .- eta .* f2)
            N[7, :] = 2.0 * (f1 .* xsi .* (1 .- zet))
            N[8, :] = 2.0 * xsi .* eta .* (1 .- zet)
            N[9, :] = 2.0 * (f1 .* eta .* (1 .- zet))
            N[10, :] = 2.0 * (f1 .* xsi .* (1 .+ zet))
            N[11, :] = 2.0 * xsi .* eta .* (1 .+ zet)
            N[12, :] = 2.0 * (f1 .* eta .* (1 .+ zet))
            N[13, :] = f1 .* f2
            N[14, :] = xsi .* f2
            N[15, :] = eta .* f2

            # derivation of shape function
            r2 = ngp * nD
            dNr[1:3:r2, 1] = (1 .- zet) .* (4 * xsi + 4 * eta + zet .- 2) / 2
            dNr[1:3:r2, 2] = (1 .- zet) .* (4 * xsi .- zet .- 2) / 2
            # dNr[1:3:r2, 3] = 0
            dNr[1:3:r2, 4] = (1 .+ zet) .* (4 * xsi + 4 * eta - zet .- 2) / 2
            dNr[1:3:r2, 5] = (1 .+ zet) .* (4 * xsi + zet .- 2) / 2
            # dNr[1:3:r2, 6] = 0
            dNr[1:3:r2, 7] = 2 * (1 .- zet) .* (1 .- 2 * xsi - eta)
            dNr[1:3:r2, 8] = 2 * eta .* (1 .- zet)
            dNr[1:3:r2, 9] = -2 * eta .* (1 .- zet)
            dNr[1:3:r2, 10] = 2 * (1 .+ zet) .* (1 .- 2 * xsi - eta)
            dNr[1:3:r2, 11] = 2 * eta .* (1 .+ zet)
            dNr[1:3:r2, 12] = -2 * eta .* (1 .+ zet)
            dNr[1:3:r2, 13] = -(1 .- zet .^ 2)
            dNr[1:3:r2, 14] = (1 .- zet .^ 2)
            # dNr[1:3:r2, 15] = 0

            dNr[2:3:r2+1, 1] = (1 .- zet) .* (4 * xsi + 4 * eta + zet .- 2) / 2
            # dNr[2:3:r2+1, 2] = 0
            dNr[2:3:r2+1, 3] = (1 .- zet) .* (4 * eta - zet .- 2) / 2
            dNr[2:3:r2+1, 4] = (1 .+ zet) .* (4 * xsi + 4 * eta - zet .- 2) / 2
            # dNr[2:3:r2+1, 5] = 0
            dNr[2:3:r2+1, 6] = (1 .+ zet) .* (4 * eta + zet .- 2) / 2
            dNr[2:3:r2+1, 7] = -2 * xsi .* (1 .- zet)
            dNr[2:3:r2+1, 8] = 2 * xsi .* (1 .- zet)
            dNr[2:3:r2+1, 9] = 2 * (1 .- zet) .* (1 .- xsi - 2 * eta)
            dNr[2:3:r2+1, 10] = -2 * xsi .* (1 .+ zet)
            dNr[2:3:r2+1, 11] = 2 * xsi .* (1 .+ zet)
            dNr[2:3:r2+1, 12] = 2 * (1 .+ zet) .* (1 .- xsi - 2 * eta)
            dNr[2:3:r2+1, 13] = -(1 .- zet .^ 2)
            # dNr[2:3:r2+1, 14] = 0
            dNr[2:3:r2+1, 15] = (1 .- zet .^ 2)

            dNr[3:3:r2+2, 1] = (1 .- xsi - eta) .* (2 * xsi + 2 * eta + 2 * zet .- 1) / 2
            dNr[3:3:r2+2, 2] = -xsi .* (2 * xsi - 2 * zet .- 1) / 2
            dNr[3:3:r2+2, 3] = -eta .* (2 * eta - 2 * zet .- 1) / 2
            dNr[3:3:r2+2, 4] = -(1 .- xsi - eta) .* (2 * xsi + 2 * eta - 2 * zet .- 1) / 2
            dNr[3:3:r2+2, 5] = xsi .* (2 * xsi + 2 * zet .- 1) / 2
            dNr[3:3:r2+2, 6] = eta .* (2 * eta + 2 * zet .- 1) / 2
            dNr[3:3:r2+2, 7] = -2 * xsi .* (1 .- xsi - eta)
            dNr[3:3:r2+2, 8] = -2 * xsi .* eta
            dNr[3:3:r2+2, 9] = -2 * eta .* (1 .- xsi - eta)
            dNr[3:3:r2+2, 10] = 2 * xsi .* (1 .- xsi - eta)
            dNr[3:3:r2+2, 11] = 2 * xsi .* eta
            dNr[3:3:r2+2, 12] = 2 * eta .* (1 .- xsi - eta)
            dNr[3:3:r2+2, 13] = -2 * zet .* (1 .- xsi - eta)
            dNr[3:3:r2+2, 14] = -2 * zet .* xsi
            dNr[3:3:r2+2, 15] = -2 * zet .* eta
        end
    else
        println("ONLY the 1, 2, and 3D functions are included here.")
    end
    return ngp, N, dNr, wp
end
