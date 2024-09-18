using Gmsh: gmsh

"""
    unithex(project_name::String, nopopup::Bool)
  
Generate the unit hexahedron element mesh by Gmsh with the output: *project_name.msh*. 
The size of element is 1 x 1 x 1. 
The name of the materials is *hex_LondonClay*, the former means the shape of the element, while the latter means the material property. 
Three boundaries are defined as the names of *top*, *sidex*, and *sidey*. 
# Arguments:
- `project_name`: name of the project, corresponding to the name of the mesh.
- `nopopup`: not show the gmsh window to visualise the mesh (*default is true*).
"""
function unithex(project_name::String, nopopup::Bool=true)
    gmsh.initialize()
    gmsh.model.add(project_name)

    p0 = gmsh.model.geo.addPoint(0, 0, 0)
    l0 = gmsh.model.geo.extrude([(0, p0)], 1, 0, 0, [1], [1], true)
    s0 = gmsh.model.geo.extrude([l0[2]], 0, 1, 0, [1], [1], true)
    v0 = gmsh.model.geo.extrude([s0[2]], 0, 0, 1, [1], [1], true)

    gmsh.model.geo.synchronize()
    # physical groups
    gmsh.model.addPhysicalGroup(3, [v0[2][2]], 30, "hex_LondonClay")
    #boundary conditions
    gmsh.model.addPhysicalGroup(2, [v0[1][2]], 21, "top")
    gmsh.model.addPhysicalGroup(2, [v0[4][2]], 22, "sidex")
    gmsh.model.addPhysicalGroup(2, [v0[5][2]], 23, "sidey")

    # generate mesh
    gmsh.option.setNumber("Mesh.ElementOrder", 0)
    gmsh.model.mesh.generate(3)
    gmsh.write(project_name * ".msh")

    if nopopup == false
        gmsh.fltk.run()
    end

    gmsh.finalize()

    return nothing
end

"""
    unittet(project_name::String, nopopup::Bool)
  
Generate the unit tetrahedron element mesh by Gmsh with the output: *project_name.msh*. 
The size of element is 1 x 1 x 1. 
The name of the materials is *tet_LondonClay*, the former means the shape of the element, while the latter means the material property. 
Three boundaries are defined as the names of *top*, *sidex*, and *sidey*. 
# Arguments:
- `project_name`: name of the project, corresponding to the name of the mesh.
- `nopopup`: not show the gmsh window to visualise the mesh (*default is true*).
"""
function unittet(project_name, nopopup::Bool=true)
    gmsh.initialize()
    gmsh.model.add(project_name)

    p0 = gmsh.model.geo.addPoint(0, 0, 0)
    l0 = gmsh.model.geo.extrude([(0, p0)], 1, 0, 0, [1], [1])
    s0 = gmsh.model.geo.extrude([l0[2]], 0, 1, 0, [1], [1])
    v0 = gmsh.model.geo.extrude([s0[2]], 0, 0, 1, [1], [1])

    gmsh.model.geo.synchronize()
    # physical groups
    gmsh.model.addPhysicalGroup(3, [v0[2][2]], 30, "tet_LondonClay")
    #boundary conditions
    gmsh.model.addPhysicalGroup(2, [v0[1][2]], 21, "top")
    gmsh.model.addPhysicalGroup(2, [v0[4][2]], 22, "sidex")
    gmsh.model.addPhysicalGroup(2, [v0[5][2]], 23, "sidey")

    # generate mesh
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)
    gmsh.model.mesh.generate(3)
    gmsh.write(project_name * ".msh")

    if nopopup == false
        gmsh.fltk.run()
    end

    gmsh.finalize()

    return nothing
end

"""
    unitwed(project_name::String, nopopup::Bool)
  
Generate the unit wedge element mesh by Gmsh with the output: *project_name.msh*. 
The size of element is 1 x 1 x 1. 
The name of the materials is *wed_LondonClay*, the former means the shape of the element, while the latter means the material property. 
Three boundaries are defined as the names of *top*, *sidex*, and *sidey*. 
# Arguments:
- `project_name`: name of the project, corresponding to the name of the mesh.
- `nopopup`: not show the gmsh window to visualise the mesh (*default is true*).
"""
function unitwed(project_name, nopopup::Bool=true)
    gmsh.initialize()
    gmsh.model.add(project_name)

    p0 = gmsh.model.geo.addPoint(0, 0, 0)
    l0 = gmsh.model.geo.extrude([(0, p0)], 1, 0, 0, [1], [1], true)
    l1 = gmsh.model.geo.extrude([(0, p0)], 0, 0, 1, [1], [1], true)
    l2 = gmsh.model.geo.extrude([(0, p0)], 1, 0, 1, [1], [1], true)
    l3 = gmsh.model.geo.extrude([l0[1]], 0, 0, 1, [1], [1], true)
    l4 = gmsh.model.geo.extrude([l1[1]], 1, 0, 0, [1], [1], true)

    c1 = gmsh.model.geo.addCurveLoop([l0[2][2], l3[2][2], -l2[2][2]])
    c2 = gmsh.model.geo.addCurveLoop([l1[2][2], l4[2][2], -l2[2][2]])
    s1 = gmsh.model.geo.addPlaneSurface([c1])
    s2 = gmsh.model.geo.addPlaneSurface([c2])

    # extrude
    v1 = gmsh.model.geo.extrude([(2, s1)], 0, 1, 0, [1], [1], true)
    v2 = gmsh.model.geo.extrude([(2, s2)], 0, 1, 0, [1], [1], true)

    gmsh.model.geo.synchronize()
    # physical groups
    gmsh.model.addPhysicalGroup(3, [v1[2][2], v2[2][2]], 30, "wed_LondonClay")
    #boundary conditions
    gmsh.model.addPhysicalGroup(2, [v2[4][2]], 21, "top")
    gmsh.model.addPhysicalGroup(2, [v1[4][2]], 22, "sidex")
    gmsh.model.addPhysicalGroup(2, [v1[1][2], v2[1][2]], 23, "sidey")

    # generate mesh
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)
    gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)  # 15-node wedge
    gmsh.model.mesh.generate(3)
    gmsh.write(project_name * ".msh")

    if nopopup == false
        gmsh.fltk.run()
    end

    gmsh.finalize()

    return nothing
end

"""
    unittetwed(project_name::String, nopopup::Bool)
  
Generate the unit tetrahedron and wedge element mesh by Gmsh with the output: *project_name.msh*. 
The size of element is 1 x 1 x 1. 
The name of the materials is *tet_LondonClay* and *wed_LondonClay*, the former means the shape of the element, while the latter means the material property. 
Three boundaries are defined as the names of *top*, *sidex*, and *sidey*. 
# Arguments:
- `project_name`: name of the project, corresponding to the name of the mesh.
- `nopopup`: not show the gmsh window to visualise the mesh (*default is true*).
"""
function unittetwed(project_name, nopopup::Bool=true)
    gmsh.initialize()
    gmsh.model.add(project_name)

    p0 = gmsh.model.geo.addPoint(0, 0, 0)
    l0 = gmsh.model.geo.extrude([(0, p0)], 1, 0, 0, [1], [1])
    s0 = gmsh.model.geo.extrude([l0[2]], 0, 0.8, 0, [1], [1])
    v0 = gmsh.model.geo.extrude([s0[2]], 0, 0, 1, [1], [1])  # tetrahedron

    # extrude the surface to be wedge element
    gmsh.model.geo.mesh.setTransfiniteSurface(v0[5][2])
    v1 = gmsh.model.geo.extrude([v0[5]], 0, 0.2, 0, [1], [1], true)  # wedge

    gmsh.model.geo.synchronize()
    # physical groups
    gmsh.model.addPhysicalGroup(3, [v0[2][2]], 30, "tet_LondonClay")
    gmsh.model.addPhysicalGroup(3, [v1[2][2]], 31, "wed_LondonClay")
    # boundary conditions
    gmsh.model.addPhysicalGroup(2, [v0[1][2], v1[5][2]], 21, "top")
    gmsh.model.addPhysicalGroup(2, [v0[4][2], v1[6][2]], 22, "sidex")
    gmsh.model.addPhysicalGroup(2, [v1[1][2]], 23, "sidey")
    # interface tag
    # gmsh.model.addPhysicalGroup(2, [v0[5][2]], 20, "inter")

    # generate mesh
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)
    gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)  # 15-node wedge
    gmsh.model.mesh.generate(3)
    gmsh.write(project_name * ".msh")

    if nopopup == false
        gmsh.fltk.run()
    end

    gmsh.finalize()
end

"""
    excavation(project_name::String, nopopup::Bool)
  
Generate the unit element mesh by Gmsh with the output: *project_name.msh*. 
The size of element is 1 x 1 x 1. 
The name of the soils is *tet_LondonClay*, the former means the shape of the element, while the latter means the material property. 
The name of the diaphragm wall is *wed_DiaphragmWall*. 
The props are added during the process of excavation. 

# Arguments:
- `project_name`: name of the project, corresponding to the name of the mesh.
- `nopopup`: not show the gmsh window to visualise the mesh (*default is true*).
"""
function excavation(project_name::String, nopopup::Bool=true)
    gmsh.initialize()
    gmsh.model.add(project_name)

    N = 4  # the division of excavated soils
    Nu = 4  # the division of unexcavated soils
    N2 = 1 # the division of the walls
    B = 0.48
    T = 0.02
    p0 = gmsh.model.geo.addPoint(0, 0, 0)

    # excavation
    l1 = gmsh.model.geo.extrude([0, p0], B, 0, 0, [N], [])
    s1 = gmsh.model.geo.extrude([l1[2]], 0, B, 0, [N], [])
    v1 = gmsh.model.geo.extrude([s1[2]], 0, 0, 0.4, [N], [])

    # For the y-direction retaining wall
    # extrude the point to obtain the surface
    p2 = gmsh.model.geo.extrude([(0, 3)], 0, T, 0, [N2], [])
    l2 = gmsh.model.geo.extrude([p2[1]], B + T, 0, 0, [N], [])
    s2 = gmsh.model.geo.extrude([l2[2]], 0, 0, 0.4, [N], [])

    # extrude the lines for the side surfaces
    s3 = gmsh.model.geo.extrude([p2[2]], 0, 0, 0.4, [N], [])
    p3 = gmsh.model.geo.extrude([(0, 4)], T, T, 0, [N2], [])
    s4 = gmsh.model.geo.extrude([p3[2]], 0, 0, 0.4, [N], [])

    s5 = gmsh.model.geo.addCurveLoop([9, 34, 30, -39])
    s6 = gmsh.model.geo.addCurveLoop([28, 29, -38, -2])
    s7 = gmsh.model.geo.addPlaneSurface([s5])
    s8 = gmsh.model.geo.addPlaneSurface([s6])
    v2 = gmsh.model.geo.addSurfaceLoop([22, 44, 43, 33, 37, 42])
    v3 = gmsh.model.geo.addVolume([v2])

    gmsh.model.geo.mesh.setTransfiniteSurface(22, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(33, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(37, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(42, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(43, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(44, "Left")
    gmsh.model.geo.mesh.setTransfiniteVolume(v3, [3, 14, 17, 15, 4, 10, 18, 16])

    # For the x-direction retaining wall
    xp2 = gmsh.model.geo.extrude([l1[1]], T, 0, 0, [N2], [])
    xs2 = gmsh.model.geo.extrude([s2[3]], 0, -B - T, 0, [N], [])
    xs3 = gmsh.model.geo.extrude([xp2[2]], 0, 0, 0.4, [N], [])
    #
    xs5 = gmsh.model.geo.addCurveLoop([8, 39, 48, -50])
    xs6 = gmsh.model.geo.addCurveLoop([4, 38, 47, -45])
    xs7 = gmsh.model.geo.addPlaneSurface([xs5])
    xs8 = gmsh.model.geo.addPlaneSurface([xs6])
    xv2 = gmsh.model.geo.addSurfaceLoop([54, 49, 53, 55, 18, 42])
    xv3 = gmsh.model.geo.addVolume([xv2])

    gmsh.model.geo.mesh.setTransfiniteSurface(18, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(49, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(53, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(54, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(55, "Left")
    gmsh.model.geo.mesh.setTransfiniteVolume(xv3, [4, 10, 18, 16, 2, 6, 21, 19])

    # The second layer
    ss0 = gmsh.model.geo.extrude([v1[1]], 0, 0, 0.2, [N], [])  # The third excavation
    ss1 = gmsh.model.geo.extrude([s3[1]], 0, 0, 0.2, [N], [])
    ss2 = gmsh.model.geo.extrude([s2[1]], 0, 0, 0.2, [N], [])
    ss3 = gmsh.model.geo.extrude([s4[1]], 0, 0, 0.2, [N], [])
    ss4 = gmsh.model.geo.extrude([xs2[3]], 0, 0, 0.2, [N], [])
    ss5 = gmsh.model.geo.extrude([xs3[1]], 0, 0, 0.2, [N], [])

    ss6 = gmsh.model.geo.addCurveLoop([59, 78, 82, -86])
    ss7 = gmsh.model.geo.addCurveLoop([86, 90, -94, 58])
    ss8 = gmsh.model.geo.addPlaneSurface([ss6])
    ss9 = gmsh.model.geo.addPlaneSurface([ss7])

    sv1 = gmsh.model.geo.addSurfaceLoop([98, 81, 85, 72, 89, 43])
    sv2 = gmsh.model.geo.addVolume([sv1])
    sv3 = gmsh.model.geo.addSurfaceLoop([99, 93, 97, 68, 89, 54])
    sv4 = gmsh.model.geo.addVolume([sv3])

    gmsh.model.geo.mesh.setTransfiniteSurface(18, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(49, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(68, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(72, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(81, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(85, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(89, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(93, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(97, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(98, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(99, "Left")
    gmsh.model.geo.mesh.setTransfiniteVolume(sv2, [14, 31, 33, 17, 10, 27, 35, 18])
    gmsh.model.geo.mesh.setTransfiniteVolume(sv4, [10, 27, 35, 18, 6, 23, 37, 21])

    # The third layer
    ts0 = gmsh.model.geo.extrude([ss0[1]], 0, 0, 0.2, [N], [])  # The second excavation
    ts1 = gmsh.model.geo.extrude([ss1[1]], 0, 0, 0.2, [N], [])
    ts2 = gmsh.model.geo.extrude([ss2[1]], 0, 0, 0.2, [N], [])
    ts3 = gmsh.model.geo.extrude([ss3[1]], 0, 0, 0.2, [N], [])
    ts4 = gmsh.model.geo.extrude([ss4[1]], 0, 0, 0.2, [N], [])
    ts5 = gmsh.model.geo.extrude([ss5[1]], 0, 0, 0.2, [N], [])

    ts6 = gmsh.model.geo.addCurveLoop([103, 122, 126, -130])
    ts7 = gmsh.model.geo.addCurveLoop([102, 130, 134, -138])
    ts8 = gmsh.model.geo.addPlaneSurface([ts6])
    ts9 = gmsh.model.geo.addPlaneSurface([ts7])

    tv1 = gmsh.model.geo.addSurfaceLoop([142, 125, 129, 98, 116, 133])
    tv2 = gmsh.model.geo.addVolume([tv1])
    tv3 = gmsh.model.geo.addSurfaceLoop([133, 143, 137, 141, 112, 99])
    tv4 = gmsh.model.geo.addVolume([tv3])

    gmsh.model.geo.mesh.setTransfiniteSurface(112, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(116, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(125, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(129, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(133, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(137, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(141, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(142, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(143, "Left")
    gmsh.model.geo.mesh.setTransfiniteVolume(tv2, [31, 47, 49, 33, 27, 43, 51, 35])
    gmsh.model.geo.mesh.setTransfiniteVolume(tv4, [27, 43, 51, 35, 23, 39, 53, 37])

    # The forth layer
    fs0 = gmsh.model.geo.extrude([ts0[1]], 0, 0, 0.2, [N], [])  # The first excavation
    fs1 = gmsh.model.geo.extrude([ts1[1]], 0, 0, 0.2, [N], [])
    fs2 = gmsh.model.geo.extrude([ts2[1]], 0, 0, 0.2, [N], [])
    fs3 = gmsh.model.geo.extrude([ts3[1]], 0, 0, 0.2, [N], [])
    fs4 = gmsh.model.geo.extrude([ts4[1]], 0, 0, 0.2, [N], [])
    fs5 = gmsh.model.geo.extrude([ts5[1]], 0, 0, 0.2, [N], [])

    fs6 = gmsh.model.geo.addCurveLoop([147, 166, 170, -174])
    fs7 = gmsh.model.geo.addCurveLoop([174, 178, -182, 146])
    fs8 = gmsh.model.geo.addPlaneSurface([fs6])
    fs9 = gmsh.model.geo.addPlaneSurface([fs7])

    fv1 = gmsh.model.geo.addSurfaceLoop([173, 186, 169, 160, 177, 142])
    fv2 = gmsh.model.geo.addVolume([fv1])
    fv3 = gmsh.model.geo.addSurfaceLoop([187, 181, 185, 177, 156, 143])
    fv4 = gmsh.model.geo.addVolume([fv3])

    gmsh.model.geo.mesh.setTransfiniteSurface(156, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(160, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(169, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(173, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(177, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(181, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(185, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(186, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(187, "Left")
    gmsh.model.geo.mesh.setTransfiniteVolume(fv2, [47, 63, 65, 49, 43, 59, 67, 51])
    gmsh.model.geo.mesh.setTransfiniteVolume(fv4, [43, 59, 67, 51, 39, 55, 69, 53])

    # unexcavated soils around the first layer
    up1 = gmsh.model.geo.extrude([xp2[1]], 0.5, 0, 0, [Nu], [])
    up2 = gmsh.model.geo.extrude([up1[1]], 0, 1, 0, [Nu], [])
    up3 = gmsh.model.geo.extrude([p2[1]], 0, 0.5, 0, [Nu], [])
    up4 = gmsh.model.geo.extrude([up3[1]], 1, 0, 0, [Nu], [])
    us1 = gmsh.model.geo.addCurveLoop([188, 189, -191, -190, 29, 47])
    us2 = gmsh.model.geo.addPlaneSurface([us1])
    uv1 = gmsh.model.geo.extrude([(2, us2)], 0, 0, 0.4, [N], [])
    uv2 = gmsh.model.geo.extrude([uv1[1]], 0, 0, 0.2, [N], [])
    uv3 = gmsh.model.geo.extrude([uv2[1]], 0, 0, 0.2, [N], [])
    uv4 = gmsh.model.geo.extrude([uv3[1]], 0, 0, 0.2, [N], [])

    # coordinates of the end points of props
    pp1 = gmsh.model.geo.addPoint(0.24, 0.48, 0.9)
    lp1 = gmsh.model.geo.extrude([(0, pp1)], 0.24, -0.24, 0, [1], [], true)

    gmsh.model.geo.synchronize()
    # embed the end points to the wall surface
    gmsh.model.mesh.embed(0, [pp1], 2, 160)
    gmsh.model.mesh.embed(0, [lp1[1][2]], 2, 156)
    # recombine the mesh to generate the wedge elements
    gmsh.model.mesh.setRecombine(2, 53)
    gmsh.model.mesh.setRecombine(2, 54)
    gmsh.model.mesh.setRecombine(2, 55)
    gmsh.model.mesh.setRecombine(2, 42)
    gmsh.model.mesh.setRecombine(2, 43)
    gmsh.model.mesh.setRecombine(2, 37)
    gmsh.model.mesh.setRecombine(2, 44)
    # recombine for the second layer
    gmsh.model.mesh.setRecombine(2, 81)
    gmsh.model.mesh.setRecombine(2, 98)
    gmsh.model.mesh.setRecombine(2, 89)
    gmsh.model.mesh.setRecombine(2, 99)
    gmsh.model.mesh.setRecombine(2, 97)
    # recombine for the third layer
    gmsh.model.mesh.setRecombine(2, 125)
    gmsh.model.mesh.setRecombine(2, 142)
    gmsh.model.mesh.setRecombine(2, 133)
    gmsh.model.mesh.setRecombine(2, 143)
    gmsh.model.mesh.setRecombine(2, 141)
    # recombine for the forth layer
    gmsh.model.mesh.setRecombine(2, 169)
    gmsh.model.mesh.setRecombine(2, 186)
    gmsh.model.mesh.setRecombine(2, 177)
    gmsh.model.mesh.setRecombine(2, 187)
    gmsh.model.mesh.setRecombine(2, 185)

    # define the physical groups
    # geostatic state for the whole model
    gmsh.model.addPhysicalGroup(3, [1, 4, 7, 10, 13, 14, 15, 16], 1, "tet_LondonClay")
    gmsh.model.addPhysicalGroup(3, [2, 3, 5, 6, 8, 9, 11, 12], 2, "wed_DiaphragmWall")
    gmsh.model.addPhysicalGroup(1, [lp1[2][2]], 101, "prop01")

    # excavation sequence
    # 1. interface element
    gmsh.model.addPhysicalGroup(2, [33, 85, 129, 173, 49, 93, 137, 181, 22, 18, 72, 68, 116, 112, 160, 156], 10, "inter0")
    gmsh.model.addPhysicalGroup(2, [33, 85, 129, 173], 11, "inter_yn")  # remaining interface in y-n
    gmsh.model.addPhysicalGroup(2, [49, 93, 137, 181], 12, "inter_xn")  # remaining interface in x-n
    gmsh.model.addPhysicalGroup(2, [22, 72], 13, "inter_yp")  # excavation side of remaining interface in y-p
    gmsh.model.addPhysicalGroup(2, [18, 68], 14, "inter_xp")  # excavation side of remaining interface in x-p
    gmsh.model.addPhysicalGroup(2, [116], 15, "inter_02_yp")  # second excavation interface in y-p
    gmsh.model.addPhysicalGroup(2, [112], 16, "inter_02_xp")  # second excavation interface in x-p 
    gmsh.model.addPhysicalGroup(2, [160], 17, "inter_01_yp")  # first excavation interface in y-p
    gmsh.model.addPhysicalGroup(2, [156], 18, "inter_01_xp")  # first excavation interface in x-p

    # 2. First excavation
    gmsh.model.addPhysicalGroup(3, [10], 101, "tet_excav01")  # excavated soils
    gmsh.model.addPhysicalGroup(2, [121], 102, "exc01bcs")  # boundary of excavation (only the bottom, ignore interface)
    gmsh.model.addPhysicalGroup(3, [1, 4, 7, 13, 14, 15, 16], 103, "tet_rem01")
    gmsh.model.addPhysicalGroup(3, [2, 3, 5, 6, 8, 9, 11, 12], 104, "wed_rem01")
    gmsh.model.addPhysicalGroup(2, [121, 156, 160], 105, "exc01bc")  # boundary of excavation (not consider interface)

    # 3. Second excavation
    gmsh.model.addPhysicalGroup(3, [7, 10], 201, "tet_excav02")  # excavated soils
    gmsh.model.addPhysicalGroup(2, [77], 202, "exc02bcs")  # boundary of excavation (only the bottom, ignore interface)
    gmsh.model.addPhysicalGroup(3, [1, 4, 13, 14, 15, 16], 203, "tet_rem02")
    gmsh.model.addPhysicalGroup(3, [2, 3, 5, 6, 8, 9, 11, 12], 204, "wed_rem02")
    gmsh.model.addPhysicalGroup(2, [77, 112, 116, 156, 160], 205, "exc02bc")  # boundary of excavation (not consider interface)

    gmsh.option.setNumber("Mesh.ElementOrder", 2)  # second order
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)
    gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)  # 15-node wedge element
    gmsh.model.mesh.generate(3)
    gmsh.write(project_name * ".msh")

    if nopopup == false
        gmsh.fltk.run()
    end

    gmsh.finalize()
end
