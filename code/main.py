import numpy as np
import open3d as o3d
import math
import sympy
from sympy import Point3D
from sympy.abc import L
from sympy.geometry import Line3D, Segment3D, Plane, Polygon

from los import sat_pos, convert_crs, tle_to_itrs
from tle import tle_json

def create_triangle_mesh(input_file):
    #Read triangle mesh from obj file
    mesh = o3d.io.read_triangle_mesh(input_file,True)

    #Print information on triangle mesh
    print(mesh)
    print('Vertices:')
    print(np.asarray(mesh.vertices))
    print('Triangles:')
    print(np.asarray(mesh.triangles))

    return mesh

def visualize(mesh):
    # Visualization preparation
    mesh.compute_vertex_normals()
    #print(np.asarray(mesh.triangle_normals))
    mesh.paint_uniform_color([0.24, 0.70, 0.44])

    # Visualize 3D model
    vis = o3d.visualization.Visualizer()
    vis.create_window(window_name='obj visualization', width=1000, height=1000)
    vis.add_geometry(mesh)
    vis.run()

def create_grid(input_mesh, cellsize_grid):
    #Create grid from mesh extents
    #first convert array of vertices to list of vertices
    list_pts_prep = []
    list_pts = []
    for coor in input_mesh.vertices:
        list_pts_prep.append([coor.tolist()])

    #list of x and y coordinates as integers
    for coor in list_pts_prep:
        for x,y,z in coor:
            list_pts.append([int(x), int(y)])

    #extents: upper x, upper y, lower x, lower y
    ux = 418676
    uy = 5653691
    lx = 418267
    ly = 5653400

    """
    print('length of list_pts_prep = ', len(list_pts_prep))
    for a in range(len(list_pts_prep)):
        if a == 0:
            ux = list_pts[a][0]
            uy = list_pts[a][1]
            lx = list_pts[a][0]
            ly = list_pts[a][1]
        else:
            if list_pts[a][0] > ux:
                ux = list_pts[a][0]
            if list_pts[a][1] > uy:
                uy = list_pts[a][1]
            if list_pts[a][0] < lx:
                lx = list_pts[a][0]
            if list_pts[a][1] < ly:
                ly = list_pts[a][1]
    """

    print('ux is', ux)
    print('uy is', uy)
    print('lx is', lx)
    print('ly is', ly)

    ncols = math.ceil(((uy - ly) / cellsize_grid))
    nrows = math.ceil(((ux - lx) / cellsize_grid))

    # Create list of centerpoints of cells
    centerpoints = []

    first_cp_y = ly + 0.5 * cellsize_grid
    first_cp_x = lx + 0.5 * cellsize_grid

    for j in range(1, nrows + 1):
        for i in range(1, ncols + 1):
            centerpoints.append([i * cellsize_grid + first_cp_x, j * cellsize_grid + first_cp_y])

    array_centerpoints = np.array(centerpoints)

    return ncols, nrows, lx, ly, centerpoints

def read_height_model(height_model_file):
    coordinates = {}
    xyz = open(height_model_file)

    for line in xyz:
        x,y,z = line.split()
        coordinates[int(float(x)),int(float(y))] = float(z)
    xyz.close()
    return coordinates

def write_raster(ncols, nrows, lx, ly, cellsize_grid, nested_satellite_values, file):
    with open(file, 'w+') as fh:
        fh.write('NCOLS ' + str(ncols) + '\n')
        fh.write('NROWS ' + str(nrows) + '\n')
        fh.write('XLLCORNER ' + str(lx) + '\n')
        fh.write('YLLCORNER ' + str(ly) + '\n')
        fh.write('CELLSIZE ' + str(cellsize_grid) + '\n')
        fh.write('NODATA_VALUE -9999')
        for i in reversed(nested_satellite_values):
            fh.write("\n")
            for point in i:
                fh.write(str(point) + " ")

    print('File written to', file)


def satellite_lines(xyz):
    satellites = sat_pos()
    GPS = []
    Galileo = []
    i = 0
    for constellation in satellites:
        for satellite in constellation:
            values = list(satellite.values())
            if i == 0:
                GPS.append(values[0][2])
            if i == 1:
                Galileo.append(values[0][2])
        i += 1
    print(GPS)
    print(Galileo)

    test_point = [418494.00, 5653466.00, 350.53]
    lines = []
    for satellite in Galileo:
        line = [test_point, [satellite[0], satellite[1], satellite[2]]]
        lines.append(line)
    return lines
'''
    test_point = Point3D(418152.00, 5653486.00, 340.56)
    lines = []
    lines_equations = []
    for satellite in Galileo:
        line = (test_point, Point3D(satellite[0], satellite[1], satellite[2]))
        lines.append(line)
        lines_equations.append(line.equation())
    print(lines_equations)

    test_point2 = [418152.00, 5653486.00, 340.56]
    points2 = []
    points2.append(test_point2)
    lines2 = []
    i = 1
    for satellite in Galileo:
        points2.append([satellite[0], satellite[1], satellite[2]])
        lines2.append([0,i])
        i+=1

    line_set = o3d.geometry.LineSet(
        points=o3d.utility.Vector3dVector(points2),
        lines=o3d.utility.Vector2iVector(lines2),
    )
    print(line_set)
'''

def mesh_to_triangles(mesh):
    list_vertices = []
    for i in mesh.vertices:
        list_vertices.append(i.tolist())
    triangles = []
    for triangle in mesh.triangles:
        point_1 = list_vertices[int(triangle[0])]
        point_2 = list_vertices[int(triangle[1])]
        point_3 = list_vertices[int(triangle[2])]
        triangles.append([point_1, point_2, point_3])

    return triangles

def sign_of_volume(a,b,c,d):
    B = np.array([b[0] - a[0], b[1] - a[1], b[2] - a[1]])
    C = np.array([c[0] - a[0], c[1] - a[1], c[2] - a[1]])
    D = np.array([d[0] - a[0], d[1] - a[1], d[2] - a[1]])
    return np.sign((1.0 / 6.0) * np.dot(np.cross(B, C), D))

def test(line, triangles):
    intersections = []
    q1 = line[0]
    q2 = line[1]
    for i in triangles:
        p1 = i[0]
        p2 = i[1]
        p3 = i[2]

        if (sign_of_volume(q1,p1,p2,p3) == sign_of_volume(q2,p1,p2,p3)) & (sign_of_volume(q1,q2,p1,p2) == sign_of_volume(q1,q2,p2,p3) == sign_of_volume(q1,q2,p3,p1)):
            return True
        else:
            continue




def main():
    #define input files
    input_file = 'LoD2_32_418_5653_1_NW.obj'
    height_model = 'dgm1_32_418_5653_1_nw.xyz'
    #define output file
    output_file = 'out.asc'
    #define cell size for grid
    cellsize_grid = 2

    mesh = create_triangle_mesh(input_file)
    visualize(mesh)
    ncols, nrows, lx, ly, array_centerpoints = create_grid(mesh, cellsize_grid)

    list_xyz = read_height_model(height_model)

    list_centerpoints_z = []
    for [x,y] in array_centerpoints:
        list_centerpoints_z.append([x,y,list_xyz[int(x),int(y)]])

    list_z_values = []
    for i in range(len(array_centerpoints)):
        list_z_values.append(list_centerpoints_z[i][2])

    nested_z_values = []
    for i in range(0, len(list_z_values), (ncols)):
        nested_z_values.append(list_z_values[i:(i+ncols)])

    lines = satellite_lines(list_centerpoints_z)
    triangles = mesh_to_triangles(mesh)

    for line in lines:
        if test(line, triangles) == True:
            print("this line has an intersection with the mesh")
        else:
            print("this line does not have intersections with the mesh")


    write_raster(ncols, nrows, lx, ly, cellsize_grid, nested_z_values, output_file)



if __name__ == "__main__":
    main()

