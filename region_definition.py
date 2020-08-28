from geometric_objects import Point, Triangle
import math
import numpy as np

import collections

import pdb


def get_region_map(NT, file_name):
    # Constants

    a = math.sqrt(2.) / 4.
    sqrt3 = math.sqrt(3.)  # square root of 3
    Nx = NT + 1  # x-direction triangle count.  1 extra needed
    Ny = int(math.ceil(NT * 2 / sqrt3))  # y-direction triangle count
    # z-direction count of tetrahedra stacked above a triangle
    Nz = int(math.ceil(float(NT) / (3. * a) + 1))

    (vertices, triangles) = parse_off_file(file_name)
    intersecting_triangles = map_grid_to_intersecting_triangles(
        NT, vertices, triangles)
    region = region_definition(intersecting_triangles, NT, Nx, Ny, Nz)
    return region
    # run_test_cases()

# Returns tuple of array of vertices, array of triangles.
# Vertices are an array of point objects.
# Triangles are three indexes, referencing points in the vertices array.


def parse_off_file(off_file):
    with open(off_file, 'r') as in_f:
        vertice_count = 0
        vertices = []
        triangles = []
        for index, line in enumerate(open(off_file).read().splitlines()):
            if index == 0:
                pass
            elif index == 1:
                count_info = line.split(' ')
                vertice_count = int(count_info[0])
            elif index < 2 + vertice_count:
                vertex_info = line.split(' ')
                vertex = Point(vertex_info[0], vertex_info[1], vertex_info[2])
                vertices.append(vertex)
            else:
                triangle_info = line.split(' ')
                vertex_1 = int(triangle_info[1])
                vertex_2 = int(triangle_info[2])
                vertex_3 = int(triangle_info[3])
                triangles.append([vertex_1, vertex_2, vertex_3])
        return (vertices, triangles)

# For each point on the grid, take the x,y value and figure out which triangles
# are intercepted by a line coming up from that point.
def map_grid_to_intersecting_triangles(nt, vertices, triangles):
    intersecting_triangles = collections.defaultdict(list)

    for triangle in triangles:
        vertex_0 = vertices[triangle[0]]
        vertex_1 = vertices[triangle[1]]
        vertex_2 = vertices[triangle[2]]

        (i_min_bound, i_max_bound, j_min_bound, j_max_bound) = bounding_rectangle(
            nt, vertex_0, vertex_1, vertex_2)

        for i in range(i_min_bound, i_max_bound + 1):
            for j in range(j_min_bound, j_max_bound + 1):
                cartesian_point = Point.create_from_grid_coordinates(
                    i, j, 0, *step_sizes(nt))
                if point_in_triangle(
                        cartesian_point,
                        vertex_0,
                        vertex_1,
                        vertex_2):
                    intersecting_triangles[(i, j)].append(
                        [vertex_0, vertex_1, vertex_2])

    for key in intersecting_triangles:
        if len(intersecting_triangles[key]) % 2 != 0:
            assert(
                "ERROR: All vertical lines from the x-y plane must intersect the surface at an even number of points")

    return intersecting_triangles


# From https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle 
def point_in_triangle(p, p0, p1, p2):
    dX = p.x-p2.x;
    dY = p.y-p2.y;
    dX21 = p2.x-p1.x;
    dY12 = p1.y-p2.y;
    D = dY12*(p0.x-p2.x) + dX21*(p0.y-p2.y);
    s = dY12*dX + dX21*dY;
    t = (p2.y-p0.y)*dX + (p0.x-p2.x)*dY;
    if (D<0):
        return s<=0 and t<=0 and s+t>=D
    return s>=0 and t>=0 and s+t<=D;


def count_triangles_above_grid_coordinate(triangles_to_check, NT, i, j, k):
    cartesian_point = Point.create_from_grid_coordinates(
        i, j, k, *step_sizes(NT))
    above_count = 0
    for triangle in triangles_to_check:
        z_intercept = z_intercept_coordinate(
            cartesian_point.x,
            cartesian_point.y,
            triangle[0],
            triangle[1],
            triangle[2])
        if (z_intercept > cartesian_point.z):
            above_count += 1
    return above_count


def z_intercept_coordinate(x0, y0, v0, v1, v2):
    return z_inter_coord(
        x0,
        y0,
        v0.x,
        v0.y,
        v0.z,
        v1.x,
        v1.y,
        v1.z,
        v2.x,
        v2.y,
        v2.z)

# Compute z coordinate where the vertical line
# above (x0, y0) crosses the triangle
# with vertices (x1,y1,z1) , (x2,y2,z2) , (x3,y3,z3)


def z_inter_coord(x0, y0, x1, y1, z1, x2, y2, z2, x3, y3, z3):
    # compute vectors of triangle sides in xy plane
    v0 = [x3 - x1, y3 - y1, 0]
    v1 = [x2 - x1, y2 - y1, 0]
# compute vector of vertical line above (x0, y0) in xy-plane
    v2 = [x0 - x1, y0 - y1, 0]
    dot00 = (x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1)
    dot01 = (x3 - x1) * (x2 - x1) + (y3 - y1) * (y2 - y1)
    dot02 = (x3 - x1) * (x0 - x1) + (y3 - y1) * (y0 - y1)
    dot11 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)
    dot12 = (x2 - x1) * (x0 - x1) + (y2 - y1) * (y0 - y1)
    dot22 = (x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1)
# revisit - if triangle is degenerate - shouldn't happen
    if ((dot00 * dot11 - dot01 * dot01) == 0):
        return z1
    invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
# with barycentric coords, v2 = u*v0 + v*v1
    u = (dot11 * dot02 - dot01 * dot12) * invDenom
    v = (dot00 * dot12 - dot01 * dot02) * invDenom
# z coordinate is computed using these barycentric coordinates.
    return z1 + u * (z3 - z1) + v * (z2 - z1)


def bounding_rectangle(nt, vertex_0, vertex_1, vertex_2):
    dx = 1. / float(nt)
    dy = dx * math.sqrt(3.) / 2.

    x_min = min(vertex_0.x, vertex_1.x, vertex_2.x)
    x_max = max(vertex_0.x, vertex_1.x, vertex_2.x)
    y_min = min(vertex_0.y, vertex_1.y, vertex_2.y)
    y_max = max(vertex_0.y, vertex_1.y, vertex_2.y)

    i_min_bound = int(math.floor((x_min / dx) - 1))
    i_max_bound = int(math.ceil((x_max / dx) + 1))
    j_min_bound = int(math.floor((y_min / dy) - 1))
    j_max_bound = int(math.ceil((y_max / dy) + 1))

    return (i_min_bound, i_max_bound, j_min_bound, j_max_bound)


def step_sizes(nt):
    a = math.sqrt(2.) / 4.

    dx = 1. / float(nt)
    dy = dx * math.sqrt(3.) / 2.
    dz = 3. * a / float(nt)

    return (dx, dy, dz)


def region_definition(intersecting_triangles, NT, i_max, j_max, k_max):
    definition = {}
    for i in range(i_max):
        for j in range(j_max):
            for k in range(k_max):
                a = (i, j)
                intersecting = intersecting_triangles[a]

                count = count_triangles_above_grid_coordinate(
                    intersecting, NT, i, j, k)
                inside_flag = -1 if (count % 2 == 1) else 1
                definition[(i, j, k)] = inside_flag
    return definition

######################################################################
# Test Cases
######################################################################

def run_test_cases():
    test_bounding_rectangle()
    test_intersecting_triangles()
    print("All test cases passed")


def test_bounding_rectangle():
    p1 = Point(0, 0, .6)
    p2 = Point(0, .01, .6)
    p3 = Point(.01, 0, .6)
    assert bounding_rectangle(100, p1, p2, p3) == (-1, 2, -1, 3)
    print("Bounding Rectangle test case passed")


def test_intersecting_triangles():
    v0 = Point(0.49, .085, .6)
    v1 = Point(0.51, .085, .6)
    v2 = Point(0.50, .087, .6)
    v3 = Point(0.49, .095, .6)
    v4 = Point(0.50, .095, .6)
    v5 = Point(0.495, .096, .6)
    vertices = [v0, v1, v2, v3, v4, v5]

    t0 = [0, 1, 2]
    t1 = [3, 4, 5]

    triangles = [t0, t1]

    intersecting_triangles = map_grid_to_intersecting_triangles(
        100, vertices, triangles)
    assert(list(intersecting_triangles.keys())) == [(50, 10), (50, 11)]
    print("Intersecting triangles test case passed")
