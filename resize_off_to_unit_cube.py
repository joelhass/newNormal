
import math

import sys
from geometric_objects import Point, Triangle
import pdb
import trimesh

UNIT_BOX_PADDING = .05

def resize_off_to_unit_cube(input_off_file, translated_off_file_name):

    mesh = trimesh.load(input_off_file)
    corners = trimesh.bounds.corners(mesh.bounding_box.bounds)
    xMin = corners[0][0]
    xMax = corners[1][0]
    yMin = corners[0][1]
    yMax = corners[2][1]
    zMin = corners[0][2]
    zMax = corners[4][2]

    xDif = xMax - xMin
    yDif = yMax - yMin
    zDif = zMax - zMin
    boxSize = max(xDif, yDif, zDif)  # cube side to enclose region
    stretch_factor = (0.9) / boxSize  # resizes to unit cube

    with open(input_off_file, 'r') as in_f, open(translated_off_file_name, 'w+') as out_f:
        for index, line in enumerate(in_f):
            if index == 0:
                out_f.write(line)
            elif index == 1:
                count_info = line.split()
                vertex_count = int(count_info[0])
                out_f.write(line)
            elif index < 2 + vertex_count:
                point = line.split(' ')
                translated_point = point
                translated_point[0] = str(
                    (float(point[0]) - xMin) * stretch_factor + UNIT_BOX_PADDING)
                translated_point[1] = str(
                    (float(point[1]) - yMin) * stretch_factor + UNIT_BOX_PADDING)
                translated_point[2] = str(
                    (float(point[2]) - zMin) * stretch_factor + UNIT_BOX_PADDING)
                out_f.write(' '.join(translated_point))
                out_f.write('\n')
            else:
                out_f.write(line)

    """
    In some cases the input off file has a carriage return between lines, in some case it doesn't. 
    To handle this, we add an extra carriage return, and then go through and delete empty lines.
    """
    remove_empty_lines(translated_off_file_name)

    return (stretch_factor, xMin, yMin, zMin)

def expand_off_to_original_dimenstions(input_off_file, output_off_file, stretch_factor, xMin, yMin, zMin):
    mesh = trimesh.load(input_off_file)

    vertices = mesh.vertices
    expanded = []
    for index, vertex in enumerate(mesh.vertices):
        expanded_x = (vertex[0] - UNIT_BOX_PADDING) * stretch_factor + xMin
        expanded_y = (vertex[1] - UNIT_BOX_PADDING) * stretch_factor + yMin
        expanded_z = (vertex[2] - UNIT_BOX_PADDING) * stretch_factor + zMin
        mesh.vertices[index] = [expanded_x, expanded_y, expanded_z]

    with open(output_off_file, "w") as output_off:
        output_off.write(trimesh.exchange.off.export_off(mesh))

def remove_empty_lines(filename):
    # Overwrite the file, removing empty lines
    with open(filename, 'r+') as f:
        lines = f.readlines()
        f.seek(0)
        f.writelines(line for line in lines if line.strip())
        f.truncate()
