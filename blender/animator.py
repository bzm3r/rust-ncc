from typing import Tuple

import bpy
from bpy.types import GreasePencil, GPencilLayer
import numpy as np

NUM_FRAMES = 30
FRAMES_SPACING = 1  # distance between frames
bpy.context.scene.frame_start = 0
bpy.context.scene.frame_end = NUM_FRAMES * FRAMES_SPACING


def create_gp(gp_data: GreasePencil, name) -> GreasePencil:
    if name not in bpy.context.scene.objects:
        gp_obj = bpy.data.objects.new(name, gp_data)
        bpy.context.scene.objects[-1].name = name
        bpy.context.scene.collection.objects.link(gp_obj)

    return bpy.context.scene.objects.get(name)


def create_gp_layer(gp_obj, layer_name, clear_layer) -> GPencilLayer:
    # Get grease pencil layer or create one if none exists
    if layer_name in gp_obj.data.layers:
        gp_layer = gp_obj.data.layers.get(layer_name)
    else:
        gp_layer = gp_obj.data.layers.new(layer_name, set_active=True)

    if clear_layer:
        gp_layer.clear()

    return gp_layer


def init_gp(name) -> Tuple[GreasePencil, GreasePencil]:
    gp_data = bpy.data.grease_pencils.new(name)
    gp_obj = create_gp(gp_data, name)

    return gp_data, gp_obj


def stroke_polyline(gp_frame: GPencilFrame, vertices, line_width) -> GPencilStroke:
    gp_stroke = gp_frame.strokes.new()
    gp_stroke.points.add(len(vertices))

    for ix, vertex in enumerate(vertices):
        gp_stroke.points[ix].co = vertex

    gp_stroke.line_width = line_width
    return gp_stroke


def create_gp_material(name: str, color) -> bpy.types.Material:
    gp_mat = bpy.data.materials.new(name)
    bpy.data.materials.create_gpencil_data(gp_mat)
    gp_mat.grease_pencil.color = color

    return gp_mat


def x_rotor(theta):
    return np.array([[1, 0, 0],
                     [0, np.cos(theta), -np.sin(theta)],
                     [0, np.sin(theta), np.cos(theta)]])


def y_rotor(theta):
    return np.array([[np.cos(theta), 0, np.sin(theta)],
                     [0, 1, 0],
                     [-np.sin(theta), 0, np.cos(theta)]])


def z_rotor(theta):
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                     [np.sin(theta), np.cos(theta), 0],
                     [0, 0, 1]])


gp_dat, gp_obj = init_gp("TestPencil")
gp_layer = create_gp_layer(gp_obj, "TestLayer", True)
gp_mat = create_gp_material("BlackLine", [0., 0., 0., 1.])
gp_dat.materials.append(gp_mat)

base_vertex = np.zeros(3)
unit_vertex = np.array([1.0, 0.0, 0.0])
for frame in range(NUM_FRAMES):
    gp_frame = gp_layer.frames.new(frame * FRAMES_SPACING)
    rotor = z_rotor((frame / NUM_FRAMES) * 2 * np.pi)
    print("rotor: {}", rotor)
    rotated_vertex = np.dot(rotor, unit_vertex)
    print("frame {}: {}", frame, rotated_vertex)
    gp_stroke = stroke_polyline(gp_frame, [base_vertex, rotated_vertex], 10)

print("Success!")
