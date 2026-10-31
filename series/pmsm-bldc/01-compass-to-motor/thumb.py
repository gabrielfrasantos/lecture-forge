"""Thumbnail hero for lesson 1: motor cutaway with a rotating-field arrow, and a compass deflected by a wire."""
import numpy as np
from engine import Camera, Scene3D, C, tube, rgb
from models3d import motor_parts, compass

D = np.pi / 180
HEADLINE = ["WHERE DOES", "TORQUE", "COME FROM?"]
HIGHLIGHT = 1
ACCENT = C["current"]
LESSON = 1
KICKER = "PMSM & BLDC motors"


def hero(ax):
    # motor cutaway, large, right side (kept clear of the bottom-right duration badge)
    cam = Camera(target=(0.05, 0, 0.0), az=18, el=22, dist=5.6, fov=30, cx=1480, cy=590)
    sc = Scene3D(cam)
    for faces, col, kind in motor_parts(rotor_angle=0.25, cut=(20 * D, 140 * D)):
        sc.mesh(faces, col)
    # rotation arrow around the shaft
    ang = np.linspace(0.3, 4.6, 50)
    P = np.stack([np.full_like(ang, 1.5), 0.33 * np.cos(ang), 0.33 * np.sin(ang)], -1)
    sc.line(P, C["motion"], lw=12, split=6, bias=-50)
    sc.arrow_head_on_path(P, len(P) - 3, C["motion"], size=58, bias=-50)
    sc.draw(ax, z0=10)
    # compass + wire, smaller, lower left of the hero area
    cam2 = Camera(target=(0, 0, 0.2), az=-100, el=55, dist=14.5, fov=30, cx=300, cy=168)
    s2 = Scene3D(cam2)
    for faces, col, kind in compass((0, 0, 0), 63 * D):
        s2.mesh(faces, col, cull=not kind.startswith("needle"),
                bias={"housing": 1.0, "dial": 1.0, "pivot": -0.1}.get(kind, 0.0))
    s2.mesh(tube(np.linspace([0, -1.35, 0.7], [0, 1.35, 0.7], 10), 0.06, n=10), C["current"])
    P = np.linspace([0, -1.35, 0.7], [0, 1.35, 0.7], 40)
    for i in (8, 30):
        s2.arrow_head_on_path(P, i, C["current"], size=36, bias=-0.05)
    s2.draw(ax, z0=400)
