"""Reusable 3D models: motor cutaway, compass, coil."""
import numpy as np
from engine import annular_sector, cylinder, box, tube, C, normalize

D = np.pi / 180
MOTOR_FRAME = np.eye(3)  # columns (u, v, w) = (x, y, z): motor axis along world x


def _check_outward(faces, center):
    bad = 0
    for f in faces:
        n = np.cross(f[1] - f[0], f[2] - f[0])
        if np.linalg.norm(n) < 1e-12:
            continue
        if n @ (f.mean(0) - center) < -1e-9:
            bad += 1
    return bad


def motor_parts(rotor_angle=0.0, cut=(15 * D, 105 * D), L=1.2, slots=12, poles=8):
    """Return list of (faces, color, kind) for a 12-slot / 8-pole surface-PM motor with a quadrant cut away."""
    parts = []
    hl = L / 2
    pitch = 2 * np.pi / slots

    def in_cut(a0, a1):
        c = 0.5 * (a0 + a1) % (2 * np.pi)
        return cut[0] - 1e-6 <= c <= cut[1] + 1e-6

    for k in range(slots):
        ac = k * pitch
        # yoke sector centred on tooth k
        y0, y1 = ac - pitch / 2, ac + pitch / 2
        if not in_cut(y0, y1):
            nb_lo = in_cut(y0 - pitch, y1 - pitch); nb_hi = in_cut(y0 + pitch, y1 + pitch)
            for f, kd in annular_sector(0.84, 1.0, y0, y1, -hl, hl, n=4, frame=MOTOR_FRAME, nz=3,
                                        sides=(True, True, nb_lo, nb_hi), tag=True):
                parts.append(([f], C["steel_cut"] if kd == "side" else C["steel"], "stator"))
            # tooth body and shoe
            tb = 6.2 * D
            for f, kd in annular_sector(0.64, 0.84, ac - tb, ac + tb, -hl, hl, n=1, frame=MOTOR_FRAME, nz=3,
                                        sides=(False, True, True, True), tag=True):
                parts.append(([f], C["steel"], "stator"))
            for f, kd in annular_sector(0.60, 0.64, ac - 13 * D, ac + 13 * D, -hl, hl, n=3, frame=MOTOR_FRAME, nz=3,
                                        tag=True):
                parts.append(([f], C["steel"], "stator"))
            # concentrated coil around tooth k: two coil sides + two end turns
            for s0, s1 in ((ac + 6.6 * D, ac + 14.2 * D), (ac - 14.2 * D, ac - 6.6 * D)):
                parts.append((annular_sector(0.655, 0.825, s0, s1, -hl, hl, n=1, frame=MOTOR_FRAME, nz=3),
                              C["copper"], "coil"))
            for z0, z1 in ((hl, hl + 0.11), (-hl - 0.11, -hl)):
                parts.append((annular_sector(0.655, 0.825, ac - 14.2 * D, ac + 14.2 * D, z0, z1, n=3,
                                             frame=MOTOR_FRAME), C["copper"], "coil"))
    # rotor: magnets + core + shaft
    pp = 2 * np.pi / poles
    for k in range(poles):
        a0 = rotor_angle + k * pp
        col = C["N"] if k % 2 == 0 else C["S"]
        parts.append((annular_sector(0.50, 0.565, a0 - 20 * D, a0 + 20 * D, -hl, hl, n=3, frame=MOTOR_FRAME, nz=3),
                      col, "magnet"))
        parts.append((annular_sector(0.16, 0.50, a0 - pp / 2, a0 + pp / 2, -hl, hl, n=2, frame=MOTOR_FRAME, nz=2,
                                     sides=(True, False, False, False)), C["core"], "core"))
    for s in (-1, 1):  # shaft stubs, subdivided for depth sorting
        zs = np.linspace(s * hl, s * 1.45, 4)
        for z0, z1 in zip(zs[:-1], zs[1:]):
            a_, b_ = (z0, z1) if s > 0 else (z1, z0)
            parts.append((cylinder([a_, 0, 0], [b_, 0, 0], 0.16, n=28, caps=(abs(z1) > 1.4 or abs(z0) > 1.4)),
                          C["shaft"], "shaft"))
    return parts


def compass(center, needle_angle, r=1.0, h=0.14):
    """Compass: housing ring, dial, needle. needle_angle measured from +y (north) toward -x (west) is positive."""
    c = np.asarray(center, float)
    parts = []
    parts.append((cylinder(c + [0, 0, 0], c + [0, 0, h], r, n=48), "#2A3240", "housing"))
    parts.append((cylinder(c + [0, 0, h], c + [0, 0, h + 0.002], 0.94 * r, n=48), "#D9DEE5", "dial"))
    # needle: flat diamond, north half red, south half light
    th = needle_angle
    fwd = np.array([-np.sin(th), np.cos(th), 0.0])
    side = np.array([np.cos(th), np.sin(th), 0.0])
    zc = c + [0, 0, h + 0.03]
    L, wd = 0.82 * r, 0.10 * r
    tipN, tipS = zc + fwd * L, zc - fwd * L
    left, right = zc + side * wd, zc - side * wd
    up = np.array([0, 0, 0.012])
    parts.append(([np.array([tipN + up, left + up, right + up])], C["N"], "needleN"))
    parts.append(([np.array([tipS + up, right + up, left + up])], "#F2F4F7", "needleS"))
    parts.append((cylinder(zc - [0, 0, 0.02], zc + [0, 0, 0.03], 0.04, n=12), "#9AA3AF", "pivot"))
    return parts


def rect_coil(center, theta, w=1.6, h=1.2, r=0.05):
    """Rectangular coil in a vertical plane; normal n=(cosθ, sinθ, 0). Returns (faces, path corners in current order)."""
    c = np.asarray(center, float)
    t = np.array([-np.sin(theta), np.cos(theta), 0.0]); z = np.array([0, 0, 1.0])
    P = [c - w / 2 * t - h / 2 * z, c + w / 2 * t - h / 2 * z, c + w / 2 * t + h / 2 * z, c - w / 2 * t + h / 2 * z]
    faces = []
    for i in range(4):
        a, b = P[i], P[(i + 1) % 4]
        d = normalize(b - a)
        faces += cylinder(a - d * r, b + d * r, r, n=10)
    return faces, P


if __name__ == "__main__":
    for f, col, kind in motor_parts()[:5]:
        pass
    # outward-normal sanity on one annular sector
    fs = annular_sector(0.5, 1.0, 0, 0.5, -0.3, 0.3, n=2, frame=MOTOR_FRAME, nz=2)
    cen = np.mean([f.mean(0) for f in fs], 0)
    print("inward faces (expect 0 for convex-ish check on outer/end/side):", _check_outward(fs, cen))
