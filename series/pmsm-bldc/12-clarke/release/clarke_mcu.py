# Lesson 12: Clarke and inverse Clarke as they run once per PWM period (amplitude-invariant).
INV_SQRT3 = 0.5773502691896258   # 1/sqrt(3), computed offline
SQRT3_2 = 0.8660254037844386     # sqrt(3)/2

def clarke(ia, ib):              # i_c = -ia - ib assumed
    return ia, (ia + 2.0*ib) * INV_SQRT3    # 2 mul, 1 add

def inv_clarke(v_al, v_be):      # for the PWM stage
    h, s = 0.5*v_al, SQRT3_2*v_be           # 2 mul
    return v_al, s - h, -s - h              # 2 add/sub

def clarke_q15(ia, ib):          # integers in and out
    return ia, ((ia + 2*ib)*18919 + 16384) >> 15
