/* Lesson 19 (PMSM & BLDC motors from first principles): one iteration of the dq current loop, single precision.
   Runs in the ADC interrupt triggered at the PWM counter zero (the middle of 000, Lesson 18), once every TS.
   Reference motor (common/motor_params.py): R = 40 mOhm, L = 0.05 mH, lambda_m = 6.5 mWb, 24 V bus, 20 kHz PWM.
   Gains by pole-zero cancellation for a 1 kHz bandwidth: KP = wc L, Ki = wc R (stored as KI_TS = Ki TS).
   The same arithmetic as CurrentController.update() in series/pmsm-bldc/common/control.py (checked by the lesson's
   physics.py, which compiles this file with a test harness). Vendor-neutral: adc_a, adc_b, rotor_angle(),
   rotor_speed() and svpwm() come from your board support code. */
#include <math.h>

#define TS        50e-6f                     /* PWM period, s (20 kHz) */
#define KP        0.3141593f                 /* wc L, V/A (wc = 2 pi 1 kHz) */
#define KI_TS     0.01256637f                /* wc R TS, V/A per sample */
#define LD        0.05e-3f                   /* H */
#define LQ        0.05e-3f                   /* H */
#define LAMBDA    6.5e-3f                    /* Wb, peak flux linkage per phase */
#define VMAX      13.85641f                  /* V, linear SVPWM limit 24 V / sqrt(3) */
#define INV_SQRT3 0.5773503f
/* Example sensing chain (an assumption, adapt to your board): 5 mOhm shunt x gain 20 = 0.1 V/A, 12-bit ADC over
   3.3 V, amplifier biased at mid-scale; a low-side shunt reads -i, hence the minus sign. */
#define K_ADC     (-3.3f / 4096.0f / 0.1f)  /* A per count */
#define OFF_A     2048.0f
#define OFF_B     2048.0f

extern float adc_a, adc_b;                   /* raw readings of the phase a and b shunts */
extern float rotor_angle(void);              /* electrical angle theta_e, rad */
extern float rotor_speed(void);              /* electrical speed omega_e, rad/s */
extern void svpwm(float v_alpha, float v_beta);   /* Lesson 18: duties -> compare registers (load at next zero) */
float id_ref, iq_ref;                        /* set by the torque or speed loop, A */
float xd, xq;                                /* integrator outputs, V */

void current_loop(void)
{
    float ia = K_ADC * (adc_a - OFF_A), ib = K_ADC * (adc_b - OFF_B);
    float th = rotor_angle(), we = rotor_speed();
    float al = ia, be = (ia + 2.0f * ib) * INV_SQRT3;      /* Clarke */
    float c = cosf(th), s = sinf(th);
    float id = c * al + s * be, iq = c * be - s * al;       /* Park */
    float ed = id_ref - id, eq = iq_ref - iq;
    float ud = KP * ed + xd - we * LQ * iq;                 /* PI + feed-forward */
    float uq = KP * eq + xq + we * (LD * id + LAMBDA);
    float vd = fminf(fmaxf(ud, -VMAX), VMAX);               /* limit: d first */
    float vqm = sqrtf(VMAX * VMAX - vd * vd);
    float vq = fminf(fmaxf(uq, -vqm), vqm);                 /* q gets the rest */
    xd += KI_TS * (ed + (vd - ud) / KP);                    /* back-calculation */
    xq += KI_TS * (eq + (vq - uq) / KP);
    float ta = th + 1.5f * we * TS;                         /* angle when v lands */
    c = cosf(ta); s = sinf(ta);
    svpwm(c * vd - s * vq, s * vd + c * vq);                /* inverse Park */
}
