#!/usr/bin/env python3
"""e-(p1) e+(p2) -> pi-(p3) pi+(p4) A(k): GoSam virtual minus vacuum polarisation.

The reference has no vacuum polarisation in the virtual photon. GoSam
(alpha(0) scheme) includes the loops of e, mu, tau, phi. At tree level
  M0 = M_ISR (photon virtuality s' = (p3+p4)^2) + M_FSR (virtuality s),
so the VP part of 2 Re(M0^* M1) is
  2 Re sum_spins M0^* [Pi(s') M_ISR + Pi(s) M_FSR],
with Pi the renormalised (Pi(q^2) - Pi(0)) complex self-energy at q^2 + i0.
M_ISR and M_FSR are computed here from the Feynman rules of scalar QED
(explicit spinors, Dirac representation); |M_ISR + M_FSR|^2 is checked
against the Born from GoSam. Spin sums over e-, e+ (as the reference).

Usage: python3 check_vp.py   (builds and runs ./check_ref first)
"""
import cmath
import math
import re
import subprocess

alpha = 0.72973525694968E-02
me = 0.510998950e-3
mmu = 105.658375e-3
mtau = 1.77684          # model default, as in the GoSam run
e = math.sqrt(4 * math.pi * alpha)

p1 = [1.8865000000000001, 0.0, 0.0, 1.8864999307924910]
p2 = [1.8865000000000001, 0.0, 0.0, -1.8864999307924910]
pm = [1.8357703467986879, -1.0085769865944880, -1.1497430849657471, -1.0057018311139998]  # pi-
pp = [1.5860223502746045, 0.83575278548617471, 0.84594586149559481, 1.0401345503233728]  # pi+
k = [0.35120730292670727, 0.17282420110831340, 0.30379722347015226, -3.4432719209372929E-002]
born_ref = 0.16283062859198236
virt_ref = -0.60017574045465627


def dot(a, b):
    return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]


def add(*vs):
    return [sum(c) for c in zip(*vs)]


def sc(x, v):
    return [x * c for c in v]


mpi = math.sqrt(dot(pm, pm))

# Dirac representation
I2 = [[1, 0], [0, 1]]
sig = [[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]]


def block(a, b, c, d):
    return [a[0] + b[0], a[1] + b[1], c[0] + d[0], c[1] + d[1]]


Z2 = [[0, 0], [0, 0]]
neg = lambda m: [[-x for x in r] for r in m]
G = [block(I2, Z2, Z2, neg(I2))] + [block(Z2, s, neg(s), Z2) for s in sig]


def mm(a, b):
    return [[sum(a[i][l] * b[l][j] for l in range(4)) for j in range(4)] for i in range(4)]


def mv(a, v):
    return [sum(a[i][l] * v[l] for l in range(4)) for i in range(4)]


def slash(p):
    # p_mu gamma^mu = p^0 g^0 - p^i g^i
    return [[p[0] * G[0][i][j] - sum(p[m] * G[m][i][j] for m in (1, 2, 3))
             for j in range(4)] for i in range(4)]


def madd(a, b, f=1):
    return [[a[i][j] + f * b[i][j] for j in range(4)] for i in range(4)]


ID = [[1 if i == j else 0 for j in range(4)] for i in range(4)]


def bar(v):
    c = [x.conjugate() for x in v]
    return [c[0], c[1], -c[2], -c[3]]


def ubar_M_u(a, M, b):
    Mb = mv(M, b)
    return sum(x * y for x, y in zip(a, Mb))


def spinor_u(p, m, s):
    chi = [0j] * 4
    chi[s] = 1
    n = math.sqrt(p[0] + m)
    return [x / n for x in mv(madd(slash(p), ID, m), chi)]


def spinor_v(p, m, s):
    """v = (-pslash + m) eta / sqrt(E + m), eta = lower basis vector."""
    eta = [0j] * 4
    eta[2 + s] = 1
    n = math.sqrt(p[0] + m)
    M = [[m * ID[i][j] - slash(p)[i][j] for j in range(4)] for i in range(4)]
    return [x / n for x in mv(M, eta)]


def polarisations(k):
    kv = k[1:]
    kn = math.sqrt(sum(x * x for x in kv))
    n = [x / kn for x in kv]
    a = [1.0, 0.0, 0.0] if abs(n[0]) < 0.9 else [0.0, 1.0, 0.0]
    d = sum(x * y for x, y in zip(a, n))
    e1 = [x - d * y for x, y in zip(a, n)]
    l1 = math.sqrt(sum(x * x for x in e1))
    e1 = [x / l1 for x in e1]
    e2 = [n[1] * e1[2] - n[2] * e1[1], n[2] * e1[0] - n[0] * e1[2], n[0] * e1[1] - n[1] * e1[0]]
    return [[0.0] + e1, [0.0] + e2]   # real, so eps* = eps


def amplitudes(s1, s2, eps):
    """M_ISR, M_FSR (common overall phase dropped)."""
    u = spinor_u(p1, me, s1)
    vb = bar(spinor_v(p2, me, s2))
    i = 1j
    epsl = slash(eps)
    # lepton currents J^mu (contravariant index mu), mu = 0..3
    q1 = add(p1, sc(-1, k))
    q2 = add(k, sc(-1, p2))
    S1 = [[i * x / (dot(q1, q1) - me ** 2) for x in r] for r in madd(slash(q1), ID, me)]
    S2 = [[i * x / (dot(q2, q2) - me ** 2) for x in r] for r in madd(slash(q2), ID, me)]
    Jisr, Jfsr = [], []
    for mu in range(4):
        g = G[mu]
        A = mm(mm(g, S1), epsl)
        B = mm(mm(epsl, S2), g)
        # vertices i e gamma (electron charge -1: -i e Q gamma = + i e gamma)
        Jisr.append((i * e) ** 2 * (ubar_M_u(vb, A, u) + ubar_M_u(vb, B, u)))
        Jfsr.append((i * e) * ubar_M_u(vb, g, u))
    # contract J^mu with -i g_{mu nu}/q^2 H^nu  -> -i (J.H)/q^2
    def cdot(J, H):
        return J[0] * H[0] - J[1] * H[1] - J[2] * H[2] - J[3] * H[3]
    sp = dot(add(pp, pm), add(pp, pm))
    s = dot(add(p1, p2), add(p1, p2))
    # scalar vertex -i e Q (p_in + p_out) along the phi+ flow (phi+ out: pp, phi- out: -pm in)
    Hisr = [(-i * e) * x for x in add(pp, sc(-1, pm))]
    Misr = -i / sp * cdot(Jisr, Hisr)
    ea = dot(add(sc(2, pp), k), eps)
    eb = dot(add(sc(2, pm), k), eps)
    pa = add(pp, k)
    pb = add(pm, k)
    Ha = [(-i * e) ** 2 * i / (dot(pa, pa) - mpi ** 2) * ea * x for x in add(pp, k, sc(-1, pm))]
    Hb = [(-i * e) ** 2 * i / (dot(pb, pb) - mpi ** 2) * (-eb) * x for x in add(pp, sc(-1, pm), sc(-1, k))]
    Hc = [2 * i * e ** 2 * x for x in eps]
    Hfsr = add(Ha, Hb, Hc)
    Mfsr = -i / s * cdot(Jfsr, Hfsr)
    return Misr, Mfsr


def G_func(q2, m):
    """beta ln((beta+1)/(beta-1)) at q2 + i0."""
    r = m * m / q2
    if 4 * r < 1:
        b = math.sqrt(1 - 4 * r)
        return b * complex(math.log((1 + b) ** 2 / (4 * r)), -math.pi)
    b = math.sqrt(4 * r - 1)
    return complex(2 * b * math.atan(1 / b), 0)


def pi_fermion(q2, m):
    r = m * m / q2
    return alpha / (3 * math.pi) * (-5.0 / 3 - 4 * r + (1 + 2 * r) * G_func(q2, m))


def pi_scalar(q2, m):
    r = m * m / q2
    return alpha / (12 * math.pi) * (-8.0 / 3 + 8 * r + (1 - 4 * r) * G_func(q2, m))


def pi_tot(q2):
    return pi_fermion(q2, me) + pi_fermion(q2, mmu) + pi_fermion(q2, mtau) + pi_scalar(q2, mpi)


def main():
    out = subprocess.run(["make", "-s", "run"], capture_output=True, text=True).stdout
    print("\n".join(l for l in out.splitlines() if l.startswith(" ") and ("GoSam" in l or "m_pi" in l)))
    gosam_virt4 = float(re.search(r"4\*VIRT\s+GoSam, ref, rel.diff:\s+(\S+)", out).group(1))
    gosam_born4 = float(re.search(r"4\*BORN\s+GoSam, ref, rel.diff:\s+(\S+)", out).group(1))

    s = dot(add(p1, p2), add(p1, p2))
    sp = dot(add(pp, pm), add(pp, pm))
    Ps, Psp = pi_tot(s), pi_tot(sp)
    born = 0.0
    vp = 0.0
    for s1 in (0, 1):
        for s2 in (0, 1):
            for eps in polarisations(k):
                Mi, Mf = amplitudes(s1, s2, eps)
                M0 = Mi + Mf
                born += abs(M0) ** 2
                vp += 2 * (M0.conjugate() * (Psp * Mi + Ps * Mf)).real
    print()
    print(" s, s'                             : %.15e %.15e" % (s, sp))
    print(" Born (this script, spin-summed)   : %.15e  GoSam*4: %.15e  rel %.2e"
          % (born, gosam_born4, born / gosam_born4 - 1))
    print(" VP part, spin-summed              : %.15e" % vp)
    corr = gosam_virt4 - vp
    print(" 4*VIRT GoSam - VP                 : %.15e  ref: %.15e  rel %.2e"
          % (corr, virt_ref, corr / virt_ref - 1))


if __name__ == "__main__":
    main()
