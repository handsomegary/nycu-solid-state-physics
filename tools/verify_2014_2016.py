"""Independent numerical checks for the 2014 and 2016 theory solutions.

Uses only Python's standard library. The checks include integer enumeration
of reciprocal directions, finite differences of the exact honeycomb band,
state counting in exact-band energy windows, and a finite-dimensional OPW
projection example. This is not a proof that every prose statement is right.
"""
from fractions import Fraction as F
from itertools import product
from math import cos, sin, hypot, sqrt, pi, isclose


def dot(x, y):
    return sum(a*b for a, b in zip(x, y))


def cross(x, y):
    return (x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2],
            x[0]*y[1]-x[1]*y[0])


def combine(coefficients, vectors):
    return tuple(sum(c*v[j] for c, v in zip(coefficients, vectors))
                 for j in range(len(vectors[0])))


def close(actual, expected, tol=1e-7):
    assert isclose(actual, expected, rel_tol=tol, abs_tol=tol), (actual, expected)


def lattice_checks():
    t = [(0, F(1, 2), F(1, 2)), (F(1, 2), 0, F(1, 2)),
         (F(1, 2), F(1, 2), 0)]
    # Reciprocal coordinates are in units of 2*pi/a, with a=1.
    g = [(-1, 1, 1), (1, -1, 1), (1, 1, -1)]
    for i, j in product(range(3), repeat=2):
        assert dot(t[i], g[j]) == int(i == j)
    volume = dot(t[0], cross(t[1], t[2]))
    assert volume == F(1, 4)
    for direction, expected_squared in [((-1, -3, 5), 35), ((0, -1, 3), 40)]:
        candidates = []
        for h in product(range(-7, 8), repeat=3):
            v = combine(h, g)
            if v != (0, 0, 0) and cross(v, direction) == (0, 0, 0):
                candidates.append(dot(v, v))
        assert min(candidates) == expected_squared
    u = combine((-2, 1, 0), t)
    v = combine((2, 0, 1), t)
    area_squared = dot(cross(u, v), cross(u, v))
    assert area_squared == F(35, 16)
    close(sqrt(float(area_squared)) / sqrt(35), float(volume))
    # Enumerating occupied-plane indices tests the common-factor issue.
    indices = {n1 + 2*n2 - 2*n3 for n1, n2, n3 in product(range(-3, 4), repeat=3)}
    assert set(range(-10, 11)).issubset(indices)
    print('PASS FCC reciprocal basis, shortest directions, occupied planes, surface cell')


def band(kx, ky):
    # Stable complex-amplitude evaluation avoids cancellation near a Dirac zero.
    x, y = kx/2, sqrt(3)*ky/2
    return hypot(1 + 2*cos(x)*cos(y), 2*cos(x)*sin(y))


def hessian(x, y, h=2e-4):
    center = band(x, y)
    xx = (band(x+h, y)-2*center+band(x-h, y))/h**2
    yy = (band(x, y+h)-2*center+band(x, y-h))/h**2
    xy = (band(x+h, y+h)-band(x+h, y-h)-band(x-h, y+h)+band(x-h, y-h))/(4*h**2)
    return xx, yy, xy


def band_checks():
    close(band(0, 0), 3)
    close(band(pi, pi/sqrt(3)), 1)
    close(band(4*pi/3, 0), 0)
    for actual, expected in zip(hessian(0, 0), (-0.5, -0.5, 0)):
        close(actual, expected, 5e-7)
    xx, yy, xy = hessian(pi, pi/sqrt(3))
    for actual, expected in zip((xx, yy, xy), (1, 0, sqrt(3)/2)):
        close(actual, expected, 5e-7)
    discriminant = sqrt((xx-yy)**2 + 4*xy**2)
    close((xx+yy+discriminant)/2, 1.5, 5e-7)
    close((xx+yy-discriminant)/2, -0.5, 5e-7)
    h = 1e-5
    for angle in (0, 0.31, 1.2, 2.9):
        for sign in (-1, 1):
            slope = band(4*pi/3 + sign*h*cos(angle), sign*h*sin(angle))/h
            close(slope, sqrt(3)/2, 5e-6)
    # Opposite one-sided x derivatives are unequal, not zero or infinite.
    right = band(4*pi/3+h, 0)/h
    left = -band(4*pi/3-h, 0)/h
    assert right > 0.86 and left < -0.86
    print('PASS exact-band finite-difference Hessians and finite Dirac directional slopes')


def count_window(center, half_width, energy_window, energy_from_band, grid=600):
    step = 2*half_width/grid
    count = 0
    lo, hi = energy_window
    for i in range(grid):
        x = center[0]-half_width+(i+0.5)*step
        for j in range(grid):
            y = center[1]-half_width+(j+0.5)*step
            energy = energy_from_band(band(x, y))
            count += lo <= energy < hi
    return count*step**2/((2*pi)**2*(hi-lo))


def dos_checks():
    gamma_dos = count_window((0, 0), 0.16, (0.001, 0.003), lambda f: 3-f)
    close(gamma_dos, 1/pi, 0.012)
    dirac_dos = count_window((4*pi/3, 0), 0.08, (0.02, 0.04), lambda f: f)
    expected = (0.02+0.04)/(4*pi*(3/4))
    close(dirac_dos, expected, 0.00015)
    print(f'PASS direct exact-band state counts: Gamma DOS={gamma_dos:.7f}; Dirac DOS={dirac_dos:.7f}')


def opw_and_spin_checks():
    # A rotated core eigenvector tests an actual non-diagonal projection.
    c = (3/5, 4/5, 0)
    v = (-4/5, 3/5, 0)
    w = (0, 0, 1)
    eigensystem = [(-20, c), (-2, v), (4, w)]
    def apply_h(x):
        return combine([energy*dot(state, x) for energy, state in eigensystem],
                       [state for _, state in eigensystem])
    phi = combine((1, 0.7), (v, c))
    psi = tuple(phi[i]-c[i]*dot(c, phi) for i in range(3))
    for actual, expected in zip(psi, v):
        close(actual, expected)
    close(dot(c, psi), 0)
    energy = -2
    repulsive = tuple((energy+20)*c[i]*dot(c, phi) for i in range(3))
    effective_phi = tuple(a+b for a, b in zip(apply_h(phi), repulsive))
    for actual, expected in zip(effective_phi, (energy*x for x in phi)):
        close(actual, expected)
    assert dot(phi, repulsive) > 0
    # Without the projection a valence/core mixture lowers the Rayleigh energy.
    mixed = combine((sqrt(0.75), 0.5), (v, c))
    assert dot(mixed, apply_h(mixed)) < energy
    for power, ratio in [(5/3, 2**(2/3)), (4/3, 2**(1/3))]:
        close((2**power+0**power)/2, ratio)
    print('PASS OPW projection direction, effective Hamiltonian, collapse example, spin factors')


if __name__ == '__main__':
    lattice_checks()
    band_checks()
    dos_checks()
    opw_and_spin_checks()
