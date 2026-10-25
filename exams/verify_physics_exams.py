"""Numerically check selected derivations in the three Physics exam solutions.

Run with Python 3.10 or newer; only the standard library is required.
Checks use lattice enumeration, quadrature of the original Fermi distribution,
and independently assembled harmonic matrices. They are mathematical checks,
not a claim to validate every conceptual statement in the written solutions.
"""
from fractions import Fraction
from itertools import product
import cmath
import json
import math


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def determinant(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def matvec(m, v):
    return [dot(row, v) for row in m]


def close(actual, expected, tolerance=1e-10):
    assert abs(actual - expected) <= tolerance, (actual, expected, tolerance)


def simpson(function, lo, hi, intervals=6000):
    assert intervals % 2 == 0
    step = (hi - lo) / intervals
    terms = [function(lo), function(hi)]
    terms.extend((4 if i % 2 else 2) * function(lo + step * i)
                 for i in range(1, intervals))
    return math.fsum(terms) * step / 3


def fermi(energy, chemical_potential, temperature):
    y = (energy - chemical_potential) / temperature
    if y > 40:
        return math.exp(-y)
    if y < -40:
        return 1.0
    return 1.0 / (1.0 + math.exp(y))


def thermal_state(temperature):
    """Units: E_F = k_B = N = 1; integrate over x = k/k_F."""
    def count(mu):
        return 3 * simpson(lambda x: x*x * fermi(x*x, mu, temperature), 0, 3)
    lo, hi = 0.5, 1.2
    assert count(lo) < 1 < count(hi)
    for _ in range(44):
        mid = (lo + hi) / 2
        if count(mid) < 1:
            lo = mid
        else:
            hi = mid
    mu = (lo + hi) / 2
    energy = 3 * simpson(lambda x: x**4 * fermi(x*x, mu, temperature), 0, 3)
    return mu, energy, count(mu)


def check_lattices():
    # Columns of T are direct vectors in units a/2; reciprocal vectors use 2*pi/a.
    t = [(0, 1, 1), (1, 0, 1), (1, 1, 0)]
    g = [(-1, 1, 1), (1, -1, 1), (1, 1, -1)]
    for i in range(3):
        for j in range(3):
            assert dot(t[i], g[j]) == (2 if i == j else 0)
    close(abs(determinant(t)) / 8, 1/4)
    assert abs(determinant(g)) == 4
    # Enumerate fcc points to detect occupied layers normal to [011].
    points = [tuple(sum(n[j] * t[j][i] for j in range(3)) for i in range(3))
              for n in product(range(-4, 5), repeat=3)]
    projections = {point[1] + point[2] for point in points}
    assert set(range(-8, 9)).issubset(projections)
    # In the central layer, x is even and z=-y: u=(2,0,0), v=(0,1,-1).
    for x, y, z in points:
        if y + z == 0:
            assert x % 2 == 0
    layer_spacing = 1/(2*math.sqrt(2))
    area = 1/math.sqrt(2)
    close(area * layer_spacing, 1/4)
    # Independent nearest-vector enumeration in reciprocal fcc coordinates.
    bcc_reciprocal = [h for h in product(range(-3, 4), repeat=3)
                      if h != (0, 0, 0) and sum(h) % 2 == 0]
    shortest_norm2 = min(dot(h, h) for h in bcc_reciprocal)
    nearest = [h for h in bcc_reciprocal if dot(h, h) == shortest_norm2]
    assert shortest_norm2 == 2 and len(nearest) == 12
    # Equal volumes alone admit a false noninteger candidate.
    volume_candidates, lattice_candidates = [], []
    for twice_alpha in range(-8, 9):
        alpha = Fraction(twice_alpha, 2)
        matrix = [[2, -1, 1], [-1, 1, alpha], [2, alpha, 1]]
        if abs(determinant(matrix)) == 1:
            volume_candidates.append(str(alpha))
            if alpha.denominator == 1:
                lattice_candidates.append(str(alpha))
    assert volume_candidates == ['-3/2', '0']
    assert lattice_candidates == ['0']
    matrix = [[2, -1, 1], [-1, 1, 0], [2, 0, 1]]
    inverse_columns = [(-1, -1, 2), (-1, 0, 2), (1, 1, -1)]
    for i, column in enumerate(inverse_columns):
        assert matvec(matrix, column) == [int(i == j) for j in range(3)]
    return {'fcc_primitive_volume_a3': 1/4,
            'fcc_occupied_011_layer_spacing_a': layer_spacing,
            'bcc_brillouin_zone_facets': len(nearest),
            'volume_only_alpha_candidates': volume_candidates,
            'same_lattice_alpha_candidates': lattice_candidates}


def check_bloch_and_sum():
    # Test a localized-orbital Bloch sum far from the truncated sum boundaries.
    k, position, shift = 0.71, 0.23, 3
    def bloch(x):
        return sum(cmath.exp(1j*k*n) * math.exp(-2*(x-n)**2)
                   for n in range(-30, 31))
    close(abs(bloch(position+shift)-cmath.exp(1j*k*shift)*bloch(position)), 0)
    terms = [x*y*z*math.exp(-x*x/2-y*y/3-z*z/5)
             for x, y, z in product(range(-8, 9), repeat=3)]
    close(math.fsum(terms), 0, 1e-12)
    return {'bloch_translation_residual':
            abs(bloch(position+shift)-cmath.exp(1j*k*shift)*bloch(position)),
            'odd_bvk_sum': math.fsum(terms)}


def check_density_of_states():
    # Anisotropic k-space volume from independent Cartesian voxel counting.
    alpha, beta, gamma = 1.3, 0.7, 2.1
    energy_above_edge = 0.8
    axes = [math.sqrt(value*energy_above_edge) for value in (alpha, beta, gamma)]
    subdivisions = 90
    inside = 0
    for i, j, k in product(range(subdivisions), repeat=3):
        x, y, z = [(2*(index+0.5)/subdivisions-1)*axis
                   for index, axis in zip((i, j, k), axes)]
        if x*x/alpha+y*y/beta+z*z/gamma <= energy_above_edge:
            inside += 1
    numerical_volume = inside * (8*math.prod(axes)/subdivisions**3)
    ellipsoid_volume = 4*math.pi/3*math.prod(axes)
    relative_error = abs(numerical_volume/ellipsoid_volume-1)
    assert relative_error < 0.003
    # After exact radial integration, the saddle DOS is a 1D indicator integral.
    cutoff, intervals = 1.0, 100000
    cusp_checks = []
    for delta in (-0.36, -0.09, 0.04):
        width = 2*cutoff/intervals
        measure = math.fsum(width for i in range(intervals)
                            if delta + (-cutoff+(i+0.5)*width)**2 >= 0)
        expected = 2*(cutoff-math.sqrt(-delta)) if delta < 0 else 2*cutoff
        close(measure, expected, 3e-5)
        cusp_checks.append({'delta': delta, 'indicator_measure': measure})
    # The ground-state mean energy is a radial Fermi-sphere integral.
    mean = simpson(lambda k: k**4, 0, 1) / simpson(lambda k: k*k, 0, 1)
    close(mean, 3/5, 1e-12)
    return {'ellipsoid_voxel_relative_error': relative_error,
            'saddle_integrals': cusp_checks, 'mean_energy_over_fermi_energy': mean}


def check_sommerfeld():
    # No Sommerfeld formula is used to solve for the numerical chemical potential.
    records = []
    for temperature in (0.04, 0.02, 0.01):
        mu, energy, number = thermal_state(temperature)
        step = temperature/100
        _, lower_energy, _ = thermal_state(temperature-step)
        _, upper_energy, _ = thermal_state(temperature+step)
        capacity = (upper_energy-lower_energy)/(2*step)
        records.append({'T_over_TF': temperature, 'N_over_N0': number,
                        'mu_over_EF': mu, 'U_over_NEF': energy,
                        'mu_quadratic_coefficient': (1-mu)/temperature**2,
                        'Cv_over_NkB_T_over_TF': capacity/temperature})
        close(number, 1, 1e-10)
    last = records[-1]
    close(last['mu_quadratic_coefficient'], math.pi**2/12, 0.0003)
    close(last['Cv_over_NkB_T_over_TF'], math.pi**2/2, 0.002)
    assert abs(last['Cv_over_NkB_T_over_TF']-1.5) > 3
    assert abs(last['mu_quadratic_coefficient']-math.pi**2/12) < abs(
        records[0]['mu_quadratic_coefficient']-math.pi**2/12)
    return records


def check_optical_vectors():
    records = []
    for m1, m2 in ((1.0, 1.0), (1.0, 4.0), (7.0, 2.0)):
        # Hessian of (spring/2)*(u1-u2)^2, then mass weighting.
        spring = 2.3
        d = [[spring/m1, -spring/math.sqrt(m1*m2)],
             [-spring/math.sqrt(m1*m2), spring/m2]]
        eigenvalue = spring*(1/m1+1/m2)
        # Solve first matrix row for the amplitude ratio, independently of the answer.
        a1 = 1.0
        a2 = -(d[0][0]-eigenvalue)*a1/d[0][1]
        residual = matvec(d, [a1, a2])
        close(residual[0], eigenvalue*a1)
        close(residual[1], eigenvalue*a2)
        close(math.sqrt(m1)*a1+math.sqrt(m2)*a2, 0)
        ratio = a1*a1/(a1*a2)
        close(ratio, -math.sqrt(m2/m1))
        close(m1*a1/math.sqrt(m1)+m2*a2/math.sqrt(m2), 0)
        records.append({'masses': [m1, m2], 'requested_ratio': ratio})
    return records


def main():
    report = {'lattices': check_lattices(),
              'bloch_translation_and_odd_sum': check_bloch_and_sum(),
              'density_of_states': check_density_of_states(),
              'sommerfeld_from_fermi_integrals': check_sommerfeld(),
              'optical_modes_from_hessian': check_optical_vectors()}
    print(json.dumps({'status': 'PASS', 'checks': report}, indent=2))


if __name__ == '__main__':
    main()
