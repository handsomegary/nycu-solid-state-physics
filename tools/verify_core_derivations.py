"""Independent numerical checks of the 2020 and undated solution derivations.

Uses only the Python standard library. Finite differences, lattice enumeration,
and direct Hamiltonian action check formulas derived in the written solutions.
These checks do not certify every explanatory claim in the books.
"""
import cmath
import math


def close(actual, expected, tolerance=1e-6):
    assert abs(actual - expected) <= tolerance * max(1, abs(expected)), (actual, expected)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def norm(a):
    return math.sqrt(dot(a, a))


def graphene(x, y):
    # Complex structure factor avoids subtractive cancellation at a Dirac zero.
    return abs(1 + 2 * math.cos(x / 2) * cmath.exp(1j * math.sqrt(3) * y / 2))


def hessian(f, x, y, h=2e-4):
    xx = (f(x+h, y)-2*f(x, y)+f(x-h, y))/h**2
    yy = (f(x, y+h)-2*f(x, y)+f(x, y-h))/h**2
    xy = (f(x+h, y+h)-f(x+h, y-h)-f(x-h, y+h)+f(x-h, y-h))/(4*h*h)
    return xx, yy, xy


def check_graphene():
    for actual, expected in zip(hessian(graphene, 0, 0), (-0.5, -0.5, 0)):
        close(actual, expected)
    xx, yy, xy = hessian(graphene, math.pi, math.pi/math.sqrt(3))
    for actual, expected in zip((xx, yy, xy), (1, 0, math.sqrt(3)/2)):
        close(actual, expected)
    delta = math.sqrt((xx-yy)**2 + 4*xy*xy)
    close((xx+yy+delta)/2, 1.5)
    close((xx+yy-delta)/2, -0.5)
    for angle in [j * math.pi / 13 for j in range(26)]:
        q = 1e-5
        slope = graphene(4*math.pi/3+q*math.cos(angle), q*math.sin(angle))/q
        close(slope, math.sqrt(3)/2, 3e-6)

    # Count exact-dispersion states near the lower-band Gamma edge by the area
    # enclosed by a radial contour, then differentiate that area numerically.
    def area(energy):
        radii_squared = []
        for j in range(120):
            angle = 2*math.pi*(j+0.5)/120
            lo, hi = 0.0, 0.2
            for _ in range(45):
                radius = (lo+hi)/2
                excess = 3-graphene(radius*math.cos(angle), radius*math.sin(angle))
                if excess < energy:
                    lo = radius
                else:
                    hi = radius
            radii_squared.append(((lo+hi)/2)**2)
        return math.pi * sum(radii_squared)/len(radii_squared)
    dos_spin_two = 2/(2*math.pi)**2*(area(2e-5)-area(1e-5))/1e-5
    close(dos_spin_two, 2/math.pi, 2e-5)


def check_fcc():
    t = [(0,.5,.5), (.5,0,.5), (.5,.5,0)]
    g = [(-1,1,1), (1,-1,1), (1,1,-1)]  # in units of 2*pi/a
    close(dot(t[0], cross(t[1], t[2])), .25)
    for i in range(3):
        for j in range(3):
            close(dot(t[i], g[j]), int(i == j))
    for normal, expected in [((-1,-3,5), math.sqrt(35)), ((0,-1,3), math.sqrt(40))]:
        vectors = []
        for i in range(-5,6):
            for j in range(-5,6):
                for k in range(-5,6):
                    vector = tuple(i*g[0][p]+j*g[1][p]+k*g[2][p] for p in range(3))
                    if norm(vector) and norm(cross(normal, vector)) < 1e-12:
                        vectors.append(norm(vector))
        close(min(vectors), expected)
    u = tuple(-2*t[0][i]+t[1][i] for i in range(3))
    v = tuple(2*t[0][i]+t[2][i] for i in range(3))
    close(norm(cross(u, v)), math.sqrt(35)/4)


def check_chain_and_gap():
    # Apply a finite periodic-chain Hamiltonian to every allowed Fourier mode.
    count, onsite, hopping = 17, 0.7, -0.31
    for j in range(count):
        k = 2*math.pi*j/count
        state = [cmath.exp(1j*k*n)/math.sqrt(count) for n in range(count)]
        predicted = onsite + 2*hopping*math.cos(k)
        for n in range(count):
            action = onsite*state[n]+hopping*(state[(n-1)%count]+state[(n+1)%count])
            close(action, predicted*state[n], 1e-12)
    # Direct two-by-two eigenvalue formula versus the lower-band mass expansion.
    reciprocal, coupling = 2*math.pi, 0.2
    def lower(q):
        left, right = (reciprocal/2+q)**2/2, (-reciprocal/2+q)**2/2
        return (left+right)/2-math.sqrt(((left-right)/2)**2+coupling**2)
    h = 1e-5
    curvature = (lower(h)-2*lower(0)+lower(-h))/h**2
    close(curvature, 1-reciprocal**2/(4*coupling), 2e-6)


def check_density_and_projector():
    # Number conservation independently fixes the polarized Fermi sphere.
    normal_k = 1.3
    density = 2 * (4*math.pi/3)*normal_k**3/(2*math.pi)**3
    polarized_k = (density*(2*math.pi)**3/(4*math.pi/3))**(1/3)
    close(polarized_k/normal_k, 2**(1/3))
    for zeta in (0, .2, .7, 1):
        ks = [normal_k*(1+zeta)**(1/3), normal_k*(1-zeta)**(1/3)]
        number = sum(k**3/(6*math.pi**2) for k in ks)
        kinetic = sum(k**5/(20*math.pi**2) for k in ks)
        close(number, density)
        expected = 3*normal_k**2/10*((1+zeta)**(5/3)+(1-zeta)**(5/3))/2
        close(kinetic/number, expected)
    for n in (.03, .2, 1, 4):
        h = n*1e-5
        density_energy = lambda value: -value**(4/3)
        derivative = (density_energy(n+h)-density_energy(n-h))/(2*h)
        close(derivative, -(4/3)*n**(1/3), 1e-9)

    # Projection must remove the core, preserve a valence component, and be
    # idempotent. The Phillips-Kleinman identity is checked on a nontrivial
    # pseudo-state containing both core and valence components.
    core = [1/math.sqrt(2), 1/math.sqrt(2), 0]
    valence = [0, 0, 1]
    pseudo = [0.7*core[i]+valence[i] for i in range(3)]
    project = lambda state: [state[i]-core[i]*dot(core,state) for i in range(3)]
    true_state = project(pseudo)
    close(dot(core, true_state), 0)
    for a, b in zip(project(true_state), true_state):
        close(a, b)
    core_energy, valence_energy = -5, 2
    for i in range(3):
        h_pseudo = core_energy*.7*core[i]+valence_energy*valence[i]
        repulsion = (valence_energy-core_energy)*core[i]*dot(core,pseudo)
        close(h_pseudo+repulsion, valence_energy*pseudo[i])


def main():
    for check in (check_graphene, check_fcc, check_chain_and_gap, check_density_and_projector):
        check()
        print('PASS', check.__name__)


if __name__ == '__main__':
    main()
