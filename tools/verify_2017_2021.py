"""Independent numerical and exact checks for the 2017 and 2021 chapters.

Standard library only. No source scans or personal information are needed.
Finite differences and direct state counting check the published expansions;
exact rational arithmetic checks the crystallographic integer-lattice claims.
"""
from fractions import Fraction as F
from itertools import product
from math import cos, sin, sqrt, pi, hypot, log, isclose
import cmath
import json

def dot(u, v):
    return sum(x*y for x, y in zip(u, v))

def cross(u, v):
    return (u[1]*v[2]-u[2]*v[1],
            u[2]*v[0]-u[0]*v[2],
            u[0]*v[1]-u[1]*v[0])

def lin(coeffs, basis):
    return tuple(sum(c*v[j] for c, v in zip(coeffs, basis)) for j in range(3))

def determinant(cols):
    return dot(cols[0], cross(cols[1], cols[2]))

def close(actual, expected, tol=1e-8):
    assert abs(actual-expected) <= tol, (actual, expected, tol)

t = [(F(0), F(1,2), F(1,2)),
     (F(1,2), F(0), F(1,2)),
     (F(1,2), F(1,2), F(0))]
g = [(-1,1,1), (1,-1,1), (1,1,-1)]  # reciprocal vectors / (2*pi)
assert determinant(t) == F(1,4)
for i, j in product(range(3), repeat=2):
    assert dot(g[i], t[j]) == int(i == j)

m = (1,2,-2)
normal = lin(m, g)
assert normal == (-1,-3,5)
u, v = lin((-2,1,0), t), lin((2,0,1), t)
assert cross(u,v) == tuple(F(x,4) for x in normal)
assert dot(cross(u,v), cross(u,v)) == F(35,16)
for n in product(range(-5,6), repeat=3):
    if dot(n,m) == 0:
        assert lin(n,t) == tuple(n[1]*u[i]+n[2]*v[i] for i in range(3))
d = 1/sqrt(35)
close(sqrt(35)/4*d, 1/4)

normal_c = lin((1,2,-2),t)
assert normal_c == (0,F(-1,2),F(3,2))
recip_c = lin((2,3,-1),g)
assert recip_c == (0,-2,6)
assert cross(normal_c,recip_c) == (0,0,0)
# Integer projections include 0 and 1, so no occupied intermediate plane is missed.
levels = {dot(n,(2,3,-1)) for n in product(range(-3,4),repeat=3)}
assert {0,1}.issubset(levels)
close(1/sqrt(dot(recip_c,recip_c)),1/(2*sqrt(10)))

u111,v111=lin((-1,1,0),t),lin((-1,0,1),t)
assert dot(u111,u111)==F(1,2) and dot(v111,v111)==F(1,2)
assert dot(u111,v111)==F(1,4)
assert dot(cross(u111,v111),cross(u111,v111))==F(3,16)
close(sqrt(3)/4 * 1/sqrt(3),1/4)

T2=((2,-1),(-1,1))
S2=((1,3),(-2,-1))
M=((4,7),(-3,-4))  # columns of the basis-change matrix
for source, coeff in zip(S2,M):
    assert tuple(sum(T2[k][j]*coeff[k] for k in range(2)) for j in range(2)) == source
assert M[0][0]*M[1][1]-M[1][0]*M[0][1]==5

def h(x,y):
    return 1+4*cos(x/2)*cos(sqrt(3)*y/2)+4*cos(x/2)**2

def band_abs(x,y):
    # Compute the magnitude before squaring, avoiding cancellation near a Dirac zero.
    return abs(1+2*cos(x/2)*cmath.exp(1j*sqrt(3)*y/2))

for i,j in product(range(-7,8),repeat=2):
    close(h(i*.37,j*.29),band_abs(i*.37,j*.29)**2,1e-13)

def hessian(func,point,step=3e-4):
    x,y=point
    f=func(x,y)
    xx=(func(x+step,y)-2*f+func(x-step,y))/step**2
    yy=(func(x,y+step)-2*f+func(x,y-step))/step**2
    xy=(func(x+step,y+step)-func(x+step,y-step)
        -func(x-step,y+step)+func(x-step,y-step))/(4*step**2)
    return xx,xy,yy

points={"Gamma":(0,0),"Q":(pi,pi/sqrt(3)),"P":(4*pi/3,0)}
targets={"Gamma":(-3,0,-3),"Q":(2,sqrt(3),0),"P":(1.5,0,1.5)}
hessians={}
for label,p in points.items():
    hs=hessian(h,p)
    for value,target in zip(hs,targets[label]):
        close(value,target,2e-6)
    hessians[label]=hs

# Check the principal directions against finite differences of the energy itself.
q0=points["Q"]
pu=(sqrt(3)/2,.5)
pv=(-.5,sqrt(3)/2)
def directional_curvature(vec,eps=5e-4):
    return (band_abs(q0[0]+eps*vec[0],q0[1]+eps*vec[1])
          +band_abs(q0[0]-eps*vec[0],q0[1]-eps*vec[1])-2)/eps**2
close(directional_curvature(pu),1.5,1e-6)
close(directional_curvature(pv),-.5,1e-6)

dirac_slopes=[]
for angle in [0,.4,1.2,2.1,3.7,5.4]:
    eps=1e-5
    slope=band_abs(points["P"][0]+eps*cos(angle),eps*sin(angle))/eps
    close(slope,sqrt(3)/2,4e-6)
    dirac_slopes.append(slope)

# Count states in the FULL dispersion near Gamma; compare N(deltaE)/deltaE
# with the predicted DOS 1/pi for a=gamma=1. A midpoint grid is deterministic.
side=.13
steps=800
dx=2*side/steps
threshold=.002
count=0
for ix in range(steps):
    x=-side+(ix+.5)*dx
    for iy in range(steps):
        y=-side+(iy+.5)*dx
        if 3-band_abs(x,y) <= threshold:
            count+=1
dos_gamma=count*dx*dx/(2*pi)**2/threshold
assert abs(dos_gamma-1/pi)/(1/pi)<.006,(dos_gamma,1/pi)

# Independent saddle integration: eliminate one quadratic coordinate, integrate
# the remaining delta-function Jacobian with a nonsingular hyperbolic change.
# Compare its exact integral to numerical Simpson integration in the new coordinate.
def simpson(func,lo,hi,n=2000):
    step=(hi-lo)/n
    total=func(lo)+func(hi)
    total+=4*sum(func(lo+(2*i-1)*step) for i in range(1,n//2+1))
    total+=2*sum(func(lo+2*i*step) for i in range(1,n//2))
    return total*step/3

alpha,beta=.75,.25
# For rectangular |X|,|Y|<L and e>0, integration over X gives
# I = 2*asinh(sqrt(L^2-e)/sqrt(e)); its log coefficient is 1.
L=.1
saddle_checks=[]
from math import asinh
for e in [1e-4,1e-5,1e-6]:
    ymax=sqrt(L*L-e)
    integ=simpson(lambda y:2/sqrt(e+y*y),0,ymax,8000)
    expected=2*asinh(ymax/sqrt(e))
    close(integ,expected,1e-8)
    coeff=integ/(4*pi*pi*sqrt(alpha*beta))
    saddle_checks.append({"energy":e,"dos":coeff})
# The ratio of differences removes the cutoff-dependent constant.
log_coefficient=(saddle_checks[-1]["dos"]-saddle_checks[-2]["dos"])/log(10)
close(log_coefficient,1/(pi*pi*sqrt(3)),3e-5)

c0=(9*pi/4)**(1/3)
def delta_hf(rs):
    kf=c0/rs
    return 3*kf*kf/10*(2**(2/3)-1)-3*kf/(4*pi)*(2**(1/3)-1)
rs_star=2*pi*c0/5*(2**(1/3)+1)
close(delta_hf(rs_star),0,1e-15)
assert delta_hf(2)>0 and delta_hf(10)<0

report={
    "status":"PASS",
    "fcc_plane_spacings_a_units":[d,1/(2*sqrt(10))],
    "fcc_planar_density_times_a_squared":4/sqrt(35),
    "h_finite_difference_hessians":hessians,
    "dirac_directional_slopes_gamma_a_units":dirac_slopes,
    "gamma_dos_full_dispersion_count":dos_gamma,
    "gamma_dos_prediction":1/pi,
    "saddle_log_coefficient_numerical":log_coefficient,
    "saddle_log_coefficient_prediction":1/(pi*pi*sqrt(3)),
    "uniform_hf_branch_crossing_rs":rs_star,
    "scope":"Lattice and band/DOS arithmetic only; this is not verification of official marking."
}
print(json.dumps(report,indent=2))
