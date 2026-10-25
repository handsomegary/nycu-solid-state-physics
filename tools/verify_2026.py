"""Independent checks for the May 11, 2026 solid-state theory exam.

Standard library only. The tests check integer charge sharing, converged Ewald
potentials, volume derivatives, charge-transfer pair counting, screening Fourier
transforms, and dispersion relations including the Drude conductivity pole.
No original question PDF or supplied answer page is required.
"""
from itertools import product
from fractions import Fraction
from math import cos, sin, exp, expm1, erfc, sqrt, pi, isclose
import cmath
import json


def close(actual, expected, atol=1e-8, rtol=1e-8):
    assert abs(actual-expected) <= atol+rtol*abs(expected), (actual, expected)


def simpson(func, left, right, intervals=20000):
    assert intervals % 2 == 0
    h=(right-left)/intervals
    total=func(left)+func(right)
    total+=4*sum(func(left+(2*j-1)*h) for j in range(1,intervals//2+1))
    total+=2*sum(func(left+2*j*h) for j in range(1,intervals//2))
    return total*h/3


# Q1: enumerate all charges of the stated Evjen cell, rather than hardcoding
# the face, edge and corner counts.
shells={1:Fraction(0),2:Fraction(0),3:Fraction(0),4:Fraction(0)}
net_charge=Fraction(0)
for site in product((-1,0,1),repeat=3):
    boundary=sum(v != 0 for v in site)
    weight=Fraction(1,2**boundary)
    charge=1 if sum(site) % 2 == 0 else -1
    net_charge+=weight*charge
    radius2=sum(v*v for v in site)
    if radius2:
        shells[radius2]-=weight*charge
assert net_charge == 0
assert list(shells.values()) == [3,-3,1,0]
evjen=sum(float(n)/sqrt(r2) for r2,n in shells.items())


# Q2: independently evaluate the two infinite-lattice Madelung constants.
# The origin ion has charge +1 and the nearest unlike separation is one.
def madelung(structure, eta, real_cut=4, reciprocal_cut=8):
    if structure == "NaCl":
        length=2.0
        basis=[(tuple(float(v) for v in p),(-1)**sum(p))
               for p in product((0,1),repeat=3)]
        def structure_factor(h):
            return 8 if all(v % 2 for v in h) else 0
    elif structure == "CsCl":
        length=2/sqrt(3)
        basis=[((0.,0.,0.),1),((length/2,)*3,-1)]
        def structure_factor(h):
            return 2 if sum(h) % 2 else 0
    else:
        raise ValueError(structure)
    real_sum=0.
    for n in product(range(-real_cut,real_cut+1),repeat=3):
        for r,q in basis:
            r2=sum((r[j]+length*n[j])**2 for j in range(3))
            if r2:
                radius=sqrt(r2)
                real_sum+=q*erfc(eta*radius)/radius
    reciprocal_sum=0.
    for h in product(range(-reciprocal_cut,reciprocal_cut+1),repeat=3):
        h2=sum(v*v for v in h)
        if not h2:
            continue
        g2=(2*pi/length)**2*h2
        reciprocal_sum+=(structure_factor(h)*exp(-g2/(4*eta*eta))/g2)
    potential=real_sum+4*pi/length**3*reciprocal_sum-2*eta/sqrt(pi)
    return -potential


alpha_na=madelung("NaCl",1.6)
alpha_cs=madelung("CsCl",1.6)
close(alpha_na,madelung("NaCl",2.1),atol=1e-10,rtol=0)
close(alpha_cs,madelung("CsCl",2.1),atol=1e-10,rtol=0)
close(alpha_na,1.74756459463318,atol=1e-10,rtol=0)
close(alpha_cs,1.76267477307099,atol=1e-10,rtol=0)
ratio_contact=(sqrt(3)+1)/2
ratio_crossover=alpha_cs/(sqrt(3)*alpha_na-alpha_cs)
assert ratio_contact < ratio_crossover < sqrt(2)+1
ratio_test=(ratio_contact+ratio_crossover)/2
r_cs=max(ratio_test+1,sqrt(3)*ratio_test,sqrt(3))
r_na=max(ratio_test+1,sqrt(2)*ratio_test,sqrt(2))
assert r_cs > ratio_test+1
assert -alpha_cs/r_cs < -alpha_na/r_na
# Equal radius result within this Coulomb/hard-sphere comparison.
assert -alpha_cs/2 < -alpha_na/2


# Q3(a,b): numerical differentiation in V per ion pair checks the bulk modulus.
r0=1.27
exponent=8.4
coulomb=alpha_na*1.31
repulsion=coulomb*r0**(exponent-1)/exponent
def pair_energy_r(r):
    return repulsion/r**exponent-coulomb/r
em=-coulomb/r0
ecoh=pair_energy_r(r0)
close(em/(em-ecoh),exponent)
v0=2*r0**3
def pair_energy_v(v):
    return pair_energy_r((v/2)**(1/3))
step=v0*2e-4
second=(-pair_energy_v(v0+2*step)+16*pair_energy_v(v0+step)
        -30*pair_energy_v(v0)+16*pair_energy_v(v0-step)
        -pair_energy_v(v0-2*step))/(12*step*step)
bulk_numerical=v0*second
bulk_formula=em*ecoh/(18*r0**3*(ecoh-em))
close(bulk_numerical,bulk_formula,atol=2e-8,rtol=2e-7)

# Q3(c): changing two ionic charges includes their mutual correction once.
sites=list(product(range(-2,3),repeat=3))
charges=[1 if sum(r) % 2 == 0 else -1 for r in sites]
i=sites.index((0,0,0))
j=sites.index((1,0,0))
def distance(a,b):
    return sqrt(sum((x-y)**2 for x,y in zip(a,b)))
def electrostatic_energy(q):
    return sum(q[a]*q[b]/distance(sites[a],sites[b])
               for a in range(len(q)) for b in range(a))
phi_i=sum(charges[k]/distance(sites[i],sites[k])
          for k in range(len(sites)) if k != i)
phi_j=sum(charges[k]/distance(sites[j],sites[k])
          for k in range(len(sites)) if k != j)
changed=charges.copy()
changed[i]=changed[j]=0
delta_u=electrostatic_energy(changed)-electrostatic_energy(charges)
delta_qi,delta_qj=-charges[i],-charges[j]
pair_change=delta_qi*phi_i+delta_qj*phi_j+delta_qi*delta_qj/distance(sites[i],sites[j])
close(delta_u,pair_change,atol=1e-11,rtol=0)
# For bulk site potentials the separated cost is 2*alpha; adjacent transfer
# has cost 2*alpha-1. Their difference is the electron-hole attraction.
close(2*alpha_na-(2*alpha_na-1),1)


# Q4: invert several hydrogenic n=1,2 excitation pairs and back-substitute.
for gap,binding in ((8.0,.2),(3.4,.007),(1.8,.6)):
    first,second_peak=gap-binding,gap-binding/4
    recovered_gap=(4*second_peak-first)/3
    recovered_binding=4*(second_peak-first)/3
    close(recovered_gap,gap)
    close(recovered_binding,binding)


# Q5,6: integrate the charge density and independently Fourier-transform it.
zt=2.3  # external charge in chosen charge units
ktf=1.7
def rho_ind(r):
    return -zt*ktf**2*exp(-ktf*r)/(4*pi*r)
limit=35/ktf
total_induced=simpson(lambda r:0 if r == 0 else 4*pi*r*r*rho_ind(r),0,limit)
close(total_induced,-zt,atol=2e-10,rtol=0)
screening_fourier=[]
for q in (0.,.3,2.,5.):
    def integrand(r):
        if r == 0:
            return 0.
        sinc=1 if q == 0 else sin(q*r)/(q*r)
        return 4*pi*r*r*rho_ind(r)*sinc
    transformed=simpson(integrand,0,limit)
    predicted=-zt*ktf**2/(q*q+ktf**2)
    close(transformed,predicted,atol=2e-10,rtol=0)
    if q:
        dielectric=1+ktf**2/q**2
        close(predicted,zt*(1/dielectric-1))
    screening_fourier.append((q,transformed,predicted))
for r in (1e-4,1e-5,1e-6):
    close(r*rho_ind(r),-zt*ktf**2/(4*pi),atol=4e-4,rtol=0)
x=1e-5
local_fraction=-expm1(-x)-x*exp(-x)
close(local_fraction/(x*x/2),1,atol=1e-5,rtol=0)


# Q7(a): zeros of the supplied truncated dielectric function.
wp=1.3
vf=.72
plasmon_results=[]
for q in (.001,.01,.1):
    freq2=wp**2/2*(1+sqrt(1+12*vf**2*q*q/(5*wp**2)))
    freq=sqrt(freq2)
    residual=1-wp**2/freq2-(3/5)*wp**2*vf**2*q*q/freq2**2
    close(residual,0,atol=3e-16,rtol=0)
    coefficient=(freq/wp-1)/(vf**2*q*q/wp**2)
    plasmon_results.append((q,freq,coefficient))
close(plasmon_results[0][2],.3,atol=1e-7,rtol=0)


# Q7(b): direct complex arithmetic, passivity, and the longitudinal loss sum rule.
gamma=.41
def drude(z):
    return 1-wp*wp/(z*(z+1j*gamma))
def loss_formula(w):
    return gamma*w*wp**2/((w*w-wp*wp)**2+gamma*gamma*w*w)
for w in (.03,.2,.8,1.3,2.7,12.):
    close(-(1/drude(w)).imag,loss_formula(w),atol=1e-13,rtol=0)
    assert loss_formula(w)>0
def weighted_loss_on_unit_interval(t):
    if t == 0: return 0.
    if t == 1: return gamma*wp*wp
    w=t/(1-t)
    return w*loss_formula(w)/(1-t)**2
loss_sum=simpson(weighted_loss_on_unit_interval,0,1,40000)
close(loss_sum,pi*wp**2/2,atol=2e-9,rtol=0)


# Q7(c,d): numerical principal-value transforms without deleting an arbitrary
# interval around the pole. The subtracted constant has zero half-line PV.
def pv_halfline(function,w,intervals=30000):
    fw=function(w)
    def regular(t):
        if t == 1:
            return -fw
        x=t/(1-t)
        if abs(x-w) < 1e-6:
            step=1e-5*max(1.,w)
            derivative=(function(w-2*step)-8*function(w-step)
                        +8*function(w+step)-function(w+2*step))/(12*step)
            quotient=derivative/(2*w)
        else:
            quotient=(function(x)-fw)/((x-w)*(x+w))
        return quotient/(1-t)**2
    return simpson(regular,0,1,intervals)

kk_results=[]
for w in (.23,.8,1.6,3.):
    # Smooth real part for nonzero real frequency.
    real_chi=lambda x:-wp*wp/(x*x+gamma*gamma)
    # x*Im(chi), extended continuously to x=0.
    ximag_chi=lambda x:wp*wp*gamma/(x*x+gamma*gamma)
    real_from_kk=2/pi*pv_halfline(ximag_chi,w)
    smooth_imag_from_kk=-2*w/pi*pv_halfline(real_chi,w)
    correct_imag=smooth_imag_from_kk+wp*wp/(gamma*w)
    close(real_from_kk,(drude(w)-1).real,atol=2e-8,rtol=0)
    close(correct_imag,drude(w).imag,atol=2e-8,rtol=0)
    assert abs(smooth_imag_from_kk-drude(w).imag)>.1
    kk_results.append((w,real_from_kk,correct_imag))

# A genuine lower-half-plane regulator obeys both ordinary KK relations.
eta=.11
def regulated(w):
    return 1-wp*wp/((w+1j*eta)*(w+1j*gamma))
for w in (.23,.8,1.6):
    real_from_kk=2/pi*pv_halfline(lambda x:x*regulated(x).imag,w)
    imag_from_kk=-2*w/pi*pv_halfline(lambda x:(regulated(x)-1).real,w)
    close(real_from_kk,(regulated(w)-1).real,atol=2e-8,rtol=0)
    close(imag_from_kk,regulated(w).imag,atol=2e-8,rtol=0)

# A causal time-domain kernel independently reproduces Drude in Im(z)>0.
# This validates causality despite the zero-frequency boundary pole.
causal_checks=[]
for z in (.7+.8j,-1.1+.3j):
    def time_integrand(t):
        kernel=wp*wp/gamma*(-expm1(-gamma*t))
        return kernel*cmath.exp(1j*z*t)
    transform=simpson(time_integrand,0,45/z.imag,40000)
    close(transform,drude(z)-1,atol=2e-9,rtol=0)
    causal_checks.append((str(z),abs(transform-(drude(z)-1))))

report={
 "status":"PASS",
 "evjen_coefficients":[int(x) for x in shells.values()],
 "evjen_small_cell_value":evjen,
 "ewald_madelung":{"NaCl":alpha_na,"CsCl":alpha_cs},
 "cscl_contact_radius_ratio":ratio_contact,
 "coulomb_model_cscl_nacl_crossover_ratio":ratio_crossover,
 "bulk_modulus":{"finite_difference":bulk_numerical,"formula":bulk_formula},
 "finite_cluster_pair_change_residual":delta_u-pair_change,
 "screening_charge_integral":total_induced,
 "screening_fourier_checks":screening_fourier,
 "plasmon_small_q_checks":plasmon_results,
 "loss_first_moment":{"numerical":loss_sum,"prediction":pi*wp**2/2},
 "conductivity_corrected_kk_checks":kk_results,
 "causal_kernel_transform_errors":causal_checks,
 "scope":"Checks the stated approximations and identifies their limits; not a claim of exact material prediction."
}
print(json.dumps(report,indent=2))
