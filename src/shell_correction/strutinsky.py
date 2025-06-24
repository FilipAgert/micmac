import numpy as np

import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import hermite, factorial

def f(x):
    val = np.exp(-x**2) / np.sqrt(np.pi)
    summer = 35/16 -x**2*35/8  + x**4*7/4 - x**6/6
    val = val * summer
    return val


def avg_lev_dens(levels, e, gamma):
    sum = 0.0
    for i in range(len(levels)):
        leve = levels[i]
        sum += f((leve-e)/gamma)*2
    return sum/gamma

def avg_e(levels, gamma, fermi):
    integrand = lambda E: avg_lev_dens(levels,E,gamma) * E
    #dens = avg_lev_dens(levels, E, gamma)
    #print(E)
    Evg = quad(integrand, -np.inf, fermi)[0]
    return Evg

def fermi(A, levels,gamma):
    rang = abs(np.max(levels)-np.min(levels))

    rolling = 0.0
    dens = lambda E: avg_lev_dens(levels,E,gamma)
    E_fsh = fermish(levels,A)
    Elower = E_fsh-2
    sum = quad(dens, -np.inf, Elower)[0]
    Eupper = E_fsh + 5

    val = sum
    nextEval = (Eupper+ Elower)/2
    func = lambda x: quad(dens,Elower,x)[0]
    lb = Elower
    while(abs(A-val) > 1e-5):

        searched = func(nextEval)
        val = sum + searched
        if(val < A):
            lb = nextEval
            nextEval = (Eupper + nextEval)/2

        else:
            Eupper = nextEval
            nextEval = (lb + nextEval)/2
    
    return nextEval


def fermish(levels,A):
    return levels[(A//2 - 1)]
    


def shell_en(levels,A):
    parts = 0
    en = 0.0
    for i in range(len(levels)):
        E = levels[i]
        en += E*2
        parts += 2
        if parts == A:
            return en
    
def getMajorShell(N):#N particle number
    NumParts = lambda shell: (shell+1)*(shell+2)
    i = 0
    particles = 0
    while(particles < N):
        parts = NumParts(i)
        particles += parts
        i+=1

    return i-1
    #which major shell are we part of?

def getMajorShellLevels(levels, majorShell, Nshells):
    NumParts = lambda shell: (shell+1)*(shell+2)//2

    i = 0
    excluded = 0
    
    for i in range(majorShell - Nshells//2): 
        excluded += NumParts(i)
    numIncluded = NumParts(majorShell -3) + NumParts(majorShell +3) + NumParts(majorShell - 1) + NumParts(majorShell - 2) + NumParts(majorShell+ 1) + NumParts(majorShell+ 2) + NumParts(majorShell)
    majors =  levels[excluded:excluded + numIncluded] 
    
    if len(majors)!= numIncluded:
        print("error. wrong number of particles")
    return majors

def getPartsUpTo(N):
    NumParts = lambda shell: (shell+1)*(shell+2)
    parts = 0
    for n in range(N + 1):
        parts += NumParts(n)
    return parts

def getNewParticleNumber(N, numShells):#
    major=getMajorShell(N)
    val = numShells//2 #3 shells => 1.5=> 1. 
    return N - getPartsUpTo(major - val)

def occupation_numbers(N, levels, gamma, fermi):
    occ_numbers = []

    integrand = lambda x: f(x)
    for i in range(len(levels)):
        E = levels[i]
        t = (fermi- E)/gamma
        occ_numbers.append(quad(integrand, -np.inf, t)[0])
    return np.array(occ_numbers)

def F(N, levels, gamma, fermi):

    F = 0.0
    integrand = lambda x: f(x)*x
    for i in range(len(levels)):
        E = levels[i]
        t = (fermi- E)/gamma
        F += (quad(integrand, -np.inf, t)[0])
    return F*2*gamma**2

A = 208
Z = 82
N = A- Z
getMajorShell(Z)
Numneutrons = 82
hbaromega = 6
Nlevels = lambda N: (N+1)*(N+2)//2
# en = lambda N: hbaromega*(N+1.5)
# levels=[]
# for N in range(16):
#     for i in range(Nlevels(N)):
#         levels.append(en(N))

levels = np.array(list(reversed([-0.08482,   -0.11710,   -0.15064,   -0.88889,   -0.98265,   -1.06116,   -1.30484,   -1.44291,   -1.49195,   -1.78015,   -2.27020,   -2.28252,   -2.50899,   -2.56979,   -3.55000,   -3.77372,   -3.77825,   -3.86038,   -3.99642,   -4.02570,   -4.18956,   -4.35498,   -4.54601,   -4.91558,   -5.13237,   -5.26893,   -5.93304,   -5.95113,   -6.43458,   -6.85591,   -6.95941,   -7.15877,   -7.60246,   -7.64318,   -7.95299,   -8.57590,   -8.65570,   -9.03908,   -9.20628,   -9.38173,   -9.78454,  -10.35152,  -10.35211,  -10.54618,  -10.85994,  -11.38454,  -11.42932,  -11.48061,  -11.82666,  -12.14451,  -12.37582,  -12.40288,  -12.56244,  -12.69512,  -13.13161,  -13.73660,  -13.76212,  -14.08518,  -14.22323,  -14.24397,  -15.46421,  -16.15328,  -16.68735,  -16.75360,  -16.92512,  -17.16079,  -17.36182,  -17.49778,  -17.92755,  -17.94355,  -18.03296,  -18.05011,  -19.12016,  -19.67618,  -19.71623,  -20.09116,  -20.37311,  -20.59331,  -20.68853,  -20.91844,  -21.48024,  -22.33032,  -22.64109,  -22.79881,  -22.89841,  -23.18351,  -23.51544,  -23.68041,  -24.29661,  -24.73720,  -25.44275,  -25.54403,  -25.77828,  -26.29156,  -26.29381,  -26.66213,  -27.46270,  -27.81508,  -28.14022,  -28.39076,  -28.92810,  -29.37788,  -29.45006,  -29.65315,  -29.66984,  -29.84097,  -31.63552,  -31.76891,  -31.89843,  -31.95714,  -33.18039,  -33.32604,  -33.70419,  -34.09381,  -34.59549,  -34.63992,  -35.35693,  -35.83291,  -36.85630,  -37.60932])))
E_sh = shell_en(levels, Numneutrons)
gamma = hbaromega
fermi_avg = fermish(levels, Numneutrons)
converged = False



while not converged:
        
    diff = lambda ef: (sum(occupation_numbers(Numneutrons,levels,gamma,ef)) - Numneutrons/2)


    print(diff(fermi_avg-15), diff(fermi_avg+5))
    ef = brentq(diff, fermi_avg - 5, fermi_avg + 5, xtol =1e-3)
    print("Fermi energy:", ef)
    print("fermi avgd", fermi_avg)
    print("num parts: ", sum(occupation_numbers(Numneutrons,levels,gamma,ef)))
    shell_energy = lambda gamma: E_sh - (2*sum(occupation_numbers(Numneutrons,levels,gamma,ef) * levels)  + F(Numneutrons,levels,gamma,ef)) 
    dg = 1e-3
    fp = shell_energy(gamma+dg)
    fm = shell_energy(gamma-dg)
    ddf = (fp - 2*shell_energy(gamma) + fm)/(dg**2)
    df = (fp - fm)/(2*dg)
    gamma -= df/ddf
    if abs(df/ddf) < 1e-3:
        converged = True
    print("gamma:", gamma)
    print("Shell energy:", shell_energy(gamma))




# gammas = np.linspace(0.5*hbaromega,2*hbaromega,30)
# fermi_energies = [fermi(newparticleNumber,newl,gamma) for gamma in gammas]
# avg_es = np.array([avg_e(newl, gamma, fermi) for gamma, fermi in zip(gammas, fermi_energies)])

# fermi_energies_2 = [fermi(Numneutrons,levels,gamma) for gamma in gammas]
# avg_es_2 = np.array([avg_e(levels, gamma, fermi) for gamma, fermi in zip(gammas, fermi_energies_2)])

# fig,ax=plt.subplots(1,1)

# ax.set_title("Shell correction")
# ax.set_ylabel(r"$E_{sh}$ (MeV)")
# ax.set_xlabel(r"$\gamma (\hbar \omega_0)$")

# #ax.plot(differences,f(differences), label="avg shell en", color='black', marker='.')
# shell_2 = shell_en(levels,Numneutrons)
# ax.plot(gammas/hbaromega, (shell_2 - avg_es_2), label="avg shell en", color='black', marker='.')

# plt.tight_layout()
# plt.show()
