import numpy as np
import matplotlib.pyplot as plt
Mh_range = np.logspace(10, 15, 10)



def moster10(N, M1, beta, gamma, Mh):
    SHM = 2 * N  * ((Mh/M1)**-beta + (Mh/M1)**gamma)**-1
    return SHM


def moster13(z, M10=11.590, M11=1.195, N10=0.0351, N11=-0.0247, beta10=1.376, beta11=-0.826, gamma10=0.608, gamma11=0.329):

    factor = z / (z+1)
    logM1 = M10 + M11 * factor
    N = N10 + N11 * factor
    beta = beta10 + beta11 * factor
    gamma = gamma10 + gamma11 * factor

    return logM1, N, beta, gamma


# ----------- Moster 13 -----------------

logM1_z0, N_z0, beta_z0, gamma_z0 = moster13(z=0)
logM1_z4, N_z4, beta_z4, gamma_z4 = moster13(z=4)
logM1_z6, N_z6, beta_z6, gamma_z6 = moster13(z=6)
SHM_moster13 = moster10(N=N_z6, M1=10**logM1_z6, beta=beta_z6, gamma=gamma_z6, Mh=3.5e11)
print(SHM_moster13)



plt.figure()
SHM_moster13_range = [moster10(N=N_z0, M1=10**logM1_z0, beta=beta_z0, gamma=gamma_z0, Mh=i) for i in Mh_range]
plt.plot(np.log10(Mh_range), np.log10(SHM_moster13_range * Mh_range), '--', label='z0')

SHM_moster13_range = [moster10(N=N_z4, M1=10**logM1_z4, beta=beta_z4, gamma=gamma_z4, Mh=i) for i in Mh_range]
plt.plot(np.log10(Mh_range), np.log10(SHM_moster13_range * Mh_range), '--', label='z4')

SHM_moster13_range = [moster10(N=N_z6, M1=10**logM1_z6, beta=beta_z6, gamma=gamma_z6, Mh=i) for i in Mh_range]
plt.plot(np.log10(Mh_range), np.log10(SHM_moster13_range * Mh_range), '--', label='z6')
plt.xlabel('log M_halo')
plt.ylabel('log M_cen, star')
plt.legend()
plt.show(block=False)


plt.figure()
SHM_moster13_range = [moster10(N=N_z0, M1=10**logM1_z0, beta=beta_z0, gamma=gamma_z0, Mh=i) for i in Mh_range]
plt.plot(np.log10(Mh_range), np.log10(SHM_moster13_range), '--', label='z0')

SHM_moster13_range = [moster10(N=N_z4, M1=10**logM1_z4, beta=beta_z4, gamma=gamma_z4, Mh=i) for i in Mh_range]
plt.plot(np.log10(Mh_range), np.log10(SHM_moster13_range), '--', label='z4')

SHM_moster13_range = [moster10(N=N_z6, M1=10**logM1_z6, beta=beta_z6, gamma=gamma_z6, Mh=i) for i in Mh_range]
plt.plot(np.log10(Mh_range), np.log10(SHM_moster13_range), '--', label='z6')
plt.xlabel('log M_halo')
plt.ylabel('log (M_cen, star / Mh)')
plt.legend()
plt.show()




# --------------------------- Behroozi+13 -----------

def z2a(z):
    a = (z+1)**-1
    return a

def nu(z):
    a = z2a(z)
    return np.exp(-4 * a**2)

def logM1(z, M10=11.514, M1a=-1.793, M1z=-0.251):
    a = z2a(z)
    nu_val = nu(z)
    return M10 + (M1a * (a - 1) + M1z * z) * nu_val


def logepi(z, epi0=-1.777, epia=-0.006, epi_z=-0, epi_a2=-0.119):
    # sec 5 of Behroozi+13
    a = z2a(z)
    nu_val = nu(z)
    return epi0 + (epia * (a - 1) + epi_z * z) * nu_val + epi_a2 * (a - 1)


def alpha(z, alpha0=-1.412, alpha_a=0.731):
    a = z2a(z)
    nu_val = nu(z)
    return alpha0 + (alpha_a * (a - 1)) * nu_val


def delta(z, delta0=3.508, delta_a=2.608, delta_z=-0.043):
    a = z2a(z)
    nu_val = nu(z)
    delta = delta0 + (delta_a * (a - 1) + delta_z * z) * nu_val
    return delta


def gamma(z, gamma0=0.316, gamma_a=1.319, gamma_z=0.279):
    a = z2a(z)
    nu_val = nu(z)
    gamma = gamma0 + ( gamma_a * (a-1) + gamma_z * z) * nu_val
    return gamma


def f_func(x, alpha, delta, gamma):
    f = -np.log10(10**(alpha * x) + 1) + delta * (np.log10(1 + np.exp(x)))**gamma / (1 + np.exp(10**-x))
    return f


def logMstar(z, Mh):

    logepi_val = logepi(z)
    epi = 10**(logepi_val)

    logM1_val = logM1(z)
    M1 = 10**logM1_val

    alpha_val = alpha(z)
    delta_val = delta(z)
    gamma_val = gamma(z)
    f0 =  f_func(0, alpha_val, delta_val, gamma_val)
    fsomething = f_func(np.log10(Mh/M1), alpha_val, delta_val, gamma_val)
    print(fsomething)
    log_mstar = np.log10(epi * M1) + fsomething - f0
    return log_mstar


logMstar = logMstar(6, 3.e11)
print(10**logMstar/3.e11)




