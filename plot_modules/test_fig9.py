'''

Reproduce Pallottini17b Fig. 9

Last mod: 10 Sept 2018

'''

import numpy as np
import matplotlib.pyplot as plt
from plot_cloud_prop import setup_plot

setup_plot()

path = '../literature/data/'


plt.close('all')
plt.figure()
# Althaea

f_in = path + 'simulations.dat'
x, y = np.loadtxt(f_in, usecols=(1, 2))
x_althaea, y_althaea = x[0], y[0]

plt.plot(np.log10(x_althaea), np.log10(y_althaea), marker='*', markersize=15, color='red', label=r'Alth{\ae}a', linestyle='None', zorder=30, markeredgecolor='k')


file_list = ['sk_kennicutt1998.dat', 'sk_bouche07a.dat',  'sk_daddi10b.dat',
             'sk_daddi10a.dat', 'sk_tacconi10a.dat', 'sk_genzel10a.dat']
label_list = ['Kennicutt+98', r"Bouch\'e+07",
              r'Daddi+10', '', 'Tacconi+10', 'Genzel+10']
list_2 = ['sk_Heiderman10.dat', 'sk_Lada10.dat']
label_list_2 = ['Heiderman+10', 'Lada+10']

colors = ('k', 'g', 'darkred', 'c', 'blue', 'orange')
markers = ('.', 's', 'D', 'v', 'h', '8', '1')
morecolors = ('darkolivegreen', 'brown')
moremarkers = ('x', '+')

for i in xrange(len(file_list)):

    f_in = path + file_list[i]
    ss = np.loadtxt(f_in, usecols=(2, 3, 4, 5, 6, 7, 8, 9))
    x = ss[:, 6]
    y = ss[:, 3]

    plt.plot(x, y, label=label_list[i], marker=markers[i], linestyle='None', color=colors[i], markersize=5)


for i in xrange(len(list_2)):
    f_in = path + list_2[i]
    ss = np.loadtxt(f_in, usecols=(5, 8))

    x = np.log10(ss[:, 1])
    y = np.log10(ss[:, 0])
    plt.plot(x, y, marker=moremarkers[i], linestyle='None', label=label_list_2[i], color=morecolors[i], markersize=8)


plt.legend(loc="best")

from matplotlib.ticker import MultipleLocator
ml = MultipleLocator(5)
plt.axes().yaxis.set_minor_locator(ml)
plt.axes().xaxis.set_minor_locator(ml)
plt.minorticks_on()

plt.ylabel(r"$\log {\left(\Sigma_{\rm SFR}\right)}$ [$M_{\odot}$ Myr$^{-1}$ pc$^{-2}$]")
plt.xlabel(r"$\log {\left(\Sigma_{\rm gas}/t_{\rm ff}\right)}$ [$M_{\odot}$  pc$^{-2}$ Myr$^{-1}$]")
plt.savefig('Test_PalloFig9.png', bbox_inches="tight")

