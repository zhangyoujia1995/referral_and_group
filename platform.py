import sympy as syp
import matplotlib.pyplot as plt
import numpy as np
from function import referral, group
from sympy.stats import Uniform, cdf


syp.init_printing(use_latex=True)
k = syp.Symbol('k')
p = syp.Symbol('p')
c = syp.Symbol('c')
delta = syp.Symbol('delta')
alpha = syp.Symbol('alpha')
beta = syp.Symbol('beta')
rr = syp.Symbol('r_{r}')
rg = syp.Symbol('r_{g}')


numerical_p = 0.6
numerical_alpha = 0.4
numerical_beta = 0.5
numerical_c = 0.05
numerical_delta = 0.9
numerical_k = 1

R, r = referral(vs=0.5, vw=0, n=1, k=numerical_k, endogenous=True)
G, g = group(endogenous=True)

#%% plot on alpha
R = R.subs([(c, numerical_c), (p, numerical_p), (beta, numerical_beta), (delta, numerical_delta)])
G = G.subs([(p, numerical_p), (delta, numerical_delta)])

r = r.subs([(c, numerical_c), (p, numerical_p), (beta, numerical_beta), (delta, numerical_delta)])
g = g.subs([(p, numerical_p), (delta, numerical_delta)])

x = np.arange(0.1, 1, 0.05)

plot_r = syp.lambdify(alpha, r, 'numpy')
plot_subsidy_r = []
X = Uniform('x', 0, 1)


for index, value in enumerate(plot_r(x)):
    if value >= numerical_p:
        if numerical_p > numerical_c/(numerical_k*x[index]):
            plot_subsidy_r.append(numerical_p)
        else:
            plot_subsidy_r.append(np.nan)
    elif value <= 0:
        if numerical_c/(numerical_k*x[index]) > numerical_p:
            plot_subsidy_r.append(np.nan)
        else:
            plot_subsidy_r.append(numerical_c/(numerical_k*x[index])) # maybe other numbers
    else:
        plot_subsidy_r.append(max(value, numerical_c/(numerical_k*x[index])))


plot_R = []

for index, value in enumerate(plot_subsidy_r):
    if np.isnan(value):
        plot_R.append(numerical_p*(1-numerical_p))
    else:
        plot_R.append(R.subs([(rr, value), (alpha, x[index])]))

#plot_R = [R.subs([(rr, value), (alpha, x[index])]) for index, value in enumerate(plot_subsidy_r)]

plot_subsidy_g = [g.subs([(p, i)]) for i in x]
plot_G = [G.subs([(rg, g), (alpha, i)]) for i in x]


fig, axs = plt.subplots(2)
axs[0].plot(x, plot_R, label='referral')
axs[0].plot(x, plot_G, label='group')
axs[0].set_ylabel('profit $\pi$')
axs[0].set_xlabel(r'$\alpha$')
axs[0].legend()

axs[1].plot(x, plot_subsidy_r, label='referral')
axs[1].plot(x, plot_subsidy_g, label='group')
axs[1].axhline(y=numerical_p, linestyle='--', color='black')
axs[1].set_ylabel('subsidy $r$')
axs[1].set_xlabel(r'$\alpha$')
axs[1].legend()
plt.show()

#%% plot on price p
R = R.subs([(c, numerical_c), (alpha, numerical_alpha), (beta, numerical_beta), (delta, numerical_delta)])
G = G.subs([(alpha, numerical_alpha), (delta, numerical_delta)])

r = r.subs([(c, numerical_c), (alpha, numerical_alpha), (beta, numerical_beta), (delta, numerical_delta)])
g = g.subs([(alpha, numerical_alpha), (delta, numerical_delta)])

x = np.arange(0.05, 1, 0.05)

plot_r = syp.lambdify(p, r, 'numpy')
plot_subsidy_r = []
X = Uniform('x', 0, 1)


for index, value in enumerate(plot_r(x)):
    if value <= 0:
        plot_subsidy_r.append(min(x[index], numerical_c/(numerical_k*numerical_alpha))) # endogenous
        #plot_subsidy_r.append(min(x[index], numerical_c / (numerical_k * (1-x[index])))) # exogenous

    elif value >= x[index]:
        plot_subsidy_r.append(x[index])
    else:
        plot_subsidy_r.append(max(value, numerical_c/(numerical_k*numerical_alpha))) # endogenous
        #plot_subsidy_r.append(max(value, numerical_c/(numerical_k*(1-x[index])))) # exogenous

plot_R = []
for index, value in enumerate(plot_subsidy_r):
    if value == 0:
        plot_R.append(x[index]*(1-x[index]))
    else:
        plot_R.append(R.subs([(rr, value), (p, x[index])]))

#plot_R = [R.subs([(rr, value), (p, x[index])]) for index, value in enumerate(plot_subsidy_r)]

plot_subsidy_g = [g.subs([(p, i)]) for i in x]
plot_G = [G.subs([(rg, g), (p, i)]) for i in x]

fig, axs = plt.subplots(2)
axs[0].plot(x, plot_R, label='referral')
axs[0].plot(x, plot_G, label='group')
axs[0].plot(x, x*(1-x), linestyle='--', color='black')
axs[0].set_ylabel('profit $\pi$')
axs[0].set_xlabel('price $p$')
axs[0].legend()

axs[1].plot(x, plot_subsidy_r, label='referral')
axs[1].plot(x, plot_subsidy_g, label='group')
axs[1].plot(x, x, color='black', linestyle='--')
axs[1].set_ylabel('subsidy $r$')
axs[1].set_xlabel('price $p$')
axs[1].legend()
plt.show()