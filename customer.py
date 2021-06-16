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
numerical_beta = 0.5
numerical_c = 0.05
numerical_delta = 0.9
numerical_k = 1

alphar = syp.Symbol('alpha_{r}')
alphag = syp.Symbol('alpha_{g}')

#%% plot on price p
R, r = referral(vs=0.5, vw=0, n=1, k=numerical_k, endogenous=False)
G, g = group(endogenous=False)

R = R.subs([(c, numerical_c), (alpha, alphar), (beta, numerical_beta), (delta, numerical_delta)])
G = G.subs([(alpha, alphag), (delta, numerical_delta)])

r = r.subs([(c, numerical_c), (alpha, alphar), (beta, numerical_beta), (delta, numerical_delta)])
g = g.subs([(alpha, alphag), (delta, numerical_delta)])

vr = p + (delta*alphag*rg + c - alphar*rr)/(1-delta)

x = np.arange(0.05, 1, 0.05)

plot_r = syp.lambdify(p, r, 'numpy')
plot_subsidy_r = []
X = Uniform('x', 0, 1)


for index, value in enumerate(plot_r(x)):
    if value <= 0:
        #plot_subsidy_r.append(min(x[index], numerical_c/(numerical_k*numerical_alpha))) # endogenous
        plot_subsidy_r.append(min(x[index], numerical_c / (numerical_k * (1-x[index])))) # exogenous

    elif value >= x[index]:
        plot_subsidy_r.append(x[index])
    else:
        #plot_subsidy_r.append(max(value, numerical_c/(numerical_k*numerical_alpha))) # endogenous
        plot_subsidy_r.append(max(value, numerical_c/(numerical_k*(1-x[index])))) # exogenous

plot_subsidy_g = [g.subs([(p, i)]) for i in x]

plot_vr = [min(vr.subs([(delta, numerical_delta), (c, numerical_c), (p, x[i]), (alphar, 1-x[i]), (alphag, 1-x[i]+plot_subsidy_g[i]), (rr, plot_subsidy_r[i]), (rg, plot_subsidy_g[i])]), 1) for i in range(len(x))]

plot_vr = np.array(plot_vr, dtype=float)

plt.plot(x, plot_vr, label='referral')
plt.plot(x, x, label='group')

plt.fill_between(x, x, plot_vr, color='orange', alpha='0.2')
plt.fill_between(x, plot_vr, 1, color='blue', alpha='0.2')


plt.legend(loc=4)
plt.xlabel('price $p$')
plt.ylabel('customer valuation $v$')
plt.show()

#%% plot on alpha
R, r = referral(vs=numerical_p, vw=0, n=1, k=numerical_k, endogenous=True)
G, g = group(endogenous=True)

R = R.subs([(c, numerical_c), (p, numerical_p), (beta, numerical_beta), (delta, numerical_delta)])
G = G.subs([(p, numerical_p), (delta, numerical_delta)])

r = r.subs([(c, numerical_c), (p, numerical_p), (beta, numerical_beta), (delta, numerical_delta)])
g = g.subs([(p, numerical_p), (delta, numerical_delta)])

vr = p + (delta*alphag*rg + c - alphar*rr)/(1-delta)

x = np.arange(0.05, 1, 0.05)

plot_r = syp.lambdify(alpha, r, 'numpy')
plot_subsidy_r = []


for index, value in enumerate(plot_r(x)):
    if value >= numerical_p:
        plot_subsidy_r.append(numerical_p)
    elif value <= 0:
        if numerical_c/(numerical_k*x[index]) > numerical_p:
            plot_subsidy_r.append(np.nan)
        else:
            plot_subsidy_r.append(numerical_c/(numerical_k*x[index])) # maybe other numbers
    else:
        plot_subsidy_r.append(max(value, numerical_c/(numerical_k*x[index])))


plot_subsidy_g = [g.subs([(p, i)]) for i in x]

plot_vr = []
for i in range(len(x)):
    if np.isnan(plot_subsidy_r[i]):
        plot_vr.append(np.nan)
    elif vr.subs([(delta, numerical_delta), (c, numerical_c), (p, x[i]), (alphar, x[i]), (alphag, syp.Min(x[i]*2, 1)), (rr, plot_subsidy_r[i]), (rg, plot_subsidy_g[i])]) < numerical_p:
        plot_vr.append(numerical_p)
    else:
        plot_vr.append(min(vr.subs([(delta, numerical_delta), (c, numerical_c), (p, x[i]), (alphar, x[i]), (alphag, syp.Min(x[i]*2, 1)), (rr, plot_subsidy_r[i]), (rg, plot_subsidy_g[i])]), 1))

#plot_vr = [min(vr.subs([(delta, numerical_delta), (c, numerical_c), (p, x[i]), (alphar, x[i]), (alphag, x[i]), (rr, plot_subsidy_r[i]), (rg, plot_subsidy_g[i])]), 1) for i in range(len(x))]

plot_vr = np.array(plot_vr, dtype=float)

plt.plot(x, plot_vr, label='referral')
plt.plot(x, [0.6]*len(x), label='group')

na_index = np.argwhere(np.isnan(plot_vr)).tolist()[-1][0]

plt.fill_between(x[:na_index+2], 0.6, 1, color='orange', alpha='0.2')
plt.fill_between(x, 0.6, plot_vr, color='orange', alpha='0.2')
plt.fill_between(x, plot_vr, 1, color='blue', alpha='0.2')


plt.legend()
plt.xlabel(r'$\alpha$')
plt.ylabel('customer valuation $v$')
plt.show()
