import sympy as syp
from sympy.stats import Uniform, cdf


def referral(vs, vw, n, k, endogenous=True):
    rr = syp.Symbol('r_{r}')
    p = syp.Symbol('p')
    c = syp.Symbol('c')
    delta = syp.Symbol('delta')
    if endogenous:
        alpha = syp.Symbol('alpha')
    else:
        alpha = 1 - p
    beta = syp.Symbol('beta')
    vs = vs
    vw = vw
    n = n
    k = k

    vr = p + c - alpha * k * rr

    g1 = (p - rr + vs) * (
            (alpha * (1 - beta) * k * delta * (1 - (alpha * (1 - beta) * k * delta) ** n)) / (
            1 - alpha * (1 - beta) * k * delta))
    g2 = (p - rr + vw) * (
            (alpha * beta * k * delta * (1 - (alpha * beta * k * delta) ** n)) / (1 - alpha * beta * k * delta))

    f1 = (1 - vr) * (g1 + g2)
    f2 = (1 - p) * p

    L = syp.simplify(f1 + f2)

    gradL = [syp.diff(L, x) for x in [rr]]

    KKT_eqs = gradL

    stationary_points = syp.solve(KKT_eqs, [rr], dict=True)

    ans_rr = stationary_points[0][rr]

    X = Uniform('x', 0, 1)
    G = syp.simplify((1 - cdf(X)(vr)) * (g1 + g2) + f2)

    return G, ans_rr


def group(endogenous=True):
    rg = syp.Symbol('r_{g}')

    p = syp.Symbol('p')
    delta = syp.Symbol('delta')
    if endogenous:
        alpha = syp.Symbol('alpha')
    else:
        alpha = 1-p+rg

    v_upp = p + (delta / (1 - delta)) * alpha * rg
    v_low = p

    f1 = (1 - v_upp) * p
    f2 = 2 * alpha * (p - rg) * (v_upp - v_low)
    f3 = (1 - alpha) * p * (v_upp - v_low)

    L = syp.simplify(f1 + f2 + f3)

    gradL = [syp.diff(L, x) for x in [rg]]

    KKT_eqs = gradL

    stationary_points = syp.solve(KKT_eqs, [rg], dict=True)

    if endogenous:
        ans_rg = stationary_points[0][rg]
    else:
        ans_rg = stationary_points[2][rg]

    X = Uniform('x', 0, 1)

    g1 = (1-cdf(X)(v_upp))*p
    g2 = (cdf(X)(v_upp)-cdf(X)(v_low)) * 2 * alpha * (p - rg)
    g3 = (cdf(X)(v_upp)-cdf(X)(v_low)) * (1 - alpha) * p
    G = syp.simplify(g1+g2+g3)

    return G, ans_rg


def mixed(vs, vw, n, k, endogenous=True):
    rr = syp.Symbol('r_{r}')
    rg = syp.Symbol('r_{g}')

    p = syp.Symbol('p')
    c = syp.Symbol('c')
    delta = syp.Symbol('delta')
    if endogenous:
        alphar = syp.Symbol('alpha_{r}')
        alphag = syp.Symbol('alpha_{g}')
    else:
        alphar = 1-p
        alphag = 1-p+rg
    beta = syp.Symbol('beta')

    vs = vs
    vw = vw
    n = n
    k = k

    vr = p + (1 / (1 - delta)) * (delta * alphag * rg + c - alphar * rr)
    v_low = p - alphag*rg

    fs1 = (p - alphar * (1 - beta) * k * rr + vs) * (
                (alphar * (1 - beta) * k * delta * (1 - (alphar * (1 - beta) * k * delta) ** n)) / (
                    1 - alphar * (1 - beta) * k * delta))
    fs2 = (p - alphar * beta * k * rr + vw) * (
                (alphar * beta * k * delta * (1 - (alphar * beta * k * delta) ** n)) / (1 - alphar * beta * k * delta))

    f1 = (1 - vr) * (fs1 + fs2)
    f2 = (1 - vr) * p

    g1 = 2 * alphag * (p - rg) * (vr - v_low)
    g2 = (1 - alphag) * p * (vr - v_low)

    L = syp.simplify(f1 + f2 + g1 + g2)

    gradL = [syp.simplify(syp.diff(L, x)) for x in [rr, rg]]

    KKT_eqs = gradL

    stationary_points = syp.solve(KKT_eqs, [rr, rg], dict=True)

    ans_rr = stationary_points[0][rr]
    ans_rg = stationary_points[0][rg]

    X = Uniform('x', 0, 1)

    gs1 = syp.Max((cdf(X)(vr)-cdf(X)(v_low)), 0) * 2 * alphag * (p - rg)
    gs2 = syp.Max((cdf(X)(vr)-cdf(X)(v_low)), 0) * (1 - alphag) * p

    Z = syp.simplify((1 - cdf(X)(vr)) * (fs1+fs2+p) + gs1+gs2)

    return Z, ans_rr, ans_rg



