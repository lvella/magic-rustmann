#!/usr/bin/env python3

import numpy as np
import random
from scipy.optimize import minimize

xs = """
-4.499990707309391553664
-3.669950373404452534729
-2.967166927905603248489
-2.325732486173857745454
-1.719992575186488932416
-1.136115585210920666319
-0.565069583255575748526
0
0.565069583255575748526
1.136115585210920666319
1.719992575186488932416
2.32573248617385774545
2.967166927905603248489
3.669950373404452534729
4.499990707309391553664
"""
xs = np.fromiter(xs.split(), float)
xs2 = xs**2

real_ws = """
1.522475804253517020161E-9
1.059115547711066635775E-6
1.00004441232499868127E-4
0.002778068842912775896079
0.03078003387254608222868
0.1584889157959357468838
0.4120286874988986270259
0.5641003087264175328526
0.4120286874988986270259
0.1584889157959357468838
0.03078003387254608222868
0.00277806884291277589608
1.00004441232499868127E-4
1.059115547711066635775E-6
1.52247580425351702016E-9
"""
real_ws = np.fromiter(real_ws.split(), float)

points = np.random.rand(100) * 20.0 - 10.0

def calc(N):
    points = np.random.rand(N+1, 100) * 20.0 - 10.0

    f = np.zeros_like(xs)
    for i in range(N+1):
        f = np.multiply.outer(points[i], xs**i) + f

    exact = np.sum(real_ws * f, 1)

    def func(ws):
        quadrature = np.sum(ws * f, 1)
        return np.sum((exact - quadrature)**2)

    def jac(ws):
        quadrature = np.sum(ws * f, 1)
        return np.sum(-2. * f.transpose() * (exact - quadrature), 1)

    def hess(ws):
        ret = np.zeros((len(xs), len(xs)))

        for fi in f:
            ret += np.multiply.outer(fi, fi)

        return ret

    x0 = np.array([0.5]*len(xs))
    #x0 = real_ws
    ans = minimize(func, x0, jac=jac, hess=hess, method='Newton-CG', tol=1e-15, options={"maxiter":2000000})

    print(ans)
    print(f"Known: {func(real_ws)}")
    print(f"Calculated: {func(ans.x)}")


for N in range(50):
    calc(N)
    print(f"##### {N} ######")
    print("=================")
