#!/usr/bin/env python3

import numpy as np
from scipy.optimize import minimize
import sympy

xs = """
-2.020182870456085632929
-0.9585724646138185071128
0
0.9585724646138185071128
2.020182870456085632929
"""
xs = np.fromiter(xs.split(), float)

real_ws = """
0.01995324205904591320774
0.3936193231522411598285
0.9453087204829418812257
0.393619323152241159828
0.01995324205904591320774
"""
real_ws = np.fromiter(real_ws.split(), float)

def calc(N):
    points = np.random.rand(N+1, 2*len(xs)) * 20.0 - 10.0

    # Exact solution
    def exact():
        x = sympy.symbols('x')
        aa = sympy.symbols(f'a0:{N+1}')
        expr = 0
        for i in range(N+1):
            expr += aa[i] * x**i
        expr = expr*sympy.exp(-x**2)

        expr = sympy.integrate(expr, (x, -sympy.oo, sympy.oo))
        print("Exact solution:", expr)

        return sympy.lambdify(aa, expr, "numpy")(*points)

    exact = exact()

    # Points for quadrature solution
    f = np.zeros_like(xs)
    for i in range(N+1):
        f = np.multiply.outer(points[i], xs**i) + f

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


for N in range(0, 2*len(xs)-1):
    print(f"##### {N} ######")
    calc(N)
    print("=================")
