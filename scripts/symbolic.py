# Solution involves (2 - erfc((2*RT)**0.5 / (2*RT))) + erfc((2*RT)**0.5 / (2*RT))

import sympy

T = 300 # kelvin
R = 8.31446261815324 # SI
RT = R*T
a, b, c, d, f, g, h = sympy.symbols("a, b, c, d, f, g, h")

rho, u, v, x, y = sympy.symbols("rho, u, v, x, y")

g = rho / (sympy.pi*2*RT) * sympy.exp(-((x-u)**2 + (y-v)**2)/(2*RT))

#f = 1
#f = u
#f = (x-u)**2 + (y-v)**2
f = a + b * x + c * y + d * x**2 + f * y**2 + g * x**3 + h * x**3

res = sympy.integrate(sympy.integrate(f * g, (y, -sympy.oo, sympy.oo)), (x, -sympy.oo, sympy.oo))

print(sympy.simplify(res))
