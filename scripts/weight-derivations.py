import sympy
import itertools

D = sympy.sympify(2)

e = [
	[0, -1,  1,  0,  0,  1, -1,  1, -1],
	[0,  0,  0, -1,  1,  1, -1, -1,  1]
]

C = sympy.symbols(f"C:{D}")
w = sympy.symbols(f"w:{3}")
w = [w[0], w[1], w[1], w[1], w[1], w[2], w[2], w[2], w[2]]
a = sympy.symbols("a")

orthogonal_tensors = [
	sympy.sympify(1),
	C[0],
	C[1],
	C[0]**2 - sympy.sympify(1)/2,
	C[0]*C[1],
	C[1]**2 - sympy.sympify(1)/2,
]

exp = sympy.exp(-(C[0]**2 + C[1]**2))

def eq(tensor):
	print(tensor)
	quadrature = 0
	for i in range(len(w)):
		sub_list = [(src, a * dest[i]) for src, dest in zip(C, e)]
		quadrature += w[i] * tensor.subs(sub_list)**2

	rhs = 1/sympy.pi**(D/2) * sympy.integrate(
		sympy.integrate(exp*tensor**2, (C[0], -sympy.oo, sympy.oo)),
		(C[1], -sympy.oo, sympy.oo))

	return sympy.Eq(quadrature, rhs)

eqs = [eq(t) for t in orthogonal_tensors]
for e in eqs:
	print(e)

print(sympy.solve(eqs, itertools.chain(w, [a])))
