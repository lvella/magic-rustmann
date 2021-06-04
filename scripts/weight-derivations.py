import sympy
import itertools

D = sympy.sympify(2)

C0 = [
	[0, -1,  1,  0,  0,  1, -1,  1, -1],
	[0,  0,  0, -1,  1,  1, -1, -1,  1]
]

C = sympy.symbols(f"C:{D}")
w = sympy.symbols(f"w:{3}")
w = [w[0], w[1], w[1], w[1], w[1], w[2], w[2], w[2], w[2]]


orthogonal_tensors = [
	sympy.sympify(1),
	2*C[0],
	2*C[1],
	2*C[0]**2 - 1,
	2*C[1]**2 - 1,
	2*C[0]*C[1]
]

exp = sympy.exp(-(C[0]**2 + C[1]**2))

def eq(tensor_a, tensor_b):
	print(tensor_a, tensor_b)
	quadrature = 0
	for i in range(len(w)):
		sub_list = [(a, b[i]) for a, b in zip(C, C0)]
		quadrature += w[i] * tensor_a.subs(sub_list) * tensor_b.subs(sub_list)

	rhs = 1/sympy.pi**(D/2) * sympy.integrate(
		sympy.integrate(exp*tensor_a*tensor_b, (C[0], -sympy.oo, sympy.oo)),
		(C[1], -sympy.oo, sympy.oo))

	return sympy.Eq(quadrature, rhs)

eqs = [eq(a,b) for a,b in itertools.combinations_with_replacement(orthogonal_tensors, 2)]
eqs = eqs[:-6]
for e in eqs:
	print(e)

print(sympy.solve(eqs, w))
