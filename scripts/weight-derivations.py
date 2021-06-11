import sympy
import sympy.logic
import mpmath
import itertools

mpmath.dps = 50

D = sympy.sympify(2)
e = [
	[0, -1,  1,  0,  0,  1, -sympy.sympify(3)/2,  1, -1],
	[0,  0,  0, -1,  1,  1, -sympy.sympify(3)/2, -1,  1]
]
w_idx = [0,  1,  2,  1,  2,  3,                   4,  5,  5]

#e = [
#	[0, -1,  1,  0,  0,  1, -1,  1, -1],
#	[0,  0,  0, -1,  1,  1, -1, -1,  1]
#]
#w_idx = [0,  1,  1,  1,  1,  2,  2,  2,  2]

C = sympy.symbols(f"C:{D}")

a = sympy.symbols("a", positive=True)
unique_w = list(sympy.symbols(f"w:{max(w_idx)+1}", positive=True))
w = [unique_w[i] for i in w_idx]

orthogonal_tensors = [
	sympy.sympify(1),
	C[0],
	C[1],
	C[0]**2 - sympy.sympify(1)/2,
	C[0]*C[1],
	C[1]**2 - sympy.sympify(1)/2,
]

exp = sympy.exp(-(C[0]**2 + C[1]**2))

def eq(tensor_a, tensor_b):
	quadrature = 0
	for i in range(len(w)):
		sub_list = [(src, a * dest[i]) for src, dest in zip(C, e)]
		quadrature += w[i] * tensor_a.subs(sub_list) * tensor_b.subs(sub_list)

	rhs = 1/sympy.pi**(D/2) * sympy.integrate(
		sympy.integrate(exp*tensor_a*tensor_b, (C[0], -sympy.oo, sympy.oo)),
		(C[1], -sympy.oo, sympy.oo))

	return sympy.Eq(quadrature, rhs)

eqs = [eq(ta, tb) for ta, tb in itertools.combinations_with_replacement(orthogonal_tensors, 2)]
eqs = set([eq for eq in eqs if not eq is sympy.S.true])

def least_squares(eqs, symbols):
	def to_func(expr):
		return sympy.lambdify(symbols, expr, "mpmath")

	eqs = [(eq.lhs - eq.rhs)**2 for eq in eqs]
	fs = [to_func(eq) for eq in eqs]
	ff = lambda *xs: [f(*xs) for f in fs]

	x0 = [0.2] * len(symbols)

	js = [
		[to_func(sympy.diff(f, x)) for x in symbols]
		for f in eqs
	]
	def jj(*xs):
		return [
			[h(*xs) for h in g]
			for g in js
		]

	print(mpmath.findroot(ff, x0, J=jj, solver='mnewton', maxsteps=100000, tol=1e-100, verbose=True))

least_squares(eqs, unique_w + [a])

print(sympy.solve(eqs, unique_w + [a]))
