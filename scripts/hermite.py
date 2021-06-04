import sympy
import itertools
from functools import reduce
import operator

def generate(symbols):
    inner_product_scale = 1/(sympy.pi**(sympy.sympify(len(symbols))/2))

    def inner(a, b):
        prod = a*b
        for s in symbols:
            prod = sympy.simplify(sympy.integrate(
                sympy.exp(-s**2) * prod,
                (s, -sympy.oo, sympy.oo)
            ))
        return inner_product_scale * prod

    def proj(v, onto):
        return inner(onto, v) / inner(onto, onto) * onto

    basis = [(sympy.sympify(1), 1)]
    yield (0, (), basis[0][0])

    for order in itertools.count(1):
        for comb in itertools.combinations_with_replacement(symbols, order):
            # Generate a polynomial in the right order
            e = sympy.sympify(1)
            for s in comb:
                e *= s

            # Gram-Schmidt it to orthonormality
            for b, b_norm in basis:
                e -= sympy.simplify(inner(b, e) / b_norm * b)
                proj(e, b)
            e_norm = inner(e, e)
            basis.append((e, e_norm))
            yield order, comb, e

if __name__ == '__main__':
    for order, comb, p in generate(sympy.symbols("x, y, z")):
        print("Order: ", order, comb)
        print(p)
