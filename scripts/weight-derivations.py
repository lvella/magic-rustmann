from sympy import symbols, Eq, solve, pi, exp, integrate

cx = [0, -1,  1,  0,  0]
cy = [0,  0,  0, -1,  1]

rho, w0, w1, w2, w3, w4, u, v, D0, RT = symbols(
'rho, w0, w1, w2, w3, w4, u, v, D0, RT')

def m(i):
	return (cx[i] - u)**2 + (cy[i] - v)**2

def f(i):
	return rho / (2*pi*RT) * exp(-m(i)/(2*RT))

#eqs = [
#	Eq(rho, w0*f(0) + w1*f(1) + w2*f(2) + w3*f(3) + w4*f(4)),
#	Eq(rho * u, w0*f(0)*cx[0] + w1*f(1)*cx[1] + w2*f(2)*cx[2] + w3*f(3)*cx[3] + w4*f(4)*cx[4]),
#	Eq(rho * v, w0*f(0)*cy[0] + w1*f(1)*cy[1] + w2*f(2)*cy[2] + w3*f(3)*cy[3] + w4*f(4)*cy[4]),
#	Eq(rho * D0 * RT, w0*f(0)*m(0) + w1*f(1)*m(1) + w2*f(2)*m(2) + w3*f(3)*m(3) + w4*f(4)*m(4)),
#	Eq(1, w0 + w1 + w2 + w3 + w4)
#]
#
#print(solve(eqs, (w0, w1, w2, w3, w4)))
#
val = (w0*f(0) + w1*f(1) + w2*f(2) + w3*f(3) + w4*f(4) - rho)**2
print(val)

val = integrate(val, u)
print(val)

val = integrate(val, v)
print(val)

val = integrate(val, rho)
print(val)

val = integrate(val, RT)
print(val)

val = integrate(val, D0)
print(val)
