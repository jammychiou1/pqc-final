#!/usr/bin/env sage

Q = 4591
P = 761

Zq = Zmod(Q)
R.<x> = PolynomialRing(Zq)

m = x ** P - x - 1
m1 = x ** P
m2 = x ** (P - 1) - 1

Zq = Zmod(Q)
R.<x> = PolynomialRing(Zq)

a_coefs = list(map(int, input().strip().split()))
b_coefs = list(map(int, input().strip().split()))

a = 0
b = 0

for i in range(P):
    a += a_coefs[i] * x ** i
    b += b_coefs[i] * x ** i

c = a * b % m

def print_poly(p):
    coefs = [p[i] for i in range(P)]
    print(' '.join(map(str, coefs)))


# print_poly(a * b % m1)
# print_poly(a * b % m2)
print_poly(a * b % m)
