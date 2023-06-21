#!/usr/bin/env sage

Q = 4591
P = 761

Zq = Zmod(Q)
R.<x> = PolynomialRing(Zq)

m = x ** 1440 - 1

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

def center_lift(val):
  val = int(val)
  return val if val <= (Q - 1) / 2 else val - Q

def print_poly(p):
    coefs = [center_lift(p[i]) for i in range(1440)]
    print(' '.join(map(str, coefs)))

print_poly(a * b % m)
