Q = 4591
Zq = Zmod(Q)
R.<x> = PolynomialRing(Zq)

m0 = x^81
m1 = x^1440 - 1
e0 = -x^1440 + 1
e1 = x^1440

M = m0 * m1
mp = x^761 - x - 1

print(M)

for i in range(81):
    print(i, (e0 * x^i % M) % mp)
for i in range(1440):
    print(i, (e1 * x^i % M) % mp)
