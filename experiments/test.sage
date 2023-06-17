Q = 4591
Zq = Zmod(Q)
R.<x> = PolynomialRing(Zq)

m = x ** 761 - x - 1
m1 = x ** 761
m2 = x ** 760 - 1

e1 = - x ** 1520 + 1
e2 = x ** 1520

# print(e1 % m1)
# print(e1 % m2)
# print(e2 % m1)
# print(e2 % m2)
#
# print('===')

# for i in range(761):
#     print((x ** i * e1) % (m1 * m2) % m)
#     print((x ** i * e2) % (m1 * m2) % m)
#
# print('===')

a = 3 * x ** 760 + 2 * x + 7
b = x ** 760 - 2 * x ** 2 - 1
c = a * b % m

c1 = a * b % m1
c2 = a * b % m2

print(a * b)
print(c)
#print((c1 * e1 + c2 * e2) % (m1 * m2))

tmp1 = [Zp(0) for i in range(761)]
tmp2 = [Zp(0) for i in range(761)]

tmp1[0] = c1[0] - c1[1]
for i in range(1, 759):
    tmp1[i] = - c1[i + 1]
tmp1[759] = - c1[0] - c1[760]
tmp1[760] = - c1[0]

tmp2[0] = c2[1]
for i in range(1, 759):
    tmp2[i] = c2[i + 1] + c2[i]
tmp2[759] = c2[0] + c2[760] + c2[759]
tmp2[760] = c2[0] + c2[760]

ans = 0
for i in range(761):
    ans += (tmp1[i] + tmp2[i]) * x ** i
print(ans)
