#!/usr/bin/env python3
import random

Q = 4591
P = 761

a = [random.randint(-(Q - 1) // 2, (Q - 1) // 2) for i in range(P)]
b = [random.randint(-(Q - 1) // 2, (Q - 1) // 2) for i in range(P)]

print(' '.join(map(str, a)))
print(' '.join(map(str, b)))
