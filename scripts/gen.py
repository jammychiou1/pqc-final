#!/usr/bin/env python3
import random

Q = 4591
P = 761

a = [random.randrange(Q) for i in range(P)]
b = [random.randrange(Q) for i in range(P)]

print(' '.join(map(str, a)))
print(' '.join(map(str, b)))
