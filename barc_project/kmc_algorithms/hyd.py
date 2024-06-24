#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

h = 200
w = 200
occ = []
ra = 1
rd = 2
N = 20

res = []
t = 0

for _ in range(2**N):
    if t > 3:
        break

    rates = [((h * w) - len(occ)) * ra, len(occ) * rd]
    cummRates = []
    c = 0
    for i, r in enumerate(rates):
        c += r
        cummRates.append(c)

    uQ = np.random.rand() * cummRates[-1]
    for i, r in enumerate(cummRates):
        if r > uQ:
            procindex = i
            break

    if procindex == 0:
        rx = np.random.rand() * w
        ry = np.random.rand() * h
        occ.append((rx, ry))
        u = np.random.rand()
        dt = -np.log(u) / cummRates[-1]
        t += dt
    else:
        if len(occ):
            i = np.random.randint(0, len(occ))
            occ.pop(i)
            u = np.random.rand()
            dt = -np.log(u) / cummRates[-1]
            t += dt

    if _ % 100 == 0:
        theta = len(occ) / (h * w)
        res.append((theta, t))


ideal_theta = lambda t: (ra / (ra + rd)) * (1 - np.exp(-(ra + rd) * t))

# plot res
# plot ideal theta
res = np.array(res)
plt.plot(res[:, 1], res[:, 0])
plt.plot(res[:, 1], [ideal_theta(t) for t in res[:, 1]])
plt.show()
