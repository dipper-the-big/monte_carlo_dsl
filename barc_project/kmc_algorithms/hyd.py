#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

h = 200
w = 200
occ = np.zeros((w, h))
ra = 1
rd = 2
Wa = 0.5
Wd = 1
N = 20

res = []
t = 0
for _ in range(2**N):
    if t > 3:
        break
    rx = np.random.randint(0, w)
    ry = np.random.randint(0, h)
    r = np.random.rand()
    if occ[rx, ry] == 0 and r < Wa:
        occ[rx, ry] = 1
        u = np.random.rand()
        dt = -np.log(u) / (h * w * ra)
        t += dt
    elif occ[rx, ry] == 1 and r < Wd:
        occ[rx, ry] = 0
        u = np.random.rand()
        dt = -np.log(u) / (h * w * rd)
        t += dt
    if _ % 100 == 0:
        theta = sum(occ.flatten()) / (h * w)
        res.append((theta, t))


ideal_theta = lambda t: (ra / (ra + rd)) * (1 - np.exp(-(ra + rd) * t))

# plot res
# plot ideal theta
res = np.array(res)
plt.plot(res[:, 1], res[:, 0])
plt.plot(res[:, 1], [ideal_theta(t) for t in res[:, 1]])
plt.show()
