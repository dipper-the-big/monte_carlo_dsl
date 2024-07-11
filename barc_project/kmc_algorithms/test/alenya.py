#!/usr/bin/env python3

from gkmc import Process, System, bkl
import numpy as np
import matplotlib.pyplot as plt
# from collections import namedtuple
from cellList import CellList

wForAll = 1e13
emJumpH = 0.9
emDestroyH = 1.9
emDestroyH2 = 0.06
jumpLen = 34.4
kb = 8.617e-5
N = 17

system = System()
system.res = []
system.temp = 600
system.h = 1000
system.w = 1000
cellSize = 20

system.particlesH = CellList(system, cellSize)
system.particlesH2 = CellList(system, cellSize)


class HParticle:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def move(self, dx, dy):
        self.x = (self.x + dx) % system.w
        self.y = (self.y + dy) % system.h

    def __repr__(self):
        return f'H({self.x}, {self.y})'


class H2Particle:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def move(self, dx, dy):
        self.x = (self.x + dx) % system.w
        self.y = (self.y + dy) % system.h

    def __repr__(self):
        return f'H2({self.x}, {self.y})'

system.eagerness = 0
system.dirty = []


def recombine(system):
    for p in system.dirty:
        for p2 in system.particlesH.neighbors(p):
            if p != p2:
                if np.sqrt((p.x - p2.x)**2 + (p.y - p2.y)**2) < recombineLen:
                    system.particlesH2.add(H2Particle((p.x + p2.x)/2, (p.y + p2.y)/2))
                    system.particlesH.remove(p)
                    system.particlesH.remove(p2)
                    break
    system.dirty = []


system.reactions = [recombine]

recombineLen = 2


def rate(v, Em, temp):
    return v * np.exp(-Em / (kb * temp))


def spawn(system):
    rx = np.random.rand() * system.w
    ry = np.random.rand() * system.h
    p = HParticle(rx, ry)
    system.particlesH.add(p)
    system.dirty.append(p)


def spawnrate(system):
    return 1e24 * 1e-20 * system.h * system.w


spawnProc = Process(spawn, spawnrate)


def jump(system):
    i = np.random.randint(len(system.particlesH))
    p = system.particlesH[i]
    direction = np.random.rand() * 2 * np.pi
    dx = np.cos(direction) * jumpLen
    dy = np.sin(direction) * jumpLen
    system.particlesH.remove(p)
    p.move(dx, dy)
    system.particlesH.add(p)
    system.dirty.append(p)


def jumpRate(system):
    return len(system.particlesH) * rate(wForAll, emJumpH, system.temp)


jumpProc = Process(jump, jumpRate)


def destroyH(system):
    if len(system.particlesH):
        p = system.particlesH[np.random.randint(len(system.particlesH))]
        system.particlesH.remove(p)
        try:
            system.dirty.remove(p)
        except:
            pass


def destroyHRate(system):
    return len(system.particlesH) * rate(wForAll, emDestroyH, system.temp)


destroyHProc = Process(destroyH, destroyHRate)


def destroyH2(system):
    if len(system.particlesH2):
        p = system.particlesH2[np.random.randint(len(system.particlesH2))]
        system.particlesH2.remove(p)


def destroyH2Rate(system):
    return len(system.particlesH2) * rate(wForAll, emDestroyH2, system.temp)


destroyH2Proc = Process(destroyH2, destroyH2Rate)

processes = [spawnProc, destroyHProc, jumpProc, destroyH2Proc]


def record(system, i):
    if i % 100 == 0:
        system.res.append((len(system.particlesH), system.time))


def viz(system):
    res = np.array(system.res)
    plt.plot(res[:, 1][100:], res[:, 0][100:])
    plt.yscale('log')
    plt.xscale('log')
    plt.show()


bkl(processes, system, 2**N, record, viz)
