#!/usr/bin/env python3

from gkmc.gkmc import Process, System, bkl
from gkmc.particles import PBCParticle2D
from gkmc.processes import jumpProcess, destroyProcess, spawnProcess
from gkmc.reactions import recombine
import numpy as np
import matplotlib.pyplot as plt
from gkmc.cellList import CellList

wForAll = 1e13
emJumpH = 0.9
emDestroyH = 1.9
emDestroyH2 = 0.06
jumpLen = 34.4
kb = 8.617e-5
N = 14

system = System()
system.res = []
system.temp = 600
system.h = 1000
system.w = 1000
cellSize = 3

system.particlesH = CellList(system, cellSize)
system.particlesH2 = CellList(system, cellSize)

HParticle = PBCParticle2D(system, system.h, system.w)
H2Particle = PBCParticle2D(system, system.h, system.w)


system.eagerness = 0
system.dirty = []


recombineLen = 2

recombineReaction = recombine(system.dirty, recombineLen, H2Particle, system.particlesH, system.particlesH2, postproc=lambda s: s.dirty.clear())

system.reactions = [recombineReaction]


def rate(v, Em, temp):
    return v * np.exp(-Em / (kb * temp))


spawn = spawnProcess(HParticle, system.particlesH, lambda s, p: s.dirty.append(p))


def spawnrate(system):
    return 1e24 * 1e-20 * system.h * system.w


spawnProc = Process(spawn, spawnrate)


jump = jumpProcess(jumpLen, system.particlesH, lambda s, p: s.dirty.append(p))


def jumpRate(system):
    return len(system.particlesH) * rate(wForAll, emJumpH, system.temp)


jumpProc = Process(jump, jumpRate)


def cleardirty(system, p):
    try:
        system.dirty.remove(p)
    except:
        pass


destroyH = destroyProcess(system.particlesH, cleardirty)


def destroyHRate(system):
    return len(system.particlesH) * rate(wForAll, emDestroyH, system.temp)


destroyHProc = Process(destroyH, destroyHRate)

destroyH2 = destroyProcess(system.particlesH2)


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
