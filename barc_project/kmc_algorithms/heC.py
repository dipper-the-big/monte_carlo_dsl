#!/usr/bin/env python3

from gkmc.gkmc import Process, System, bkl, rate_calc
from gkmc.particles import PBCParticle3D, SizeParticle
from gkmc.processes import emit3DProcess, jump3DProcess, spawn3DProcess
from gkmc.reactions import recombine3D, absorb, absorbInto
from gkmc.cellList import CellList3D
import numpy as np
import matplotlib.pyplot as plt

spawnRate = 22
jumpEm = 1.257
jumpW = 2.0e12
jumpLen = 9.0
emitDist = 10.0
minDist = 5.0
boxDim = 500
temperature = 800
kB = 8.617e-5
w = 2e12

system = System()
system.temp = temperature
system.w = boxDim
system.h = boxDim
system.d = boxDim

cellSize = 100


system.particlesHe = CellList3D(system, cellSize)
system.particlesHeC = CellList3D(system, cellSize)


HeParticle = PBCParticle3D('He', system.h, system.w, system.d)
HeCluster = SizeParticle(PBCParticle3D('HeC', system.h, system.w, system.d))

recombineHe = recombine3D(system.particlesHe, minDist, HeCluster, system.particlesHe, system.particlesHeC)

recombineHeC = absorbInto(system.particlesHeC, lambda p: minDist + p.radius(), system.particlesHe)
recombineHeC2 = absorb(system.particlesHeC, lambda p, p2: minDist + p.radius() + p2.radius(), system.particlesHeC)

system.reactions = [recombineHe, recombineHeC, recombineHeC2]

spawnHe = spawn3DProcess(HeParticle, system.particlesHe)
spawnProc = Process(spawnHe, 22)

jumpHe = jump3DProcess(jumpLen, system.particlesHe)


def jumpHeRate(system):
    return len(system.particlesHe) * rate_calc(w, jumpEm, system.temp)


jumpProc = Process(jumpHe, jumpHeRate)

emitHeC = emit3DProcess(emitDist, HeParticle, system.particlesHeC, system.particlesHe)

nrg = lambda i: 1.863 * pow(i, 0.8998425)
tm = lambda i: 5e10 * np.exp(-(nrg(i) - nrg(i - 1)) / (kB * system.temp))


def emitHecRate(system):
    return sum(1.0/tm(p.size) for p in system.particlesHeC)


emitProc = Process(emitHeC, emitHecRate)

processes = [spawnProc, jumpProc, emitProc]

def viz(system):
    print([p.size for p in system.particlesHeC], system.time)

bkl(processes, system, 2**20, viz=viz)
