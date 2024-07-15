#!/usr/bin/env python3

from gkmc.gkmc import Process, System, bkl, rate
from gkmc.processes import emit3DProcess, jump3DProcess, spawn3DProcess
from gkmc.reactions import recombine3D, absorb, absorbInto
from gkmc.cellList import CellList3D
import numpy as np

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


class HeParticle:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def move(self, dx, dy, dz):
        self.x = (self.x + dx) % system.w
        self.y = (self.y + dy) % system.h
        self.z = (self.z + dz) % system.d

    def __repr__(self):
        return f'He({self.x}, {self.y}, {self.z})'


class HeCluster:
    def __init__(self, x, y, z, size=2):
        self.x = x
        self.y = y
        self.z = z
        self.size = size

    def move(self, dx, dy, dz):
        self.x = (self.x + dx) % system.w
        self.y = (self.y + dy) % system.h
        self.z = (self.z + dz) % system.d

    def changeSize(self, size):
        self.size = size

    def radius(self):
        return pow(8.37 * pow(self.size, 1.02), 0.333333);

    def __repr__(self):
        return f'HeC({self.x}, {self.y}, {self.z}, {self.size})'


def recombineHe(system):
    recombine3D(system, system.particlesHe, minDist, HeCluster, system.particlesHe, system.particlesHeC)



def recombineHeC(system):
    absorbInto(system, system.particlesHeC, lambda p: minDist + p.radius(), system.particlesHe)


def recombineHeC2(system):
    absorb(system, system.particlesHeC, lambda p, p2: minDist + p.radius() + p2.radius(), system.particlesHeC)


system.reactions = [recombineHe, recombineHeC, recombineHeC2]


def spawnHe(system):
    spawn3DProcess(system, HeParticle, system.particlesHe)


spawnProc = Process(spawnHe, 22)


def jumpHe(system):
    jump3DProcess(system, jumpLen, system.particlesHe)


def jumpHeRate(system):
    return len(system.particlesHe) * rate(w, jumpEm, system.temp)


jumpProc = Process(jumpHe, jumpHeRate)


def emitHeC(system):
    emit3DProcess(system, emitDist, HeParticle, system.particlesHeC)


nrg = lambda i: 1.863 * pow(i, 0.8998425)
tm = lambda i: 5e10 * np.exp(-(nrg(i) - nrg(i - 1)) / (kB * system.temp))


def emitHecRate(system):
    return sum(1.0/tm(p.size) for p in system.particlesHeC)


emitProc = Process(emitHeC, emitHecRate)

processes = [spawnProc, jumpProc, emitProc]

bkl(processes, system, 2**16, viz=lambda s: print(len(s.particlesHe), len(s.particlesHeC), s.time))
