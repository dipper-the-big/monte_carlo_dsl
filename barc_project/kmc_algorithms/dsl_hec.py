#!/usr/bin/env python3

from gkmc.gkmc import Process, System, bkl, rate_calc, PBCParticle3D, SizeParticle
from gkmc.processes import emit3DProcess, jump3DProcess, spawn3DProcess
from gkmc.reactions import recombine3D, absorb, absorbInto
from gkmc.cellList import CellList3D
from gkmc.simulation import Simulation
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

sim = Simulation()
sim.system.temp = temperature
sim.system.w = boxDim
sim.system.h = boxDim
sim.system.d = boxDim

cellSize = 100

sim.add_species('He', 'PBCParticle3D', store_type='celllist', cellSize=cellSize)
sim.add_species('HeC', ['PBCParticle3D', 'SizeParticle'], store_type='celllist', cellSize=cellSize)

# system.particlesHe = CellList3D(system, cellSize)
# system.particlesHeC = CellList3D(system, cellSize)


# HeParticle = PBCParticle3D('He', system.h, system.w, system.d)
# HeCluster = SizeParticle(PBCParticle3D('HeC', system.h, system.w, system.d))

sim.add_reaction('recombine3D', combineFrom='He', recombineLen=minDist, combineTo='HeC')
sim.add_reaction('absorbInto', absorb_into='HeC', absorb='He', absorbDist=lambda p: minDist + p.radius())
sim.add_reaction('absorb', species='HeC', absorbDist=lambda p, p2: minDist + p.radius() + p2.radius())

# recombineHe = recombine3D(sim.system.species['He'].store, minDist, sim.system.species['HeC'].particle, sim.system.species['He'].store, sim.system.species['HeC'].store)

# recombineHeC = absorbInto(sim.system.species['HeC'].store, lambda p: minDist + p.radius(), sim.system.species['He'].store)
# recombineHeC2 = absorb(sim.system.species['HeC'].store, lambda p, p2: minDist + p.radius() + p2.radius(), sim.system.species['HeC'].store)

# sim.system.reactions = [recombineHe, recombineHeC, recombineHeC2]

spawnHe = spawn3DProcess(sim.system.species['He'].particle, sim.system.species['He'].store)
spawnProc = Process(spawnHe, 22)

jumpHe = jump3DProcess(jumpLen, sim.system.species['He'].store)


def jumpHeRate(system):
    return len(sim.system.species['He'].store) * rate_calc(w, jumpEm, sim.system.temp)


jumpProc = Process(jumpHe, jumpHeRate)

emitHeC = emit3DProcess(emitDist, sim.system.species['He'].particle, sim.system.species['HeC'].store, sim.system.species['He'].store)

nrg = lambda i: 1.863 * pow(i, 0.8998425)
tm = lambda i: 5e10 * np.exp(-(nrg(i) - nrg(i - 1)) / (kB * sim.system.temp))


def emitHecRate(system):
    return sum(1.0/tm(p.size) for p in sim.system.species['HeC'].store)


emitProc = Process(emitHeC, emitHecRate)

processes = [spawnProc, jumpProc, emitProc]
print(sim.system.reactions)
bkl(processes, sim.system, 2**18, viz=lambda s: print(len(s.species['He'].store), len(s.species['HeC'].store), s.time))
