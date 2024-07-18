#!/usr/bin/env python3

from gkmc.simulation import Simulation
import numpy as np
import matplotlib.pyplot as plt

N = 20
cellSize = 3

sim = Simulation()
sim.system.h = 1000
sim.system.w = 1000
sim.system.temp = 600
sim.system.res = []


sim.add_species('H', 'PBCParticle2D', store_type='celllist', cellSize=cellSize)
sim.add_species('H2', 'PBCParticle2D', store_type='celllist', cellSize=cellSize)
sim.add_process('spawn', 'H', 1e24 * 1e-20 * sim.system.h * sim.system.w)
sim.add_process('jump', 'H', w=1e13, em=0.9, jumpLen=34.4, nParticles=True)
sim.add_process('destroy', 'H', w=1e13, em=1.9, nParticles=True)
sim.add_process('destroy', 'H2', w=1e13, em=0.06, nParticles=True)
sim.add_reaction('recombine', combineFrom='H', combineTo='H2', recombineLen=2)


def record(system, i):
    if i % 100 == 0:
        system.res.append((len(system.species['H'].store), system.time))


def viz(system):
    res = np.array(system.res)
    plt.plot(res[:, 1][100:], res[:, 0][100:])
    # x axis time, y axis number of H particles
    plt.xlabel('Time')
    plt.ylabel('Number of H particles')
    plt.yscale('log')
    plt.xscale('log')
    plt.show()


sim.bkl(2**N, record=record, viz=viz)
