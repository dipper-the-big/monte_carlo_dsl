#!/usr/bin/env python3

from .gkmc import System, Process, bkl, rate_calc, Species
from .particles import PBCParticle2D, PBCParticle3D, SizeParticle
from .processes import spawnProcess, jumpProcess, destroyProcess, emit3DProcess, jump3DProcess, spawn3DProcess
from .reactions import recombine, recombine3D, absorb, absorbInto
from .cellList import CellList, CellList3D


class Simulation:
    def __init__(self):
        self.system = System()
        self.system.time = 0
        self.system.species = {}
        self.system.reactions = []
        self.processes = []

    def add_species(self, species, particle_type, store_type='list', **params):
        if not isinstance(particle_type, list):
            particle_type = [particle_type]
        if 'PBCParticle2D' in particle_type:
            particle = PBCParticle2D(species, self.system.h, self.system.w)
        if 'PBCParticle3D' in particle_type:
            particle = PBCParticle3D(species, self.system.h, self.system.w, self.system.d)
        if 'SizeParticle' in particle_type:
            particle = SizeParticle(particle)

        if store_type == 'list':
            self.system.species[species] = Species(particle, [])
        if store_type == 'celllist':
            cl = CellList(self.system, params['cellSize'])
            self.system.species[species] = Species(particle, cl)
        if store_type == 'celllist3D':
            cl = CellList3D(self.system, params['cellSize'])
            self.system.species[species] = Species(particle, cl)

    def add_process(self, process, species, rate=None, nParticles=False, **params):
        if not rate:
            rate = rate_calc(params['w'], params['em'], self.system.temp)

        if nParticles:
            rate_single = rate

            def rate_all(sys):
                return len(self.system.species[species].store) * rate_single

            rate = rate_all

        if process == 'spawn':
            postproc = lambda sys, p: sys.species[species].dirty.append(p)
            proc = spawnProcess(self.system.species[species].particle, self.system.species[species].store, postproc=postproc)
            self.processes.append(Process(proc, rate))
        if process == 'jump':

            postproc = lambda sys, p: sys.species[species].dirty.append(p)
            proc = jumpProcess(params['jumpLen'], self.system.species[species].store, postproc=postproc)
            self.processes.append(Process(proc, rate))
        if process == 'destroy':
            proc = destroyProcess(self.system.species[species].store)
            self.processes.append(Process(proc, rate))
            # self.processes.append(Process(process, rate))
        if process == 'emit3D':
            postproc = lambda sys, p: sys.species[species].dirty.append(p)
            proc = emit3DProcess(params['emitDist'], self.system.species[species].particle, self.system.species[species].store, self.system.species[species].store, postproc=postproc)
            self.processes.append(Process(proc, rate))
        if process == 'jump3D':
            postproc = lambda sys, p: sys.species[species].dirty.append(p)
            proc = jump3DProcess(params['jumpLen'], self.system.species[species].store, postproc=postproc)
            self.processes.append(Process(proc, rate))
        if process == 'spawn3D':
            postproc = lambda sys, p: sys.species[species].dirty.append(p)
            proc = spawn3DProcess(self.system.species[species].particle, self.system.species[species].store, postproc=postproc)
            self.processes.append(Process(proc, rate))

    def add_reaction(self, reaction, **params):
        if reaction == 'recombine':
            postproc = lambda sys: sys.species[params['combineFrom']].dirty.clear()
            reac = recombine(self.system.species[params['combineFrom']].dirty,
                             params['recombineLen'],
                             self.system.species[params['combineFrom']].particle,
                             self.system.species[params['combineFrom']].store,
                             self.system.species[params['combineTo']].store,
                             postproc=postproc)
            self.system.reactions.append(reac)
        if reaction == 'recombine3D':
            postproc = lambda sys: sys.species[params['combineFrom']].dirty.clear()
            reac = recombine3D(self.system.species[params['combineFrom']].dirty,
                               params['recombineLen'],
                               self.system.species[params['combineFrom']].particle,
                               self.system.species[params['combineFrom']].store,
                               self.system.species[params['combineTo']].store,
                               postproc=postproc)
            self.system.reactions.append(reac)
        if reaction == 'absorb':
            postproc = lambda sys: sys.species[params['species']].dirty.clear()
            reac = absorb(self.system.species[params['species']].dirty,
                          params['absorbDist'],
                          self.system.species[params['species']].store,
                          postproc=postproc)
            self.system.reactions.append(reac)
        if reaction == 'absorbInto':
            postproc = lambda sys: sys.species[params['absorb_into']].dirty.clear()
            reac = absorbInto(self.system.species[params['absorb_into']].dirty,
                              params['absorbDist'],
                              self.system.species[params['absorb']].store,
                              postproc=postproc)
            self.system.reactions.append(reac)

    def bkl(self, nSteps, record=None, viz=None):
        bkl(self.processes, self.system, nSteps, record=record, viz=viz)
