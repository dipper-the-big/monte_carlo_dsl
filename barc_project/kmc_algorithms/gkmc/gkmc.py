#!/usr/bin/env python3

import numpy as np
from tqdm import tqdm

kb = 8.617e-5


def rate_calc(v, Em, temp):
    return v * np.exp(-Em / (kb * temp))


class Species:
    def __init__(self, particle, store):
        self.particle = particle
        self.store = store
        self.dirty = []


class Process:
    def __init__(self, procfunc, rate):
        self.procfunc = procfunc
        self.ratefunc = rate if callable(rate) else lambda x: rate

    def rate(self, system):
        return self.ratefunc(system)

    def execute(self, system):
        self.procfunc(system)


class System(dict):
    def __init__(self, reactions=[], eagerness=0):
        self.time = 0
        self.eagerness = eagerness
        self.elapsed = 0
        self.reactions = reactions

    def __getattr__(self, key):
        return self[key]

    def __setattr__(self, key, value):
        self[key] = value

    def update(self, dtime):
        self.time += dtime
        self.elapsed += dtime
        if self.elapsed > self.eagerness:
            self.elapsed = 0
            for reaction in self.reactions:
                reaction(self)


def bkl(processes, system, steps, record=None, viz=None):
    for i in tqdm(range(steps)):
        rates = [proc.rate(system) for proc in processes]
        # print(rates)
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

        processes[procindex].execute(system)

        u = np.random.rand()
        dt = -np.log(u) / cummRates[-1]
        system.update(dt)

        if record:
            record(system, i)

    if viz:
        viz(system)


