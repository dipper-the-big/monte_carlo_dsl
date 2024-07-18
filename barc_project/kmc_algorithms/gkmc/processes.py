#!/usr/bin/env python3

import numpy as np


def spawnProcess(spawn_particle, store, postproc=None):
    def proc(system):
        rx = np.random.rand() * system.w
        ry = np.random.rand() * system.h
        p = spawn_particle(rx, ry)
        store.add(p)
        if postproc:
            postproc(system, p)
    return proc


def spawn3DProcess(spawn_particle, store, postproc=None):
    def proc(system):
        rx = np.random.rand() * system.w
        ry = np.random.rand() * system.h
        rz = np.random.rand() * system.d
        p = spawn_particle(rx, ry, rz)
        store.add(p)
        if postproc:
            postproc(system, p)
    return proc


def destroyProcess(store, postproc=None):
    def proc(system):
        if len(store):
            p = store[np.random.randint(len(store))]
            store.remove(p)
            if postproc:
                postproc(system, p)
    return proc


def jumpProcess(jumpLen, store, postproc=None):
    def proc(system):
        i = np.random.randint(len(store))
        p = store[i]
        direction = np.random.rand() * 2 * np.pi
        dx = np.cos(direction) * jumpLen
        dy = np.sin(direction) * jumpLen
        store.remove(p)
        p.move(dx, dy)
        store.add(p)
        if postproc:
            postproc(system, p)
    return proc


def jump3DProcess(jumpLen, store, postproc=None):
    def proc(system):
        phi = np.random.uniform(0, np.pi*2)
        costheta = np.random.uniform(-1, 1)

        theta = np.arccos(costheta)
        dx = np.sin(theta) * np.cos(phi) * jumpLen
        dy = np.sin(theta) * np.sin(phi) * jumpLen
        dz = np.cos(theta) * jumpLen
        p = store[np.random.randint(len(store))]
        store.remove(p)
        p.move(dx, dy, dz)
        store.add(p)
        if postproc:
            postproc(system, p)
    return proc


def emit3DProcess(emitDist, emitParticle, emmiter_store, emmitee_store, postproc=None):
    def proc(system):
        p = emmiter_store[np.random.randint(len(emmiter_store))]
        phi = np.random.uniform(0, np.pi*2)
        costheta = np.random.uniform(-1, 1)

        theta = np.arccos(costheta)
        dx = np.sin(theta) * np.cos(phi) * emitDist
        dy = np.sin(theta) * np.sin(phi) * emitDist
        dz = np.cos(theta) * emitDist
        p2 = emitParticle(p.x + dx, p.y + dy, p.z + dz)
        emmitee_store.add(p2)
        p.changeSize(p.size - 1)
        if postproc:
            postproc(system, p, p2)
    return proc
