#!/usr/bin/env python3

import numpy as np


def spawnProcess(system, spawn_particle, store, postproc=None):
    rx = np.random.rand() * system.w
    ry = np.random.rand() * system.h
    p = spawn_particle(rx, ry)
    store.add(p)
    if postproc:
        postproc(system, p)


def spawn3DProcess(system, spawn_particle, store, postproc=None):
    rx = np.random.rand() * system.w
    ry = np.random.rand() * system.h
    rz = np.random.rand() * system.d
    p = spawn_particle(rx, ry, rz)
    store.add(p)
    if postproc:
        postproc(system, p)


def destroyProcess(system, store, postproc=None):
    if len(system.particlesH):
        p = store[np.random.randint(len(store))]
        store.remove(p)
        if postproc:
            postproc(system, p)


def jumpProcess(system, jumpLen, store, postproc=None):
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


def jump3DProcess(system, jumpLen, store, postproc=None):
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


def emit3DProcess(system, emitDist, emitParticle, store, postproc=None):
    p = store[np.random.randint(len(store))]
    phi = np.random.uniform(0, np.pi*2)
    costheta = np.random.uniform(-1, 1)

    theta = np.arccos(costheta)
    dx = np.sin(theta) * np.cos(phi) * emitDist
    dy = np.sin(theta) * np.sin(phi) * emitDist
    dz = np.cos(theta) * emitDist
    p2 = emitParticle(p.x + dx, p.y + dy, p.z + dz)
    store.add(p2)
    p.changeSize(p.size - 1)
    if postproc:
        postproc(system, p, p2)
