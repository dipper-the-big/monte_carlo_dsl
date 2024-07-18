#!/usr/bin/env python3

import numpy as np


def recombine(candidates, recombineLen, combineParticle, store_combiner, store_combinee, postproc=None):
    def reac(system):
        for p in candidates:
            for p2 in store_combiner.neighbors(p):
                if p != p2:
                    if np.sqrt((p.x - p2.x)**2 + (p.y - p2.y)**2) < recombineLen:
                        store_combinee.add(combineParticle((p.x + p2.x)/2, (p.y + p2.y)/2))
                        store_combiner.remove(p)
                        store_combiner.remove(p2)
                        break
        if postproc:
            postproc(system)
    return reac


def recombine3D(candidates, recombineLen, combineParticle, store_combiner, store_combinee, postproc=None):
    def reac(system):
        for p in candidates:
            for p2 in store_combiner.neighbors(p):
                if p != p2:
                    if np.sqrt((p.x - p2.x)**2 + (p.y - p2.y)**2 + (p.z - p2.z)**2) < recombineLen:
                        store_combinee.add(combineParticle((p.x + p2.x)/2, (p.y + p2.y)/2, (p.z + p2.z)/2))
                        store_combiner.remove(p)
                        store_combiner.remove(p2)
                        break
        if postproc:
            postproc(system)
    return reac


def absorb(candidates, abs_rad, store, postproc=None):
    def reac(system):
        for p in candidates:
            for p2 in store.neighbors(p):
                if p != p2:
                    if np.sqrt((p.x - p2.x)**2 + (p.y - p2.y)**2 + (p.z - p2.z)**2) < abs_rad(p, p2):
                        if p.size > p2.size:
                            p.changeSize(p.size + p2.size)
                            store.remove(p2)
                        else:
                            p2.changeSize(p.size + p2.size)
                            store.remove(p)
                        break
        if postproc:
            postproc(system)
    return reac


def absorbInto(candidates, abs_rad, store, postproc=None):
    def reac(system):
        for p in candidates:
            for p2 in store.neighbors(p):
                if np.sqrt((p.x - p2.x)**2 + (p.y - p2.y)**2 + (p.z - p2.z)**2) < abs_rad(p):
                    p.changeSize(p.size + 1)
                    store.remove(p2)
                    break
        if postproc:
            postproc(system)
    return reac
