#!/usr/bin/env python3

import numpy as np


class CellList:
    def __init__(self, system, cellSize):
        self.cellSize = cellSize
        self.h = system.h
        self.w = system.w
        self.pList = []
        self.cells = np.empty((int(np.ceil(self.h/cellSize)), int(np.ceil(self.w/cellSize))), dtype=object)
        for i in range(int(np.ceil(self.h/cellSize))):
            for j in range(int(np.ceil(self.w/cellSize))):
                self.cells[i, j] = []

    def add(self, p):
        self.cells[int(p.y/self.cellSize) % self.h, int(p.x/self.cellSize) % self.w].append(p)
        self.pList.append(p)

    def remove(self, p):
        try:
            self.cells[int(p.y/self.cellSize), int(p.x/self.cellSize)].remove(p)
            self.pList.remove(p)
        except ValueError:
            print(f'particle {p} not found')

    def neighbors(self, p):
        # periodic boundary conditions
        i = int(p.y/self.cellSize)
        j = int(p.x/self.cellSize)
        for k in range(-1, 2):
            for l in range(-1, 2):
                for p2 in self.cells[(i + k) % len(self.cells), (j + l) % len(self.cells[0])]:
                    if p != p2:
                        yield p2

    def __iter__(self):
        for row in self.cells:
            for cell in row:
                for p in cell:
                    yield p

    def __getitem__(self, key):
        return self.pList[key]

    def __len__(self):
        return len(self.pList)


class CellList3D:
    def __init__(self, system, cellSize):
        self.cellSize = cellSize
        self.pList = []
        self.h = system.h
        self.w = system.w
        self.d = system.d
        self.cells = np.empty((int(np.ceil(self.h/cellSize)), int(np.ceil(self.w/cellSize)), int(np.ceil(self.d/cellSize))), dtype=object)
        for i in range(int(np.ceil(self.h/cellSize))):
            for j in range(int(np.ceil(self.w/cellSize))):
                for k in range(int(np.ceil(self.d/cellSize))):
                    self.cells[i, j, k] = []

    def add(self, p):
        self.cells[int(p.y/self.cellSize) % self.h, int(p.x/self.cellSize) % self.w, int(p.z/self.cellSize) % self.d].append(p)
        self.pList.append(p)

    def remove(self, p):
        try:
            self.cells[int(p.y/self.cellSize), int(p.x/self.cellSize), int(p.z/self.cellSize)].remove(p)
            self.pList.remove(p)
        except ValueError:
            print(f'particle {p} not found')

    def neighbors(self, p):
        # periodic boundary conditions
        i = int(p.y/self.cellSize)
        j = int(p.x/self.cellSize)
        k = int(p.z/self.cellSize)
        for l in range(-1, 2):
            for m in range(-1, 2):
                for n in range(-1, 2):
                    for p2 in self.cells[(i + l) % len(self.cells), (j + m) % len(self.cells[0]), (k + n) % len(self.cells[0][0])]:
                        if p != p2:
                            yield p2

    def __iter__(self):
        for i in self.cells:
            for j in i:
                for k in j:
                    for p in k:
                        yield p

    def __getitem__(self, key):
        return self.pList[key]

    def __len__(self):
        return len(self.pList)
