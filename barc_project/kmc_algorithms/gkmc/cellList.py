#!/usr/bin/env python3

import numpy as np


class CellList:
    def __init__(self, system, cellSize):
        self.cellSize = cellSize
        self.pList = []
        self.cells = np.empty((system.h//cellSize, system.w//cellSize), dtype=object)
        for i in range(system.h//cellSize):
            for j in range(system.w//cellSize):
                self.cells[i, j] = []

    def add(self, p):
        self.cells[int(p.y//self.cellSize), int(p.x//self.cellSize)].append(p)
        self.pList.append(p)

    def remove(self, p):
        try:
            self.cells[int(p.y//self.cellSize), int(p.x//self.cellSize)].remove(p)
            self.pList.remove(p)
        except ValueError:
            print(f'particle {p} not found')
            for i, row in enumerate(self.cells):
                for j, cell in enumerate(row):
                    for pp in cell:
                        if p == pp:
                            print(f'found in cell {(i, j)} in stead of {(int(p.y//self.cellSize), int(p.x//self.cellSize))}')
                            break
            else:
                print('not found anywhere')

    def neighbors(self, p):
        # periodic boundary conditions
        i = int(p.y//self.cellSize)
        j = int(p.x//self.cellSize)
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
