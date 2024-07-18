#!/usr/bin/env python3

def PBCParticle2D(species, h, w):
    class PBCParticle2D:
        def __init__(self, x, y):
            self.x = x
            self.y = y

        def move(self, dx, dy):
            self.x = (self.x + dx) % w
            self.y = (self.y + dy) % h

        def __repr__(self):
            return f'{species} Particle({self.x}, {self.y})'
    return PBCParticle2D


def PBCParticle3D(species, h, w, d):
    class PBCParticle3D:
        def __init__(self, x, y, z):
            self.x = x
            self.y = y
            self.z = z

        def move(self, dx, dy, dz):
            self.x = (self.x + dx) % w
            self.y = (self.y + dy) % h
            self.z = (self.z + dz) % d

        def __repr__(self):
            return f'{species} Particle({self.x}, {self.y})'
    return PBCParticle3D


def SizeParticle(Particle):
    class SizeParticle(Particle):
        def __init__(self, *args, size=2):
            super().__init__(*args)
            self.size = size

        def changeSize(self, size):
            self.size = size

        def radius(self):
            return pow(8.37 * pow(self.size, 1.02), 0.333333);

        def __repr__(self):
            return super().__repr__() + f' size={self.size}'
    return SizeParticle
