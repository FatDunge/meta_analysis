from dataclasses import dataclass

import numpy as np

@dataclass
class Study(object):
    name: str
    method: str
    m1: float
    s1: float
    n1: int
    m2: float
    s2: float
    n2: int

    def __post_init__(self):
        self.s = None
        if self.method == 'cohen_d':
            self.func = self.cohen_d

    def cohen_d(self):
        m1 = self.m1
        s1 = self.s1
        n1 = self.n1
        m2 = self.m2
        s2 = self.s2
        n2 = self.n2
        self.s = np.sqrt(((n1-1)*(s1**2)+(n2-1)*(s2**2))/(n1+n2-2))
        d = (m1 - m2) / s
        return d

    def get_effect_size(self):
        return self.func()
    
    def get_variance(self):
        if not self.s:
            self.func()
        return self.s
