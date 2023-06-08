#!/usr/bin/env python3


from math import sqrt, ceil, log


Relams = [15.0 + 20.0 * i for i in range(16)]


L = 1.0
ltarget = 0.2 * L
etaodx = 0.5


composites = [
    2**i * 3**j * 5**k
    for i in range(ceil(log(1024) / log(2)))
    for j in range(3)
    for k in range(1)
]


for Rel in Relams:
    eta = ltarget * (Rel / sqrt(15.0))**(-3.0 / 2)
    dxmax = eta / etaodx
    Nmin = L / dxmax
    Nuse = min(N for N in composites if N > Nmin)
    print(f"Rel = {Rel} requires N >= {Nmin:.2f}, use N = {Nuse}")


