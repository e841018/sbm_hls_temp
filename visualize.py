# %%
import os
import numpy as np
from matplotlib import pyplot as plt

def parse_line(line):
    return [float(num) for num in line.split(' | ')[1:-1]]

# %%
os.system('g++ sbm.cpp sbm_tb.cpp && a')

x, p = [], []
i = 0
with open('log.txt') as f:
    while f:
        line_x = f.readline()
        if line_x == '':
            break
        line_p = f.readline()
        x.append(parse_line(line_x))
        p.append(parse_line(line_p))
x = np.array(x).T
p = np.array(p).T
N, T = x.shape

fig, axs = plt.subplots(N, 1, figsize=(6, 9))
for a, ax in enumerate(axs):
    ax.plot(x[a])
    ax.plot(p[a])
    ax.set_ylim(-1, 1)
plt.legend(['x', 'p'])

# %%
