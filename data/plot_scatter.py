#!/usr/bin/env python
#
# Created by Tai, Yuan-yen on 11/17/17.
# All rights reserved.
#
import matplotlib.pyplot as plt
import sys

fig, ax = plt.subplots()

colors = ['green', 'blue', 'brown', 'black', 'orange', 'purple', 'gray', 'cyan', 'DarkGreen', 'DarkBlue', 'magenta']
infile = open(sys.argv[1], 'r')
count = 0
for line in infile.readlines():
    line_split = line.split()
    if line_split[0] == 'C':
        print(line_split)
        x = float(line_split[1])
        y = float(line_split[2])
        print(x, ' ',y)
        ax.scatter(x, y, c='red', s=80, alpha=1, edgecolors='none')
    elif count % 40 == 0:
        cluster_id = hash(line_split[0]) % 10000
        x = float(line_split[2])
        y = float(line_split[3])
        ax.scatter(x, y, c=colors[cluster_id % len(colors)], s=10, label=cluster_id, alpha=0.5, edgecolors='none')
    count += 1
infile.close()

ax.grid(True)

plt.axis('equal')
plt.show()