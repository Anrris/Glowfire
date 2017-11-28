#!/usr/bin/env python
#
# Created by Tai, Yuan-yen on 11/17/17.
# All rights reserved.
#
import matplotlib.pyplot as plt
import sys

fig, ax = plt.subplots()

colors = ['yellow', 'red', 'green', 'blue', 'brown', 'black', 'orange', 'purple', 'gray']
infile = open(sys.argv[1], 'r')
count = 0
for line in infile.readlines():
    if count % 100 == 0:
        line_split = line.split()
        cluster_id = int(line_split[0])
        x = float(line_split[2])
        y = float(line_split[3])

        ax.scatter(x, y, c=colors[cluster_id % len(colors)], s=10, label=cluster_id, alpha=0.3, edgecolors='none')
    count += 1
infile.close()

ax.grid(True)

plt.axis('equal')
plt.show()