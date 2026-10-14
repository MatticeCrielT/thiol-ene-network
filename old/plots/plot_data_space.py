import sys, os, re
import numpy as np
import signac
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm
import gsd, gsd.hoomd 
import networkx as nx 
from collections import defaultdict
from matplotlib.pyplot import cm
from mpl_toolkits.mplot3d import axes3d

ax = plt.figure().add_subplot()

project = signac.get_project()

u = project.detect_schema() 
print(u)


for job in project:
    c = job.sp['chain_side_reaction_probability']
    e = job.sp['crosslinker_percent']

    if job.isfile('equi.gsd'):
        color='red'
    else: 
        color='grey'

    if job.isfile('polymerize.gsd'):
        color='green'

    if job.isfile('run.gsd'):
        color='blue'
    if c==0:
        print(job.id,job.sp)

    ax.scatter(c,e,color=color)

ax.set_xlabel('chain_side_reaction_probability')
ax.set_ylabel('crosslinker_percent')

plt.show()