import scipy as sp
import numpy as np
from scipy import stats
import networkx as nx
from matplotlib import pyplot as plt
import matplotlib as mpl
import random
import matplotlib.image as mpimg
import pygraphviz
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
import matplotlib.patches as mpatches

Dm = np.loadtxt('dataXu.csv', delimiter=' ')
gensNames = np.loadtxt('dataXu.geneNames', dtype=np.str)


font = {'weight' : 'normal',
        'size'   : 8}
mpl.rc('font', **font)

## first you need to define your color map and value name as a dic
t = 1 ## alpha value
cmap = {0:[1,1,0.95,t], 1:[0.2,0.2,0.4,t], 3:[0.8,0.5,0.5,t]}
labels = {0:'0', 1:'1', 3:'missed'}
arrayShow = np.array([[cmap[i] for i in j] for j in Dm])    
## create patches as legend
patches =[mpatches.Patch(color=cmap[i],label=labels[i]) for i in cmap]

plt.imshow(arrayShow, interpolation="nearest")
plt.legend(handles=patches, loc=2, borderaxespad=-5)
plt.yticks(range(Dm.shape[0]), ['%s'%n for n in gensNames])
plt.xticks(range(Dm.shape[1]), ['cell %d'%i for i in range(Dm.shape[1])])
plt.xticks(rotation=60)
plt.xlabel('Noisy Genes-Cells Matrix with Missed Data ($D_m$)')
plt.title("Xu's Data")
plt.savefig('Xu_Dm.png')
plt.close()
