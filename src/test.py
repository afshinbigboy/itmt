# import community
import numpy as np
import networkx as nx
import matplotlib as mpl
from matplotlib.pyplot import imshow
from matplotlib import pyplot as plt
import matplotlib.image as mpimg
import pygraphviz
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
import random
import pydoc
from ds import McmcTree as Tree
from utils import ColorPrint as _

import sys
sys.path.append("..")
from datasets.synthetic.generator import TreeGenerator



font = {'weight' : 'normal',
        'size'   : 24}

mpl.rc('font', **font)






### load random data
M = 20
N = 30
ZETA = 1
Gamma = 0.15
alpha = 0.00
beta = 0.00
MR = 0.00

tg = TreeGenerator(
    M = M,
    N = N,
    ZETA = ZETA,
    Gamma = Gamma,
    alpha = alpha,
    beta = beta,
    MR = MR,
)
(gt_E, gt_D, D, gt_T) = tg.generate()
gensNames = list( str(i) for i in range(M) )
print(gensNames)



C_num = D.shape[1]
G_num = D.shape[0]
_.print_warn( 'There is {} cells and {} mutations at {} genes in this dataset.'.format(C_num, G_num, len(gensNames)) )



### fill missed data
def tf(m,c):
    os = len(np.where(D[:,c]==1.))*1.
    zs = len(np.where(D[:,c]==0.))*1.
    return 1. if np.random.rand() < os/(os+zs) else 0.

for m in range(G_num):
    for c in range(C_num):
        if D[m,c] == 3.:
            D[m,c] = tf(m,c)


### Run
dl = list(d for d in D)
root = [n for n,d in gt_T.in_degree() if d==0][0]
T = Tree(gensNames, data_list=dl, root=str(root))
T.set_ground_truth(gt_D, gt_E, gt_T=gt_T)

T.randomize()
# T.plot('T0')

for i in range(3000):
    if T.next():
        break


T.plot_all_results()
