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



CELL_RANGE = np.linspace(5, 30, num=3, dtype=np.int)

### load random data
M = 20
N = 40
ZETA = 1
Gamma = 0.15
alpha = 0.01
beta = 0.01
MR = 0.005

prefix_dir = '../outputs/logs/'
for N in CELL_RANGE:
    for M in np.linspace(N, N+N//2, min(10, N-1), dtype=np.int):

        step_num = 10
        _.print_info( 'There is {} cells and {} mutations.'.format(N, M) )

        while True:    
            # try:
                # with open('n{}_m{}_z{}_g{}_a{}_b{}_mr{}.txt'.format(N, M, ZETA, Gamma, alpha, beta, MR), 'w') as file:
                logfile = open('{}n{}_m{}_z{}_g{}_a{}_b{}_mr{}.txt'.format(prefix_dir, N, M, ZETA, Gamma, alpha, beta, MR), 'w')
                logfile.write( 'There is {} cells and {} mutations.\n'.format(N, M) )
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
                
                ### Run
                dl = list(d for d in D)
                root = [n for n,d in gt_T.in_degree() if d==0][0]
                print('ROOT:', root)
                logfile.write( 'ROOT:{}\n'.format(root) )
                T = Tree(gensNames, D, data_list=dl, root=str(root), alpha=alpha, beta=beta, logfile=logfile)
                T.set_ground_truth(gt_D, gt_E, gt_T=gt_T)
                T.randomize()

                for i in range(step_num):
                    if T.next():
                        break
                # T.plot_all_results()
                T.save_mats(prefix_dir)
                logfile.close()
                break

            # except:
            #     logfile.close()
            #     print('Ops...')
            #     print('try again!')



