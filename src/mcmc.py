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


font = {'weight' : 'normal',
        'size'   : 24}

mpl.rc('font', **font)




### load Navin's data
D = np.loadtxt('../../dataset/real/dataNavin.csv', delimiter=' ')
gensNames = np.loadtxt('../../dataset/real/dataNavin.geneNames', dtype=np.str)

D.shape
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



### SCITE Tree with Navin's data
SCITE_Navin_Tree = nx.DiGraph()
edges = [
    ('DNM3','ITGAD'), ('ITGAD','BTLA'), ('BTLA','PAMK3'), ('PAMK3', 'FCHSD2'), ('FCHSD2','LSG1'), 
    ('LSG1','DCAF8L1'), ('DCAF8L1','PIK3CA'), ('PIK3CA','CASP3'), ('CASP3','TRIM58'), ('TRIM58','TCP11'),
    ('TCP11','MARCH11'), ('MARCH11','DUSP12'), ('DUSP12','PPP2RE'), ('PPP2RE','ROPN1B'), ('ROPN1B','PITRM1'),
    ('PITRM1','FBN2'), ('FBN2','PRDM9'), ('FBN2','GPR64'), ('PRDM9','CABP2'), ('PRDM9','ZEHX4'), 
    ('PRDM9','H1ENT'), ('PRDM9', 'WDR16'), ('CABP2', 'TRIB2'), ('ZEHX4','DKEZ'), ('WDR16','GLCE'), 
    ('GLCE','CALD1'), ('CABP2','C15orf23'), ('CABP2','CNDP1'), ('CNDP1','CXXX1'), ('CNDP1','c1orf223'), 
    ('CXXX1','FUBP3'), ('c1orf223','TECTA'), ('GPR64','MUTHY'), ('MUTHY','SEC11A'), ('SEC11A','KIAA1539'), 
    ('SEC11A','RABGAP1L'), ('RABGAP1L','ZNE318'), ('KIAA1539','FGFR2'), ('FGFR2','PLXNA2')
]

dl = list(d for d in D)
SNT = Tree(gensNames, data_list=dl, name='Paper Tree')
SNT.set_edges(edges, remove_edges=True)

_.print_bold( 'SCITE Navis\'s Tree Error:', SNT.get_best_error() )
# SNT.plot('SCITE')




### Run
dl = list(d for d in D)
T = Tree(gensNames, data_list=dl)
T.randomize()
edges = [('PIK3CA', 'c1orf223'),('PIK3CA', 'TCP11'),('DNM3', 'ITGAD'),('TRIM58', 'DUSP12'),('DCAF8L1', 'FCHSD2'),('DCAF8L1', 'GLCE'),('FBN2', 'PPP2RE'),('FCHSD2', 'LSG1'),('CASP3', 'PITRM1'),('CASP3', 'RABGAP1L'),('ITGAD', 'DCAF8L1'),('PPP2RE', 'ROPN1B'),('LSG1', 'PIK3CA'),('ROPN1B', 'MARCH11'),('BTLA', 'FBN2'),('DUSP12', 'BTLA'),('MARCH11', 'CASP3'),('MARCH11', 'CALD1'),('TCP11', 'TRIM58'),('TCP11', 'CNDP1'),('PITRM1', 'PRDM9'),('PITRM1', 'ZEHX4'),('PITRM1', 'MUTHY'),('PITRM1', 'CXXX1'),('PRDM9', 'CABP2'),('PRDM9', 'PAMK3'),('MUTHY', 'FGFR2'),('GPR64', 'TRIB2'),('SEC11A', 'C15orf23'),('SEC11A', 'PLXNA2'),('C15orf23', 'H1ENT'),('DKEZ', 'FUBP3'),('FUBP3', 'WDR16'),('WDR16', 'GPR64'),('RABGAP1L', 'TECTA'),('RABGAP1L', 'ZNE318'),('RABGAP1L', 'KIAA1539'),('RABGAP1L', 'SEC11A'),('GLCE', 'DKEZ')]
T.set_edges(edges, remove_edges=True)
# T.plot('T0')


# T.set_edges(edges, remove_edges=True)
for i in range(1):
    T.next()
# T.plot_run('energy_chart')
T.plot_best_T('best_tree')
SNT.plot_best_T('paper_tree')


T.plot_all_results()
# T.calc_tree_likelihood()
# SNT.calc_tree_likelihood()
# T.plot_best_T()
