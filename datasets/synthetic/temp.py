import scipy as sp
import numpy as np
from scipy import stats
import networkx as nx
from matplotlib import pyplot as plt
import matplotlib as mpl
from graph_plot import plot_graph as plg
import random
import matplotlib.image as mpimg
import pygraphviz
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
import matplotlib.patches as mpatches


T = nx.read_edgelist(path=, delimiter=":")
T = nx.read_gpickle("Tree_N30_M20_Z1_G0.15.edgelist.gpickle")
pdot = nx.drawing.nx_pydot.to_pydot(T)
for i, node in enumerate(pdot.get_nodes()):
    node_name = str(node)[:-1]
    if 'cell' in node_name:
        node.set_label('s%s'%node_name.split()[-1][:-1])
        node.set_shape('egg')
        node.set_fillcolor('#db8625')
        node.set_color('red')
pdot.write_png('reloaded_tree.png')

