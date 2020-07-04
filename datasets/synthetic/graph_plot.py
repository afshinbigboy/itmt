import numpy as np
import networkx as nx
import matplotlib as mpl
from matplotlib import pyplot as plt


font = {'weight' : 'normal',
        'size'   : 24}

mpl.rc('font', **font)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






# ~~~~~~~~~~~~~~~~~~~~~~~~~~ plot function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def plot_graph(G, pos=2, node_size=500, alpha=.5, with_node_labels=True, 
               with_edge_labels=False, show=True, figsize=(10,6)):
    
    plt.figure(figsize=figsize)
    
    if pos == 0: pos = nx.random_layout(G)
    elif pos == 1: pos = nx.spectral_layout(G)
    elif pos == 2: pos = nx.spring_layout(G)
    elif pos == 3: pos = nx.circular_layout(G)
        
    elabels = nx.get_edge_attributes(G,'weight')
    
    try:
        weights = [d['weight'] for _,_,d in G.edges(data=True)]
        weights = np.log(weights/(np.max(weights)))
    except:
        weights = 1
    
    nx.draw(
        G, 
        pos, 
        with_labels=with_node_labels,
        node_size=node_size,
        width=weights,
        font_size=14,
        font_weight='bold',
        font_color='black',
        alpha=alpha,
        edge_color='gray',
        node_color='cyan'
    )
    if with_edge_labels:
        nx.draw_networkx_edge_labels(G,pos,edge_labels=elabels)
    plt.axis('off')
    if show: 
        plt.show()

# ==================================================================================





# ~~~~~~~~~~~~~~~~~~~~~~ generate < n > distinc colors function ~~~~~~~~~~~~~~~~~~~~

import random
def get_random_color(pastel_factor = 0.5):
    return [(x+pastel_factor)/(1.0+pastel_factor) for x in [random.uniform(0,1.0) for i in [1,2,3]]]
def color_distance(c1,c2):
    return sum([abs(x[0]-x[1]) for x in zip(c1,c2)])
def generate_new_color(existing_colors,pastel_factor = 0.5):
    max_distance = None
    best_color = None
    for i in range(0,100):
        color = get_random_color(pastel_factor = pastel_factor)
        if not existing_colors:
            return color
        best_distance = min([color_distance(color,c) for c in existing_colors])
        if not max_distance or best_distance > max_distance:
            max_distance = best_distance
            best_color = color
    return best_color

def get_colors(n):
    ret = ['red', 'cyan', 'yellow', 'magenta', 'orange', 'black', 'green', 'gray', 'w', 'pink', 'violet']
    ret = []
    colors = []
    for i in range(n):
        color = generate_new_color(colors,pastel_factor = 0.9)
        colors.append(color)
        r, g, b = color
        r,g,b = r*256,g*256,b*256
        r,g,b = int(r), int(g), int(b)
        ret.append("#{0:02x}{1:02x}{2:02x}".format(r,g,b))
    return ret

# ==================================================================================
