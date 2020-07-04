import scipy as sp
import numpy as np
from scipy import stats
import networkx as nx
from matplotlib import pyplot as plt
import matplotlib as mpl
from graph_plot import plot_graph as plg
import random



def rand_pmf(xk, pk, size=1):
    custm = stats.rv_discrete(name='custm', values=(xk, pk))
    cnt = 0
    while True:
        rs = custm.rvs(size = size)
        if len(set(rs)) == len(rs):
            break
        cnt+=1
    return rs


def weighted_drand(xk, wk, size=1):
    pk = wk/np.sum(wk, dtype=np.float128)
    return rand_pmf(xk, pk, size)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tree = dict()
N = 20
ZETA = 1
cnt = 1


xk = np.arange(N)
name_k = [str(i) for i in xk]
wk = np.ones(N, dtype=np.float128)



def do_next(xk, wk, name_k):

    u, v = weighted_drand(xk, wk, size=2)

    idx_u = np.where(xk==u)[0]
    idx_v = np.where(xk==v)[0]

    w_u = wk[idx_u]
    w_v = wk[idx_v]
    w_uv = (w_u+w_v)/(ZETA**0.25)

    nu = name_k[int(idx_u)]
    nv = name_k[int(idx_v)]

    nuv = '{}.{}'.format(nu, nv)
    Tree[nuv] = [nu, nv]

    xk = np.delete(xk, [idx_u, idx_v])
    name_k = np.delete(name_k, [idx_u, idx_v])
    wk = np.delete(wk, [idx_u, idx_v])

    xk = np.append(xk, N+cnt)
    name_k = np.append(name_k, nuv)
    wk = np.append(wk, w_uv)

    # print(w_u, w_v)
    # print(uv, w_uv)

    return (xk, wk, name_k, u, v)

    



while True:
    xk, wk, name_k, u, v = do_next(xk, wk, name_k)
    cnt+=1
    # print(u,v)
    if len(xk) < 2:
        break


T = nx.DiGraph(Tree)
# plg(T)

# pos = graphviz_layout(T, prog='twopi', args='-Nfontsize=10 -Nwidth=".2" -Nheight=".2" -Nmargin=0 -Gfontsize=8')
# plt.figure(figsize=(10,5))
# nx.draw(T, pos, node_size=1000, alpha=0.5, node_color="cyan", font_size=9, font_weight='bold',with_labels=True)
# # plt.axis('equal')
# plt.show()


import matplotlib.image as mpimg
import pygraphviz
from networkx.drawing.nx_agraph import write_dot, graphviz_layout

shapes = ['box', 'polygon', 'ellipse', 'oval', 'circle', 'egg', 'triangle', 'exagon', 'star', ]
colors = ['blue', 'black', 'red', '#db8625', 'green', 'gray', 'cyan', '#ed125b']
styles = ['filled', 'rounded', 'rounded, filled', 'dashed', 'dotted, bold']

pdot = nx.drawing.nx_pydot.to_pydot(T)
for i, node in enumerate(pdot.get_nodes()):
    try:
        cell_num = int(str(node)[:-1])
        node.set_label("cell %d" % cell_num)
        node.set_shape('egg')
        # node.set_fontcolor('green')
        node.set_fillcolor('#db8625')
        # node.set_style('dashed')
        node.set_color('red')
    except:
        node.set_height('".25"')
        node.set_width('".25"')
        node.set_label('')

pdot.write_png('example.png')
# img = mpimg.imread('example.png')
# plt.imshow(img)
# plt.axis('off')
# plt.show()






## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def next_gen_mut(T, n):
    global mutations

    childs = T.successors(n)
    if not childs:
        return 

    parent_mutations = T.nodes[n]['mutations']
    for child in childs:
        rand_idx = random.randrange(len(mutations))
        t = parent_mutations.copy()
        if T.out_degree(child):
            m = mutations[rand_idx]
            mutations = np.delete(mutations, rand_idx)
            t.append(m)
        T.nodes[child]['mutations'] = t.copy()    

        next_gen_mut(T, child)




M = max(int(np.ceil(len(T.nodes())/2)), 10)

root = [n for n,d in T.in_degree() if d==0][0]
mutations = np.arange(M)

rand_idx = random.randrange(len(mutations))
m = mutations[rand_idx]
mutations = np.delete(mutations, rand_idx)

T.nodes[root]['mutations'] = [m]

next_gen_mut(T, root)



E = np.zeros([M, N])
for n,d in T.out_degree():
    if d == 0: #leaves
        for ms in T.nodes[n]['mutations']:
            E[int(ms), int(n)] = 1

print(mutations)
font = {'weight' : 'normal',
        'size'   : 8}

mpl.rc('font', **font)

plt.imshow(E, cmap='GnBu', interpolation="nearest")
plt.yticks(range(E.shape[0]), ['gene %d'%i for i in range(M)])
plt.xticks(range(E.shape[1]), ['cell %d'%i for i in range(N)])
plt.xticks(rotation=60)
plt.xlabel('Genes-Cells Matrix (E)')

plt.show()

print(E)




