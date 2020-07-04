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



class TreeGenerator():
    '''
    Inputs: (M, N, ZETA, Gamma, alpha, beta, MR, save_images, save_trees, save_mats)
        --------------------------------------
        M       :   num of genes (mutations)
        --------------------------------------
        N       :   num of samples (cells)
        --------------------------------------
        ZETA    :   homogeness of tree
        --------------------------------------
        Gamma   :   merge genes
        --------------------------------------
        alpha   :   ~ P(D=1|E=0)
        --------------------------------------
        beta    :   ~ P(D=0|E=1)
        --------------------------------------
        MR      :   missing ratio
        --------------------------------------

    Outputs: (E, D, Dm)
        ------------------------------------------------
        E       :   Mutation-cell matrix without errors
        ------------------------------------------------
        D       :   Mutation-cell matrix with errors
        ------------------------------------------------
        Dm      :   <D> with missed data
        ------------------------------------------------
    '''

    def __init__(self,
        M, N, 
        ZETA=1,Gamma=0.15,
        alpha=0.1,beta=0.08,MR=0.05,
        save_images=False, save_trees=False, save_mats=False):

        self.M=M
        self.N=max(M, N)
        self.ZETA=ZETA
        self.Gamma=Gamma
        self.alpha=alpha
        self.beta=beta
        self.MR=MR
        self.save_images=save_images
        self.save_trees=save_trees
        self.save_mats=save_mats


    def generate(self,):
        self.Tree = dict()
        self.cnt = 2

        xk = np.arange(self.M+1)
        name_k = [str(i) for i in xk]
        wk = np.ones(self.M+1, dtype=np.float128)

        while True:
            xk, wk, name_k, u, v = self.do_next(xk, wk, name_k)
            self.cnt+=1
            if len(xk) < 2:
                break

        T = nx.DiGraph(self.Tree)
        T_leaves = [x for x in T.nodes() if T.out_degree(x)==0 and T.in_degree(x)==1]
        T.remove_nodes_from(T_leaves)

        t = np.arange(self.M)
        np.random.shuffle(t)
        t = dict((i, j) for i,j in zip(T.nodes(), t))
        T = nx.relabel_nodes(T, t)
        raw_T = T.copy()


        if self.save_images:
            pdot = nx.drawing.nx_pydot.to_pydot(T)
            pdot.write_png('mt_N{}_M{}_Z{}_G{}_step1.png'.format(self.N,self.M,self.ZETA,self.Gamma))

        A = int(np.floor(self.Gamma*self.M))
        if A:
            for i in range(A):
                p, c = random.sample(T.edges(),1)[0]
                for child in T.successors(c):
                    T.add_edge(p,child)        
                T.remove_node(c)
                T = nx.relabel_nodes(T, {p: '{} . {}'.format(p,c)})

        if self.save_images:
            pdot = nx.drawing.nx_pydot.to_pydot(T)
            pdot.write_png('mt_N{}_M{}_Z{}_G{}_step2.png'.format(self.N,self.M,self.ZETA,self.Gamma))


        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ## Add cells

        Mutaions_T = T.copy()
        mutaion_nodes = Mutaions_T.nodes()

        cells = np.array(['cell %d'%i for i in range(self.N)])
        np.random.shuffle(cells)

        for n in mutaion_nodes:
            T.add_edge(n, cells[0])
            cells = cells[1:]

        for cell in cells:
            node = random.sample(mutaion_nodes, 1)[0]
            T.add_edge(node, cell)

        if self.save_images:
            pdot = nx.drawing.nx_pydot.to_pydot(T)
            for i, node in enumerate(pdot.get_nodes()):
                node_name = str(node)[:-1]
                if 'cell' in node_name:
                    node.set_label('s%s'%node_name.split()[-1][:-1])
                    node.set_shape('egg')
                    node.set_fillcolor('#db8625')
                    node.set_color('red')
            pdot.write_png('mt_N{}_M{}_Z{}_G{}_step3.png'.format(self.N,self.M,self.ZETA,self.Gamma))
        if self.save_trees:
            nx.write_gpickle(T, 'Tree_N{}_M{}_Z{}_G{}.gpickle'.format(self.N,self.M,self.ZETA,self.Gamma))

        # nx.write_graphml_lxml(T, 'Tree_N{}_M{}_Z{}_G{}.graphml'.format(N,M,ZETA,Gamma))
        
        
        ## ========================================================
        ## ~~~~~~~~~~~~~~~~~~~~~ Tree to E ~~~~~~~~~~~~~~~~~~~~~~~~
        ## ========================================================
        E = np.zeros([self.M, self.N])
        root = [n for n,d in T.in_degree() if d==0][0]
        E[root, :] = 1
        for n in range(self.N):
            path = list(nx.all_simple_paths(T, root, 'cell %d'%n))[0]
            for g in path[:-1]:
                try:
                    E[int(g),n] = 1
                except:
                    gs = g.split(' . ')
                    for g in gs:
                        E[int(g),n] = 1

        if self.save_images:
            font = {'weight' : 'normal',
                    'size'   : 8}
            mpl.rc('font', **font)

            plt.imshow(E, cmap='GnBu', interpolation="nearest")
            plt.yticks(range(E.shape[0]), ['gene %d'%i for i in range(self.M)])
            plt.xticks(range(E.shape[1]), ['cell %d'%i for i in range(self.N)])
            plt.xticks(rotation=60)
            plt.xlabel('Genes-Cells Matrix (E)')
            plt.title(r'Parameters: ($N={},M={},\zeta={},\gamma={}$)'.format(self.N,self.M,self.ZETA,self.Gamma))
            plt.savefig('mt_N{}_M{}_Z{}_G{}_E.png'.format(self.N,self.M,self.ZETA,self.Gamma))
            plt.close()


        ## ========================================================
        ## ~~~~~~~~~~~~~~~~~~~~~ E to D ~~~~~~~~~~~~~~~~~~~~~~~~
        ## ========================================================
        D = E.copy()

        nz_idxs = np.nonzero(E)
        z_idxs = np.nonzero(E-1)

        z_rnds = np.random.rand(len(z_idxs[0]))
        nz_rnds = np.random.rand(len(nz_idxs[0]))

        z_rnds  = [1 if i < self.alpha  else 0 for i in z_rnds ]
        nz_rnds = [0 if i < self.beta else 1 for i in nz_rnds]

        D[nz_idxs] = nz_rnds
        D[z_idxs]  = z_rnds

        if self.save_images:
            plt.imshow(D, cmap='GnBu', interpolation="nearest")
            plt.yticks(range(E.shape[0]), ['gene %d'%i for i in range(self.M)])
            plt.xticks(range(E.shape[1]), ['cell %d'%i for i in range(self.N)])
            plt.xticks(rotation=60)
            plt.xlabel('Noisy Genes-Cells Matrix (D)')
            plt.title(r'Parameters: ($N={},M={},\zeta={},\gamma={},\alpha={},\beta={}$)'.format(self.N,self.M,self.ZETA,self.Gamma,self.alpha,self.beta))
            plt.savefig('mt_N{}_M{}_Z{}_G{}_a{}_b{}_D.png'.format(self.N,self.M,self.ZETA,self.Gamma,self.alpha,self.beta))
            plt.close()

            ## first you need to define your color map and value name as a dic
            t = 1 ## alpha value
            cmap = {0:[1,1,0.95,t], 1:[0.5,0.5,0.8,t], -1:[0.8,0.5,0.5,t]}
            labels = {0:'true', 1:'false positive', -1:'false negetive'}
            arrayShow = np.array([[cmap[i] for i in j] for j in D-E])    
            ## create patches as legend
            patches =[mpatches.Patch(color=cmap[i],label=labels[i]) for i in cmap]

            plt.imshow(arrayShow, interpolation="nearest")
            plt.legend(handles=patches, loc=2, borderaxespad=-6)
            plt.yticks(range(E.shape[0]), ['gene %d'%i for i in range(self.M)])
            plt.xticks(range(E.shape[1]), ['cell %d'%i for i in range(self.N)])
            plt.xticks(rotation=60)
            plt.xlabel('D-E')
            plt.title(r'Parameters: ($N={},M={},\zeta={},\gamma={},\alpha={},\beta={},m_r={}$)'.format(self.N,self.M,self.ZETA,self.Gamma,self.alpha,self.beta,self.MR))
            plt.savefig('mt_N{}_M{}_Z{}_G{}_a{}_b{}_DmE.png'.format(self.N,self.M,self.ZETA,self.Gamma,self.alpha,self.beta))
            plt.close()



        ## ========================================================
        ## ~~~~~~~~~~~~~~~~~ add missing data ~~~~~~~~~~~~~~~~~~~~~
        ## ========================================================

        Dm = D.copy()

        idxs = np.nonzero(D+1)
        rnds = np.random.rand(self.M,self.N)
        for m in range(self.M):
            for n in range(self.N):
                if rnds[m,n] < self.MR:
                    Dm[m,n] = 3

        if self.save_images:
            ## first you need to define your color map and value name as a dic
            t = 1 ## alpha value
            cmap = {0:[1,1,0.95,t], 1:[0.2,0.2,0.4,t], 3:[0.8,0.5,0.5,t]}
            labels = {0:'0', 1:'1', 3:'missed'}
            arrayShow = np.array([[cmap[i] for i in j] for j in Dm])    
            ## create patches as legend
            patches =[mpatches.Patch(color=cmap[i],label=labels[i]) for i in cmap]

            plt.imshow(arrayShow, interpolation="nearest")
            plt.legend(handles=patches, loc=2, borderaxespad=-6)
            plt.yticks(range(E.shape[0]), ['gene %d'%i for i in range(self.M)])
            plt.xticks(range(E.shape[1]), ['cell %d'%i for i in range(self.N)])
            plt.xticks(rotation=60)
            plt.xlabel('Noisy Genes-Cells Matrix with Missed Data ($D_m$)')
            plt.title(r'Parameters: ($N={},M={},\zeta={},\gamma={},\alpha={},\beta={},m_r={}$)'.format(self.N,self.M,self.ZETA,self.Gamma,self.alpha,self.beta,self.MR))
            plt.savefig('mt_N{}_M{}_Z{}_G{}_a{}_b{}_MR_{}_Dm.png'.format(self.N,self.M,self.ZETA,self.Gamma,self.alpha,self.beta,self.MR))
            plt.close()



        # --------------- save mats ------------------
        if self.save_mats:
            p = 'Parameters: (N={},M={},zeta={},gamma={},alpha={},beta={},m_r={})\n'.format(self.N,self.M,self.ZETA,self.Gamma,self.alpha,self.beta,self.MR)
            np.savetxt('E.csv', E, fmt='%.0f', delimiter=',', header=p)
            np.savetxt('D.csv', D, fmt='%.0f', delimiter=',', header=p)
            np.savetxt('DmE.csv', D-E, fmt='%.0f', delimiter=',', header=p)
            np.savetxt('Dm.csv', Dm, fmt='%.0f', delimiter=',', header=p)


        return (E.astype(int), D.astype(int), Dm.astype(int), raw_T)



    def do_next(self, xk, wk, name_k):
        u, v = weighted_drand(xk, wk, size=2)

        idx_u = np.where(xk==u)[0]
        idx_v = np.where(xk==v)[0]

        w_u = wk[idx_u]
        w_v = wk[idx_v]
        w_uv = (w_u+w_v)/(self.ZETA**0.25)

        nu = name_k[int(idx_u)]
        nv = name_k[int(idx_v)]

        nuv = '{}.{}'.format(nu, nv)
        self.Tree[nuv] = [nu, nv]

        xk = np.delete(xk, [idx_u, idx_v])
        name_k = np.delete(name_k, [idx_u, idx_v])
        wk = np.delete(wk, [idx_u, idx_v])

        xk = np.append(xk, self.M+self.cnt)
        name_k = np.append(name_k, nuv)
        wk = np.append(wk, w_uv)

        return (xk, wk, name_k, u, v)



