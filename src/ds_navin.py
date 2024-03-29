# import community
import numpy as np
import networkx as nx
import matplotlib as mpl
from matplotlib.pyplot import imshow
from matplotlib import pyplot as plt
import matplotlib.image as mpimg
import graphviz
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
import random
import pydoc

from decimal import Decimal
import json

from utils import ColorPrint as _
import matplotlib.patches as mpatches


class McmcTree():
    def __set_plt(self):
        font = {
            'weight' : 'normal',
            'size'   : 8
        }
        mpl.rc('font', **font)
        

    def __init__(self, nodes, D, data_list=None, alpha=0.001, beta=0.1, root='DNM3', name='My Tree', logfile=None):
        self.__T = nx.DiGraph()
        
        if data_list:
            NWD = list( (n, dict(data=d)) for n, d in zip(nodes, data_list) )
            self.__T.add_nodes_from(NWD)
            self.num_cells = len(data_list[0])
            self.num_genes = len(nodes)
        else:
            self.__T.add_nodes_from(nodes)
        self.__best_T = self.__T.copy()
        self.alpha = alpha
        self.beta = beta
        self.root = root
        self.gene_names = nodes
        self.name = name
        self.__rho = 6
        self.D = D
        self.step = 0
        self.logfile = logfile

        self.last_A_hash = None
        self.last_A = None
        self.last_E_hash = None
        self.last_E = None

        self.gt_D = None
        self.gt_E = None
        self.gt_T = None

        self.__random_errors = []
        self.__best_errors = []
        self.__errors = []
        
        self.run_data = []

        self.SSNODES = set(['DNM3','ITGAD','ITGAD','BTLA','BTLA','PAMK3','PAMK3', 'FCHSD2','FCHSD2','LSG1', 
                            'LSG1','DCAF8L1','DCAF8L1','PIK3CA','PIK3CA','CASP3','CASP3','TRIM58','TRIM58','TCP11',
                            ])
    
#         _.print_info('\nNew mcmc tree named:', name)

    
    def set_rho(self, rho):
        self.__rho = rho
    
    def set_ground_truth(self, gt_D, gt_E, gt_T=None):
        self.gt_D = gt_D
        self.gt_E = gt_E
        if gt_T:
            self.gt_T = gt_T
        return



    def __plot_matrix(self, M, title, cmap=None, labels=None, xlabel='cells', ylabel='genes', xticks=None, yticks=None, filename=None):
        t = 1 ## alpha value
        if not cmap:
            cmap = {0:[1,1,0.95,t], 1:[0.5,0.5,0.8,t], -1:[0.8,0.5,0.5,t]}
        if not labels:
            labels = {0:'0', 1:'1', -1:'-1'}
        arrayShow = np.array([[cmap[i] for i in j] for j in M])    
        ## create patches as legend
        patches =[mpatches.Patch(color=cmap[i],label=labels[i]) for i in cmap]
        plt.imshow(arrayShow, interpolation="nearest")
        plt.legend(handles=patches, loc=2, borderaxespad=-6)

        m, n = M.shape
        plt.yticks(range(m), [i for i in self.gene_names] if not yticks else yticks)
        plt.xticks(range(n), ['cell %d'%i for i in range(n)] if not xticks else xticks)
        plt.xticks(rotation=60)

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if title:
            plt.title(title)
        if filename:
            plt.savefig(filename)
            plt.close()
        else:
            return plt

    def __plot_E(self, E, title, filename=None):
        self.__plot_matrix(E, title, filename=filename)

    def __plot_D(self, D, title, filename=None):
        labels = {0:'0', 1:'1', 3:'missed'}
        cmap = {0:[1,1,0.95,1], 1:[0.5,0.5,0.8,1], 3:[0.75,0.75,0.875,0.5]}
        self.__plot_matrix(D, title, filename=filename, labels=labels, cmap=cmap)

    def __plot_DmE(self, DmE, title, filename=None):
        labels = {0:'true', 1:r'false positive (\alpha)', -1:r'false negetive (\beta)', 2:'accepted (1)', 3:'accepted (0)'}
        cmap = {0:[1,1,0.95,1], 1:[0.5,0.5,0.8,1], -1:[0.8,0.5,0.5,1], 2:[0.5,0.5,1,0.25], 3:[1,0.5,0.5,0.25]}
        self.__plot_matrix(DmE, title, filename=filename, labels=labels, cmap=cmap)
    
    def __plot_EmGtE(self, DmE, title, filename=None):
        labels = {0:'true', 1:r'false positive (\alpha)', -1:r'false negetive (\beta)'}
        self.__plot_matrix(DmE, title, filename=filename, labels=labels)

    def __plot_DmGtE(self, DmE, title, filename=None):
        labels = {0:'true', 1:r'false positive (\alpha)', -1:r'false negetive (\beta)', 2:'accepted (1)', 3:'accepted (0)'}
        cmap = {0:[1,1,0.95,1], 1:[0.5,0.5,0.8,1], -1:[0.8,0.5,0.5,1], 2:[0.5,0.5,1,0.25], 3:[1,0.5,0.5,0.25]}
        self.__plot_matrix(DmE, title, filename=filename, labels=labels, cmap=cmap)

    def __plot_A(self, A, title, filename=None):
        xlabel = 'Attached cells to every node(gene) in the Mutation Tree'
        xticks = ['cell %s'%i for i in self.gene_names]
        self.__plot_matrix(A, title, filename=filename, xlabel=xlabel, xticks=xticks)

    def __plot_charts(self, filename=None):
        plt.plot(self.__errors, 'r', label='Accepted Error') # accepted errors
        plt.plot(self.__random_errors, 'k', label='Random Error') # random errors
        # plt.plot(self.enrgs) # best errors
        plt.legend()
        plt.xlabel('Iteration')
        plt.ylabel('Loss')
        plt.title('Changing loss after {} step'.format(self.step))
        if filename:
            plt.savefig(filename)
        # plt.show()
        return
    
    def __plot_pm_charts(self, filename=None, p=5):
        new_acc_errors = []
        new_random_errors = []
        for i, t in enumerate(np.array(self.run_data)):
            rnd = np.random.rand()
            if t[-1] > rnd/2:
                new_acc_errors.append(self.__errors[i])
                new_random_errors.append(t[-2])
                
        plt.plot(new_acc_errors[1:], 'r', label='Accepted Error') # accepted errors
        plt.plot(new_random_errors, 'k', label='Random Error') # random errors
        # plt.plot(self.enrgs) # best errors
        plt.legend()
        plt.xlabel('Iteration')
        plt.ylabel('Loss')
        plt.title('Changing loss after {} step'.format(len(new_random_errors)))
        if filename:
            plt.savefig(filename+"_pm")
        # plt.show()
        return
    
    def plot_eng_chart(self):
        self.__plot_pm_charts()

    def save_mats(self, prefix_dir=None):
        D = self.__get_D()
        A = self.__get_A()
        E = self.__get_E()
        best_T = self.get_best_tree()
        gt_D = self.gt_D
        gt_E = self.gt_E
        gt_T = self.gt_T
        np.savetxt('{}D_n{}_m{}_a{}_b{}.mat'.format(prefix_dir, self.num_cells, self.num_genes, self.alpha, self.beta), D)
        np.savetxt('{}A_n{}_m{}_a{}_b{}.mat'.format(prefix_dir, self.num_cells, self.num_genes, self.alpha, self.beta), A)
        np.savetxt('{}E_n{}_m{}_a{}_b{}.mat'.format(prefix_dir, self.num_cells, self.num_genes, self.alpha, self.beta), E)
        np.savetxt('{}gtD_n{}_m{}_a{}_b{}.mat'.format(prefix_dir, self.num_cells, self.num_genes, self.alpha, self.beta), gt_D)
        np.savetxt('{}gtE_n{}_m{}_a{}_b{}.mat'.format(prefix_dir, self.num_cells, self.num_genes, self.alpha, self.beta), gt_E)


    def plot_all_results(self, plot_pm=False, p=5):
        D = self.__get_D()
        best_T = self.get_best_tree()
        A = self.__get_A(T=best_T)
        E = self.__get_E(T=best_T)
        gt_D = self.gt_D
        gt_E = self.gt_E
        gt_T = self.gt_T
        
        self.__set_plt()
        if gt_D is not None and gt_E is not None:
            
            plt.figure(figsize=(30, 50))

            plt.subplot2grid((8, 2), (0, 0), rowspan=3)
            pdot = nx.drawing.nx_pydot.to_pydot(best_T)
            pdot.write_png('example.png')
            img = mpimg.imread('example.png')
            # plt.figure(figsize=(30,40))
            plt.imshow(img)
            plt.title('best tree with error:{:0.3f}'.format(self.get_best_error()))
            plt.axis('off')

            plt.subplot2grid((8, 2), (0, 1), rowspan=3)
            pdot = nx.drawing.nx_pydot.to_pydot(gt_T)
            pdot.write_png('example.png')
            img = mpimg.imread('example.png')
            # plt.figure(figsize=(40,40))
            plt.imshow(img)
            plt.title('truth tree')
            plt.axis('off')


            plt.subplot2grid((8, 3), (3, 0))
            self.__plot_D(D, 'D')

            plt.subplot2grid((8, 3), (3, 1))
            self.__plot_D(gt_D, 'Ground Truth D (GT_D)')

            plt.subplot2grid((8, 3), (3, 2))
            self.__plot_DmE(D - gt_D, 'D - GT_D')


            plt.subplot2grid((8, 3), (4, 0))
            self.__plot_E(E, 'E')

            plt.subplot2grid((8, 3), (4, 1))
            self.__plot_E(gt_E, 'Ground Truth E (GT_E)')

            plt.subplot2grid((8, 3), (4, 2))
            self.__plot_EmGtE(E - gt_E, 'E - GT_E')


            plt.subplot2grid((8, 3), (5, 0))
            self.__plot_DmE(D-E, 'D - E')

            plt.subplot2grid((8, 3), (5, 1))
            self.__plot_EmGtE(gt_D-gt_E, 'GT_D - GT_E')

            plt.subplot2grid((8, 3), (5, 2))
            self.__plot_DmE(D - gt_E, 'D - GT_E')


            plt.subplot2grid((8, 3), (6, 0), colspan=3)
            if plot_pm:
                self.__plot_pm_charts(p=p)
            else:
                self.__plot_charts()

        else:
            plt.figure(figsize=(20, 20))
            plt.subplot(211)
            self.__plot_charts()
            plt.subplot(245)
            self.__plot_D(D, 'D')
            plt.subplot(246)
            self.__plot_E(E, 'E')
            plt.subplot(247)
            self.__plot_DmE(D-E, 'D-E')
            plt.subplot(248)
            self.__plot_A(A, 'A')
            plt.title('Sharing x per column, y per row')

        plt.savefig('benchmark')
        plt.close()
        



    def __data(self, n):
        if Decimal(nx.__version__) < 2.4:
            return self.__T.node[n]['data']
        else:
            return self.__T.nodes[n]['data']
            
    def parent(self, n):
        return list(self.__T.predecessors(n))[0]


    
    def __swap_nodes(self,):
        n1, n2, p1, p2 = '', '', '', ''
        T = self.__T.copy()
        while True:
            n1 = random.choice(list(T.nodes))
            n2 = random.choice(list(T.nodes))

            if n1 in self.SSNODES or n2 in self.SSNODES:
                if n1 in self.SSNODES and n2 in self.SSNODES:
                    continue

            if not n1 == n2:
                if n1 != self.root and n2 != self.root:
                    break
        p1 , p2 = self.parent(n1), self.parent(n2)
        n1_childs = list(T.neighbors(n1))
        n2_childs = list(T.neighbors(n2))

        T.remove_edge(p1, n1)
        T.remove_edge(p2, n2)

        if p1 == n2 :
            T.add_edge(p2, n1)
            T.add_edge(n1, n2)
            # _.print_info('n1:{}, n2_childs:{}, n2:{}'.format(n1, n2_childs, n2))
        elif p2 == n1 :
            T.add_edge(p1, n2)
            T.add_edge(n2, n1)
            # _.print_info('n2:{}, n1_childs:{}, n1:{}'.format(n2, n1_childs, n1))
        else:
            T.add_edge(p1, n2)
            T.add_edge(p2, n1)

        for n1_child in n1_childs:
            if not n1_child == n2:
                T.remove_edge(n1, n1_child)
                T.add_edge(n2, n1_child)

        for n2_child in n2_childs:
            if not n2_child == n1:
                T.remove_edge(n2, n2_child)
                T.add_edge(n1, n2_child)

        self.__finilizing_step(T, 'swap_nodes')

    def __swap_subtrees(self,):
        n1, p1 , n2, p2 = '', '', '', ''
        PTN = []
        while True:
            T = self.__T.copy()
            n1 = random.choice(list(T.nodes))
            if n1 in self.SSNODES:
                continue
            if not n1 == self.root:
                p1 = self.parent(n1)
                if not p1 == self.root:
                    T.remove_edge(p1, n1)
                    WCC = list(list(wc) for wc in nx.algorithms.components.weakly_connected_components(T))
                    PTN = WCC[0] if p1 in WCC[0] else WCC[-1]
                    if len(PTN) > 2:
                        break
       
        while True:
            n2 = random.choice(PTN)
            if not n2 == self.root and not n2 == p1:
                break

        p2 = self.parent(n2)
        # T = self.__T.copy()
        # T.remove_edge(p1, n1)      
        T.remove_edge(p2, n2)

        T.add_edge(p2, n1)
        if n2 in list(nx.algorithms.ancestors(self.__T, p1)):
            T.add_edge(n1, n2)
        else:
            T.add_edge(p1, n2)

                    
        self.__finilizing_step(T, 'swap_subtrees')

    def __prune_reattach(self,):
        T = self.__T.copy()
        n1, p , n2 = '', '', ''
        while True:
            n1 = random.choice(list(T.nodes))
            if n1 in self.SSNODES:
                continue
            if not n1 == self.root:
                p = self.parent(n1)
                break    
        T.remove_edge(p, n1)
        WCC = list(list(wc) for wc in nx.algorithms.components.weakly_connected_components(T))
        PTN = WCC[0] if p in WCC[0] else WCC[-1]
        while True:
            n2 = random.choice(PTN)
            if n1 in self.SSNODES:
                continue
            break
        T.add_edge(n2, n1)
        
        # print(n1, n2, )
        self.__finilizing_step(T.copy(), 'prune_reattach')
        
    def __finilizing_step(self, new_T, method):
        # _.print_warn(hash(new_T))
        new_error = self.__calc_tree_error(new_T)
        if new_error > 0:
            acc_prob = min(1, (self.__errors[-1]/new_error))
        else:
            acc_prob = 1

        if method == 'swap_nodes':
            if acc_prob<1: acc_prob = acc_prob/100
        else:
            if acc_prob<1: acc_prob = acc_prob/100

        self.__random_errors.append(new_error)
        
        if np.random.rand() <= acc_prob:
            self.__T = new_T.copy()
            self.__errors.append(new_error)
        else:
            self.__errors.append(self.__errors[-1])
    
        if new_error < self.__best_errors[-1]:
            self.__best_errors.append(new_error)
            self.__best_T = new_T.copy()
        else:
            self.__best_errors.append(self.__best_errors[-1])
        
        self.run_data.append(
            [self.step, new_error, acc_prob]
        )
        
        if self.step % 25 == 0:
            print(
                ',\t'.join([
                    'step:{:3d}'.format(self.step),
                    'mode:{}'.format(method),
                    'new_error:{:.2f}'.format(new_error),
                    'last_error:{:.2f}'.format(self.__errors[-1]),
                    'acc_prob:{:0.3f}'.format(acc_prob),
                ])
            )
        if self.logfile:
            self.logfile.write(
                ',\t'.join([
                'step:{:3d}'.format(self.step),
                'mode:{}'.format(method),
                'new_error:{:.2f}'.format(new_error),
                'last_error:{:.2f}'.format(self.__errors[-1]),
                'acc_prob:{:0.3f}\n'.format(acc_prob),
                ])
            )


       
    def get_tree(self,):
        return self.__T.copy()

    def get_best_tree(self,):
        return self.__best_T.copy()

    def get_best_error(self):
        return self.__best_errors[-1]

    def get_errors(self):
        return self.__errors




    def __initialize_params(self):
        self.step = 0
        error = self.__calc_tree_error()
        self.__errors = [error]
        self.__random_errors = [error]
        self.__best_errors = [error]
        self.__best_T = self.__T.copy()
 
    def set_edges(self, edges, remove_edges=False):
        if remove_edges:
            self.__T.remove_edges_from(list(self.__T.edges))
        self.__T.add_edges_from(edges)
        self.__initialize_params()
     
    def randomize(self,):
        
        self.__errors = []
        self.__random_errors = []
        self.__best_errors = []
        
        self.__T.remove_edges_from(list(self.__T.edges))
        nodes = list(self.__T.nodes())
        np.random.shuffle(nodes)
        nodes.remove(self.root)
        nodes.append(self.root)
        for i, n in enumerate(nodes[:-1]):
            p = random.choice(nodes[i+1:])
            self.__T.add_edge(p, n)
        self.__initialize_params()
       
    
    
    def plot_probs(self):
        plt.figure(figsize=(30,10))
        plt.hist(self.probs)
        plt.title('after {} step'.format(self.step))
        # plt.show()
        return
    

    def plot_best_T_just(self, filename='pure_best_T.png'):
        pdot = nx.drawing.nx_pydot.to_pydot(best_T)
        pdot.write_png(filename)
        img = mpimg.imread(filename)
        # plt.figure(figsize=(30,40))
        plt.imshow(img)
        plt.title('best tree with error:{:0.3f}'.format(self.get_best_error()))
        plt.axis('off')
        plt.savefig(filename)


    def plot_best_T(self, filename=None):
        T = self.__best_T
        
        # p = nx.drawing.nx_pydot.to_pydot(T)
        # p.write_png('example.png')
        # img = mpimg.imread('example.png')
        # plt.figure(figsize=(20,20))
        # plt.imshow(img)
        # plt.title('best tree with error:{}'.format(self.get_best_error()))
        # plt.axis('off')

        # tree with attached cells...
        D = self.__get_D()
        A = self.__get_A(T)

        attachments = dict((g, []) for g in self.gene_names)

        for j in range(D.shape[1]):
            d = D[:,j]
            probs = []
            for e in A.T:
                prob = self.__calc_sample_distance(e, d)
                probs.append(prob)
            max_idx = np.argmax(probs)
            attachments[self.gene_names[max_idx]].append('{}|{:0.1f}'.format(j, probs[max_idx]))
        
        
        attached_T = T.copy()
        for gene, samples in attachments.items():
            for s in samples:
                attached_T.add_edge(gene, s)

        pdot = nx.drawing.nx_pydot.to_pydot(attached_T)
        for i, node in enumerate(pdot.get_nodes()):
            node_name = str(node)[:-1]
            print(node_name)
            if '|' in node_name:
                # node.set_label('s%s'%node_name.split()[-1][:-1])
                node.set_shape('egg')
                # node.set_fillcolor('#db8625')
                node.set_color('red')

        pdot.write_png('example.png')
        img = mpimg.imread('example.png')
        plt.figure(figsize=(20,20))
        plt.imshow(img)
        plt.title('Navin best tree with error:{:0.3f}'.format(self.__calc_tree_error(T=T)))
        plt.axis('off')

        if filename:
            plt.savefig(filename)
            edges = tuple(T.edges())
            with open(filename+'.edges', 'w') as f:
                f.write('['+','.join(['({}, {})'.format(u, v) for u, v in edges])+']')
        # plt.show()
        return



    def __get_D(self):
        return self.D
        # if self.D:
        #     return D
        # else:
        #     genes = self.__best_T.nodes(data=True)
        #     D = np.zeros([self.num_genes, self.num_cells], dtype=np.int)
        #     for i, g in enumerate(genes):
        #         D[i, :] = np.array(g[1]['data'])
        #     return D

    def __get_A(self, T=None):
        if not T:
            T = self.__T.copy()

        if hash(tuple(T.edges)) == self.last_A_hash:
            return self.last_A
        
        A = np.eye(self.num_genes)
        for j, gene_name in enumerate(self.gene_names):
            ancestors = nx.algorithms.ancestors(T, gene_name)
            for i in range(len(self.gene_names)):
                if self.gene_names[i] in ancestors:
                    A[i, j] = 1
        
        self.last_A_edges = hash(tuple(T.edges))
        self.last_A = A
        return A

    def __calc_sample_distance(self, from_sample, to_sample):
        ''' FromSamle(E) ~> ToSample(D) '''
        prob = 0
        a, b = self.alpha, self.beta
        for f,t in zip(from_sample, to_sample):
            if int(t) == 3: continue # missed data
            p = f*t*(1-b) + (1-f)*(1-t)*(1-a) + f*(1-t)*b + (1-f)*t*a
            prob += np.log(p)
        return prob
    
    def __get_E(self, T=None):

        if not T:
            T = self.__T.copy()
        
        if hash(tuple(T.edges)) == self.last_E_hash:
            return self.last_E

        D = self.__get_D()
        A = self.__get_A(T=T)
        E = np.zeros_like(D)
        for j in range(D.shape[1]):
            d = D[:,j]
            probs = []
            for e in A.T:
                prob = self.__calc_sample_distance(e, d)
                probs.append(prob)
            max_idx = np.argmax(probs)
            E[:,j] = A[:, max_idx]
            # print(probs[max_idx])

        self.last_E_hash = hash(tuple(T.edges))
        self.last_E = E
        return E

    def __calc_tree_error(self, T=None):
        D = self.__get_D()
        E = self.__get_E(T=T)
        
        DPE = D+E
        DME = D-E
        
        tn_cnt = self.num_cells*self.num_genes - np.count_nonzero(DPE)
        tp_cnt = self.num_cells*self.num_genes - np.count_nonzero(DPE-2)
        
        fn_cnt = self.num_cells*self.num_genes - np.count_nonzero(DME+1)
        fp_cnt = self.num_cells*self.num_genes - np.count_nonzero(DME-1)

#         DmE = D-E
#         DmE[np.where(np.abs(DmE) > 1)] = 0
#         D_and_E = D*E
#         D_and_E[np.where(np.abs(D_and_E) > 1)] = 1

#         ze_cnt = np.ones(self.num_cells)*self.num_genes - np.count_nonzero(DmE, 0)
#         tp_cnt = np.count_nonzero(D_and_E-1, 0)
#         tn_cnt = ze_cnt - tp_cnt
#         fp_cnt = np.count_nonzero(DmE-1, 0) - ze_cnt
#         fn_cnt = np.count_nonzero(DmE+1, 0) - ze_cnt
        ep = 0
        prob =  np.log((1-self.beta )**(tp_cnt+1)/1 + ep) +\
                np.log((1-self.alpha)**(tn_cnt+1)/1 + ep) +\
                np.log(self.beta **(fn_cnt+1)/1 + ep) +\
                np.log(self.alpha**(fp_cnt+1)/1 + ep)

#         error = np.sum(self.beta*fn_cnt + self.alpha*fp_cnt)*10
        # error = np.abs(np.linalg.norm(E) - np.linalg.norm(D))
        # error = np.sum(np.abs( (DmE**1)[:] ) )

        # _.print_info('Error:', error)
        error = -1*prob/400
        return error**1.3
    
    
    def calc_curr_t_error(self, T=None):
        self.__plot_charts()
        return self.__calc_tree_error(T)
    


    base_step = 8
    face_step = 8
    def next(self,):
        if not self.__best_errors[-1]:
            return True

        self.step += 1

        if self.step % (self.base_step+self.face_step) < self.base_step:
            rnd = np.random.rand()
            if rnd < 0.5:
                self.__prune_reattach()
            else:
                self.__swap_subtrees()

        else:
            self.__swap_nodes()
        
        return False
            


#         pos = graphviz_layout(T, prog='twopi', args='-Nfontsize=10 -Nwidth=".2" -Nheight=".2" -Nmargin=0 -Gfontsize=8')
#         plt.figure(figsize=(20,20))
#         nx.draw_networkx_nodes(T, pos, T.nodes, node_size=9000, alpha=0.5, 
#                                node_color="k", font_size=34, font_weight='bold')
#         nx.draw_networkx_labels(T, pos, font_size=14, font_weight='bold', font_color='w')
#         nx.draw_networkx_edges(T, pos, edge_color='gray', alpha=0.925, width=1)
#         plt.axis('equal')
#         plt.axis('off')
#         plt.show()