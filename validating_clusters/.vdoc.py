# type: ignore
# flake8: noqa
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
import os
# local vs UGER
if os.path.expanduser('~') in ["/Users/youyun", "/Users/youyunzheng"]: 
    # in a local mac, the home directory is usuaully at '/Users/[username]'
    workdir = os.path.expanduser('~')+"/Documents/HMS/PhD/beroukhimlab/broad_mount/"
else:
    # in dipg or uger, the home directory is usuaully at '/home/unix/[username]'
    workdir = "/xchip/beroukhimlab/"

import pandas as pd
import numpy as np
import networkx as nx
# import torch_geometric
import sklearn
import matplotlib.pyplot as plt
from tqdm import tqdm
from joblib import Parallel, delayed
import multiprocessing
import pickle
import re
import scipy.sparse as sp
from sklearn.model_selection import StratifiedKFold
import umap.umap_ as umap
import csv
import math
import torch
import torch.nn as nn
import torch.nn.functional as F
torch.manual_seed(123)
import numpy as np
np.random.seed(123)
import time
#
#
#
#
#

"""Adapted from https://github.com/weihua916/powerful-gnns/blob/master/util.py"""

class S2VGraph(object):
    def __init__(self, g, label, node_tags=None, node_features=None):
        '''
            g: a networkx graph
            label: an integer graph label
            node_tags: a list of integer node tags
            node_features: a torch float tensor, one-hot representation of the tag that is used as input to neural nets
            edge_mat: a torch long tensor, contain edge list, will be used to create torch sparse tensor
            neighbors: list of neighbors (without self-loop)
        '''
        self.label = label
        self.g = g
        self.node_tags = node_tags
        self.neighbors = []
        self.node_features = 0
        self.edge_mat = 0
        self.max_neighbor = 0

def create_s2v(graph_file):
    G = pickle.load(open(graph_file, 'rb'))
    graph_label = re.sub('_[0-9]*$','',list(G.nodes)[0])
    # relabel nodes to integer
    G = nx.relabel_nodes(G, {n: int(re.sub('.*_', '',n)) for n in G.nodes})
    # use fake label of 0 for everything
    return S2VGraph(G, label = 0), graph_label

def load_data_from_nx(graph_files):
    num_cores = multiprocessing.cpu_count()
    g_list, label_list = zip(*Parallel(n_jobs = num_cores)(
        delayed(create_s2v)(graph_file) for graph_file in graph_files if \
            # pickle.load(open(graph_file, 'rb')).number_of_nodes() > 9 and
            # pickle.load(open(graph_file, 'rb')).number_of_nodes() <= 30 and
            pickle.load(open(graph_file, 'rb')).number_of_edges() > 0
    ))
        
    return g_list, label_list
#
#
#
#
#
#
#
# graph_embeddings = pd.read_csv(workdir + "youyun/complexSV/data/TCGA/graph_embedding/UGformerV2_UnSup_lr0.0005_bs16_ep30_ff1024_nn5_sd512_do0.5_nl2_nt2_graph_embeddings.csv",index_col=0)
graph_embeddings = pd.read_csv(workdir + "youyun/complexSV/data/TCGA/graph_embedding/UGformerV2_UnSup_lr0.0005_bs16_ep60_ff1024_nn20_sd512_do0.5_nl3_nt2_graph_embeddings.csv",index_col=0)
ct_overlap = pd.read_csv(workdir + "youyun/complexSV/data/validation/cluster_overlaps_w_ct.csv")
#
#
#
#
#
graph_dir = workdir + "youyun/complexSV/data/TCGA/SV/nx/"
nx_files = [graph_dir+f+'_nx_graph.txt' for f in graph_embeddings.index]
print(len(nx_files))
graphs, labels = load_data_from_nx(nx_files)
#
#
#
#
#
# UMAP dimensionality reduction
reducer = umap.UMAP(
    n_neighbors = 50, min_dist=0.05, n_components=2, n_epochs = 500
)
embedding = reducer.fit_transform(graph_embeddings)
#
#
#
# Plotting
plt.figure(figsize=(10, 6))
plt.scatter(
    embedding[:, 0], embedding[:, 1], 
    c=[np.log10(graph.g.number_of_nodes()) for graph in graphs], 
    cmap='Spectral', s=5, alpha = 0.5
)
plt.title('UMAP Visualization of Graph Embeddings \n Colored by log10 SV Count')
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.colorbar()
plt.show()
#
#
#
feature_names = [
    'svclass_DEL', 'svclass_DUP','svclass_TRA', 'svclass_h2hINV', 'svclass_t2tINV',
    'dist_to_line_breakend1','dist_to_sine_breakend1','dist_to_ltr_breakend1','dist_to_dna_breakend1',
    'dist_to_TAD_breakend1', 'dist_to_telomere_breakend1', 'dist_to_centromere_breakend1',
    'total_cn_breakend1', 'major_cn_breakend1', 'minor_cn_breakend1',
    'total_cn_breakend2', 'major_cn_breakend2', 'minor_cn_breakend2'
    
]
fig, axs = plt.subplots(
    math.ceil(len(feature_names) / 3), 3,
    figsize=(18, math.ceil(len(feature_names) / 2) * 2.5)
)
# need to adjust for the fact that there are three columns
for i, feat_name in enumerate(feature_names):
    axs[i//3,i % 3].scatter(
        embedding[:, 0], embedding[:, 1], 
        c=[
            sum(nx.get_node_attributes(graph.g, feat_name).values()) / graph.g.number_of_nodes()
#             np.std(list(nx.get_node_attributes(graph.g, feat_name).values()))
            for graph in graphs
        ], 
        cmap='Spectral', s=5, alpha = 0.5
    )
    axs[i//3,i % 3].set_title(feat_name)
# colorbar for each panel
fig.colorbar(
    plt.cm.ScalarMappable(
        norm=plt.Normalize(0, 1), cmap='Spectral'
    ), ax=axs, orientation='vertical',
    shrink = 1/math.ceil(len(feature_names) / 3)
)
plt.show()
#
#
#
# Plotting
plt.figure(figsize=(10, 6))
plt.scatter(
    embedding[:, 0], embedding[:, 1], 
    c=[ct_overlap.loc[(ct_overlap['cluster_ID'] == x), 'ct_overlap'] for x in graph_embeddings.index], 
    cmap='Spectral', s=5, alpha = 0.5
)
plt.title('UMAP Visualization of Graph Embeddings \n Colored by log10 SV Count')
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.colorbar()
plt.show()
#
#
#
#
#
#
#
#
#
