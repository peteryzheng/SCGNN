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
import os
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import re
import pickle
from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()
# local vs UGER
if os.path.expanduser('~') in ["/Users/youyun", "/Users/youyunzheng"]: 
    # in a local mac, the home directory is usuaully at '/Users/[username]'
    workdir = os.path.expanduser('~')+"/Documents/HMS/PhD/beroukhimlab/broad_mount/"
else:
    # in dipg or uger, the home directory is usuaully at '/home/unix/[username]'
    workdir = "/xchip/beroukhimlab/"
#
#
#
#
#
#
#
#
#
graph_files = glob.glob(workdir+"/youyun/complexSV/data/TCGA/SV/nx/*.txt")
print("Number of graph files: ", len(graph_files))
# multicore read in graph_files
def read_in_graph(graph_file):
    G = pickle.load(open(graph_file, 'rb'))
    return G
g_list = Parallel(n_jobs = num_cores)(
    delayed(read_in_graph)(graph_file) for graph_file in graph_files
)
#
#
#
#
#
#
#
#
#
n_svs = [G.number_of_nodes() for G in g_list]
print('Quantiles of number of SVs: ', 
    np.quantile(n_svs, np.arange(0, 1.1, 0.1))
)
# summary statistics
pd.Series(n_svs).describe()
# plot both the histogram and the cumulative distribution of the number of SVs on the same graph
plt.figure(figsize=(6, 6))
plt.hist(n_svs, bins=2000, color='gray', alpha=0.5, cumulative=True, label='Cumulative Distribution')
plt.hist(n_svs, bins=2000, color='#80DDF9', alpha=0.9, label='Histogram')
plt.xlabel('Number of SVs')
plt.ylabel('Cumulative Counts')
plt.title('Number of SVs in SV Clusters')
plt.xscale('log')
breaks = [2, 5, 10, 100, 1000, 2074]
plt.xticks(breaks, breaks)
plt.legend()
plt.show()
#
#
#
n_svs_larger_than_5 = [x for x in n_svs if x > 5]
# plot both the histogram and the cumulative distribution of the number of SVs on the same graph
plt.figure(figsize=(6, 6))
plt.hist(n_svs_larger_than_5, bins=2000, color='gray', alpha=0.5, cumulative=True, label='Cumulative Distribution')
plt.hist(n_svs_larger_than_5, bins=2000, color='#80DDF9', alpha=0.9, label='Histogram')
plt.xlabel('Number of SVs')
plt.ylabel('Cumulative Counts')
plt.title('Number of SVs in SV Clusters (SV > 5)')
plt.xscale('log')
breaks = [5, 10, 100, 1000, 2074]
plt.xticks(breaks, breaks)
plt.legend()
plt.show()
#
#
#
#
#
# get the degree of each node in the graph
node_degrees = [[v for k,v in G.degree()] for G in g_list]

# plot the distribution of node degrees in three panels -- 
# 2 SVs, 2-10 SVs, > 10 SVs
fig, ax = plt.subplots(1, 4, figsize=(24, 6))
cutoffs = [5, 50, 500, 2074]
for i, cutoff in enumerate(cutoffs):
    if i == 0:
        plot_list = [x for a, y in enumerate(node_degrees) if n_svs[a] <= cutoff for x in y]
    else:
        plot_list = [x for a, y in enumerate(node_degrees) if n_svs[a] <= cutoff and n_svs[a] > cutoffs[i-1] for x in y]
    ax[i].hist(plot_list, bins=100, color='#80DDF9', alpha=0.9, label='Histogram')
    ax[i].set_xlabel('Node Degree')
    ax[i].set_ylabel('Cumulative Counts')
    ax[i].set_title('SV Clusters with <= {} SVs'.format(cutoff))
    ax[i].legend()
plt.show()
#
#
#
#
#
nd_graph = pd.DataFrame({
    'n_svs': n_svs,
    'node_degrees_avg': [np.mean(x) for x in node_degrees],
    'node_degrees_med': [np.median(x) for x in node_degrees],
    'node_degrees_std': [np.std(x) for x in node_degrees],
    'node_degrees_max': [max(x) for x in node_degrees],
    'node_degrees_min': [min(x) for x in node_degrees]
})

fig, ax = plt.subplots(1, 2, figsize=(12, 6))
# plot the average node degree vs the number of SVs, with error bars
ax[0].errorbar(
    nd_graph['n_svs'], nd_graph['node_degrees_avg'], 
    yerr=nd_graph['node_degrees_std'], fmt='o', color='#80DDF9'
)
ax[0].set_xlabel('Number of SVs')
ax[0].set_ylabel('Average Node Degree')
ax[0].set_title('Average Node Degree vs Number of SVs')
breaks = [5, 10, 100, 1000, 2074]
ax[0].set_xticks(breaks, breaks)
ax[0].set_xscale('log')
ax[0].set_yscale('log')
# plot the median node degree vs the number of SVs
ax[1].plot(
    nd_graph['n_svs'], nd_graph['node_degrees_med'], 'o', color='#80DDF9'
)
ax[1].set_xlabel('Number of SVs')
ax[1].set_ylabel('Median Node Degree')
ax[1].set_title('Median Node Degree vs Number of SVs')
ax[1].set_xticks(breaks, breaks)
ax[1].set_xscale('log')
ax[1].set_yscale('log')
plt.show()

#
#
#
#
#
#
#
