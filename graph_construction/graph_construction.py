# take in all bedpe in a directory, concatenate them
# normalize continuous features, one hot encode categorical features
# output a graph dataset for each cluster of every patient

# set up the paths =====================================================================================================
import os
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from tqdm import tqdm
from joblib import Parallel, delayed
import multiprocessing
import pickle
num_cores = multiprocessing.cpu_count()
# local vs UGER
if os.path.expanduser('~') in ["/Users/youyun", "/Users/youyunzheng"]: 
    # in a local mac, the home directory is usuaully at '/Users/[username]'
    workdir = os.path.expanduser('~')+"/Documents/HMS/PhD/beroukhimlab/broad_mount/"
else:
    # in dipg or uger, the home directory is usuaully at '/home/unix/[username]'
    workdir = "/xchip/beroukhimlab/"

bedpe_dir = workdir + "youyun/complexSV/data/TCGA/SV/bedpe/"
output_bedpe_dir = workdir + "youyun/complexSV/data/TCGA/SV/total_bedpe/"
output_networkx_dir = workdir + "youyun/complexSV/data/TCGA/SV/nx/"

# load data ==========================================================================================================
print("Loading data")
bedpe_files = [bedpe_dir +  f for f in os.listdir(bedpe_dir) if f.endswith('annotated.bedpe')]
# read in all the bedpe files and concatenate them into one variable
bedpe = pd.concat([pd.read_csv(f, sep='\t') for f in bedpe_files])
print('Finished loading data')
print('Number of SVs:', bedpe.shape[0])
bedpe[bedpe['cluster_size'] > 1].to_csv(output_bedpe_dir + 'cluster_bedpe_concat.csv', index=False)

# normalize and one hot encode ========================================================================================
print("Normalizing and one hot encoding")

# one hot encode chromosome, cnt type and sv class
scaled_bedpe = pd.get_dummies(bedpe[['seqnames_breakend1', 'seqnames_breakend2', 'cnt_type_breakend1', 'cnt_type_breakend2', 'svclass']])

# normalize continuous features

# start positions
scaled_bedpe['start_breakend1'] = bedpe['start_breakend1']/249250621
scaled_bedpe['start_breakend2'] = bedpe['start_breakend2']/249250621

# breakpoint characteristics
scaled_bedpe['ins_len'] = np.log10(bedpe['ins_len']+1)/np.log10(bedpe['ins_len']+1).max()
scaled_bedpe['mh_len'] = np.log10(bedpe['mh_len']+1)/np.log10(bedpe['mh_len']+1).max()
scaled_bedpe['ins_len'] = scaled_bedpe['ins_len'].fillna(0)
scaled_bedpe['mh_len'] = scaled_bedpe['mh_len'].fillna(0)
scaled_bedpe['N_ALT'] = np.log10(bedpe['N_ALT']+1)/np.log10(bedpe['N_ALT']+1).max()
scaled_bedpe['N_ALT_RP'] = np.log10(bedpe['N_ALT_RP']+1)/np.log10(bedpe['N_ALT_RP']+1).max()
scaled_bedpe['N_ALT_SR'] = np.log10(bedpe['N_ALT_SR']+1)/np.log10(bedpe['N_ALT_SR']+1).max()

# epigenetic features
scaled_bedpe['rep_time_breakend1'] = bedpe['rep_time_breakend1'].fillna(0)
scaled_bedpe['rep_time_breakend2'] = bedpe['rep_time_breakend2'].fillna(0)
scaled_bedpe['rep_time_breakend1'] = (scaled_bedpe['rep_time_breakend1']-scaled_bedpe['rep_time_breakend1'].min())/(scaled_bedpe['rep_time_breakend1'].max()-scaled_bedpe['rep_time_breakend1'].min())
scaled_bedpe['rep_time_breakend2'] = (scaled_bedpe['rep_time_breakend2']-scaled_bedpe['rep_time_breakend2'].min())/(scaled_bedpe['rep_time_breakend2'].max()-scaled_bedpe['rep_time_breakend2'].min())

# sub out '.' with 0 first
scaled_bedpe['gc_breakend1'] = bedpe['gc_breakend1'].replace('.', 0).fillna(0).astype(np.float64)/100
scaled_bedpe['gc_breakend2'] = bedpe['gc_breakend2'].replace('.', 0).fillna(0).astype(np.float64)/100

scaled_bedpe['gene_density_breakend1'] = np.log10(bedpe['gene_density_breakend1']+1)/np.log10(bedpe['gene_density_breakend1']+1).max()
scaled_bedpe['gene_density_breakend2'] = np.log10(bedpe['gene_density_breakend2']+1)/np.log10(bedpe['gene_density_breakend2']+1).max()

scaled_bedpe['dist_to_line_breakend1'] = bedpe['dist_to_line_breakend1']/bedpe['dist_to_line_breakend1'].max()
scaled_bedpe['dist_to_line_breakend2'] = bedpe['dist_to_line_breakend2']/bedpe['dist_to_line_breakend2'].max()
scaled_bedpe['dist_to_sine_breakend1'] = bedpe['dist_to_sine_breakend1']/bedpe['dist_to_sine_breakend1'].max()
scaled_bedpe['dist_to_sine_breakend2'] = bedpe['dist_to_sine_breakend2']/bedpe['dist_to_sine_breakend2'].max()
scaled_bedpe['dist_to_ltr_breakend1'] = bedpe['dist_to_ltr_breakend1']/bedpe['dist_to_ltr_breakend1'].max()
scaled_bedpe['dist_to_ltr_breakend2'] = bedpe['dist_to_ltr_breakend2']/bedpe['dist_to_ltr_breakend2'].max()
scaled_bedpe['dist_to_dna_breakend1'] = bedpe['dist_to_dna_breakend1']/bedpe['dist_to_dna_breakend1'].max()
scaled_bedpe['dist_to_dna_breakend2'] = bedpe['dist_to_dna_breakend2']/bedpe['dist_to_dna_breakend2'].max()
scaled_bedpe['dist_to_sr_breakend1'] = bedpe['dist_to_sr_breakend1']/bedpe['dist_to_sr_breakend1'].max()
scaled_bedpe['dist_to_sr_breakend2'] = bedpe['dist_to_sr_breakend2']/bedpe['dist_to_sr_breakend2'].max()
scaled_bedpe['dist_to_TAD_breakend1'] = bedpe['dist_to_TAD_breakend1']/bedpe['dist_to_TAD_breakend1'].max()
scaled_bedpe['dist_to_TAD_breakend2'] = bedpe['dist_to_TAD_breakend2']/bedpe['dist_to_TAD_breakend2'].max()
scaled_bedpe['dist_to_telomere_breakend1'] = bedpe['dist_to_telomere_breakend1']/bedpe['dist_to_telomere_breakend1'].max()
scaled_bedpe['dist_to_telomere_breakend2'] = bedpe['dist_to_telomere_breakend2']/bedpe['dist_to_telomere_breakend2'].max()
scaled_bedpe['dist_to_centromere_breakend1'] = bedpe['dist_to_centromere_breakend1']/bedpe['dist_to_centromere_breakend1'].max()
scaled_bedpe['dist_to_centromere_breakend2'] = bedpe['dist_to_centromere_breakend2']/bedpe['dist_to_centromere_breakend2'].max()

# copy number
scaled_bedpe['total_cn_breakend1'] = np.log2(bedpe['total_cn_breakend1']+1)/np.log2(bedpe['total_cn_breakend1']+1).max()
scaled_bedpe['major_cn_breakend1'] = np.log2(bedpe['major_cn_breakend1']+1)/np.log2(bedpe['major_cn_breakend1']+1).max()
scaled_bedpe['minor_cn_breakend1'] = np.log2(bedpe['minor_cn_breakend1']+1)/np.log2(bedpe['minor_cn_breakend1']+1).max()
scaled_bedpe['total_cn_breakend2'] = np.log2(bedpe['total_cn_breakend2']+1)/np.log2(bedpe['total_cn_breakend2']+1).max()
scaled_bedpe['major_cn_breakend2'] = np.log2(bedpe['major_cn_breakend2']+1)/np.log2(bedpe['major_cn_breakend2']+1).max()
scaled_bedpe['minor_cn_breakend2'] = np.log2(bedpe['minor_cn_breakend2']+1)/np.log2(bedpe['minor_cn_breakend2']+1).max()
scaled_bedpe['total_cn_breakend1'] = scaled_bedpe['total_cn_breakend1'].fillna(0)
scaled_bedpe['major_cn_breakend1'] = scaled_bedpe['major_cn_breakend1'].fillna(0)
scaled_bedpe['minor_cn_breakend1'] = scaled_bedpe['minor_cn_breakend1'].fillna(0)
scaled_bedpe['total_cn_breakend2'] = scaled_bedpe['total_cn_breakend2'].fillna(0)
scaled_bedpe['major_cn_breakend2'] = scaled_bedpe['major_cn_breakend2'].fillna(0)
scaled_bedpe['minor_cn_breakend2'] = scaled_bedpe['minor_cn_breakend2'].fillna(0)

scaled_bedpe['cluster_ID'] = bedpe['cluster_ID']
scaled_bedpe['cluster_size'] = bedpe['cluster_size']
scaled_bedpe['Sample'] = bedpe['Sample']

# Reviewing the features ===============================================================================================
print(scaled_bedpe.columns)
print('Starting Positions')
print(scaled_bedpe[['start_breakend1', 'start_breakend2']].describe())
print('Insertion Lengths, Microhomology Lengths, Number of ALTs, Number of ALTs from RP, Number of ALTs from SR')
print(scaled_bedpe[['ins_len', 'mh_len', 'N_ALT', 'N_ALT_RP', 'N_ALT_SR']].describe())
print('Replication Time')
print(scaled_bedpe[['rep_time_breakend1', 'rep_time_breakend2']].describe())
print('GC Content')
print(scaled_bedpe[['gc_breakend1', 'gc_breakend2']].describe())
print('Gene Density')
print(scaled_bedpe[['gene_density_breakend1', 'gene_density_breakend2']].describe())
print('Distances to Genomic Features')
print(scaled_bedpe[[
    'dist_to_line_breakend1', 'dist_to_line_breakend2', 'dist_to_sine_breakend1', 'dist_to_sine_breakend2', 
    'dist_to_dna_breakend1', 'dist_to_dna_breakend2', 'dist_to_ltr_breakend1', 'dist_to_ltr_breakend2'
]].describe())
print(scaled_bedpe[[
    'dist_to_sr_breakend1', 'dist_to_sr_breakend2', 'dist_to_TAD_breakend1', 'dist_to_TAD_breakend2', 
    'dist_to_telomere_breakend1', 'dist_to_telomere_breakend2', 'dist_to_centromere_breakend1', 'dist_to_centromere_breakend2'
]].describe())
print('Copy Number')
print(scaled_bedpe[['total_cn_breakend1', 'major_cn_breakend1', 'minor_cn_breakend1', 'total_cn_breakend2', 'major_cn_breakend2', 'minor_cn_breakend2']].describe())

# clustered SVs only ==================================================================================================
clustered_bedpe = scaled_bedpe[scaled_bedpe['cluster_size']>1]
print(clustered_bedpe[['Sample','cluster_ID','cluster_size']].drop_duplicates())
print('Number of clusters: ',clustered_bedpe[['Sample','cluster_ID','cluster_size']].drop_duplicates().shape[0])
print('Cluster size distribution: ', clustered_bedpe[['Sample','cluster_ID','cluster_size']].drop_duplicates()['cluster_size'].describe())

# save the clustered bedpe
clustered_bedpe.to_csv(output_bedpe_dir + 'cluster_bedpe_scaled.csv', index=False)

# # Graph construction ==================================================================================================
# print('Constructing graph')
# # parallelize the for loop
# # for each sample
# #     for each cluster
# #         do parallel for each row in the bedpe
# #             add the node to the graph
# #             for each other row in the bedpe
# #                 if the two nodes are connected
# #                     add the edge to the graph

# # construct graph
# G = nx.Graph()

# def construct_cluster_graph(sample, cluster, cluster_bedpe):
#     # print('Constructing graph for', sample, cluster)
#     # construct graph
#     G_cluster = nx.Graph()
#     # for each row in the bedpe
#     for i in range(cluster_bedpe.shape[0]):
#         # adding in all info aside from the Sample, cluster_ID, and cluster_size
#         nodal_features = cluster_bedpe.iloc[i].drop(['Sample', 'cluster_ID', 'cluster_size']).to_dict()
#         # add the node to the graph
#         G_cluster.add_node(
#             cluster_bedpe.iloc[i]['Sample'] + '_' + str(cluster_bedpe.iloc[i]['cluster_ID']) + '_' + str(i), 
#             **nodal_features
#         )
#         # for each other row in the bedpe
#         for j in range(i+1, cluster_bedpe.shape[0]):
#             # checking if either of the breakends of node A is downstream of either of the breakends of node B
#             # if the two nodes are connected
#             if (
#                 cluster_bedpe.iloc[i]['cnt_type_breakend1_+'] and 
#                 cluster_bedpe.iloc[j]['cnt_type_breakend2_-'] and 
#                 cluster_bedpe.iloc[i]['start_breakend1']<cluster_bedpe.iloc[j]['start_breakend2']
#             ) or (
#                 cluster_bedpe.iloc[i]['cnt_type_breakend1_-'] and 
#                 cluster_bedpe.iloc[j]['cnt_type_breakend2_+'] and 
#                 cluster_bedpe.iloc[i]['start_breakend1']>cluster_bedpe.iloc[j]['start_breakend2']
#             ) or (
#                 cluster_bedpe.iloc[i]['cnt_type_breakend2_+'] and 
#                 cluster_bedpe.iloc[j]['cnt_type_breakend1_-'] and 
#                 cluster_bedpe.iloc[i]['start_breakend2']<cluster_bedpe.iloc[j]['start_breakend1']
#             ) or (
#                 cluster_bedpe.iloc[i]['cnt_type_breakend2_-'] and 
#                 cluster_bedpe.iloc[j]['cnt_type_breakend1_+'] and 
#                 cluster_bedpe.iloc[i]['start_breakend2']>cluster_bedpe.iloc[j]['start_breakend1']
#             ) or (
#                 cluster_bedpe.iloc[i]['cnt_type_breakend1_+'] and 
#                 cluster_bedpe.iloc[j]['cnt_type_breakend1_-'] and 
#                 cluster_bedpe.iloc[i]['start_breakend1']<cluster_bedpe.iloc[j]['start_breakend1']
#             ) or (
#                 cluster_bedpe.iloc[i]['cnt_type_breakend1_-'] and 
#                 cluster_bedpe.iloc[j]['cnt_type_breakend1_+'] and
#                 cluster_bedpe.iloc[i]['start_breakend1']>cluster_bedpe.iloc[j]['start_breakend1']
#             ) or (
#                 cluster_bedpe.iloc[i]['cnt_type_breakend2_+'] and 
#                 cluster_bedpe.iloc[j]['cnt_type_breakend2_-'] and 
#                 cluster_bedpe.iloc[i]['start_breakend2']<cluster_bedpe.iloc[j]['start_breakend2']
#             ) or (
#                 cluster_bedpe.iloc[i]['cnt_type_breakend2_-'] and 
#                 cluster_bedpe.iloc[j]['cnt_type_breakend2_+'] and 
#                 cluster_bedpe.iloc[i]['start_breakend2']>cluster_bedpe.iloc[j]['start_breakend2']
#             ):
#                 G_cluster.add_edge(
#                     cluster_bedpe.iloc[i]['Sample'] + '_' + str(cluster_bedpe.iloc[i]['cluster_ID']) + '_' + str(i),
#                     cluster_bedpe.iloc[j]['Sample'] + '_' + str(cluster_bedpe.iloc[j]['cluster_ID']) + '_' + str(j)
#                 )
#     # save cluster graph
#     pickle.dump(G_cluster, open(output_networkx_dir + sample + '_' + str(cluster) + "_nx_graph.txt", 'wb'))
#     return G_cluster

# for sample in tqdm(clustered_bedpe['Sample'].unique()):
#     sample_bedpe = clustered_bedpe[clustered_bedpe['Sample']==sample]
#     G_sample = nx.Graph()
#     G_sample = Parallel(n_jobs=num_cores)(
#         delayed(construct_cluster_graph)(
#             sample, cluster, sample_bedpe[sample_bedpe['cluster_ID']==cluster]
#         ) for cluster in sample_bedpe['cluster_ID'].unique()
#     )
#     # G_sample = nx.compose_all(G_sample)
#     # G = nx.compose(G, G_sample)
    
# print('Finished constructing graph')
# print('Graphs are saved in ', output_networkx_dir)


# # # save the graph =====================================================================================================
# # print('Saving the graph')
# # import pickle
# # pickle.dump(G, open(output_networkx_dir + "total_nx_graph.txt", 'wb'))
# # print('Done!')
