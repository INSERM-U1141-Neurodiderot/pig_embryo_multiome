import sys
import os
import scanpy as sc
from pycisTopic.cistopic_class import *
import pickle
from pycisTopic.lda_models import *
from pycisTopic.clust_vis import *
from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates

from scenicplus.wrappers.run_pycistarget import run_pycistarget
from scenicplus.TF_to_gene import *
from scenicplus.grn_builder.gsea_approach import build_grn
import scenicplus

from scenicplus.scenicplus_class import SCENICPLUS, create_SCENICPLUS_object
from scenicplus.preprocessing.filtering import *

import dill
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import pyranges as py
import numpy as np

from scenicplus.cistromes import *
import time
from scenicplus.enhancer_to_gene import get_search_space, calculate_regions_to_genes_relationships, GBM_KWARGS

tmp_dir = "/tmp"
work_dir = '/home/adufour/work/scenic_omics/all_v4'
adata = sc.read_h5ad(os.path.join(work_dir, 'scRNA/adata.h5ad'))
scRNA_bc = adata.obs_names
cell_data = adata.obs
cell_data['sample_id'] = cell_data['Sample']
cell_data['celltype'] = cell_data['Clusters'].astype(str) # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.

fragments_dict = {'J7': "/home/adufour/work/mat_multiome/atac_fragments_j7.tsv.gz",
                  'J9-1': "/home/adufour/work/mat_multiome/atac_fragments_j9-1.tsv.gz",
                  'J9-2': "/home/adufour/work/mat_multiome/atac_fragments_j9-2.tsv.gz",
                  'J11-1': "/home/adufour/work/mat_multiome/atac_fragments_j11-1.tsv.gz",
                  'J11-2': "/home/adufour/work/mat_multiome/atac_fragments_j11-2.tsv.gz",
                  'J11-3': "/home/adufour/work/mat_multiome/atac_fragments_j11-3.tsv.gz",
                  'J11-4': "/home/adufour/work/mat_multiome/atac_fragments_j11-4.tsv.gz",
                  'lw7': "/home/adufour/work/mat_multiome/atac_fragments_lw7.tsv.gz",
                  'lw9': "/home/adufour/work/mat_multiome/atac_fragments_lw9.tsv.gz",
                 }
path_to_regions = {'J7':os.path.join(work_dir, 'scATAC/consensus_peak_calling/consensus_regions.bed'),
                  'J9-1':os.path.join(work_dir, 'scATAC/consensus_peak_calling/consensus_regions.bed'),
                  'J9-2':os.path.join(work_dir, 'scATAC/consensus_peak_calling/consensus_regions.bed'),
                  'J11-1':os.path.join(work_dir, 'scATAC/consensus_peak_calling/consensus_regions.bed'),
                  'J11-2':os.path.join(work_dir, 'scATAC/consensus_peak_calling/consensus_regions.bed'),
                  'J11-3':os.path.join(work_dir, 'scATAC/consensus_peak_calling/consensus_regions.bed'),
                  'lw7':os.path.join(work_dir, 'scATAC/consensus_peak_calling/consensus_regions.bed'),
                  'lw9':os.path.join(work_dir, 'scATAC/consensus_peak_calling/consensus_regions.bed'),
                  'J11-4':os.path.join(work_dir, 'scATAC/consensus_peak_calling/consensus_regions.bed')}
metadata_bc = pickle.load(open(os.path.join(work_dir, 'scATAC/quality_control/metadata_bc.pkl'), 'rb'))

bc_passing_filters = {"J7": metadata_bc['J7'].index.values.categories.tolist(),
                     "J9-1": metadata_bc['J9-1'].index.values.categories.tolist(),
                     "J9-2": metadata_bc['J9-2'].index.values.categories.tolist(),
                     "J11-1": metadata_bc['J11-1'].index.values.categories.tolist(),
                     "J11-2": metadata_bc['J11-2'].index.values.categories.tolist(),
                     "J11-3": metadata_bc['J11-3'].index.values.categories.tolist(),
                     "lw7": metadata_bc['lw7'].index.values.categories.tolist(),
                     "lw9": metadata_bc['lw9'].index.values.categories.tolist(),
                     "J11-4": metadata_bc['J11-4'].index.values.categories.tolist()}

#print(f"{len(list(set(bc_passing_filters['merge']) & set(scRNA_bc)))} cell barcodes pass both scATAC-seq and scRNA-seq based filtering")

cistopic_obj_list = [create_cistopic_object_from_fragments(
                            path_to_fragments=fragments_dict[key],
                            path_to_regions=path_to_regions[key],
                            metrics=metadata_bc[key],
                            valid_bc=list(set(bc_passing_filters[key]) & set(scRNA_bc)),
                            n_cpu=40,
                            project=key,
                            split_pattern='-') for key in fragments_dict.keys()]

cistopic_obj = merge(cistopic_obj_list)

cell_data['barcode'] = [x for x in cell_data.index.tolist()]

cell_data.index = cell_data['barcode'] + '-' + cell_data['sample_id'].astype(str) + '___' + cell_data['sample_id'].astype(str)

pickle.dump(cistopic_obj, open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))

cistopic_obj.add_cell_data(cell_data, split_pattern='-')
print(cistopic_obj)

models = run_cgs_models(cistopic_obj,
                        n_topics=[16,32],
                        n_cpu=40,
                        n_iter=500,
                        random_state=555,
                        alpha=50,
                        alpha_by_topic=True,
                        eta=0.1,
                        eta_by_topic=False,
                        save_path=None,
                        _temp_dir=tmp_dir)

if not os.path.exists(os.path.join(work_dir, 'scATAC/models')):
    os.makedirs(os.path.join(work_dir, 'scATAC/models'))

pickle.dump(models,
            open(os.path.join(work_dir, 'scATAC/models/10x_pbmc_models_500_iter_LDA.pkl'), 'wb'))

model = evaluate_models(models,
                       select_model=16,
                       return_model=True,
                       metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                       plot_metrics=False)

cistopic_obj.add_LDA_model(model)
pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))

#run_umap(cistopic_obj, target  = 'cell', scale=True)

region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)

imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable='celltype', var_features=variable_regions, split_pattern = '-')

if not os.path.exists(os.path.join(work_dir, 'scATAC/candidate_enhancers')):
    os.makedirs(os.path.join(work_dir, 'scATAC/candidate_enhancers'))
pickle.dump(region_bin_topics_otsu, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
pickle.dump(region_bin_topics_top3k, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))
pickle.dump(markers_dict, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict.pkl'), 'wb'))

####
#
# PYCISTARGET
#
####

region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top_3'] = {}
region_sets['DARs'] = {}
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for topic in region_bin_topics_top3k.keys():
    regions = region_bin_topics_top3k[topic].index
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index
    region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))

for key in region_sets.keys():
    print(f'{key}: {region_sets[key].keys()}')

rankings_db = "/home/adufour/work/cistargetdb/feather_v3_updown/Sus_scrofa.genes_vs_motifs.rankings.feather"
scores_db =  "/home/adufour/work/cistargetdb/feather_v3_updown/Sus_scrofa.genes_vs_motifs.scores.feather"
motif_annotation = '/home/adufour/work/cistargetdb/motif2tf_v3.tbl'

if not os.path.exists(os.path.join(work_dir, 'motifs')):
    os.makedirs(os.path.join(work_dir, 'motifs'))

tss = pd.read_csv("/home/adufour/work/scenic_omics/tss.csv")

run_pycistarget(
    region_sets = region_sets,
    species = 'custom',
    custom_annot = tss,
    save_path = os.path.join(work_dir, 'motifs'),
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    n_cpu = 40,
    _temp_dir = tmp_dir
    )

menr = dill.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))

####
#
# SCENICPLUS
#
####

_stderr = sys.stderr
null = open(os.devnull,'wb')

if not os.path.exists(os.path.join(work_dir, 'SCENIC')):
    os.makedirs(os.path.join(work_dir, 'SCENIC'))

adata = sc.read_h5ad(os.path.join(work_dir, 'scRNA/adata.h5ad'))
cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))

adata.obs_names = adata.obs_names + '-' + adata.obs['Sample'].astype(str) + '___' + adata.obs['Sample'].astype(str)

scplus_obj = create_SCENICPLUS_object(
    GEX_anndata = adata.raw.to_adata(),
    cisTopic_obj = cistopic_obj,
    menr = menr
    )
scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())
scplus_obj

start_time = time.time()
merge_cistromes(scplus_obj)
time = time.time()-start_time
print(time/60)

tss_df = pd.read_csv("/home/adufour/work/scenic_omics/tss_scp.csv")
chromsize_df = pd.read_csv("/home/adufour/work/scenic_omics/chromsize.csv")

tss = pr.PyRanges(tss_df)
chromsize = pr.PyRanges(chromsize_df)

get_search_space(scplus_obj,
                 pr_annot = tss,
                 pr_chromsizes = chromsize,
                 upstream = [1000, 150000],
                 downstream = [1000, 150000])

calculate_regions_to_genes_relationships(scplus_obj,
                    ray_n_cpu = 40,
                    _temp_dir = tmp_dir,
                    importance_scoring_method = 'GBM',
                    importance_scoring_kwargs = GBM_KWARGS)

tf_file = '/home/adufour/work/cistargetdb/tf_list.txt'

calculate_TFs_to_genes_relationships(scplus_obj,
                    tf_file = tf_file,
                    ray_n_cpu = 40,
                    method = 'GBM',
                    _temp_dir = tmp_dir,
                    key= 'TF2G_adj')

build_grn(scplus_obj,
         min_target_genes = 10,
         adj_pval_thr = 1,
         min_regions_per_gene = 0,
         quantiles = (0.85, 0.90, 0.95),
         top_n_regionTogenes_per_gene = (5, 10, 15),
         top_n_regionTogenes_per_region = (),
         binarize_using_basc = True,
         rho_dichotomize_tf2g = True,
         rho_dichotomize_r2g = True,
         rho_dichotomize_eregulon = True,
         rho_threshold = 0.05,
         keep_extended_motif_annot = True,
         merge_eRegulons = True,
         order_regions_to_genes_by = 'importance',
         order_TFs_to_genes_by = 'importance',
         key_added = 'eRegulons_importance',
         cistromes_key = 'Unfiltered',
         disable_tqdm = False, #If running in notebook, set to True
         ray_n_cpu = 40,
         _temp_dir = tmp_dir)

with open(os.path.join(work_dir, 'SCENIC/scplus_obj.pkl'), 'wb') as f:
  dill.dump(scplus_obj, f)
