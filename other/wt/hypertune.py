#!/usr/bin/env python
# coding: utf-8

# In[9]:


import scanpy as sc
import torch
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt
import numpy as np
import argparse
import anndata


# In[3]:


adata_whole=sc.read("../../imporatant_process_data/4.12_wt_integrate_SCRAN_log1p.h5ad")


# In[6]:


parser = argparse.ArgumentParser()
parser.add_argument('--hvg_n' , required=True) # int 
parser.add_argument('--nlatent' , required=True) 

args = parser.parse_args()
print(args)
hvg_n=int(args.hvg_n)
nlatent=int(args.nlatent)

print('hvg:',hvg_n,"nlatent:",nlatent)

unlabeled = adata_whole.obs.index[[ coarse_label != "Osteoblasts" for coarse_label in adata_whole.obs.coarse_label]]
adata_whole.obs['scanvi_label'] = adata_whole.obs.coarse_label.tolist()
adata_whole.obs.loc[unlabeled,'scanvi_label'] = "unlabeled"

condition_key = 'batch'
cell_type_key = 'scanvi_label'
unlabeled_category = "unlabeled"


# In[7]:



sc.pp.highly_variable_genes(adata_whole, batch_key="batch",min_mean=0.035, flavor="cell_ranger",n_top_genes=hvg_n)
adata = adata_whole[:,adata_whole.var.highly_variable].copy()
adata.X = adata.layers['counts']
adata.raw = adata
raw = adata.raw.to_adata()
raw.X = adata.layers['counts']
adata.raw = raw
sca.models.SCVI.setup_anndata(adata, batch_key=condition_key, labels_key=cell_type_key)
vae = sca.models.SCVI(
    adata,
    n_latent=nlatent,
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
    gene_likelihood="zinb",
    n_layers=1,
    dropout_rate=0.1
    )
vae.train(max_epochs=50)
scanvae = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category =unlabeled_category)
scanvae.train(max_epochs=10)
reference_latent = sc.AnnData(scanvae.get_latent_representation())
reference_latent.obs["cell_type"] = adata.obs[cell_type_key].tolist()
reference_latent.obs["batch"] = adata.obs[condition_key].tolist()
reference_latent.obs["Project"] = adata.obs["Project"].tolist()
reference_latent.obs["coarse_label"] = adata.obs["coarse_label"].tolist()
sc.pp.neighbors(reference_latent, n_neighbors=8)
sc.tl.leiden(reference_latent)
sc.tl.umap(reference_latent)
sc.pl.umap(reference_latent,
    color=['batch', 'coarse_label'],
    frameon=False,
    wspace=0.6,save="gene_{}_latent_{}_umap".format(hvg_n, nlatent)
    )
reference_latent.write("../../unimportant_processed_data/hypertune/gene_{}_latent_{}_umap.h5ad".format(hvg_n, nlatent))

