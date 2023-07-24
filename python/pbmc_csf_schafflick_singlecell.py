#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import numpy as np
import pandas as pd
import os


# In[2]:


os.chdir('/home/shao11/data_kfitzg13/msgwas_shao11')


# In[3]:


schafflick = sc.read_h5ad('data/sc_data/pbmc/schafflick/processed/MS_CSF.h5ad')
# copy of data exists in 'data/sc_data/csf/MS_CSF.h5ad'


# In[4]:


schafflick.obs


# # Filter to PBMC/CSF only

# In[5]:


# remove doublets and contamination
filt = schafflick[~((schafflick.obs.labels=='Mono Doublet')|(schafflick.obs.labels=='B cell doublets')|(schafflick.obs.labels=='doublet')|(schafflick.obs.labels=='contamination1'))]
# filter out pbmc
pbmc = filt[filt.obs.CSF == 'False']
# labels includes RBC
# repeat for csf
csf = filt[filt.obs.CSF =='True']
# no CD8n in CSF


# In[6]:


pbmc.shape


# In[7]:


csf.shape


# # MAGMA input

# ## PBMC

# In[8]:


gene_coordinates = pd.read_csv('data/NCBI37.3.gene.loc', sep = '\t', names = ['ENTREZ', 'chr', 'gene_start', 'gene_end', 'strand', 'symbol'])


# In[9]:


# filter out extended MHC locus1
gene_coordinates = gene_coordinates[~((gene_coordinates['chr']!=6) & (gene_coordinates['gene_start']>25000000) & (gene_coordinates['gene_end']<34000000))]


# In[10]:


gene_coordinates['start'] = gene_coordinates['gene_start'].apply(lambda x: 0 if (x-100000)<0 else (x-100000))
gene_coordinates['end'] = gene_coordinates['gene_end'].apply(lambda x: x+100000)


# In[13]:


# convert to dataframe for pandas operations
pbmc_counts = pbmc.to_df()
csf_counts = csf.to_df()


# ## Run this for PBMC

# In[15]:


aggr_exp = pbmc_counts.merge(pbmc.obs['labels'], right_index=True, left_index=True)
cell_pop = "pbmc"


# ## Run this for CSF

# In[16]:


aggr_exp = csf_counts.merge(csf.obs['labels'], right_index=True, left_index=True)
cell_pop = "csf"


# ## get total raw counts for each gene for each cell type

# In[17]:


aggr_exp = aggr_exp.groupby('labels').sum()


# In[18]:


# calculate TPM for each gene in each cell type
TPM = aggr_exp.apply(lambda i: (i*(10**6))/i.sum(), axis=1)


# In[19]:


# calculate specificity = proportion of gene expressed by the cell type
# drop any genes that are not expressed
TPM = TPM.loc[:,(TPM!=0).any(axis=0)]
specificity = TPM.apply(lambda j: j/j.sum(), axis=0)


# In[20]:


# filter for MAGMA genes
magma_genes = gene_coordinates[gene_coordinates['symbol'].isin(specificity.columns)]
specificity = specificity.loc[:,magma_genes['symbol']]


# In[22]:


# write specificity matrix to file
specificity.to_csv(f'data/processed/{cell_pop}_MAGMA_specificity.csv', index=True, header=True)


# In[23]:


# convert gene symbols to ENTREZ IDs
rename_dict = dict(zip(magma_genes['symbol'], magma_genes['ENTREZ']))
ENTREZ_specificity = specificity.rename(columns=rename_dict)


# In[24]:


# get top 10 specific genes
n_genes = round(0.1*ENTREZ_specificity.shape[1])
top10 = ENTREZ_specificity.apply(lambda i: ENTREZ_specificity.sort_values(by=i.name, axis=1).columns[:n_genes], axis=1)


# In[25]:


df = pd.DataFrame(dict(zip(top10.index, top10)))


# In[26]:


df.transpose().to_csv(f'data/processed/{cell_pop}_MAGMA_genesets.bed', index=True, header=False, sep='\t')


# In[27]:


# symbol version
symb_top10 = specificity.apply(lambda i: specificity.sort_values(by=i.name, axis=1).columns[:n_genes], axis=1)
df = pd.DataFrame(dict(zip(symb_top10.index, symb_top10)))
df.transpose().to_csv(f'data/processed/{cell_pop}_MAGMA_genesets_symbols.tsv', index=True, header=True, sep='\t')


# In[33]:


# sanity check: FOXP3 should be relatively specific for Tregs
specificity['FOXP3']

