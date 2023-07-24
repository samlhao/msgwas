#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd


# In[6]:


cross_sect = pd.read_csv('../../../data/oct/cross_sectional_results.csv')
longitudinal = pd.read_csv('../../../data/oct/longitudinal_results.csv')


# In[4]:


cross_sect


# In[10]:


cross_sect = cross_sect.rename({'rsid':'SNP', 'chr':'CHR', 'pos':'BP'}, axis=1)
longitudinal = longitudinal.rename({'rsid':'SNP', 'chr':'CHR', 'pos':'BP'}, axis=1)


# In[12]:


cross_sect[['SNP', 'CHR', 'BP']].to_csv('../../../data/oct/cross_sectional_snploc.tsv', sep='\t', index=False)
longitudinal[['SNP', 'CHR', 'BP']].to_csv('../../../data/oct/longitudinal_snploc.tsv', sep='\t', index=False)


# In[13]:


cross_sect_pval = cross_sect[['SNP', 'pval.HIPOD1']]
cross_sect_pval = cross_sect_pval.rename({'pval.HIPOD1':'P'}, axis=1)
cross_sect_pval.to_csv('../../../data/oct/cross_sectional_composite_pval.tsv', sep='\t', index=False)
longitudinal_pval = longitudinal[['SNP', 'pval.HIPOD1']]
longitudinal_pval = longitudinal_pval.rename({'pval.HIPOD1':'P'}, axis=1)
longitudinal_pval.to_csv('../../../data/oct/longitudinal_composite_pval.tsv', sep='\t', index=False)


# In[ ]:




