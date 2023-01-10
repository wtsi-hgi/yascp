#!/usr/bin/env python3

__date__ = '2020-03-13'
__version__ = '0.0.1'

import sys
import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
from distutils.version import LooseVersion
import scipy
import scipy.io
import gzip
import pandas
import scanpy
import anndata
## split h5ad file by chromium channel

# for testing use
# infnam = "/lustre/scratch123/hgi/mdt1/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/franke_data/work/2a/ebdc5079ae949777263ce3b1aca510/BF61CE54F4603C9F-adata.h5ad"

def write_h5_out_for_ct(ad,oufn_list_AZ,oufnam,oufn_list,samples,samples_AZ,bl,count,anndata_compression_opts,batch_labels=None):
    samples[count]={'experiment_id':bl,'h5ad_filepath':f"{os.getcwd()}/{oufnam}"}
    samples_AZ[count]={'experiment_id':bl,'h5ad_filepath':f"{os.getcwd()}/AZ_{oufnam}"}
    if (bl=='full'):
        adb =ad
    else:
        adb = ad[batch_labels == bl,:]
        
    # strip unneccessary meta data - this seems to aid subsequent transformation to Seurat h5 format
    # assumes that cells failing qc already were stripped out
    # ad.obs = ad.obs[['convoluted_samplename', 'cell_passes_qc']]
    # ad.var = ad.var[['feature_types', 'genome', 'gene_symbols']]
    adb_AZ = adb.copy()
    # disable
    adb_AZ.obs = pandas.DataFrame(adb_AZ.obs.index, index = adb_AZ.obs.index, columns = ["cell_barcode"])

    adb.layers['counts'] = adb.X.copy()

    # Total-count normalize (library-size correct) the data matrix X to
    # counts per million, so that counts become comparable among cells.
    scanpy.pp.normalize_total(
        adb,
        target_sum=1e4,
        exclude_highly_expressed=False,
        key_added='normalization_factor',  # add to adata.obs
        inplace=True
    )
    # Logarithmize the data: X = log(X + 1) where log = natural logorithm.
    # Numpy has a nice function to undo this np.expm1(adata.X).
    scanpy.pp.log1p(adb)
    # Delete automatically added uns - UPDATE: bad idea to delete as this slot
    # is used in _highly_variable_genes_single_batch.
    # del adata.uns['log1p']
    # Add record of this operation.
    # adata.layers['log1p_cpm'] = adata.X.copy()
    # adata.uns['log1p_cpm'] = {'transformation': 'ln(CPM+1)'}
    adb.layers['log1p_cp10k'] = adb.X.copy()
    adb.uns['log1p_cp10k'] = {'transformation': 'ln(CP10k+1)'}

    # Reset X to counts
    adb.X = adb.layers['counts'].copy()  

    try:
        vdf = adb_AZ.var[["feature_types", "genome"]]
    except:
        adb_AZ.var['genome']='GRCh38'
        vdf = adb_AZ.var[["feature_types", "genome"]]
    vdf.insert(1,"gene_ids", vdf.index)
    vdf.index = pandas.Index(adb_AZ.var['gene_symbols'].astype('str'))
    #ad.var = vdf.set_index("gene_symbols", drop = True, verify_integrity = False)
    adb_AZ.var = vdf

    del adb_AZ.uns
    adata = scanpy.AnnData(adb_AZ.X)
    adata.obs= adb_AZ.obs
    adata.var = adb_AZ.var
    adb_AZ = adata
    # disable
    if anndata_compression_opts is None:
        adb.write(oufnam)
        adb_AZ.write(f"AZ_{oufnam}")
    else:
        adb.write(
            oufnam,
            compression='gzip',
            compression_opts=anndata_compression_opts
        )
        adb_AZ.write(
            f"AZ_{oufnam}",
            compression='gzip',
            compression_opts=anndata_compression_opts
        )                
        
    oufn_list_AZ.append(f'AZ_{oufnam}')
    oufn_list.append(oufnam)
    return oufn_list_AZ,oufn_list,samples,samples_AZ

def split_h5ad_by_batch(ad, oufnprfx, colnam_batch = 'batch', anndata_compression_opts=None, option='true'):
    oufn_list_fnam = '{}_files.txt'.format(oufnprfx)
    oufn_list = []
    oufn_list_AZ = []
    batch_labels = pandas.Categorical(ad.obs[colnam_batch].apply(lambda a: a.split('__')[0])) # <class 'pandas.core.series.Series'>
    samples = {}
    samples_AZ = {}
    count=0
    if (option == 'true'):
        print("Here we split the h5ad we just normalise for celltype assignmet")
        for bl in batch_labels.categories:
            print(bl)
            count+=1
            oufnam = '{0}_{1}.h5ad'.format(oufnprfx, bl)
            oufn_list_AZ,oufn_list,samples,samples_AZ = write_h5_out_for_ct(ad,oufn_list_AZ,oufnam,oufn_list,samples,samples_AZ,bl,count,anndata_compression_opts,batch_labels=batch_labels)

    else:
        print("Here we dont split the h5ad we just normalise for celltype assignmet")
        bl ='full'
        count+=1
        oufnam = '{0}_{1}.h5ad'.format(oufnprfx, bl)
        oufn_list_AZ,oufn_list,samples,samples_AZ = write_h5_out_for_ct(ad,oufn_list_AZ,oufnam,oufn_list,samples,samples_AZ,bl,count,anndata_compression_opts)

    Samples = pandas.DataFrame(samples).T
    Samples.to_csv('Samples.tsv',sep='\t',index=False)
    
    Samples_AZ = pandas.DataFrame(samples_AZ).T
    Samples_AZ.to_csv('AZ_Samples.tsv',sep='\t',index=False)
    with open(oufn_list_fnam, 'w') as oufh:
        for fn in oufn_list:
            oufh.write(fn + '\n')
    return oufn_list_fnam

if __name__ == '__main__':
    nargs = len(sys.argv)
    # if nargs != 3:
    #     sys.exit("usage: %s ./<input_file *.h5ad> <output_file_prefix>".format(sys.argv[0]))
    infnam = sys.argv[1]
    oufnprfx = sys.argv[2]
    option = sys.argv[3]

    # print(infnam)
    ad = scanpy.read(infnam)
    split_h5ad_by_batch(ad, oufnprfx, colnam_batch = 'batch', anndata_compression_opts = None,option=option)
    sys.exit()
