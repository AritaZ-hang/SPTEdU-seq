import scvelo as scv
import pandas as pd
import numpy as np
import os

loom_data = scv.read("./Hippocampus_Top5W_degen_V3CE0.loom", cache=False)

meta_path="./Velocyto_data/"
sample_obs=pd.read_csv(os.path.join(meta_path, "cellID_obs.csv"))
barcode_coords=pd.read_csv(os.path.join(meta_path, "spatial_coords.csv"), header=0, names=["CellID", "x_coords", "y_coords"])
barcode_clusters=pd.read_csv(os.path.join(meta_path, "cell_clusters.csv"), header=0, names=["CellID", "cluster"])
barcode_umap=pd.read_csv(os.path.join(meta_path, "cell_embeddings.csv"), header=0, names=["CellID", "UMAP_1", "UMAP_2"])

hippocampus=loom_data[np.isin(loom_data.obs.index, sample_obs)]
hippocampus_index=pd.DataFrame(hippocampus.obs.index)
coords_ordered=pd.merge(hippocampus_index, barcode_coords, on = ["CellID"])
clusters_ordered=pd.merge(hippocampus_index, barcode_clusters, on = ["CellID"])
umap_ordered=pd.merge(hippocampus_index, barcode_umap, on = ["CellID"])

coords_ordered=coords_ordered.iloc[:,1:]
clusters_ordered=clusters_ordered.iloc[:,1:]
umap_ordered=umap_ordered.iloc[:,1:]

hippocampus.obsm["X_coords"] = coords_ordered.values
hippocampus.obsm["X_umap"] = umap_ordered.values
hippocampus.uns["clusters"] = clusters_ordered.values
hippocampus.obs["clusters"] = clusters_ordered.values

adata = hippocampus
adata.var_names_make_unique()
# save model to file
adata.write("./Velocyto_data/Hippocampus_Velocyto.h5ad", compression='gzip')

# for cortex

cortex_barcodes=pd.read_csv(os.path.join("./", "cortex_barcodes_obs.csv"), header=0, names=["CellID"])
adata_cortex = adata[np.isin(adata.obs.index, cortex_barcodes["CellID"])]
os.mkdir("./Cortex")
adata_cortex.write("./Cortex/cortex_velocyto.h5ad")

# scVelo
scv.pp.filter_and_normalize(adata_cortex, min_shared_counts =5, n_top_genes=2000)
scv.pp.moments(adata_cortex, n_pcs=25, n_neighbors=30)
scv.tl.recover_dynamics(adata_cortex, n_jobs = 12)
scv.tl.velocity(adata_cortex, mode='stochastic')
scv.tl.velocity_graph(adata_cortex, basis = "X_umap", n_jobs = 12)

colors=["#246CA5", "#5395CC", "#C2D8EE"]
scv.pl.velocity_embedding_grid(adata_cortex, basis="umap", color="clusters", scale=0.6, arrow_size=1, dpi = 300, title = "", arrow_length=2, size = 3, alpha = 1, figsize=(4, 4), 
                               save="scvelo_velocity_embedding_grid_umap_cortex_only.pdf", 
                               palette = colors)

scv.pl.velocity_embedding_grid(adata_cortex, basis="coords", color="clusters", scale=0.6, arrow_size=1, dpi = 300, title = "", arrow_length=2, size = 3, alpha = 1, figsize=(4, 4), 
                               save="scvelo_velocity_embedding_grid_coords_cortex_only.pdf", 
                               palette = colors)

adata_cortex.write("./Cortex/Cortex_processed.h5ad")