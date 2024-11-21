import sys
import numpy as np
import cellrank as cr
import scanpy as sc
import scvelo as scv
from colormap import Colormap
import matplotlib.pyplot as pl
from matplotlib.pyplot import rc_context
c = Colormap()
rnacmap = c.cmap_linear("#2971b1", "#F7F7F7","#bb2933")

scv.settings.verbosity = 3
cr.settings.verbosity = 2

adata = scv.read("KB_python_adata.h5ad")
scv.pp.filter_and_normalize(adata,
                            min_counts = 10, 
                            n_top_genes=1000, enforce = False, subset_highly_variable = False)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs = 30, n_neighbors = 30, random_state=0)
scv.pp.moments(adata, n_pcs = None, n_neighbors=None)

scv.tl.recover_dynamics(adata, n_jobs = 8)
scv.tl.velocity(adata, mode = "dynamical")

vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()

merge_colors=["#009392", "#39B185", "#E9E29C", "#E88471", "#CF597E"]
vk.plot_projection(basis = "X_coords", size=15, color="merge_clusters", palette = merge_colors, alpha = 1, title = "", legend_loc = "none", stream=False, dpi=300, save="CellRank_1st_velocity.pdf", figsize=(4,4))

# select root cells
sc.tl.diffmap(adata)
root_ixs = adata.obsm["X_diffmap"][:,2].argmax()
scv.pl.scatter(
adata, 
basis = ["X_diffmap"], color = ["merge_clusters", root_ixs], 
    legend_loc = "right", 
    title = "", 
    dpi = 300, 
    save = "CellRank_Diffmap.pdf", 
    figsize=(8, 4)
)
adata.uns["iroot"] = root_ixs

sc.tl.dpt(adata)
scv.pl.scatter(
adata, 
basis = ["X_coords"], color = ["dpt_pseudotime"],  
    title = "", 
    dpi = 300, 
    figsize=(4, 4), 
    size=10, 
    color_map = "RdBu_r", 
    save="CellRank_DPT_Pseudotime.pdf"
)

traj = [3, 1, 5]
mask = np.in1d(adata.obs["merge_clusters"], traj)

with rc_context({"figure.figsize":(4,4)}):
    sc.pl.violin(adata[mask], 
            keys = ["dpt_pseudotime"], 
            groupby = "merge_clusters", 
            rotation = -90, 
            order=traj, 
            inner = "box", jitter=False, 
            stripplot = False, 
            alpha=1, 
            show=False)
pl.savefig("./figures/CellRank_DPT_Pseudotime.pdf")

pk = cr.kernels.PseudotimeKernel(adata, time_key = "dpt_pseudotime")
pk.compute_transition_matrix()
print(pk)

pk.plot_projection(basis = "X_coords", 
                   size=15, color="merge_clusters", 
                   palette = merge_colors, alpha = 1, title = "", legend_loc = "none", stream=False, recompute=True, 
                    dpi = 300, save = "CellRank_2nd_Pseudotime.pdf")

# genes trend
g = cr.estimators.GPCCA(pk)
print(g)
g.fit(cluster_key = "merge_clusters", n_states=2)
g.plot_macrostates(which = "all", discrete = True, legend_loc = "right", s= 100, basis = "X_merge_umap")
g.predict_initial_states()
g.plot_macrostates(which = "initial", discrete=False, basis="X_merge_umap")
g.set_terminal_states(states=["5"])
g.plot_macrostates(which = "terminal", discrete = False, basis = "X_merge_umap")
g.compute_fate_probabilities()
model = cr.models.GAMR(adata, n_knots = 6, smoothing_penalty = 10.0)

# all ENSEMBL IDs
all_genes = ["ENSMUSG00000020932", "ENSMUSG00000026728", "ENSMUSG00000046805", "ENSMUSG00000002985", "ENSMUSG00000027797", "ENSMUSG00000057738", "ENSMUSG00000029309", "ENSMUSG00000060924", "ENSMUSG00000024109", "ENSMUSG00000035864", "ENSMUSG00000041014", "ENSMUSG00000005089", "ENSMUSG00000002107", "ENSMUSG0000006299", "ENSMUSG00000027210", "ENSMUSG00000061911", "ENSMUSG00000030209", "ENSMUSG00000066392"]
all_genes_intersection = list(set(all_genes) & set(adata.var_names))
all_genes_intersection
cr.pl.gene_trends(adata, 
                 model = model, 
                 genes = all_genes_intersection, 
                 same_plot = True, 
                 ncols = 2, 
                 time_key = "dpt_pseudotime", 
                 hide_cells=True)
cr.pl.heatmap(adata, 
             model = model, 
             lineages = "5", 
             cluster_key = "merge_clusters", 
             show_fate_probabilities=False, 
             genes = all_genes_intersection, 
             time_key = "dpt_pseudotime", 
             show_all_genes=True, 
             weight_threshold=(1e-3, 1e-3), 
             figsize=(7,10), 
             cmap = "RdBu_r", 
             cluster_genes = False, 
             dpi=300, 
             save="CellRank_DPT_Pseudotime_Heatmap.pdf"
             )

# PAGA
scv.tl.velocity_graph(adata)
scv.tl.paga(adata, groups = "merge_clusters", use_time_prior = "dpt_pseudotime")
scv.pl.paga(adata, basis = "merge_umap", size=10, alpha=1, color="merge_clusters", min_edge_width=2, node_size_scale=1.5, dpi=300, 
            save = "DPT_Pseudotime_paga.pdf", 
            figsize=(4, 4))
scv.pl.paga(adata, basis = "coords", size=10, alpha=1, color="merge_clusters", min_edge_width = 2, node_size_scale=1.5, dpi = 300, 
            save = "DPT_Pseudotime_paga_coords.pdf", 
            figsize=(4,4))