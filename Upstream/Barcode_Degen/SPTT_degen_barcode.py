import gzip
import logging
from collections import Counter
from pathlib import Path
import click
import networkx as nx
import numpy as np
import scipy.sparse
from sklearn.neighbors import radius_neighbors_graph
import pandas as pd
import pickle
import argparse

'''
inspired by Slide-seqV2
'''

def args_parse():
    parser = argparse.ArgumentParser(description="Degen beads barcodes.", add_help=True)
    parser.add_argument("-bc_loc", type = str, help = "The Barcode Location file.")
    parser.add_argument("-seq_barcodes", type = str, help = "The .txt file containing Top Beads Barcodes of raw reads.")
    parser.add_argument("-workdir", type = str, help = "The working aka the saving directory.")
    return parser.parse_args()

def degen_barcode(barcode_set):
    return "".join(s.pop() if len(s) == 1 else "N" for s in map(set, zip(*barcode_set)))

def initial_h_set(*barcodes):
    base_d = {"A": 0, "C":1, "G": 2, "T":3, "N":4}
    return {tuple(base_d[c] for c in barcode) for barcode in barcodes}
def hamming_set(h_set,  d = 1, include_N = False):
    h_set = h_set.copy()
    bc_len = len(next(iter(h_set)))

    new_base = [i * np.eye(bc_len, dtype = np.uint8) for i in range(4 + include_N) ]
    other_bases = 1- np.eye(bc_len, dtype = np.uint8)

    for _ in range(d):
        for a in list(map(np.array, h_set)):
            h_set.update(t for i in new_base for t in map(tuple, a * other_bases + i))
    h_set = {"".join("ACGTN"[i] for i in h) for h in h_set}
    return h_set

def hamming1_adjacency(barcodes):
    assert len(barcodes) == len(set(barcodes)), "barcodes should be unique"
    bc_to_i = {bc: i for i, bc in enumerate(barcodes)}
    adjacency = scipy.sparse.dok_matrix((len(barcodes), len(barcodes)), dtype = int)

    for i, bead_bc in enumerate(barcodes):
        for bc1 in hamming_set(initial_h_set(bead_bc), d = 1, include_N = True):
            if bc1 in bc_to_i and bead_bc != bc1:
                adjacency[i, bc_to_i[bc1]] = 1

    return adjacency.tocsr()

def bipartite_matching(
    bead_barcodes, 
    degen_barcodes, 
    bead_groups, seq_barcodes): 
    bead_lens = set(map(len, bead_barcodes))
    seq_lens = set(map(len, seq_barcodes))
    assert ( bead_lens == seq_lens), f"Beads have length {bead_lens} but sequenced lengths {seq_lens} "

    assert {c for bc in bead_barcodes for c in bc}.issubset("ACGTN")
    seq_nset = {c for bc in seq_barcodes for c in bc}
    assert seq_nset.issubset("ACGTN")

    include_N = "N" in seq_nset

    ## construct bipartite graph: barcode to hamming set
    matching_graph = nx.Graph()
    for i, bead_bg in enumerate(bead_groups):
        matching_graph.add_node(i, bipartite = 0)
        h_set = initial_h_set(*(bead_barcodes[j] for j in bead_bg))
        for bc1 in hamming_set(h_set, d = 1, include_N = include_N):
            matching_graph.add_node(bc1, bipartite = 1)
            matching_graph.add_edge(i, bc1)

    barcode_mapping = dict()
    for seq_bc in seq_barcodes:
        if seq_bc in matching_graph:
            sample_set = set(matching_graph[seq_bc])
            if len(sample_set) == 1:
                barcode_mapping[seq_bc] = degen_barcodes[sample_set.pop()]
    return barcode_mapping

def write_barcode_mapping(barcode_mapping, bead_xy, output_file):
    with gzip.open(output_file, "wt") as out:
        for seq_bc, bead_bc in barcode_mapping.items():
            x, y = bead_xy[bead_bc]
            print(f"{seq_bc}\t{bead_bc}\t{x:.1f}\t{y:.1f}", file = out)

def write_barcode_xy(bead_list, bead_xy, output_file):
    with gzip.open(output_file, "wt") as out:
        for bead_bc in bead_list:
            x, y = bead_xy[bead_bc]
            print(f"{bead_bc}\t{x:.1f}\t{y:.1f}", file = out)

def main(beads_location_file: str, seq_barcode_file:str, workdir:str):
    bead_location_file = beads_location_file
    bead_location = pd.read_csv(bead_location_file, header = None, names = ["new_x", "new_y", "barcode"])
    x_line = np.array(bead_location[bead_location.columns[0:1]].values)
    y_line = np.array(bead_location[bead_location.columns[1:2]].values)
    bead_barcodes = list(bead_location.iloc[:, 2])

    xy = np.hstack((x_line, y_line))
    if xy.shape[0] == len(bead_barcodes):
        print("Read " + str(xy.shape[0]) + " bead locations for " + str(len(bead_barcodes)) + " beads")
    ## pre-emptively remove poly-T/N sequences
    ok_barcodes = [not set(bc).issubset({"T", "N"}) for bc in bead_barcodes]
    xy = xy[ok_barcodes, :]
    bead_barcodes = [bc for ok, bc in zip (ok_barcodes, bead_barcodes) if ok]
    print(len(ok_barcodes)) 
    #### adjacency matrix for all beads within radius of each other
    radius = 10
    radius_matrix = radius_neighbors_graph(xy, radius = radius)
    ## adjacency matrix for all barcodes within hamming distance 1
    hamming_matrix = hamming1_adjacency(bead_barcodes)
    # just multiply together to get the combined adjacency matrix
    combined_graph = nx.from_scipy_sparse_array(radius_matrix.multiply(hamming_matrix))
    ## add xy coordinates to graph so we can analyze later
    for n, (x, y) in zip(combined_graph.nodes, xy):
        combined_graph.nodes[n]["x"] = x
        combined_graph.nodes[n]["y"] = y
    
    ## get connected components to find groups of similar/close barcodes
    bead_groups = list(nx.connected_components(combined_graph))

    # calculate degenerate (ambiguous bases -> N) barcodes
    degen_bead_barcodes = [degen_barcode({bead_barcodes[j] for j in bg}) for bg in bead_groups]
    print("Collapsed " + str(len(bead_groups)) + " bead groups into ")
    print(str(len(set(degen_bead_barcodes))) + " barcodes" )

    ## just in case, we'll add integer tags to each one so they are unique
    barcode_counter = Counter()
    for i, barcode in enumerate(degen_bead_barcodes):
        barcode_counter[barcode] += 1
        degen_bead_barcodes[i] = f"{barcode}-{barcode_counter[barcode]}"
    
    # average xy for grouped beads to get centroids
    bead_xy = dict()
    for bg, degen_bc in zip(bead_groups, degen_bead_barcodes): 
        bg_graph = combined_graph.subgraph(bg)
        mean_x, mean_y = np.array(
        [[nd["x"], nd["y"]] for _, nd in bg_graph.nodes(data = True)]).mean(0)
        bead_xy[degen_bc] = (mean_x, mean_y)

    seq_barcode_file = seq_barcode_file
    seq_barcodes_f = open(seq_barcode_file)
    seq_barcodes = list(map(lambda x:x.strip(), seq_barcodes_f.readlines()))

    barcode_matching = bipartite_matching(bead_barcodes, degen_bead_barcodes, bead_groups, seq_barcodes)

    ### filter degen_bead_barcodes to only barcodes that were observed
    matched_set = set(barcode_matching.values())
    degen_bead_barcodes = [bc for bc in degen_bead_barcodes if bc in matched_set]

    barcode_list = degen_bead_barcodes
    barcode_mapping = barcode_matching
    bead_xy = bead_xy
    bead_graph = combined_graph

    output_file = workdir + "/barcode_mapping.txt.gz"
    write_barcode_mapping(barcode_mapping, bead_xy, output_file)

    with open(workdir + "/bead_graph.gpickle", "wb") as f:
        pickle.dump(bead_graph, f)

    return

if __name__ == "__main__":
    args = args_parse()
    main(beads_location_file=args.bc_loc, 
         seq_barcode_file=args.seq_barcodes, 
         workdir=args.workdir)