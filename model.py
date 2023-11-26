#!/usr/bin/env python
# coding: utf-8

import itertools as it
import json
import pickle

import markov_clustering as mc
import networkx as nx
from networkx.readwrite import json_graph

with open('../data/axen_loclass_input.pickle', 'rb') as f:
    domains = pickle.load(f)

pfams = {}
for loci in domains:
    pfams[loci] = set(item for lst in domains[loci] for item in lst if ('PF00936' not in lst and 'PF03319' not in lst))

identifying_pfams = {'PF12288',
                     'PF08936',
                     'PF00132',
                     'PF06751',
                     'PF05985',
                     'PF02286',
                     'PF02287',
                     'PF02288',
                     'PF00596'}

weights = {}
for loci in domains:

    for index, gene in enumerate(domains[loci]):
        if 'PF00936' in gene or 'PF03319' in gene:
            envelope_start = index
            break
    for index, gene in reversed(list(enumerate(domains[loci]))):
        if 'PF00936' in gene or 'PF03319' in gene:
            envelope_end = index
            break

    weights[loci] = {}
    for idx, o in enumerate(domains[loci]):
        if idx > envelope_end:
            dis = idx - envelope_end
        elif idx < envelope_start:
            dis = envelope_start - idx
        else:
            dis = 0
        for domain in o:
            weights[loci][domain] = dis

rarity = {}

for loci in domains:
    for index, pf in enumerate(domains[loci]):
        for pfam in pf:
            if pfam not in rarity.keys():
                rarity[pfam] = len([i for i, P in pfams.items() if pfam not in P]) / len(pfams.keys())


def identifying(p):
    if p in identifying_pfams:
        return 3
    else:
        return 1


def distance(d):
    return max(1 - 0.1 * d, 0.1)


def coocurrence(p, c):
    return len([i for i, P in pfams.items() if c <= P and p in P]) / len([i for i, P in pfams.items() if c <= P])


def similarity(i, j):
    common = pfams[i] & pfams[j]
    c_score = sum(identifying(p)
                  * rarity[p]
                  * distance(min(weights[i][p], weights[j][p])) for p in common)
    d_score = sum(identifying(p)
                  * rarity[p]
                  * distance(weights[i][p])
                  * coocurrence(p, common) for p in pfams[i] - pfams[j])

    d1_score = sum(identifying(p)
                   * rarity[p]
                   * distance(weights[j][p])
                   * coocurrence(p, common) for p in pfams[j] - pfams[i])

    return c_score - .5 * (d_score + d1_score)


def main():
    graph = nx.Graph()

    for pair in list(it.combinations(domains.keys(), 2)):
        score = similarity(pair[0], pair[1])
        if score > 1:
            graph.add_edge(pair[0], pair[1], weight=score)
    return graph


if __name__ == '__main__':
    G = main()

matrix = nx.to_scipy_sparse_matrix(G)
result = mc.run_mcl(matrix)
clusters = mc.get_clusters(result)

cmap = {}
for i, key in enumerate(G.nodes):
    for j, value in enumerate(clusters):
        if i in value:
            cmap[key] = "group" + str(j)
            break

nx.set_node_attributes(G, cmap, "group")

with open('../data.json', 'w') as outfile1:
    outfile1.write(json.dumps(json_graph.node_link_data(G)))
