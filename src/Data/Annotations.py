"""
Created on Jun 23, 2013

@author: agross
"""

from Stats.Scipy import chi2_cont_test
import pandas as pd


def read_in_pathways(mSigDbFile):
    """
    Reads in mSigDb pathways.
    File format should be that downloaded from website, tsv with format:
    pathway \t gene1 \t gene2 ... geneX \n
    input:
      mSigDbFile:        pathway file, can be downloaded from the web-site
    output:
      geneSets:          dict mapping pathway name to set of genes in the pathway
      geneLookup:        dict mapping genes to the pathways they occur in
    """
    f = open(mSigDbFile, 'r')
    geneSets = {}
    genes = []
    for line in f:
        line = line.replace('\"', '')
        tmp = line.strip().split("\t")
        setName = tmp[0]
        geneSets[setName] = set(tmp[2:])
        genes.extend(tmp[2:])
    f.close()
    genes = set(genes)

    geneLookup = dict([(gene, set()) for gene in genes])
    for pathway in geneSets:
        for gene in geneSets[pathway]:
            geneLookup[gene].add(pathway)
    return geneSets, geneLookup


def filter_pathway_hits(hits, gs, cutoff=.00001):
    """
    Returns a filtered list of p-values with redundant pathways
    removed.

    hits:
        Series of p-values.
    gs:
        DataFrame of gene-set asignments encoded as a binary matrix.
    cutoff:
        p-value cutoff for overlab between gene-sets
    """
    hits = hits.order()
    l = [hits.index[0]]
    for gg in hits.index:
        flag = 0
        for g2 in l:
            if gg in l:
                flag = 1
                break
            elif chi2_cont_test(gs[gg], gs[g2])['p'] < cutoff:
                flag = 1
                break
        if flag == 0:
            l.append(gg)
    hits_filtered = hits.ix[l]
    return hits_filtered


def unstack_geneset_csv(filename):
    """
    File should be a stacked binary table of gene set assignments.
    """
    s = pd.read_csv(filename, header=None, index_col=[0, 1], squeeze=True)
    df = s.unstack().T.fillna(0)
    df.columns.name = 'Gene_Set'
    df.index.name = 'Gene'
    return df

