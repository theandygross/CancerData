"""
Created on Jun 23, 2013

@author: agross
"""


def read_in_pathways(mSigDbFile):
    '''
    Reads in mSigDb pathways.
    File format should be that downloaded from website, tsv with format:
    pathway \t gene1 \t gene2 ... geneX \n
    input:
      mSigDbFile:        pathway file, can be downloaded from the web-site
    output:
      geneSets:          dict mapping pathway name to set of genes in the pathway
      geneLookup:        dict mapping genes to the pathways they occur in
    '''    
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

