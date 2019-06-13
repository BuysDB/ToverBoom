import networkx as nx
import pandas as pd
import numpy as np
import itertools

def cnvDistance(cnvStates, cnvStateA, cnvStateB ):
    """
    Returns the total distance and edits per chromosome from cnvStateA to cnvStateB
    Returns None if it is not possible to go from  cnvStateA to cnvStateB because of
    an allele existing in B but not in A.

    Parameters
    ----------
    cnvStates : pandas.DataFrame
        DataFrame containing the columns 'binIndex' (int),'chromosome' (str), 'cluster' (int),'copyNumber' (int),'endCoordinate' (int),'startCoordinate' (int)
    cnvStateA : int
        Integer value which points to a cluster id in cnvStates['cluster']
    cnvStateB : int
        Integer value which points to a cluster id in cnvStates['cluster']
    """

    a = cnvStates.loc[ cnvStates['cluster']==cnvStateA,:][['binIndex','chromosome', 'cluster','copyNumber','endCoordinate','startCoordinate']].set_index('binIndex')
    b =  cnvStates.loc[ cnvStates['cluster']==cnvStateB,:][['binIndex','cluster','copyNumber']].set_index('binIndex')
    merge = a.join(b, lsuffix='_A', rsuffix='_B')
    distancesPerChromosome = collections.Counter()
    for i,row in merge.iterrows():
        distance = row.copyNumber_B - row.copyNumber_A
        # We cannot make a chromosome out of nothing
        if row.copyNumber_B>0 and  row.copyNumber_A==0:
            return None,None
        if distance != 0:
            distancesPerChromosome[(row.chromosome,distance)]+=1

    return len(distancesPerChromosome), distancesPerChromosome


def create_ancestry_graph(cnvStates, maxCNVdistance= 3, select_clusters=None)
    """
    Create a directional graph where every node is an unique copy number state.
    Edges indicate a possibility to traverse between copy number state with a
    cost. The cost of traversing between states is defined by the cnvDistance
    function.

    Parameters
    ----------
    cnvStates : pandas.DataFrame
        DataFrame containing the columns 'binIndex' (int),'chromosome' (str), 'cluster' (int),'copyNumber' (int),'endCoordinate' (int),'startCoordinate' (int)
        'chromosome' values suffixed with _A or _B are seen as specific alleles of chromosome.

    maxCNVdistance : int
        maximum edit distance to keep in the ancestry graph

    select_clusters : iterable or None
        copy number states/clusters to only use. None to use all clusters in cnvStates
    """

    clusters_to_use = select_clusters if select_clusters is not None else cnvStates['cluster'].unique()
    ancestry = nx.MultiDiGraph()
    for cnvStateA,cnvStateB in itertools.product(

        (
            clusters_to_use
        )
        , repeat=2):
        if cnvStateA!=cnvStateB:
            d, changes = cnvDistance(cnvStates, cnvStateA, cnvStateB)
            if d is not None and d<=maxCNVdistance:
                # Convert edge into readable format
                desc = ', '.join([f'{chrom} loss' if change<0 else
                                  f'{chrom} gain' for (chrom, change),obs in changes.items() ])
                ancestry.add_edge( cnvStateA,cnvStateB, distance=d, dtype='cnv', desc=desc )
    return ancestry
    #####

def visualize_ancestry(G, target_dot_path=None, target_png_path=None):

    #nx.write_gpickle(G,f'{replicate}_cnv_distance.gpickle')

    agraph = nx.drawing.nx_agraph.to_agraph(ancestry)
    for i,(a,b) in enumerate(agraph.edges_iter()):
        for mindex in G[int(a)][int(b)]:
            if G[int(a)][int(b)][mindex]['dtype']=='ssnv':
                d = G[int(a)][int(b)][mindex]['distance']
                agraph.get_edge( int(a), int(b)).attr['color'] = 'blue'
                agraph.get_edge( int(a), int(b), mindex).attr['label'] = '%.2f' % d
            elif G[int(a)][int(b)][mindex]['dtype']=='cnv':
                d = G[int(a)][int(b)][mindex]['distance']
                desc = G[int(a)][int(b)][mindex]['desc']
                agraph.get_edge( int(a), int(b), mindex).attr['label'] = f'{desc}: {int(d)}'

    agraph.layout(prog='dot')
    if target_dot_path is not None:
        agraph.write(target_dot_path)
    if target_png_path is not None:
        agraph.draw(target_png_path)
