import networkx as nx
import pandas as pd
import numpy as np

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
