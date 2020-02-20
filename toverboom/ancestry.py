import networkx as nx
import pandas as pd
import numpy as np
import itertools
import collections

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

    a = cnvStates.loc[ cnvStates['cluster']==cnvStateA,:][['binIndex','chromosome', 'cluster','copyNumber','endCoordinate','startCoordinate']]
    a.index = pd.MultiIndex.from_frame(a[['chromosome','binIndex']])
    b = cnvStates.loc[ cnvStates['cluster']==cnvStateB,:][['binIndex','chromosome', 'cluster','copyNumber','endCoordinate','startCoordinate']]
    b.index = pd.MultiIndex.from_frame(b[['chromosome','binIndex']])

    merge = a.join(b,lsuffix='_A',rsuffix='_B')

    distancesPerChromosome = collections.Counter()
    for i,row in merge.iterrows():
        distance = row.copyNumber_B - row.copyNumber_A
        # We cannot make a chromosome out of nothing
        if row.copyNumber_B>0 and  row.copyNumber_A==0:
            return None,None
        if distance != 0:
            distancesPerChromosome[(row.chromosome_A,distance)]+=1

    return len(distancesPerChromosome), distancesPerChromosome


def create_ancestry_graph(cnvStates, maxCNVdistance= 3, select_clusters=None):
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

    agraph = nx.drawing.nx_agraph.to_agraph(G)
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

def visualize_distance_graph(ancestry, prefix):
    """
    visualize ancestry graph using pygraphviz
    """

    agraph = nx.drawing.nx_agraph.to_agraph(ancestry)
    G = ancestry
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
    agraph.write(f'{prefix}_cnv_distance.dot')
    agraph.draw(f'{prefix}_cnv_distance.png')


def time_expand(reducedGraph, replicate, replicateToStateCounter,timePunish=0.001):
    
    expandedGraph = nx.DiGraph()

    allTimePoints = set([-5,0])
    for node in reducedGraph:
        timePoints = [
            pas for (rep, pas), obs in replicateToStateCounter[node].items()
            if rep==replicate or replicate=='ALL']
        allTimePoints.update(set(timePoints))

    allTimePoints = sorted(list(allTimePoints))
    cnvToTimePointMapping = {}
    cloneExtinction = {} # CNV -> time of extinction
    for node in reducedGraph:
        #Obtain timepoints for node at x
        timePoints = [ pas for (rep, pas), obs in replicateToStateCounter[node].items() if rep==replicate or replicate=='ALL']
        if node==0: # Prior
            #inOtherReplicates=None
            earliestTimepoint =-5
            lastTimePoint=0
        else:
            inOtherReplicates = [ pas for (rep, pas), obs in replicateToStateCounter[node].items() if rep!=replicate and obs>1]
            earliestTimepoint = min(timePoints)
            lastTimePoint = max(timePoints)

            if len(inOtherReplicates)>1 and createSharedTimepoint:
                earliestTimepoint = 0 # Already present

        cnvToTimePointMapping[node] =[timePoint for timePoint in allTimePoints
            if timePoint>=earliestTimepoint and timePoint<=lastTimePoint]

        if lastTimePoint<max(allTimePoints):
            cloneExtinction[node] = lastTimePoint

    for cnvState, timePoints in cnvToTimePointMapping.items():
        prev = None # Previous cnv state and timepoint
        for timePoint in timePoints:
            current = (cnvState, timePoint)
            expandedGraph.add_node( current )
            if prev:
                expandedGraph.add_edge(prev, current, distance=0)
            prev = current

        for tp in allTimePoints:
            if (u,tp) in expandedGraph:
                fromNode = (u,tp)
                break
        prevTimePoint = None
        for tp in allTimePoints:
            if (v,tp) in expandedGraph:
                toNode = (v, tp)
                break
            prevTimePoint=tp

        if (u, prevTimePoint) in expandedGraph:
            fromNode = (u, prevTimePoint)

    ## Add ancestry edges:
    edgesToPut = reducedGraph.edges()
    for copyStateA,copyStateB,d in reducedGraph.edges(data=True):


        for fromNode in [ node for node in expandedGraph if node[0]==copyStateA  ]:
            # To node should be the smallest timepoint
            toNode = (copyStateB, min(cnvToTimePointMapping[copyStateB]))

            # Fromnode needs to be earlier than to-node
            if fromNode[1]<=toNode[1]:
                 expandedGraph.add_edge(fromNode, toNode, desc=d.get('desc'), distance=d['distance'] - (timePunish*fromNode[1]) + 10*(fromNode[1]==toNode[1]))

    for copyState, timePoint in expandedGraph:
        if cloneExtinction.get(copyState, None)==timePoint:
            expandedGraph.nodes[(copyState, timePoint)]['extinct'] = True
        if max(cnvToTimePointMapping[copyState])==timePoint:
            expandedGraph.nodes[(copyState, timePoint)]['leaf'] = True

    return expandedGraph


def add_wildtype_state(df):
    """
    Add a wiltype state to the dataframe

    For chromosomes with an "_" in the name the wt copy number is expected to be 1
    2 otherwise

    Args:
        df (pd.Dataframe) Expected dataframe columns:
        chromosome, binIndex,startCoordinate, endCoordinate, cluster, copyNumber

    Returns:
        DataFrame
    """
    zero_state = []
    for (chromosome, binIndex,startCoordinate, endCoordinate ),group in df.groupby(['chromosome','binIndex','startCoordinate', 'endCoordinate']):

        zero_state.append({'chromosome':chromosome,
                           'binIndex':binIndex,
                           'startCoordinate':startCoordinate,
                           'endCoordinate':endCoordinate,
                            'cluster':0,
                           'copyNumber':2 if chromosome.count('_')==0 else 1
                          })

    return pd.concat( (df,pd.DataFrame(zero_state) ), sort=True)
