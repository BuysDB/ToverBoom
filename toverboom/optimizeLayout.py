#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import itertools

def swap(indexA,indexB, listToSwap):
    return [
            listToSwap[indexA] if i==indexB else
            (listToSwap[indexB] if i==indexA else x)
        for i,x in enumerate(listToSwap)]

def optimize_layout(
            lineage_graph,
            plotArgs = {
                        'xDistance':1,
                        'verticalSpacing':0.1,
                        'defaultRadius' :0.3,
                        'radiusAttribute' : 'ratio',
                        'radiusMultiplier':8,
                        'minRadius':0.3,
                        'initRadius':0.3
                        } ,
                        visualize_progress_fig = None,
                        visualize_progress_ax = None
                    ):
        if visualize_progress_ax is not None and visualize_progress_ax is None:
            raise ValueError("Supply the figure")

        if visualize_progress_ax is None:
            plot=False
            fig = None
            ax = None
        else:
            ax = visualize_progress_ax
            fig = visualize_progress_fig
            plot=True

        # Load the topology:
        graph = lineage_graph.graph

        # The order of the clones (trellis) is stored here:
        trellisOrder = list(set([lineage for lineage,timepoint in lineage_graph.graph]))

        # Obtain the score for the first trellis:
        result = lg.optimizeGraph( graph, tuple(trellisOrder), plot=plot, ax=ax, **plotArgs)

        curScore, intersectingClones = result


    prevBest = None
    totalSwaps = 0
    didSwap=True
    while didSwap:
        didSwap=False
        prevBest=best
        for iA, iB in itertools.combinations(
            [ best.index(clone) for clone in intersectingClones ],2):

            if iA!=iB:
                newTrellisOrder = swap(iA,iB, best )
                totalSwaps+=1
                result = lg.optimizeGraph( graph, trellisOrder=tuple(newTrellisOrder), plot=False, quitAtScore=curScore, **plotArgs)
                if result is not None:
                    score, intersectingClones = result
                if result is not None and (curScore is None or score<curScore):
                    best = newTrellisOrder
                    if plot:
                        ax.clear()
                    didSwap=True
                    curScore, intersectingClones = lg.optimizeGraph( graph, trellisOrder=tuple(newTrellisOrder), plot=plot, ax=ax, **plotArgs)

                    if plot:
                        ax.set_title(f'{curScore} iter:{totalSwaps}' )
                        fig.canvas.draw()
                    #plt.pause(0.001)
                    break

        for iA, iB in itertools.product(
            [ best.index(clone) for clone in intersectingClones ],
            list(range(len(best)))):
            if iA!=iB:
                newTrellisOrder = swap(iA,iB, best )
                totalSwaps+=1
                result = lg.optimizeGraph( graph, trellisOrder=tuple(newTrellisOrder), plot=False, quitAtScore=curScore, **plotArgs)
                if result is not None:
                    score, intersectingClones = result
                if result is not None and (curScore is None or score<curScore):
                    best = newTrellisOrder
                    if plot:
                        ax.clear()
                    didSwap=True
                    curScore, intersectingClones = lg.optimizeGraph( graph, trellisOrder=tuple(newTrellisOrder), plot=plot, ax=ax, **plotArgs)

                    if plot:
                        ax.set_title(f'{curScore} iter:{totalSwaps}' )
                        fig.canvas.draw()
                    #plt.pause(0.001)
                    break

        for iA, iB in itertools.combinations( list(range(len(best))),r=2 ):
            if iA!=iB:
                newTrellisOrder = swap(iA,iB, best )
                totalSwaps+=1
                result = lg.optimizeGraph( graph, trellisOrder=tuple(newTrellisOrder), plot=False, quitAtScore=curScore, **plotArgs)
                if result is not None:
                    score, intersectingClones = result
                if result is not None and (curScore is None or score<curScore):
                    best = newTrellisOrder
                    if plot:
                        ax.clear()
                    didSwap=True
                    curScore, intersectingClones = lg.optimizeGraph( graph, trellisOrder=tuple(newTrellisOrder), plot=plot, ax=ax, **plotArgs)

                    if plot:
                        ax.set_title(f'{curScore} iter:{totalSwaps}' )
                        fig.canvas.draw()
                    #plt.pause(0.001)
                    break

    if plot:
        ax.clear()
    didSwap=True
    curScore = lg.optimizeGraph( graph, trellisOrder=tuple(best), plot=plot, ax=ax, storeParams=True, **plotArgs)
    if plot:
        ax.set_title(f'{curScore} iter:{totalSwaps}' )
        ax.set_aspect('equal')
        fig.canvas.draw()
    return best
