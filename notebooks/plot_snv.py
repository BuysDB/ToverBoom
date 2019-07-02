import PIL.Image as Image
import argparse
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import os
import networkx as nx
import toverboom
import toverboom.lineageGraph
import toverboom.optimizeLayout
import toverboom.preprocessing
# import importlib
# importlib.reload(toverboom)
# importlib.reload(toverboom.lineageGraph)
# importlib.reload(toverboom.preprocessing)
# importlib.reload(toverboom.optimizeLayout)


def create_topfolder(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def gen_imputed_matrix_for_visualization(raw_matrix, imputed_matrix, transparency=0.45):
    '''
    raw values -> 0,1
    imputed values -> 0.45, 0.55
    return post_impute_matrix
    = a matrix with graded values representing the evidence of values (can be even more! different round of imputation!)
    '''
    # 0. Reduce to the same matrix size
    try:
        raw_matrix = raw_matrix.loc[imputed_matrix.index]
    except Exceptions:
        print(Exceptions)


    # 1. check if only include 0, 1. If -1: convert to np.nan
    raw_matrix[(raw_matrix == -1)] = np.nan
    raw_matrix[(raw_matrix > 0.5)] = 1
    raw_matrix[(raw_matrix < 0.5)] = 0

    # 2. enhance the transparency of the imputed values
    transparent_matrix = imputed_matrix.copy()
    transparent_matrix[(transparent_matrix == 0)] = transparency
    transparent_matrix[(transparent_matrix == 1)] = (1 - transparency)

    # 3. combine the measured values and imputed values in a graded matrix
    combined_matrix = transparent_matrix.copy()
    combined_matrix.update(raw_matrix)

    return combined_matrix


def test_visual():
    raw = pd.DataFrame.from_dict({'col1': [0, np.nan, 0, np.nan], 'col2': [1, np.nan, np.nan, 1]})
    imputed = pd.DataFrame.from_dict({'col1': [0.45, 0.45, 0.45, 0.45], 'col2': [0.55, 0.55, np.nan, 0.55]})
    combined = pd.DataFrame.from_dict(
        {'col1': {0: 0.0, 1: 0.45, 2: 0.0, 3: 0.45}, 'col2': {0: 1.0, 1: 0.55, 2: nan, 3: 1.0}})

    if combined == gen_imputed_matrix_for_visualization(raw, imputed):
        pass
    else:
        raise "Error"
    return None


def construct_df_per_snv(cellData, snvData, column):
    cellData['color'] = [{0: '#1d2bf7', 0.45: '#5b94ff',
                          1: 'r', 0.55: '#ff7575'}.get(cluster, 'grey') for cluster in snvData[column]]

    # Assign markers
    cellData['marker'] = [{0: 'o', 0.45: 'o', 1: 's', 0.55: 's'}.get(cluster, '.') for cluster in snvData[column]]
    cellData['size'] = [{0: 28, 0.45: 20,
                         1: 28, 0.55: 20}.get(cluster, 3) for cluster in snvData[column]]
    # drawing order
    cellData['drawing order'] = [{1: 1, 0: 2, 0.55: 3, 0.45: 4}.get(cluster, 5) for cluster in snvData[column]]
    cellData.sort_values(by='drawing order', ascending=False, inplace=True, na_position='last')
    return cellData


def plot_per_snv(lg, cellData, replicate, column, output):
    fig, ax = lg.getEmptyPlot()
    lg.plotPatches(ax, facecolor=(0.8, 0.8, 0.8, 1))
    lg.plotSingleCells(cellData, ax=ax, fig=fig, plotPatches=False, enableShadow=True)

    # Add labels to the clones:
    lg.annotateNodes(ax, plotArgs={'size': 9},
                     # Use the nodesToAnnotate argument to select which nodes to annotate
                     nodesToAnnotate=[
                         (cluster, tp)
                         for cluster, tp in lg.graph if cluster in [1, 4, 3, 5, 2, 8]],
                     x_offset=5  # How much to move the labels to the right
                     )

    # Add vertical lines to indicate sampled timepoints
    lg.plot_vertical_lines(ax, cellData['tp'].unique(), c='black')

    lg.plot_xticks(ax, cellData['tp'].unique())

    ax.set_xlabel('Time (weeks)')
    plt.title(f"{replicate} {column[0]}:{column[1]}")
    plt.savefig(f"{output}/{replicate}/{replicate}_{column[0]}_{column[1]}.png", dpi=300)
    plt.close(fig)

def per_replicate(replicate, cnv_tree, raw_snv_matrix, imputed_snv_matrix, cellCnv, output):
    create_topfolder(f"{output}/{replicate}")
    # Instantiate the lineage graph object
    lg = toverboom.lineageGraph.LineageGraph(cnv_tree)
    # Create a figure for our plot:
    fig, ax = plt.subplots()
    # Find the best layout
    trellisOrder = toverboom.optimizeLayout.optimize_layout(lg,
                                             visualize_progress_ax=ax,
                                             visualize_progress_fig=fig,
                                             initial_order=None)
    # Plot the polygons of the tree
    fig, ax = plt.subplots()
    # wavyness controls how wavy the segments of the tree are
    wavyness = 0.4
    # xDistance controls the stretch in the x direction
    lg.xDistance = 10
    lg.verticalSpacing = 0.1
    lg.plotEdges(ax, bezier=True, wavyness=wavyness, stepCount=30, plotArgs={'linewidth': 1}, offsetCentroid=True)
    lg.plotPatches(ax=ax, wavyness=wavyness)
    # Remove plot spines:
    toverboom.lineageGraph.despine(ax)
    # Scale labels to plot size:
    toverboom.lineageGraph.format_x_axis_labels(ax)
    # Add labels to the clones:
    lg.annotateNodes(ax, plotArgs={'size': 8})
    fig.canvas.draw()

    # combine imputed and raw snv matrix
    snvData = gen_imputed_matrix_for_visualization(raw_snv_matrix, imputed_snv_matrix)
    snvData = snvData.loc[replicate]
    # start constructing plotting dataframe by copying cellCnv
    cellData = cellCnv.loc[replicate]
    cellData['tp'] = [passage for passage, plate, cell in list(cellData.index)]
    # Overlapped df
    snvData = snvData.loc[cellData.index]
    # plot snv on cnv tree per snv
    for column in snvData.loc[list(cellData.index)].columns:
        plotData = construct_df_per_snv(cellData, snvData, column)
        plot_per_snv(lg, plotData, replicate, column, output)
        break


def pillow_grid(images, max_horiz=np.iinfo(int).max):
    n_images = len(images)
    n_horiz = min(n_images, max_horiz)
    h_sizes, v_sizes = [0] * n_horiz, [0] * ((n_images // n_horiz) + (1 if n_images % n_horiz > 0 else 0))
    for i, im in enumerate(images):
        h, v = i % n_horiz, i // n_horiz
        h_sizes[h] = max(h_sizes[h], im.size[0])
        v_sizes[v] = max(v_sizes[v], im.size[1])
    h_sizes, v_sizes = np.cumsum([0] + h_sizes), np.cumsum([0] + v_sizes)
    im_grid = Image.new('RGB', (h_sizes[-1], v_sizes[-1]), color='white')
    for i, im in enumerate(images):
        im_grid.paste(im, (h_sizes[i % n_horiz], v_sizes[i // n_horiz]))
    return im_grid


def combine_images(output, combine_all= False):
    # find all file with the same snv name replicate here doesn't matter.
    replicate = 'APKS1'
    fullpiclist = []

    for root, dirs, files in os.walk(f"{output}/{replicate}/", topdown=True):
        for i, name in enumerate(files):
            piclist = []
            rep, snvname = name.split(sep='_', maxsplit=1)
            if rep[:-1] != "APKS":
                pass
            else:
                print(snvname)
                apkslist = ['APKS1', 'APKS2', 'APKS3']
                for replicate in apkslist:
                    im = Image.open(f"{output}/{replicate}/{replicate}_{snvname}")
                    piclist.append(im)
                    fullpiclist.append(im)
                    combined_image = pillow_grid(piclist, 3)
                    create_topfolder(f"{output}/combined/")
                    combined_image.save(f"{output}/combined/{snvname}")

    if combine_all == True:
        full_image = pillow_grid(fullpiclist, 3)
        full_image.save(f"{output}/combined/all_snv.png")




parser = argparse.ArgumentParser(description='Plot SNV on CNV tree. Provide path and select if all replicate should be plotted or only one')
parser.add_argument('-r', '--replicate', default=None, type = str)
parser.add_argument('-t', '--tree_graphs', required=True,  nargs='+') # list of trees in pickle format
parser.add_argument('-s', '--raw_snv_matrix', required=True, type = str)
parser.add_argument('-si', '--imputed_snv_matrix', required=True, type = str) # imputed snv matrix
parser.add_argument('-cnv', '--cellCnv', required=True, type = str) # cnv mapping
parser.add_argument('-o', '--output', required=True, type = str)  # output directory
parser.add_argument('-co', '--combine_all', default = False, type = bool)  # output directory
args = parser.parse_args()




def main():
    imputed_snv_matrix = pd.read_pickle(args.imputed_snv_matrix)
    raw_snv_matrix = pd.read_pickle(args.raw_snv_matrix).loc[imputed_snv_matrix.index]
    cellCnv = pd.read_pickle(args.cellCnv)


    if args.replicate == None:
        replicates = ['APKS1', 'APKS2', 'APKS3']
        for i, replicate in enumerate(replicates):
            cnv_tree_graph = nx.read_gpickle(args.tree_graphs[i])
            per_replicate(replicate, cnv_tree_graph, raw_snv_matrix, imputed_snv_matrix, cellCnv, args.output)

        print("Generating comparison between replicates")
        combine_images(args.output, combine_all = args.combine_all)

    else:
        if len(args.tree_graphs) == 1:
            cnv_tree_graph = nx.read_gpickle(args.tree_graphs[0])
            per_replicate(args.replicate, cnv_tree_graph, raw_snv_matrix, imputed_snv_matrix, cellCnv, args.output)
            print(f"{args.replicate} figure generation finished!")
        else:
            print("Error: -t received >1 tree graphs. Please modify the parameter.")
            exit()

if __name__ == "__main__":
    main()




# ./data/APKS2_CNV_tree_33.pickle ./data/APKS3_CNV_tree_33.pickle