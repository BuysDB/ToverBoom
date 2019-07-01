import PIL
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import os



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


def construct_df_per_snv(snvData, cellCnv):
    for column in snvData.loc[list(cellCnv.index)].columns:
        # color snvs
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


def plot_per_snv(lg, cellData):
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
    plt.savefig(f"./output/{replicate}/{replicate}_{column[0]}_{column[1]}.png", dpi=300)


def per_replicate(replicate):
    create_topfolder(f"./output/{replicate}")

    # Instantiate the lineage graph object
    lg = toverboom.lineageGraph.LineageGraph(graph)
    # Create a figure for our plot:
    fig, ax = plt.subplots()
    # Find the best layout
    toverboom.optimizeLayout.optimize_layout(lg,
                                             visualize_progress_ax=ax,
                                             visualize_progress_fig=fig)
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
    imputed_visualize = gen_imputed_matrix_for_visualization(args.raw_snv_matrix, args.imputed_snv_matrix)
    snvData = imputed_visualize.loc[args.replicate].loc[args.cellCnv.index]
    cellData = cellCnv.copy()
    cellData['tp'] = [passage for passage, plate, cell in list(cellData.index)]
    # plot snv on cnv tree per snv
    for column in snvData.loc[list(cellCnv.index)].columns:
        cellData = plot_per_snv(snvData, cellCnv, column)
        plot_per_snv(lg, cellData)


def combine_images(combine_all= False):
    # find all file with the same snv name replicate here doesn't matter.
    replicate = 'APKS1'
    fullpiclist = []

    for root, dirs, files in os.walk(f"./output/{replicate}/", topdown=True):
        for i, name in enumerate(files):
            piclist = []
            rep, snvname = name.split(sep='_', maxsplit=1)
            if rep[:-1] != "APKS":
                pass
            else:
                print(snvname)
                apkslist = ['APKS1', 'APKS2', 'APKS3']
                for replicate in apkslist:
                    im = Image.open(f"./output/{replicate}/{replicate}_{snvname}")
                    piclist.append(im)
                    fullpiclist.append(im)
                    combined_image = pil_grid(piclist, 3)
                    create_topfolder(f"./output/combined/")
                    combined_image.save(f"./output/combined/{snvname}")

    if combine_all == True:
        full_image = pil_grid(fullpiclist, 3)
        full_image.save(f"./output/combined/all_snv.png")




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
    print(args.tree_graphs)
    if args.replicate == None:
        replicates = ['APKS1', 'APKS2', 'APKS3']
        for replicate in replicates:
            per_replicate(replicate)

        print("Generating comparison between replicates")
        combine_images(combine_all = args.combine_all)

    else:
        main(args.replicate)
        print(f"{args.replicate} figure generation finished!")


if __name__ == "__main__":
    main()




