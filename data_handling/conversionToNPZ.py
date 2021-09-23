import pandas as pd
from collections import namedtuple
import numpy as np
from recordtype import recordtype
from collections import OrderedDict
import os

base_prefix = 'testGraph'
graph_dir = '../graphs/'.format(base_prefix)
print(graph_dir)
print(os.listdir(graph_dir))
graph_paths = [os.path.join(graph_dir, i) for i in os.listdir(graph_dir) if (i != "graphsNPZ" and not i.startswith('.'))]
output_dir = '../graphs/graphsNPZ/'
for graph_num, graph_path in enumerate(graph_paths):
    print(graph_path)
    df = pd.read_csv(graph_path, delimiter=',')
    dfNP = df.to_numpy()
    if(len(dfNP) == 0):
        continue
    column_names = ['x_0', 'x_1', 'x_2', 'edge_attr_0', 'edge_attr_1', 'edge_attr_2', 'edge_attr_3', 'edge_index_1', 'edge_index_2', 'y', 'pid']
    column_names_string = ' '.join(column_names[i] for i in range(11))
    column_dict = {i:name for i, name in enumerate(column_names)}
    dictNoNans = OrderedDict((i,i) for i in column_names)
    dfNoNans = pd.DataFrame()
    for i in range(11):
        newcol = dfNP[:,i]
        newcolNoNaNs = newcol[~np.isnan(newcol)]
        if(i >= 7):
            if(i == 9):
                newcolNoNaNs = newcolNoNaNs.astype(bool)
            else:
                newcolNoNaNs = newcolNoNaNs.astype(int)
        if(i != 9):
            newcolNoNaNs = newcolNoNaNs[..., None]

        dictNoNans[i] = newcolNoNaNs
    graph = dict(x = np.column_stack((dictNoNans[0],dictNoNans[1],dictNoNans[2])),
    edge_attr = np.column_stack((dictNoNans[3],dictNoNans[4],dictNoNans[5],dictNoNans[6])),
    edge_index = np.column_stack((dictNoNans[7],dictNoNans[8])), y = dictNoNans[9], pid = dictNoNans[10])

    output_filename = os.path.join(output_dir, '%s%03i' % (base_prefix, graph_num))
    #print(output_filename)
    np.savez(output_filename, **graph)

newGraph = dict(np.load('../graphs/graphsNPZ/testGraph000.npz').items())
