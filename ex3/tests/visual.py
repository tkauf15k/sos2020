# %%

import pandas as pd
import numpy as np
import gzip
import sklearn.metrics
import pandas as pd
import minisom as som
from sklearn import datasets, preprocessing
import matplotlib.pyplot as plt
import seaborn as sbs
from matplotlib.collections import LineCollection


class SOMToolBox_Parse:

    def __init__(self, filename):
        self.filename = filename

    def read_weight_file(self, ):
        df = pd.DataFrame()
        if self.filename[-3:len(self.filename)] == '.gz':
            with gzip.open(self.filename, 'rb') as file:
                df, vec_dim, xdim, ydim = self._read_vector_file_to_df(df, file)
        else:
            with open(self.filename, 'rb') as file:
                df, vec_dim, xdim, ydim = self._read_vector_file_to_df(df, file)

        file.close()
        return df.astype('float64'), vec_dim, xdim, ydim

    def _read_vector_file_to_df(self, df, file):
        xdim, ydim, vec_dim, position = 0, 0, 0, 0
        for byte in file:
            line = byte.decode('UTF-8')
            if line.startswith('$'):
                xdim, ydim, vec_dim = self._parse_vector_file_metadata(line,
                                                                       xdim,
                                                                       ydim,
                                                                       vec_dim)
                if xdim > 0 and ydim > 0 and len(df.columns) == 0:
                    df = pd.DataFrame(index=range(0, ydim * xdim),
                                      columns=range(0, vec_dim))
            else:
                if len(df.columns) == 0 or vec_dim == 0:
                    raise ValueError(
                        'Weight file has no correct Dimensional information.')
                position = self._parse_weight_file_data(line, position, vec_dim,
                                                        df)
        return df, vec_dim, xdim, ydim

    def _parse_weight_file_data(self, line, position, vec_dim, df):
        splitted = line.split(' ')
        try:
            df.values[position] = list(
                np.array(splitted[0:vec_dim]).astype(float))
            position += 1
        except:
            raise ValueError(
                'The input-vector file does not match its unit-dimension.')
        return position

    def _parse_vector_file_metadata(self, line, xdim, ydim, vec_dim):
        splitted = line.split(' ')
        if splitted[0] == '$XDIM':
            xdim = int(splitted[1])
        elif splitted[0] == '$YDIM':
            ydim = int(splitted[1])
        elif splitted[0] == '$VEC_DIM':
            vec_dim = int(splitted[1])
        return xdim, ydim, vec_dim

    # %%


import numpy as np
from scipy.spatial import distance_matrix, distance
from ipywidgets import Layout, HBox, Box, widgets, interact
import plotly.graph_objects as go
from sklearn.neighbors import KDTree
from scipy.ndimage.filters import gaussian_filter

from IPython.core.display import display, HTML

display(HTML("<style>.container { width:100% !important; }</style>"))


class SomViz:

    def __init__(self, weights=[], m=None, n=None):
        self.weights = weights
        self.m = m
        self.n = n

    def umatrix(self, som_map=None, color="Viridis", interp="best", title=""):
        um = np.zeros((self.m * self.n, 1))
        neuron_locs = list()
        for i in range(self.m):
            for j in range(self.n):
                neuron_locs.append(np.array([i, j]))
        neuron_distmat = distance_matrix(neuron_locs, neuron_locs)

        for i in range(self.m * self.n):
            neighbor_idxs = neuron_distmat[i] <= 1
            neighbor_weights = self.weights[neighbor_idxs]
            um[i] = distance_matrix(np.expand_dims(self.weights[i], 0),
                                    neighbor_weights).mean()

        if som_map == None:
            return self.plot(um.reshape(self.m, self.n), color=color,
                             interp=interp, title=title)
        else:
            som_map.data[0].z = um.reshape(self.m, self.n)

    def hithist(self, som_map=None, idata=[], color='RdBu', interp="best",
                title=""):
        hist = [0] * self.n * self.m
        for v in idata:
            position = np.argmin(
                np.sqrt(np.sum(np.power(self.weights - v, 2), axis=1)))
            hist[position] += 1

        if som_map == None:
            return self.plot(np.array(hist).reshape(self.m, self.n),
                             color=color, interp=interp, title=title)
        else:
            som_map.data[0].z = np.array(hist).reshape(self.m, self.n)

    def component_plane(self, som_map=None, component=0, color="Viridis",
                        interp="best", title=""):
        if som_map == None:
            return self.plot(self.weights[:, component].reshape(-1, self.n),
                             color=color, interp=interp, title=title)
        else:
            som_map.data[0].z = self.weights[:, component].reshape(-1, n)

    def sdh(self, som_map=None, idata=[], sdh_type=1, factor=1, draw=True,
            color="Cividis", interp="best", title=""):

        import heapq
        sdh_m = [0] * self.m * self.n

        cs = 0
        for i in range(0, factor): cs += factor - i

        for vector in idata:
            dist = np.sqrt(np.sum(np.power(self.weights - vector, 2), axis=1))
            c = heapq.nsmallest(factor, range(len(dist)), key=dist.__getitem__)
            if (sdh_type == 1):
                for j in range(0, factor):  sdh_m[c[j]] += (
                                                                   factor - j) / cs  # normalized
            if (sdh_type == 2):
                for j in range(0, factor): sdh_m[c[j]] += 1.0 / dist[
                    c[j]]  # based on distance
            if (sdh_type == 3):
                dmin = min(dist)
                for j in range(0, factor): sdh_m[c[j]] += 1.0 - (
                        dist[c[j]] - dmin) / (max(dist) - dmin)

        if som_map == None:
            return self.plot(np.array(sdh_m).reshape(-1, self.n), color=color,
                             interp=interp, title=title)
        else:
            som_map.data[0].z = np.array(sdh_m).reshape(-1, self.n)

    def project_data(self, som_m=None, idata=[], title=""):

        data_y = []
        data_x = []
        for v in idata:
            position = np.argmin(
                np.sqrt(np.sum(np.power(self.weights - v, 2), axis=1)))
            x, y = position % self.n, position // self.n
            data_x.extend([x])
            data_y.extend([y])

        if som_m != None: som_m.add_trace(
            go.Scatter(x=data_x, y=data_y, mode="markers",
                       marker_color='rgba(255, 255, 255, 0.8)', ))

    def time_series(self, som_m=None, idata=[], wsize=50,
                    title=""):  # not tested

        data_y = []
        data_x = [i for i in range(0, len(idata))]

        data_x2 = []
        data_y2 = []

        qmin = np.Inf
        qmax = 0

        step = 1

        ps = []
        for v in idata:
            matrix = np.sqrt(np.sum(np.power(self.weights - v, 2), axis=1))
            position = np.argmin(matrix)
            qerror = matrix[position]
            if qmin > qerror: qmin = qerror
            if qmax < qerror: qmax = qerror
            ps.append((position, qerror))

        markerc = []
        for v in ps:
            data_y.extend([v[0]])
            rez = v[1] / qmax

            markerc.append('rgba(0, 0, 0, ' + str(rez) + ')')

            x, y = v[0] % self.n, v[0] // self.n
            if x == 0:
                y = np.random.uniform(low=y, high=y + .1)
            elif x == self.m - 1:
                y = np.random.uniform(low=y - .1, high=y)
            elif y == 0:
                x = np.random.uniform(low=x, high=x + .1)
            elif y == self.n - 1:
                x = np.random.uniform(low=x - .1, high=x)
            else:
                x, y = np.random.uniform(low=x - .1,
                                         high=x + .1), np.random.uniform(
                    low=y - .1, high=y + .1)

            data_x2.extend([x])
            data_y2.extend([y])

        ts_plot = go.FigureWidget(
            go.Scatter(x=[], y=[], mode="markers", marker_color=markerc,
                       marker=dict(colorscale='Viridis', showscale=True,
                                   color=np.random.randn(500))))
        ts_plot.update_xaxes(range=[0, wsize])

        ts_plot.data[0].x, ts_plot.data[0].y = data_x, data_y
        som_m.add_trace(go.Scatter(x=data_x2, y=data_y2, mode="markers", ))

        som_m.layout.height = 500
        ts_plot.layout.height = 500
        som_m.layout.width = 500
        ts_plot.layout.width = 1300

        return HBox([go.FigureWidget(som_m), go.FigureWidget(ts_plot)])

    def plot(self, matrix, color="Viridis", interp="best", title=""):
        return go.FigureWidget(
            go.Heatmap(z=matrix, zsmooth=interp, showscale=False,
                       colorscale=color),
            layout=go.Layout(width=1400, height=700, title=title,
                             title_x=0.5, ))

    # helper function for drawing from som unit (x1,y1) to (x2,y2)
    def draw_line(self, x1, y1, x2, y2, figure, color='red'):
        figure.add_scatter(x=[x1, x2], y=[y1, y2], line_color=color,
                           mode='lines', showlegend=False)

    # helper function for getting corrected (x,y) indices for weight array indexes
    def get_reshapesindex(self, position):
        return position % self.n, position // self.n

    def prepare_um_figure(self, color="Viridis", interp=False, title=""):
        # First compute U-matrix values
        um = np.zeros((self.m * self.n, 1))
        neuron_locs = list()
        for i in range(self.m):
            for j in range(self.n):
                neuron_locs.append(np.array([i, j]))
        neuron_distmat = distance_matrix(neuron_locs, neuron_locs)

        for i in range(self.m * self.n):
            neighbor_idxs = neuron_distmat[i] <= 1
            neighbor_weights = self.weights[neighbor_idxs]
            um[i] = distance_matrix(np.expand_dims(self.weights[i], 0),
                                    neighbor_weights).mean()

        fig, ax1 = plt.subplots(figsize=(15, 10))
        dat = um.reshape(self.m, self.n)

        if interp is not False:
            dat = gaussian_filter(dat, sigma=0.5)

        sbs.heatmap(data=dat, ax=ax1, cmap=color, cbar=False)

        return fig

        # Create U-matrix plot
        SCALE = 20
        layout = go.Layout(width=self.n * SCALE, height=self.m * SCALE,
                           title=title, title_x=0.5, )
        figure = go.FigureWidget(
            go.Heatmap(z=um.reshape(self.m, self.n), zsmooth=interp,
                       showscale=False, colorscale=color), layout=layout)
        return figure

    def neighbourhood_knn(self, idata, k=1, color="Viridis", interp=False,
                          title=""):
        # First compute U-matrix values
        um = np.zeros((self.m * self.n, 1))
        neuron_locs = list()
        for i in range(self.m):
            for j in range(self.n):
                neuron_locs.append(np.array([i, j]))
        neuron_distmat = distance_matrix(neuron_locs, neuron_locs)

        for i in range(self.m * self.n):
            neighbor_idxs = neuron_distmat[i] <= 1
            neighbor_weights = self.weights[neighbor_idxs]
            um[i] = distance_matrix(np.expand_dims(self.weights[i], 0),
                                    neighbor_weights).mean()

        # Create U-matrix plot
        SCALE = 20
        layout = go.Layout(width=self.n * SCALE, height=self.m * SCALE,
                           title=title, title_x=0.5, )
        figure = go.FigureWidget(
            go.Heatmap(z=um.reshape(self.m, self.n), zsmooth=interp,
                       showscale=False, colorscale=color), layout=layout)

        # Start k-NN computation

        idata = idata.to_numpy()
        # build kd-tree on input vectors
        tree = KDTree(idata)  # euclidean metric is already used here
        # use cache for best matching unit computation
        positionchache = {}

        # for each input vector do knn computation
        for ind_orig, v in enumerate(idata):
            if tuple(v) in positionchache:
                position1 = positionchache[tuple(v)]
            else:
                position1 = np.argmin(
                    np.sqrt(np.sum(np.power(self.weights - v, 2), axis=1)))

            nearest_dist, nearest_ind = tree.query([v], k=(
                    k + 1))  # k+1 because we also get the query point
            inds = nearest_ind[0][:]
            for ind in inds:
                if ind != ind_orig:
                    position2 = np.argmin(np.sqrt(
                        np.sum(np.power(self.weights - idata[ind], 2), axis=1)))
                    if tuple(idata[ind]) in positionchache:
                        position2 = positionchache[tuple(idata[ind])]
                    else:
                        position2 = np.argmin(np.sqrt(
                            np.sum(np.power(self.weights - idata[ind], 2),
                                   axis=1)))

                    if position1 != position2:
                        # different units, draw line
                        x1, y1 = self.get_reshapesindex(position1)
                        x2, y2 = self.get_reshapesindex(position2)
                        self.draw_line(x1, y1, x2, y2, figure)

        return figure

    def neighbourhood_radius(self, idata, radius=0.2, color="Viridis",
                             interp=False,
                             title="", highlight_longest_n: int = None):

        figure = self.prepare_um_figure(color, interp, title)

        num_nodes = idata.shape[0]
        feature_dim = idata.shape[1]
        input = idata.to_numpy()

        input_assigned_units = np.apply_along_axis(lambda x: np.argmin(
            np.linalg.norm(self.weights - x.reshape((1, feature_dim)), axis=1)),
                                                   1, input)

        assigned_unit_coords = np.apply_along_axis(
            lambda x: self.get_reshapesindex(x),
            axis=0, arr=input_assigned_units)

        assignment_x = assigned_unit_coords[0]
        assignment_y = assigned_unit_coords[1]

        distances = sklearn.metrics.pairwise_distances(input)

        tmp = distances < radius
        np.fill_diagonal(tmp, False)

        tmp2 = np.tril(tmp)
        tmp3 = tmp2.astype(np.int)
        index_matrix = np.array(
            [list(range(0, num_nodes)), ] * num_nodes).transpose()
        tmp4 = np.multiply(tmp3, index_matrix)

        tmp5 = np.where(tmp4 > 0, tmp4, -1)

        lines = set()

        for i in range(0, num_nodes):
            my_coords = (assignment_x[i], assignment_y[i])

            my_partners_filtered = np.where(tmp5[:, i] > -1)

            if len(my_partners_filtered[0]) == 0:
                continue

            partner_x_coords = np.vectorize(lambda x: assignment_x[x])(
                my_partners_filtered)
            partner_y_coords = np.vectorize(lambda y: assignment_y[y])(
                my_partners_filtered)

            coords = np.concatenate([partner_x_coords, partner_y_coords],
                                    axis=0).transpose()

            array_of_tuples = list(map(list, coords))

            neighbors = {tuple(val) for val in array_of_tuples}

            neighbors = {t for t in neighbors if t != my_coords}

            for n in neighbors:
                lines.add((my_coords, n))

        longest_lines = []
        if highlight_longest_n is not None:
            line_lengths = [(((x1, y1), (x2, y2)), np.linalg.norm(
                np.array((x1, y1) - np.array((x2, y2))))) for (x1, y1), (x2, y2)
                            in lines]
            longest_lines = [x[0] for x in
                             sorted(line_lengths, key=lambda x: x[1],
                                    reverse=True)[0:highlight_longest_n]]

        # for line in lines:
        #     if highlight_longest_n is not None and line in longest_lines:
        #         continue
        #     (x1, y1), (x2, y2) = line
        #     self.draw_line(x1, y1, x2, y2, figure)
        #
        # if highlight_longest_n is not None:
        #     for line in longest_lines:
        #         (x1, y1), (x2, y2) = line
        #         self.draw_line(x1, y1, x2, y2, figure, color='black')

        if highlight_longest_n is None:
            lc = LineCollection([l for l in lines], color="red", lw=2)
            plt.gca().add_collection(lc)
        else:
            lc = LineCollection([l for l in lines if l not in longest_lines],
                                color="red", lw=2)
            lc2 = LineCollection([l for l in lines if l in longest_lines],
                                 color="black", lw=3)
            plt.gca().add_collection(lc)
            plt.gca().add_collection(lc2)

        return figure


# %%

# interp: False, 'best', 'fast',
# color = 'viridis': https://plotly.com/python/builtin-colorscales/

##########################################
######## read from SOMToolBox ############
##########################################

def chainlink():
    trainedmap = SOMToolBox_Parse('data/chainlink_input.vec')
    idata, idim, idata_x, idata_y = trainedmap.read_weight_file()

    smap = SOMToolBox_Parse('data/chainlink_100x60.wgt.gz')
    smap, sdim, smap_x, smap_y = smap.read_weight_file()

    # Visualizaton
    viz_SOMToolBox = SomViz(smap.values.reshape(-1, sdim), smap_y, smap_x)
    # um = viz_SOMToolBox.neighbourhood_knn(k = 5, idata = idata, color='viridis', interp=False, title='U-matrix SOMToolBox')
    # um.show()

    um = viz_SOMToolBox.neighbourhood_radius(radius=0.1, idata=idata,
                                             color='viridis', interp=False,
                                             title='U-matrix SOMToolBox 222',
                                             highlight_longest_n=None)
    um.show()

chainlink()


def tencluster():
    trainedmap = SOMToolBox_Parse('data/10clusters_input.vec')
    idata, idim, idata_x, idata_y = trainedmap.read_weight_file()

    smap = SOMToolBox_Parse('data/10clusters_40x20.wgt.gz')
    smap, sdim, smap_x, smap_y = smap.read_weight_file()

    # Visualizaton
    viz_SOMToolBox = SomViz(smap.values.reshape(-1, sdim), smap_y, smap_x)

    # um = viz_SOMToolBox.neighbourhood_knn(k = 5, idata = idata, color='viridis', interp=False, title='U-matrix SOMToolBox')

    # todo: increase from radius 0.1 ... 3, then we can see the scarce cluster
    um = viz_SOMToolBox.neighbourhood_radius(radius=0.01, idata=idata,
                                             color='viridis', interp=True,
                                             title='U-matrix SOMToolBox 222',
                                             highlight_longest_n=10)

    um.show()


# tencluster()


# interp: False, 'best', 'fast',
# color = 'viridis': https://plotly.com/python/builtin-colorscales/


#############################
######## miniSOM ############1/0
#############################
def miniSom():
    m = 10
    n = 10

    # Pre-processing
    iris = datasets.load_iris().data
    min_max_scaler = preprocessing.MinMaxScaler()
    iris = min_max_scaler.fit_transform(iris)

    # Train
    s = som.MiniSom(m, n, iris.shape[1], sigma=0.8, learning_rate=0.7)
    s.train_random(iris, 10000, verbose=False)

    # Visualizaton
    viz_miniSOM = SomViz(s._weights.reshape(-1, 4), m, n)
    um1 = viz_miniSOM.umatrix(color='magma', interp='best',
                              title='U-matrix miniSOM')

    um1.show()

    return um1


##########################################
######## read from SOMToolBox ############
##########################################

def iris():
    trainedmap = SOMToolBox_Parse('iris.vec')
    idata, idim, idata_x, idata_y = trainedmap.read_weight_file()

    smap = SOMToolBox_Parse('iris.wgt.gz')
    smap, sdim, smap_x, smap_y = smap.read_weight_file()

    # Visualizaton
    viz_SOMToolBox = SomViz(smap.values.reshape(-1, sdim), smap_y, smap_x)
    # um2 = viz_SOMToolBox.umatrix(color='viridis', interp='fast', title='U-matrix SOMToolBox')
    # um2.show()
    #
    # um = viz_SOMToolBox.neighbourhood_knn(k=5, idata=idata, color='viridis',
    #                                       interp=False,
    #                                       title='U-matrix SOMToolBox')
    # um.show()

    um = viz_SOMToolBox.neighbourhood_radius(radius=0.3, idata=idata,
                                             color='viridis', interp=False,
                                             title='U-matrix SOMToolBox 222')
    um.show()
