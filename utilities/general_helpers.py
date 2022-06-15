import os
import numpy as np
import yaml
import os
import pandas
from collections import Counter
import matplotlib.pyplot as plt
# from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import sys
from shutil import copyfile
from matplotlib.colors import LinearSegmentedColormap
from shap.plots.colors import red_white_blue
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy


def flatten_list(l):
    return [item for sublist in l for item in sublist]


def plot_counter_list(l):
    letter_counts = Counter(l)
    df = pandas.DataFrame.from_dict(letter_counts, orient='index')
    df.plot(kind='bar')


def search_in_list(count_list, key):
    d = {v[0]: v[1] for v in count_list}
    return d.get(key, 0)


def is_there_overlap_in_lists(l1 , l2):
    return len([f for f in l1 if f in l2])!=0


def are_the_lists_identical(l1, l2):
    """
    NOT IDENTICAL IN ORDER BUT IDENTICAL IN CONTENT.
    :param l1:
    :param l2:
    :return:
    """
    l1 = sorted(l1)
    l2 = sorted(l2)
    return sum([l1[ii] != l2[ii] for ii in range(len(l1))]) == 0


def intersection_of_lists(l1, l2):
    """
    Returns unique appearance of each value of both lists.
    :param l1: list
    :param l2: list
    """
    s1 = list(set(l1))
    s2 = list(set(l2))
    inter_l = [v for v in s1 if v in s2]
    return inter_l


def avg(l):
    return sum(l)/len(l)


class Experiments_manager:
    def __init__(self, experiment_name, experiment_folder, config_path=None):
        # Builds experiment folder
        exp_folder_path = os.path.join(experiment_folder, experiment_name)
        if not os.path.isdir(exp_folder_path):
            os.mkdir(exp_folder_path)
        # copy config
        if config_path:
            copyfile(config_path, os.path.join(exp_folder_path, 'config.yaml'))
        self.print_file_path = os.path.join(exp_folder_path, 'prints.txt')
        self.out_file = open(self.print_file_path, 'w')
        self.orig_stdout = None

    def activate_prints_to_file(self):
        self.orig_stdout = sys.stdout
        sys.stdout = self.out_file
        return self

    def finish_run(self):
        sys.stdout = self.orig_stdout
        self.out_file.close()
        with open(self.print_file_path, 'r+') as f:
            lines = f.readlines()
            for line in lines:
                print(line, end='')

    def print(self, txt, end='\n'):
        print(txt, end=end)
        self.out_file.write(txt+end)


def experiment_manager(experiment_name, experiment_folder, experiment_config=None):
    def experiment_manager_wrapper(func):
        def inner(*args, **kwargs):
            print(f"Experiment \'{experiment_name}\' has started, prints will be saved in \'{experiment_folder}\'")
            em = Experiments_manager(experiment_name, experiment_folder, experiment_config).activate_prints_to_file()
            output = func(*args, **kwargs)
            em.finish_run()
            return output
        return inner
    return experiment_manager_wrapper


def load_yml(yml_path):
    with open(yml_path, 'r') as f:
        config = yaml.safe_load(f)
    return config['EXPERIMENT']['experiment_name'], config['EXPERIMENT']['experiments_folder'], config


def binary_search(lst, target):
    def _binary_search_loop(lst, i, j, target):
        avg = int((i + j) / 2)
        if j <= i:
            return -1
        if lst[avg] == target:
            return avg
        elif lst[avg] > target:
            return _binary_search_loop(lst, i, avg, target)
        else:
            return _binary_search_loop(lst, avg + 1, j, target)
    i = 0
    j = len(lst)
    return _binary_search_loop(lst, i, j, target)


def visualization_confusion_matrix(labels, predictions, title=None, save_path=None, display_labels=None):
    from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
    import matplotlib

    cm = confusion_matrix(labels, predictions)
    if display_labels:
        disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                                      display_labels=display_labels)
    else:
        disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                                      display_labels=['non-response', 'response'])
    _, ax = plt.subplots()
    if title:
        ax.set_title(title)
    disp.plot(include_values=True,
              cmap='viridis', ax=ax, xticks_rotation='horizontal',
              values_format=None)
    if save_path:
        plt.savefig(os.path.join(save_path+".png"))
    plt.show()

    return cm


def binary_list_to_indices_list(lst):
    new_lst = [ii for ii in range(len(lst)) if lst(ii)]
    return new_lst


def indices_list_to_binary_list(lst):
    new_lst = [ii in lst for ii in range(len(lst))]
    return new_lst


def create_folder(folder):
    if not os.path.isdir(folder):
        os.mkdir(folder)
        print(f'folder: {folder} has been created')


def get_common_indices_of_boolean_lists(*argv):
    return [sum([not argv[jj][ii] for jj in range(len(argv))]) == 0 for ii in range(len(argv[0]))]


def flip_sign_of_boolean_list(l1):
    return [not aa for aa in l1]


def sort_dic(dic, by_key=True, descending=False):
    if by_key:
        return {k: v for k, v in sorted(dic.items(), key=lambda item: item[0], reverse=descending)}
    return {k: v for k, v in sorted(dic.items(), key=lambda item: item[1], reverse=descending)}


def transpose(l1):
    numpy_array = np.array(l1)#.astype(np.float32)
    transposed_list = numpy_array.T
    return transposed_list


def annotate_boxplot(bpdict, data, annotate_params=None,
                     x_offset=0.05, x_loc=0,
                     text_offset_x=35,
                     text_offset_y=20):
    """Annotates a matplotlib boxplot with labels marking various centile levels.

    Parameters:
    - bpdict: The dict returned from the matplotlib `boxplot` function. If you're using pandas you can
    get this dict by setting `return_type='dict'` when calling `df.boxplot()`.
    - annotate_params: Extra parameters for the plt.annotate function. The default setting uses standard arrows
    and offsets the text based on other parameters passed to the function
    - x_offset: The offset from the centre of the boxplot to place the heads of the arrows, in x axis
    units (normally just 0-n for n boxplots). Values between around -0.15 and 0.15 seem to work well
    - x_loc: The x axis location of the boxplot to annotate. Usually just the number of the boxplot, counting
    from the left and starting at zero.
    text_offset_x: The x offset from the arrow head location to place the associated text, in 'figure points' units
    text_offset_y: The y offset from the arrow head location to place the associated text, in 'figure points' units

    Example of use:
    annotate_boxplot(pd.DataFrame(arr).boxplot(whis=[5, 95], return_type='dict'), arr)
    """

    data.sort()

    if annotate_params is None:
        annotate_params = dict(xytext=(text_offset_x, text_offset_y), textcoords='offset points',
                               arrowprops={'arrowstyle': '->'})
    v_95 = round(data[int(len(data) * 0.95) - 1], 2)
    v_05 = round(data[int(len(data) * 0.05) - 1], 2)
    v_50 = round(data[int(len(data) * 0.5) - 1], 2)
    v_25 = round(data[int(len(data) * 0.25) - 1], 2)
    v_75 = round(data[int(len(data) * 0.75) - 1], 2)
    plt.annotate(f'{v_50} (median)', (x_loc + 1 + x_offset, bpdict['medians'][x_loc].get_ydata()[0]), **annotate_params)
    plt.annotate(f'{v_25} (25%)', (x_loc + 1 + x_offset, bpdict['boxes'][x_loc].get_ydata()[0]), **annotate_params)
    plt.annotate(f'{v_75} (75%)', (x_loc + 1 + x_offset, bpdict['boxes'][x_loc].get_ydata()[2]), **annotate_params)
    plt.annotate(f'{v_05} (5%)', (x_loc + 1 + x_offset, bpdict['caps'][x_loc * 2].get_ydata()[0]), **annotate_params)
    plt.annotate(f'{v_95} (95%)', (x_loc + 1 + x_offset, bpdict['caps'][(x_loc * 2) + 1].get_ydata()[0]),
                 **annotate_params)


def transpose_list(l):
    return list(map(list, zip(*l)))


def bold(string):
    """
    Return Bold string for print
    :param string:
    :return:
    """
    return '\033[1m' + string + '\033[0m'


def plot_stackedbar_p(df, labels, title, subtitle, file_path=None, colors=[]):
    if not len(colors):
        # colors = ['darkred', 'yellow', 'tomato', 'orange', 'bisque', 'aqua', 'powderblue', 'violet', 'purple',
        #   '#1D2F6F', '#8390FA', '#6EAF46', '#FAC748', '#1D2F6F', '#8390FA', '#6EAF46', '#FAC748']
        colormap = plt.cm.gist_ncar  # nipy_spectral, Set1,Paired           [1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12]
        colors = np.array([colormap(i) for i in np.linspace(0, 0.9, 100)])[
            [84, 99, 15, 70, 60, 53, 92, 23, 30, 1, 40, 80, 40, 12, 13, 14]]
    fields = df.columns.tolist()

    # figure and axis
    fig, ax = plt.subplots(1, figsize=(12, 10))
    # plot bars
    left = len(df) * [0]
    for idx, name in enumerate(fields):
        plt.barh(df.index, df[name], left=left, color=colors[idx])
        left = left + df[name]
        # title and subtitle
        plt.title(title, loc='left')
        plt.text(0, ax.get_yticks()[-1] + 0.75, subtitle)
        # legend
        # plt.legend(labels, bbox_to_anchor=([0.58, 1, 0, 0]), ncol=4, frameon=False)
        plt.legend(labels, bbox_to_anchor=([-0.1, 1, 0, 0]), ncol=1, frameon=False)

        # remove spines
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        # format x ticks
        xticks = np.arange(0, 1.1, 0.1)
        xlabels = ['{}%'.format(i) for i in np.arange(0, 101, 10)]
        plt.xticks(xticks, xlabels)
        # adjust limits and draw grid lines
        plt.ylim(-0.5, ax.get_yticks()[-1] + 0.5)
        ax.xaxis.grid(color='gray', linestyle='dashed')
    plt.show()
    if file_path:
        fig.savefig(file_path)

def annotate_seaborn_barplot(axs):
    for p in axs.patches:
        axs.annotate(str(np.round(p.get_height(),3)), (p.get_x() * 1.005, p.get_height() * 1.005))


def show_values_on_bars(axs):
    """
    Seaborn
    :param axs:
    :return:
    """
    def _show_on_single_plot(ax):
        for p in ax.patches:
            _x = p.get_x() + p.get_width() / 2
            _y = p.get_y() + p.get_height()
            value = '{:.2f}%'.format(p.get_height())
            ax.text(_x, _y, value, ha="center")

    if isinstance(axs, np.ndarray):
        for idx, ax in np.ndenumerate(axs):
            _show_on_single_plot(ax)
    else:
        _show_on_single_plot(axs)


def get_hierarchical_cells_order(values):
    linked = linkage(values, 'single')
    order = hierarchy.leaves_list(hierarchy.optimal_leaf_ordering(linked, values))
    return order


def shapely_heatmap(shap_values, feature_names, instance_order, cell_values,
                    max_display=10, cmap=red_white_blue, show=True):
    """ Create a heatmap plot of a set of SHAP values.

    This plot is designed to show the population substructure of a dataset using supervised
    clustering and a heatmap. Supervised clustering involves clustering data points not by their original
    feature values but by their explanations. By default we cluster using shap.utils.hclust_ordering
    but any clustering can be used to order the samples.

    Parameters
    ----------
    shap_values : shap.Explanation
        A multi-row Explanation object that we want to visualize in a cluster ordering.

    instance_order : OpChain or numpy.ndarray
        A function that returns a sort ordering given a matrix of SHAP values and an axis, or
        a direct sample ordering given as an numpy.ndarray.

    feature_values : OpChain or numpy.ndarray
        A function that returns a global summary value for each input feature, or an array of such values.

    feature_order : None, OpChain, or numpy.ndarray
        A function that returns a sort ordering given a matrix of SHAP values and an axis, or
        a direct input feature ordering given as an numpy.ndarray. If None then we use
        feature_values.argsort

    max_display : int
        The maximum number of features to display.

    show : bool
        If show is set to False then we don't call the matplotlib.pyplot.show() function. This allows
        further customization of the plot by the caller after the bar() function is finished.

    """
    feature_values = np.abs(cell_values).mean(0)
    feature_order = np.argsort(np.abs(shap_values).sum(axis=0))[::-1]
    xlabel = "Instances"
    values = shap_values[:, feature_order]
    feature_values = feature_values[feature_order]
    feature_names = feature_names[feature_order]


    if not instance_order is None:
        values = values[instance_order]

    # collapse
    if values.shape[1] > max_display:
        new_values = np.zeros((values.shape[0], max_display))
        new_values[:, :max_display - 1] = values[:, :max_display - 1]
        new_values[:, max_display - 1] = values[:, max_display - 1:].sum(1)
        new_feature_values = np.zeros(max_display)
        new_feature_values[:max_display - 1] = feature_values[:max_display - 1]
        new_feature_values[max_display - 1] = feature_values[max_display - 1:].sum()
        feature_names = list(feature_names[:max_display])
        feature_names[-1] = "Sum of %d other features" % (values.shape[1] - max_display + 1)
        values = new_values
        feature_values = new_feature_values

    # define the plot size
    row_height = 0.5
    plt.gcf().set_size_inches(8, values.shape[1] * row_height + 2.5)

    # plot the matrix of SHAP values as a heat map
    vmin = np.nanpercentile(values.flatten(), 1)
    vmax = np.nanpercentile(values.flatten(), 99)
    plt.imshow(
        values.T, aspect=0.7 * values.shape[0] / values.shape[1], interpolation="nearest", vmin=min(vmin, -vmax),
        vmax=max(-vmin, vmax),
        cmap=cmap
    )
    yticks_pos = np.arange(values.shape[1])
    yticks_labels = feature_names

    plt.yticks([-1.5] + list(yticks_pos), ["f(x)"] + list(yticks_labels), fontsize=13)

    plt.ylim(values.shape[1] - 0.5, -3)

    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')
    plt.gca().spines['right'].set_visible(True)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.axhline(-1.5, color="#aaaaaa", linestyle="--", linewidth=0.5)
    fx = values.T.mean(0)
    plt.plot(-fx / np.abs(fx).max() - 1.5, color="#000000", linewidth=1)
    # pl.colorbar()
    plt.gca().spines['left'].set_bounds(values.shape[1] - 0.5, -0.5)
    plt.gca().spines['right'].set_bounds(values.shape[1] - 0.5, -0.5)
    b = plt.barh(
        yticks_pos, (feature_values / np.abs(feature_values).max()) * values.shape[0] / 20,
        0.7, align='center', color="#000000", left=values.shape[0] * 1.0 - 0.5
        # color=[colors.red_rgb if shap_values[feature_inds[i]] > 0 else colors.blue_rgb for i in range(len(y_pos))]
    )
    for v in b:
        v.set_clip_on(False)
    plt.xlim(-0.5, values.shape[0] - 0.5)
    plt.xlabel(xlabel)

    if True:
        import matplotlib.cm as cm
        m = cm.ScalarMappable(cmap=cmap)
        m.set_array([min(vmin, -vmax), max(-vmin, vmax)])
        cb = plt.colorbar(m, ticks=[min(vmin, -vmax), max(-vmin, vmax)], aspect=1000, fraction=0.0090, pad=0.10,
                          panchor=(0, 0.05))
        # cb.set_ticklabels([min(vmin,-vmax), max(-vmin,vmax)])
        cb.set_label("SHAP value", size=12, labelpad=-10)
        cb.ax.tick_params(labelsize=11, length=0)
        cb.set_alpha(1)
        cb.outline.set_visible(False)
        bbox = cb.ax.get_window_extent().transformed(plt.gcf().dpi_scale_trans.inverted())
        cb.ax.set_aspect((bbox.height - 0.9) * 15)
        cb.ax.set_anchor((1, 0.2))
        # cb.draw_all()

    for i in [0]:
        plt.gca().get_yticklines()[i].set_visible(False)

    if show:
        plt.show()