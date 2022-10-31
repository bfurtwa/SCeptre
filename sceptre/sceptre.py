from typing import Tuple, Sequence, Mapping, Optional, Union
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

from anndata import AnnData
from copy import copy
import pandas as pd
import ntpath
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import logging
import seaborn as sns
import warnings

figwd = 7.2  # standard figure width
cellsize = 20  # size to plot cells
wspace = 1  # space between scanpy plots to make room for legends
hspace = 0.5  # space between scanpy plots to make room for legends


def create_meta_data(input_dir: str, output_dir: str):
    """Create a meta data table from PD output and mapping tables.
    Requires the following tables in the `input_dir`:

    * PD Protein
    * PD InputFiles
    * file_sample_mapping
    * plate_layout_mapping
    * sort_layout
    * sample_layout
    * facs_data

    See the example files to understand the structure of these tables.
    Alternatively, create meta data on you own and use :func:`~sceptre.load_dataset`

    Parameters
    ----------
    input_dir
        The path to the input directory.
    output_dir
        The path to the output directory.
    Returns
    -------
    :obj:`None`

    Saves the meta table in `output_dir`.
    """
    import os

    # PD tables
    for file in os.listdir(input_dir):
        if "_Proteins.txt" in file:
            prot = pd.read_table(input_dir + file, low_memory=False)
        if "_InputFiles.txt" in file:
            files = pd.read_table(input_dir + file)
            files["File Name"] = files["File Name"].apply(lambda x: ntpath.basename(x))

    # mapping tables
    file_sample_mapping = pd.read_table("{}file_sample_mapping.txt".format(input_dir))
    plate_layout_mapping = pd.read_table(
        "{}plate_layout_mapping.txt".format(input_dir)
    ).set_index("Plate")

    # plate data tables
    plate_data = {k: {} for k in plate_layout_mapping.index.unique()}
    for plate in plate_data.keys():
        plate_data[plate]["sort_layout"] = pd.read_table(
            "{}{}".format(input_dir, plate_layout_mapping.loc[plate, "Sort Layout"]),
            index_col=0,
        )
        plate_data[plate]["label_layout"] = pd.read_table(
            "{}{}".format(input_dir, plate_layout_mapping.loc[plate, "Label Layout"]),
            index_col=0,
        ).fillna("")
        plate_data[plate]["sample_layout"] = pd.read_table(
            "{}{}".format(input_dir, plate_layout_mapping.loc[plate, "Sample Layout"]),
            index_col=0,
        )
        plate_data[plate]["facs_data"] = pd.read_table(
            "{}{}".format(input_dir, plate_layout_mapping.loc[plate, "Facs Data"])
        )
        plate_data[plate]["facs_data"] = plate_data[plate]["facs_data"].drop(
            ["Row", "Column"], axis=1
        )
        plate_data[plate]["facs_data"] = plate_data[plate]["facs_data"].set_index(
            "Well"
        )

    # create cell metadata
    # add each channel from each file to the rows
    meta = pd.DataFrame(
        [
            x.split(" ")[2:]
            for x in prot.columns[prot.columns.str.contains("Abundance")]
        ],
        columns=["File ID", "Channel"],
    )
    # add the file name
    meta = meta.merge(
        files.set_index("File ID")["File Name"],
        left_on="File ID",
        right_index=True,
        validate="many_to_one",
    )
    # add the plate and sample
    _ = len(meta)
    meta = meta.merge(file_sample_mapping, on="File Name", validate="many_to_one")
    if len(meta) < _:
        raise ValueError("Error in file_sample_mapping.txt")
    # add the well information via the sample_layout and label_layout for each plate
    for i in meta.index:
        p, s, c = meta.loc[i, ["Plate", "Sample", "Channel"]]
        p_d = plate_data[p]
        well = (
            p_d["sample_layout"][
                (p_d["sample_layout"] == s) & (p_d["label_layout"] == c)
            ]
            .stack()
            .index.tolist()
        )
        if len(well) > 1:
            raise ValueError(
                "Error in plate layout data: Plate {}, Sample {}, Channel {}".format(
                    p, s, c
                )
            )
        elif len(well) == 0:
            row, col, well = pd.NA, pd.NA, pd.NA
        else:
            row = well[0][0]
            col = well[0][1]
            well = "".join(well[0])

        meta.loc[i, ["Row", "Column", "Well"]] = row, col, well

        # use the sort layout to map the sorted population and add the facs data
        if not pd.isna(well):
            meta.loc[i, "Sorted Population"] = plate_data[p]["sort_layout"].loc[
                row, col
            ]
            # add the facs data
            # meta.loc[i] = meta.loc[i].append(plate_data[p]['facs_data'].loc[well, :])
        else:
            meta.loc[i, "Sorted Population"] = pd.NA

    # add the facs data for each plate
    _ = []
    for p in meta["Plate"].unique():
        _.append(
            meta.loc[meta["Plate"] == p].merge(
                plate_data[p]["facs_data"],
                left_on="Well",
                right_index=True,
                how="left",
            )
        )
    meta = pd.concat(_)
    meta = meta.rename(columns={"Population": "Gated Population"})
    meta.to_csv(output_dir + "meta.txt", sep="\t", index=False)


def load_dataset(proteins: str, psms: str, msms: str, files: str, meta: str):
    """Load the dataset from specified paths.

    Parameters
    ----------
    proteins
        The path to the PD protein table.
    psms
        The path to the PD PSMs table.
    msms
        The path to the PD MSMS table.
    files
        The path to the PD InputFiles table.
    meta
        The path to the Meta table.
    Returns
    -------
    A dict containing all required tables.
    """
    prot = pd.read_table(proteins, low_memory=False)
    # To use the Gene Symbol as index:
    # Set nan Gene Symbol to protein accession
    # and if Gene Symbol is not unique, add the protein accession to make duplicates unique.
    nans = prot["Gene Symbol"].isna()
    for i in nans.index:
        if nans[i]:
            prot.loc[i, "Gene Symbol"] = prot.loc[i, "Accession"]
    duplicates = prot["Gene Symbol"].duplicated(keep=False)
    for i in duplicates.index:
        if duplicates[i]:
            prot.loc[i, "Gene Symbol"] = (
                prot.loc[i, "Gene Symbol"] + "_" + prot.loc[i, "Accession"]
            )

    psms = pd.read_table(psms, low_memory=False)
    msms = pd.read_table(msms, low_memory=False)
    files = pd.read_table(files, low_memory=False)
    files["File Name"] = files["File Name"].apply(lambda x: ntpath.basename(x))
    meta = pd.read_table(meta, low_memory=False)

    return {"proteins": prot, "psms": psms, "msms": msms, "files": files, "meta": meta}


def plot_psms_msms(dataset: Mapping, figsize: Tuple[float, float] = (figwd, 3.5)):
    """Plot the number of MSMS spectra and PSMs per file.

    Parameters
    ----------
    dataset
        Dict containing PD tables and meta information.
    figsize
        Size of the plotted figure.
    Returns
    -------
    :obj:`None`
    """
    ms_perf = pd.concat(
        [
            dataset["msms"].groupby("File ID").size(),
            dataset["psms"].groupby("File ID").size(),
        ],
        axis=1,
    )
    ms_perf.columns = ["MSMS", "PSMs"]
    ordered_file_id = sorted(
        ms_perf.index.tolist(), key=lambda x: int(x.split("F")[-1])
    )

    fig, ax = plt.subplots(figsize=figsize)
    ms_perf.loc[ordered_file_id, :].plot.barh(
        title="Number of MSMS spectra and PSMs per File", grid=True, ax=ax
    )
    # legend:
    # shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    # put a legend below current axis
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=2)
    fig.tight_layout()


def plot_avg_sn(dataset: Mapping, figsize: Tuple[float, float] = (figwd, 3)):
    """Boxplot of the average S/N of each PSM for each file.

    Parameters
    ----------
    dataset
        Dict containing PD tables and meta information.
    figsize
        Size of the plotted figure.
    Returns
    -------
    :obj:`None`
    """
    ms_perf = pd.concat(
        [
            dataset["msms"].groupby("File ID").size(),
            dataset["psms"].groupby("File ID").size(),
        ],
        axis=1,
    )
    ms_perf.columns = ["MSMS", "PSMs"]
    ordered_file_id = sorted(
        ms_perf.index.tolist(), key=lambda x: int(x.split("F")[-1])
    )
    # set the File ID column to a ordered categorical
    cat_dtype = pd.api.types.CategoricalDtype(categories=ordered_file_id, ordered=True)
    dataset["psms"]["File ID"] = dataset["psms"]["File ID"].astype(cat_dtype)

    fig, ax = plt.subplots(figsize=figsize)
    dataset["psms"].boxplot(
        "Average Reporter SN", by="File ID", vert=False, grid=True, ax=ax
    )
    ax.xaxis.label.set_visible(False)
    ax.set_ylabel("File ID")
    plt.title("Average Reporter S/N of PSMs per File")
    plt.suptitle("")
    fig.tight_layout()


def plot_set_overview(dataset: Mapping, figsize: Tuple[float, float] = (figwd, 10)):
    """Barplot of the median log10(S/N) of each quantification channel in each file.

    Parameters
    ----------
    dataset
        Dict containing PD tables and meta information.
    figsize
        Size of the plotted figure.
    Returns
    -------
    :obj:`None`
    """
    grouped = dataset["psms"].groupby("File ID")
    group_data = (
        grouped[
            dataset["psms"].columns[dataset["psms"].columns.str.contains("Abundance")]
        ]
        .median()
        .apply(np.log10)
    )
    group_data.columns = [x.split(" ")[1] for x in group_data.columns]
    rowlength = int(np.ceil(grouped.ngroups / 3))
    fig, axs = plt.subplots(
        nrows=rowlength,
        ncols=3,
        gridspec_kw=dict(hspace=0.35, wspace=0.1),
        sharex=True,
        sharey=True,
        figsize=figsize,
    )

    targets = zip(grouped.groups.keys(), axs.flatten())
    for i, (key, ax) in enumerate(targets):
        group_data.loc[key, :].plot.bar(ax=ax, grid=True, title=key)

    # remove empty axes
    empty_axis = len(grouped.groups.keys()) - len(axs.flatten())
    for i in range(empty_axis, 0):
        axs.flat[i].set_visible(False)


def print_ms_stats(dataset: Mapping, s_c_channels: Sequence[str]):
    """Print statistics about the dataset.

    Parameters
    ----------
    dataset
        Dict containing PD tables and meta information.
    s_c_channels
        The channels containing single-cells.
    Returns
    -------
    :obj:`None`
    """

    print("Protein IDs: {}".format(len(dataset["proteins"])))
    print("Peptide IDs: {}".format(len(dataset["psms"]["Annotated Sequence"].unique())))
    print("PSMs: {}".format(len(dataset["psms"])))
    print("PSM rate: {}".format(round(len(dataset["psms"]) / len(dataset["msms"]), 3)))
    print(
        "Median of median S/N in single-cell channels: {}".format(
            round(
                dataset["psms"][
                    dataset["psms"].columns[
                        dataset["psms"].columns.str.contains("Abundance")
                        & dataset["psms"].columns.str.contains("|".join(s_c_channels))
                    ]
                ]
                .median()
                .median(),
                3,
            )
        )
    )
    print(
        "Median of mean S/N in single-cell channels: {}".format(
            round(
                dataset["psms"][
                    dataset["psms"].columns[
                        dataset["psms"].columns.str.contains("Abundance")
                        & dataset["psms"].columns.str.contains("|".join(s_c_channels))
                    ]
                ]
                .mean()
                .median(),
                3,
            )
        )
    )
    print(
        "Median S/N of booster channel: {}".format(
            round(dataset["psms"]["Abundance 126"].median(), 3)
        )
    )
    # counting number of protein IDs per file depends on version of Proteome Discoverer
    # in version 2 missing proteins are "Not Found", in version 3 NaN
    nums = []
    for f in dataset["psms"]["File ID"].unique():
        df = dataset["proteins"][
            dataset["proteins"].columns[
                dataset["proteins"].columns.str.contains("Found in Sample")
                & dataset["proteins"].columns.str.contains(f)
            ]
        ]
        nums.append(((df != "Not Found") & (~pd.isna(df))).apply(any, axis=1).sum())
    print("Mean protein IDs per file: {}".format(round(np.mean(nums), 3)))


def plot_interference(dataset: Mapping):
    """Violin plot of the isolation interference of all PSMs.

    Parameters
    ----------
    dataset
        Dict containing PD tables and meta information.
    Returns
    -------
    :obj:`None`
    """
    import matplotlib.ticker as ticker

    fig, ax = plt.subplots(figsize=(2, 2.5))
    sns.violinplot(
        y="Isolation Interference in Percent",
        data=dataset["psms"],
        inner="quartile",
        bw=0.1,
        cut=0,
        ax=ax,
    )
    ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
    # linestyles = [':', '-', ':']
    plt.title("Interference of all PSMs")
    plt.grid()
    fig.tight_layout()


def dataset_to_scanpy(dataset: Mapping, temp_dir: str = "../results/tmp"):
    """Load the dataset dict into a scanpy AnnData object.

    Parameters
    ----------
    dataset
        Dict containing PD tables and meta information.
    temp_dir
        A directory to save temporary data in.
    Returns
    -------
    :class:`~anndata.AnnData`
    """
    quant = dataset["proteins"].set_index("Gene Symbol").copy()
    if not quant.index.is_unique:
        raise IndexError("Protein index not unique")
    quant = quant[
        quant.columns[quant.columns.str.contains("Abundance")]
    ]  # only quantification columns
    file_id = [x.split(" ")[-2] for x in quant.columns]
    channel = [x.split(" ")[-1] for x in quant.columns]
    quant.columns = pd.MultiIndex.from_tuples(
        zip(file_id, channel), names=["File ID", "Channel"]
    )
    # quant[quant < 1.1] = pd.NA  # set S/N values below 1.1 to NA
    quant = quant.dropna(how="all").fillna(
        0
    )  # remove all NA proteins and fill remaining NA with 0
    quant_meta = (
        dataset["meta"]
        .set_index(["File ID", "Channel"])
        .loc[quant.columns, :]
        .copy()
        .reset_index()
    )

    # save to file and load it in scanpy
    quant.to_csv(
        "{}/scanpy_data.txt".format(temp_dir), sep="\t", header=False, index=True
    )
    adata = sc.read_text(
        "{}/scanpy_data.txt".format(temp_dir), delimiter="\t", first_column_names=False
    ).T

    adata.obs = quant_meta
    prots = dataset["proteins"].copy()
    prot_anno = prots[
        [
            "Accession",
            "Gene Symbol",
            "Description",
            "Biological Process",
            "Cellular Component",
            "Molecular Function",
            "KEGG Pathways",
            "Reactome Pathways",
            "WikiPathways",
        ]
    ]
    # object columns to category for .var
    prot_anno = pd.concat(
        [
            prot_anno.select_dtypes([], ["object"]),
            prot_anno.select_dtypes(["object"]).apply(
                pd.Series.astype, dtype="category"
            ),
        ],
        axis=1,
    ).reindex(prot_anno.columns, axis=1)

    adata.var = adata.var.merge(
        prot_anno,
        how="left",
        left_index=True,
        right_on="Gene Symbol",
    ).set_index("Gene Symbol")

    return adata


def normalize(
    adata: AnnData,
    method: Literal["file_channel", "file", "channel"] = "file_channel",
    iter_thresh: float = 1.1,
    na_thresh: Optional[float] = 1.1,
    drop_na: bool = True,
):
    """Normalize expression values in AnnData object.
    The median expression of each protein is equalized across files/channels
    by applying correction factors.
    If 'file_channel' is selected, the factors are applied iteratively for
    file and channel until the highest change in values compared to the
    previous iteration is below 'iter_thresh'.

    Parameters
    ----------
    adata
        The annotated data matrix.
    method : {``'file_channel'``, ``'file'``, ``'channel'``}
        The normalization method. Options are:
        * 'file_channel': Equalize medians across files and channels.
        * 'file': Equalize medians across files.
        * 'channel': Equalize medians across channels.
    iter_thresh
        Only for 'file_channel'. Stop iterating when the highest
        change in the expression matrix is below this value.
    na_thresh
        If not None, replace expression values below this value with 0.
    drop_na
        Drop proteins with no expression values (i.e. 0)
    Returns
    -------
    :obj:`None`

    Updates `adata` with the normalized expression values.
    """
    # get a df from adata
    quant = pd.DataFrame(
        adata.X.T.copy(),
        columns=adata.obs[["File ID", "Channel"]]
        .set_index(["File ID", "Channel"])
        .index,
    ).replace(0, np.nan)

    if method == "file_channel":
        for i in range(100):  # iterate to converge to normalized channel and file
            quant_0 = quant.copy()
            # file bias normalization
            # calculate median for each protein in each sample
            med = quant.T.reset_index().groupby("File ID").median().T
            # calculate the factors needed for a median shift
            med_tot = med.median(axis=1)
            factors = med.divide(med_tot, axis=0)

            quant = quant.groupby(axis=1, level=0).apply(
                lambda x: x.divide(factors.loc[:, x.name], axis=0)
            )

            # channel bias normalization
            # calculate median for each protein in each channel
            med = quant.T.reset_index().groupby("Channel").median().T
            # calculate the factors needed for a median shift
            med_tot = med.median(axis=1)
            factors = med.divide(med_tot, axis=0)

            quant = quant.groupby(axis=1, level=1).apply(
                lambda x: x.divide(factors.loc[:, x.name], axis=0)
            )

            # stop iterating when the change in quant to the previous iteration is below iter_thresh
            if (abs(quant - quant_0).max().max()) <= iter_thresh:
                break
        print("performed {} iterations".format(i + 1))

    elif method == "file":

        # file bias normalization
        # calculate median for each protein in each sample
        med = quant.T.reset_index().groupby("File ID").median().T
        # calculate the factors needed for a median shift
        med_tot = med.median(axis=1)
        factors = med.divide(med_tot, axis=0)

        quant = quant.groupby(axis=1, level=0).apply(
            lambda x: x.divide(factors.loc[:, x.name], axis=0)
        )

    elif method == "channel":

        # channel bias normalization
        # calculate median for each protein in each channel
        med = quant.T.reset_index().groupby("Channel").median().T
        # calculate the factors needed for a median shift
        med_tot = med.median(axis=1)
        factors = med.divide(med_tot, axis=0)

        quant = quant.groupby(axis=1, level=1).apply(
            lambda x: x.divide(factors.loc[:, x.name], axis=0)
        )

    # apply na_thresh and remove all NA proteins
    if na_thresh:
        print(
            "{} values below {} were set to 0".format(
                (quant < na_thresh).sum().sum(), na_thresh
            )
        )
        quant[quant < na_thresh] = pd.NA
    adata.X = quant.fillna(0).values.T
    if drop_na:
        sc.pp.filter_genes(adata, min_cells=1)


def calculate_cell_filter(
    adata: AnnData,
    min_proteins: int = 500,
    thresh_sum: float = 3,
    plot_qc: bool = True,
    scatter_labels: Tuple[str, str] = ("Channel", "Sorted Population"),
    figsizes: Tuple[Tuple[float, float], Tuple[float, float]] = (
        (figwd, 4.5),
        (figwd, 3),
    ),
):
    """Calculate the cell filtering based on defined thresholds.
    Thresholds control MAD-based outlier detection.
    Optionally shows diagnostic plots.

    Parameters
    ----------
    adata
        The annotated data matrix.
    min_proteins
        The minimum number of proteins per cell
    thresh_sum
        The threshold for the Log2 Sum S/N per cell
    plot_qc
        Whether to show quality control plots
    scatter_labels
        The labels to use for the two scatterplots.
    figsizes
        Figure sizes for both figures
    Returns
    -------
    A list of two :class:`~matplotlib.figure` objects if plot_qc==True.

    Updates `adata` with the following fields.

    Num Proteins : :class:`~pandas.Series` (``adata.obs``, dtype ``int``)
    Array of dim (number of samples) that stores the number of proteins
    identified in each cell.
    Log2 Sum S/N : :class:`~pandas.Series` (``adata.obs``, dtype ``float``)
    Array of dim (number of samples) that stores the Log2 Sum S/N of each cell.
    Pass Cell Filter : :class:`~pandas.Series` (``adata.obs``, dtype ``bool``)
    Array of dim (number of samples) that stores if the cell passed the filter.
    """

    def mad_based_outlier(points, thresh):
        if len(points.shape) == 1:
            points = points.to_numpy()[:, None]
        median = np.median(points, axis=0)
        diff = np.sum((points - median) ** 2, axis=-1)
        diff = np.sqrt(diff)
        med_abs_deviation = np.median(diff)

        modified_z_score = 0.6745 * diff / med_abs_deviation
        df = pd.DataFrame(
            {"points": points.flatten(), "modified_z_score": modified_z_score.flatten()}
        )
        outlier_up = df[
            (df["modified_z_score"] <= thresh) & (df["points"] > median[0])
        ]["points"].max()
        outlier_down = df[
            (df["modified_z_score"] <= thresh) & (df["points"] < median[0])
        ]["points"].min()

        return outlier_up, outlier_down

    # annotate cells with qc parameters
    adata.obs["Num Proteins"] = (adata.X > 0).sum(axis=1)
    adata.obs["Log2 Sum S/N"] = np.log2((adata.X).sum(axis=1))

    sum_sn_max, sum_sn_min = mad_based_outlier(
        adata.obs["Log2 Sum S/N"], thresh=thresh_sum
    )

    # annotate cells with filter pass
    adata.obs["Pass Cell Filter"] = (
        (adata.obs["Log2 Sum S/N"] >= sum_sn_min)
        & (adata.obs["Log2 Sum S/N"] <= sum_sn_max)
        & (adata.obs["Num Proteins"] >= min_proteins)
    )

    print(
        "{} of {} cells do not pass filter".format(
            (len(adata.obs) - adata.obs["Pass Cell Filter"].sum()), len(adata.obs)
        )
    )
    
    if plot_qc:
        fig_objects = []
        # plot qc params for each cell
        gridspec = dict(
            hspace=0.0, height_ratios=[1, 1, 0.5, 1, 1]
        )  # invisible axis for title between axes
        fig, axs = plt.subplots(nrows=5, ncols=1, figsize=figsizes[0], gridspec_kw=gridspec)
        labels = ["Log2 Sum S/N", "Num Proteins"]
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        for i, l in enumerate(labels):
            ax = axs[i]
            adata.obs[l].plot(ax=ax, grid=True, color=colors[i])
            ax.set_ylabel(l)
            ax.legend().remove()
            if i == 0:
                ax.set_title("Cells before filtering")
                ax.axhline(sum_sn_max, color="black", linestyle="--")
                ax.axhline(sum_sn_min, color="black", linestyle="--")
                ax.xaxis.set_ticklabels([])
            if i == 1:
                ax.axhline(min_proteins, color="black", linestyle="--")

        # make sure that x-axis is shared. Needed for files with 0 proteins
        axs[0].get_shared_x_axes().join(axs[0], axs[1])
        axs[1].autoscale()

        axs[2].set_visible(False)

        # plot after filter
        for i, l in enumerate(labels):
            ax = axs[i + 3]
            adata.obs.loc[adata.obs["Pass Cell Filter"], l].plot(
                ax=ax, grid=True, color=colors[i]
            )
            if i == 0:
                ax.set_title("Cells after filtering")
                ax.xaxis.set_ticklabels([])
            ax.set_ylabel(labels[i])
            ax.legend().remove()
        ax.set_xlabel("Cell index")
        fig.tight_layout()
        fig_objects.append(fig)

        # scatterplots
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=figsizes[1])
        sc.pl.scatter(
            adata,
            x="Log2 Sum S/N",
            y="Num Proteins",
            color=scatter_labels[0],
            size=cellsize,
            show=False,
            title="Cell filter by {}".format(scatter_labels[0]),
            ax=axs[0],
        )
        axs[0].axvline(sum_sn_max, color="black", linestyle="--")
        axs[0].axvline(sum_sn_min, color="black", linestyle="--")
        axs[0].axhline(min_proteins, color="black", linestyle="--")

        sc.pl.scatter(
            adata,
            x="Log2 Sum S/N",
            y="Num Proteins",
            color=scatter_labels[1],
            size=cellsize,
            show=False,
            title="Cell filter by {}".format(scatter_labels[1]),
            ax=axs[1],
        )
        axs[1].axvline(sum_sn_max, color="black", linestyle="--")
        axs[1].axvline(sum_sn_min, color="black", linestyle="--")
        axs[1].axhline(min_proteins, color="black", linestyle="--")
        fig.tight_layout()

        fig_objects.append(fig)

        return fig_objects


def apply_cell_filter(adata: AnnData, min_cells: int = 3):
    """Remove cells that do not pass cell filter from `adata`.
    Also removes Proteins that are detected in less than `min_cells`.
    Requires having ran :func:`~sceptre.calculate_cell_filter` first.

    Parameters
    ----------
    adata
        The annotated data matrix.
    min_cells
        Minimum number of cells a protein needs to be detected to
        be retained.
    Returns
    -------
    :obj:`None`

    Filters and updates `adata`.
    """
    # apply the filter
    print(
        "removed {} cells".format(int((adata.obs["Pass Cell Filter"] == False).sum()))
    )
    adata._inplace_subset_obs(adata.obs[adata.obs["Pass Cell Filter"]].index)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    # recalculate cell stats after some genes were removed
    adata.obs["Num Proteins"] = (adata.X > 0).sum(axis=1)
    adata.obs["Log2 Sum S/N"] = np.log2(adata.X.sum(axis=1))


def plot_batch_qc(
    adata: AnnData,
    labels: Sequence[str] = ("Log2 Sum S/N", "Num Proteins"),
    groupby: Sequence[str] = (
        "File ID",
        "Channel",
        "Row",
        "Column",
        "Sorted Population",
        "Gated Population",
        "Plate",
    ),
    figsize: Tuple[float, float] = (8, 15),
):
    """Plot QC parameters grouped by possible sources of batch effects.

    Parameters
    ----------
    adata
        The annotated data matrix.
    labels
        QC parameters to plot.
    groupby
        Batch variables to group by.
    figsize
        Size of the plotted figure.
    Returns
    -------
    :obj:`None`
    """
    fig, axs = plt.subplots(len(groupby), len(labels), figsize=figsize)
    for r in range(len(groupby)):
        for c in range(len(labels)):
            ax = axs[r, c]
            sc.pl.violin(
                adata,
                labels[c],
                jitter=1,
                size=2,
                groupby=groupby[r],
                show=False,
                rotation=45,
                ax=ax,
            )
            ax.set_ylabel(labels[c])
            ax.set_xlabel(groupby[r])
            ax.grid(alpha=0.5)
    fig.tight_layout()


def plot_plate_qc(
    adata: AnnData,
    labels: Sequence[str] = ("Log2 Sum S/N", "Num Proteins"),
    figsize: Tuple[float, float] = (10, 3),
):
    """Plot QC parameters for each cell on the plate layout.

    Parameters
    ----------
    adata
        The annotated data matrix.
    labels
        QC parameters to plot.
    figsize
        Size of the plotted figure.
    Returns
    -------
    :obj:`None`
    """
    rows = [
        "A",
        "B",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "J",
        "K",
        "L",
        "M",
        "N",
        "O",
        "P",
    ]
    columns = range(1, 25)
    plate = pd.DataFrame(
        index=pd.MultiIndex.from_tuples(
            [x for subl in [[(r, c) for c in columns] for r in rows] for x in subl],
            names=("Row", "Column"),
        )
    )
    plates = adata.obs["Plate"].unique()
    color_map = copy(plt.get_cmap("copper"))
    color_map.set_bad("lightgray")
    fig, axs = plt.subplots(len(plates), len(labels), figsize=figsize)
    for i, p in enumerate(plates):
        for c in range(len(labels)):
            if len(plates) == 1:
                ax = axs[c]
            else:
                ax = axs[i, c]
            df = adata.obs[adata.obs["Plate"] == p].copy()
            df["Column"] = df["Column"].astype(int)
            df = plate.merge(
                df.set_index(["Row", "Column"])[labels[c]],
                left_index=True,
                right_index=True,
                how="left",
            )[labels[c]].unstack()
            sns.heatmap(df, ax=ax, linewidths=0.3, square=True, cmap=color_map)
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
            ax.xaxis.label.set_visible(False)
            ax.set_ylabel("")
            ax.xaxis.tick_top()
            # ax.set_facecolor("lightgrey")

    if len(plates) == 1:
        for ax, col in zip(axs, labels):
            ax.set_title(col)
        for ax, row in zip(axs, plates):
            ax.set_ylabel(row, rotation=90)
    else:
        for ax, col in zip(axs[0], labels):
            ax.set_title(col)
        for ax, row in zip(axs[:, 0], plates):
            ax.set_ylabel(row, rotation=90)

    fig.tight_layout()


def plot_data_completeness(adata: AnnData, figsize: Tuple[float, float] = (2.9, 2)):
    """Plot the data completeness of the expression matrix.

    Parameters
    ----------
    adata
        The annotated data matrix.
    figsize
        Size of the plotted figure.
    Returns
    -------
    :obj:`None`
    """
    print(
        "mean protein IDs per cell: {}".format(
            round((~pd.DataFrame(adata.X.T).replace(0, pd.NA).isna()).sum().mean(), 1)
        )
    )
    print(
        "median protein IDs per cell: {}".format(
            round((~pd.DataFrame(adata.X.T).replace(0, pd.NA).isna()).sum().median(), 1)
        )
    )
    print(
        "percent missing values: {}".format(
            round(
                100
                * (pd.DataFrame(adata.X.T).replace(0, pd.NA).isna()).sum().sum()
                / (adata.X.shape[0] * adata.X.shape[1]),
                2,
            )
        )
    )

    res = {}
    frac_missing = (
        (adata.X.shape[0] - np.count_nonzero(adata.X, axis=0)) / adata.X.shape[0] * 100
    )
    for i in [0, 20, 40, 60, 80, 99]:
        res[i] = sum(frac_missing <= i)
    fig, ax = plt.subplots(figsize=figsize)
    plt.bar(res.keys(), res.values(), width=5)
    plt.axhline(adata.X.shape[1], color="black", linestyle="--")
    plt.xticks(list(res.keys()), [100, 80, 60, 40, 20, 1])
    plt.text(
        10,
        adata.X.shape[1] * 0.97,
        str(adata.X.shape[1]),
        horizontalalignment="center",
        verticalalignment="top",
    )
    # yt = ax.get_yticks()
    # yt = np.append(yt, adata.X.shape[1])
    # ax.set_yticks(yt)
    plt.xlabel("Percent data completeness")
    plt.ylabel("Number of proteins")
    plt.grid()
    plt.title("Data completeness across {} cells".format(adata.shape[0]))
    fig.tight_layout()


def plot_facs_qc(
    adata: AnnData,
    marker_x: str,
    marker_y: str,
    labels: Sequence[str] = (
        "Log2 Sum S/N",
        "Num Proteins",
        "Sorted Population",
        "Gated Population",
    ),
    figsize: Tuple[float, float] = (2.7, 1.66),
):
    """Plot selected facs parameter of each cell and overlay with annotation.

    Parameters
    ----------
    adata
        The annotated data matrix.
    marker_x
        The marker from `adata.obs` to plot on x-axis.
    marker_y
        The marker from `adata.obs` to plot on y-axis.
    labels
        The annotations to overlay the cells with.
    figsize
        Size of the plotted figure.
    Returns
    -------
    A list of two :class:`~matplotlib.figure` objects.
    """
    # temporarily change figsize
    curr = plt.rcParams["figure.figsize"]
    plt.rcParams["figure.figsize"] = figsize

    figs = []
    for l in labels:
        axs = sc.pl.scatter(
            adata, x=marker_x, y=marker_y, color=l, title=l, size=25, show=False
        )
        plt.tight_layout()
        figs.append(plt.gcf())

    plt.rcParams["figure.figsize"] = curr

    return figs


def impute(adata: AnnData, **kwargs):
    """Impute missing values in `adata`

    Parameters
    ----------
    adata
        The annotated data matrix.
    **kwargs
        Parameters are passed to :class:`sklearn.impute.KNNImputer`.
    Returns
    -------
    :obj:`None`

    Updates `adata` with the imputed expression values.
    """
    if 0 in adata.X:
        from sklearn.impute import KNNImputer

        imputer = KNNImputer(missing_values=0, **kwargs)
        if (adata.X == 0).all(axis=0).any():
            print("remove all-zero proteins:")
            sc.pp.filter_genes(adata, min_cells=1)
        adata.X = imputer.fit_transform(adata.X)
    else:
        warnings.warn("No zeros in adata")


def find_embedding_params(
    adata: AnnData,
    labels: str = "Sorted Population",
    use_rep: Literal["umap", "fa", "diffmap"] = "umap",
    impute_knn: Sequence[int] = (5,),
    embedding_knn: Sequence[int] = (15,),
    figsize: Tuple[float, float] = (4, 1.3),
    save_all_plots: bool = True,
    save_path: Union[str, None] = None,
    size: float = 20,
):
    """Calculate the silhouette score of given labels across embeddings with different parameters.
    A high silhouette score indicates good cell separation based on labels.
    Expects log2 transformed values with 0 denoting missing values.
    Test the following parameters:
    * fraction of valid values protein filter
    * knn for imputation
    * knn for neighborhood graph

    Parameters
    ----------
    adata
        The annotated data matrix.
    labels
        The cell label to optimize separation on.
    use_rep
        The representation to optimize on.
    impute_knn
        The numbers of nearest neighbors to test for KNN imputation.
    embedding_knn
        The number of nearest neighbors to test for the neighborhood graph for the embedding.
    figsize
        Size of the plotted figure.
    save_all_plots
        Saves plots of all generated embeddings in 'save_path'.
    save_path
        The directory to save all the embedding plots.
    size
        Size of the plotted cells in all plots.
    Returns
    -------
     A :class:`~pandas.DataFrame` with the scores for each test.
    """

    from sklearn import metrics

    current_verbosity = sc.settings.verbosity
    sc.settings.verbosity = 1

    if save_all_plots:
        if save_path is None:
            raise ValueError("Specify 'save_path' to save all plots.")

    scores = []
    for i, i_knn in enumerate(impute_knn):
        print("impute_knn: {}/{}".format(i + 1, len(impute_knn)))
        fracs = np.linspace(0.0, 1, 11)
        for j, frac in enumerate(fracs):
            print("\tfraction filter: {}/{}".format(j + 1, len(fracs)))
            ad = adata[
                :, (adata.X != 0).sum(axis=0) >= adata.obs.shape[0] * frac
            ].copy()
            if 0 in ad.X:
                impute(ad, n_neighbors=i_knn)
            sc.pp.scale(ad)
            # scale can introduce nan genes, remove if there
            if np.isnan(ad.X).any():
                retain = [p for p in ad.var.index if not (p in ad.var.index[list(set(np.where(np.isnan(ad.X))[1]))])]
                ad = ad[:, retain].copy()
            sc.pp.pca(ad)
            for e_knn in embedding_knn:
                sc.pp.neighbors(ad, n_neighbors=e_knn)
                if use_rep == "umap":
                    sc.tl.umap(ad)
                    score = metrics.silhouette_score(ad.obsm["X_umap"], ad.obs[labels])
                    scores.append(
                        {
                            "impute_knn": i_knn,
                            "frac_valid": frac,
                            "embedding_knn": e_knn,
                            "silhouette_score": score,
                        }
                    )
                elif use_rep == "fa":
                    sc.tl.draw_graph(ad)
                    score = metrics.silhouette_score(
                        ad.obsm["X_draw_graph_fa"], ad.obs[labels]
                    )
                    scores.append(
                        {
                            "impute_knn": i_knn,
                            "frac_valid": frac,
                            "embedding_knn": e_knn,
                            "silhouette_score": score,
                        }
                    )
                elif use_rep == "diffmap":
                    sc.tl.diffmap(ad, n_comps=3)
                    score = metrics.silhouette_score(
                        ad.obsm["X_diffmap"][:, 1:3], ad.obs[labels]
                    )
                    scores.append(
                        {
                            "impute_knn": i_knn,
                            "frac_valid": frac,
                            "embedding_knn": e_knn,
                            "silhouette_score": score,
                        }
                    )
                if save_all_plots:
                    title = f"sc:{score:.2f} i_n:{i_knn} frac:{frac:.1f} e_n:{e_knn}"
                    fig, ax = plt.subplots(figsize=(2, 1.5))
                    if use_rep == "umap":
                        sc.pl.umap(
                            ad, color=labels, ax=ax, title=title, size=size, show=False
                        )
                    elif use_rep == "fa":
                        sc.pl.draw_graph(
                            ad, color=labels, ax=ax, title=title, size=size, show=False
                        )
                    elif use_rep == "diffmap":
                        sc.pl.diffmap(
                            ad, color=labels, ax=ax, title=title, size=size, show=False
                        )
                    ax.get_legend().remove()
                    plt.savefig(
                        "{}/{}.pdf".format(
                            save_path, title.replace(":", "_").replace(" ", "_")
                        )
                    )
                    plt.close()

    sc.settings.verbosity = current_verbosity
    scores = pd.DataFrame(scores).sort_values("silhouette_score", ascending=False)

    # plot the scores
    scores["embedding_knn"] = scores["embedding_knn"].astype("category")
    scores["impute_knn"] = scores["impute_knn"].astype("category")
    scores = scores.rename(
        columns={
            "frac_valid": "Fraction of valid values filter",
            "silhouette_score": "Silhouette Score",
        }
    )
    fig, ax = plt.subplots(figsize=figsize)
    # if there are only two dimensions, plot simpler plot
    if len(impute_knn) == 1 and len(embedding_knn) == 1:
        g = sns.scatterplot(
            x="Fraction of valid values filter",
            y="Silhouette Score",
            data=scores,
            s=20,
            legend=False,
            ax=ax,
        )
        ax.grid()
    else:
        g = sns.scatterplot(
            x="Fraction of valid values filter",
            y="Silhouette Score",
            data=scores,
            hue="embedding_knn",
            style="impute_knn",
            s=20,
            legend="full",
            ax=ax,
        )
        g.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        for handle in ax.get_legend().legendHandles:
            handle.set_sizes([20])

    return scores


def find_appropriate_neighbors(
    adata: AnnData,
    color: str = "Gated Population",
    test: Sequence[int] = (5, 10, 15, 20, 30, 40, 50),
    diffmap_denoise: bool = False,
    figsize: Tuple[float, float] = (4, 12),
):
    """Plot UMAP and FA embeddings of the data using different numbers of neighbors.
    Helps to find a number of neighbors that is approriate for the specific dataset.

    Parameters
    ----------
    adata
        The annotated data matrix.
    color
        Key for annotation of observation.
    test
        The different numbers of neighbors to test
    diffmap_denoise
        Apply diffusion map on data to denoise the data and subsequently build
        neighborhood graph from diffusion map.
    figsize
        Size of the plotted figure.
    Returns
    -------
    :obj:`None`
    """

    # temporarily change verbosity
    current_verbosity = sc.settings.verbosity
    sc.settings.verbosity = 1
    neighbors = test

    fig, axs = plt.subplots(
        nrows=len(neighbors),
        ncols=2,
        gridspec_kw=dict(hspace=0.35, wspace=0.2),
        figsize=figsize,
    )

    for i, n in enumerate(neighbors):
        sc.pp.neighbors(adata, n_neighbors=n)
        if diffmap_denoise:
            sc.tl.diffmap(adata, n_comps=n)
            sc.pp.neighbors(adata, n_neighbors=n, use_rep="X_diffmap")
        sc.tl.umap(adata)
        sc.pl.umap(
            adata, color=color, title=[n], show=False, ax=axs[i][0], size=cellsize
        )
        axs[i][0].get_legend().remove()
        sc.tl.draw_graph(adata)
        sc.pl.draw_graph(
            adata, color=color, title=[n], show=False, ax=axs[i][1], size=cellsize
        )
        axs[i][1].get_legend().remove()

    sc.settings.verbosity = current_verbosity


def de_test(
    adata: AnnData,
    by: str,
    group1: str,
    group2: Union[str, None] = None,
    test: Literal["mannwhitneyu", "welch"] = "welch",
    key: str = "de_test",
    use_raw: bool = True,
    ignore: Union[float, None] = 0,
    is_log: bool = True,
):
    """Test differential protein expression between groups of cells.
    Performs Mann-Whitney rank test or Welch test and corrects p-values with Benjamini/Hochberg.
    Ignore missing values by setting `ignore` to `0` or `np.nan`

    Parameters
    ----------
    adata
        The annotated data matrix.
    by
        The variable in `adata.obs` containing the group annotation.
    group1
        The name of the first group.
    group2
        The name of the second group. If None, test `group1` against rest.
    test
        The test method to use.
    key
        The key under which the results are stored in `adata.uns`.
    use_raw
        If the raw data in `adata.raw` should be used.
    ignore
        A value that should be ignored during testing. Used to ignore
        missing values (0).
    is_log
        If the data is log transformed or not. Relevant for the fold change calculation.
        If False, normal foldchanges instead of log2foldchanges are calculated.

    Returns
    -------
    :obj:`None`

    Results are stored in `adata.uns[key]`.
    """
    from scipy.stats import mannwhitneyu
    from scipy.stats import ttest_ind
    from statsmodels.stats.multitest import multipletests

    results = []
    if use_raw:
        genes = adata.raw.var_names
        group1_vals = adata.raw[adata.obs[by] == group1].X
        if group2 is not None:
            group2_vals = adata.raw[adata.obs[by] == group2].X
        else:
            group2_vals = adata.raw[adata.obs[by] != group1].X
    else:
        genes = adata.var_names
        group1_vals = adata[adata.obs[by] == group1].X
        if group2 is not None:
            group2_vals = adata[adata.obs[by] == group2].X
        else:
            group2_vals = adata[adata.obs[by] != group1].X
    for i, gene in enumerate(genes):
        vals1 = group1_vals[:, i]
        vals2 = group2_vals[:, i]
        if ignore is not None:
            if np.isnan(ignore):
                vals1 = vals1[~np.isnan(vals1)]
                vals2 = vals2[~np.isnan(vals2)]
            else:
                vals1 = vals1[vals1 != ignore]
                vals2 = vals2[vals2 != ignore]
        if len(vals1) > 2 and len(vals2) > 2:  # at least 3 values per group
            if test == "welch":
                test_res = ttest_ind(vals1, vals2, equal_var=False)
            if test == "mannwhitneyu":
                test_res = mannwhitneyu(vals1, vals2, alternative="two-sided")
            if is_log:
                fold_change = np.mean(vals1) - np.mean(vals2)
                # fold_change = np.log2(2**(np.mean(vals1) + 1e-9) / 2**(np.mean(vals2) + 1e-9))
            else:
                fold_change = np.mean(vals1) / np.mean(vals2)
            results.append(
                {
                    "gene": gene,
                    "pval": test_res[1],
                    "log2foldchange": fold_change,
                    "size1": len(vals1),
                    "size2": len(vals2),
                }
            )
        else:
            results.append(
                {
                    "gene": gene,
                    "pval": 1.0,
                    "log2foldchange": np.nan,
                    "size1": len(vals1),
                    "size2": len(vals2),
                }
            )

    results = pd.DataFrame.from_records(results)
    results["pval"] = results["pval"].fillna(1.0)
    _, pvals_adj, _, _ = multipletests(results["pval"], alpha=0.05, method="fdr_bh")
    results["pval_adj"] = pvals_adj
    results = results.sort_values("log2foldchange")
    adata.uns[key] = {
        "by": by,
        "group1": group1,
        "group2": group2 if group2 is not None else "Rest",
        "results": results,
    }


def plot_volcano(
    adata: AnnData,
    test_key: str,
    adj_pval_cutoff: float = 0.05,
    log2foldchange_cutoff: float = 0.5,
    title: Union[str, None] = None,
    figsize: Tuple[float, float] = (4, 4),
):
    """Plot the results of :func:`~sceptre.de_test` as volcano plot.

    Parameters
    ----------
    adata
        The annotated data matrix.
    test_key
        The key under which the results are stored in `adata.uns`.
    adj_pval_cutoff
        The significance cutoff of the adjusted p-values.
    log2foldchange_cutoff
        The effect size cutoff of the log2 fold change.
    title
        The title of the plot. If `None`, sets automatically from the group names.
    key
        The key under which the results are stored in `adata.uns`.
    figsize
        Size of the plotted figure.
    Returns
    -------
    :obj:`None`
    """
    from adjustText import adjust_text
    from matplotlib import patches

    # get data
    df = adata.uns[test_key]["results"].dropna(subset=["log2foldchange"], axis=0).copy()
    df["nlog_pval_adj"] = -np.log10(df["pval_adj"])
    df["up"] = (df["pval_adj"] <= adj_pval_cutoff) & (
        df["log2foldchange"] >= log2foldchange_cutoff
    )
    df["down"] = (df["pval_adj"] <= adj_pval_cutoff) & (
        df["log2foldchange"] <= -log2foldchange_cutoff
    )
    df["unchanged"] = (df["up"] == False) & (df["down"] == False)
    # plot
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(
        x=df[df["up"]]["log2foldchange"], y=df[df["up"]]["nlog_pval_adj"], color="green"
    )
    ax.scatter(
        x=df[df["down"]]["log2foldchange"],
        y=df[df["down"]]["nlog_pval_adj"],
        color="red",
    )
    ax.scatter(
        x=df[df["unchanged"]]["log2foldchange"],
        y=df[df["unchanged"]]["nlog_pval_adj"],
        color="grey",
    )
    hl = plt.axhline(
        -np.log10(adj_pval_cutoff), color="black", linestyle="--", alpha=0.5
    )
    vl1 = plt.axvline(-log2foldchange_cutoff, color="black", linestyle="--", alpha=0.5)
    vl2 = plt.axvline(log2foldchange_cutoff, color="black", linestyle="--", alpha=0.5)
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10 P-value")
    # expand the axis limits
    x_min, x_max = df["log2foldchange"].min(), df["log2foldchange"].max()
    x_len = x_max - x_min
    plt.xlim(x_min - x_len * 0.1, x_max + x_len * 0.1)
    plt.ylim(0, df["nlog_pval_adj"].max() * 1.15)
    # add gene labels
    texts = [
        ax.text(
            df.loc[i, "log2foldchange"],
            df.loc[i, "nlog_pval_adj"],
            df.loc[i, "gene"],
            ha="center",
            va="center",
        )
        for i in df[~df["unchanged"]].index
    ]
    # add boxes to be avoided by the text labels during adjustment
    patch_ymax = patches.Rectangle(
        (-3, df["nlog_pval_adj"].max() * 1.10), 15, 3, fill=False, alpha=0.0
    )
    ax.add_patch(patch_ymax)
    patch_hl = patches.Rectangle(
        (-3, 0), 15, -np.log10(adj_pval_cutoff), fill=False, alpha=0.0
    )
    ax.add_patch(patch_hl)
    # adjust the text labels
    adjust_text(
        texts,
        x=df["log2foldchange"].values,
        y=df["nlog_pval_adj"].values,
        arrowprops=dict(arrowstyle="-", color="grey", alpha=0.5, lw=0.5),
        add_objects=[patch_ymax, patch_hl],
        expand_objects=(1, 1),
        lim=100,
    )
    if title:
        plt.title(title)
    else:
        plt.title(
            "{} vs. {}".format(
                adata.uns[test_key]["group1"], adata.uns[test_key]["group2"]
            )
        )
    fig.tight_layout()


def plot_de_heatmap(
    adata: AnnData,
    test_key: str,
    adj_pval_cutoff: float = 0.05,
    log2foldchange_cutoff: float = 0.5,
    use_raw: bool = True,
    dendrogram: bool = False,
    standard_scale: Literal["var", "obs"] = "var",
    figsize: Tuple[float, float] = (4, 4),
    **kwargs,
):
    """Plot the results of :func:`~sceptre.de_test` as heatmap.

    Parameters
    ----------
    adata
        The annotated data matrix.
    test_key
        The key under which the results are stored in `adata.uns`.
    adj_pval_cutoff
        The significance cutoff of the adjusted p-values.
    log2foldchange_cutoff
        The effect size cutoff of the log2 fold change.
    use_raw
        Use `raw` attribute of `adata` if present.
    dendrogram
        If True or a valid dendrogram key, a dendrogram based on the hierarchical
        clustering between the `groupby` categories is added.
        The dendrogram information is computed using :func:`scanpy.tl.dendrogram`.
        If `tl.dendrogram` has not been called previously the function is called
        with default parameters.
    standard_scale
        Whether or not to standardize that dimension between 0 and 1,
        meaning for each variable or group,
        subtract the minimum and divide each by its maximum.
    figsize
        Size of the plotted figure.
    kwargs
        Parameters are passed to :func:`scanpy.pl.heatmap`.

    Returns
    -------
    :obj:`None`
    """
    res = adata.uns[test_key]["results"]
    group1_genes = res.loc[
        (res["pval_adj"] <= adj_pval_cutoff)
        & (res["log2foldchange"] >= log2foldchange_cutoff),
        "gene",
    ]
    group2_genes = res.loc[
        (res["pval_adj"] <= adj_pval_cutoff)
        & (res["log2foldchange"] <= -log2foldchange_cutoff),
        "gene",
    ]
    axs = sc.pl.heatmap(
        adata,
        groupby=adata.uns[test_key]["by"],
        var_names={
            adata.uns[test_key]["group2"]: group2_genes,
            adata.uns[test_key]["group1"]: group1_genes,
        },
        dendrogram=dendrogram,
        standard_scale=standard_scale,
        use_raw=use_raw,
        show=False,
        figsize=figsize,
        **kwargs,
    )


def enrichment_test(
    adata: AnnData,
    gene_set: Sequence[str],
    categories: Sequence[str],
    background: Union[pd.DataFrame, None] = None,
    sep: str = ";",
    key: str = "enrichment_test",
    pval_thesh: float = 0.05,
    use_raw: bool = True,
):
    """Perform a term enrichment test on the selected genes in selected term category.
    If no background is provided, all genes in the matrix are used as background.

    Parameters
    ----------
    adata
        The annotated data matrix.
    gene_set
        Genes to be tested against the background.
    background
        Dataframe of genes in the rows and categories in the columns.
    categories
        The names of the variables in `adata.var` that store the term annotations.
        E.g. `'Biological Process'`.
    sep
        The separator used in the `categories` variables.
    key
        The key under which the results are stored in `adata.uns`.
    pval_thesh
        Filter the output table based on the adjusted p-value.
    use_raw
        Use genes in 'adata.raw'.
    Returns
    -------
    :obj:`None`

    Results are stored in `adata.uns[key]`.
    """

    from scipy.stats import hypergeom
    from statsmodels.stats.multitest import multipletests

    if use_raw:
        ad = adata.raw.to_adata()
    else:
        ad = adata

    results = pd.DataFrame()
    for cat in categories:
        gene_terms = (
            ad.var[cat]
            .astype(str)[~(ad.var[cat] == "nan")]
            .apply(lambda x: x.split(sep))
        )
        # remove 'nan' term
        gene_terms = gene_terms.apply(lambda x: [t for t in x if t != 'nan'])

        all_gene_set_terms = pd.DataFrame(
            index=set([l for subl in gene_terms.loc[gene_set].values for l in subl])
        )
        # skip if all_gene_set_terms is empty
        if len(all_gene_set_terms) == 0:
            continue
        for term in all_gene_set_terms.index:
            x = (
                gene_terms[gene_set].apply(lambda x: term in x).sum()
            )  # number of test proteins with the term
            if background is not None:
                M = len(background) # number of background proteins
            else:
                M = len(gene_terms)  # number of background proteins
            N = len(gene_set)  # number of test proteins
            if background is not None:
                n = background[cat].str.contains(term, regex=False).sum() # number of background proteins with the term
            else:
                n = gene_terms.apply(
                    lambda x: term in x
                ).sum()  # number of background proteins with the term
            p = hypergeom.sf(x - 1, M, n, N)
            expected = n / M * N
            enrichment = x / expected
            all_gene_set_terms.loc[
                term,
                [
                    "size background",
                    "# in background",
                    "size subset",
                    "# in subset",
                    "expected",
                    "enrichment",
                    "pval",
                ],
            ] = (M, n, N, x, expected, enrichment, p)

        _, pvals_adj, _, _ = multipletests(
            all_gene_set_terms["pval"], alpha=0.05, method="fdr_bh"
        )
        all_gene_set_terms["pvals_adj"] = pvals_adj
        if pval_thesh:
            all_gene_set_terms = all_gene_set_terms[
                all_gene_set_terms["pvals_adj"] <= pval_thesh
            ]
        all_gene_set_terms["Category"] = cat
        results = results.append(all_gene_set_terms)

    adata.uns[key] = results.sort_values("pvals_adj")


def ordered_clustermap(
    adata: AnnData,
    genes: Sequence[str],
    n_cluster: int = 5,
    annotations: Sequence[str] = ("CD34 APC-Cy7-A", "CD38 PE-A", "dpt_pseudotime"),
    annotations_cmap: Mapping = {"CD34 APC-Cy7-A": "cividis", "CD38 PE-A": "cividis", "dpt_pseudotime": "copper"},
    moving_avg: int = 10,
    use_raw: bool = True,
    plot_imputed: bool = False,
    show_gene_names="auto",
    figsize: Tuple[float, float] = (4, 4),
):
    """Plot a heatmap of cells and `genes`.
    Cells are ordered by `dpt_pseudotime`, genes are clustered hierarchically.

    Parameters
    ----------
    adata
        The annotated data matrix.
    genes
        Genes to be plotted.
    n_cluster
        The number of clusters to form from the hierarchical clustering.
    annotations
        The cell annotations to plot from `adata.obs`.
    annotations_cmap
        The colormaps for the annotations.
    moving_avg
        Apply moving average along the ordered cells. Applied on genes and annotations.
        If set to 1, no averaging is applied.
    use_raw
        Use geneexpression from 'adata.raw'. Missing values are imputed for the clustering.
    plot_imputed
        Plot the imputed values instead of the raw values. Only applies if `use_raw` is `True`.
    show_gene_names
        Plot the gene names.
    figsize
        Size of the plotted figure.
    Returns
    -------
    DataFrame of gene-cluster mapping
    Figure of heatmap
    Figure of cluster legend
    """

    from scipy.cluster import hierarchy
    from matplotlib import cm
    import fastcluster

    sorted_cell_idx = adata.obs["dpt_pseudotime"].sort_values().index
    if use_raw:
        x = pd.DataFrame(
            adata.raw[sorted_cell_idx, genes].X.T, index=genes, columns=sorted_cell_idx
        )
    else:
        x = pd.DataFrame(
            adata[sorted_cell_idx, genes].X.T, index=genes, columns=sorted_cell_idx
        )
    x.index.name = None

    # cluster the genes
    if use_raw:
        # using an imputed matrix
        adata_imp = adata.raw.to_adata()
        impute(adata_imp)
        sc.pp.scale(adata_imp)
        if plot_imputed:
            x = pd.DataFrame(
                adata_imp[sorted_cell_idx, genes].X.T,
                index=genes,
                columns=sorted_cell_idx,
            )
    else:
        adata_imp = adata.copy()
        sc.pp.scale(adata_imp)
    x_imp = adata_imp[:, genes].X.T
    # linkage = hierarchy.linkage(x_imp, method='ward')
    linkage = fastcluster.linkage(x_imp, method="ward")
    clusters = hierarchy.fcluster(linkage, criterion="maxclust", t=n_cluster)
    cluster_cols = pd.Series(plt.cm.tab10(clusters).tolist(), index=genes)

    anno = adata.obs.loc[sorted_cell_idx, annotations].copy()

    # moving average across cells
    if moving_avg > 1:
        x = x.apply(
            lambda x: pd.Series(sc._utils.moving_average(x.values, n=moving_avg)),
            axis=1,
        )
        anno = adata.obs.loc[sorted_cell_idx, annotations].apply(
            lambda x: pd.Series(sc._utils.moving_average(x.values, n=moving_avg)),
            axis=0,
        )

    # create colors of annotations
    for a in anno.columns:
        if a in annotations_cmap.keys():
            cmap = plt.get_cmap(annotations_cmap[a])
        else:
            cmap = "copper"
        anno[a] = (
            cm.ScalarMappable(
                norm=plt.Normalize(vmin=anno[a].min(), vmax=anno[a].max()), cmap=cmap
            )
            .to_rgba(anno[a])
            .tolist()
        )

    anno = anno.rename(columns={"dpt_pseudotime": "Pseudotime"})
    cl_map = sns.clustermap(
        x,
        col_cluster=False,
        col_colors=anno,
        row_colors=cluster_cols,
        row_linkage=linkage,
        standard_scale=0,
        cmap="viridis",
        cbar_pos=None,
        xticklabels=False,
        yticklabels=show_gene_names,
        figsize=figsize,
    )

    cls = pd.DataFrame(
        zip(genes, clusters, cluster_cols), columns=["gene", "cluster", "cluster_col"]
    )
    cl = cls.drop_duplicates(subset="cluster").sort_values("cluster")
    cl_legend, ax = plt.subplots(figsize=(1, 1))
    ax.imshow(np.array([cl["cluster_col"].tolist()]))
    ticks = [0] + cl["cluster"].tolist()
    ax.set_xticks(ticks[:-1])
    ax.set_xticklabels(ticks[1:])
    ax.get_yaxis().set_visible(False)
    plt.title("Clusters")

    return cls[["gene", "cluster"]], cl_map, cl_legend
