import time

import matplotlib
import numpy as np
from scipy.stats import zscore

import tsmd.competitors.competitors_tools.motiflets_tools as ml

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


class Motiflets:

    def __init__(
            self,
            k_max,
            min_wlen,
            max_wlen,
            elbow_deviation=1.00,
            slack=0.5,
    ):
        """Computes the AU_EF plot to extract the best motif lengths

            This is the method to find and plot the characteristic motif-lengths, for k in
            [2...k_max], using the area AU-EF plot.

            Details are given within the paper 5.2 Learning Motif Length l.

            Parameters
            ----------
            ds_name: String
                Name of the time series for displaying
            series: array-like
                the TS
            ground_truth: pd.Series
                Ground-truth information as pd.Series.
            elbow_deviation : float, default=1.00
                The minimal absolute deviation needed to detect an elbow.
                It measures the absolute change in deviation from k to k+1.
                1.05 corresponds to 5% increase in deviation.
            slack: float
                Defines an exclusion zone around each subsequence to avoid trivial matches.
                Defined as percentage of m. E.g. 0.5 is equal to half the window length.
            n_jobs : int
                Number of jobs to be used.

            Returns
            -------
            best_motif_length: int
                The motif length that maximizes the AU-EF.

            """
    
        self.elbow_deviation = elbow_deviation
        self.slack = slack

        self.motif_length_range = np.arange(min_wlen,max_wlen+1)
        self.motif_length = 0
    

        self.k_max = k_max +1 

    def fit(self,signal):
        self.signal=signal
        self.fit_motif_length()
        self.fit_k_elbow()

    def fit_motif_length(
            self,
            subsample=2,
            # plot=True,
            # plot_elbows=False,
            # plot_motifs_as_grid=True,
            # plot_best_only=True
    ):
        """Computes the AU_EF plot to extract the best motif lengths

            This is the method to find and plot the characteristic motif-lengths, for k in
            [2...k_max], using the area AU-EF plot.

            Details are given within the paper 5.2 Learning Motif Length l.

            Parameters
            ----------
            k_max: int
                use [2...k_max] to compute the elbow plot.
            motif_length_range: array-like
                the interval of lengths

            Returns
            -------
            best_motif_length: int
                The motif length that maximizes the AU-EF.

            """


        self.motif_length = ml.find_au_ef_motif_length(self.signal,self.k_max,self.motif_length_range,exclusion=None,elbow_deviation=self.elbow_deviation,slack=self.slack,subsample=subsample)[0]

        return self.motif_length

    def fit_k_elbow(
            self,
            motif_length=None,  # if None, use best_motif_length
            exclusion=None,
    ):
        """Plots the elbow-plot for k-Motiflets.

            This is the method to find and plot the characteristic k-Motiflets within range
            [2...k_max] for given a `motif_length` using elbow-plots.

            Details are given within the paper Section 5.1 Learning meaningful k.

            Parameters
            ----------
            k_max: int
                use [2...k_max] to compute the elbow plot (user parameter).
            motif_length: int
                the length of the motif (user parameter)
            exclusion: 2d-array
                exclusion zone - use when searching for the TOP-2 motiflets
            filter: bool, default=True
                filters overlapping motiflets from the result,
            plot_elbows: bool, default=False
                plots the elbow ploints into the plot

            Returns
            -------
            Tuple
                dists:          distances for each k in [2...k_max]
                candidates:     motifset-candidates for each k
                elbow_points:   elbow-points

            """

        if motif_length is None:
            motif_length = self.motif_length
        else:
            self.motif_length = motif_length
        self.dists, self.motiflets, self.elbow_points, _ = ml.search_k_motiflets_elbow(
        self.k_max,
        self.signal,
        motif_length,
        exclusion=exclusion,
        elbow_deviation=self.elbow_deviation,
        slack=self.slack)

        return self.dists, self.motiflets, self.elbow_points

    @property
    def prediction_mask_(self)->np.ndarray: 
        n_motifs=self.elbow_points.shape[0]
        mask=np.zeros((n_motifs,self.signal.shape[0]))
        for i in range(self.elbow_points.shape[0]):
            elbow=self.elbow_points[i]
            motif_starts=self.motiflets[elbow]
            for j in range(motif_starts.shape[0]):
                mask[i,motif_starts[j]:motif_starts[j]+self.motif_length]=1
        return mask