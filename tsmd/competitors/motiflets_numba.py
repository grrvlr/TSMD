import numpy as np

import motiflets.motiflets as ml 

class Motiflets_numba:
    """k-Motiflets algorithm for motif discovery.

    Parameters
    ----------
    k_max : int 
        Maximum number of occurences of a single motif.
    min_wlen : int
        Minimium window length.
    max_wlen : int 
        Maximum window length.
    elbow_deviation : float, optional (default=1.0)
        The minimal absolute deviation needed to detect an elbow.
        It measures the absolute change in deviation from k to k+1.
        1.05 corresponds to 5% increase in deviation.
    slack : float, optional (default=0.5)
        Defines an exclusion zone around each subsequence to avoid trivial matches.
        Defined as percentage of m. E.g. 0.5 is equal to half the window length.
    
    Attributes
    ----------
    prediction_mask_ : np.ndarray of shape (n_patterns, n_samples)
        Binary mask indicating the presence of motifs across the signal.  
        Each row corresponds to one discovered motif, and each column to a time step.  
        A value of 1 means the motif is present at that time step, and 0 means it is not.
        """
    def __init__(
            self,
            k_max,
            min_wlen,
            max_wlen,
            elbow_deviation=1.00,
            slack=0.5,
    ):
        self.elbow_deviation = elbow_deviation
        self.slack = slack

        self.motif_length_range = np.arange(min_wlen,max_wlen+1)
        self.motif_length = 0
    

        self.k_max = k_max +1 

    def fit(self,signal):
        """Fit Motiflets
        
        Parameters
        ----------
        signal : numpy array of shape (n_samples, )
            The input samples (time series length).
        
        Returns
        -------
        self : object
            Fitted estimator.
        """
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
