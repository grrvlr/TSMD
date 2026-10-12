
import copy
import functools
from typing import Literal

import numba as nb
import numpy as np

type_motif_set = nb.types.Tuple(
    (
        nb.types.Array(nb.types.boolean, 1, "C"),  # The mask
        nb.types.Array(nb.types.float32, 3, "C"),  # The motif locations
    )
)


class SubTSMD:
    """SubTSMD algorithm for subspace motif discovery.

    Parameters
    ----------
    univariate_motif_discovery : object
        An instance of a motif discovery algorithm of this package (e.g., LocoMotif), which
        implements ``fit(signal)`` and ``prediction_mask_``. The instance is copied for each
        attribute, so it is never fitted itself.
    delta : float, optional (default=0.5)
        The threshold for merging motifs across dimensions. A value between 0.5 and 1.
    linkage : str, optional (default="average")
        The linkage method to use when merging motifs across subspaces. Can be either "complete" or "average".
    size_constraint : int, optional (default=1)
        The minimum number of attributes in which a motif set must occur.
    normalize : bool, optional (default=True)
        Whether to z-normalize each attribute before discovering the univariate motifs.

    Attributes
    ----------
    prediction_mask_ : np.ndarray of shape (n_patterns, n_samples)
        Binary mask indicating the presence of motifs across the signal.
        Each row corresponds to one discovered motif, and each column to a time step.
        A value of 1 means the motif is present at that time step (in any attribute of its subspace),
        and 0 means it is not.
    subspaces_ : np.ndarray of shape (n_patterns, n_channels)
        Boolean mask indicating, for each discovered motif, the attributes in which it occurs.
    motif_sets_ : list of tuples (np.ndarray of shape (n_channels,), np.ndarray of shape (n_occurrences, 2, n_attributes))
        The subspace motif sets, each consisting of the subspace mask and the start and end
        of every occurrence in each attribute of the subspace.

    Notes
    -----
    The `prediction_mask_` variable loses information with respect to the start and end of
    each motif on each attribute: it assumes that every motif occurs on the same interval
    across all attributes. SubTSMD, however, finds motifs that are not perfectly aligned
    across attributes. This information can be obtained through the `motif_sets_` variable,
    which contains the start and end of every subspace motif.
    """

    def __init__(
            self,
            univariate_motif_discovery,
            delta: float = 0.5,
            linkage: Literal["complete", "average"] = "average",
            size_constraint: int = 1,
            normalize: bool = True,
    ):
        self.univariate_motif_discovery = univariate_motif_discovery
        self.delta = delta
        self.linkage = linkage
        self.size_constraint = size_constraint
        self.normalize = normalize

    def fit(
            self,
            signal: np.ndarray,
            inclusion_constraint: list[int] | None = None,
            exclusion_constraint: list[int] | None = None,
            co_occurrence_constraint: list[list[int]] | None = None,
    ) -> "SubTSMD":
        """Fit SubTSMD

        Parameters
        ----------
        signal : numpy array of shape (n_samples, n_channels)
            The input samples (time series length, number of channels).
        inclusion_constraint : list of int, optional
            A list of attribute indices that must be included in each discovered motif.
        exclusion_constraint : list of int, optional
            A list of attribute indices in which no motif may appear.
        co_occurrence_constraint : list of list of int, optional
            A list of lists, where each inner list contains attribute indices that must co-occur in the same motif.

        Constraints with multiple attributes are discovered on the multivariate signal of those
        attributes, which requires a motif discovery method that supports multivariate time series.
        The constraints may not share attributes.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        assert len(signal.shape) == 2 and signal.shape[1] > 1, "Input signal must be multivariate (2D array)."

        self.signal = signal
        (self.n, self.d) = self.signal.shape

        # Mine motifs in the attributes independently
        independent_motif_sets = []
        all_masks = self._get_motif_discovery_masks(
            inclusion_constraint=inclusion_constraint,
            exclusion_constraint=exclusion_constraint,
            co_occurrence_constraint=co_occurrence_constraint
        )
        for mask in all_masks:
            sub_signal = signal[:, mask]
            if self.normalize:
                sub_signal = _z_normalize(sub_signal)
            if mask.sum() == 1:
                sub_signal = sub_signal[:, 0]

            # Discover the motif sets
            motif_discovery = copy.deepcopy(self.univariate_motif_discovery)
            try:
                motif_discovery.fit(sub_signal)
            except Exception as e:
                if mask.sum() == 1:
                    raise
                raise ValueError(
                    f"{motif_discovery.__class__.__name__} failed on the multivariate signal of attributes "
                    f"{np.where(mask)[0].tolist()}. A co-occurrence constraint or an inclusion constraint with "
                    f"multiple attributes requires a motif discovery method that supports multivariate time "
                    f"series (e.g., LocoMotif)."
                ) from e
            motif_sets = [_mask_to_intervals(line) for line in motif_discovery.prediction_mask_]

            # Format the motif sets
            independent_motif_sets.append(
                [
                    (
                        mask,
                        np.repeat(
                            motif_set.reshape(len(motif_set), 2, 1),
                            mask.sum(),
                            axis=-1,
                        ),
                    )
                    for motif_set in motif_sets
                    if len(motif_set) >= 2
                ]
            )

        # Apply SubTSMD
        inclusion_constraint_set = (
            inclusion_constraint is not None and len(inclusion_constraint) >= 1
        )
        subspace_motif_sets = apply_sub_tsmd(
            independent_motif_sets,
            delta=self.delta,
            linkage=self.linkage,
            inclusion_constraint_set=inclusion_constraint_set,
        )

        # Store the motif sets that satisfy the size constraint
        self.motif_sets_ = [
            motif_set
            for motif_set in subspace_motif_sets
            if motif_set[0].sum() >= self.size_constraint
        ]

        return self

    def _get_motif_discovery_masks(
            self,
            inclusion_constraint: list[int] | None = None,
            exclusion_constraint: list[int] | None = None,
            co_occurrence_constraint: list[list[int]] | None = None,
    ) -> list[np.ndarray]:
        attributes = np.arange(self.d)
        used = np.zeros(self.d, dtype=bool)
        masks = []

        # The constraints may not share attributes, as overlapping masks cannot be merged
        groups = [inclusion_constraint or [], exclusion_constraint or []] + list(co_occurrence_constraint or [])
        all_constrained = [attribute for group in groups for attribute in group]
        if len(all_constrained) != len(set(all_constrained)):
            raise ValueError(
                "The inclusion, exclusion and co-occurrence constraints may not share attributes, "
                f"but got inclusion={inclusion_constraint}, exclusion={exclusion_constraint}, "
                f"co-occurrence={co_occurrence_constraint}."
            )

        # The inclusion mask must come first, as apply_sub_tsmd only merges into the first motif sets
        if inclusion_constraint:
            masks.append(np.isin(attributes, inclusion_constraint))
            used[inclusion_constraint] = True

        if exclusion_constraint:
            used[exclusion_constraint] = True

        for constraint in co_occurrence_constraint or []:
            masks.append(np.isin(attributes, constraint))
            used[constraint] = True

        # Remaining attributes each get their own univariate mask
        masks.extend(np.eye(self.d, dtype=bool)[~used])

        return masks

    @property
    def subspaces_(self) -> np.ndarray:
        return np.array([subspace for subspace, _ in self.motif_sets_], dtype=bool).reshape(-1, self.d)

    @property
    def prediction_mask_(self) -> np.ndarray:
        mask = np.zeros((len(self.motif_sets_), self.n))
        for i, (_, motifs) in enumerate(self.motif_sets_):
            for motif in motifs:
                mask[i, int(motif[0].min()):int(motif[1].max())] = 1
        return mask


def _z_normalize(X: np.ndarray) -> np.ndarray:
    std = X.std(axis=0)
    std[std == 0] = 1
    return (X - X.mean(axis=0)) / std


def _mask_to_intervals(line: np.ndarray) -> np.ndarray:
    # Convert a binary mask to an array of shape (n_occurrences, 2) with [start, end) intervals
    diff = np.diff(np.concatenate(([0], (line > 0).astype(int), [0])))
    return np.stack((np.where(diff == 1)[0], np.where(diff == -1)[0]), axis=1)


def apply_sub_tsmd(
    independent_motif_sets,
    delta: float = 0.5,
    linkage: Literal["complete", "average"] = "average",
    inclusion_constraint_set: bool = False,
):
    # No motif can satisfy the inclusion constraint if the first (inclusion) motif sets are empty
    if inclusion_constraint_set and (len(independent_motif_sets) == 0 or len(independent_motif_sets[0]) == 0):
        return []

    # Filter empty motif sets
    independent_motif_sets = list(filter(lambda s: len(s) > 0, independent_motif_sets))

    # Nothing to do if there are no motif sets
    if len(independent_motif_sets) == 0:
        return []

    # Set the type of the motifs, as required by the signature of _match
    independent_motif_sets = [
        [_as_motif_set(mask, locations) for mask, locations in motif_sets]
        for motif_sets in independent_motif_sets
    ]

    def _merge_safe(motif_sets, motif_sets_prime):
        # Numba cannot type empty lists, so handle them in python
        if len(motif_sets) == 0:
            return [] if inclusion_constraint_set else motif_sets_prime
        if len(motif_sets_prime) == 0:
            return motif_sets
        return _merge(
            motif_sets,
            motif_sets_prime,
            delta=delta,
            linkage=linkage,
            inclusion_constraint_set=inclusion_constraint_set,
        )

    # Iteratively merge the motif sets
    return functools.reduce(_merge_safe, independent_motif_sets)


def _as_motif_set(mask: np.ndarray, locations: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    # Copy into writeable, C-contiguous arrays of the types in type_motif_set
    return (
        np.array(mask, dtype=np.bool_, order="C", copy=True),
        np.array(locations, dtype=np.float32, order="C", copy=True),
    )


@nb.njit(
    nb.float32(nb.float32, nb.float32, nb.float32, nb.float32),
    fastmath=True,
    cache=True,
)
def _overlap_rate(start_1, end_1, start_2, end_2):
    return max(
        0,
        (min(end_1, end_2) - max(start_1, start_2))
        / (max(end_1, end_2) - min(start_1, start_2)),
    )


@nb.njit(
    nb.types.Tuple((type_motif_set, type_motif_set, type_motif_set))(
        type_motif_set, type_motif_set, nb.float32, nb.types.unicode_type
    ),
    fastmath=True,
    cache=True,
)
def _match(motif_set, motif_set_prime, delta, linkage):

    matches = np.full(shape=motif_set[1].shape[0], fill_value=-1)
    matches_prime = np.full(shape=motif_set_prime[1].shape[0], fill_value=-1)

    for i in range(motif_set[1].shape[0]):
        for j in range(motif_set_prime[1].shape[0]):

            # We can skip this due to the property
            if matches_prime[j] != -1:
                continue

            # Between the overlap rates across each attribute
            overlap_rates = np.ones(
                shape=(motif_set[1][i].shape[1], motif_set_prime[1][j].shape[1])
            )
            for a in range(motif_set[1][i].shape[1]):
                for b in range(motif_set_prime[1][j].shape[1]):
                    overlap_rates[a, b] = _overlap_rate(
                        start_1=motif_set[1][i, 0, a],
                        end_1=motif_set[1][i, 1, a],
                        start_2=motif_set_prime[1][j, 0, b],
                        end_2=motif_set_prime[1][j, 1, b],
                    )

            # The overlap rate between subspace motifs
            if linkage == "complete":
                overlap = (
                    overlap_rates.min()
                )  # min instead of max because overlap rate is a similarity and not a distance
            elif linkage == "average":
                overlap = overlap_rates.mean()
            else:
                raise ValueError(
                    f"Unknown linkage method: '{linkage}'. Expected one of 'complete' or 'average'."
                )

            # Check if there is a match
            if overlap > delta or (
                delta == 1 and overlap == 1
            ):  # '>' and not '>=', to ensure single motif overlaps
                matches[i] = j
                matches_prime[j] = i
                break  # We can break the inner for-loop due to the property

    # Combine the matched motifs
    matched_mask = motif_set[0] | motif_set_prime[0]
    matched_motifs = np.zeros(
        shape=(np.sum(matches_prime != -1), 2, matched_mask.shape[0]), dtype=nb.float32
    )
    matched_motifs[:, :, motif_set[0]] = motif_set[1][matches != -1]
    matched_motifs[:, :, motif_set_prime[0]] = motif_set_prime[1][
        matches[matches != -1]
    ]
    matched_motifs = matched_motifs[:, :, matched_mask]

    # Extract the motifs that were not matched
    unmatched_motifs = (motif_set[0], motif_set[1][matches == -1])
    unmatched_motifs_prime = (
        motif_set_prime[0],
        motif_set_prime[1][matches_prime == -1],
    )

    # Return all the motifs
    return (matched_mask, matched_motifs), unmatched_motifs, unmatched_motifs_prime


@nb.njit(fastmath=True, cache=True)
def _merge(motif_sets, motif_sets_prime, delta, linkage, inclusion_constraint_set):

    merged_motif_sets = []

    for i in range(len(motif_sets)):
        for j in range(len(motif_sets_prime)):

            # Match the motifs (and also obtain the unmatched motifs)
            matched_motifs, unmatched_motifs, unmatched_motifs_prime = _match(
                motif_sets[i], motif_sets_prime[j], delta=delta, linkage=linkage
            )

            if matched_motifs[1].shape[0] >= 2:
                # Add the matched motif sets
                merged_motif_sets.append(matched_motifs)

                # Update the motif sets that have not been matched
                motif_sets[i] = unmatched_motifs
                motif_sets_prime[j] = unmatched_motifs_prime

                # Discard if only one motif set remains
                if unmatched_motifs[1].shape[0] <= 1:
                    motif_sets[i] = (
                        np.zeros(shape=matched_motifs[0].shape[0], dtype=nb.boolean),
                        np.empty(shape=(0, 2, 0), dtype=nb.float32),
                    )
                if unmatched_motifs_prime[1].shape[0] <= 1:
                    motif_sets_prime[j] = (
                        np.zeros(shape=matched_motifs[0].shape[0], dtype=nb.boolean),
                        np.empty(shape=(0, 2, 0), dtype=nb.float32),
                    )

    # Add all the motif sets
    merged_motif_sets += [
        motif_set for motif_set in motif_sets if motif_set[1].shape[0] > 0
    ]
    if not inclusion_constraint_set:
        # Only include the second list if there was no inclusion constraint set.
        merged_motif_sets += [
            motif_set_prime
            for motif_set_prime in motif_sets_prime
            if motif_set_prime[1].shape[0] > 0
        ]

    # Return all the motif sets
    return merged_motif_sets
