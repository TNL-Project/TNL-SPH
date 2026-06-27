"""
Vectorized normalization of surface normal arrays.

Replaces the per-particle Python loop found in several ``init.py`` scripts
with a vectorized numpy operation (~100x faster for large particle sets).
"""

import numpy as np


def normalize_vectors(vectors, tol=1e-12):
   """
   Normalize each row vector to unit length.

   Parameters
   ----------
   vectors : ndarray, shape (n, 3)
      Input vectors (one per row).
   tol : float
      Vectors with magnitude below this threshold are left unchanged
      and counted as zero-length.

   Returns
   -------
   normalized : ndarray, shape (n, 3), dtype float64
      Unit-length vectors.  Zero-length inputs are left as-is (not NaN).
   zero_count : int
      Number of vectors whose magnitude was below ``tol``.
   """
   vectors = np.asarray(vectors, dtype=float)
   magnitudes = np.sqrt(np.sum(vectors**2, axis=1))
   zero_mask = magnitudes < tol
   zero_count = int(np.count_nonzero(zero_mask))

   # avoid division by zero; zero-length vectors keep their original values
   safe_mags = np.where(zero_mask, 1.0, magnitudes)
   normalized = vectors / safe_mags[:, np.newaxis]

   return normalized, zero_count
