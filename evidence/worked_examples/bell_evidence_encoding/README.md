# Bell Evidence Encoding

This example shows the narrow job of `evidence`.

- Full pairwise Gaussian law data is encoded as exact visible evidence.
- The normalized correlator summary is preserved separately as a structured
  table.
- The correlator summary is not promoted into the full-law `ProblemSpec`.
- `nomodescent` then reproduces the known compatibility split.

The point is not automatic discovery. The point is to prevent the encoding layer
from silently flattening full-law evidence into a coarser summary.

