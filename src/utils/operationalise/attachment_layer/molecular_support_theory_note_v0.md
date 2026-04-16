# Molecular Support Theory Note v0

Status: working doctrine for the Hessian QM9 support/nullspace layer.

## Why The First Atlas Projected The Hessian

The first molecular reader used `visible_precision`. That surface requires a strictly positive definite finite object.

A raw Cartesian molecular Hessian at an equilibrium is not strictly positive definite in the full Cartesian coordinate space. It contains rigid translation and rotation directions. In ideal nonlinear molecules these are six zero modes; in linear molecules the rotational rank can differ. In numerical finite-difference Hessians, those directions may appear as small positive or negative residual directions.

So the first atlas used the honest SPD route:

`Cartesian Hessian -> mass weighting -> rigid-mode removal -> vibrational SPD block -> transformed observer -> visible_precision`

This is not a loss of interest. It is a route discipline. The visible-precision atlas reads the vibrational support, while the removed rigid/null/near-null structure belongs to a separate support provenance layer.

## What The Support Layer Can Use

The support layer is allowed to record:

- rigid translation/rotation rank;
- raw and projected inertia;
- rigid-block residual;
- rigid-vibrational cross residual;
- lowest positive vibrational eigenvalues;
- negative projected vibrational eigenvalues;
- observer rank loss after projection;
- soft-mode participation by atom group or element;
- environment-dependent support status across separately optimised Hessians.

These are not hidden-load values. They are support-stratum, compiler-provenance, and near-boundary diagnostics.

## What The Support Layer Must Not Claim Yet

The support layer must not claim:

- hidden load without a declared ceiling/reference;
- a birth/death event from two unrelated static records;
- a continuous solvent or reaction path without declared path data;
- a chemical transition-state classification from a negative mode alone;
- a theorem-level support clock without the appropriate support-event coefficients or path derivative data.

The module has support-event machinery, but those public surfaces require explicit support coordinates, finite-order event jets, or path derivative data. A static Hessian database can suggest support-boundary candidates; it does not by itself provide a support-event theorem instance.

## Current Interpretation

The molecular support atlas turns the information that was projected away in the first atlas into a second scientific surface.

The clean division is:

- vibrational SPD atlas: what a declared atom observer sees inside the positive vibrational support;
- support atlas: how the record sits relative to the boundary of that positive support.

This is potentially valuable because near-boundary molecular curvature, soft modes, and environment-specific projected indefiniteness are not noise by default. They may be numerical artifacts, but they may also mark flexible, weakly stabilised, or environment-sensitive structures. The support layer lets us keep those cases as disciplined leads rather than discarding them.

## Current Gate

The full four-environment Hessian QM9 support pass processed `166,580` records. It emitted support rows for all records with zero parser failures.

The first strong pattern is environment-dependent support admission:

- vacuum SPD rate: `0.9745`;
- THF SPD rate: `0.9699`;
- toluene SPD rate: `0.9909`;
- water SPD rate: `0.9560`.

This is an empirical support-gate pattern, not yet a solvent-physics theorem. It is worth investigating because it is large, systematic, and aligned with soft-mode drift.
