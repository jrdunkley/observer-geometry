# Micro Case Studies

- `2604_06078_observability`
  - We took the paper's explicit partial-observation examples and re-encoded them as observer relations.
  - Result: exact coarsening structure; one case stays underdetermined instead of collapsing to a guessed source.
- `2604_06108_amp_protocol`
  - We took an exact ACS/WFC table slice and tested shared-observer versus richer-blue observer hypotheses.
  - Result: the richer blue observer reduces contradiction sharply, but the slice still does not close to common descent.
- `2604_06064_spmi_degeneracy`
  - We took the published summary table and split power-only, power-plus-escape, and power-plus-geometry observers.
  - Result: power-only is a coarsening of both richer observers; the richer pair is non-nested and underdetermined.
- `2604_06099_robustness_views`
  - We took exact benchmark tables and separated realistic-shift from adversarial-stress observers.
  - Result: they are non-nested; the table does not justify one scalar called robustness.

If you mention one of these case studies, cite both this workspace and the
source paper. See [CITATIONS.md](../../CITATIONS.md).

Run one case:

```powershell
python tests\micro_case_studies\<case>\run_main.py
```

