# Replication Fragility As Observer Mismatch

Run:

```bash
python -m worked_examples.replication_fragility.run_main
```

This worked example uses one latent Gaussian object and two different protocol
observers. The point is not sociology. The point is that observer mismatch
alone can generate visibly different protocol outputs.

What it shows:

- one latent precision can project to visibly different protocol laws
- the difference can be quantified with an exact Gaussian mismatch score
- a finite refinement search can identify a richer observer that reconciles the
  two protocols
