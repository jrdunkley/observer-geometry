# Free Gaussian RG As Quotient Descent

Run:

```bash
python -m worked_examples.free_gaussian_rg.run_main
```

This worked example stays inside the exact free Gaussian regime. It does not
claim an interacting RG theorem.

What it shows:

- mode retention is an exact quotient descent on the precision object
- staged mode elimination agrees with direct elimination by the tower law
- block averaging is a different coarse observer and is generally non-nested
  with retained-mode observation
- the first boundary is exactly where interacting terms would enter, which this
  example does not model
