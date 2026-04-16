# Incoming Data

Place manual static datasets here before analysis.

Expected first folders:

```text
incoming_data/rc_decay_trace/
incoming_data/beer_lambert_calibration/
```

Each folder should contain:

```text
data.csv
metadata.json
```

Use the shapes in `../dataset_templates/`.

Do not treat a dataset as admitted because the CSV is present. Admission requires the metadata declarations in `../empirical_dataset_protocol_v0.md`, especially units, noise/uncertainty convention, sensor regime, and refusal conditions.
