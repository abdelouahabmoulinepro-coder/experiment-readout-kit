# Experiment Readout Kit

A small, auditable command-line tool for the final readout of a two-arm online experiment with a binary conversion outcome. It reports conversion rates, absolute and relative lift, a Newcombe score interval for the difference, a pooled two-proportion z test, and a sample-ratio mismatch (SRM) check.

## Business use

Product and growth teams need a reproducible answer to: did the treatment change conversion, and does the assignment data look trustworthy enough to interpret? The tool produces a compact Markdown report plus machine-readable JSON. An SRM failure blocks a positive or negative conclusion. A result whose interval crosses zero is labelled `inconclusive`.

This is a final analysis utility. It does not select metrics, monitor experiments, or authorize a product launch.

## Quick start

Python 3.10 or later; no third-party dependencies.

```bash
python experiment_readout.py demo --output demo-output
python -m unittest -v
```

For your own experiment:

```bash
python experiment_readout.py analyze assignments.csv --output readout
python experiment_readout.py analyze assignments.csv --output readout \
  --control holdout --treatment new-flow \
  --expected-control-share 0.5 --confidence 0.95 --srm-alpha 0.001
```

The output folder receives `readout.json` and `report.md`. `demo-output/` in this repository contains a deterministic example input and its generated readout.

## Input contract

UTF-8 CSV (UTF-8 BOM accepted) with one row per randomized unit:

| Column | Meaning |
|---|---|
| `unit_id` | Unique user/account/session identifier at the actual randomization unit |
| `variant` | Exactly one of the configured control or treatment labels |
| `converted` | Binary outcome encoded as `0` or `1` |

Duplicate units, missing IDs, unknown variants, and non-binary outcomes are rejected with a row-specific error. The input should contain the predeclared analysis population and a final outcome window; the tool cannot determine exposure or attribution correctness.

## Statistical choices

- Conversion rates use all assigned units in each arm as denominators (intention-to-treat form, provided the input population was constructed correctly).
- The difference interval is Newcombe's hybrid score interval built from separate Wilson intervals. It behaves better near zero conversions than the plain Wald interval.
- The reported p-value is a two-sided pooled two-proportion z test under equal conversion rates.
- The two-arm SRM check compares observed arm counts with the configured expected control share using a one-degree-of-freedom chi-square goodness-of-fit statistic. The default block threshold is `p < 0.001`.
- A decision of `positive` or `negative` requires the difference interval to exclude zero and the SRM check to pass. Otherwise it is `inconclusive` or `invalid_srm`.

The interval and p-value use related but not identical procedures, so they need not match exactly at a significance boundary. Treat both as evidence, not as a substitute for a pre-analysis plan.

## Reproducible demo

The fixed-seed generator creates 12,000 fictional assignments, 6,000 per arm. The observed sample has 572 conversions in control (9.53%) and 676 in treatment (11.27%). Its estimated absolute lift is **+1.73 percentage points**, with a 95% Newcombe interval of **[+0.64, +2.83] points**. The pooled z-test p-value is **0.00187**, the SRM check passes, and the tool returns `positive` under its stated rule. These are synthetic demo results, not evidence about a real product or audience.

## Limits and responsible use

The calculations assume independent randomized units, a binary metric chosen before examining results, a valid final observation window, and no repeated peeking or optional stopping. Clustering, multiple comparisons, covariate adjustment, sequential designs, noncompliance, and missing outcomes are not modelled. The SRM check is only a diagnostic for aggregate arm counts and cannot explain the cause of an imbalance. Review instrumentation, exposure, guardrail metrics, experiment duration, and practical effect size before making a product decision.

All included data are synthetic and generated locally with a fixed seed. No external or private data are included.

## License

MIT. See [LICENSE](LICENSE).
