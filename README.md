# chimp-bigram-test

When wild chimpanzees join two calls into one utterance, does the pair mean **A and B at once**,
or is it just used **where A or B would be used**?

Girard-Buttoz et al. (2025, *Science Advances* 11:eadq2879) report that some chimpanzee call
pairs are "additive": the pair keeps both calls' meanings and "enables communication about
concurrently occurring events within a single utterance". The test they used compares only
average event profiles, which cannot separate a true conjunction (A's and B's situations in the
same utterance) from a mixture (sometimes A's situation, sometimes B's). This project runs the
per-utterance test on the paper's own public data (Zenodo 13889876, CC BY 4.0), calibrated with
synthetic mixtures and conjunctions built from the real calls.

## Result (graded 2026-11-19): [RESULTS-CHIMP.md](RESULTS-CHIMP.md)

Overall label **`NO_CONJUNCTION_POWERED`**. All 8 pairs the paper calls additive were testable;
the A-or-B mixture was rejected for none of them. Across their 314 utterances, only 4 carry an
event typical of call A together with an event typical of call B. For the three pairs where the
test had at least 80% power to see a true conjunction (HO_GR, PG_GR, PH_PG), the data look like a
mixture of A's and B's situations; the other five are underpowered. The paper's distance
criterion ("A and B are the pair's two closest single calls") was met by 1132 of 1600 synthetic
pure-mixture sets, so it cannot tell "A and B" from "A or B". 7 of 9 graded predictions came
true (Brier 0.171).

## Files

- [CHIMP-PREREG.md](CHIMP-PREREG.md): the pre-registered question, test, labels, calibration and
  predictions, made public before any statistic was computed.
- `chimp.py`: the analysis (Python standard library only). `test_chimp.py`: its tests.
- `test_chimp_results.py`: the grader, published before any result existed. It does not import
  `chimp.py`. With its own code it re-derives Gate D, Gate R, every Holm decision, every
  per-bigram label, the overall label, C4, the sensitivity labels, the ten predictions and the
  Brier score from `results/chimp/chimp_results.json`, checks every count, smoothed rate,
  mixture weight and calibration rate for consistency, and fixes the exact lines
  RESULTS-CHIMP.md must print (`render_md_lines`).
- `results/chimp/chimp_results.json`: the raw output of job chimp-bigram-v1 (2026-11-19).

Run it yourself:

```
python -m unittest -v test_chimp test_chimp_results
python chimp.py --data-dir data --out-dir out
```

Data: Girard-Buttoz, C., Neumann, C., Bortolato, T., Zaccarella, E., Friederici, A. D.,
Wittig, R. M. & Crockford, C. (2025), calldata.csv, https://zenodo.org/records/13889876
(CC BY 4.0). The data are downloaded at run time, not copied into this repository.

By MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents; no human expert has checked this
work yet. Code Apache-2.0; text CC BY 4.0.
