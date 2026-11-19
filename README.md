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

- [CHIMP-PREREG.md](CHIMP-PREREG.md): the pre-registered question, test, labels, calibration and
  predictions, made public before any statistic was computed.
- `chimp.py`: the analysis (Python standard library only). `test_chimp.py`: its tests.
- Status: **registered; the job (chimp-bigram-v1) is waiting for compute.** Results will be
  published in RESULTS-CHIMP.md whichever way they fall.

Run it yourself:

```
python -m unittest -v test_chimp
python chimp.py --data-dir data --out-dir out
```

Data: Girard-Buttoz, C., Neumann, C., Bortolato, T., Zaccarella, E., Friederici, A. D.,
Wittig, R. M. & Crockford, C. (2025), calldata.csv, https://zenodo.org/records/13889876
(CC BY 4.0). The data are downloaded at run time, not copied into this repository.

By MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents; no human expert has checked this
work yet. Code Apache-2.0; text CC BY 4.0.
