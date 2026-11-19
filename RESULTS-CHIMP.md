# RESULTS-CHIMP: do chimpanzee 'additive' call combinations mean A AND B, or A OR B?

Every line below is printed by render_md_lines in test_chimp_results.py from results/chimp/chimp_results.json (job chimp-bigram-v1), as registered in CHIMP-PREREG.md. MV2, Munim, Inc. Made with AI agents; no human expert has checked it yet.

Data: zenodo, calldata.csv sha256 7cfe322d68c88c5c56bd2be3abce94d6f86e446c2a85f2728b0c54033c0bdf0d
Gate D: 4,322 rows read, 4,322 used, 0 length mismatches; missing bigrams none; missing single calls none: **pass**
Gate R: our closest single call matches Table 1 for 14 of 16 bigrams (need 12): **pass**

| bigram | Table 1 closest | our closest | our second | distance | match |
|---|---|---|---|---|---|
| GR_BK | HO | SC | PS | 0.391 | no |
| GR_HO | HO | HO | PH | 0.296 | yes |
| GR_PG | PG | PG | PN | 0.197 | yes |
| GR_PN | PG | PN | PG | 0.301 | no |
| HO_GR | HO | HO | GR | 0.260 | yes |
| HO_PG | PG | PG | HO | 0.763 | yes |
| HO_PH | PH | PH | HO | 0.232 | yes |
| HO_PN | PG | PG | HO | 0.626 | yes |
| PG_GR | PG | PG | PN | 0.408 | yes |
| PG_PN | PN | PN | PG | 0.449 | yes |
| PH_HO | HO | HO | PH | 0.190 | yes |
| PH_PB | PH | PH | PB | 0.271 | yes |
| PH_PG | PG | PG | PN | 0.263 | yes |
| PH_PS | PH | PH | HO | 0.141 | yes |
| PH_SC | PH | PH | HO | 0.360 | yes |
| PN_PG | PG | PG | PN | 0.413 | yes |

## T1: conjunction versus best-fitting mixture

Family F1: 8 of 8 claimed additive bigrams testable; Holm at 0.05 over them; calibration level 0.0063

| bigram | family | n | markers of A | markers of B | k (both) | J | w | pi (mixture) | J - pi | p T1 | p T2 | label |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GR_BK | - | 21 | AP, FE | BA, EN, IG, IP | 0 | 0.000 | 0.627 | 0.033 | -0.033 | 1.0000 | 1.0000 | exploratory: NOT_REJECT |
| GR_HO | F1 | 26 | AP | RE, TR | 0 | 0.000 | 0.177 | 0.018 | -0.018 | 1.0000 | 1.0000 | `UNDERPOWERED` |
| GR_PG | - | 124 | FE | AF, NS, TR | 1 | 0.008 | 0.064 | 0.014 | -0.006 | 0.7324 | 0.7712 | exploratory: NOT_REJECT |
| GR_PN | - | 54 | FE | AF, PL | 0 | 0.000 | 0.095 | 0.007 | -0.007 | 1.0000 | 1.0000 | exploratory: NOT_REJECT |
| HO_GR | F1, strict | 53 | RE, TR | AP | 0 | 0.000 | 0.522 | 0.024 | -0.024 | 1.0000 | 1.0000 | `MIXTURE_NOT_REJECTED` |
| HO_PG | - | 62 | FE, RE | AF, AP, NS | 0 | 0.000 | 0.139 | 0.018 | -0.018 | 1.0000 | 1.0000 | exploratory: NOT_REJECT |
| HO_PH | F1, strict | 136 | AP | IP | 0 | 0.000 | 0.294 | 0.002 | -0.002 | 1.0000 | 1.0000 | `UNDERPOWERED` |
| HO_PN | - | 11 | FE, RE, TR | AF, AP, PL | 0 | 0.000 | 0.712 | 0.032 | -0.032 | 1.0000 | n/a | exploratory: NOT_REJECT |
| PG_GR | F1 | 27 | AF, NS, TR | FE | 1 | 0.037 | 0.693 | 0.011 | +0.027 | 0.2391 | 0.3665 | `MIXTURE_NOT_REJECTED` |
| PG_PN | F1, strict | 13 | NS, TR | PL | 0 | 0.000 | 0.322 | 0.003 | -0.003 | 1.0000 | 1.0000 | `UNDERPOWERED` |
| PH_HO | F1 | 23 | IP | AP | 1 | 0.043 | 0.095 | 0.002 | +0.041 | 0.0523 | 0.1072 | `UNDERPOWERED` |
| PH_PB | - | 69 | FE, TR | BA, IG, IP | 17 | 0.246 | 0.644 | 0.165 | +0.082 | 0.0735 | 0.3041 | exploratory: NOT_REJECT |
| PH_PG | F1 | 15 | FE, IP, RE, TR | AF, AP, NS | 2 | 0.133 | 0.030 | 0.081 | +0.052 | 0.3394 | 0.4489 | `MIXTURE_NOT_REJECTED` |
| PH_PS | - | 62 | FU, IP, TR | BE, RA | 0 | 0.000 | 1.000 | 0.002 | -0.002 | 1.0000 | 1.0000 | exploratory: NOT_REJECT |
| PH_SC | - | 17 | FE, FU, IP, RE, TR | GA, RA | 0 | 0.000 | 0.864 | 0.031 | -0.031 | 1.0000 | 1.0000 | exploratory: NOT_REJECT |
| PN_PG | F1, strict | 21 | PL | NS, TR | 0 | 0.000 | 0.219 | 0.003 | -0.003 | 1.0000 | 1.0000 | `UNDERPOWERED` |

HO_PH: **`UNDERPOWERED`** (C2 label power 0.00; C1 false-positive rate 0.00)
PN_PG: **`UNDERPOWERED`** (C2 label power 0.37; C1 false-positive rate 0.00)
HO_GR: **`MIXTURE_NOT_REJECTED`** (C2 label power 0.81; C1 false-positive rate 0.03)
PG_PN: **`UNDERPOWERED`** (C2 label power 0.14; C1 false-positive rate 0.00)
PH_HO: **`UNDERPOWERED`** (C2 label power 0.07; C1 false-positive rate 0.00)
GR_HO: **`UNDERPOWERED`** (C2 label power 0.61; C1 false-positive rate 0.01)
PG_GR: **`MIXTURE_NOT_REJECTED`** (C2 label power 1.00; C1 false-positive rate 0.01)
PH_PG: **`MIXTURE_NOT_REJECTED`** (C2 label power 1.00; C1 false-positive rate 0.04)

Overall label: **`NO_CONJUNCTION_POWERED`**
Strict subset (distance criterion): HO_PH UNDERPOWERED, PN_PG UNDERPOWERED, HO_GR MIXTURE_NOT_REJECTED, PG_PN UNDERPOWERED
T2 Holm (length-matched, Holm over F1): HO_PH no, PN_PG no, HO_GR no, PG_PN no, PH_HO no, GR_HO no, PG_GR no, PH_PG no
Exploratory, all 16 bigrams: 0 of 16 testable rejected (none)

## Calibration (C1 mixture, C2 conjunction, C4 distance criterion)

| bigram | level | C1 nominal FPR | C1 label rate | C2 label power | C4 mixture: A, B closest | C4 conjunction: A, B closest | real: A, B closest |
|---|---|---|---|---|---|---|---|
| HO_PH | 0.0063 | 0.000 | 0.000 | 0.000 | 200 of 200 | 200 of 200 | yes |
| PN_PG | 0.0063 | 0.000 | 0.000 | 0.370 | 193 of 200 | 200 of 200 | yes |
| HO_GR | 0.0063 | 0.030 | 0.000 | 0.815 | 190 of 200 | 199 of 200 | yes |
| PG_PN | 0.0063 | 0.000 | 0.000 | 0.145 | 190 of 200 | 199 of 200 | yes |
| PH_HO | 0.0063 | 0.000 | 0.000 | 0.065 | 157 of 200 | 199 of 200 | yes |
| GR_HO | 0.0063 | 0.010 | 0.000 | 0.615 | 101 of 200 | 194 of 200 | no |
| PG_GR | 0.0063 | 0.015 | 0.000 | 1.000 | 99 of 200 | 193 of 200 | no |
| PH_PG | 0.0063 | 0.035 | 0.000 | 1.000 | 2 of 200 | 100 of 200 | no |

C4: A and B are the two closest single calls in 1132 of 1600 mixture sets and 1484 of 1600 conjunction sets; real bigrams 5 of 8
D1: 0.354 of 734 bigram utterances carry 2 or more events, versus 0.188 of 3,588 single-call utterances

## Sensitivity

S1_dedup (markers p >= 0.10, 2.0x; R 2,000): GR_HO NOT_SUPPORTED, HO_GR NOT_SUPPORTED, HO_PH NOT_SUPPORTED, PG_GR NOT_SUPPORTED, PG_PN NOT_SUPPORTED, PH_HO NOT_SUPPORTED, PH_PG NOT_SUPPORTED, PN_PG NOT_SUPPORTED
S2_thresholds_1 (markers p >= 0.05, 1.5x; R 2,000): GR_HO NOT_SUPPORTED, HO_GR NOT_SUPPORTED, HO_PH NOT_SUPPORTED, PG_GR NOT_SUPPORTED, PG_PN NOT_SUPPORTED, PH_HO NOT_SUPPORTED, PH_PG NOT_SUPPORTED, PN_PG NOT_SUPPORTED
S2_thresholds_2 (markers p >= 0.20, 3.0x; R 2,000): GR_HO NOT_SUPPORTED, HO_GR NOT_SUPPORTED, HO_PH NOT_TESTABLE, PG_GR NOT_SUPPORTED, PG_PN NOT_SUPPORTED, PH_HO NOT_TESTABLE, PH_PG NOT_SUPPORTED, PN_PG NOT_SUPPORTED

## Predictions (CHIMP-PREREG.md section 7)

| | Prediction | Prior | Outcome |
|---|---|---|---|
| P1 | Gate D passes | 0.90 | **true** |
| P2 | Gate R passes (closest single matches Table 1 for >= 12 of 16) | 0.70 | **true** |
| P3 | At least 4 of the 8 F1 bigrams are testable | 0.50 | **true** |
| P4 | HO_GR (the paper's lead additive example) is testable | 0.30 | **true** |
| P5 | Overall label is CONJUNCTION_FOUND | 0.35 | **false** |
| P6 | C1 nominal false-positive rate <= 0.10 for every testable F1 bigram | 0.80 | **true** |
| P7 | Median C2 label power over testable F1 bigrams >= 0.80 | 0.50 | **false** |
| P8 | C4: mixture sets put A and B as the two closest singles at least half as often as conjunction sets | 0.60 | **true** |
| P9 | D1: bigrams carry 2+ events more often than single calls | 0.65 | **true** |
| P10 | Every CONJUNCTION_SUPPORTED bigram also has T2 p < 0.05 (graded only if one exists) | 0.60 | **not graded** |

7 of 9 graded predictions true
Brier score **0.171**

## Reading (written by hand; not printed by render_md_lines)

Everything in this section is a plain-words reading of the numbers above. The rules come from
CHIMP-PREREG.md sections 4 to 7 and were applied exactly as coded.

**What the pre-registered test says.** Girard-Buttoz et al. (2025) call 8 two-call combinations
"additive": the pair is said to keep both calls' meanings and to let a chimpanzee signal two
concurrent events in one utterance. All 8 were testable. Holm rejected the A-or-B mixture for
none of them (smallest p 0.0523, for PH_HO, against a first Holm threshold of 0.05 / 8). The
overall label is `NO_CONJUNCTION_POWERED`: for the three pairs where the test had at least 80%
power to see a true conjunction (HO_GR 0.815, PG_GR 1.000, PH_PG 1.000), the data look like a
mixture of A's and B's situations. For the other five, the test had too little power (0.000 to
0.615), so it says nothing either way.

**The count behind it.** Across the 314 utterances of the 8 "additive" pairs, only 4 carry at
least one event typical of call A together with at least one event typical of call B
(PG_GR 1 of 27, PH_HO 1 of 23, PH_PG 2 of 15, and 0 for the other five). The paper's lead
example, HO_GR, has 0 of 53 such utterances; under a pure conjunction built from real HO and GR
utterances the same rule would have labelled it a conjunction in 81.5% of synthetic sets.

**Why the published criterion could not tell.** The paper's criterion asks whether A and B are
the two single calls whose average event profiles lie closest to the pair's profile. In our
synthetic sets this was true in 1132 of 1600 pure-mixture sets (A OR B, no added meaning) and
1484 of 1600 pure-conjunction sets (A AND B). So the criterion is satisfied by mixtures most of
the time; it is a test of "between A and B", not of "A and B at once" (P8 true).

**What we did not expect.** Pairs do carry two or more events more often than single calls
(0.354 of 734 bigram utterances versus 0.188 of 3,588 single calls; P9 true), but in the
"additive" pairs those extra events are not A's events plus B's events. The pair with the most
A-plus-B utterances, PH_PB (17 of 69), is not one the paper calls additive, and it too stays
within the mixture null in the exploratory test (p 0.0735, not rejected). The two false
predictions are P5 (a conjunction found; prior 0.35) and P7 (median power at least 0.80; prior
0.50; it was 0.4925, the mean of 0.370 and 0.615). P4 came true against a prior of only 0.30.

**Limits.** Utterances from the same caller and day are not independent; after collapsing
identical rows (S1) no pair is supported as a conjunction. Marker events are defined by our
thresholds; under two other threshold settings (S2) no pair is supported either, and HO_PH and
PH_HO become untestable at the strictest one.
Our event vectors are raw shares, not the paper's model-based vectors; they still pick the
paper's closest single call for 14 of 16 pairs (Gate R). The test addresses only the paper's
"additive" class and only its per-utterance reading ("concurrently occurring events within a
single utterance"); it says nothing about the other three mechanisms, about acoustics, or about
what the calls mean to listeners. Five of the eight pairs are underpowered with these data. The
event coding was made by the original observers; we used it as released. This work was carried
out with AI agents working for MV2, Munim, Inc.; no human expert has checked it yet.
