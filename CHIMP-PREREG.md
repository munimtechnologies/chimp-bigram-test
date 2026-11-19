# CHIMP-PREREG: do chimpanzee "additive" call combinations mean A AND B, or A OR B?

Pre-registration, written 2026-11-19 before any statistic was computed on the data.
MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents; no human expert has reviewed it yet.
Text CC BY 4.0; code Apache-2.0.

## 1. The claim and the gap

Girard-Buttoz, Neumann, Bortolato, Zaccarella, Friederici, Wittig & Crockford (2025),
"Versatile use of chimpanzee call combinations promotes meaning expansion", *Science Advances*
11:eadq2879, recorded 4,323 utterances (3,589 single calls, 734 two-call "bigrams") from 53
wild chimpanzees in Taï National Park, each tagged with one or more of 22 co-occurring events.
They report four ways bigrams change meaning. One is **additive** (mechanism 3, "A + B = A + B"):
the bigram "retains the meaning of each of the units composing it and those meanings are added
together", and it "enables communication about concurrently occurring events within a single
utterance". Their example: A = feeding, B = resting, AB = feeding + resting.

That is a **conjunction**: AB should be used when A's situation and B's situation hold at the
same time. But the test used for it compares only the *average* event profile of AB with those
of the single calls (Euclidean distance; "P.support.pair" = posterior share in which A and B
are AB's two closest units). A bigram that is used sometimes where A is used and sometimes where
B is used (a **mixture**, A OR B, with no added meaning) has an average profile between A and B
too, so it can earn the same high pair support. The per-utterance event coding, released with
the paper, can tell the two apart: under conjunction, single utterances of AB carry A's events
and B's events together; under a mixture they rarely do, no more often than A's or B's own
utterances.

Prior art (searched 2026-11-19, Europe PMC full text and web): the paper reports no such test
("we do not strictly demonstrate that bigrams are compositional"); Klinedinst (2025, *Animal
Cognition* 28:90) discusses its subset/intersection readings and calls for "finer-grained
analyses" but runs none; Wartel et al. (2026, *PeerJ* 14:e21651) re-test the bonobo study
(Berthet et al. 2025) and explicitly set this chimpanzee study aside. We found no per-utterance
conjunction-versus-mixture test of any animal call combination.

## 2. Reading gate (what we saw before registering)

Seen: the paper's text (Europe PMC PMC12063654), its Table 1 classifications, the GitHub
README and a summary of the R/Stan code, the Zenodo record metadata (13889876, CC BY 4.0,
md5 96e88c4dff3ceb7f647ad5bd4ca213c0), and the CSV header plus its first 8 rows and the list of
call and event codes (to write the parser). No count, rate, marker set or statistic of the data
was computed by us before this file was public. The job waits until this file is public on
GitHub before it downloads anything.

## 3. Data

`calldata.csv` from the Zenodo zip (md5 checked; fallback: GitHub tag `chimpbigrams`, sha256
recorded). One row per utterance: date, caller, call.name (single "HO" or bigram "HO_GR",
first call first), context (event codes joined by "_"), context.length. A row is used only if its
event codes are distinct and as many as context.length.

**Gate D** (else stop, status gate_failed): rows read in 4,300-4,350; length mismatches at most
1%; all 16 bigrams of Table 1 and every component single call have at least 10 utterances.

**Gate R** (replication check, reported, does not stop T1): with our raw vectors (share of
utterances carrying each event), the single call closest to each bigram matches the paper's
Table 1 for at least 12 of 16 bigrams. The paper's vectors come from a Bayesian model for an
average caller, so exact agreement is not expected.

## 4. Test T1 (primary): conjunction versus best-fitting mixture

For a bigram AB with single calls A and B, and p_X(e) the share of X's utterances carrying
event e:

- **Markers.** M_A = events with p_A(e) >= 0.10 and p_A(e) >= 2 p_B(e); M_B the same with A and B
  swapped (disjoint by construction). A bigram is **testable** if M_A and M_B are both non-empty
  and AB, A and B each have at least 10 utterances. Otherwise it is NOT_TESTABLE: then A and B
  have no distinct situations for an "A + B" reading to add (itself a reportable outcome).
- **Statistic.** k = number of AB utterances whose event set contains at least one M_A event
  and at least one M_B event; J = k / n.
- **Mixture null.** Each AB utterance is a draw of a real A utterance (probability w) or a real
  B utterance (1 - w). w is the least-squares fit of AB's event profile by w p_A + (1 - w) p_B,
  clipped to [0, 1]. Under it K ~ Binomial(n, pi), pi = w q_A + (1 - w) q_B, where q_X =
  (k_X + 0.5) / (n_X + 1) is the smoothed share of X's own utterances that carry both marker kinds.
- **p value.** P(K >= k) averaged over 10,000 bootstrap resamples of the A and B utterance pools
  (w, q_A, q_B refitted in each), one-sided.
- **Family F1** = the 8 bigrams the paper calls additive by either criterion: HO_PH, PN_PG,
  HO_GR, PG_PN, PH_HO, GR_HO, PG_GR, PH_PG (strict subset by distances: HO_PH, PN_PG, HO_GR,
  PG_PN, reported separately). Holm at 0.05 over the testable bigrams of F1.

**Labels** (per bigram of F1, in this order):
NOT_TESTABLE; NULL_MISCALIBRATED if C1's nominal false-positive rate > 0.10;
CONJUNCTION_SUPPORTED if Holm rejects, J - pi_hat >= 0.10 and k >= 3; SMALL_EXCESS if Holm
rejects but not those two; MIXTURE_NOT_REJECTED if not rejected and C2 label power >= 0.80;
UNDERPOWERED otherwise.

**Overall:** CONJUNCTION_FOUND if any F1 bigram is CONJUNCTION_SUPPORTED;
NO_CONJUNCTION_POWERED if at least 2 are MIXTURE_NOT_REJECTED; NO_TESTABLE_BIGRAMS if all are
NOT_TESTABLE; else UNDECIDED.

## 5. Calibration (per testable F1 bigram; 200 synthetic sets each, same n, 500 bootstraps)

- **C1, negative control (pure mixture):** each synthetic utterance is a random real A utterance
  with probability w_hat, else a random real B utterance. Reports the nominal rate p <= 0.05 and
  the label rate.
- **C2, positive control (pure conjunction):** each synthetic utterance is the union of the event
  sets of a random real A utterance and a random real B utterance. Label power = share of sets
  meeting the CONJUNCTION_SUPPORTED rule at level 0.05 / (number of testable F1 bigrams).
- **C4, does the paper's distance criterion tell them apart?** For every synthetic set, whether
  A and B are the two singles closest (Euclidean, raw vectors) to the set's profile; counts are
  pooled over F1 for mixture and conjunction sets.

## 6. Secondary and sensitivity (reported, not used for the overall label)

- **T2 (length-matched):** the same test restricted to AB utterances with at least 2 events,
  against pools restricted the same way (q_X over X's multi-event utterances; w fitted on all AB
  utterances as in T1; same bootstrap); Holm over F1;
  needs at least 5 such AB utterances.
- **All 16 bigrams:** T1 with Holm over all testable bigrams (exploratory).
- **S1:** T1 after collapsing identical rows (same date, caller, call and event set), 2,000
  bootstraps. **S2:** marker thresholds (0.05, 1.5x) and (0.20, 3x).
- **D1:** share of utterances with at least 2 events, all bigrams versus all single calls.

Known limits, stated now: utterances from the same caller and day are not independent (S1
addresses exact repeats only); the mixture null draws from pooled single-call utterances of all
callers; marker thresholds are ours (S2 varies them); our event vectors are raw shares, not the
paper's model-based vectors (Gate R checks how close they are).

## 7. Predictions (probabilities for the Brier score)

| | Prediction | Prior |
|---|---|---|
| P1 | Gate D passes | 0.90 |
| P2 | Gate R passes (closest single matches Table 1 for >= 12 of 16) | 0.70 |
| P3 | At least 4 of the 8 F1 bigrams are testable | 0.50 |
| P4 | HO_GR (the paper's lead additive example) is testable | 0.30 |
| P5 | Overall label is CONJUNCTION_FOUND | 0.35 |
| P6 | C1 nominal false-positive rate <= 0.10 for every testable F1 bigram | 0.80 |
| P7 | Median C2 label power over testable F1 bigrams >= 0.80 | 0.50 |
| P8 | C4: mixture sets put A and B as the two closest singles at least half as often as conjunction sets | 0.60 |
| P9 | D1: bigrams carry 2+ events more often than single calls | 0.65 |
| P10 | Every CONJUNCTION_SUPPORTED bigram also has T2 p < 0.05 (graded only if one exists) | 0.60 |

The job (chimp.py, tested by test_chimp.py) computes every label, prediction and the Brier score
itself and writes them to chimp_results.json; RESULTS-CHIMP.md will print them unchanged, with
power next to every non-rejection, whichever way they fall. Any change after this file is
public goes in a dated CHIMP-AMENDMENT file; this file is never edited.

What would be remarkable either way: a clear CONJUNCTION_FOUND would be the first per-utterance
evidence that a wild ape call combination carries two concurrent meanings at once; a powered
NO_CONJUNCTION result would show the published "additive" class is indistinguishable from
calls used in either part's situations.
