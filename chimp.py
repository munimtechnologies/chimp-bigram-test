# SPDX-License-Identifier: Apache-2.0
# chimp.py (chimp-bigram-test): when wild chimpanzees combine two calls A and B into a bigram AB
# with "additive" meaning (Girard-Buttoz et al. 2025, Sci. Adv. 11:eadq2879), is AB used when A's
# situation and B's situation happen together in the same utterance (conjunction, A AND B), or is
# it used sometimes where A is used and sometimes where B is used (mixture, A OR B)? The paper's
# distance criterion compares only average event profiles, which cannot tell these apart; the
# per-utterance event coding can. Registered in CHIMP-PREREG.md before any statistic was computed.
# Data: calldata.csv, Zenodo record 13889876 (CC BY 4.0), the paper's own data release.
# Stdlib only. MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents. Apache-2.0.
import argparse
import csv
import hashlib
import io
import json
import math
import os
import random
import statistics
import sys
import time
import urllib.request
import zipfile
from collections import Counter

VERSION = "chimp-1"
SEED = 20261119
UA = {"User-Agent": "chimp-bigram-test (MV2, Munim, Inc.; mv2@munimtech.com)"}
ZENODO_URL = ("https://zenodo.org/api/records/13889876/files/"
              "tozbu/Chimpanzee_bigram_meaning-chimpbigrams.zip/content")
ZENODO_MD5 = "96e88c4dff3ceb7f647ad5bd4ca213c0"
GITHUB_URL = ("https://raw.githubusercontent.com/tozbu/Chimpanzee_bigram_meaning/"
              "chimpbigrams/calldata.csv")
PREREG_URL = ("https://raw.githubusercontent.com/munimtechnologies/chimp-bigram-test/main/"
              "CHIMP-PREREG.md")
PREREG_TAG = b"CHIMP-PREREG"

ALPHA = 0.05
EXCESS = 0.10        # J - pi_mix must be at least 10 percentage points for CONJUNCTION_SUPPORTED
MIN_K = 3            # and at least 3 bigram utterances must carry both marker kinds
MIN_N = 10           # tokens needed for a bigram and for each of its single calls
MIN_N2 = 5           # multi-event bigram tokens needed for the length-matched test T2
FPR_MAX = 0.10       # C1: mixture false-positive rate above this -> NULL_MISCALIBRATED
POWER_OK = 0.80      # C2: label power at or above this -> a non-rejection is "powered"

# Table 1 of Girard-Buttoz et al. 2025: the single unit closest to each bigram (Euclidean).
PAPER_CLOSEST = {
    "HO_PN": "PG", "GR_BK": "HO", "GR_PN": "PG", "HO_PG": "PG", "PH_SC": "PH", "PH_PB": "PH",
    "PH_HO": "HO", "HO_PH": "PH", "PN_PG": "PG", "HO_GR": "HO", "PG_PN": "PN", "GR_PG": "PG",
    "PH_PS": "PH", "GR_HO": "HO", "PG_GR": "PG", "PH_PG": "PG"}
# Bigrams the paper classifies as additive (mechanism 3, alone or with 2) by either criterion.
CLAIMED = ("HO_PH", "PN_PG", "HO_GR", "PG_PN", "PH_HO", "GR_HO", "PG_GR", "PH_PG")
# Additive by the Euclidean-distance criterion (Table 1, "Mech. (Euclid.)" contains 3).
STRICT = ("HO_PH", "PN_PG", "HO_GR", "PG_PN")

CFG = {"seed": SEED, "bigrams": sorted(PAPER_CLOSEST), "claimed": list(CLAIMED),
       "strict": list(STRICT), "paper_closest": PAPER_CLOSEST, "gate_r_need": 12,
       "rows_lo": 4300, "rows_hi": 4350, "min_p": 0.10, "ratio": 2.0,
       "sens_thresholds": [[0.05, 1.5], [0.20, 3.0]],
       "r_main": 10000, "r_cal": 500, "s_cal": 200, "r_sens": 2000}

PRIORS = {"P1": 0.90, "P2": 0.70, "P3": 0.50, "P4": 0.30, "P5": 0.35,
          "P6": 0.80, "P7": 0.50, "P8": 0.60, "P9": 0.65, "P10": 0.60}


# ---------------------------------------------------------------- data
def fetch(url, timeout=120):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def wait_for_registration(max_minutes, fetcher=None):
    fetcher = fetcher or fetch
    t0 = time.time()
    tries = 0
    last = None
    while True:
        tries += 1
        try:
            blob = fetcher(PREREG_URL, timeout=60)
            if PREREG_TAG in blob:
                return {"found": True, "tries": tries, "url": PREREG_URL,
                        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                        "sha256": hashlib.sha256(blob).hexdigest()}
        except Exception as e:                       # 404 until published
            last = repr(e)[:200]
        if time.time() - t0 > max_minutes * 60:
            return {"found": False, "tries": tries, "url": PREREG_URL, "last_error": last}
        time.sleep(60)


def get_data(data_dir, fetcher=None):
    """calldata.csv from the Zenodo zip (md5-checked, cached); GitHub tag as fallback."""
    fetcher = fetcher or fetch
    os.makedirs(data_dir, exist_ok=True)
    cache = os.path.join(data_dir, "chimpbigrams.zip")
    info = {}
    try:
        if os.path.exists(cache):
            with open(cache, "rb") as f:
                blob = f.read()
        else:
            blob = fetcher(ZENODO_URL, timeout=120)
        md5 = hashlib.md5(blob).hexdigest()
        if md5 != ZENODO_MD5:
            raise ValueError("md5 mismatch " + md5)
        if not os.path.exists(cache):
            with open(cache, "wb") as f:
                f.write(blob)
        with zipfile.ZipFile(io.BytesIO(blob)) as zf:
            names = sorted(n for n in zf.namelist()
                           if n.endswith("calldata.csv") and "/._" not in n
                           and not n.startswith("._"))
            raw = zf.read(names[0])
        info = {"source": "zenodo", "url": ZENODO_URL, "zip_md5": md5, "member": names[0]}
    except Exception as e:
        info = {"zenodo_error": repr(e)[:300]}
        raw = fetcher(GITHUB_URL, timeout=120)
        info.update({"source": "github-tag-chimpbigrams", "url": GITHUB_URL})
    info["csv_sha256"] = hashlib.sha256(raw).hexdigest()
    info["csv_bytes"] = len(raw)
    return raw.decode("utf-8-sig"), info


def parse_rows(text):
    """-> (rows, n_length_mismatch, n_read). A row is kept only if its context codes, split on
    '_', are distinct and as many as context.length."""
    BOM = chr(0xFEFF)
    text = text[1:] if text.startswith(BOM) else text
    # BOM removed above; the line below is kept only as a comment: ("﻿") else text
    lines = text.splitlines()  # was: lstrip("﻿").splitlines()
    if lines and lines[0].startswith('"') and lines[0].endswith('"') and lines[0].count('"') == 2:
        lines = [l[1:-1] if len(l) >= 2 and l.startswith('"') and l.endswith('"') else l
                 for l in lines]
    rows, mismatch, n = [], 0, 0
    for r in csv.DictReader(lines):
        n += 1
        ev = [e.strip() for e in (r.get("context") or "").split("_") if e.strip()]
        try:
            L = int(r.get("context.length") or -1)
        except ValueError:
            L = -1
        if not ev or len(set(ev)) != L:
            mismatch += 1
            continue
        rows.append({"date": (r.get("date") or "").strip(),
                     "caller": (r.get("caller") or "").strip(),
                     "call": (r.get("call.name") or "").strip(),
                     "events": frozenset(ev)})
    return rows, mismatch, n


def call_tokens(rows):
    tok = {}
    for r in rows:
        tok.setdefault(r["call"], []).append(r["events"])
    return tok


def all_events(tok):
    return sorted(set(e for ts in tok.values() for t in ts for e in t))


def dedup(rows):
    """S1: one row per (date, caller, call, event set)."""
    seen = set()
    out = []
    for r in rows:
        key = (r["date"], r["caller"], r["call"], r["events"])
        if key not in seen:
            seen.add(key)
            out.append(r)
    return out


def gate_d(rows, mismatch, n_read, cfg):
    tok = call_tokens(rows)
    missing = [b for b in cfg["bigrams"] if len(tok.get(b, [])) < MIN_N]
    comp = sorted(set(c for b in cfg["bigrams"] for c in b.split("_")
                      if len(tok.get(c, [])) < MIN_N))
    checks = {"rows_in_range": cfg["rows_lo"] <= n_read <= cfg["rows_hi"],
              "mismatch_at_most_1pct": mismatch <= 0.01 * n_read,
              "all_bigrams_present": not missing,
              "all_components_present": not comp}
    return {"n_read": n_read, "n_used": len(rows), "length_mismatch": mismatch,
            "missing_bigrams": missing, "missing_components": comp, "checks": checks,
            "pass": all(checks.values())}


# ---------------------------------------------------------------- statistics
def pvec(tokens, events):
    """Share of tokens whose event set contains each event (an expected-count vector)."""
    n = len(tokens)
    c = Counter(e for t in tokens for e in t)
    return [c[e] / n for e in events]


def markers(pA, pB, events, min_p=0.10, ratio=2.0):
    MA = frozenset(e for e, a, b in zip(events, pA, pB) if a >= min_p and a >= ratio * b)
    MB = frozenset(e for e, a, b in zip(events, pA, pB) if b >= min_p and b >= ratio * a)
    return MA, MB


def both(t, MA, MB):
    return bool(t & MA) and bool(t & MB)


def qsmooth(k, n):
    return (k + 0.5) / (n + 1)


def fit_w(m, pA, pB):
    """Least-squares mixture weight of A: argmin_w |m - (w pA + (1-w) pB)|^2, clipped to [0,1]."""
    d = [a - b for a, b in zip(pA, pB)]
    dd = sum(x * x for x in d)
    if dd == 0:
        return 0.5
    w = sum((mi - b) * di for mi, b, di in zip(m, pB, d)) / dd
    return min(1.0, max(0.0, w))


def binom_sf(k, n, p):
    """P(X >= k), X ~ Binomial(n, p)."""
    if k <= 0:
        return 1.0
    if k > n or p <= 0:
        return 0.0
    if p >= 1:
        return 1.0
    lp = (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
          + k * math.log(p) + (n - k) * math.log1p(-p))
    term = math.exp(lp)
    s = term
    r = p / (1 - p)
    for j in range(k, n):
        term *= (n - j) / (j + 1) * r
        s += term
    return min(1.0, s)


def holm(pdict, alpha=ALPHA):
    items = sorted(pdict.items(), key=lambda kv: (kv[1], kv[0]))
    m = len(items)
    rej = {}
    going = True
    for i, (name, p) in enumerate(items):
        if going and p <= alpha / (m - i):
            rej[name] = True
        else:
            going = False
            rej[name] = False
    return rej


def pool_patterns(tokens):
    c = Counter(tokens)
    pats = sorted(c, key=lambda s: sorted(s))
    return pats, [c[p] for p in pats]


def resample_counts(counts, rng):
    n = sum(counts)
    draws = rng.choices(range(len(counts)), weights=counts, k=n)
    c = Counter(draws)
    return [c[i] for i in range(len(counts))]


def summarize(pats, counts, events, MA, MB):
    idx = {e: i for i, e in enumerate(events)}
    v = [0.0] * len(events)
    n = kb = n2 = 0
    for p, c in zip(pats, counts):
        if not c:
            continue
        n += c
        for e in p:
            v[idx[e]] += c
        if both(p, MA, MB):
            kb += c
        if len(p) >= 2:
            n2 += c
    return {"n": n, "p": [x / n for x in v], "k": kb, "n2": n2}


def boots_pair(tA, tB, events, MA, MB, R, rng):
    """R bootstrap resamples of the single-call pools A and B (the mixture null's ingredients)."""
    pa, ca = pool_patterns(tA)
    pb, cb = pool_patterns(tB)
    out = []
    for _ in range(R):
        a = summarize(pa, resample_counts(ca, rng), events, MA, MB)
        b = summarize(pb, resample_counts(cb, rng), events, MA, MB)
        out.append((a["p"], b["p"], qsmooth(a["k"], a["n"]), qsmooth(b["k"], b["n"]),
                    qsmooth(a["k"], a["n2"]), qsmooth(b["k"], b["n2"])))
    return out


def pvals(m, k, n, n2, boots):
    """T1: P(K >= k) under the best-fitting mixture, averaged over the pool bootstrap.
    T2: the same among multi-event tokens only (k counts only multi-event tokens anyway,
    because the two marker sets are disjoint)."""
    s1 = s2 = 0.0
    for pa, pb, qa, qb, qa2, qb2 in boots:
        w = fit_w(m, pa, pb)
        s1 += binom_sf(k, n, w * qa + (1 - w) * qb)
        if n2:
            s2 += binom_sf(k, n2, w * qa2 + (1 - w) * qb2)
    R = len(boots)
    return s1 / R, (s2 / R if n2 >= MIN_N2 else None)


def ab_test(tokens, events, MA, MB, boots):
    n = len(tokens)
    k = sum(1 for t in tokens if both(t, MA, MB))
    n2 = sum(1 for t in tokens if len(t) >= 2)
    m = pvec(tokens, events)
    p1, p2 = pvals(m, k, n, n2, boots)
    return {"n": n, "k": k, "n2": n2, "J": k / n, "p_t1": p1, "p_t2": p2, "m": m}


def point(m, pA, pB, sA, sB):
    w = fit_w(m, pA, pB)
    return w, w * qsmooth(sA["k"], sA["n"]) + (1 - w) * qsmooth(sB["k"], sB["n"])


def closest(m, sp):
    return sorted((math.dist(m, v), s) for s, v in sp.items())


def pair_closest(m, sp, A, B):
    d = closest(m, sp)
    return len(d) >= 2 and {d[0][1], d[1][1]} == {A, B}


def seed_for(*parts):
    return int(hashlib.sha256(repr(parts).encode()).hexdigest()[:12], 16)


def testability(name, tok, events, th):
    A, B = name.split("_")
    for c in (name, A, B):
        if len(tok.get(c, [])) < MIN_N:
            return False, "%s has fewer than %d tokens" % (c, MIN_N), None, None
    pA, pB = pvec(tok[A], events), pvec(tok[B], events)
    MA, MB = markers(pA, pB, events, th[0], th[1])
    if not MA or not MB:
        who = [c for c, M in ((A, MA), (B, MB)) if not M]
        return False, "no distinctive event for " + " and ".join(who), MA, MB
    return True, "", MA, MB


def calibrate(n, tA, tB, pA, pB, sA, sB, events, MA, MB, w_hat, boots, sp, A, B, S, level,
              rng):
    """C1: S synthetic bigram sets drawn as a pure mixture (each token is a real A token with
    probability w_hat, else a real B token). C2: S sets drawn as a pure conjunction (each token
    is the union of a real A token and a real B token). Same n, same test, same label rule."""
    res = {"S": S, "R_cal": len(boots), "level": level}
    for kind in ("mixture", "conjunction"):
        nom = lab = pair = 0
        for _ in range(S):
            if kind == "mixture":
                syn = [rng.choice(tA) if rng.random() < w_hat else rng.choice(tB)
                       for _ in range(n)]
            else:
                syn = [rng.choice(tA) | rng.choice(tB) for _ in range(n)]
            t = ab_test(syn, events, MA, MB, boots)
            w, pi = point(t["m"], pA, pB, sA, sB)
            if t["p_t1"] <= ALPHA:
                nom += 1
            if t["p_t1"] <= level and t["J"] - pi >= EXCESS and t["k"] >= MIN_K:
                lab += 1
            if pair_closest(t["m"], sp, A, B):
                pair += 1
        res[kind] = {"nominal_rate": nom / S, "label_rate": lab / S,
                     "pair_closest_rate": pair / S, "pair_closest_count": pair}
    return res


def analyze_bigram(name, tok, events, cfg, th, R, level=None, sp=None, tag="main"):
    A, B = name.split("_")
    ok, why, MA, MB = testability(name, tok, events, th)
    rec = {"bigram": name, "A": A, "B": B, "n": len(tok.get(name, [])), "testable": ok}
    if MA is not None:
        rec["markers_A"], rec["markers_B"] = sorted(MA), sorted(MB)
    if not ok:
        rec["why"] = why
        return rec
    tA, tB, tab = tok[A], tok[B], tok[name]
    pA, pB = pvec(tA, events), pvec(tB, events)
    rng = random.Random(seed_for(cfg["seed"], tag, name, tuple(th)))
    boots = boots_pair(tA, tB, events, MA, MB, R, rng)
    sA = summarize(*pool_patterns(tA), events, MA, MB)
    sB = summarize(*pool_patterns(tB), events, MA, MB)
    t = ab_test(tab, events, MA, MB, boots)
    w, pi = point(t["m"], pA, pB, sA, sB)
    rec.update({"k": t["k"], "n2": t["n2"], "J": t["J"], "w": w, "pi_mix": pi,
                "excess": t["J"] - pi, "p_t1": t["p_t1"], "p_t2": t["p_t2"], "R": R,
                "n_A": sA["n"], "n_B": sB["n"],
                "q_A": qsmooth(sA["k"], sA["n"]), "q_B": qsmooth(sB["k"], sB["n"]),
                "multi_rate": {"A": sA["n2"] / sA["n"], "B": sB["n2"] / sB["n"],
                               "AB": t["n2"] / t["n"]}})
    if sp:
        rec["pair_closest_real"] = pair_closest(t["m"], sp, A, B)
    if level is not None:
        crng = random.Random(seed_for(cfg["seed"], "cal", tag, name, tuple(th)))
        rec["calibration"] = calibrate(len(tab), tA, tB, pA, pB, sA, sB, events, MA, MB, w,
                                       boots[:cfg["r_cal"]], sp or {}, A, B, cfg["s_cal"],
                                       level, crng)
    return rec


def label_bigram(rec, rejected):
    if not rec.get("testable"):
        return "NOT_TESTABLE"
    cal = rec.get("calibration")
    if cal and cal["mixture"]["nominal_rate"] > FPR_MAX:
        return "NULL_MISCALIBRATED"
    if rejected:
        if rec["excess"] >= EXCESS and rec["k"] >= MIN_K:
            return "CONJUNCTION_SUPPORTED"
        return "SMALL_EXCESS"
    if cal and cal["conjunction"]["label_rate"] >= POWER_OK:
        return "MIXTURE_NOT_REJECTED"
    return "UNDERPOWERED"


def overall_label(labels):
    v = list(labels.values())
    if "CONJUNCTION_SUPPORTED" in v:
        return "CONJUNCTION_FOUND"
    if v.count("MIXTURE_NOT_REJECTED") >= 2:
        return "NO_CONJUNCTION_POWERED"
    if all(x == "NOT_TESTABLE" for x in v):
        return "NO_TESTABLE_BIGRAMS"
    return "UNDECIDED"


def sens_run(tok, events, cfg, th, R, tag):
    recs = {b: analyze_bigram(b, tok, events, cfg, th, R, tag=tag) for b in cfg["claimed"]}
    rej = holm({b: r["p_t1"] for b, r in recs.items() if r.get("testable")})
    lab = {}
    for b, r in recs.items():
        if not r.get("testable"):
            lab[b] = "NOT_TESTABLE"
        elif rej[b] and r["excess"] >= EXCESS and r["k"] >= MIN_K:
            lab[b] = "CONJUNCTION_SUPPORTED"
        else:
            lab[b] = "NOT_SUPPORTED"
    keep = ("n", "k", "J", "pi_mix", "excess", "p_t1", "testable", "why", "markers_A",
            "markers_B")
    return {"thresholds": list(th), "R": R, "labels": lab,
            "per_bigram": {b: {x: r.get(x) for x in keep} for b, r in recs.items()}}


def run_analysis(rows, cfg):
    tok = call_tokens(rows)
    events = all_events(tok)
    singles = sorted(c for c in tok if "_" not in c and len(tok[c]) >= MIN_N)
    sp = {s: pvec(tok[s], events) for s in singles}
    out = {"events": events, "singles": {s: len(tok[s]) for s in singles},
           "family": list(cfg["claimed"]), "strict": list(cfg["strict"])}
    gr = {}
    for b in cfg["bigrams"]:
        if tok.get(b) and len(sp) >= 2:
            d = closest(pvec(tok[b], events), sp)
            gr[b] = {"ours": d[0][1], "second": d[1][1], "dist": d[0][0],
                     "paper": cfg["paper_closest"].get(b)}
    agree = sum(1 for v in gr.values() if v["ours"] == v["paper"])
    out["gate_r"] = {"per_bigram": gr, "agree": agree, "of": len(cfg["bigrams"]),
                     "need": cfg["gate_r_need"], "pass": agree >= cfg["gate_r_need"]}
    th = (cfg["min_p"], cfg["ratio"])
    fam_t = [b for b in cfg["claimed"] if testability(b, tok, events, th)[0]]
    level = ALPHA / max(1, len(fam_t))
    recs = {}
    for b in cfg["bigrams"]:
        t1 = time.time()
        recs[b] = analyze_bigram(b, tok, events, cfg, th, cfg["r_main"],
                                 level if b in fam_t else None, sp)
        print("bigram", b, "testable" if recs[b]["testable"] else "not testable",
              "%.0f s" % (time.time() - t1), flush=True)
    rej = holm({b: recs[b]["p_t1"] for b in fam_t})
    for b in cfg["claimed"]:
        recs[b]["label"] = label_bigram(recs[b], rej.get(b, False))
    rej2 = holm({b: recs[b]["p_t2"] for b in fam_t if recs[b].get("p_t2") is not None})
    for b in fam_t:
        recs[b]["t2_holm_reject"] = rej2.get(b)
    rej_all = holm({b: r["p_t1"] for b, r in recs.items() if r.get("testable")})
    for b, r in recs.items():
        if not r.get("testable"):
            r["exploratory_all"] = "NOT_TESTABLE"
        elif rej_all[b] and r["excess"] >= EXCESS and r["k"] >= MIN_K:
            r["exploratory_all"] = "REJECT"
        else:
            r["exploratory_all"] = "NOT_REJECT"
    out["bigrams"] = recs
    out["labels"] = {b: recs[b]["label"] for b in cfg["claimed"]}
    out["strict_labels"] = {b: recs[b]["label"] for b in cfg["strict"]}
    out["overall"] = overall_label(out["labels"])
    out["c4"] = {
        "real_pair_closest": {b: recs[b].get("pair_closest_real") for b in fam_t},
        "mixture_pair_closest": sum(recs[b]["calibration"]["mixture"]["pair_closest_count"]
                                    for b in fam_t),
        "conjunction_pair_closest": sum(
            recs[b]["calibration"]["conjunction"]["pair_closest_count"] for b in fam_t),
        "sets_per_kind": len(fam_t) * cfg["s_cal"]}
    s_tok = [t for c, ts in tok.items() if "_" not in c for t in ts]
    b_tok = [t for c, ts in tok.items() if "_" in c for t in ts]
    out["d1"] = {"single_tokens": len(s_tok), "bigram_tokens": len(b_tok),
                 "single_multi_rate": (sum(len(t) >= 2 for t in s_tok) / len(s_tok)
                                       if s_tok else None),
                 "bigram_multi_rate": (sum(len(t) >= 2 for t in b_tok) / len(b_tok)
                                       if b_tok else None)}
    sens = {"S1_dedup": sens_run(call_tokens(dedup(rows)), events, cfg, th, cfg["r_sens"],
                                 "s1")}
    for i, t2 in enumerate(cfg["sens_thresholds"]):
        sens["S2_thresholds_%d" % (i + 1)] = sens_run(tok, events, cfg, tuple(t2),
                                                      cfg["r_sens"], "s2_%d" % i)
    out["sensitivity"] = sens
    return out


def grade(out):
    t = {k: None for k in PRIORS}
    if out.get("gate_d") is not None:
        t["P1"] = bool(out["gate_d"]["pass"])
    if "bigrams" in out:
        B = out["bigrams"]
        tf = [b for b in out["family"] if B[b].get("testable")]
        t["P2"] = bool(out["gate_r"]["pass"])
        t["P3"] = len(tf) >= 4
        t["P4"] = bool(B["HO_GR"].get("testable")) if "HO_GR" in B else None
        t["P5"] = out["overall"] == "CONJUNCTION_FOUND"
        if tf:
            cal = [B[b]["calibration"] for b in tf]
            t["P6"] = all(c["mixture"]["nominal_rate"] <= FPR_MAX for c in cal)
            t["P7"] = statistics.median(c["conjunction"]["label_rate"] for c in cal) >= POWER_OK
            mix = out["c4"]["mixture_pair_closest"]
            conj = out["c4"]["conjunction_pair_closest"]
            t["P8"] = (mix >= 0.5 * conj) if conj > 0 else None
        d1 = out["d1"]
        if d1["single_multi_rate"] is not None and d1["bigram_multi_rate"] is not None:
            t["P9"] = d1["bigram_multi_rate"] > d1["single_multi_rate"]
        sup = [b for b in tf if B[b].get("label") == "CONJUNCTION_SUPPORTED"]
        if sup:
            t["P10"] = all(B[b]["p_t2"] is not None and B[b]["p_t2"] < ALPHA for b in sup)
    graded = [k for k in PRIORS if t[k] is not None]
    brier = (sum((PRIORS[k] - (1.0 if t[k] else 0.0)) ** 2 for k in graded) / len(graded)
             if graded else None)
    return {"truth": t, "priors": PRIORS, "graded": graded, "brier": brier}


def dump(obj, path):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=1, sort_keys=True)


def main(argv=None, fetcher=None, cfg=None):
    cfg = cfg or CFG
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--wait-minutes", type=float, default=90)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    out_path = os.path.join(a.out_dir, "chimp_results.json")
    t0 = time.time()
    with open(os.path.abspath(__file__), "rb") as f:
        code_sha = hashlib.sha256(f.read()).hexdigest()
    out = {"version": VERSION, "python": sys.version.split()[0], "code_sha256": code_sha,
           "config": cfg, "status": "started"}
    reg = wait_for_registration(a.wait_minutes, fetcher)
    out["registration"] = reg
    if not reg["found"]:
        out["status"] = "not_run_unregistered"
        dump(out, out_path)
        print("CHIMP-PREREG.md not public; nothing analysed")
        return 4
    text, info = get_data(a.data_dir, fetcher)
    out["data"] = info
    rows, mismatch, n_read = parse_rows(text)
    out["gate_d"] = gate_d(rows, mismatch, n_read, cfg)
    print("gate D", json.dumps(out["gate_d"]["checks"]), flush=True)
    if not out["gate_d"]["pass"]:
        out["status"] = "gate_failed"
        out["grading"] = grade(out)
        dump(out, out_path)
        return 3
    out.update(run_analysis(rows, cfg))
    out["grading"] = grade(out)
    out["status"] = "ok"
    out["seconds"] = time.time() - t0
    dump(out, out_path)
    print("overall", out["overall"], json.dumps(out["labels"]), flush=True)
    print("grading", json.dumps(out["grading"]["truth"]), "brier", out["grading"]["brier"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
