# SPDX-License-Identifier: Apache-2.0
# test_chimp_results.py: the grader for chimp-bigram-test (CHIMP-PREREG.md, job chimp-bigram-v1),
# written and published BEFORE the job has produced any result. It does not import chimp.py.
# With its own code it re-derives Gate D, Gate R, every Holm decision (T1 over F1, T2, the
# exploratory all-bigram family and the three sensitivity runs), every per-bigram label, the
# overall label, C4, the ten graded predictions and the Brier score from the raw output; checks
# every count, share, smoothed rate, mixture weight and calibration rate for internal consistency;
# checks them against the job's own labels and grading; and fixes now the exact lines
# RESULTS-CHIMP.md must print (render_md_lines).
# Until results/chimp/chimp_results.json exists the result tests are skipped; the toy-data tests
# always run. Stdlib only. MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents.
import copy
import hashlib
import json
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results", "chimp", "chimp_results.json")
MD = os.path.join(HERE, "RESULTS-CHIMP.md")
PREREG = os.path.join(HERE, "CHIMP-PREREG.md")
HAVE = os.path.exists(RES)

VERSION = "chimp-1"
ZIP_MD5 = "96e88c4dff3ceb7f647ad5bd4ca213c0"
ALPHA = 0.05
EXCESS = 0.10
MIN_K = 3
MIN_N = 10
MIN_N2 = 5
FPR_MAX = 0.10
POWER_OK = 0.80

# CHIMP-PREREG.md sections 3-5, typed in again from the registration (not imported).
TABLE1 = {"HO_PN": "PG", "GR_BK": "HO", "GR_PN": "PG", "HO_PG": "PG", "PH_SC": "PH",
          "PH_PB": "PH", "PH_HO": "HO", "HO_PH": "PH", "PN_PG": "PG", "HO_GR": "HO",
          "PG_PN": "PN", "GR_PG": "PG", "PH_PS": "PH", "GR_HO": "HO", "PG_GR": "PG",
          "PH_PG": "PG"}
F1 = ["HO_PH", "PN_PG", "HO_GR", "PG_PN", "PH_HO", "GR_HO", "PG_GR", "PH_PG"]
STRICT = ["HO_PH", "PN_PG", "HO_GR", "PG_PN"]
REGISTERED = {"seed": 20261119, "bigrams": sorted(TABLE1), "claimed": F1, "strict": STRICT,
              "paper_closest": TABLE1, "gate_r_need": 12, "rows_lo": 4300, "rows_hi": 4350,
              "min_p": 0.10, "ratio": 2.0, "sens_thresholds": [[0.05, 1.5], [0.20, 3.0]],
              "r_main": 10000, "r_cal": 500, "s_cal": 200, "r_sens": 2000}

ORDER = ("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10")
PRIORS = {"P1": 0.90, "P2": 0.70, "P3": 0.50, "P4": 0.30, "P5": 0.35,
          "P6": 0.80, "P7": 0.50, "P8": 0.60, "P9": 0.65, "P10": 0.60}
PRED_TEXT = {
    "P1": "Gate D passes",
    "P2": "Gate R passes (closest single matches Table 1 for >= 12 of 16)",
    "P3": "At least 4 of the 8 F1 bigrams are testable",
    "P4": "HO_GR (the paper's lead additive example) is testable",
    "P5": "Overall label is CONJUNCTION_FOUND",
    "P6": "C1 nominal false-positive rate <= 0.10 for every testable F1 bigram",
    "P7": "Median C2 label power over testable F1 bigrams >= 0.80",
    "P8": "C4: mixture sets put A and B as the two closest singles at least half as often as "
          "conjunction sets",
    "P9": "D1: bigrams carry 2+ events more often than single calls",
    "P10": "Every CONJUNCTION_SUPPORTED bigram also has T2 p < 0.05 (graded only if one exists)"}
CONJ, SMALL = "CONJUNCTION_SUPPORTED", "SMALL_EXCESS"
MNR, UNDER = "MIXTURE_NOT_REJECTED", "UNDERPOWERED"
NT, MISCAL = "NOT_TESTABLE", "NULL_MISCALIBRATED"


# --------------------------------------------------------------------------
# independent re-implementation of CHIMP-PREREG.md sections 3 to 7
# --------------------------------------------------------------------------
def holm_set(pdict, alpha=ALPHA):
    """Names rejected by Holm's step-down procedure (ties cannot change the outcome)."""
    items = sorted(pdict.items(), key=lambda kv: (kv[1], kv[0]))
    m = len(items)
    out = set()
    for j, (name, p) in enumerate(items):
        if p > alpha / (m - j):
            break
        out.add(name)
    return out


def median(xs):
    s = sorted(xs)
    n = len(s)
    if n == 0:
        return None
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


def near_int(x, tol=1e-6):
    return abs(x - round(x)) < tol


def label_of(rec, rejected):
    """Section 4 labels, in the registered order."""
    if not rec.get("testable"):
        return NT
    cal = rec.get("calibration")
    if cal and cal["mixture"]["nominal_rate"] > FPR_MAX:
        return MISCAL
    if rejected:
        return CONJ if rec["excess"] >= EXCESS and rec["k"] >= MIN_K else SMALL
    if cal and cal["conjunction"]["label_rate"] >= POWER_OK:
        return MNR
    return UNDER


def overall_of(labels):
    v = list(labels.values())
    if CONJ in v:
        return "CONJUNCTION_FOUND"
    if v.count(MNR) >= 2:
        return "NO_CONJUNCTION_POWERED"
    if v and all(x == NT for x in v):
        return "NO_TESTABLE_BIGRAMS"
    return "UNDECIDED"


def gate_d_checks(gd, cfg):
    n = gd["n_read"]
    return {"rows_in_range": cfg["rows_lo"] <= n <= cfg["rows_hi"],
            "mismatch_at_most_1pct": gd["length_mismatch"] <= 0.01 * n,
            "all_bigrams_present": not gd["missing_bigrams"],
            "all_components_present": not gd["missing_components"]}


def sens_names(cfg):
    return ["S1_dedup"] + ["S2_thresholds_%d" % (i + 1)
                           for i in range(len(cfg["sens_thresholds"]))]


def sens_labels(run):
    pb = run["per_bigram"]
    t = [b for b in pb if pb[b].get("testable")]
    rej = holm_set({b: pb[b]["p_t1"] for b in t})
    out = {}
    for b in pb:
        if not pb[b].get("testable"):
            out[b] = NT
        elif b in rej and pb[b]["excess"] >= EXCESS and pb[b]["k"] >= MIN_K:
            out[b] = CONJ
        else:
            out[b] = "NOT_SUPPORTED"
    return out


def regrade(r):
    cfg = r.get("config") or {}
    g = {k: None for k in ORDER}
    rg = {"g": g}
    gd = r.get("gate_d")
    if gd is not None:
        rg["gate_d_pass"] = all(gate_d_checks(gd, cfg).values())
        g["P1"] = rg["gate_d_pass"]
    B = r.get("bigrams")
    if B is not None:
        fam = list(cfg["claimed"])
        tf = [b for b in fam if B[b].get("testable")]
        rej = holm_set({b: B[b]["p_t1"] for b in tf})
        labels = {b: label_of(B[b], b in rej) for b in fam}
        overall = overall_of(labels)
        t2 = holm_set({b: B[b]["p_t2"] for b in tf if B[b].get("p_t2") is not None})
        t2_rej = {b: ((b in t2) if B[b].get("p_t2") is not None else None) for b in tf}
        t_all = [b for b in cfg["bigrams"] if B[b].get("testable")]
        ra = holm_set({b: B[b]["p_t1"] for b in t_all})
        explo = {}
        for b in cfg["bigrams"]:
            if not B[b].get("testable"):
                explo[b] = NT
            elif b in ra and B[b]["excess"] >= EXCESS and B[b]["k"] >= MIN_K:
                explo[b] = "REJECT"
            else:
                explo[b] = "NOT_REJECT"
        agree = sum(1 for b, v in r["gate_r"]["per_bigram"].items()
                    if v["ours"] == cfg["paper_closest"].get(b))
        mix = sum(B[b]["calibration"]["mixture"]["pair_closest_count"] for b in tf)
        conj = sum(B[b]["calibration"]["conjunction"]["pair_closest_count"] for b in tf)
        rg.update({"fam_t": tf, "level": ALPHA / max(1, len(tf)), "rejected": rej,
                   "labels": labels, "overall": overall, "t2_rej": t2_rej, "explo": explo,
                   "agree": agree, "gate_r_pass": agree >= cfg["gate_r_need"],
                   "c4": {"mix": mix, "conj": conj, "sets": len(tf) * cfg["s_cal"],
                          "real": sum(1 for b in tf if B[b].get("pair_closest_real"))},
                   "sens": {name: sens_labels(run)
                            for name, run in (r.get("sensitivity") or {}).items()}})
        g["P2"] = rg["gate_r_pass"]
        g["P3"] = len(tf) >= 4
        g["P4"] = bool(B["HO_GR"].get("testable")) if "HO_GR" in B else None
        g["P5"] = overall == "CONJUNCTION_FOUND"
        if tf:
            cal = [B[b]["calibration"] for b in tf]
            g["P6"] = all(c["mixture"]["nominal_rate"] <= FPR_MAX for c in cal)
            g["P7"] = median([c["conjunction"]["label_rate"] for c in cal]) >= POWER_OK
            g["P8"] = (mix >= 0.5 * conj) if conj > 0 else None
        d1 = r["d1"]
        if d1["single_multi_rate"] is not None and d1["bigram_multi_rate"] is not None:
            g["P9"] = d1["bigram_multi_rate"] > d1["single_multi_rate"]
        sup = [b for b in tf if labels[b] == CONJ]
        if sup:
            g["P10"] = all(B[b]["p_t2"] is not None and B[b]["p_t2"] < ALPHA for b in sup)
    graded = [k for k in ORDER if g[k] is not None]
    rg["graded"] = graded
    rg["brier"] = (sum((PRIORS[k] - (1.0 if g[k] else 0.0)) ** 2 for k in graded) / len(graded)
                   if graded else None)
    return rg


# --------------------------------------------------------------------------
# configuration and internal consistency of the whole output
# --------------------------------------------------------------------------
def config_problems(r):
    bad = []
    if r.get("version") != VERSION:
        bad.append("version")
    c = r.get("config") or {}
    for k, v in REGISTERED.items():
        if c.get(k) != v:
            bad.append("config " + k)
    return bad


def grading_problems(r, rg):
    bad = []
    gr = r.get("grading") or {}
    if gr.get("truth") != rg["g"]:
        bad.append("grading truth")
    if list(gr.get("graded") or []) != rg["graded"]:
        bad.append("grading graded list")
    if (gr.get("brier") is None) != (rg["brier"] is None) or \
            (rg["brier"] is not None and abs(gr["brier"] - rg["brier"]) > 1e-12):
        bad.append("grading brier")
    return bad


def problems(r):
    bad = []
    cfg = r.get("config") or {}
    rg = regrade(r)
    if not (r.get("registration") or {}).get("found"):
        bad.append("registration")
    gd = r["gate_d"]
    chk = gate_d_checks(gd, cfg)
    if gd["checks"] != chk or gd["pass"] != all(chk.values()):
        bad.append("gate D flags")
    if gd["n_used"] + gd["length_mismatch"] != gd["n_read"]:
        bad.append("gate D counts")
    if r.get("status") == "gate_failed":
        if gd["pass"] or "bigrams" in r:
            bad.append("gate_failed status with a passing gate or a test")
        return bad + grading_problems(r, rg)
    if r.get("status") != "ok" or r.get("version") != VERSION:
        bad.append("status or version")
    if not gd["pass"]:
        bad.append("status ok with a failed Gate D")
    d = r.get("data") or {}
    if d.get("source") == "zenodo" and d.get("zip_md5") != ZIP_MD5:
        bad.append("data md5")
    if len(d.get("csv_sha256", "")) != 64:
        bad.append("data sha256")

    B = r["bigrams"]
    events = set(r["events"])
    singles = r["singles"]
    tf = rg["fam_t"]
    if sorted(B) != sorted(cfg["bigrams"]):
        bad.append("bigram set")
    if any(n < MIN_N for n in singles.values()):
        bad.append("singles below 10")
    if r["family"] != list(cfg["claimed"]) or r["strict"] != list(cfg["strict"]):
        bad.append("family")
    for b in cfg["bigrams"]:
        x = B[b]
        A, Bc = b.split("_")
        tag = " of " + b
        if (x["bigram"], x["A"], x["B"]) != (b, A, Bc):
            bad.append("names" + tag)
        MA, MB = x.get("markers_A"), x.get("markers_B")
        if MA is not None:
            if (set(MA) & set(MB)) or not (set(MA) | set(MB)) <= events:
                bad.append("markers" + tag)
        if x.get("exploratory_all") != rg["explo"][b]:
            bad.append("exploratory label" + tag)
        if not x["testable"]:
            if not x.get("why"):
                bad.append("why" + tag)
            if MA and MB and x["n"] >= MIN_N and singles.get(A, 0) >= MIN_N and \
                    singles.get(Bc, 0) >= MIN_N:
                bad.append("testable flag" + tag)
            if "calibration" in x:
                bad.append("calibration presence" + tag)
            continue
        if not MA or not MB or x["n"] < MIN_N:
            bad.append("testable flag" + tag)
        if x["n_A"] != singles.get(A) or x["n_B"] != singles.get(Bc):
            bad.append("pool sizes" + tag)
        if not 0 <= x["k"] <= x["n2"] <= x["n"]:
            bad.append("k, n2, n" + tag)
        if abs(x["J"] - x["k"] / x["n"]) > 1e-12:
            bad.append("J" + tag)
        if not 0 <= x["w"] <= 1:
            bad.append("w" + tag)
        for q, n in ((x["q_A"], x["n_A"]), (x["q_B"], x["n_B"])):
            kk = q * (n + 1) - 0.5
            if not near_int(kk) or not 0 <= round(kk) <= n:
                bad.append("q" + tag)
        if abs(x["pi_mix"] - (x["w"] * x["q_A"] + (1 - x["w"]) * x["q_B"])) > 1e-9:
            bad.append("pi" + tag)
        if abs(x["excess"] - (x["J"] - x["pi_mix"])) > 1e-12:
            bad.append("excess" + tag)
        if x["R"] != cfg["r_main"] or not 0 <= x["p_t1"] <= 1:
            bad.append("p_t1" + tag)
        if x["k"] == 0 and x["p_t1"] != 1.0:
            bad.append("p_t1 at k = 0" + tag)
        p2 = x.get("p_t2")
        if (p2 is None) != (x["n2"] < MIN_N2) or (p2 is not None and not 0 <= p2 <= 1):
            bad.append("p_t2" + tag)
        mr = x["multi_rate"]
        for key, n in (("A", x["n_A"]), ("B", x["n_B"]), ("AB", x["n"])):
            if not near_int(mr[key] * n) or not 0 <= mr[key] <= 1:
                bad.append("multi rate" + tag)
        if abs(mr["AB"] - x["n2"] / x["n"]) > 1e-12:
            bad.append("multi rate" + tag)
        if not isinstance(x.get("pair_closest_real"), bool):
            bad.append("real pair" + tag)
        if (b in tf) != ("calibration" in x):
            bad.append("calibration presence" + tag)
        if b in tf:
            if x.get("t2_holm_reject") != rg["t2_rej"][b]:
                bad.append("T2 Holm" + tag)
            if x.get("label") != rg["labels"][b]:
                bad.append("job label" + tag)
            c = x["calibration"]
            if c["S"] != cfg["s_cal"] or c["R_cal"] != min(cfg["r_cal"], cfg["r_main"]) or \
                    abs(c["level"] - rg["level"]) > 1e-15:
                bad.append("calibration setup" + tag)
            for kind in ("mixture", "conjunction"):
                ck = c[kind]
                for key in ("nominal_rate", "label_rate", "pair_closest_rate"):
                    if not near_int(ck[key] * c["S"]) or not 0 <= ck[key] <= 1:
                        bad.append("calibration %s %s%s" % (kind, key, tag))
                if abs(ck["pair_closest_rate"] - ck["pair_closest_count"] / c["S"]) > 1e-12:
                    bad.append("calibration %s pair count%s" % (kind, tag))
                if ck["label_rate"] > ck["nominal_rate"] + 1e-12:
                    bad.append("calibration %s label above nominal%s" % (kind, tag))
        elif b in cfg["claimed"] and x.get("label") != rg["labels"][b]:
            bad.append("job label" + tag)
    for b in cfg["claimed"]:
        if not B[b]["testable"] and B[b].get("label") != NT:
            bad.append("job label of " + b)
    if r["labels"] != rg["labels"]:
        bad.append("labels")
    if r["strict_labels"] != {b: rg["labels"][b] for b in cfg["strict"]}:
        bad.append("strict labels")
    if r["overall"] != rg["overall"]:
        bad.append("overall")
    c4 = r["c4"]
    if c4["real_pair_closest"] != {b: B[b]["pair_closest_real"] for b in tf}:
        bad.append("C4 real")
    if (c4["mixture_pair_closest"], c4["conjunction_pair_closest"], c4["sets_per_kind"]) != \
            (rg["c4"]["mix"], rg["c4"]["conj"], rg["c4"]["sets"]):
        bad.append("C4 counts")
    gr = r["gate_r"]
    if (gr["agree"], gr["of"], gr["need"], gr["pass"]) != \
            (rg["agree"], len(cfg["bigrams"]), cfg["gate_r_need"], rg["gate_r_pass"]):
        bad.append("gate R")
    if sorted(gr["per_bigram"]) != sorted(cfg["bigrams"]):
        bad.append("gate R bigram set")
    for b, v in gr["per_bigram"].items():
        if v["paper"] != cfg["paper_closest"].get(b) or v["ours"] == v["second"] or \
                v["ours"] not in singles or v["second"] not in singles or v["dist"] < 0:
            bad.append("gate R row of " + b)
    d1 = r["d1"]
    for rate, n in ((d1["single_multi_rate"], d1["single_tokens"]),
                    (d1["bigram_multi_rate"], d1["bigram_tokens"])):
        if rate is not None and (not near_int(rate * n) or not 0 <= rate <= 1):
            bad.append("D1 rates")
    if d1["single_tokens"] < sum(singles.values()) or \
            d1["bigram_tokens"] < sum(B[b]["n"] for b in cfg["bigrams"]):
        bad.append("D1 token counts")
    sens = r.get("sensitivity") or {}
    if sorted(sens) != sorted(sens_names(cfg)):
        bad.append("sensitivity runs")
    for name in sens_names(cfg):
        run = sens.get(name)
        if run is None:
            continue
        th = [cfg["min_p"], cfg["ratio"]] if name == "S1_dedup" else \
            list(cfg["sens_thresholds"][int(name[-1]) - 1])
        if run["thresholds"] != th or run["R"] != cfg["r_sens"]:
            bad.append("sensitivity setup " + name)
        if sorted(run["per_bigram"]) != sorted(cfg["claimed"]):
            bad.append("sensitivity bigrams " + name)
        if run["labels"] != rg["sens"][name]:
            bad.append("sensitivity labels " + name)
        for b, x in run["per_bigram"].items():
            if not x.get("testable"):
                continue
            if abs(x["J"] - x["k"] / x["n"]) > 1e-12 or \
                    abs(x["excess"] - (x["J"] - x["pi_mix"])) > 1e-12:
                bad.append("sensitivity numbers %s %s" % (name, b))
            if name == "S1_dedup" and x["n"] > B[b]["n"]:
                bad.append("S1 larger than main for " + b)
    return bad + grading_problems(r, rg)


# --------------------------------------------------------------------------
# formats fixed before the result exists (RESULTS-CHIMP.md must contain these lines)
# --------------------------------------------------------------------------
def fp(p):
    if p is None:
        return "n/a"
    return "%.1e" % p if p < 1e-4 else "%.4f" % p


def mk(m):
    return "-" if m is None else (", ".join(m) if m else "none")


def fam_tag(b, cfg):
    return "F1, strict" if b in cfg["strict"] else ("F1" if b in cfg["claimed"] else "-")


def t1_row(b, x, cfg, rg):
    """| bigram | family | n | markers of A | markers of B | k | J | w | pi | J - pi | p T1 |
    p T2 | label |"""
    head = "| {} | {} | {} | {} | {} |".format(b, fam_tag(b, cfg), x["n"],
                                               mk(x.get("markers_A")), mk(x.get("markers_B")))
    if not x["testable"]:
        return head + " - | - | - | - | - | - | - | NOT_TESTABLE ({}) |".format(x.get("why", ""))
    lab = ("`%s`" % rg["labels"][b]) if b in rg["labels"] else \
        "exploratory: " + rg["explo"][b]
    return head + " {} | {:.3f} | {:.3f} | {:.3f} | {:+.3f} | {} | {} | {} |".format(
        x["k"], x["J"], x["w"], x["pi_mix"], x["excess"], fp(x["p_t1"]), fp(x.get("p_t2")), lab)


def cal_row(b, x):
    """| bigram | level | C1 nominal FPR | C1 label rate | C2 label power | C4 mixture |
    C4 conjunction | real |"""
    c = x["calibration"]
    return "| {} | {:.4f} | {:.3f} | {:.3f} | {:.3f} | {} of {} | {} of {} | {} |".format(
        b, c["level"], c["mixture"]["nominal_rate"], c["mixture"]["label_rate"],
        c["conjunction"]["label_rate"], c["mixture"]["pair_closest_count"], c["S"],
        c["conjunction"]["pair_closest_count"], c["S"],
        "yes" if x["pair_closest_real"] else "no")


def label_line(b, x, lab):
    if lab == NT:
        return "{}: **`{}`** ({})".format(b, NT, x.get("why", ""))
    c = x["calibration"]
    return "{}: **`{}`** (C2 label power {:.2f}; C1 false-positive rate {:.2f})".format(
        b, lab, c["conjunction"]["label_rate"], c["mixture"]["nominal_rate"])


def gate_r_row(b, v):
    """| bigram | Table 1 closest | our closest | our second | distance | match |"""
    return "| {} | {} | {} | {} | {:.3f} | {} |".format(
        b, v["paper"], v["ours"], v["second"], v["dist"], "yes" if v["ours"] == v["paper"]
        else "no")


def sens_line(name, run, labs):
    th = run["thresholds"]
    return "{} (markers p >= {:.2f}, {:.1f}x; R {:,}): {}".format(
        name, th[0], th[1], run["R"], ", ".join("%s %s" % (b, labs[b]) for b in sorted(labs)))


def pred_end(k, v):
    return "| {:.2f} | **{}** |".format(
        PRIORS[k], "not graded" if v is None else ("true" if v else "false"))


def _doc(r, rg):
    """(line, required) pairs, in reading order."""
    D = []

    def add(line, req=True):
        D.append((line, req))

    cfg = r.get("config") or {}
    add("# RESULTS-CHIMP: do chimpanzee 'additive' call combinations mean A AND B, or A OR B?",
        False)
    add("", False)
    add("Every line below is printed by render_md_lines in test_chimp_results.py from "
        "results/chimp/chimp_results.json (job chimp-bigram-v1), as registered in "
        "CHIMP-PREREG.md. MV2, Munim, Inc. Made with AI agents; no human expert has checked it "
        "yet.", False)
    add("", False)
    d = r.get("data") or {}
    add("Data: {}, calldata.csv sha256 {}".format(d.get("source", "n/a"),
                                                   d.get("csv_sha256", "n/a")))
    gd = r.get("gate_d")
    if gd is not None:
        add("Gate D: {:,} rows read, {:,} used, {} length mismatches; missing bigrams {}; "
            "missing single calls {}: {}".format(
                gd["n_read"], gd["n_used"], gd["length_mismatch"],
                ", ".join(gd["missing_bigrams"]) or "none",
                ", ".join(gd["missing_components"]) or "none",
                "**pass**" if rg["gate_d_pass"] else "**fail**"))
    B = r.get("bigrams")
    if B is None:
        add("Status: **{}** (no test was run)".format(r.get("status")))
        return D
    gr = r["gate_r"]
    add("Gate R: our closest single call matches Table 1 for {} of {} bigrams (need {}): {}"
        .format(rg["agree"], len(cfg["bigrams"]), cfg["gate_r_need"],
                "**pass**" if rg["gate_r_pass"] else "**fail**"))
    add("", False)
    add("| bigram | Table 1 closest | our closest | our second | distance | match |", False)
    add("|---|---|---|---|---|---|", False)
    for b in cfg["bigrams"]:
        if b in gr["per_bigram"]:
            add(gate_r_row(b, gr["per_bigram"][b]))
    add("", False)
    add("## T1: conjunction versus best-fitting mixture", False)
    add("", False)
    add("Family F1: {} of {} claimed additive bigrams testable; Holm at 0.05 over them; "
        "calibration level {:.4f}".format(len(rg["fam_t"]), len(cfg["claimed"]), rg["level"]))
    add("", False)
    add("| bigram | family | n | markers of A | markers of B | k (both) | J | w | pi (mixture) "
        "| J - pi | p T1 | p T2 | label |", False)
    add("|---|---|---|---|---|---|---|---|---|---|---|---|---|", False)
    for b in cfg["bigrams"]:
        add(t1_row(b, B[b], cfg, rg))
    add("", False)
    for b in cfg["claimed"]:
        add(label_line(b, B[b], rg["labels"][b]))
    add("", False)
    add("Overall label: **`{}`**".format(rg["overall"]))
    add("Strict subset (distance criterion): " +
        ", ".join("%s %s" % (b, rg["labels"][b]) for b in cfg["strict"]))
    add("T2 Holm (length-matched, Holm over F1): " +
        (", ".join("%s %s" % (b, "n/a" if rg["t2_rej"][b] is None else
                              ("yes" if rg["t2_rej"][b] else "no")) for b in rg["fam_t"])
         or "none testable"))
    rej_all = sorted(b for b, v in rg["explo"].items() if v == "REJECT")
    add("Exploratory, all {} bigrams: {} of {} testable rejected ({})".format(
        len(cfg["bigrams"]), len(rej_all), sum(1 for v in rg["explo"].values() if v != NT),
        ", ".join(rej_all) or "none"))
    add("", False)
    add("## Calibration (C1 mixture, C2 conjunction, C4 distance criterion)", False)
    add("", False)
    add("| bigram | level | C1 nominal FPR | C1 label rate | C2 label power | C4 mixture: A, B "
        "closest | C4 conjunction: A, B closest | real: A, B closest |", False)
    add("|---|---|---|---|---|---|---|---|", False)
    for b in rg["fam_t"]:
        add(cal_row(b, B[b]))
    add("", False)
    c4 = rg["c4"]
    add("C4: A and B are the two closest single calls in {} of {} mixture sets and {} of {} "
        "conjunction sets; real bigrams {} of {}".format(
            c4["mix"], c4["sets"], c4["conj"], c4["sets"], c4["real"], len(rg["fam_t"])))
    d1 = r["d1"]
    f3 = (lambda v: "n/a" if v is None else "%.3f" % v)
    add("D1: {} of {:,} bigram utterances carry 2 or more events, versus {} of {:,} single-call "
        "utterances".format(f3(d1["bigram_multi_rate"]), d1["bigram_tokens"],
                            f3(d1["single_multi_rate"]), d1["single_tokens"]))
    add("", False)
    add("## Sensitivity", False)
    add("", False)
    for name in sens_names(cfg):
        if name in (r.get("sensitivity") or {}):
            add(sens_line(name, r["sensitivity"][name], rg["sens"][name]))
    return D


def required_strings(r, rg):
    lines = [ln for ln, req in _doc(r, rg) if req]
    ntrue = sum(1 for k in rg["graded"] if rg["g"][k])
    lines.append("%d of %d graded predictions true" % (ntrue, len(rg["graded"])))
    if rg["graded"]:
        lines.append("Brier score **%.3f**" % rg["brier"])
    return lines


def render_md_lines(r, rg):
    """The whole fixed part of RESULTS-CHIMP.md, prediction table included."""
    out = [ln for ln, _ in _doc(r, rg)]
    out += ["", "## Predictions (CHIMP-PREREG.md section 7)", "",
            "| | Prediction | Prior | Outcome |", "|---|---|---|---|"]
    out += ["| %s | %s %s" % (k, PRED_TEXT[k], pred_end(k, rg["g"][k])) for k in ORDER]
    out += [""] + required_strings(r, rg)[-2:] if rg["graded"] else \
        [""] + required_strings(r, rg)[-1:]
    return out


def check_md(md, r, rg):
    missing = [x for x in required_strings(r, rg) if x not in md]
    lines = md.splitlines()
    for k in ORDER:
        hits = [ln for ln in lines if ln.startswith("| %s |" % k)]
        if len(hits) != 1 or not hits[0].rstrip().endswith(pred_end(k, rg["g"][k])):
            missing.append("prediction row " + k)
    return missing


def sha256_file(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


# --------------------------------------------------------------------------
# toy output (hand-built; the expected values are worked out in the comments)
# --------------------------------------------------------------------------
def _cal(m_nom, m_lab, m_pair, c_nom, c_lab, c_pair, S=200, level=0.025):
    """Calibration block from counts out of S synthetic sets."""
    return {"S": S, "R_cal": 500, "level": level,
            "mixture": {"nominal_rate": m_nom / S, "label_rate": m_lab / S,
                        "pair_closest_rate": m_pair / S, "pair_closest_count": m_pair},
            "conjunction": {"nominal_rate": c_nom / S, "label_rate": c_lab / S,
                            "pair_closest_rate": c_pair / S, "pair_closest_count": c_pair}}


def _tb(name, n, k, n2, nA, nB, kA, kB, w, MA, MB, p1, p2, mrA, mrB, real, cal):
    A, Bc = name.split("_")
    qA = (kA + 0.5) / (nA + 1)
    qB = (kB + 0.5) / (nB + 1)
    J = k / n
    pi = w * qA + (1 - w) * qB
    rec = {"bigram": name, "A": A, "B": Bc, "n": n, "testable": True, "markers_A": MA,
           "markers_B": MB, "k": k, "n2": n2, "J": J, "w": w, "pi_mix": pi, "excess": J - pi,
           "p_t1": p1, "p_t2": p2, "R": 10000, "n_A": nA, "n_B": nB, "q_A": qA, "q_B": qB,
           "multi_rate": {"A": mrA, "B": mrB, "AB": n2 / n}, "pair_closest_real": real}
    if cal:
        rec["calibration"] = cal
    return rec


def _sp(n, k, p, excess, MA, MB):
    J = k / n
    return {"n": n, "k": k, "J": J, "pi_mix": J - excess, "excess": excess, "p_t1": p,
            "testable": True, "why": None, "markers_A": MA, "markers_B": MB}


def _sp_nt(n, why, MA, MB):
    return {"n": n, "k": None, "J": None, "pi_mix": None, "excess": None, "p_t1": None,
            "testable": False, "why": why, "markers_A": MA, "markers_B": MB}


def toy():
    cfg = {"seed": 20261119, "bigrams": ["GR_BK", "HO_GR", "PH_HO"],
           "claimed": ["HO_GR", "PH_HO"], "strict": ["HO_GR"],
           "paper_closest": {"GR_BK": "HO", "HO_GR": "HO", "PH_HO": "HO"}, "gate_r_need": 2,
           "rows_lo": 4300, "rows_hi": 4350, "min_p": 0.1, "ratio": 2.0,
           "sens_thresholds": [[0.05, 1.5], [0.2, 3.0]], "r_main": 10000, "r_cal": 500,
           "s_cal": 200, "r_sens": 2000}
    # HO_GR: n 20, k 8 -> J 0.4; q_A = 2.5/51, q_B = 1.5/31; w 0.6 -> pi 0.04877, excess 0.351
    ho_gr = _tb("HO_GR", 20, 8, 10, 50, 30, 2, 1, 0.6, ["FE"], ["AG", "TR"], 1e-6, 2e-4,
                0.2, 0.2, True, _cal(8, 0, 100, 200, 190, 80))
    # PH_HO: n 15, k 1 -> J 0.0667; pi = 0.5 * 3.5/41 + 0.5 * 2.5/51 = 0.0672; n2 4 < 5 -> no T2
    ph_ho = _tb("PH_HO", 15, 1, 4, 40, 50, 3, 2, 0.5, ["RE"], ["FE"], 0.3, None,
                0.2, 0.2, False, _cal(10, 0, 60, 180, 170, 100))
    gr_bk = {"bigram": "GR_BK", "A": "GR", "B": "BK", "n": 12, "testable": False,
             "markers_A": ["AG"], "markers_B": [], "why": "no distinctive event for BK"}
    # Holm over F1 {HO_GR 1e-6, PH_HO 0.3}: 1e-6 <= 0.025 rejects, 0.3 > 0.05 stops.
    # HO_GR: excess 0.351 >= 0.10 and k 8 >= 3 -> CONJUNCTION_SUPPORTED.
    # PH_HO: C1 0.05 <= 0.10, not rejected, C2 power 0.85 >= 0.80 -> MIXTURE_NOT_REJECTED.
    ho_gr.update({"label": CONJ, "t2_holm_reject": True, "exploratory_all": "REJECT"})
    ph_ho.update({"label": MNR, "t2_holm_reject": None, "exploratory_all": "NOT_REJECT"})
    gr_bk["exploratory_all"] = NT
    sens = {
        # Holm {0.001, 0.4}: 0.001 <= 0.025 rejects; excess 0.3, k 7 -> CONJ
        "S1_dedup": {"thresholds": [0.1, 2.0], "R": 2000,
                     "labels": {"HO_GR": CONJ, "PH_HO": "NOT_SUPPORTED"},
                     "per_bigram": {"HO_GR": _sp(18, 7, 0.001, 0.3, ["FE"], ["AG", "TR"]),
                                    "PH_HO": _sp(14, 1, 0.4, 0.0, ["RE"], ["FE"])}},
        # Holm {0.04, 0.4}: 0.04 > 0.025 -> nothing rejected
        "S2_thresholds_1": {"thresholds": [0.05, 1.5], "R": 2000,
                            "labels": {"HO_GR": "NOT_SUPPORTED", "PH_HO": "NOT_SUPPORTED"},
                            "per_bigram": {"HO_GR": _sp(20, 8, 0.04, 0.3, ["FE"], ["AG", "TR"]),
                                           "PH_HO": _sp(15, 1, 0.4, 0.0, ["RE"], ["FE"])}},
        # Holm {0.01}: rejects, but excess 0.05 < 0.10 -> NOT_SUPPORTED
        "S2_thresholds_2": {"thresholds": [0.2, 3.0], "R": 2000,
                            "labels": {"HO_GR": NT, "PH_HO": "NOT_SUPPORTED"},
                            "per_bigram": {"HO_GR": _sp_nt(20, "no distinctive event for GR",
                                                           ["FE"], []),
                                           "PH_HO": _sp(15, 2, 0.01, 0.05, ["RE"], ["FE"])}}}
    truth = {"P1": True, "P2": True, "P3": False, "P4": True, "P5": True, "P6": True,
             "P7": True, "P8": True, "P9": True, "P10": True}
    # Brier: (.01 + .09 + .25 + .49 + .4225 + .04 + .25 + .16 + .1225 + .16) / 10 = 0.1995
    return {"version": VERSION, "status": "ok", "config": cfg,
            "registration": {"found": True, "sha256": "0" * 64},
            "data": {"source": "zenodo", "zip_md5": ZIP_MD5, "csv_sha256": "a" * 64},
            "gate_d": {"n_read": 4320, "n_used": 4310, "length_mismatch": 10,
                       "missing_bigrams": [], "missing_components": [],
                       "checks": {"rows_in_range": True, "mismatch_at_most_1pct": True,
                                  "all_bigrams_present": True, "all_components_present": True},
                       "pass": True},
            "events": ["AG", "FE", "RE", "TR"],
            "singles": {"BK": 20, "GR": 30, "HO": 50, "PH": 40},
            "family": ["HO_GR", "PH_HO"], "strict": ["HO_GR"],
            "gate_r": {"per_bigram": {
                "GR_BK": {"ours": "GR", "second": "BK", "dist": 0.3, "paper": "HO"},
                "HO_GR": {"ours": "HO", "second": "GR", "dist": 0.1, "paper": "HO"},
                "PH_HO": {"ours": "HO", "second": "PH", "dist": 0.2, "paper": "HO"}},
                "agree": 2, "of": 3, "need": 2, "pass": True},
            "bigrams": {"GR_BK": gr_bk, "HO_GR": ho_gr, "PH_HO": ph_ho},
            "labels": {"HO_GR": CONJ, "PH_HO": MNR}, "strict_labels": {"HO_GR": CONJ},
            "overall": "CONJUNCTION_FOUND",
            "c4": {"real_pair_closest": {"HO_GR": True, "PH_HO": False},
                   "mixture_pair_closest": 160, "conjunction_pair_closest": 180,
                   "sets_per_kind": 400},
            "d1": {"single_tokens": 200, "bigram_tokens": 50, "single_multi_rate": 0.2,
                   "bigram_multi_rate": 0.4},
            "sensitivity": sens,
            "grading": {"truth": truth, "priors": dict(PRIORS), "graded": list(ORDER),
                        "brier": 0.1995}}


def toy_gate_failed():
    r = toy()
    gd = {"n_read": 4000, "n_used": 3990, "length_mismatch": 10, "missing_bigrams": [],
          "missing_components": [],
          "checks": {"rows_in_range": False, "mismatch_at_most_1pct": True,
                     "all_bigrams_present": True, "all_components_present": True},
          "pass": False}
    truth = {k: None for k in ORDER}
    truth["P1"] = False
    return {"version": VERSION, "status": "gate_failed", "config": r["config"],
            "registration": {"found": True}, "data": r["data"], "gate_d": gd,
            "grading": {"truth": truth, "priors": dict(PRIORS), "graded": ["P1"],
                        "brier": 0.81}}


class TestToy(unittest.TestCase):
    def test_helpers(self):
        self.assertEqual(holm_set({"a": 0.01, "b": 0.04, "c": 0.03}), {"a"})
        self.assertEqual(holm_set({"a": 0.01, "b": 0.02, "c": 0.03}), {"a", "b", "c"})
        self.assertEqual(holm_set({}), set())
        self.assertAlmostEqual(median([0.95, 0.85]), 0.9)
        self.assertEqual(median([3, 1, 2]), 2)
        self.assertIsNone(median([]))
        self.assertTrue(near_int(2.0000000001) and not near_int(2.05))
        self.assertEqual(fp(None), "n/a")
        self.assertEqual(fp(1e-6), "1.0e-06")
        self.assertEqual(fp(0.3), "0.3000")
        self.assertEqual(mk(None), "-")
        self.assertEqual(mk([]), "none")
        self.assertEqual(mk(["AG", "TR"]), "AG, TR")

    def test_label_rules(self):
        ok = _cal(10, 0, 0, 200, 170, 0)                 # C1 0.05, C2 power 0.85
        weak = _cal(10, 0, 0, 200, 100, 0)               # C2 power 0.50
        bad = _cal(21, 0, 0, 200, 170, 0)                # C1 0.105 > 0.10
        r = {"testable": True, "excess": 0.10, "k": 3, "calibration": ok}
        self.assertEqual(label_of(r, True), CONJ)       # 0.10 >= 0.10 and 3 >= 3
        self.assertEqual(label_of(dict(r, k=2), True), SMALL)
        self.assertEqual(label_of(dict(r, excess=0.09), True), SMALL)
        self.assertEqual(label_of(r, False), MNR)
        self.assertEqual(label_of(dict(r, calibration=weak), False), UNDER)
        self.assertEqual(label_of(dict(r, calibration=bad), True), MISCAL)
        self.assertEqual(label_of({"testable": False}, True), NT)
        self.assertEqual(overall_of({"a": MNR, "b": MNR}), "NO_CONJUNCTION_POWERED")
        self.assertEqual(overall_of({"a": NT, "b": NT}), "NO_TESTABLE_BIGRAMS")
        self.assertEqual(overall_of({"a": UNDER, "b": NT}), "UNDECIDED")
        self.assertEqual(overall_of({"a": CONJ, "b": MISCAL}), "CONJUNCTION_FOUND")

    def test_toy_regrade(self):
        r = toy()
        rg = regrade(r)
        self.assertEqual(rg["g"], r["grading"]["truth"])
        self.assertAlmostEqual(rg["brier"], 0.1995, places=12)
        self.assertEqual(rg["labels"], r["labels"])
        self.assertEqual(rg["overall"], "CONJUNCTION_FOUND")
        self.assertEqual(rg["t2_rej"], {"HO_GR": True, "PH_HO": None})
        self.assertEqual(rg["explo"], {"GR_BK": NT, "HO_GR": "REJECT", "PH_HO": "NOT_REJECT"})
        self.assertEqual(rg["c4"], {"mix": 160, "conj": 180, "sets": 400, "real": 1})
        self.assertEqual(problems(r), [])

    def test_toy_problems_detected(self):
        r = toy()
        r["bigrams"]["HO_GR"]["label"] = SMALL
        self.assertIn("job label of HO_GR", problems(r))
        r = toy()
        r["bigrams"]["PH_HO"]["q_A"] = 0.1
        self.assertIn("q of PH_HO", problems(r))
        r = toy()
        r["c4"]["mixture_pair_closest"] = 161
        self.assertIn("C4 counts", problems(r))
        r = toy()
        r["sensitivity"]["S2_thresholds_1"]["labels"]["HO_GR"] = CONJ
        self.assertIn("sensitivity labels S2_thresholds_1", problems(r))
        r = toy()
        r["gate_d"]["n_read"] = 4400
        p = problems(r)
        self.assertIn("gate D flags", p)
        self.assertIn("gate D counts", p)
        self.assertFalse(regrade(r)["g"]["P1"])
        r = toy()
        r["bigrams"]["PH_HO"]["calibration"]["mixture"]["nominal_rate"] = 0.2
        p = problems(r)
        self.assertIn("labels", p)                       # PH_HO is now NULL_MISCALIBRATED
        self.assertIn("grading truth", p)                # and P6 is now false
        r = toy()
        r["d1"]["bigram_multi_rate"] = 0.401
        self.assertIn("D1 rates", problems(r))
        r = toy()
        r["bigrams"]["PH_HO"]["p_t2"] = 0.5               # n2 = 4 < 5 must have no T2 p
        self.assertIn("p_t2 of PH_HO", problems(r))

    def test_flip_to_no_conjunction(self):
        r = toy()
        r["bigrams"]["HO_GR"]["p_t1"] = 0.2
        rg = regrade(r)
        # Holm {0.2, 0.3}: 0.2 > 0.025 -> none; both have C2 power >= 0.80 -> MNR, MNR
        self.assertEqual(rg["labels"], {"HO_GR": MNR, "PH_HO": MNR})
        self.assertEqual(rg["overall"], "NO_CONJUNCTION_POWERED")
        self.assertFalse(rg["g"]["P5"])
        self.assertIsNone(rg["g"]["P10"])
        self.assertEqual(len(rg["graded"]), 9)

    def test_config(self):
        self.assertIn("config bigrams", config_problems(toy()))
        full = {"version": VERSION, "config": json.loads(json.dumps(REGISTERED))}
        self.assertEqual(config_problems(full), [])
        self.assertEqual(len(REGISTERED["bigrams"]), 16)
        self.assertTrue(set(F1) <= set(TABLE1) and set(STRICT) <= set(F1))

    def test_gate_failed(self):
        r = toy_gate_failed()
        rg = regrade(r)
        self.assertEqual(rg["graded"], ["P1"])
        self.assertAlmostEqual(rg["brier"], 0.81, places=12)
        self.assertEqual(problems(r), [])
        lines = render_md_lines(r, rg)
        self.assertIn("Status: **gate_failed** (no test was run)", lines)
        self.assertIn("| P1 | Gate D passes | 0.90 | **false** |", lines)
        self.assertIn("| P5 | Overall label is CONJUNCTION_FOUND | 0.35 | **not graded** |", lines)
        self.assertEqual(check_md("\n".join(lines), r, rg), [])

    def test_render_and_check(self):
        r = toy()
        rg = regrade(r)
        lines = render_md_lines(r, rg)
        self.assertIn("Gate D: 4,320 rows read, 4,310 used, 10 length mismatches; missing "
                      "bigrams none; missing single calls none: **pass**", lines)
        self.assertIn("Gate R: our closest single call matches Table 1 for 2 of 3 bigrams "
                      "(need 2): **pass**", lines)
        self.assertIn("| GR_BK | HO | GR | BK | 0.300 | no |", lines)
        self.assertIn("Family F1: 2 of 2 claimed additive bigrams testable; Holm at 0.05 over "
                      "them; calibration level 0.0250", lines)
        self.assertIn("| HO_GR | F1, strict | 20 | FE | AG, TR | 8 | 0.400 | 0.600 | 0.049 | "
                      "+0.351 | 1.0e-06 | 0.0002 | `CONJUNCTION_SUPPORTED` |", lines)
        self.assertIn("| GR_BK | - | 12 | AG | none | - | - | - | - | - | - | - | NOT_TESTABLE "
                      "(no distinctive event for BK) |", lines)
        self.assertIn("HO_GR: **`CONJUNCTION_SUPPORTED`** (C2 label power 0.95; C1 "
                      "false-positive rate 0.04)", lines)
        self.assertIn("Overall label: **`CONJUNCTION_FOUND`**", lines)
        self.assertIn("Strict subset (distance criterion): HO_GR CONJUNCTION_SUPPORTED", lines)
        self.assertIn("T2 Holm (length-matched, Holm over F1): HO_GR yes, PH_HO n/a", lines)
        self.assertIn("Exploratory, all 3 bigrams: 1 of 2 testable rejected (HO_GR)", lines)
        self.assertIn("| HO_GR | 0.0250 | 0.040 | 0.000 | 0.950 | 100 of 200 | 80 of 200 | yes |",
                      lines)
        self.assertIn("C4: A and B are the two closest single calls in 160 of 400 mixture sets "
                      "and 180 of 400 conjunction sets; real bigrams 1 of 2", lines)
        self.assertIn("D1: 0.400 of 50 bigram utterances carry 2 or more events, versus 0.200 "
                      "of 200 single-call utterances", lines)
        self.assertIn("S1_dedup (markers p >= 0.10, 2.0x; R 2,000): HO_GR CONJUNCTION_SUPPORTED, "
                      "PH_HO NOT_SUPPORTED", lines)
        self.assertIn("S2_thresholds_2 (markers p >= 0.20, 3.0x; R 2,000): HO_GR NOT_TESTABLE, "
                      "PH_HO NOT_SUPPORTED", lines)
        self.assertIn("| P3 | At least 4 of the 8 F1 bigrams are testable | 0.50 | **false** |",
                      lines)
        self.assertIn("9 of 10 graded predictions true", lines)
        self.assertTrue(any(ln.startswith("Brier score **") for ln in lines))
        md = "\n".join(lines)
        self.assertEqual(check_md(md, r, rg), [])
        self.assertIn("prediction row P3", check_md(md.replace("| P3 |", "| P 3 |"), r, rg))
        self.assertTrue(check_md(md.replace("Overall label", "Overall"), r, rg))


@unittest.skipUnless(HAVE, "results/chimp/chimp_results.json not published yet")
class TestChimpResults(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(RES, encoding="utf-8") as f:
            cls.r = json.load(f)
        cls.rg = regrade(cls.r)

    def test_registration(self):
        if os.path.exists(PREREG):
            self.assertEqual(self.r["registration"]["sha256"], sha256_file(PREREG))

    def test_config(self):
        self.assertEqual(config_problems(self.r), [])

    def test_consistency(self):
        self.assertEqual(problems(self.r), [])

    def test_results_md(self):
        self.assertTrue(os.path.exists(MD), "RESULTS-CHIMP.md missing")
        with open(MD, encoding="utf-8") as f:
            md = f.read()
        self.assertEqual(check_md(md, self.r, self.rg), [])


if __name__ == "__main__":
    unittest.main()
