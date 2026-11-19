# SPDX-License-Identifier: Apache-2.0
# Tests for chimp.py (stdlib only). Expected values are worked out by hand in the comments.
# MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents. Apache-2.0.
import hashlib
import io
import json
import math
import os
import random
import shutil
import tempfile
import unittest
import zipfile

import chimp as C


def fs(*e):
    return frozenset(e)


def toy_tok(conj=True, n_ab=20):
    # A: 30 x {a} + 10 x {a,c}; B: 40 x {b}.  events a, b, c.
    # pA = (a 1.0, b 0, c 0.25); pB = (0, 1, 0)  ->  markers A = {a, c}, B = {b}.
    tok = {"A": [fs("a")] * 30 + [fs("a", "c")] * 10, "B": [fs("b")] * 40}
    if conj:
        tok["A_B"] = [fs("a", "b")] * n_ab
    else:
        tok["A_B"] = [fs("a")] * (n_ab // 2) + [fs("b")] * (n_ab - n_ab // 2)
    return tok


def rows_from(tok):
    rows = []
    for c in sorted(tok):
        for i, t in enumerate(tok[c]):
            rows.append({"date": "d%d" % (i % 3), "caller": "id%d" % (i % 4), "call": c,
                         "events": t})
    return rows


def toy_cfg():
    return {"seed": 7, "bigrams": ["A_B"], "claimed": ["A_B"], "strict": ["A_B"],
            "paper_closest": {"A_B": "B"}, "gate_r_need": 1, "rows_lo": 90, "rows_hi": 110,
            "min_p": 0.10, "ratio": 2.0, "sens_thresholds": [[0.05, 1.5], [0.20, 3.0]],
            "r_main": 40, "r_cal": 20, "s_cal": 10, "r_sens": 20}


def csv_text(rows):
    lines = ["date,caller,call.name,combi.length,context,context.length,demosubset"]
    for r in rows:
        ev = sorted(r["events"])
        lines.append("%s,%s,%s,%d,%s,%d,0" % (r["date"], r["caller"], r["call"],
                                              r["call"].count("_") + 1, "_".join(ev), len(ev)))
    return "\n".join(lines) + "\n"


class TestBasics(unittest.TestCase):
    def test_binom_sf_small(self):
        self.assertAlmostEqual(C.binom_sf(1, 2, 0.5), 0.75)     # 1 - 0.25
        self.assertAlmostEqual(C.binom_sf(2, 2, 0.5), 0.25)
        self.assertEqual(C.binom_sf(0, 5, 0.3), 1.0)
        self.assertEqual(C.binom_sf(3, 2, 0.5), 0.0)
        self.assertEqual(C.binom_sf(1, 5, 0.0), 0.0)

    def test_binom_sf_brute(self):
        n, p = 7, 0.3
        for k in range(0, 9):
            want = sum(math.comb(n, j) * p ** j * (1 - p) ** (n - j) for j in range(k, n + 1))
            self.assertAlmostEqual(C.binom_sf(k, n, p), want, places=12)

    def test_parse_rows(self):
        text = ("date,caller,call.name,combi.length,context,context.length,demosubset\n"
                "2019-01-03,id10,GR_PG,2,FU_RA_RE,3,0\n"
                "2019-01-03,id10,GR,1,AP,1,1\n"
                "2019-01-03,id10,BK,1,BA_FE,3,0\n")            # 2 codes, length 3: mismatch
        rows, mm, n = C.parse_rows(text)
        self.assertEqual((len(rows), mm, n), (2, 1, 3))
        self.assertEqual(rows[0]["call"], "GR_PG")
        self.assertEqual(rows[0]["events"], fs("FU", "RA", "RE"))
        self.assertEqual(rows[1]["caller"], "id10")
        quoted = "\n".join('"%s"' % l for l in text.splitlines())
        rows2, mm2, n2 = C.parse_rows(quoted)
        self.assertEqual((len(rows2), mm2, n2), (2, 1, 3))
        rows3, mm3, n3 = C.parse_rows(chr(0xFEFF) + text)
        self.assertEqual((len(rows3), mm3, n3), (2, 1, 3))

    def test_markers(self):
        # a: 0.5 >= 0.1 and >= 2*0.1 -> A.  b: 0.3 >= 0.1 and >= 2*0.1 -> B.  c: below 0.1.
        MA, MB = C.markers([0.5, 0.1, 0.0], [0.1, 0.3, 0.05], ["a", "b", "c"])
        self.assertEqual((MA, MB), (fs("a"), fs("b")))

    def test_fit_w(self):
        # d = (1,-1), |d|^2 = 2; (m - pB).d = 0.3 + 0.3 = 0.6 -> w = 0.3
        self.assertAlmostEqual(C.fit_w([0.3, 0.7], [1, 0], [0, 1]), 0.3)
        self.assertEqual(C.fit_w([2, 0], [1, 0], [0, 1]), 1.0)          # 3/2 clipped
        self.assertEqual(C.fit_w([0.2, 0.2], [0.5, 0.5], [0.5, 0.5]), 0.5)

    def test_patterns_and_summarize(self):
        toks = [fs("a"), fs("a", "b"), fs("b"), fs("a", "b")]
        pats, cnt = C.pool_patterns(toks)
        self.assertEqual(pats, [fs("a"), fs("a", "b"), fs("b")])
        self.assertEqual(cnt, [1, 2, 1])
        s = C.summarize(pats, cnt, ["a", "b"], fs("a"), fs("b"))
        self.assertEqual((s["n"], s["k"], s["n2"]), (4, 2, 2))
        self.assertEqual(s["p"], [0.75, 0.75])

    def test_resample_counts(self):
        rng = random.Random(1)
        for _ in range(20):
            c = C.resample_counts([3, 0, 5], rng)
            self.assertEqual(sum(c), 8)
            self.assertEqual(c[1], 0)

    def test_holm(self):
        self.assertEqual(C.holm({"x": 0.01, "y": 0.02, "z": 0.04}),
                         {"x": True, "y": True, "z": True})
        self.assertEqual(C.holm({"x": 0.02, "y": 0.2, "z": 0.001}),
                         {"z": True, "x": True, "y": False})
        self.assertEqual(C.holm({"x": 0.03, "y": 0.03}), {"x": False, "y": False})
        self.assertEqual(C.holm({}), {})

    def test_dedup(self):
        r = {"date": "d", "caller": "c", "call": "A", "events": fs("a")}
        self.assertEqual(len(C.dedup([r, dict(r), dict(r, caller="e")])), 2)


class TestNull(unittest.TestCase):
    def setUp(self):
        self.events = ["a", "b", "c"]
        self.MA, self.MB = fs("a", "c"), fs("b")

    def boots(self, tok, R=30):
        return C.boots_pair(tok["A"], tok["B"], self.events, self.MA, self.MB, R,
                            random.Random(3))

    def test_conjunction_rejected(self):
        tok = toy_tok(True)
        t = C.ab_test(tok["A_B"], self.events, self.MA, self.MB, self.boots(tok))
        # no A or B token has both marker kinds: q = 0.5/41 = 0.0122 for every resample;
        # all 20 bigram tokens have both -> p = 0.0122^20, far below 1e-6
        self.assertEqual((t["n"], t["k"], t["n2"]), (20, 20, 20))
        self.assertLess(t["p_t1"], 1e-6)
        self.assertLess(t["p_t2"], 0.05)

    def test_mixture_not_rejected(self):
        tok = toy_tok(False)
        t = C.ab_test(tok["A_B"], self.events, self.MA, self.MB, self.boots(tok))
        self.assertEqual((t["k"], t["p_t1"]), (0, 1.0))      # P(K >= 0) = 1
        self.assertIsNone(t["p_t2"])                        # no multi-event tokens (< 5)

    def test_point(self):
        tok = toy_tok(True)
        pA, pB = C.pvec(tok["A"], self.events), C.pvec(tok["B"], self.events)
        self.assertEqual(pA, [1.0, 0.0, 0.25])
        sA = C.summarize(*C.pool_patterns(tok["A"]), self.events, self.MA, self.MB)
        sB = C.summarize(*C.pool_patterns(tok["B"]), self.events, self.MA, self.MB)
        # m = (1,1,0); d = (1,-1,0.25), |d|^2 = 2.0625; (m-pB).d = 1 -> w = 1/2.0625
        w, pi = C.point([1.0, 1.0, 0.0], pA, pB, sA, sB)
        self.assertAlmostEqual(w, 1 / 2.0625)
        self.assertAlmostEqual(pi, 0.5 / 41)

    def test_testability(self):
        tok = toy_tok(True)
        ok, why, MA, MB = C.testability("A_B", tok, self.events, (0.10, 2.0))
        self.assertTrue(ok)
        self.assertEqual((MA, MB), (self.MA, self.MB))
        tok2 = dict(tok, B=[fs("a")] * 40)          # B identical in meaning to A's main event
        ok2, why2, _, _ = C.testability("A_B", tok2, self.events, (0.10, 2.0))
        self.assertFalse(ok2)
        self.assertIn("no distinctive event", why2)
        tok3 = dict(tok, A_B=[fs("a", "b")] * 9)
        self.assertFalse(C.testability("A_B", tok3, self.events, (0.10, 2.0))[0])


class TestPipeline(unittest.TestCase):
    def test_conjunction_end_to_end(self):
        out = C.run_analysis(rows_from(toy_tok(True)), toy_cfg())
        b = out["bigrams"]["A_B"]
        # distances from m = (1,1,0): to A sqrt(1.0625) = 1.031, to B 1.0 -> closest B
        self.assertEqual(out["gate_r"]["per_bigram"]["A_B"]["ours"], "B")
        self.assertTrue(out["gate_r"]["pass"])
        self.assertTrue(b["testable"])
        self.assertEqual((b["k"], b["n"]), (20, 20))
        self.assertGreater(b["excess"], 0.9)
        self.assertEqual(b["calibration"]["mixture"]["nominal_rate"], 0.0)
        self.assertEqual(b["calibration"]["conjunction"]["label_rate"], 1.0)
        self.assertEqual(b["label"], "CONJUNCTION_SUPPORTED")
        self.assertEqual(out["overall"], "CONJUNCTION_FOUND")
        self.assertTrue(b["t2_holm_reject"])
        self.assertEqual(b["exploratory_all"], "REJECT")
        # d1: singles 10 of 80 multi-event; bigram 20 of 20
        self.assertAlmostEqual(out["d1"]["single_multi_rate"], 0.125)
        self.assertEqual(out["d1"]["bigram_multi_rate"], 1.0)
        # S1 keeps 12 bigram rows (12 distinct (date, caller) pairs), still testable
        s1 = out["sensitivity"]["S1_dedup"]
        self.assertEqual(s1["per_bigram"]["A_B"]["n"], 12)
        self.assertEqual(s1["labels"]["A_B"], "CONJUNCTION_SUPPORTED")
        self.assertIn("S2_thresholds_2", out["sensitivity"])
        g = C.grade(out)["truth"]
        self.assertEqual((g["P1"], g["P2"], g["P3"], g["P4"], g["P5"]),
                         (None, True, False, None, True))
        self.assertEqual((g["P6"], g["P7"], g["P8"], g["P9"], g["P10"]),
                         (True, True, True, True, True))
        json.dumps(out)                              # serialisable

    def test_mixture_end_to_end(self):
        out = C.run_analysis(rows_from(toy_tok(False)), toy_cfg())
        b = out["bigrams"]["A_B"]
        self.assertEqual((b["k"], b["p_t1"]), (0, 1.0))
        self.assertEqual(b["label"], "MIXTURE_NOT_REJECTED")    # power 1.0 from C2
        self.assertEqual(out["overall"], "UNDECIDED")           # only one powered bigram
        g = C.grade(out)
        self.assertFalse(g["truth"]["P5"])
        self.assertIsNone(g["truth"]["P10"])

    def test_labels(self):
        cal = {"mixture": {"nominal_rate": 0.2}, "conjunction": {"label_rate": 0.9}}
        self.assertEqual(C.label_bigram({"testable": True, "calibration": cal}, True),
                         "NULL_MISCALIBRATED")
        cal["mixture"]["nominal_rate"] = 0.05
        rec = {"testable": True, "calibration": cal, "excess": 0.05, "k": 9}
        self.assertEqual(C.label_bigram(rec, True), "SMALL_EXCESS")
        self.assertEqual(C.label_bigram(rec, False), "MIXTURE_NOT_REJECTED")
        cal["conjunction"]["label_rate"] = 0.5
        self.assertEqual(C.label_bigram(rec, False), "UNDERPOWERED")
        self.assertEqual(C.label_bigram({"testable": False}, True), "NOT_TESTABLE")
        self.assertEqual(C.overall_label({"x": "NOT_TESTABLE"}), "NO_TESTABLE_BIGRAMS")
        self.assertEqual(C.overall_label({"x": "MIXTURE_NOT_REJECTED",
                                          "y": "MIXTURE_NOT_REJECTED"}),
                         "NO_CONJUNCTION_POWERED")

    def test_gate_d(self):
        rows = rows_from(toy_tok(True))                 # 100 rows
        cfg = toy_cfg()
        self.assertTrue(C.gate_d(rows, 0, 100, cfg)["pass"])
        self.assertFalse(C.gate_d(rows, 2, 100, cfg)["pass"])      # 2 > 1% of 100
        self.assertFalse(C.gate_d(rows, 0, 120, cfg)["pass"])
        cfg["bigrams"] = ["A_B", "A_C"]
        g = C.gate_d(rows, 0, 100, cfg)
        self.assertEqual((g["missing_bigrams"], g["missing_components"]), (["A_C"], ["C"]))


class TestIO(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d)

    def zipblob(self, csv_bytes):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as zf:
            zf.writestr("tozbu-Chimpanzee_bigram_meaning-x/._calldata.csv", b"\x00junk")
            zf.writestr("tozbu-Chimpanzee_bigram_meaning-x/calldata.csv", csv_bytes)
        return buf.getvalue()

    def test_get_data_fallback_and_zenodo(self):
        csvb = csv_text(rows_from(toy_tok(True))).encode()
        blob = self.zipblob(csvb)
        calls = []

        def fake(url, timeout=60):
            calls.append(url)
            return blob if url == C.ZENODO_URL else csvb
        text, info = C.get_data(self.d, fake)                   # md5 differs -> GitHub
        self.assertEqual(info["source"], "github-tag-chimpbigrams")
        self.assertIn("md5 mismatch", info["zenodo_error"])
        self.assertEqual(text, csvb.decode())
        old = C.ZENODO_MD5
        try:
            C.ZENODO_MD5 = hashlib.md5(blob).hexdigest()
            text2, info2 = C.get_data(self.d, fake)
            self.assertEqual(info2["source"], "zenodo")
            self.assertEqual(info2["member"], "tozbu-Chimpanzee_bigram_meaning-x/calldata.csv")
            self.assertEqual(text2, csvb.decode())
            self.assertTrue(os.path.exists(os.path.join(self.d, "chimpbigrams.zip")))
        finally:
            C.ZENODO_MD5 = old

    def test_wait_for_registration(self):
        ok = C.wait_for_registration(-1, fetcher=lambda u, timeout=60: b"# CHIMP-PREREG ...")
        self.assertTrue(ok["found"])

        def missing(u, timeout=60):
            raise OSError("HTTP Error 404")
        r = C.wait_for_registration(-1, fetcher=missing)
        self.assertFalse(r["found"])
        self.assertIn("404", r["last_error"])

    def test_main_ok_and_gate_failed(self):
        csvb = csv_text(rows_from(toy_tok(True))).encode()

        def fake(url, timeout=60):
            if url == C.PREREG_URL:
                return b"# CHIMP-PREREG"
            if url == C.ZENODO_URL:
                return b"not a zip"
            return csvb
        od = os.path.join(self.d, "out")
        rc = C.main(["--data-dir", os.path.join(self.d, "data"), "--out-dir", od,
                     "--wait-minutes", "0"], fetcher=fake, cfg=toy_cfg())
        self.assertEqual(rc, 0)
        with open(os.path.join(od, "chimp_results.json"), encoding="utf-8") as f:
            res = json.load(f)
        self.assertEqual(res["status"], "ok")
        self.assertEqual(res["overall"], "CONJUNCTION_FOUND")
        self.assertTrue(res["grading"]["truth"]["P1"])
        rc2 = C.main(["--data-dir", os.path.join(self.d, "data"), "--out-dir", od,
                      "--wait-minutes", "0"], fetcher=fake)           # real CFG: gate D fails
        self.assertEqual(rc2, 3)
        with open(os.path.join(od, "chimp_results.json"), encoding="utf-8") as f:
            res2 = json.load(f)
        self.assertEqual(res2["status"], "gate_failed")
        self.assertFalse(res2["grading"]["truth"]["P1"])


if __name__ == "__main__":
    unittest.main()
