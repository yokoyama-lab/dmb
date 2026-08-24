#!/usr/bin/env python3
"""コマンドラインと JSON 入出力の検査。

自分の複体を持ち込めること（load_json）と，機械可読な出力（--json）が
壊れていないことを固定する。
"""

import contextlib
import io
import json
import os
import tempfile
import unittest

import complexes as X
import dmb_core as D


def run(argv):
    """dmb_core.main を走らせ，(終了コード, 標準出力) を返す。"""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        code = D.main(argv)
    return code, buf.getvalue()


class TestJsonRoundTrip(unittest.TestCase):
    def test_simplicial(self):
        for name, build in (("S^2", lambda: X.sphere(2)),
                            ("RP^2", X.projective_plane),
                            ("T(4,4)", lambda: D.torus(4, 4))):
            K = build()
            f = D.canonical_dmf(K)
            with self.subTest(name):
                K2, f2 = X.load_json(json.loads(X.dump_json(K, f)))
                self.assertEqual(sorted(K2.cells), sorted(K.cells))
                self.assertEqual(f2, f)
                self.assertEqual(D.MorseBott(K2, f2).report(),
                                 D.MorseBott(K, f).report())

    def test_cw(self):
        for build in (D.cw_sphere2_minimal, D.cw_torus_minimal,
                      D.cw_projective_plane_minimal):
            K = build()
            f = {c: K.dim(c) for c in K.cells}
            with self.subTest(K.name):
                K2, f2 = X.load_json(json.loads(X.dump_json(K, f)))
                self.assertEqual(K2.counts(), K.counts())
                self.assertEqual(K2.betti(), K.betti())
                self.assertEqual(K2.homology_z(), K.homology_z())
                self.assertEqual(f2, f)
                self.assertEqual(sorted(K2.faces), sorted(K.faces))
                self.assertEqual(K2.regular_pairs, K.regular_pairs)

    def test_regularity_survives(self):
        """regular / irregular の区別が往復で保たれる（(M1)(M3) に効く）。"""
        K = X.sphere(2)
        C = D.cw_from_simplicial(K)
        C2, _ = X.load_json(json.loads(X.dump_json(C)))
        for a in C.cells:
            for b in C.cells:
                self.assertEqual(C.is_regular(a, b), C2.is_regular(a, b))
        self.assertEqual(C2.irregular_above(C2.cells[0]), [])

    def test_writes_to_file(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "k.json")
            K = X.moebius()
            X.dump_json(K, D.constant_fn(K), path)
            K2, f2 = X.load_json(path)
            self.assertEqual(sorted(K2.cells), sorted(K.cells))
            self.assertEqual(set(f2.values()), {0})

    def test_rejects_bad_input(self):
        with self.assertRaises(ValueError):
            X.load_json({"type": "no-such"})
        with self.assertRaises(ValueError):          # 複体に無いセルの値
            X.load_json({"type": "simplicial", "facets": [[0, 1]], "f": {"9": 0}})
        with self.assertRaises(ValueError):          # 値が足りない
            X.load_json({"type": "simplicial", "facets": [[0, 1]], "f": {"0": 0}})


class TestCommandLine(unittest.TestCase):
    def test_default_report(self):
        code, out = run([])
        self.assertEqual(code, 0)
        self.assertIn("トーラス T(4, 4)", out)
        self.assertIn("すべての検査を通過", out)

    def test_positional_sizes(self):
        code, out = run(["5", "4"])
        self.assertEqual(code, 0)
        self.assertIn("トーラス T(5, 4)", out)

    def test_table(self):
        code, out = run(["--table"])
        self.assertEqual(code, 0)
        self.assertIn("inv.DMF", out)

    def test_json_output(self):
        code, out = run(["--ni", "4", "--nj", "4", "--json"])
        self.assertEqual(code, 0)
        payload = json.loads(out)
        self.assertEqual(payload["torus"], [4, 4])
        for name, r in payload["functions"].items():
            with self.subTest(name):
                self.assertIn("is_dmb", r)
                self.assertIn("P_K", r)
                self.assertEqual(r["P_K"], [1, 2, 1])
                self.assertTrue(r["is_dmb"])

    def test_field_option(self):
        code, out = run(["--ni", "4", "--nj", "4", "--field", "2", "--json"])
        self.assertEqual(code, 0)
        for r in json.loads(out)["functions"].values():
            self.assertEqual(r["field"], "F_2")

    def test_complex_from_json(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "rp2.json")
            K = D.cw_projective_plane_minimal()
            X.dump_json(K, {c: K.dim(c) for c in K.cells}, path)

            code, out = run(["--complex", path, "--json"])
            self.assertEqual(code, 0)
            r = json.loads(out)
            self.assertEqual(r["P_K"], [1])            # Q 上は 1
            self.assertEqual(r["M"], [1, 1, 1])
            self.assertEqual(r["R_M"], [0, 1])         # R(t) = t

            code, out = run(["--complex", path, "--field", "2", "--json"])
            r = json.loads(out)
            self.assertEqual(r["P_K"], [1, 1, 1])      # F_2 上は 1 + t + t^2
            self.assertEqual(D.poly_trim(r["R_M"]), [])
            self.assertTrue(r["sharp"])

    def test_complex_without_values_is_an_error(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "k.json")
            X.dump_json(X.sphere(2), None, path)
            with self.assertRaises(SystemExit):
                run(["--complex", path])

    def test_rejects_small_grid(self):
        with self.assertRaises(SystemExit):
            run(["2", "2"])

    def test_accepts_sys_argv_style(self):
        """main(sys.argv) の形（先頭がスクリプト名）でも動く。"""
        code, out = run(["dmb_core.py", "--table"])
        self.assertEqual(code, 0)
        self.assertIn("inv.DMF", out)


class TestReportJson(unittest.TestCase):
    def test_matches_the_object_api(self):
        K = D.torus(4, 4)
        f = D.height_fn(K, 4, 4)
        r = D.report_json(K, f)
        M = D.MorseBott(K, f)
        self.assertEqual(r["MB_sum"], M.morse_bott_polynomial())
        self.assertEqual(r["is_dmb"], M.is_dmb())
        self.assertEqual(r["critical"], len(M.critical()))
        self.assertTrue(r["sharp"])
        self.assertEqual(len(r["contributing_collections"]), 2)

    def test_is_json_serialisable(self):
        for build, f in ((lambda: D.torus(3, 3), None),
                         (D.cw_projective_plane_minimal, None)):
            K = build()
            g = f or {c: K.dim(c) for c in K.cells}
            json.dumps(D.report_json(K, g), ensure_ascii=False)


if __name__ == "__main__":
    unittest.main(verbosity=2)
