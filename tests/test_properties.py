#!/usr/bin/env python3
"""定義・定理が満たすべき性質を，多数のランダムな例で確かめる（property-based test）。

個々の値を決め打ちする回帰テストと違い，「どんな入力でも成り立つはずのこと」を
たくさんの例で殴る。乱数は固定シードなので再現する（PYTHONHASHSEED にも依存しない）。

    python3 -m unittest test_properties -v
    DMB_TRIALS=50 python3 -m unittest test_properties   # 試行回数を増やす
"""

import os
import random
import unittest

import complexes as X
import dmb_core as D

TRIALS = int(os.environ.get("DMB_TRIALS", "8"))

# 小さめの複体を使う（ベッチ数の計算が O(n^3) なので）
SPACES = [
    ("S^1(5)", lambda: X.circle(5)),
    ("Δ^2", lambda: X.simplex(2)),
    ("S^2", lambda: X.sphere(2)),
    ("メビウス", X.moebius),
    ("円筒(3,2)", lambda: X.annulus(3, 2)),
    ("RP^2", X.projective_plane),
    ("T(3,3)", lambda: X.torus(3, 3)),
    ("S^1⊔S^1", lambda: X.two_circles(3)),
]


GENERATORS = [D.random_function, D.random_max_extension, D.random_matching_dmf]
hash_seed = D.hash_seed


def instances(seed0=0):
    """(名前, 複体, 関数) を順に返す。"""
    for name, build in SPACES:
        K = build()
        for gen in GENERATORS:
            for t in range(TRIALS):
                rng = random.Random(hash_seed(name, gen.__name__, t, seed0))
                yield f"{name}/{gen.__name__}/{t}", K, gen(K, rng)


# ------------------------------------------------------------------ 不変条件


class TestStructuralInvariants(unittest.TestCase):
    """f が何であっても成り立つべきこと。"""

    def test_collections_partition_the_complex(self):
        for label, K, f in instances():
            with self.subTest(label):
                M = D.MorseBott(K, f)
                Ls = M.collections()
                seen = [c for L in Ls for c in L]
                self.assertEqual(sorted(seen), sorted(K.cells))
                self.assertEqual(len(seen), len(set(seen)))
                for L in Ls:
                    self.assertEqual(len({f[c] for c in L}), 1)

    def test_collections_are_r_path_connected_and_maximal(self):
        for label, K, f in instances():
            with self.subTest(label):
                M = D.MorseBott(K, f)
                for L in M.collections():
                    members = set(L)
                    # 連結性: 1 点から facet 関係（同値）で全体に届く
                    seen, stack = {L[0]}, [L[0]]
                    while stack:
                        c = stack.pop()
                        for d in K.above[c] + K.below[c]:
                            if d in members and d not in seen and f[d] == f[c]:
                                seen.add(d)
                                stack.append(d)
                    self.assertEqual(seen, members)
                    # 極大性: 外に，同じ値で余次元 1 で隣接するセルはない
                    for c in L:
                        for d in K.above[c] + K.below[c]:
                            if f[d] == f[c]:
                                self.assertIn(d, members)

    def test_critical_implies_weakly_critical(self):
        for label, K, f in instances():
            with self.subTest(label):
                M = D.MorseBott(K, f)
                self.assertTrue(set(M.critical()) <= set(M.weakly_critical()))

    def test_arrows_are_exactly_the_snc_pairs(self):
        for label, K, f in instances():
            with self.subTest(label):
                M = D.MorseBott(K, f)
                pairs = {(s, t) for s in K.cells for t in K.above[s] if f[t] < f[s]}
                self.assertEqual(set(M.arrows()), pairs)
                # DMBF ⇔ 各セルの上向き・下向きの矢印が高々 1 本
                ok = all(len(M.up_snc(c)) <= 1 and len(M.down_snc(c)) <= 1 for c in K.cells)
                self.assertEqual(ok, M.is_dmb())

    def test_reduced_collections_are_the_weakly_critical_part(self):
        for label, K, f in instances():
            with self.subTest(label):
                M = D.MorseBott(K, f)
                wc = set(M.weakly_critical())
                got = {frozenset(C) for C in M.reduced_collections()}
                want = {frozenset(c for c in L if c in wc)
                        for L in M.collections()
                        if any(c in wc for c in L)}
                self.assertEqual(got, want)

    def test_morsification_critical_cells(self):
        """Morsification の臨界セル = ⋃ reduced collection（f が DMBF でなくても成立）。"""
        for label, K, f in instances():
            with self.subTest(label):
                M = D.MorseBott(K, f)
                Y = D.MorseBott(K, D.morsify(K, f))
                self.assertEqual(set(Y.critical()),
                                 {c for C in M.reduced_collections() for c in C})

    def test_deterministic(self):
        """同じ入力なら何度計算しても同じ答え。"""
        for label, K, f in instances():
            with self.subTest(label):
                a = D.MorseBott(K, f).report()
                b = D.MorseBott(K, dict(f)).report()
                for key in ("is_dmb", "is_dmf", "P_K", "MB_sum", "R_MB", "collections"):
                    self.assertEqual(a[key], b[key])


class TestInvarianceUnderReparametrisation(unittest.TestCase):
    """理論は f の値の「順序」しか見ないはずである。"""

    def test_strictly_increasing_reparametrisation(self):
        for label, K, f in instances():
            with self.subTest(label):
                for g in (lambda x: 3 * x + 1, lambda x: x - 7, lambda x: 2 ** x):
                    a = D.MorseBott(K, f)
                    b = D.MorseBott(K, {c: g(v) for c, v in f.items()})
                    self.assertEqual(a.is_dmb(), b.is_dmb())
                    self.assertEqual(a.is_dmf(), b.is_dmf())
                    self.assertEqual(sorted(map(sorted, a.collections())),
                                     sorted(map(sorted, b.collections())))
                    self.assertEqual(a.weakly_critical(), b.weakly_critical())
                    self.assertEqual(a.morse_bott_polynomial(), b.morse_bott_polynomial())

    def test_relabelling_vertices(self):
        """頂点に付け替え（複体の同型）をしても結論は変わらない。"""
        rng = random.Random(hash_seed("relabel"))
        for label, K, f in instances():
            with self.subTest(label):
                verts = [v[0] for v in K.cells_of_dim(0)]
                perm = list(verts)
                rng.shuffle(perm)
                sub = dict(zip(verts, perm))
                K2 = D.Complex([tuple(sub[v] for v in c) for c in K.cells_of_dim(2)]
                               or [tuple(sub[v] for v in c) for c in K.cells_of_dim(1)])
                if len(K2.cells) != len(K.cells):
                    continue                      # 1 次元と 2 次元が混ざる複体は飛ばす
                f2 = {tuple(sorted(sub[v] for v in c)): val for c, val in f.items()}
                a, b = D.MorseBott(K, f).report(), D.MorseBott(K2, f2).report()
                for key in ("is_dmb", "is_dmf", "P_K", "MB_sum", "R_MB", "collections"):
                    self.assertEqual(a[key], b[key], key)


class TestTheorems(unittest.TestCase):
    """論文の定理が成り立つこと。"""

    def test_theorem_412_whenever_dmb(self):
        """DMBF なら Σ_C P_t(C) = P_t(K) + (1+t)R(t) かつ R(t) ≥ 0。"""
        checked = 0
        for label, K, f in instances():
            M = D.MorseBott(K, f)
            if not M.is_dmb():
                continue
            with self.subTest(label):
                r = M.report()
                self.assertIsNotNone(r["R_MB"], label)
                self.assertTrue(all(c >= 0 for c in r["R_MB"]), (label, r["R_MB"]))
                checked += 1
        self.assertGreater(checked, 20, "DMBF の例が少なすぎる")

    def test_euler_identity_whenever_dmb(self):
        """χ(K) = Σ_{σ ∈ ⋃C_red} (-1)^{dim σ}（Lean の IsDMB.euler_eq_reduced）。"""
        for label, K, f in instances():
            M = D.MorseBott(K, f)
            if not M.is_dmb():
                continue
            with self.subTest(label):
                cells = {c for C in M.reduced_collections() for c in C}
                self.assertEqual(X.euler(K), sum((-1) ** K.dim(c) for c in cells))

    def test_morsification_gives_a_dmf(self):
        for label, K, f in instances():
            M = D.MorseBott(K, f)
            if not M.is_dmb():
                continue
            with self.subTest(label):
                self.assertTrue(D.MorseBott(K, D.morsify(K, f)).is_dmf())

    def test_dmf_implies_dmb(self):
        """Lemma 3.1（修正版）: 単体的複体では DMF ⇒ DMBF。"""
        checked = 0
        for label, K, f in instances():
            M = D.MorseBott(K, f)
            if not M.is_dmf():
                continue
            with self.subTest(label):
                self.assertTrue(M.is_dmb(), label)
                checked += 1
        self.assertGreater(checked, 20, "DMF の例が少なすぎる")

    def test_forman_lemma_u_plus_d(self):
        """Forman Lemma 2.5: DMF なら U(σ) + D(σ) ≤ 1。"""
        for label, K, f in instances():
            M = D.MorseBott(K, f)
            if not M.is_dmf():
                continue
            with self.subTest(label):
                for c in K.cells:
                    self.assertLessEqual(len(M.up_nc(c)) + len(M.down_nc(c)), 1, (label, c))

    def test_theorem_32_both_directions(self):
        """Theorem 3.2（両方向）: f が DMF ⇔ f が DMBF で，かつどの collection L も
            (i) #L = 1，または (ii) #L = 2 かつ L が reduced
        を満たす。

        (i) では reduced を要求しないところが非対称で，実際に必要:
        単射な DMF では matching の相手と狭義の大小がつくので，
        singleton の collection は weakly critical とは限らない。"""
        for label, K, f in instances():
            with self.subTest(label):
                M = D.MorseBott(K, f)
                wc = set(M.weakly_critical())
                rhs = (M.is_dmb()
                       and all(len(L) == 1
                               or (len(L) == 2 and all(c in wc for c in L))
                               for L in M.collections()))
                self.assertEqual(M.is_dmf(), rhs, label)

    def test_theorem_32_reducedness_of_pairs_cannot_be_dropped(self):
        """(ii) の「reduced」は落とせない。

        パス v0 - v1 - v2 上で f(v0)=0, f(v1)=1, f(v2)=0, f(e01)=1, f(e12)=0 と置くと
        L = {v1, e01} は #L = 2 だが v1 は U^snc(v1) = 1 で weakly critical でない。
        このとき f は DMBF だが DMF ではない（U(v1) = 2）。"""
        K = D.Complex([(0, 1), (1, 2)])
        v0, v1, v2 = (0,), (1,), (2,)
        e01, e12 = (0, 1), (1, 2)
        M = D.MorseBott(K, {v0: 0, v1: 1, v2: 0, e01: 1, e12: 0})
        self.assertTrue(M.is_dmb())
        self.assertFalse(M.is_dmf())
        self.assertIn("M2", {k for k, _, _ in M.dmf_violations()})
        pair = [L for L in M.collections() if len(L) == 2 and v1 in L]
        self.assertEqual(len(pair), 1)
        self.assertNotIn(v1, M.weakly_critical())   # だから reduced ではない

    def test_theorem_32_pair_that_is_reduced_gives_a_dmf(self):
        """逆に #L = 2 で reduced なら DMF になる（matching の対を同じ値に潰した形）。"""
        K = D.Complex([(0, 1)])
        M = D.MorseBott(K, {(0,): 0, (1,): 1, (0, 1): 1})
        self.assertTrue(M.is_dmb())
        self.assertTrue(M.is_dmf())
        self.assertEqual(sorted(len(L) for L in M.collections()), [1, 2])
        self.assertEqual(len(M.weakly_critical()), len(K.cells))

    def test_strong_morse_inequalities(self):
        """DMF の強いモース不等式と，DMBF 版の強い不等式。"""
        for label, K, f in instances():
            M = D.MorseBott(K, f)
            if not M.is_dmb():
                continue
            with self.subTest(label):
                p = M.report()["P_K"]
                for poly in ([M.morse_polynomial()] if M.is_dmf() else []) + \
                            [M.morse_bott_polynomial()]:
                    top = max(len(poly), len(p))
                    for k in range(top):
                        lhs = sum((-1) ** (k - i) * (poly[i] if i < len(poly) else 0)
                                  for i in range(k + 1))
                        rhs = sum((-1) ** (k - i) * (p[i] if i < len(p) else 0)
                                  for i in range(k + 1))
                        self.assertGreaterEqual(lhs, rhs, (label, k))
                    # 最上位では等号（オイラー標数）
                    self.assertEqual(
                        sum((-1) ** i * (poly[i] if i < len(poly) else 0) for i in range(top)),
                        sum((-1) ** i * (p[i] if i < len(p) else 0) for i in range(top)))


class TestPolynomialHelpers(unittest.TestCase):
    def test_div_by_one_plus_t_round_trip(self):
        rng = random.Random(hash_seed("poly"))
        for _ in range(200):
            q = [rng.randrange(-5, 6) for _ in range(rng.randrange(1, 6))]
            p = D.poly_add(q, [0] + q)                     # p = (1+t) q
            self.assertEqual(D.poly_trim(D.poly_div_1pt(p)), D.poly_trim(q))

    def test_div_rejects_non_multiples(self):
        rng = random.Random(hash_seed("poly2"))
        rejected = 0
        for _ in range(200):
            p = [rng.randrange(-3, 4) for _ in range(rng.randrange(1, 5))]
            q = D.poly_div_1pt(p)
            if q is None:
                rejected += 1
            else:
                self.assertEqual(D.poly_trim(D.poly_add(q, [0] + q)), D.poly_trim(p))
        self.assertGreater(rejected, 0, "割り切れない例が 1 つも出ていない")


if __name__ == "__main__":
    unittest.main(verbosity=2)
