# Boyd's conductor-11 Mahler measure conjecture — certification repository

This repository accompanies the paper

> **Boyd's conductor-11 Mahler measure conjecture: proof of the
> split-integral identity (C3), with an exact structural analysis of the
> family S_k**, Huimin Zheng (2026).

It contains the full computational certification archive for the proof:
all scripts, exact-arithmetic certificates, referee reports and responses,
and the manuscript sources.

## Layout

- `report/paper.tex` — the manuscript (compiles with `pdflatex`, 34 pp).
- `report/main.tex` — a detailed Chinese companion report (`xelatex`, 27 pp).
- `report/referee_report_round*.md`, `report/response_round*.md` — the
  complete adversarial-review trail (six rounds) and the responses.
- `code/` — all certification scripts (Python 3 + mpmath/sympy/python-flint,
  and PARI/GP). The key scripts print explicit PASS/FAIL summaries and are
  listed with their purposes in Appendix A of the paper.
- `notes/` — frozen certificates and archives (e.g.
  `notes/attack10-interval.txt` for the Arb period certification,
  `notes/attack16-siegel-anchor.txt` and `notes/attack17-*.txt` for the
  Siegel-unit regulator anchor).
- `requirements.txt` — Python dependencies.

## Reproduction

```sh
python -m pip install -r requirements.txt
cd code && python b11.py && python attack1.py && python attack2.py \
  && python attack3.py && python torsion.py && python endpoint_torsion3.py \
  && python kappa_exact.py \
  && python boundary_torsion.py && python closedness_check.py \
  && python ntilde_family.py && python b_family.py && python winding.py \
  && python dilog.py && python k53_attack.py && python kneg_m.py \
  && python n1_certify.py
python n1_interval.py       # Arb certification, k=0
python branch_certify.py    # branch assignments + modulus ordering
python sign_certify.py      # all-Arb sign certificate, k=0
python k1_interval.py       # Arb certification, conductor 17
python k1_branch_certify.py # k=1 branch assignments + ordering
python k1_sign_certify.py   # all-Arb sign certificate, k=1
python k1_certify.py && python k1_diamond.py && python k1_dilog.py
gp -q verify_family.gp && gp -q verify_ratios.gp
gp -q winding.gp && gp -q dilog.gp
gp -q k53.gp && gp -q k53b.gp && gp -q kfamily_torsion.gp
gp -q k1_pari.gp && gp -q k1_points.gp && gp -q k1_zvals.gp
gp -q verify_coinvariant.gp  # quotient vs subgroup generator integrals
python siegel_anchor_step11.py  # Siegel anchor: final value (Thm. anchor)
python siegel_anchor_step12.py  # membership M_2(Gamma_1(121)) + Sturm
python siegel_anchor_step13.py  # primitivity of the cycle (Smith nf)
python siegel_anchor_step14.py  # exact constants, D_U/D_V, cusp table
# full anchor chain rebuild: siegel_anchor_step4.py -> step5.py -> step6.py
#   -> step7.gp -> step9.py -> step8.py   (step8 takes 10-20 min)
#   (step10 is a discarded defective experiment and is not part of the chain)
```

The revision history is preserved in git tags `rev2`–`rev9`, matching the
referee rounds documented in `report/`.

## Literature

Third-party papers used in the proof are not redistributed here; they are
publicly available (Brunault, *Regulators of Siegel units and applications*,
J. Number Theory 163 (2016), arXiv:1504.08127; Brunault's thesis, 2005;
Bertin, CRM Proc. 36 (2004); Lalín–Ramamonjisoa, IJNT 13 (2017); and
Boyd's PNWNT 2015 talk slides). See the bibliography of the paper and
`notes/literature-notes.md`.

## Declaration on the use of AI tools

As declared in the paper, the research (computation, proof-strategy
exploration, certified verification, manuscript preparation) was carried
out by the author with the assistance of the AI system Kimi (Moonshot AI);
all mathematical content has been checked and verified by the author, who
takes full responsibility.

## License

MIT (see `LICENSE`). Contact: Huimin Zheng, zhhm@ahstu.edu.cn.
