"""Small-sample-safe(ish), transparent two-arm conversion experiment readout."""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import NormalDist


def load_assignments(path: Path, control: str, treatment: str) -> list[dict]:
    if control == treatment:
        raise ValueError("Les noms des variantes doivent être différents")
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"unit_id", "variant", "converted"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Colonnes requises: unit_id, variant, converted")
        rows, seen = [], set()
        for line, raw in enumerate(reader, 2):
            unit = (raw.get("unit_id") or "").strip()
            variant = (raw.get("variant") or "").strip()
            value = (raw.get("converted") or "").strip()
            if not unit:
                raise ValueError(f"Ligne {line}: unit_id vide")
            if unit in seen:
                raise ValueError(f"Ligne {line}: unité dupliquée {unit}")
            seen.add(unit)
            if variant not in (control, treatment):
                raise ValueError(f"Ligne {line}: variante inconnue {variant!r}")
            if value not in ("0", "1"):
                raise ValueError(f"Ligne {line}: converted doit valoir 0 ou 1")
            rows.append({"unit_id": unit, "variant": variant, "converted": int(value)})
    return rows


def wilson(successes: int, n: int, z: float) -> tuple[float, float]:
    if n <= 0:
        raise ValueError("Chaque variante doit contenir au moins une unité")
    p = successes / n
    z2 = z * z
    denom = 1 + z2 / n
    center = (p + z2 / (2 * n)) / denom
    radius = z * math.sqrt(p * (1 - p) / n + z2 / (4 * n * n)) / denom
    return max(0.0, center - radius), min(1.0, center + radius)


def analyze(rows: list[dict], control: str = "control", treatment: str = "treatment",
            expected_control_share: float = .5, confidence: float = .95,
            srm_alpha: float = .001) -> dict:
    if control == treatment:
        raise ValueError("Les noms des variantes doivent être différents")
    if not 0 < expected_control_share < 1:
        raise ValueError("expected_control_share doit être strictement entre 0 et 1")
    if not .5 < confidence < 1:
        raise ValueError("confidence doit être entre 0.5 et 1")
    if not 0 < srm_alpha < 1:
        raise ValueError("srm_alpha doit être strictement entre 0 et 1")
    grouped: dict[str, list[int]] = defaultdict(list)
    seen = set()
    for row in rows:
        unit, variant, converted = row["unit_id"], row["variant"], row["converted"]
        if unit in seen:
            raise ValueError(f"Unité dupliquée {unit}")
        seen.add(unit)
        if variant not in (control, treatment):
            raise ValueError(f"Variante inconnue {variant!r}")
        if converted not in (0, 1):
            raise ValueError("converted doit valoir 0 ou 1")
        grouped[variant].append(converted)
    if not grouped[control] or not grouped[treatment]:
        raise ValueError("Les deux variantes doivent avoir au moins une unité")

    nc, nt = len(grouped[control]), len(grouped[treatment])
    xc, xt = sum(grouped[control]), sum(grouped[treatment])
    pc, pt = xc / nc, xt / nt
    delta = pt - pc
    zcrit = NormalDist().inv_cdf((1 + confidence) / 2)
    lc, uc = wilson(xc, nc, zcrit)
    lt, ut = wilson(xt, nt, zcrit)
    # Newcombe hybrid score interval for a difference of independent proportions.
    ci_low = delta - math.sqrt((pt - lt) ** 2 + (uc - pc) ** 2)
    ci_high = delta + math.sqrt((ut - pt) ** 2 + (pc - lc) ** 2)

    pooled = (xc + xt) / (nc + nt)
    se = math.sqrt(pooled * (1 - pooled) * (1 / nc + 1 / nt))
    z_stat = delta / se if se else 0.0
    p_value = math.erfc(abs(z_stat) / math.sqrt(2))

    total = nc + nt
    expected_c = total * expected_control_share
    expected_t = total * (1 - expected_control_share)
    srm_chi2 = (nc - expected_c) ** 2 / expected_c + (nt - expected_t) ** 2 / expected_t
    srm_p = math.erfc(math.sqrt(srm_chi2 / 2))  # chi-square survival, df=1
    assignment_ok = srm_p >= srm_alpha
    if not assignment_ok:
        decision = "invalid_srm"
    elif ci_low > 0:
        decision = "positive"
    elif ci_high < 0:
        decision = "negative"
    else:
        decision = "inconclusive"

    return {
        "control": {"name": control, "assigned": nc, "conversions": xc, "conversion_rate": pc,
                    "wilson_interval": [lc, uc]},
        "treatment": {"name": treatment, "assigned": nt, "conversions": xt, "conversion_rate": pt,
                      "wilson_interval": [lt, ut]},
        "absolute_lift": delta, "relative_lift": (delta / pc if pc else None),
        "difference_interval": [max(-1.0, ci_low), min(1.0, ci_high)],
        "confidence": confidence, "pooled_z_statistic": z_stat, "pooled_z_p_value": p_value,
        "sample_ratio_check": {"expected_control_share": expected_control_share, "chi_square": srm_chi2,
                                "p_value": srm_p, "alpha": srm_alpha, "passed": assignment_ok},
        "decision": decision,
    }


def write_report(result: dict, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    (output / "readout.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    c, t = result["control"], result["treatment"]
    low, high = result["difference_interval"]
    lines = ["# Experiment readout", "",
        f"- Contrôle: **{c['conversions']}/{c['assigned']}** ({c['conversion_rate']:.2%}).",
        f"- Traitement: **{t['conversions']}/{t['assigned']}** ({t['conversion_rate']:.2%}).",
        f"- Lift absolu: **{result['absolute_lift']:+.2%}**; intervalle de différence à {result['confidence']:.0%}: [{low:+.2%}, {high:+.2%}].",
        f"- Test z bilatéral des proportions: p = {result['pooled_z_p_value']:.4g}.",
        f"- Contrôle du ratio d'affectation: p = {result['sample_ratio_check']['p_value']:.4g} "
        f"(seuil {result['sample_ratio_check']['alpha']:.3g}; {'OK' if result['sample_ratio_check']['passed'] else 'ÉCHEC'}).",
        f"- Décision: **{result['decision']}**.", "",
        "## Lecture", "",
        "La décision positive ou négative exige que l'intervalle de différence exclue zéro et que le contrôle du ratio d'affectation passe. "
        "Une conclusion `inconclusive` n'est pas une preuve d'absence d'effet. `invalid_srm` bloque l'interprétation.", "",
        "## Limites", "",
        "Le test z et l'intervalle supposent des unités indépendantes, une seule métrique binaire pré-déclarée et une analyse finale sans arrêt anticipé. "
        "Le test SRM détecte un déséquilibre global entre deux bras selon la part de contrôle configurée; il ne localise pas sa cause. "
        "Vérifiez aussi la randomisation, l'exposition, les métriques de garde et la durée de l'expérience avant toute décision produit.", ""]
    (output / "report.md").write_text("\n".join(lines), encoding="utf-8")


def generate_demo(path: Path, seed: int = 20261106) -> None:
    rng = random.Random(seed)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["unit_id", "variant", "converted"])
        for arm, probability in (("control", .10), ("treatment", .112)):
            for index in range(6000):
                writer.writerow([f"{arm}-{index:05d}", arm, int(rng.random() < probability)])


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyse finale reproductible d'un test A/B binaire à deux bras")
    subs = parser.add_subparsers(dest="command", required=True)
    demo = subs.add_parser("demo", help="génère un jeu de données fictif et son readout")
    demo.add_argument("--output", type=Path, default=Path("demo-output"))
    run = subs.add_parser("analyze", help="analyse un CSV d'affectations")
    run.add_argument("csv", type=Path); run.add_argument("--output", type=Path, default=Path("readout"))
    run.add_argument("--control", default="control"); run.add_argument("--treatment", default="treatment")
    run.add_argument("--expected-control-share", type=float, default=.5)
    run.add_argument("--confidence", type=float, default=.95)
    run.add_argument("--srm-alpha", type=float, default=.001)
    args = parser.parse_args()
    if args.command == "demo":
        args.output.mkdir(parents=True, exist_ok=True)
        data = args.output / "assignments.csv"
        generate_demo(data)
        result = analyze(load_assignments(data, "control", "treatment"))
    else:
        result = analyze(load_assignments(args.csv, args.control, args.treatment), args.control, args.treatment,
                         args.expected_control_share, args.confidence, args.srm_alpha)
    write_report(result, args.output)
    print(f"Décision: {result['decision']} | lift {result['absolute_lift']:+.2%} | p={result['pooled_z_p_value']:.4g} | rapport {args.output / 'report.md'}")


if __name__ == "__main__":
    main()
