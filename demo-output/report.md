# Experiment readout

- Contrôle: **572/6000** (9.53%).
- Traitement: **676/6000** (11.27%).
- Lift absolu: **+1.73%**; intervalle de différence à 95%: [+0.64%, +2.83%].
- Test z bilatéral des proportions: p = 0.00187.
- Contrôle du ratio d'affectation: p = 1 (seuil 0.001; OK).
- Décision: **positive**.

## Lecture

La décision positive ou négative exige que l'intervalle de différence exclue zéro et que le contrôle du ratio d'affectation passe. Une conclusion `inconclusive` n'est pas une preuve d'absence d'effet. `invalid_srm` bloque l'interprétation.

## Limites

Le test z et l'intervalle supposent des unités indépendantes, une seule métrique binaire pré-déclarée et une analyse finale sans arrêt anticipé. Le test SRM détecte un déséquilibre global entre deux bras selon la part de contrôle configurée; il ne localise pas sa cause. Vérifiez aussi la randomisation, l'exposition, les métriques de garde et la durée de l'expérience avant toute décision produit.
