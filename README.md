# LikelihoodFreePhylogenetics (LFP)

Simulation et inférence sans vraisemblance (approche likelihood-free / ABC) pour
des modèles de mutation-sélection codons (mut-sel) à partir de chaînes MCMC de
PhyloBayes.

LFP peut :
- **simuler** des séquences codantes sous un modèle mut-sel en tirant au hasard
  les paramètres (priors) et des points de chaîne PhyloBayes, puis calculer des
  statistiques descriptives ;
- produire la **table de référence** (`.simu` : paramètres + statistiques)
  utilisée ensuite pour une inférence basée sur simulation (NRE / SNRE par ex.) ;
- faire de la **posterior prediction** (`PPred`) et de l'**ABC** par rejet.

---

## Compilation

```
cd sources
make              # ou make <cible>
```

Les binaires sont produits dans `data/`. Cibles principales : `LFP`,
`SimuMutSelPhylobayesMPI`, `ratesFromLeaves`, `codemlM7M8`, `codemlM7M8TpA`,
`codemlFMutSel`, `BayescodeSimuMutSel{AA,C,AAC}`. Il n'y a **pas** d'accélération
OpenMP significative sur la boucle de simulation : chaque point reste mono-fil.

---

## Modèles (`-m`)

- `AACodonMutSelPhylobayesMPI` — modèle codon mutation-sélection AA dont les
  chaînes sont des chaînes **PhyloBayes-MPI** d'emblème `AACODONMUTSELSBDP`
  (fichiers `*MSAAGTRW-A`). La sélection s'applique sur les **61 vecteurs de
  fitness codon** (61 codons sens du code utilisé, sans les stops). Le lecteur
  `readChainAACodonMutSelPhylobayesMPI` lit directement la `.chain`
  (arbre, nucp, nucrr, codonprofile, omega, Ncat, ssprofiles AA, alloc) —
  **pas de fichier `.pvalues`**.
  Variantes : `AACodonMutSelPhylobayesMPIABC` (ABC par rejet),
  `AACodonMutSelPhylobayesMPIPPred` (posterior prediction),
  `AACodonMutSelPhylobayesMPISeq` (association de séquences).
- `CodonMutSelFinite` — simulation depuis M[GTR+λCpG]-[1CatAA] ou
  M[GTR+λCpG]-[1CatCodon]. Variantes : `CodonMutSelFiniteABC`,
  `CodonMutSelFinitePPred`.
- `stats`, `show` — aide / contrôle.

---

## Fichier de contrôle (`*.conf`)

Directives lues (ordre libre dans le fichier) :

| Directive | Rôle |
|---|---|
| `#SUMMARIES  col1 col2 …` | Statistiques descriptives à calculer |
| `#ANCSUMMARIES` | Statistiques d'ancêtres (optionnel) |
| `#ACCSUMMARIES` | Statistiques supplémentaires du point (optionnel) |
| `#PARAM  chainID root  …` | Paramètres simulés en tête des `.simu` |
| `#SSMAP`, `#MAP`, `#ANCESTRALMAP` | Cartes des statistiques / substitutions extraites |
| `#DIST`, `#TRANS`, `#SPEUDODATA` | Options de sortie |
| `#SAMPLING  start every until` | Bornes de l'échantillonnage de la chaîne |
| `#NRUN  <N> <T>` | Dans `SimuMutSelPhylobayesMPI` : `N` simulats par point, au plus `T` simulats au total |
| `#NTHREADS  n` | Mono-fil en pratique : `1` recommandé |
| `#OUTPUT  /chemin/vers/base` | Préfixe des fichiers écrits (`+ .realdata`, `+ .simu`, `+ .ppp`) |
| `#LOCALPARAM  …` | Paramètres du modèle (voir ci-dessous) |

### `#LOCALPARAM` — principales options

```
-d <align>        alignement .puz            -chain <base_chain>   chaîne PhyloBayes (sans .chain)
-iscodon          données codantes           -code <Code>          Universal | MtMam | MtInv | …
-rootlength <r>   longueur racine
-freeroot <t1 t2> racine libre (2 taxa sortants)     -fixroot <t1 t2> <v> | -freek
-freelambdaTBL    libre (branches)           -freegtr | -freehky | -freegtnr   libre (modèle nuc)
-freelambdaomega  libre (ω)                  -freelambdaCpG        libre (hyper-mut. CpG)
-priorlambdaTBL log2Unif  -priorlambdaomega log2Unif  -priorlambdaCpG log10Unif
-start <s> -every <e> -until <u>   plage de points de chaîne
```

Les paramètres « multiplicatifs » : **λCpG** module le taux de transition en
contexte CpG, **λTBL** module toutes les longueurs de branches (taille de
l'arbre), **λω** module la valeur d'ω.

---

## Exemple 1 — générer une table de référence (driver)

```
data/SimuMutSelPhylobayesMPI -m AACodonMutSelPhylobayesMPI NHS-4353-MSAAGTRW-A-GTRWCpG.conf
```

### `NHS-4353-MSAAGTRW-A-GTRWCpG.conf`

```
#SUMMARIES  pwAC pwAG pwAT pwCG pwCT pwGT K80nuc pwSynAC … nucA nucC nucG nucT dinuc12AA … TTT TTC …
#PARAM  chainID root lambda_CpG lambda_TBL lambda_omega nucsA nucsC nucsG nucsT nucrrAC nucrrAG nucrrAT nucrrCG nucrrCT nucrrGT
#MAP  Nsub Nsynsub
#SAMPLING 600 4 999
#NRUN 100 10000
#NTHREADS 1
#OUTPUT /data_antv/NHS-4353-MSAAGTRW-A-GTRWCpG
#LOCALPARAM -iscodon -code Universal -freelambdaTBL -freeroot Echinops Procavia -freegtr \
            -freelambdaomega -rootlength 100 -freelambdaCpG \
            -priorlambdaomega log2Unif -priorlambdaCpG log10Unif -priorlambdaTBL log2Unif \
            -d NHS.puz -chain /step_1/NHS-4353-MSAAGTRW-A
```

Ici `#NRUN 100 10000` tire 100 réplicats sur chacun des ~100 points (600→999,
pas de 4) soit **10 000 simulats** couvrant toute la plage de chaîne.

Sorties :
- **`.realdata`** — statistiques observées (lues de la chaîne / de l'alignement) ;
- **`.simu`** — une ligne par simul : `chainID root <13 params> <stats descriptives> T_Nsub T_Nsynsub`
  (en-tête 41 colonnes ; lignes de données 42 champs, tabulation finale) ;
- **`.ppp`** — posterior predictive p-values.

### Assemblage vers une table SBI

Le `.simu` (colonnes dans l'ordre défini par le conf) se réordonne vers le layout
canonique `chainID root + 13 params + 132 stats` avec
`bioinfo/workflow/scripts/sbi_table_from_simu.py` du dépôt
`simulation-based-inference`.

---

## Exemple 2 — ABC par rejet (binaire `LFP`)

```
./LFP -m AACodonMutSelPhylobayesMPIABC ACO1-2673-MSAAGTRW-A-GTRWCpG.conf
```

---

## Licence

GPLv3 (voir l'en-tête des sources).
