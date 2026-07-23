# Phandango files

The files in this directory reproduce the group 2 and group 3 K-locus
visualisations associated with the published EC-K-Typing study. They are
based on the 90-locus reference set used for that analysis and are retained
unchanged for reproducibility.

These files are not the current EC-K-Typing database. The current Kaptive
database is provided by `EC-K-typing_group2and3.gbk` and
`EC-K-typing_group2and3.toml` in the repository root.

The historical Phandango metadata uses the locus identifiers current at the
time of the analysis. In database v4.0.0:

- KL110 was renamed KL51.
- KL116 was renamed KL97.
- KL124 was renamed KL6.

The `Order` column in `group2and3_phandango_meta.csv` identifies loci
classified as having atypical group 2 locus architecture in the published
90-locus dataset. For the current classification, see
[`DB/panaroo_refset/atypical_loci.csv`](../DB/panaroo_refset/atypical_loci.csv).

The directory contains:

- group 2 and group 3 gene presence–absence matrices;
- group 2 and group 3 midpoint-rooted neighbour-joining trees; and
- the combined metadata used to colour and label the visualisations.

Any future visualisation based on a later database release should be added
as a separately versioned dataset rather than overwriting these files.
