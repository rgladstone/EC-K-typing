# Changelog: EC-K-typing_group2and3

## [v4.0.0] - 2026-07-23
* Expanded the database from 90 to 109 group 2 and group 3 K-locus records.
* Added genotype loci KL176-KL191.
* Added phenotype-associated reference loci K18-K22 (KL18), K94 (KL94), and K95 (KL95).
* Renamed KL124 to KL6, KL110 to KL51, and KL116 to KL97 to match their established K phenotypes.
* Added complete phenotype metadata for K6, K51, K97, K18-K22, K94, and K95.
* Added the unconditional `Capsule null predicted` phenotype for KL158.
* Preserved complete composite and descriptive phenotype values in Kaptive by formatting GenBank type notes as `K type: VALUE`.
* Excluded the deprecated K74 designation and did not import kTYPr allele names.
* Reconciled annotations across all loci with Bakta 1.10.x, Panaroo 1.5.2, and panaroo-generate-gffs 1.6.0.
* Annotated the divergent terminal transport clusters in KL183, KL186, and KL187 as KpsT and KpsM variants.
* Added current atypical-order locus metadata and documented the Phandango files as frozen publication artifacts.
* Removed IS-element annotations and newly assigned hypothetical CDS annotations shorter than 200 bp without altering locus sequences.
* Added reproducible scripts, environments, metadata mappings, and validation guidance for future database updates.
* Added a containerised, offline private-database builder that leaves the public repository unchanged and validates every derivative reference by self-typing.
* Confirmed 109/109 perfect reference self-matches with Kaptive 3.1.0, with no untypeable references.
* Added a provenance crosswalk between the v4 EC-K-Typing loci, kTYPr catalogue labels, and source genome accessions.

## [v3.0.1] - 2026-07-02
* Update keyword for EC-K typing configuration (9cf1146)
* Update genbank filename in TOML configuration (cde6ea7)
* Update EC-K-typing_group2and3.gbk (e8397bb)
* Rename EC-K-typing_group2and3_v3.0.0.gbk to EC-K-typing_group2and3.gbk (b6595e0)
* Update name and pathway in EC-K-typing config (ac335f5)
