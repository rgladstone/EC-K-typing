# Main database

EC-K-typing_group2and3_v1_170724.gbk is the main database, it is an _E. coli_ Group 2 and Group 3 K-typing database formatted for use with the tool Kaptive version 3 https://kaptive.readthedocs.io/en/latest/

Phenotypes, if known, are given in K-type and serotype fields. 18/26 group 2 phenotypes are represented, and 5/6 group 3 phenotypes are represented.

K-antigen loci (with and without phenotype). This database is expected to cover the majority of invasive isolates; group 2 and 3 capsule prevalence is near complete in phylogroups B2 and D, whilst phylogroups A, B1, and C have a lower prevalence. ~50,000 E. coli genomes were screened from Bloodstream Infections (human), carriage (human and animal), from Europe, Africa and Asia, and all kpsF-positive assemblies in Blackwell _et al_ PLOS 2018.

## Skeleton database
EC-K-typing_group1and4_v1_180724.gbk is a skeleton Group 1 and Group 4 K-typing database formatted for use with the tool Kaptive version 3 https://kaptive.readthedocs.io/en/latest/
Phenotypes, if known, are given in K-type and serotype fields. Only 3 reference loci with phenotypes are represented, in group 1 and group 4 (K38, K84, K102).

Annotated using prokka 1.14.5 and panaroo 3.1.4 from https://github.com/tseemann/prokka and https://github.com/gtonkinhill/panaroo to give consistent annotation across the DB.

Currently, IS-element-associated annotations, but not sequence, have been removed. Kaptive V3 is, therefore, required for this version, and length discrepancies will result when the IS-element is not present in the query assembly.
Other known non-capsule genes (beta-lactamase) annotations within the K-locus have also been removed, but not the sequence.
