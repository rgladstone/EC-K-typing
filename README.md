<img align="right" src="https://github.com/klebgenomics/Kaptive/blob/master/docs/assets/logo.png?raw=true" alt="Kaptive" width="200">

# _E. coli_ group 2 and group 3 capsular K-typing database

[![Streamlit App](https://img.shields.io/badge/Streamlit-%23FE4B4B.svg?logo=streamlit&logoColor=white)](https://kaptive-database-validator.streamlit.app/)
[![Release Database](https://github.com/rgladstone/EC-K-typing/actions/workflows/release.yml/badge.svg)](https://github.com/rgladstone/EC-K-typing/actions/workflows/release.yml)
[![DOI:10.1038/s41564-026-02283-w](https://zenodo.org/badge/DOI/10.1038/s41564-026-02283-w.svg)](https://doi.org/10.1038/s41564-026-02283-w)

## Database update summary: v4.0.0

Database v4.0.0 expands the group 2 and group 3 reference set from 90 to
109 K-locus records. It adds 16 sequence-defined loci, assigned the new
EC-K-Typing identifiers KL176–KL191 following comparison with the kTYPr
capsule catalogue. These loci do not currently have assigned K phenotypes.
A crosswalk between the EC-K-Typing identifiers, kTYPr catalogue labels,
and source genome accessions is provided in the
[v4 locus provenance](DB/panaroo_refset/v4_locus_provenance.md). The update
also adds three references with known phenotype associations: K18-K22
(KL18), K94 (KL94), and K95 (KL95).

Three existing locus records have been renamed to match their established
K phenotypes: the records previously designated KL124, KL110, and KL116
are now designated KL6, KL51, and KL97, respectively. K74 is excluded
because it is no longer recognised as a valid K serotype.

KL158 now has the descriptive type `Capsule null predicted`. When comparing
results generated with a database version before v4.0.0, we recommend
rerunning Kaptive with the current database or applying the nomenclature
crosswalk above.

## About the database

This Kaptive database supports _in silico_ typing of _E. coli_ group 2 and
group 3 capsular K loci, which use ABC transporter-dependent capsule
assembly systems <sup>[1]</sup>.

The reference set was developed by screening approximately 50,000
_E. coli_ genomes from bloodstream infection, human and animal carriage,
and collections spanning Europe, North America, Africa, and Asia. This was
supplemented with group 2 `kpsF`-positive and group 3 `kpsM`-positive
assemblies identified at 90% k-mer identity in a collection of 661,000
bacterial assemblies <sup>[2–12]</sup>.

The database is expected to cover most group 2 and group 3 loci found among
invasive _E. coli_. These capsule groups are especially prevalent in
phylogroups B2 and D and occur less frequently in phylogroups A, B1, and C.

The database includes all currently recognised group 2 and group 3 K
phenotypes with an established locus association: 33 K phenotypes
represented by 30 locus records. Three phenotype pairs cannot currently
be distinguished reliably from their K-locus sequences and are reported
as composite types: K13-K23 (KL13), K18-K22 (KL18), and K96-K54 (KL96).
KL158 additionally has the descriptive type `Capsule null predicted`.

### Atypical loci and reference boundaries

Some K loci have atypical architectures in which capsule-specific genes
occur outside the conserved export regions. A current list of loci
classified as having atypical architecture is provided in the
[atypical-locus metadata](DB/panaroo_refset/atypical_loci.csv).

A close Kaptive match may fail to capture additional capsule genes located
beyond the reference boundaries. For these loci, inspect the assembly
context and gene-level Kaptive results. See the
[Kaptive documentation](https://kaptive.readthedocs.io/en/latest/) for
general guidance on interpreting results.

### Database construction

Reference genes were annotated using Bakta 1.10.x and reconciled using
Panaroo 1.5.2. The final reference GFFs were generated using Panaroo 1.6.0
<sup>[13,14]</sup>.

Non-capsular IS-element-associated annotations and newly assigned
hypothetical CDS annotations shorter than 200 bp were excluded from the
final gene annotations. Their underlying nucleotide sequences were not
removed from the locus references.

The reproducible database-construction workflow and curated reference GFFs
are provided in [`DB/panaroo_refset`](DB/panaroo_refset/).

### Database additions and corrections

To have a novel locus considered for inclusion in the public database,
contact rebeccgl@uio.no or open an issue in this repository.

Requests to add or correct a K-phenotype association for an existing locus
are also welcome. Please provide the locus designation, proposed phenotype,
reference strain or sequence accession where available, and the supporting
serological, experimental, or published evidence. Phenotype associations
will only be added when their provenance and relationship to the locus can
be evaluated.

For sequences that cannot yet be shared, the containerised
[private database builder](private_builder/README.md) can add one local
locus to a frozen copy of the public reference set and produce a validated
private Kaptive database. The sequence is processed locally, the repository
is not modified, and nothing is uploaded. Local identifiers of KL9000 or
greater are deliberately kept separate from official EC-K-Typing locus
assignments.

## Using the database

The database is formatted for Kaptive 3 and has been validated using
Kaptive 3.1.0. After cloning the repository or downloading a release, run:

```bash
kaptive assembly EC-K-typing_group2and3.gbk your_assembly.fasta
```

See the [Kaptive documentation](https://kaptive.readthedocs.io/en/latest/)
for installation, output interpretation, and additional options.

## Citation

If you use this database, please cite:

Gladstone, R. A., Pesonen, M., Pöntinen, A. K. _et al._
Identification of transporter-dependent capsular loci associated with the
invasive potential of _Escherichia coli_. _Nature Microbiology_ **11**,
1205–1216 (2026).
https://doi.org/10.1038/s41564-026-02283-w

If you use Kaptive, please also cite:

Stanton, T. D., Hetland, M. A. K., Löhr, I. H., Holt, K. E. & Wyres, K. L.
Fast and accurate in silico antigen typing with Kaptive 3.
_Microbial Genomics_ **11**, 001428 (2025).
https://doi.org/10.1099/mgen.0.001428

## References

1. Stanton TD, Hetland MAK, Löhr IH, Holt KE, Wyres KL. Fast and accurate
   in silico antigen typing with Kaptive 3. _Microbial Genomics_.
   2025;11:001428.
2. Gladstone RA, McNally A, Pöntinen AK, Tonkin-Hill G, Lees JA, et al.
   Emergence and dissemination of antimicrobial resistance in
   _Escherichia coli_ causing bloodstream infections in Norway in 2002–17:
   a nationwide, longitudinal, microbial population genomic study.
   _The Lancet Microbe_. 2021;2:e331–e341.
3. Arredondo-Alonso S, Pöntinen AK, Gama JA, Gladstone RA, Harms K, et al.
   Plasmid-driven strategies for clone success in _Escherichia coli_.
   _Nature Communications_. 2025;16:2921.
4. Kallonen T, Brodrick HJ, Harris SR, Corander J, Brown NM, et al.
   Systematic longitudinal survey of invasive _Escherichia coli_ in England
   demonstrates a stable population structure only transiently disturbed
   by the emergence of ST131. _Genome Research_. 2017;27:1437–1449.
5. Pöntinen AK, Gladstone RA, Pesonen H, Pesonen M, Cléon F, et al.
   Modulation of multidrug-resistant clone success in _Escherichia coli_
   populations: a longitudinal, multi-country, genomic and antibiotic
   usage cohort study. _The Lancet Microbe_. 2024;5:e142–e150.
6. Shao Y, Garcia-Mauriño C, Clare S, Dawson NJR, Mu A, et al. Primary
   succession of Bifidobacteria drives pathogen resistance in neonatal
   microbiota assembly. _Nature Microbiology_. 2024;9:2570–2582.
7. Mäklin T, Thorpe HA, Pöntinen AK, Gladstone RA, Shao Y, et al. Strong
   pathogen competition in neonatal gut colonisation.
   _Nature Communications_. 2022;13:7417.
8. Liu CM, Aziz M, Park DE, Wu Z, Stegger M, et al. Using
   source-associated mobile genetic elements to identify zoonotic
   extraintestinal _E. coli_ infections. _One Health_. 2023;16:100518.
9. Ludden C, Raven KE, Jamrozy D, Gouliouris T, Blane B, et al. One Health
   genomic surveillance of _Escherichia coli_ demonstrates distinct
   lineages and mobile genetic elements in isolates from humans versus
   livestock. _mBio_. 2019;10:e02693-18.
10. Sands K, Carvalho MJ, Portal E, Thomson K, Dyer C, et al.
    Characterization of antimicrobial-resistant Gram-negative bacteria
    that cause neonatal sepsis in seven low- and middle-income countries.
    _Nature Microbiology_. 2021;6:512–523.
11. Dicks J, Fazal M-A, Oliver K, Grayson NE, Turnbull JD, et al.
    NCTC3000: a century of bacterial strain collecting leads to a rich
    genomic data resource. _Microbial Genomics_. 2023;9:mgen000976.
12. Blackwell GA, Hunt M, Malone KM, Lima L, Horesh G, et al. Exploring
    bacterial diversity via a curated and searchable snapshot of archived
    DNA sequences. _PLoS Biology_. 2021;19:e3001421.
13. Schwengers O, Jelonek L, Dieckmann MA, Beyvers S, Blom J, Goesmann A.
    Bakta: rapid and standardized annotation of bacterial genomes via
    alignment-free sequence identification. _Microbial Genomics_.
    2021;7:000685.
14. Tonkin-Hill G, MacAlasdair N, Ruis C, Weimann A, Horesh G, et al.
    Producing polished prokaryotic pangenomes with the Panaroo pipeline.
    _Genome Biology_. 2020;21:180.
