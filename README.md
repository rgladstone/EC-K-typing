# Main Group 2 and Group 3 database

EC-K-typing_group2and3_v1_170724.gbk is a _E. coli_ Group 2 and Group 3 (ABC-transporter dependant) K-typing database formatted for use with the tool Kaptive version 3 https://kaptive.readthedocs.io/en/latest/

Phenotypes, if known, are given in the genbank K-type fields, when the phenotype is not known and the locus shares multiple region 2 genes with a known phenotype K-like is denoted in the K-type field. According to the classifcation reported on https://www.iith.ac.in/EK3D/classification.php, 18/26 group 2 phenotypes are represented, and 7/7 group 3 phenotypes are represented in a database of 101 loci according the classfication here https://www.iith.ac.in/EK3D/classification.php. These loci represent unique gene presence-absence patterns excluding IS-elements.

K-antigen loci (with the ex cept of K19 which we reclassified as group 3). This database is expected to cover the majority of invasive _E. coli_ isolates; group 2 and 3 capsule prevalence is near complete in phylogroups B2 and D, whilst phylogroups A, B1, and C have a lower prevalence. ~50,000 _E. coli_ genomes were screened from bloodstream infections (human), carriage (human and animal), from Europe, North America, Africa and Asia, and all kpsF-positive assemblies (90% ID) in Blackwell _et al_ PLOS 2018.

Some capsules have the same gene presence absence pattern and differ only in sequence, the genetic determinants need validation before these are included as phenotypic logic for Kaptive to distinguish these.
###	K96 and K54
1 simple non-synomous change in K96-07 and one 1 non-synomous change that results in missense for 9 AA and a truncation of 8 AA in K96-07 an alginate o-acetyltransferase.
###	K13 and K23
These capsular types shared the same genes and are highly conserved, they have previosuly been reported to belong to a serogroup along with K20 (which has a slighly different gene set) https://journals.asm.org/doi/10.1128/iai.39.2.623-629.1983
2 non-synomous changes in vatD acetyltransferase.
2 insertions of 10 and 12bp in a IS-element and 1 non-synomous change.
###	K1 and K92
Despite having the same gene content these two capsular types have divergent sequences allowing them to be typed sepearately, both reference loci are included in the database. There is a particularly high density of changes in neuS resulting in 60 NS AA changes, known to be the genetic determinant.

Currently, IS-element-associated annotations, but not sequence, have been removed, 33 of the G2 and G3 K-loci have an IS and of those 17 have TIRs in capsule gene CDS, and 9 have putative capsular genes as cargo.
Another known non-capsule gene (a beta-lactamase) annotation within KL143 has also been removed, but not the sequence.

The K-loci genes were Annotated using prokka 1.14.5 and panaroo 3.1.4 from https://github.com/tseemann/prokka and https://github.com/gtonkinhill/panaroo to give consistent annotation across the DB.

# References
1) Blackwell, Grace A., Martin Hunt, Kerri M. Malone, Leandro Lima, Gal Horesh, Blaise T. F. Alako, Nicholas R. Thomson, and Zamin Iqbal. 2021. “Exploring Bacterial Diversity via a Curated and Searchable Snapshot of Archived DNA Sequences.” PLoS Biology 19 (11): e3001421. https://doi.org/10.1371/journal.pbio.3001421.
2) Gladstone, Rebecca A., Alan McNally, Anna K. Pöntinen, Gerry Tonkin-Hill, John A. Lees, Kusti Skytén, François Cléon, et al. 2021. “Emergence and Dissemination of Antimicrobial Resistance in Escherichia Coli Causing Bloodstream Infections in Norway in 2002–17: A Nationwide, Longitudinal, Microbial Population Genomic Study.” The Lancet Microbe 2 (7): e331–41. https://doi.org/10.1016/S2666-5247(21)00031-8.
3) Horesh, Gal, Grace A. Blackwell, Gerry Tonkin-Hill, Jukka Corander, Eva Heinz, and Nicholas R. Thomson. 2021. “A Comprehensive and High-Quality Collection of Escherichia Coli Genomes and Their Genes.” Microbial Genomics 7 (2). https://doi.org/10.1099/mgen.0.000499.
4) Lam, Margaret M. C., Ryan R. Wick, Louise M. Judd, Kathryn E. Holt, and Kelly L. Wyres. 2022. “Kaptive 2.0: Updated Capsule and Lipopolysaccharide Locus Typing for the Klebsiella Pneumoniae Species Complex.” Microbial Genomics 8 (3). https://doi.org/10.1099/mgen.0.000800.
5) Liu, Cindy M., Maliha Aziz, Daniel E. Park, Zhenke Wu, Marc Stegger, Mengbing Li, Yashan Wang, et al. 2023. “Using Source-Associated Mobile Genetic Elements to Identify Zoonotic Extraintestinal E. Coli Infections.” One Health (Amsterdam, Netherlands) 16 (June): 100518. https://doi.org/10.1016/j.onehlt.2023.100518.
6) Mäklin, Tommi, Harry A. Thorpe, Anna K. Pöntinen, Rebecca A. Gladstone, Yan Shao, Maiju Pesonen, Alan McNally, et al. 2022. “Strong Pathogen Competition in Neonatal Gut Colonisation.” Nature Communications 13 (1): 7417. http://dx.doi.org/10.1038/s41467-022-35178-5.
7) Pöntinen, Anna K., Rebecca A. Gladstone, Henri Pesonen, Maiju Pesonen, François Cléon, Benjamin J. Parcell, Teemu Kallonen, et al. 2023. “Β-Lactam Antibiotic Use Modulates Multi-Drug Resistant Clone Success in Escherichia Coli Populations: A Longitudinal Multi-Country Genomic Cohort Study.” https://doi.org/10.2139/ssrn.4364886.
8) Wyres, Kelly L., Ryan R. Wick, Claire Gorrie, Adam Jenney, Rainer Follador, Nicholas R. Thomson, and Kathryn E. Holt. 2016. “Identification of Klebsiella Capsule Synthesis Loci from Whole Genome Data.” Microbial Genomics 2 (12): e000102. https://doi.org/10.1099/mgen.0.000102.
