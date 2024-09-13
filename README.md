# Main Group 2 and Group 3 database

EC-K-typing_group2and3_v1_170724.gbk is the main database, it is an _E. coli_ Group 2 and Group 3 (ABC-transporter dependant) K-typing database formatted for use with the tool Kaptive version 3 https://kaptive.readthedocs.io/en/latest/

Phenotypes, if known, are given in K-type and serotype fields. 18/26 group 2 phenotypes are represented, and 6/7 group 3 phenotypes are represented in a database of 102 loci. These loci represent unique gene presence-absence patterns excluding IS-elements.

K-antigen loci (with and without phenotype). This database is expected to cover the majority of invasive isolates; group 2 and 3 capsule prevalence is near complete in phylogroups B2 and D, whilst phylogroups A, B1, and C have a lower prevalence. ~50,000 E. coli genomes were screened from Bloodstream Infections (human), carriage (human and animal), from Europe, Africa and Asia, and all kpsF-positive assemblies in Blackwell _et al_ PLOS 2018.

Some capsules have the same gene presence absence pattern and different only in sequence
###	K96 and K54
1 non-synomous change in kpsC
1 simple non-synomous change in K96-07 and one 1 non-synomous change that results in missense for 9 AA and a truncation of 8 AA in K96-07 an alginate o-acetyltransferase.
###	K13 and K20 
Reported to belong to a serogroup https://journals.asm.org/doi/10.1128/iai.39.2.623-629.1983
2 non-synomous changes in vatD acetyltransferase
2 insertions of 10 and 12bp in IS-element and 1 non-synomous change.
	K1 and K92
These two loci are divergent, there is a particularly high density of changes in neuS resulting in 60 NS AA changes.

## Skeleton database
EC-K-typing_group1and4_v1_180724.gbk is a skeleton Group 1 and Group 4 K-typing database formatted for use with the tool Kaptive version 3 https://kaptive.readthedocs.io/en/latest/
Phenotypes, if known, are given in K-type and serotype fields. Only 3 reference loci with phenotypes are represented, in group 1 and group 4 (K38, K84, K102).

Annotated using prokka 1.14.5 and panaroo 3.1.4 from https://github.com/tseemann/prokka and https://github.com/gtonkinhill/panaroo to give consistent annotation across the DB.

Currently, IS-element-associated annotations, but not sequence, have been removed, 33 of the K-loci have an IS and of those 17 have TIRs in capsule gene CDS, 8/17 carry a capsular gene in the IS. Kaptive V3 is, therefore, required for this version, and length discrepancies will result when the IS-element is not present in the query assembly.
Another known non-capsule gene (a beta-lactamase) annotation within KL143 has also been removed, but not the sequence.

# References
1) Blackwell, Grace A., Martin Hunt, Kerri M. Malone, Leandro Lima, Gal Horesh, Blaise T. F. Alako, Nicholas R. Thomson, and Zamin Iqbal. 2021. “Exploring Bacterial Diversity via a Curated and Searchable Snapshot of Archived DNA Sequences.” PLoS Biology 19 (11): e3001421. https://doi.org/10.1371/journal.pbio.3001421.
2) Gladstone, Rebecca A., Alan McNally, Anna K. Pöntinen, Gerry Tonkin-Hill, John A. Lees, Kusti Skytén, François Cléon, et al. 2021. “Emergence and Dissemination of Antimicrobial Resistance in Escherichia Coli Causing Bloodstream Infections in Norway in 2002–17: A Nationwide, Longitudinal, Microbial Population Genomic Study.” The Lancet Microbe 2 (7): e331–41. https://doi.org/10.1016/S2666-5247(21)00031-8.
3) Horesh, Gal, Grace A. Blackwell, Gerry Tonkin-Hill, Jukka Corander, Eva Heinz, and Nicholas R. Thomson. 2021. “A Comprehensive and High-Quality Collection of Escherichia Coli Genomes and Their Genes.” Microbial Genomics 7 (2). https://doi.org/10.1099/mgen.0.000499.
4) Lam, Margaret M. C., Ryan R. Wick, Louise M. Judd, Kathryn E. Holt, and Kelly L. Wyres. 2022. “Kaptive 2.0: Updated Capsule and Lipopolysaccharide Locus Typing for the Klebsiella Pneumoniae Species Complex.” Microbial Genomics 8 (3). https://doi.org/10.1099/mgen.0.000800.
5) Liu, Cindy M., Maliha Aziz, Daniel E. Park, Zhenke Wu, Marc Stegger, Mengbing Li, Yashan Wang, et al. 2023. “Using Source-Associated Mobile Genetic Elements to Identify Zoonotic Extraintestinal E. Coli Infections.” One Health (Amsterdam, Netherlands) 16 (June): 100518. https://doi.org/10.1016/j.onehlt.2023.100518.
6) Mäklin, Tommi, Harry A. Thorpe, Anna K. Pöntinen, Rebecca A. Gladstone, Yan Shao, Maiju Pesonen, Alan McNally, et al. 2022. “Strong Pathogen Competition in Neonatal Gut Colonisation.” Nature Communications 13 (1): 7417. http://dx.doi.org/10.1038/s41467-022-35178-5.
7) Pöntinen, Anna K., Rebecca A. Gladstone, Henri Pesonen, Maiju Pesonen, François Cléon, Benjamin J. Parcell, Teemu Kallonen, et al. 2023. “Β-Lactam Antibiotic Use Modulates Multi-Drug Resistant Clone Success in Escherichia Coli Populations: A Longitudinal Multi-Country Genomic Cohort Study.” https://doi.org/10.2139/ssrn.4364886.
8) Wyres, Kelly L., Ryan R. Wick, Claire Gorrie, Adam Jenney, Rainer Follador, Nicholas R. Thomson, and Kathryn E. Holt. 2016. “Identification of Klebsiella Capsule Synthesis Loci from Whole Genome Data.” Microbial Genomics 2 (12): e000102. https://doi.org/10.1099/mgen.0.000102.
