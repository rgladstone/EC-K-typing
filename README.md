# _E. coli_ Group 2 and Group 3 capsular K-typing database

This _E. coli_ Group 2 and Group 3 (ABC-transporter dependent) K-typing database is formatted for use with the tool Kaptive<sup>1</sup> version 3 https://kaptive.readthedocs.io/en/latest/

This database includes 91 loci with unique capsular gene presence-absence patterns (excluding IS elements). It is expected to cover the majority of invasive _E. coli_ isolates; group 2 and 3 capsule prevalence is near complete in phylogroups B2 and D, whilst phylogroups A, B1, and C have a lower prevalence. ~50,000 _E. coli_ genomes were screened from bloodstream infections (human), carriage (human and animal), from Europe, North America, Africa and Asia, and all G2 kpsF-positive assemblies (90% ID) in a collection of 661k bacterial assemblies <sup>3–12</sup>

Phenotypes, if known, are given in the database Genbank K-type fields. When the phenotype is unknown and the locus shares multiple region 2 genes with a known phenotype, K-like is denoted in the K-type field. According to the classification reported on https://www.iith.ac.in/EK3D/classification.php<sup>2</sup>, 21/26 group 2 phenotypes and 7/7 group 3 phenotypes are represented. We reclassified K19 as group 3.

Some capsules have the same gene presence-absence pattern and differ only in sequence, in this cases multiple K-types are reported for a given K-locus; the genetic determinants need validation before these are included as phenotypic logic for Kaptive to distinguish.
###	K96 and K54
These capsular types share the same genes and are highly conserved.
Putative determinant - one non-synonymous change results in missense for 9 AA and a truncation of 8 AA in K96-07, an alginate o-acetyltransferase.
###	K13 and K23
These capsular types share the same genes and are highly conserved. They have previously been reported to belong to a serogroup along with K20<sup>13</sup> (which has a slightly different gene set).
Putative determinant - 2 non-synonymous changes in the VatD acetyltransferase.
###	K1 and K92
Despite having the same gene content, these two capsular types have divergent sequences, allowing them to be typed separately. Therefore, both reference loci are included in the database under KL8 allele 1 (K1) and KL8 allele 2 (K92). A divergent _neuS_ gene is known to be the genetic determinant<sup>14</sup>.
###	Atypical K-loci
These K-loci have an atypical locus architecture, with the capsular-specific region 2 genes outside of regions 1 and 3. If similar loci exist with additional novel region 2 genes after the end of the loci included in the DB, they will look typeable, but the extra region 2 genes may not be captured. We recommend additional inspection of atypical loci marked as atypical in the DB K-type field that will appear in the Kaptive output. 

Currently, non-capsular IS-element-associated annotations, but not sequence, have been removed; 33 of the G2 and G3 K-loci have an IS, 17 of those have TIRs in capsule gene CDS, and 9 have putative capsular genes as cargo.

The K-loci genes were annotated using prokka 1.14.5 and panaroo 3.1.4 from https://github.com/tseemann/prokka<sup>15</sup> and https://github.com/gtonkinhill/panaroo<sup>16</sup> to give consistent annotation across the DB.

#### Cite
Gladstone, R. A. et al. Groups 2 and 3 ABC-transporter dependent capsular K-loci contribute significantly to variation in the invasive potential of _Escherichia coli_. medRxiv 2024.11.22.24317484 (2024) doi:10.1101/2024.11.22.24317484.
https://medrxiv.org/cgi/content/short/2024.11.22.24317484v1
Stanton TD, et al. Fast and Accurate in silico Antigen Typing with Kaptive 3. bioRxiv; 2025. p. 2025.02.05.636613. https://www.biorxiv.org/content/10.1101/2025.02.05.636613v1 

#### References
1.	Stanton TD, et al. Fast and Accurate in silico Antigen Typing with Kaptive 3. bioRxiv; 2025. p. 2025.02.05.636613. https://www.biorxiv.org/content/10.1101/2025.02.05.636613v1
2.	Kunduru, B. R., Nair, S. A. & Rathinavelan, T. EK3D: an E. coli K antigen 3-dimensional structure database. Nucleic Acids Res. 44, D675–81 (2016).
3.	Gladstone, R. A. et al. Emergence and dissemination of antimicrobial resistance in Escherichia coli causing bloodstream infections in Norway in 2002–17: a nationwide, longitudinal, microbial population genomic study. The Lancet Microbe 2, e331–e341 (2021).
4.	Arredondo-Alonso, S. et al. Plasmid-driven strategies for clone success in Escherichia coli. bioRxiv 2023.10.14.562336 (2023) doi:10.1101/2023.10.14.562336.
5.	Kallonen, T. et al. Systematic longitudinal survey of invasive Escherichia coli in England demonstrates a stable population structure only transiently disturbed by the emergence of ST131. Genome Res. (2017) doi:10.1101/gr.216606.116.
6.	Pöntinen, A. K. et al. Β-Lactam Antibiotic Use Modulates Multi-Drug Resistant Clone Success in Escherichia Coli Populations: A Longitudinal Multi-Country Genomic Cohort Study. (2023) doi:10.2139/ssrn.4364886.
7.	Mäklin, T. et al. Strong pathogen competition in neonatal gut colonisation. Nat. Commun. 13, 7417 (2022).
8.	Shao, Y. et al. Primary succession of Bifidobacteria drives pathogen resistance in neonatal microbiota assembly. Nat. Microbiol. 9, 2570–2582 (2024).
9.	Liu, C. M. et al. Using source-associated mobile genetic elements to identify zoonotic extraintestinal E. coli infections. One Health 16, 100518 (2023).
10.	Ludden, C. et al. One Health Genomic Surveillance of Escherichia coli Demonstrates Distinct Lineages and Mobile Genetic Elements in Isolates from Humans versus Livestock. MBio 10, (2019).
11.	Sands, K. et al. Characterization of antimicrobial-resistant Gram-negative bacteria that cause neonatal sepsis in seven low- and middle-income countries. Nat Microbiol 6, 512–523 (2021).
12.	Dicks, J. et al. NCTC3000: a century of bacterial strain collecting leads to a rich genomic data resource. Microb. Genom. 9, mgen000976 (2023).
13.	Vann, W. F. et al. Serological, chemical, and structural analyses of the Escherichia coli cross-reactive capsular polysaccharides K13, K20, and K23. Infect. Immun. 39, 623–629 (1983).
14.	Roberts, I. S. The Expression of Polysaccharide Capsules in Escherichia coli. in Glycomicrobiology 441–464 (Kluwer Academic Publishers, Boston, 2005).
15.	Seemann, T. Prokka: rapid prokaryotic genome annotation. Bioinformatics (2014).
16.	Tonkin-Hill, G. et al. Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biol. 21, 180 (2020).
