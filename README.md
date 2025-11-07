# _E. coli_ Group 2 and Group 3 capsular K-typing database

###	Database update summary V3
KL nomenclature now reflects K-phenotype numbering with unknown phenotypes as KL110+. 
For database v1 and v2 downloaded prior to 6/11/25, we recommend rerunning kaptive with the lastest version of the DB or using the lookup table provided in EC-K-typing/supplementary.

###	Information
This database allows the _in-silico_ K-typing of _E. coli_ Group 2 and Group 3 (ABC-transporter dependent) K-loci with the tool Kaptive<sup>[1]</sup>.

This database includes K-loci with unique capsular gene presence-absence patterns (excluding IS elements). It is expected to cover the majority of invasive _E. coli_ isolates; group 2 and 3 capsule prevalence is near complete in phylogroups B2 and D, whilst phylogroups A, B1, and C have a lower prevalence. ~50,000 _E. coli_ genomes were screened from bloodstream infections (human), carriage (human and animal), from Europe, North America, Africa and Asia, and all G2 kpsF-positive and G3 kpsM-positive assemblies (90% kmerID) in a collection of 661k bacterial assemblies <sup>[2–12]</sup>.

Phenotypes, if known, are provided in the GenBank K-type fields. Most (20/26) group 2 phenotypes and all (6/6) group 3 phenotypes are represented by the database <sup>[13]</sup>.

Some capsules have the same gene presence-absence pattern and differ only in sequence; in these cases, multiple K-types are reported for a given K-locus. The genetic determinants need validation before they are included as phenotypic logic for Kaptive to distinguish.

####	K96 and K54
These capsular types share the same genes and are highly conserved. K54 has been reported to be a modified K96 polysaccharide where transfer of threonine and serine replaces some of the glucuronic acid. Such modifications are thought to be encoded outside of the kps <sup>[14, 15]</sup>
####	K13 and K23
These capsular types share the same genes and are highly conserved. They have previously been reported to belong to a serogroup along with K20<sup>[16]</sup> that has replaced one gene relative to K13 and K23.
####	K1 and K92
Despite having the same gene content, these two capsular types have divergent sequences, allowing them to be typed separately. Therefore, both reference loci are included in the database under KL1 (K1) and KL92 (K92). A divergent _neuS_ gene is known to be the genetic determinant<sup>[17]</sup>.
####	Atypical K-loci
These K-loci have an atypical locus architecture, with the capsular-specific region 2 genes outside of regions 1 and 3. If similar loci exist with additional novel region 2 genes after the end of the loci included in the DB, they will look typeable, but the extra region 2 genes may not be captured. We recommend additional inspection of atypical loci marked as atypical in the DB K-type field that will appear in the Kaptive output. Please refer to the Kaptive documentation for general advice on interpreting Kaptive outputs. Kaptive<sup>1</sup> version 3 https://kaptive.readthedocs.io/en/latest/

The K-loci genes were annotated using bakta 1.10.4 and panaroo 1.5.2 from https://github.com/oschwengers/bakta<sup>[18]</sup> and https://github.com/gtonkinhill/panaroo<sup>[19]</sup> to give consistent annotation across the DB. Non-capsular IS-element-associated annotations, but not sequence, have been removed.

#### Adding novel K-loci
Please contact rebeccgl@uio.no or log an issue on this GitHub repository to have novel alleles added to the database. Alternatively, use the provided EC-K-typing/DB/panaroo_refset/gffs and panaroo-generate-gff to annotate your locus with the database gene cluster names and use the K-gff_to_gbk.py script to generate an in-house Kaptive database from the gffs.

###	How to use the database
K-typing database is formatted for use with the tool Kaptive<sup>[1]</sup> version 3.0.0b5

Quick start: kaptive assembly EC-K-typing_group2and3_v3.0.0.gbk your_assembly.fasta

For more information read https://kaptive.readthedocs.io/en/latest/

### Cite
Gladstone, R. A. et al. Groups 2 and 3 ABC-transporter-dependent capsular K-loci contribute significantly to variation in the invasive potential of _Escherichia coli_. 2025 medRxiv https://doi.org/10.1101/2024.11.22.24317484

Stanton TD, Hetland MAK, Löhr IH, Holt KE, Wyres KL. Fast and accurate in silico antigen typing with Kaptive 3. Microb Genom 2025;11:001428.

#### References
1.	Stanton TD, Hetland MAK, Löhr IH, Holt KE, Wyres KL. Fast and accurate in silico antigen typing with Kaptive 3. Microb Genom 2025;11:001428.
2.	Gladstone RA, McNally A, Pöntinen AK, Tonkin-Hill G, Lees JA, et al. Emergence and dissemination of antimicrobial resistance in Escherichia coli causing bloodstream infections in Norway in 2002–17: a nationwide, longitudinal, microbial population genomic study. The Lancet Microbe 2021;2:e331–e341.
3.	Arredondo-Alonso S, Pöntinen AK, Gama JA, Gladstone RA, Harms K, et al. Plasmid-driven strategies for clone success in Escherichia coli. Nat Commun 2025;16:2921.
4.	Kallonen T, Brodrick HJ, Harris SR, Corander J, Brown NM, et al. Systematic longitudinal survey of invasive Escherichia coli in England demonstrates a stable population structure only transiently disturbed by the emergence of ST131. Genome Res. Epub ahead of print 18 July 2017. DOI: 10.1101/gr.216606.116.
5.	Pöntinen AK, Gladstone RA, Pesonen H, Pesonen M, Cléon F, et al. Modulation of multidrug-resistant clone success in Escherichia coli populations: a longitudinal, multi-country, genomic and antibiotic usage cohort study. Lancet Microbe 2024;5:e142–e150.
6.	Shao Y, Garcia-Mauriño C, Clare S, Dawson NJR, Mu A, et al. Primary succession of Bifidobacteria drives pathogen resistance in neonatal microbiota assembly. Nat Microbiol 2024;9:2570–2582.
7.	Mäklin T, Thorpe HA, Pöntinen AK, Gladstone RA, Shao Y, et al. Strong pathogen competition in neonatal gut colonisation. Nat Commun 2022;13:7417.
8.	Liu CM, Aziz M, Park DE, Wu Z, Stegger M, et al. Using source-associated mobile genetic elements to identify zoonotic extraintestinal E. coli infections. One Health 2023;16:100518.
9.	Ludden C, Raven KE, Jamrozy D, Gouliouris T, Blane B, et al. One Health Genomic Surveillance of Escherichia coli Demonstrates Distinct Lineages and Mobile Genetic Elements in Isolates from Humans versus Livestock. MBio;10. Epub ahead of print 22 January 2019. DOI: 10.1128/mBio.02693-18.
10.	Sands K, Carvalho MJ, Portal E, Thomson K, Dyer C, et al. Characterization of antimicrobial-resistant Gram-negative bacteria that cause neonatal sepsis in seven low- and middle-income countries. Nat Microbiol 2021;6:512–523.
11.	Dicks J, Fazal M-A, Oliver K, Grayson NE, Turnbull JD, et al. NCTC3000: a century of bacterial strain collecting leads to a rich genomic data resource. Microb Genom 2023;9:mgen000976.
12.	Blackwell GA, Hunt M, Malone KM, Lima L, Horesh G, et al. Exploring bacterial diversity via a curated and searchable snapshot of archived DNA sequences. PLoS Biol 2021;19:e3001421.
13.	Kunduru BR, Nair SA, Rathinavelan T. EK3D: an E. coli K antigen 3-dimensional structure database. Nucleic Acids Res 2016;44:D675–81.
14.	Jann B, Kochanowski H, Jann K. Structure of the capsular K96 polysaccharide (K96 antigen) from Escherichia coli O77:K96:H- and comparison with the capsular K54 polysaccharide (K54 antigen) from Escherichia coli O6:K54:H10. Carbohydr Res 1994;253:323–327.
15.	Whitfield C. Biosynthesis and assembly of capsular polysaccharides in Escherichia coli. Annu Rev Biochem 2006;75:39–68.
16.	Vann WF, Soderstrom T, Egan W, Tsui FP, Schneerson R, et al. Serological, chemical, and structural analyses of the Escherichia coli cross-reactive capsular polysaccharides K13, K20, and K23. Infect Immun 1983;39:623–629.
17.	Roberts IS. The Expression of Polysaccharide Capsules in Escherichia coli. In: Glycomicrobiology. Boston: Kluwer Academic Publishers; 2005. pp. 441–464.
18.	Schwengers O, Jelonek L, Dieckmann MA, Beyvers S, Blom J, et al. Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification: Find out more about Bakta, the motivation, challenges and applications, here. Microb Genom 2021;7:000685.
19.	Tonkin-Hill G, MacAlasdair N, Ruis C, Weimann A, Horesh G, et al. Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biol 2020;21:180.
