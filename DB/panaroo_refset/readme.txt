###Adding novel loci to an in house database

### 1. Extract novel K-locus and annotate
Provided the K-locus assembles in one piece, and with the exception of atypical K-loci, the locus can be extracted with kaptive assembly.

kaptive assembly --fasta EC-K-typing/DB/EC-K-typing_group2and3_v3.0.0.gbk novel.fa

It is important to confirm that the K-locus sequence is free from Ns.

Rename your file and fasta header to look like GX_KLXXX_Accession. Where the correct group (G2 or G3) is first and the KL number is greater than KL175 and the accession is from the genome you extracted the K-locus from. Then you can annotate the K-locus sequence and inspect your K-locus boundaries. Do they make sense, e.g. are the first and last genes kpsF and kpsM for G2? Flip the sequence so the orientation is kpsF-kpsM for G2 and kpsD/M-kpsS for G3. 

Add your annotated gff to your local EC-K-typing/DB/panaroo_refset/gff folder so that it is included in the clustering, tested with (Bakta v1.10.4).

Fix for Panaroo 1.5.2 with Batka input: Add empty panaroo_id=;prepanaroo_inference=Unknown_inference; fields to the novel gff using edit_novel_gff.sh

### Panaroo 1.5.2 command 
panaroo --clean-mode sensitive --search_radius 30000 --trailing_recursive 0 --min_trailing_support 1 -t 36 -i EC-K-typing/DB/panaroo_refset/gffs/*gff -o panaroo_add_novel

### Panaroo 1.5.2 commandPanaroo generate gffs command
panaroo-generate-gffs -i EC-K-typing/DB/panaroo_refset/gffs/*.gff -o panaroo_add_novel

#NB
Some gene clusters may have been merged.You can check this check this in gene_presence_absence.csv
awk -F "," '{print $1}' gene_presence_absence.csv | grep "~"

### Convert gffs to a single kaptive gbk database
EC-K-typing/DB/panaroo_refset/K-gff_to_gbk.py
