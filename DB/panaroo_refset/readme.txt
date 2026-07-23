# Adding novel loci to the EC-K-Typing group 2/3 database

This workflow adds a newly identified K locus to the Panaroo reference set,
reconciles its gene annotations with the existing loci, generates updated GFFs,
and converts the reference set into a Kaptive GenBank database.

Run the commands below from the root of the EC-K-typing repository.


## Software requirements

The reference workflow has been tested with:

- Bakta 1.10.x
- Panaroo 1.5.2 for clustering
- Panaroo 1.6.0 for panaroo-generate-gffs
- Biopython 1.84
- NumPy 1.26.4

Record the actual Bakta software and database versions used for each release.

Panaroo 1.5.2 should be run with Biopython 1.84. Biopython 1.85 and later
use a stricter FASTA parser and can cause Panaroo 1.5.2 to fail with:

    This FASTA file contains comments at the beginning of the file.

Panaroo 1.5.2 also requires NumPy 1.x. NumPy 1.26.4 is used for this
workflow. Newer NumPy releases remove the deprecated ndarray.tostring()
method used by Panaroo 1.5.2 and can cause:

    AttributeError: 'numpy.ndarray' object has no attribute 'tostring'

Confirm the active versions before running:

    bakta --version
    panaroo --version
    python -c 'import Bio, numpy; print("Biopython:", Bio.__version__); print("NumPy:", numpy.__version__)'

The repository contains separate scheduler scripts for clustering and GFF
generation because the validated environments use different Panaroo versions.

## 1. Extract and annotate the novel K locus

Provided that the K locus assembles in one piece, it can be extracted from an
assembly using Kaptive. Atypical or fragmented loci may require manual
extraction and boundary inspection.

Example using Kaptive 3:

    kaptive assembly \
        EC-K-typing_group2and3.gbk \
        novel_assembly.fasta \
        -f novel_k_locus.fasta

Confirm that the extracted locus:

- is a single sequence;
- contains no Ns;
- has biologically sensible boundaries;
- has the expected orientation where possible.

For a conventional group 2 locus, the expected orientation is generally
kpsF to kpsM. For a conventional group 3 locus, it is generally kpsD/kpsM
to kpsS. Some valid atypical loci do not have conventional boundaries, so
do not trim a locus solely to force this pattern.

Annotate the locus with Bakta and retain the complete GFF FASTA section.


## 2. Name and prepare the Bakta GFF

Use one GFF per K locus. The filename and embedded sequence identifier must
both be unique.

For an EC K locus without a known phenotype, use:

    GX_KLXXX_Accession.gff

For a locus with an established K phenotype, use:

    KXX_GX_KLXX_Accession.gff

Unknown novel EC KL numbers should be assigned after the highest identifier
in the current database; for v4.0.0 the next available identifier is KL192.
Established K phenotypes should follow the current database naming policy.

Do not use kTYPr allele or genotype names as EC-K-Typing locus names,
phenotypes, notes, or metadata.

The filename must end exactly in .gff, not .gff3. Panaroo 1.5.2
panaroo-generate-gffs removes the .gff suffix when matching input filenames
to isolate identifiers stored in gene_data.csv. A .gff3 extension can cause:

    StopIteration

Copy the annotated GFF into:

    DB/panaroo_refset/gffs/

Run the preparation script on every newly added Bakta GFF:

    DB/panaroo_refset/edit_novel_gff.sh \
        DB/panaroo_refset/gffs/GX_KLXXX_Accession.gff

The preparation script:

1. adds empty panaroo_id and prepanaroo_inference fields to every CDS;
2. converts Bakta's standard single-# metadata lines to ## comments;
3. validates that every CDS was updated;
4. is safe to run more than once.

Bakta normally adds these five single-# metadata lines:

    # Annotated with Bakta
    # Software:
    # Database:
    # DOI:
    # URL:

These must be changed to ## comments before panaroo-generate-gffs. Otherwise,
the lines may be interpreted as GFF feature rows and cause:

    IndexError: list index out of range

Do not remove or alter the ##FASTA directive.


## 3. Validate the complete Panaroo input set

The input directory must contain exactly one GFF for every reference locus.

Do not include a combined database GFF alongside the individual locus GFFs.
This creates duplicated loci and non-unique sequence identifiers.

Check the number of inputs:

    find DB/panaroo_refset/gffs \
        -maxdepth 1 \
        -type f \
        -name '*.gff' |
    wc -l

Check that all filenames are unique:

    find DB/panaroo_refset/gffs \
        -maxdepth 1 \
        -type f \
        -name '*.gff' \
        -printf '%f\n' |
    sort |
    uniq -d

The duplicate-filename check should produce no output.

Check that every GFF has a FASTA sequence:

    for gff in DB/panaroo_refset/gffs/*.gff; do
        grep -q '^##FASTA$' "$gff" ||
            echo "Missing FASTA section: $gff"
    done

Check that the first non-empty line after ##FASTA is a FASTA header:

    for gff in DB/panaroo_refset/gffs/*.gff; do
        first=$(awk '
            found && NF {print; exit}
            /^##FASTA$/ {found=1}
        ' "$gff")

        case "$first" in
            ">"*) ;;
            *) printf '%s\tFIRST_AFTER_FASTA=%s\n' "$gff" "$first" ;;
        esac
    done

Both FASTA checks should produce no output.


## 4. Run Panaroo 1.5.2

Run Panaroo through the scheduler when requesting multiple CPUs.

The reference command is:

    panaroo \
        --clean-mode sensitive \
        --search_radius 30000 \
        --trailing_recursive 0 \
        --min_trailing_support 1 \
        -t 36 \
        -i DB/panaroo_refset/gffs/*.gff \
        -o panaroo_add_novel

The output directory must not already contain a failed or incomplete run.
Move a failed directory aside before resubmitting so that its contents cannot
be mistaken for a successful result.

Confirm that the Panaroo log reports the expected executable, Panaroo 1.5.2,
Biopython 1.84, and the expected number of input GFFs.


## 5. Inspect the Panaroo clustering

Panaroo may merge gene clusters. Inspect the first column of
gene_presence_absence.csv for merged cluster names containing a tilde:

    awk -F ',' 'NR > 1 {print $1}' \
        panaroo_add_novel/gene_presence_absence.csv |
    grep '~'

No output means that no cluster name contains a tilde. Any merged cluster must
be reviewed before generating the final database.

Also inspect the new loci for:

- unexpected gene splitting or merging;
- missing conserved export genes;
- duplicated genes;
- inappropriate annotation transfer between unrelated loci;
- changes to established reference loci.

## Curate merged cluster annotations

Adding loci can cause Panaroo to merge previously distinct gene clusters.
Every cluster name containing ~~~ must be reviewed before generating GFFs.

Mappings must be recorded in:

    DB/panaroo_refset/merged_cluster_curation.csv

Where a new Bakta name joins one established reference cluster, retain the
established curated cluster name.

Where Panaroo joins two established suffixed clusters, do not replace them
with a generic name such as kpsM, kpsT or kpsC because other distinct clusters
with those base names remain in the database. Use a unique composite name,
for example:

    kpsM_5_6
    kpsT_4_5
    kpsC_6_7

Apply the reviewed mapping using:

    python DB/panaroo_refset/curate_merged_clusters.py \
        panaroo_add_novel/final_graph.gml \
        DB/panaroo_refset/merged_cluster_curation.csv \
        panaroo_add_novel/final_graph_curated.gml

The script must edit only the name, annotation and description node
attributes. Node identifiers, labels, membership, sequences and edges must
remain unchanged.

When writing Panaroo GML files with NetworkX, use the default read_gml()
label handling. Do not write a graph loaded using label=None because
NetworkX will replace the original GML label attributes with numeric node
identifiers.

Keep an unedited copy of final_graph.gml and validate the curated graph
before making it active.

## Restore all historical and new node metadata

Panaroo-generated group names are transient and must not be matched between
runs by name. The comprehensive metadata mapping in
node_metadata_curation.csv was produced by exact protein/DNA sequence
membership comparison with the historical curated graph.

It:

- restores 269 exact historical clusters;
- retains five reviewed composite clusters;
- assigns stable DB_group identifiers to 89 genuinely new clusters.

Apply the complete mapping after the merged-cluster review:

    python DB/panaroo_refset/apply_node_metadata_curation.py \
        panaroo_add_novel/final_graph_curated.gml \
        DB/panaroo_refset/node_metadata_curation.csv \
        panaroo_add_novel/final_graph_fully_curated.gml

The script requires exactly one mapping for every graph node and validates
that only name, annotation and description metadata change. A future Panaroo
run that produces a different node set must receive a newly reviewed mapping;
do not reuse a partial or name-based mapping.

After validation, keep a backup and make the fully curated graph active:

    cp panaroo_add_novel/final_graph.gml \
        panaroo_add_novel/final_graph_before_full_metadata_curation.gml

    cp panaroo_add_novel/final_graph_fully_curated.gml \
        panaroo_add_novel/final_graph.gml

## 6. Generate the updated GFFs

panaroo-generate-gffs must receive the same set of input GFFs used for the
corresponding Panaroo run.

Panaroo 1.5.2 and 1.6.0 can fail when historical GFFs contain the auxiliary
attribute panaroo_ID with an uppercase ID:

    IndexError: list index out of range

Do not alter the source GFFs to work around this parser issue. Create a
validated copy in which only this auxiliary key is changed to lowercase:

    python DB/panaroo_refset/prepare_gffs_for_generate.py \
        DB/panaroo_refset/gffs \
        panaroo_add_novel/gffs_for_generate

The script refuses to overwrite an existing output directory and confirms
that filenames, FASTA identifiers and sequences are unchanged.

Generate GFFs with Panaroo 1.6.0 using:

    panaroo-generate-gffs \
        -i panaroo_add_novel/gffs_for_generate/*.gff \
        -o panaroo_add_novel \
        -t 36

The equivalent scheduler command is:

    sbatch DB/panaroo_refset/run_panaroo_generate_gffs.sbatch

Do not mix GFFs from another Panaroo run or add or remove inputs between the
Panaroo and panaroo-generate-gffs commands.

Confirm that one updated GFF was produced for every input locus.

## Filter reviewed short annotations and stage canonical GFFs

Do not apply a general minimum gene-length filter. Existing curated loci
contain legitimate CDS features shorter than 200 bp and these historical
features are retained.

For this release:

- remove Panaroo candidate_gene refounds shorter than 200 bp;
- remove hypothetical CDS features shorter than 200 bp only when their
  DB_group identifier was newly assigned during this database expansion
  (DB_group_257 or later);
- retain all named or functionally annotated short CDS features;
- retain historical DB_group annotations, including those shorter than
  200 bp;
- retain the historical 629-bp vatD candidate annotation and three complete
  1086-bp rmlB_2 candidates in the curated GFFs.

The newly assigned short hypothetical CDS filter removes Bakta-predicted
ORFs shorter than 67 amino acids without altering any locus sequence. These
features caused full or partial off-target Kaptive gene matches and made five
perfect self-references untypeable. The established short CDS annotations
were preserved for backward compatibility.

Stage the generated GFFs with canonical reference filenames using:

    python DB/panaroo_refset/finalise_generated_gffs.py \
        panaroo_add_novel/postpanaroo_gffs \
        DB/panaroo_refset/gffs \
        panaroo_add_novel/final_curated_gffs

The script:

- removes candidate_gene features shorter than 200 bp;
- removes newly assigned hypothetical CDS features shorter than 200 bp;
- never removes an established, named or functionally annotated CDS;
- rejects IS annotations, merged names and unnamed genes;
- verifies all FASTA identifiers and sequences;
- uses the existing reference directory only as the source of canonical
  filenames;
- refuses to overwrite an existing output directory.

Review the reported retained candidates before replacing the reference GFF
directory.


## 7. Convert the curated GFF reference set to a Kaptive database

After reviewing the generated GFFs and completing any intended annotation or
IS-element curation, convert the reference set using:

    python DB/panaroo_refset/K-gff_to_gbk.py \
        DB/panaroo_refset/gffs/*.gff \
        --metadata DB/panaroo_refset/locus_metadata.csv \
        --output new_DB/EC-K-typing_group2and3_vX.X.X.gbk

The converter derives each accession from the embedded FASTA identifier, not
from the GFF filename. This preserves GCA/GCF accessions when a historical
filename contains an ERZ identifier.

Candidate genes with valid lengths are converted to GenBank CDS features
because Kaptive ignores non-CDS gene features. The three complete 1086-bp
rmlB_2 refounds are therefore included. The historical KL114 vatD candidate
contains a one-base deletion, is 629 bp long and is not promoted to a CDS.
The converter reports every invalid candidate it skips.

Known K phenotypes are derived from canonical GFF filename prefixes. Explicit
phenotype overrides and the single additional note for a locus are stored in:

    DB/panaroo_refset/locus_metadata.csv

For example, the KL158 capsule-null prediction is an unconditional locus
phenotype and is recorded there. Add a reviewed K-locus family note to the
same row when one is required.

The converter writes phenotype notes as `K type: VALUE`, including the space
after the colon. Kaptive 3.1.0's default type parser requires this space to
retain punctuation and additional words. Without it, `K18-K22` is truncated
to `K18` and `Capsule null predicted` is truncated to `Capsule`.

Before accepting the generated GenBank database, verify:

- the release version and date;
- locus names;
- accessions;
- K phenotype metadata;
- special phenotype notes, including capsule-null predictions;
- locus-family notes;
- the number of records;
- unique locus names and locus tags;
- successful loading by Kaptive.

The conversion script derives locus names, accessions, and known K phenotypes
from GFF filenames. Therefore, generated filenames must be checked before
conversion.


## 8. Final validation

Before release:

1. Type every database reference sequence against the new database.
2. Confirm that each reference returns its intended EC KL.
3. Confirm the intended phenotype for loci with established K types.
4. Compare results for all pre-existing references with the previous release.
5. Apply the documented KL rename crosswalk where identifiers have changed.
6. Test representative genomes and previously untypeable isolates.
7. Confirm that deprecated or kTYPr-only designations have not entered the
   EC-K-Typing metadata.
8. Review git diff before committing the updated database and reference GFFs.

For the 109-locus database expansion, the reference self-test produced:

- 109/109 perfect locus matches;
- zero untypeable references;
- 20 references with Kaptive problem flags.

The previous 90-locus database produced 12 problem-flagged self-matches and
zero untypeable references. Problem flags alone are therefore not a release
failure. Review every newly introduced flag, but require all self-references
to remain typeable and to match their own locus at 100% identity and
coverage.
