# EC-K-Typing private database builder

This container adds one private group 2 or group 3 K-locus to a frozen
EC-K-Typing reference set and produces a new, local Kaptive database. It does
not alter the Git repository, upload the private sequence, or add the private
locus to the public database.

## Requirements

Use Docker or Podman on a laptop or workstation, or Apptainer/Singularity on
an HPC system. The primary input is a single-locus Bakta-style GFF3 containing
one embedded FASTA sequence. Annotating a whole genome is deliberately outside
this container.

The private locus must use a local-only identifier of `KL9000` or greater.
These numbers are not public EC-K-Typing assignments and must not be reported
as official KL designations.

The builder checks the expected locus orientation: kpsF-to-kpsM for group 2,
and kpsD/kpsM-to-kpsS for group 3. A deliberately reviewed atypical locus can
be accepted with `--allow-atypical-boundaries`; this decision is recorded in
the output manifest.

## Build the image

From the root of the EC-K-typing repository:

```bash
docker build \
  --platform linux/amd64 \
  -f private_builder/Containerfile \
  -t ec-k-typing-private-builder:4.0.0 \
  .
```

Podman accepts the same command with `podman build`. The explicit platform
also works on Apple Silicon through container emulation and matches typical
x86-64 HPC systems.

The image and wrapper support rootless Podman installations that have only a
single UID mapping. Output created by rootless Podman remains owned by the
user running Podman; no `/etc/subuid` or `/etc/subgid` change is required for
this builder.

On an HPC system with a single UID mapping, build with
`--storage-opt ignore_chown_errors=true`. Podman's temporary build directory
must be on a filesystem with sufficient space and support for SELinux extended
attributes. The Podman graph-root filesystem is normally a suitable choice;
an NFS data filesystem may not be.

For an Apptainer-only system, build the Docker image elsewhere and either
publish it to a private registry or convert it to a SIF:

```bash
apptainer build ec-k-typing-private-builder_4.0.0.sif \
  docker-daemon://ec-k-typing-private-builder:4.0.0
```

## Build a private database

The wrapper detects an installed container runtime:

```bash
private_builder/ec-k-private-db \
  --input-gff /absolute/path/private_locus.gff3 \
  --output-dir /absolute/path/KL9001_private_DB \
  --group 2 \
  --locus KL9001 \
  --accession PRIVATE
```

An optional phenotype and one note can be included:

```bash
private_builder/ec-k-private-db \
  --input-gff /absolute/path/private_locus.gff3 \
  --output-dir /absolute/path/KL9001_private_DB \
  --group 2 \
  --locus KL9001 \
  --accession PRIVATE \
  --phenotype "Private phenotype" \
  --note "Local private locus; not an official KL assignment"
```

To select Apptainer and a SIF explicitly:

```bash
private_builder/ec-k-private-db \
  --runtime apptainer \
  --image /absolute/path/ec-k-typing-private-builder_4.0.0.sif \
  --input-gff /absolute/path/private_locus.gff3 \
  --output-dir /absolute/path/KL9001_private_DB \
  --group 2 \
  --locus KL9001 \
  --accession PRIVATE
```

## Output

The output directory contains:

- `EC-K-typing_private_KL9001.gbk`: the new local Kaptive database;
- `validation.tsv`: Kaptive self-typing of every reference record;
- `G2_KL9001_PRIVATE.gff`: the final curated private-locus GFF;
- `reconciliation.tsv`: the Panaroo metadata decision for every cluster; and
- `manifest.json`: versions, checksums, counts, and validation summary.

The build succeeds only when every reference self-matches at 100% identity and
100% coverage and none is untypeable. Kaptive problem flags are retained in
`validation.tsv`; as in the public database, a flag can reflect a valid
off-target homolog rather than a failed self-match.

On failure, `FAILED.txt` and the `.work` directory are retained for diagnosis.
On success, intermediates are deleted unless `--keep-work` is supplied.

## Reproducibility and privacy

The image contains three isolated environments:

- Panaroo 1.5.2 with NumPy 1.26.4 for clustering;
- Panaroo 1.6.0 with NumPy 1.26.4 for GFF generation; and
- Kaptive 3.1.0 plus the conversion and validation scripts.

Micromamba is used only while constructing the image. Runtime commands invoke
the pinned environment executables directly, avoiding dependence on shell
activation or micromamba runtime wrappers.

The official reference files are baked into the image and made read-only.
Docker and Podman runs use `--network none`; the input GFF is mounted
read-only, and only the selected output parent is writable. The manifest
records hashes of both the private input and the frozen official reference.
The repository `.dockerignore` allow-lists only the files required by the
builder, so unrelated working files and private data are not added to the
container build context.

Panaroo-generated annotations are retained only for the private record. The
109 frozen official GFFs are copied byte-for-byte into the derivative
database after clustering, so a Panaroo merge or split cannot rewrite a
public reference. Ambiguous private annotations receive deterministic
`private_composite_` or `private_split_` names and are listed explicitly in
`reconciliation.tsv`.

The result is a private derivative database. Running the builder does not
create a commit, branch, pull request, or public KL assignment.
