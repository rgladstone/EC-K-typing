import sys
import datetime
import os
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from io import StringIO

# --- PYTHON VERSION COMPATIBILITY NOTE ---
# This script is designed for use with Python 3.6+ and has been tested 
# successfully with Python 3.12.0.
# ----------------------------------------

def sanitize_string(s):
    """
    Aggressively cleans a string by replacing all non-standard whitespace 
    and removing non-printable ASCII characters to prevent parser errors.
    """
    if not isinstance(s, str):
        return s
    
    # 1. Replace all non-standard whitespace with a single standard space
    s = re.sub(r'\s+', ' ', s).strip()
    
    # 2. Remove any remaining control or non-printable ASCII characters
    s = ''.join(char for char in s if 32 <= ord(char) <= 126 or ord(char) == 10)
    
    return s

def extract_accession_and_kl_type(filename):
    """
    Extracts the Accession ID and KL type from the filename. 
    Applies sanitization to ensure clean output.
    """
    full_base_name = os.path.splitext(filename)[0]
    kl_type = "KL_DEFAULT"

    default_accession = full_base_name.replace('_panarooupdated', '').replace('_noIS', '')
    accession = default_accession

    kl_match = re.search(r'(KL[A-Z0-9]+)', full_base_name, re.IGNORECASE)
    panaroo_match = re.search(r'_panarooupdated', full_base_name)

    if kl_match:
        kl_type = sanitize_string(kl_match.group(1).upper())
    
    if kl_match and panaroo_match:
        kl_end_index = kl_match.end()
        panaroo_start_index = panaroo_match.start()

        if kl_end_index < panaroo_start_index:
            accession_raw = full_base_name[kl_end_index:panaroo_start_index]
            extracted_accession = accession_raw.strip('_')

            if extracted_accession:
                accession = sanitize_string(extracted_accession)

    return accession, kl_type

def parse_k_type(filename):
    """
    Extracts the K type (e.g., K2Ab or K96-K54) from the filename.
    """
    match = re.search(r'(K[A-Z0-9-]+)_G', filename, re.IGNORECASE)
    if match:
        return sanitize_string(match.group(1).upper())
    return None

def apply_post_processing_fixes(genbank_content, locus_name, new_version, new_date, comment_text):
    """
    Applies all header fixes (LOCUS, VERSION, COMMENT) and critical feature alignment fixes.
    This function processes the raw GenBank string for a single record.
    """
    
    # CRITICAL FIX 1: Explicitly replace non-breaking spaces (\xa0)
    genbank_content = genbank_content.replace('\xa0', ' ')
    
    # CRITICAL FIX 4: Replace ANY non-ASCII character with a standard space
    genbank_content = re.sub(r'[^\x00-\x7F\s]+', ' ', genbank_content)

    # --- CRITICAL FIX 3 (ULTIMATE): FORCE PURE ASCII INDENTATION OVERWRITE ---
    
    # Feature line fix: Overwrites the start of 'source' and 'CDS' lines with 5 clean spaces.
    genbank_content = re.sub(
        r'^\s*(source|CDS)(.*)', 
        r'     \1\2', 
        genbank_content, 
        flags=re.MULTILINE
    )

    # Qualifier line fix: Overwrites the start of qualifier lines with 21 clean spaces.
    genbank_content = re.sub(
        r'^\s*(/[a-zA-Z0-9_]+="?.*)', 
        r'                     \1', 
        genbank_content, 
        flags=re.MULTILINE
    )
    # --------------------------------------------------------------------------

    lines = genbank_content.splitlines()
    output_lines = []

    # Format the comment text (using sanitized lines)
    genbank_comment_lines = []
    comment_prefix = "COMMENT     "
    continuation_prefix = " " * 12
    for i, line in enumerate(comment_text.split('\n')):
        clean_line = sanitize_string(line)
        if i == 0:
            genbank_comment_lines.append(f"{comment_prefix}{clean_line}")
        else:
            genbank_comment_lines.append(f"{continuation_prefix}{clean_line}")

    version_inserted = False
    comment_inserted = False

    for line in lines:

        if line.startswith('ACCESSION'):
            output_lines.append(line)
            # FORCE correct VERSION insertion
            output_lines.append(f"VERSION     {new_version}")
            version_inserted = True

        elif line.startswith('VERSION'):
            if version_inserted:
                continue
            else:
                output_lines.append(line)

        elif line.startswith('LOCUS'):
            match = re.search(r'LOCUS\s+.*?\s+(\d+)\s+bp', line) 

            if match:
                length = match.group(1)
                required_padding = 45 - (9 + len(locus_name))
                if required_padding < 1:
                    required_padding = 1

                new_locus_line = (
                    f"LOCUS       {locus_name}"
                    f"{' ' * required_padding}{length:>6} bp    DNA     linear   UNK {new_date}"
                )
                output_lines.append(new_locus_line)
            else:
                output_lines.append(line)

        elif line.startswith('KEYWORDS') or line.startswith('SOURCE') or line.startswith('FEATURES'):
            if not comment_inserted:
                 output_lines.extend(genbank_comment_lines)
                 comment_inserted = True
            output_lines.append(line)
        
        elif line.startswith('//'):
            # This is the separator line, let Biopython handle it, but trim surrounding whitespace
            output_lines.append(line.strip())

        else:
            output_lines.append(line)

    return "\n".join(output_lines)

def process_gff_to_cleaned_string(gff_filename, base_version='v2'):
    """
    Converts a single GFF file to a fully-cleaned GenBank record string.
    """
    # --- Custom Derivation & Formatting ---
    ACCESSION_ID, KL_TYPE = extract_accession_and_kl_type(gff_filename)
    
    FULL_LOCUS_NAME = os.path.splitext(gff_filename)[0]
    FULL_LOCUS_NAME = FULL_LOCUS_NAME.replace('_panarooupdated', '').replace('_noIS', '')
    FULL_LOCUS_NAME = sanitize_string(FULL_LOCUS_NAME)

    LOCUS_DISPLAY_NAME = sanitize_string(f"{KL_TYPE}_{ACCESSION_ID}")

    K_TYPE_FULL = parse_k_type(gff_filename) 
    K_TYPE_SEROTYPE = K_TYPE_FULL.split('-')[0] if K_TYPE_FULL and '-' in K_TYPE_FULL else K_TYPE_FULL
    
    GB_VERSION = sanitize_string(f"{KL_TYPE}_{base_version}")
    
    CURRENT_DATE = datetime.date.today().strftime("%d-%b-%Y").upper()

    CUSTOM_COMMENT = (
        "Annotated using Bakta v1.10.4 and panaroo v1.5.2 from\n"
        "https://github.com/oschwengers/bakta and\n"
        "https://github.com/gtonkinhill/panaroo"
    )

    base_def = "Escherichia coli DNA, capsular polysaccharide synthesis gene cluster"
    full_def = f"{base_def}, serotype: {K_TYPE_SEROTYPE}".rstrip('. ') if K_TYPE_SEROTYPE else base_def.rstrip('. ')
    full_def = sanitize_string(full_def)

    TAX_LINEAGE = (
        "Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales;"
        " Enterobacteriaceae; Escherichia."
    )

    # --- GFF Parsing and SeqRecord Creation ---
    features = []
    sequence = ""
    in_fasta_block = False
    cds_counter = 1

    try:
        with open(gff_filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('##FASTA') or (line.startswith('#') and '>' in line):
                    in_fasta_block = True
                    continue
                if in_fasta_block:
                    if not line.startswith('>'):
                        sequence += line.upper().replace(' ', '')
                    continue
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 9:
                    continue

                seqid, source_gff, feature_type, start, end, score, strand, phase, attributes = parts
                if feature_type != 'CDS':
                    continue

                qualifiers = {}
                attrs = dict(item.split('=', 1) for item in attributes.split(';') if '=' in item)

                clean_kl_type = sanitize_string(KL_TYPE)
                new_locus_tag = f"{clean_kl_type}_{cds_counter:05d}"
                qualifiers['locus_tag'] = new_locus_tag
                cds_counter += 1

                if 'name' in attrs:
                    qualifiers['gene'] = sanitize_string(attrs['name'])
                if 'description' in attrs:
                    qualifiers['product'] = sanitize_string(attrs['description'])

                loc_start = int(start) - 1
                loc_end = int(end)
                loc_strand = 1 if strand == '+' else -1

                feature = SeqFeature(
                    FeatureLocation(loc_start, loc_end, strand=loc_strand),
                    type="CDS",
                    qualifiers=qualifiers
                )
                features.append(feature)
    except FileNotFoundError:
        print(f"\nError: File '{gff_filename}' not found. Skipping.")
        return None
    except Exception as e:
        print(f"\nError processing file '{gff_filename}': {e}. Skipping.")
        return None

    if not sequence:
        print(f"\nError: Could not extract DNA sequence from the GFF file's FASTA block in '{gff_filename}'. Skipping.")
        return None

    seq_object = Seq(sequence)

    record = SeqRecord(
        seq_object,
        id=sanitize_string(ACCESSION_ID),
        name=sanitize_string(LOCUS_DISPLAY_NAME),
        description=full_def 
    )

    record.annotations['sequence_version'] = base_version
    record.annotations['source'] = "Escherichia coli"
    record.annotations['organism'] = "Escherichia coli"
    record.annotations['taxonomy'] = TAX_LINEAGE.split('; ')
    record.annotations['molecule_type'] = 'DNA'

    note_list = [f"K locus:{KL_TYPE}"]
    source_qualifiers = {
        "organism": sanitize_string("Escherichia coli"),
        "mol_type": sanitize_string("genomic DNA"),
        "db_xref": sanitize_string("taxon:562"),
    }

    if K_TYPE_SEROTYPE:
        source_qualifiers["serotype"] = sanitize_string(K_TYPE_SEROTYPE)
    
    if K_TYPE_FULL:
        note_list.append(f"K type:{K_TYPE_FULL}")

    note_list.append("IS-element related annotations removed")
    
    source_qualifiers["note"] = [sanitize_string(n) for n in note_list]

    source_feature = SeqFeature(
        FeatureLocation(0, len(sequence), strand=1),
        type="source",
        qualifiers=source_qualifiers
    )
    record.features.append(source_feature)
    record.features.extend(features)

    # 4. Write to buffer
    handle = StringIO()
    # Write a single record. SeqIO adds the final '//'
    SeqIO.write(record, handle, "genbank")
    raw_genbank_content = handle.getvalue()

    # 5. POST-PROCESS: Apply all the critical fixes to the single record string
    final_genbank_content = apply_post_processing_fixes(
        raw_genbank_content,
        LOCUS_DISPLAY_NAME,
        GB_VERSION,
        CURRENT_DATE,
        CUSTOM_COMMENT
    )
    
    return final_genbank_content


# --- Command Line Execution ---
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python K-gff_to_gbk_v5.py <GFF_FILENAME_1> [<GFF_FILENAME_2> ...] or python K-gff_to_gbk_v5.py *.gff")
        sys.exit(1)

    gff_files = sys.argv[1:]
    all_genbank_content = ""
    
    # Process each GFF file one by one
    for gff_file in gff_files:
        print(f"\n--- Processing file: {gff_file} ---")
        cleaned_content = process_gff_to_cleaned_string(gff_file, base_version='v2')
        
        if cleaned_content:
            if all_genbank_content:
                # CRITICAL FIX for //LOCUS: Ensure a newline separates the previous "//" from the next "LOCUS"
                # The GenBank standard requires the separator to be on its own line: //\nLOCUS
                all_genbank_content += "\n" + cleaned_content
            else:
                all_genbank_content = cleaned_content
        
    
    # 6. Write the final master GenBank database file
    if all_genbank_content:
        # Set the custom output path
        output_filename = "new_DB/EC-K-typing_group2and3_vX.X.X.gbk"
        
        # Ensure the output directory exists
        output_dir = os.path.dirname(output_filename)
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
                print(f"Created directory: {output_dir}")
            except OSError as e:
                print(f"Error creating directory {output_dir}: {e}")
                sys.exit(1)
        
        # Ensure the final output file does not have an extra newline at the end of the last record's "//"
        final_output = all_genbank_content.strip()

        try:
            # CRITICAL FIX 2: Force ASCII encoding on write to eliminate ALL non-standard characters
            with open(output_filename, "w", encoding='ascii', errors='replace') as out_handle:
                out_handle.write(final_output)
            
            print(f"\n\n✅ Final Database Successfully Generated!")
            print(f"The multi-record GenBank file '{output_filename}' has been created with corrected separators and whitespace alignment.")
        except Exception as e:
            print(f"\nError writing Final GenBank file: {e}")
