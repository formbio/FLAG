import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser(description="Process gene annotation and related files.")

parser.add_argument('--annotation_file', default=None, help='Path to the annotation GTF file')
parser.add_argument('--functional_simplified', default=None, help='Path to the simplified functional file')
parser.add_argument('--output_file', default=None, help='Output file')

args = parser.parse_args()

annotation_file_path = args.annotation_file
functional_simplified = args.functional_simplified
output_file = args.output_file

# Load simplified functional annotations into a dictionary
simplified_functional_df = pd.read_csv(functional_simplified, sep='\t', header=None)
simplified_functional_dict = pd.Series(simplified_functional_df.iloc[:, 1:].values.tolist(), index=simplified_functional_df[0]).to_dict()


def adjust_annotations_for_tRNA_with_line_numbers_and_unique_ids(input_file_path, output_file_path):
    annotations = []
    dummy_gene_counter = 1
    last_gene_biotype = None  # Add a variable to track the last non-exon gene_biotype

    with open(input_file_path, 'r') as infile:
        annotations = infile.readlines()

    corrected_annotations = []
    for i, line in enumerate(annotations):
        if line.startswith('#') or not line.strip():
            corrected_annotations.append(line)
            continue

        parts = line.strip().split('\t')
        feature_type = parts[2]
        annotation_tool = parts[1]
        start_site = parts[4]
        
        # If it's a gene, tRNA, or pseudogene, update the last_gene_biotype variable
        if feature_type in ['gene', 'pseudogene', 'tRNA']:
            attributes = parse_attributes(parts[8])
            last_gene_biotype = attributes.get('gene_biotype', None)

        if feature_type in ['tRNA', 'exon', 'RNA']:
            attributes = parse_attributes(parts[8])
            if feature_type == 'exon' and annotation_tool == 'tRNAscan-SE' and last_gene_biotype:
                # print(parts)
                # print(last_gene_biotype)
                attributes['gene_biotype'] = last_gene_biotype  # Set the exon's biotype based on the last gene/pseudogene/tRNA
            
                # Update only if we have previously inserted a dummy gene
                if corrected_annotations and 'dummy_gene_' in corrected_annotations[-1]:
                    attributes['gene_id'] = last_dummy_gene_id
                    attributes['transcript_id'] = f"{last_dummy_gene_id}.t1"
                    if feature_type == 'exon':
                        currentID = attributes['ID']
                        newExonNum = currentID.split('.exon', 1)[-1]
                        attributes['ID'] = f"{last_dummy_gene_id}.t1.exon{newExonNum}"
                        attributes['Parent'] = f"{last_dummy_gene_id}.t1"
                    else:
                        attributes['ID'] = f"{last_dummy_gene_id}.t1"
                        attributes['Parent'] = f"{last_dummy_gene_id}"
                    parts[8] = '; '.join([f'{k} "{v}"' for k, v in attributes.items()]) + ';'
                    line = '\t'.join(parts) + '\n'
                else:
                    parts[8] = '; '.join([f'{k} "{v}"' for k, v in attributes.items()]) + ';'
                    line = '\t'.join(parts) + '\n'
            elif feature_type == 'RNA':
                attributes['gene_biotype'] = last_gene_biotype  # Set the exon's biotype based on the last gene/pseudogene/tRNA
                parts[2] = 'tRNA'  # Ensure this line remains a tRNA entry

                # Update only if we have previously inserted a dummy gene
                if corrected_annotations and 'dummy_gene_' in corrected_annotations[-1]:
                    attributes['gene_id'] = last_dummy_gene_id
                    attributes['transcript_id'] = f"{last_dummy_gene_id}.t1"
                    if feature_type == 'exon':
                        currentID = attributes['ID']
                        newExonNum = currentID.split('.exon', 1)[-1]
                        attributes['ID'] = f"{last_dummy_gene_id}.t1.exon{newExonNum}"
                        attributes['Parent'] = f"{last_dummy_gene_id}.t1"
                    else:
                        attributes['ID'] = f"{last_dummy_gene_id}.t1"
                        attributes['Parent'] = f"{last_dummy_gene_id}"
                    parts[8] = '; '.join([f'{k} "{v}"' for k, v in attributes.items()]) + ';'
                    line = '\t'.join(parts) + '\n'
                else:
                    parts[8] = '; '.join([f'{k} "{v}"' for k, v in attributes.items()]) + ';'
                    line = '\t'.join(parts) + '\n'
            else:
                # Update only if we have previously inserted a dummy gene
                if corrected_annotations and 'dummy_gene_' in corrected_annotations[-1]:
                    attributes['gene_id'] = last_dummy_gene_id
                    attributes['transcript_id'] = f"{last_dummy_gene_id}.t1"
                    if feature_type == 'exon':
                        currentID = attributes['ID']
                        newExonNum = currentID.split('.exon', 1)[-1]
                        attributes['ID'] = f"{last_dummy_gene_id}.t1.exon{newExonNum}"
                        attributes['Parent'] = f"{last_dummy_gene_id}.t1"
                    else:
                        attributes['ID'] = f"{last_dummy_gene_id}.t1"
                        attributes['Parent'] = f"{last_dummy_gene_id}"
                    parts[8] = '; '.join([f'{k} "{v}"' for k, v in attributes.items()]) + ';'
                    line = '\t'.join(parts) + '\n'
            

        if feature_type == 'tRNA' or feature_type == 'RNA':
            # Check if the previous line is not a gene or pseudogene
            if i == 0 or not (annotations[i-1].split('\t')[2] in ['gene', 'pseudogene']):
                dummy_gene_line, dummy_gene_id = generate_dummy_gene_line_based_on_tRNA(parts, dummy_gene_counter)
                corrected_annotations.append(dummy_gene_line)
                last_dummy_gene_id = dummy_gene_id  # Keep the last inserted dummy gene ID
                dummy_gene_counter += 1
                # Correct the feature_type for the tRNA line to preserve the original 'tRNA'
                parts[2] = 'tRNA'  # Ensure this line remains a tRNA entry
                attributes = parse_attributes(parts[8])
                attributes['gene_id'] = f"dummy_gene_{dummy_gene_counter - 1}"  # Update gene_id to match the last inserted dummy gene
                parts[8] = '; '.join([f'{k} "{v}"' for k, v in attributes.items()]) + ';'
                line = '\t'.join(parts) + '\n'

        corrected_annotations.append(line)
        
        if attributes['gene_id'] == 'agat-pseudogene-1':
            print('yes')
            print(line)

    with open(output_file_path, 'w') as outfile:
        for annotation in corrected_annotations:
            outfile.write(annotation)

def generate_dummy_gene_line_based_on_tRNA(tRNA_parts, dummy_gene_counter):
    """
    Generate a dummy gene line with unique gene_id and attributes based on tRNA parts provided.
    Returns both the dummy gene line and the new dummy gene ID.
    """
    tRNA_attributes_str = tRNA_parts[8]
    tRNA_attributes = parse_attributes(tRNA_attributes_str)
    
    # Extract attributes from tRNA needed for dummy gene
    anticodon = tRNA_attributes.get('anticodon', 'unknown')
    gene_biotype = tRNA_attributes.get('gene_biotype', 'tRNA')
    isotype = tRNA_attributes.get('isotype', 'unknown')
    
    # Create a unique gene_id for the dummy gene
    dummy_gene_id = f"dummy_gene_{dummy_gene_counter}"

    dummy_gene_attributes_str = f'\tgene_id "{dummy_gene_id}"; anticodon "{anticodon}"; gene_biotype "{gene_biotype}"; isotype "{isotype}";'
    tRNA_parts[2] = 'gene'
    dummy_gene_line = '\t'.join(tRNA_parts[:8]) + dummy_gene_attributes_str + '\n'
    return dummy_gene_line, dummy_gene_id

def parse_attributes(attributes_str):
    """
    Parse a string of GTF attributes into a dictionary.
    """
    attributes = {}
    for attr in attributes_str.split(';'):
        if not attr.strip():
            continue
        key_value_pair = attr.strip().split(' ', 1)
        if len(key_value_pair) == 2:
            key, value = key_value_pair
            attributes[key] = value.strip('"')
    return attributes

# Adjust annotations and write to a new file
adjust_annotations_for_tRNA_with_line_numbers_and_unique_ids(annotation_file_path, 'CorrectedAnnotation.gtf')


def load_annotations_with_transcripts(annotation_file_path):
    # Initialize containers for genes (including pseudogenes) and their related features
    gene_related_entries = {}

    with open(annotation_file_path, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue  # Skip header lines or empty lines
            parts = line.strip().split('\t')
            feature_type = parts[2]
            # print(parts)
            attributes = parse_attributes(parts[8])


            # Check if the entry is a gene, pseudogene, mRNA, or related to either (not 'five_prime_utr', 'three_prime_utr' for now)
            if feature_type in ['gene', 'RNA', 'tRNA', 'pseudogene', 'mRNA', 'exon', 'CDS'] and ('gene_id' in attributes or 'transcript_id' in attributes):
                gene_id = attributes.get('gene_id', attributes.get('transcript_id', ''))

                if gene_id not in gene_related_entries:
                    gene_related_entries[gene_id] = []

                gene_related_entries[gene_id].append({
                    'type': feature_type,
                    'start': parts[3],
                    'end': parts[4],
                    'strand': parts[6],
                    'attributes': attributes,
                    'raw_line': line.strip()  # Keeping the raw line for further processing
                })

    return gene_related_entries

# Updated to parse both gene_id and transcript_id
# def parse_attributes(attributes_str):
#     attributes = {}
#     for attr in attributes_str.split(';'):
#         if attr.strip():  # Ensure non-empty strings
#             print(attr.strip())
#             key, value = attr.strip().split(' ', 1)
#             attributes[key.strip()] = value.strip('"')
#     return attributes
def parse_attributes(attributes_str):
    """
    Parse a string of GTF attributes into a dictionary using regex to handle
    semicolons and spaces within quoted values correctly.
    """
    attributes = {}
    # Regex pattern to match key-value pairs where value is wrapped in quotes
    # and may contain semicolons or spaces.
    pattern = re.compile(r'([^;]+) "(.*?)"(?=;|$)')
    
    for match in pattern.finditer(attributes_str):
        key, value = match.groups()
        key = key.strip()
        if key:  # Ensure key is not empty
            attributes[key] = value.strip()

    return attributes

# Use the updated function
annotations_with_transcripts = load_annotations_with_transcripts('CorrectedAnnotation.gtf')

def add_functional_annotations_optimized(annotations, functional_dict):
    for gene_id, entries in annotations.items():
        # Initialize a flag to identify if the gene instance has been updated
        updated = False
        # Store the updates here to apply to all entries once a match is found
        updates = {}
        
        # First, check if any of the entries match the keys in the functional_dict
        for entry in entries:
            key = entry['attributes'].get('ID', entry['attributes'].get('transcript_id', ''))
            if key in functional_dict and not updated:
                functional_annot = functional_dict[key]
                # Process the gene_name to include only text after the final '|'
                gene_name = functional_annot[-1]  # Assuming last element is gene_name
                # if '|' in gene_name:
                #     gene_name = gene_name.split('|')[-1]
                
                eggnog_taxscope = functional_annot[1]  # Assuming last element is gene_name
                if '|' in eggnog_taxscope:
                    eggnog_taxscope = eggnog_taxscope.split('|')[-1]

                updates = {
                    'eggnog_ortholog': functional_annot[0],
                    'gene_name': gene_name,
                    'eggnog_taxscope': eggnog_taxscope,
                    'gene_description': functional_annot[2].replace(';','')
                }

                updated = True  # Set flag to true to indicate updates are ready to apply
                
        # If updates are found, apply them to all entries of the gene instance
        if updated:
            for entry in entries:
                entry['attributes'].update(updates)
                # Reconstruct the raw line with updated attributes
                attributes_str = '; '.join([f'{k} "{v}"' for k, v in entry['attributes'].items()])
                entry['raw_line'] = '\t'.join(entry['raw_line'].split('\t')[:8] + [attributes_str + ';'])

    return annotations
    
annotations_with_functional_info = add_functional_annotations_optimized(annotations_with_transcripts, simplified_functional_dict)

def reformat_and_number_exons(annotations):
    updated_annotations = {}
    formbioid_counter = 1  # Start counter for FormBio IDs

    for gene_id, entries in annotations.items():
        formbioid = f"FormBio{formbioid_counter:06d}"  # Format the FormBio ID
        exon_number = None
        last_exon_number = None
        increment_direction = None  # Will determine if we increment or decrement
        firstExon = -1
        CDSFound = False
        for entry in entries:
            fields = entry['raw_line'].split('\t')
            annotationTool = fields[1].replace(" ", "")
            # Universal changes to source and formatting based on entry type
            fields[1] = "FormBio"  # Change source to "FormBio"
            entry['attributes']['gene_id'] = formbioid  # Update gene_id
            # Delete the 'Name' attribute
            if 'Name' in entry['attributes']:
                del entry['attributes']['Name']
            if entry['type'] == "RNA":
                fields[2] = "tRNA"  # Change type RNA to tRNA
                entry['type'] == "tRNA"
            if 'transcript_id' in entry['attributes']:
                entry['attributes']['transcript_id'] = f"{formbioid}.t1"
            if entry['type'] == "gene":
                entry['attributes']['ID'] = f"{formbioid}"
                if annotationTool == "tRNAscan-SE":
                    entry['attributes']['gene_biotype'] = "tRNA"
            elif entry['type'] == "pseudogene":
                entry['attributes']['ID'] = f"{formbioid}"
                entry['type'] = "pseudogene"
                fields[2] = "pseudogene"
            elif entry['type'] == "RNA" or entry['type'] == "tRNA":
                entry['attributes']['ID'] = f"{formbioid}.t1"
                entry['attributes']['Parent'] = f"{formbioid}"
            elif entry['type'] == "mRNA":
                entry['attributes']['ID'] = f"{formbioid}.t1"
                entry['attributes']['Parent'] = f"{formbioid}"
            elif entry['type'] == "Transcript":
                entry['attributes']['ID'] = f"{formbioid}.t1"
                entry['attributes']['Parent'] = f"{formbioid}"
            elif entry['type'] == "transcript":
                entry['attributes']['ID'] = f"{formbioid}.t1"
                entry['attributes']['Parent'] = f"{formbioid}"
            if 'gene_biotype' in entry['attributes']:
                if entry['attributes']['gene_biotype'] == "pseudogene":
                    entry['attributes']['gene_biotype'] = "pseudogene"
                elif 'anticodon' in entry['attributes']:
                    if entry['attributes']['gene_biotype'] != "pseudogene":
                        entry['attributes']['gene_biotype'] = "tRNA"
                elif annotationTool == "tRNAscan-SE":
                    if entry['attributes']['gene_biotype'] != "pseudogene":
                        entry['attributes']['gene_biotype'] = "tRNA"
                elif entry['attributes']['gene_biotype'] != "tRNA":
                    if entry['attributes']['gene_biotype'] != "pseudogene":
                        entry['attributes']['gene_biotype'] = "protein_coding"
            else:
                if annotationTool == "tRNAscan-SE":
                    print(entry)
                    if entry['attributes']['gene_biotype'] != "pseudogene":
                        entry['attributes']['gene_biotype'] = "tRNA"
                else:
                    entry['attributes']['gene_biotype'] = 'protein_coding'

            if 'gene_name' not in entry['attributes']:
                entry['attributes']['gene_name'] = f"{formbioid}"
            if 'gene_description' not in entry['attributes']:
                entry['attributes']['gene_description'] = ""
            if 'eggnog_ortholog' not in entry['attributes']:
                entry['attributes']['eggnog_ortholog'] = ""
            if 'eggnog_taxscope' not in entry['attributes']:
                entry['attributes']['eggnog_taxscope'] = ""
                    
            # Handling exons with numbering
            # if annotationTool != "tRNAscan-SE":
            if entry['type'] == "exon":
                entry['attributes']['Parent'] = f"{formbioid}.t1"
                currentID = entry['attributes']['ID']
                # print(f"{annotationTool}: {currentID}")
                # print(annotationTool)
                if annotationTool == "Helixer":
                    # Splitting the string by '.exon.' and taking the part after it
                    exonNum = currentID.split('.exon.', 1)[-1]
                    # Then, finding the substring before the first quote
                    # exonNum = currentID_exon.split('";', 1)[0]
                elif annotationTool == "Liftoff":
                    # Splitting the string by '.exon.' and taking the part after it
                    currentID_exon = currentID.split('exon-', 1)[-1]
                    currentID_exon = currentID_exon.split('"', 1)[0]
                    # Then, finding the substring before the first quote
                    exonNum = currentID_exon.split('-', 1)[-1]
                    if 'exon' in exonNum:
                        exonNum = currentID.split('.exon', 1)[-1]
                elif annotationTool == "tRNAscan-SE":
                    # Splitting the string by '.exon.' and taking the part after it
                    exonNum = currentID.split('.exon', 1)[-1]
                else:
                    # print(currentID)
                    # Splitting the string by '.exon.' and taking the part after it
                    exonNum = currentID.split('.exon', 1)[-1]
                    # Then, finding the substring before the first quote
                    # exonNum = currentID_exon.split('";', 1)[0]
                    # print(currentID_exon)
                    # print(exonNum)

                if firstExon == -1:
                    firstExon = int(exonNum)
                    if firstExon == 1:
                        incrementUp = True
                    else:
                        incrementUp = False
                entry['attributes']['exon_number'] = str(exonNum)  # Updating the exon number
                entry['attributes']['ID'] = f"{formbioid}.t1.exon{exonNum}"
                # last_exon_number = current_exon_number  # Keep track of the last seen exon number for direction determination
            elif entry['type'] == "CDS":
                if CDSFound == False:
                    CDSNum = firstExon
                    CDSFound = True
                    if CDSNum > 1:
                        incrementUp = False
                    else:
                        incrementUp = True
                else:
                    if incrementUp == True:
                        CDSNum += 1
                    else:
                        CDSNum += -1
                entry['attributes']['Parent'] = f"{formbioid}.t1"
                entry['attributes']['cds_number'] = str(CDSNum)  # Updating the exon number
                entry['attributes']['ID'] = f"{formbioid}.t1.cds{CDSNum}"
            # Update attributes string and reconstruct the raw line
            attributes_str = '; '.join([f'{k} "{v}"' for k, v in entry['attributes'].items()])
            fields[8] = attributes_str + ';'
            entry['raw_line'] = '\t'.join(fields)

        updated_annotations[formbioid] = entries
        formbioid_counter += 1  # Increment for the next gene

    return updated_annotations

updated_annotations_2 = reformat_and_number_exons(annotations_with_functional_info)

# Correcting the error in the code and directly writing the lines from the 'raw_line' attribute since it already contains the updated information.
with open(output_file, 'w') as file:
    for gene_entries in updated_annotations_2.values():
        for entry in gene_entries:
            file.write(entry['raw_line'] + '\n')  # Writing the raw_line directly, it's already updated in previous steps


