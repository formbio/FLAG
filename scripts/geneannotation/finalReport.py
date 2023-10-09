import datetime
import re
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--input_stats',dest='input_stats',required=True)
parser.add_argument('--lineage',dest='lineage')
parser.add_argument('--scientific_name',dest='scientific_name')
args = parser.parse_args()

lineage=str(args.lineage)
scientific_name=str(args.scientific_name)
file_path=str(args.input_stats)

# Get the current date and time
current_datetime = datetime.datetime.now()
formatted_datetime = current_datetime.strftime("%Y-%m-%d %H:%M:%S")
with open('busco.txt', 'r') as busco_file:
    busco_data = busco_file.readlines()

busco_results = []
parsing_results = False

for line in busco_data:
    if "***** Results: *****" in line:
        parsing_results = True
        continue
    if parsing_results and line.strip():
        busco_results.append(line.strip())

# Extract specific data
busco_version = busco_results[0]
complete_buscos = busco_results[1].split('\t')[0]
complete_single_copy_buscos = busco_results[2].split('\t')[0]
complete_duplicated_buscos = busco_results[3].split('\t')[0]
fragmented_buscos = busco_results[4].split('\t')[0]
missing_buscos = busco_results[5].split('\t')[0]
total_busco_groups = busco_results[6].split('\t')[0]


# Calculate percentages for Busco categories
complete_single_copy_percentage = (int(complete_single_copy_buscos) / int(total_busco_groups)) * 100
complete_duplicated_percentage = (int(complete_duplicated_buscos) / int(total_busco_groups)) * 100
fragmented_percentage = (int(fragmented_buscos) / int(total_busco_groups)) * 100
missing_percentage = (int(missing_buscos) / int(total_busco_groups)) * 100

# Format the percentages for HTML style attribute
complete_single_copy_percentage_style = f"{complete_single_copy_percentage}%"
complete_duplicated_percentage_style = f"{complete_duplicated_percentage}%"
fragmented_percentage_style = f"{fragmented_percentage}%"
missing_percentage_style = f"{missing_percentage}%"

#stats_data = [line.strip() for line in stats_data if line.strip()]  # Remove empty lines

# Extract specific data from stats.txt
#number_of_genes = stats_data[2].split()[-1]
#number_of_transcripts = stats_data[3].split()[-1]
#number_of_cds = stats_data[4].split()[-1]

# Create a regex pattern to match lines with numbers
number_pattern = r'Number of (.+?)\s+(\d+)'
mean_pattern = r'mean (.+?)\s+([\d.]+)'
total_pattern = r'Total (.+?)\s+(\d+)'

# Initialize dictionaries to store the results
results1_number = {}
results1_mean = {}
results1_total = {}
results2_number = {}
results2_mean = {}
results2_total = {}
results3_number = {}
results3_mean = {}
results3_total = {}
results4_number = {}
results4_mean = {}
results4_total = {}

# Open the file and read its contents
with open(file_path, 'r') as file:
    text = file.read()
# Split the text into sections using the "-----" separator
sections = text.split('ompute ')

# Process each section
for section in sections:
    lines = section.strip().split('\n')
    section_name = lines[0].strip()
    section_text = '\n'.join(lines[1:])

    # Use regex to extract numbers and store them in the appropriate dictionary
    if "mrna with isoforms" in section_name:
        print(section_name)
        print('1111111111')
        results1_number[section_name] = dict(re.findall(number_pattern, section_text))
        results1_mean[section_name] = dict(re.findall(mean_pattern, section_text))
        results1_total[section_name] = dict(re.findall(total_pattern, section_text))
    elif "mrna without isoforms" in section_name:
        print(section_name)
        print('2222222222')
        results2_number[section_name] = dict(re.findall(number_pattern, section_text))
        results2_mean[section_name] = dict(re.findall(mean_pattern, section_text))
        results2_total[section_name] = dict(re.findall(total_pattern, section_text))
    elif "trna with isoforms" in section_name:
        print(section_name)
        print('3333333333')
        results3_number[section_name] = dict(re.findall(number_pattern, section_text))
        results3_mean[section_name] = dict(re.findall(mean_pattern, section_text))
        results3_total[section_name] = dict(re.findall(total_pattern, section_text))
    elif "trna without isoforms" in section_name:
        print(section_name)
        print('4444444444')
        results4_number[section_name] = dict(re.findall(number_pattern, section_text))
        results4_mean[section_name] = dict(re.findall(mean_pattern, section_text))
        results4_total[section_name] = dict(re.findall(total_pattern, section_text))

# test if the AGAT Stats file is formatted differently
if results1_number == {}:
    sections = text.split('---------- ')
    # Process each section
    for section in sections:
        lines = section.strip().split('\n')
        section_name = lines[0].strip()
        section_text = '\n'.join(lines[1:])
        print(section_name)
        print(section_text)
        # Use regex to extract numbers and store them in the appropriate dictionary
        if "mrna" in section_name:
            results1_number[section_name] = dict(re.findall(number_pattern, section_text))
            results1_mean[section_name] = dict(re.findall(mean_pattern, section_text))
            results1_total[section_name] = dict(re.findall(total_pattern, section_text))
        elif "trna" in section_name:
            results3_number[section_name] = dict(re.findall(number_pattern, section_text))
            results3_mean[section_name] = dict(re.findall(mean_pattern, section_text))
            results3_total[section_name] = dict(re.findall(total_pattern, section_text))
        

# Print the results
for section_name, data in results1_number.items():
    for key, value in data.items():
        #print(f"{key}: {value}")
        if key == 'gene':
            number_of_isoform_genes = value
        elif key == 'mrna':
            number_of_isoform_mrna = value
        elif key == 'single exon gene':
            number_of_single_exon_isoform_mrna = value
        elif key == 'single exon mrna':
            number_of_single_exon_isoform_mrna = value
            
for section_name, data in results1_mean.items():
    for key, value in data.items():
        #print(f"{key}: {value}")
        if key == 'mrnas per gene':
            mean_mrnas_per_gene = value
        elif key == 'exons per mrna':
            mean_exons_per_mrna = value
        elif key == 'gene length' or key == 'gene length (bp)':
            mean_gene_length = value
            
#print("\nResults for mrna without isoforms:")
#for section_name, data in results2.items():
#    print(section_name)
#    for key, value in data.items():
#        print(f"{key}: {value}")

print("\nResults for trna with isoforms:")
for section_name, data in results3_number.items():
    for key, value in data.items():
        #print(f"{key}: {value}")
        if key == 'trna':
            number_of_isoform_trna = value

#print("\nResults for trna without isoforms:")
#for section_name, data in results4.items():
#    print(section_name)
#    for key, value in data.items():
#        print(f"{key}: {value}")



html_template = """
<!DOCTYPE html>
<html>
<head>
    <title>FLAG Annotation Report</title>
    <style>
        h1 {{
            font-family: Arial, sans-serif;
            font-size: 36px;
            text-align: center;
            padding: 20px;
            background-color: #3498db; /* Background color */
            color: #fff; /* Text color */
            border-radius: 10px; /* Rounded corners */
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1); /* Shadow */
        }}
       .percentage-bar-container {{
           width: 550px;
           display: flex; /* or display: inline-block; */
       }}

       .percentage-bar {{
           height: 20px;
           background-color: lightgray;
           margin-right: 0px; /* Add some spacing between bars */
       }}
       
      .percentage-number {{
           position: absolute;
           top: 0;
           left: 50%;
           transform: translateX(-50%);
       }}
       .complete-single-copy {{
           background-color: blue;
       }}

       .complete-duplicated {{
           background-color: #00008B;
       }}

       .fragmented {{
           background-color: yellow;
       }}

       .missing {{
           background-color: red;
       }}
        .inline-div {{
            display: flex;
            justify-content: center;
           align-items: center;
           width: 550px;
           height: 25px;
        }}
        .label-box {{
           width: 20px;
           height: 20px;
           margin-right: 5px;
           margin-left: 5px;
           display: inline-block;
       }}

       .complete-single-copy-label {{
           background-color: blue;
       }}

       .complete-duplicated-label {{
           background-color: #00008B;
       }}

       .fragmented-label {{
           background-color: yellow;
       }}

       .missing-label {{
           background-color: red;
       }}
   </style>
</head>
<body>
   <h1>FLAG: Find, Label, Annotate, Genes Annotation Report</h1>
   <h2>Annotation of {7}</h2>
   <p>Completed on {8}</p>
   <h2>Busco Results</h2>
   <p>Busco Version: 5.3.2</p>
   <div class="inline-div">
    <div class="label-box complete-single-copy-label"></div>Complete (C) and single-copy (S)
   <div class="label-box complete-duplicated-label"></div>Complete (C) and duplicatex (D)
   </div>
   <div class="inline-div">
    <div class="label-box fragmented-label"></div>Fragmented (F)
   <div class="label-box missing-label"></div>Missing (M)
   </div>
   <div class="percentage-bar-container">
       <div class="percentage-bar complete-single-copy" style="width: {1};">
       </div>
       <div class="percentage-bar complete-duplicated" style="width: {2};">
       </div>
       <div class="percentage-bar fragmented" style="width: {3};">
       </div>
       <div class="percentage-bar missing" style="width: {4};">
       </div>
   </div>
      <p>Busco Assessment: {0}</p>

   <p>Total BUSCO Groups: {5} for the {6} lineage</p>

    <h2>AGAT Statistics</h2>
    <p>Number of Protein-coding Genes: {9}</p>
    <p>Number of mRNAs: {10}</p>
    <p>Number of Single-exon Protein-coding Genes: {11}</p>
    <p>Number of Single-exon mRNA: {12}</p>
    <p>Mean Number of mRNAs per Gene: {13}</p>
    <p>Mean Number of Exons per mRNA: {14}</p>
    <p>Mean Protein-coding Gene Length: {15}</p>
    <p>Number of tRNAs: {16}</p>
    <!-- Add more placeholders for other stats data -->
</body>
</html>
""".format(
    busco_version, complete_single_copy_percentage_style,
    complete_duplicated_percentage_style, fragmented_percentage_style,
    missing_percentage_style, total_busco_groups, lineage, scientific_name, formatted_datetime, number_of_isoform_genes, number_of_isoform_mrna, number_of_single_exon_isoform_mrna, number_of_single_exon_isoform_mrna, mean_mrnas_per_gene, mean_exons_per_mrna, mean_gene_length, number_of_isoform_trna
)

with open('results.html', 'w') as output_file:
    output_file.write(html_template)
