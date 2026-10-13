import csv
from pathlib import Path

# -------------------------PATH OF INTERPROSCAN RESULT--------------------------
interpro_folder = Path(r'C:\Users\DELL\Desktop\interproscan')
input_file = interpro_folder / 'tp53_human.tsv'

#----------------------file to store our reults of tsv file in clean form in csv file-------------------
output_file = interpro_folder / 'tp53_human.csv'

# ---------------------read interproscan result----------------------------

domains = []

with open(input_file, 'r', encoding='utf-8') as file:

    reader = csv.reader(file, delimiter='\t')
#------------------go through all rows in file------------------
    for row in reader:
# -------------------------Skip empty rows----------------------
        if not row:
            continue
# InterProScan columns
        protein_id = row[0]
        sequence_md5 = row[1]
        sequence_length = row[2]
        analysis = row[3]
        signature_accession = row[4]
        signature_description = row[5]
        start = row[6]
        end = row[7]
        score = row[8]
        status = row[9]
        date = row[10]
        interpro_accession = row[11]
        interpro_description = row[12]
# Only keep rows that actually have an InterPro domain
        if interpro_accession != '-':

            domains.append({
                'protein_id': protein_id,
                'sequence_length': sequence_length,
                'analysis': analysis,
                'signature_accession': signature_accession,
                'signature_description': signature_description,
                'start': start,
                'end': end,
                'interpro_accession': interpro_accession,
                'interpro_description': interpro_description
            })
# ---------------------- PRINT DOMAIN INFORMATION------------------

print("=" * 70)
print("TP53 PROTEIN DOMAIN RESULTS")
print("=" * 70)

print(f"Total domain hits found: {len(domains)}")

for domain in domains:

    print("\nProtein:", domain['protein_id'])
    print("Analysis:", domain['analysis'])
    print("Signature:", domain['signature_accession'])
    print("Description:", domain['signature_description'])
    print("Position:", domain['start'], "-", domain['end'])
    print("InterPro ID:", domain['interpro_accession'])
    print("InterPro description:", domain['interpro_description'])

# --------------------------SAVE CLEAN DOMAIN RESULTS AS CSV-------------------------

with open(output_file, 'w', newline='', encoding='utf-8') as file:

    fieldnames = [
        'protein_id',
        'sequence_length',
        'analysis',
        'signature_accession',
        'signature_description',
        'start',
        'end',
        'interpro_accession',
        'interpro_description'
    ]

    writer = csv.DictWriter(
        file,
        fieldnames=fieldnames
    )

    writer.writeheader()

    for domain in domains:
        writer.writerow(domain)


print("\n" + "-" * 70)
print("DOMAIN ANALYSIS COMPLETE")
print("-" * 70)

print(f"\nClean domain file saved at:")
print(output_file)