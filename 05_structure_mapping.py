import csv
from pathlib import Path

from Bio.PDB import PDBParser

# CHANGED: relative folders instead of C:\Users\DELL\Desktop\...
data_folder = results_folder = Path(r'C:\Users\DELL\Desktop\MSA') # the folder this script is in

pdb_file = data_folder / 'TP53_structure_3D.pdb'        # AlphaFold model (UniProt P04637)
interpro_file = data_folder / 'tp53_human.csv'          # output of the InterProScan script
conserved_file = results_folder / 'TP53_conserved_regions_human.csv'  # output of 02_msa_conservation.py


def to_ranges(numbers):
    """Turn [1,2,3,7,8] into '1-3, 7-8'."""
    if not numbers:
        return ''
    ranges = []
    start = prev = numbers[0]
    for n in numbers[1:]:
        if n == prev + 1:
            prev = n
            continue
        ranges.append(f'{start}-{prev}')
        start = prev = n
    ranges.append(f'{start}-{prev}')
    return ', '.join(ranges)


# ---------------- 1. read the 3D structure ----------------
parser = PDBParser(QUIET=True)
structure = parser.get_structure('TP53', pdb_file)
chain = structure[0]['A']  # TP53 has one model and one chain
residues = [r for r in chain if r.id[0] == ' ']  # standard residues only
print('PDB structure read successfully')
print('Number of residues in chain A:', len(residues))  # should be 393 for human TP53

# ---------------- 2. read conserved regions (CHANGED: from CSV, not typed by hand) ----------------
conserved_human_regions = []
with open(conserved_file, 'r', newline='', encoding='utf-8') as in_file:
    for row in csv.DictReader(in_file):
        conserved_human_regions.append((int(row['human_start']), int(row['human_end'])))

conserved_residues = set()
for start, end in conserved_human_regions:
    conserved_residues.update(range(start, end + 1))
print('\nTotal conserved human residues:', len(conserved_residues))

# ---------------- 3. read functional domains (CHANGED: one list, taken from InterProScan) ----------------
DOMAIN_LABELS = {
    'p53, transactivation domain': 'Transactivation domain',
    'Cellular tumor antigen p53, transactivation domain 2': 'Transactivation domain 2',
    'p53, DNA-binding domain': 'DNA-binding domain',
    'p53, central conserved site': 'Central conserved site',
    'p53, tetramerisation domain': 'Tetramerisation domain',
}

domain_ranges = {}
with open(interpro_file, 'r', newline='', encoding='utf-8') as in_file:
    for row in csv.DictReader(in_file):
        label = DOMAIN_LABELS.get(row['interpro_description'])
        if label is None:
            continue
        start, end = int(row['start']), int(row['end'])
        if label in domain_ranges:  # several signatures can hit the same domain: keep the widest range
            old_start, old_end = domain_ranges[label]
            domain_ranges[label] = (min(start, old_start), max(end, old_end))
        else:
            domain_ranges[label] = (start, end)

functional_domains = [(name, s, e) for name, (s, e) in domain_ranges.items()]

print('\nFUNCTIONAL DOMAINS FROM INTERPROSCAN')
print('-' * 60)
for name, start, end in functional_domains:
    print(f'{name}: {start}-{end}')

# ---------------- 4. conservation coverage per domain (CHANGED: computed, not typed in) ----------------
summary_rows = []
for name, start, end in functional_domains:
    total_residues = end - start + 1
    conserved_in_domain = sorted(r for r in conserved_residues if start <= r <= end)
    percentage = len(conserved_in_domain) / total_residues * 100
    summary_rows.append([
        name,
        f'{start}-{end}',
        to_ranges(conserved_in_domain),
        len(conserved_in_domain),
        total_residues,
        round(percentage, 2),
    ])
    print(f'\n{name}')
    print(f'Domain range: {start}-{end}')
    print(f'Total residues: {total_residues}')
    print(f'Conserved residues: {len(conserved_in_domain)}')
    print(f'Conservation coverage: {percentage:.2f}%')

summary_file = results_folder / 'TP53_final_summary.csv'
with open(summary_file, 'w', newline='', encoding='utf-8') as out_file:
    writer = csv.writer(out_file)
    writer.writerow([
        'Functional Region',
        'Functional Range',
        'Conserved Residues/Ranges',
        'Conserved Residues (n)',
        'Domain Length',
        'Conservation Coverage (%)',
    ])
    writer.writerows(summary_rows)
print(f'\nFinal summary saved to: {summary_file}')

# ---------------- 5. find the conserved residues in the 3D structure ----------------
conserved_pdb_residues = [
    (r.id[1], r.resname) for r in residues if r.id[1] in conserved_residues
]
print('\nConserved residues found in AlphaFold structure:', len(conserved_pdb_residues))

# ---------------- 6. which domain does each conserved residue belong to? ----------------
conserved_with_domains = []
for residue_number, amino_acid in conserved_pdb_residues:
    matched = [
        name for name, start, end in functional_domains
        if start <= residue_number <= end
    ]
    conserved_with_domains.append((residue_number, amino_acid, matched))

residue_file = results_folder / 'TP53_conserved_residues_3D_functional.csv'
with open(residue_file, 'w', newline='', encoding='utf-8') as out_file:
    writer = csv.writer(out_file)
    writer.writerow(['Residue', 'Amino_acid', 'Functional_domains'])
    for residue_number, amino_acid, matched in conserved_with_domains:
        writer.writerow([residue_number, amino_acid, '; '.join(matched)])
print(f'3D functional residue results saved to: {residue_file}')

inside_residue_numbers = [r for r, aa, m in conserved_with_domains if m]
outside_residue_numbers = [r for r, aa, m in conserved_with_domains if not m]

print('\nResidues to highlight inside functional domains:')
print(to_ranges(inside_residue_numbers))
print('\nConserved residues outside functional domains:')
print(to_ranges(outside_residue_numbers))

# ---------------- 7. overlap between conserved regions and functional domains ----------------
overlap_results = []
for conserved_start, conserved_end in conserved_human_regions:
    for domain_name, domain_start, domain_end in functional_domains:
        if conserved_start <= domain_end and conserved_end >= domain_start:
            overlap_start = max(conserved_start, domain_start)
            overlap_end = min(conserved_end, domain_end)
            overlap_results.append([
                conserved_start, conserved_end,
                domain_name, domain_start, domain_end,
                overlap_start, overlap_end,
            ])

print('\n' + '-' * 80)
print('CONSERVED TP53 REGIONS OVERLAPPING FUNCTIONAL REGIONS')
print('-' * 80)

# CHANGED: every row is now unpacked from its own result, so each row prints its own values
for result in overlap_results:
    c_start, c_end, d_name, d_start, d_end, o_start, o_end = result
    print(f'\nConserved region: {c_start}-{c_end}')
    print(f'Functional region: {d_name}')
    print(f'Functional range: {d_start}-{d_end}')
    print(f'Exact overlap: {o_start}-{o_end}')

overlap_file = results_folder / 'TP53_conserved_functional_overlap.csv'
with open(overlap_file, 'w', newline='', encoding='utf-8') as out_file:
    writer = csv.writer(out_file)
    writer.writerow(['Conserved_region', 'Functional_region', 'Functional_range', 'Overlap', 'Overlap_length'])
    for result in overlap_results:  # CHANGED: each CSV row uses its own values
        c_start, c_end, d_name, d_start, d_end, o_start, o_end = result
        writer.writerow([
            f'{c_start}-{c_end}',
            d_name,
            f'{d_start}-{d_end}',
            f'{o_start}-{o_end}',
            o_end - o_start + 1,
        ])
print(f'\nOverlap results saved to: {overlap_file}')
