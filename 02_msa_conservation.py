import csv
from collections import Counter
from pathlib import Path

from Bio import AlignIO

# CHANGED: relative folders instead of C:\Users\DELL\Desktop\MSA
data_folder = results_folder = Path(r'C:\Users\DELL\Desktop\MSA')  # the folder this script is in

msa_file = data_folder / 'TP53_MSA.fasta'

# ---------------- read MSA file ----------------
alignment = AlignIO.read(msa_file, 'fasta')

MIN_COUNT = 3  # a position is kept only if at least 3 sequences share the same amino acid


# CHANGED: new helper, uses >= instead of exact == checks, so gaps cannot break the label
def conservation_label(percentage):
    if percentage >= 100:
        return 'Highly conserved'
    elif percentage >= 80:
        return 'Moderately conserved'
    else:
        return 'Slightly conserved'


# ---------------- find conserved positions ----------------
def find_conserved_positions(alignment):
    conserved_positions = []
    seq_length = alignment.get_alignment_length()
    for i in range(seq_length):
        amino_acids = [record.seq[i] for record in alignment]
        non_gap_aa = [aa for aa in amino_acids if aa != '-']
        if not non_gap_aa:
            continue
        aa_counts = Counter(non_gap_aa)
        most_common_aa, count = aa_counts.most_common(1)[0]
        percentage = (count / len(non_gap_aa)) * 100
        if count >= MIN_COUNT:
            conserved_positions.append({
                'position': i + 1,
                'amino_acid': most_common_aa,
                'percentage': round(percentage, 2),
                'conservation': conservation_label(percentage)  # CHANGED: always defined
            })
    return conserved_positions


# ---------------- find conserved regions ----------------
def find_conserved_regions(conserved_positions, min_conserved=5):
    conserved_regions = []
    current_region = []
    for pos in conserved_positions:
        if pos['percentage'] < 80:
            if len(current_region) >= min_conserved:
                conserved_regions.append(current_region)
            current_region = []
            continue
        if not current_region:
            current_region.append(pos)
        elif pos['position'] == current_region[-1]['position'] + 1:
            current_region.append(pos)
        else:
            if len(current_region) >= min_conserved:
                conserved_regions.append(current_region)
            current_region = [pos]
    if len(current_region) >= min_conserved:
        conserved_regions.append(current_region)
    return conserved_regions


# ---------------- map MSA columns to human TP53 residue numbers ----------------
def map_msa_to_human(alignment):
    human_record = None
    for record in alignment:
        if 'human' in record.id.lower():
            human_record = record
            break
    if human_record is None:
        raise ValueError('Human TP53 sequence was not found in the MSA.')

    msa_to_human = {}
    human_residue_number = 0
    for msa_position, amino_acid in enumerate(human_record.seq, start=1):
        if amino_acid != '-':
            human_residue_number += 1
            msa_to_human[msa_position] = human_residue_number
    return msa_to_human


def map_conserved_regions_to_human(conserved_regions, msa_to_human):
    mapped_regions = []
    for region in conserved_regions:
        human_positions = []
        for pos in region:
            if pos['position'] in msa_to_human:
                human_positions.append(msa_to_human[pos['position']])
        if human_positions:
            mapped_regions.append({
                'msa_start': region[0]['position'],
                'msa_end': region[-1]['position'],
                'human_start': human_positions[0],
                'human_end': human_positions[-1],
                'length': len(human_positions),
                'sequence': ''.join(pos['amino_acid'] for pos in region)
            })
    return mapped_regions


# ---------------- run (CHANGED: each function is called once) ----------------
conserved_positions = find_conserved_positions(alignment)
conserved_regions = find_conserved_regions(conserved_positions, min_conserved=5)
msa_to_human = map_msa_to_human(alignment)
mapped_regions = map_conserved_regions_to_human(conserved_regions, msa_to_human)

print('Conserved positions')
for pos in conserved_positions:
    print(
        f"Position: {pos['position']} | "
        f"Amino acid: {pos['amino_acid']} | "
        f"Percentage: {pos['percentage']} | "
        f"Conservation: {pos['conservation']}"
    )

print('\nConserved regions (MSA numbering)')
for region in conserved_regions:
    start = region[0]['position']
    end = region[-1]['position']
    sequence = ''.join(pos['amino_acid'] for pos in region)
    print(f'Region: {start}-{end} | Length: {len(region)} | Sequence: {sequence}')

print('\nCONSERVED REGIONS -> HUMAN TP53')
print('=' * 60)
for region in mapped_regions:
    print(
        f"MSA positions: {region['msa_start']}-{region['msa_end']} "
        f"-> Human TP53 residues: {region['human_start']}-{region['human_end']} "
        f"({region['length']} residues)"
    )

# ---------------- save CSV files ----------------
positions_file = results_folder / 'TP53_conserved_positions.csv'
with open(positions_file, 'w', newline='') as out_file:  # CHANGED: out_file instead of file
    writer = csv.DictWriter(out_file, fieldnames=['position', 'amino_acid', 'percentage', 'conservation'])
    writer.writeheader()
    writer.writerows(conserved_positions)
print(f'\nConserved positions saved to: {positions_file}')

regions_file = results_folder / 'TP53_conserved_regions.csv'
with open(regions_file, 'w', newline='') as out_file:
    writer = csv.DictWriter(out_file, fieldnames=['start', 'end', 'length', 'sequence'])
    writer.writeheader()
    for region in conserved_regions:
        writer.writerow({
            'start': region[0]['position'],
            'end': region[-1]['position'],
            'length': len(region),
            'sequence': ''.join(pos['amino_acid'] for pos in region)
        })
print(f'Conserved regions (MSA numbering) saved to: {regions_file}')

# CHANGED: NEW file with human TP53 numbering, read by the structure script
human_regions_file = results_folder / 'TP53_conserved_regions_human.csv'
with open(human_regions_file, 'w', newline='') as out_file:
    writer = csv.DictWriter(
        out_file,
        fieldnames=['msa_start', 'msa_end', 'human_start', 'human_end', 'length', 'sequence']
    )
    writer.writeheader()
    writer.writerows(mapped_regions)
print(f'Conserved regions (human numbering) saved to: {human_regions_file}')
