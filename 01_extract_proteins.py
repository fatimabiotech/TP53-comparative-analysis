import csv
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2 
from itertools import combinations
from Bio import AlignIO
def find_orf_nucleotides(sequence):
    orf_list=[]
    stop_codons = ['TAA', 'TGA', 'TAG']
    seq_str = str(sequence)
    start_pos=0
    while True:
        start_pos = seq_str.find('ATG', start_pos)
        if start_pos == -1:
            break
        for stop_pos in range(start_pos + 3, len(seq_str) - 2, 3):
            codon = seq_str[stop_pos:stop_pos + 3]
            if codon in stop_codons:
                orf_list.append((start_pos, stop_pos + 3))
                break
        start_pos += 1
    return (orf_list)
def find_longest_tp53_orf(sequence):
    orf_list=find_orf_nucleotides(sequence)
    if not orf_list:
        return None,'', 0, '', 0
    longest_coordinates = None
    longest_dna = ''
    longest_dna_length = 0
    longest_protein = ''
    longest_protein_length = 0
    all_orfs=[]
    for start, end in orf_list:
        orf_dna = sequence[start:end]
        current_length = len(orf_dna)
        protein = str(Seq(orf_dna).transcribe().translate(to_stop=True))
        protein_len = len(protein)
        all_orfs.append({
            'start': start,
            'end': end,
            'dna_length': current_length,
            'protein_length': protein_len,
            'protein': protein + '...' if len(protein) > 50 else protein
        })
        if current_length > longest_dna_length:
                longest_dna_length= current_length
                longest_coordinates= (start, end)
                longest_dna = orf_dna
                longest_protein= protein
                longest_protein_length = protein_len
    return(longest_coordinates,longest_dna, longest_dna_length, longest_protein, longest_protein_length, all_orfs)
def calculate_protein_summary(sequence):
    coordinates, dna, dna_length, protein, protein_length, all_orfs = find_longest_tp53_orf(sequence)
    return {
        'coordinates': coordinates,
        'dna_length': dna_length,
        'protein_length': protein_length,
        'protein': protein,
        'all_orfs': all_orfs,
        'total_orfs': len(all_orfs)
    }
def calculate_identity(seq1, seq2):
    alignments=pairwise2.align.globalxx(seq1, seq2)
    best_alignment=alignments[0]
    matches = 0
    for i in range(len(best_alignment[0])):
        if best_alignment[0][i] == best_alignment[1][i]:
            matches += 1
    identity=(matches/len((best_alignment[0]))) *100
    return(identity)
folder=Path(r'C:\Users\DELL\Desktop\project')
results=[]
for file in folder.glob('**/*.fasta'):
    parts = file.stem.split('_')
    gene = parts[0]      # "TP53"
    species = parts[1]
    if gene != 'TP53':
        continue 
    print(f"Processing: {species}")
    record=SeqIO.read(file, 'fasta')
    summary=calculate_protein_summary(record.seq)
    
    results.append({
        'species': species,
        'protein_length': summary['protein_length'],
        'total_orfs': summary['total_orfs'],
        'coordinates': summary['coordinates'],
        'protein': summary['protein']
    })

for result in results:
    print(
            result['species'],
            result['protein_length'],
            result['total_orfs'],
           # result['protein']
        )
human_protein = None

for result in results:
    if result['species'].lower() == 'human':
        human_protein = result['protein']
        break
if human_protein is not None:
    for result in results:
        species=result['species']
        #other_protein=result['protein']
        if species.lower() == 'human':
            identity=100.0
        else:
            identity = calculate_identity(human_protein, result['protein'])
        print(f"Human vs {species} TP53 identity: {identity:.2f}%")
proteins={}
for result in results:
    proteins[result['species']]=result['protein']
for species1, species2 in combinations(proteins.keys(), 2):
    identity = calculate_identity(
        proteins[species1],
        proteins[species2]
    )
    print(f"{species1} vs {species2}: " f"{identity:.2f}%")
#with open('TP53_proteins.fasta', 'w') as file:
    #for result in results:
        #file.write(f">{result['species']}\n")
        #file.write(f"{result['protein']}\n")