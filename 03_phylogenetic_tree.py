from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
import matplotlib.pyplot as plt
import csv
import random
from collections import Counter
#------------path of our file---------------
msa_folder=Path(r'C:\Users\DELL\Desktop\MSA')
msa_file=msa_folder/'TP53_MSA.fasta'
#----------------read MSA file---------------------
alignment = AlignIO.read(msa_file, 'fasta')
#-----------------calculate distance between multiple species protein sequences------------------
calculator=DistanceCalculator('identity') #identity method to calculate sequence distance
distance_matrix=calculator.get_distance(alignment) #to find pairwise distances b/w species
print(distance_matrix)
def bootstrap_analysis(alignment, n_replicates=100):
    calculator=DistanceCalculator('identity') 
    distance_matrix=calculator.get_distance(alignment)
    constructor=DistanceTreeConstructor()
    constructor.nj(distance_matrix)
    sequence_length=alignment.get_alignment_length()
    bootstrap_trees=[]
    for replicate in range(n_replicates):
        random_columns=[random.randint(0, sequence_length-1)
        for _ in range(sequence_length)]
        bootstrap_sequences=[]
        for record in alignment:
            sequence=''.join(
                str(record.seq[i])
                for i in random_columns
            )
            bootstrap_sequences.append(sequence)
        bootstrap_alignment = MultipleSeqAlignment([
    SeqRecord(
        Seq(sequence),
        id=alignment[j].id
    )
    for j, sequence in enumerate(bootstrap_sequences)
])
        distance_matrix=calculator.get_distance(bootstrap_alignment)
        tree=constructor.nj(distance_matrix)
        bootstrap_trees.append(tree)
        print(
                    f"Bootstrap replicate {replicate + 1}/{n_replicates}"
                )
        
    return bootstrap_trees
#-------------------call function-----------------
#random.seed()
bootstrap_trees = bootstrap_analysis(
    alignment,
    n_replicates=100
)
#------------now create tree constructor--------------------
constructor=DistanceTreeConstructor()
tree=constructor.nj(distance_matrix) #actually construct tree using neighbor joining(nj) alogritham
print(tree)
tree_file=msa_folder/ 'TP53_phylogenetictree.nwk'
Phylo.write(tree, tree_file, 'newick')
#---------------------------connecting bootstrap to original tree--------------
def get_clades(tree):
    clades = []
    for clade in tree.get_nonterminals():
        species = frozenset(
            terminal.name
            for terminal in clade.get_terminals())
        clades.append(species)
    return clades
original_clades = get_clades(tree)
print("Original tree clades:")
def calculate_bootstrap_support(original_tree, bootstrap_trees):
    original_clades = get_clades(original_tree)
    support_values = {}
    for clade in original_clades:
        # Ignore the clade containing all species
        if len(clade) == len(original_tree.get_terminals()):
            continue
        support_count = 0
        for bootstrap_tree in bootstrap_trees:
            bootstrap_clades = get_clades(bootstrap_tree)
            if clade in bootstrap_clades:
                support_count += 1
        support_percentage = (
            support_count / len(bootstrap_trees)
        ) * 100
        support_values[clade] = support_percentage
    return support_values
support_values = calculate_bootstrap_support(
    tree,
    bootstrap_trees
)
for clade, support in support_values.items():
    print(
        f"{set(clade)} → {support:.1f}%"
    )
#---------------save the as png-------------------
fig, ax= plt.subplots(figsize=(10, 7))
Phylo.draw(tree, 
           axes=ax,
           do_show=False)
ax.set_title('TP53_phylogenetic_tree.png')
plt.tight_layout()
png_file=msa_folder/ 'TP53_phylogenetic_tree.png'
plt.savefig(png_file,
            dpi=300,
            bbox_inches='tight')
plt.show()
#-------------csv file for distance matrix storage-----------------------
matrix_file=msa_folder/'TP53_distance_matrix.csv'
with open(matrix_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Species'] + distance_matrix.names)
    # Distance values
    for i, name in enumerate(distance_matrix.names):
        row = [name]

        for j in range(len(distance_matrix.names)):
            row.append(round(distance_matrix[i, j], 4))

        writer.writerow(row)

print(f"Distance matrix saved to: {matrix_file}")