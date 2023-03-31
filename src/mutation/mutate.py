import sys
import uuid

import pandas as pd
from Bio.Seq import Seq

from fasta_reader import read_wuhan_1
from utils import get_mutated_phenotype


def mutation(original_phenotype, original_genotype, genome_df):
    next_gen = pd.DataFrame(
        columns=[
            "id",
            "genome",
            "generation",
            "parent",
            "mut_pos",
            "mutated_protein_index",
        ]
    )
    nucleotides = ["A", "C", "G", "T"]
    df_loc = 0
    total_revert = 0
    total_dupes = 0
    full_revert = 0
    iters = 0
    total_muts = 0

    for row in genome_df.itertuples(index=False):
        # get the parent genome to mutate
        genome = row[1]
        iters += 1
        # get where the parent genome has been previously mutated
        pos = row[4]
        # get the parent genome's generatation
        generation = row[2]
        # get the parent genome's id
        parent_id = row[0]

        # mutation
        for i in range(len(genome)):
            if genome[i] not in nucleotides:
                print("ERROR: Unexpected nucleotide value: ", genome[i])
                sys.exit(1)

            if not pos:
                mut_pos = [i]
            else:
                mut_pos = pos[:]
                mut_pos.append(i)

            for nuc in nucleotides:
                # only create a mutation if it is different from original
                if genome[i] is not nuc:
                    total_muts += 1
                    # creates the (mut)ated genome
                    if i == 0:
                        mut = nuc + genome[1:]
                    elif i == len(genome):
                        mut = genome[: len(genome)] + nuc
                    else:
                        mut = genome[:i] + nuc + genome[i + 1:]

                     # check if mutation is reversion and don't add if it is
                    if mut == original_genotype:
                        full_revert += 1
                        continue

                    if mut in genome_df["genome"].values:
                        total_revert += 1
                        continue

                    if mut in next_gen["genome"].values:
                        parent_list = (
                            next_gen[
                                next_gen.genome == mut].parent.item() + ", "
                                                                        "" + parent_id
                        )
                        next_gen.loc[
                            next_gen.genome == mut, "parent"] = parent_list
                        total_dupes += 1
                        continue

                    # id for this genome, the genome, the generation of the
                    # genome, and its parent's id
                    temp = [
                        str(uuid.uuid1().int),
                        mut,
                        generation + 1,
                        str(parent_id),
                        mut_pos,
                        None,
                    ]
                    temp = get_mutated_phenotype(original_phenotype, temp)
                    next_gen.loc[df_loc] = temp
                    df_loc += 1
    print("total partial reversions: ", total_revert)
    print("total duplicates: ", total_dupes)
    print("total reversion to original: ", full_revert)
    print("total iterations: ", iters)
    print("total mutations: ", total_muts)
    return next_gen


def test_mutation(orig):

    orig_id = str(uuid.uuid1().int)

    original_phenotype = Seq(orig).translate()
    print(f"\n\noriginal_phenotype: {original_phenotype}")
    gen0 = pd.DataFrame(
        columns=[
            "id",
            "genome",
            "generation",
            "parent",
            "mut_pos",
            "mutated_protein_index",
        ]
    )
    gen0.loc[0] = [orig_id, orig, 0, "none", [], [-1]]

    gen1 = mutation(original_phenotype, orig, gen0)

    # test_gen2 = mutation(original_phenotype, orig, gen1)
    # print("full mutation of gen 2 is: ", test_gen2.shape)
    # test_gen3 = mutation(original_phenotype,orig, test_gen2)
    # print("full mutation of gen 3 is: ", test_gen3.shape)


    print("neutral network to follow")
    neutral1 = gen1.loc[gen1["mutated_protein_index"].isin([[-1]])]

    print("GEN1 ", gen1.shape)
    print(neutral1.shape)

    gen2 = mutation(original_phenotype, orig, neutral1)
    # print(tabulate(gen2, headers="keys", tablefmt="psql"))
    print("GEN2 ", gen2.shape)

    neutral2 = gen2.loc[gen2["mutated_protein_index"].isin([[-1]])]
    print(neutral2.shape)

    gen3 = mutation(original_phenotype, orig, neutral2)
    # print(tabulate(gen3, headers="keys", tablefmt="psql"))
    print("GEN3 ", gen3.shape)

    neutral3 = gen3.loc[gen3["mutated_protein_index"].isin([[-1]])]
    # print(tabulate(neutral3, headers="keys", tablefmt="psql"))
    print(neutral3.shape)

    gen4 = mutation(original_phenotype, orig, neutral3)
    # print(tabulate(gen3, headers="keys", tablefmt="psql"))
    print("GEN4 ", gen4.shape)

    neutral4 = gen4.loc[gen4["mutated_protein_index"].isin([[-1]])]
    # print(tabulate(neutral3, headers="keys", tablefmt="psql"))
    print(neutral4.shape)



if __name__ == "__main__":
    # orig = "TACCATGGAATTACTGCG"  # Phenotype: YHGITA
    # test_mutation(orig)
    _, rbd = read_wuhan_1(r"/genomes/wuhan-hu-1.txt")

    # cleaning up RBD here
    # TODO: fix fasta_reader to remove whitespace/newlines
    newrbd = []
    for i in range(len(rbd)):
        if rbd[i] != "\n":
          newrbd.append(rbd[i])

    rbd_final = ''.join(str(i) for i in newrbd)
    print(rbd_final)
    rbd_loop = rbd_final[372:399]
    print(len(rbd_loop))
    # print(len(rbd_loop))
    print(rbd_loop)


    test_mutation(rbd_loop)
