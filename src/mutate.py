import sys
import uuid

import pandas as pd

from tabulate import tabulate


def mutation(genome, generation, id, pos):
    """
    Takes a string genome and returns a dataframe containing its mutations.
    generation=0 will build a dataframe containing generation 1 with 1
    mutation per string and mut_pos no nucleotides because pos will be an empty
    array. generation=1 will build a df containing generation 2 with 2
    mutations per string but will skip the nucleotide site listed in pos.
    And generation=2  will generate the third generation with three mutations
    per string.
    :param genome: a string of the genome
    :param generation: int representation of generation of the genome being
    passed with 0 being original
    strain and each subsequent iteration indicating number of mutations in
    the generation. [0,2]
    :param id: id associated with this parent genome
    :param pos: array containing int val(s) to skip when mutating as they
    are already mutated
    generation 0 should be [], gen 1 should be
    :return: dataframe with cols: id, genome, generation, parent (id), mut_pos
    """

    df = pd.DataFrame(columns=["id", "genome", "generation", "parent", "mut_pos"])
    nucleotides = ["A", "C", "G", "T"]
    df_loc = 0

    for i in range(len(genome)):
        if genome[i] not in nucleotides:
            print("ERROR: Unexpected nucleotide value: ", genome[i])
            sys.exit(1)

        if i not in pos:
            if pos == []:
                mut_pos = [i]
            else:
                mut_pos = pos[:]
                mut_pos.append(i)

            for nuc in nucleotides:
                # only create a mutation if it is different from original
                if genome[i] is not nuc:
                    # creates the (mut)ated genome
                    if i == 0:
                        mut = nuc + genome[1:]
                    elif i == len(genome):
                        mut = genome[: len(genome)] + nuc
                    else:
                        mut = genome[:i] + nuc + genome[i + 1 :]
                    # id for this genome, the genome, the generation of the
                    # genome, and its parent's id
                    temp = [uuid.uuid1().int, mut, generation + 1, id, mut_pos]
                    df.loc[df_loc] = temp
                    df_loc += 1
    return df


def test_mutation():
    orig = "AAATGAC"
    orig_id = uuid.uuid1().int

    gen1 = mutation(orig, 0, orig_id, [])
    print(tabulate(gen1, headers="keys", tablefmt="psql"))
    print(gen1.shape)

    gen2_list = []

    for row in gen1.itertuples(index=False):
        gen2_list.append(mutation(row[1], row[2], row[0], row[4]))

    gen2 = pd.concat(gen2_list)
    gen2 = gen2.drop_duplicates(subset=["genome"])
    print(tabulate(gen2, headers="keys", tablefmt="psql"))
    print(gen2.shape)

    gen3_list = []

    for row in gen2.itertuples(index=False):
        gen3_list.append(mutation(row[1], row[2], row[0], row[4]))

    gen3 = pd.concat(gen3_list)
    gen3 = gen3.drop_duplicates(subset=["genome"])
    print(tabulate(gen3, headers="keys", tablefmt="psql"))
    print(gen3.shape)
    print("done")


if __name__ == "__main__":
    test_mutation()
