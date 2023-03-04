import uuid

import pandas as pd
from tabulate import tabulate
from genomes import fasta_reader


def mutation(genome, generation, id, skip):
    """
    Takes a string genome and returns a dataframe containing its mutations.
    generation=0 will build a dataframe containing generation 1 with 1
    mutation per string and skip no nucleotides because skip will be an empty
    array. generation=1 will build a df containing generation 2 with 2
    mutations per string but will skip the nucleotide site listed in skip.
    And generation=2  will generate the third generation with three mutations
    per string.
    :param genome: a string of the genome
    :param generation: int representation of generation of the genome being
    passed with 0 being original
    strain and each subsequent iteration indicating number of mutations in
    the generation. [0,2]
    :param id: id associated with this parent genome
    :param skip: array containing int val(s) to skip when mutating.
    generation 0 should be [], gen 1 should be
    :return: dataframe with cols: id, genome, generation, parent (id), skip
    """

    df = pd.DataFrame(columns=["id", "genome", "generation", "parent", "skip"])
    nucleotides = ["A", "C", "G", "T"]
    df_loc = 0

    for i in range(len(genome)):
        if i not in skip:
            if skip == []:
                mut_skip = [i]
            else:
                mut_skip = skip[:]
                mut_skip.append(i)

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
                    if mut in df["genome"].unique():
                        print("duplicate", mut)
                        continue
                    # id for this genome, the genome, the generation of the
                    # genome, and its parent's id
                    temp = [uuid.uuid1().int, mut, generation + 1, id, mut_skip]
                    # print(temp)
                    df.loc[df_loc] = temp
                    df_loc += 1
    return df

def test_mutation():
    orig = "AAAT"
    orig_id = uuid.uuid1().int

    gen1 = mutation(orig, 0, orig_id, [])
    # print(tabulate(gen1, headers='keys', tablefmt='psql'))

    gen2_list = []

    for row in gen1.itertuples(index=False):
        gen2_list.append(mutation(row[1], row[2],
                                  row[0], row[4]))

    gen2 = pd.concat(gen2_list)
    # gen2 = gen2.drop_duplicates(subset=["genome", "parent"])
    print(tabulate(gen2, headers='keys', tablefmt='psql'))

    gen3_list = []

    for row in gen2.itertuples(index=False):
        gen3_list.append(mutation(row[1], row[2],
                                  row[0], row[4]))

    gen3 = pd.concat(gen3_list)
    # gen3 = gen3.drop_duplicates(subset=["genome", "parent"])
    print(tabulate(gen3, headers='keys', tablefmt='psql'))
    print("done")

if __name__ == "__main__":
    # test_mutation()
    hu1_full_genome, utr5, hu1_gene = fasta_reader.read_wuhan_1(
        r"D:\CS523\CAS-523-Proj-2\genomes\wuhan-hu-1.txt")

    orig_id = uuid.uuid1().int

    gen1 = mutation(hu1_gene, 0, orig_id, [])

    print(gen1.shape)
    print("done")