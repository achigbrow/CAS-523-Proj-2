import uuid

import pandas as pd

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
    nucleotides = ['A', 'C', 'G', 'T']
    df_loc = 0

    for i in range(len(genome)):

        if i not in skip:
            for nuc in nucleotides:
                # only create a mutation if it is different from original
                if genome[i] is not nuc:
                    # adds id for the mutated genome
                    temp = [uuid.uuid1()]
                    # creates the (mut)ated genome
                    if i == 0:
                        mut = nuc + genome[1:]
                    elif i == len(genome):
                        mut = genome[:len(genome)] + nuc
                    else:
                        mut = genome[:i] + nuc + genome[i+1:]
                    temp.append(mut)
                    # adds the generation according to how many mutations are
                    # accumulated
                    temp.append(generation+1)
                    # adds the id of the parent genome
                    temp.append(id)
                    # adds the pos that is mutated so that it is skipped in
                    # future generations to accumulate mutations
                    skip.append(i)
                    temp.append(skip)
                    df[df_loc] = temp
                    df_loc += 1

    return df


