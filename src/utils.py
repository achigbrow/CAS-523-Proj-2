"""
Utility functions
"""
from Bio.Seq import Seq


def get_mutated_phenotype(original_phenotype, values):
    """
    :param original_phenotype: string, the original phenotype of the genome sqeuence of A, C, T, and Gs
    :param values: a list with the following values, in order:
        Expects columns: "id", "genome", "generation", "parent", "mut_pos", "mutated_protein_index"
        These are columns in the final genome mutation dataframe.

    :return: a list that contains an updated "mutated_protein_index" column.
    """

    genome = Seq(values[1])
    phenotype = genome.translate()
    # print(f"genome: {genome} phenotype: {phenotype}")
    if phenotype == original_phenotype:
        values[5] = [-1]
        return values

    diff = [i for i in range(len(original_phenotype)) if original_phenotype[i] != phenotype[i]]

    values[5] = diff

    return values
