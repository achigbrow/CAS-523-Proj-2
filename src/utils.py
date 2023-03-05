"""
Utility functions
"""


def get_mutated_phenotypes(original_genome, mutations):
    """
    :param original_genome: string, the original genome sqeuence of A, C, T, and Gs
    :param mutations: a dataframe of mutations to the original genome.
        Expects columns: "id", "genome", "generation", "parent"

    :return: dataframe that contains extra colums indicating whether the phenotype has
        changed and which protein was changed (indexed from 0, starting at the beginning
        of the original genome)
        Columns: "id", "genome", "generation", "parent", "mut_pheno", "protein_index"
    """

    raise NotImplementedError
