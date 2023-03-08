import argparse


def get_genome_string(filepath):
    """
    Reads a fasta txt file, prints its header, and returns a string
    containing the genome
    :param filepath: real string containing the filepath to the fasta txt file
    :return: a string containing the genome from the fasta txt file
    """
    file = open(filepath, "r")
    all = file.readlines()
    print(all[0])
    genome = all[1:]
    return "".join(genome)


def read_wuhan_1(filepath):
    """
    returns the wuhan-hu-1 genome broken into three parts
    :param filepath:
    :return: hu1_full_genome, utr5, hu1_gene
    """
    hu1_full_genome = get_genome_string(filepath).upper()

    # spike protein gene is 21563..25384
    # RBD may be 22517-23185 start point is 1 less b/c of 0 index, end is same
    # b/c it is not included
    hu1_rbd = hu1_full_genome[22516:23185]

    return hu1_full_genome, hu1_rbd


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-filepath",
        type=str,
        required=False,
        default=r"D:\CS523\CAS-523-Proj-2\genomes\wuhan-hu-1.txt",
        help="fully qualified filepath to the fasta file you are reading",
    )
    parser.add_argument(
        "-opt",
        type=int,
        default=0,
        help="Which virus are " "you reading? 0 is wuhan-hu-1",
    )
    return parser


if __name__ == "__main__":
    parser = build_parser()
    options, _ = parser.parse_known_args()

    if options.opt == 0:
        hu1_full_genome, hu1_rbd = read_wuhan_1(options.filepath)
        # test confirmed that reading the wuhan-hu-1 file returns the
        # appropriate length strings
        print(len(hu1_full_genome))
        print(len(hu1_rbd))
