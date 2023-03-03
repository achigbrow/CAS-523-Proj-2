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
    return ''.join(genome)


def seperate_genome(genome, start, end):
    """
    This function allows you to seperate out specific parts of the genome for
    analysis.
    Make sure that your start/end account for the fact that the string is
    indexed from 0
    :param genome: a string of the genome
    :param start: the first nucleotide of the section of the genome you want to
    analyze
    (will be included in the result)
    :param end: the last nucleotide in the section of the genome that you
    want to analyze (will be included in the result)
    :return: a string containing a fraction of the genome
    """
    # TODO: figure out if I need to do anything else or if it dumb to have a
    #  one line function omg
    return genome[start : end + 1]


def read_wuhan_1(filepath):
    """
    returns the wuhan-hu-1 genome broken into three parts
    :param filepath:
    :return: hu1_full_genome, utr5, hu1_gene
    """
    hu1_full_genome = get_genome_string(filepath)

    utr5 = seperate_genome(hu1_full_genome, 0, 264)
    hu1_gene = seperate_genome(hu1_full_genome, 265, 21554)

    return hu1_full_genome, utr5, hu1_gene


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-filepath",
        type=str,
        required=False,
        default=r"D:\CS523\CAS-523-Proj-2\original\wuhan-hu-1.txt",
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
        hu1_full_genome, utr5, hu1_gene = read_wuhan_1(options.filepath)
        # test confirmed that reading the wuhan-hu-1 file returns the
        # appropriate length strings
        print(len(hu1_full_genome))
        print(len(utr5))
        print(len(hu1_gene))
