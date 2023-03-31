from collections import Counter

def count_aas(protein):
    """
    creates a histogram of aa's present in the protein
    :param protein: a string form of a protein
    :return:
    """
    aas = ["Ala", "Arg", "Asn", "Asp", "Cys", "Gin", "Glu", "Gly", "His",
           "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp",
           "Tyr", "Val"]

    aa_abbrevs = ["A", "R", "N", "D", "C", "E", "Q", "G", "H","I",
                  "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

    count = Counter(protein)
    length = len(protein)
    freqs = {}
    print(count)
    print(sum(count.values()))
    print(length)
    for aa in aa_abbrevs:
       freqs[aa] = round(count[aa]/length, 5)

    print(freqs)
    print(freqs.values())
    print(sum(freqs.values()))
    return count, freqs

def get_s_substring(start, end):
    """
    returns some subsection of the Spike protein. For example if you want the
    RBD, pass 319, 541
    :param start: first residue number in the spike protein section you want
    :param end: final residue number in the spike protein section you want
    :return:
    """
    s = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVS" \
        "GTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLG" \
        "VYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLV" \
        "RDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGT" \
        "ITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNR" \
        "KRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLP" \
        "DDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGF" \
        "QPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRD" \
        "IADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTG" \
        "SNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNN" \
        "SIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEV" \
        "FAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICA" \
        "QKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLI" \
        "ANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQID" \
        "RLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHV" \
        "TYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNT" \
        "VYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYE" \
        "QYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"

    return s[start - 1 : end]

if __name__ == "__main__":
    # rbd = get_s_substring(319, 541)
    # print(rbd)
    # spike = get_s_substring(1, 1275)
    # count_aas(spike)
    part_rbd = get_s_substring(443, 451)
    print("\n" in part_rbd)
    print(part_rbd)
