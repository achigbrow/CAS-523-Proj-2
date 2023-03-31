import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from rbd import *


def plot_aas():
    aas = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
           "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    cols = ["Source", "Amino Acid", "Frequency"]

    df = pd.DataFrame(columns=cols, index=range(0, 100))

    # get RBD and spike strings
    rbd = get_s_substring(319, 541)
    spike = get_s_substring(1, 1274)
    omicron_spike = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFF" \
                    "SNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLL" \
                    "IVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLM" \
                    "DLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINIT" \
                    "RFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDP" \
                    "LSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWN" \
                    "RKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPG" \
                    "QTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEI" \
                    "YQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKST" \
                    "NLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPC" \
                    "SFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGC" \
                    "LIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYS" \
                    "NNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRAL" \
                    "TGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKV" \
                    "TLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSG" \
                    "WTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASA" \
                    "LGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQS" \
                    "LQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVF" \
                    "LHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTF" \
                    "VSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQ" \
                    "KEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTS" \
                    "CCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"

    # get RBD and spike amino acid frequencies
    _, rfreqs = count_aas(rbd)
    _, sfreqs = count_aas(spike)
    _, osfreqs = count_aas(omicron_spike)

    rbd_freqs = list(rfreqs.values())
    spike_freqs = list(sfreqs.values())
    os_freqs = list(osfreqs.values())

    # Hardcoded from paper referenced in README
    human_host_freqs = [0.05054, 0.04901, 0.05869, 0.05532, 0.02181, 0.03646,
                        0.0565, 0.05148, 0.02222, 0.07077, 0.09374, 0.06532,
                        0.02208, 0.04273, 0.04821, 0.07559, 0.0667, 0.01248,
                        0.04059, 0.05976]
    vertebrate_freqs = [0.0681465, 0.0583065, 0.0482765, 0.0563865, 0.0208265,
                        0.0320565, 0.0560365, 0.0583265, 0.0217265, 0.0592565,
                        0.0904965, 0.0565965, 0.0251565, 0.0426065, 0.0521865,
                        0.0744765, 0.0601965, 0.0119865, 0.0390565, 0.0678965]

    sources = ["Human-Host Viral Proteome", "Wuhan-Hu-1 Spike Protein",
               "Wuhan-Hu-1 RBD", "Vertebrate-Host Viral Proteome",
               "Omicron Spike Protein"]
    all_freqs = [human_host_freqs, spike_freqs, rbd_freqs, vertebrate_freqs,
                 os_freqs]
    counter = 0

    for i in range(0, 100):
        temp = [sources[int(i / 20)], aas[counter],
                all_freqs[int(i / 20)][counter]]

        counter += 1
        if counter > 19:
            counter = 0

        df.iloc[i] = temp

    # plot data
    sns.set_theme(style="whitegrid")


    # g = sns.catplot(
    #     data=df, kind="bar",
    #     x="Amino Acid", y="Frequency", hue="Source",
    #     errorbar="sd", palette="dark", alpha=.6, height=6
    # )
    # g.despine(left=True)
    # g.set_axis_labels("", "Frequency")
    # g.legend.set_title("Frequency Source")
    # plt.show()

    # no_spike = df[df["Source"] != "Wuhan-Hu-1 Spike Protein"]
    #
    # g1 = sns.catplot(
    #     data=no_spike, kind="bar",
    #     x="Amino Acid", y="Frequency", hue="Source",
    #     errorbar="sd", palette="dark", alpha=.6, height=6
    # )
    # g1.despine(left=True)
    # g1.set_axis_labels("", "Frequency")
    # g1.legend.set_title("Frequency Source")
    # plt.show()
    #
    no_rbd = df[df["Source"] != "Wuhan-Hu-1 RBD"]

    g = sns.lineplot(data=no_rbd,
                     x="Amino Acid", y="Frequency", hue="Source",
                     style="Source"
                     )
    # g2 = sns.catplot(
    #     data=no_rbd, kind="bar",
    #     x="Amino Acid", y="Frequency", hue="Source",
    #     errorbar="sd", palette="dark", alpha=.6, height=6
    # )
    # g2.despine(left=True)
    # g2.set_axis_labels("", "Frequency")
    # g2.legend.set_title("Frequency Source")
    plt.show()

def plot_probs():

    df = pd.DataFrame()

    df["Amino Acid"] = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
                        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    # Hardcoded from paper referenced in README
    df["Probability"] = [0.333, 0.333, 0.111, 0.111, 0.111, 0.111, 0.111, 0.333, 0.111,
                         0.222, 0.333, 0.111, 0, 0.111, 0.333, 0.259, 0.333, 0, 0.111,
                         0.333]

    sns.set_theme(style="whitegrid")

    g = sns.catplot(
        data=df, kind="bar",
        x="Amino Acid", y="Probability",
        errorbar="sd", palette="dark", alpha=.6, height=6
    )
    g.despine(left=True)
    g.set_axis_labels("", "Probability of Self-Preservation")
    plt.show()



if __name__ == "__main__":
    # plot_probs()
    plot_aas()
