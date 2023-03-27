"""
Driver file to create an antigenically neutral network
"""
from antigenic_neutral_network import AntigenicNeutralNetwork


def main():
    # Neutral network params
    tolerance = 0.95
    size = 2000

    # Init the antigenically neutral network
    ann = AntigenicNeutralNetwork(size=size, tolerance=tolerance)

    ann.build(to_print=False, only_mutate_neutral=False)
    ann.make_titer_and_distance_tables()
    path = "C:/Users/jason/AppData/Local/R/win-library/4.2/Racmacs/extdata/titer.csv"
    ann.save_titer_table(path)

    #ann.create_figure()

    # for node in ann.nodes.keys():
    #     print(ann.nodes[node])


if __name__ == "__main__":
    main()
