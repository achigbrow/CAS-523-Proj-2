"""
Driver file to create an antigenically neutral network
"""
from antigenic_neutral_network import AntigenicNeutralNetwork


def main():
    # Neutral network params
    tolerance = 0.95
    size = 25

    # Init the antigenically neutral network
    ann = AntigenicNeutralNetwork(size=size, tolerance=tolerance)

    ann.build(to_print=True)
    # for node in ann.nodes.keys():
    #     print(ann.nodes[node])

    # ann.create_figure()


if __name__ == "__main__":
    main()
