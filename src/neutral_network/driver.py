"""
Driver file to create an antigenically neutral network
"""
from antigenic_neutral_network import AntigenicNeutralNetwork


def main():
    # Neutral network params
    tolerance = 0.7
    size = 540

    # Init the antigenically neutral network
    ann = AntigenicNeutralNetwork(size=size, tolerance=tolerance)

    ann.build()

    ann.create_figure()


if __name__ == "__main__":
    main()
