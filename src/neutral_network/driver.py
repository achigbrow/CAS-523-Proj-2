"""
Driver file to create an antigenically neutral network
"""
from antigenic_neutral_network import AntigenicNeutralNetwork, get_epistatic_statistics
from src.neutral_network.bindingcalculator import BindingCalculator


def main():
    # Neutral network params
    tolerance = 0.95
    size = 200
    print_network = True
    only_mutate_neutral = True
    save_titer_table = True

    # Params for epistatic calcs
    # If consider_epistatic_change is True, only epistatic change will be considered and the above parameters will
    # be ignored.
    consider_epistatic_change = False  # If set to False, the following values will not be used.
    starting_nodes = 10
    times_to_mutate = 300
    parent_escape = 0.6  # Between 0 and 1

    #################################################
    # Nothing needs to be edited below this comment #
    #################################################
    # Init the antigenically neutral network
    print("Initializing binding calculator...", end="")
    bc = BindingCalculator()
    print("Done")
    ann = AntigenicNeutralNetwork(size=size, tolerance=tolerance, binding_calculator=bc)

    if consider_epistatic_change:
        get_epistatic_statistics(y=starting_nodes,
                                 x=times_to_mutate,
                                 tolerance=tolerance,
                                 binding_calc=bc,
                                 parent_escape=parent_escape)
        return

    else:
        ann.build(to_print=True, only_mutate_neutral=only_mutate_neutral)

    if save_titer_table:
        ann.make_titer_and_distance_tables()
        path = "C:/Users/jason/AppData/Local/R/win-library/4.2/Racmacs/extdata/titer.csv"
        ann.save_titer_table(path)

    if print_network:
        ann.create_figure()

    # for node in ann.nodes.keys():
    #     print(ann.nodes[node])


if __name__ == "__main__":
    main()
