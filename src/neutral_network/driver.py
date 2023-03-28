"""
Driver file to create an antigenically neutral network
"""
from antigenic_neutral_network import AntigenicNeutralNetwork, get_epistatic_statistics
from src.neutral_network.bindingcalculator import BindingCalculator


def main():
    # Neutral network params
    tolerance = 0.95
    size = 10
    print_network = False
    only_mutate_neutral = True
    consider_epistatic_change = True

    # Params for epistatic calcs
    num_epistatic_mutations = 10
    times_to_mutate = 300
    parent_escape = 0.6

    save_titer_table = False

    # Init the antigenically neutral network
    print("Initializing binding calculator...", end="")
    bc = BindingCalculator()
    print("Done")
    ann = AntigenicNeutralNetwork(size=size, tolerance=tolerance, binding_calculator=bc)

    if consider_epistatic_change:
        get_epistatic_statistics(y=num_epistatic_mutations,
                                 x=times_to_mutate,
                                 tolerance=tolerance,
                                 binding_calc=bc,
                                 parent_escape=parent_escape)
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
