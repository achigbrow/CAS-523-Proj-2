"""
Python code for building an antigenic neutral network
"""
import random

import networkx as nx
from bindingcalculator import *


def main():
    # Global params
    tolerance = 0.8
    mutations = []  # Root node has no mutations

    # Make a new escape calculator
    bd = BindingCalculator()

    # Create the neutral network
    nn = nx.DiGraph()

    # Create the root node of the neutral network where there are no mutations
    node = Node(bd, mutations, tolerance)
    nn.add_node(node)

    print(nn)


def get_escape_remaining(bd, mutations):
    """
    :param bd: the binding calculator
    :param mutations: a list of sites on the spike protein to mutate. Sites must be between 331 and 531.
    :return: the amount of escape remaining, as a decimal between 0 and 1
    """
    return round(bd.binding_retained(mutations), 3)


class Node:
    """
    A node in the neutral network
    """

    def __init__(self, bd, mutations, tolerance):
        """
        :param mutations: a list of mutations performed on the mode.
            Sites must be between 331 and 531.
        """
        self.mutations = mutations
        self.escape = get_escape_remaining(bd, mutations)
        self.tolerance = tolerance
        self.is_neutral = True if self.escape > tolerance else False

    def create_child(self, bd):
        """
        Adds a random mutation to the self node and returns a new node
        :return: a child node
        """
        mutations = self.mutations

        site = 0

        # Do not mutate the same site twice
        while site not in mutations:
            site = random.randrange(331, 531)
        mutations.append(site)

        child_node = Node(bd, mutations, self.tolerance)

        return child_node


if __name__ == "__main__":
    main()
