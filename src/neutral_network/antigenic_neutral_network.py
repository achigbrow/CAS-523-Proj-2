"""
Python code for building an antigenic neutral network
"""
import copy
import random
import numpy as np
import pandas as pd
from math import cos, sin, pi
import networkx as nx
from bindingcalculator import *
import plotly.graph_objects as go


class AntigenicNeutralNetwork:

    def __init__(self, size, tolerance, binding_calculator):
        self.size = size
        self.tolerance = tolerance
        self.bc = binding_calculator
        self.sites = list(self.bc.sites)  # Used to track the sites we can mutate

        self.nodes = {}  # Holds the node data
        self.neutral_nodes = 0

        # The distance table for this neutral network
        self.distance_table = None  # Type: numpy array

        # The titer table for this neutral network
        self.titer_table = None  # Type: numpy array

        # This is the neutral network, which is a directed graph
        self.nn = nx.DiGraph()

    def build(self, to_print=False, only_mutate_neutral=True, epistatic=False, test_epistatic=False, parent_escape=0.75):
        """
        Builds the neutral network self.nn to have self.size nodes
        If a mutation has a binding ability above self.tolerance it is considered neutral

        :param test_epistatic: Used to ensure the non-epistatic node has the same escape threshold for comparing
            with an epistatic node.
        :param parent_escape: For epistatic nodes, the epistatic node will have at most this much antibody escape
        without accounting for it being epistatic
        :param epistatic: If the effects of epistatic nodes should be considered.
        :param only_mutate_neutral: Only mutate neutral nodes. Set as false to make a
            more interesting antigenic map
        :param to_print: True/False, whether node coordinates should be generated for printing
        :return:
        """
        # Create the root node of the neutral network where there are no mutations
        node = Node(node_id=0,
                    xy_pos=(0, 0, 0),
                    bd=self.bc,
                    tolerance=self.tolerance,
                    mutations=[],
                    to_print=to_print,
                    epistatic=epistatic,
                    root_epistatic_node=epistatic)
        if epistatic:
            # Create the root node of the neutral network where there is one epistatic mutation
            # Mutate the node until it is no longer neutral
            while node.escape > parent_escape:
                node = node.mutate(child_node_id=0, bd=self.bc, sites=self.sites, root_epistatic_node=True)
            print(f"found non-neutral epistatic node with properties:\n{node}")
            print(f"escape adjustment: {node.escape_adjustment}")
        else:
            print(f"found neutral node with properties:\n{node}")
            print(f"escape adjustment: {node.escape_adjustment}")
        self.nodes[0] = node
        self.nn.add_node(node.id, pos=(node.xy_pos[0], node.xy_pos[1]))

        # Set counters
        current_size = 1

        while current_size < self.size:

            # Randomly choose a node that is neutral
            rand_node = self.nodes[random.sample(list(self.nn.nodes()), 1)[0]]
            # print(f"Mutating node {rand_node.id} to make new node {current_size}")
            counter = 0
            while not rand_node.is_neutral and only_mutate_neutral:
                print(rand_node)
                rand_node = self.nodes[random.sample(list(self.nn.nodes()), 1)[0]]
                counter += 1
                if counter == 25:
                    exit(1)

            # Create a child node
            child = rand_node.mutate(child_node_id=current_size, bd=self.bc, sites=self.sites)
            if child.is_neutral:
                self.neutral_nodes += 1

            self.nodes[current_size] = child

            # Add the child to the network
            if to_print:
                self.nn.add_node(child.id, pos=(child.xy_pos[0], child.xy_pos[1]))
                self.nn.add_edge(rand_node.id, child.id)

            current_size += 1

        print(f"Total nodes: {self.size}\nNeutral nodes: {self.neutral_nodes}")

    def make_titer_and_distance_tables(self):
        """
        Tested with tableDistances(map, optimization_number = 1)

        This function computes the distances table for the neutral network and then converts it to a "fake"
        titer table because the titer table cannot be calculated directly
        The formula used is:
        distance_(i, j) = escape_i / escape_j
        titer_(i, j) = 10 * 2^(distance_(i, j))
        This formula was derived from the titer table to log titer table formula at
        https://acorg.github.io/Racmacs/articles/intro-to-antigenic-cartography.html
        but without including a column base.

        This conversion is necessary because the Racmacs program to generate a antigenic map only takes
        titer tables as input

        :return: populates self.distance_table and self.titer_table
        """
        # Check if a neutral network has been built
        node_ids = self.nodes.keys()
        if not len(node_ids):
            return

        # Make the distance table
        self.distance_table = np.ndarray(shape=(self.size, self.size))
        for i in range(self.size):
            for j in range(self.size):
                value = round(abs(1 - self.nodes[i].escape / self.nodes[j].escape), 2)

                # Avoid divide-by-zero errors
                if value < 0.09:
                    self.distance_table[i, j] = 0.09
                else:
                    self.distance_table[i, j] = value

        # Make the titer table
        column_base = 10
        self.titer_table = np.ndarray(shape=(self.size, self.size))
        for i in range(self.size):
            for j in range(self.size):
                self.titer_table[i, j] = round(pow(2, column_base - self.distance_table[i, j]), 0)

    def save_titer_table(self, path):
        # Convert titer table to a pandas dataframe and get it in the format Racmacs wants
        headers = range(self.size)
        titer_df = pd.DataFrame(self.titer_table, columns=headers, index=headers)

        titer_df.to_csv(path)
        print(f"Saved titer.csv to {path}")

    def create_figure(self):
        """
        Code copied from: https://plotly.com/python/network-graphs/
        :return:
        """
        G = self.nn
        # G = nx.random_geometric_graph(200, 0.125)
        edge_x = []
        edge_y = []
        for edge in G.edges():
            x0, y0 = G.nodes[edge[0]]['pos']
            x1, y1 = G.nodes[edge[1]]['pos']
            edge_x.append(x0)
            edge_x.append(x1)
            edge_x.append(None)
            edge_y.append(y0)
            edge_y.append(y1)
            edge_y.append(None)

        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5, color='#888'),
            hoverinfo='none',
            mode='lines')

        node_x = []
        node_y = []
        for node in G.nodes():
            x, y = G.nodes[node]['pos']
            node_x.append(x)
            node_y.append(y)

        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers',
            hoverinfo='text',
            marker=dict(
                showscale=True,
                # colorscale options
                # 'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
                # 'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
                # 'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
                colorscale='YlGnBu',
                reversescale=True,
                color=['#ff0000'],
                size=10,
                cmin=0,
                cmax=1,
                colorbar=dict(
                    thickness=15,
                    title='Node Connections',
                    xanchor='left',
                    titleside='right'
                ),
                line_width=2))

        node_neutral = []
        node_text = []
        for node, adjacencies in enumerate(G.adjacency()):
            val = 1 if self.nodes[node].is_neutral else 0
            node_neutral.append(val)
            node_text.append(f'id: {self.nodes[node].id} neutral: {self.nodes[node].is_neutral}')

        node_trace.marker.color = node_neutral
        node_trace.text = node_text

        fig = go.Figure(data=[edge_trace, node_trace],
                        layout=go.Layout(
                            title='<br>Network graph made with Python',
                            titlefont_size=16,
                            showlegend=False,
                            hovermode='closest',
                            margin=dict(b=20, l=5, r=5, t=40),
                            annotations=[dict(
                                text="Python code: <a href='https://plotly.com/ipython-notebooks/network-graphs/'> https://plotly.com/ipython-notebooks/network-graphs/</a>",
                                showarrow=False,
                                xref="paper", yref="paper",
                                x=0.005, y=-0.002)],
                            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                        )
        fig.show()


def get_epistatic_statistics(y, x, tolerance, parent_escape, binding_calc):
    """
    Calculates the chances of antibody escape from genomes with epistatic interactions,
    according to the following logic:

    There is a 0.001% chance that a non-neutral mutation will be considered neutral because
    it is "epistatic and allows for mutations in antigenically active sites to persist".
    For y of these genomes, we run "x" more mutations on it and calculate how many of
    these x mutations are not neutral (i.e. they escape antibodies better than the parent)
    with respect to the parent genome.
    :param: parent_escape: The epistatic node will have at most this much antibody escape
        without accounting for it being epistatic
    :param: empty_ann: an AntigenicNeutralNetwork object to be used to calculate this question
    :param: y, x: The number of non-neutral mutations to be considered neutral, y, (i.e., epistatic),
        each to be mutated x times

    :return: The probability of antibody escape from genomes with epistatic interations
    """
    # Lists of neutral_nodes numbers
    results_epistatic = []
    results_normal = []
    for i in range(y):
        ann_ep = AntigenicNeutralNetwork(size=x, tolerance=tolerance, binding_calculator=binding_calc)
        ann_norm = AntigenicNeutralNetwork(size=x, tolerance=tolerance, binding_calculator=binding_calc)
        ann_ep.build(epistatic=True, parent_escape=parent_escape)
        #print("finished ann ep")
        ann_norm.build(epistatic=False, test_epistatic=True)
        #print("finished ann norm")
        results_epistatic.append(ann_ep.neutral_nodes)
        results_normal.append(ann_norm.neutral_nodes)
        print(f"Completed for epistatic node {i+1}/{y}...")

    # Calculate the neutral-to-total-size ratios
    results_epistatic = [item / x for item in results_epistatic]
    results_normal = [item / x for item in results_normal]

    avg_epi = sum(results_epistatic) / len(results_epistatic)
    avg_norm = sum(results_normal) / len(results_normal)

    result = avg_epi - avg_norm

    print(f"If {y} epistatic genomes are each randomly mutated {x} more times\nthere will be {result * 100:0.2f}% more"
          f" antigenically-neutral genomes than\nif {y} neutral genomes are each mutated {x} more times")

    return result


class Node:
    """
    A node in the neutral network
    """

    def __init__(self, node_id, xy_pos, bd, tolerance, mutations, parent_escape_adjustment=0, to_print=False,
                 epistatic=False, root_epistatic_node=False):
        """
        :param xy_pos: tuple(xy coords of parent, rotation)
        :param mutations: a list of mutations performed on the mode.
            Sites must be between 331 and 531.
        """
        self.id = node_id
        self.mutations = mutations
        self.child_mutations = []  # Used to ensure no mutations is made twice
        self.tolerance = tolerance
        self.epistatic = epistatic
        self.escape_adjustment = parent_escape_adjustment  # Will be updated if the node is epistatic, default 0

        self.escape = self.get_escape_remaining(bd)
        #print(f"node {self.id} escape:{self.escape + self.escape_adjustment}")

        # Epistatic nodes have their escape adjusted back to 1 after mutation.
        # Used for the get_epistatic_statistics function
        if self.epistatic:
            if root_epistatic_node:
                self.escape_adjustment = 1 - self.escape

        self.xy_pos = (0, 0, 0)
        self.to_print = to_print

        # Calculate if the node is neutral
        self.is_neutral = True if self.escape + self.escape_adjustment > tolerance else False

        # If the network is to be printed, set the xy position of the node depending on whether it is neutral
        if self.to_print:
            self.generate_coords(xy_pos, mutations)

    def generate_coords(self, xy_pos, mutations):
        if len(mutations):
            neutral_radius = random.randint(15, 25)
            mut_radius = 5

            if self.is_neutral:
                radius = neutral_radius
            else:
                radius = mut_radius

            x, y, z = xy_pos
            dx = radius * cos(pi * xy_pos[2])
            dy = radius * sin(pi * xy_pos[2])
            self.xy_pos = (x + dx, y + dy, z)

        else:
            self.xy_pos = xy_pos

    def get_escape_remaining(self, bd):
        """
        :param bd: the binding calculator
        :return: the amount of escape remaining, as a decimal between 0 and 1
        """
        return round(bd.binding_retained(self.mutations), 3)

    def mutate(self, child_node_id, bd, sites, root_epistatic_node=False):
        """
        Adds a random mutation to the self node and returns a new node
        :return: a child node
        """
        mutations = copy.deepcopy(self.mutations)

        # Do not mutate the same site twice
        site = random.choice(sites)
        while site in mutations and site not in self.child_mutations:
            site = random.choice(sites)
        mutations.append(site)
        self.child_mutations.append(site)

        # Set the position for the child node if printing
        if self.to_print:
            rotation_inc = random.random() - 0.5
            x, y, z = self.xy_pos
            z += rotation_inc
            self.xy_pos = (x, y, z)
        child_node = Node(bd=bd,
                          xy_pos=self.xy_pos,
                          node_id=child_node_id,
                          tolerance=self.tolerance,
                          mutations=mutations,
                          to_print=self.to_print,
                          epistatic=self.epistatic,
                          parent_escape_adjustment=self.escape_adjustment,
                          root_epistatic_node=root_epistatic_node)

        return child_node

    def __repr__(self):
        string = f"NODE {self.id}:\n    mutations: {len(self.mutations)}\n"
        string += f"    escape: {self.escape}\n    neutral: {self.is_neutral}\n"
        string += f"    escape adjustment: {self.escape_adjustment}\n"
        return string
