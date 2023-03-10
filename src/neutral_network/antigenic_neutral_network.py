"""
Python code for building an antigenic neutral network
"""
import random
from math import cos, sin, pi, dist
import networkx as nx
from bindingcalculator import *
import plotly.graph_objects as go


class AntigenicNeutralNetwork:

    def __init__(self, size, tolerance):
        self.size = size
        self.tolerance = tolerance
        print("Initializing binding calculator...", end="")
        self.bd = BindingCalculator()
        print("Done")
        self.sites = list(self.bd.sites)  # Used to track the sites we can mutate

        self.nodes = {}  # Holds the node data

        # This is the neutral network, which is a directed graph
        self.nn = nx.DiGraph()

    def build(self):
        """
        Builds the neutral network self.nn to have self.size nodes
        If a mutation has a binding ability above self.tolerance it is considered neutral

        :return:
        """

        # Create the root node of the neutral network where there are no mutations
        node = Node(node_id=0, xy_pos=(0, 0, 0), bd=self.bd, tolerance=self.tolerance, mutations=[])
        self.nodes[0] = node
        self.nn.add_node(node.id, pos=(node.xy_pos[0], node.xy_pos[1]))
        current_size = 1
        neutral_nodes = 1
        while current_size < self.size:

            # Randomly choose a node that is neutral
            index = random.sample(list(self.nn.nodes()), 1)[0]
            rand_node = self.nodes[random.sample(list(self.nn.nodes()), 1)[0]]
            print(f"Mutating node {rand_node.id} to make new node {current_size}")
            while not rand_node.is_neutral:
                rand_node = self.nodes[random.sample(list(self.nn.nodes()), 1)[0]]

            # Create a child node
            child = rand_node.mutate(child_node_id=current_size, bd=self.bd, sites=self.sites)
            if child.is_neutral:
                neutral_nodes += 1
            self.nodes[current_size] = child

            # Add the child to the network
            self.nn.add_node(child.id, pos=(child.xy_pos[0], child.xy_pos[1]))
            self.nn.add_edge(rand_node.id, child.id)

            current_size += 1

        print(f"Number of neutral nodes: {neutral_nodes}\nTotal nodes: {self.size}")

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
            node_text.append(f'neutral: {self.nodes[node].is_neutral}')

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


class Node:
    """
    A node in the neutral network
    """

    def __init__(self, node_id, xy_pos, bd, tolerance, mutations):
        """
        :param: xy_pos: tuple(xy coords of parent, rotation)
        :param mutations: a list of mutations performed on the mode.
            Sites must be between 331 and 531.
        """
        self.id = node_id
        self.mutations = mutations
        self.child_mutations = []  # Used to ensure no mutations is made twice
        self.tolerance = tolerance
        self.escape = self.get_escape_remaining(bd)
        self.is_neutral = True if self.escape > tolerance else False

        # Set the xy position of the node depending on whether it is neutral
        if len(mutations):
            neutral_radius = 25
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

    def mutate(self, child_node_id, bd, sites):
        """
        Adds a random mutation to the self node and returns a new node
        :return: a child node
        """
        mutations = self.mutations

        # Do not mutate the same site twice
        site = random.choice(sites)
        while site in mutations and site not in self.child_mutations:
            site = random.choice(sites)
        mutations.append(site)
        self.child_mutations.append(site)

        # Set the position for the child node
        rotation_inc = 0.1
        x, y, z = self.xy_pos
        z += rotation_inc
        self.xy_pos = (x, y, z)
        child_node = Node(bd=bd, xy_pos=self.xy_pos, node_id=child_node_id, tolerance=self.tolerance,
                          mutations=mutations)

        return child_node

    def __repr__(self):
        string = f"NODE:\nmutations: {len(self.mutations)}\n"
        string += f"escape: {self.escape}\nneutral: {self.is_neutral}"
        return string
