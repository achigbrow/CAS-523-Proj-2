"""
Python code for building an antigenic neutral network
"""
import random

import networkx as nx
from bindingcalculator import *
import plotly.graph_objects as go


class AntigenicNeutralNetwork:

    def __init__(self, size, tolerance):
        self.size = size
        self.tolerance = tolerance
        self.bd = BindingCalculator()
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
        node = Node(node_id=0, bd=self.bd, tolerance=self.tolerance, mutations=[])
        self.nodes[0] = node
        self.nn.add_node(node.id, pos=(0, 0))
        current_size = 1
        while current_size < self.size:

            # Randomly choose a node that is neutral
            index = random.sample(list(self.nn.nodes()), 1)[0]
            rand_node = self.nodes[random.sample(list(self.nn.nodes()), 1)[0]]
            print(f"rand node id: {rand_node.id}")
            while not rand_node.is_neutral:
                rand_node = self.nodes[random.sample(list(self.nn.nodes()), 1)[0]]

            # Create a child node
            child = rand_node.mutate(child_node_id=current_size, bd=self.bd, sites=self.sites)
            self.nodes[current_size] = child

            # Add the child to the network
            self.nn.add_node(child.id, pos=(child.id % 3, child.id))
            self.nn.add_edge(rand_node.id, child.id)

            print(f"made node {current_size}")
            current_size += 1

    def create_figure(self):
        """
        Code copied from: https://plotly.com/python/network-graphs/
        :return:
        """
        G = self.nn
        #G = nx.random_geometric_graph(200, 0.125)
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
                colorbar=dict(
                    thickness=15,
                    title='Node Connections',
                    xanchor='left',
                    titleside='right'
                ),
                line_width=2))

        node_adjacencies = []
        node_text = []
        for node, adjacencies in enumerate(G.adjacency()):
            node_adjacencies.append(len(adjacencies[1]))
            node_text.append('# of connections: ' + str(len(adjacencies[1])))

        node_trace.marker.color = node_adjacencies
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

    def __init__(self, node_id, bd, tolerance, mutations):
        """
        :param mutations: a list of mutations performed on the mode.
            Sites must be between 331 and 531.
        """
        self.id = node_id
        self.mutations = mutations
        self.child_mutations = []  # Used to ensure no mutations is made twice
        self.tolerance = tolerance
        self.escape = self.get_escape_remaining(bd)
        self.is_neutral = True if self.escape > tolerance else False

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

        child_node = Node(bd=bd, node_id=child_node_id, tolerance=self.tolerance, mutations=mutations)

        return child_node

    def __repr__(self):
        string = f"NODE:\nmutations: {len(self.mutations)}\n"
        string += f"escape: {self.escape}\nneutral: {self.is_neutral}"
        return string
