import numpy as np


class Fine2Coarse:
    def __init__(self, adjacency_matrix):
        self.adjacency_matrix = adjacency_matrix
        self.number_of_nodes = self.adjacency_matrix.shape[0]

    def beck(self):
        """
        We sequentially iterate over all nodes. If selected node
        is not in coarse set or fine set, we append the node to
        the coarse set.
        Then we append all the unlabled connected nodes to fine set.
        """
        fine = []
        coarse = []
        for node in range(self.number_of_nodes):
            if node not in fine and node not in coarse:
                # Append to coarse set
                coarse.append(node)
                # Find the connected nodes to this node
                # Connected nodes in row of adj. matrix
                connected_nodes_in_row = set(np.argwhere(
                    self.adjacency_matrix[node] != 0)[:, 1]) - set([node])
                # Connected nodes in column of adj. matrix
                connected_nodes_in_column = set(np.argwhere(
                    self.adjacency_matrix[:, node] != 0)[:, 0]) - set([node])
                # All connected nodes
                connected_nodes = connected_nodes_in_column | connected_nodes_in_row
                # We assign all the unlabled connected nodes to fine set
                for neighbor_node in connected_nodes:
                    if neighbor_node not in fine and neighbor_node not in coarse:
                        fine.append(neighbor_node)
        return fine, coarse
