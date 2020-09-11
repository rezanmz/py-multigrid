class Prolongator:
    def __init__(self, adjacency_matrix, fine, coarse):
        self.adjacency_matrix = adjacency_matrix
        self.fine = fine
        self.coarse = coarse
        self.number_of_nodes = len(fine) + len(coarse)

    def beck(self):
        from scipy import sparse
        import numpy as np

        p = sparse.lil_matrix((self.number_of_nodes, len(self.coarse)))
        for node in self.coarse:
            # Identity for coarse nodes in prolongation matrix
            p[node, np.argwhere(np.array(self.coarse) == node).reshape(-1)] = 1

        for node in self.fine:
            # Find all connected nodes to the current fine node
            # Connected nodes in row of adj. matrix
            connected_nodes_in_row = set(np.argwhere(
                self.adjacency_matrix[node] != 0)[:, 1]) - set([node])
            # Connected nodes in column of adj. matrix
            connected_nodes_in_column = set(np.argwhere(
                self.adjacency_matrix[:, node] != 0)[:, 0]) - set([node])
            # All connected nodes
            connected_nodes = connected_nodes_in_column | connected_nodes_in_row
            # From all the connected nodes, find coarse nodes
            connected_coarse = connected_nodes & set(self.coarse)
            # Calculate sigma (number of coarse nodes in all connected nodes)
            sigma = len(connected_coarse)
            # Find index of each connected coarse node in the list of all coarse nodes
            connected_coarse = list(map(
                lambda con_coarse: True if con_coarse in connected_coarse else False, self.coarse))
            # According to back algorithm, fill the prolongator matrix
            p[node, connected_coarse] = 1 / sigma
        return p.tocsr()
