import numpy as np
from copy import deepcopy
import networkx as nx
import rdkit
from rdkit import Chem
# def find_fitting_lists(lists, adj_matrix):
#     # return the indices of the lists that can be appended together to create a longer list
#     # with no overlapping values in those lists

#     all_values = set()
#     for l in lists:
#         for value in l:
#             all_values.add(value)
#     max_val = max(all_values)
#     min_val = min(all_values)
#     matrix_data = []
#     for l in lists:
#         row = []
#         for i in range(0, max_val+1):
#             if i in l:
#                 row.append(1)
#             else:
#                 row.append(0)
#         matrix_data.append(row)
#     matrix = np.matrix(matrix_data)
#     n_rows, n_cols = matrix.shape
#     # first, find unique rows to start the recursive algorithm
#     unique_row_ids = list(range(0, len(lists)))
#     # we only want one unique_col for each unique_row_ids
#     unique_mapping_groups = []
#     while len(unique_row_ids) != 0:
#         unique_row = unique_row_ids[0]
#         row = np.squeeze(np.asarray(matrix[unique_row, :]))
#         search_scope, = np.where(row==1)
#         list_groups = reduce_matrix(matrix, adj_matrix, search_scope.tolist())
#         print(list_groups)
#         if list_groups: # if list group found
#             # get the largest list group with the fewest number of rows
#             for list_group in list_groups:
#                 matrix[list_group, :] = np.zeros((len(list_group), n_cols))
#                 group_length = matrix[list_group, :].sum()
#                 unique_mapping_groups.append(tuple([group_length, list_group]))
#         else: # returns False if there is a violation
#             continue
#         # remove the found unique mapping from unique_row_ids
#         found_row_ids = []
#         for list_group in list_groups:
#             found_row_ids += list_group
#         found_row_ids = list(set(found_row_ids))
#         new_unique_row_ids = deepcopy(unique_row_ids)
#         for row_id in unique_row_ids:
#             if row_id in found_row_ids:
#                 new_unique_row_ids.remove(row_id)
#         unique_row_ids = new_unique_row_ids
#     # finally pick from unique_mapping_groups the largest list
#     largest_group_id = 0
#     largest_group_length = 0
#     last_group_list_size = 0
#     i = 0
#     for group_length, mapping_group in unique_mapping_groups:
#         if group_length >= largest_group_length and len(mapping_group) < last_group_list_size:
#             largest_group_id = i
#             largest_group_length = group_length
#             last_group_list_size = len(mapping_group)
#         i += 1
#     l, largest_group = unique_mapping_groups[largest_group_id]
#     return largest_group
# bla = 1
# def reduce_matrix(matrix, adj_matrix, search_scope=[]):
#     def get_adjacencies(adj_matrix, col_ids):
#         # returns the column ids of adjacent nodes to the given col_ids
#         all_connections = np.asarray(adj_matrix[col_ids, :].sum(axis=0)).reshape(-1)
#         all_connections[col_ids] = np.zeros(len(col_ids))
#         adjacent_connections, = np.nonzero(all_connections)
#         return adjacent_connections.tolist()
#     # reduces the matrix and returns the rows that were cut plus a smaller matrix for recursion

#     n_rows, n_cols = matrix.shape

#     if search_scope == []:
#         search_scope = list(range(0, n_cols))
#     new_search_scope = []
#     for idx in search_scope:
#         if idx < n_cols:
#             new_search_scope.append(idx)
#     search_scope = new_search_scope
            
#     sums = np.asarray(matrix.sum(axis=0)).reshape(-1)
#     print(sums.sum())
#     # required_unique_cols, = np.where(sums==1)
#     # required_unique_row_ids, = np.nonzero(np.asarray(matrix[:, required_unique_cols].sum(axis=1)).reshape(-1))
#     # required_unique_row_ids = set(required_unique_row_ids.tolist())
#     allowed_sums = sums[search_scope]
#     allowed_sums_ids, = np.where(allowed_sums > 0)
#     allowed_nonzero_sums = allowed_sums[allowed_sums_ids]
#     if allowed_nonzero_sums.size == 0:
#         return [[]]
#     min_allowed_value = min(allowed_nonzero_sums)
#     if min_allowed_value > 1:
#         test = 1
#     min_ids, = np.where(allowed_sums==min_allowed_value)
#     last_min_value_id = search_scope[min_ids[0]]
    
#     # identify columns that have only single elements, and find the unique rows
#     fcolumn = np.squeeze(np.asarray(matrix[:, last_min_value_id]))
#     unique_row_lists = []

#     for value, frow_id in zip(fcolumn, range(0, n_rows)):
#         if value == 1:
#             # first, eliminate all rows with overlap to the current row
#             # then call reduce_matrix recursively for each occurance of 1
#             chosen_cols = np.asarray(matrix[frow_id, :].sum(axis=0)).reshape(-1).astype(bool)
#             chosen_col_ids, = np.where(chosen_cols==True)
#             new_search_scope = get_adjacencies(adj_matrix, chosen_col_ids) + search_scope
#             overlap_col_ids = np.array(range(0,len(sums)))[np.logical_and(sums > 1, chosen_cols)]
#             # find rows that have overlapping values (and therefore overlapping substructures)
#             overlap_rows = set()
#             for col_id in overlap_col_ids:
#                 column = np.squeeze(np.asarray(matrix[:, col_id]))
#                 for row_val, row_id in zip(column, range(0, n_rows)):
#                     # if row_val == 1 and row_id not in unique_rows:
#                     #     overlap_rows.add(row_id)
#                     if row_val == 1:
#                         overlap_rows.add(row_id)
#             # if the overlap_rows (which are discarded) share any rows with the required_unique_row_ids
#             # then, this path of adjacent isomorphisms is invalid 
#             # if len(overlap_rows & (required_unique_row_ids - {frow_id})) > 0:
#             #     return False
    
#             cut_rows = list(set([frow_id]) | overlap_rows)
#             matrix_copy = deepcopy(matrix)
#             matrix_copy[list(cut_rows), :] = np.zeros((len(cut_rows), n_cols))
#             new_unique_rows = reduce_matrix(matrix_copy, adj_matrix, new_search_scope)
#             if new_unique_rows == False:
#                 return False
#             for unique_row in new_unique_rows:
#                 unique_row_lists.append([frow_id, *unique_row])
#             break
#     if len(unique_row_lists) > 15:
#         test = 1
#     return unique_row_lists

# def find_fitting_lists(isomorphism_info):
#     # return the indices of the lists that can be appended together to create a longer list
#     # with no overlapping values in those lists

#     all_values = set()
#     for isomorphism in isomorphism_info:
#         iso_ids, adjacencies = isomorphism
#         for value in iso_ids:
#             all_values.add(value)
#     max_val = max(all_values)
#     min_val = min(all_values)
#     # create a matrix to find atom overlaps
#     matrix_data = []
#     adjacency_list = []
#     for isomorphism in isomorphism_info:
#         iso_ids, adjacencies = isomorphism
#         row = []
#         for i in range(0, max_val+1):
#             if i in iso_ids:
#                 row.append(1)
#             else:
#                 row.append(0)
#         matrix_data.append(row)
#         adjacencies_tupled = [] # adjacencies in tuple form for easier searching
#         for bond, bond_type in adjacencies.items():
#             edge_atom, cap_atom = bond
#             adjacencies_tupled.append(tuple([edge_atom, cap_atom, bond_type]))
#         adjacency_list.append(adjacencies_tupled)

#     # matrix_data and adjacency_list are parallel to eachother
#     # each item in adjacency_list corresponds to the isomorphism found in the 
#     # corresponding row of matrix_data
#     matrix = np.matrix(matrix_data)
#     n_rows, n_cols = matrix.shape
#     # first, find unique rows to start the recursive algorithm
#     unique_row_ids = list(range(0, n_rows))
#     # we only want one unique_col for each unique_row_ids
#     mapping_groups = []
#     while len(unique_row_ids) != 0:
#         unique_row = unique_row_ids[0]
#         row = np.squeeze(np.asarray(matrix[unique_row, :]))
#         search_scope, = np.where(row==1)
#         list_groups = reduce_matrix([unique_row], adjacency_list, found_isomorphisms = [], found_edges = [])
#         if list_groups: # if list group found
#             # get the largest list group with the fewest number of rows
#             group_length = matrix[list_groups, :].sum()
#             matrix[list_groups, :] = np.zeros((len(list_groups), n_cols))
#             mapping_groups.append(tuple([group_length, list_groups]))
#         else:
#             continue
#         # remove the found unique mapping from unique_row_ids
#         new_unique_row_ids = deepcopy(unique_row_ids)
#         for row_id in unique_row_ids:
#             if row_id in list_groups:
#                 new_unique_row_ids.remove(row_id)
#         unique_row_ids = new_unique_row_ids
#     # finally pick from unique_mapping_groups the largest list
#     largest_group_id = 0
#     largest_group_length = 0
#     last_group_list_size = 0
#     i = 0
#     for group_length, mapping_group in mapping_groups:
#         if group_length >= largest_group_length and len(mapping_group) < last_group_list_size:
#             largest_group_id = i
#             largest_group_length = group_length
#             last_group_list_size = len(mapping_group)
#         i += 1
#     l, largest_group = mapping_groups[largest_group_id]
#     return largest_group

# def reduce_matrix(queue, adjacency_list, found_isomorphisms = [], found_edges = []):
#     # searches for isomorphisms that fit together based on isomorphism adjacency info: 

#     if queue == []:
#         return found_isomorphisms
#     v = queue.pop(0)
#     found_isomorphisms.append(v)
#     # search v for adjacent isomorphisms
#     adjacency_info = adjacency_list[v]
#     for edge_atom, cap_atom, bond_type in adjacency_info:
#         neighbor_bond = (cap_atom, edge_atom, bond_type) # adjacent isomorphism will contain a bond with reversed indices
#         # attempt to find the first instance of the cap_atom in neighboring isomorphisms
#         for neighbor_bond_info, id in zip(adjacency_list, range(0, len(adjacency_list))):
#             if neighbor_bond in neighbor_bond_info and id not in found_isomorphisms and {edge_atom, cap_atom} not in found_edges:
#                 queue.append(id)
#                 found_edges.append(set([edge_atom, cap_atom]))
#                 break
#     return reduce_matrix(queue, adjacency_list, found_isomorphisms, found_edges)

# def find_fitting_lists(isomorphism_info):
#     # return the indices of the lists that can be appended together to create a longer list
#     # with no overlapping values in those lists
#     def _reduce_matrix(queue, adjacency_list, found_isomorphisms = [], found_edges = [], favor_isomorphisms = []):
#         # searches for isomorphisms that fit together based on isomorphism adjacency info: 

#         if queue == []:
#             return found_isomorphisms
#         v = queue.pop(0)
#         found_isomorphisms.append(v)
#         # search v for adjacent isomorphisms
#         adjacency_info = adjacency_list[v]
#         for edge_atom, cap_atom, bond_type in adjacency_info:
#             neighbor_bond = (cap_atom, edge_atom, bond_type) # adjacent isomorphism will contain a bond with reversed indices
#             # attempt to find the first instance of the cap_atom in neighboring isomorphisms
#             found_iso_id = -1
#             for neighbor_bond_info, id in zip(adjacency_list, range(0, len(adjacency_list))):
#                 if neighbor_bond in neighbor_bond_info and id not in found_isomorphisms and {edge_atom, cap_atom} not in found_edges:
#                     found_iso_id = id
#                     if id in favor_isomorphisms:
#                         break
#             if found_iso_id > 0:
#                 queue.append(found_iso_id)
#                 found_edges.append(set([edge_atom, cap_atom]))
#         return _reduce_matrix(queue, adjacency_list, found_isomorphisms, found_edges)

#     all_values = set()
#     for isomorphism in isomorphism_info:
#         iso_ids, adjacencies = isomorphism
#         for value in iso_ids:
#             all_values.add(value)
#     max_val = max(all_values)
#     min_val = min(all_values)
#     # create a matrix to find atom overlaps
#     matrix_data = []
#     adjacency_list = []
#     for isomorphism in isomorphism_info:
#         iso_ids, adjacencies = isomorphism
#         row = []
#         for i in range(0, max_val+1):
#             if i in iso_ids:
#                 row.append(1)
#             else:
#                 row.append(0)
#         matrix_data.append(row)
#         adjacencies_tupled = [] # adjacencies in tuple form for easier searching
#         for bond, bond_type in adjacencies.items():
#             edge_atom, cap_atom = bond
#             adjacencies_tupled.append(tuple([edge_atom, cap_atom, bond_type]))
#         adjacency_list.append(adjacencies_tupled)

#     # matrix_data and adjacency_list are parallel to eachother
#     # each item in adjacency_list corresponds to the isomorphism found in the 
#     # corresponding row of matrix_data
#     matrix = np.matrix(matrix_data)
#     full_matrix = deepcopy(matrix)
#     n_rows, n_cols = matrix.shape
#     # first, find unique rows to start the recursive algorithm
#     # we only want one unique_col for each unique_row_ids
#     largest_mapping_group = []
#     largest_group_length = 0
#     for unique_row in range(0, n_rows):
#         if unique_row in largest_mapping_group:
#             continue
#         row = np.squeeze(np.asarray(matrix[unique_row, :]))
#         search_scope, = np.where(row==1)
#         list_groups = _reduce_matrix([unique_row], adjacency_list, found_isomorphisms = [], found_edges = [], favor_isomorphisms=largest_mapping_group)
#         if list_groups: # if list group found
#             # get the largest list group with the fewest number of rows
#             group_length = full_matrix[list_groups, :].sum()
#             matrix[list_groups, :] = np.zeros((len(list_groups), n_cols))
#             if (group_length > largest_group_length) or (group_length == largest_group_length and len(list_groups) < len(largest_mapping_group)):
#                 largest_group_length = group_length
#                 largest_mapping_group = list_groups
#         else:
#             continue
#     # finally pick from unique_mapping_groups the largest list
#     return largest_mapping_group

        # def _find_fitting_lists(isomorphism_info):
        #     # return the indices of the lists that can be appended together to create a longer list
        #     # with no overlapping values in those lists
        #     def _reduce_matrix(queue, matrix, adjacency_list, found_isomorphisms = [], found_edges = [], favor_isomorphisms = []):
        #         # searches for isomorphisms that fit together based on isomorphism adjacency info: 
        #         if queue == []:
        #             return found_isomorphisms
        #         v = queue.pop(0)
        #         found_isomorphisms.append(v)
        #         # search v for adjacent isomorphisms
        #         adjacency_info = adjacency_list[v]
        #         for edge_atom, cap_atom, bond_type in adjacency_info:
        #             neighbor_bond = (cap_atom, edge_atom, bond_type) # adjacent isomorphism will contain a bond with reversed indices
        #             # attempt to find the first instance of the cap_atom in neighboring isomorphisms
        #             found_iso_id = -1
        #             biggest_neighbor_size = 0
        #             for neighbor_bond_info, id in zip(adjacency_list, range(0, len(adjacency_list))):
        #                 if neighbor_bond in neighbor_bond_info and id not in (found_isomorphisms + queue) and {edge_atom, cap_atom} not in found_edges:
        #                     # if neighbor_bond in [(32, 31, bond_type), (34,35, bond_type), (39,40, bond_type), (39,44, bond_type)]:
        #                     #     test = 1
        #                     if np.any(np.squeeze(np.asarray(matrix[found_isomorphisms + [id], :].sum(axis=0))) > 1): # if there is an overlap with previous isomorphisms
        #                         continue
        #                     # this does solve the problem of overlapping isomorphisms, but sometimes you
        #                     # actually want isomorphisms to replace previous isomorphisms that didn't work so well 
        #                     # if id in favor_isomorphisms:
        #                     #     found_iso_id = id
        #                     #     break
        #                     if matrix[id, :].sum() > biggest_neighbor_size:
        #                         found_iso_id = id
        #                         biggest_neighbor_size = matrix[id, :].sum()
                            
        #                     #     ####
        #                     #     found_isomorphisms += favor_isomorphisms
        #                     #     found_iso_id = -1
        #                     #     break
        #                     #     #####
        #                     # else:
        #                     #     found_iso_id = id
        #                     # break
                            
        #             if found_iso_id > 0:
        #                 queue.append(found_iso_id)
        #                 found_edges.append(set([edge_atom, cap_atom]))
        #         return _reduce_matrix(queue, matrix, adjacency_list, found_isomorphisms, found_edges)

        #     all_values = set()
        #     for isomorphism in isomorphism_info:
        #         iso_ids, adjacencies = isomorphism
        #         for value in iso_ids:
        #             all_values.add(value)
        #     max_val = max(all_values)
        #     min_val = min(all_values)
        #     # create a matrix to find atom overlaps
        #     matrix_data = []
        #     adjacency_list = []
        #     for isomorphism in isomorphism_info:
        #         iso_ids, adjacencies = isomorphism
        #         row = []
        #         for i in range(0, max_val+1):
        #             if i in iso_ids:
        #                 row.append(1)
        #             else:
        #                 row.append(0)
        #         matrix_data.append(row)
        #         adjacencies_tupled = [] # adjacencies in tuple form for easier searching
        #         for bond, bond_type in adjacencies.items():
        #             edge_atom, cap_atom = bond
        #             adjacencies_tupled.append(tuple([edge_atom, cap_atom, bond_type]))
        #         adjacency_list.append(adjacencies_tupled)

        #     # matrix_data and adjacency_list are parallel to eachother
        #     # each item in adjacency_list corresponds to the isomorphism found in the 
        #     # corresponding row of matrix_data
        #     matrix = np.matrix(matrix_data)
        #     full_matrix = deepcopy(matrix)
        #     n_rows, n_cols = matrix.shape
        #     # first, find unique rows to start the recursive algorithm
        #     # we only want one unique_col for each unique_row_ids
        #     largest_mapping_group = []
        #     largest_group_length = 0
        #     for unique_row in range(0, n_rows):
        #         if unique_row in largest_mapping_group:
        #             continue
        #         row = np.squeeze(np.asarray(matrix[unique_row, :]))
        #         search_scope, = np.where(row==1)
        #         list_groups = _reduce_matrix([unique_row], full_matrix, adjacency_list, found_isomorphisms = [], found_edges = [], favor_isomorphisms=largest_mapping_group)
        #         if list_groups: # if list group found
                    
        #             # get the largest list group with the fewest number of rows
        #             group_length = full_matrix[list_groups, :].sum()
        #             matrix[list_groups, :] = np.zeros((len(list_groups), n_cols))
        #             if (group_length > largest_group_length) or (group_length == largest_group_length and len(list_groups) < len(largest_mapping_group)):
        #                 largest_group_length = group_length
        #                 largest_mapping_group = list_groups
        #         else:
        #             continue
        #     # finally pick from unique_mapping_groups the largest list
        #     return largest_mapping_group

def find_fitting_lists(isomorphism_info):
    # return the indices of the lists that can be appended together to create a longer list
    # with no overlapping values in those lists
    def _create_choice_graph(isomorphism_info):
        all_values = set()
        for isomorphism in isomorphism_info:
            iso_ids, adjacencies = isomorphism
            for value in iso_ids:
                all_values.add(value)
        max_val = max(all_values)
        min_val = min(all_values)
        # create a matrix to find atom overlaps
        matrix_data = []
        adjacency_list = []
        for isomorphism in isomorphism_info:
            iso_ids, adjacencies = isomorphism
            cap_ids = []
            adjacencies_tupled = [] # adjacencies in tuple form for easier searching
            for bond, bond_type in adjacencies.items():
                edge_atom, cap_atom = bond
                cap_ids.append(cap_atom)
                adjacencies_tupled.append(tuple([edge_atom, cap_atom, bond_type]))
            adjacency_list.append(adjacencies_tupled)

            row = []
            for i in range(0, max_val+1):
                if i in iso_ids and i not in cap_ids:
                    row.append(1)
                else:
                    row.append(0)
            matrix_data.append(row)
            

        # matrix_data and adjacency_list are parallel to eachother
        # each item in adjacency_list corresponds to the isomorphism found in the 
        # corresponding row of matrix_data
        matrix = np.matrix(matrix_data)
        n_rows, n_cols = matrix.shape
        G = nx.Graph()
        for iso_id, isomorphism in enumerate(adjacency_list):
            n_atoms = len(isomorphism_info[iso_id][0])
            cap_atoms = []
            neighbors = []
            for bond in isomorphism:
                edge_a, cap_a, bond_type = bond
                cap_atoms.append(cap_a)
                neighbor_bond = (cap_a, edge_a, bond_type)
                for neighbor_id, neighbor_isomorphism in enumerate(adjacency_list):
                    for bond in neighbor_isomorphism:
                        if neighbor_bond == bond:
                            neighbors.append(tuple([neighbor_id, bond_type, edge_a]))

            row, = np.nonzero(np.squeeze(np.asarray(matrix[iso_id, :])))
            row = list(row)
            # find overlaps using matrix 
            overlaps, = np.nonzero(np.asarray(matrix[:, row].sum(axis=1)).reshape(-1))
            overlaps = list(overlaps)
            overlaps.remove(iso_id)
            G.add_node(
                iso_id,
                selected = False,
                overlaps = overlaps,
                bonds = isomorphism,
                n_atoms = n_atoms
            )
            for n_id, bond_type, edge_a in neighbors:
                G.add_edge(
                    iso_id,
                    n_id,
                    order = bond_type
                )
        return G

    def _traverse(queue, graph, found_nodes, found_edges):
        # searches for isomorphisms that fit together based on isomorphism adjacency info: 
        if queue == []:
            return found_nodes
        v = queue.pop(0)
        found_nodes.append(v)
        for bond in graph.nodes[v]['bonds']:
            edge_a, cap_a, bond_type = bond
            if {edge_a, cap_a} in found_edges:
                continue
            selected_neighbor = -1
            for neighbor in graph.neighbors(v):
                # make sure this is the neighbor with the correct bond info
                neighbor_node = graph.nodes[neighbor]
                if (cap_a, edge_a, bond_type) not in neighbor_node['bonds']:
                    continue
                if neighbor_node['selected'] == True:
                    selected_neighbor = -1
                    break
                selected_neighbor = neighbor
            if selected_neighbor >= 0:
                queue.append(selected_neighbor)
                found_edges.append({edge_a, cap_a})
    
        return _traverse(queue, graph, found_nodes, found_edges)


    # first, find unique rows to start the recursive algorithm
    # we only want one unique_col for each unique_row_ids
    choice_G = _create_choice_graph(isomorphism_info)
    connected_component_counts = []
    biggest_chain = None
    biggest_chain_length = 0
    for chain in nx.connected_components(choice_G):
        subgraph = choice_G.subgraph(chain)
        not_searched_nodes = list(subgraph.nodes)
        while len(not_searched_nodes) != 0:
            unique_group = _traverse([not_searched_nodes[0]], subgraph, [], [])
            new_tally = 0
            old_tally = 0
            overlapping_nodes = []
            for node in unique_group:
                new_tally += subgraph.nodes[node]['n_atoms']
                for overlapping_node in subgraph.nodes[node]['overlaps']:
                    if overlapping_node in subgraph.nodes and subgraph.nodes[overlapping_node]['selected'] == True:
                        overlapping_nodes.append(overlapping_node)
                        old_tally += subgraph.nodes[overlapping_node]['n_atoms']
            if new_tally > old_tally:
                # exchange new for old choices 
                for node in unique_group:
                    subgraph.nodes[node]['selected'] = True
                for node in overlapping_nodes:
                    subgraph.nodes[node]['selected'] = False
            [not_searched_nodes.remove(i) for i in unique_group]
        tally = 0
        for node in subgraph.nodes:
            if subgraph.nodes[node]['selected']:
                tally += subgraph.nodes[node]['n_atoms']
        if tally > biggest_chain_length:
            biggest_chain = subgraph
            biggest_chain_length = tally
        # finally pick from unique_mapping_groups the largest list
    largest_mapping_group = []
    for node in biggest_chain.nodes:
        if biggest_chain.nodes[node]['selected']:
            largest_mapping_group.append(node)
    return largest_mapping_group


if __name__ == "__main__":
    lists = [([4, 5, 6, 7, 8, 9, 10, 0, 11], {(4, 0): rdkit.Chem.rdchem.BondType.SINGLE, (8, 11): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([11, 8, 9, 10, 5, 6, 7, 12, 4], {(11, 12): rdkit.Chem.rdchem.BondType.SINGLE, (5, 4): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([26, 27, 28, 29, 30, 31, 32, 23, 33], {(26, 23): rdkit.Chem.rdchem.BondType.SINGLE, (30, 33): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([33, 34, 35, 36, 37, 38, 39, 30, 40], {(33, 30): rdkit.Chem.rdchem.BondType.SINGLE, (37, 40): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([33, 30, 31, 32, 27, 28, 29, 34, 26], {(33, 34): rdkit.Chem.rdchem.BondType.SINGLE, (27, 26): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([40, 37, 38, 39, 34, 35, 36, 41, 33], {(40, 41): rdkit.Chem.rdchem.BondType.SINGLE, (34, 33): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([55, 56, 57, 58, 59, 60, 61, 52, 62], {(55, 52): rdkit.Chem.rdchem.BondType.SINGLE, (59, 62): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([62, 59, 60, 61, 56, 57, 58, 63, 55], {(62, 63): rdkit.Chem.rdchem.BondType.SINGLE, (56, 55): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([101, 102, 103, 104, 105, 106, 107, 95, 108], {(101, 95): rdkit.Chem.rdchem.BondType.SINGLE, (105, 108): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([108, 105, 106, 107, 102, 103, 104, 109, 101], {(108, 109): rdkit.Chem.rdchem.BondType.SINGLE, (102, 101): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([108, 109, 110, 111, 112, 113, 114, 105, 115], {(108, 105): rdkit.Chem.rdchem.BondType.SINGLE, (112, 115): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([115, 112, 113, 114, 109, 110, 111, 116, 108], {(115, 116): rdkit.Chem.rdchem.BondType.SINGLE, (109, 108): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([108, 109, 110, 111, 112, 113, 114, 105, 115, 116], {(108, 105): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, 11], {(8, 11): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([11, 12, 13, 14, 15, 16, 17, 18, 19, 8, 20], {(11, 8): rdkit.Chem.rdchem.BondType.SINGLE, (14, 20): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([40, 41, 42, 43, 44, 45, 46, 47, 48, 37, 49], {(40, 37): rdkit.Chem.rdchem.BondType.SINGLE, (43, 49): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([62, 63, 64, 65, 66, 67, 68, 69, 70, 59, 71], {(62, 59): rdkit.Chem.rdchem.BondType.SINGLE, (65, 71): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([83, 84, 85, 86, 87, 88, 89, 90, 91, 80, 92], {(83, 80): rdkit.Chem.rdchem.BondType.SINGLE, (86, 92): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([92, 93, 94, 95, 96, 97, 98, 99, 100, 86, 101], {(92, 86): rdkit.Chem.rdchem.BondType.SINGLE, (95, 101): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([20, 21, 22, 23, 24, 25, 14, 26], {(20, 14): rdkit.Chem.rdchem.BondType.SINGLE, (23, 26): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([49, 50, 51, 52, 53, 54, 43, 55], {(49, 43): rdkit.Chem.rdchem.BondType.SINGLE, (52, 55): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([71, 72, 73, 74, 75, 76, 65, 77], {(71, 65): rdkit.Chem.rdchem.BondType.SINGLE, (74, 77): rdkit.Chem.rdchem.BondType.SINGLE}), 
            ([77, 78, 79, 80, 81, 82, 74, 83], {(77, 74): rdkit.Chem.rdchem.BondType.SINGLE, (80, 83): rdkit.Chem.rdchem.BondType.SINGLE})]
    
    # graph = nx.path_graph(20)
    # adj_matrix = nx.to_numpy_matrix(graph).astype(int)

    lists = find_fitting_lists(lists)
    print(lists)