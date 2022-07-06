import numpy as np
from copy import deepcopy
import networkx as nx

# def reduce_matrix(matrix):
#     # reduces the matrix and returns the rows that were cut plus a smaller matrix for recursion
#     if matrix.size == 0:
#         return [[]]
#     n_rows, n_cols = matrix.shape
#     sums = np.asarray(matrix.sum(axis=0)).reshape(-1)
#     min_value = min(sums)
#     min_col_ids, = np.where(sums==min_value)
#     last_min_value_id = min_col_ids[-1]
    
#     # identify columns that have only single elements, and find the unique rows
#     fcolumn = np.squeeze(np.asarray(matrix[:, last_min_value_id]))
#     unique_row_lists = []
#     for value, frow_id in zip(fcolumn, range(0, n_rows)):
#         if value == 1:
#             # first, eliminate all rows with overlap to the current row
#             # then call reduce_matrix recursively for each occurance of 1
#             chosen_cols = np.asarray(matrix[frow_id, :].sum(axis=0)).reshape(-1).astype(bool)
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
    
#             cut_rows = list(set([frow_id]) | overlap_rows)
#             matrix_copy = deepcopy(matrix)
#             matrix_copy[list(cut_rows), :] = np.zeros((len(cut_rows), n_cols))
#             sub_matrix = matrix_copy[:, ~np.asarray(np.all(matrix_copy == 0, axis = 0)).reshape(-1)]
#             new_unique_rows = reduce_matrix(sub_matrix)
#             for unique_row in new_unique_rows:
#                 unique_row_lists.append([frow_id, *unique_row])
#     return unique_row_lists

def find_fitting_lists(lists, adj_matrix):
    # return the indices of the lists that can be appended together to create a longer list
    # with no overlapping values in those lists

    all_values = set()
    for l in lists:
        for value in l:
            all_values.add(value)
    max_val = max(all_values)
    min_val = min(all_values)
    matrix_data = []
    for l in lists:
        row = []
        for i in range(0, max_val+1):
            if i in l:
                row.append(1)
            else:
                row.append(0)
        matrix_data.append(row)
    matrix = np.matrix(matrix_data)
    n_rows, n_cols = matrix.shape
    # first, find unique rows to start the recursive algorithm
    unique_row_ids = list(range(0, len(lists)))
    # we only want one unique_col for each unique_row_ids
    unique_mapping_groups = []
    while len(unique_row_ids) != 0:
        unique_row = unique_row_ids[0]
        row = np.squeeze(np.asarray(matrix[unique_row, :]))
        search_scope, = np.where(row==1)
        list_groups = reduce_matrix(matrix, adj_matrix, search_scope.tolist())
        print(list_groups)
        if list_groups: # if list group found
            # get the largest list group with the fewest number of rows
            for list_group in list_groups:
                matrix[list_group, :] = np.zeros((len(list_group), n_cols))
                group_length = matrix[list_group, :].sum()
                unique_mapping_groups.append(tuple([group_length, list_group]))
        else: # returns False if there is a violation
            continue
        # remove the found unique mapping from unique_row_ids
        found_row_ids = []
        for list_group in list_groups:
            found_row_ids += list_group
        found_row_ids = list(set(found_row_ids))
        new_unique_row_ids = deepcopy(unique_row_ids)
        for row_id in unique_row_ids:
            if row_id in found_row_ids:
                new_unique_row_ids.remove(row_id)
        unique_row_ids = new_unique_row_ids
    # finally pick from unique_mapping_groups the largest list
    largest_group_id = 0
    largest_group_length = 0
    last_group_list_size = 0
    i = 0
    for group_length, mapping_group in unique_mapping_groups:
        if group_length >= largest_group_length and len(mapping_group) < last_group_list_size:
            largest_group_id = i
            largest_group_length = group_length
            last_group_list_size = len(mapping_group)
        i += 1
    l, largest_group = unique_mapping_groups[largest_group_id]
    return largest_group

def reduce_matrix(matrix, adj_matrix, search_scope=[]):
    def get_adjacencies(adj_matrix, col_ids):
        # returns the column ids of adjacent nodes to the given col_ids
        all_connections = np.asarray(adj_matrix[col_ids, :].sum(axis=0)).reshape(-1)
        all_connections[col_ids] = np.zeros(len(col_ids))
        adjacent_connections, = np.nonzero(all_connections)
        return adjacent_connections.tolist()
    # reduces the matrix and returns the rows that were cut plus a smaller matrix for recursion

    n_rows, n_cols = matrix.shape

    if search_scope == []:
        search_scope = list(range(0, n_cols))
    new_search_scope = []
    for idx in search_scope:
        if idx < n_cols:
            new_search_scope.append(idx)
    search_scope = new_search_scope
            
    sums = np.asarray(matrix.sum(axis=0)).reshape(-1)
    # required_unique_cols, = np.where(sums==1)
    # required_unique_row_ids, = np.nonzero(np.asarray(matrix[:, required_unique_cols].sum(axis=1)).reshape(-1))
    # required_unique_row_ids = set(required_unique_row_ids.tolist())
    allowed_sums = sums[search_scope]
    allowed_sums_ids, = np.where(allowed_sums > 0)
    allowed_nonzero_sums = allowed_sums[allowed_sums_ids]
    if allowed_nonzero_sums.size == 0:
        return [[]]
    min_allowed_value = min(allowed_nonzero_sums)
    min_ids, = np.where(allowed_sums==min_allowed_value)
    last_min_value_id = search_scope[min_ids[0]]
    
    # identify columns that have only single elements, and find the unique rows
    fcolumn = np.squeeze(np.asarray(matrix[:, last_min_value_id]))
    unique_row_lists = []

    for value, frow_id in zip(fcolumn, range(0, n_rows)):
        if value == 1:
            # first, eliminate all rows with overlap to the current row
            # then call reduce_matrix recursively for each occurance of 1
            chosen_cols = np.asarray(matrix[frow_id, :].sum(axis=0)).reshape(-1).astype(bool)
            chosen_col_ids, = np.where(chosen_cols==True)
            new_search_scope = get_adjacencies(adj_matrix, chosen_col_ids) + search_scope
            overlap_col_ids = np.array(range(0,len(sums)))[np.logical_and(sums > 1, chosen_cols)]
            # find rows that have overlapping values (and therefore overlapping substructures)
            overlap_rows = set()
            for col_id in overlap_col_ids:
                column = np.squeeze(np.asarray(matrix[:, col_id]))
                for row_val, row_id in zip(column, range(0, n_rows)):
                    # if row_val == 1 and row_id not in unique_rows:
                    #     overlap_rows.add(row_id)
                    if row_val == 1:
                        overlap_rows.add(row_id)
            # if the overlap_rows (which are discarded) share any rows with the required_unique_row_ids
            # then, this path of adjacent isomorphisms is invalid 
            # if len(overlap_rows & (required_unique_row_ids - {frow_id})) > 0:
            #     return False
    
            cut_rows = list(set([frow_id]) | overlap_rows)
            matrix_copy = deepcopy(matrix)
            matrix_copy[list(cut_rows), :] = np.zeros((len(cut_rows), n_cols))
            new_unique_rows = reduce_matrix(matrix_copy, adj_matrix, new_search_scope)
            if new_unique_rows == False:
                return False
            for unique_row in new_unique_rows:
                unique_row_lists.append([frow_id, *unique_row])
    return unique_row_lists

if __name__ == "__main__":
    lists = [
        [0,1,2,3,4],
        [0,1,2,3],
        [0,1],
        [2,3],
        [3,4],
        [4,5],
        [5,6]
        ]
    
    graph = nx.path_graph(20)
    adj_matrix = nx.to_numpy_matrix(graph).astype(int)

    lists = find_fitting_lists(lists, adj_matrix)
    print(lists)