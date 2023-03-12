import sys
import numpy as np
import pandas as pd
import mykmeanssp

np.random.seed(0)

def parse_datapoints(file_name):
    file_data = pd.read_csv(file_name, header=None)
    datapoints = file_data.to_numpy().tolist()
    return datapoints, data.shape[0], data.shape[1]

def print_row(row):
    print(",".join('{:.4f}'.format(np.round(num, 4) for num in row)))

def print_matrix(matrix):
    for row in matrix:
        print_row(row)

def handle_jacobi(datapoints, num_points, point_size, k):
    eigen_vals, eigen_vecs = mykmeanssp.get_jacobi(datapoints, num_points, point_size)
    if eigen_vals == None or eigen_vecs == None:
        exit(1)
    print_row(eigen_vals)
    eigen_vecs = np.array(eigen_vecs).reshape(num_points, num_points) # todo - why?
    print_matrix(eigen_vecs) # todo print eigen_vecs.T - why?

def init_clusters_pp(k, datapoints, indices, num_datapoints):
    np.random.seed(0)
    chosen_cent_indexes = []
    probabilities = None
    point_min_distances = np.inf # cause everything is smaller than +infinity

    for i in range(k):
        chosen = np.random.choice(num_datapoints, p=probabilities) # wont choose a chosen one cause probability will be zero
        chosen_cent_indexes.append(chosen)
        calc_dist_from_chosen = lambda p: get_dist(p, datapoints[chosen])
        distances = np.array(list(map(calc_dist_from_chosen, datapoints)))  # activates lambda on each element
        point_min_distances = np.minimum(point_min_distances, distances)
        probabilities = np.divide(point_min_distances, sum(point_min_distances))

    chosen_centroids = datapoints[chosen_cent_indexes].tolist()
    return chosen_centroids, chosen_cent_indexes

def print_indices_chosen(indices):
    print(",".join([str(i) for i in indices]))


def handle_spk(datapoints, num_points, point_size, k):
    u_matrix, tmp_k = mykmeanssp.get_u(datapoints, num_points, point_size)
    if u_matrix == None:
        exit(1)
    if(k == -1):
        k = tmp_k
    u_matrix = np.array(u_matrix).reshape(n,k)
    num_points, point_size = u_matrix.shape
    
    chosen_clusters, chosen_indices = init_clusters_pp(k, u_matrix, indices, num_points)
    print_indices_chosen(chosen_indices)
    clusters = mykmeanssp.fit(chosen_clusters, datapoints, k, num_points, point_size)
    if clusters == None:
        exit(1)
    print_matrix(clusters)

def handle_matrix_output(result_matrix):
    if result_matrix == None:
        exit(1)
    result_matrix = np.array(result_matrix).reshape(num_points, num_points) #todo - needed?
    print_matrix(result_matrix)
        
def run_goal(goal, datapoints, num_points, point_size, k):
    print(goal)
    if goal == "spk":
        handle_spk(datapoints, num_points, point_size, k)
    elif goal == "wam":
        handle_matrix_output(mykmeanssp.get_wam(datapoints, num_points, point_size))
    elif goal == "ddg":
        handle_matrix_output(mykmeanssp.get_ddg(datapoints, num_points, point_size))
    elif goal =="gl":
        handle_matrix_output(mykmeanssp.get_gl(datapoints, num_points, point_size))
    elif goal == "jacobi":
        handle_jacobi(datapoints, num_points, point_size)
    else:
        exit(1)


def main(args):
    try:
        if(len(args) == 4):
            k = int(args[1])
            goal = args[2]
            file_name = args[3]
        else:
            k = -1
            goal = args[1]
            file_name = args[2]

        datapoints, num_points, point_size = parse_datapoints(file_name)
        run_goal(goal, datapoints, num_points, point_size, k)
    except Exception as e:
        print("An Error Has Occurred")
        exit(1)


if __name__ == '__main__':
    main(sys.argv)