import sys
import numpy as np
import traceback
import pandas as pd
import math
import mykmeanssp

np.random.seed(0)

def parse_datapoints(file_name):
    file_data = pd.read_csv(file_name, header=None)
    datapoints = file_data.to_numpy().tolist()
    return datapoints, file_data.shape[0], file_data.shape[1]

def print_row(row):
    print(",".join('{:.4f}'.format(np.round(num, 4)) for num in row))

def print_matrix(matrix):
    for row in matrix:
        print_row(row)

def handle_jacobi(datapoints, num_points, point_size):
    eigen_vals, eigen_vecs = mykmeanssp.get_jacobi(datapoints, num_points, point_size)
    if eigen_vals == None or eigen_vecs == None:
        exit(1)
    print_row(eigen_vals)
    eigen_vecs = np.array(eigen_vecs).reshape(num_points, num_points) # todo - why?
    print_matrix(eigen_vecs) # todo print eigen_vecs.T - why?

def get_dist(p, q):
    sum_part = 0
    for i in range(len(p)):
        sum_part += (p[i] - q[i]) ** 2
    return math.sqrt(sum_part)

def init_clusters_pp(k, datapoints, num_datapoints):
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

def kmeans_pp(vectors, K):
    """
    run k means on vectors, select initial centroids using
    k means++ strategy
    """
    N, d = vectors.shape
    selected = []
    min_dist = None
    p = None
    for Z in range(K):
        cur_select = np.random.choice(N, p=p)
        selected.append(cur_select)
        cur_dist = ((vectors - vectors[cur_select]) ** 2).sum(axis=1)
        min_dist = cur_dist if min_dist is None else np.minimum(min_dist, cur_dist)
        p = min_dist / min_dist.sum()
    c_res = mykmeanssp.fit(vectors[selected].tolist(),vectors.tolist(),K, d, K)
    final_centroids = np.array(c_res).reshape(K, d)
    print(','.join([str(i) for i in selected]))
    print_matrix(final_centroids)


def print_indices_chosen(indices):
    print(",".join([str(i) for i in indices]))


def handle_spk(datapoints, num_points, point_size, k):
    u_matrix, tmp_k = mykmeanssp.get_u(datapoints, num_points, point_size)
    if u_matrix == None:
        exit(1)
    if(k == -1):
        k = tmp_k

    u_matrix = np.array(u_matrix).reshape(num_points,tmp_k)
    num_points, point_size = u_matrix.shape
    
    chosen_clusters, chosen_indices = init_clusters_pp(k, u_matrix, num_points)
    print_indices_chosen(chosen_indices)
    clusters = mykmeanssp.fit(chosen_clusters, datapoints, k, num_points, point_size)
    if clusters == None:
        exit(1)
    print_matrix(clusters)
    #kmeans_pp(u_matrix, k)# todo del func maybe


def handle_matrix_output(result_matrix, num_points):
    if result_matrix == None:
        exit(1)
    result_matrix = np.array(result_matrix).reshape(num_points, num_points) #todo - needed?
    print_matrix(result_matrix)
        
def run_goal(goal, datapoints, num_points, point_size, k):
    if goal == "spk":
        handle_spk(datapoints, num_points, point_size, k)
    elif goal == "wam":
        handle_matrix_output(mykmeanssp.get_wam(datapoints, num_points, point_size), num_points)
    elif goal == "ddg":
        handle_matrix_output(mykmeanssp.get_ddg(datapoints, num_points, point_size), num_points)
    elif goal =="gl":
        handle_matrix_output(mykmeanssp.get_gl(datapoints, num_points, point_size), num_points)
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
        elif len(args) == 3:
            k = -1
            goal = args[1]
            file_name = args[2]
        else:
            exit(1)

        datapoints, num_points, point_size = parse_datapoints(file_name)
        run_goal(goal, datapoints, num_points, point_size, k)
    except Exception as e:
        print("An Error Has Occurred")
        print(e)
        print(traceback.format_exc())  
        exit(1)


if __name__ == '__main__':
    main(sys.argv)