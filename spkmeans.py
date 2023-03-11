import sys
import numpy as np
import pandas as pd
import mykmeanssp

np.random.seed(0)

def print_and_exit():
    print("An Error Has Occurred")
    exit(1)

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
    print_row(eigen_vals)
    eigen_vecs = np.array(eigen_vecs).reshape(num_points, num_points) # todo - why?
    print_matrix(eigen_vecs) # print eigen_vecs.T - why?

def handle_spk(datapoints, num_points, point_size, k):
    u_matrix, tmp_k = mykmeanssp.get_u(datapoints, num_points, point_size)
    if(k == -1):
        k = tmp_k
    u_matrix = np.array(u_matrix).reshape(n,k)
    # todo print chosen indices?
    # todo kmeans

def run_goal(goal, datapoints, num_points, point_size, k):
    match goal:
        case "spk":
            handle_spk(datapoints, num_points, point_size, k)
        case "wam":
            result_matrix = mykmeanssp.get_wam(datapoints, num_points, point_size)
            print_matrix(result_matrix)
        case "ddg":
            result_matrix = mykmeanssp.get_ddg(datapoints, num_points, point_size)
            result_matrix = np.array(result_matrix).reshape(num_points, num_points) #todo - needed?
            print_matrix(result_matrix)
        case "gl":
            result_matrix = mykmeanssp.get_gl(datapoints, num_points, point_size)
            print_matrix(result_matrix)
        case "jacobi":
            handle_jacobi(datapoints, num_points, point_size)
        case _:
            print_and_exit() 


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


if __name__ == '__main__':
    main(sys.argv)