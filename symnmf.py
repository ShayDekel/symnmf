import numpy as np
import pandas as pd
import sys
import symnmfmodule

beta = 0.5
max_iters = 300
eps = 1e-4
frobenius_norm = lambda x: np.linalg.norm(x, ord='fro')

def sym(X):
    return symnmfmodule.sym(X, len(X), len(X[0]))

def ddg(X):
    return symnmfmodule.ddg(X, len(X), len(X[0]))

def norm(X):
    return symnmfmodule.norm(X, len(X), len(X[0]))

def symnmf(W, H):
    n = len(W)
    k = len(H[0])
    ret = symnmfmodule.symnmf(W, H, n, k)
    return ret

def symnmf_from_X(X, k):
    W = norm(X)
    H = initialize_H(W, k).tolist()
    return symnmf(W, H)

def initialize_H(W, k):
    np.random.seed(0)
    m = np.mean(W)
    H = np.random.uniform(low=0, high=2 * np.sqrt(m / k), size=(len(W), k))
    return H

def print_vector(vec):
    d = len(vec)
    string = ''
    for i in range(d - 1):
        string += '{:.4f}'.format(vec[i]) + ','
    string += '{:.4f}'.format(vec[d - 1])
    print(string)

def print_matrix(mat):
    for i in range(len(mat)):
        print_vector(mat[i])

def main():
    if len(sys.argv) < 4:
        print("An Error Has Occurred")
        sys.exit(1)

    # Save user input arugemtns
    k = int(sys.argv[1])
    goal = sys.argv[-2]
    file_name = sys.argv[-1]

    # Read the input file
    X = pd.read_csv(file_name, header=None).values.tolist()
    n = len(X)

    # Choose function
    if goal == 'sym':
        A = sym(X)
        print_matrix(A)

    elif goal == 'ddg':
        D = ddg(X)
        print_matrix(D)

    elif goal == 'norm':
        W = norm(X)
        print_matrix(W)

    elif goal == 'symnmf':
        final_H = symnmf_from_X(X,k)
        print_matrix(final_H)

    else:
        print("An Error Has Occurred")


if __name__ == "__main__":
    main()
