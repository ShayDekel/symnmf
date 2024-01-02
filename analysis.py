import numpy as np
from sklearn.metrics import silhouette_score
import sys
from symnmf import symnmf_from_X
import kmeans

ITERS = 200
EPS = 0.001

# Calculate clusters labels from the kmeans algorithm
def get_clusters_from_kmeans(X, clusters):
    n = len(X)
    ret = np.zeros(n)
    for i in range(n):
        # calculate list of euclidean distances squared between X[i] and clusters' rows
        distances = np.sum((X[i]-clusters)**2, axis=1)
        # X[i]'s cluster is the one closest to it in euclidean distance, either squared or not
        ret[i] = np.argmin(distances)
    return ret

def main():
    if len(sys.argv) < 3:
        pass
    
    k = int(sys.argv[1])
    file_name = sys.argv[2]

    X = np.loadtxt(file_name, delimiter=',')
    n = len(X)
    if k >= n:
        pass  # error here
    
    x_list = X.tolist()
    
    d = len(x_list[0])

    # Silhouette score of symnmf algorithm
    result_symnmf = symnmf_from_X(x_list, k)
    # X[i]'s cluster is j iff argmax(H[i])=j
    symnmf_clustering = np.argmax(result_symnmf, axis=1)
    symnmf_ss = silhouette_score(X, symnmf_clustering)
    print("nmf:", '{:.4f}'.format(symnmf_ss))

    # Silhouette score of kmeans algorithm
    kmeans_clusters = np.array(kmeans.kmeans(x_list, d, k, ITERS, EPS))
    kmeans_clustering = get_clusters_from_kmeans(X, kmeans_clusters)
    kmeans_ss = silhouette_score(X, kmeans_clustering)
    print("kmeans:", '{:.4f}'.format(kmeans_ss))

if __name__ == "__main__":
    main()
