import sys


def distance(vec1, vec2):
    lst = [(vec1[i] - vec2[i]) ** 2 for i in range(len(vec1))]
    ret = (sum(lst)) ** 0.5
    return ret


def sum_vectors(vec1, vec2):
    ret = [vec1[i] + vec2[i] for i in range(len(vec1))]
    return ret

def kmeans(vectors, d, k, iters, eps):
    centroids = [vectors[i].copy() for i in range(k)]

    for iteration in range(iters):
        clusters = [[] for i in range(k)]
        for vector in vectors:
            dist = distance(vector, centroids[0])
            min_dist_centroid_index = 0
            for j in range(k):
                new_dist = distance(vector, centroids[j])
                if new_dist < dist:
                    dist = new_dist
                    min_dist_centroid_index = j
            clusters[min_dist_centroid_index].append(vector)

        has_converged = True
        # update centroids
        for j in range(k):
            new_centroid = [0 for i in range(d)]
            for vector in clusters[j]:
                new_centroid = sum_vectors(new_centroid, vector)
            for i in range(d):
                new_centroid[i] /= (len(clusters[j]) if len(clusters[j]) else 1)
            if distance(new_centroid, centroids[j]) >= eps:
                has_converged = False
            centroids[j] = new_centroid

        if has_converged:  # check if lower than eps
            break

    return centroids


def main():
    args = sys.argv
    file_name = args[-1]
    eps = 0.001
    file = open(file_name, "r")
    data = file.readlines()
    vectors = [[float(value) for value in line.strip().split(",")] for line in data]
    n = len(vectors)
    d = len(vectors[0]) if n else 0
    err = 0
    try:
        k = int(args[1])
        if k <= 1 or k >= n:
            raise ValueError
    except Exception as e:
        print("Invalid number of clusters!")
        err = 1
    try:
        iters = int(args[2]) if len(args) == 4 else 200
        if iters <= 1 or iters >= 1000:
            raise ValueError
    except Exception as e:
        print("Invalid maximum iteration!")
        err = 1
    if err:
        return
    centroids = kmeans(vectors, d, k, iters, eps)
    for centroid in centroids:
        print(",".join(["%.4f" % centroid[i] for i in range(d)]))


if __name__ == "__main__":
    main()
