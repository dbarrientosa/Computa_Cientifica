import numpy as np
import scipy as sp
from numba import jit, float64

## Criterio de convergencia para update¿?
# This functions were taken from https://github.com/cmackenziek/tsfl

@jit(float64[:, :](float64[:, :], float64[:], float64[:], float64[:],
    float64[:], float64, float64), nopython=True)
def pairwise_tweds(TWED, A, A_times, B, B_times, lam=0.5, nu=1e-5):
    for i in range(1, len(A)):
        for j in range(1, len(B)):
            TWED[i, j] = min(
                # insertion
                (TWED[i - 1, j] + abs(A[i - 1] - A[i]) +
                 nu*(A_times[i] - A_times[i - 1]) + lam),
                # deletion
                (TWED[i, j - 1] + abs(B[j - 1] - B[j]) +
                 nu*(B_times[j] - B_times[j - 1]) + lam),
                # match
                (TWED[i - 1, j - 1] + abs(A[i] - B[j]) +
                 nu*(A_times[i] - B_times[j]) +
                 abs(A[i - 1] - B[j - 1]) +
                 nu*(A_times[i - 1] - B_times[j - 1]))
            )
    return TWED


def twed(A, A_times, B, B_times, lam=0.5, nu=1e-5):
    n, m = len(A), len(B)

    A, A_times = np.append(0.0, A), np.append(0.0, A_times)
    B, B_times = np.append(0.0, B), np.append(0.0, B_times)

    TWED = np.zeros((n + 1, m + 1))
    TWED[:, 0] = np.finfo(np.float).max
    TWED[0, :] = np.finfo(np.float).max
    TWED[0, 0] = 0.0

    TWED = pairwise_tweds(TWED, A, A_times, B, B_times, lam=lam, nu=nu)
return TWED[n, m]


def twed_distance(n):
    """
    Calculate the distance between i and j lightcurves.
    i = int(n/ncolumns).
    j = n - int(n/ncolumns).
    """
    return  twed(kmedoids[i], lightcurves[j])


def clasiffier(dist_matrix, lightcurves):
    """
    Returns an array with the index of the cluster of each lightcurve
    using the distance matrix.
    dist_matrix : Distance matrix, (cluster, lc).
    lightcurves : Array with the lightcurves (time, magnitude).
    """
    n_clusters = len(dist_matrix) # Number of clusters
    n_lc = len(lightcurves)  # Number of lightcurves
    cluster_index = sp.zeros(n_lc)
    for i in range(n_lc)
        lc = lightcurves[i]
        index = np.argmin(dist_matrix[:,i])
        cluster_index[i] = index
    return cluster_index


def update_center(lightcurves, medoid):
    """
    Update the representative of a cluster and return it.
    lightcurves : Array with the lightcurves (time, magnitude).
    medoid : Actual center of the cluster.
    """
    n_lc = len(lightcurves) ## NUmber of lightcurves in the cluster.
    distance_matrix = sp.zeros((n_lc,n_lc))
    print('Calculating distance matrix in the cluster...')
    for i in range(n_lc):
        for j in range(i+1, nlc):
            distance_matrix[i,j] = TWED(lightcurves[i,1], lightcurves[i,0],
                                        lightcurves[j,1], lightcurves[j,0])
            distance_matrix[j,i] = distance_matrix[i,j]
    distances = sp.sum(distance_matrix**2)
    index_min = sp.argmin(sp.sum(distances,axis=0))
    return(lightcurves[index_min])



def k_medoids_clustering(lightcurves, medoid_center, index=-1):
    """
    Clusters a sample of lightcurves using K-medoids algorithm with TWED metric.
    Returns the sample of lightcurves of index-cluster or all of them.
    lightcurves : Array with the lightcurves (time, magnitude)..
    medoids : Centers of the clusters (time, magnitude).
    index : index of the cluster returned.
    """
    n_clusters = len(dist_matrix) # Number of clusters
    n_lc = len(lightcurves)  # Number of lightcurves
    dist_matrix = sp.zeros((n_clusters,n_lc))
    ### Build the distance matrix
    print('Calculating distance matrix...')
    for i in range(n_clusters):
        for j in range(n_lc):
            dist_matrix = TWED(lightcurve[j,1], lightcurve[j,0],
                               medoids[i,1], medoids[i,0])
    lc_indexes = classifiers(dist_matrix, lightcurves) # lightcurves Clusters
    if index>=0 and index<n_clusters:
        lc_index = lightcurves[classifiers==index]
        return lc_index, medoids[index], dist_matrix
    else:
        return lc_indexes, medoids, dist_matrix
