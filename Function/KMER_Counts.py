import os
import numpy as np

alphabet = {'A': 0, 'T': 1, 'C': 2, 'G': 3}

def count_kmers(sequence, numFeature):
    k = 6
    feature = np.zeros(numFeature)
    for i in range(len(sequence) - k + 1):
        kmer = [alphabet[char] for char in sequence[i:i+k]]
        index = 0
        index = int(index)
        for digit in kmer:
            index = index * 4 + digit
        feature[index] += 1
    return feature
