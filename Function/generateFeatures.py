import os
import numpy as np
import pandas as pd
from Bio import SeqIO


seq_len = 150
alphabet = {'A': 0, 'U': 1, 'C': 2, 'G': 3}


def count_kmers(sequence, numFeature):
    k = 6
    feature = np.zeros(numFeature)
    for i in range(len(sequence) - k + 1):
        kmer = [alphabet[char] for char in sequence[i:i+k]]
        index = 0
        for digit in kmer:
            index = index * 4 + digit
        feature[index] += 1

    return feature


def calFeature_1D(fasta_pos, fasta_neg):
    numFeature = pow(4, 6)
    fasta_positives = SeqIO.parse(open(fasta_pos), 'fasta')
    fasta_negatives = SeqIO.parse(open(fasta_neg), 'fasta')

    dataset = []
    for fasta in fasta_positives:
        name, sequence = fasta.id, str(fasta.seq).upper()
        sequence = sequence.replace('T', 'U')
        sequence = [item for item in sequence if item.isupper()]
        sequence = ''.join(sequence)
        feature = count_kmers(sequence, numFeature)
        feature = np.append(feature, 1)
        dataset.append(feature)
    total = 2464977
    count = 0
    for fasta in fasta_negatives:
        if count % 10000 == 0:
            print(count, '/', total)
        count +=1
        # if count >= 5000000:
        name, sequence = fasta.id, str(fasta.seq).upper()
        sequence = sequence.replace('T', 'U')
        sequence = [item for item in sequence if item.isupper()]
        sequence = ''.join(sequence)
        feature = count_kmers(sequence, numFeature)
        feature = np.append(feature, 0)
        dataset.append(feature)

    dataset = np.array(dataset)

    # X = dataset[:, :-1]
    # Y = dataset[:, -1]
    return dataset


def calFeature_2D(fasta_pos, fasta_neg):
    with open('/mount/cheng-scratch/ydong/m6a/dict2D.txt') as file:
        motifList = file.readlines()
    motifList = [motif.strip() for motif in motifList]
    motifDict = {key: value for value, key in enumerate(motifList)}
    numFeature = len(motifDict)

    fasta_positives = SeqIO.parse(open(fasta_pos), 'fasta')
    fasta_negatives = SeqIO.parse(open(fasta_neg), 'fasta')

    dataset = []
    for fasta in fasta_positives:
        name, sequence = fasta.id, str(fasta.seq)
        feature = [0] * numFeature
        for i in range(len(sequence) - 8 + 1):
            motif = sequence[i: i + 8]
            if motif in motifDict:
                feature[motifDict[motif]] += 1

        feature = np.append(feature, 1)
        dataset.append(feature)

    for fasta in fasta_negatives:
        name, sequence = fasta.id, str(fasta.seq)
        feature = [0] * numFeature
        for i in range(len(sequence) - 8 + 1):
            motif = sequence[i: i + 8]
            if motif in motifDict:
                feature[motifDict[motif]] += 1

        feature = np.append(feature, 0)
        dataset.append(feature)

    dataset = np.array(dataset)

    # X = dataset[:, :-1]
    # Y = dataset[:, -1]
    return dataset


def calFeature_3D(fasta_pos, fasta_neg):
    fasta_positives = SeqIO.parse(open(fasta_pos), 'fasta')
    fasta_negatives = SeqIO.parse(open(fasta_neg), 'fasta')

    dataset = []
    for fasta in fasta_positives:
        name, sequence = fasta.id, str(fasta.seq)
        feature = sequence.split(',')
        feature = np.array([float(i) for i in feature])
        if len(feature) != 592:
            print('positive', name)

        feature[feature > 0] = 1
        feature[feature < 0] = 0
        feature = np.append(feature, 1)

        dataset.append(feature)

    for fasta in fasta_negatives:
        name, sequence = fasta.id, str(fasta.seq)
        feature = sequence.split(',')
        feature = np.array([float(i) for i in feature])
        if len(feature) != 592:
            print('negative', name)

        feature[feature > 0] = 1
        feature[feature < 0] = 0
        feature = np.append(feature, 0)
        dataset.append(feature)

    dataset = np.array(dataset)

    # X = dataset[:, :-1]
    # Y = dataset[:, -1]
    return dataset


