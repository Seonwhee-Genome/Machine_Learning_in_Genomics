#!/usr/bin/python
##############################################################################
## The following python code is converted from Pseudocode written by Heng Li and Richard Durbin
## in Bioinformatics, Volume 25, Issue 14, 15 July 2009, Pages 1754–1760
## https://doi.org/10.1093/bioinformatics/btp324

##############################################################################
import os, re
import numpy as np


class BWA(object):

    def __init__(self, X="googol_"):
        self.X = X
        self.S = []  # Suffix array
        self.B = []  # BWT string
        self.B_inverse = []
        self.C = []
        self.O = []  # Occurrence array
        self.C_inverse = []
        self.O_inverse = []

    def Burrows_Wheeler_Transform(self, X, direction="Forward"):
        seqs = []
        i = 0

        for i in range(0, len(X)):
            if i > 0:
                suffix = X[0:i]
            else:
                suffix = ""
            seqs.append([X[i:] + suffix, i])
        array = np.array(seqs)
        array = array[array[:, 0].argsort()[::1]]

        for i in range(0, len(X)):
            self.S.append(array[i][1])
            if direction == "Forward":
                self.B.append(array[i][0][-1])
            else:
                self.B_inverse.append(array[i][0][-1])

    def Suffix_array_interval(self, W="LOL"):
        if W in self.X:
            R_lower = self.X.find(W)
            k = self.X.find(W[-1])
            R_upper = R_lower
            while W in self.X[k+1:]:
                R_upper = k + self.X[k+1:].find(W)
                k = k + self.X[k+1:].find(W[-1])
                print(k)

    def CalculateD(self, W="LOL"):
        k = 1
        l = len(self.X) - 1
        z = 0
        O_inverse = 0
        C = 0
        D = []
        for i in range(0, len(W)):
            for a in self.X[0:len(self.X)-2]:
                if W[i] > a:
                    C = C + 1

            if W[i] in self.B_inverse[0:k-1]:
                O_inverse = self.B_inverse[0:k-1].count(W[i])
            k = C + O_inverse + 1

            if W[i] in self.B_inverse[0:l]:
                O_inverse = self.B_inverse[0:l].count(W[i])
            l = C + O_inverse
            C = 0
            O_inverse = 0
            if k > l:
                k = l
                l = len(self.X) -1
                z = z + 1
            D.append(z)
        I = self.InexRecur(W, len(W)-1, z, 1, len(self.X)-1, D)
        print(I)

    def InexRecur(self, W, i, z, k, l, D):
        C = 0
        O = 0

        if z < D[i]:
            return [0]
        if i < 0:
            return [k, l]
        I = [0]
        I = I + self.InexRecur(W, i-1, z-1, k, l, D)  # insertions to self.X
        for b in ["A", "C", "G", "T"]:
            for a in self.X[0:len(self.X)-2]:
                if b > a:
                    C = C + 1
            if b in self.B[0:k-1]:
                O = self.B[0:k-1].count(b)
            k = C + O + 1

            if b in self.B[0:l]:
                O = self.B[0:l].count(b)
            l = C + O
            C = 0
            O = 0
            if k <= l:
                I = I + self.InexRecur(W, i, z-1, k, l, D)  # deletions from self.X
                if b == W[i]:
                    I = I + self.InexRecur(W, i-1, z, k, l, D)
                else:
                    I = I + self.InexRecur(W, i-1, z-1, k, l, D)
        return I



if __name__=="__main__":

    bwa = BWA("ATGCGCCA_")
    bwa.Burrows_Wheeler_Transform(X= bwa.X, direction="Forward")
    bwa.Burrows_Wheeler_Transform(X= ''.join(reversed(bwa.X)), direction="Reverse")
    bwa.CalculateD("CTGC")
