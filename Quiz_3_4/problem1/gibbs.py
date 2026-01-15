#!/usr/bin/env python
import sys
import string
import random

alphabet = ['A', 'G', 'C', 'T']

def build_pwm(S, positions, L, exclude_idx):
    counts = [[1 for _ in range(L)] for _ in range(len(alphabet))]
    
    for i, seq in enumerate(S):
        if i == exclude_idx:
            continue
        pos = positions[i]
        motif = seq[pos:pos+L]
        for j, base in enumerate(motif):
            if base in alphabet:
                counts[alphabet.index(base)][j] += 1
    
    pwm = [[0.0 for _ in range(L)] for _ in range(len(alphabet))]
    for j in range(L):
        col_sum = sum(counts[k][j] for k in range(len(alphabet)))
        for k in range(len(alphabet)):
            pwm[k][j] = float(counts[k][j]) / col_sum
    
    return pwm

def score_kmer(kmer, pwm):
    score = 1.0
    for i, base in enumerate(kmer):
        if base in alphabet:
            score *= pwm[alphabet.index(base)][i]
    return score

def sample_position(seq, L, pwm):
    scores = []
    for i in range(len(seq) - L + 1):
        kmer = seq[i:i+L]
        scores.append(score_kmer(kmer, pwm))
    
    total = sum(scores)
    if total == 0:
        return random.randint(0, len(seq) - L)
    
    probs = [s/total for s in scores]
    rand_val = random.random()
    cumsum = 0
    for i, p in enumerate(probs):
        cumsum += p
        if rand_val <= cumsum:
            return i
    return len(seq) - L

def GibbsSampler(S, L):
    positions = [random.randint(0, len(seq) - L) for seq in S]
    
    iterations = 1000
    for iteration in range(iterations):
        exclude_idx = random.randint(0, len(S) - 1)
        
        pwm = build_pwm(S, positions, L, exclude_idx)
        
        new_pos = sample_position(S[exclude_idx], L, pwm)
        positions[exclude_idx] = new_pos
    
    final_pwm = build_pwm(S, positions, L, -1)
    
    return final_pwm

def main():
    L = int(sys.argv[1])
    datafile = sys.argv[2]
    S = readdata(datafile)
	
    P = GibbsSampler(S,L)
	
    print "    ", 
    for i in range(L):
        print "%-5d " % (i+1),
    print ""
	
    for j in range(len(alphabet)):
        print " %s " % alphabet[j], 
        for i in range(L):
            print " %5.3f" % P[j][i],
        print ""
	
def readdata(file):
    data = [];
    for line in open(file,'r'):
        data.append(line[0:-1])
    return data

main()

