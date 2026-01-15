#!/usr/bin/env python2
import sys
from collections import defaultdict

def analyze_kmers(seq_file, cons_file, k):
    with open(seq_file, 'r') as f:
        sequence = f.read().strip()
    with open(cons_file, 'r') as f:
        conservation = f.read().strip()
    
    kmer_count = defaultdict(int)
    kmer_conserved = defaultdict(int)
    
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k].upper()
        
        valid = True
        for base in kmer:
            if base not in ['A', 'C', 'G', 'T']:
                valid = False
                break
        
        if not valid:
            continue
        
        kmer_count[kmer] += 1
        
        cons_region = conservation[i:i+k]
        if cons_region.count('*') == k:
            kmer_conserved[kmer] += 1
    
    kmer_conservation_ratio = {}
    for kmer in kmer_count:
        if kmer_count[kmer] > 0:
            kmer_conservation_ratio[kmer] = float(kmer_conserved[kmer]) / kmer_count[kmer]
    
    return kmer_count, kmer_conserved, kmer_conservation_ratio

def main():
    seq_file = 'allinter'
    cons_file = 'allintercons'
    k = 6
    
    print "Analyzing %d-mers..." % k
    kmer_count, kmer_conserved, kmer_conservation_ratio = analyze_kmers(seq_file, cons_file, k)
    
    sorted_by_freq = sorted(kmer_count.items(), key=lambda x: x[1], reverse=True)
    sorted_by_cons = sorted(kmer_conservation_ratio.items(), key=lambda x: x[1], reverse=True)
    
    with open('top50_frequent.txt', 'w') as f:
        f.write('Rank\tKmer\tFrequency\tConserved\tConservation_Ratio\n')
        for i, (kmer, freq) in enumerate(sorted_by_freq[:50]):
            cons_count = kmer_conserved[kmer]
            cons_ratio = kmer_conservation_ratio[kmer]
            f.write('%d\t%s\t%d\t%d\t%.4f\n' % (i+1, kmer, freq, cons_count, cons_ratio))
    
    with open('top50_conserved.txt', 'w') as f:
        f.write('Rank\tKmer\tConservation_Ratio\tFrequency\tConserved\n')
        for i, (kmer, cons_ratio) in enumerate(sorted_by_cons[:50]):
            freq = kmer_count[kmer]
            cons_count = kmer_conserved[kmer]
            if freq >= 10:
                f.write('%d\t%s\t%.4f\t%d\t%d\n' % (i+1, kmer, cons_ratio, freq, cons_count))
    
    print "Done! Created top50_frequent.txt and top50_conserved.txt"
    print "Total unique %d-mers: %d" % (k, len(kmer_count))
    print "Top 10 most frequent:"
    for i, (kmer, freq) in enumerate(sorted_by_freq[:10]):
        print "  %d. %s: %d times (%.2f%% conserved)" % (i+1, kmer, freq, kmer_conservation_ratio[kmer]*100)

if __name__ == '__main__':
    main()
