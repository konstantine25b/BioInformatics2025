#!/usr/bin/env python2
import sys

def read_top_kmers(filename):
    kmers = []
    with open(filename, 'r') as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                kmers.append(parts[1])
    return kmers

def gc_content(seq):
    gc_count = seq.count('G') + seq.count('C')
    return float(gc_count) / len(seq)

def is_repeat(seq):
    if len(seq) % 2 == 0:
        half = len(seq) // 2
        if seq[:half] == seq[half:]:
            return True
    if len(seq) % 3 == 0:
        third = len(seq) // 3
        if seq[:third] == seq[third:2*third] == seq[2*third:]:
            return True
    return False

def read_known_motifs(filename):
    motifs = {}
    with open(filename, 'r') as f:
        next(f)
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                name = parts[0]
                motif = parts[1]
                motifs[name] = motif
    return motifs

def extract_6mers_from_motif(motif):
    sixmers = set()
    motif = motif.upper()
    clean_motif = ''
    for c in motif:
        if c in ['A', 'C', 'G', 'T']:
            clean_motif += c
    
    for i in range(len(clean_motif) - 5):
        sixmers.add(clean_motif[i:i+6])
    
    return sixmers

def main():
    frequent_kmers = read_top_kmers('top50_frequent.txt')
    conserved_kmers = read_top_kmers('top50_conserved.txt')
    
    print "=== ANALYSIS OF TOP 50 FREQUENT vs CONSERVED 6-MERS ==="
    print ""
    
    print "FREQUENT 6-MERS:"
    gc_values = [gc_content(k) for k in frequent_kmers]
    avg_gc = sum(gc_values) / len(gc_values)
    repeats = sum(1 for k in frequent_kmers if is_repeat(k))
    print "  Average GC content: %.2f%%" % (avg_gc * 100)
    print "  Number of repeats: %d / %d (%.1f%%)" % (repeats, len(frequent_kmers), 100.0*repeats/len(frequent_kmers))
    
    print ""
    print "CONSERVED 6-MERS:"
    gc_values_cons = [gc_content(k) for k in conserved_kmers]
    avg_gc_cons = sum(gc_values_cons) / len(gc_values_cons)
    repeats_cons = sum(1 for k in conserved_kmers if is_repeat(k))
    print "  Average GC content: %.2f%%" % (avg_gc_cons * 100)
    print "  Number of repeats: %d / %d (%.1f%%)" % (repeats_cons, len(conserved_kmers), 100.0*repeats_cons/len(conserved_kmers))
    
    print ""
    print "COMPARISON:"
    print "  GC difference: %.2f%% (conserved has higher GC)" % ((avg_gc_cons - avg_gc) * 100)
    print "  Repeat difference: %.1f%% (frequent has more repeats)" % (100.0*(repeats - repeats_cons)/len(frequent_kmers))
    
    print ""
    print "=== MATCHING WITH KNOWN YEAST MOTIFS ==="
    known_motifs = read_known_motifs('yeast_motifs.txt')
    
    frequent_set = set(frequent_kmers)
    conserved_set = set(conserved_kmers)
    
    matches_frequent = {}
    matches_conserved = {}
    
    for name, motif in known_motifs.items():
        motif_6mers = extract_6mers_from_motif(motif)
        
        freq_overlap = motif_6mers & frequent_set
        cons_overlap = motif_6mers & conserved_set
        
        if freq_overlap:
            matches_frequent[name] = freq_overlap
        if cons_overlap:
            matches_conserved[name] = cons_overlap
    
    print ""
    print "FOUND IN FREQUENT LIST: %d motifs" % len(matches_frequent)
    for name in sorted(matches_frequent.keys()):
        print "  %s: %s" % (name, ', '.join(sorted(matches_frequent[name])))
    
    print ""
    print "FOUND IN CONSERVED LIST: %d motifs" % len(matches_conserved)
    for name in sorted(matches_conserved.keys()):
        print "  %s: %s" % (name, ', '.join(sorted(matches_conserved[name])))

if __name__ == '__main__':
    main()
