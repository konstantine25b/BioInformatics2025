#!/usr/bin/env python2
import sys
import random

alphabet = ['A', 'G', 'C', 'T']

def read_pwm(pwm_file):
    pwm = []
    with open(pwm_file, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if i == 0:
                continue
            parts = line.strip().split()
            if len(parts) > 1 and parts[0] in alphabet:
                values = [float(x) for x in parts[1:]]
                pwm.append(values)
    return pwm

def generate_sequence_from_pwm(pwm, method='sample'):
    seq = ''
    motif_length = len(pwm[0])
    
    for pos in range(motif_length):
        probs = [pwm[base_idx][pos] for base_idx in range(len(alphabet))]
        
        if method == 'consensus':
            max_idx = probs.index(max(probs))
            seq += alphabet[max_idx]
        else:
            rand_val = random.random()
            cumsum = 0
            for i, p in enumerate(probs):
                cumsum += p
                if rand_val <= cumsum:
                    seq += alphabet[i]
                    break
    
    return seq

def main():
    if len(sys.argv) < 4:
        print "Usage: python pwm_to_fasta.py <pwm_file> <output_fasta> <num_sequences>"
        sys.exit(1)
    
    pwm_file = sys.argv[1]
    output_file = sys.argv[2]
    num_seqs = int(sys.argv[3])
    
    pwm = read_pwm(pwm_file)
    
    with open(output_file, 'w') as f:
        consensus = generate_sequence_from_pwm(pwm, method='consensus')
        f.write('>consensus\n')
        f.write(consensus + '\n')
        
        for i in range(num_seqs - 1):
            seq = generate_sequence_from_pwm(pwm, method='sample')
            f.write('>seq' + str(i+1) + '\n')
            f.write(seq + '\n')
    
    print "Generated " + str(num_seqs) + " sequences in " + output_file

if __name__ == '__main__':
    main()
