#!/usr/bin/env python2
import random
from nussinov import nussinov

def generate_random_sequence(length, gc_content=0.5):
    bases = ['A', 'U', 'G', 'C']
    
    if gc_content == 0.5:
        return ''.join(random.choice(bases) for _ in range(length))
    
    gc_count = int(length * gc_content)
    au_count = length - gc_count
    
    gc_bases = [random.choice(['G', 'C']) for _ in range(gc_count)]
    au_bases = [random.choice(['A', 'U']) for _ in range(au_count)]
    
    sequence = gc_bases + au_bases
    random.shuffle(sequence)
    
    return ''.join(sequence)

def test_length_100():
    print "=== Testing 100 random sequences of length 100 ==="
    scores = []
    for i in range(100):
        if i % 10 == 0:
            print "  Progress: %d/100" % i
        seq = generate_random_sequence(100)
        score = nussinov(seq)
        scores.append(score)
    
    avg_score = sum(scores) / float(len(scores))
    print "Average score: %.2f" % avg_score
    print "Min score: %d" % min(scores)
    print "Max score: %d" % max(scores)
    return avg_score

def test_varying_length():
    print "\n=== Testing score vs length ==="
    lengths = [20, 40, 60, 80, 100, 120, 140]
    
    with open('length_results.txt', 'w') as f:
        f.write('Length\tAverage_Score\tScore_per_base\n')
        for length in lengths:
            scores = []
            num_tests = 50
            for i in range(num_tests):
                seq = generate_random_sequence(length)
                score = nussinov(seq)
                scores.append(score)
            
            avg_score = sum(scores) / float(len(scores))
            score_per_base = avg_score / length
            print "Length %d: avg score = %.2f (%.4f per base)" % (length, avg_score, score_per_base)
            f.write('%d\t%.2f\t%.4f\n' % (length, avg_score, score_per_base))

def test_varying_gc():
    print "\n=== Testing score vs GC content ==="
    gc_contents = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    length = 100
    
    with open('gc_results.txt', 'w') as f:
        f.write('GC_Content\tAverage_Score\n')
        for gc in gc_contents:
            scores = []
            for i in range(50):
                seq = generate_random_sequence(length, gc_content=gc)
                score = nussinov(seq)
                scores.append(score)
            
            avg_score = sum(scores) / float(len(scores))
            print "GC content %.1f: avg score = %.2f" % (gc, avg_score)
            f.write('%.1f\t%.2f\n' % (gc, avg_score))

def main():
    random.seed(42)
    
    test_length_100()
    test_varying_length()
    test_varying_gc()
    
    print "\nResults saved to length_results.txt and gc_results.txt"

if __name__ == '__main__':
    main()
