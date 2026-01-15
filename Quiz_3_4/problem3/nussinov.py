#!/usr/bin/env python2
import sys

def can_pair(base1, base2):
    pairs = [('A', 'U'), ('U', 'A'), ('G', 'U'), ('U', 'G'), ('C', 'G'), ('G', 'C')]
    return (base1, base2) in pairs

def nussinov(sequence):
    n = len(sequence)
    dp = [[0 for _ in range(n)] for _ in range(n)]
    
    for length in range(2, n + 1):
        for i in range(n - length + 1):
            j = i + length - 1
            
            dp[i][j] = dp[i][j-1]
            
            for k in range(i, j):
                if can_pair(sequence[k], sequence[j]):
                    score = dp[i][k-1] if k > i else 0
                    score += dp[k+1][j-1] if k + 1 < j else 0
                    score -= 1
                    dp[i][j] = min(dp[i][j], score)
                
                bifurc = dp[i][k] + dp[k+1][j]
                dp[i][j] = min(dp[i][j], bifurc)
    
    return dp[0][n-1]

def main():
    if len(sys.argv) < 2:
        print "Usage: python nussinov.py <RNA_sequence>"
        sys.exit(1)
    
    sequence = sys.argv[1].upper()
    score = nussinov(sequence)
    print "Sequence: %s" % sequence
    print "Score: %d" % score

if __name__ == '__main__':
    main()
