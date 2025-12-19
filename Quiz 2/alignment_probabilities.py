#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter
from math import comb, log, exp
from typing import Dict, List, Tuple
import random


def compute_column_score(column_chars: List[str]) -> int:
	"""
	Compute the alignment score for a column of 4 characters (A/C/G/T/-).
	Score is the number of unique pairs (among 6) that share the same symbol.
	Valid for 4 sequences: possible scores are {0, 1, 2, 3, 6}.
	"""
	if len(column_chars) != 4:
		raise ValueError("Expected exactly 4 characters per column.")
	counts = Counter(column_chars)
	# Sum of nC2 for each symbol group
	score = sum(comb(n, 2) for n in counts.values())
	return score


def compute_scores_for_alignment(alignment: List[str]) -> List[int]:
	"""
	Given an alignment with 4 equal-length sequences, compute the score per column.
	"""
	if len(alignment) != 4:
		raise ValueError("Alignment must contain exactly 4 sequences.")
	seq_lengths = {len(s) for s in alignment}
	if len(seq_lengths) != 1:
		raise ValueError("All sequences in an alignment must have equal length.")
	L = seq_lengths.pop()
	scores: List[int] = []
	for i in range(L):
		column = [alignment[r][i] for r in range(4)]
		scores.append(compute_column_score(column))
	return scores


def model_log_likelihood(scores: List[int], model_probs: Dict[int, float]) -> float:
	"""
	Compute log-likelihood of scores under given model probabilities.
	model_probs maps score -> P(score | model).
	"""
	ll = 0.0
	for s in scores:
		if s not in model_probs:
			raise ValueError(f"Score {s} not present in model probability table.")
		p = model_probs[s]
		if p <= 0.0:
			raise ValueError(f"Non-positive probability {p} for score {s}.")
		ll += log(p)
	return ll


def summarize_alignment(
	name: str,
	alignment: List[str],
	model_N: Dict[int, float],
	model_C: Dict[int, float],
) -> Tuple[List[int], float, float]:
	"""
	Compute per-column scores and log-likelihoods under models N and C.
	"""
	scores = compute_scores_for_alignment(alignment)
	ll_N = model_log_likelihood(scores, model_N)
	ll_C = model_log_likelihood(scores, model_C)

	print(f"\n=== {name} ===")
	print("Alignment:")
	for row in alignment:
		print(row)
	print(f"Scores per column: {scores}")
	print(f"log P(alignment | N): {ll_N:.6f}")
	print(f"log P(alignment | C): {ll_C:.6f}")
	try:
		print(f"P(alignment | N): {exp(ll_N):.6e}")
		print(f"P(alignment | C): {exp(ll_C):.6e}")
	except OverflowError:
		# If numbers are too small/large, skip exp
		pass
	print(f"Delta log-likelihood (C - N): {ll_C - ll_N:.6f}")
	return scores, ll_N, ll_C


def simulate_decision_rate_under_N(
	model_N: Dict[int, float],
	model_C: Dict[int, float],
	length: int = 10,
	num_sequences: int = 10_000,
	seed: int | None = 42,
) -> Tuple[float, float, int]:
	"""
	Simulate sequences of scores from N and compute how often P(S|C) > P(S|N).
	Returns (proportion, standard_error, count_C_greater_than_N).
	"""
	rng = random.Random(seed)
	support = sorted(model_N.keys())
	weights_N = [model_N[s] for s in support]

	count_favor_C = 0
	for _ in range(num_sequences):
		S = rng.choices(support, weights=weights_N, k=length)
		ll_N = model_log_likelihood(S, model_N)
		ll_C = model_log_likelihood(S, model_C)
		if ll_C > ll_N:
			count_favor_C += 1
	prop = count_favor_C / num_sequences
	se = (prop * (1.0 - prop) / num_sequences) ** 0.5
	return prop, se, count_favor_C


def simulate_decision_rate_under_C(
	model_N: Dict[int, float],
	model_C: Dict[int, float],
	length: int = 10,
	num_sequences: int = 10_000,
	seed: int | None = 123,
) -> Tuple[float, float, int]:
	"""
	Simulate sequences of scores from C and compute how often P(S|N) > P(S|C).
	Returns (proportion, standard_error, count_N_greater_than_C).
	"""
	rng = random.Random(seed)
	support = sorted(model_C.keys())
	weights_C = [model_C[s] for s in support]

	count_favor_N = 0
	for _ in range(num_sequences):
		S = rng.choices(support, weights=weights_C, k=length)
		ll_N = model_log_likelihood(S, model_N)
		ll_C = model_log_likelihood(S, model_C)
		if ll_N > ll_C:
			count_favor_N += 1
	prop = count_favor_N / num_sequences
	se = (prop * (1.0 - prop) / num_sequences) ** 0.5
	return prop, se, count_favor_N

def main() -> None:
	# Probability tables for scores given models N and C
	# Scores: 0, 1, 2, 3, 6
	model_N = {
		0: 0.1,
		1: 0.35,
		2: 0.25,
		3: 0.2,
		6: 0.1,
	}
	model_C = {
		0: 0.05,
		1: 0.15,
		2: 0.2,
		3: 0.3,
		6: 0.3,
	}

	# Two separate 4-sequence alignments of length 10 each
	alignment_1 = [
		"ACGACGACTA",
		"CAGACGCTGA",
		"TTCCTCTGAT",
		"AGATGTGACT",
	]
	alignment_2 = [
		"ACAACGAGTA",
		"AAAACGAATA",
		"TCATCGAGTT",
		"ACATCTAACT",
	]

	summarize_alignment("Alignment 1", alignment_1, model_N, model_C)
	summarize_alignment("Alignment 2", alignment_2, model_N, model_C)

	# Part (b): simulation under N
	prop, se, count = simulate_decision_rate_under_N(model_N, model_C, length=10, num_sequences=10_000, seed=42)
	print("\n=== Simulation (Part b) ===")
	print("Generated 10,000 sequences of length 10 from model N.")
	print(f"Count where P(S | C) > P(S | N): {count} / 10000")
	print(f"Proportion: {prop:.4f}")
	print(f"Approx 95% CI: [{max(0.0, prop - 1.96*se):.4f}, {min(1.0, prop + 1.96*se):.4f}]")

	# Part (c): simulation under C
	prop_c, se_c, count_c = simulate_decision_rate_under_C(model_N, model_C, length=10, num_sequences=10_000, seed=123)
	print("\n=== Simulation (Part c) ===")
	print("Generated 10,000 sequences of length 10 from model C.")
	print(f"Count where P(S | N) > P(S | C): {count_c} / 10000")
	print(f"Proportion: {prop_c:.4f}")
	print(f"Approx 95% CI: [{max(0.0, prop_c - 1.96*se_c):.4f}, {min(1.0, prop_c + 1.96*se_c):.4f}]")


if __name__ == "__main__":
	main()


