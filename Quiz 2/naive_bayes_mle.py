#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter, defaultdict
from typing import Dict, List, Tuple

Record = Tuple[str, str, str, str]  # (GC, Length, Complexity, Class)


def load_training() -> List[Record]:
	"""
	Return the training dataset as (GC, Length, Complexity, Class) tuples.
	"""
	return [
		("Low", "Long", "High", "Gene"),
		("Low", "Long", "Low", "Gene"),
		("High", "Long", "High", "Repeat"),
		("Medium", "Short", "High", "Motif"),
		("Medium", "Short", "Low", "Motif"),
		("High", "Long", "Low", "Repeat"),
		("High", "Short", "High", "Motif"),
		("Medium", "Long", "High", "Gene"),
		("High", "Long", "Low", "Repeat"),
		("High", "Short", "High", "Motif"),
	]


def compute_priors(data: List[Record]) -> Dict[str, float]:
	class_counts = Counter(label for (_, _, _, label) in data)
	total = len(data)
	return {y: class_counts[y] / total for y in class_counts}


def compute_conditionals(
	data: List[Record],
	feature_index: int,
	feature_values: List[str],
) -> Dict[str, Dict[str, float]]:
	"""
	Compute P(X=val | Y=y) for a single feature given by feature_index in [0,1,2].
	Returns nested dict: {class: {feature_value: probability}} without smoothing (pure MLE).
	"""
	class_to_value_counts: Dict[str, Counter] = defaultdict(Counter)
	class_counts: Counter = Counter()
	for rec in data:
		value = rec[feature_index]
		y = rec[3]
		class_to_value_counts[y][value] += 1
		class_counts[y] += 1
	result: Dict[str, Dict[str, float]] = {}
	for y in class_counts:
		total_in_class = class_counts[y]
		result[y] = {v: class_to_value_counts[y][v] / total_in_class for v in feature_values}
	return result


def pretty_print_table(title: str, table: Dict[str, Dict[str, float]]) -> None:
	print(title)
	classes = sorted(table.keys())
	values = sorted({val for y in table for val in table[y].keys()})
	header = ["Class"] + values
	print("\t".join(header))
	for y in classes:
		row = [y] + [f"{table[y].get(v, 0.0):.3f}" for v in values]
		print("\t".join(row))
	print()


def main() -> None:
	data = load_training()

	priors = compute_priors(data)
	print("P(Y) (priors)")
	for y in sorted(priors.keys()):
		print(f"{y}\t{priors[y]:.3f}")
	print()

	# Define discrete domains
	gc_values = ["Low", "Medium", "High"]
	length_values = ["Long", "Short"]
	complexity_values = ["High", "Low"]

	p_gc_given_y = compute_conditionals(data, feature_index=0, feature_values=gc_values)
	p_len_given_y = compute_conditionals(data, feature_index=1, feature_values=length_values)
	p_comp_given_y = compute_conditionals(data, feature_index=2, feature_values=complexity_values)

	pretty_print_table("P(GC | Y)", p_gc_given_y)
	pretty_print_table("P(Length | Y)", p_len_given_y)
	pretty_print_table("P(Complexity | Y)", p_comp_given_y)


if __name__ == "__main__":
	main()


