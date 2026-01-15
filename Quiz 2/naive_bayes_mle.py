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


def classify_map(
	priors: Dict[str, float],
	p_gc_given_y: Dict[str, Dict[str, float]],
	p_len_given_y: Dict[str, Dict[str, float]],
	p_comp_given_y: Dict[str, Dict[str, float]],
	x_gc: str,
	x_len: str,
	x_comp: str,
) -> Tuple[str, Dict[str, float], Dict[str, float]]:
	"""
	Compute unnormalized posteriors and normalized posteriors for the MAP decision.
	Returns (map_class, unnormalized, normalized).
	"""
	unnormalized: Dict[str, float] = {}
	for y in priors:
		unnormalized[y] = (
			priors[y]
			* p_gc_given_y[y].get(x_gc, 0.0)
			* p_len_given_y[y].get(x_len, 0.0)
			* p_comp_given_y[y].get(x_comp, 0.0)
		)
	total = sum(unnormalized.values())
	if total > 0.0:
		normalized = {y: unnormalized[y] / total for y in unnormalized}
	else:
		normalized = {y: 0.0 for y in unnormalized}
	map_class = max(unnormalized, key=unnormalized.get)
	return map_class, unnormalized, normalized


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

	# Part (c): MAP for x = (GC=Medium, Length=Long, Complexity=Low)
	x_gc, x_len, x_comp = "Medium", "Long", "Low"
	map_class, unnorm, norm = classify_map(
		priors, p_gc_given_y, p_len_given_y, p_comp_given_y, x_gc, x_len, x_comp
	)
	print("MAP classification for x = (GC=Medium, Length=Long, Complexity=Low)")
	print("Unnormalized posteriors (P(Y)*P(GC|Y)*P(Length|Y)*P(Complexity|Y)):")
	for y in sorted(unnorm.keys()):
		print(f"{y}\t{unnorm[y]:.6f}")
	print("Normalized posteriors:")
	for y in sorted(norm.keys()):
		print(f"{y}\t{norm[y]:.6f}")
	print(f"MAP class: {map_class}")
	print()

if __name__ == "__main__":
	main()


