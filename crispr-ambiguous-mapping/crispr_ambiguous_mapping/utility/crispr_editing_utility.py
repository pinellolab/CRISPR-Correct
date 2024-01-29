from typing import Counter as CounterType

def calculate_average_editing_frequency(total_mutation_counter: CounterType) -> float:
    total_editing_counts = (sum(total_mutation_counter.values()) - total_mutation_counter[0])
    total_counts = sum(total_mutation_counter.values())
    if total_counts > 0:
        editing_efficiency = total_editing_counts/total_counts
        return editing_efficiency
    else:
        return 0