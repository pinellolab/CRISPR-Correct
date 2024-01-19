from typing import Counter as CounterType

def calculate_average_editing_frequency(total_mutation_counter: CounterType) -> float:
    editing_efficiency = (sum(total_mutation_counter.values()) - total_mutation_counter[0])/sum(total_mutation_counter.values())
    return editing_efficiency