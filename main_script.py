from typing import List, Union, Dict, Tuple
from utils.helpers import (
     is_valid_sequence, transcribe, reverse,
     complement, reverse_complement, gc_content
)

# --- Основные функции: run_dna_rna_tools и filter_fastq ---


def run_dna_rna_tools(*args: str) -> Union[str, List[str]]:
    """
    Запускает указанные процедуры над последовательностями.

    :param args: Последовательности (ДНК или РНК) и операция.
                 Операция должна быть последним аргументом.
    :return: Результат выполнения процедуры над последовательностью(ями).
    :raises ValueError: Если не указаны последовательности или процедура.
    """
    if len(args) < 2:
        raise ValueError(
            "Необходимо предоставить как минимум одну "
            "последовательность и операцию."
        )

    # Последний аргумент - это операция
    *sequences, operation = args

    valid_bases = {'A', 'T', 'C', 'G', 'U', 'a', 't', 'c', 'g', 'u'}
    valid_procedures = {'transcribe', 'reverse', 'complement',
                        'reverse_complement', 'gc_content'}

    if operation not in valid_procedures:
        raise ValueError(f"Недействительная процедура: {operation}")

    results = []

    for seq in sequences:
        if not is_valid_sequence(seq, valid_bases):
            raise ValueError(f"Недействительная последовательность: {seq}")
        if operation == 'transcribe':
            results.append(transcribe(seq))
        elif operation == 'reverse':
            results.append(reverse(seq))
        elif operation == 'complement':
            results.append(complement(seq))
        elif operation == 'reverse_complement':
            results.append(reverse_complement(seq))
        elif operation == 'gc_content':
            results.append(gc_content(seq))

    return results if len(results) > 1 else results[0]

# Добавим эрзац для функции run_dna_rna_tools как rdrt,
# чтобы соответствовать тестам
rdrt = run_dna_rna_tools


def filter_fastq(
    seqs: Dict[str, Tuple[str, str]],
    gc_bounds: Tuple[int, int] = (0, 100),
    length_bounds: Tuple[int, int] = (0, 2**32),
    quality_threshold: int = 0
) -> Dict[str, Tuple[str, str]]:
    """
    Фильтрует FASTQ последовательности по GC-составу,
    длине и среднему качеству.

    :param seqs: Словарь с последовательностями. Ключ — имя,
                 значение — кортеж (последовательность, качество).
    :param gc_bounds: Интервал GC-состава для фильтрации (в процентах).
    :param length_bounds: Интервал длины последовательности для фильтрации.
    :param quality_threshold: Пороговое значение среднего качества.
    :return: Отфильтрованный словарь с последовательностями.
    """

    # Если передано одно число, считаем его верхней границей
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)

    def is_within_bounds(value: float, bounds: Tuple[int, int]) -> bool:
        """Проверяет, находится ли значение в указанных пределах."""
        lower, upper = bounds
        return lower <= value <= upper

    filtered_seqs = {}

    for name, (sequence, quality) in seqs.items():
        # Проверяем длину последовательности
        seq_length = len(sequence)
        if not is_within_bounds(seq_length, length_bounds):
            continue

        # Рассчитываем процент GC
        gc_count = sequence.upper().count('G') + sequence.upper().count('C')
        gc_content = (
            (gc_count / len(sequence)) * 100
            if len(sequence) > 0 else 0.0
        )
        if not is_within_bounds(gc_content, gc_bounds):
            continue

        # Рассчитываем среднее качество
        total_quality = sum(ord(char) - 33 for char in quality)
        avg_quality = total_quality / len(quality)
        if avg_quality < quality_threshold:
            continue

        # Если все условия выполнены, добавляем последовательность в результат
        filtered_seqs[name] = (sequence, quality)

    return filtered_seqs


if __name__ == "__main__":
    pass
