def is_valid_sequence(seq: str, valid_bases: set) -> bool:
    """
    Проверяет, является ли последовательность допустимой.
    """
    return all(base in valid_bases for base in seq)


def transcribe(dna: str) -> str:
    """
    Транскрибирует ДНК в РНК, заменяя T на U и t на u.
    """
    return dna.replace('T', 'U').replace('t', 'u')


def reverse(seq: str) -> str:
    """
    Возвращает развернутую последовательность.
    """
    return seq[::-1]


def complement(dna: str) -> str:
    """
    Возвращает комплементарную последовательность ДНК или РНК.
    """
    complement_dict = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
        'U': 'A', 'u': 'a'
    }
    return ''.join(complement_dict[base] for base in dna)


def reverse_complement(dna: str) -> str:
    """
    Возвращает обратно комплементарную последовательность ДНК или РНК.
    """
    return complement(dna)[::-1]


def gc_content(sequence: str) -> float:
    """
    Рассчитывает процент содержания G и C в последовательности.
    """
    if not sequence:
        raise ValueError("Последовательность не может быть пустой.")

    sequence_upper = sequence.upper()
    gc_count = sequence_upper.count('G') + sequence_upper.count('C')
    return (gc_count / len(sequence)) * 100
