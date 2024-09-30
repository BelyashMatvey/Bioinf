from typing import Dict, Union


def filter_fastq(seqs: Dict[str, tuple[str, str]], gc_bounds: Union[int, tuple[int, int]],
                 length_bounds: Union[int, tuple[int, int]] = [0, 2 ** 32], quality_threshold: int = 0) -> Dict[
    str, tuple[str, str]]:
    result_seqs: Dict[str, tuple[str, str]] = {}
    for name, fastq in seqs.items():
        if isinstance(gc_bounds, int):
            gc_bounds = (0, gc_bounds)
        if isinstance(length_bounds, int):
            length_bounds = (0, length_bounds)
        gc_content: float = (fastq[0].count('G') + fastq[0].count('C')) / len(fastq[0]) * 100
        threshold: int = 0
        for char in fastq[1]:
            threshold += ord(char)
        threshold /= len(fastq[0])
        if gc_bounds[0] <= gc_content <= gc_bounds[1] and length_bounds[0] <= len(fastq[0]) <= length_bounds[
            1] and threshold >= quality_threshold:
            result_seqs[name] = fastq

    return result_seqs
