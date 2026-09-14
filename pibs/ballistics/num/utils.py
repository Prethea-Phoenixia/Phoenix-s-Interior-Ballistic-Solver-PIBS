from __future__ import annotations


def find_last_le(record: list[tuple], target: float) -> int:
    """
    Find the last index where record[i][0] <= target using binary search.

    Assumes record is sorted by the first element.
    Returns 0 if no element satisfies the condition.
    """
    lo, hi = 0, len(record) - 1
    result = 0
    while lo <= hi:
        mid = (lo + hi) // 2
        if record[mid][0] <= target:
            result = mid
            lo = mid + 1
        else:
            hi = mid - 1
    return result
