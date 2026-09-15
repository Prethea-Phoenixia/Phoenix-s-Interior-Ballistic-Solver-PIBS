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


def merge_sorted_records(record: list[tuple], new_points: list[tuple]) -> None:
    """
    Merge new_points into record in place, both sorted by the first element.
    Skips duplicates (same first element).
    """
    if not new_points:
        return

    existing_xs = {v[0] for v in record}
    unique_new = [v for v in new_points if v[0] not in existing_xs]

    if not unique_new:
        return

    # Merge two sorted lists
    i, j = 0, 0
    n, m = len(record), len(unique_new)
    merged = []

    while i < n and j < m:
        if record[i][0] <= unique_new[j][0]:
            merged.append(record[i])
            i += 1
        else:
            merged.append(unique_new[j])
            j += 1

    merged.extend(record[i:])
    merged.extend(unique_new[j:])

    # Update in place
    record.clear()
    record.extend(merged)
