from typing import Iterable


def all_equal(iterable: Iterable):
    iterator = iter(iterable)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(item == first for item in iterator)


def contains_duplicates(iterable: Iterable):
    iterator = iter(iterable)
    items = set()
    for item in iterator:
        if item in items:
            return True
        items.add(item)
    return False
