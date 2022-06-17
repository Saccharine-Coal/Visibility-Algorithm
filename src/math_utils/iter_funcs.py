"""Collection of functions that perform operations on
ordered iterables: list, tuple ..."""


from math_utils import numerical

def pair(iterable, closed=False) -> list:
    """Arrange an iterable in pairs in the order of the original iterable.
    @param: closed = True pair last element with first.
    """
    # N = 3 => [1, 2, 3] => [(1, 2), (2, 3)] => Even size
    # N = 4 => [1, 2, 3, 4] => [(1, 2), (2, 3), (3, 4)] => odd size
    size = len(iterable)
    if size > 2:
        pairs = []
        for i, element in enumerate(iterable):
            if i < size - 1:
                pair_i = (element, iterable[i+1])
                pairs.append(pair_i)
            else:
                break
        if closed:
            # pair last with first
            pairs.append((iterable[-1], iterable[0]))
        return pairs
    elif size == 2:
        # smallest size
        return [(iterable[0], iterable[1])]
    else:
        # a iter of size 1 should not be allowed
        raise ValueError

def get_unique_ordered(iterable_of_arr) -> list:
    """Get a unique list of vals that preserves the order of the original iterable"""
    unique_list = []
    for arr in iterable_of_arr:
        if not numerical.is_array_in(arr, unique_list):
            # arr is not in unique
            unique_list.append(arr)
    return unique_list


def complementary_slice(iterable, start, stop, inclusive=False):
    """Get a slice of an iterable that is the complement of the
    given slice [start:stop:step=1].
    complementary_slice + slice = full iterable
    """
    if inclusive:
        return iterable[:start+1] + iterable[stop:]
    return iterable[:start] + iterable[stop:-1]
