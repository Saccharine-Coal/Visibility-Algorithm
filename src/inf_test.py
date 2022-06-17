a = (-1, 0, 1, float("inf"))
b = max(a)
print(b)


def three_switch_generator() -> bool:
    """A binary like switch that takes 3 iterations for the state to flip."""
    counter = 1
    while True:
        mod = counter % 3
        is_mod_zero = mod == 0
        if is_mod_zero:
            counter = 0
        counter += 1
        yield is_mod_zero


def four_switch_generator() -> bool:
    """A binary like switch that takes 4 iterations for the state to flip."""
    counter = 1
    while True:
        mod = counter % 4
        is_mod_zero = mod == 0
        if is_mod_zero:
            counter = 0
        counter += 1
        yield is_mod_zero


def binary_four_switch_generator(starting_state=False) -> bool:
    """A binary like switch that takes 2 iterations for the state to flip."""
    # True, True, False, False
    counter = 0
    while True:
        if counter < 2:
            yield starting_state
        else:
            if counter == 3:
                counter = -1
            yield not starting_state
        counter += 1

gener = binary_four_switch_generator()
for i in range(0, 20):
    result = next(gener)
    print(i, result)

for i in range(0, 10):
    print(i, i % 6)
