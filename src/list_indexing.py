def list_to_dict(values) -> dict:
    # make sure vals are sorted
    vals_dict = {}
    for i, val in enumerate(values):
        vals_dict[i] = val
    return vals_dict



if __name__ == "__main__":
    l = [4, 3, 6, 7]
    l.sort()
    print(list_to_dict(l))
