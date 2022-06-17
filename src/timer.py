from time import time

def timer(func):
    # This function shows the execution time of
    # the function object passed
    def wrap_func(*args, **kwargs):
        t1 = time()
        result = func(*args, **kwargs)
        t2 = time()
        dt = t2 - t1
        print(f'Function {func.__name__!r} executed in {dt:.1f}s = {(dt*1000):.1f}ms = {dt*1000000:.1f}us')
        return result
    return wrap_func

@timer
def test():
    print('hi')


if __name__ == "__main__":
    test()
