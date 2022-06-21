import io
import pstats
import cProfile

import main

# FROM: https://gist.github.com/ralfstx/a173a7e4c37afa105a66f371a09aa83e
def prof_to_csv(prof):
    out_stream = io.StringIO()
    pstats.Stats(prof, stream=out_stream).print_stats()
    result = out_stream.getvalue()
    # chop off header lines
    result = 'ncalls' + result.split('ncalls')[-1]
    lines = [','.join(line.rstrip().split(None, 5)) for line in result.split('\n')]
    return '\n'.join(lines)


pr = cProfile.Profile()
pr.enable()
main.main(DEBUG=False, PRETTY=True, PROFILING=True)
pr.disable()
csv = prof_to_csv(pr)
with open("prof.csv", "w+") as file:
    file.write(csv)
