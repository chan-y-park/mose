from api import analysis, save
from cmath import pi

def phase_scan(config, n_steps):
    for i in range(n_steps):
        phase = (pi / n_steps) * i
        print("\n\t-----------------------\n" \
            + "\n\t(%s) Analyzing phase %s\n" \
            + "\n\t-----------------------\n") % (i, phase)
        d = analysis(config, phase=phase)
        save(config, d)


