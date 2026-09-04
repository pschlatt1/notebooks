"""Minimal Matlab-style tic/toc timer.

Usage:
    from tictoc import tic, toc
    tic()
    ...work...
    toc()                  # prints elapsed time
    dt = toc(quiet=True)   # returns seconds without printing

Independent timers:
    tic('solve'); ...; toc('solve')
"""
import time

_start = {}

def tic(name='default'):
    """Start (or restart) the timer `name`."""
    _start[name] = time.perf_counter()

def toc(name='default', quiet=False):
    """Seconds elapsed since tic(name). Prints unless quiet=True."""
    if name not in _start:
        print(f"toc: timer '{name}' was never started")
        return None
    dt = time.perf_counter() - _start[name]
    if not quiet:
        print(f"Elapsed time{'' if name == 'default' else f' [{name}]'}: {dt:.3f} s")
    return dt
