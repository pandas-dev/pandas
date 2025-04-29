from numba.testing._runtests import _main


if __name__ == '__main__':
    import sys
    # For parallel testing under Windows
    from multiprocessing import freeze_support
    freeze_support()
    sys.exit(0 if _main(sys.argv) else 1)
