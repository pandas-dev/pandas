def log_unhandled_exceptions(logger):
    def outer_wrapper(main):
        def wrapper(*args, **kwargs):
            try:
                main(*args, **kwargs)
            except Exception as e:
                logger.exception(f"ERROR: {e}")
                raise e

        return wrapper

    return outer_wrapper
