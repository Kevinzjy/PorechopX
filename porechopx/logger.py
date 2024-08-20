import sys
import logging


def get_logger(logger_name='logger', fname=None, debug=False):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    level = logging.DEBUG if debug else logging.INFO

    # LOG format
    fmt = "%(asctime)-15s %(levelname)-5s %(message)s"
    datefmt = "[%a %Y-%m-%d %H:%M:%S]"
    formatter = logging.Formatter(fmt, datefmt)

    # LOG file
    if fname is not None:
        fh = open(fname, 'w')
        fh.close()
        file_handler = logging.FileHandler(fname)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # LOG console
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    return logger
