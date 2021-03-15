import logging
import sys

# "DEBUG" < "INFO" < "WARNING" < "ERROR" < "CRITICAL"
LOW_PASS_STREAM = sys.stdout
HIGH_PASS_STREAM = sys.stderr
LOG_LVL = logging.DEBUG
# Min inclusive log level for which messages will be sent to stderr
FILTER_LVL = logging.WARNING

class LessThanFilter(logging.Filter):
    def __init__(self, max_level):
        super(LessThanFilter).__init__()
        self.max_level = max_level

    def filter(self, record):
        # Log if severity level lower than set level (filter_lvl).
        return record.levelno < self.max_level

def set_logging_level(log_lvl=LOG_LVL, filter_lvl=FILTER_LVL,
                      low_pass_stream=LOW_PASS_STREAM,
                      high_pass_stream=HIGH_PASS_STREAM):
    """
    Allow output from logs below a certain level to be sent to stdout, and the
    rest to stderr (default for logging).
    Changing stdout/stderr allows for sending to alternative StreamHandler.
    Adapted from SO.
    """
    # Responsible for logs below filter_lvl
    handler_low = logging.StreamHandler(low_pass_stream)
    handler_low.addFilter(LessThanFilter(filter_lvl))
    handler_low.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))

    # Responsible for logs filter_lvl and above
    handler_high = logging.StreamHandler(high_pass_stream)
    handler_high.setLevel(filter_lvl)
    handler_high.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))

    logger = logging.getLogger()  # Get the root logger
    logger.handlers = []      # Remove any previous handlers
    logger.setLevel(log_lvl)  # Defaults to logging.WARNING
    logger.addHandler(handler_low)
    logger.addHandler(handler_high)
    return