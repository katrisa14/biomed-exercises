from .cli import cli
import logging
import sys

try:
    cli()
except SystemExit:
    raise
except:
    logging.error("Exception", exc_info=sys.exc_info())
    sys.exit(1)
