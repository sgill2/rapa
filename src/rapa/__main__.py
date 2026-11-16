"""Allows `python -m rapa ...`."""
import sys

from .cli import main

main(sys.argv[1:])
