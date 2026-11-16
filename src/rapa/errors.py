"""
Exceptions raised by RAPA.

Copyright (c) 2025 Ventus Therapeutics U.S., Inc.
Released under the MIT License (see LICENSE).
"""


class RapaInputError(ValueError):
    """Invalid input or options (missing file, unknown residue, bad --fix entry, ...)."""
