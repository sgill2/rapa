"""
RAPA (Rotamer And Protonation Assignment): assigns protonation and rotamer states of HIS/ASN/GLN/ASP/GLU
and outputs all energetically accessible models. See docs/readme.md.

    import rapa
    result = rapa.run('protein.pdb', ligands=['LIG'], around=8)

Copyright (c) 2025 Ventus Therapeutics U.S., Inc.
Released under the MIT License (see LICENSE).
"""

__version__ = "1.2.0"

from .errors import RapaInputError  # noqa: E402
from .cli import run, RapaResult, main  # noqa: E402

__all__ = ['run', 'RapaResult', 'RapaInputError', 'main', '__version__']
