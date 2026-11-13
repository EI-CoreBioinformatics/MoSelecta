"""MoSelecta package.

MoSelecta selects protein-coding gene models from multiple annotation sources
using CDS-structure-aware scoring and conservative pruning.
"""

from .constants import PROGRAM_NAME, PROGRAM_VERSION

__all__ = ["PROGRAM_NAME", "PROGRAM_VERSION"]
__version__ = PROGRAM_VERSION
