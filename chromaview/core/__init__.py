"""Core modules: data models and file parsers."""
from .models import Chromatogram, ChromatogramMetadata, BASE_COLORS, BASE_COLORS_DARK
from .ab1_parser import parse_ab1
from .scf_parser import parse_scf
from .sequence_ops import reverse_complement, align_to_reference, search_sequence
