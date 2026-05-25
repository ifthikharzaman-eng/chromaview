"""Analysis modules for quality, peaks, and mutation detection."""
from .quality import trim_low_quality, quality_summary
from .peaks import peak_heights, signal_to_noise
from .mutation import detect_ambiguous_bases, call_mutations
