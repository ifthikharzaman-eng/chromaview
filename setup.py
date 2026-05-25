"""ChromaView — Sanger sequencing chromatogram analyzer."""
from setuptools import setup, find_packages

setup(
    name="chromaview",
    version="0.1.0",
    description="Modern chromatogram viewer for Sanger DNA sequencing",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="ChromaView Contributors",
    license="MIT",
    packages=find_packages(),
    python_requires=">=3.10",
    install_requires=[
        "PyQt6>=6.5",
        "pyqtgraph>=0.13.3",
        "biopython>=1.83",
        "numpy>=1.24",
    ],
    extras_require={
        "dev": ["pytest>=7.4"],
    },
    entry_points={
        "console_scripts": [
            "chromaview=chromaview.app:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
)
