"""
Shared pytest fixtures for the DAD publication package test suite.

Tests in this directory import dad.* from the co-located packaged code at
``Publication/code/dad/``. We add only that path so pytest does not silently
resolve imports against the development tree at ``06_Report/Mr_Pipeline``.

Codex Round 8 P1 fix: previous version inserted the development tree at
``sys.path[0]`` after the publication path, which made imports resolve to the
unpackaged source. The development-tree fallback has been removed.
"""

import sys
from pathlib import Path

# Add Publication/code/ to sys.path so `import dad` resolves to the
# packaged copy under Publication/code/dad/, not the development tree.
_CODE_DIR = Path(__file__).resolve().parent.parent
if str(_CODE_DIR) not in sys.path:
    sys.path.insert(0, str(_CODE_DIR))
