"""Package contract: version, import surface, no personal paths."""
from pathlib import Path

import wgdlife as wd

def test_version_and_surface():
    assert wd.__version__ == "1.0.0"
    assert wd.__all__ == ["cohort", "selection", "timing", "metastasis"]
    for m in wd.__all__:
        assert hasattr(wd, m)

def test_no_personal_paths_or_private_imports():
    src = Path(wd.__file__).parent
    for f in src.rglob("*.py"):
        txt = f.read_text(encoding="utf-8")
        assert "/home/" not in txt, f"{f}: absolute home path"
