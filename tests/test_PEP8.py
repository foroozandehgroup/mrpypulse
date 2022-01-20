from pathlib import Path

import pycodestyle


def test_style_mrpypulse():
    s = pycodestyle.StyleGuide()
    p = (Path(__file__).parents[1].resolve()) / "mrpypulse"
    res = s.check_files([f for f in p.iterdir() if f.suffix == ".py"])
    assert res.total_errors == 0


def test_style_examples():
    s = pycodestyle.StyleGuide()
    p = (Path(__file__).parents[1].resolve()) / "examples"
    res = s.check_files([f for f in p.iterdir() if f.suffix == ".py"])
    assert res.total_errors == 0


def test_style_tests():
    s = pycodestyle.StyleGuide()
    p = Path(__file__).parent
    res = s.check_files([f for f in p.iterdir() if f.suffix == ".py"])
    assert res.total_errors == 0
