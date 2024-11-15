from main_script import rdrt


def test_transcribe():
    assert rdrt("ATG", "transcribe") == "AUG"
    assert rdrt("AGt", "transcribe") == "AGu"


def test_reverse():
    assert rdrt("ATG", "reverse") == "GTA"
    assert rdrt("cUG", "reverse") == "GUc"


def test_complement():
    assert rdrt("AtG", "complement") == "TaC"
    assert rdrt("CUG", "complement") == "GAC"


def test_reverse_complement():
    assert rdrt("ATg", "reverse_complement") == "cAT"
    assert rdrt("CUG", "reverse_complement") == "CAG"


def test_multiple_args():
    assert rdrt("ATG", "aT", "reverse") == ["GTA", "Ta"]
    assert rdrt("ttG", "AT", "ATc", "complement") == ["aaC", "TA", "TAg"]
