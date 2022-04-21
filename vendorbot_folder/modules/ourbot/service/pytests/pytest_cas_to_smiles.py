
from modules.ourbot.service.cas_to_smiles import pubchempy_smiles_resolve, \
    cirpy_smiles_resolve, cas_to_smiles, get_cas_smiles

from unittest.mock import patch


class PubChempyComponent:
    def __init__(self, isomeric_smiles):
        self.isomeric_smiles = isomeric_smiles


@patch("pubchempy.get_compounds", side_effect=[[PubChempyComponent("1")], [PubChempyComponent("2")]])
def test_pubchempy_smiles_resolve(pubchempy_smiles_resolve_patched):

    assert "1" == pubchempy_smiles_resolve("1-1-1")
    assert "2" == pubchempy_smiles_resolve("1-1-1")

    assert pubchempy_smiles_resolve_patched.call_count == 2


@patch("cirpy.resolve", side_effect=["1"])
def test_cirpy_smiles_resolve(cirpy_resolve_patched):

    assert "1" == cirpy_smiles_resolve("1")
    assert cirpy_resolve_patched.call_count == 1


def test_cas_to_smiles():

    with patch("cirpy.resolve", side_effect=["1"]) as cirpy_patched:
        with patch("pubchempy.get_compounds",
                   side_effect=[[PubChempyComponent("2")]]) as pubchempy_patched:
            assert "1" == cas_to_smiles("1")
            assert cirpy_patched.call_count == 1
            assert pubchempy_patched.call_count == 0

    with patch("cirpy.resolve", side_effect=[None]) as cirpy_patched:
        with patch("pubchempy.get_compounds",
                   side_effect=[[PubChempyComponent("2")]]) as pubchempy_patched:

            assert "2" == cas_to_smiles("1")
            assert cirpy_patched.call_count == 1
            assert pubchempy_patched.call_count == 1

    with patch("cirpy.resolve", side_effect=[Exception("test")]) as cirpy_patched:
        with patch("pubchempy.get_compounds",
                   side_effect=[[PubChempyComponent("2")]]) as pubchempy_patched:

            assert "2" == cas_to_smiles("1")
            assert cirpy_patched.call_count == 1
            assert pubchempy_patched.call_count == 1

    with patch("cirpy.resolve", side_effect=[Exception("test")]) as cirpy_patched:
        with patch("pubchempy.get_compounds",
                   side_effect=[Exception("test")]) as pubchempy_patched:

            assert not cas_to_smiles("1")
            assert cirpy_patched.call_count == 1
            assert pubchempy_patched.call_count == 1


def test_get_cas_smiles():

    with patch("modules.ourbot.service.cas_to_smiles.cas_to_smiles",
               side_effect=[None, "C", Exception("test")]) as patched:

        cas, smiles = get_cas_smiles("1")
        assert "1" == cas
        assert not smiles

        cas, smiles = get_cas_smiles("1")
        assert "1" == cas
        assert "C" == smiles

        cas, smiles = get_cas_smiles("1")
        assert "1" == cas
        assert not smiles
