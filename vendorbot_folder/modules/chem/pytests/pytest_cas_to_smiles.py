
from unittest.mock import patch
from modules.chem.cas_to_smiles import pubchempy_smiles_resolve, \
    cirpy_smiles_resolve, cas_to_smiles, get_reagent_cas_smiles, is_cas_number, what_reagent

from modules.reagent import Reagent
from . import pubchempy_smiles_return, PubChempyComponent


@patch("pubchempy.get_compounds", side_effect=[pubchempy_smiles_return("1"), pubchempy_smiles_return("2")])
def test_pubchempy_smiles_resolve(pubchempy_smiles_resolve_patched):

    assert "1" == pubchempy_smiles_resolve("1-1-1")
    assert "2" == pubchempy_smiles_resolve("1-1-1")

    assert pubchempy_smiles_resolve_patched.call_count == 2


@patch("cirpy.resolve", side_effect=["1"])
def test_cirpy_smiles_resolve(cirpy_resolve_patched):

    assert "1" == cirpy_smiles_resolve("1")
    assert cirpy_resolve_patched.call_count == 1


def test_cas_to_smiles(mock_pubchempy_get_compounds,
                       mock_cirpy_resolve):

    mock_cirpy_resolve.side_effect = ["1"]
    mock_pubchempy_get_compounds.side_effect = [[PubChempyComponent("2")]]

    assert "1" == cas_to_smiles("1")
    assert mock_cirpy_resolve.call_count == 1
    assert mock_pubchempy_get_compounds.call_count == 0

    mock_cirpy_resolve.side_effect = [None]
    mock_pubchempy_get_compounds.side_effect = [[PubChempyComponent("2")]]

    assert "2" == cas_to_smiles("1")
    assert mock_cirpy_resolve.call_count == 2
    assert mock_pubchempy_get_compounds.call_count == 1

    mock_cirpy_resolve.side_effect = [Exception("it's OK test exception")]
    mock_pubchempy_get_compounds.side_effect = [[PubChempyComponent("2")]]

    assert "2" == cas_to_smiles("1")
    assert mock_cirpy_resolve.call_count == 3
    assert mock_pubchempy_get_compounds.call_count == 2

    mock_cirpy_resolve.side_effect = [Exception("it's OK test exception")]
    mock_pubchempy_get_compounds.side_effect = [Exception("it's OK test exception")]

    assert not cas_to_smiles("1")
    assert mock_cirpy_resolve.call_count == 4
    assert mock_pubchempy_get_compounds.call_count == 3

    # другой вариант patch прям в функции
    #with patch("cirpy.resolve", side_effect=[Exception("it's OK test exception")]) as cirpy_patched:
    #    with patch("pubchempy.get_compounds",
    #               side_effect=[Exception("it's OK test exception")]) as pubchempy_patched:


def test_get_cas_smiles():

    with patch("modules.chem.cas_to_smiles.cas_to_smiles",
               side_effect=[None, "C", Exception("it's OK test exception")]) as patched:

        r = get_reagent_cas_smiles(Reagent(cas="1"))
        assert "1" == r.cas
        assert not r.smiles

        r = get_reagent_cas_smiles(Reagent(cas="1"))
        assert "1" == r.cas
        assert "C" == r.smiles

        r = get_reagent_cas_smiles(Reagent(cas="1"))
        assert "1" == r.cas
        assert not r.smiles


def test_is_cas_number():
    assert is_cas_number('75-64-9')
    assert is_cas_number('120-46-7')

    assert not is_cas_number('+79168681111')
    assert not is_cas_number('110--86-1')
    assert not is_cas_number('102-95-5')
    assert not is_cas_number('@d12412')
    assert not is_cas_number('50 g')


def test_what_reagent(mock_cirpy_resolve,
                      mock_pubchempy_get_compounds):

    mock_cirpy_resolve.side_effect = [None]
    mock_pubchempy_get_compounds.side_effect = [[PubChempyComponent(None)]]

    cas_list, smiles_list = what_reagent("1-1-1")

    assert not cas_list
    assert not smiles_list

    #
    mock_cirpy_resolve.side_effect = ["COC(=O)CC#N"]

    cas_list, smiles_list = what_reagent("105-34-0")

    assert cas_list == ["105-34-0"]
    assert smiles_list == ["COC(=O)CC#N"]

    #
    mock_cirpy_resolve.side_effect = ["C1=CC=C(C=C1)CC(=O)CC2=CC=CC=C2", ["102-04-5"], ["102-04-5"], ["102-04-5"]]
    mock_pubchempy_get_compounds.side_effect = [[PubChempyComponent("O=C(Cc1ccccc1)Cc2ccccc2")]]

    cas_list, smiles_list = what_reagent("Benzyl ketone")

    assert cas_list == ["102-04-5"]
    assert "C1=CC=C(C=C1)CC(=O)CC2=CC=CC=C2" in smiles_list
    #assert "O=C(Cc1ccccc1)Cc2ccccc2" in smiles_list  TODO а куда делся?

    #
    mock_cirpy_resolve.side_effect = ["O=C(Cc1ccccc1)Cc2ccccc2", ["61346-73-4", "120-46-7"],
                                      ["61346-73-4", "120-46-7"], ["61346-73-4", "120-46-7"]]
    mock_pubchempy_get_compounds.side_effect = [[PubChempyComponent("C1=CC=C(C=C1)C(=O)CC(=O)C2=CC=CC=C2")]]

    cas_list, smiles_list = what_reagent("Dibenzoylmethane")

    assert "120-46-7" in cas_list
    assert "61346-73-4" in cas_list
    #assert "C1=CC=C(C=C1)C(=O)CC(=O)C2=CC=CC=C2" in smiles_list  TODO а куда делся?
    assert "O=C(Cc1ccccc1)Cc2ccccc2" in smiles_list
