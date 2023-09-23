import pytest
from unittest import mock
from unittest.mock import MagicMock


@pytest.fixture(autouse=True)
def mock_batch_reagent_cas_to_smiles() -> None:
    with mock.patch("modules.chem.batch.batch_reagent_cas_to_smiles") as mocker:
        yield mocker


# тестовая БД для функции cirpy.resolve
cirpy_db = {
    "2749-11-3": "C[C@H]([NH3+])CO",
    "120-46-7": "O=C(CC(=O)c1ccccc1)c2ccccc2",
    "13896-65-6": "[Ru+3].[I-].[I-].[I-]",
    "15243-33-1": "C",
    "917-64-6": "CC",
    "94-02-0": "CCC"
}


@pytest.fixture(autouse=True)
def mock_cirpy_resolve() -> None:
    with mock.patch("cirpy.resolve") as mocker:
        def _mocked_function(cas):
            if cas in cirpy_db:
                return cirpy_db[cas]
            return None
        yield MagicMock(side_effect=_mocked_function)


@pytest.fixture(autouse=True)
def mock_is_similar() -> None:
    with mock.patch("modules.db.blacklist.BlackList.is_similar") as mocker:
        yield mocker


@pytest.fixture(autouse=True)
def mock_pubchempy_get_compounds() -> None:
    with mock.patch("pubchempy.get_compounds") as mocker:
        yield mocker


@pytest.fixture(autouse=True)
def mock_neutralize_atoms() -> None:
    with mock.patch("modules.chem.helpers.neutralize_atoms") as mocker:
        yield mocker


@pytest.fixture(autouse=True)
def mock_filter_smiles_by_neutralize_atoms() -> None:
    with mock.patch("modules.chem.helpers.filter_smiles_by_neutralize_atoms") as mocker:
        yield mocker


@pytest.fixture(autouse=True)
def mock_smiles_to_inchikey() -> None:
    with mock.patch("modules.chem.helpers.smiles_to_inchikey") as mocker:
        yield mocker
