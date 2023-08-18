import pytest
from unittest import mock


@pytest.fixture(autouse=True)
def mock_batch_reagent_cas_to_smiles() -> None:
    with mock.patch("modules.chem.batch.batch_reagent_cas_to_smiles") as mocker:
        yield mocker


@pytest.fixture(autouse=True)
def mock_cirpy_resolve() -> None:
    with mock.patch("cirpy.resolve") as mocker:
        yield mocker


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
