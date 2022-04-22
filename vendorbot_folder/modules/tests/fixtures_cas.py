import pytest
from unittest import mock


@pytest.fixture(autouse=True)
def mock_batch_cas_to_smiles() -> None:
    with mock.patch("modules.ourbot.service.batch.batch_cas_to_smiles") as mocker:
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
