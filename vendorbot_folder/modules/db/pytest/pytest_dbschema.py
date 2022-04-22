
import logging
from modules.db.dbschema import parse_cas_list
from datetime import datetime

logger = logging.getLogger(__name__)

"""
without smiles: [{'reagent_internal_id': '11c83c110a904ca59b171b6dd1e9f61c', 'CAS': '2749-11-3'},
{'reagent_internal_id': '7ebcbc3923524e998df1ac503ed1eb14', 'CAS': '75-64-9'},
{'reagent_internal_id': '3b16b8a33c06436caf23cc5c6221e601', 'CAS': '120-46-7'}]

2022-04-07 01:32:51,404 - root - INFO - resolved_list:
[{'reagent_internal_id': '11c83c110a904ca59b171b6dd1e9f61c', 'CAS': '2749-11-3', 'SMILES': 'C[C@H]([NH3+])CO'},
{'reagent_internal_id': '7ebcbc3923524e998df1ac503ed1eb14', 'CAS': '75-64-9', 'SMILES': 'CC(C)(C)N'},
{'reagent_internal_id': '3b16b8a33c06436caf23cc5c6221e601', 'CAS': '120-46-7', 'SMILES': 'O=C(CC(=O)c1ccccc1)c2ccccc2'}],
not_found_list: []

"""


def test_parse_cas_list(mock_batch_cas_to_smiles,
                        mock_is_similar):

    cas_list = ["15243-33-1"]

    # side_effect - список возвращаемых значений. 3 раза будем вызывать, значит нам надо 3 значения подставить
    mock_batch_cas_to_smiles.side_effect = [[("15243-33-1", "C")]]
    mock_is_similar.side_effect = [False]

    reagents, text_report = parse_cas_list(cas_list)

    assert len(reagents) == 1

    reagent = reagents[0]

    for f in ["reagent_internal_id", "CAS", "SMILES", "sharing_status", "timestamp"]:
        assert f in reagent

    assert reagent["sharing_status"] == "shared"
    assert reagent["CAS"] == "15243-33-1"
    assert reagent["SMILES"] == "C"
    assert datetime.strptime(reagent["timestamp"], "%d.%m.%Y %H:%M")

    assert "Строк в вашем списке <b>1</b>" in text_report
    assert "Правильных CAS-номеров <b>1</b>" in text_report
    assert "Опечатка в CAS: <b></b>" in text_report
    assert "Не найдено SMILES для: <b>0</b> позиций" in text_report
    assert "Ошибка обработки SMILES <b>0</b> позиций" in text_report
    assert "Найдено SMILES для: <b>1</b> реагентов" in text_report
    assert "Прекурсоров найдено и вычеркнуто: <b>0</b>" in text_report


def test_parse_cas_list_bad_cas(mock_batch_cas_to_smiles,
                                mock_is_similar):

    cas_list = ["1-1-1"]

    reagents, text_report = parse_cas_list(cas_list)

    assert not len(reagents)

    assert "Строк в вашем списке <b>1</b>" in text_report
    assert "Правильных CAS-номеров <b>0</b>" in text_report
    assert "Опечатка в CAS: <b>1-1-1</b>" in text_report
    assert "Не найдено SMILES для: <b>0</b> позиций" in text_report
    assert "Ошибка обработки SMILES <b>0</b> позиций" in text_report
    assert "Найдено SMILES для: <b>0</b> реагентов" in text_report
    assert "Прекурсоров найдено и вычеркнуто: <b>0</b>" in text_report


def test_parse_cas_list_fail_smiles(mock_batch_cas_to_smiles,
                                    mock_is_similar):

    cas_list = ["15243-33-1"]

    mock_batch_cas_to_smiles.side_effect = [[("15243-33-1", None)]]

    reagents, text_report = parse_cas_list(cas_list)

    assert not len(reagents)

    assert "Строк в вашем списке <b>1</b>" in text_report
    assert "Правильных CAS-номеров <b>1</b>" in text_report
    assert "Опечатка в CAS: <b></b>" in text_report
    assert "Не найдено SMILES для: <b>1</b> позиций\n15243-33-1" in text_report
    assert "Ошибка обработки SMILES <b>0</b> позиций" in text_report
    assert "Найдено SMILES для: <b>0</b> реагентов" in text_report
    assert "Прекурсоров найдено и вычеркнуто: <b>0</b>" in text_report


def test_parse_cas_list_precursor(mock_batch_cas_to_smiles,
                                  mock_is_similar):

    cas_list = ["15243-33-1"]

    mock_batch_cas_to_smiles.side_effect = [[("15243-33-1", "C")]]
    mock_is_similar.side_effect = [True]

    reagents, text_report = parse_cas_list(cas_list)

    assert not len(reagents)

    assert "Строк в вашем списке <b>1</b>" in text_report
    assert "Правильных CAS-номеров <b>1</b>" in text_report
    assert "Опечатка в CAS: <b></b>" in text_report
    assert "Не найдено SMILES для: <b>0</b> позиций" in text_report
    assert "Ошибка обработки SMILES <b>0</b> позиций" in text_report
    assert "Найдено SMILES для: <b>1</b> реагентов" in text_report
    assert "Прекурсоров найдено и вычеркнуто: <b>1</b>" in text_report


def test_parse_cas_list3(mock_batch_cas_to_smiles,
                         mock_is_similar):

    mock_batch_cas_to_smiles.side_effect = [[("15243-33-1", "C"), ("917-64-6", "CC"), ("94-02-0", "CCC")]]
    mock_is_similar.side_effect = [False, False, False]

    cas_list = ["15243-33-1", "917-64-6", "94-02-0"]

    reagents, text_report = parse_cas_list(cas_list)

    assert len(reagents) == 3

    reagent = reagents[0]

    for f in ["reagent_internal_id", "CAS", "SMILES", "sharing_status", "timestamp"]:
        assert f in reagent

    assert reagent["sharing_status"] == "shared"
    assert reagent["CAS"] == "15243-33-1"
    assert reagent["SMILES"] == "C"
    assert datetime.strptime(reagent["timestamp"], "%d.%m.%Y %H:%M")

    assert "Строк в вашем списке <b>3</b>" in text_report
    assert "Правильных CAS-номеров <b>3</b>" in text_report
    assert "Опечатка в CAS: <b></b>" in text_report
    assert "Не найдено SMILES для: <b>0</b> позиций" in text_report
    assert "Ошибка обработки SMILES <b>0</b> позиций" in text_report
    assert "Найдено SMILES для: <b>3</b> реагентов" in text_report
    assert "Прекурсоров найдено и вычеркнуто: <b>0</b>" in text_report
