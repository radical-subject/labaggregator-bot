
import logging
from modules.reagent import Reagent
from modules.db.parser import parse_reagent_list
from datetime import datetime

logger = logging.getLogger(__name__)

"""
without smiles: [{'reagent_id': '11c83c110a904ca59b171b6dd1e9f61c', 'CAS': '2749-11-3'},
{'reagent_id': '7ebcbc3923524e998df1ac503ed1eb14', 'CAS': '75-64-9'},
{'reagent_id': '3b16b8a33c06436caf23cc5c6221e601', 'CAS': '120-46-7'}]

2022-04-07 01:32:51,404 - root - INFO - resolved_list:
[{'reagent_id': '11c83c110a904ca59b171b6dd1e9f61c', 'CAS': '2749-11-3', 'SMILES': 'C[C@H]([NH3+])CO'},
{'reagent_id': '7ebcbc3923524e998df1ac503ed1eb14', 'CAS': '75-64-9', 'SMILES': 'CC(C)(C)N'},
{'reagent_id': '3b16b8a33c06436caf23cc5c6221e601', 'CAS': '120-46-7', 'SMILES': 'O=C(CC(=O)c1ccccc1)c2ccccc2'}],
not_found_list: []

"""


def test_parse_cas_list(mock_batch_reagent_cas_to_smiles,
                        mock_is_similar):

    reagents_in = [Reagent(cas="15243-33-1")]

    # side_effect - список возвращаемых значений. 3 раза будем вызывать, значит нам надо 3 значения подставить
    mock_batch_reagent_cas_to_smiles.side_effect = [[Reagent(cas="15243-33-1",
                                                             smiles="C")]]
    mock_is_similar.side_effect = [False]

    reagents, text_report = parse_reagent_list(reagents_in)

    assert len(reagents) == 1

    r = reagents[0]

    assert r.sharing_status == "shared"  # REAGENT_SHARED
    assert r.cas == "15243-33-1"
    assert r.smiles == "C"
    assert datetime.strptime(r.timestamp, "%d.%m.%Y %H:%M:%S")  # DATETIME_FMT

    assert "Строк в вашем списке <b>1</b>" in text_report, text_report
    assert "Правильных CAS-номеров <b>1</b>" in text_report, text_report
    assert "Опечатка в CAS: <b>0</b>" in text_report, text_report
    assert "Не найдено SMILES для: <b>0</b> позиций" in text_report, text_report
    assert "Ошибка обработки SMILES <b>0</b> позиций" in text_report, text_report
    assert "Найдено SMILES для: <b>1</b> реагентов" in text_report, text_report
    assert "Прекурсоров найдено и вычеркнуто: <b>0</b>" in text_report, text_report


def test_parse_cas_list_bad_cas(mock_batch_reagent_cas_to_smiles,
                                mock_is_similar):

    reagents_in = [Reagent(cas="1-1-1")]

    reagents, text_report = parse_reagent_list(reagents_in)

    assert not len(reagents)

    assert "Строк в вашем списке <b>1</b>" in text_report, text_report
    assert "Правильных CAS-номеров <b>0</b>" in text_report, text_report
    assert "Опечатка в CAS: <b>1</b>" in text_report, text_report
    assert "Не найдено SMILES для: <b>0</b> позиций" in text_report, text_report
    assert "Ошибка обработки SMILES <b>0</b> позиций" in text_report, text_report
    assert "Найдено SMILES для: <b>0</b> реагентов" in text_report, text_report
    assert "Прекурсоров найдено и вычеркнуто: <b>0</b>" in text_report, text_report


def test_parse_cas_list_fail_smiles(mock_batch_reagent_cas_to_smiles,
                                    mock_is_similar):

    reagents_in = [Reagent(cas="15243-33-1")]

    mock_batch_reagent_cas_to_smiles.side_effect = [[Reagent(cas="15243-33-1",
                                                             smiles="")]]

    reagents, text_report = parse_reagent_list(reagents_in)

    assert not len(reagents)

    assert "Строк в вашем списке <b>1</b>" in text_report, text_report
    assert "Правильных CAS-номеров <b>1</b>" in text_report, text_report
    assert "Опечатка в CAS: <b>0</b>" in text_report, text_report
    assert "Не найдено SMILES для: <b>1</b> позиций\n" \
           "(cas=15243-33-1)" in text_report, text_report
    assert "Ошибка обработки SMILES <b>0</b> позиций" in text_report, text_report
    assert "Найдено SMILES для: <b>0</b> реагентов" in text_report, text_report
    assert "Прекурсоров найдено и вычеркнуто: <b>0</b>" in text_report, text_report


def test_parse_cas_list_precursor(mock_batch_reagent_cas_to_smiles,
                                  mock_is_similar):

    reagents_in = [Reagent(cas="15243-33-1")]

    mock_batch_reagent_cas_to_smiles.side_effect = [[Reagent(cas="15243-33-1",
                                                             smiles="C")]]
    mock_is_similar.side_effect = [True]

    reagents, text_report = parse_reagent_list(reagents_in)

    assert not len(reagents)

    assert "Строк в вашем списке <b>1</b>" in text_report, text_report
    assert "Правильных CAS-номеров <b>1</b>" in text_report, text_report
    assert "Опечатка в CAS: <b>0</b>" in text_report, text_report
    assert "Не найдено SMILES для: <b>0</b> позиций" in text_report, text_report
    assert "Ошибка обработки SMILES <b>0</b> позиций" in text_report, text_report
    assert "Найдено SMILES для: <b>1</b> реагентов" in text_report, text_report
    assert "Прекурсоров найдено и вычеркнуто: <b>1</b>" in text_report, text_report


def test_parse_cas_list3(mock_batch_reagent_cas_to_smiles,
                         mock_is_similar):

    mock_batch_reagent_cas_to_smiles.side_effect = [[Reagent(cas="15243-33-1",
                                                             smiles="C"),
                                                     Reagent(cas="917-64-6",
                                                             smiles="CC"),
                                                     Reagent(cas="94-02-0",
                                                             smiles="CCC")]]
    mock_is_similar.side_effect = [False, False, False]

    reagents_in = [Reagent(cas="15243-33-1"),
                   Reagent(cas="917-64-6"),
                   Reagent(cas="94-02-0")]

    reagents, text_report = parse_reagent_list(reagents_in)

    assert len(reagents) == 3

    assert "Строк в вашем списке <b>3</b>" in text_report, text_report
    assert "Правильных CAS-номеров <b>3</b>" in text_report, text_report
    assert "Опечатка в CAS: <b>0</b>" in text_report, text_report
    assert "Не найдено SMILES для: <b>0</b> позиций" in text_report, text_report
    assert "Ошибка обработки SMILES <b>0</b> позиций" in text_report, text_report
    assert "Найдено SMILES для: <b>3</b> реагентов" in text_report, text_report
    assert "Прекурсоров найдено и вычеркнуто: <b>0</b>" in text_report, text_report
