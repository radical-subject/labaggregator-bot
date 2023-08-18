from typing import List
import time
import traceback
import logging
import uuid

from modules.db.blacklist import blacklist_engine
from modules.chem.cas_to_smiles import is_cas_number
from modules.chem import batch
from modules.chem import helpers
from modules.reagent import Reagent, REAGENT_SHARED, DATETIME_FMT

logger = logging.getLogger(__name__)


def parse_reagent_list(reagents_in: List[Reagent], contact: str = ''):
    """
    Фильтруем список CAS, ищем SMILES, удаляем прекурсоры,
    возвращаем список компонентов для БД и статистику
    :param cas_file:
    :param contact:
    :return:
    """
    now = time.strftime(DATETIME_FMT, time.localtime())

    good_cas = [r for r in reagents_in if is_cas_number(r.cas)]
    bad_cas_count = len(reagents_in) - len(good_cas)

    good_cas = batch.batch_reagent_cas_to_smiles(good_cas)
    no_smiles_list = [r for r in good_cas if not r.smiles]
    cas_smiles_list = [r for r in good_cas if r.smiles]

    cas_smiles_whitelist = []
    errors = []

    for r in cas_smiles_list:
        try:
            if not blacklist_engine.is_similar(r.smiles):
                cas_smiles_whitelist.append(r)

        except Exception as err:
            tb = traceback.format_exc()
            logger.error(f"is_similar failed for ({r}). Error: {tb}")
            errors.append(r)

    for r in cas_smiles_whitelist:
        """
        исправляем ошибки в зарядах: делаем все реагенты электронейтральными. 
        этот модуль не оттестирован на прочие ошибки, возможно он приведет 
        к ошибкам при импорте больших баз.
        """
        r.filter_smiles = helpers.filter_smiles_by_neutralize_atoms(r.smiles)

        r.inchikey_standard = helpers.smiles_to_inchikey(r.filter_smiles)

        r.reagent_id = uuid.uuid4().hex
        r.sharing_status = REAGENT_SHARED
        r.timestamp = now
        r.contact = contact

        #for option_values in ['location', 'name']:
        #    if option_values in valid_cas_list.columns:
        #        values = [i for i in set(valid_cas_list.loc[valid_cas_list['CAS'] == cas][option_values].dropna().values) if i]
        #        if values:
        #            r[option_values] = '\n'.join(values)

    message = f"file was successfully parsed and uploaded.\n"
    message += f"<b>import results</b>:\n"
    message += f"Строк в вашем списке <b>{len(reagents_in)}</b>\n"
    message += f"Правильных CAS-номеров <b>{len(good_cas)}</b>\n"
    message += f"Опечатка в CAS: <b>{bad_cas_count}</b>\n"
    message += f"Не найдено SMILES для: <b>{len(no_smiles_list)}</b> позиций\n"
    if no_smiles_list:
        message += "\n".join([str(r) for r in no_smiles_list]) + "\n"
    message += f"Ошибка обработки SMILES <b>{len(errors)}</b> позиций\n"
    if errors:
        message += "\n".join([str(r) for r in errors]) + "\n"
    message += f"Найдено SMILES для: <b>{len(cas_smiles_list)}</b> реагентов\n"
    message += f"Прекурсоров найдено и вычеркнуто: <b>{len(cas_smiles_list) - len(cas_smiles_whitelist)}</b>\n"

    return cas_smiles_whitelist, message
