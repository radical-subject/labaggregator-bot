
import cirpy
import pubchempy
import traceback
import time
import re
from modules.reagent import Reagent

import logging
logger = logging.getLogger(__name__)

# Chemical Abstract Service Registry Number regular expression pattern
CAS_REGEXP = "^[1-9][0-9]{1,6}\\-[0-9]{2}\\-[0-9]$"


def is_cas_number(cas: str) -> bool:
    """
    Функция работает грамотно, проверено кровью, ее, сука, не трогать!
    """
    pattern = re.compile(CAS_REGEXP)
    if not pattern.match(cas):
        return False

    only_digits = cas.replace('-', '')
    control_digit = only_digits[-1]
    only_digits = only_digits[::-1]
    list_digits = list(only_digits[1:])

    sum_digits = 0
    index = 1
    for digit in list_digits:
        sum_digits += int(digit) * index
        index += 1

    return sum_digits % 10 == int(control_digit)


def get_reagent_cas_smiles(r: Reagent, delay: float = 0.2):
    try:
        r.smiles = cas_to_smiles(r.cas)
        time.sleep(delay)  # чтобы не грузить сервер
    except Exception as err:
        logger.error(err)
    return r


def cirpy_smiles_resolve(cas_or_name: str):
    """
    param: 106-95-6
    :return: C=CCBr  # TODO а может быть несколько SMILES?
    """
    try:
        return cirpy.resolve(cas_or_name, "smiles")
    except Exception as err:
        logger.warning(traceback.format_exc())   # посмотрим сетевые ошибки


def cirpy_cas_resolve(smiles_or_name: str):
    """
    :param smiles: C=CCBr
    :return: [106-95-6]
    """
    try:
        ret = cirpy.resolve(smiles_or_name, "cas")
        if ret:
            return [ret] if isinstance(ret, str) else ret
    except Exception as err:
        logger.warning(traceback.format_exc())   # посмотрим сетевые ошибки
    return []


def pubchempy_smiles_resolve(cas_or_name: str):
    """
    подходит для cas и для name.
    param: cas name. example: "1-1-1"
    TODO: кажется функция не работает. написал свою: pubchempy_get_smiles().
    """
    try:
        pubchem_response = pubchempy.get_compounds(cas_or_name, "name")
        return pubchem_response[0].isomeric_smiles if pubchem_response else None
    except Exception as err:
        logger.warning(traceback.format_exc())   # посмотрим сетевые ошибки


def cas_to_smiles(cas):
    """
    {
        "reagent_internal_id": uuid.uuid4().hex,
        "CAS": CAS_number
    }
    """
    res = cirpy_smiles_resolve(cas)
    if not res:
        logger.warning(f"cirpy: CAS({cas}) not found")
        res = pubchempy_smiles_resolve(cas)
        if not res:
            logger.warning(f"pubchempy: CAS({cas}) not found")
    return res


def what_reagent(text):
    """
    :param text: CAS or SMILES or name
    :return: [cas1, cas2], [smiles1, smiles2]
    """
    cas_list = []
    smiles_list = []

    if is_cas_number(text):
        cas_list.append(text)
        smiles_list.append(cas_to_smiles(text))
    else:
        for fun in (cirpy_smiles_resolve, pubchempy_smiles_resolve):
            smiles = fun(text)
            if smiles:
                smiles_list.append(smiles)      
                break  

        # smiles = cirpy_smiles_resolve(text)
        # if smiles == None:
        #     smiles = pubchempy_smiles_resolve(text)
        #     if smiles =! None
        # smiles_list.append(smiles)  # name2smiles
        

        # smiles_list = list(filter(None, smiles_list))
        for smiles in smiles_list:  # TODO если передали SMILES, нужен ли поиск по CAS?
            cas_list.extend(cirpy_cas_resolve(smiles))  # smiles2cas

        if smiles_list and text not in smiles_list:
            cas_list.extend(cirpy_cas_resolve(text))  # name2cas

    if cas_list:
        cas_list = list(set(cas_list))   # remove dublicates
        cas_list = list(filter(None, cas_list))  # remove None

    if smiles_list:
        smiles_list = list(set(smiles_list))
        smiles_list = list(filter(None, smiles_list))

    return cas_list, smiles_list
