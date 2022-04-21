
import cirpy
import pubchempy
import traceback
import time

import logging
logger = logging.getLogger(__name__)


def get_cas_smiles(cas: str, delay: float = 0.2):
    try:
        smiles = cas_to_smiles(cas)
        time.sleep(delay)  # чтобы не грузить сервер
        return cas, smiles
    except Exception as err:
        logger.error(err)

    return cas, None


def cirpy_smiles_resolve(cas: str):
    """
    param: 106-95-6
    :return: C=CCBr
    """
    return cirpy.resolve(cas, 'smiles')


def cirpy_cas_resolve(smiles: str):
    """
    :param smiles: C=CCBr
    :return: 106-95-6
    """
    return cirpy.resolve(smiles, 'cas')


def pubchempy_smiles_resolve(cas: str):
    """
    подходит для cas и для name.
    param: cas name. example: "1-1-1"
    TODO: кажется функция не работает. написал свою: pubchempy_get_smiles().
    """
    pubchem_response = pubchempy.get_compounds(cas, "name")
    return pubchem_response[0].isomeric_smiles if pubchem_response else None


def cas_to_smiles(cas):
    """
    {
        "reagent_internal_id": uuid.uuid4().hex,
        "CAS": CAS_number
    }
    """
    res = None
    try:
        res = cirpy.resolve(cas, "smiles")
    except Exception as err:
        logger.warning(traceback.format_exc())   # посмотрим сетевые ошибки
        logger.warning(f"cirpy: CAS({cas}) not found")

    try:
        if not res:
            res = pubchempy_smiles_resolve(cas)
    except Exception as err:
        logger.warning(traceback.format_exc())   # посмотрим сетевые ошибки
        logger.warning(f"pubchempy_smiles_resolve: CAS({cas}) not found")

    return res
