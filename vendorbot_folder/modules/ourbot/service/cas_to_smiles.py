from typing import List
import sys
import argparse
import cirpy
import pubchempy
import requests
from multiprocessing import Pool
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


def banch_cas_to_smiles(cas_list: List[str]):
    """
    :param cas_list: list of CAS numbers ['1-2-1', '14-1-5']
    :return: list of tuples [('1-2-1', 'COC'), ('14-1-5', None)]
    """
    n = 50
    with Pool(processes=n) as pool:
        out = pool.map(get_cas_smiles, cas_list)
    return out


def cirpy_smiles_resolve(cas: str):
    """
    param: cas - line 1-1-1
    """
    return cirpy.resolve(cas, 'smiles')


def pubchempy_get_smileses(name: str):
    """
    dontsovcmc: я взял и накатал функцию... но возможно она не то возвращает)
    :param name: reagent name (eng)
    :return: list of smiles values
    """
    r = requests.post('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/JSON',
                      data=f'name={name}',
                      timeout=10.0)
    assert r.ok, f"pubchempy_get_smileses({name}) http error"

    j = r.json()

    smiles = []
    if 'PC_Compounds' in j:
        for compound in j['PC_Compounds']:
            if 'props' in compound:
                for prop in compound['props']:
                    if 'urn' in prop and 'label' in prop['urn']:
                        if prop['urn']['label'] == 'SMILES':
                            if 'value' in prop and 'sval' in prop['value']:
                                smiles.append(prop['value']['sval'])
                            else:
                                logger.warning(f"pubchem: incorrect SMILES value {name}")
                    else:
                        logger.warning(f"pubchem: incorrect prop {name}")
            else:
                logger.warning(f"pubchem: no props for {name}")

    return list(set(smiles))  # удалим дубликаты


def pubchempy_get_smiles(name: str):
    """
    :param name: reagent name (eng)
    :return: first smiles value found
    """
    smileses = pubchempy_get_smileses(name)
    if smileses:
        return smileses[0]


def pubchempy_smiles_resolve(cas: str):
    """
    подходит для cas и для name.
    param: cas name. example: "1-1-1"
    TODO: кажется функция не работает. написал свою: pubchempy_get_smiles().
    """
    pubchem_response = pubchempy.get_compounds(cas, "name")
    return pubchem_response[0].isomeric_smiles


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
        if res is None:
            pubchem_response = pubchempy.get_compounds(cas, "name")
            res = pubchem_response[0].isomeric_smiles
    except:
        try:
            pubchem_response = pubchempy.get_compounds(cas, "name")
            res = pubchem_response[0].isomeric_smiles
        except:
            logger.warning(f"CAS: {cas} not found")

    return res


def parse_lines(lines: List[str]):

    result = []
    for cas in lines:
        cas = cas.strip()
        try:
            smiles = cas_to_smiles(cas)
            result.append((cas, smiles))
            logger.info(f'{cas} - {smiles}')
        except Exception as err:
            logger.error(f'{cas}: {err}')
            result.append((cas, 'resolver_error'))

    result_object_list = [{"CAS": i[0], "SMILES": i[1]} for i in result] # if i[0]!="resolver_error"
    
    errors_CAS_list = [i[0] for i in result if i[1]=="resolver_error"]

    # # remove all error indications - this gives clean SMILES list
    # SMILES_list = list(filter(("resolver_error").__ne__, SMILES_list))
    
    return result_object_list, errors_CAS_list


if __name__ == "__main__":

    argparser = argparse.ArgumentParser()
    argparser.add_argument('-i', dest='inpath', help='txt file one line one cas')
    argparser.add_argument('-o', dest='outpath', help='outputfile')

    args = argparser.parse_args()

    if not args.inpath:
        logger.error('no infile argument')
        sys.exit(0)
    else:
        logger.info(f'Parse input file: {args.inpath}')

    with open(args.inpath, 'r') as infile:
        lines = infile.readlines()
        out = parse_lines(lines)
        print(out)
        if args.outpath:
            logger.info(f'Write output file: {args.outpath}')
            with open(args.outpath, 'w') as outfile:
                for cas, smiles in out:
                    outfile.write(f'{cas}\t{smiles}')
