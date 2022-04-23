from typing import List
from multiprocessing import Pool
from .cas_to_smiles import get_cas_smiles


# нужна отдельно, чтобы мокать в pytest
def batch_cas_to_smiles(cas_list: List[str]):
    """
    :param cas_list: list of CAS numbers ['1-2-1', '14-1-5']
    :return: list of tuples [('1-2-1', 'COC'), ('14-1-5', None)]
    """

    #out = [get_cas_smiles(cas) for cas in cas_list]
    n = 50
    with Pool(processes=n) as pool:
        out = pool.map(get_cas_smiles, cas_list)
    return out
