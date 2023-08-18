from typing import List
from multiprocessing import Pool
from .cas_to_smiles import get_reagent_cas_smiles
from modules.reagent import Reagent


# нужна отдельно, чтобы мокать в pytest
def batch_reagent_cas_to_smiles(reagents: List[Reagent]) -> List[Reagent]:
    """
    :param reagents: list of reagents with CAS numbers
            [Reagent(cas='1-2-1', smiles=None), Reagent(cas='14-1-5', smiles=None)]
    :return: list of reagents with filled SMILES if possible
            [Reagent(cas='1-2-1', smiles='COC'), Reagent(cas='14-1-5', smiles=None)]
    """
    with Pool(processes=50) as pool:
        out = pool.map(get_reagent_cas_smiles, reagents)
    return out
