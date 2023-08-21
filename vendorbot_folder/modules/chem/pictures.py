
from typing import List, Tuple
import os
import traceback
import logging

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import SimilarityMaps
except Exception:
    pass

PICTURES_PATH = os.getenv('PICTURES_PATH')

def create_smiles_picture(smiles: str, path: str = PICTURES_PATH) -> str:
    """
    :param smiles: TODO опять глупый вопрос нужен filter_smiles или обычный
    :return:
    """
    smiles = smiles.replace("|", "")  # вот этого не было и поэтому поиск падал

    fpath = os.path.join(path, smiles + '.png')
    if not os.path.exists(fpath):  # TODO нужно ли перегенерировать или картинки всегда одинаковые?
        mol = Chem.MolFromSmiles(smiles)
        refmol = Chem.MolFromSmiles(smiles)

        fp = SimilarityMaps.GetAPFingerprint(mol, fpType='normal')
        fp = SimilarityMaps.GetTTFingerprint(mol, fpType='normal')
        fp = SimilarityMaps.GetMorganFingerprint(mol, fpType='bv')
        fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(refmol, mol, SimilarityMaps.GetMorganFingerprint)

        fig.savefig(fpath, bbox_inches="tight")
    return fpath


def create_similar_smiles_grid_picture(request_smiles: str, molecules: List[Tuple[str, str, float]],
                                       path: str = PICTURES_PATH):
    """
    TODO smiles - или всё же filter_smiles?
    TODO добавить контакты реактивов. Нужно вынести 1 фукнцию, которая ищет по БД все реактивы и возвращает нам и её результат
    передавать сюда в виде  smile, similarity, contact, location
    """
    def mol_from_smiles(s: str):
        # вертикальная черта в SMILES - непонятно что несёт, и RDKIT ее не понимает, убираем ее
        s = s.replace("|", "")
        return Chem.MolFromSmiles(s)

    ms = [mol_from_smiles(request_smiles), ]
    legends = ["requested structure", ]
    for molecule in molecules:
        same_inchikey, same_smiles, similarity = molecule
        ms.append(mol_from_smiles(same_smiles))
        legends.append(f"Similarity={similarity * 100:.2f}% SMILES:{same_smiles}\n")

    img = Draw.MolsToGridImage(ms, molsPerRow=3, subImgSize=(400, 400), legends=legends)

    fname = request_smiles.replace("|", "")
    fpath = os.path.join(path, fname + '_grid.png')
    img.save(fpath, bbox_inches="tight")
    return fpath
