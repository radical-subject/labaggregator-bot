
import traceback
import logging

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem.rdmolops import GetFormalCharge
    from mongordkit.Database import registration
    from mongordkit.Database import write
    from mongordkit.Search import similarity, substructure
except:
    pass


# Импортировать только так: from modules.chem import helpers


def neutralize_atoms(mol):
    """
    atom neutralizer
    """
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


def filter_smiles_by_neutralize_atoms(smiles: str):
    """
    TODO: добавить обработку ошибок.
    Почему в black_list ищется без этой фильтрации smiles,
    а в БД уникальных молекул добавляется после неё?
    Добавлять в БД если эта функция завершиться с ошибкой?

    Neutralize molecules atom by atom
    """
    try:
        smiles = smiles.replace("|", "")  # for RDKIT
        # Create RDKit molecular objects
        mol = Chem.MolFromSmiles(smiles)
        before = GetFormalCharge(mol)
        if before != 0:
            mol = neutralize_atoms(mol)
            after = GetFormalCharge(mol)
            if before != after:
                logger.info(f'Formal change: {before} -> {after}')
        return Chem.MolToSmiles(mol)
    except Exception as err:
        tb = traceback.format_exc()
        logger.error(f"filter_smiles_by_neutralize_atoms error: {tb}")
    return ''


def smiles_to_inchikey(smiles: str):
    try:
        molfile = Chem.MolFromSmiles(smiles)
        scheme = registration.MolDocScheme()
        moldoc = scheme.generate_mol_doc(molfile)
        return moldoc['index']
    except Exception as err:
        tb = traceback.format_exc()
        logger.error(f"smiles_to_inchikey error: {tb}")
    return ''
