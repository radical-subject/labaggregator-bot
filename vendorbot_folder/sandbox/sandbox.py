
# место для мусора и игрищ


from typing import List, Tuple


def banch_cas_to_smiles(cas_list: List[str]):
    from multiprocessing import Pool
    from .cas_to_smiles import cas_to_smiles

    n = 10
    with Pool(processes=n) as pool:
        def get_smiles(cas: str):
            try:
                smiles = cas_to_smiles(cas)
                return cas, smiles
            except Exception as err:
                pass

        return pool.map(get_smiles, cas_list)


def demo_insert_cas_to_user_functions(CAS_list):
    # импорт листа реагентов с фильтрациями
    from modules.ourbot.service.helpers import is_CAS_number, Reagent

    valid_cas_list = [r for r in CAS_list if is_CAS_number(r)]
    cas_smiles_list = banch_cas_to_smiles(valid_cas_list)

    whitelist = [cas_smile for cas_smile in cas_smiles_list if blacklist.is_similar(cas_smile[1])]

    user_reagents_object.replace_reagents(contact, whitelist)

    return {
        "input_lines_number": len(CAS_list),
        "valid_CAS_numbers": len(valid_cas_list),
        "failed_CAS_check_number": len(CAS_list)-len(valid_cas_list),
        "SMILES_not_found": len(valid_cas_list) - len(cas_smiles_list),
        "SMILES_found": len(cas_smiles_list),
        "blacklist_filter_result": len(cas_smiles_list) - len(whitelist),
        "total_reagents_imported": len(whitelist),
        #"total_reagents_count_in_DB": len(self.user_reagents)
    }


def replace_reagents(contact: str, cas_smiles: Tuple[str, str]):
    # эта функция у юзера
    reagents = []

    now = time.strftime("%d.%m.%Y %H:%M", time.localtime())

    for cas, smiles in cas_smiles:
        reagents.append({
            "reagent_internal_id": uuid.uuid4().hex,
            "CAS": cas,
            "contact": contact,
            "sharing_status": "shared",
            "timestamp": now
        })


class BlackList:

    def __init__(self, client, db):  # blacklist_rdkit_db
        self.client = client
        self.db = db

        self.molecules = client[db].molecules
        self.mfp_counts = client[db].mfp_counts
        self.permutations = client[db].permutations

    def drop_database(self):
        self.client.drop_database(self.db)

    def update_blacklist(self, path_sdf: str = './srs/Narkotiki_test.sdf'):
        from rdkit import RDLogger
        from mongordkit.Database import write
        from mongordkit.Search import similarity, substructure

        # Disable rdkit warnings
        RDLogger.DisableLog('rdApp.*')

        # clean previous version of blacklist_rdkit_db
        self.drop_database()

        result = write.WriteFromSDF(self.molecules, path_sdf)

        # Search.PrepareForSearch(rdkit_db, rdkit_db.molecules, rdkit_db.mfp_counts, rdkit_db.permutations)

        substructure.AddPatternFingerprints(self.molecules)
        similarity.AddMorganFingerprints(self.molecules, self.mfp_counts)

        # Generate 100 different permutations of length 2048 and save them in demo_db.permutations as separate documents.
        similarity.AddRandPermutations(self.permutations)

        # Add locality-sensitive hash values to each documents in demo_db.molecules by splitting the 100 different permutations
        # in demo_db.permutations into 25 different buckets.
        similarity.AddLocalityHashes(self.molecules, self.permutations, 25)

        # Create 25 different collections in db_demo each store a subset of hash values for molecules in demo_db.molecules.
        similarity.AddHashCollections(self.client[self.db], self.molecules)

    def is_similar(self, smiles: str, threshold: float = 0.75):
        from rdkit import Chem
        from mongordkit.Search import similarity
        from operator import itemgetter

        smiles = smiles.replace("|", "") # вертикальная черта в SMILES - непонятно что несёт, и RDKIT ее не понимает, убираем ее
        mol = Chem.MolFromSmiles(smiles)

        res = similarity.SimSearchAggregate(mol, self.molecules, self.mfp_counts, 0.1)
        res = sorted(res, key=itemgetter(0), reverse=True)

        return res[0][0] > threshold

#blacklist = BlackList(db_client, 'blacklist_rdkit_db')