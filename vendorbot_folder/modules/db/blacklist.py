from typing import Optional, List

try:
    from rdkit import RDLogger, Chem
    from mongordkit.Database import write
    from mongordkit.Search import similarity, substructure

    from rdkit.Chem import PandasTools, SanitizeMol
except:
    pass  # for debugging bot without import rdkit

import logging

from operator import itemgetter
from modules.db.dbconfig import root_client

logger = logging.getLogger(__name__)


class BlackList:

    def __init__(self, client, db):
        self.client = client
        self.db = db

        self.molecules = client[db].molecules
        self.mfp_counts = client[db].mfp_counts
        self.permutations = client[db].permutations

        self.pandas = client[db].pandas

    def drop_database(self):
        logger.info('BlackList: drop database')
        self.client.drop_database(self.db)

    def reload_rdkit(self, path_sdf: str = './srs/Narkotiki_test.sdf'):
        """
        :param path_sdf:
        :return:
        """
        logger.info('reload_rdkit')

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

        return result

    def reload_pandas(self, path_sdf: str = './srs/Narkotiki_test.sdf'):
        molecules = PandasTools.LoadSDF(path_sdf)

        molecules_dict = molecules.to_dict()

        logger.info(f"reload_pandas: found {len(molecules_dict['ID'].keys())} molecules")
        data = []

        for ID in molecules_dict['ID'].keys():
            SanitizeMol(molecules_dict["ROMol"][ID])
            data.append({
                "ID": ID,
                "references": molecules_dict['references'][ID],
                "UPAC_name": molecules_dict['UPAC_name'][ID],
                "CAS": molecules_dict['CAS'][ID],
                "List": molecules_dict['List'][ID],
                "category": molecules_dict['category'][ID],
                "NameRUS": molecules_dict['NameRUS'][ID],
                "NameTRIVIAL": molecules_dict['NameTRIVIAL'][ID],
                "Synonym1": molecules_dict['Synonym1'][ID],
                "precursor": molecules_dict['precursor'][ID],
                "Primechanie": molecules_dict['Primechanie'][ID],
                "Comment": molecules_dict['Comment'][ID],
                "DateADD": molecules_dict['DateADD'][ID],
                "NumberPOST": molecules_dict['NumberPOST'][ID],
                "SMILES": Chem.rdmolfiles.MolToSmiles(molecules_dict["ROMol"][ID])
            })

        self.pandas.insert_many(data)

        return True

    def similarity_search(self, smiles: str) -> Optional[List[List]]:
        """
        Ищем похожие реагенты умными функциями
        :param smiles:
        :return: отсортированный по похожести список реагентов
        """
        # вертикальная черта в SMILES - непонятно что несёт,
        # и RDKIT ее не понимает, убираем ее
        smiles = smiles.replace("|", "")
        mol = Chem.MolFromSmiles(smiles)

        if not mol:
            raise Exception("MolFromSmiles returned None")

        res = similarity.SimSearchAggregate(mol, self.molecules, self.mfp_counts, 0.1)

        logger.debug(f"{smiles} similarity ret = {len(res)}")

        if not res:
            return False  # TODO что мы возвращаем?

        res = sorted(res, key=itemgetter(0), reverse=True)

        if res[0] is None:
            raise Exception(f"incorrect result: {str(res)}")

        return res

    def is_similar(self, smiles: str, threshold: float = 0.75):
        """
        :param smiles: SMILES
        :param threshold: порог похожести
        :return: Есть ли подходящий на threshold реагент в blacklist?
        """
        res = self.similarity_search(smiles)
        if not res:
            return False

        logger.info(f"is_similar: {smiles} = {res[0][0]}")
        return res[0][0] > threshold

    def find_pandas_data(self, smiles: str):
        return self.pandas.find({"SMILES": smiles})

    def find_similar_pandas_data(self, smiles: str):
        """
        :param smiles: SMILES
        :return: наиболее подходящий реагент из blacklist_pandas или информацию он нём
        """
        res = self.similarity_search(smiles)
        if res:
            result = self.molecules.find({"index": f"{res[0][1]}"})
            if result:
                return self.find_pandas_data(result["smiles"])


blacklist_engine = BlackList(root_client, 'blacklist_rdkit_db')
