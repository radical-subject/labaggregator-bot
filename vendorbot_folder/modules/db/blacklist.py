try:
    from rdkit import RDLogger, Chem
    from mongordkit.Database import write
    from mongordkit.Search import similarity, substructure
except:
    pass  # for debugging bot without import rdkit

import logging
logger = logging.getLogger(__name__)

from operator import itemgetter
from modules.db.dbconfig import db_client


class BlackList:

    def __init__(self, client, db):
        self.client = client
        self.db = db

        self.molecules = client[db].molecules
        self.mfp_counts = client[db].mfp_counts
        self.permutations = client[db].permutations

    def drop_database(self):
        logger.info('BlackList: drop database')
        self.client.drop_database(self.db)

    def update_blacklist(self, path_sdf: str = './srs/Narkotiki_test.sdf'):

        logger.info('BlackList: update_blacklist')

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

        smiles = smiles.replace("|", "")  # вертикальная черта в SMILES - непонятно что несёт, и RDKIT ее не понимает, убираем ее
        mol = Chem.MolFromSmiles(smiles)

        res = similarity.SimSearchAggregate(mol, self.molecules, self.mfp_counts, 0.1)

        if res is None or res[0] is None:
            raise Exception(f"incorrect result: {str(res)}")

        else:
            res = sorted(res, key=itemgetter(0), reverse=True)
            if res is None or res[0] is None:
                raise Exception(f"incorrect sorted result: {str(res)}")

            logger.info(f"is_similar: {smiles} = {res[0][0]}")

        return res[0][0] > threshold


blacklist_engine = BlackList(db_client, 'blacklist_rdkit_db')