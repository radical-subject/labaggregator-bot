from operator import itemgetter

from rdkit import RDLogger, Chem
from mongordkit.Database import write
from mongordkit.Search import similarity, substructure
from mongordkit.Database import registration

from rdkit.Chem import PandasTools, SanitizeMol

from modules.db.dbconfig import db_client, MONGO_VENDORBOT_DATABASE, MOLECULES_DATABASE
import logging
logger = logging.getLogger(__name__)


class UniqueMolecules:

    def __init__(self, client, db):
        self.name = 'unique_molecules_collection'
        self.collection = client[db][self.name]
        self.client = client
        self.db = client[db]

    def get_molecule(self, index: str):
        return self.collection.find_one({"index": index})

    def add_molecule(self, data):
        result = self.collection.insert_one(data)
        #assert result.modified_count == 1

    # def get_reagents(self, user_id: int):
    #     user = self.get_user(user_id)
    #     if "user_reagents" in user:
    #         return user["user_reagents"]
    #     return []

    # def get_all_users(self):
    #     return list(self.collection.find({}))

    def update_molecule(self, index: int, moldoc):
        logger.info(f"molecule entry is being updated")
        query = {'index': index}
        # if '_id' in user_data:
        #     del user_data['_id']  # Performing an update on the path '_id' would modify the immutable field '_id'
        insertion_result = self.collection.update_one(query, {"$set": moldoc}, upsert=True) 
        logger.info(f"molecule entry updated: {insertion_result.acknowledged}")
        #return result
        #assert result.modified_count == 1

    # def get_users_by_cas(self, cas: str):
    #     return list(self.collection.find({"user_reagents": {'$elemMatch': {'CAS': cas}}}))

    # def get_users_by_smiles(self, smiles: str):
    #     return list(self.collection.find({"user_reagents": {'$elemMatch': {'SMILES': smiles}}}))
    
    
    def reagent_registration(self, SMILES):
        """
        регистрация уникальной молекулы в коллекции уникальных молекул в отдельной базе. 
        набивка ссылками записи уникальной молекулы на конкретные айдишники банок с реагентами
        """
        
        molfile = Chem.MolFromSmiles(SMILES)
        scheme = registration.MolDocScheme()
        # scheme.add_value_field('reagent_internal_id_list', [reagent_internal_id])
        moldoc = scheme.generate_mol_doc(molfile)

        result = write.WriteFromMolList(self.collection, [molfile], scheme=scheme) 
        print (result)
        # print(query)
        # if result == 0:

        #     reagent_internal_id_list = self.get_molecule(moldoc['index'])["value_data"]['reagent_internal_id_list']
            
        #     # if reagent_internal_id not in reagent_internal_id_list:
        #     #     reagent_internal_id_list.append(reagent_internal_id)

        #         # scheme.add_value_field('reagent_internal_id_list', reagent_internal_id_list)

        #         # moldoc = scheme.generate_mol_doc(molfile)
        #     self.update_molecule(moldoc['index'], moldoc)

        inchikey_standard = moldoc['index']
        return inchikey_standard
    

    def calculate_hashes(self):
        """
        необходимая подготовка к поиску, обсчет структурных данных
        """
        # во избежание дубликатов
        db_client[MOLECULES_DATABASE].mfp_counts.drop()
        db_client[MOLECULES_DATABASE].permutations.drop()

        # Search.PrepareForSearch(rdkit_db, rdkit_db.molecules, rdkit_db.mfp_counts, rdkit_db.permutations)
        substructure.AddPatternFingerprints(self.collection)
        similarity.AddMorganFingerprints(self.collection, self.db.mfp_counts) # db_client[MOLECULES_DATABASE]

        # Generate 100 different permutations of length 2048 and save them in demo_db.permutations as separate documents.
        similarity.AddRandPermutations(self.db.permutations)

        # Add locality-sensitive hash values to each documents in demo_db.molecules by splitting the 100 different permutations
        # in demo_db.permutations into 25 different buckets.
        similarity.AddLocalityHashes(self.collection, self.db.permutations, 25)

        # Create 25 different collections in db_demo each store a subset of hash values for molecules in demo_db.molecules.
        similarity.AddHashCollections(self.db, self.collection)


    def similarity_search(self, database, smiles: str):
        """
        Ищем похожие реагенты умными функциями
        :param smiles:
        :return: отсортированный по похожести список реагентов
        """
        smiles = smiles.replace("|", "")  # вертикальная черта в SMILES - непонятно что несёт, и RDKIT ее не понимает, убираем ее
        mol = Chem.MolFromSmiles(smiles)

        if not mol:
            raise Exception("MolFromSmiles returned None")

        res = similarity.SimSearchAggregate(mol, database.unique_molecules_collection, database.mfp_counts, 0.1) # 0.5 = similarity threshold

        if not res:
            return

        res = sorted(res, key=itemgetter(0), reverse=True)

        if res[0] is None:
            raise Exception(f"incorrect result: {str(res)}")

        return res

    def get_most_similar_reagent(self, REQUESTED_SMILES: str):
        if self.similarity_search(self.db, REQUESTED_SMILES) != None:
            best_similarity_result_probability = self.similarity_search(self.db, REQUESTED_SMILES)[0][0]

            best_similarity_result_id = self.similarity_search(self.db, REQUESTED_SMILES)[0][1]

            inchi_key = self.get_molecule(best_similarity_result_id)["index"]
            best_match_smiles = self.get_molecule(best_similarity_result_id)["smiles"]
            return (inchi_key, best_match_smiles, best_similarity_result_probability)
        
        else: 
            return None
        
unique_molecules_collection = UniqueMolecules(db_client, MOLECULES_DATABASE)
