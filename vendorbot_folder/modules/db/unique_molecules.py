from typing import List
from operator import itemgetter

try:
    from rdkit import RDLogger, Chem
    from rdkit.Chem import PandasTools, SanitizeMol
    from mongordkit.Database import write
    from mongordkit.Search import similarity, substructure
    from mongordkit.Database import registration
except:
    pass

from modules.reagent import Reagent
from modules.db.dbconfig import db_client, MOLECULES_DATABASE
import logging
logger = logging.getLogger(__name__)


class UniqueMolecules:

    def __init__(self, client, db):
        self.client = client
        self.db = client[db]
        self.unique_molecules_collection = client[db]['unique_molecules_collection']
        self.mfp_counts = client[db]['mfp_counts']
        self.permutations = client[db]['permutations']

    def get_molecule(self, index: str):
        return self.unique_molecules_collection.find_one({"index": index})

    def add_molecule(self, data):
        result = self.unique_molecules_collection.insert_one(data)
        if not result.acknowledged:
            logger.error(f"add_molecule error: {result}")

    # def get_reagents(self, user_id: int):
    #     user = self.get_user(user_id)
    #     if "user_reagents" in user:
    #         return user["user_reagents"]
    #     return []

    # def get_all_users(self):
    #     return list(self.unique_molecules_collection.find({}))

    def update_molecule(self, index: int, moldoc):
        logger.info(f"molecule entry is being updated")
        query = {'index': index}
        # if '_id' in user_data:
        #     del user_data['_id']  # Performing an update on the path '_id' would modify the immutable field '_id'
        result = self.unique_molecules_collection.update_one(query, {"$set": moldoc}, upsert=True)
        if not result.acknowledged:
            logger.error(f"update_molecule error: {result}")

    # def get_users_by_cas(self, cas: str):
    #     return list(self.unique_molecules_collection.find({"user_reagents": {'$elemMatch': {'CAS': cas}}}))

    # def get_users_by_smiles(self, smiles: str):
    #     return list(self.unique_molecules_collection.find({"user_reagents": {'$elemMatch': {'SMILES': smiles}}}))

    def register_reagents(self, reagents: List[Reagent]) -> None:
        """
        регистрируем реагент в базе молекул, в коллекции уникальных молекул.
        если молекула уже зарегистрирована, то reagent_internal_id реагента из user_reagents

        TODO: А где поверяется что молекула уже зарегистиррована?

        вписывается в лист value_data в уникальной записи уникальной молекулы.

        TODO: а что за лист value_data?

        таким образом, пропуская дубликаты, создаются ссылки из сущности "банка реактивов"
        на уникальную молекулу, которая содержится в банке и обратно.

        функция регистрации возвращает уникальный индекс - inchikey_standard,
        который добавляется в запись о реагенте у пользователя.
        """
        try:
            for r in reagents:
                self.reagent_registration(r.filter_smiles)  # TODO может нефильтованную smiles ?
        except Exception as err:
            logger.error(err)

    def reagent_registration(self, smiles: str):
        """
        регистрация уникальной молекулы в коллекции уникальных молекул в отдельной базе. 
        набивка ссылками записи уникальной молекулы на конкретные айдишники банок с реагентами
        """
        molfile = Chem.MolFromSmiles(smiles)
        scheme = registration.MolDocScheme()
        # scheme.add_value_field('reagent_internal_id_list', [reagent_internal_id])
        moldoc = scheme.generate_mol_doc(molfile)

        result = write.WriteFromMolList(self.unique_molecules_collection, [molfile], scheme=scheme) 
        logger.info(f"reagent_registration (smiles={smiles}) result: {result}")
        # print(query)
        # if result == 0:
        #     reagent_internal_id_list = self.get_molecule(moldoc['index'])["value_data"]['reagent_internal_id_list']
        #     # if reagent_internal_id not in reagent_internal_id_list:
        #     #     reagent_internal_id_list.append(reagent_internal_id)
        #         # scheme.add_value_field('reagent_internal_id_list', reagent_internal_id_list)
        #         # moldoc = scheme.generate_mol_doc(molfile)
        #     self.update_molecule(moldoc['index'], moldoc)

    def calculate_hashes(self):
        """
        необходимая подготовка к поиску, обсчет структурных данных
        """
        # во избежание дубликатов
        db_client[MOLECULES_DATABASE].mfp_counts.drop()
        db_client[MOLECULES_DATABASE].permutations.drop()

        # Search.PrepareForSearch(rdkit_db, rdkit_db.molecules, rdkit_db.mfp_counts, rdkit_db.permutations)
        substructure.AddPatternFingerprints(self.unique_molecules_collection)
        similarity.AddMorganFingerprints(self.unique_molecules_collection, self.mfp_counts) # db_client[MOLECULES_DATABASE]

        # Generate 100 different permutations of length 2048 and save them in demo_db.permutations as separate documents.
        similarity.AddRandPermutations(self.permutations)

        # Add locality-sensitive hash values to each documents in demo_db.molecules by splitting the 100 different permutations
        # in demo_db.permutations into 25 different buckets.
        similarity.AddLocalityHashes(self.unique_molecules_collection, self.permutations, 25)

        # Create 25 different collections in db_demo each store a subset of hash values for molecules in demo_db.molecules.
        similarity.AddHashCollections(self.db, self.unique_molecules_collection)


    def similarity_search(self, smiles: str):
        """
        Ищем похожие реагенты умными функциями
        :param smiles:
        :return: отсортированный по похожести список реагентов
        """
        smiles = smiles.replace("|", "")  # вертикальная черта в SMILES - непонятно что несёт, и RDKIT ее не понимает, убираем ее
        mol = Chem.MolFromSmiles(smiles)

        if not mol:
            raise Exception("MolFromSmiles returned None")

        # 0.1 or 0.5 = similarity threshold
        res = similarity.SimSearchAggregate(mol, self.unique_molecules_collection, self.mfp_counts, 0.1)

        if not res: 
            return

        res = sorted(res, key=itemgetter(0), reverse=True)

        if res[0] is None:
            raise Exception(f"incorrect result: {str(res)}")

        return res

    def get_most_similar_reagent(self, REQUESTED_SMILES: str):
        searchResult = self.similarity_search(REQUESTED_SMILES)
        if searchResult != None:

            best_similarity_result_probability = searchResult[0][0]
            best_similarity_result_id = searchResult[0][1]

            molecule = self.get_molecule(best_similarity_result_id)
            inchi_key = molecule["index"]
            best_match_smiles = molecule["smiles"]

            return (inchi_key, best_match_smiles, best_similarity_result_probability)
        
        else: 
            return None

    def get_5_most_similar_reagents(self, REQUESTED_SMILES: str):

        """
        returns list of tuples with max length of 5
        [(inchi_key, best_match_smiles, best_similarity_result_probability), (), ()..]

        """

        searchResult = self.similarity_search(REQUESTED_SMILES)
        if searchResult != None:
            number_of_hits = min(len(searchResult), 5)

            return [(self.get_molecule(searchResult[i][1])["index"], self.get_molecule(searchResult[i][1])["smiles"], searchResult[i][0]) for i in range(number_of_hits)]
        
        else: 
            return None


unique_molecules_collection = UniqueMolecules(db_client, MOLECULES_DATABASE)
