
import os
from modules.ourbot.handlers.decorators import log_errors
import logging
logger = logging.getLogger(__name__)


# ----------------------------------------------------------
# ----------------------------------------------------------
# VENDOR CATALOGS RDKIT FUNCTIONS
# ----------------------------------------------------------
# ----------------------------------------------------------

try:
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import PandasTools
    # for similarity search need mongo-rdkit
    from mongordkit.Search import similarity, substructure
    from mongordkit.Database import write
except:
    logger.error('rdkit not found')


@log_errors
def update_rdkit_with_sialdrich (client, db_instance):
    # Disable rdkit warnings
    rdkit.RDLogger.DisableLog('rdApp.*')

    db_name = db_instance.DATABASE_NAME # sialdrich_rdkit_db
    molecules_collection = client[db_name].molecules
    mfp_counts_collection = client[db_name].mfp_counts
    permutations_collection = client[db_name].permutations

    # clean previous version 
    # DB_rdkit_connection.connection.db.command("dropDatabase")
	# or like this
    # clean previous version of sialdrich_rdkit_db
    client.drop_database(db_name) # sialdrich_rdkit_db

    # смотрим на разбитые на 10000 молекул куски SDF базы данных sialdrich
    parts=os.listdir("./jupyter-scripts/sdf_sial/")
    
    for i in parts:
        result = write.WriteFromSDF(molecules_collection, f'./jupyter-scripts/sdf_sial/{i}')
    # result = write.WriteFromSDF(vendors_DB_rdkit_connection.molecules, f'./jupyter-scripts/sdf_sial/output0rename.sdf')
    # Search.PrepareForSearch(DB_rdkit_connection, DB_rdkit_connection.molecules, DB_rdkit_connection.mfp_counts, DB_rdkit_connection.permutations)

    substructure.AddPatternFingerprints(molecules_collection)
    similarity.AddMorganFingerprints(molecules_collection, mfp_counts_collection)

    # Generate 100 different permutations of length 2048 and save them in demo_db.permutations as separate documents.
    similarity.AddRandPermutations(permutations_collection)

    # Add locality-sensitive hash values to each documents in demo_db.molecules by splitting the 100 different permutations
    # in demo_db.permutations into 25 different buckets. 
    similarity.AddLocalityHashes(molecules_collection, permutations_collection, 25)

    # Create 25 different collections in db_demo each store a subset of hash values for molecules in demo_db.molecules.
    similarity.AddHashCollections(client[db_name], molecules_collection)

    return result


def similarity_search(DB_rdkit_connection, SMILES_input):

    mol_input = Chem.MolFromSmiles(SMILES_input)
    # similarity_results = similarity.SimSearch(mol_input, DB_rdkit_connection.molecules, DB_rdkit_connection.mfp_counts, 0.1)
    similarity_results = similarity.SimSearchAggregate(mol_input, DB_rdkit_connection.molecules, DB_rdkit_connection.mfp_counts, 0.1)
    # results_substructure = substructure.SubSearch(mol_input, DB_rdkit_connection.molecules, chirality=False)
    from operator import itemgetter
    similarity_results = sorted(similarity_results, key=itemgetter(0), reverse=True)
    # print(similarity_results)
    return similarity_results


# ----------------------------------------------------------
# ----------------------------------------------------------
# BLACKLIST RDKIT FUNCTIONS
# ----------------------------------------------------------
# ----------------------------------------------------------


def convert_to_smiles(client, db_instance, similarity_results):
    """
    Конвертирует результат поиска (хеш) вытаскивая структурную информацию о результате в виде SMILES
    #TODO не пойму для чего нужна... find_similar_pandas_data
    """
    db_name = db_instance.DATABASE_NAME # blacklist_rdkit_db
    molecules_collection = client[db_name].molecules

    similarity_search_compound_dict = []

    for i in range(len(similarity_results)):
        if similarity_results[i][0] > 0.70:
            myquery = { "index": "{}".format(similarity_results[i][1]) }
            search_result = molecules_collection.find(myquery)
            for x in search_result:
                SMILES = x["smiles"]
            similarity_search_compound_dict.append({ "SMILES": SMILES })
        else: 
            continue

    return similarity_search_compound_dict


def convert_to_smiles_and_get_additional_data(client, db_instance, similarity_results):
    """
    TODO: Не понимаю смысла этой фнукции
    Переписал find_similar_pandas_data
    """
    db_name = db_instance.DATABASE_NAME # blacklist_rdkit_db
    molecules_collection = client[db_name].molecules
    blacklist_pandas_collection = client[db_name].blacklist_pandas

    myquery = { "index": "{}".format(similarity_results[0][1]) }
    search_result = molecules_collection.find(myquery)
    try:
        for x in search_result:
            SMILES = x["smiles"]
            # print(SMILES)
    except:
        return False

    myquery = { "SMILES": SMILES }

    search_result = blacklist_pandas_collection.find(myquery)
    try:
        for x in search_result:
            # print(f"pandas_blacklist results: {x}")
            return SMILES, x
    except:
        return SMILES