import logging
from modules.db import dbschema
from modules.ourbot.service.decorators import log_errors
logger = logging.getLogger(__name__)


@log_errors
def purge(client, db_instance):
    db_name = db_instance.DATABASE_NAME
    db = client[db_name]
    db.command("dropDatabase")
    logger.info(f"database {db_name} dropped.")

# @log_errors
def add_records(client, db_instance, collection_name: str, data: dict):
    db_name = db_instance.DATABASE_NAME
    # logger.info(f"{client}, {db_instance}, {collection_name}, {data}")
    collection = client[db_name][collection_name]
    result = collection.insert_one(data)
    logger.info(f"data inserted into {db_name}, {collection_name}")
    return result

@log_errors
def update_record(client, db_instance, collection_name: str, query: dict, data: dict):
    db_name = db_instance.DATABASE_NAME
    collection = client[db_name][collection_name]
    # result = collection.insert_one(data)
    result = collection.update(query, {"$set":data}, upsert=True)
    logger.info(f"data upserted into {db_name}, {collection_name}")
    return result

@log_errors
def get_records(client, db_instance, collection_name: str, query: dict, *args):
    db_name = db_instance.DATABASE_NAME
    collection = client[db_name][collection_name]
    search = list(collection.find(query, *args))
    # logger.info(search)
    return search


@log_errors
def get_timerdata_object(client, db_instance, collection_name: str, query: dict, user_id):
        
    # достаем ее из бд
    previous_records=get_records(client, db_instance, collection_name, query)
    
    # результат поиска может оказаться пустым
    if previous_records == [] or previous_records == None:
        timer_object = dbschema.TimerData(
            **{
                "user_id": user_id
            }
        )
        
    else:
        # если раньше у пользователя были записи то импортируем данные пользователя в объект таймера
        timer_object = dbschema.TimerData(
            **previous_records[0]
        )
    
    return timer_object



    # --------------------------------------------------
    # BLACKLIST REBORN
    # ---------------------------------------------------

#RDKIT
import rdkit
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import rdmolfiles, AllChem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
# for similarity search need mongo-rdkit
from mongordkit.Search import similarity, substructure, utils
from mongordkit import Search
from mongordkit.Database import create, write


@log_errors    
def update_rdkit_db_blacklist ():
    # echo $PYTHONPATH
    # home_path = "/home/oikura/github/reagent_checker_bot"
    # os.chdir(home_path)
    # sys.path.append(home_path)
    # os.system('cd ./mongo-rdkit && export PYTHONPATH="$PWD\n')
    # Disable rdkit warnings
    rdkit.RDLogger.DisableLog('rdApp.*')

    # clean previous version of rdkit_db
    connection.drop_database('rdkit_db')     
    # mfp_counts.drop()
    # permutations.drop()
    # molecules.drop()

    result = write.WriteFromSDF(rdkit_db.molecules, './srs/Narkotiki_test.sdf')
    # Search.PrepareForSearch(rdkit_db, rdkit_db.molecules, rdkit_db.mfp_counts, rdkit_db.permutations)

    substructure.AddPatternFingerprints(rdkit_db.molecules)
    similarity.AddMorganFingerprints(rdkit_db.molecules, rdkit_db.mfp_counts)

    # Generate 100 different permutations of length 2048 and save them in demo_db.permutations as separate documents.
    similarity.AddRandPermutations(rdkit_db.permutations)

    # Add locality-sensitive hash values to each documents in demo_db.molecules by splitting the 100 different permutations
    # in demo_db.permutations into 25 different buckets. 
    similarity.AddLocalityHashes(rdkit_db.molecules, rdkit_db.permutations, 25)

    # Create 25 different collections in db_demo each store a subset of hash values for molecules in demo_db.molecules.
    similarity.AddHashCollections(rdkit_db, rdkit_db.molecules)

    return result


def similarity_search(SMILES_input):
    # echo $PYTHONPATH

    # home_path = "/home/oikura/github/reagent_checker_bot"
    # os.chdir(home_path)
    # sys.path.append(home_path)
    # os.system('cd ./mongo-rdkit && export PYTHONPATH="$PWD')
    mol_input = Chem.MolFromSmiles(SMILES_input)
   # results_similarity = similarity.SimSearch(mol_input, rdkit_db.molecules, rdkit_db.mfp_counts, 0.1)
    results_similarity = similarity.SimSearchAggregate(mol_input, rdkit_db.molecules, rdkit_db.mfp_counts, 0.1)
    # results_substructure = substructure.SubSearch(mol_input, rdkit_db.molecules, chirality=False)
    from operator import itemgetter
    results_similarity = sorted(results_similarity, key=itemgetter(0), reverse=True)

    return results_similarity



def update_blacklist_with_pandas ():
    SDFFile = "./srs/Narkotiki_test.sdf"
    molecules = PandasTools.LoadSDF(SDFFile)
    molecules_dict = molecules.to_dict()
    blacklist_data = []
    for ID in molecules_dict['ID'].keys():
        rdkit.Chem.SanitizeMol(molecules_dict["ROMol"][ID])
        blacklist_data.append({ "ID" : ID,
                                "references" : molecules_dict['references'][ID],
                                "UPAC_name" : molecules_dict['UPAC_name'][ID],
                                "CAS" : molecules_dict['CAS'][ID],
                                "List" : molecules_dict['List'][ID],
                                "category" : molecules_dict['category'][ID],
                                "NameRUS" : molecules_dict['NameRUS'][ID],
                                "NameTRIVIAL" : molecules_dict['NameTRIVIAL'][ID],
                                "Synonym1" : molecules_dict['Synonym1'][ID],
                                "precursor" : molecules_dict['precursor'][ID],
                                "Primechanie" : molecules_dict['Primechanie'][ID],
                                "Comment" : molecules_dict['Comment'][ID],
                                "DateADD" : molecules_dict['DateADD'][ID],
                                "NumberPOST" : molecules_dict['NumberPOST'][ID],
                                "SMILES" : Chem.rdmolfiles.MolToSmiles(molecules_dict["ROMol"][ID])
                              }) 

    blacklist_pandas.insert_many(blacklist_data)
    return True

def convert_to_smiles_and_get_additional_data(results_similarity):
    myquery = { "index": "{}".format(results_similarity[0][1]) }
    print (myquery)
    search_result = molecules.find(myquery)
    try:
        for x in search_result:
            SMILES = x["smiles"]
            print (SMILES)
    except:
        return False
    myquery = { "SMILES": SMILES }
    print(myquery)
    search_result = blacklist_pandas.find(myquery)
    try:
        for x in search_result:
            return x
    except:
        return False