from multiprocessing import Pool
import pandas as pd
import cirpy, pubchempy
import re
from modules.ourbot.service.timer import Timer


def get_SMILES(request_query):
    try:
        pubchem_response = pubchempy.get_compounds(request_query, "name")
        return pubchem_response[0].isomeric_smiles
    except:
        return cirpy.resolve("{}".format(request_query), 'smiles')


def CIRPY_resolve(reagent_without_SMILES):
    """
    {
        "reagent_internal_id": uuid.uuid4().hex, 
        "CAS": CAS_number
    }
    """
    print(reagent_without_SMILES)
    try:
        res = cirpy.resolve(reagent_without_SMILES['CAS'], 'smiles')
        reagent_without_SMILES["SMILES"] = res
        if res == None:
            try:
                pubchem_response = pubchempy.get_compounds(reagent_without_SMILES['CAS'], "name")
                res = pubchem_response[0].isomeric_smiles
                reagent_without_SMILES["SMILES"] = res
            except:
                reagent_without_SMILES["SMILES"] = "resolver_error"
    except:
        try:
            pubchem_response = pubchempy.get_compounds(reagent_without_SMILES['CAS'], "name")
            res = pubchem_response[0].isomeric_smiles
            reagent_without_SMILES["SMILES"] = res
        except:
            reagent_without_SMILES["SMILES"] = "resolver_error"

    return reagent_without_SMILES


def batch_SMILES_resolve(reagents_without_SMILES_list):
    # import_CAS_df = pd.read_csv(input_txt_file_path, header = None)
    # CAS_list = import_CAS_df[0].tolist()
    
    timer = Timer()
    timer.start()
    n = 50
    with Pool(processes=n) as pool:
        result = pool.map(CIRPY_resolve, reagents_without_SMILES_list)
    timer.stop()

    result_object_list = result # if i[0]!="resolver_error"
    
    errors_list = [i for i in result if i["SMILES"]=="resolver_error"]

    # # remove all error indications - this gives clean SMILES list
    # SMILES_list = list(filter(("resolver_error").__ne__, SMILES_list))
    
    return (result_object_list, errors_list)



def get_IUPAC(request_query):

    try:
        pubchem_response = pubchempy.get_compounds(request_query, 'name')
        return pubchem_response[0].iupac_name
    except:
        IndexError
        try:
            return cirpy.resolve("{}".format(request_query), 'iupac_name') #This is alternative, but it is bad due to the fact that list of synonyms in cirpy is not ranked by quality
        except:
            return None


def get_CAS(request_query):
    #PATTERN FOR CAS MATCHING
    pattern = re.compile("^\d+-\d+-\d+$")
    i = 0
    cas_list = []

    cirpy_response = cirpy.resolve("{}".format(request_query), 'cas')
    case = type(cirpy_response)

    if "{}".format(case) != "<class 'list'>":
        cas_list.append(cirpy_response)
    else:
        cas_list = cirpy_response
    try:
        while True:
            if pattern.match(cas_list[i]):
                return cas_list[i]
            i += 1
    except:
        Exception
        if pattern.match(request_query) != False:
            try:
                pubchem_response = pubchempy.get_compounds(request_query, 'name')
                pubchem_response[0].iupac_name ## EXAMPLE
                return None
            except:
                return None
        return None



def get_SYNONYMS(request_query):
        pubchem_response = pubchempy.get_compounds(request_query, 'name')
        try:
            return pubchem_response[0].synonyms
        except:
            return None
