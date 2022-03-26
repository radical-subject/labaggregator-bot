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


def CIRPY_resolve(str_input):
    try:
        res = cirpy.resolve(str_input, 'smiles')
        (smiles, failed_CAS) = (res, "NaN")
        if res == None:
            try:
                pubchem_response = pubchempy.get_compounds(str_input, "name")
                (smiles, failed_CAS) = (pubchem_response[0].isomeric_smiles, "NaN")
            except:
                (smiles, failed_CAS) = ("resolver_error", str_input)
    except:
        try:
            pubchem_response = pubchempy.get_compounds(str_input, "name")
            res = pubchem_response[0].isomeric_smiles
            (smiles, failed_CAS) = (res, "NaN")
        except:
            print(str_input)
            (smiles, failed_CAS) = ("resolver_error", str_input)

    return (smiles, failed_CAS)

def batch_SMILES_resolve(input_txt_file_path):
    import_CAS_df = pd.read_csv(input_txt_file_path, header = None)
    CAS_list = import_CAS_df[0].tolist()
    
    timer = Timer()
    timer.start()
    n = 50
    with Pool(processes=n) as pool:
        result = pool.map(CIRPY_resolve, CAS_list)
    timer.stop()

    return result



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
