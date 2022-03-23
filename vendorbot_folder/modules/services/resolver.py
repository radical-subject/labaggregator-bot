import cirpy, pubchempy
import re


def get_SMILES(request_query):
    try:
        pubchem_response = pubchempy.get_compounds(request_query, "name")
        return pubchem_response[0].isomeric_smiles
    except:
        return cirpy.resolve("{}".format(request_query), 'smiles')

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
