
# место для мусора и игрищ


from typing import List, Tuple
import requests
import pubchempy
import pprint
import cirpy

pp = pprint.PrettyPrinter(indent=4)

from modules.ourbot.service.cas_to_smiles import pubchempy_smiles_resolve

#
# convert_to_smiles_and_get_additional_data(client, db_instance, result)[1]['NameRUS']


#res = pubchempy_smiles_resolve('Allyl bromide')
#pubchem_response = pubchempy.get_compounds('3-bromoprop-1-ene', "name")
#smiles = pubchem_response[0].isomeric_smiles

#pp.pprint(smiles)

#pp.pprint(pubchem_response[0].to_dict())

#ret = cirpy.resolve(smiles, 'cas')

#pp.pprint(ret)

#ret = cirpy.resolve('Allyl bromide', 'names')

#pp.pprint(ret)

#ret = cirpy.resolve('Allyl bromide', 'iupac_name')

#pp.pprint('iupac_name:' + ret)  # 3-bromoprop-1-ene

#ret = cirpy.resolve('Allyl bromide', 'formula')

#pp.pprint('formula:' + ret) #  C3H5Br

#ret = cirpy.resolve('3-bromoprop-1-ene', 'names')

#pp.pprint(ret)