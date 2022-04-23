
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
#ret = cirpy.resolve('120-46-7', 'names')

#pp.pprint(ret)

#ret = pubchempy_smiles_resolve('O=C(CC(=O)c1ccccc1)c2ccccc2')
#

def cas_to_name(cas):
    smiles = cirpy.resolve(cas, 'smiles')
    pp.pprint(smiles)
    ret = cirpy.resolve(smiles, 'names')
    pp.pprint(ret)

#cas_to_name('75-64-9')

ret = cirpy.resolve('Dibenzoylmethane', 'cas')
pp.pprint(ret)

#smiles = pubchempy_smiles_resolve('Dibenzoylmethane')

#pp.pprint(smiles)
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