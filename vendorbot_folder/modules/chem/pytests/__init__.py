

class PubChempyComponent:
    def __init__(self, isomeric_smiles):
        self.isomeric_smiles = isomeric_smiles


def pubchempy_smiles_return(ret):
    return [PubChempyComponent(ret)]
