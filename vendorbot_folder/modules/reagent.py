
class Reagent:
    cas: str = ''
    smiles: str = ''
    contact: str = ''
    filter_smiles: str = ''   # для регистрации в unique_molecules DB

    reagent_id: str = ''      # почему не id ?
    inchikey_standard: str = ''
    sharing_status: str = ''  # какие есть кроме shared?
    timestamp: str = ''   # DATETIME_FMT

    def __init__(self, cas='', smiles=''):
        self.cas = cas
        self.smiles = smiles

    def __str__(self):
        if self.smiles:
            return f"(cas={self.cas}, smiles={self.smiles})"
        else:
            return f"(cas={self.cas})"

    def __repr__(self):
        return str(self)

    def to_dict(self):
        return {
            "reagent_id": self.reagent_id,
            "inchikey_standard": self.inchikey_standard,
            "CAS": self.cas,
            "SMILES": self.smiles,
            "sharing_status": self.sharing_status,
            "timestamp": self.timestamp,
            "contact": self.contact
        }


REAGENT_SHARED = "shared"
DATETIME_FMT = "%d.%m.%Y %H:%M:%S"
