
class Reagent:
    cas: str = ''
    smiles: str = ''
    filter_smiles: str = ''   # для регистрации в unique_molecules DB

    reagent_id: str = ''      # почему не id ?
    inchikey_standard: str = ''
    sharing_status: str = ''  # какие есть кроме shared?
    timestamp: str = ''   # DATETIME_FMT

    contact: str = ''     # у кого, если было заполнено в csv
    location: str = ''    # где лежит, если было заполнено в csv

    user_id:  int = 0     # заполняется только когда ищем реагент по базе

    name: str = ''  # TODO не понятно откуда взялся, а по нему уже поиск есть ))

    def __init__(self, cas='', smiles=''):
        self.cas = cas
        self.smiles = smiles

    def __str__(self):
        if self.smiles:
            return f"(CAS:{self.cas}, SMILES:{self.smiles})"
        else:
            return f"(CAS:{self.cas})"

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
            "contact": self.contact,
            "location": self.location,
            "name": self.name
        }

    def from_dict(self, d):
        self.reagent_id = d.get("reagent_id", "")
        self.inchikey_standard = d.get("inchikey_standard", "")
        self.cas = d.get("CAS", "")
        self.smiles = d.get("SMILES", "")
        self.sharing_status = d.get("sharing_status", "")
        self.timestamp = d.get("timestamp", "")
        self.contact = d.get("contact", "")
        self.location = d.get("location", "")
        self.name = d.get("name", "")


REAGENT_SHARED = "shared"
DATETIME_FMT = "%d.%m.%Y %H:%M:%S"
