
class Reagent:
    cas: str = ""
    smiles: str = ""
    filter_smiles: str = ""   # для регистрации в unique_molecules DB

    reagent_id: str = ""     # почему не id ?
    inchikey_standard: str = ""
    sharing_status: str = ""  # какие есть кроме shared?
    timestamp: str = ""   # DATETIME_FMT

    contact: str = ""     # у кого есть: из 1й строки csv или самого юзера
    location: str = ""    # из excel столбца comment

    user_id:  int = 0     # заполняется только когда ищем реагент по базе

    name: str = ""        # из excel столбца name
    comment: str = ""     # из excel столбца comment

    def __init__(self, cas="", smiles=""):
        self.cas = cas
        self.smiles = smiles

    def __eq__(self, other):
        """
        Чтобы можно было сравнивать Reagent ы
        :param other:
        :return:
        """
        if isinstance(other, Reagent):
            return self.cas == other.cas\
                   and self.smiles == other.smiles \
                   and self.location == other.location \
                   and self.contact == other.contact
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

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
            "name": self.name,
            "comment": self.comment,
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
        self.comment = d.get("comment", "")


REAGENT_SHARED = "shared"
DATETIME_FMT = "%d.%m.%Y %H:%M:%S"
