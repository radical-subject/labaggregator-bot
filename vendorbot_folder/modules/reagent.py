
class Reagent:
    """
    Реагент = пробирка. Список Reagent хранится в БД пользователей.
    """
    reagent_id: str = ""         # Уникальный в БД
    inchikey_standard: str = ""  # Молекула вещества. Может быть 2 reagent с одинаковыми inchikey, разными location.

    # поля для поиска
    name: str = ""        # Название Реагента из XLSX пользователя. Может быть любым. Не обязательно.
    cas: str = ""         # CAS реагента из XLSX. Реагенты без CAS не добавляются. Обязательно.
                          # TODO: Сделать CAS не обязательным. могут добавится, если заполнено NAME
    smiles: str = ""      # из XLSX пользователя
    filter_smiles: str = ""  # нейтрализованный SMILES. вдруг по ним тоже будут искать.
    #

    contact: str = ""     # если указан в 1й строке. Если нет, то пустое поле.
    location: str = ""    # из excel столбца
    comment: str = ""     # из excel столбца comment

    sharing_status: str = ""  # какие есть кроме shared? TODO убрать, если не нужен.
    timestamp: str = ""   # DATETIME_FMT TODO убрать, если не нужен.

    user_id:  int = 0     # заполняется только когда ищем реагент по базе

    def __init__(self, cas="", smiles=""):
        self.cas = cas
        self.smiles = smiles

    def __eq__(self, other):
        """
        TODO: кажется нужно переписать или дополнить. Может кстати только reagent_id сравнивать и всё.
        Чтобы можно было сравнивать Reagent ы
        :param other:
        :return:
        """
        if isinstance(other, Reagent):
            return self.reagent_id == other.reagent_id
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
            "SMILES_filtered": self.filter_smiles,
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
