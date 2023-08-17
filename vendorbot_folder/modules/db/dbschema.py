
from typing import List
import uuid
import time
import logging
import traceback
import pandas as pd

from modules.db.blacklist import blacklist_engine
from modules.chem.cas_to_smiles import is_cas_number
from modules.chem import batch
from modules.db.unique_molecules import *

try:
    from rdkit import Chem
    from rdkit.Chem.rdmolops import GetFormalCharge
except:
    pass

logger = logging.getLogger(__name__)


def get_shared_reagents(user):
    if 'user_reagents' in user:
        return filter(lambda r: r['sharing_status'] == 'shared', user['user_reagents'])
    return []


def reagent_name(r):
    """
    У нас может быть что-то не заполнено (?)
    """
    if 'CAS' in r and r['CAS']:
        return r['CAS']
    if 'reagent_name' in r and r['reagent_name']:
        return r['reagent_name']


def reagent_CAS(r):
    """
    У нас может быть что-то не заполнено (?)
    """
    if 'CAS' in r and r['CAS']:
        return r['CAS']


def get_contact(user):
    if user["username"]:
        return f"@{user['username']}"
    elif user["phone_number"]:
        return user["phone_number"]  # TODO добавить + если 7
    else:
        return "user haven't specified phone or username!"


def reagent_contact(user, reagent):
    if 'contact' in reagent and reagent['contact']:
        return reagent['contact']
    else:
        return get_contact(user)


def find_reagent(user, value):
    reagents = []
    if 'user_reagents' in user:
        for reagent in user['user_reagents']:
            if value in reagent.values():
                reagents.append(reagent)
    return reagents


def get_reagent_contacts(users, text):
    """
    :param users: список объектов пользователей
    :param text: CAS или SMILES
    :return: список контактов
    """
    ret = []
    if users:
        for user in users:
            reagents = find_reagent(user, text)
            for r in reagents:
                contact = reagent_contact(user, r)
                if contact:
                    ret.append(contact)
    return ret


def neutralize_atoms(mol):
    '''
    atom neutralizer
    '''
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol



def charge_fixer(SMILES):
    '''
    Neutralize molecules atom by atom
    '''
    SMILES = SMILES.replace("|", "")
    # Create RDKit molecular objects
    mol = Chem.MolFromSmiles(SMILES)
    if GetFormalCharge(mol) != 0:
        charge=GetFormalCharge(mol)
        print(f"FORMAL CHARGE WAS:{charge}")
        neutralize_atoms(mol)
        charge=GetFormalCharge(mol)
        print(f"FORMAL CHARGE BECAME:{charge}")
    return Chem.MolToSmiles(mol)
    

def parse_cas_list(cas_file: pd.DataFrame, contact: str = ''):
    """
    Фильтруем список CAS, ищем SMILES, удаляем прекурсоры, возвращаем список компонентов для БД и статистику
    :param cas_file:
    :param contact:
    :return:
    """
    condition = cas_file['CAS'].fillna('Missing').map(lambda x: is_cas_number(x))
    valid_cas_list = cas_file.loc[condition]
    failed_cas = cas_file.loc[~condition]

    valid_cas_list['cas_smiles'] = batch.batch_cas_to_smiles(valid_cas_list['CAS'])

    no_smiles_list = valid_cas_list.loc[valid_cas_list['cas_smiles'] == None]['CAS'].to_list()
    cas_smiles_list = valid_cas_list.loc[valid_cas_list['cas_smiles'] != None]

    cas_smiles_whitelist = []
    errors = []
    for cas, smiles in zip(valid_cas_list['CAS'], valid_cas_list['cas_smiles']):
        try:
            if not blacklist_engine.is_similar(smiles):
                cas_smiles_whitelist.append((cas, smiles))
        except Exception as err:
            tb = traceback.format_exc()
            logger.error(f"is_similar failed for ({cas}, {smiles}). Error: {tb}")
            #errors.append(f"{cas}, {smiles}")
            errors.append(f"{cas}")

    reagents = []

    now = time.strftime("%d.%m.%Y %H:%M", time.localtime())
    
    unique_molecules_collection_instance = UniqueMolecules(db_client, MOLECULES_DATABASE)

    for cas, smiles in cas_smiles_whitelist:

        """
        исправляем ошибки в зарядах: делаем все реагенты электронейтральными. 
        этот модуль не оттестирован на прочие ошибки, возможно он приведет к ошибкам при импорте больших баз.
        """
        SMILES = charge_fixer(smiles)
        reagent_internal_id = uuid.uuid4().hex
        """
        регистрируем реагент в базе молекул, в коллекции уникальных молекул. 
        если молекула уже зарегистрирована, то reagent_internal_id реагента из user_reagents вписывается в лист value_data в уникальной записи уникальной молекулы. 
        таким образом, пропуская дубликаты, создаются ссылки из сущности "банка реактивов" на уникальную молекулу, которая содержится в банке и обратно.

        функция регистрации возвращает уникальный индекс - inchikey_standard, который добавляется в запись о реагенте у пользователя. 

        """
        inchikey_standard = unique_molecules_collection_instance.reagent_registration(SMILES)
        r = {
            "reagent_internal_id": reagent_internal_id,
            "inchikey_standard" : inchikey_standard,
            "CAS": cas,
            "SMILES": SMILES,
            "sharing_status": "shared",
            "timestamp": now,
        }
        for option_values in ['location', 'name']:
            if option_values in valid_cas_list.columns:
                values = [i for i in set(valid_cas_list.loc[valid_cas_list['CAS'] == cas][option_values].dropna().values) if i]
                if values:
                    r[option_values] = '\n'.join(values)

        if contact:
            r["contact"] = contact
        reagents.append(r)

    message = f"file was successfully parsed and uploaded.\n"
    message += f"<b>import results</b>:\n"
    message += f"Строк в вашем списке <b>{len(cas_file.index)}</b>\n"
    message += f"Правильных CAS-номеров <b>{len(valid_cas_list)}</b>\n"
    message += f"Опечатка в CAS: <b>{len(failed_cas.index)}</b>\n"
    message += f"Не найдено SMILES для: <b>{len(no_smiles_list)}</b> позиций\n"
    if no_smiles_list:
        message += "\n".join(no_smiles_list) + "\n"
    message += f"Ошибка обработки SMILES <b>{len(errors)}</b> позиций\n"
    if errors:
        message += "\n".join(errors) + "\n"
    message += f"Найдено SMILES для: <b>{len(cas_smiles_list)}</b> реагентов\n"
    message += f"Прекурсоров найдено и вычеркнуто: <b>{len(cas_smiles_list) - len(cas_smiles_whitelist)}</b>\n"

    return reagents, message


class UserReagents:
    """
    This is object for manipulating easily with user 
    attributes before sending updated data to database. 

    test_record = {
        _id: "980159954",
        user_id: "980159954",
        username: "@None",
        firstname: "Alex",
        lastname: "Fedorov"
        laboratory: [
            {
                laboratory_object
            },
            {
                laboratory_object
            }
        ]
        reagent_requests: [
            {
                requested_CAS: "50-00-0"
            }
        ],
        user_reagents: [
            {
                CAS: "50-00-0",
                SMILES: "???",
                #reagent_name: "something 4-something"
                sharing_status: "shared",
                contact: "" # если админ добавил
            }
        ]
    }

    """
    def __init__(self, *args, **kwargs):
        if args:
            self.args = args
        if kwargs:
            for key, value in kwargs.items():
                setattr(self, key, value)

        # присвоить id если его не было 
        if "_id" not in kwargs.keys():
            self._id = uuid.uuid4().hex
    
    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value

    def get_contacts_for_reagent(self, value):
        """
        find value in reagents:
        {
            CAS: "75-64-9",
            SMILES: "CC(C)(C)N",
            reagent_name: "something 4-something"
            sharing_status: "shared"
        }
        :param value:
        :return:
        """
        contacts = []
        for reagent in self.user_reagents:
            if value in reagent.values():
                if reagent["contact"] not in contacts:
                    contacts.append(reagent["contact"])
        return contacts

    def export(self):
        """
        UserReagents_export = {
            "_id": self._id
            "user_id": self.user_id,
            "username": self.username,  
            "time": self.time,
            "firstname": self.firstname,
            "lastname": self.lastname,
            "laboratory": self.laboratory,
            "reagent_requests": self.reagent_requests,
            "user_reagents": self.user_reagents
        }
        """
        return {**{"_id": self._id}, **dict(self)}
        # json.dumps(UserReagents_export) # exports json string (to use it as python object you should convert it by json.loads())
