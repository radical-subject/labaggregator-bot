
import logging
import traceback
logger = logging.getLogger(__name__)
from modules.db.dbschema import parse_cas_list
from unittest.mock import patch


def parse_cas_list():
    #TODO пока не победил
    #with patch("modules.ourbot.service.cas_to_smiles.banch_cas_to_smiles",
    #           side_effect=[[("15243-33-1", "C"), ("917-64-6", "CC"), ("94-02-0", "CCC")]]) as cas_to_smiles_patched:

    with patch("cirpy.resolve", side_effect=["C", "CC",  "CCC"]) as cirpy_patched:

            with patch("modules.db.blacklist.blacklist_engine.is_similar",
                       side_effect=[False, False, False]) as is_similar_patched:

                contact = None
                cas_list = ["15243-33-1", "917-64-6", "94-02-0"]

                reagents, text_report = parse_cas_list(cas_list, contact)

                assert len(reagents) == 3

                logger.info(text_report)
