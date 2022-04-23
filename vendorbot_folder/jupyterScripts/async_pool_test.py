from multiprocessing import Pool
import sys

sys.path.append("..")
import pandas as pd

from sandbox.resolver import resolve
from modules.service import Timer

timer = Timer()

import_CAS_df = pd.read_csv('../srs/user_reagent_lists_import/radical_subject_reagents.txt', header = None)
radical_subject_CAS_list = import_CAS_df[0].tolist()

data = pd.read_csv('../srs/user_reagent_lists_import/Chusov_1.txt', header = None)
chusov_reagents_list = data[0].tolist()

# radical_subject_CAS_list = [row[1][0] for row in import_CAS_df.iterrows()]

# print(radical_subject_CAS_list)

# print(chusov_reagents_list)

if __name__ == '__main__':
    
    timer.start()
    n = 50
    with Pool(processes=n) as pool:
        smiles_list = pool.map(resolve, chusov_reagents_list)
        
    print(smiles_list)
    timer.stop()
    # freeze_support()      