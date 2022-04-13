import sys
import argparse
# from logger import log
import csv
import cirpy, pubchempy
import re


def cirpy_smiles_resolve(cas: str):
    """
    param: cas - line 1-1-1
    """
    return cirpy.resolve(cas, 'smiles')


def pubchempy_smiles_resolve(cas: str):
    """
    param: cas - line 1-1-1
    """
    pubchem_response = pubchempy.get_compounds(cas, "name")
    return pubchem_response[0].isomeric_smiles


def cas_to_smiles(cas: str):
    """
    param: cas - line 1-1-1
    """
    try:
        # log.info(f'Processing... {cas}')
        res = cirpy_smiles_resolve(cas)
        if not res:
            return pubchempy_smiles_resolve(cas)

        return res
    except Exception as err:
        # log.error(f'{cas}: cirpy_smiles_resolve error: {err}')
        return pubchempy_smiles_resolve(cas)


def parse_lines(lines: str):

    result = []
    for cas in lines:
        cas = cas.strip()
        try:
            smiles = cas_to_smiles(cas)
            result.append((cas, smiles))
            # log.info(f'{cas} - {smiles}')
        except Exception as err:
            # log.error(f'{cas}: {err}')
            result.append((cas, 'resolver_error'))

    result_object_list = [{"CAS": i[0], "SMILES": i[1]} for i in result] # if i[0]!="resolver_error"
    
    errors_CAS_list = [i[0] for i in result if i[1]=="resolver_error"]

    # # remove all error indications - this gives clean SMILES list
    # SMILES_list = list(filter(("resolver_error").__ne__, SMILES_list))
    
    return (result_object_list, errors_CAS_list)


if __name__ == "__main__":

    argparser = argparse.ArgumentParser()
    argparser.add_argument('-i', dest='inpath', help='txt file one line one cas')
    argparser.add_argument('-o', dest='outpath', help='outputfile')

    args = argparser.parse_args()

    if not args.inpath:
        # log.error('no infile argument')
        sys.exit(0)
    else:
        # log.info(f'Parse input file: {args.inpath}')
        pass

    with open(args.inpath, 'r') as infile:
        lines = infile.readlines()
        out = parse_lines(lines)
        print(out)
        if args.outpath:
            # log.info(f'Write output file: {args.outpath}')
            with open(args.outpath, 'w') as outfile:
                for cas, smiles in out:
                    outfile.write(f'{cas}\t{smiles}')
