import os, sys
import pandas as pd
from validator.validator import Validator
from validator.residue_index import ResidueIndexes


def validate(pwd,filename):
    df=pd.read_csv(f'{pwd}/{filename}',sep='\t',header='infer')

    pdb_list=list();pdb_skip=list()
    
    for i in df.index:
        pdb=df.at[i,'PDBid'][0:4].lower()     #filename is in lower case
    
        if str(df.at[i,'Author_Aspnum']).lower() == 'nan':   #Do not include the structures where Phe is not resolved; Pandas read null as nan
                pdb_skip.append(pdb)
                continue
    
        if pdb in pdb_list:
            continue
        pdb_list.append(pdb)
        dir_name=pdb[1:3].lower()
    
        validator = Validator("Kincore") # Same as in the JSON
        validator.load_schema(f"{pwd}/funpdbe-validator/data/funpdbe_schema.json")
        validator.load_json(f"{pwd}/JSON/{dir_name}/{pdb}.json")
        print(pdb)
        if validator.basic_checks() and validator.validate_against_schema():
              # Passed data validations
              print('Passed schema validation')
              residue_indexes = ResidueIndexes(validator.json_data)
              if residue_indexes.check_every_residue():
                  # Passed the index validation
                  print('Passed index validation')
              else:
                  print(residue_indexes.mismatches)
        else:
              print('Failed schema validation')
              print(validator.error_log)


if __name__=='__main__':
    pwd='/home/vivek/Applications/Flask/Kincore'     #location in workhorse
    filename=sys.argv[1]   #Input csv file
    validate(pwd,filename)