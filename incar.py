from pathlib import Path
from typing import Union
import numpy as np
from vasprunio import read_vasprun
import yaml

class Incar:

    def __init__(self, path: Path):
        '''Creates an Incar object from a vasprun.xml file.'''
        self.path = path
        self.root = read_vasprun(path)
        self.incar = self.create_incar_dictionary()


    def create_incar_dictionary(self) -> dict:
        '''Creates a dictionary of the INCAR file from vasprun.xml'''
        incar_children = self.root.find('incar').getchildren()
        incar_dict = {child.attrib['name']                      : child.text for child in incar_children}
        incar_dict = {k: v.strip() for k, v in incar_dict.items()}
        #remove KPOINT_BSE
        incar_dict.pop('KPOINT_BSE', None)
        #update T with True and F with False
        incar_dict = {k: True if v == 'T' else False if v == 'F' else v for k, v in incar_dict.items()}
        return incar_dict
    
    def __getitem__(self, key: str) -> str:
        return self.incar[key]
    
    def __setitem__(self, key: str, value: str):
        
        try: 
            
            self.incar[key] = value

        except KeyError:

            print(f'ERROR: INCAR key {key} does not exist')
            print('ERROR: Check the INCAR file for typos')


    def __str__(self) -> str:
        incar_str = ''
        for key, value in self.incar.items():
            incar_str += f'{key} = {value} \n'

        return incar_str

    def write_incar(self, path: Path, dry_run: bool = False):
        '''Writes an identical INCAR file for the current calculation'''

        if dry_run:
            print(self.__str__() )
        else:
            with open(path, 'w') as f:
                f.write(self.__str__())


    def increment(self, increment_dict: dict[str, list], dry_run: bool = False):
        '''Increments the INCAR file by a list of values for a given INCAR key'''

        for key, values in increment_dict.items():
            for value in values:
                self.incar[key] = value
                self.write_incar(f'INCAR_{key}_{value}', dry_run = dry_run)

    
    def to_yaml(self, path: Path):
        '''Writes an identical INCAR file for the current calculation'''

        with open(path, 'w') as f:
            yaml.dump(self.incar, f)

 




