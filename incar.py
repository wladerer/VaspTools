from pathlib import Path
from vasprunio import read_vasprun

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
        # trim leading and trailing whitespace from the values
        incar_dict = {k: v.strip() for k, v in incar_dict.items()}

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

    


 




