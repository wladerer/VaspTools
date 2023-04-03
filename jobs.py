from incar import Incar
from structure import Structure
from kspace import ElectronicStructure
from pathlib import Path
import yaml

class Job:

    def __init__(self, incar: Incar, structure: Structure, electronics: ElectronicStructure, file: str = None):
        self.incar: Incar = incar
        self.structure: Incar = structure
        self.electronics: ElectronicStructure = electronics
        self.incar_dict: dict = incar.incar
        self.kpoint_dict: dict = electronics.as_dict()
        self.structure_dict: dict = structure.as_dict()
        self.title: str = file

    @classmethod
    def from_xml(cls, xml_file: str | Path): 
        '''Reads in an XML file and returns a Job object'''
        incar = Incar(xml_file)
        structure = Structure(xml_file)
        electronics = ElectronicStructure(xml_file)
        return cls(incar, structure, electronics)


    @property
    def yaml_filename(self):
        
        if self.title:
            return self.title
        else:
            prefix = 'gen_'
            suffix = '.yaml'
            if 'SYSTEM' not in self.incar_dict:
                system = self.structure.formula
            else:
                system = self.incar_dict['SYSTEM']
                
            filename = prefix + system + suffix
            return filename


    def __repr__(self):
        return f'Job({self.incar}, {self.structure})'
    
    def __str__(self):
        return f'VASP Job for {self.structure.formula}'
    
    def __eq__(self, other):
        return self.incar == other.incar and self.structure == other.structure
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def to_yaml(self): # fix this
        job_dict = {'incar': self.incar_dict, 'structure': self.structure_dict, 'kpoints': self.kpoint_dict}
        with open(self.yaml_filename, 'w') as file:
            yaml.dump(job_dict, file)





