from pathlib import Path
from incar import Incar

def test_incar():

    incar_test_file = Incar(Path('vasprun.xml'))
    
    #make sure that the incar dictionary is not empty
    assert incar_test_file.incar.keys() != []

    #make sure that the incar dictionary keys only has int, float, or str values
    for key, value in incar_test_file.incar.items():
        assert type(key) == str
        assert type(value) == str or type(value) == int or type(value) == float

    #add a new key to the incar dictionary
    incar_test_file['test_key'] = 'test_value'

    #make sure that the new key is in the incar dictionary
    assert 'test_key' in incar_test_file.incar.keys()

    


if __name__ == '__main__':
    test_incar()
