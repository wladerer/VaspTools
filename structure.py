

class Structure:

    def __init__(self, vasprun_xml_path: Path | str):
        '''Creates a Structure object from a vasprun.xml file.'''
        self.path = path
        self.root = read_vasprun(path)
        self.vasprun = Vasprun(path)

    @property
    def initial_stucture_element(self):
        initial_pos_element = self.root.find('structure[@name="initialpos"]')
        initial_structure_element = initial_pos_element.find('crystal')
        return initial_structure_element

    @property
    def initial_basis(self) -> np.ndarray:
        initial_basis_varray = self.initial_stucture_element.find(
            'varray[@name="basis"]')
        initial_basis = unpack_varray(initial_basis_varray)

        return initial_basis

    @property
    def initial_reciprocal_basis(self) -> np.ndarray:
        initial_reciprocal_basis_varray = self.initial_stucture_element.find(
            'varray[@name="rec_basis"]')
        initial_reciprocal_basis = unpack_varray(
            initial_reciprocal_basis_varray)

        return initial_reciprocal_basis

    @property
    def initial_positions(self) -> np.ndarray:
        initial_positions_varray = self.root.find(
            'structure[@name="initialpos"]').find('varray[@name="positions"]')
        initial_positions = unpack_varray(initial_positions_varray)

        return initial_positions

    @property
    def final_stucture_element(self):
        final_pos_element = self.root.find('structure[@name="finalpos"]')
        final_structure_element = final_pos_element.find('crystal')

        return final_structure_element

    @property
    def final_basis(self) -> np.ndarray:
        final_basis_varray = self.final_stucture_element.find(
            'varray[@name="basis"]')
        final_basis = unpack_varray(final_basis_varray)

        return final_basis

    @property
    def final_reciprocal_basis(self) -> np.ndarray:
        final_reciprocal_basis_varray = self.final_stucture_element.find(
            'varray[@name="rec_basis"]')
        final_reciprocal_basis = unpack_varray(final_reciprocal_basis_varray)

        return final_reciprocal_basis

    @property
    def final_positions(self) -> np.ndarray:
        final_positions_varray = self.root.find(
            'structure[@name="finalpos"]').find('varray[@name="positions"]')
        final_positions = unpack_varray(final_positions_varray)

        return final_positions

    @property
    def formula(self) -> str:
        formula_dict = {atom_type: self.vasprun.atom_types.count(
            atom_type) for atom_type in self.vasprun.atom_types}
        formula = ''.join(
            [f'{atom_type}{formula_dict[atom_type]}' for atom_type in formula_dict])
        # remove 1s
        formula = formula.replace('1', '')

        return formula

    @property
    def selective_dynamics(self) -> list[list[bool]]:
        '''Returns a list of lists of booleans indicating whether the atom is fixed in that direction'''
        selective_dynamics_element = self.root.find(
            'structure[@name="initialpos"]').find('varray[@name="selective"]')
        selective_dynamics = [v.text.split()
                              for v in selective_dynamics_element.findall('v')]

        return selective_dynamics


def write_poscar(self, filename: str = 'POSCAR', final: bool = False, scale: float = 1.0):
    '''Writes a POSCAR file with the final or initial structure'''
    basis = self.final_basis if final else self.initial_basis
    positions = self.final_positions if final else self.initial_positions

    with open(filename, 'w') as f:
        f.write(f'{self.formula}\n{scale}\n')
        for basis_vector in basis:
            f.write(
                f'{basis_vector[0]:8f} {basis_vector[1]:8f} {basis_vector[2]:8f}\n')

        atom_dict = {atom_type: self.vasprun.atom_types.count(
            atom_type) for atom_type in self.vasprun.atom_types}
        atom_line = ' '.join(atom_dict.keys())
        atom_counts = ' '.join(str(count) for count in atom_dict.values())

        f.write(f'{atom_line}\n{atom_counts}\nDirect\n')
        selective_dynamics = self.selective_dynamics
        for atom, position in enumerate(positions):
            selective = ' '.join(selective_dynamics[atom])
            f.write(
                f'{position[0]:.8f} {position[1]:.8f} {position[2]:.8f} {selective} {self.vasprun.atom_types[atom]}\n')