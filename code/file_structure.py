from pathlib import Path

class Simulation:
    """
    * The analyses use `Simulation` objects to specify a simulation to be loaded.
    * The Simulation class knows about folder structure, filenames etc., not about
      how those files should be processed.
    * Any of the H-REMD simulations with quadratic or exponential scheme can be
      chosen by passing a Simulation object to the analysis functions. To change
      settings or use a different kind of simulation setup, sub-class Simulaton or
      make an equivalent class.
    * Optimally, all the low-level functions should work identically with a different
      Simulation class. However, some of the higher-level functions (plotting etc.)
      might have to be adapted.
      
    **Usage:** `Simulation("compound", "forcefield", "simulation_method", "solvent", state=0)`
    """
    def __init__(self, compound, forcefield, method, solvent, state=0):
        self.compound = compound
        self.forcefield = forcefield
        self.method = method
        self.solvent = solvent
        self.state = state

    @classmethod
    def fields(cls):
        return ["compound", "forcefield", "method", "solvent", "state"]

    def __repr__(self):
        args = ", ".join(repr(getattr(self, field)) for field in self.fields())
        return f"{self.__class__.__name__}({args})"

    def trajectory_path(self):
        return str(self.basepath() / f'remd{self.state}' / 'traj_comp.xtc')

    def basepath(self):
        return Path(f'../data/{self.compound}/{self.solvent}/{self.method}/{self.forcefield}/outputs')

    def logfile_path(self):
        return self.basepath() / f'remd{self.state}' / 'logfile.log.gz'

    def n_states(self):
        return 12
        
    def with_param(self, key, value):
        out = type(self)(self.compound, self.forcefield, self.method, self.solvent)
        setattr(out, key, value)
        return out

    def with_compound(self, value):
        return self.with_param('compound', value)

    def with_forcefield(self, value):
        return self.with_param('forcefield', value)

    def with_method(self, value):
        return self.with_param('method', value)

    def with_solvent(self, value):
        return self.with_param('solvent', value)

    def pdb_filename(self):
        return str(self.basepath() / 'npt-bonds.pdb')

    def equilibration_length(self):
        return 10_000

    def as_tuple(self):
        return tuple([getattr(self, name) for name in self.fields()])

    def __hash__(self):
        return hash(self.as_tuple())

    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()


def molfile_path(compound: str):
    return f"../data/{compound}/reference-atom-order.mol"


def atom_assignment_path(compound, solvent):
    return f"../data/{compound}/{solvent}/atom-assignment.csv"


def restraintfile_path(compound, solvent):
    return f"../data/{compound}/{solvent}/noe-bounds.csv"
