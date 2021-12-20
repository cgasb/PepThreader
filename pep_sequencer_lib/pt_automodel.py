import modeller
from modeller import *
from modeller.automodel import *


class PT_automodel(automodel):
    """
    'automodel' class for pep_threader.
    TODO:
        - test: 7.0, 8.0, 9.0 radiuses.
        - test a shorter MDSA protocol.
        - test the extraction of MDSA trajectories.
        - add DOPE or SOAP in the objective functionself.
        - tweak the restraints.
    """

    def set_defaults(self):
        automodel.set_defaults(self)

        # Selectes a sphere of 'self.pt_select_sphere_radius' Angstroms around the peptide atoms
        # and optimize only the atoms inside it.
        self.pt_select_sphere_radius = 8.0 # 9.0
        # If set to 2, the atoms in the selection will not feel other atoms in the optimization process.
        self.pt_nonbonded_sel_atoms = 2


    def select_atoms(self):
        # Select residues neighbouring the peptide chain.
        sel = selection(self.chains['B']).select_sphere(self.pt_select_sphere_radius).by_residue()
        # .extend_by_residue(1)
        return sel
    

    def special_restraints(self, aln):
        '''
        self.env.edat.nonbonded_sel_atoms = self.pt_nonbonded_sel_atoms
        return None
        '''

        return None

        self.env.schedule_scale = physical.values(default=1,
                                                  nonbond_spline=0.5,
                                                  ca_distance=1)

        # Allow calculation of statistical (dynamic_modeller) potential
        edat = self.env.edat
        edat.contact_shell = 7.0
        # edat.dynamic_sphere=False
        # edat.dynamic_lennard=True
        # edat.dynamic_coulomb=False
        # edat.relative_dielectric=1.0

        edat.dynamic_modeller = True # Allow calculation of statistical (dynamic_modeller) potential

        #----------------
        # Energy terms. -
        #----------------

        # GBSA.
        # edat.energy_terms.append(gbsa.Scorer(cutoff=edat.contact_shell))

        # SOAP.
        # from modeller import soap_pp
        # edat.energy_terms.append(soap_pp.PairScorer())

        #--------------------
        # Group restraints. -
        #--------------------

        # Read Fiser/Melo loop modeling potential
        # gprsr = group_restraints(self.env, classes='$(LIB)/atmcls-melo.lib', parameters='$(LIB)/melo1-dist.lib')

        # Read DOPE loop modeling potential (the same one used in assess_dope).
        gprsr = group_restraints(self.env, classes='$(LIB)/atmcls-mf.lib', parameters='$(LIB)/dist-mf.lib')

        # Read DOPE-HR loop modeling potential
        # gprsr = group_restraints(self.env, classes='$(LIB)/atmcls-mf.lib', parameters='$(LIB)/dist-mfhr.lib')

        # Use this potential for the 1fas model
        self.group_restraints = gprsr
