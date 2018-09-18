#!/usr/bin/python

from Bio.PDB import *
p = PDBParser()
ub = p.get_structure('A_1d3z', 'D_aligned_20000/1d3z.pdb')
ub_model = ub[0]
ub_chain = ub_model['A'] 

##
dimer = p.get_structure('dimer','D_aligned_20000/mnodes001.pdb')
dimer_model = dimer[0]
dimer_chain = dimer_model['B']

##
ub[0]['A'] = dimer_chain
io=PDBIO()
io.set_structure(ub)
io.save("out.pdb")

