#!/usr/bin/python
from pymol.cgo import *    
from pymol import cmd
from glob import glob


#chain z
#cmd.load('mobile_pdb/b_1d.pdb','b_1d3z')
#chain y
curr_dir = os.getcwd()
#for files in glob("D_aligned_20000/mnodes*.pdb"):
i = 1
for files in glob("D_aligned_20000/mnodes*.pdb"):
	print(files)
	cmd.load('test_D/new_B.pdb','new_dimer')
	cmd.load(files,"mnodes%s"%str(i),state=1)
	#cmd.align('a_1d3z','/mnodes//A')
	cmd.align('/mnodes%s//B'%str(i),'/new_dimer//a', mobile_state=1, target_state=1)
	cmd.save('B_aln/B_aln_%s.pdb'%str(i))
	cmd.delete('new_dimer')
	cmd.delete('/mnodes%s//A'%str(i))
	i += 1
