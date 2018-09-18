#!/usr/bin/python
from pymol.cgo import *    
from pymol import cmd
from glob import glob

cmd.load('a_1d.pdb','a_1d3z')
#chain z
cmd.load('b_1d.pdb','b_1d3z')
#chain y
curr_dir = os.getcwd()
i = 1
for files in glob("D_aligned_20000/mnodes*.pdb"):
	cmd.load(files,"mnodes")
	cmd.align('a_1d3z','/mnodes//c')
	cmd.align('b_1d3z','/mnodes//d')
	cmd.save('myalign/ali_' + str(i) +'.pdb')
	cmd.delete(all)
	i += 1
