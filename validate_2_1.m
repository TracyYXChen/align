%find paramagnetic center in 2+1 pdb file.
%variables
pdb_file = 'test_data/chainA_001.pdb';
pcs_raw = 'data_PCS/Mono_pcs_exp.txt';
out_pcs_file = 'out/pred_pcs.txt';
out_chi_file = 'out/chi_calc.txt';
%chain z aligned with chain A 
%chain y aligned with chain B
chain_zA = {'z'};
chain_zB = {'y'};
which_chi = 'xx';
%get coordinates
pdb_coor_A = dimer_preprocess(pdb_file, chain_zA);
pdb_coor_B = dimer_preprocess(pdb_file, chain_zB);
%remove NA and zero
pcs_table = readtable(pcs_raw, 'HeaderLines', 0);
pcs_res = pcs_table.Var1;
pcs_exp = pcs_table.Var3;
pdb_x = pdb_coor_A(:,1);
pdb_x = pdb_x([pcs_res]);
pdb_y = pdb_coor_A(:,2);
pdb_y = pdb_y([pcs_res]);
pdb_z = pdb_coor_A(:,3);
pdb_z = pdb_z([pcs_res]);
sele_pdb_coor = [pdb_x, pdb_y, pdb_z];
%simplex
guess = [57,30, -20]*10^-10;
options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',1000000,'MaxIter',100000);
fprintf('Start searching process...\n')
[position, Chi2]=fminsearch(@(guess) svd_solver(guess,sele_pdb_coor,pcs_exp,out_pcs_file, out_chi_file, which_chi),guess,options);







