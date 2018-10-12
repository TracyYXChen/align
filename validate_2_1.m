%find paramagnetic center in 2+1 pdb file.
%variables
%chain A and chain B are the same
pdb_file = 'reprod/1AAR_hydr.pdb';
pcs_raw = 'data_PCS/Mono_pcs_exp.txt';
out_pcs_file = 'out/pred_pcs.txt';
out_chi_file = 'out/chi_calc.txt';
%chain z aligned with chain A 
%chain y aligned with chain B
chain = {'A'};
%chain = {'B'};
which_chi = 'xx';
sele_atom = {'H'}
% new_3ns8 chain A 
%guess = [-8.5705   -4.7741   11.3346];
% new_3ns8 chain B 
%guess = [-19.8951    6.6872   34.3609];
% 2bgf chain A
%guess = [10.8231   16.6988    9.6177];
% 2bgf chain B
%guess = [1.1403  -15.1452  -15.4806];
% 1AARhydro A
guess =[33.5705   25.9841    7.8993];
% 1AARhydro B
%guess =[8.8028   -8.2468   -8.6791];
%get coordinates
pdb_coor_A = dimer_preprocess(pdb_file, chain,sele_atom);
%remove NA and zero
pcs_table = readtable(pcs_raw, 'HeaderLines', 0);
pcs_res = pcs_table.Var1;
pcs_exp = pcs_table.Var3*10^-6;
pdb_x = pdb_coor_A(:,1);
pdb_x = pdb_x([pcs_res]);
pdb_y = pdb_coor_A(:,2);
pdb_y = pdb_y([pcs_res]);
pdb_z = pdb_coor_A(:,3);
pdb_z = pdb_z([pcs_res]);
sele_pdb_coor = [pdb_x, pdb_y, pdb_z];
tic
best_Chi2 = 1;
options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',1000000,'MaxIter',100000);
fprintf('Start searching process...\n')
[position, Chi2]=fminsearch(@(guess) svd_solver(guess,sele_pdb_coor,pcs_exp,out_pcs_file, out_chi_file, which_chi),guess,options);
position
Chi2
toc






