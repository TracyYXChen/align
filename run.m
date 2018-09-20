%% variables
pdb_file = 'yxc_aln/ali_1.pdb';
pcs_raw = 'data_PCS/Ub2b6.8.txt';
out_pcs_file = 'out/pred_pcs.txt';
out_chi_file = 'out/pred_chi.txt';
%para center on chain A
para_center = [32.92,26.63,9.19] * 10^-10;
%chain z aligned with chain A 
%chain y aligned with chain B
now_chain = {'z'};
which_chi = 'xx';
%remove NA and zero
t = readtable(pcs_raw, 'HeaderLines', 2);
res_pcs = [t.Var1,t.Var4];
pcs_exp = res_pcs(:,2);
res_pcs(any(isnan(res_pcs),2)|any(res_pcs==0,2),:) = [];
%% read pdb file
num_pdb = 200;
A = zeros(length(pcs_exp),num_pdb)
for ii=1:200
    prefix = 'yxc_aln/ali';
    pdb_coor_B = dimer_preprocess(prefix '_' ii, now_chain);
    sele_B_x = pdb_coor_B(:,1);
    sele_B_y = pdb_coor_B(:,2);
    sele_B_z = pdb_coor_B(:,3);
    res_num = res_pcs(:,1);
    pcs_exp = res_pcs(:,2);
    %remove residue 75,76 ...
    nozero = res_num(1:length(res_num)-3);
    nonzero_B_x = sele_B_x([nozero]);
    nonzero_B_y = sele_B_y([nozero]);
    nonzero_B_z = sele_B_z([nozero]);
    pcs_exp = pcs_exp(1:length(res_num)-3);
    pdb_coor_B = [nonzero_B_x, nonzero_B_y, nonzero_B_z];
    %build A matrix on chain B
    A = build_Amat(pcs_exp, para_center, pdb_coor_B);
    [U, S, V] = svd(A);
    A_inv = V * pinv(S) * U';
    chi_mat = A_inv * pcs_exp;
    pcs_calc = A * chi_mat;
    pcs_A(:,ii) = pcs_calc;










