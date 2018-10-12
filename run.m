%% variables
type='Ub2a6.8';
pcs_raw = ['data_PCS/' type '.txt'];
out_pcs_file = 'out/pred_pcs.txt';
out_chi_file = 'out/pred_chi.txt';
%para_center on chain A
name = 'reprod/new_3ns8.pdb';
%3ns8 needs dimer preprocess H02
pdbname = name(end-7:end-4);
%pdbname = '1AARhydr';
% new_3ns8 chain A 
%para_center = [-8.5705   -4.7741   11.3346];
% new_3ns8 chain B 
para_center = [-19.8951    6.6872   34.3609];
% 2bgf chain A
%para_center = [10.8231   16.6988    9.6177];
% 2bgf chain B
%para_center = [1.1403  -15.1452  -15.4806];
% 1AAR chain A
%para_center = [33.5705   25.9841    7.8993];
% 1AAR chain B
%para_center = [8.8028   -8.2468   -8.6791];
%chi_mat for new_3ns8 chain A
%chi_mat=[2.4412,-0.099032,7.6665,-2.5611,1.78]'.*10^-32;
%chi_mat for new_3ns8 chain B
chi_mat=[0.23225,-1.1807,10.431,-2.707,-1.8449]'.*10^-32;
%chi_mat for 2bgf chain A
%chi_mat=[-6.3249,-0.75892,-2.3032,-0.70688,3.2808]'.* 10^-32;
%chi_mat for 2bgf chain B
%chi_mat=[3.1354,3.1607,2.5597,-6.1534,-0.76776]'.*10^-32;
%chi_mat for 1AAR hydr chain A
%chi_mat = [-9.179,1.6056,-1.5625,4.3597,-0.14309]'.*10^-32;
%chi_mat for 1AAR hydr chain B
chi_mat = [-2.7439,1.3022,2.7592,-8.4111,-5.7667]' .*10^-32;
%chain z, y aligned with chain A,B
now_chain = {'A'};
which_chi = 'xx';
%remove NA and zero
t = readtable(pcs_raw, 'HeaderLines', 2);
res_pcs = [t.Var1,t.Var4];
pcs_exp = res_pcs(:,2)*10^-6;
res_pcs(any(isnan(res_pcs),2)|any(res_pcs==0,2),:) = [];
%% read pdb file
num_pdb = 1;
pcs_A = [];
%% 
for ii=1:num_pdb
    %ii
    %prefix = 'single_aln_A/A_aln';
    %pdb_coor_B = dimer_preprocess([prefix '_' num2str(ii) '.pdb'], now_chain);
    pdb_coor_B = dimer_preprocess(name, now_chain);
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
    %[U, S, V] = svd(A);
    %A_inv = V * pinv(S) * U';
    %chi_mat = A_inv * pcs_exp;
    %pcs_calc = A * chi_mat;
    %pcs_A = [pcs_A, pcs_calc];
end
pcs_calc = A * chi_mat * 10^6;
dlmwrite(['/Users/yuexichen/Desktop/nmr/experiment/SES/data/pcs/Amat/' pdbname 'Amat_' type '.txt'],pcs_A,'\t')
err = ones(length(pcs_exp),1.0000);
pcs_err = [pcs_exp, err];
dlmwrite(['/Users/yuexichen/Desktop/nmr/experiment/SES/data/pcs/y/' type '_y_err.txt'], pcs_err,'\t')
scatter(pcs_exp, pcs_calc);
xlabel("PCS experiment");
ylabel("PCS returned by SES");
title([pdbname '  PCS pred vs exp single conformation ' type]);
hold on;
x = -0.2 : 0.01 : 0;
y = -0.2 : 0.01 : 0;
plot(x,y);









