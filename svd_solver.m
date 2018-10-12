%% -----Info-----
%Solve Pseudocontact Shift(PCS) tensors 
%Yuexi (Tracy) Chen
%July 31, 2018
%% -----solve svd -----
function [Chi_2] = svd_solver(guess, pdb_coor, pcs_exp, out_pcs_file, out_chi_file, which_chi)
%guess: initialization of paramagnetic center
%pbd_coor: a list of pdb coordinates
%pcs_exp: pcs experiment
%pcs_file: file to store exp and predicted pcs for comparison
%construct A matrix
para_center = guess;
num_res = length(pcs_exp);
A = zeros(num_res,5);
pdb_x = pdb_coor(:,1);
x = (pdb_x - para_center(1))*10^-10;
pdb_y = pdb_coor(:,2);
y = (pdb_y - para_center(2))*10^-10;
pdb_z = pdb_coor(:,3);
z = (pdb_z - para_center(3))*10^-10;
r_sqr = x.^2 + y.^2 +z.^2;
if which_chi == 'zz'
    A(:,1)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* [x.^2 - z.^2];
    A(:,2)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* [2.*x.*y];
    A(:,3)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* [2.*x.*z];
    A(:,4)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* [y.^2 - z.^2];
    A(:,5)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* [2.*y.*z];
elseif which_chi == 'xx'
    A(:,1)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* (2*x.*y); 
    A(:,2)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* (2*x.*z); 
    A(:,3)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* (y.^2 - x.^2); 
    A(:,4)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* (2*y.*z); 
    A(:,5)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* (z.^2 - x.^2); 
elseif which_chi == 'yy'
    fprintf('underconstruction');
else
    fprintf('which_chi could only be xx,yy or zz, others are not supported');
%SVD
end
fprintf('condition number of A is %f\n',cond(A));
[U, S, V] = svd(A);
%x = A-1 * PCS_exp
%S is not square, use pinv
A_inv = V * pinv(S) * U';
chi_mat = A_inv * pcs_exp;
if which_chi == 'zz'
    [chi_xx,chi_xy,chi_xz,chi_yy,chi_yz] = deal(chi_mat(1),chi_mat(2),chi_mat(3),chi_mat(4),chi_mat(5));
    chi_zz = -(chi_xx + chi_yy);
elseif which_chi == 'xx'
    [chi_xy,chi_xz,chi_yy,chi_yz,chi_zz] = deal(chi_mat(1),chi_mat(2),chi_mat(3),chi_mat(4),chi_mat(5));
    chi_xx = -(chi_yy + chi_zz);
elseif which_chi == 'yy'
    [chi_xx,chi_xy,chi_xz,chi_yz,chi_zz] = deal(chi_mat(1),chi_mat(2),chi_mat(3),chi_mat(4),chi_mat(5));
    chi_yy = -(chi_xx + chi_zz);
else
    fprintf('which_chi could only be xx,yy or zz, others are not supported');
end
pcs_calc = A * chi_mat;
dlmwrite(out_chi_file, [chi_xx, chi_xy,chi_xz,chi_yy,chi_yz, chi_zz]);
norm_pcs_calc = 10^6 * pcs_calc;
dlmwrite(out_pcs_file, [norm_pcs_calc]);
Chi_2 = sum((pcs_calc - pcs_exp).^2);
return


