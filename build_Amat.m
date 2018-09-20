%-----
%Build A matrix for PCS data
%Yuexi (Tracy) Chen, Sept 2018
%-----
function [A]=build_Amat(pcs_exp, para_center,pdb_coor)
num_res = length(pcs_exp);
pcs_exp = pcs_exp;
A = zeros(num_res,5);
pdb_x = pdb_coor(:,1);
x = pdb_x - para_center(1);
pdb_y = pdb_coor(:,2);
y = pdb_y - para_center(2);
pdb_z = pdb_coor(:,3);
z = pdb_z - para_center(3);
r_sqr = x.^2 + y.^2 +z.^2;
A(:,1)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* (2*x.*y); 
A(:,2)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* (2*x.*z); 
A(:,3)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* (y.^2 - x.^2); 
A(:,4)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* (2*y.*z); 
A(:,5)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* (z.^2 - x.^2); 
