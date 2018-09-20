%-----
%extract selected residues from a chain of a pdb file
%Yuexi (Tracy) Chen Sept-2018
%-----
function [pdb_coor] = dimer_preprocess(pdb_file, mychain)
%---read pdb file ---
PDBdata = pdb2mat(pdb_file);
fprintf('Finished reading PDB file.\n')
%num of all atoms and selected atoms
num_all = length(PDBdata.resNum);
sele_num = 76;
%create an array to store coordinates of selected atoms
pdb_coor = zeros(sele_num,3);
for ii= 1:sele_num
    for jj=1: num_all
      if (PDBdata.resNum(jj)==ii) && isequal(PDBdata.atomName(jj),{'H'})...
              && isequal(PDBdata.chainID(jj),mychain)
          pdb_coor(ii,1) = PDBdata.X(jj)*10^-10; 
          pdb_coor(ii,2) = PDBdata.Y(jj)*10^-10;
          pdb_coor(ii,3) = PDBdata.Z(jj)*10^-10; 
      end
    end
end
