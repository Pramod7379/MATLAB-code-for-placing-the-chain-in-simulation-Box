function [] = Sort_Molecules()
% Syntax : Sort_Molecules()
% This function sorts the data according to molecule type, and stores it in the 'lammpsfile.mat'. It also creates a sample for each type of molecule, storing their atom, bond, angle and sihedral sequences. 
% For more details, refer to the manual provided with the package.
if exist('lammpsfile.mat')
    load 'lammpsfile.mat';
else
    disp('Please upload lammps data before attempting to set parameters');
end
%----------------------------------------------------------------------------
% initialize MoleculeData
MoleculeData = cell(1,1);
%-----------------------------------------------------------------------------------------------
% sorting algorithm 
AtomDataMod = [];
AtomDataMod = AtomData;
Pattern = cell(1,1);
w = 1;
while isempty(AtomDataMod) ~= 1
    % storing the data for the first molecule of a type
    x = 1;
    size_AtomDataMod = size(AtomDataMod);
    z = AtomDataMod(1,2);
    MoleculeData{w} = [];
    while ((x<=size_AtomDataMod(1)) && (AtomDataMod(x, 2) == z))
        MoleculeData{w} = [MoleculeData{w}; AtomDataMod(x, :)];
        x=x+1;
    end
    % storing the atom sequences for a type of molecule
    y = 1;
    zz = AtomDataMod(1,2);
    Pattern{w} = [];
    while ((y<=size_AtomDataMod(1)) && (AtomDataMod(y, 2) == zz))
        Pattern{w} = [Pattern{w}, AtomDataMod(y, 3)];
        y=y+1;
    end
    
    size_patt = size(Pattern{w});
    % finding the all the molecules of that type
    patt_rec =  strfind(AtomDataMod(:, 3)', Pattern{w});
    patt_rec(1) = [];
    size_patt_rec = size(patt_rec);
    ind_mat = [];
    for ii=1:size_patt_rec(2)
        if AtomDataMod(patt_rec(ii), 2) ~= AtomDataMod((patt_rec(ii))-1, 2)
            indices = find(AtomDataMod(:,2) == AtomDataMod(patt_rec(ii), 2));
            ind_mat = [ind_mat; indices];
            size_indices = size(indices);
            if size_indices(1) == size_patt(2)
                MoleculeData{w} = [MoleculeData{w}; AtomDataMod(indices, :)];
            end
        end
    end
    ind_mat = [[1:(x-1)]'; ind_mat]; 
    AtomDataMod(ind_mat, :) = [];
    MoleculeData = [MoleculeData; []];
    w=w+1;
end
size_in_type = [];
disp(['There is/are ' num2str(w-1) ' type(s) of molecules.']);
disp('The atomic sequences (by atom type) in that/each type -');
for jj=1:(w-1)
    disp(['Type ' num2str(jj) ' - ' num2str(Pattern{jj}) '.']);
    dum_size = size(Pattern{jj});
    size_in_type = [size_in_type; dum_size(2)];
end

disp('No. of molecules in that/each type - ');
NumMolInType = [];
for kk=1:(w-1)
    num_in_type = size(MoleculeData{kk});
    disp(['Type ' num2str(kk) ' - ' num2str((num_in_type(1))/size_in_type(kk, 1)) '.']);
    NumMolInType = [NumMolInType; ((num_in_type(1))/(size_in_type(kk, 1)))];
end

NumAllMolecules = unique(AtomData(:, 2));
NumTypeMolecules = w-1;
%--------------------------------------------------------------------------------------------
% sample AtomData, AngleData, BondData, DihedralData
sample_BondData = cell(NumTypeMolecules, 1);
sample_AngleData = cell(NumTypeMolecules, 1);
sample_DihedralData = cell(NumTypeMolecules, 1);
sample_ImproperData = cell(NumTypeMolecules, 1);

sample_Molecules = cell(NumTypeMolecules, 1);
molecule1_of_each_type = [];
% taking the first molecule of each type, to be used as the sample for all
% molecules of that type
for ii = 1:NumTypeMolecules
    dummy = unique(MoleculeData{ii}(:, 2));
    molecule1_of_each_type = [molecule1_of_each_type; dummy(1,1)];
end    
% filling up the sample Data matrices
for jj = 1:NumTypeMolecules
    indices = find(AtomData(:, 2) == molecule1_of_each_type(jj, 1));
    
    % sample_Molecules
    sample_Molecules{jj} = [sample_Molecules{jj}, AtomData(indices, :)];
    
    % sample_BondData
    num_bonds = 1;
    isel1 = [];
    isel2 = [];
    for nn=1:length(AtomData)
       isel1 = [isel1; find(BondData(:,3) == indices(nn, 1))];
       isel2 = [isel2; find(BondData(:,4) == indices(nn, 1))];
    end
    isel = intersect(isel1(:, 1), isel2(:, 1));
    size_isel = size(isel);
    % Entering data into the sample_BondData 
    for kk = 1:size_isel(1)
        atom1 = BondData(isel(kk, 1), 3);
        atom2 = BondData(isel(kk, 1), 4);
        sample_BondData{jj} = [sample_BondData{jj}; [(num_bonds), (BondData(isel(kk, 1), 2)), (atom1), (atom2)]];
        num_bonds = num_bonds + 1;
    end

    % sample_AngleData
    num_angles = 1;
    isel_ang1 = [];
    isel_ang2 = [];
    isel_ang3 = [];
    for nn=1:length(AtomData)
       isel_ang1 = [isel_ang1; find(AngleData(:,3) == indices(nn, 1))];
       isel_ang2 = [isel_ang2; find(AngleData(:,4) == indices(nn, 1))];
       isel_ang3 = [isel_ang3; find(AngleData(:,5) == indices(nn, 1))];
    end
    isel_ang_dum = intersect(isel_ang1(:, 1), isel_ang2(:, 1));
    isel_ang = intersect(isel_ang_dum, isel_ang3(:, 1));
    size_isel_ang = size(isel_ang);
    % Entering data into the sample_AngleData 
    for kk = 1:size_isel_ang(1, 1)
        atom1 = AngleData(isel_ang(kk, 1), 3);
        atom2 = AngleData(isel_ang(kk, 1), 4);
        atom3 = AngleData(isel_ang(kk, 1), 5);
        sample_AngleData{jj} = [sample_AngleData{jj}; [(num_angles), (AngleData(isel_ang(kk, 1), 2)), (atom1), (atom2), (atom3)]];
        num_angles = num_angles + 1;
    end
    
    % sample_DihedralData
    num_dihedrals = 1;
    isel_dih1 = [];
    isel_dih2 = [];
    isel_dih3 = [];
    isel_dih4 = [];
    for nn=1:length(AtomData)
       isel_dih1 = [isel_dih1; find(DihedralData(:,3) == indices(nn, 1))];
       isel_dih2 = [isel_dih2; find(DihedralData(:,4) == indices(nn, 1))];
       isel_dih3 = [isel_dih3; find(DihedralData(:,5) == indices(nn, 1))];
       isel_dih4 = [isel_dih4; find(DihedralData(:,6) == indices(nn, 1))];
    end
    isel_dih_dum1 = intersect(isel_dih1(:, 1), isel_dih2(:, 1));
    isel_dih_dum2 = intersect(isel_dih_dum1, isel_dih3(:, 1));
    isel_dih = intersect(isel_dih_dum2, isel_dih4(:, 1));
    size_isel_dih = size(isel_dih);
    % Entering data into the sample_DihedralData 
    for kk = 1:size_isel_dih(1, 1)
        atom1 = DihedralData(isel_dih(kk, 1), 3);
        atom2 = DihedralData(isel_dih(kk, 1), 4);
        atom3 = DihedralData(isel_dih(kk, 1), 5);
        atom4 = DihedralData(isel_dih(kk, 1), 6);
        sample_DihedralData{jj} = [sample_DihedralData{jj}; [(num_dihedrals), (DihedralData(isel_dih(kk, 1), 2)), (atom1), (atom2), (atom3), (atom4)]];
        num_dihedrals = num_dihedrals + 1;
    end

 




% sample_ImproperData
num_impropers = 1;

isel_imp1 = [];
isel_imp2 = [];
isel_imp3 = [];
isel_imp4 = [];

for nn = 1:size(indices,1)

    isel_imp1 = [isel_imp1; find(ImproperData(:,3) == indices(nn,1))];
    isel_imp2 = [isel_imp2; find(ImproperData(:,4) == indices(nn,1))];
    isel_imp3 = [isel_imp3; find(ImproperData(:,5) == indices(nn,1))];
    isel_imp4 = [isel_imp4; find(ImproperData(:,6) == indices(nn,1))];

end

isel_imp_dum1 = intersect(isel_imp1,isel_imp2);
isel_imp_dum2 = intersect(isel_imp_dum1,isel_imp3);
isel_imp = intersect(isel_imp_dum2,isel_imp4);

for kk = 1:length(isel_imp)

    atom1 = ImproperData(isel_imp(kk),3);
    atom2 = ImproperData(isel_imp(kk),4);
    atom3 = ImproperData(isel_imp(kk),5);
    atom4 = ImproperData(isel_imp(kk),6);


sample_ImproperData{jj} = [sample_ImproperData{jj}; [(num_impropers), (ImproperData(isel_dih(kk, 1), 2)), (atom1), (atom2), (atom3), (atom4)]];



    % sample_ImproperData{jj} = ...
    %     [sample_ImproperData{jj}; ...
    %     num_impropers, ...
    %     ImproperData(isel_imp(kk),2), ...
    %     atom1, atom2, atom3, atom4];

    num_impropers = num_impropers + 1;

end

end





%---------------------------------------------------------------
% saving new data along with older data in old .mat file
save('lammpsfile.mat','NumAllAtoms', 'NumAllBonds', 'NumAllDihedrals','NumAllImpropers', 'NumAllAngles', 'NumTypeAtoms', 'NumTypeBonds', 'NumTypeDihedrals','NumTypeImpropers',...
    'NumTypeAngles', 'MassData', 'AtomData', 'AtomDataVirgin', 'MoleculeData', 'NumAllMolecules', 'NumTypeMolecules',  'BondData', 'AngleData', 'DihedralData','ImproperData', ...
    'DihedCoeff', 'NumMolInType', 'BondCoeff','PairCoeff', 'AngleCoeff','ImpCoeff',  'Lx', 'Ly', 'Lz',...
    'sample_Molecules', 'sample_BondData', 'sample_AngleData', 'sample_DihedralData', 'sample_ImproperData', 'Pattern');
end

