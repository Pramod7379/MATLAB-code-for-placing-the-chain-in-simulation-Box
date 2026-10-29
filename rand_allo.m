function [] = rand_allo(box_size, num_mol_of_type)
% Syntax : rand_allo(box_size, num_mol_of_type)
% This function randomly allocates a position and orientation to the molecules selected
% from each type, given by 'num_mol_of_type', within a new box, whose
% dimensions are given by 'box_size'.
% 'box_size' should be a 3 - element array, with inputs as new box's dimensions.
% 'num_mol_of_type' should contain as many no. of elements, as there are
% types of molecules, with 0 as an entry if no molecules are to be selected
% from that type.
% For more details, refer to the manual provided with the package.
if exist('lammpsfile.mat')
    load 'lammpsfile.mat';
else
    disp('Please upload lammps data before attempting to set parameters');
end
%----------------------------------------------------------------------------
% Conditions imposed on input arguments
size_box_size = size(box_size);
if (size_box_size(2) ~= 3)
    disp('Error in dimensions of box_size. Input should be a 3-element array, with format [Lx Ly Lz].');
    return;
end
size_num_mol_of_type = size(num_mol_of_type);
if (size_num_mol_of_type(2) ~= NumTypeMolecules)
    disp('Error in dimensions of num_mol_of_type. Input should be an array of size NumTypeMolecules, with 0 as entry if no molecules of that type are required.');
    return;
end
%-------------------------------------------------------------------------------
% Random selection of molecules from MoleculeData
rand_sel_mol = [];
for ii = 1:NumTypeMolecules
    % Selecting random indices
    ind_rand = [];
    if (num_mol_of_type(ii) < 0)
        disp(['The ' num2str(ii) 'th element of the input argument "num_mol_of_type" is invalid.']);
        return;
    elseif (num_mol_of_type(ii) > 0) && (num_mol_of_type(ii) <= NumMolInType(ii, 1))
        ind_rand = randperm(NumMolInType(ii, 1), num_mol_of_type(ii));
    elseif (num_mol_of_type(ii) > NumMolInType(ii, 1))
        a = (num_mol_of_type(ii))/(NumMolInType(ii, 1));
        b = floor(a);
        c = floor((NumMolInType(ii, 1))*(a-b));
        for jj = 1:b
            ind_rand = [ind_rand, randperm(NumMolInType(ii, 1), NumMolInType(ii, 1))];
        end    
        ind_rand = [ind_rand, randperm(NumMolInType(ii, 1), c)];
    end
    % Selecting the molecules from MoleculeData{ii} according to those
    % random indices
    B = unique(MoleculeData{ii}(:, 2));
    size_ind_rand = size(ind_rand);
    if (size_ind_rand(2) ~= 0)
        for kk = 1:size_ind_rand(2)
            rand_sel_mol = [rand_sel_mol; B(ind_rand(kk))];
        end
    end    
end
%--------------------------------------------------------------------------------------------
% Checking whether any molecule selected or not
if (isempty(rand_sel_mol) == 1)
    disp('No molecules selected. Quitting.');
    return;
end
%---------------------------------------------------------------------------------------------
% Taking the old AtomData matrix and modifying it to 'no PBC' condition in
% all directions
Pos_new = [AtomDataVirgin(:,5) AtomDataVirgin(:,6) AtomDataVirgin(:,7)];
L_size = [Lx(2)-Lx(1) Ly(2)-Ly(1) Lz(2)-Lz(1)];
for k=1:3
    for i = 1:(NumAllAtoms - 1)
        if AtomDataVirgin(i, 2) == AtomDataVirgin(i+1, 2)
            dw = AtomDataVirgin(i+1, k+4) - Pos_new(i, k);
            dw = dw - round(dw/L_size(k)) * L_size(k);
            Pos_new(i+1, k) = Pos_new(i, k) + dw;
        elseif AtomDataVirgin(i, 2) ~= AtomDataVirgin(i+1, 2)
            Pos_new(i+1, k) = AtomDataVirgin(i+1, k+4);
        end
    end           
end





% No PBC AtomData
AtomData = [AtomDataVirgin(:, 1:4) Pos_new];
%----------------------------------------------------------------------------------------------
% Random rotation and translation; and creating a new AtomData matrix
size_rand_sel_mol = size(rand_sel_mol);
num_atoms = 1;
new_AtomData = [];
for ii = 1:size_rand_sel_mol(1)
    dum_indices = find(AtomData(:, 2) == rand_sel_mol(ii, 1));
    size_dum_indices = size(dum_indices);
    dum_Pos = AtomData(dum_indices, 5:7);
    Pos = (dum_Pos)';
    % Translating the first atom of every molecule to origin
    Pos_origin = [Pos(1, :)-Pos(1, 1); Pos(2, :)-Pos(2, 1); Pos(3, :)-Pos(3, 1)];
    % Applying rotation matrix; random rotation
    angle_phi = (rand*180) - 90;
    angle_theta = rand*360;
    R_theta = [cosd((angle_theta)), sind((angle_theta)), 0; -sind((angle_theta)), cosd((angle_theta)), 0; 0, 0, 1];
    R_phi = [cosd(angle_phi), 0, sind(angle_phi); 0, 1, 0; -sind(angle_phi), 0, cosd(angle_phi)];
    mod_Pos_origin = (R_theta)*(R_phi)*(Pos_origin);
    % Applying random translation in new box defined by 'box_size'
    rand_x = rand*(box_size(1));
    rand_y = rand*(box_size(2));
    rand_z = rand*(box_size(3)/3);
  % rand_z=0;
    dum_mod_Pos = [(mod_Pos_origin(1, :)+rand_x); (mod_Pos_origin(2, :)+rand_y); (mod_Pos_origin(3, :)+rand_z)];
    mod_Pos = (dum_mod_Pos)';
    % Filling up new_AtomData matrix
    for jj = 1:size_dum_indices(1)
        new_AtomData = [new_AtomData; [num_atoms, (ii-1), AtomData((dum_indices(jj, 1)), 3), AtomData((dum_indices(jj, 1)), 4), mod_Pos(jj, :)]];
        num_atoms = num_atoms + 1;
    end
end
%---------------------------------------------------------------------------------------------------------------------------------------------------
% New NumAllAtoms
new_NumAllAtoms = num_atoms - 1;
%---------------------------------------------------------------------------------------------------------------------------------------------------
% Applying PBC to new_AtomData

% Applying PBC to new_AtomData
  % Len=83.78;
  % for kk=1:3
  %      for ii = 1:new_NumAllAtoms
  % 
  %             new_AtomData(ii, kk+4)=new_AtomData(ii, kk+4)-floor(( new_AtomData(ii, kk+4))/Len)*Len;       
  %      end
  %  end


%---------------------------------------------------------------------------------------------------------------------------------------------------
% New BondData, AngleData and DihedralData
new_BondData = [];
new_AngleData = [];
new_DihedralData = [];
new_ImproperData = [];
num_bonds = 1;
num_angles = 1;
num_dihedrals = 1;
num_impropers=1;
for ii = 1:size_rand_sel_mol(1)
    
    % Finding the type of the old molecule
    old_indices = find(AtomData(:, 2) == rand_sel_mol(ii, 1));
    dum_Pattern = (AtomData(old_indices, 3))';
    for jj = 1:NumTypeMolecules
        if (size(dum_Pattern) == size(Pattern{jj}))
            if dum_Pattern == Pattern{jj}
                Type_of_molecule = jj;
            end
        end
    end

    % Finding the molecule in the new_AtomData
    new_indices = find(new_AtomData(:, 2) == (ii-1));
    
    % new_BondData
    size_sample_BondData = size(sample_BondData{Type_of_molecule});
    for kk = 1:size_sample_BondData(1)
        dum_atom1 = sample_BondData{Type_of_molecule}(kk, 3);
        dum_atom2 = sample_BondData{Type_of_molecule}(kk, 4);
        dum_index1 = find(sample_Molecules{Type_of_molecule}(:, 1) == dum_atom1);
        dum_index2 = find(sample_Molecules{Type_of_molecule}(:, 1) == dum_atom2);
        atom1 = new_AtomData(new_indices(dum_index1), 1);
        atom2 = new_AtomData(new_indices(dum_index2), 1);
        new_BondData = [new_BondData; [(num_bonds), sample_BondData{Type_of_molecule}(kk, 2), (atom1), (atom2)]];
        num_bonds = num_bonds + 1;
    end
    
    % new_AngleData
    size_sample_AngleData = size(sample_AngleData{Type_of_molecule});
    for kk = 1:size_sample_AngleData(1)
        dum_atom1 = sample_AngleData{Type_of_molecule}(kk, 3);
        dum_atom2 = sample_AngleData{Type_of_molecule}(kk, 4);
        dum_atom3 = sample_AngleData{Type_of_molecule}(kk, 5);
        dum_index1 = find(sample_Molecules{Type_of_molecule}(:, 1) == dum_atom1);
        dum_index2 = find(sample_Molecules{Type_of_molecule}(:, 1) == dum_atom2);
        dum_index3 = find(sample_Molecules{Type_of_molecule}(:, 1) == dum_atom3);
        atom1 = new_AtomData(new_indices(dum_index1), 1);
        atom2 = new_AtomData(new_indices(dum_index2), 1);
        atom3 = new_AtomData(new_indices(dum_index3), 1);
        new_AngleData = [new_AngleData; [(num_angles), sample_AngleData{Type_of_molecule}(kk, 2), (atom1), (atom2), (atom3)]];
        num_angles = num_angles + 1;
    end
    
    % new_DihedralData
    size_sample_DihedralData = size(sample_DihedralData{Type_of_molecule});
    for kk = 1:size_sample_DihedralData(1)
        dum_atom1 = sample_DihedralData{Type_of_molecule}(kk, 3);
        dum_atom2 = sample_DihedralData{Type_of_molecule}(kk, 4);
        dum_atom3 = sample_DihedralData{Type_of_molecule}(kk, 5);
        dum_atom4 = sample_DihedralData{Type_of_molecule}(kk, 6);
        dum_index1 = find(sample_Molecules{Type_of_molecule}(:, 1) == dum_atom1);
        dum_index2 = find(sample_Molecules{Type_of_molecule}(:, 1) == dum_atom2);
        dum_index3 = find(sample_Molecules{Type_of_molecule}(:, 1) == dum_atom3);
        dum_index4 = find(sample_Molecules{Type_of_molecule}(:, 1) == dum_atom4);
        atom1 = new_AtomData(new_indices(dum_index1), 1);
        atom2 = new_AtomData(new_indices(dum_index2), 1);
        atom3 = new_AtomData(new_indices(dum_index3), 1);
        atom4 = new_AtomData(new_indices(dum_index4), 1);
        new_DihedralData = [new_DihedralData; [(num_dihedrals), sample_DihedralData{Type_of_molecule}(kk, 2), (atom1), (atom2), (atom3), (atom4)]];
        num_dihedrals = num_dihedrals + 1;
    end

% new_ImproperData
size_sample_ImproperData = size(sample_ImproperData{Type_of_molecule});
for kk = 1:size_sample_ImproperData(1)
    dum_atom1 = sample_ImproperData{Type_of_molecule}(kk, 3);
    dum_atom2 = sample_ImproperData{Type_of_molecule}(kk, 4);
    dum_atom3 = sample_ImproperData{Type_of_molecule}(kk, 5);
    dum_atom4 = sample_ImproperData{Type_of_molecule}(kk, 6);
    dum_index1 = find(sample_Molecules{Type_of_molecule}(:, 1) == dum_atom1);
    dum_index2 = find(sample_Molecules{Type_of_molecule}(:, 1) == dum_atom2);
    dum_index3 = find(sample_Molecules{Type_of_molecule}(:, 1) == dum_atom3);
    dum_index4 = find(sample_Molecules{Type_of_molecule}(:, 1) == dum_atom4);
    atom1 = new_AtomData(new_indices(dum_index1), 1);
    atom2 = new_AtomData(new_indices(dum_index2), 1);
    atom3 = new_AtomData(new_indices(dum_index3), 1);
    atom4 = new_AtomData(new_indices(dum_index4), 1);
    new_ImproperData = [new_ImproperData; [(num_impropers), sample_ImproperData{Type_of_molecule}(kk, 2), (atom1), (atom2), (atom3), (atom4)]];
    num_impropers = num_impropers + 1;
end
end







%-------------------------------------------------------------------------------------------------------------------------------------------
% New NumAllBonds, NumAllAngles, NumAllDihedrals
new_NumAllBonds = num_bonds - 1;
new_NumAllAngles = num_angles - 1;
new_NumAllDihedrals = num_dihedrals - 1;
new_NumAllImpropers = num_impropers - 1;

%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Storing all other relevant data 
new_Lx = [0 (box_size(1))];
new_Ly = [0 (box_size(2))];
new_Lz = [0 (box_size(3))];
new_NumTypeAngles = NumTypeAngles;
new_NumTypeAtoms = NumTypeAtoms; 
new_NumTypeBonds = NumTypeBonds; 
new_NumTypeDihedrals = NumTypeDihedrals; 
new_NumTypeImpropers = NumTypeImpropers; 
new_MassData = MassData;
new_DihedCoeff = DihedCoeff; 
new_BondCoeff = BondCoeff; 
new_AngleCoeff = AngleCoeff; 
new_ImpCoeff = ImpCoeff;
new_PairCoeff=PairCoeff;
%---------------------------------------------------------------------------------------------------------------------------------------------------
% Saving new data
save('newlammpsfile.mat','new_NumAllAtoms', 'new_AtomData', 'new_NumAllBonds', 'new_NumAllImpropers','new_BondData',...
    'new_NumAllAngles', 'new_AngleData', 'new_NumAllDihedrals', 'new_DihedralData','new_ImproperData',...
    'new_NumTypeAngles', 'new_NumTypeAtoms', 'new_NumTypeBonds', 'new_NumTypeDihedrals', 'new_NumTypeImpropers', 'new_MassData',...
    'new_DihedCoeff', 'new_BondCoeff','new_ImpCoeff','new_PairCoeff', 'new_AngleCoeff', 'new_Lx', 'new_Ly', 'new_Lz');
end

