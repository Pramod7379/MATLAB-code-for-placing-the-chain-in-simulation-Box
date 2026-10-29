function [] = print_new_data_file()
% Syntax : print_new_data_file()
% This function prints new data generated from the 'rand_allo()' function
% to create a new input (text) file.
% CAUTION : This function is to be used only after rand_allo().
% For more details, refer to the manual provided with the package.
if exist('newlammpsfile.mat')
    load 'newlammpsfile.mat';
else
    disp('Run Sort_molecules(), then rand_allo(), and after that use this function.');
end
%----------------------------------------------------------------------------
% Printing code
fileID = fopen('new_Data.txt','w');
% Name of file
fprintf(fileID,'%s\r\n','NEW RANDOMIZED DATA');
fprintf(fileID,'%s\r\n',' ');
% Number of each item
fprintf(fileID,'%6d %12s\r\n',new_NumAllAtoms,'atoms');
fprintf(fileID,'%6d %12s\r\n',new_NumAllBonds,'bonds');
fprintf(fileID,'%6d %12s\r\n',new_NumAllAngles,'angles');
fprintf(fileID,'%6d %12s\r\n',new_NumAllDihedrals,'dihedrals');
fprintf(fileID,'%6d %12s\r\n',new_NumAllImpropers,'impropers');

fprintf(fileID,'%s\r\n',' ');
% Types
fprintf(fileID,'%6d %12s\r\n',new_NumTypeAtoms,'atom types');
fprintf(fileID,'%6d %12s\r\n',new_NumTypeBonds,'bond types');
fprintf(fileID,'%6d %12s\r\n',new_NumTypeAngles,'angle types');
fprintf(fileID,'%6d %12s\r\n',new_NumTypeDihedrals,'dihedral types');
fprintf(fileID,'%6d %12s\r\n',new_NumTypeImpropers,'improper types');

fprintf(fileID,'%s\r\n',' ');
% Box size
fprintf(fileID,'%8f %8f %6s %6s\r\n',new_Lx,'xlo','xhi');
fprintf(fileID,'%8f %8f %6s %6s\r\n',new_Ly,'ylo','yhi');
fprintf(fileID,'%8f %8f %6s %6s\r\n',new_Lz,'zlo','zhi');
fprintf(fileID,'%s\r\n',' ');
% MassData
if isempty(new_MassData) ~= 1
    fprintf(fileID,'%s\r\n','Masses');
    fprintf(fileID,'%s\r\n',' ');
    fprintf(fileID,'%3d %8.5f\r\n',(new_MassData)');
    fprintf(fileID,'%s\r\n',' ');
end

% Angle Coeffs
if isempty(new_PairCoeff) ~= 1    
    fprintf(fileID,'%s\r\n','Pair Coeffs');
    fprintf(fileID,'%s\r\n',' ');
    fprintf(fileID,'%3d %13.9f %13.9f\r\n',(new_PairCoeff)');
    fprintf(fileID,'%s\r\n',' ');
end



% Bond Coeffs
if isempty(new_BondCoeff) ~= 1    
    fprintf(fileID,'%s\r\n','Bond Coeffs');
    fprintf(fileID,'%s\r\n',' ');
    fprintf(fileID,'%3d %14.9f %14.9f\r\n',(new_BondCoeff)');
    fprintf(fileID,'%s\r\n',' ');
end
% Angle Coeffs
if isempty(new_AngleCoeff) ~= 1    
    fprintf(fileID,'%s\r\n','Angle Coeffs');
    fprintf(fileID,'%s\r\n',' ');
    fprintf(fileID,'%3d %13.9f %13.9f\r\n',(new_AngleCoeff)');
    fprintf(fileID,'%s\r\n',' ');
end
% Dihedral Coeffs
if isempty(new_DihedCoeff) ~= 1    
    fprintf(fileID,'%s\r\n','Dihedral Coeffs');
    fprintf(fileID,'%s\r\n',' ');
    fprintf(fileID,'%3d %13.9f %6d %6d\r\n',(new_DihedCoeff)');
    fprintf(fileID,'%s\r\n',' ');
    
    % Angle Coeffs
    if isempty(new_ImpCoeff) ~= 1    
        fprintf(fileID,'%s\r\n','Improper Coeffs');
        fprintf(fileID,'%s\r\n',' ');
        fprintf(fileID,'%3d %13.9f %13.9f\r\n',(new_ImpCoeff)');
        fprintf(fileID,'%s\r\n',' ');
    end
    
    
end
% Atoms
if isempty(new_AtomData) ~= 1    
    fprintf(fileID,'%s\r\n','Atoms');
    fprintf(fileID,'%s\r\n',' ');
    fprintf(fileID,'%6d %6d %6d %8.6f %14.9f %14.9f %14.9f\r\n',(new_AtomData)');
    fprintf(fileID,'%s\r\n',' ');
end
% Bonds
if isempty(new_BondData) ~= 1    
    fprintf(fileID,'%s\r\n','Bonds');
    fprintf(fileID,'%s\r\n',' ');
    fprintf(fileID,'%6d %6d %6d %6d\r\n',(new_BondData)');
    fprintf(fileID,'%s\r\n',' ');
end
% Angles
if isempty(new_AngleData) ~= 1    
    fprintf(fileID,'%s\r\n','Angles');
    fprintf(fileID,'%s\r\n',' ');
    fprintf(fileID,'%6d %6d %6d %6d %6d\r\n',(new_AngleData)');
    fprintf(fileID,'%s\r\n',' ');
end
% Dihedrals
if isempty(new_DihedralData) ~= 1    
    fprintf(fileID,'%s\r\n','Dihedrals');
    fprintf(fileID,'%s\r\n',' ');
    fprintf(fileID,'%6d %6d %6d %6d %6d %6d\r\n',(new_DihedralData)');
    fprintf(fileID,'%s\r\n',' ');
end

% Dihedrals
if isempty(new_ImproperData) ~= 1    
    fprintf(fileID,'%s\r\n','Impropers');
    fprintf(fileID,'%s\r\n',' ');
    fprintf(fileID,'%6d %6d %6d %6d %6d %6d\r\n',(new_ImproperData)');
    fprintf(fileID,'%s\r\n',' ');
end




fclose(fileID);
end

