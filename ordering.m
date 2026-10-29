% sort atom bond,angle.... improper data according to id #

function [] = ordering()
if exist('lammpsfile.mat')
    load 'lammpsfile.mat';
else
    disp('Please upload lammps data before attempting to set parameters');
end

%%%%%%%%%%%%%%%%%%%%%%%%
AtomData = sortrows(AtomData,1);
AtomDataVirgin = sortrows(AtomDataVirgin,1);
BondData = sortrows(BondData,1);
AngleData = sortrows(AngleData,1);
DihedralData=sortrows(DihedralData,1);
ImproperData=sortrows(ImproperData,1);

save('lammpsfile.mat','NumAllAtoms', 'NumAllBonds', 'NumAllDihedrals', 'NumAllImpropers', 'NumAllAngles', 'NumTypeAtoms', 'NumTypeBonds', 'NumTypeDihedrals','NumTypeImpropers',...
    'NumTypeAngles', 'MassData', 'AtomData', 'AtomDataVirgin', 'BondData', 'AngleData', 'DihedralData','ImproperData',...
    'DihedCoeff', 'BondCoeff', 'AngleCoeff','ImpCoeff', 'PairCoeff', 'Lx', 'Ly', 'Lz');
end