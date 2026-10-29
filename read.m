function[] = read(Filenam)
%read_lammps("")
AtomDataVirgin=[];    
MassData=[];AtomData=[];BondData=[];AngleData=[];DihedralData=[];
DihedCoeff=[]; BondCoeff=[]; AngleCoeff=[]; PairCoeff=[]; ImproperData=[];ImpCoeff=[];
NumTypeAtoms=[];
fid=fopen(Filenam);

tline=fgets(fid);
while tline~=-1
 
    if contains(tline,' atoms') && ~contains(tline,'atom types')
        disp(tline)
        NumAllAtoms = sscanf(tline,'%d',1);
    end
    
    if contains(tline,'atom types')
        NumTypeAtoms = sscanf(tline,'%d',1);
        fprintf('NumTypeAtoms = %d\n',NumTypeAtoms);
    end



    if contains(tline,' bonds') && ~contains(tline,'bond types')
        NumAllBonds = sscanf(tline,'%d',1);
    end

    if contains(tline,'bond types')
        NumTypeBonds = sscanf(tline,'%d',1);
    end


    if contains(tline,' angles') && ~contains(tline,'angle types')
        NumAllAngles = sscanf(tline,'%d',1);
    end

    if contains(tline,'angle types')
        NumTypeAngles = sscanf(tline,'%d',1);
    end


    if contains(tline,' dihedrals') && ~contains(tline,'dihedral types')
        NumAllDihedrals = sscanf(tline,'%d',1);
    end



    if contains(tline,'dihedral types')
        NumTypeDihedrals = sscanf(tline,'%d',1);
    end


    if contains(tline,' impropers') && ~contains(tline,'improper types')
        NumAllImpropers = sscanf(tline,'%d',1);
    end

    if contains(tline,'improper types')
        NumTypeImpropers = sscanf(tline,'%d',1);
    end


    if ismember('xlo',tline)==1
        k = findstr(tline,'xlo');
        if size(k)>0
        
        Lx=str2num(tline(1:k-1));
        end
    end
    if ismember('ylo',tline)==1
        k = findstr(tline,'ylo');
        if size(k)>0
        
        Ly=str2num(tline(1:k-1))
        end
    end
    if ismember('zlo',tline)==1
        k = findstr(tline,'zlo');
        if size(k)>0
        
        Lz=str2num(tline(1:k-1))
        end
    end


    if ismember('Masses',tline)==1
        tline=fgets(fid);
        for i=1:NumTypeAtoms
            tline=fgets(fid);
             s=str2num(tline);
            MassData=[MassData;[s]];
        end
    end

    if ismember('Pair Coeffs',tline)==1
        tline=fgets(fid);
        for i=1:NumTypeAtoms
            tline=fgets(fid);
            s=str2num(tline);
            PairCoeff=[PairCoeff;[str2num(tline)]];
        end
    end

    if ismember('Bond Coeffs',tline)
        tline=fgets(fid);

        for i=1:NumTypeBonds
            tline=fgets(fid);
            s=str2num(tline);
            BondCoeff=[BondCoeff;[s]];
        end
    end

    if ismember('Angle Coeffs',tline)
        tline=fgets(fid);

        for i=1:NumTypeAngles
            tline=fgets(fid);
            s=str2num(tline);
            AngleCoeff=[AngleCoeff;[s]];
        end

    end

    if ismember('Dihed Coeffs',tline)
        tline=fgets(fid);
        for i=1:NumTypeDihedrals
            tline=fgets(fid);
            s=str2num(tline);
            DihedCoeff=[DihedCoeff;[s]];
        end
    end

    if ismember('Improper Coeffs',tline)
        tline=fgets(fid);
        for i=1:NumTypeImpropers
            tline=fgets(fid);
            s=str2num(tline);
            ImpCoeff=[ImpCoeff;[s]];
        end
    end


    if ismember('Atoms',tline)
        tline=fgets(fid);

        for i=1:NumAllAtoms
            tline=fgets(fid);
            AtomDataVirgin=[AtomDataVirgin;[str2num(tline)]];
            AtomData=[AtomData;[str2num(tline)]];
        end
    end
    % 
    % 
    if ismember('Bonds',tline)
        tline=fgets(fid);

        for i=1:NumAllBonds
            tline=fgets(fid);
            BondData=[BondData;[str2num(tline)]];
        end
    end

    if ismember('Angles',tline)
        tline=fgets(fid);
        for i=1:NumAllAngles
            tline=fgets(fid);
            AngleData=[AngleData;[str2num(tline)]];
        end
    end

    if ismember('Dihedrals',tline)
        tline=fgets(fid);
        for i=1:NumAllDihedrals
            tline=fgets(fid);
            DihedralData=[DihedralData;[str2num(tline)]];
        end
    end
    
    if ismember('Impropers',tline)
        tline=fgets(fid);
        for i=1:NumAllImpropers
            tline=fgets(fid);
            ImproperData=[ImproperData;[str2num(tline)]];
        end
    end





    tline=fgets(fid);
    
end


disp1=['Read ' num2str(NumAllAtoms) ' atoms of ' num2str(NumTypeAtoms) ' types'];
disp2=['Read ' num2str(NumAllBonds) ' bondss of ' num2str(NumTypeBonds) ' types'];
disp3=['Read ' num2str(NumAllAngles) ' angles of ' num2str(NumTypeAngles) ' types'];
disp4=['Read ' num2str(NumAllDihedrals) ' dihedrals of ' num2str(NumTypeDihedrals) ' types'];
disp5=['Read ' num2str(NumAllImpropers) ' Impropers of ' num2str(NumTypeImpropers) ' types'];
disp(disp1);
disp(disp2);
disp(disp3);
disp(disp4);
disp(disp5);
disp6=['Box size goes from (' num2str([Lx(1) Ly(1) Lz(1)]) ') to (' num2str([Lx(2) Ly(2) Lz(2)]) ')'];
disp(disp6);
disp('Invoke set_parameters_plot_initial();');
save('lammpsfile.mat','NumAllAtoms', 'NumAllBonds', 'NumAllDihedrals', 'NumAllImpropers', 'NumAllAngles', 'NumTypeAtoms', 'NumTypeBonds', 'NumTypeDihedrals','NumTypeImpropers',...
    'NumTypeAngles', 'MassData', 'AtomData', 'AtomDataVirgin', 'BondData', 'AngleData', 'DihedralData','ImproperData',...
    'DihedCoeff', 'BondCoeff', 'AngleCoeff','ImpCoeff', 'PairCoeff', 'Lx', 'Ly', 'Lz');





%===========================================================================
% END OF READING INPUT FILE
%===========================================================================
