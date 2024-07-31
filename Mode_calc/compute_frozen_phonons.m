% Get the current directory path
currentDirPath = pwd;

% Split the path into parts
pathParts = strsplit(currentDirPath, '/');

% Extract the last part of the path, which is the current directory name
IRREP = pathParts{end};
% Construct the working directory path
workingDirectory = strcat(currentDirPath, '/');
% Get the current directory
%workingDirectory = strcat(pwd,'/',IRREP,'/');
% Load the forces on the unperturbed atoms
A_ifcs_0 = load(strcat(workingDirectory,'forces_0.dat'));
% Load the displacement magnitude (Cartesian)
load(strcat(workingDirectory,'Udis.dat'));
SGN=+1; %sign of displacement
% Load the forces on the perturbed atoms
A_ifcs_1 = load(strcat(workingDirectory,'forces_',IRREP,'.dat'));
% Load the masses associated with each SAM mode
load(strcat(workingDirectory,'mode_masses'));
MASS = mode_masses;
BASIS = load(strcat(workingDirectory,'SAMs_',IRREP,'.txt'));
% Read the number of SAMs
NumIons=length(MASS);
% Get the number of atoms from the length of the SAM basis set
num_atoms = length(BASIS)/NumIons;
% Format the basis correctly
BASIS = permute(reshape(BASIS', 3, num_atoms, NumIons), [2, 1, 3]);
F00=A_ifcs_0(:,[4 5 6]); F0=F00;
for i=1:(NumIons-1),
  F0=[F00;F0];
end
F11=A_ifcs_1(:,[4 5 6])-F0;
for i=1:NumIons,
  F1(:,:,i)=F11(1+(i-1)*num_atoms:num_atoms+(i-1)*num_atoms,:);
end
for i=1:NumIons,
  for j=1:NumIons,
    FFF(j,i)= sum(sum(BASIS(:,:,j).*F1(:,:,i)));
  end
end
e=1.602177E-19;         % Coloumbs
perm=8.854187E-12;   % C^2/(J m)
c = 3e10;                % cm/s speed of light
AMU = 1.66053E-27;   % kg/amu
CONST=e/(1E-10)^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%a=lattPARvasp(lattFILE);
%[a,Vol]=lattPARvasp(lattFILE);
a=3.7827485118871818;
%a=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumIons=length(MASS);
NumMods=1*NumIons;
M=ones(length(MASS))*diag(sqrt(MASS));
MM=M.*M';
D=-FFF./(SGN*a*Udis);
D=(D+D')./2;%force contants matrix         % D(1:2:5,1:2:5)
DTO=D./MM; %dynamical matrix
[vTO,dTO]=eig(DTO);
Veig_TO=vTO;
Vreal_AU_TO=vTO./M';
wTO= (dTO.*CONST./AMU).^.5;
muTO=diag(wTO/(2*pi*c));
VTO=reshape(Vreal_AU_TO,1,NumIons,NumMods);
VTOeig=reshape(Veig_TO,1,NumIons,NumMods);
[vf,vd]=eig(D);
%vEu1 are eigenvectors of force constants matrix, rEu1 are eigenvectors of
%dynamical matrix
vEu1=0*ones(size(BASIS));
rEu1=0*ones(size(BASIS));
eigEu1=0*ones(size(BASIS));
for j=1:NumIons
  for i=1:NumIons,
    vEu1(:,:,j)=vEu1(:,:,j)+vf(i,j)*BASIS(:,:,i);
    rEu1(:,:,j)=rEu1(:,:,j)+VTO(:,i,j)*BASIS(:,:,i);
    eigEu1(:,:,j)=eigEu1(:,:,j)+VTOeig(:,i,j)*BASIS(:,:,i);
  end
end
% Get wTO in units cm-1
wTO_thz = wTO/(2*pi*10^12);
wTO_cm = wTO_thz * 33.35641;
dlmwrite(strcat(workingDirectory,'phonon_freqs_THz'), wTO_thz, ' ')
dlmwrite(strcat(workingDirectory,'phonon_freqs_cm-1'), wTO_cm, ' ')
dlmwrite(strcat(workingDirectory,'force_consts'),vd, ' ')
% Remove the old 'force_const_eigenvecs' file if it exists
forceConstEigenvecsFile = strcat(workingDirectory, 'force_const_eigenvecs');
if exist(forceConstEigenvecsFile, 'file')
    delete(forceConstEigenvecsFile);
end
for i = 1:NumIons
  dlmwrite(strcat(workingDirectory,'force_const_eigenvecs'), vEu1(:,:,i), ' ', "-append", "newline", "unix", "precision", 16)
  dlmwrite(strcat(workingDirectory,'force_const_eigenvecs'), [], ' ', "-append", "newline", "unix", "precision", 16)
end
