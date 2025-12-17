clear
clc
close all

%% Parameters
h    = 0.325e-4;          % spacing (cm)
dh   = 0.325e-4;          % kept for reference
zeta = 1.0;               % SBM interface thickness (nondim here)

sets = '3D_disk_full';

% --- Rectangular mesh (feel free to tweak) ---
me = 80;   % y-elements
ne = 120;  % x-elements (wider -> rectangular)
le = 6;    % z-elements

mn = me + 1;
nn = ne + 1;
ln = le + 1;

% --- Allocate distance fields for A and C (same size as old psi/dsf) ---
dsf_A = zeros(mn, nn, ln);   % anode distance field
dsf_C = zeros(mn, nn, ln);   % cathode distance field

% --- Disk parameters (same shape as before) ---
radius = 25;                 % in index units (matches your earlier code)
k0 = round(mn/2);            % vertical center (y-like index)
jA = 22.5;                   % left-side center for anode
jC = nn - 22.5;              % right-side center for cathode (mirrored)

% Build distance fields: + inside, - outside (radius - r)
for iY = 1:mn
    for jX = 1:nn
        rA = sqrt((iY - k0)^2 + (jX - jA)^2);
        rC = sqrt((iY - k0)^2 + (jX - jC)^2);
        dsf_A(iY, jX, :) = radius - rA;
        dsf_C(iY, jX, :) = radius - rC;
    end
end

% --- Smoothed indicator functions (SBM) ---
psA = 0.5 * (1 + tanh(dsf_A(:,:,1) / zeta));  % anode mask
psC = 0.5 * (1 + tanh(dsf_C(:,:,1) / zeta));  % cathode mask

figure; imagesc(psA); axis equal tight; colorbar; title('psA (anode)');
figure; imagesc(psC); axis equal tight; colorbar; title('psC (cathode)');

%% Element connectivity (CUBE elements like before)
eid = 0;
for k = 1:le
    for i = 1:me
        for j = 1:ne
            eid = eid + 1;
            ERC(eid, 1:3) = [i j k];
            ERC(eid, 4) = eid - 1;
        end
    end
end
nE = size(ERC, 1);

%% Vertex table (x-fastest when writing later)
nid = 0;
for k = 1:ln
    for i = 1:mn
        for j = 1:nn
            nid = nid + 1;
            N_ps(nid, 1) = (j - 1) * h;  % x (column)
            N_ps(nid, 2) = (i - 1) * h;  % y (row)
            N_ps(nid, 3) = (k - 1) * h;  % z (height)
            N_ps(nid, 4) = nid - 1;      % label
            NRC(nid, 1:3) = [i j k];
            NRC(nid, 4) = nid - 1;
        end
    end
end
nV = size(N_ps, 1);

%% Hex element node ordering (same as before)
for el = 1:nE
    nd1 = (ERC(el,1) - 1) * nn + (ERC(el,2) - 1) + (ERC(el,3) - 1) * (mn * nn);
    nd2 = nd1 + 1;
    nd3 = nd2 + nn;
    nd4 = nd1 + nn;
    nd5 = nd1 + mn * nn;
    nd6 = nd2 + mn * nn;
    nd7 = nd3 + mn * nn;
    nd8 = nd4 + mn * nn;
    IEN(el, 1:8) = [nd1 nd2 nd3 nd4 nd5 nd6 nd7 nd8];
end

%% element table
for el = 1:nE
    Efile(el, 1) = 1;         % attribute
    Efile(el, 2) = 5;         % 5 = CUBE
    Efile(el, 3:10) = IEN(el, 1:8);
end

%% boundary table (west/east/south/north/top/bottom)
cnt = 1;
for el = 1: nE
    % south (i == 1)
    if ERC(el,1) == 1
        Bndry(cnt,1:6) = [2 3 IEN(el,1) IEN(el,2) IEN(el,6) IEN(el,5)]; cnt = cnt+1;
    end
    % north (i == me)
    if ERC(el,1) == me
        Bndry(cnt,1:6) = [4 3 IEN(el,4) IEN(el,3) IEN(el,7) IEN(el,8)]; cnt = cnt+1;
    end
    % west  (j == 1)
    if ERC(el,2) == 1
        Bndry(cnt,1:6) = [1 3 IEN(el,1) IEN(el,5) IEN(el,8) IEN(el,4)]; cnt = cnt+1;
    end
    % east  (j == ne)
    if ERC(el,2) == ne
        Bndry(cnt,1:6) = [3 3 IEN(el,2) IEN(el,6) IEN(el,7) IEN(el,3)]; cnt = cnt+1;
    end
    % top   (k == 1)
    if ERC(el,3) == 1
        Bndry(cnt,1:6) = [5 3 IEN(el,1) IEN(el,2) IEN(el,3) IEN(el,4)]; cnt = cnt+1;
    end
    % bottom (k == le)
    if ERC(el,3) == le
        Bndry(cnt,1:6) = [6 3 IEN(el,5) IEN(el,6) IEN(el,7) IEN(el,8)]; cnt = cnt+1;
    end
end
Bndry = sortrows(Bndry,1);

%% MFEM mesh output
dim = 3;
fname = ['Mesh_' num2str(me) 'x' num2str(ne) 'x' num2str(le) '_' sets '.mesh'];

stn = {
'MFEM mesh v1.0'
''
'#'
'# MFEM Geometry Types (see mesh/geom.hpp):'
'#'
'# POINT       = 0'
'# SEGMENT     = 1'
'# TRIANGLE    = 2'
'# SQUARE      = 3'
'# TETRAHEDRON = 4'
'# CUBE        = 5'
'# PRISM       = 6'
'# PYRAMID     = 7'
'#'
''};

fid = fopen(fname, 'w');
for k = 1:numel(stn)
    fprintf(fid, '%s\n', stn{k});
end
fprintf(fid, 'dimension\n%d\n\n', dim);

fprintf(fid, 'elements\n%d\n', size(Efile, 1));
for k = 1:size(Efile, 1)
    fprintf(fid, '%d %d %d %d %d %d %d %d %d %d\n', Efile(k, :));
end
fprintf(fid, '\n');

fprintf(fid, 'boundary\n%d\n', size(Bndry, 1));
for k = 1:size(Bndry, 1)
    fprintf(fid, '%d %d %d %d %d %d\n', Bndry(k, :));
end
fprintf(fid, '\n');

fprintf(fid, 'vertices\n%d\n%d\n', size(N_ps, 1), dim);
for k = 1:size(N_ps, 1)
    fprintf(fid, '%.7f %.7f %.7f\n', N_ps(k, 1:3));
end
fprintf(fid, '\n');
fclose(fid);

%% Distance field exports (x-fastest: i -> x, j -> y, k -> z)
fnameA = ['dsF_A_' num2str(me) 'x' num2str(ne) 'x' num2str(le) '_' sets '.txt'];
fid = fopen(fnameA, 'w');
for k = 1:ln          % Z
    for j = 1:mn      % Y
        for i = 1:nn  % X (fastest)
            fprintf(fid, '%.16e\n', dsf_A(j, i, k));  % note swap to (y,x,z)
        end
    end
end
fclose(fid);

fnameC = ['dsF_C_' num2str(me) 'x' num2str(ne) 'x' num2str(le) '_' sets '.txt'];
fid = fopen(fnameC, 'w');
for k = 1:ln
    for j = 1:mn
        for i = 1:nn
            fprintf(fid, '%.16e\n', dsf_C(j, i, k));  % note swap to (y,x,z)
        end
    end
end
fclose(fid);

%% (Optional) quick sanity check of ranges
fprintf('[psA] min/max: %.6f / %.6f\n', min(psA,[],'all'), max(psA,[],'all'));
fprintf('[psC] min/max: %.6f / %.6f\n', min(psC,[],'all'), max(psC,[],'all'));


% --- Overlay psA & psC ---
figure;
imagesc(psA); 
hold on;
h = imagesc(psC);
set(h, 'AlphaData', 0.5);   % 50% transparency for cathode
axis equal tight; colorbar;
title('Overlay: psA + psC');
xlabel('X'); ylabel('Y');
