clear
clc
close all

%% plot FEM element

% row and column of elements
me = 40;
ne = 30;

sets = 2;

h = 0.2e-4;

% lines of points
mn = me+1;      % # of rows
nn = ne+1;      % # of columns
% 
for i = 1: mn
    plot([1,nn],[i,i],'k-','LineWidth',2)
    hold on
end

for i = 1: nn
    plot([i,i],[1,mn],'k-','LineWidth',2)
end

eid = 0;
for i = 1: me
    for j = 1: ne

        v = num2str(eid);
        text(j+0.3,i+0.4,['$\fbox{' v '}$'], ...
            'Interpreter', 'latex','FontSize',16,'Color','r')

        eid = eid + 1;
        ERC(eid,1:2) = [i j];
        ERC(eid,3) = eid-1;
    end
end

% row and column of nodes
nid = 0;
for i = 1: mn
    for j = 1: nn

        v = num2str(nid);
        text(j-0.25,i+0.1,['$\textcircled{' v '}$'], ...
            'Interpreter', 'latex','FontSize',18,'Color','b')

        nid = nid + 1;

        N_ps(nid,2) = (j-1)*h;  % column
        N_ps(nid,1) = (i-1)*h;  % row
        N_ps(nid,3) = nid - 1;

        % N_ps(nid,2) = (j-1);  % column
        % N_ps(nid,1) = (i-1);  % row
        % N_ps(nid,3) = nid - 1;

        NRC(nid,1:2) = [i,j];
        NRC(nid,3) = nid - 1;

    end
end

% ze = 0.2*h;
% ze = 5*h;
ze = 1.0e-5;

% xc = 12.1e-4;
% xc = ((nn - 1) * h) / 2;
% xc = h*nn/2 + 15*h;
xc = 3.5e-4;


for i = 1: size(N_ps,1)

    dis(i,1) = +(N_ps(i,2) - xc);
    % dis(i,1) = +(N_ps(i,2));


    psi(i) = 0.5*(1+tanh(dis(i,1)/ze));

    dis2(NRC(i,1),NRC(i,2)) = dis(i,1);
    psi2(NRC(i,1),NRC(i,2)) = psi(i);

end

figure()

imagesc(psi2);
colorbar;
title('Psi Distribution');

figure()

imagesc(dis2);
colorbar;
title('distance');

% figure()
% 
% hold on
% set(gcf,'color','w')
% hold off

for i = 1: me*ne

    nd1 = (ERC(i,1)-1)*nn + (ERC(i,2)-1) ;

    nd2 = nd1 + nn;
    nd3 = nd2 + 1;
    nd4 = nd1 + 1;

    IEN(i,1:4) = [nd1 nd4 nd3 nd2]; % changed this to show ccw
    [i-1 nd1 nd4 nd3 nd2];
end

% IEN

%% element table
for el = 1: me*ne

    Efile(el,1) = 1;

    Efile(el,2) = 3;
    Efile(el,3:6) = IEN(el,1:4);

end


%% boundary table

cnt = 1;
for el = 1: me*ne

    % first row of elements (south)
    if ERC(el,1) == 1
        Bndry(cnt,1:4) = [1 1 IEN(el,1) IEN(el,2)];
        cnt = cnt+1;
    end

    % last row of elements (north)
    if ERC(el,1) == me
        Bndry(cnt,1:4) = [3 1 IEN(el,4) IEN(el,3)];
        cnt = cnt+1;
    end

    % first column of elements (west)
    if ERC(el,2) == 1
        Bndry(cnt,1:4) = [4 1 IEN(el,1) IEN(el,4)];
        cnt = cnt+1;
    end

    % last column of elements (east)
    if ERC(el,2) == ne
        Bndry(cnt,1:4) = [2 1 IEN(el,2) IEN(el,3)];
        cnt = cnt+1;
    end
end
size(Bndry)
Bndry = sortrows(Bndry,1);



%% output

dim = 2;
fname = ['Mesh_' num2str(me) 'x' num2str(ne) '_' num2str(sets) '.mesh'];

stn{1} =  'MFEM mesh v1.0';
stn{2} =  '';
stn{3} =  '#';
stn{4} =  '# MFEM Geometry Types (see mesh/geom.hpp):';
stn{5} =  '#';
stn{6} =  '# POINT       = 0';
stn{7} =  '# SEGMENT     = 1';
stn{8} =  '# TRIANGLE    = 2';
stn{9} =  '# SQUARE      = 3';
stn{10} = '# TETRAHEDRON = 4';
stn{11} = '# CUBE        = 5';
stn{12} = '# PRISM       = 6';
stn{13} = '# PYRAMID     = 7';
stn{14} = '#';
stn{15} = '';

fid = fopen(fname,'w');

for k = 1: 15
    fprintf(fid,'%s\n',stn{k});
end
fprintf(fid,'%s\n','dimension');    
fprintf(fid,'%s\n',num2str(dim));
fprintf(fid,'%s\n','');

fprintf(fid,'%s\n','elements');    
fprintf(fid,'%s\n',num2str(size(Efile,1)));
for k = 1: size(Efile,1)
    fprintf(fid,'%d %d %d %d %d %d\n',Efile(k,:));
end
fprintf(fid,'%s\n','');

fprintf(fid,'%s\n','boundary');    
fprintf(fid,'%s\n',num2str(size(Bndry,1)));
for k = 1: size(Bndry,1)
    fprintf(fid,'%d %d %d %d\n',Bndry(k,:))  ;
end
fprintf(fid,'%s\n','');

fprintf(fid,'%s\n','vertices');  
fprintf(fid,'%s\n',num2str(size(N_ps,1)));
fprintf(fid,'%s\n',num2str(dim));
for k = 1: size(N_ps,1)
    fprintf(fid,'%.7f %.7f\n',N_ps(k,1:2));
end
fprintf(fid,'%s\n','');

% fprintf(fid,'%s\n','vertices');  
% fprintf(fid,'%s\n',num2str(size(N_ps,1)));
% fprintf(fid,'%s\n',num2str(dim));
% for el = 1: nid
%     fprintf(fid, '%d %d\n', NRC(el,1)-1, NRC(el,3));
% end
% fprintf(fid,'%s\n','');

fclose(fid);



fname = ['dsF_' num2str(me) 'x' num2str(ne) '_' num2str(sets) '.txt'];
fid = fopen(fname,'w');
for k = 1: size(dis,1)
     fprintf(fid,'%.7f\n',dis(k,1))    ;
end 
fclose(fid);
