%
% exflow.m
%

% Exercise (p): We find the conductivity

clear all; 
clf;
% First , find the backbone
% Generate spanning cluster (l-r spanning) 
limit = 0.00000000001;

pc = 0.59275;
pv = pc;%(pc:0.05:1.0);
L = [10 20 30 40 50 60 70 80 90 100];
%lx = 20;
%ly = 20;

nx = size(pv,2);
nsample = 10;
nL = size(L,2);

conductivity = zeros(nL,1);%zeros(nx,1);
M = zeros(nL,1);
Msingle = zeros(nL,1);
Mdangle = zeros(nL,1);
Mback = zeros(nL,1);

for lv = 1:nL
%for vp = 1:nx
    for sample = 1:nsample
    p = pv;
    lx = L(lv);
    ly = lx;
    
    ncount = 0;
    perc = [];

    while (size(perc ,1)==0)
        ncount = ncount + 1; 
        if (ncount >1000)
            return 
        end
        z = rand(lx,ly) < p; 
        [lw,num]=bwlabel(z,4);
        perc_x = intersect(lw(1,:),lw(lx,:)); 
        perc = find(perc_x >0);
    end

    s = regionprops(lw,'Area'); 
    clusterareas = cat(1,s.Area);
    maxarea = max(clusterareas);
    i = find(clusterareas==maxarea);
    zz = lw == i;

    % zz now contains the spanning cluster 

    [l,r] = walk(zz);
    single_connected_bonds = l.*r;
    
    % Transpose
    zzz = zz';

    % Generate bond lattice from this
    g = sitetobond(zzz);

    % Generate conductivity matrix
    [p c_eff] = FIND_COND(g,lx,ly);
    %conductivity(vp) = conductivity(vp) + c_eff; % for different p-values
    conductivity(lv) = conductivity(lv) + c_eff;% for different L-values
    
   
    % Transform this onto a nx x ny lattice 
    x = coltomat(full(p),lx,ly);
    P = x.*zzz;
    g1 = g(:,1);
    g2 = g(:,2);
    z1 = coltomat(g1,lx,ly);
    z2 = coltomat(g2,lx,ly);

    % Plotting
    subplot(2,2,1), imagesc(zzz); 
    title('Spanning cluster')
    axis equal
    subplot(2,2,2), imagesc(P); 
    title('Pressure');
    axis equal

    f2 = zeros(lx,ly);
    for iy = 1:ly-1
        f2(:,iy) = (P(:,iy) - P(:,iy+1)).*z2(:,iy); 
    end
    
    f1 = zeros(lx,ly);
    for ix = 1:lx-1
        f1(ix,:) = (P(ix,:) - P(ix+1,:)).*z1(ix,:);
    end
    
    % Find the sum of absolute fluxes into each site
    fn = zeros(lx,ly);
    fn = fn + abs(f1);
    fn = fn + abs(f2);
    fn(:,2:ly) = fn(:,2:ly) + abs(f2(:,1:ly-1));
    fn(:,1) = fn(:,1) + abs((P(:,1) - 1.0).*(zzz(:,1))); 
    fn(2:lx,:) = fn(2:lx,:) + abs(f1(1:lx-1,:)); 

    %plots
    subplot(2,2,3), imagesc(fn);
    title('Flux');
    axis equal
    zfn = fn>limit;
    zbb = (zzz + 2*zfn);
    zbb = zbb/max(max(zbb)); 
    subplot(2,2,4), imagesc(zbb); 
    title('BB and DE');
    axis equal
    hold off;
    
    
    figure(10)
    for i=1:length(single_connected_bonds)
        for j=1:length(single_connected_bonds)
            if (single_connected_bonds(i,j) > 1)
                single_connected_bonds(i,j) = 1;
            end
        end
    end
    
    subplot (2 ,2 ,1) , imagesc ( single_connected_bonds );
    title ('Single connected bonds');

    zJJ = zeros(length(zbb), length(zbb));

    for i=1:length(zbb)
        for j=1:length(zbb)
            if (zbb(i,j) == 1)
                zJJ(i,j) = 1;
            end
        end
    end

    zJK = zbb - zJJ;

    subplot (2 ,2 ,2) , imagesc ( zJJ );
    title ('Backbone');
    subplot (2 ,2 ,3) , imagesc ( zJK );
    title ('Dangling ends');
    subplot (2 ,2 ,4) , imagesc (zJJ+zJK+2*single_connected_bonds);
    title ('BB, DE og SCB in the same plot');
    
    end
    %conductivity(vp) = conductivity(vp)/nsample; %p-values
    conductivity(lv) = conductivity(lv)/nsample; %L-values   
end

% plot for conductivity as a function of p-pc 
% figure(20);
% for i=2:length(conductivity)
%     XX(i-1) = pv(i)-pc;
%     YY(i-1) = conductivity(i);
% end
% 
% plot(log10(XX), log10(YY), '-o')

% plot for conductivity as a function of L
figure(30);
for i=2:length(conductivity)
    XX(i-1) = L(i);
    YY(i-1) = conductivity(i);
end

plot(log10(XX), log10(YY), '-o')

