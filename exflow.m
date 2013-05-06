%
% exflow.m
%

% Exercise (o): We find the dimensionality of singly connected bonds,
% the backbone and the dangling ends.


clear all; 
clf;
% First , find the backbone
% Generate spanning cluster (l-r spanning) 
limit = 0.00000000001;

pc = 0.59275;
L = [10 20 30 40 50 60 70 80 90 100];

nsample = 10;
nL = size(L,2);

conductivity = zeros(nL,1);
M = zeros(nL,1);
Msingle = zeros(nL,1);
Mdangle = zeros(nL,1);
Mback = zeros(nL,1);



for k = 1:nL
    lx = L(k);
    ly = lx;
    for sample = 1:nsample
        ncount = 0;
        perc = [];
        
        while (size(perc ,1)==0)
            ncount = ncount + 1; 
            if (ncount >1000)
                return 
            end
            z = rand(lx,ly) < pc; 
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
        conductivity(k) = conductivity(k) + c_eff;
    
        M(k) = M(k) + sum(sum(zz));
        Msingle(k) = Msingle(k) + sum(sum(single_connected_bonds));
    
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


        for i=1:length(zbb)
            for j=1:length(zbb)
                if(zbb(i,j) < 1)
                    if(zbb(i,j) > 0)
                        Mdangle(k) = Mdangle(k) + 1;
                    end
                end
            end
        end

        for i=1:length(zbb)
            for j=1:length(zbb)
                if(zbb(i,j) == 1)
                    Mback(k) = Mback(k) + 1;
                end
            end
        end
        
    end

    conductivity(k) = conductivity(k)/nsample; 
    M(k) = M(k)/nsample;
    Msingle(k) = Msingle(k)/nsample;
    Mdangle(k) = Mdangle(k)/nsample;
    Mback(k) = Mback(k)/nsample;
end
hold off;
figure(2);
plot(log10(L), log10(M), '-o')
xlabel('log10(L)');
ylabel('log10(M_t_o_t_a_l)');

figure(3);
plot(log10(L), log10(Msingle), '-o')
xlabel('log10(L)');
ylabel('log10(M_S_C)');

figure(4);
plot(log10(L), log10(Mdangle), '-o')
xlabel('log10(L)');
ylabel('log10(M_D_E)');

figure(5);
plot(log10(L), log10(Mback), '-o')
xlabel('log10(L)');
ylabel('log10(M_B_B)');
