%
% exwalk.m %
% Example of use of the walk routine
% Generate spanning cluster (l-r spanning)

close all;

% Here we want to find the mass Msc and Dsc

pc = 0.59275;

L = [16 32 64 128 256 512];

sizeOfL = size(L,2);
nsample = 30;

M = zeros(1,sizeOfL);

for k = 1:sizeOfL
    lx = L(k);
    ly = lx;

    for j = 1:nsample
        ncount = 0;
        perc = [];
        while (size(perc ,1)==0)
            ncount = ncount + 1; 
            if (ncount >1000)
                return 
            end
            z=rand(lx,ly)<pc; 
            [lw,num]=bwlabel(z,4);
            perc_x = intersect(lw(1,:),lw(lx,:)); 
            perc = find(perc_x >0)
        end

        s = regionprops(lw,'Area');
        clusterareas = cat(1,s.Area);
        maxarea = max(clusterareas);
        i = find(clusterareas==maxarea);
        zz = lw == i;
    
        % zz now contains the spanning cluster 
        %imagesc(zz); % Display spanning cluster

        % Run walk on this cluster
        [l,r] = walk(zz);
        zzz = l.*r; % Find points where both l and r are non-zero 
        zadd = zz + zzz;
        
        M(k) = M(k) + sum(sum(zzz));
    end
    M(k) = M(k) / nsample;
end

figure(1);
plot(log10(L), log10(M),'-o');
xlabel('log10(L)')
ylabel('log10(M_S_C)')