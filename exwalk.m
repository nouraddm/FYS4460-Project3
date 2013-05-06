%
% exwalk.m %
% Example of use of the walk routine
% Generate spanning cluster (l-r spanning)

close all;

% Here we want to find Psc as a function of p-pc for different L-values

pc = 0.59275;
p =[0.6 0.61 0.62 0.63 0.64 0.68 0.72 0.76 0.84 0.92];
L = [25 50 100 200 400];

sizeOfp = size(p,2);
sizeOfL = size(L,2);
nsample = 30;


M = zeros(1,sizeOfp);

for k = 1:sizeOfL
    lx = L(k);
    ly = lx;
    for j = 1:sizeOfp
        for sample = 1:nsample
            ncount = 0;
            perc = [];
            while (size(perc ,1)==0)
                ncount = ncount + 1; 
                if (ncount >1000)
                    return 
                end
                z = rand(lx,ly)<p(j); 
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
        
            M(j) = M(j) + sum(sum(zzz));
        end
        M(j) = M(j) / nsample;
    end
    figure(1);
    cstring='rgbcmyk'; % color string
    plot((p-pc), (M./(lx*lx)),cstring(mod(k,7)+1));
    xlabel('p-p_c')
    ylabel('P_S_C')
    hold on;
end
hold off;

figure(2);
subplot(2,2,1), imagesc(zz); 
subplot(2,2,2), imagesc(zadd); 
subplot(2,2,3), imagesc(zzz>0); 
subplot(2,2,4), imagesc(l+r>0);
