%-------------------------------------------------------------------------
%                          
%                           Project 3
%
%                        Exercise (a) & (b)
%
%-------------------------------------------------------------------------

pc = 0.59275;

p = (0.6:0.001:1.0);
sizeOfp = size(p,2);

Lx = 100;
Ly = Lx;
A = Lx * Ly;

numberOfSamples = 20;

Pi = zeros(1,sizeOfp);
P = zeros(1,sizeOfp);
M = zeros(1,sizeOfp);


for i = 1:numberOfSamples
    r = rand(Lx,Ly);
    for j = 1:sizeOfp
        z = r < p(j);
        [lw, num] = bwlabel(z,4);
        
        perc_x = intersect(lw(:,1),lw(:,Ly));
        perc_y = intersect(lw(1,:),lw(Lx,:));
        perc_xy = intersect(perc_y, perc_x);
        nonZeros = find(perc_xy>0);
        perc = perc_xy(nonZeros);
        
        if(length(nonZeros)>0)
            s = regionprops(lw,'Area');
            area = cat(1,s.Area);
            ar = sum(area(perc));
            P(1,j) = P(1,j) + ar/A;
        end
         M(1,j) = A * P(1,j);
    end
end
Pi(1,:) = Pi(1,:)/numberOfSamples;
P(1,:)  = P(1,:)/numberOfSamples;
M(1,:)  = M(1,:)/numberOfSamples;

figure(1);
plot(p,P(1,:));

figure(2);
plot(log10(p-pc),log10(P(1,:)))

