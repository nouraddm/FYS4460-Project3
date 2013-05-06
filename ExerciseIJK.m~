%-------------------------------------------------------------------------
%                          
%                           Project 3
%                           
%                       Exercise (i), (j) & (k)
%
%                         Finit Size Scaling
%-------------------------------------------------------------------------
close all;

pc = 0.59275;
numberOfSamples = 100;

L = [25 50222222 100 200 400 800];
p = (0.55:0.001:.65);

lengthL = size(L,2);
lengthP = size(p,2);

x1 = 0.8;
x2 = 0.3;

%Pi = zeros(lengthP,1);
PiMatrix = zeros(lengthL,lengthP);
pPi8 = zeros(lengthL, 1);
pPi3 = zeros(lengthL, 1);

for i = 1:lengthL
    Lx = L(i);
    Ly = Lx;
    A = Lx*Ly;
    
    Pi = zeros(lengthP,1);
    
    for j = 1 : lengthP
        for k = 1 : numberOfSamples
            z = rand(Lx,Ly);
            zz = z < p(j);
            [lw,num] = bwlabel(zz,4);
            perc_y = intersect(lw(:,1),lw(:,Ly)); 
            perc_x = intersect(lw(1,:),lw(Lx,:)); 
            perc_xy = union(perc_x,perc_y);
            perc = find(perc_xy >0);
            s = regionprops(lw,'Area'); 
            clusterareas = cat(1,s.Area);
            if (length(perc)>0)
                Pi(j) = Pi(j) + 1;
            end
        end
        Pi(j) = Pi(j)/numberOfSamples;
        PiMatrix(i,j) = Pi(j);
    end
    % Exersice (k)
    % Data collapse plot for Pi:
%     figure(1);
%     cstring='rgbcmyk'; % color string
%     plot(((p-pc)*Lx^(3/4)), Pi,cstring(mod(i,7)+1))
%     xlabel('(p-p_c)L^1^/^\ni')
%     ylabel('\Pi(p,L)')
%     hold on;
    
end
%hold off;

% The method to find the closest value to p:

for i = 1:lengthL
    tempConst1 = 100;
    for j = 1:lengthP
        tempConst2 = x1 - PiMatrix(i, j);
        tempConst2 = abs(tempConst2);
        if (tempConst2 < tempConst1)
            tempConst1 = tempConst2;
            pPi8(i) = p(j);
        end
    end
end

for i = 1:lengthL
    tempConst1 = 100;
    for j=1:lengthP
        tempConst2 = x2 - PiMatrix(i, j);
        tempConst2 = abs(tempConst2);
        if (tempConst2 < tempConst1)
            tempConst1 = tempConst2;
            pPi3(i) = p(j);
        end
    end
end

% Exercise i

figure(2)
plot(L, pPi8,'-o')
xlabel('L');
ylabel('p_\Pi_=_0_._8');

figure(3)
plot(L, pPi3,'-o')
xlabel('L');
ylabel('p_\Pi_=_0_._3');


% Exercise j

% Here we find that nu er 1/0.78359 = 1.2762, which is very close
% to the exact value nu = 4/3 = 1.333....

pPix = pPi8 - pPi3;

figure(4)
plot(log10(L), log10(pPix), '-o')
xlabel('log10(L)');
ylabel('log10(p_\Pi_=_0_._8 - p_\Pi_=_0_._3)');



% Exercise k

figure(5);

plot(L.^(-3/4), pPi8,'-or')
hold on
plot(L.^(-3/4), pPi3,'-o')
xyz = [0 0.1];
xyzk = [0.59275 0.59275];
hold on;
plot(xyz,xyzk, 'g');
xlabel('L^-^1^/^\nu');
ylabel('p_\Pi_=_x');
hold off;
