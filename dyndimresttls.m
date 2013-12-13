%Same as dyndimrest but we solve it with the TLS.
% it is a guitar string in tension
% for which suddeny damage is added at the center.
% and a brittle damage model Y = Yc is used
clear all;
close all;
%faire converger le cas mineaire

ts_refine = 2;
end_t = 2;
h_refine = 46.41;

Ntim = 200*h_refine*ts_refine*end_t;
Nelt = 200*h_refine;
Nnod = Nelt+1;
E = 1; %(beton)
rho = 1; %(beton)
A = 1; % barre de 10cm sur 10cm
c = sqrt(E/rho);
L = 1;
h = 1/Nelt; %
cfl = 1./ts_refine;
dt = cfl * h/c;%
Yc = E/10;
ec = sqrt(2 * Yc / E);
sigc = E * ec;
lc = 0.02;
% two gauss point on the element
pg(1) = (1-sqrt(3)/3)/2;
pg(2) = (1+sqrt(3)/3)/2;
intorder = 6;

d = zeros(Ntim, Nelt);
s = zeros(Ntim, Nelt);
e = zeros(Ntim, Nelt);
energy = zeros(Ntim, Nelt);
Y = zeros(Ntim, Nelt);
strain_energy = zeros(Ntim, 1);
kinetic_energy = zeros(Ntim, 1);
dissip_energy = zeros(Ntim, 1);
ext_energy = zeros(Ntim, 1);
tot_energy = zeros(Ntim, 1);

bpos = zeros(Ntim, Nelt);
grad = zeros(Ntim, Nelt);

%s(1, 1:Nelt) = ones(1,Nelt);
%e(1, 1:Nelt) = ones(1,Nelt);

u = zeros(Ntim, Nnod);
ustat = zeros(Ntim, Nnod);
Ystat = zeros(Ntim, Nelt);
v = zeros(Ntim, Nnod);
a = zeros(Ntim, Nnod);
phi = zeros(Ntim, Nnod);
YmYc= zeros(Ntim, Nelt);
phidot = zeros(Ntim,1);
ddotbar = zeros(Ntim,1);
dissip = zeros(Ntim,1);
nbiter = zeros(Ntim,1);
nfrags = zeros(Ntim,1);


m = rho * h * A * ones(Nnod);
m(1) = m(1)/2;
m(Nnod) = m(Nnod)/2;

%initialization
for j=1:Nnod;
    x(j) = (j-1) * h;
end;
for j=1:Nelt;
    xe(j) = 0.5 * (x(j) + x(j+1));
end;
for i=1:Ntim;
    t(i) = (i-1) * dt;
end;

% Initially the bar is loaded and all elements are at Yc
% A tls is placed on the first element, obviously it satifies
%   Ybar = Yc
% This extra damage will create new stress that are less than in the next element
% and ill thus give an unloading wave.

%constant strain rate applied
csr = 0.25;
vbc = 1;
v(1,:) = vbc*csr*x;

for j=1:Nnod;
    u(1,j) = x(j) * ec * L * 0.999*(1-vbc);
    ustat(1,j) = x(j) * ec * L * 0.999*(1-vbc);
    phi(1,j) = (2*h-x(j))*(1-vbc)-vbc;%*0-1;
% %     if (x(j) > 0.3)
% %        phi(1,j) = x(j)-0.6+2*h;
% %     end
% % %         if (x(j) > 0.6)
% % %             phi(1,j) = 0.6-x(j)+2*h;
% % %         end
% %     if (x(j) > 0.7)
% %       phi(1,j) = 0.8-x(j)+2*h;
% %     end
%     if (x(j) > 0.1)
%        phi(1,j) = x(j)-0.2+2*h;
%     end
%     if (x(j) > 0.3)
%       phi(1,j) = 0.4-x(j)+2*h;
%     end

end

for j=1:Nelt;
    e(1,j) = (u(1,j+1)-u(1,j))/h;
    dloc = [0;0];
    if (phi(1,j) > 0  && phi(1,j+1) > 0)
        for k=1:2
            philoc = pg(k)*phi(1,j)+ (1-pg(k))*phi(1,j+1);
            dloc(k) = dval(philoc,lc);
            s(1,j) = s(1,j) + 0.5 * (1-dloc(k)) * E * e(1,j);
        end
    elseif  (phi(1,j) <= 0 && phi(1,j+1) <= 0)
        s(1,j) = E * e(1,j);
    elseif  (phi(1,j) > 0 && phi(1,j+1) <= 0)
        delta = abs(phi(1,j))/(abs(phi(1,j))+abs(phi(1,j+1)));
        sloc = 0;
        for k=1:2
            philoc = pg(k)*phi(1,j);
            dloc(k) = dval(philoc,lc);
            sloc = sloc + 0.5 * (1-dloc(k)) * E * e(1,j);
        end
        s(1,j) = delta *  sloc +  (1-delta) * E * e(1,j);
    elseif  (phi(1,j) <= 0 && phi(1,j+1) > 0)
        delta = abs(phi(1,j+1))/(abs(phi(1,j))+abs(phi(1,j+1)));
        sloc = 0;
        for k=1:2
            philoc = pg(k)*phi(1,j+1);
            dloc(k) = dval(philoc,lc);
            sloc = sloc + 0.5 * (1-dloc(k)) * E * e(1,j);
        end
        s(1,j) = delta *  sloc +  (1-delta) * E * e(1,j);
    end
    d(1,j) = 0.5*(dloc(1)+dloc(2));
end

%acceleration
a(1,1) = 0;
a(1,Nnod) = 0;
for j=2:Nnod-1;
    a(1,j) = A*(s(1,j) - s(1,j-1)) /m(j);
end

nbiter(1) = 0;
[phi(1,:),segments]=analyzeDamage(x,phi(1,:),h);
len = 0;
for l = 1:length(segments)
        if (size(segments{l},2)==0)
            continue;
        end
        len = len +1;
end
phidot = zeros(1,len);

for i=2:Ntim;
    
    % prediction
    v(i,:)= v(i-1,:) + 0.5*dt*a(i-1,:);
    u(i,:)= u(i-1,:) + dt*v(i-1,:) + 0.5*dt*dt*a(i-1,:);
    
    %def computation and Y update.
    
    for j=1:Nelt;
        e(i,j) = (u(i,j+1)-u(i,j))/(h);
        %b=0.5*E*e(i,j)*e(i,j)-Yc;
    end
    
    % moving the localization front
    % we compute a = integral (Yn+1 - Yc) d' in the current non-local zone
    % then we compute b = (Yn+1-Yc) d' on the front
    % the shift in level set if the ratio of the two.
    
    for j=1:Nnod;
        phi(i,j) = phi(i-1,j);
    end
    
    for l = 1:length(segments)
        if (size(segments{l},2)==0)
            continue;
        end

        sbegin = min(segments{l});
        send = max(segments{l});
        
        %skip if all negative
        if (sum(phi(i,sbegin:send)<0) == send-sbegin+1)
            continue;
        end
        
        err_crit = 1e15;
        nbiter(i) = 0;
        residu = 0;
        while (err_crit > 1.e-6)
            nbiter(i) = nbiter(i) + 1;
            residu_Y = 0; tangent_Y = 0;
            loop_residu = 0;
            loop_tangent = 0;
            for j=sbegin:min(send-1,Nelt)
                if (phi(i,j) > 0 && phi(i,j+1) > 0)
                    for k=1:2
                        philoc = pg(k)*phi(i,j) + (1-pg(k))*phi(i,j+1);
                        residu_Y = residu_Y + h * 0.5 * (0.5 * E * e(i,j) * e(i,j) - Yc) * dp(philoc,lc);
                        tangent_Y = tangent_Y + h * 0.5 * (0.5 * E * e(i,j) * e(i,j) - Yc) * dpp(philoc,lc);
                    end
                    loop_residu = loop_residu + 1;
                elseif  (phi(i,j) > 0 && phi(i,j+1) <= 0)
                    delta = h * abs(phi(i,j))/(abs(phi(i,j))+abs(phi(i,j+1))); %phi>0 portion
                    for k=1:2
                        philoc = pg(k)*phi(i,j);
                        residu_Y = residu_Y + delta *  0.5 * (0.5 * E * e(i,j) * e(i,j) - Yc) * dp(philoc,lc);
                        tangent_Y = tangent_Y + delta *  0.5 * (0.5 * E * e(i,j) * e(i,j) - Yc) * dpp(philoc,lc);
                    end
                    loop_residu = loop_residu + 1;
                    if (delta < h) tangent_Y = tangent_Y + (0.5 * E * e(i,j) * e(i,j) - Yc)* dp(0,lc);
                    else tangent_Y = tangent_Y + (0.5 * E * e(i,j+1) * e(i,j+1) - Yc)  * dp(0,lc);
                    end
                    loop_tangent = loop_tangent + 1;
                elseif  (phi(i,j) <= 0 && phi(i,j+1) > 0)
                    delta = h * abs(phi(i,j+1))/(abs(phi(i,j))+abs(phi(i,j+1))); %phi>0 portion
                    for k=1:2
                        philoc = pg(k)*phi(i,j+1);
                        residu_Y = residu_Y + delta *  0.5 * (0.5 * E * e(i,j) * e(i,j) - Yc) * dp(philoc,lc);
                        tangent_Y = tangent_Y + delta *  0.5 * (0.5 * E * e(i,j) * e(i,j) - Yc) * dpp(philoc,lc);
                    end
                    loop_residu = loop_residu + 1;
                    if (delta < h) tangent_Y = tangent_Y + (0.5 * E * e(i,j) * e(i,j) - Yc)* dp(0,lc);  %todo-doublecheck this
                    else tangent_Y = tangent_Y + (0.5 * E * e(i,j-1) * e(i,j-1) - Yc)  * dp(0,lc);              %todo-doublecheck this
                    end
                    loop_tangent = loop_tangent + 1;
                    
                end
            end
            
            % law Ybar = Yc
            if 1
                %i
                %loop_residu
                %residu_Y
                %loop_tangent
                %tangent_Y
                %i
                %dval(phi(i-1,1),lc)
                %err_crit = abs(residu_Y)/(Yc*dval(phi(i-1,1),lc));
                %if (abs(tangent_Y) <= 1.e-10) tangent_Y = 0.02; end
                %tangent_Y
                %dphi = -residu_Y/tangent_Y;
                %dphi
                %i
                %dval(phi(i,1),lc)
                flag = 1; %0 is exterior (damage centered on edge of domain); 1 is interior
                if (l == 1 && phi(i,1) > phi(i,2))
                    flag = 0;
                end
                if (l == length(segments) && phi(i,end) > phi(i,end-1))
                    flag = 0;
                end
                phimax = max(phi(i,sbegin:send));
                iphimax = find(phi(i,sbegin:send) == phimax);
                if (iphimax == Nnod)
                    phimaxY = 1/2*E*e(i,Nelt).^2; %1/2*s(i,Nelt)*e(i,Nelt); %;
                else
                    phimaxY = 1/2*E*e(i,iphimax).^2;% 1/2*s(i,iphimax)*e(i,iphimax); %
                end
                
                phimaxY = phimaxY(1);

                %phimax=phi(i,sbegin);
                YbarmYc = residu_Y/(dval(phimax,lc));
                oldresidu = residu;
                residu = YbarmYc/Yc;
                err_crit = abs(residu-oldresidu);
                %tangent = (tangent_Y + 0*(phimaxY-Yc)*dp(phimax,lc)/2)/(Yc*dval(phimax,lc)) - ((1+0/2)*dp(phimax,lc)/dval(phimax,lc)^2) * (YbarmYc/Yc);
                tangent = (tangent_Y + flag*(phimaxY-Yc)*dp(phimax,lc)/2)/(Yc*dval(phimax,lc)) - ((1+flag/2)*dp(phimax,lc)/dval(phimax,lc)^2) * (YbarmYc/Yc);
                if (abs(tangent) <= 1.e-10) err_crit = 0.; dphi = 0;
                else
                dphi = - residu/tangent;
                end
                %dphi = max(min(dphi,L/2),-L/2);
                
            end
            
            %delay effect law.
            if 0
                tauc=0.4;
                coef_a=1;
                l = max(phi(i,1),lc);
                dpbar = dval(phi(i,1),lc)/l;
                dphi_iter = phi(i,1) - phi(i-1,1);
                YbarmYc = residu_Y/(dval(phi(i,1),lc));
                if (YbarmYc < 0) YbarmYc = 0.; end
                residu_delay = - dpbar * dphi_iter * tauc/dt + (1-exp(-coef_a*YbarmYc/Yc));
                err_crit = abs(residu_delay);
                
                if (YbarmYc < 0) tangent = 0;
                else
                    tangent = tangent_Y/(Yc*dval(phi(i,1),lc)) - (dp(phi(i,1),lc)/dval(phi(i,1),lc)^2) * (YbarmYc/Yc);
                end
                
                tangent_delay = ((-dp(phi(i,1),lc)/l+dval(phi(i,1),lc)/l^2) * dphi_iter - dval(phi(i,1),lc)/l ) * tauc/dt + coef_a*exp(-coef_a*YbarmYc/Yc)*tangent;
                
                dphi = -residu_delay/tangent_delay;
            end
            

            if (nbiter(i)>50)
                dphi = 0;
                %assert(1==0);
            end
            
            if (isnan(dphi))
                dphi = 0;
                assert(1==0);                
            end
            
            for j=sbegin:send;
                if (i > 6 && intorder >= 6)
                    w = [60/147 360/147 -450/147 400/147 -225/147 72/147 -10/147];
                    phihist = [dphi; phi(i-1,j); phi(i-2,j); phi(i-3,j); phi(i-4,j); phi(i-5,j); phi(i-6,j)];                    
                    phi(i,j) = dphi*2/3 + 4/3*phi(i-1,j) - 1/3*phi(i-2,j);
                elseif (i > 5 && intorder >= 5)
                    w = [60/137 300/137 -300/137 200/137 -75/137 12/137];
                    phihist = [dphi; phi(i-1,j); phi(i-2,j); phi(i-3,j); phi(i-4,j); phi(i-5,j)];
                    phi(i,j) = w*phihist;  
                elseif (i > 4 && intorder >= 4)
                    w = [12/25 48/25 -36/25 16/25 -3/25];
                    phihist = [dphi; phi(i-1,j); phi(i-2,j); phi(i-3,j); phi(i-4,j)];
                    phi(i,j) = w*phihist;                    
                elseif (i > 3 && intorder >= 3)
                    w = [6/11 18/11 -9/11 2/11];
                    phihist = [dphi; phi(i-1,j); phi(i-2,j); phi(i-3,j)];
                    phi(i,j) = w*phihist;
                elseif (i > 2 && intorder >= 2)
                    w = [2/3 4/3 -1/3];
                    phihist = [dphi; phi(i-1,j); phi(i-2,j)];
                    phi(i,j) = w*phihist;
                else
                    phi(i,j) = phi(i-1,j) + dphi;
                end
                phi(i,j) = max(phi(i,j),phi(i-1,j)); %constraint: dphi >= 0
                %enforcing limit of level-set motion
                %phi(i,j) = min(phi(i-1,j)+h,phi(i,j));
            end
        end %while
        
        %err_crit = 0.
        
    end %for segments
    
    %check for nucleation
    phi(i,:) = checkFailureCriteria(t(i),x,phi(i,:),Yc*(randn(Nelt,1)*.00+1),'elem',0.5*E*e(i,:).^2,0,0,2*h);
    
    %enforce phi constraints
    [phi(i,:),segments]=analyzeDamage(x,phi(i,:),h);
    
    for l=1:length(segments)
        if (size(segments{l},2)==0)
            continue;
        end
        median(segments{l});
        smid = floor(median(segments{l}));
        phidot(i,l) = (phi(i,smid) - phi(i-1,smid))/dt;
        if phidot(i,1)*dt > h*1.01
            sprintf('level-set front advancing more than one element per time-step: t=%f, segment %u , dphi/h = %f',t(i),l,phidot(i,1)*dt/h)
        end
    end
        
    %updating the stress
    for j=1:Nelt;
        if (phi(i,j) > 0 && phi(i,j+1) > 0)
            for k=1:2
                philoc = pg(k)*phi(i,j)+ (1-pg(k))*phi(i,j+1);
                dlocg(k) = dval(philoc,lc);
                s(i,j) = s(i,j) + 0.5 * (1-dlocg(k)) * E * e(i,j);
            end
        elseif  (phi(i,j) <= 0 && phi(i,j+1) <= 0)
            s(i,j) = E * e(i,j);
            dlocg(1) = 0; dlocg(2) = 0;
        elseif  (phi(i,j) > 0 && phi(i,j+1) <= 0)
            delta = abs(phi(i,j))/(abs(phi(i,j))+abs(phi(i,j+1)));
            sloc = 0;
            for k=1:2
                philoc = pg(k)*phi(i,j);
                dlocg(k) = dval(philoc,lc);
                sloc = sloc + 0.5 * (1-dlocg(k)) * E * e(i,j);
            end
            s(i,j) = delta * sloc +  (1-delta) * E * e(i,j);
        elseif  (phi(i,j) <= 0 && phi(i,j+1) > 0)
            delta = abs(phi(i,j+1))/(abs(phi(i,j))+abs(phi(i,j+1)));
            sloc = 0;
            for k=1:2
                philoc = pg(k)*phi(i,j+1);
                dlocg(k) = dval(philoc,lc);
                sloc = sloc + 0.5 * (1-dlocg(k)) * E * e(i,j);
            end
            s(i,j) = delta * sloc +  (1-delta) * E * e(i,j);
        end
        d(i,j) = 0.5*(dlocg(1)+dlocg(2));
    end
    
    %acceleration
    %if (phi(i,1)/lc < 1)
    if (phi(i,1) <= lc) a(i,1) = 0;
    else a(i,1) =  A*s(i,1) /m(1);
    end
    a(i,Nnod) = 0; %note: this is a constraint which was hidden. if free, it would be -s(i,Nnod)/m(j);
    for j=2:Nnod-1;
        a(i,j) = A*(s(i,j) - s(i,j-1)) /m(j);
    end
    
    
    %correction
    v(i,:)= v(i,:) + 0.5*dt*a(i,:);
    
    %record number of fragments and quantities per fragment
    [nfrags(i),fragment_list] = findFragments(x,phi(i,:),d(i,:));
    
    
end

for i = 1:Ntim;
    for j=1:Nelt;
        Y(i,j) = 0.5*E*e(i,j)^2;
        YmYc(i,j) = Y(i,j)/Yc - 1;
        if (i > 1)
            dissip(i) = dissip(i) + h * Y(i,j) * (d(i,j)-d(i-1,j));
            kinetic_energy(i) = kinetic_energy(i) + 0.5 * h * rho * 0.5 * (((u(i,j) - u(i-1,j))^2)+((u(i,j+1) - u(i-1,j+1))^2))/dt^2;
            %ustat(i,j+1) = ustat(i,j) + h*s(1,Nelt)/(E*(1-d(i,j)));
        end;
        energy(i,j)=h*Y(i,j)*(1-d(i,j));
    end
    %ustat(i,:) *= u(1,Nnod)/ustat(i,Nnod);
    for j=1:Nelt;
        Ystat(i,j) = 0.5*E*((ustat(i,j+1)-ustat(i,j))/h)^2;
    end;
    if (i > 1) ext_energy(i) = (a(i,Nnod)*m(j)+s(i,Nnod-1))*v(i,Nnod)*dt + ext_energy(i-1); end
    if (i > 1) dissip_energy(i) = dissip_energy(i-1) + dissip(i); end
    strain_energy(i) = sum (energy(i,:));
    tot_energy(i) = strain_energy(i) + kinetic_energy(i) + dissip_energy(i) - ext_energy(i);
end


minfrag = L*2;
for i=1:length(round(nfrags(end)/2))
    
    fragLen = (fragment_list{i}(2)-fragment_list{i}(1))*h;
    if ((mod(nfrags(end),2) == 1) && (i == 1))
        fragLen = fragLen * 2;
    end

    if (fragLen < minfrag)
        minfrag = fragLen;
    end
end
sprintf('Final number of fragments: %i \nMinimum fragment length: %f \nFinal dissipated energy: %f \n',nfrags(end),minfrag,dissip_energy(end))


%%
fig=figure;
%plot(x+u(end,:),0*x,'x-');
%plot(x+u(end,:),[d(end,1) d(end,:)])
%plot(x+u(end,:),[s(end,1) s(end,:)])
%plot(x+u(end,:),phi(end,:),'x-')
%axis([min(min(ones(size(u,1),1)*x+u)),max(max(ones(size(u,1),1)*x+u)),min(min(phi))*1.1,max(max(d))*1.1])
set(gca,'NextPlot','replaceChildren');


% Preallocate the struct array for the struct returned by getframe
F(size(u,1)) = struct('cdata',[],'colormap',[]);
% Record the movie
for j = 1:size(u,1)
    clf(fig)
    col = [s(j,1) s(j,:)];
    Y = (e(j,:).*s(j,:)*0.5)'; %this equals Y above *(1-d)
    X = x+u(j,:);
    Z = x*0;
    %plot(x,u(j,:)+x,'x-')
    %plot(x+u(j,:),v(j,:),'x-')
    %surface(X,Y,Z,col);
    %plot(x+u(j,:),x*0,'x-');
    hold on
    %plot(x+0*u(j,:),[s(end,1) s(j,:)],'x-')
    plot(x+0*u(j,:),[d(j,1) d(j,:)],'x-')
    plot(x+0*u(j,:),phi(j,:),'rx-')
    plot(x(1:end-1)*0.5+x(2:end)*0.5,Y/Yc,'xg-')
    %plot(x+0*u(j,:),[s(j,1) s(j,:)],'gx-')
        xlabel('Position, x')
        ylabel('d,\phi,Y/Yc')
    legend('d','\phi','Y/Yc',0)
    title(sprintf('t = %f',t(j)));
    %plot(x,phi(j,:),'x-')
    F(j) = getframe(fig);    
end
writerObj = VideoWriter('video.avi'); writerObj.Quality = 100;
open(writerObj);writeVideo(writerObj,F);close(writerObj);
% use 1st frame to get dimensions
%[h, w, p] = size(F(1).cdata);
%hf = figure; 
% resize figure based on frame's w x h, and place at (150, 150)
%set(hf,'Position', [150 150 w h]);
%axis off
% Place frames at bottom left
%movie(hf,F,4,30,[0 0 0 0]);

% figure
% surf(t,x,u','EdgeAlpha',0.25)
% xlabel('time')
% ylabel('x')
% zlabel('displacement')
% surf(t,x,[s(1:s) s]','EdgeAlpha',0.5)
% zlabel('stress')
%%

figure
plot(t,phidot)
xlabel('Time, t')
ylabel('d\phi/dt')
title('Level-Set Movement by Segment')

figure
plot(x,phi)
xlabel('Position, x')
ylabel('\phi')
title('Level-Set Movement')

figure
plot(t,[nfrags,dissip_energy/dissip_energy(end)*nfrags(end)])
xlabel('Time, t')
ylabel('# of Fragments')
title('Fragment Count')
legend('# frags','scaled dissip energy',0)


figure
plot(t,[strain_energy, kinetic_energy, ext_energy, dissip_energy,  tot_energy])
xlabel('Time, t')
ylabel('Energy')
title('Energy Plot')
legend('str_E', 'kin_E', 'ext_E', 'dissip_E','total_E',0)

