function [ktildez0,noptz0,cz0,vz0,kmax,kbar,ksd,kssM,kssno,nbar,nsd,cbar,csd,ibar,isd,ybar,ysd,v1,v2,v3,v4,iter,trashit,trashmax,v,ktildeposition,ktilde,n,c,z,kmin,sigma,sdz] = GRiskGHH2RCNfocfun(sdz,sigma,A,alpha,rho,beta,mu,delta,theta,psi,gp,gpz,T,Ti,itermax,k,v,UsekssTh,kssTh,countsdz,Reps)
% Solves the Markow growth problem with GHH-utility by iterating on a
% discrete grid.
% Resource constraint
% Uncertain productivity z, n=gpz different states
% Endog work effort solved from FOC

iter=0;
trashit=0;
trashmax=0;

% This part is for regular
%[P, z]=rouwen(rho, 0, sdz, gpz); % P(n,m) = Prob(z(t+1)=z_n|z(t)=z_m)

[P,z,~]=markovapprcj(rho,sdz*sqrt(1-rho^2),4,gpz);
z=z'; % Required with markovapprcj
P=P'; % Required with markovapprcj

% v=zeros(gp,gpz); % v(i,m) value function before iteration, initially all zeros
nextv=ones(gp,gpz); %value function after iteration, initially all ones (must differ from initial v)
U0=nan(gp,gp,gpz); %U(i,j,m) stores present period utility of going from k(i),z(m) to k(j)
n=nan(gp,gpz);

%toc
disp('Setting up U0...');
% Compute u0 first because it will not change as we iterate, so only do once
for i=1:gp
    for m=1:gpz

        nn=(((1-alpha)*(1+z(m))*A*(k(i))^alpha)/psi)^(1/(theta+alpha-1));
        n(i,m)=nn;

        for j=1:gp

            cc=(1+z(m))*A*(k(i))^alpha*nn^(1-alpha)+(1-delta)*k(i)+mu-k(j);
            u0=((cc-psi*(nn^theta)/theta)^(1-sigma)-1)/(1-sigma);

            if sigma == 0
                if cc >= 0 && nn >= 0
                    U0(i,j,m)=u0;
                else
                    U0(i,j,m)=-inf;
                end

            else
                if cc >= 0 && nn >= 0 && isreal(u0) && cc-psi*(nn^theta)/theta > 0
                    U0(i,j,m)=u0;
                else
                    U0(i,j,m)=-inf;
                end
            end

            % if isnan(nn)
            %     U0(i,j,m)=-inf;
            % elseif nn < 0  % non-pos work
            %     U0(i,j,m)=-inf;
            % else
            %     cc=(1+z(m))*A*(k(i))^alpha*nn^(1-alpha)+(1-delta)*k(i)+mu-k(j);
            %     if cc < 0 % Non-pos consumption
            %         U0(i,j,m)=-inf;
            %     else
            %         u0=((cc-psi*(nn^theta)/theta)^(1-sigma)-1)/(1-sigma);
            %         % u0=((cc-psi*((n(i,m))^theta)/theta)^(1-sigma))/(1-sigma);
            %         if isreal(u0)
            %             U0(i,j,m)=u0;
            %         else
            %             U0(i,j,m)=-inf;
            %         end
            %     end
            % end

        end
    end
end
clear nn cc u0;
disp('Done setting up U0!');
%toc
disp('Starting iterations...');

% keep iterating until the two values are close enough
% while max(abs(nextv-v),[],'all')/min(abs(nextv),[],'all')>tol
% while max(max(abs(nextv-v)))/min(min(abs(nextv)))>tol
while max(max(abs(nextv-v)))>0
    
    V=nan(gp,gp,gpz);
    
    if iter~=0 %update v, except in first iteration
        v=nextv;
    end
    
    iter=iter+1; %counter
    
    if iter>itermax
        break; %stop if max iterations reached
    end
    
    for m=1:gpz
        M=zeros(gp,gpz);
        M(:,m)=ones(gp,1);
        V(:,:,m)=U0(:,:,m)+beta*M*P'*v';
    end
    
    clear M;
    V(isnan(V))=-inf;
    V=permute(V,[1 3 2]);
    [nextv, ktildeposition]=max(V,[],3);
    clear V;
end

clear U0;
iter
%toc

while iter >= itermax
    disp('XXXXXXXXXXXX DID NOT CONVERGE XXXXXXXXXXXX');
    trashit=1;
    break;
end

while max(max(ktildeposition)) >=  gp
    disp('XXXXXXXXXXXX Kmax too low, increase kmaxgp XXXXXXXXXXXX');
    trashmax=1;
    break;
end


% Things below this line are only for plots and displaying results


ktilde=k(ktildeposition);
c=nan(gp,gpz);

for i=1:gp
    for m=1:gpz
        c(i,m)=(1+z(m))*A*(k(i))^alpha*(n(i,m))^(1-alpha)+(1-delta)*k(i)+mu-ktilde(i,m);
    end
end


% Finding the steady-state from computed rules, and feeding it into the
% simulation, or using theoretical ss
kssMpoints=find((1:gp)'==ktildeposition(:,round(gpz/2)));
kssM=sum(k(kssMpoints))/(size(kssMpoints,1)-1);
kssno=size(kssMpoints,1)-1;
if UsekssTh==1
    Dist1=abs(k-kssTh);
    idk1=round(find(Dist1==min(Dist1))); % theoretical kss
else
    Dist1=abs(k-kssM);
    idk1=find(Dist1==min(Dist1)); % computed kss, median of all ss points exluding 0
end
% if kssno>0
%     Kpositionsim(1,1)=idk1(1,1); %initial capital stock position
% else
%     Dist1=abs(k-kssTh);
%     idk1=round(find(Dist1==min(Dist1))); % theoretical kss
%     warning('Had to use theoretical ss, make k grid finer')
% end


kmax_s=nan(Reps,1);
kmin_s=nan(Reps,1);
kbar_s=nan(Reps,1);
ksd_s=nan(Reps,1);
nbar_s=nan(Reps,1);
nsd_s=nan(Reps,1);
cbar_s=nan(Reps,1);
csd_s=nan(Reps,1);
ibar_s=nan(Reps,1);
isd_s=nan(Reps,1);
ybar_s=nan(Reps,1);
ysd_s=nan(Reps,1);
v3_s=nan(Reps,1);


disp('Starting simulations...');

for r = 1:Reps
    
    sim=rand(Ti+T,1);
    Zsim=hitm_zn(z,P,sim); %simulated values of z, initial always z=E(z)
    Zpositionsim=hitm_sn(P,sim); %simulated positions in z-vector, consistent with sims from previous line
    clear sim;
    
    Kpositionsim=nan(Ti+T+1,1); %will hold simulated values for k
    Kpositionsim(1,1)=idk1(1,1); %initial capital stock position
    Nsim=nan(Ti+T,1);
    

    % Applying optimal decision rule
    for t=1:Ti+T
        Kpositionsim(t+1,1)=ktildeposition(Kpositionsim(t,1),Zpositionsim(t));
        Nsim(t,1)=n(Kpositionsim(t,1),Zpositionsim(t));
    end
    

    Kpositionsim=Kpositionsim(Ti+1:Ti+T);
    Ksim=k(Kpositionsim(1:T));
    clear Kpositionsim Zpositionsim;
    Zsim=Zsim(Ti+1:Ti+T);
    Nsim=Nsim(Ti+1:Ti+T);
        
    Ysim=A*(1+Zsim).*(Ksim.^alpha).*(Nsim.^(1-alpha));
    Csim=Ysim+(1-delta)*Ksim+mu-lagmatrix(Ksim,-1);
    Csim=Csim(1:T-1);
    Isim=Ysim(1:T-1)-Csim;
    
    kmax_s(r)=max(Ksim);
    kmin_s(r)=min(Ksim);
    kbar_s(r)=mean(Ksim);
    ksd_s(r)=std(Ksim);
    nbar_s(r)=mean(Nsim);
    nsd_s(r)=std(Nsim);
    cbar_s(r)=mean(Csim);
    csd_s(r)=std(Csim);
    ibar_s(r)=mean(Isim);
    isd_s(r)=std(Isim);
    ybar_s(r)=mean(Ysim);
    ysd_s(r)=std(Ysim);
    
    clear Ksim Isim Ysim Csim Nsim;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % HP-filtering and % deviaiton from trend before computing stats, as in
    % % data
    
    % % [Ksimtrend, Ksimcycle]=hpfilter(Ksim,1600);
    % % Ksimhat=Ksimcycle./Ksimtrend*100;
    % [Nsimtrend, Nsimcycle]=hpfilter(Nsim,1600);
    % Nsimhat=Nsimcycle./Nsimtrend*100;    
    % [Csimtrend, Csimcycle]=hpfilter(Csim,1600);
    % Csimhat=Csimcycle./Csimtrend*100;
    % % [Isimtrend, Isimcycle]=hpfilter(Isim,1600);
    % % Isimhat=Isimcycle./Isimtrend*100;
    % % [Ysimtrend, Ysimcycle]=hpfilter(Ysim,1600);
    % % Ysimhat=Ysimcycle./Ysimtrend*100;
    % clear Csimtrend Csimcycle Csim Nsim Nsimtrend Nsimcycle
    % 
    % nsdcompdata_s(r) = std(Nsimhat);
    % csdcompdata_s(r) = std(Csimhat);
    % acNsimhat = autocorr(Nsimhat,'NumLags',1);
    % acCsimhat = autocorr(Csimhat,'NumLags',1);
    % naccompdata_s(r) = acNsimhat(2);
    % caccompdata_s(r) = acCsimhat(2);
    % 
    % clear Csimhat Nsimhat;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % Dist2=abs(k-kbar);
    Dist2=abs(k-kbar_s(r));
    idk2=find(Dist2==min(Dist2)); % kbar mean simulated k
    % v3=v(idk2(1,1),round(gpz/2));  % v(kbar,E(z))
    v3_s(r)=v(idk2(1,1),round(gpz/2));  % v(kbar,E(z))
    % v4=v(idk2(1,1),:)*diag(P10000); % E(v(kbar,z))
    % v4_s(r)=v(idk2(1,1),:)*diag(P10000); % E(v(kbar,z))
    
end


disp('Done with simulations!');
%toc

kmax=max(kmax_s);
kmin=min(kmin_s);
kbar=mean(kbar_s);
ksd=mean(ksd_s);
nbar=mean(nbar_s);
nsd=mean(nsd_s);
cbar=mean(cbar_s);
csd=mean(csd_s);
ibar=mean(ibar_s);
isd=mean(isd_s);
ybar=mean(ybar_s);
ysd=mean(ysd_s);
v3=mean(v3_s);


v1=v(idk1(1,1),round(gpz/2)); % v(kss,E(z))
v2=nan; % E(v(kss,z))
v4=nan; % E(v(kbar,z))
%P10000=P^10000;
%v2=v(idk1(1,1),:)*diag(P10000); % E(v(kss,z))

ktildez0=ktilde(:,round(gpz/2));
noptz0=n(:,round(gpz/2));
cz0=c(:,round(gpz/2));
vz0=v(:,round(gpz/2));

%toc