clear;
clc;
tic;

% tv = datetime('now', 'Format','yyyy_MM_dd_HH_mm');
% logname = sprintf('GRiskCESRCNfocPC_%s.log',tv);
% diary(logname); % Saves log

display(datetime);

% %cpus = feature('numcores') % cpus = maxNumCompThreads
% cpus = 10 % 6 12 
% t.TimeZone = 'America/New_York';
% patchJobStorageLocation; % Needed for 2020 Matlab
% delete(gcp('nocreate'));
% mypool = parpool(cpus);
% disp('Pool loaded!');

alpha=1/3 %1/3;
rho=.9 %.9
beta=.99; %.99 .989-.995
mu=0; %safe income
delta=.025; %.025 also yields safe income .01-.04 quarterly, .025 popular
theta=1.5 % Frisch elast = 1/(theta-1), 1.25-2.25 since otherwise makes MU leisure non-decreasing in leisure
%sdzsteps=4
%sdzstep=.04/sdzsteps; %.005
sigmavect=[0]; %[20; 15; 10; 5; 4; 3; 2; 1.01; .5; 0]; CHANGE FIG LABELS IF CHANGE, [4, 3, 2, 1.01, .5, 0] % sigma can't equal 1
labels={[char(949) '=4']; [char(949) '=3']; [char(949) '=2']; [char(949) '=1']; [char(949) '=.5']; [char(949) '=0']};
colors={[0.6859 0.4035 0.2412]; [0.2941 0.5447 0.7494]; [0.3718 0.7176 0.3612]; [0.8650 0.8110 0.4330]; [1.0000 0.5482 0.1000]; [0.9047 0.1918 0.1988]}; % colors generated by linspecer.m

gp=8000 % 4000 even number of gridpoints k-vector
gpz=51; % 51 19,21 odd number>=5 Stops running properly if too large, and optimal k-rule goes out of bounds with so high shock values

T=1e8 % 1e8 (if HP filtering) or 2e8, simulation length. 4e8 is too much for 16GB RAM
Reps=1e3 % 1e3 Repetitions in simulations. Faster with larger T, less reps.
Ti=1e6; % 1e6


kssTh=100; % still needed even if not used in simulations
nssTh=1/3; %1/3 40/(24*7);



kmaxgp=kssTh*13 % 5 13 (sdz=.2);  (sdz=.25)
%k=linspace(0,kmaxgp,gp)'; %makes k be a discrete grid
k1=linspace(1e-8,kssTh,.32*gp+1)';
k2=linspace(kssTh+k1(2),kssTh*3,.64*gp)'; % Adjust so kmax from sim is within k2 (fine grid)
k3=linspace(max(k2)+k1(2),kmaxgp,gp-size(k1,1)-size(k2,1))'; % adjust kmaxgp until no complaints about this
k=[k1;k2;k3];
clear k1 k2 k3;
gporig=gp;
gp=size(k,1);

% Changes in utility or production function requires modifying equations below:
UsekssTh=1; % in simulations, = 1 to use theoretical in simulations, 0 to use nummerical approx from solution
A=(1/beta-1+delta)/(alpha*kssTh^(alpha-1)*nssTh^(1-alpha));
psi=(1-alpha)*A*kssTh^alpha*nssTh^(1-alpha-theta);
cssTh=A*kssTh^alpha*nssTh^(1-alpha) - delta*kssTh;

sigmavect=sigmavect/(cssTh/(cssTh-psi*nssTh^theta/theta))

GBRAMreq=2*(2*gp^2*gpz+4*gp*gpz+gpz^2)*8/1024^3 % 8 bit with double precision

itermax=10000; %5000 maximum iterations allowed
%tol=1e-8; % 1e-8 tolerance

%GBRAMreq=((gp^2*gpz*(sdzsteps+1))+2*(gp^2*gpz))*8/1024^3 % 8 bit with double precision

% if UsekssTh==1
%     Ti=0;
% else
%     Ti=1e4; % Drops first Ti simulated. Starts at SS, but multiple SS
% end

sdzvect = (.015); %(.0325);
lsdzvect = length(sdzvect);
lsigmavect = length(sigmavect);

iterations = [lsdzvect, lsigmavect];
numberiterations = prod(iterations);

S_KTILDEZ0n=nan(gp,numberiterations);
S_NOPTZ0n=nan(gp,numberiterations);
S_CZ0n=nan(gp,numberiterations);
S_VZ0n=nan(gp,numberiterations);
S_KMAXn=nan(1,numberiterations);
S_KBARn=nan(1,numberiterations);
S_KSDn=nan(1,numberiterations);
S_KSSMn=nan(1,numberiterations);
S_KSSNOn=nan(1,numberiterations);
S_NBARn=nan(1,numberiterations);
S_NSDn=nan(1,numberiterations);
S_CBARn=nan(1,numberiterations);
S_CSDn=nan(1,numberiterations);
S_IBARn=nan(1,numberiterations);
S_ISDn=nan(1,numberiterations);
S_YBARn=nan(1,numberiterations);
S_YSDn=nan(1,numberiterations);
S_V1n=nan(1,numberiterations);
S_V2n=nan(1,numberiterations);
S_V3n=nan(1,numberiterations);
S_V4n=nan(1,numberiterations);
S_ITERn=nan(1,numberiterations);
S_TRASHITn=nan(1,numberiterations);
S_TRASHMAXn=nan(1,numberiterations);
Vn=nan(gp,gpz,numberiterations);
Ktildepositionn=nan(gp,gpz,numberiterations);
Ktilden=nan(gp,gpz,numberiterations);
Nn=nan(gp,gpz,numberiterations);
Cn=nan(gp,gpz,numberiterations);
Zn=nan(gpz,numberiterations);
S_KMINn=nan(1,numberiterations);
S_SIGMAn=nan(1,numberiterations);
S_SDZn=nan(1,numberiterations);

v=zeros(gp,gpz); % v(i,m) value function before iteration. Below is recycled from previous 

%par
for ix = 1:numberiterations
    [countsdz, countsigma] = ind2sub(iterations,ix);
    sigma=sigmavect(countsigma);
    sdz=sdzvect(countsdz);
    
    %v=zeros(gp,gpz); % v(i,m) value function before iteration. Below is recycled from previous
    
    [S_KTILDEZ0n(:,ix),S_NOPTZ0n(:,ix),S_CZ0n(:,ix),S_VZ0n(:,ix),S_KMAXn(:,ix),S_KBARn(:,ix),S_KSDn(:,ix),S_KSSMn(:,ix),S_KSSNOn(:,ix),S_NBARn(:,ix),S_NSDn(:,ix),S_CBARn(:,ix),S_CSDn(:,ix),S_IBARn(:,ix),S_ISDn(:,ix),S_YBARn(:,ix),S_YSDn(:,ix),S_V1n(:,ix),S_V2n(:,ix),S_V3n(:,ix),S_V4n(:,ix),S_ITERn(:,ix),S_TRASHITn(:,ix),S_TRASHMAXn(:,ix),Vn(:,:,ix),Ktildepositionn(:,:,ix),Ktilden(:,:,ix),Nn(:,:,ix),Cn(:,:,ix),Zn(:,ix),S_KMINn(:,ix),S_SIGMAn(:,ix),S_SDZn(:,ix)] = GRiskGHH2RCNfocfun(sdz,sigma,A,alpha,rho,beta,mu,delta,theta,psi,gp,gpz,T,Ti,itermax,k,v,UsekssTh,kssTh,countsdz,Reps);
end

clear v;
Runtime=toc
display(datetime);

while max(S_TRASHITn,[],'all') > 0
    error('XXXXXX NOT SAVED - ALL DID NOT CONVERGE XXXXXX');
end


% Reconverting from linear index to matrix
S_KTILDEZ0=nan(gp,lsdzvect,lsigmavect);
S_NOPTZ0=nan(gp,lsdzvect,lsigmavect);
S_CZ0=nan(gp,lsdzvect,lsigmavect);
S_VZ0=nan(gp,lsdzvect,lsigmavect);
S_KMAX=nan(1,lsdzvect,lsigmavect);
S_KBAR=nan(1,lsdzvect,lsigmavect);
S_KSD=nan(1,lsdzvect,lsigmavect);
S_KSSM=nan(1,lsdzvect,lsigmavect);
S_KSSNO=nan(1,lsdzvect,lsigmavect);
S_NBAR=nan(1,lsdzvect,lsigmavect);
S_NSD=nan(1,lsdzvect,lsigmavect);
S_CBAR=nan(1,lsdzvect,lsigmavect);
S_CSD=nan(1,lsdzvect,lsigmavect);
S_IBAR=nan(1,lsdzvect,lsigmavect);
S_ISD=nan(1,lsdzvect,lsigmavect);
S_YBAR=nan(1,lsdzvect,lsigmavect);
S_YSD=nan(1,lsdzvect,lsigmavect);
S_V1=nan(1,lsdzvect,lsigmavect);
S_V2=nan(1,lsdzvect,lsigmavect);
S_V3=nan(1,lsdzvect,lsigmavect);
S_V4=nan(1,lsdzvect,lsigmavect);
S_ITER=nan(1,lsdzvect,lsigmavect);
S_TRASHIT=nan(1,lsdzvect,lsigmavect);
S_TRASHMAX=nan(1,lsdzvect,lsigmavect);
V=nan(gp,gpz,lsdzvect,lsigmavect);
Ktildeposition=nan(gp,gpz,lsdzvect,lsigmavect);
Ktilde=nan(gp,gpz,lsdzvect,lsigmavect);
N=nan(gp,gpz,lsdzvect,lsigmavect);
C=nan(gp,gpz,lsdzvect,lsigmavect);
Z=nan(gpz,lsdzvect,lsigmavect);
S_KMIN=nan(1,lsdzvect,lsigmavect);
S_SIGMA=nan(1,lsdzvect,lsigmavect);
S_SDZ=nan(1,lsdzvect,lsigmavect);


for csigma=1:lsigmavect
    for csdz=1:lsdzvect
        S_KTILDEZ0(:,csdz,csigma)=S_KTILDEZ0n(:,sub2ind(iterations,csdz,csigma));
        S_NOPTZ0(:,csdz,csigma)=S_NOPTZ0n(:,sub2ind(iterations,csdz,csigma));
        S_CZ0(:,csdz,csigma)=S_CZ0n(:,sub2ind(iterations,csdz,csigma));
        S_VZ0(:,csdz,csigma)=S_VZ0n(:,sub2ind(iterations,csdz,csigma));
        S_KMAX(:,csdz,csigma)=S_KMAXn(:,sub2ind(iterations,csdz,csigma));
        S_KBAR(:,csdz,csigma)=S_KBARn(:,sub2ind(iterations,csdz,csigma));
        S_KSD(:,csdz,csigma)=S_KSDn(:,sub2ind(iterations,csdz,csigma));
        S_KSSM(:,csdz,csigma)=S_KSSMn(:,sub2ind(iterations,csdz,csigma));
        S_KSSNO(:,csdz,csigma)=S_KSSNOn(:,sub2ind(iterations,csdz,csigma));
        S_NBAR(:,csdz,csigma)=S_NBARn(:,sub2ind(iterations,csdz,csigma));
        S_NSD(:,csdz,csigma)=S_NSDn(:,sub2ind(iterations,csdz,csigma));
        S_CBAR(:,csdz,csigma)=S_CBARn(:,sub2ind(iterations,csdz,csigma));
        S_CSD(:,csdz,csigma)=S_CSDn(:,sub2ind(iterations,csdz,csigma));
        S_IBAR(:,csdz,csigma)=S_IBARn(:,sub2ind(iterations,csdz,csigma));
        S_ISD(:,csdz,csigma)=S_ISDn(:,sub2ind(iterations,csdz,csigma));
        S_YBAR(:,csdz,csigma)=S_YBARn(:,sub2ind(iterations,csdz,csigma));
        S_YSD(:,csdz,csigma)=S_YSDn(:,sub2ind(iterations,csdz,csigma));
        S_V1(:,csdz,csigma)=S_V1n(:,sub2ind(iterations,csdz,csigma));
        S_V2(:,csdz,csigma)=S_V2n(:,sub2ind(iterations,csdz,csigma));
        S_V3(:,csdz,csigma)=S_V3n(:,sub2ind(iterations,csdz,csigma));
        S_V4(:,csdz,csigma)=S_V4n(:,sub2ind(iterations,csdz,csigma));
        S_ITER(:,csdz,csigma)=S_ITERn(:,sub2ind(iterations,csdz,csigma));
        S_TRASHIT(:,csdz,csigma)=S_TRASHITn(:,sub2ind(iterations,csdz,csigma));
        S_TRASHMAX(:,csdz,csigma)=S_TRASHMAXn(:,sub2ind(iterations,csdz,csigma));
        V(:,:,csdz,csigma)=Vn(:,:,sub2ind(iterations,csdz,csigma));
        Ktildeposition(:,:,csdz,csigma)=Ktildepositionn(:,:,sub2ind(iterations,csdz,csigma));
        Ktilde(:,:,csdz,csigma)=Ktilden(:,:,sub2ind(iterations,csdz,csigma));
        N(:,:,csdz,csigma)=Nn(:,:,sub2ind(iterations,csdz,csigma));
        C(:,:,csdz,csigma)=Cn(:,:,sub2ind(iterations,csdz,csigma));
        Z(:,csdz,csigma)=Zn(:,sub2ind(iterations,csdz,csigma));
        S_KMIN(:,csdz,csigma)=S_KMINn(:,sub2ind(iterations,csdz,csigma));
        S_SDZ(:,csdz,csigma)=S_SDZn(:,sub2ind(iterations,csdz,csigma));
        S_SIGMA(:,csdz,csigma)=S_SIGMAn(:,sub2ind(iterations,csdz,csigma));
    end
end

clear S_KTILDEZ0n S_NOPTZ0n S_CZ0n S_VZ0n S_KMAXn S_KBARn S_KSDn S_KSSMn S_KSSNOn S_NBARn S_NSDn S_CBARn S_CSDn S_IBARn S_ISDn S_YBARn S_YSDn S_V1n S_V2n S_V3n S_V4n S_ITERn S_TRASHITn S_TRASHMAXn Vn Ktildepositionn Ktilden Nn Cn Zn S_KMINn S_SDZn S_SIGMAn

tv = datetime('now', 'Format','yyyy_MM_dd_HH_mm');
filename = sprintf('GRiskGHH4RCNfoc_%s.mat',tv);
save(fullfile(filename), '-v7.3'); % Saves results with timestamp

disp('SAVED & DONE!');

while max(S_TRASHMAX,[],'all') >  0
    error('XXXXXX Kmax too low, increase kmaxgp XXXXXX');
end

delete(mypool);





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plots below this point
% 
% % figure;
% % hold on;
% % plot(k, k, 'black-');
% % for i=1:lsigmavect
% %     plot(k, S_KTILDEZ0(1:gp,:,i),'color',colors{i,1});
% %     text(k(gp)*1.01, S_KTILDEZ0(gp,1,i), labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % %axis([0 2*kssTh 0 2*kssTh]);
% % title('Optimal capital rule, Z=E(Z)','FontSize',12);
% % xlabel('\it k_t','FontSize',12);
% % ylabel('\it k_{t+1}','FontSize',12,'rotation',0);
% % centeraxes(gca);
% % hold off;
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     plot(k,  S_CZ0(1:gp,:,i),'color',colors{i,1});
% %     text(k(gp)*1.01,  S_CZ0(gp,1,i), labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % %axis([0 2*kssTh 0 2.1]);
% % title('Optimal consumption rule, Z=E(Z)','FontSize',12);
% % xlabel('\it k_t','FontSize',12);
% % ylabel('\it c_t','FontSize',12,'rotation',0);
% % centeraxes(gca);
% % hold off;
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     p{i}=plot(k(1:gp), S_NOPTZ0(1:gp,:,i),'color',colors{i,1},'linewidth',.9);
% % end
% % % axis([-1 1.45*kssTh -.007 1.001]);
% % title('Optimal work effort rule, Z=E(Z)','FontSize',12);
% % xlabel('\it k_t','FontSize',12);
% % ylabel('\it n_t','FontSize',12,'rotation',0);
% % centeraxes(gca);
% % legend([p{6}(1,1), p{5}(1,1), p{4}(1,1), p{3}(1,1), p{2}(1,1), p{1}(1,1)], {[char(949) '=0'], [char(949) '=.5'], [char(949) '=1'], [char(949) '=2'], [char(949) '=3'], [char(949) '=4']});
% % hold off;
% 
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     plot(k, S_VZ0(1:gp,:,i),'color',colors{i,1});
% %     text(k(gp)*1.01, S_VZ0(gp,1,i), labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % title('Expected life-time Utility v(k,E(z)','FontSize',12);
% % xlabel('\it k_t','FontSize',12);
% % ylabel('\it v_t','FontSize',12,'rotation',0);
% % centeraxes(gca);
% % hold off;
% 
% 
% 
% % %Optimal decison rules and v, for different risk levels
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     plot(k, (S_KTILDEZ0(1:gp,:,i)./S_KTILDEZ0(1:gp,lsdzvect,i)-1)*100,'color',colors{i,1});
% %     text(k(gp)*1.01, (S_KTILDEZ0(gp,1,i)./S_KTILDEZ0(gp,lsdzvect,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % %axis([0 95 -1 5]);
% % title('Optimal capital rule, Z=E(Z), relative to no risk','FontSize',12);
% % xlabel('\it k_t','FontSize',12);
% % centeraxes(gca);
% % hold off;
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     plot(k, (S_CZ0(1:gp,:,i)./S_CZ0(1:gp,lsdzvect,i)-1)*100,'color',colors{i,1});
% %     text(k(gp)*1.01, (S_CZ0(gp,1,i)./S_CZ0(gp,lsdzvect,i)-1)*100 ,labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % %axis([0 95 -10 10]);
% % title('Optimal consumption rule, Z=E(Z), relative to no risk','FontSize',12);
% % xlabel('\it k_t','FontSize',12);
% % centeraxes(gca);
% % hold off;
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     plot(k, (S_VZ0(1:gp,:,i)-S_VZ0(1:gp,lsdzvect,i))./abs(S_VZ0(1:gp,lsdzvect,i))*100,'color',colors{i,1});
% %     text(k(gp)*1.01,(S_VZ0(gp,1,i)-S_VZ0(gp,lsdzvect,i))./abs(S_VZ0(gp,lsdzvect,i))*100, labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % ytickformat('%.1f')
% % %axis([0 475 -.7 1.3]);
% % title('Expected life-time Utility v(k,z), relative to no risk','FontSize',12);
% % xlabel('\it k_t','FontSize',12);
% % centeraxes(gca);
% % hold off;
% 
% % i=6; %{1, 2,..., 6) varepsilon levels
% % figure;
% % hold on;
% % for j=1:lsdzvect
% %     colors{j,1}=[.1*j .45 .74-.1*j];
% %     plot(k, (S_VZ0(1:gp,j,i)-S_VZ0(1:gp,lsdzvect,i))./abs(S_VZ0(1:gp,lsdzvect,i))*100,'color',colors{j,1});
% % end
% % %ytickformat('%.0f')
% % %axis([-.3 50 -8 12.7]);
% % title(['Expected life-time Utility v(k,E(z)), relative to no risk, ' char(949) '=' num2str(sigmavect(i),1)],'FontSize',12);
% % xlabel('\it k_t','FontSize',12);
% % centeraxes(gca);
% % hold off;
% 
% 
% 
% % Simulated variables vs risk, different sigma values
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     plot(S_SDZ(:,:,i), (S_KSSM(:,:,i)/S_KSSM(:,lsdzvect,i)-1)*100,'color',colors{i,1});
% %     text(S_SDZ(:,1,i)*1.01, (S_KSSM(:,1,i)/S_KSSM(:,lsdzvect,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % %axis([-.001 .1 -.05 2.9]);
% % ytickformat('%.1f')
% % title('SS capital computational, percentage deviation from no risk','FontSize',12);
% % xlabel('\sigma_z','FontSize',14);
% % centeraxes(gca);
% % hold off;
% 
% figure;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), (S_KBAR(:,:,i)/S_KBAR(:,lsdzvect,i)-1)*100, 'color',colors{i,1},'linewidth',1.5);
%     text(S_SDZ(:,1,i)*1.01, (S_KBAR(:,1,i)/S_KBAR(:,lsdzvect,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',12);
% end
% %axis([-.02 .25 -.0 50]);
% %ytickformat('%.0f')
% title('Mean capital $\bar k$ and investment, percentage deviation from no risk','Interpreter','Latex','FontSize',14);
% xlabel('\sigma','FontSize',14);
% centeraxes(gca);
% hold off;
% 
% figure;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), (S_CBAR(:,:,i)/S_CBAR(:,lsdzvect,i)-1)*100, 'color',colors{i,1},'linewidth',1.5);
%     text(S_SDZ(:,1,i)*1.01, (S_CBAR(:,1,i)/S_CBAR(:,lsdzvect,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',12);
% end
% %axis([-.002 .25 -.075 8]);
% %ytickformat('%.0f')
% title('Mean consumption, percentage deviation from no risk','Interpreter','Latex','FontSize',14);
% xlabel('\sigma','FontSize',14);
% centeraxes(gca);
% hold off;
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     plot(S_SDZ(:,:,i), (S_IBAR(:,:,i)/S_IBAR(:,lsdzvect,i)-1)*100, 'color',colors{i,1});
% %     text(S_SDZ(:,1,i)*1.01, (S_IBAR(:,1,i)/S_IBAR(:,lsdzvect,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % %axis([-.002 .25 -.075 8]);
% % ytickformat('%.0f')
% % title('Mean investment, percentage deviation from no risk (\sigma=0)','FontSize',12);
% % xlabel('\sigma','FontSize',14);
% % centeraxes(gca);
% % hold off;
% 
% figure;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), (S_YBAR(:,:,i)/S_YBAR(:,lsdzvect,i)-1)*100, 'color',colors{i,1},'linewidth',1.5);
%     text(S_SDZ(:,1,i)*1.01, (S_YBAR(:,1,i)/S_YBAR(:,lsdzvect,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',12);
% end
% %axis([-.002 .25 -.075 8]);
% %ytickformat('%.0f')
% title('Mean production, percentage deviation from no risk','Interpreter','Latex','FontSize',14);
% xlabel('\sigma','FontSize',14);
% centeraxes(gca);
% hold off;
% 
% figure;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), (S_NBAR(:,:,i)/S_NBAR(:,lsdzvect,i)-1)*100, 'color',colors{i,1},'linewidth',1.5);
%     text(S_SDZ(:,1,i)*1.01, (S_NBAR(:,1,i)/S_NBAR(:,lsdzvect,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',12);
% end
% %axis([-.002 .2 -.9 1.6]);
% %ytickformat('%.1f')
% title('Mean labor, percentage deviation from no risk','Interpreter','Latex','FontSize',14);
% xlabel('\sigma','FontSize',14);
% centeraxes(gca);
% hold off;
% 
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     plot(S_SDZ(:,:,i), (S_V1(:,:,i)-S_V1(:,lsdzvect,i))/abs(S_V1(:,lsdzvect,i))*100, 'color',colors{i,1});
% %     text(S_SDZ(:,1,i)*1.01, (S_V1(:,1,i)-S_V1(:,lsdzvect,i))/abs(S_V1(:,lsdzvect,i))*100, labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % %axis([-.002 .2 -8.7 1.5]);
% % ytickformat('%.0f')
% % title('v(kss,E(z)), percentage deviation from no risk','FontSize',12);
% % xlabel('\sigma','FontSize',14);
% % centeraxes(gca);
% % hold off;
% 
% 
% figure;
% hold on;
% cssGainPercentVkss=nan(lsigmavect,lsdzvect);
% for i=1:lsigmavect
%     for j=1:lsdzvect
%         cssGainPercentVkss(i,j)=((((1-beta)*(1-sigmavect(i))*(S_V1(1,j,i)-S_V1(1,lsdzvect,i))+(cssTh-psi*nssTh^theta/theta)^(1-sigmavect(i)))^(1/(1-sigmavect(i)))+psi*nssTh^theta/theta)/cssTh-1)*100;
%     end
% end
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), cssGainPercentVkss(i,:), 'color',colors{i,1},'linewidth',1.5);
%     text(S_SDZ(:,1,i)*1.01, cssGainPercentVkss(i,1), labels{i,1},'color',colors{i,1},'FontSize',12);
% end
% %axis([-.001 .2 -20 3]);
% %ytickformat('%.0f')
% title('$x\%, c_{ss}$ equivalent gains, $v(k_{ss},E(z))$','Interpreter','Latex','FontSize',14);
% xlabel('\sigma','FontSize',14);
% centeraxes(gca);
% hold off;
% 
% % % css equivalent of going from sdz=.0 to .1 with captial level k and z0
% % xss=.5;
% % figure;
% % hold on;
% % cssGainPercentVk=nan(gp, lsigmavect);
% % Dist1=abs(k-kssTh);
% % idk1=round(find(Dist1==min(Dist1))); % theoretical kss
% % dVZ0=S_VZ0(:,round(lsdzvect/2),:)-S_VZ0(:,lsdzvect,:);
% % dVZ0=permute(dVZ0,[1 3 2]);
% % for i=1:lsigmavect
% %         cssGainPercentVk(:,i)=((((1-beta)*(1-sigmavect(i))*(dVZ0(:,i))+(cssTh-psi*nssTh^theta/theta).^(1-sigmavect(i))).^(1/(1-sigmavect(i)))+psi*nssTh^theta/theta)/cssTh-1)*100;
% % end
% % for i=1:lsigmavect
% %     plot(k(3:round(idk1*xss)), cssGainPercentVk(3:round(idk1*xss),i), 'color',colors{i,1},'linewidth',1.5);
% %     text(k(round(idk1*xss))*1.0135, cssGainPercentVk(round(idk1*xss),i), labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % %axis([-.25 50 -11 1.8]);
% % %ytickformat('%.0f')
% % xlabel('\it k','FontSize',14);
% % %xlabel('k','FontName','mwa_cmmi10','FontSize',14);
% % centeraxes(gca);
% % title('$x\%, c_{ss}$ equivalent gains, $v(k,E(z))$','Interpreter','Latex','FontSize',14);
% % hold off;
% 
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     plot(S_SDZ(:,:,i), (S_V2(:,:,i)-S_V2(:,lsdzvect,i))/abs(S_V2(:,lsdzvect,i))*100, 'color',colors{i,1});
% %     text(S_SDZ(:,1,i)*1.01, (S_V2(:,1,i)-S_V2(:,lsdzvect,i))/abs(S_V2(:,lsdzvect,i))*100, labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % %axis([-.001 .1 -.05 2.9]);
% % ytickformat('%.1f')
% % title('E(v(kss,z)), percentage deviation from no risk','FontSize',12);
% % xlabel('\sigma','FontSize',14);
% % centeraxes(gca);
% % hold off;
% 
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     plot(S_SDZ(:,:,i), (S_V3(:,:,i)-S_V3(:,lsdzvect,i))/abs(S_V3(:,lsdzvect,i))*100, 'color',colors{i,1});
% %     text(S_SDZ(:,1,i)*1.01, (S_V3(:,1,i)-S_V3(:,lsdzvect,i))/abs(S_V3(:,lsdzvect,i))*100, labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % %axis([-.001 .1 -.05 2.9]);
% % ytickformat('%.1f')
% % title('v(kbar,E(z)), percentage deviation from no risk','FontSize',12);
% % xlabel('\sigma','FontSize',14);
% % centeraxes(gca);
% % hold off;
% 
% 
% figure;
% hold on;
% cssGainPercentVkbar=nan(lsigmavect,lsdzvect);
% for i=1:lsigmavect
%     for j=1:lsdzvect
%         cssGainPercentVkbar(i,j)=((((1-beta)*(1-sigmavect(i))*(S_V3(1,j,i)-S_V3(1,lsdzvect,i))+(cssTh-psi*nssTh^theta/theta)^(1-sigmavect(i)))^(1/(1-sigmavect(i)))+psi*nssTh^theta/theta)/cssTh-1)*100;
%     end
% end
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), cssGainPercentVkbar(i,:), 'color',colors{i,1},'linewidth',1.5);
%     text(S_SDZ(:,1,i)*1.01, cssGainPercentVkbar(i,1), labels{i,1},'color',colors{i,1},'FontSize',12);
% end
% %axis([-.001 .2 -20 3]);
% %ytickformat('%.0f')
% title('$x\%, c_{ss}$ equivalent gains, $v(\bar{k},E(z))$','Interpreter','Latex','FontSize',14);
% xlabel('\sigma','FontSize',14);
% centeraxes(gca);
% hold off;
% 
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     plot(S_SDZ(:,:,i), (S_V4(:,:,i)-S_V4(:,lsdzvect,i))/abs(S_V4(:,lsdzvect,i))*100, 'color',colors{i,1},'linewidth',1.5);
% %     text(S_SDZ(:,1,i)*1.01, (S_V4(:,1,i)-S_V4(:,lsdzvect,i))/abs(S_V4(:,lsdzvect,i))*100, labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % %axis([-.001 .1 -.05 2.9]);
% % ytickformat('%.1f')
% % title('E(v(kbar,z)), percentage deviation from no risk','FontSize',12);
% % xlabel('\sigma','FontSize',14);
% % centeraxes(gca);
% % hold off;
% 
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     plot(S_SDZ(:,:,i), (S_CSD(:,:,i)./S_CBAR(:,:,i))*100, 'color',colors{i,1},'linewidth',1.5);
% %     text(S_SDZ(:,1,i)*1.01, (S_CSD(:,1,i)./S_CBAR(:,1,i))*100, labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % %axis([-.001 .25 -.05 3.15]);
% % %ytickformat('%.1f')
% % title('SD consumption, percentage deviation from mean','Interpreter','Latex','FontSize',14);
% % xlabel('\sigma','FontSize',14);
% % centeraxes(gca);
% % hold off;
% 
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     plot(S_SDZ(:,:,i), (S_ISD(:,:,i)./S_IBAR(:,:,i))*100, 'color',colors{i,1},'linewidth',1.5);
% %     text(S_SDZ(:,1,i)*1.01, (S_ISD(:,1,i)./S_IBAR(:,1,i))*100, labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % %axis([-.001 .25 -.05 3.15]);
% % ytickformat('%.1f')
% % title('SD investment, percentage deviation from trend','FontSize',12);
% % xlabel('\sigma','FontSize',14);
% % centeraxes(gca);
% % hold off;
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     plot(S_SDZ(:,:,i), (S_YSD(:,:,i)./S_YBAR(:,:,i))*100, 'color',colors{i,1},'linewidth',1.5);
% %     text(S_SDZ(:,1,i)*1.01, (S_YSD(:,1,i)./S_YBAR(:,1,i))*100, labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % %axis([-.001 .25 -.05 3.15]);
% % ytickformat('%.1f')
% % title('SD GDP, percentage deviation from trend','FontSize',12);
% % xlabel('\sigma','FontSize',14);
% % centeraxes(gca);
% % hold off;
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     plot(S_SDZ(:,:,i), (S_KSD(:,:,i)./S_KBAR(:,:,i))*100, 'color',colors{i,1},'linewidth',1.5);
% %     text(S_SDZ(:,1,i)*1.01, (S_KSD(:,1,i)./S_KBAR(:,1,i))*100, labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % %axis([-.001 .25 -.05 3.15]);
% % ytickformat('%.1f')
% % title('SD k, percentage deviation from trend','FontSize',12);
% % xlabel('\sigma','FontSize',14);
% % centeraxes(gca);
% % hold off;
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     plot(S_SDZ(:,:,i), (S_NSD(:,:,i)./S_NBAR(:,:,i))*100, 'color',colors{i,1},'linewidth',1.5);
% %     text(S_SDZ(:,1,i)*1.01, (S_NSD(:,1,i)./S_NBAR(:,1,i))*100, labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % %axis([-.001 .25 -.05 3.15]);
% % %ytickformat('%.1f')
% % title('SD labor, percentage deviation from mean','Interpreter','Latex','FontSize',14);
% % xlabel('\sigma','FontSize',14);
% % centeraxes(gca);
% % hold off;
% 
% 
% figure;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), S_CSDCOMPDATA(:,:,i), 'color',colors{i,1},'linewidth',1.5);
%     text(S_SDZ(:,1,i)*1.01, S_CSDCOMPDATA(:,1,i), labels{i,1},'color',colors{i,1},'FontSize',12);
% end
% %text(.01,3.25, labels{i,1},'color',colors{i,1},'FontSize',12);
% plot(S_SDZ(:,:,i), 2.755*ones(lsdzvect),'color',[.5 .5 .5],'LineStyle','--');
% plot(S_SDZ(:,:,i), 1.3754*ones(lsdzvect),'color',[.5 .5 .5],'LineStyle',':');
% plot(S_SDZ(:,:,i), .2912*ones(lsdzvect),'color',[.5 .5 .5],'LineStyle','--');
% axis([-.001 .1 -.07 3.5]);
% %ytickformat('%.0f')
% title('SD consumption, HP-filtered','FontSize',12);
% xlabel('\sigma','FontSize',14);
% centeraxes(gca);
% hold off;
% % In data max 1.79 (by decade), 2.755 for rolling window (5), 2 for 10y
% % For whole period 1.3754 
% 
% 
% figure;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), S_NSDCOMPDATA(:,:,i), 'color',colors{i,1},'linewidth',1.5);
%     text(S_SDZ(:,1,i)*1.01, S_NSDCOMPDATA(:,1,i), labels{i,1},'color',colors{i,1},'FontSize',12);
% end
% plot(S_SDZ(:,:,i), 3.9455*ones(lsdzvect),'color',[.5 .5 .5],'LineStyle','--');
% plot(S_SDZ(:,:,i), 2.0870*ones(lsdzvect),'color',[.5 .5 .5],'LineStyle',':');
% plot(S_SDZ(:,:,i), .3388*ones(lsdzvect),'color',[.5 .5 .5],'LineStyle','--');
% %axis([-.001 .1 -.05 4.4]);
% %ytickformat('%.1f')
% title('SD work effort, HP-filtered','FontSize',12);
% xlabel('\sigma','FontSize',14);
% centeraxes(gca);
% hold off;
% 
% 
% figure;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), S_YSDCOMPDATA(:,:,i), 'color',colors{i,1},'linewidth',1.5);
%     text(S_SDZ(:,1,i)*1.01, S_YSDCOMPDATA(:,1,i), labels{i,1},'color',colors{i,1},'FontSize',12);
% end
% plot(S_SDZ(:,:,i), 2.7632*ones(lsdzvect),'color',[.5 .5 .5],'LineStyle','--');
% plot(S_SDZ(:,:,i), 1.6411*ones(lsdzvect),'color',[.5 .5 .5],'LineStyle',':');
% plot(S_SDZ(:,:,i), .3983*ones(lsdzvect),'color',[.5 .5 .5],'LineStyle','--');
% %axis([-.001 .1 -.05 4.4]);
% %ytickformat('%.1f')
% title('SD output, HP-filtered','FontSize',12);
% xlabel('\sigma','FontSize',14);
% centeraxes(gca);
% hold off;
% 
% 
% % figure;
% % hold on;
% % for i=1:lsigmavect
% %     plot(S_SDZ(:,:,i), S_KMAX(:,:,i), 'color',colors{i,1},'linewidth',1.5);    
% %     %plot(S_SDZ(:,:,i), S_KMIN(:,:,i), 'color',colors{i,1});
% %     text(S_SDZ(:,1,i)*1.01, S_KMAX(:,1,i), labels{i,1},'color',colors{i,1},'FontSize',12);
% % end
% % %plot(S_SDZ(:,:,i), max(k2)*ones(lsdzvect), 'color',[.5 .5 .5]);
% % %axis([-.001 .25 -.05 3.15]);
% % ytickformat('%.1f')
% % title('Max K in simulations','FontSize',12);
% % xlabel('\sigma','FontSize',14);
% % centeraxes(gca);
% % hold off;
