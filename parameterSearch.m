% Clear all data and close all plots
clear
close all

rng(1)

% Time this script
tic


% Number of neurons in network along
% each dimension (e.g., the exc network
% is Ne1xNe1)
Ne1=200;
Ni1=100;
NF1=75;
Ne=Ne1*Ne1;
Ni=Ni1*Ni1;
NF=NF1*NF1;
N=Ne+Ni;

% Proportion of neurons that are excitatory
q=Ne/N;

% Recurrent connection probabilities averaged over 
% all distances (scale the Gaussians with these numbers)
% pab0 is the proportion of postsynaptic neurons in pop.
% a that receive input from a randomly chosen presynaptic
% neuron in pop. b
pee0=.05;
pei0=.05;
pie0=.05;
pii0=.05;

% Recurrent connection strength scalings
jee=40;   % same firing rates
jei=345;
jie=100 
jii=380;


%jee=33; %35   33 best CV isi
%jei=344;%342, 350
%jie=75;%75   nicht ändern
%jii=450;%380;   nicht ändern

jee=40;%40;   % test
jei=345;
jie=100;%120; 
jii=350;%400;

jee=40;   % old (paper)
jei=250;
jie=150; 
jii=300;
 

% Effective connection strengths normalized by sqrt(N)
% and averaged over all distances
wee0=jee*pee0*q;
wei0=jei*pei0*(1-q);
wie0=jie*pie0*q;
wii0=jii*pii0*(1-q);

% Same for feedforward connections
jeF=100;%120;
jiF=100;%120;
rF=.005;%.0035;%.01;%.005
peF0=.25;
piF0=.08;
piF0=.25
weF0=jeF*peF0*(NF/N);
wiF0=jiF*piF0*(NF/N);

% Feedforward and recurrent
% connection widths
sigmaffwd=.5;
sigmarec=.05; %.05;

% Time constants
%taujitter=0;
tausyne=6;
tausyni=5;
tausynF=tausyne+1;

tausyne = 5;
tausyni = 5;
tausynF = 5;
% Random Excitatory neurons to record 
% currents and voltages from
nrecord0=200;
Irecord=randi(Ne1,2,nrecord0);


% For balanced state to exist this vector should be decreasing
disp(sprintf('\nThis list should be decreasing for\n  a balanced state to exist: %.2f %.2f %.2f\n',weF0/wiF0,wei0/wii0,wee0/wie0));

% and these values should be >1
disp(sprintf('\nAlso, this number should be greater than 1: %.2f\n',wii0/wee0));

% Avg. firing rates in the balanced network limit
disp(sprintf('\nFiring rates for large N: re=%.2f, ri=%.2f Hz',rF*1000*(weF0*wii0-wiF0*wei0)/(wei0*wie0-wee0*wii0),rF*1000*(weF0*wie0-wiF0*wee0)/(wei0*wie0-wee0*wii0)))

% Neuron params
gl=[1/15 1/10];
Cm=[1 1];
vlb=[-100 -100];
vth=[-10 -10];
DeltaT=[2 .5]; 
DeltaT=[0 0];
vT=[-50 -50]; %mV
vl=[-60 -60]; % mV
vre=[-65 -65];
tref=[1.5 .5];

% Width of recurrent connections
sigmaee=sigmarec;
sigmaie=sigmarec;
sigmaei=sigmarec;
sigmaii=sigmarec;
sigmaeF=sigmaffwd;
sigmaiF=sigmaffwd;

% Length, time bin, and burn-in period
% for simulation
T=22000;
dt=.1;
Nt=round(T/dt);
Tburn=2000;
nburn=round(Tburn/dt);

% Widths of connections in terms of number
% of postsynaptic neurons
betaee=sigmaee*(Ne1);
betaei=sigmaei*(Ne1);
betaie=sigmaie*(Ni1);
betaii=sigmaii*(Ni1);
betaeF=sigmaeF*(Ne1);
betaiF=sigmaiF*(Ni1);

% Number of outgoing connections
Kee=pee0*Ne;
Kei=pei0*Ne;
Kie=pie0*Ni;
Kii=pii0*Ni;
KeF=peF0*Ne;
KiF=piF0*Ni;

% Actual connection strengths
Jee=(jee/sqrt(N));
Jei=-(jei/sqrt(N));
Jie=(jie/sqrt(N));
Jii=-(jii/sqrt(N));
JeF=jeF/sqrt(N);
JiF=jiF/sqrt(N);

% Generate Poisson spike times for the ffwd network
tempspikes=sort(T*rand(poissrnd(rF*NF*T),1));

% Store in the data format expected by the C code,
% where neuron indices are assigned randomly
sF=zeros(3,numel(tempspikes));
sF(1,:)=tempspikes;
clear tempspikes;
sF(2,:)=ceil(rand(1,size(sF,2))*NF1);
sF(3,:)=ceil(rand(1,size(sF,2))*NF1);

% Maximum number of spikes.
% Simulation will terminate with a warning if this is exceeded
maxns=N*T*.1;

% Random initial membrane potentials
V0min=vre(1);
V0max=vT(1);
V0=(V0max-V0min).*rand(N,1)+V0min;

% Number of neurons in network along each dimension
Ne1 = 200;
Ni1 = 100;
NF1 = 75;
Ne = Ne1 * Ne1;
Ni = Ni1 * Ni1;
NF = NF1 * NF1;
N = Ne + Ni;

% Proportion of neurons that are excitatory
q = Ne / N;

% Recurrent connection probabilities averaged over all distances
pee0 = .05;
pei0 = .05;
pie0 = .05;
pii0 = .05;

% Feedforward connection strength scalings
jeF = 100;
jiF = 100;

rF = .005;
peF0 = .25;
piF0 = .08;
piF0 = .25;

sigmaffwd = .5;
sigmarec = .05;

tausyne = 5;
tausyni = 5;
tausynF = 5;

% Random Excitatory neurons to record
nrecord0 = 200;
Irecord = randi(Ne1, 2, nrecord0);

% Length, time bin, and burn-in period for simulation
T = 22000;
dt = .1;
Nt = round(T / dt);
Tburn = 2000;
nburn = round(Tburn / dt);

% Generate Poisson spike times for the ffwd network
tempspikes = sort(T * rand(poissrnd(rF * NF * T), 1));

% Store in the data format expected by the C code
sF = zeros(3, numel(tempspikes));
sF(1, :) = tempspikes;
clear tempspikes;
sF(2, :) = ceil(rand(1, size(sF, 2)) * NF1);
sF(3, :) = ceil(rand(1, size(sF, 2)) * NF1);

maxns = N * T * .1;

V0min = -65;
V0max = -50;
V0 = (V0max - V0min) .* rand(N, 1) + V0min;

% Initialize variables for the grid search
best_total_input = Inf; % Start with a large value
best_params = [];

min = 50
max = 500

min = 70
max = 71

% Perform grid search over Kee, Kei, Kie, Kii, KeF, KiF
for Kee = min:1:max
    for Kei = min:1:max
        for Kie = min:1:max
            for Kii = min:1:max
                for KeF = min:1:max
                    for KiF = min:1:max
                        
                        % Actual connection strengths for current combination
                        Jee = Kee / sqrt(N);
                        Jei = -Kei / sqrt(N);
                        Jie = Kie / sqrt(N);
                        Jii = -Kii / sqrt(N);
                        JeF = jeF / sqrt(N);
                        JiF = jiF / sqrt(N);
                        
                        % Run simulation
                        [s,IF,Ie,Ii,~,f_total,e_total,i_total,total,count,input] = IF2DSpatialNetworkNonlinearities(sF,NF1,Ne1,Ni1,JeF,JiF,Jee,Jei,Jie,Jii,KeF,KiF,Kee,Kei,Kie,Kii,betaeF,betaiF,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynF,tausyne,tausyni,V0,T,dt,maxns,Irecord);

                        
                        % Compute firing rates from simulation
                        reSim = 1000 * nnz(s(1, :) > Tburn & s(2, :) > 0) / (Ne * (T - Tburn));
                        riSim = 1000 * nnz(s(1, :) > Tburn & s(2, :) < 0) / (Ni * (T - Tburn));
                        
                        % Evaluate the total input
                        current_total_input = abs(sum(input)); % Example: using total input
                        
                        % Check the additional condition: firing rates, excitatory, and inhibitory inputs must be non-zero
                        if reSim > 0 && riSim > 0 && e_total > 0 && i_total > 0
                            % Check if this is the best combination
                            if current_total_input < best_total_input
                                best_total_input = current_total_input;
                                best_params = [Kee, Kei, Kie, Kii, KeF, KiF];
                            end
                        end
                    end
                end
            end
        end
    end
end

% Display the best combination of parameters and the corresponding total input
disp('Best combination of parameters:');
disp(best_params);
disp('Minimum total input:');
disp(best_total_input);

% Time taken for the grid search
t0 = toc;
disp(sprintf('\nTotal time for grid search: %.1f min', t0 / 60));

% Save results
save BestParameters.mat best_params best_total_input
