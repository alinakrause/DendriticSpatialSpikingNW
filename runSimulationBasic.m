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

jee=38;
jei=400; %350
jie=130; 
jii=450;

jee=18; %25  30 %35 niedriger?
jei=550; %360 350   357  360   370   376 380 %400   höher für niedriger isi aber zufällig? 
jei=360;
jie=140;%119  120  130  140  150    130; %120  niedriger?
jii=450;%500 50  500  600;
jii=400;

jee=40;   % balanced LIF
jei=370; % 370 300 370  470 445  345
jie=120; %105 110 90  niedriger mehr clustering 
jii=30; %100

jee=55;   % same firing rates
jei=420; % 370 300 370  470 445  345
jie=170; %105 110 90  niedriger mehr clustering 
jii=20; %100



jee=55;   % same firing rates
jei=320; % 370 300 370  470 445  345
jie=170; %105 110 90  niedriger mehr clustering 
jii=40; %100

% balanced IF
jee=40;%40
jei=300; %290 
jie=120; 
jii=400;%400;


% balanced LIF
jee=120;%40
jei=900; %290 
jie=360; 
jii=1200;%400;
jeF=240;%80 120;
jiF=240;%120;
%ffwd 0.1



% 
% % % % EIF
% jee=40;
% jei=400;
% jie=120; 
% jii=400;
% 


% Effective connection strengths normalized by sqrt(N)
% and averaged over all distances
wee0=jee*pee0*q;
wei0=jei*pei0*(1-q);
wie0=jie*pie0*q;
wii0=jii*pii0*(1-q);

% Same for feedforward connections
jeF=240;%80 120;
jiF=240;%120;
% jeF = 70;
% jiF = 70;
rF=.01;%.006; .007;%.01;%.005;%.0035;%.01;%.005     um cv ändern aber mehr synchron
peF0=.25;
piF0=.08;
%piF0=.25;
weF0=jeF*peF0*(NF/N);
wiF0=jiF*piF0*(NF/N);

% Feedforward and recurrent
% connection widths
sigmaffwd=.1;%.1;%.5;%5;
sigmarec=.25;%.05;%.05; %.05; %.05;

% Time constants
%taujitter=0;
tausyne=6;
tausyni=5;
tausynF=tausyne+1;
% 
% tausyne = 5;
% tausyni = 5;
% tausynF = 5;
% Random Excitatory neurons to record 
% currents and voltages from
nrecord_inhibitory = 100;
record_min = Ne1+1
record_max = Ne1+Ni1-1
Irecord_i = randi([record_min,record_max],2,nrecord_inhibitory)

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
vT=[-50 -50]; %mV
vl=[-60 -60]; % mV
vre=[-65 -65];
tref=[1.5 .5];

% Define a range of sigmarec values to iterate over
sigmarec_values = [0.05, 0.10,0.15,0.2,0.25]; % Example values
sigmarec_values = [0.05];
sigmarec_values = [0.05,0.1,0.2,0.3]
sigmarec_values = [0.05,0.1,0.2,0.3,0.4,0.5];
sigmarec_values = [0.05];
sigmarec_inhib = [1,2,3,4];
centers = [0,-2/4,-8/4,-14/4];
index = 1;
sigmarec_inhib = [1];
%inhib = 1;
for sigmarec = sigmarec_values
    for inhib= sigmarec_inhib
        NL_center = 0;%centers(index);
        %index = index + 1;
        
        % Width of recurrent connections
        sigmaee=sigmarec;
        sigmaie=sigmarec;
        sigmaei=sigmarec*inhib;
        sigmaii=sigmarec*inhib;
        sigmaeF=sigmaffwd;
        sigmaiF=sigmaffwd;
        
        % Length, time bin, and burn-in period
        % for simulation
        T=22000;
        dt=.1;
        Nt=round(T/dt);
        Tburn=2000;
        nburn=round(Tburn/dt);
        T = 4000;
        
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
        
        
        % Run simulation
        %[s,IF,Ie,Ii,~]=EIF2DSpatialNetworkNoJitter(sF,NF1,Ne1,Ni1,JeF,JiF,Jee,Jei,Jie,Jii,KeF,KiF,Kee,Kei,Kie,Kii,betaeF,betaiF,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynF,tausyne,tausyni,V0,T,dt,maxns,Irecord);
        %[s,IF,Ie,Ii,~]=EIF2DSpatialNetworkNonlinearities(sF,NF1,Ne1,Ni1,JeF,JiF,Jee,Jei,Jie,Jii,KeF,KiF,Kee,Kei,Kie,Kii,betaeF,betaiF,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynF,tausyne,tausyni,V0,T,dt,maxns,Irecord);
        %[s,IF,Ie,Ii,~,e,i,ce,ci] = IF2DSpatialNetwork(sF,NF1,Ne1,Ni1,JeF,JiF,Jee,Jei,Jie,Jii,KeF,KiF,Kee,Kei,Kie,Kii,betaeF,betaiF,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynF,tausyne,tausyni,V0,T,dt,maxns,Irecord,Irecord2);
        %[s,IF,Ie,Ii,~,f_total,e_total,i_total,total,count,input,ae1,ai1,ax1,ae2,ai2,ax2,ae3,ai3,ax3,ae4,ai4,ax4] = IF2DSpatialNetworkNonlinearities(sF,NF1,Ne1,Ni1,JeF,JiF,Jee,Jei,Jie,Jii,KeF,KiF,Kee,Kei,Kie,Kii,betaeF,betaiF,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynF,tausyne,tausyni,V0,T,dt,maxns,Irecord,Irecord_i);
        %[s,IF,Ie,Ii,~,f_total,e_total,i_total,total,count,input,ae1,ai1,ax1,ae2,ai2,ax2,ae3,ai3,ax3,ae4,ai4,ax4] = IF2D_NL_test(sF,NF1,Ne1,Ni1,JeF,JiF,Jee,Jei,Jie,Jii,KeF,KiF,Kee,Kei,Kie,Kii,betaeF,betaiF,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynF,tausyne,tausyni,V0,T,dt,maxns,Irecord,Irecord2);
        %[s,IF,Ie,Ii,~,f_total,e_total,i_total,total,count,input,ae] = IF2D_NL_test(sF,NF1,Ne1,Ni1,JeF,JiF,Jee,Jei,Jie,Jii,KeF,KiF,Kee,Kei,Kie,Kii,betaeF,betaiF,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynF,tausyne,tausyni,V0,T,dt,maxns,Irecord,Irecord2);
        %[s,IF,Ie,Ii,~,e,i,IeI,IiI,IfI] = IF2DSpatialNetwork(sF,NF1,Ne1,Ni1,JeF,JiF,Jee,Jei,Jie,Jii,KeF,KiF,Kee,Kei,Kie,Kii,betaeF,betaiF,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynF,tausyne,tausyni,V0,T,dt,maxns,Irecord,Irecord_i);
        
        
        dendrites = true;
        NL_center = 0;
        NL_slope = 0.5;
        NL_height = 1.0;
    
        slopes = [0.1,0.7];
        slopes = [10];
        %slopes = [1, 5, 10, 50] %0.5,1, 5,
        %slopes = [5,10,50,100]
        slopes = [0.5, 1,2, 5];
        slopes = [10,20];
        slopes= [20];
        heights = [0.9,1.0,2.0,3.75,4.0];
        heights = [2.5, 3.75,5];
        heights = [7,10,15];
        heights = [100];
        %heights = [0.5,1]; %, 0.8
        %heights = [0.8, 1.2, 1.6, 2.0, 2.4];
        %heights = [50.0];
    
       
        for NL_height = heights
            for NL_slope = slopes
            
                if dendrites
                
                    NL_center = 0;
                    % NL_slope = 0.5;
                    % NL_height = 1.0;
                
                    
                    %last wotking version with dendrites: 
                    %[s,IF,Ie,Ii,v,f_total,e_total,i_total,total,count,input,ae1,ai1,ax1,ae2,ai2,ax2,ae3,ai3,ax3,ae4,ai4,ax4,NL,NL_after] = IF2DSpatialNetworkNonlinearities(sF,NF1,Ne1,Ni1,JeF,JiF,Jee,Jei,Jie,Jii,KeF,KiF,Kee,Kei,Kie,Kii,betaeF,betaiF,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynF,tausyne,tausyni,V0,T,dt,maxns,Irecord,Irecord_i);
                    
                    %[s,IF,Ie,Ii,v,f_total,e_total,i_total,total,count,inpFut,ae1,ai1,ax1,ae2,ai2,ax2,ae3,ai3,ax3,ae4,ai4,ax4,NL,NL_after] = IF2DSpatialNetworkNonlinearities(sF,NF1,Ne1,Ni1,JeF,JiF,Jee,Jei,Jie,Jii,KeF,KiF,Kee,Kei,Kie,Kii,betaeF,betaiF,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynF,tausyne,tausyni,V0,T,dt,maxns,Irecord,Irecord_i,NL_center,NL_height,NL_slope);
                     
                    % exponential IF
                    %[s,IF,Ie,Ii,v,f_total,e_total,i_total,total,count,input,ae1,ai1,ax1,ae2,ai2,ax2,ae3,ai3,ax3,ae4,ai4,ax4,NL,NL_after] = EIF2DSpatialNetworkNonlinearities(sF,NF1,Ne1,Ni1,JeF,JiF,Jee,Jei,Jie,Jii,KeF,KiF,Kee,Kei,Kie,Kii,betaeF,betaiF,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynF,tausyne,tausyni,V0,T,dt,maxns,Irecord,Irecord_i,NL_center,NL_height,NL_slope);
                        
                
                
                    %only dendrites without NLs
                    [s,IF,Ie,Ii,v,f_total,e_total,i_total,total,count,input,ae1,ai1,ax1,ae2,ai2,ax2,ae3,ai3,ax3,ae4,ai4,ax4,NL,NL_after] = IF2DSpatialNetworkDendrites(sF,NF1,Ne1,Ni1,JeF,JiF,Jee,Jei,Jie,Jii,KeF,KiF,Kee,Kei,Kie,Kii,betaeF,betaiF,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynF,tausyne,tausyni,V0,T,dt,maxns,Irecord,Irecord_i);
                
                
                    % attempt indexing:
                    %[s,IF,Ie,Ii,~,f_total,e_total,i_total,total,count,input,ae1,ai1,ax1,ae2,ai2,ax2,ae3,ai3,ax3,ae4,ai4,ax4] = IF2D_NL_test(sF,NF1,Ne1,Ni1,JeF,JiF,Jee,Jei,Jie,Jii,KeF,KiF,Kee,Kei,Kie,Kii,betaeF,betaiF,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynF,tausyne,tausyni,V0,T,dt,maxns,Irecord,Irecord_i);
                else
                    %no dendrites
                    [s,IF,Ie,Ii,~]=EIF2DSpatialNetworkNoJitter(sF,NF1,Ne1,Ni1,JeF,JiF,Jee,Jei,Jie,Jii,KeF,KiF,Kee,Kei,Kie,Kii,betaeF,betaiF,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynF,tausyne,tausyni,V0,T,dt,maxns,Irecord);
                
                    %[s,IF,Ie,Ii,v,e,i,ce,ci] = IF2DSpatialNetwork(sF,NF1,Ne1,Ni1,JeF,JiF,Jee,Jei,Jie,Jii,KeF,KiF,Kee,Kei,Kie,Kii,betaeF,betaiF,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynF,tausyne,tausyni,V0,T,dt,maxns,Irecord,Irecord_i);
                
                end
                
                % Get rid of padded zeros
                s=s(:,s(2,:)~=0);
                
                % Average firig rates from simulation
                reSim=1000*nnz(s(1,:)>Tburn & s(2,:)>0)/(Ne*(T-Tburn));
                riSim=1000*nnz(s(1,:)>Tburn & s(2,:)<0)/(Ni*(T-Tburn));
                
                disp(sprintf('\nFiring rates from sim: re=%.2f, ri=%.2f Hz',reSim,riSim))
                
                % How long did the simulation take
                t0=toc;
                disp(sprintf('\nSimulation time for Fig3: %.1f min',t0/60))
                
                
                % Consolidate all variables into a structure
                simulation_data = struct();
                
                % Core data
                simulation_data.s = s;
                simulation_data.IF = IF;
                simulation_data.Ie = Ie;
                simulation_data.Ii = Ii;
                % simulation_data.v = v;
                % 
                % % % Dendritic data
                % 
                % simulation_data.ae1 = ae1;
                % simulation_data.ai1 = ai1;
                % simulation_data.ax1 = ax1;
                % simulation_data.ae2 = ae2;
                % simulation_data.ai2 = ai2;
                % simulation_data.ax2 = ax2;
                % simulation_data.ae3 = ae3;
                % simulation_data.ai3 = ai3;
                % simulation_data.ax3 = ax3;
                % simulation_data.ae4 = ae4;
                % simulation_data.ai4 = ai4;
                % simulation_data.ax4 = ax4;
                % 
                % simulation_data.NL = NL;
                % simulation_data.NL_after = NL_after;
            
                
                filename = sprintf('slopes_simulation_sigmarec_%.3f_height_%.3f_slope_%.3f_inhib_%.3f.mat', sigmarec,NL_height,NL_slope,inhib);
                filename = sprintf('no_dendrites_simulation_sigmarec_%.3f_inhib_%.3f_.mat', sigmarec,inhib);
                filename = sprintf('test_simulation_sigmarec_%.3f_inhib_%.3f__regular.mat', sigmarec,inhib);
    
                filename = sprintf('/Volumes/TOSHIBA EXT/new_NL_dendrites_recurrent_largeNL/LIF_no_Dendrites_simulation_sigmarec_%.3f_inhib_%.3f_height_%.3f_slope_%.3f_.mat', sigmarec,inhib,NL_height,NL_slope);
               
                % Save the structure to the generated filename
                %save(filename, '-struct', 'simulation_data');
                % Save the structure to a single .mat file
                save('Simulation_rec_inib_height_slope.mat', '-struct', 'simulation_data');
            end
        end

    end
end



% 
% %disp(i_total)
% %disp(e_total)
% %disp(f_total)
% %disp(e)
% %disp(i)
% %disp(ce)
% %disp(ci)
% % Save simulation data
% %disp(count/N)
% save NetworkSimForFigure3.mat
% save IF.mat 
% save Ie.mat 
% save Ii.mat
% save v.mat
% save NL.mat
% save NL_after.mat
% %save ae1.mat
% % save IfI.mat
% % save IeI.mat
% % save IiI.mat
% save ae1.mat
% save ae2.mat
% save ae3.mat
% save ae4.mat
% save ai1.mat
% save ai2.mat
% save ai3.mat
% save ai4.mat
% save ax1.mat
% save ax2.mat
% save ax3.mat
% save ax4.mat


