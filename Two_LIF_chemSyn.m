
%This program is generating two LIF neurons which are connected through a 
%chemical synapse. It solves LIF equation as follows:
%     tau * v_dot = -(v - v_rest) + I_ext - I_syn,
%where I_ext is the external input, which can be a sensory stimulus or any
%input from the population activity of other neurons. v_rest is the
%membrane potential of resting state. I_syn is accumulative synaptic inputs 
%arrived from all the pre-synaptic neurons to the given post-synaptic one.
%v_th is the threshold potential. tau is the time constant of the neurons.
%alpha, beta, and N_duration characterize the post-synaptic potential. The 
%first elements of alpha and beta correspond to negative I_syn, 
%therefore to EPSP, and the second ones correspond to positive I_syn, and so to
%IPSP. N is equivalent to concentration of neurotransmitters in synaptic 
%cleft which rises during synaptic activity (Action Potential). It is a 
%rectangular pulse. S represents the fraction of bound receptors. 
%alpha and beta and N specify the shape of S.
%g is conductivity matrix; its elements demonstrate maximal conductance of 
%the synapses. Whether a synapse is excitatory or inhibitory is 
%specified using E_syn, the synaptic reversal potential (E_syn).
%For inhibitory synapses, E_syn equals to E_syn_inh, and for excitatory 
%ones the value is set to zero.
%This program returns 4 graphs. Figure(1) depicts the sub-threshold
%dynamics of neurons' membrane potential. Figure(2) is raster plot of the
%two neurons. Figure(3) represents the dynamics of N and S,
%which determine the input I_syn. You can produce EPSP and IPSP
%corresponding to NMDA, GABAa, GABAb, etc. by setting relevant values for 
%for alpha and beta and therefore for the shape the S profile.
%Figure(4) demonstrates the accumulative synaptic inputs from pre-synaptic 
%neurons to the post-synaptic one. Note that negative I_syn corresponds 
%to EPSP and positive I_syn corresponds to IPSP (Consider the negative 
%sign before I_syn in the LIF equation). 
%


close all; clc; clear all

no_neurons = 2;
v_th = -65;
v_a0 = -70;
v_b0 = -65.1;
dt = .01;
v_rest = -70;
spk_pick= 40;
N_duration = 1;
t_f = 150;

%tau = rand (no_neurons,1)*10+10; 
tau = [10,12];
alpha = [.6,1.5];
beta = [.3,.3];


%E_syn = zeros(2,2);
E_syn_inh = -80;
E_syn = [0,0; E_syn_inh,0];

g=[0,.01; .02,0];





n_tSteps = t_f/dt +1;
V = zeros(n_tSteps,no_neurons);
V(1,1:no_neurons) = rand(1,no_neurons)*4 +v_rest;
V(1,1:2) = [v_a0,v_b0];


spike_train = zeros(n_tSteps,no_neurons);
T = zeros(n_tSteps,1);

S = zeros(no_neurons, no_neurons,n_tSteps);
N = zeros(n_tSteps,no_neurons);
I_synps = zeros(n_tSteps,no_neurons);

I_ext = [v_th - v_rest + .01,v_th - v_rest + .01];
t = 0:dt:t_f;


for tStep=1:n_tSteps -1
    
    for j=1:no_neurons
        v_a1 = V(tStep,j);
        
        [i_synps,s] = I_chem_synps(j,tStep,g,S,N,alpha,beta,dt,E_syn,v_a1,E_syn_inh);
        I_synps (tStep,j) = i_synps;
        S(j,:,tStep+1) = s;
        
        %i_synps = 0;
        [v_a2,spk] = LIF_ODE(v_th, v_rest, tau(j), dt, I_ext(j), i_synps, v_a1 );
        V(tStep+1,j) = v_a2;
        
        if spk == true
            spike_train(tStep,j) = 1;
            n = rectpuls(t-T(tStep)-.5*N_duration,N_duration);
            n = n';
            N(:,j) = N(:,j)+n;
        end
         
    end    
    T(tStep+1) = T(tStep)+dt;
  
end



figure(1);  
plot(T,V)
title(' Dynamics of Membrane Potential')
xlabel('Time')
ylabel('V')
legend ('Neuron 1' , 'Neuron 2', 'location', 'east')


figure(2)
rasterPlot(spike_train,T,no_neurons)


figure (3)
p1 = plot (T,N(:,1));
hold on;
p2 = plot (T,N(:,2));
ss(:,1) = S(1,2,:);
p3 = plot (T,ss);
hold on
ss(:,1) = S(2,1,:);
p4 = plot (T,ss);
title('Model Prameters: Concentration of Transmitters, [N], & Fraction of Bound Receptors, S')
xlabel('Time')
ylabel('[N] & S')
legend([p1 p2 p3 p4],{'[N_1]', '[N_2]', 'S_{1,2}', 'S_{2,1}'})


figure(4)
plot(T,I_synps)
hold on 
%plot(T,I_ext)
% plot(T,-I_synps+I_ext)
title('Synaptic Input') 
xlabel('Time')
ylabel('I_{syn}')
legend('I_{syn,1}', 'I_{syn,2}')



% figure (5)
% plot (T,spike_train)
% title('Spike Train')