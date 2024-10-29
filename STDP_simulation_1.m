close all
clear
% the number of excitatory and inhibitory neurons in a neural network
N = 100;
p = 0.2;
Ni = ceil(N*p);
Ne = N - Ni;

% parameter setting
unit = 1e-9;
delta_t = 1e-3;
neuron.index_e = 1:Ne;
neuron.index_i = Ne+1:N;
neuron.N = N;
neuron.Ne = Ne;
neuron.Ni = Ni;
neuron.C = 281 * 1e-12; %pF. membrane capacitance.
neuron.VT = -50.4 *1e-3; %mV. votage threshold for exponential function.
neuron.Vslop = 2* 1e-3; %mV, slop of the exponential.
neuron.Vcut = neuron.VT + 10e-3; %mV, votage threshold for spike generation.
neuron.Vpeak = 20 * 1e-3; %mV, peak votage of a spike.

neuron.Jstrength = 1e-7; 
neuron.I0 = 0.8;         

neuron.gL = [30*1e-9*ones(neuron.Ne,1);...
             30*1e-9*ones(neuron.Ni,1) + 3e-9 *rand(neuron.Ni,1)];
neuron.EL = [-65*1e-3*ones(Ne,1);-70.6 *1e-3*ones(neuron.Ni,1)+1e-2 *rand(neuron.Ni,1)];
neuron.a = [4*1e-9 *ones(neuron.Ne,1);...
            2*neuron.C/(144*1e-3)-0.4*1e-9 *rand(neuron.Ni,1)];  
neuron.b = [0.0805*1e-9 *ones(neuron.Ne,1);...
            zeros(neuron.Ni,1)];
neuron.tau_adapt =[144*1e-3*ones(neuron.Ne,1);...
                   144*1e-3 + 1e-2 *rand(neuron.Ni,1)];
neuron.Vreset = [-65*1e-3*ones(neuron.Ne, 1);...
                 neuron.EL(neuron.index_i)];

W  = 0.5*ones(neuron.N,neuron.N)+0.5*rand(neuron.N,neuron.N); %synaptic weights
% edge connection probability
connectivity = 0.3; 

Je = 1;
Ji = 1.5;
Jie = 1.1; 
W(neuron.index_e,neuron.index_e) = (rand(neuron.Ne,neuron.Ne)<connectivity).*W(neuron.index_e,neuron.index_e);
W(neuron.index_e,neuron.index_i) = (rand(neuron.Ne,neuron.Ni)<0.9).*W(neuron.index_e,neuron.index_i);
W(neuron.index_i,neuron.index_e) = (rand(neuron.Ni,neuron.Ne)<0.9).*W(neuron.index_i,neuron.index_e);
W(neuron.index_i,neuron.index_i) = (rand(neuron.Ni,neuron.Ni)<0.9).*W(neuron.index_i,neuron.index_i);
neuron.Je = Je /neuron.Ne/connectivity;
neuron.Ji = -Ji/neuron.Ni/0.9;
neuron.Jie = Jie/neuron.Ne/0.9;

neuron.S =...
    neuron.Jstrength * [neuron.Je*W(neuron.index_e,neuron.index_e),...
    neuron.Ji*W(neuron.index_e,neuron.index_i);...
    neuron.Jie*W(neuron.index_i,neuron.index_e), 1*neuron.Ji*rand(neuron.Ni,neuron.Ni)];

refrac_t =2*delta_t;  %refrac_t=0.002


%neuron.T =20.9; 
neuron.T =14.8; 
%neuron.T =4.5; 
%neuron.T =7.8; 
learning_start = 1;
learning_end = 21000;


u = zeros(neuron.N,1);
tSpan = 0:delta_t:neuron.T;
v = neuron.EL + 10*rand(N,1)*1e-3;
firings = [];
vt_neuron = zeros(N,length(tSpan));
last_spike = zeros(N,1);

his_current = zeros(size(vt_neuron));
% synaptic parameter
f = 0.5;
f_ampa = f;
f_gaba = f;
tau_amap = 4*delta_t;
tau_gaba = 10*delta_t;
connect = (neuron.S~=0);
% g_ampa = zeros(neuron.N,neuron.Ne);
% g_gaba = zeros(neuron.N,neuron.Ni);
r_e = zeros(neuron.N,length(tSpan));
r_i = zeros(neuron.N,length(tSpan));
current = zeros(neuron.N,length(tSpan));
i_step = 1; 

e_spikes = 0;
i_spikes = 0;
%pattern_input = zeros(Ne,1);
%pattern_input(rand(size(pattern_input))<0.5) = 1;
%p_index = find(pattern_input);
%p_index=[10:10:80];
p_index=[3:3:72];


avg_firing_rate = zeros(N,length(tSpan));
firing_time = false(N,length(tSpan));
delta_w = zeros(1,length(tSpan));
dt1 = [];
dt2 = [];
wij = zeros(5,length(tSpan));
wi = zeros(length(p_index),length(tSpan));
cv_w = zeros(1,length(tSpan));
mf_S = 15*1e-9;
ext_ee_S = 16*1e-9;
ext_ei_S = 15*1e-9;
ext_ie_S = 12*1e-9;
ext_ii_S = 15*1e-9;
rate_ee = 0.8;
rate_ei = 0.8;
rate_ie = 0.8;
rate_ii = 0.8;
% scale factor
g_ampa_ee = rand(neuron.Ne,1)*1e-3;
g_gaba_ei = rand(neuron.Ne,1)*1e-3;
g_ampa_ie = rand(neuron.Ni,1)*1e-3;
g_gaba_ii = rand(neuron.Ni,1)*1e-3;

g_ampa_mf = rand(neuron.Ne,1)*1e-3;
f_ext = 0.5;

ll=1.0;
mm=1.0;

for t = tSpan   %tSpan = 0:0.001:20.9

    I_mf = zeros(neuron.Ne,1); 
    ref_ind = find((t - last_spike)<refrac_t & (last_spike>0));
    v(ref_ind) = neuron.Vreset(ref_ind);%refractoriness 

    vt_neuron(:,i_step) = v;
    fired=find(v>=neuron.Vcut); % indices of cells spikes 
    v(fired) = neuron.Vreset(fired);
    u(fired) = u(fired)+ neuron.b(fired);
    last_spike(fired) = t;
    firings=[firings; t+0*fired,fired];
    ind_e_fired = intersect(fired, neuron.index_e);
    ind_i_fired = intersect(fired, neuron.index_i);
    firing_time(fired,i_step) = true;   

% apply the stimulus

    if i_step>=1 && i_step<=21000

        I_mf(p_index) = random('poisson',1,size(p_index));
        
    end  



        current_spike = zeros(Ne,1);
        current_spike(ind_e_fired) = t; 
        if i_step>=learning_start
            w = neuron.S(neuron.index_e,neuron.index_e);
            [tmp1]=STDP(w,last_spike(neuron.index_e),current_spike,neuron);
            neuron.S(neuron.index_e,neuron.index_e) = w+tmp1;
        end

    e_spikes = e_spikes + length(ind_e_fired); 
    i_spikes = i_spikes + length(ind_i_fired);

    wi(:,i_step) = sum(neuron.S(p_index,p_index),2);
    delta_w(i_step) = sum(wi(:,i_step));
    tmp = neuron.S(p_index,p_index);
    index = connect(p_index,p_index)==1;
   % cv_w(i_step) = std(tmp(index),0,'all')/mean(tmp(index),'all');
    wij(:,i_step) = tmp(2:2:10);
    
    spikes_e = [random('poisson',rate_ee,neuron.Ne,1),random('poisson',rate_ei,neuron.Ne,1)];
    spikes_i = [random('poisson',rate_ie,neuron.Ni,1),random('poisson',rate_ii,neuron.Ni,1)];
    ext_ee_fired = find(spikes_e(:,1)>0);
    ext_ei_fired = find(spikes_e(:,2)>0);
    ext_ie_fired = find(spikes_i(:,1)>0);
    ext_ii_fired = find(spikes_i(:,2)>0);
    spikes_e(spikes_e==0) = 1;
    spikes_i(spikes_i==0) = 1;
     
    spikes_mf = I_mf;
    mf_fired = find(spikes_mf>0);
    spikes_mf(spikes_mf==0) = 1;
    

    g_ampa_mf = g_ampa_mf - delta_t*g_ampa_mf/tau_amap;
    g_ampa_mf(mf_fired) = g_ampa_mf(mf_fired) + f_ext*(1-g_ampa_mf(mf_fired));%¸ø¶¨´Ì¼¤ÊäÈë
   
    g_ampa_ee = g_ampa_ee - delta_t*g_ampa_ee/tau_amap;
    g_ampa_ee(ext_ee_fired) = g_ampa_ee(ext_ee_fired) + f_ext*(1-g_ampa_ee(ext_ee_fired));
    
    g_ampa_ie = g_ampa_ie -  delta_t*g_ampa_ie/tau_amap;
    g_ampa_ie(ext_ie_fired) = g_ampa_ie(ext_ie_fired) + f_ext*(1-g_ampa_ie(ext_ie_fired));
  
    g_gaba_ei = g_gaba_ei - delta_t*g_gaba_ei/tau_gaba;
    g_gaba_ei(ext_ei_fired) = g_gaba_ei(ext_ei_fired) + f_ext*(1-g_gaba_ei(ext_ei_fired));
    
    g_gaba_ii = g_gaba_ii  - delta_t*g_gaba_ii/tau_gaba;
    g_gaba_ii(ext_ii_fired) = g_gaba_ii(ext_ii_fired) + f_ext*(1-g_gaba_ii(ext_ii_fired));
    

    diff_ampa = 0-v;
    diff_gaba = -75*1e-3-v;
 
    

     I_ext = [ext_ee_S.*diff_ampa(neuron.index_e).*g_ampa_ee.*spikes_e(:,1) + ext_ei_S.*diff_gaba(neuron.index_e).*g_gaba_ei.*spikes_e(:,2);...
              ext_ie_S.*diff_ampa(neuron.index_i).*g_ampa_ie.*spikes_i(:,1) + ext_ii_S.*diff_gaba(neuron.index_i).*g_gaba_ii.*spikes_i(:,2)];
     I_ext(neuron.index_e) =  I_ext(neuron.index_e)+ mf_S.*diff_ampa(neuron.index_e).*g_ampa_mf.*spikes_mf;
   
     I =  I_ext + randn(neuron.N,1)*0.1*1e-9;

     his_current(:,i_step) = I;
%    I(neuron.index_e) = I(neuron.index_e) + recurrent_e(neuron.index_e) + recurrent_i(neuron.index_e);
%    I(neuron.index_e) = I(neuron.index_e) + sum(neuron.S(neuron.index_e,fired),2);
     I(neuron.index_i) = I(neuron.index_i) + sum(neuron.S(neuron.index_i,fired),2);

     his_current(:,i_step) = I-his_current(:,i_step);
     
     v = v + delta_t*(-neuron.gL .* (v-neuron.EL) +...
          neuron.Vslop * neuron.gL .*  exp((v-neuron.VT)/neuron.Vslop) + I - u)/neuron.C;
     u = u + delta_t*(neuron.a .* (v-neuron.EL) - u) ./ neuron.tau_adapt;
     current(:,i_step) = I;
     i_step = i_step + 1;
     %vv(t)=v;
     %plot(t,vv(t))
     if mod(i_step,100) == 0
        fprintf('program execution schedule--> %.2f%%\n',i_step/length(tSpan)*100);
     end
end


f_range = firings(:,1)>=0 & firings(:,1)<=neuron.T;
f_exi = ismember(firings(:,2),neuron.index_e);
spike_time_e = firings(f_exi&f_range, 1);
spike_cell_e = firings(f_exi&f_range, 2);

f_inhib = ismember(firings(:,2),neuron.index_i);
spike_time_i = firings(f_inhib&f_range, 1);
spike_cell_i = firings(f_inhib&f_range, 2);

figure;
hold on;
plot(spike_time_e,spike_cell_e,'b.');
plot(spike_time_i,spike_cell_i,'k.');
legend('inhibitory neuron','excitatory neuron');
hold off;
disp(['excitatory population firing rate:',num2str(e_spikes/Ne/neuron.T)]);
disp(['inhibitory population firing rate:',num2str(i_spikes/Ni/neuron.T)]);
 


for i = 250:length(tSpan)-250
       avg_firing_rate(:,i) = sum(firing_time(:,i-249:i+250),2)*2.0;
end


% mean firing rate
figure;
indices = 5:100:15000;
hold on;
h1=plot(tSpan,sum(avg_firing_rate(p_index,:))/length(p_index),'r','linewidth',0.8);
h2=plot(tSpan,sum(avg_firing_rate(neuron.index_e,:))/Ne,'k','linewidth',0.8);
h3=plot(tSpan,sum(avg_firing_rate(neuron.index_i,:))/Ni,'b','linewidth',0.8);
ylabel('Firing rate of a group of neurons£¨Hz£©');
hold off;
xlabel('Time£¨s£©');
legend([h1 h2 h3],'target neuron population','excitatory neuron population','inhibitory neuron population');

