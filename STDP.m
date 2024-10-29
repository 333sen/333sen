%Before running, you need to change the file name to STDP because it is called by the S1A and S1B files
function [delta_w] = STDP(w,last_spike,current_spike,net)
% % STDP learning rule;
delta_t = 1e-3;
tau_p = 15*delta_t;
tau_p1 = 30*delta_t;
tau_n = 30*delta_t;
alpha = 0.8;
A = 1;
u = 0.3;
LR = 0.5;
w0 = 3*net.Jstrength*net.Je;
% size = net.Ne - net.Nb;
% range = net.Nb+1:net.Ne;
Size = net.Ne ;
% range = 1:net.Ne;
tmp = find(~eye(Size,Size)& (w~=0));
delta_w = zeros(size(w));

  last_spike = last_spike';
  cs = current_spike(:,ones(1,Size));
  ls = last_spike(ones(Size,1),:);
  is_valid1 = intersect(find(cs ~= 0),find(ls ~= 0)); 
%   is_zero = find(cs==ls);
  is_valid1 = intersect(is_valid1,tmp);
%   is_valid1 = setdiff(is_valid1,is_zero); 
  d_t = cs - ls;
  pos = find(d_t>0);
  is_valid1 = intersect(is_valid1,pos);
 

  delta_w(is_valid1) =  LR*w0^(1-u)*(w(is_valid1).^u).*(A*exp(-d_t(is_valid1)/tau_p)-(A-1)*exp(-d_t(is_valid1)/tau_p1));

  current_spike = current_spike'; 
  last_spike = last_spike';
  ls = current_spike(ones(Size,1),:);
  cs = last_spike(:,ones(1,Size));
  is_valid2 = intersect(find(cs ~= 0),find(ls ~= 0));
  d_t = cs - ls;
  neg = find(d_t<0);
  is_valid2 = intersect(is_valid2,tmp); 
  is_valid2 = intersect(is_valid2,neg); 
  delta_w(is_valid2) =  - LR*alpha*w(is_valid2).*exp(d_t(is_valid2)/tau_n);


 tmp = delta_w + w ;
 Neg = (tmp<0);
 delta_w(Neg) =  0;
 


 if isempty(find(delta_w>1)) == 0
     disp('STDP bug!!!!!!!!!')
 end
end
