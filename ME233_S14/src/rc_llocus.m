% 2014-01-31
% Xu Chen 
% maxchen@berkeley.edu 
N = 4;
q = 2;
p = 0.8;
G = tf([q],conv([1 0 0 0 -1],[1 -p]),1);
rlocus(G)