function F = box_gauss_funct(X,A,B,C,D,E) % From /afs/ipp/home/s/sprd/XXX_DIAG/LIB
gam   = double(D);
width = double(E);
rl    = abs(0.5d0*width./gam);
Z     = abs((double(X)-double(C))./gam);
F     = double(B)*(0.5d0./width.*(erf(Z+rl) - erf(Z-rl)))+double(A);

% Normalization and cutoff as implemented in fplot
%F(F<1e-5) = 0;
F = F./sum(F,1);

end