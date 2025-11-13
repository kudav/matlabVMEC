
function R = tb_zyx(a, b, g)
% Calculates Tait-Bryan z-y'-x" active rotation matrix given rotation angles `alpha`,`beta`,`gamma` in radians
% Arguments:
%     a: rotation angle about z [radians]
%     b: rotation angle about y' [radians]
%     g: rotation angle about x" [radians]
% Return Value:
%     Rotation Matrix
%Function converted from D3D FIDASIM idl routines

sa = sin(a); ca = cos(a);
sb = sin(b); cb = cos(b);
sg = sin(g); cg = cos(g);

R = zeros(3, 3);
R(1, 1) = ca * cb; R(1, 2) = ca * sb * sg - cg * sa; R(1, 3) = sa * sg + ca * cg * sb;
R(2, 1) = cb * sa; R(2, 2) = ca * cg + sa * sb * sg; R(2, 3) = cg * sa * sb - ca * sg;
R(3, 1) = -sb;     R(3, 2) = cb * sg;              R(3, 3) = cb * cg;

% If you prefer returning a transposed matrix
% R = R';

% If you want to convert the result to single precision (float)
% R = single(R);
end
