function [A, B, c] = linearize_SRB(Zbar, Ubar, p, Ts)
% linearize_SRB  外环提升模型（双线性 Koopman）离散线性化（用于SCP/SQP）
%
% 提升状态:
%   Z = [x; y; z; vx; vy; vz; s; c; 1],  s=sin(psi), c=cos(psi)
% 输入:
%   U = [v; xi; Sfb; r]
%
% 离散化: Euler
%   Z_{k+1} = f(Z_k, U_k)
% 线性化:  Z_{k+1} ≈ A Z_k + B U_k + c
%
% 用法:
%   p = get_params();
%   Ts = p.Ts_outer;
%   [A,B,c] = linearize_SRB(Zbar, Ubar, p, Ts);

if nargin < 4 || isempty(Ts)
    Ts = p.Ts_outer;
end

m  = p.mass;
g  = p.g;
kx = p.kx;  ky = p.ky;  kz = p.kz;

% unpack Zbar
vx = Zbar(4);
vy = Zbar(5);
vz = Zbar(6);
s  = Zbar(7);
cc = Zbar(8);

% unpack Ubar
v   = Ubar(1);
xi  = Ubar(2);
Sfb = Ubar(3);
r   = Ubar(4);

% ---- f(Z,U) (Euler discrete) ----
f = zeros(9,1);
f(1) = Zbar(1) + Ts*vx;
f(2) = Zbar(2) + Ts*vy;
f(3) = Zbar(3) + Ts*vz;
f(4) = vx + Ts*(-(kx/m)*vx + (1/m)*cc*v);
f(5) = vy + Ts*(-(ky/m)*vy + (1/m)*s*v);
f(6) = vz + Ts*(-(kz/m)*vz + (1/m)*(xi+Sfb) - g);
f(7) = s  + Ts*(cc*r);
f(8) = cc + Ts*(-s*r);
f(9) = 1;

% ---- Jacobians ----
A = eye(9);
B = zeros(9,4);

% rows 1..3 (positions)
A(1,4) = Ts;
A(2,5) = Ts;
A(3,6) = Ts;

% row 4: vx
A(4,4) = 1 + Ts*(-kx/m);
A(4,8) = Ts*(1/m)*v;     % d/dc of (c*v)
B(4,1) = Ts*(1/m)*cc;    % d/dv

% row 5: vy
A(5,5) = 1 + Ts*(-ky/m);
A(5,7) = Ts*(1/m)*v;     % d/ds of (s*v)
B(5,1) = Ts*(1/m)*s;     % d/dv

% row 6: vz
A(6,6) = 1 + Ts*(-kz/m);
B(6,2) = Ts*(1/m);       % d/dxi
B(6,3) = Ts*(1/m);       % d/dSfb

% row 7: s
A(7,8) = Ts*r;           % d/dc of (c*r)
B(7,4) = Ts*cc;          % d/dr

% row 8: c
A(8,7) = Ts*(-r);        % d/ds of (-s*r)
B(8,4) = Ts*(-s);        % d/dr

% affine term
c = f - A*Zbar - B*Ubar;
end
