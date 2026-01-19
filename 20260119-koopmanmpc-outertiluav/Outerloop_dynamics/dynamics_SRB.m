function dXdt = dynamics_SRB(~, X, U, p)
% dynamics_SRB  外环连续时间动力学（按纸面外环 + 航向状态化）
%
% 状态: X = [x; y; z; vx; vy; vz; psi]   (7x1)
% 输入: U = [v; xi; Sfb; r]              (4x1)
%
% 说明:
%   v   : 水平合力分量（沿偏航方向）
%   xi  : 左右对合力的竖直分量（等效）
%   Sfb : 前后对合力（竖直）
%   r   : 偏航角速度（外环输入），psi_dot = r
%
% 模型（纸面）:
%   x_ddot = (cos(psi)/m)*v - (kx/m)*vx
%   y_ddot = (sin(psi)/m)*v - (ky/m)*vy
%   z_ddot = (xi+Sfb)/m - g - (kz/m)*vz

m  = p.mass;
g  = p.g;
kx = p.kx;  ky = p.ky;  kz = p.kz;

% unpack state
vx  = X(4);
vy  = X(5);
vz  = X(6);
psi = X(7);

% unpack input
v   = U(1);
xi  = U(2);
Sfb = U(3);
r   = U(4);

cpsi = cos(psi);
spsi = sin(psi);

% kinematics
xdot = vx;
ydot = vy;
zdot = vz;

% translational dynamics
vxdot = (cpsi/m)*v - (kx/m)*vx;
vydot = (spsi/m)*v - (ky/m)*vy;
vzdot = (1/m)*(xi + Sfb) - g - (kz/m)*vz;

% yaw
psidot = r;

dXdt = [xdot; ydot; zdot; vxdot; vydot; vzdot; psidot];
end
