function p = get_params()
% get_params  外环/内环多速率 + 纸面约束参数
%
% 约定:
%   外环周期 Ts_outer = 0.10 s  (MPC)
%   内环周期 Ts_inner = 0.05 s  (姿态/抗扰)
%
% 外环模型参数（用于 dynamics_SRB 和 linearize_SRB）

%% 物理参数
p.mass = 0.48;          % kg
p.g    = 9.81;          % m/s^2

% 一阶线性阻尼（N·s/m）—— 若不确定，宁可偏大一些，外环更保守
p.kx = 0.25;
p.ky = 0.25;
p.kz = 0.25;

%% 多速率采样
p.Ts_inner = 0.05;      % s
p.Ts_outer = 0.10;      % s

%% 外环输入与约束（纸面）
% 左右对合力上限（推力圆半径）: sqrt(v^2 + xi^2) <= Tf_max
p.Tf_max = 12.0;        % N

% 前后对合力上限: 0 <= Sfb <= Ts_max
p.Ts_max = 12.0;        % N

% 倾转角上限 |a| <= a_max  (用 |v| <= xi*tan(a_max), xi>=0 实现)
p.a_max  = deg2rad(90); % rad

% 变化率约束（外环每步）
% |u_k - u_{k-1}| <= u_dot_max * Ts_outer
p.v_dot_max   = 60.0;   % (N/s)  对 v
p.xi_dot_max  = 60.0;   % (N/s)  对 xi
p.Sfb_dot_max = 60.0;   % (N/s)  对 Sfb

% 航向约束（外环输入 r）
p.r_max      = 3.0;     % rad/s
p.r_dot_max  = 30.0;    % rad/s^2  可选（若做更平滑的偏航）

%% 约束收紧（给内环抗扰/分配留余量）
% 建议初期先设为 0，跑通后再逐渐调大
p.margin_Tf = 0.0;      % N
p.margin_Ts = 0.0;      % N

%% 兼容旧命名（如你的工程里其他地方还在用）
p.Slr_max   = p.Tf_max;
p.Sfb_max   = p.Ts_max;
p.alpha_max = p.a_max;
p.Ts        = p.Ts_outer;
end
