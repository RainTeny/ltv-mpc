function [X, U, X1, X2, U1, U2] = get_rnd_trajectories_outer_pe( ...
    X0, n_traj, t_traj, show_plot, flag_mode)
% get_rnd_trajectories_outer_pe
% 用于双线性 Koopman 训练的外环随机轨迹数据（PE导向）
%
% 状态: x = [x;y;z;vx;vy;vz;psi]  (7)
% 输入: u = [v; xi; Sfb; r]        (4)
%
% 核心：多级PRBS + 随机驻留时间 + 变化率限制 + 约束投影
%      竖直通道围绕 mg 做零均值激励，避免z发散

P = get_params();

%% ===== 固定步长 Ts = 0.1 s =====
Ts = 0.1;

% 时间轴
if nargin < 3 || isempty(t_traj)
    t_traj = 0:Ts:25;
elseif isscalar(t_traj)
    t_traj = 0:Ts:t_traj;
else
    t_traj = t_traj(:).';
end
if nargin < 4 || isempty(show_plot), show_plot = false; end
if nargin < 5 || isempty(flag_mode), flag_mode = 'train'; end

t_traj = t_traj(:)';
T      = numel(t_traj);
Tm1    = T-1;

%% 参数与约束
m  = P.mass;
g  = P.g;

Tf_max = P.Tf_max - getfieldwithdef(P,'margin_Tf',0.0);   % 推力圆半径 on (v,xi)
Ts_max = P.Ts_max - getfieldwithdef(P,'margin_Ts',0.0);   % Sfb in [0,Ts_max]

v_dot_max   = P.v_dot_max;
xi_dot_max  = P.xi_dot_max;
Sfb_dot_max = P.Sfb_dot_max;

r_max      = P.r_max;
r_dot_max  = getfieldwithdef(P,'r_dot_max',30.0);

pos_lim = getfieldwithdef(P,'pos_lim',200);
vel_lim = getfieldwithdef(P,'vel_lim',10);

% 训练/验证/上线：激励幅度与驻留时间范围（越小越“像MPC采集”）
switch lower(flag_mode)
    case 'mpc'
        amp_u   = 0.45;         % 输入幅度比例
        amp_fz  = 0.20;         % 竖直激励相对 mg
        hold_rng = [6, 14];     % 每段驻留(步) 0.6~1.4s
        imp_prob = 0.01;        % 脉冲概率
    case 'val'
        amp_u   = 0.60;
        amp_fz  = 0.30;
        hold_rng = [5, 12];
        imp_prob = 0.015;
    otherwise % train
        amp_u   = 0.80;
        amp_fz  = 0.45;
        hold_rng = [4, 10];     % 0.4~1.0s
        imp_prob = 0.02;
end

%% 初始状态模板
x0 = zeros(7,1);
if nargin >= 1 && ~isempty(X0)
    x0(1:min(7,numel(X0))) = X0(1:min(7,numel(X0)));
end

%% 输出累计
X = []; U = []; X1 = []; X2 = []; U1 = []; U2 = [];
ok_cnt = 0; bad_cnt = 0;
bad = struct('finite',0,'pos',0,'vel',0,'sat',0);

if show_plot
    figure(92); clf; hold on; grid on; axis equal;
    xlabel('x'); ylabel('y'); zlabel('z');
    title('Outer-loop PE trajectories (Ts=0.1s)');
end

for epi = 1:n_traj

    % --- 初值：小扰动 + 随机psi（避免退化） ---
    xk = x0;
    xk(7) = (2*rand-1)*deg2rad(15);
    xk(4:6) = clip(0.5*(2*rand(3,1)-1), -vel_lim, vel_lim);

    x_rec = zeros(T,7);
    u_rec = zeros(4,Tm1);
    x_rec(1,:) = xk.';

    % ===== 生成“目标输入序列”(未做rate/约束) =====
    % 多级PRBS：levels = {-1,0,1}，不同通道独立驻留
    v_tgt  = amp_u * Tf_max * gen_multilevel_prbs(Tm1, hold_rng, [-1 0 1]);
    xi_tgt = amp_u * Tf_max * gen_multilevel_prbs(Tm1, hold_rng, [-1 0 1]);
    r_tgt  = amp_u * r_max  * gen_multilevel_prbs(Tm1, hold_rng, [-1 0 1]);

    % 竖直总力目标：Fz = mg + delta(t)，delta零均值
    delta  = amp_fz * (m*g) * gen_multilevel_prbs(Tm1, hold_rng, [-1 0 1]);
    Fz_tgt = m*g + delta;

    % 偶发“脉冲/阶跃”注入：增加频谱覆盖（受限幅约束）
    for k = 1:Tm1
        if rand < imp_prob
            v_tgt(k)  = sign(2*rand-1) * amp_u * Tf_max;
        end
        if rand < imp_prob
            xi_tgt(k) = sign(2*rand-1) * amp_u * Tf_max;
        end
        if rand < imp_prob
            r_tgt(k)  = sign(2*rand-1) * amp_u * r_max;
        end
        if rand < imp_prob
            Fz_tgt(k) = m*g + sign(2*rand-1)*amp_fz*(m*g);
        end
    end

    % ===== 在线施加变化率 + 约束投影，并推进动力学 =====
    u_prev = zeros(4,1);
    sat_ratio_acc = 0;

    ok = true;
    for k = 1:Tm1

        % 1) rate limit
        v  = limit_rate(v_tgt(k),  u_prev(1), v_dot_max,   Ts);
        xi = limit_rate(xi_tgt(k), u_prev(2), xi_dot_max,  Ts);
        r  = limit_rate(r_tgt(k),  u_prev(4), r_dot_max,   Ts);

        % 2) 推力圆投影：sqrt(v^2+xi^2) <= Tf_max
        [v, xi] = project_to_disc(v, xi, Tf_max);

        % 3) 竖直分配：Sfb = Fz - xi，然后限幅+rate
        Sfb_raw = Fz_tgt(k) - xi;
        Sfb_raw = clip(Sfb_raw, 0, Ts_max);
        Sfb     = limit_rate(Sfb_raw, u_prev(3), Sfb_dot_max, Ts);

        % 4) 最终限幅
        r   = clip(r, -r_max, r_max);
        Sfb = clip(Sfb, 0, Ts_max);

        % 饱和率统计（用于剔除“几乎全饱和”的坏数据）
        sat_ratio_acc = sat_ratio_acc + ...
            (abs(hypot(v,xi) - Tf_max) < 1e-6) + ...
            (abs(Sfb - 0) < 1e-6) + (abs(Sfb - Ts_max) < 1e-6);

        u = [v; xi; Sfb; r];
        u_rec(:,k) = u;

        % 5) RK4 积分外环动力学
        xk = rk4(@dynamics_SRB, xk, u, P, Ts);

        % 6) 检查
        if any(~isfinite(xk))
            bad.finite = bad.finite + 1; ok = false; break;
        end
        if max(abs(xk(1:3))) > pos_lim
            bad.pos = bad.pos + 1; ok = false; break;
        end
        if max(abs(xk(4:6))) > vel_lim
            bad.vel = bad.vel + 1; ok = false; break;
        end

        % 安全剪裁
        xk(1:3) = clip(xk(1:3), -pos_lim, pos_lim);
        xk(4:6) = clip(xk(4:6), -vel_lim, vel_lim);

        x_rec(k+1,:) = xk.';
        u_prev = u;
    end

    if ~ok
        bad_cnt = bad_cnt + 1;
        continue;
    end

    % 若长期饱和（数据对识别没帮助），丢弃
    sat_ratio = sat_ratio_acc / max(1, 3*Tm1);
    if sat_ratio > 0.35
        bad.sat = bad.sat + 1;
        bad_cnt = bad_cnt + 1;
        continue;
    end

    ok_cnt = ok_cnt + 1;

    % 输出拼装（与原管线一致）
    X  = [X,  x_rec.'];                       %#ok<AGROW>
    U  = [U,  [u_rec, u_rec(:,end)]];         %#ok<AGROW>
    X1 = [X1, x_rec(1:end-1,:).'];            %#ok<AGROW>
    X2 = [X2, x_rec(2:end,:).'];              %#ok<AGROW>
    U1 = [U1, u_rec];                         %#ok<AGROW>
    U2 = [U2, [u_rec(:,2:end), u_rec(:,end)]];%#ok<AGROW>

    if show_plot
        plot3(x_rec(:,1), x_rec(:,2), x_rec(:,3), 'LineWidth', 1.0);
        drawnow limitrate;
    end
end

fprintf('\n== Outer-loop PE trajectories (Ts=%.2fs) ==\n', Ts);
fprintf('  成功: %d / %d (%.1f%%)\n', ok_cnt, n_traj, 100*ok_cnt/max(1,n_traj));
fprintf('  失败: %d  [finite=%d, pos=%d, vel=%d, sat=%d]\n', ...
    bad_cnt, bad.finite, bad.pos, bad.vel, bad.sat);

end

%% ===== 工具函数 =====
function s = gen_multilevel_prbs(N, hold_rng, levels)
% 多级“PRBS风格”序列：随机选 level，保持随机步数
    s = zeros(N,1);
    i = 1;
    while i <= N
        hold = randi([hold_rng(1), hold_rng(2)]);
        val  = levels(randi(numel(levels)));
        j = min(N, i+hold-1);
        s(i:j) = val;
        i = j + 1;
    end
    % 去均值（避免偏置）
    s = s - mean(s);
    % 归一化到 [-1,1]（防止levels不对称时幅值漂）
    mx = max(abs(s));
    if mx > 1e-9
        s = s / mx;
    end
end

function [v,xi] = project_to_disc(v, xi, Rmax)
    n = hypot(v,xi);
    if n > Rmax && n > 1e-12
        s = Rmax/n;
        v = v*s; xi = xi*s;
    end
end

function y = clip(x, lo, hi)
    y = min(max(x, lo), hi);
end

function y = limit_rate(v, pv, rmax, dt)
    y = pv + clip(v - pv, -rmax*dt, rmax*dt);
end

function x1 = rk4(f, x0, u, P, h)
    k1 = f(0, x0, u, P);
    k2 = f(0, x0 + 0.5*h*k1, u, P);
    k3 = f(0, x0 + 0.5*h*k2, u, P);
    k4 = f(0, x0 + h*k3, u, P);
    x1 = x0 + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function val = getfieldwithdef(s, name, default)
    if isstruct(s) && isfield(s, name)
        val = s.(name);
    else
        val = default;
    end
end
