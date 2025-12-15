function [K_TEh, K_TMh, F_out, theta_deg_air, eta_out, eta_wg_sub, eta_sp, eta_abs, eta_emitter_eff, eta_out_h, eta_wg_sub_h, eta_sp_h, eta_abs_h, eta_emitter_eff_h, eta_out_TMv, eta_wg_sub_TMv, eta_sp_TMv, eta_abs_TMv, eta_emitter_eff_TMv, u_crit, K_h, uK_h, K_h_prime, uK_h_prime, K_TMv, uK_TMv, K_TMv_prime, uK_TMv_prime, K, uK, K_prime, uK_prime, Purcell_h, Purcell_TMv, Purcell] = Func_CPS_clean(thickness_upper, thickness_below, z_plus, z_minus, n, lambda0, eta_emitter, du, u_upper)
k0 = 2*pi/lambda0;
k(1,:) = k0*n(1,:);

thickness_upper_prime = [thickness_upper 0]; %thickness_init = 0
thickness_below_prime = [0 thickness_below]; %thickness_init = 0
thickness_EML = z_plus + z_minus;
thickness = [thickness_upper thickness_EML thickness_below];

n_upper = n(1:length(thickness_upper));
n_below = n(length(thickness_upper)+2:length(thickness));
n_upper_prime = n(1:(length(thickness_upper)+1)); %n_init = n_init
n_below_prime = n((length(thickness_upper)+1):length(thickness)); % Func_TMM : 나가는 매질이 1이 아닐때도 되는가 확인

floor_EML = length(thickness_upper)+1;

%% Normalized in-plane wave vector
u = 0 : du : u_upper;
u_crit = n(1)/n(floor_EML);
du_count_below1 = [1:(1/du)];%ex 1~10000
du_count_upper1 = [(1/du)+2:(u_upper/du)];%ex 10002~20000

%% Plus
count = 0;
for u_temp = u
    count = count + 1;
    %%[r_fin_TE,R_fin_TE,t_fin_TE,T_fin_TE,r_fin_TM,R_fin_TM,t_fin_TM,T_fin_TM] = Func_TMM(thickness, n, lambda, u)
    [r_TE_plus_temp,R_TE_plus_temp,t_TE_plus_temp,T_TE_plus_temp,r_TM_plus_temp,R_TM_plus_temp,t_TM_plus_temp,T_TM_plus_temp] = Func_TMM(flip(thickness_upper_prime), flip(n_upper_prime), lambda0, u_temp);
    r_TE_plus(count) = r_TE_plus_temp;
    R_TE_plus(count) = R_TE_plus_temp;
    t_TE_plus(count) = t_TE_plus_temp;
    T_TE_plus(count) = T_TE_plus_temp;
    
    r_TM_plus(count) = r_TM_plus_temp;
    R_TM_plus(count) = R_TM_plus_temp;
    t_TM_plus(count) = t_TM_plus_temp;
    T_TM_plus(count) = T_TM_plus_temp;
    
end

%% Minus
count = 0;
for u_temp = u
    count = count + 1;
    %%[r_fin_TE,R_fin_TE,t_fin_TE,T_fin_TE,r_fin_TM,R_fin_TM,t_fin_TM,T_fin_TM] = Func_TMM(thickness, n, lambda, u)
    [r_TE_minus_temp,R_TE_minus_temp,t_TE_minus_temp,T_TE_minus_temp,r_TM_minus_temp,R_TM_minus_temp,t_TM_minus_temp,T_TM_minus_temp] = Func_TMM(thickness_below_prime, n_below_prime, lambda0, u_temp);
    r_TE_minus(count) = r_TE_minus_temp;
    R_TE_minus(count) = R_TE_minus_temp;
    t_TE_minus(count) = t_TE_minus_temp;
    T_TE_minus(count) = T_TE_minus_temp;
    
    r_TM_minus(count) = r_TM_minus_temp;
    R_TM_minus(count) = R_TM_minus_temp;
    t_TM_minus(count) = t_TM_minus_temp;
    T_TM_minus(count) = T_TM_minus_temp;
    
end

%% a_TE+-, a_TM+-
a_TE_plus = r_TE_plus.*exp(2*i*k0*n(floor_EML)*sqrt((1-u.^2))*z_plus);
a_TE_minus = r_TE_minus.*exp(2*i*k0*n(floor_EML)*sqrt((1-u.^2))*z_minus);
a_TE = a_TE_plus.*a_TE_minus;

a_TM_plus = r_TM_plus.*exp(2*i*k0*n(floor_EML)*sqrt((1-u.^2))*z_plus);
a_TM_minus = r_TM_minus.*exp(2*i*k0*n(floor_EML)*sqrt((1-u.^2))*z_minus);
a_TM = a_TM_plus.*a_TM_minus;

%% K_TEh and uK_TEh
K_TEh = (3/8)*real((1./sqrt(1-u.^2)).* ((1+a_TE_plus).*(1+a_TE_minus))./(1-a_TE));
uK_TEh = u.*K_TEh;

%% K_TMh, uK_TMh, K_TMv, uK_TMv, K, and uK
K_TMh = (3/8)*real(sqrt(1-u.^2).*((1-a_TM_plus).*(1-a_TM_minus))./(1-a_TM));
uK_TMh = u.*K_TMh;

K_TMv = (3/4)*real((u.^2./sqrt(1-u.^2)).*((1+a_TM_plus).*(1+a_TM_minus))./(1-a_TM));
uK_TMv = u.*K_TMv;
uK_v = uK_TMv;%%

K = (1/3)*(K_TMv+2*K_TMh+2*K_TEh);
uK = u.*K;

K_h = K_TMh+K_TEh;
uK_h = u.*(K_TMh+K_TEh);

%% K_prime
K_TEh_prime = (3/16) * 1./sqrt(1-u.^2) .* (abs(1+a_TE_minus)).^2 ./ (abs(1-a_TE)).^2 .* T_TE_plus;
K_TMh_prime = (3/16) * sqrt(1-u.^2) .* (abs(1-a_TM_minus)).^2 ./ (abs(1-a_TM)).^2 .* T_TM_plus;
K_TMv_prime = (3/8) * u.^2 ./ sqrt(1-u.^2) .* (abs(1+a_TM_minus)).^2 ./ (abs(1-a_TM)).^2 .* T_TM_plus;
K_prime = (1/3)*(K_TMv_prime+2*K_TMh_prime+2*K_TEh_prime);
K_h_prime = (K_TMh_prime+K_TEh_prime);

uK_TEh_prime = u.*K_TEh_prime;
uK_TMh_prime = u.*K_TMh_prime;
uK_TMv_prime = u.*K_TMv_prime;
uK_prime = u.*K_prime;
uK_h_prime = u.*K_h_prime;

%% Value of Integration
F_TEh_out = 2*sum((uK_TEh(1:floor(u_crit*1/du)))*du, "omitnan");
F_TMh_out = 2*sum((uK_TMh(1:floor(u_crit*1/du)))*du, "omitnan");
F_TMv_out = 2*sum((uK_TMv(1:floor(u_crit*1/du)))*du, "omitnan");

F_TEh_prime_out = 2*sum((uK_TEh_prime(1:floor(u_crit*1/du)))*du, "omitnan");
F_TMh_prime_out = 2*sum((uK_TMh_prime(1:floor(u_crit*1/du)))*du, "omitnan");
F_TMv_prime_out = 2*sum((uK_TMv_prime(1:floor(u_crit*1/du)))*du, "omitnan");

%% Total(vertical+horizontal)
F_out_direction = 2*sum((uK(1:floor(u_crit*1/du)))*du, "omitnan");%
F_wg_sub = 2*sum((uK(floor(u_crit*1/du):1/du))*du, "omitnan");
F_sp = 2*sum((uK(1/du + 1:u_upper*1/du))*du, "omitnan");
F_out = 2*sum(uK_prime(1:floor(u_crit*1/du))*du, "omitnan");%
F_abs = F_out_direction - F_out;

eta_out = F_out / (F_out+F_wg_sub+F_sp+F_abs);
eta_wg_sub = F_wg_sub / (F_out+F_wg_sub+F_sp+F_abs);
eta_sp = F_sp / (F_out+F_wg_sub+F_sp+F_abs);
eta_abs = F_abs / (F_out+F_wg_sub+F_sp+F_abs);

Purcell = 2*sum(uK*du, "omitnan");
eta_emitter_eff = eta_emitter*Purcell / (eta_emitter*Purcell + (1-eta_emitter));

%% Horizontal
F_out_direction_TEh = 2*sum((uK_TEh(1:floor(u_crit*1/du)))*du, "omitnan");
F_wg_sub_TEh = 2*sum((uK_TEh(floor(u_crit*1/du):1/du))*du, "omitnan");
F_sp_TEh = 2*sum((uK_TEh(1/du + 1:u_upper*1/du))*du, "omitnan");
F_out_TEh = 2*sum(uK_TEh_prime(1:floor(u_crit*1/du))*du, "omitnan");
F_abs_TEh = F_out_direction_TEh - F_out_TEh;

eta_out_TEh = F_out_TEh / (F_out_TEh+F_wg_sub_TEh+F_sp_TEh);
eta_wg_sub_TEh = F_wg_sub_TEh / (F_out_TEh+F_wg_sub_TEh+F_sp_TEh);
eta_sp_TEh = F_sp_TEh / (F_out_TEh+F_wg_sub_TEh+F_sp_TEh);
eta_abs_TEh = F_abs_TEh / (F_out_TEh+F_wg_sub_TEh+F_sp_TEh);

Purcell_TEh = 2*sum(uK_TEh*du, "omitnan");
eta_emitter_eff_TEh = eta_emitter*Purcell_TEh / (eta_emitter*Purcell_TEh + (1-eta_emitter));

%%
F_out_direction_TMh = 2*sum((uK_TMh(1:floor(u_crit*1/du)))*du, "omitnan");
F_wg_sub_TMh = 2*sum((uK_TMh(floor(u_crit*1/du):1/du))*du, "omitnan");
F_sp_TMh = 2*sum((uK_TMh(1/du + 1:u_upper*1/du))*du, "omitnan");
F_out_TMh = 2*sum(uK_TMh_prime(1:floor(u_crit*1/du))*du, "omitnan");
F_abs_TMh = F_out_direction_TMh - F_out_TMh;

eta_out_TMh = F_out_TMh / (F_out_TMh+F_wg_sub_TMh+F_sp_TMh);
eta_wg_sub_TMh = F_wg_sub_TMh / (F_out_TMh+F_wg_sub_TMh+F_sp_TMh);
eta_sp_TMh = F_sp_TMh / (F_out_TMh+F_wg_sub_TMh+F_sp_TMh);
eta_abs_TMh = F_abs_TMh / (F_out_TMh+F_wg_sub_TMh+F_sp_TMh);

Purcell_TMh = 2*sum(uK_TMh*du, "omitnan");
eta_emitter_eff_TMh = eta_emitter*Purcell_TMh / (eta_emitter*Purcell_TMh + (1-eta_emitter));

%%
F_out_direction_h = (F_out_direction_TEh+F_out_direction_TMh);
F_wg_sub_h = (F_wg_sub_TEh+F_wg_sub_TMh);
F_sp_h = (F_sp_TEh+F_sp_TMh);
F_out_h = (F_out_TEh+F_out_TMh);
F_abs_h = (F_abs_TEh+F_abs_TMh);

eta_out_h = F_out_h / (F_out_h+F_wg_sub_h+F_sp_h+F_abs_h);
eta_wg_sub_h = F_wg_sub_h / (F_out_h+F_wg_sub_h+F_sp_h+F_abs_h);
eta_sp_h = F_sp_h / (F_out_h+F_wg_sub_h+F_sp_h+F_abs_h);
eta_abs_h = F_abs_h / (F_out_h+F_wg_sub_h+F_sp_h+F_abs_h);

Purcell_h = 2*sum(uK_h*du, "omitnan");
eta_emitter_eff_h = eta_emitter*Purcell_h / (eta_emitter*Purcell_h + (1-eta_emitter));

%% Vertical
F_out_direction_TMv = (1/3)*2*sum((uK_TMv(1:floor(u_crit*1/du)))*du, "omitnan");
F_wg_sub_TMv = (1/3)*2*sum((uK_TMv(floor(u_crit*1/du):1/du))*du, "omitnan");
F_sp_TMv = (1/3)*2*sum((uK_TMv(1/du + 1:u_upper*1/du))*du, "omitnan");
F_out_TMv = (1/3)*2*sum(uK_TMv_prime(1:floor(u_crit*1/du))*du, "omitnan");
F_abs_TMv = F_out_direction_TMv - F_out_TMv;

eta_out_TMv = F_out_TMv / (F_out_TMv+F_wg_sub_TMv+F_sp_TMv+F_abs_TMv);
eta_wg_sub_TMv = F_wg_sub_TMv / (F_out_TMv+F_wg_sub_TMv+F_sp_TMv+F_abs_TMv);
eta_sp_TMv = F_sp_TMv / (F_out_TMv+F_wg_sub_TMv+F_sp_TMv+F_abs_TMv);
eta_abs_TMv = F_abs_TMv / (F_out_TMv+F_wg_sub_TMv+F_sp_TMv+F_abs_TMv);

Purcell_TMv = 2*sum(uK_TMv*du, "omitnan");
eta_emitter_eff_TMv = eta_emitter*Purcell_TMv / (eta_emitter*Purcell_TMv + (1-eta_emitter));

theta_rad_air = asin(n(floor_EML)/n(1)*u);
theta_deg_air = theta_rad_air*180/pi;
% save('sample_uK.mat')
end