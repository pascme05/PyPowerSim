%% Derive the coefficients for the transfer function of the simplified equivalent circuit with a RLe load at the output
% Leveraging symbolic math toolbox

%% Init
clear all
syms s
syms r1 r2
syms Lm1 Lm2
syms Ll1 Ll2
syms M
syms w
syms Rc1 Rc2
syms Is
syms RL
syms LL

%% Transfer function for transformer with RL+e load

%% Iw1 (from (14))
% Iw1 = G11 * Iw2 + G12 *E
N_G11 = Rc2*(-RL-s*LL) - r2*(Rc2+RL+s*LL) - s*Ll2*(Rc2+RL+s*LL) - s*Lm2*(Rc2+RL+s*LL);
D_G11 = s*M*(Rc2+RL+s*LL);

N_G12 = Rc2;
D_G12 = D_G11;

%% Iw2: (from rewriting (2) and plugging in the previous)
% Iw2 = G21*V1 + G22*E
N_G21 = D_G11;
D_G21 = (r1+s*Ll1+s*Lm1)*N_G11 + s*M*D_G11;

N_G22 = N_G21 * (-1)*(r1+s*Ll1+s*Lm1) * Rc2;
D_G22 = D_G21 * (s*M*(Rc2+RL+s*LL));

%% Coefficients
N_G11 = coeffs(N_G11,s, 'All')
D_G11 = coeffs(D_G11,s, 'All')
N_G12 = coeffs(N_G12,s, 'All')
D_G12 = coeffs(D_G12,s, 'All')

N_G21 = coeffs(N_G21,s, 'All')
D_G21 = coeffs(D_G21,s, 'All')
N_G22 = coeffs(N_G22,s, 'All')
D_G22 = coeffs(D_G22,s, 'All')
