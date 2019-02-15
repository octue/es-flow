% Script to get results to compare between matlab and c++ implementations of
% adem.

format long

pi_coles = 0.42;
kappa = 0.41;
delta_c = 1000.0;
u_inf = 20.0;
s = 23.6;
beta = 0.0;
zeta = 0.0;

ademResult = adem(delta_c, u_inf, pi_coles, s, zeta, beta);

% Translate to data compatible with C++ AdemData class
clearvars results
results.eddy_types = 'A, B1+B2+B3+B4';
results.z = ademResult.z;
results.eta = ademResult.eta;
results.k1z = ademResult.k1z;
results.lambda_e = ademResult.lambdaE;
results.delta_c = ademResult.deltac;
results.u_inf = ademResult.U1;
results.kappa = kappa;
results.pi_coles = ademResult.Pi;
results.shear_ratio = ademResult.S;
results.u_tau = ademResult.Utau;
results.zeta = ademResult.zeta;
results.beta = ademResult.beta;
results.u_horizontal = ademResult.Ux;
results.t2wa = ademResult.T2wA;
results.t2wb = ademResult.T2wB;
results.residual_a = ademResult.residualA;
results.residual_b = ademResult.residualB;
results.reynolds_stress = ademResult.R;
results.reynolds_stress_a = ademResult.RA;
results.reynolds_stress_b = ademResult.RB;
results.psi = ademResult.Psi;
results.psi_a = ademResult.PsiA;
results.psi_b = ademResult.PsiB;

save('adem_data_check.mat', '-struct', 'results')
