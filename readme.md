# [CRB Boundary] MUSIC, ESPRIT, and CRB with CLAMEROS LAO [with MATLAB Code]

References: MUSIC_maximum_likelihood_and_Cramer-Rao_bound_further_results_and_comparisons.  
CRB main derivation process:
MATLABCode%% \nclear\nclc\naddpath('.\func')\nrng(666);\n%% Parameter settings\nN               =           10;                % Antenna at base station\nfc              =           100e6;              %100GHz\nlambda_c        =           3e8/fc;\nd               =           lambda_c/2;\nT               =           100;\nNit             =           1e2;\nSNRdBs          =           -10:5:20;\n%% \ntheta   = 0.5;%pi/6\nH   = genSteerVector(theta, N, d, lambda_c);\nY   = zeros(N, T);\nA   = eye(N);\na1  = genPartialSteerVector(theta, N, d, lambda_c, 1);\na2  = genPartialSteerVector(theta, N, d, lambda_c, 2);\nMseMUSIC    = zeros(length(SNRdBs), 1);\nMseESPRIT   = zeros(length(SNRdBs), 1);\nCRLB        = zeros(length(SNRdBs), 1);\nD           = a1;\nCst         = D'*(eye(N)-H*inv(H'*H)*H')*D;\nfor i_SNR = 1:length(SNRdBs)\n    fprintf('==============================================\n');\n    fprintf('SNR        = %ddB\n', SNRdBs(i_SNR));\n    for i_it = 1:Nit\n        \n        if mod(i_it,10) == 0\n            fprintf('i_it/Nit   = %3d/%3d\n', i_it, Nit);\n        end\n        X = (sqrt(2)/2*randn(1,T)+1j*sqrt(2)/2*randn(1,T));\n        for i_T = 1:T\n            Y(:, i_T) = awgn(H*X(i_T), SNRdBs(i_SNR), 'measured');\n        end\n        \n        theta_MUISC     = MUSIC(Y);\n        psi = TLS_ESPRIT_Algorithm(Y, 1);\n        theta__ESPRIT = log(psi)/(1j*pi); \n        MseMUSIC(i_SNR)     = MseMUSIC(i_SNR) + abs(asin(theta_MUISC) - asin(theta))^2;\n        MseESPRIT(i_SNR)    = MseESPRIT(i_SNR) + abs(asin(theta__ESPRIT) - asin(theta))^2;\n        \n        sigma2  = 10^(-SNRdBs(i_SNR)/10);\n% This sigma2 is progressively consistent with the above value\n%         sigma2  = (norm(Y, 'fro')^2-norm(H*X, 'fro')^2)/...\n%                 (size(Y, 1)*size(Y, 2));\n%         X_bar   = kron(X.', eye(N));\n%         y       = reshape(Y,[],1);\n%         CRLB(i_SNR)     = CRLB(i_SNR) + sigma2/2./...\n%                             real(-y'*X_bar*a2 + a2'*(X_bar'*X_bar)*H + a1'*(X_bar'*X_bar)*a1);\n        CRLB(i_SNR)     = CRLB(i_SNR) + sigma2/2/real((Cst*(X*X')));\n    end\nend\nMseMUSIC     = MseMUSIC/Nit;\nMseESPRIT    = MseESPRIT/Nit;\nCRLB         = CRLB/Nit;\n​\n% fprintf("Calculated SNRdB = %f",\n", 10*log10((norm(H*X, 'fro')^2)/(norm(Y, 'fro')^2-norm(H*X, 'fro')^2)))\n%% Drawing\nset(0,'defaultfigurecolor','w') \nfigure; hold on; grid on; box on;\nxlabel('SNR/dB');\nylabel('MSE(rad^2)');\nset(gca, 'YScale', 'log');\n​\nplot(SNRdBs, MseESPRIT,'r:s', 'LineWidth', 2);\nplot(SNRdBs, MseMUSIC,'g:o', 'LineWidth', 2);\nplot(SNRdBs, CRLB,'m:<', 'LineWidth', 2);\n​\nl = legend({  \n        'ESPRIT',...    \n        'MUSIC',...\n        'CRLB'},...    \n        'Interpreter','latex', 'Box','off'); \nl.FontSize = 12;\n​\nMATLABSimulation results  

## Images

![[CRB Boundary] MUSIC, ESPRIT, and CRB with CLAMEROS [with MATLAB code attached]](01.jpg)
![[CRB Boundary] MUSIC, ESPRIT, and CRB with CLAMEROS [with MATLAB code attached]](02.jpg)
![[CRB Boundary] MUSIC, ESPRIT, and CRB with CLAMEROS [with MATLAB code attached]](03.jpg)
![[CRB Boundary] MUSIC, ESPRIT, and CRB with CLAMEROS [with MATLAB code attached]](04.jpg)
![[CRB Boundary] MUSIC, ESPRIT, and CRB with CLAMEROS [with MATLAB code attached]](05.jpg)
![[CRB Boundary] MUSIC, ESPRIT, and CRB with CLAMEROS [with MATLAB code attached]](06.jpg)
![[CRB Boundary] MUSIC, ESPRIT, and CRB with CLAMEROS [with MATLAB code attached]](07.jpg)

<!-- PAYMENT_BLOCK -->

## Get this dataset

- Payment link (Stripe): [Buy now](https://buy.stripe.com/3cs8yP7sY87d0vu9AB)
- Payment QR code: ![Payment QR code](pay.png)

<!-- /PAYMENT_BLOCK -->
