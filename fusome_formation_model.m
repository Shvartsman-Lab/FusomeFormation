
%make range of parameter values from model
alph = linspace(0.25, 0.75, 500);
beta = linspace(0.25, 0.75, 500);

%use mean and st dev values of measurements for comparing with model
measfrac2 = [0.700 0.300];
measfrac4 = [0.492 0.283 0.128 0.096];
measfrac8 = [0.327 0.255 0.142 0.095 0.057 0.040 0.040 0.045];
measfrac16 = [0.205 0.166 0.106 0.109 0.063 0.053 0.052 0.049 0.031 0.027 0.025 0.022 0.022 0.026 0.023 0.022];

sterr2 = [0.024 0.024];
sterr4 = [0.058 0.035 0.028 0.035];
sterr8 = [0.047 0.032 0.025 0.016 0.030 0.023 0.026 0.018];
sterr16 = [0.023 0.018 0.018 0.019 0.012 0.010 0.014 0.012 0.013 0.007 0.009 0.009 0.008 0.007 0.010 0.009];

%set error for calculating for each value of parmater
error = zeros(length(alpha), length(beta));

for i = 1:length(alpha)
    for j = 1:length(beta)
        
        a = alph(i);
        b = beta(j);
        
        if a > b
        
        %derived values of fusome pieces added at each division
        f0 = 1;
        f1 = (1-a)/(a-b);
        f2 = (1-b)*(4-11*b+8*b^2)/2/(a-b)/(8-31*b+42*b^2-20*b^3);
        f3 = (1-b)^2*(2-3*b)/2/(a-b)/(8-31*b+42*b^2-20*b^3);
        f4 = (1-b)^3/2/(a-b)/(8-31*b+42*b^2-20*b^3);
        
        %calculate error for each cell fraction at each division
        cellfrac2 = [(f0+b*f1)/(f0+f1) ((1-b)*f1)/(f0+f1)];
        cellfrac4 = [(f0+b*f1+b*f2)/(f0+f1+2*f2) ((1-b)*f1+b*f2)/(f0+f1+2*f2) ((1-b)*f2)/(f0+f1+2*f2) ((1-b)*f2)/(f0+f1+2*f2)];
        cellfrac8 = [(f0+b*f1+b*f2+b*f3)/(f0+f1+2*f2+4*f3) ((1-b)*f1+b*f2+b*f3)/(f0+f1+2*f2+4*f3) ((1-b)*f2+b*f3)/(f0+f1+2*f2+4*f3) ((1-b)*f2+b*f3)/(f0+f1+2*f2+4*f3) ((1-b)*f3)/(f0+f1+2*f2+4*f3) ((1-b)*f3)/(f0+f1+2*f2+4*f3) ((1-b)*f3)/(f0+f1+2*f2+4*f3) ((1-b)*f3)/(f0+f1+2*f2+4*f3)];
        cellfrac16 = [(f0+b*f1+b*f2+b*f3+b*f4)/(f0+f1+2*f2+4*f3+8*f4) ((1-b)*f1+b*f2+b*f3+b*f4)/(f0+f1+2*f2+4*f3+8*f4) ((1-b)*f2+b*f3+b*f4)/(f0+f1+2*f2+4*f3+8*f4) ((1-b)*f2+b*f3+b*f4)/(f0+f1+2*f2+4*f3+8*f4) ((1-b)*f3+b*f4)/(f0+f1+2*f2+4*f3+8*f4) ((1-b)*f3+b*f4)/(f0+f1+2*f2+4*f3+8*f4) ((1-b)*f3+b*f4)/(f0+f1+2*f2+4*f3+8*f4) ((1-b)*f3+b*f4)/(f0+f1+2*f2+4*f3+8*f4) ((1-b)*f4)/(f0+f1+2*f2+4*f3+8*f4) ((1-b)*f4)/(f0+f1+2*f2+4*f3+8*f4) ((1-b)*f4)/(f0+f1+2*f2+4*f3+8*f4) ((1-b)*f4)/(f0+f1+2*f2+4*f3+8*f4) ((1-b)*f4)/(f0+f1+2*f2+4*f3+8*f4) ((1-b)*f4)/(f0+f1+2*f2+4*f3+8*f4) ((1-b)*f4)/(f0+f1+2*f2+4*f3+8*f4) ((1-b)*f4)/(f0+f1+2*f2+4*f3+8*f4)];
        
        %calculate abs z-score for each cell fraction compared with data
        z_error2 = sum(abs(cellfrac2-measfrac2)./sterr2);
        z_error4 = sum(abs(cellfrac4-measfrac4)./sterr4);
        z_error8 = sum(abs(cellfrac8-measfrac8)./sterr8);
        z_error16 = sum(abs(cellfrac16-measfrac16)./sterr16);
        
        %take sum of z-scores to get total error
        z_error = z_error2 + z_error4 + z_error8 + z_error16;
        error(i,j) = z_error;
        else
            error(i,j) = -Inf;
        end
    end
end

%normalize the error to the maximum found value
error_norm = error./max(error(:));

%identify the minimum error and parameter value at this min
min_error = min(error_norm(error_norm > 0));
[row, col] = find(error_norm == min_error);
alpha_min = alpha(row);
beta_min = beta(col);

%plot results (heatmap of error)
figure; hold on; box on; grid off;
h = gca;
h.FontSize = 20;
surf(beta, alpha, error_norm, 'EdgeColor', 'none'); colorbar; colormap(bone); alpha 0.8;
c = colorbar;
set(c,'TickLabelInterpreter','latex')
plot3(beta_min, alpha_min, 1, 'rx','MarkerSize',14);
ylabel('2-Cell Volume Assymetry, $\alpha$','interpreter','latex','FontSize',20)
xlabel('Fusome Division Asymmetry, $\beta$','interpreter','latex','FontSize',20)
axis square;

axis([min(beta) max(beta) min(alpha) max(alpha) 0 1])

