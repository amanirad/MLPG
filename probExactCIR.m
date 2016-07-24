function prob = probExactCIR(theta_in, data,t, FIGS)


% Theta 
theta = (theta_in)



%prob = zeros(length(data(1:end-1)),1);
prob = 0;
for di = 2:length(data(1:end-1))
    
% Nsqnonlin case:
% Return f_i such that min_theta ( f_i^2 ) 
% evaluate at y, t+1 starting at x at time s
% a = theta2, b = theta1, sigma = theta3
if(FIGS==1)
figure(44)
y = linspace(0, 2*max(data), 100);
end
dens = cirpdf(data(di),t(di), data(di-1),t(di-1), theta(2),theta(1),theta(3));

prob = prob  -log(dens + 1e-300);
%     prob(di) = sqrt(-log( dens ));
 if(FIGS==1)  
     if(mod(di, 100)==0)
plot(y, cirpdf(y,t(di), data(di-1),t(di-1), theta(2),theta(1),theta(3)))
hold on
plot(data(di), cirpdf(data(di),t(di), data(di-1),t(di-1), theta(2),theta(1),theta(3)), 'ro')
     end
 end
end


hold off
prob = prob./length(data(1:end-1))


end