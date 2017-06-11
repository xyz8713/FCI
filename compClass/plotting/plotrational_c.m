%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%          Test rational approximations                   %%%%%%%%
%%%%%%%            Yuanzhe Xi, 02/01/2016                       %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
%% Part 3: Setup params for the circle \Gamma
%%-------- radius
r = 1;
%%-------- intial center
c = 0.0;

%% Part 4: Compute the weights and nodes for rational approximation to h
%%-------- number of poles on the upper half plane
nC  = 8;
%%-------- Get nodes (z) and weights (om) on a unit circle
[z, om] = contQuad(nC,2);
%%-------- Get rational approximation to 1/z outside circle
%% 1/z \approx \sum omega(k)/(sigma(k)+c\sigma(k)^2) *(1/(1/sigma(k)+c-z))
%% 1/z \approx \sum coeff(k)* 1/(shift - z)
sigma = z/r;
omega = om/r;
%c = 0-max(real(z))*r
coefs = omega./(sigma + c*sigma.*sigma);
shift = 1./sigma +c;
coefs = [coefs conj(coefs)];
shift = [shift conj(shift)];
%%-------- Plot the approximation of 1/z outside the circle
figure(1)
comph_c(coefs,shift,r,c);
%%-------- Get rational approximation to 1 inside circle
%%-------- 1_[c-r,c+r] = \sum coef2(k)*1/(shift(k) -z)
sigma2 = z*r;
omega2 = om*r;
coefs2 = -conj(omega2);
coefs2 = [coefs2 conj(coefs2)];
shift2 = shift;
figure(2)
compp_c(coefs2, shift2, r, c);
%%------- Plot the nodes on the circle
figure(3)
plotnodes_c(shift,r,c);
figure(4)
compp2_c(coefs2, shift2, r, c);



