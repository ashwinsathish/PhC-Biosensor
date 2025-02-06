%%%%%%%%%% for tm binary grating glass metal  metalgrating dielectric
%%%%%%%%%% wavelength modulation
clear;
close all;
clc;
clf;
hold all;

c=3e8;

count=1;
w=60:0.01:70;

kkk_range = 1:1:2;
iii_range = w;
total_iterations = length(iii_range);
td = 2500e-9; % Blood sample thickness
IR12_nd_133 = [];
% Initialize arrays to hold FOM and QF values

nsio2 = 1.444;
nPbS = 4.25;
nptse2 = 3.7803 - 1i*0.1694;
ntio2 = 2.4328;
nMgF2 = 1.3709;
nZnO = 2.0034;
nmos2 = 3.647;
nmxene = 2.1095 - 1.1338i;
nmeta = -3.0984 - 0.001i;

% Defining the thicknesses
t_PbS = 350e-9;
t_sio2 = 400e-9;
t_MgF2 = 250e-9; 
t_si = 40e-9;
t_ptse2 = 4e-9;
t_mos2 = 0.65e-9;
t_gr = 0.34e-9;
t_tio2 = 2e-9;
t_Al = 30e-9;
t_ZnO = 2e-9;

totalOuterIterations = length(t_si);
completedOuterIterations = 0;

N=6;
sensitivities = zeros(1, N);
FOM_values = zeros(1, N);

startTime = tic;

for stack = 1:N
    %stack=1;
    completedOuterIterations = completedOuterIterations + 1;
    figure
    title(sprintf('SPR Curve for stacks = %.2f', stack));
    hold all;
    for kkk=kkk_range
        for iii=iii_range
            lambda0=1.550e-6;
            lamdac = 2.4511e-5;
            lamdap = 1.0657e-7;
            epsilonreal = 1-((lambda0.^2.*(lamdac)^2)./(lamdap^2.*(lambda0.^2+lamdac^2)));
            epsilonim = ((lambda0.^3.*(lamdac))./(lamdap^2.*(lambda0.^2+lamdac^2)));
    
            n2 = sqrt((sqrt(epsilonreal.^2+epsilonim.^2)+epsilonreal)./2);
            k2 = sqrt((sqrt(epsilonreal.^2+epsilonim.^2)-epsilonreal)./2);
    
            Numords = 101;  % number of diffractive orders maintained
            nc = 1.426;   % CaF2 prism
            %nc = 1.3194; % NaF prism
            ns = 1;  % region 3 substrate refractive index
            period = 100e-9;  % grating period in microns
    
            if kkk==1
                nd = 1.33;
            else
                nd = 1.34;
            end
            nm = n2-1i*k2; % Aluminium refractive index
            A = 3.44904;
            A1 = 2271.88813;
            A2 = 3.39538;
            t1 = 0.058304;
            t2 = 0.30384;
    
            neth = 1.36;
            nALO = 1.746;
            teth = 0e-9;
            tALO = kkk*5e-9;
            mos2_layers = 0;
            j = sqrt(-1);
            c1 = 5.446e6;
    
            nsi = A + A1*exp(-1.55/t1) + A2*exp(-1.55/t2);
            nx = 1; % dummy variable
            ngr =  3.0 - 1.1491i;
            nTi = 2.7043 - 3.7657i;
            nAg = 0.0562 - 4.2776i;
            nBFO = 2.968;
            nBTO = 2.30;
            %nmeta = -sqrt((-4.0 + 1i*0.001)*(-2.4 + 1i*0.001));
    
            % Initialize nr and depth arrays with Aluminium and ZnO layers
            nr = [nm];  
            depth = [t_Al];
            
            % Add Si/SiO2 stacks based on the value of N
            for i = 1:stack
                nr = [nr, nsi, nptse2];  % Adding Si and SiO2 to the refractive index array
                depth = [depth, t_si, t_ptse2];  % Adding Si and SiO2 thicknesses to the depth array
            end
            
            % Add the substrate as the last element of the arrays
            nr = [nr, nd];  
            depth = [depth, td];
    
            Ngrat = length(nr);   % number of grating slices
                   
            ng = nr; % Groove
            Filfac = repmat(0.5, 1, length(nr));  % Fill factor for ridges, default value of 0.5 for all
            Disp = zeros(1, length(nr));  % Ridge displacement, default value of 0 for all
    
            theta0=iii;  %   angle of incidence
            phi0=0;    % azimuthal angle of incidence
            deg=pi/180;
            Nmax=(Numords-1)/2; % highest order number retained
            I=(-Nmax:Nmax)';  % I is the order index
            p=(Numords+1)/2;  % index of zeroth order
            theta0=theta0*deg;
            phi0=phi0*deg;  % converting in radians
            epsc=nc^2;    % relative permittivity
            epss=ns^2;
    
            k0=2*pi/lambda0;  % free space vector
            K=2*pi/period;    % grating vector
    
            kc=k0*nc;
            kx=kc*sin(theta0)*cos(phi0)-I*K;  % region1 wave vector components
    
            ky=kc*sin(theta0)*sin(phi0)*ones(size(kx));
            kzc=sqrt(kc^2-kx.^2-ky.^2);
            bad_indices=find((real(kzc)-imag(kzc))<0);  
            kzc(bad_indices)=-kzc(bad_indices);              
            ks=k0*ns;
            kzs=sqrt(ks^2-kx.^2-ky.^2);
            bad_indices=find((real(kzs)-imag(kzs))<0);  % region3 wavevector
            kzs(bad_indices)=-kzs(bad_indices); 
    
            %%%%%% define some  matrices and vectors %%%%%%%
            Zv=zeros(Numords,1);
            Dv=Zv;
            Dv(p)=1;
            Zv2=[Zv;Zv];
            Eye=eye(Numords);  % identity matrix
            Kx=diag(kx/k0);
            Kxsq=Kx.^2;
            Kzc=diag(kzc/k0);
            Kzcsq= Kzc.^2;
            Kzs=diag(kzs/k0);
            Kzssq= Kzs.^2;
            M=Numords-1;
            temp1=Eye/ns;
            fmat=Eye;
            gmat=j*Kzs/ns^2;
    
            for ii=Ngrat:-1:1
                epsg=ng(ii).^2;  % groove permittivity
                epsr=nr(ii).^2;   % ridge permittivity
                epsG=(1-Filfac(ii))*epsg+Filfac(ii)*epsr;  % average grating
                iepsG=(1-Filfac(ii))/epsg+Filfac(ii)/epsr; 
                Sinc=sin(pi*Filfac(ii)*(1:M))./(pi*(1:M));
                vm=(epsr-epsg)*fliplr(Sinc);
                v0=epsG;
                vp=(epsr-epsg)*Sinc;
                v=[vm v0 vp].*exp(+j*2*pi*Disp(ii)*(-M:M));
                ivm=(1/epsr-1/epsg)*fliplr(Sinc);
                iv0=iepsG;
                ivp=(1/epsr-1/epsg)*Sinc;
                iv=[ivm iv0 ivp].*exp(+j*2*pi*Disp(ii)*(-M:M));
                Epsilon=toeplitz(fliplr(v(1:Numords)),v(Numords:2*Numords-1));
                Alpha=toeplitz(fliplr(iv(1:Numords)),iv(Numords:2*Numords-1));
    
                clear Sinc v  vm  v0  vp 
                B=Kx*(Epsilon\Kx)-Eye;  % cofficient matrix
                [W,V]=eig(Alpha\B);    % W is the eigen vector and V are the eigen values
                Q=sqrt(V);
                M0=Alpha*W*Q;
                E=expm(-k0*Q*depth(ii));
                v=[W W;M0,-M0]\[fmat;gmat];
                temp2=v(1:Numords,:)\E;
                temp3=E*v(Numords+1:2*Numords,:)*temp2;
                temp1=temp1*temp2;
                fmat=W+W*temp3;
                gmat=M0-M0*temp3;
    
            end
    
            gfi=gmat/fmat;
            RHS=-gfi(:,p);
            RHS(p)=RHS(p)+j*kzc(p)/k0/epsc;
            Rs=(gfi+j*Kzc/nc^2)\RHS;
            Ts=(temp1/fmat)*(Rs+Dv)*nc;
         
    
           IR1=(abs(Rs).^2).*real(kzc./kzc(p));
           IT1=(abs(Ts).^2).*real(kzs./kzc(p));
    
           e=sum(IT1);
           f=sum(IR1);
           g=1-e-f;
           IT12(count)=e;
           IR12(count)=f; 
           loss(count)=g;
    
           if kkk == 1
               IR12_nd_133 = IR12;
           elseif kkk == 2
               IR12_nd_134 = IR12;
           end
    
           count=count+1;
     
           progress = (count/total_iterations)*100;
           clc;
           fprintf('nd = %.3f\t stacks=%d', nd, stack);
           fprintf('\t\t');
           fprintf('\nProgress: [%s%s] %.2f%%\r', repmat('=',1,floor(progress/2)), repmat(' ',1,50-floor(progress/2)), floor(progress));
           elapsedTime = toc(startTime);
           remainingTime = elapsedTime * (total_iterations / count - 1);
            
           % Convert to minutes and seconds
           elapsedMinutes = floor(elapsedTime / 60);
           elapsedSeconds = rem(elapsedTime, 60);
           remainingMinutes = floor(remainingTime / 60);
           remainingSeconds = rem(remainingTime, 60);
    
           % Display the elapsed time and estimated time left for the current inner loop in minutes and seconds
           fprintf('\nElapsed Time: %d minutes and %.2f seconds', elapsedMinutes, elapsedSeconds);
           fprintf('\nEstimated Time Left for Current Thickness: %d minutes and %.2f seconds', remainingMinutes, remainingSeconds);
    
           % Calculate and display the total estimated time left for all iterations
           % Adjust the remaining time estimate to account for the remaining outer loop iterations
           totalRemainingIterations = totalOuterIterations - completedOuterIterations + 1;
           totalRemainingTime = remainingTime * totalRemainingIterations;
           totalRemainingMinutes = floor(totalRemainingTime / 60);
           totalRemainingSeconds = rem(totalRemainingTime, 60);
           fprintf('\nTotal Estimated Time Left for All Thicknesses: %d minutes and %.2f seconds', totalRemainingMinutes, totalRemainingSeconds);
    
        end
    
        %%%%% Plotting results %%%%%
       
        % Finding the minimum y-value and its corresponding x-value
        [minVal, minIdx] = min(IR12);
        minW = w(minIdx);
        
        % Storing min values and respective x values for each kkk
        minVals(kkk) = minVal;
        minWs(kkk) = minW;
        
        % Plotting the curve for current nd value
        plot(w, IR12, 'DisplayName', ['nd = ', num2str(nd)]);
        
        % Resetting count for next iteration
        count=1;
        
        % Storing IR12 values for sensitivity calculation in next step
        IR12_store{kkk} = IR12;
        
    end    
    
    % Calculating sensitivity
    sensitivity = abs(minWs(1)-minWs(2)) / 0.01;
    sensitivities(stack) = sensitivity;
    
     % Finding all the local minima (dips) in the IR12 curve
    [mins, minIndices] = findpeaks(-IR12);
    
    % Finding the main dip in the IR12 curve for nd=1.33
    [minValue1, minIndex1] = min(IR12_nd_133);
    
    % Calculating the FWHM for the SPR dip
    half_max_value = (max(IR12_nd_133) + min(IR12_nd_133))/2;
    
    left_cross_index = find(IR12_nd_133(1:minIndex1) >= half_max_value, 1, 'last');
    right_cross_index = find(IR12_nd_133(minIndex1:end) >= half_max_value, 1, 'first') + minIndex1 - 1;
    FWHM = w(right_cross_index) - w(left_cross_index);
    
    % Calculating FOM, DA and QF
    FOM = sensitivity / FWHM;
    FOM_values(stack) = FOM; 
    DA = (1-min(IR12_nd_133))/FWHM;
    QF =sensitivity*DA; 
    
    scatter(w(left_cross_index), IR12_nd_133(left_cross_index), 'rx'); % Marking the left cross point with a red 'x'
    scatter(w(right_cross_index), IR12_nd_133(right_cross_index), 'bx'); % Marking the right cross point with a blue 'x'
    scatter(w(minIndex1), IR12_nd_133(minIndex1), 'bx')
    
    % sensitivity_values(n) = sensitivity;
    % QF_values(n) = QF; 
    
    % Displaying sensitivity on the plot
    dim = [0.65 0.2 0.25 0.25];
    str = {['Sensitivity: ', num2str(sensitivity)], ['FWHM: ', num2str(FWHM)],['FOM: ', num2str(FOM)], ['Min Angle: ', num2str(minWs(1))], ['Rmin (nd=1.34): ', num2str(min(IR12_nd_133))]...
        ['DA: ', num2str(DA)], ['QF: ', num2str(QF)]}; 
    annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','white');
    
    % Setting plot title and labels
    title(sprintf('SPR Curve for stacks = %d', stack));
    xlabel('Angle');
    ylabel('Reflectivity');
    legend;
    
    % Resetting minVals and minWs for next td value
    minVals = [];
    minWs = [];
end

% Plotting stacks vs sensitivity
figure;
plot(1:N, FOM_values, '-o');
xlabel('Number of Stacks');
ylabel('FOM');
title('Stacks vs FOM');
grid on;
