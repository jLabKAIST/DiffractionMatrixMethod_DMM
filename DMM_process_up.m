function [result] = DMM_process_up(param,up,rd,ru,tu)
%% AFTER PROCESSING PARAMETERS

small = 1e-12*1i;      % very small value to improve calculation speed
E_Coefficient = 1;
L_coefficient = 3/8;   % Factor at power inside
I_Coefficient = 3/16;  % Factor at power outside

loss = 0;
Dox = param.Dox;
Doy = param.Doy;
Doz = param.Doz;
x_dipole = param.x_dipole; % x-position of dipole

uxsize = param.uxsize;
uysize = param.uysize;
resolution = param.resolution;

k0 = param.k0;
n_org = up.n(1);
n_top = up.n(end);

Org_thickness = param.Org_thickness;
Dipole_pos = param.Dipole_pos;


% Matrices

lengthofux = length(uxsize);
Diffsize = ceil(lengthofux/resolution);

theta = zeros(resolution,Diffsize*2,length(uysize));
theta_air = zeros(resolution,Diffsize*2,length(uysize));
bplus = zeros(resolution,2*Diffsize,2*Diffsize,length(uysize));
bminus = zeros(resolution,2*Diffsize,2*Diffsize,length(uysize));
e12 = zeros(resolution,2*Diffsize,2*Diffsize,length(uysize));
e21 = zeros(resolution,2*Diffsize,2*Diffsize,length(uysize));

r1 = zeros(resolution,2*Diffsize,2*Diffsize,length(uysize));
t1 = zeros(resolution,2*Diffsize,2*Diffsize,length(uysize));
r2 = zeros(resolution,2*Diffsize,2*Diffsize,length(uysize));
a2 = zeros(resolution,2*Diffsize,2*Diffsize,length(uysize));
a1a2 = zeros(resolution,2*Diffsize,2*Diffsize,length(uysize));
inva1a2 = zeros(resolution,2*Diffsize,2*Diffsize,length(uysize));
cplus = zeros(resolution,2*Diffsize,2*Diffsize,length(uysize));
cminus = zeros(resolution,2*Diffsize,2*Diffsize,length(uysize));
invcpm = zeros(resolution,2*Diffsize,2*Diffsize,length(uysize));
invcmp = zeros(resolution,2*Diffsize,2*Diffsize,length(uysize));

Eplusx = zeros(resolution,Diffsize*2,length(uysize));
Eminusx = zeros(resolution,Diffsize*2,length(uysize));
Eplusy = zeros(resolution,Diffsize*2,length(uysize));
Eminusy = zeros(resolution,Diffsize*2,length(uysize));
Eplusz = zeros(resolution,Diffsize*2,length(uysize));
Eminusz = zeros(resolution,Diffsize*2,length(uysize));

EinTEx = zeros(resolution,Diffsize*2,length(uysize));
EinTMx = zeros(resolution,Diffsize*2,length(uysize));
Linx = zeros(resolution,Diffsize*2,length(uysize));

EinTEy = zeros(resolution,Diffsize*2,length(uysize));
EinTMy = zeros(resolution,Diffsize*2,length(uysize));
Liny = zeros(resolution,Diffsize*2,length(uysize));

EinTEz = zeros(resolution,Diffsize*2,length(uysize));
EinTMz = zeros(resolution,Diffsize*2,length(uysize));
Linz = zeros(resolution,Diffsize*2,length(uysize));

Eout1x = zeros(resolution,Diffsize*2,length(uysize));
Eout2x = zeros(resolution,Diffsize*2,length(uysize));
Eoutx = zeros(resolution,Diffsize*2,length(uysize));
Ioutx = zeros(resolution,Diffsize*2,length(uysize));
Koutx = zeros(resolution,Diffsize*2,length(uysize));

Eout1y = zeros(resolution,Diffsize*2,length(uysize));
Eout2y = zeros(resolution,Diffsize*2,length(uysize));
Eouty = zeros(resolution,Diffsize*2,length(uysize));
Iouty = zeros(resolution,Diffsize*2,length(uysize));
Kouty = zeros(resolution,Diffsize*2,length(uysize));

Eout1z = zeros(resolution,Diffsize*2,length(uysize));
Eout2z = zeros(resolution,Diffsize*2,length(uysize));
Eoutz = zeros(resolution,Diffsize*2,length(uysize));
Ioutz = zeros(resolution,Diffsize*2,length(uysize));
Koutz = zeros(resolution,Diffsize*2,length(uysize));

Kouttempx = zeros(resolution*Diffsize*2,length(uysize));
Iouttempx = zeros(resolution*Diffsize*2,length(uysize));

Lxytempx = zeros(resolution*Diffsize*2,length(uysize));
Lxyx = zeros(length(uxsize),length(uysize));

Kouttempy = zeros(resolution*Diffsize*2,length(uysize));
Iouttempy = zeros(resolution*Diffsize*2,length(uysize));

Lxytempy = zeros(resolution*Diffsize*2,length(uysize));
Lxyy = zeros(length(uxsize),length(uysize));

Kouttempz = zeros(resolution*Diffsize*2,length(uysize));
Iouttempz = zeros(resolution*Diffsize*2,length(uysize));

Lxytempz = zeros(resolution*Diffsize*2,length(uysize));
Lxyz = zeros(length(uxsize),length(uysize));


for z = 1:resolution
for uyindex = 1:length(uysize)
for uxindex = 1:Diffsize

    if (uxindex-1)*resolution+z<=lengthofux

        ux = uxsize((uxindex-1)*resolution+z);
        uy = uysize(uyindex);

        if uy == 0
            uy = 1e-8;
        end
        ull = sqrt(ux^2+uy^2);
        uz = sqrt(1-ull^2);
        if abs(uz) < 1e-6
            uz = 1e-6;
        end

        theta(z,uxindex,uyindex) = asind(sqrt(ux^2+uy^2));
        theta(z,uxindex+Diffsize,uyindex) = asind(sqrt(ux^2+uy^2));

        % propagation matrices
        % bplus : from dipole to top layer
        % bmins : from dipole to bottom layer
        % e12 = e21 : from top to bottom or bottom to top
        bplus(z,uxindex,uxindex,uyindex) = small+exp(1i*k0*(n_org*cosd(theta(z,uxindex,uyindex)+loss))*(Org_thickness-Dipole_pos));
        bplus(z,uxindex+Diffsize,uxindex+Diffsize,uyindex) = small+exp(1i*k0*(n_org*cosd(theta(z,uxindex,uyindex)+loss))*(Org_thickness-Dipole_pos));
        bminus(z,uxindex,uxindex,uyindex) = small+exp(1i*k0*(n_org*cosd(theta(z,uxindex,uyindex)+loss))*(Dipole_pos));
        bminus(z,uxindex+Diffsize,uxindex+Diffsize,uyindex) = small+exp(1i*k0*(n_org*cosd(theta(z,uxindex,uyindex)+loss))*(Dipole_pos));

        e12(z,uxindex,uxindex,uyindex) = small+exp(1i*k0*(n_org*cosd(theta(z,uxindex,uyindex)+loss))*(Org_thickness));
        e12(z,uxindex+Diffsize,uxindex+Diffsize,uyindex) = small+exp(1i*k0*(n_org*cosd(theta(z,uxindex,uyindex)+loss))*(Org_thickness));
        e21(z,uxindex,uxindex,uyindex) = small+exp(1i*k0*(n_org*cosd(theta(z,uxindex,uyindex)+loss))*(Org_thickness));
        e21(z,uxindex+Diffsize,uxindex+Diffsize,uyindex) = small+exp(1i*k0*(n_org*cosd(theta(z,uxindex,uyindex)+loss))*(Org_thickness));

        % reflectance by RCWA

        theta_air(z,uxindex,uyindex) = asind(sind(theta(z,uxindex,uyindex))*n_org/n_top);
        theta_air(z,uxindex+Diffsize,uyindex) = asind(sind(theta(z,uxindex,uyindex))*n_org/n_top);

        % r1(ru), r2(rd), t1(tu)
        % |  TETE0(ux)     TETE1(ux)     TETE2(ux)     TETE3(ux)      TETM
        % |  TETEm1(ux+G)  TETE0(ux+G)   TETE1(ux+G)   TETE2(ux+G)    TETM
        % |  TETEm2(ux+2G) TETEm1(ux+2G) TETE0(ux+2G)  TETE1(ux+2G)   TETM
        % |  TETEm3(ux+3G) TETEm2(ux+3G) TETEm1(ux+3G) TETE0(ux+3G)   TETM
        % |
        % |  TMTE0(ux)     TMTE1(ux)     TMTE2(ux)     TMTE3(ux)      TMTM
        % |  TMTEm1(ux+G)  TMTE0(ux+G)   TMTE1(ux+G)   TMTE2(ux+G)    TMTM
        % |  TMTEm2(ux+2G) TMTEm1(ux+2G) TMTE0(ux+2G)  TMTE1(ux+2G)   TMTM
        % |  TMTEm3(ux+3G) TMTEm2(ux+3G) TMTEm1(ux+3G) TMTE0(ux+3G)   TMTM
        %
        % small : make righten side complex number -> faster filling matrix
        r1(z,uxindex,uxindex,uyindex) = small+ru.TETE((uxindex-1)*resolution+z,uyindex,11);
        r1(z,uxindex+Diffsize,uxindex+Diffsize,uyindex) = small+ru.TMTM((uxindex-1)*resolution+z,uyindex,11);
        r1(z,uxindex,uxindex+Diffsize,uyindex) = small+ru.TETM((uxindex-1)*resolution+z,uyindex,11);
        r1(z,uxindex+Diffsize,uxindex,uyindex) = small+ru.TMTE((uxindex-1)*resolution+z,uyindex,11);

        for kk = 0:9
            if uxindex< Diffsize-kk
                r1(z,uxindex,uxindex+1+kk,uyindex) = small+ru.TETE((uxindex-1)*resolution+z,uyindex,12+kk);
                r1(z,uxindex+Diffsize,uxindex+Diffsize+1+kk,uyindex) = small+ru.TMTM((uxindex-1)*resolution+z,uyindex,12+kk);
                r1(z,uxindex,uxindex+Diffsize+1+kk,uyindex) = small+ru.TETM((uxindex-1)*resolution+z,uyindex,12+kk);
                r1(z,uxindex+Diffsize,uxindex+1+kk,uyindex) = small+ru.TMTE((uxindex-1)*resolution+z,uyindex,12+kk);
            end

            if uxindex>1+kk
                r1(z,uxindex,uxindex-1-kk,uyindex) = small+ru.TETE((uxindex-1)*resolution+z,uyindex,10-kk);
                r1(z,uxindex+Diffsize,uxindex+Diffsize-1-kk,uyindex) = small+ru.TMTM((uxindex-1)*resolution+z,uyindex,10-kk);
                r1(z,uxindex,uxindex+Diffsize-1-kk,uyindex) = small+ru.TETM((uxindex-1)*resolution+z,uyindex,10-kk);
                r1(z,uxindex+Diffsize,uxindex-1-kk,uyindex) = small+ru.TMTE((uxindex-1)*resolution+z,uyindex,10-kk);
            end
        end % for kk = 0:9 
        
        % t1 = t_up
        t1(z,uxindex,uxindex,uyindex) = small+tu.TETE((uxindex-1)*resolution+z,uyindex,11);
        t1(z,uxindex+Diffsize,uxindex+Diffsize,uyindex) = small+tu.TMTM((uxindex-1)*resolution+z,uyindex,11);
        t1(z,uxindex,uxindex+Diffsize,uyindex) = small+tu.TETM((uxindex-1)*resolution+z,uyindex,11);
        t1(z,uxindex+Diffsize,uxindex,uyindex) = small+tu.TMTE((uxindex-1)*resolution+z,uyindex,11);

        for kk = 0:9
            if uxindex< Diffsize-kk
                t1(z,uxindex,uxindex+1+kk,uyindex) = small+tu.TETE((uxindex-1)*resolution+z,uyindex,12+kk);
                t1(z,uxindex+Diffsize,uxindex+Diffsize+1+kk,uyindex) = small+tu.TMTM((uxindex-1)*resolution+z,uyindex,12+kk);
                t1(z,uxindex,uxindex+Diffsize+1+kk,uyindex) = small+tu.TETM((uxindex-1)*resolution+z,uyindex,12+kk);
                t1(z,uxindex+Diffsize,uxindex+1+kk,uyindex) = small+tu.TMTE((uxindex-1)*resolution+z,uyindex,12+kk);
            end

            if uxindex>1+kk
                t1(z,uxindex,uxindex-1-kk,uyindex) = small+tu.TETE((uxindex-1)*resolution+z,uyindex,10-kk);
                t1(z,uxindex+Diffsize,uxindex+Diffsize-1-kk,uyindex) = small+tu.TMTM((uxindex-1)*resolution+z,uyindex,10-kk);
                t1(z,uxindex,uxindex+Diffsize-1-kk,uyindex) = small+tu.TETM((uxindex-1)*resolution+z,uyindex,10-kk);
                t1(z,uxindex+Diffsize,uxindex-1-kk,uyindex) = small+tu.TMTE((uxindex-1)*resolution+z,uyindex,10-kk);
            end
        end % for kk = 0:9 
        
        % r2 = diffraction at bottom
        r2(z,uxindex,uxindex,uyindex) = small+rd.TETE((uxindex-1)*resolution+z,uyindex,11);
        r2(z,uxindex+Diffsize,uxindex+Diffsize,uyindex) = small+rd.TMTM((uxindex-1)*resolution+z,uyindex,11);
        r2(z,uxindex,uxindex+Diffsize,uyindex) = small+rd.TETM((uxindex-1)*resolution+z,uyindex,11);
        r2(z,uxindex+Diffsize,uxindex,uyindex) = small+rd.TMTE((uxindex-1)*resolution+z,uyindex,11);

        for kk = 0:9
            if uxindex< Diffsize-kk
                r2(z,uxindex,uxindex+1+kk,uyindex) = small+rd.TETE((uxindex-1)*resolution+z,uyindex,12+kk);
                r2(z,uxindex+Diffsize,uxindex+Diffsize+1+kk,uyindex) = small+rd.TMTM((uxindex-1)*resolution+z,uyindex,12+kk);
                r2(z,uxindex,uxindex+Diffsize+1+kk,uyindex) = small+rd.TETM((uxindex-1)*resolution+z,uyindex,12+kk);
                r2(z,uxindex+Diffsize,uxindex+1+kk,uyindex) = small+rd.TMTE((uxindex-1)*resolution+z,uyindex,12+kk);
            end

            if uxindex>1+kk
                r2(z,uxindex,uxindex-1-kk,uyindex) = small+rd.TETE((uxindex-1)*resolution+z,uyindex,10-kk);
                r2(z,uxindex+Diffsize,uxindex+Diffsize-1-kk,uyindex) = small+rd.TMTM((uxindex-1)*resolution+z,uyindex,10-kk);
                r2(z,uxindex,uxindex+Diffsize-1-kk,uyindex) = small+rd.TETM((uxindex-1)*resolution+z,uyindex,10-kk);
                r2(z,uxindex+Diffsize,uxindex-1-kk,uyindex) = small+rd.TMTE((uxindex-1)*resolution+z,uyindex,10-kk);
            end
        end % for kk = 0:9
    end % if uxsize()<
end % uxsize

a2(z,:,:,uyindex) = squeeze(r2(z,:,:,uyindex))*squeeze(e21(z,:,:,uyindex));
a1a2(z,:,:,uyindex) = squeeze(r1(z,:,:,uyindex))*squeeze(e12(z,:,:,uyindex))*squeeze(r2(z,:,:,uyindex))*squeeze(e21(z,:,:,uyindex));
inva1a2(z,:,:,uyindex) = inv(eye(2*Diffsize)-squeeze(a1a2(z,:,:,uyindex)));

% c's -> for calculation of electric field at dipole
cplus(z,:,:,uyindex) = squeeze(bplus(z,:,:,uyindex))*squeeze(r1(z,:,:,uyindex))*squeeze(bplus(z,:,:,uyindex));
cminus(z,:,:,uyindex) = squeeze(bminus(z,:,:,uyindex))*squeeze(r2(z,:,:,uyindex))*squeeze(bminus(z,:,:,uyindex));
invcpm(z,:,:,uyindex) = inv(eye(2*Diffsize)-squeeze(cplus(z,:,:,uyindex))*squeeze(cminus(z,:,:,uyindex)));
invcmp(z,:,:,uyindex) = inv(eye(2*Diffsize)-squeeze(cminus(z,:,:,uyindex))*squeeze(cplus(z,:,:,uyindex)));

end % uysize
end % z

timeMatrix = toc



for z = 1:resolution
for uyindex = 1:length(uysize)
for uxindex = 1:Diffsize
    
    if (uxindex-1)*resolution+z<=lengthofux
        
        ux = uxsize((uxindex-1)*resolution+z);
        uy = uysize(uyindex);

        if uy == 0
            uy = 1e-8;
        end
        ull = sqrt(ux^2+uy^2);
        uz = sqrt(1-ull^2);
        if abs(uz) < 1e-6
            uz = 1e-6;
        end
        
        Expe = exp(-1i*k0*n_org*x_dipole*ux);
        
        % E field ( x-dipole )
        if Dox == 1
            Eplusx(z,uxindex,uyindex) = E_Coefficient * uy/uz/ull * Expe; %% Up TE
            Eplusx(z,uxindex+Diffsize,uyindex) = E_Coefficient * ux/ull * Expe; %% Up TM
            Eminusx(z,uxindex,uyindex) = E_Coefficient * uy/uz/ull * Expe; %% Down TE
            Eminusx(z,uxindex+Diffsize,uyindex) = -E_Coefficient * ux/ull * Expe; %% Down TM
        end
        % E field ( y-dipole )
        if Doy == 1
            Eplusy(z,uxindex,uyindex) = E_Coefficient * -ux/uz/ull * Expe;
            Eplusy(z,uxindex+Diffsize,uyindex) = E_Coefficient * uy/ull * Expe;
            Eminusy(z,uxindex,uyindex) = E_Coefficient * -ux/uz/ull * Expe;
            Eminusy(z,uxindex+Diffsize,uyindex) = -E_Coefficient * uy/ull * Expe;
        end
        % E field ( z-dipole )
        if Doz ==1
            Eplusz(z,uxindex,uyindex) = E_Coefficient * 0 * Expe;
            Eplusz(z,uxindex+Diffsize,uyindex) = E_Coefficient * ull/uz * Expe;
            Eminusz(z,uxindex,uyindex) = E_Coefficient * 0 * Expe;
            Eminusz(z,uxindex+Diffsize,uyindex) = E_Coefficient * ull/uz * Expe;
        end
    end
end

A_ITO = squeeze(bplus(z,:,:,uyindex))*squeeze(inva1a2(z,:,:,uyindex))*squeeze((t1(z,:,:,uyindex)));
B_ITO = squeeze(bminus(z,:,:,uyindex))*squeeze(a2(z,:,:,uyindex))*squeeze(inva1a2(z,:,:,uyindex))*squeeze(t1(z,:,:,uyindex));

Ap_ITOTE = eye(2*Diffsize) + squeeze(cplus(z,:,:,uyindex))*squeeze(cminus(z,:,:,uyindex))*squeeze(invcpm(z,:,:,uyindex)) + squeeze(cplus(z,:,:,uyindex))*squeeze(invcmp(z,:,:,uyindex));
Bp_ITOTE = squeeze(cminus(z,:,:,uyindex))*squeeze(invcpm(z,:,:,uyindex)) + squeeze(cminus(z,:,:,uyindex))*squeeze(cplus(z,:,:,uyindex))*squeeze(invcmp(z,:,:,uyindex));
Ap_ITOTM = eye(2*Diffsize) + squeeze(cplus(z,:,:,uyindex))*squeeze(cminus(z,:,:,uyindex))*squeeze(invcpm(z,:,:,uyindex)) - squeeze(cplus(z,:,:,uyindex))*squeeze(invcmp(z,:,:,uyindex));
Bp_ITOTM = squeeze(cminus(z,:,:,uyindex))*squeeze(invcpm(z,:,:,uyindex)) - squeeze(cminus(z,:,:,uyindex))*squeeze(cplus(z,:,:,uyindex))*squeeze(invcmp(z,:,:,uyindex));

Az_ITOTM = eye(2*Diffsize) + squeeze(cplus(z,:,:,uyindex))*squeeze(cminus(z,:,:,uyindex))*squeeze(invcpm(z,:,:,uyindex)) + squeeze(cplus(z,:,:,uyindex))*squeeze(invcmp(z,:,:,uyindex));
Bz_ITOTM = squeeze(cminus(z,:,:,uyindex))*squeeze(invcpm(z,:,:,uyindex)) + squeeze(cminus(z,:,:,uyindex))*squeeze(cplus(z,:,:,uyindex))*squeeze(invcmp(z,:,:,uyindex));

% Eout
if Dox == 1
    Eout1x(z,:,uyindex) = squeeze(Eplusx(z,:,uyindex))*A_ITO;
    Eout2x(z,:,uyindex) = squeeze(Eminusx(z,:,uyindex))*B_ITO;
    Eoutx(z,:,uyindex) = Eout1x(z,:,uyindex) + Eout2x(z,:,uyindex);
    % Eoutx(z,:,uyindex) = Eout1x(z,:,uyindex) + Eout2x(z,:,uyindex);

    C1_ITO = conj(squeeze(Eplusx(z,:,uyindex)))'.*A_ITO;
    C2_ITO = conj(squeeze(Eminusx(z,:,uyindex)))'.*B_ITO;
    %                 C_ITO = abs(C1_ITO + C2_ITO).^2;
    C_ITO = zeros(Diffsize,Diffsize*2);
    for iss = 1:Diffsize
        C_ITO(iss,:) = abs(C1_ITO(iss,:) + C1_ITO(iss+Diffsize,:) + C2_ITO(iss,:) + C2_ITO(iss+Diffsize,:)).^2;
    end
    if param.singledip ==1
        Ioutx(z,:,uyindex) = abs(Eoutx(z,:,uyindex)).^2 .* n_top/n_org .* cosd(theta_air(z,:,uyindex)) ./ cosd(theta(z,:,uyindex))*  I_Coefficient .* max(0,cosd(theta_air(z,:,uyindex)).^2);     
    else
        Ioutx(z,:,uyindex) = sum(C_ITO,1) .* n_top/n_org .* cosd(theta_air(z,:,uyindex)) ./ cosd(theta(z,:,uyindex)).* I_Coefficient .* max(0,cosd(theta_air(z,:,uyindex)).^2);
    end
    Koutx(z,:,uyindex) = Ioutx(z,:,uyindex)./ max(0,cosd(theta_air(z,:,uyindex)).^2);
    if param.singledip == 1
        EinTEx(z,:,uyindex) = squeeze(Eplusx(z,:,uyindex)) + ( squeeze(Eminusx(z,:,uyindex))*squeeze(cminus(z,:,:,uyindex)) + squeeze(Eplusx(z,:,uyindex))*squeeze(cplus(z,:,:,uyindex))*squeeze(cminus(z,:,:,uyindex))  )*squeeze(invcpm(z,:,:,uyindex));
        EinTEx(z,:,uyindex) = EinTEx(z,:,uyindex) + ( squeeze(Eplusx(z,:,uyindex))*squeeze(cplus(z,:,:,uyindex)) + squeeze(Eminusx(z,:,uyindex))*squeeze(cminus(z,:,:,uyindex))*squeeze(cplus(z,:,:,uyindex))  )*squeeze(invcmp(z,:,:,uyindex));
        EinTMx(z,:,uyindex) = +squeeze(Eplusx(z,:,uyindex)) + ( squeeze(Eminusx(z,:,uyindex))*squeeze(cminus(z,:,:,uyindex)) + squeeze(Eplusx(z,:,uyindex))*squeeze(cplus(z,:,:,uyindex))*squeeze(cminus(z,:,:,uyindex))  )*squeeze(invcpm(z,:,:,uyindex));
        EinTMx(z,:,uyindex) = EinTMx(z,:,uyindex) - ( squeeze(Eplusx(z,:,uyindex))*squeeze(cplus(z,:,:,uyindex)) + squeeze(Eminusx(z,:,uyindex))*squeeze(cminus(z,:,:,uyindex))*squeeze(cplus(z,:,:,uyindex))  )*squeeze(invcmp(z,:,:,uyindex));
    else
        for iss = 1:Diffsize
            EinTEx(z,iss,uyindex) = Eplusx(z,iss,uyindex) * Ap_ITOTE(iss,iss) + Eplusx(z,iss+Diffsize,uyindex) * Ap_ITOTE(iss+Diffsize,iss) + Eminusx(z,iss,uyindex) * Bp_ITOTE(iss,iss) + Eminusx(z,iss+Diffsize,uyindex) * Bp_ITOTE(iss+Diffsize,iss);
            EinTEx(z,iss+Diffsize,uyindex) = Eplusx(z,iss,uyindex) * Ap_ITOTE(iss,iss+Diffsize) + Eplusx(z,iss+Diffsize,uyindex) * Ap_ITOTE(iss+Diffsize,iss+Diffsize) + Eminusx(z,iss,uyindex) * Bp_ITOTE(iss,iss+Diffsize) + Eminusx(z,iss+Diffsize,uyindex) * Bp_ITOTE(iss+Diffsize,iss+Diffsize);
            EinTMx(z,iss,uyindex) = Eplusx(z,iss,uyindex) * Ap_ITOTM(iss,iss) + Eplusx(z,iss+Diffsize,uyindex) * Ap_ITOTM(iss+Diffsize,iss) + Eminusx(z,iss,uyindex) * Bp_ITOTM(iss,iss) + Eminusx(z,iss+Diffsize,uyindex) * Bp_ITOTM(iss+Diffsize,iss);
            EinTMx(z,iss+Diffsize,uyindex) = Eplusx(z,iss,uyindex) * Ap_ITOTM(iss,iss+Diffsize) + Eplusx(z,iss+Diffsize,uyindex) * Ap_ITOTM(iss+Diffsize,iss+Diffsize) + Eminusx(z,iss,uyindex) * Bp_ITOTM(iss,iss+Diffsize) + Eminusx(z,iss+Diffsize,uyindex) * Bp_ITOTM(iss+Diffsize,iss+Diffsize);
        end
    end
end

if Doy == 1
    Eout1y(z,:,uyindex) = squeeze(Eplusy(z,:,uyindex))*A_ITO;
    Eout2y(z,:,uyindex) = squeeze(Eminusy(z,:,uyindex))*B_ITO;
    Eouty(z,:,uyindex) = Eout1y(z,:,uyindex) + Eout2y(z,:,uyindex);
   
    C1_ITO = conj(squeeze(Eplusy(z,:,uyindex)))'.*A_ITO;
    C2_ITO = conj(squeeze(Eminusy(z,:,uyindex)))'.*B_ITO;
    C_ITO = zeros(Diffsize,Diffsize*2);
    for iss = 1:Diffsize
        C_ITO(iss,:) = abs(C1_ITO(iss,:) + C1_ITO(iss+Diffsize,:) + C2_ITO(iss,:) + C2_ITO(iss+Diffsize,:)).^2;
    end
    if param.singledip ==1
        Iouty(z,:,uyindex) = abs(Eouty(z,:,uyindex)).^2 .* n_top/n_org .* cosd(theta_air(z,:,uyindex)) ./ cosd(theta(z,:,uyindex))*  I_Coefficient .* max(0,cosd(theta_air(z,:,uyindex)).^2);     
    else
        Iouty(z,:,uyindex) = sum(C_ITO,1) .* n_top/n_org .* cosd(theta_air(z,:,uyindex)) ./ cosd(theta(z,:,uyindex)).* I_Coefficient .* max(0,cosd(theta_air(z,:,uyindex)).^2);
    end
    Kouty(z,:,uyindex) = Iouty(z,:,uyindex)./ max(0,cosd(theta_air(z,:,uyindex)).^2);
    if param.singledip == 1
        EinTEy(z,:,uyindex) = squeeze(Eplusy(z,:,uyindex)) + ( squeeze(Eminusy(z,:,uyindex))*squeeze(cminus(z,:,:,uyindex)) + squeeze(Eplusy(z,:,uyindex))*squeeze(cplus(z,:,:,uyindex))*squeeze(cminus(z,:,:,uyindex))  )*squeeze(invcpm(z,:,:,uyindex));
        EinTEy(z,:,uyindex) = EinTEy(z,:,uyindex) + ( squeeze(Eplusy(z,:,uyindex))*squeeze(cplus(z,:,:,uyindex)) + squeeze(Eminusy(z,:,uyindex))*squeeze(cminus(z,:,:,uyindex))*squeeze(cplus(z,:,:,uyindex))  )*squeeze(invcmp(z,:,:,uyindex));
        EinTMy(z,:,uyindex) = +squeeze(Eplusy(z,:,uyindex)) + ( squeeze(Eminusy(z,:,uyindex))*squeeze(cminus(z,:,:,uyindex)) + squeeze(Eplusy(z,:,uyindex))*squeeze(cplus(z,:,:,uyindex))*squeeze(cminus(z,:,:,uyindex))  )*squeeze(invcpm(z,:,:,uyindex));
        EinTMy(z,:,uyindex) = EinTMy(z,:,uyindex) - ( squeeze(Eplusy(z,:,uyindex))*squeeze(cplus(z,:,:,uyindex)) + squeeze(Eminusy(z,:,uyindex))*squeeze(cminus(z,:,:,uyindex))*squeeze(cplus(z,:,:,uyindex))  )*squeeze(invcmp(z,:,:,uyindex));

    else

        for iss = 1:Diffsize
            EinTEy(z,iss,uyindex) = Eplusy(z,iss,uyindex) * Ap_ITOTE(iss,iss) + Eplusy(z,iss+Diffsize,uyindex) * Ap_ITOTE(iss+Diffsize,iss) + Eminusy(z,iss,uyindex) * Bp_ITOTE(iss,iss) + Eminusy(z,iss+Diffsize,uyindex) * Bp_ITOTE(iss+Diffsize,iss);
            EinTEy(z,iss+Diffsize,uyindex) = Eplusy(z,iss,uyindex) * Ap_ITOTE(iss,iss+Diffsize) + Eplusy(z,iss+Diffsize,uyindex) * Ap_ITOTE(iss+Diffsize,iss+Diffsize) + Eminusy(z,iss,uyindex) * Bp_ITOTE(iss,iss+Diffsize) + Eminusy(z,iss+Diffsize,uyindex) * Bp_ITOTE(iss+Diffsize,iss+Diffsize);
            EinTMy(z,iss,uyindex) = Eplusy(z,iss,uyindex) * Ap_ITOTM(iss,iss) + Eplusy(z,iss+Diffsize,uyindex) * Ap_ITOTM(iss+Diffsize,iss) + Eminusy(z,iss,uyindex) * Bp_ITOTM(iss,iss) + Eminusy(z,iss+Diffsize,uyindex) * Bp_ITOTM(iss+Diffsize,iss);
            EinTMy(z,iss+Diffsize,uyindex) = Eplusy(z,iss,uyindex) * Ap_ITOTM(iss,iss+Diffsize) + Eplusy(z,iss+Diffsize,uyindex) * Ap_ITOTM(iss+Diffsize,iss+Diffsize) + Eminusy(z,iss,uyindex) * Bp_ITOTM(iss,iss+Diffsize) + Eminusy(z,iss+Diffsize,uyindex) * Bp_ITOTM(iss+Diffsize,iss+Diffsize);
        end
    end
end

if Doz == 1
    Eout1z(z,:,uyindex) = squeeze(Eplusz(z,:,uyindex))*A_ITO;
    Eout2z(z,:,uyindex) = squeeze(Eminusz(z,:,uyindex))*B_ITO;
    Eoutz(z,:,uyindex) = Eout1z(z,:,uyindex) + Eout2z(z,:,uyindex);
   
    C1_ITO = conj(squeeze(Eplusz(z,:,uyindex)))'.*A_ITO;
    C2_ITO = conj(squeeze(Eminusz(z,:,uyindex)))'.*B_ITO;
    C_ITO = zeros(Diffsize,Diffsize*2);
    for iss = 1:Diffsize
        C_ITO(iss,:) = abs(C1_ITO(iss,:) + C1_ITO(iss+Diffsize,:) + C2_ITO(iss,:) + C2_ITO(iss+Diffsize,:)).^2;
    end
    if param.singledip ==1
        Ioutz(z,:,uyindex) = abs(Eoutz(z,:,uyindex)).^2 .* n_top/n_org .* cosd(theta_air(z,:,uyindex)) ./ cosd(theta(z,:,uyindex))*  I_Coefficient .* max(0,cosd(theta_air(z,:,uyindex)).^2);     
    else
        Ioutz(z,:,uyindex) = sum(C_ITO,1) .* n_top/n_org .* cosd(theta_air(z,:,uyindex)) ./ cosd(theta(z,:,uyindex)).* I_Coefficient .* max(0,cosd(theta_air(z,:,uyindex)).^2);
    end
    Koutz(z,:,uyindex) = Ioutz(z,:,uyindex)./ max(0,cosd(theta_air(z,:,uyindex)).^2);
    if param.singledip == 1

        EinTEz(z,:,uyindex) = squeeze(Eplusz(z,:,uyindex)) + ( squeeze(Eminusz(z,:,uyindex))*squeeze(cminus(z,:,:,uyindex)) + squeeze(Eplusz(z,:,uyindex))*squeeze(cplus(z,:,:,uyindex))*squeeze(cminus(z,:,:,uyindex))  )*squeeze(invcpm(z,:,:,uyindex));
        EinTEz(z,:,uyindex) = EinTEz(z,:,uyindex) + ( squeeze(Eplusz(z,:,uyindex))*squeeze(cplus(z,:,:,uyindex)) + squeeze(Eminusz(z,:,uyindex))*squeeze(cminus(z,:,:,uyindex))*squeeze(cplus(z,:,:,uyindex))  )*squeeze(invcmp(z,:,:,uyindex));
        EinTMz(z,:,uyindex) = squeeze(Eplusz(z,:,uyindex)) + ( squeeze(Eminusz(z,:,uyindex))*squeeze(cminus(z,:,:,uyindex)) + squeeze(Eplusz(z,:,uyindex))*squeeze(cplus(z,:,:,uyindex))*squeeze(cminus(z,:,:,uyindex))  )*squeeze(invcpm(z,:,:,uyindex));
        EinTMz(z,:,uyindex) = EinTMz(z,:,uyindex) + ( squeeze(Eplusz(z,:,uyindex))*squeeze(cplus(z,:,:,uyindex)) + squeeze(Eminusz(z,:,uyindex))*squeeze(cminus(z,:,:,uyindex))*squeeze(cplus(z,:,:,uyindex))  )*squeeze(invcmp(z,:,:,uyindex));
    else
        for iss = 1:Diffsize
            EinTEz(z,iss,uyindex) = Eplusz(z,iss,uyindex) * Ap_ITOTE(iss,iss) + Eplusz(z,iss+Diffsize,uyindex) * Ap_ITOTE(iss+Diffsize,iss) + Eminusz(z,iss,uyindex) * Bp_ITOTE(iss,iss) + Eminusz(z,iss+Diffsize,uyindex) * Bp_ITOTE(iss+Diffsize,iss);
            EinTEz(z,iss+Diffsize,uyindex) = Eplusz(z,iss,uyindex) * Ap_ITOTE(iss,iss+Diffsize) + Eplusz(z,iss+Diffsize,uyindex) * Ap_ITOTE(iss+Diffsize,iss+Diffsize) + Eminusz(z,iss,uyindex) * Bp_ITOTE(iss,iss+Diffsize) + Eminusz(z,iss+Diffsize,uyindex) * Bp_ITOTE(iss+Diffsize,iss+Diffsize);
            EinTMz(z,iss,uyindex) = Eplusz(z,iss,uyindex) * Az_ITOTM(iss,iss) + Eplusz(z,iss+Diffsize,uyindex) * Az_ITOTM(iss+Diffsize,iss) + Eminusz(z,iss,uyindex) * Bz_ITOTM(iss,iss) + Eminusz(z,iss+Diffsize,uyindex) * Bz_ITOTM(iss+Diffsize,iss);
            EinTMz(z,iss+Diffsize,uyindex) = Eplusz(z,iss,uyindex) * Az_ITOTM(iss,iss+Diffsize) + Eplusz(z,iss+Diffsize,uyindex) * Az_ITOTM(iss+Diffsize,iss+Diffsize) + Eminusz(z,iss,uyindex) * Bz_ITOTM(iss,iss+Diffsize) + Eminusz(z,iss+Diffsize,uyindex) * Bz_ITOTM(iss+Diffsize,iss+Diffsize);
        end
    end
end
end
end

% Lin generator

for z = 1:resolution
for uxindex = 1:Diffsize
    for uyindex = 1:length(uysize)
        if (uxindex-1)*resolution+z<=lengthofux

            ux = uxsize((uxindex-1)*resolution+z);
            uy = uysize(uyindex);
           
            if uy == 0
                uy = 1e-8;
            end
            ull = sqrt(ux^2+uy^2);
            uz = sqrt(1-ull^2);
            if abs(uz) < 1e-6
                uz = 1e-6;
            end
            Expe2 = exp(1i*k0*n_org*x_dipole*ux); %%% Expe 무효화

            if Dox == 1
                Linx(z,uxindex,uyindex) = real(EinTEx(z,uxindex,uyindex) * uy/ull * Expe2); % x_dipole TE
                Linx(z,uxindex+Diffsize,uyindex) = real(EinTMx(z,uxindex+Diffsize,uyindex) * (uz*ux/ull) * Expe2); % x_dipole TM
            end

            if Doy ==1
                Liny(z,uxindex,uyindex) = real(EinTEy(z,uxindex,uyindex) * -ux/ull * Expe2); % y_dipole TE
                Liny(z,uxindex+Diffsize,uyindex) = real(EinTMy(z,uxindex+Diffsize,uyindex) * (uz*uy/ull) * Expe2); % y_dipole TM
            end

            if Doz ==1
                Linz(z,uxindex+Diffsize,uyindex) = real(EinTMz(z,uxindex+Diffsize,uyindex) * ull * Expe2); % z_dipole TM
            end
        end
    end
end
end


% Make L,K matrix
Koutx(isnan(Koutx))=0;
Kouty(isnan(Kouty))=0;
Koutz(isnan(Koutz))=0;

for z = 1:resolution
for uxindex = 1:Diffsize
    for uyindex = 1:length(uysize)
        if (uxindex-1)*resolution+z<=lengthofux

            ux = uxsize((uxindex-1)*resolution+z);
            uy = uysize(uyindex);
            ull = sqrt(ux^2+uy^2);
            uz = sqrt(1-ull^2);
            uz_air = sqrt(1- (n_org/n_top)^2 * ull^2);
            if abs(uz) < 1e-8
                uz = 0;
            end
            
            if Dox ==1
                Iouttempx((uxindex-1)*resolution+z,uyindex) = Ioutx(z,uxindex,uyindex)*real(uz/uz_air);
                Iouttempx((uxindex-1)*resolution+z+Diffsize*resolution,uyindex) = Ioutx(z,uxindex+Diffsize,uyindex)*real(uz/uz_air);

                Kouttempx((uxindex-1)*resolution+z,uyindex) = Koutx(z,uxindex,uyindex)*uz;
                Kouttempx((uxindex-1)*resolution+z+Diffsize*resolution,uyindex) = Koutx(z,uxindex+Diffsize,uyindex)*uz;

                Lxytempx((uxindex-1)*resolution+z,uyindex) = Linx(z,uxindex,uyindex);
                Lxytempx((uxindex-1)*resolution+z+Diffsize*resolution,uyindex) = Linx(z,uxindex+Diffsize,uyindex);
            end

            if Doy ==1
                Iouttempy((uxindex-1)*resolution+z,uyindex) = Iouty(z,uxindex,uyindex)*real(uz/uz_air);
                Iouttempy((uxindex-1)*resolution+z+Diffsize*resolution,uyindex) = Iouty(z,uxindex+Diffsize,uyindex)*real(uz/uz_air);

                Kouttempy((uxindex-1)*resolution+z,uyindex) = Kouty(z,uxindex,uyindex)*uz;
                Kouttempy((uxindex-1)*resolution+z+Diffsize*resolution,uyindex) = Kouty(z,uxindex+Diffsize,uyindex)*uz;

                Lxytempy((uxindex-1)*resolution+z,uyindex) = Liny(z,uxindex,uyindex);
                Lxytempy((uxindex-1)*resolution+z+Diffsize*resolution,uyindex) = Liny(z,uxindex+Diffsize,uyindex);
            end

            if Doz ==1
                Iouttempz((uxindex-1)*resolution+z,uyindex) = Ioutz(z,uxindex,uyindex)*real(uz/uz_air);
                Iouttempz((uxindex-1)*resolution+z+Diffsize*resolution,uyindex) = Ioutz(z,uxindex+Diffsize,uyindex)*real(uz/uz_air);

                Kouttempz((uxindex-1)*resolution+z,uyindex) = Koutz(z,uxindex,uyindex)*uz;
                Kouttempz((uxindex-1)*resolution+z+Diffsize*resolution,uyindex) = Koutz(z,uxindex+Diffsize,uyindex)*uz;

                Lxytempz((uxindex-1)*resolution+z,uyindex) = Linz(z,uxindex,uyindex);
                Lxytempz((uxindex-1)*resolution+z+Diffsize*resolution,uyindex) = Linz(z,uxindex+Diffsize,uyindex);
            end

        end
    end
end
end

if Dox == 1
Lxyx = ( Lxytempx(1:length(uxsize),:) +Lxytempx(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),:) )* L_coefficient ;

Ioutuxuyx = Iouttempx(1:length(uxsize),:);
Ioutuxuy2x = Iouttempx(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),:);
Ioutnormx = (Ioutuxuyx'+Ioutuxuy2x');

Koutuxuyx = Kouttempx(1:length(uxsize),:);
Koutuxuy2x = Kouttempx(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),:);
Koutnormx = (Koutuxuyx'+Koutuxuy2x');
end

if Doy == 1
Lxyy = ( Lxytempy(1:length(uxsize),:) +Lxytempy(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),:) )* L_coefficient ;

Ioutuxuyy = Iouttempy(1:length(uxsize),:);
Ioutuxuy2y = Iouttempy(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),:);
Ioutnormy = (Ioutuxuyy'+Ioutuxuy2y');

Koutuxuyy = Kouttempy(1:length(uxsize),:);
Koutuxuy2y = Kouttempy(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),:);
Koutnormy = (Koutuxuyy'+Koutuxuy2y');
end

if Doz == 1
Lxyz = ( Lxytempz(1:length(uxsize),:) +Lxytempz(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),:) )* L_coefficient ;

Ioutuxuyz = Iouttempz(1:length(uxsize),:);
Ioutuxuy2z = Iouttempz(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),:);
Ioutnormz = (Ioutuxuyz'+Ioutuxuy2z');

Koutuxuyz = Kouttempz(1:length(uxsize),:);
Koutuxuy2z = Kouttempz(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),:);
Koutnormz = (Koutuxuyz'+Koutuxuy2z');
end

if Dox == 1
LEEx = sum(sum(Koutnormx))/sum(sum(Lxyx));
Purcellx = sum(sum(Lxyx))*(uxsize(2)-uxsize(1))*(uysize(2)-uysize(1))*2/pi;
Transx = LEEx * Purcellx;
end

if Doy == 1
LEEy = sum(sum(Koutnormy))/sum(sum(Lxyy));
Purcelly = sum(sum(Lxyy))*(uxsize(2)-uxsize(1))*(uysize(2)-uysize(1))*2/pi;
Transy = LEEy * Purcelly;
end

if Doz == 1
LEEz = sum(sum(Koutnormz))/sum(sum(Lxyz));
Purcellz = sum(sum(Lxyz))*(uxsize(2)-uxsize(1))*(uysize(2)-uysize(1))*2/pi;
Transz = LEEz * Purcellz;
end

[a,b] = min(abs(uxsize-n_top/n_org));
[c,d] = min(abs(uxsize+n_top/n_org));

if Dox*Doy*Doz == 1

Trans = (Transx+Transy+Transz)/3;
PurcellA = (Purcellx+Purcelly+Purcellz)/3;
LEE = Trans / PurcellA;

Ksum = (Koutnormx+Koutnormy+Koutnormz)/3;
Isum = (Ioutnormx+Ioutnormy+Ioutnormz)/3;
Lsum = (Lxyx+Lxyy+Lxyz)/3;

c_DMM = sum(sum(Ksum))*n_org^2*(uxsize(2)-uxsize(1))*(uysize(2)-uysize(1));

Ifin = Isum * LEE / c_DMM;
end

if Dox == 1
result.xdip.LEE = LEEx;
result.xdip.Purcell = Purcellx;
result.xdip.T = Transx;
result.xdip.K = Koutnormx;
result.xdip.I = Ioutnormx;
result.xdip.L = Lxyx;
end

if Doy == 1
result.ydip.LEE = LEEy;
result.ydip.Purcell = Purcelly;
result.ydip.T = Transy;
result.ydip.K = Koutnormy;
result.ydip.I = Ioutnormy;
result.ydip.L = Lxyy;
end

if Doz == 1
result.zdip.LEE = LEEz;
result.zdip.Purcell = Purcellz;
result.zdip.T = Transz;
result.zdip.K = Koutnormz;
result.zdip.I = Ioutnormz;
result.zdip.L = Lxyz;
end

if Dox*Doy*Doz == 1
result.tdip.LEE = LEE;
result.tdip.Purcell = PurcellA;
result.tdip.T = Trans;
result.tdip.K = Ksum;
result.tdip.I = Isum;
result.tdip.L = Lsum;
result.Ifin = Ifin;

end

result.d = d;
result.b = b;

% save('DMM_datas.mat')
end