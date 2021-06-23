function [tm, rm, prob] = PDEM_solve(asgn_prob, d, v, t, dt, dt_ratio, thres)
%function [tm, rm, prob, pdf_rps] = PDEM_solve(asgn_prob, d, v, t, dt, dt_ratio, thres)

%PDEM_solve - Description
%solve the GDEE by finite difference method with TVD scheme

%*******************************************
%code by J. S. Yang
%date: 2019-11-05

%modified by J. S. Yang
%date: 2020-06-16
%interploate the velocity history inside this function

%modified by J. S. Yang
%date: 2020-09-22
%add an option to output the probability density function at each representative points
%at the end of time history if equivalent extreme value distribution approach is considered

%modified by J.S. Yang
%data: 2020-11-18
%vectorize the FDM with TVD scheme

%Reference:
%Chen JB, Yang JS, Jensen H. (2020) Structural optimization considering dynamic 
%   reliability constraints via probability density evolution method and change
%   of probability measure. Structural and Multidisciplinary Optimization 62(5): 2499-2516.

%*******************************************

%====================================================================
    %INUPT:
    %----------------------------
    %asgn_prob: assigned probabilities
    %d        : displacement history (each column according to a representative point)
    %v        : velocity history (each column according to a representative point)
    %t        : time instance series
    %dt       : time increasement
    %dt_ratio : ratio between the time increasement of GDEE solution and structural response
    %thres    : absorbing boundary condition

    %OUTPUT:
    %----------------------------
    %prob     : matrix of time variant probability density
    %tm       : mesh of time
    %rm       : mesh of displacement
    %pdf_rps  : PDF at each representative points at the end of time history
%=====================================================================
    vectorize = 'on';
    %number of representative points
    np    = length(asgn_prob);
    
    %---------------------------------------
    %2020-06-16
    %interpolated time series
    dt    = dt_ratio*dt;
    t_i     = (t(1):dt:t(end))';
    %---------------------------------------
    
    %number of time step
    nt    = length(t_i);
    d_max = max(max(abs(d)));   %maximum displacement
    v_max = max(max(abs(v)));   %maximum velocity
    
    %---------------------------------------
    %2020-06-16
    %To reduce memeory consumption, only keep the first column in memory
    d     = d(1,:);
    %---------------------------------------
    
    %range of the mesh
    f_r = 2;    %parameters for the range of the mesh
    if length(thres(:)) == 1
        % circle boundary, i.e., the upper threshold is equal to
        % the lower threshold in terms of absolute value
        if isinf(thres)
            thres0 = f_r*abs(d_max);
        else
            thres0 = thres;
        end
        bd_absorb(1) = -thres0;
        bd_absorb(2) =  thres0;
    elseif length(thres(:)) == 2 && thres(1) ~= thres(2)
        thres0 = thres;
        %two-sided boundary
        for k = 1:1:2
            if isinf(thres(k))
                thres0(k) = sign(thres(k))*f_r*abs(d_max);
            end
        end
        bd_absorb(1) = min(thres0);
        bd_absorb(2) = max(thres0);
    elseif length(thres(:)) == 2 && thres(1) == thres(2)
        thres0 = thres(1);
        if isinf(thres0)
            thres0 = f_r*abs(d_max);   
        end
        bd_absorb(1) = -thres0;
        bd_absorb(2) =  thres0;
    else
        n_in = length(thres);
        error('The dimension of <thres> must be 1 or 2 rather than %d!', n_in);
    end

    %factor for mesh adjustion
    f_msh_adj = 0.95;
    %mesh size in response
    dr        = dt/(f_msh_adj/abs(v_max));
    %number of mesh
    nr1       = abs(bd_absorb(1))/dr;
    nr2       = abs(bd_absorb(2))/dr;
    %adjusted number of mesh in response
    if bd_absorb(1) * bd_absorb(2) <= 0
        nr = ceil(nr1) + ceil(nr2);
    elseif bd_absorb(1) > 0 && bd_absorb(2) > 0
        nr = ceil(nr2) - floor(nr1);
    elseif bd_absorb(1) < 0 && bd_absorb(2) < 0
        nr = ceil(nr1) - floor(nr2);
    end
    %adjusted mesh size in response
    dr = (bd_absorb(2) - bd_absorb(1))/nr;
    %adjusted mesh ratio
    mesh_ratio = dt / dr;
    if bd_absorb(1) <= 0 && bd_absorb(2) >= 0
        grid_r1 = -(ceil(nr1):-1:0)' * dr;
        grid_r2 =  (1:1:ceil(nr2))' * dr;
        grid_r = [grid_r1; grid_r2];
    elseif bd_absorb(1) <= 0 && bd_absorb(2) < 0
        grid_r = -(ceil(nr1):-1:floor(nr2))' * dr;
    elseif bd_absorb(1) > 0 && bd_absorb(2) > 0
        grid_r =  (floor(nr1):1:ceil(nr2))' * dr;
    end
    nr = length(grid_r);
    %mesh in time
    grid_t = t_i;
    %grid of mesh
    [tm, rm] = meshgrid(grid_t, grid_r);
    % %information display
    % fprintf('The number of the representative points is %15d.\n', np);
    % fprintf('The maximum value of the response is       %15.6f.\n', d_max);
    % fprintf('The lower boundary of the mesh domain is   %15.6f.\n', bd_absorb(1));
    % fprintf('The upper boundary of the mesh domain is   %15.6f.\n', bd_absorb(2));
    % fprintf('The mesh size for response is              %15.6f.\n', dr);
    % fprintf('The mesh size for time is                  %15.6f.\n', dt);
    % fprintf('The mesh ratio is                          %15.6f.\n', mesh_ratio);
    % fprintf('The number of grid for response is         %15d.\n', length(grid_r));
    % fprintf('The number of grid for time is             %15d.\n', length(grid_t));

    prob    = zeros(size(tm));
    % pdf_rps = zeros(length(rm(:, end)), np);
    % for each of the representative point
    parfor nn = 1:np
        %---------------------------------------------
        %interpolate the time history of velocity
        vv0   = v(:,nn);
        vv    = interp1(t, vv0, t_i, 'linear');
        %---------------------------------------------
        % fprintf('Representative point:                  %15d/%6d.\n',nn,np);
        %initial condition
        prob0 = zeros(size(tm));
        %location of the initial response
        grid_r1 = grid_r;
        bd_absorb1 = bd_absorb;
        %---------------------------------------
        %2020-06-16
        %n_ini = round(abs(d(1,nn)-bd_absorb1(1))/dr);
        %---------------------------------------
        n_ini = round(abs(d(nn)-bd_absorb1(1))/dr);
        prob0(n_ini+1, 1) = 1/dr;
        %---------------------------------------
        %2020-06-16
        %if d(1,nn) > bd_absorb1(2) || d(1,nn) < bd_absorb1(1)
        if d(nn) > bd_absorb1(2) || d(nn) < bd_absorb1(1)
            error('The initial response is out of the mesh domain!');
        end
        %---------------------------------------
        %for each time step
        for kk = 1:1:nt-1
            % the velocity at current time step
            %a = vv(kk,nn);
            a = vv(kk);
            %for each response grid
            if strcmpi(vectorize, 'off')
                for jj = 1:1:nr
                    %flux limiter
                    %-------------------------------------------------------------------------
                    %r_j+1/2^+, r3
                    if jj == nr
                        r11 = 1;
                    elseif jj == nr-1
                        if prob0(jj+1, kk) == prob0(jj, kk)
                            r11 = 1;
                        else
                            r11 = -prob0(jj+1, kk)/(prob0(jj+1, kk)-prob0(jj, kk));
                        end
                    else
                        if prob0(jj+1, kk) == prob0(jj, kk)
                            r11 = 1;
                        else
                            r11 = (prob0(jj+2, kk)-prob0(jj+1, kk))/(prob0(jj+1, kk)-prob0(jj, kk));
                        end
                    end
                    %-------------------------------------------------------------------------
                    %r_j+1/2^-, r1
                    if jj == nr
                        if prob0(jj, kk) == 0
                            r12 = 1;
                        else
                            r12 = (prob0(jj, kk)-prob0(jj-1, kk))/(-prob0(jj, kk));
                        end
                    elseif jj == 1
                        if prob0(jj+1, kk) == prob0(jj, kk)
                            r12 = 1;
                        else
                            r12 = prob0(jj, kk)/(prob0(jj+1, kk)-prob0(jj, kk));
                        end
                    else
                        if prob0(jj+1, kk) == prob0(jj, kk)
                            r12 = 1;
                        else
                            r12 = (prob0(jj, kk)-prob0(jj-1, kk))/(prob0(jj+1, kk)-prob0(jj, kk));
                        end
                    end
                    %-------------------------------------------------------------------------
                    %r_j-1/2^+, r4
                    if jj == nr
                        if prob0(jj, kk) == prob0(jj-1, kk)
                            r21 = 1;
                        else
                            r21 = (-prob0(jj, kk))/(prob0(jj, kk)-prob0(jj-1, kk));
                        end
                    elseif jj == 1
                        if prob0(jj, kk) == 0
                            r21 = 1;
                        else
                            r21 = (prob0(jj+1, kk)-prob0(jj, kk))/(prob0(jj, kk));
                        end
                    else
                        if prob0(jj, kk) == prob0(jj-1, kk)
                            r21 = 1;
                        else
                            r21 = (prob0(jj+1, kk)-prob0(jj, kk))/(prob0(jj, kk)-prob0(jj-1, kk));
                        end
                    end
                    %-------------------------------------------------------------------------
                    %r_j-1/2^-, r2
                    if jj == 1
                        r22 = 1;
                    elseif jj == 2
                        if prob0(jj, kk) == prob0(jj-1, kk)
                            r22 = 1;
                        else
                            r22 = prob0(jj-1, kk)/(prob0(jj, kk)-prob0(jj-1, kk));
                        end
                    else
                        if prob0(jj, kk) == prob0(jj-1, kk)
                            r22 = 1;
                        else
                            r22 = (prob0(jj-1, kk)-prob0(jj-2, kk))/(prob0(jj, kk)-prob0(jj-1, kk));
                        end
                    end
                    %phi
                    phi11    = max(max(0,min(2*r11,1)),min(r11,2)); %r3
                    phi12    = max(max(0,min(2*r12,1)),min(r12,2)); %r1
                    phi21    = max(max(0,min(2*r21,1)),min(r21,2)); %r4
                    phi22    = max(max(0,min(2*r22,1)),min(r22,2)); %r2
                    phi1     = heaviside(-a)*phi11 + heaviside(a)*phi12;
                    phi2     = heaviside(-a)*phi21 + heaviside(a)*phi22;
                    %finite difference update with tvd scheme
                    if jj == 1
                        d_prob_1 =  prob0(jj+1, kk) - prob0(jj, kk);
                        d_prob_2 =  prob0(jj, kk);
                    elseif jj == nr
                        d_prob_1 = -prob0(jj, kk);
                        d_prob_2 =  prob0(jj, kk)   - prob0(jj-1, kk);
                    else
                        d_prob_1 =  prob0(jj+1, kk) - prob0(jj, kk);
                        d_prob_2 =  prob0(jj, kk)   - prob0(jj-1, kk);
                    end
                    tmp      = mesh_ratio*a;
                    %check the CFL conditon
                    if abs(tmp) > 1
                        warning('The CFL condition is not verified at current time step!')
                    end
                    prob0(jj, kk+1) = prob0(jj, kk) - 1/2*(tmp-abs(tmp))*d_prob_1...
                                                    - 1/2*(tmp+abs(tmp))*d_prob_2...
                                                    - 1/2*(abs(tmp)-tmp^2)*(phi1*d_prob_1-phi2*d_prob_2);
                    if grid_r1(jj)>bd_absorb1(2) || grid_r1(jj)<bd_absorb1(1)
                        prob0(jj, kk+1) = 0;
                    end
                end
            elseif strcmpi(vectorize, 'on')
                probk      = prob0(:, kk);

                %flux limter
                dp10       = probk(3:end) - probk(2:end-1);
                dp20       = probk(2:end-1) - probk(1:end-2);
                abnormInd1 = (dp10 == 0);
                abnormInd2 = (dp20 == 0);

                dp1 = dp10; dp2 = dp20;
                r11v            = zeros(nr,1);
                r11v(end-1:end) = 1;
                dp1(abnormInd2) = 1;
                dp2(abnormInd2) = 1;
                r11v(1:end-2)   = dp1./dp2;

                dp1 = dp10; dp2 = dp20;
                r12v            = zeros(nr,1);
                r12v(end) = 1; r12v(1) = 1;
                dp1(abnormInd1) = 1;
                dp2(abnormInd1) = 1;
                r12v(2:end-1)   = dp2./dp1;

                dp1 = dp10; dp2 = dp20;
                r21v            = zeros(nr,1);
                r21v(end) = 1; r21v(1) = 1;
                dp1(abnormInd2) = 1;
                dp2(abnormInd2) = 1;
                r21v(2:end-1)   = dp1./dp2;

                dp1 = dp10; dp2 = dp20;
                r22v            = zeros(nr,1);
                r22v(1:2)       = 1;
                dp1(abnormInd1) = 1;
                dp2(abnormInd1) = 1;
                r22v(3:end)     = dp2./dp1;

                %phi
                zerosv   = zeros(nr,1);
                onesv    = ones(nr,1);
                twosv    = 2*onesv;
                phi11    = max(max(zerosv,min(2*r11v,onesv)),min(r11v,twosv)); %r3
                phi12    = max(max(zerosv,min(2*r12v,onesv)),min(r12v,twosv)); %r1
                phi21    = max(max(zerosv,min(2*r21v,onesv)),min(r21v,twosv)); %r4
                phi22    = max(max(zerosv,min(2*r22v,onesv)),min(r22v,twosv)); %r2
                phi1     = heaviside(-a)*phi11 + heaviside(a)*phi12;
                phi2     = heaviside(-a)*phi21 + heaviside(a)*phi22;

                %FDM with TVD
                ddp1          = zeros(nr, 1);
                ddp2          = zeros(nr, 1);
                ddp1(end)     = -probk(end);
                ddp2(1)       = probk(1);
                ddp1(1:end-1) = probk(2:end)-probk(1:end-1);
                ddp2(2:end)   = probk(2:end)-probk(1:end-1);

                tmp      = mesh_ratio*a;
                %check the CFL conditon
                if abs(tmp) > 1
                    warning('The CFL condition is not verified at current time step!')
                end
                probk1 = probk  - 1/2*(tmp-abs(tmp))*ddp1...
                                - 1/2*(tmp+abs(tmp))*ddp2...
                                - 1/2*(abs(tmp)-tmp^2)*(phi1.*ddp1-phi2.*ddp2);
                iiind1 = grid_r1>bd_absorb1(2);
                iiind2 = grid_r1<bd_absorb1(1);
                iiind  = iiind1 | iiind2;
                probk1(iiind) = 0;
                prob0(:, kk+1) = probk1;
            end
        end
        prob           = prob + asgn_prob(nn)*prob0;
        %pdf_rps(:, nn) = asgn_prob(nn)*prob0(:,end);
    end

end
