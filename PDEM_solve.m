function [tm, rm, prob] = PDEM_solve(asgn_prob, d, v, t, dt, thres)
%PDEM_solve - Description
%solve the GDEE by finite difference method with TVD scheme

%by J. S. Yang
%date: 2019-11-05
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
%=====================================================================

    %number of representative points
    np    = length(asgn_prob);
    %number of time step
    nt    = length(t);
    d_max = max(max(abs(d)));   %maximum displacement
    v_max = max(max(abs(v)));   %maximum velocity
    
    %range of the mesh
    f_r = 2;    %parameters for the range of the mesh
    if length(thres(:)) == 1
        % circle boundary, i.e., the upper threshold is equal to
        % the lower threshold in terms of absolute value
        if isinf(thres)
            thres0 = f_r*abs(d_max);
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
    %mesh in time
    grid_t = t;
    %grid of mesh
    [tm, rm] = meshgrid(grid_t, grid_r);
    %information display
    fprintf('The number of the representative points is %15d.\n', np);
    fprintf('The maximum value of the response is       %15.6f.\n', d_max);
    fprintf('The lower boundary of the mesh domain is   %15.6f.\n', bd_absorb(1));
    fprintf('The upper boundary of the mesh domain is   %15.6f.\n', bd_absorb(2));
    fprintf('The mesh size for response is              %15.6f.\n', dr);
    fprintf('The mesh size for time is                  %15.6f.\n', dt);
    fprintf('The mesh ratio is                          %15.6f.\n', mesh_ratio);
    fprintf('The number of grid for response is         %15d.\n', length(grid_r));
    fprintf('The number of grid for time is             %15d.\n', length(grid_t));

    prob = zeros(size(tm));
    % for each of the representative point
    parfor nn = 1:1:np
        fprintf('Representative point:                  %15d/%6d.\n',nn,np);
        %initial condition
        prob0 = zeros(size(tm));
        %location of the initial response
        grid_r1 = grid_r;
        bd_absorb1 = bd_absorb;
        n_ini = round(abs(d(1,nn)-bd_absorb1(1))/dr);
        prob0(n_ini+1, 1) = 1/dr;
        if d(1,nn) > bd_absorb1(2) || d(1,nn) < bd_absorb1(1)
            error('The initial response is out of the mesh domain!');
        end
        %for each time step
        for kk = 1:1:nt-1
            % the velocity at current time step
            a = v(kk,nn);
            %for each response grid
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
                    warining('The CFL condition is not verified at current time step!')
                end
                prob0(jj, kk+1) = prob0(jj, kk) - 1/2*(tmp-abs(tmp))*d_prob_1...
                                                - 1/2*(tmp+abs(tmp))*d_prob_2...
                                                - 1/2*(abs(tmp)-tmp^2)*(phi1*d_prob_1-phi2*d_prob_2);
                if grid_r1(jj)>bd_absorb1(2) || grid_r1(jj)<bd_absorb1(1)
                    prob0(jj, kk+1) = 0;
                end
            end
        end
        prob = prob + asgn_prob(nn)*prob0;
    end
end