clear;
clc;

//geomertical properties
L1 = 1; L2 = 1;
imax = 42;  jmax = 42;
dx = L1/(imax-2); dy = L2/(jmax-2);
Dx = [dx/2 dx*ones(1,imax-2) dx/2]; Dy = [dy/2 dy*ones(1,jmax-2) dy/2];
del_x = [dx/2 dx*ones(1,imax-3) dx/2]; del_y = [dy/2 dy*ones(1,jmax-3) dy/2];  
//velocity conditons
u0 = 0; v0 = 0; P0 = 0;
uc = 1;
Un = 1; //top surface

// thermophysical properties
Re = 100;   rho = 1;  mu = (rho*uc*L1)/Re;    nu = mu/rho;

// IC and BC
U(1:imax,jmax) = Un;
U(1:imax,1:jmax-1) = u0/uc;
V(1:imax,1:jmax) = v0/uc;
P(1:imax,1:jmax) = P0;

//convergence tolerances
epsilon_st = 1e-3;  epsilon = 1e-8;
unsteadiness_nd = 1;
n = 0;

U_star = U;  V_star = V;
//Quick scheme

w = [3/8,6/8,-1/8];   //[D,U,UU] for QUICK Scheme
function u_q = quick(w,ue,up,uuu)
    u_q = w(1)*uuu + w(2)*up + w(3)*ue;
endfunction
// Main function starts, outer loop
while unsteadiness_nd >= epsilon_st
    n = n+1;    //increasing count of time steps taken
    U_old = U; V_old = V; P_old = P;   //storing old iteration velocities and pressure
    
    //time stamp 
    dt_adv = min(0.99./(abs(U)/dx + abs(V)/dy));
    dt_diff = 0.99*0.5*(1/dx^2 + 1/dy^2)^(-1)/mu;
    dt = min(dt_adv, dt_diff);
    
    // u-CV ,u-velocity and mass flux prediction 
    for j = 2:jmax-1
        for i = 1:imax-2
            m_ux_old(i,j) = 0.5*rho*(U_old(i+1,j) + U_old(i,j));  
            if (i==1)
                uq=U_old(i,j);
            else
                uq=quick(w,U_old(i+1,j),U_old(i,j),U_old(i-1,j));   //quick scheme
            end
            
            
            
            a_ux_old(i,j) = max(m_ux_old(i,j),0)*uq ;
            d_ux_old(i,j) = mu*(U_old(i+1,j) - U_old(i,j))/dx;
        end
    end
    
    for i = 2:imax-2
        for j = 1:jmax-1
            m_uy_old(i,j) = 0.5*rho*(V_old(i+1,j) + V_old(i,j));
            if (j==1)
                uq=U_old(i,j);
            else
                [uq]=quick(w,U_old(i,j+1),U_old(i,j),U_old(i,j-1));   //quick scheme
            end
            
            a_uy_old(i,j) = max(m_uy_old(i,j),0)*uq ;
            if (j == 1)|| (j == jmax-1) then
                d_uy_old(i,j) = mu*(U_old(i,j+1) - U_old(i,j))/(dy/2);
            else
                d_uy_old(i,j) = mu*(U_old(i,j+1) - U_old(i,j))/(dy);
            end
        end
    end
    
    for i = 2:imax-2
        for j = 2:jmax-1
            _A_ux_old(i,j) = (a_ux_old(i,j) - a_ux_old(i-1,j))*dy + (a_uy_old(i,j) - a_uy_old(i,j-1))*dx;
            _D_ux_old(i,j) = (d_ux_old(i,j) - d_ux_old(i-1,j))*dy + (d_uy_old(i,j) - d_uy_old(i,j-1))*dx;
            S_u_old(i,j) = (P_old(i,j) - P_old(i+1,j))*dy;
            U_star(i,j) = U_old(i,j) + (dt/(rho*dx*dy))*(_D_ux_old(i,j) - _A_ux_old(i,j) + S_u_old(i,j));
            m_star_x(i,j) = rho*U_star(i,j);
        end
    end
    
    m_star_x(1,2:imax-1) = rho*U_old(1,2:imax-1);
    m_star_x(imax-1,2:imax-1) = rho*U_old(imax-1,2:imax-1);
    
    //v-CV ,v-velocity and mass flux prediction 
    for i = 1:imax-1
        for j = 2:jmax-2
            m_vx_old(i,j) = 0.5*rho*(U_old(i,j+1) + U_old(i,j));
            if (i==1)
                vq=V_old(i,j);
            else
                [vq]=quick(w,V_old(i+1,j),V_old(i+1,j),V_old(i-1,j));  //quick scheme
            end
            
            a_vx_old(i,j) = max(m_vx_old(i,j),0)*vq;
            if (i == 1)|| (i == imax-1) then
                d_vx_old(i,j) = mu*(V_old(i,j+1) - V_old(i,j))/(dx/2);
            else
                d_vx_old(i,j) = mu*(V_old(i,j+1) - V_old(i,j))/(dx);
            end
        end
    end
    
    for i = 2:imax-1
        for j = 1:jmax-2
            m_vy_old(i,j) = 0.5*rho*(V_old(i,j+1) + V_old(i,j));
            if (j==1)
                vq=V_old(i,j);
            else
                [vq]=quick(w,V_old(i,j+1),V_old(i,j),V_old(i,j-1));
            end
            
            a_vy_old(i,j) = max(m_vy_old(i,j),0)*V_old(i,j) +  min(m_vy_old(i,j),0)*V_old(i,j+1);
            d_vy_old(i,j) = mu*(V_old(i,j+1) - V_old(i,j))/(dy);
        end
    end
    
    for i = 2:imax-1
        for j = 2:jmax-2
            _A_vx_old(i,j) = (a_vx_old(i,j) - a_vx_old(i-1,j))*dy + (a_vy_old(i,j) - a_vy_old(i,j-1))*dx;
            _D_vx_old(i,j) = (d_vx_old(i,j) - d_vx_old(i-1,j))*dy + (d_vy_old(i,j) - d_vy_old(i,j-1))*dx;
            S_v_old(i,j) = (P_old(i,j) - P_old(i,j+1))*dx;
            V_star(i,j) = V_old(i,j) + (dt/(rho*dx*dy))*(_D_vx_old(i,j) - _A_vx_old(i,j) + S_v_old(i,j));
            m_star_y(i,j) = rho*V_star(i,j);
        end
    end
    
    m_star_y(2:imax-1,1) = rho*V_old(2:imax-1,1);
    m_star_y(2:imax-1,jmax-1) = rho*V_old(2:imax-1,jmax-1);
    
    //predicted value of mass source
    S_star_m(2:imax-1,2:jmax-1) = (m_star_x(2:imax-1,2:jmax-1) - m_star_x(1:imax-2,2:jmax-1))*dy + (m_star_y(2:imax-1,2:jmax-1) - m_star_y(2:imax-1,1:jmax-2))*dx;
    
    P_corr = zeros(imax,jmax);
    
    N = 0;
    
    //inner loop
    while max(S_star_m) >= epsilon
        
        N = N+1;
        S_star_m_old = S_star_m;
        m_star_x_old = m_star_x;
        mstar_y_old = m_star_y;
        P_corr_old = P_corr;
        P_old = P;
        
        for i = 2:imax-1
            for j = 2:jmax-1
                Sm_corr(i,j) = -dt*(dy*(((P_corr_old(i+1,j) - P_corr_old(i,j))/del_x(i)) - ((P_corr_old(i,j) - P_corr(i-1,j))/del_x(i-1))) - dx*(((P_corr_old(i,j+1) - P_corr_old(i,j))/del_y(i)) - ((P_corr_old(i,j) - P_corr(i,j-1))/del_y(i-1))));
                ap(i,j) = dt*(Dy(i)*(1/del_x(i) + 1/del_x(i-1)) + Dx(i)*(1/del_y(j) + 1/del_y(j-1)));
                P_corr(i,j) = P_corr_old(i,j) - (S_star_m_old(i,j) + Sm_corr(i,j))/ap(i,j);
            end
        end
        
        for i = 1:imax-1
            for j = 2:jmax-1
                mx_corr(i,j) = -dt*((P_corr(i+1,j) - P_corr(i,j))/del_x(i));
                m_star_x(i,j) = m_star_x_old(i,j) + mx_corr(i,j);
                U_star(i,j) = m_star_x(i,j)/rho;
            end
        end
        
        for i = 1:imax-1
            for j = 2:jmax-1
                my_corr(i,j) = -dt*((P_corr(i,j+1) - P_corr(i,j))/del_y(i));
                m_star_y(i,j) = mstar_y_old(i,j) + my_corr(i,j);
                V_star(i,j) = m_star_y(i,j)/rho;
            end
        end
        
        for i = 2:imax-1
            for j = 2:imax-1
                S_star_m(i,j) = (m_star_x(i,j) - m_star_x(i-1,j))*dy + (m_star_y(i,j) - m_star_y(i,j-1))*dx;
                P(i,j) = P_old(i,j) + P_corr(i,j);
            end
        end
        U(2:imax-1,2:jmax-1) = U_star(2:imax-1,2:jmax-1);
        V(2:imax-1,2:jmax-1) = V_star(2:imax-1,2:jmax-1);
        m_star_x(2:imax-1,2:jmax-1) = m_star_x(2:imax-1,2:jmax-1);
        m_star_y(2:imax-1,2:jmax-1) = m_star_y(2:imax-1,2:jmax-1);

        P(1,1:jmax) = P(2,1:jmax);  P(imax,1:jmax) = P(imax-1,1:jmax);  P(1:jmax,imax) = P(1:jmax,imax-1);
        P_corr(1,1:jmax) = P_corr(2,1:jmax);  P_corr(imax,1:jmax) = P_corr(imax-1,1:jmax);  P_corr(1:jmax,imax) = P(1:jmax,imax-1);
    end
    
    //updating non dimensional unsteadiness
    dtau = uc*dt/L1;
    unsteadiness_nd = max(max(abs(U-U_old))/dtau,max(abs(V-V_old))/dtau);;
    disp(unsteadiness_nd);
end

disp(U,V,P);
x = [0 dx/2:dx:L1-dx/2 L1];
y = [0 dy/2:dy:L2-dy/2 L2];

scf(1)
champ(x,y,U,V,'k');
title('Velocity vectors');
xlabel('X axis')
ylabel('Y axis')

figure(2)
contourf(x,y,P);
title('Pressure contour overlapped with the streamlines');
xlabel('X axis')
ylabel('Y axis')
colorbar

