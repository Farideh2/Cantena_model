clc
close all
clear

epsilon = (8.854.*10^-12);
mu = 4.*pi.*10^-7;
m_rara = input("Ingresa el valor de m: ");
n_rara = input("Ingresa el valor de n: ");
A =  input("Ingresa el valor de a(cm): ")*10^-2;
B = input("Ingresa el valor de b(cm): ")*10^-2;
salto = 10;

while true
    modo = int16(input("Modo Tm(1) o modo Te(2): "));
    if modo == 1 || modo == 2
        break
    elseif modo == 3
        break
    else
        disp("Porfavor escoga un modo valido")
    end
end

h2 = (((m_rara*pi)/A)^2)+(((n_rara*pi)/B)^2);
Fc = (1/(2*pi.*sqrt(mu.*epsilon))).*sqrt(h2).*10^-9;

disp("La frecuencia de corte es " + Fc)

while true
    frecuencia = (input("Ingrese una frecuencia: "));
    if frecuencia <= Fc
        disp("La frecuencia tiene que ser mayor a " + Fc)
    else
        break
    end
end

omega = frecuencia*10^9*2*pi;
k2 = (omega^2)*mu*epsilon;
beta = sqrt(k2-h2);
gamma = 1i*beta;
Rg = (2*pi)/beta;

if modo == 1
    TM_mode(A,B,m_rara,n_rara,Rg,gamma,omega,h2,k2,salto)
elseif modo == 2
    TE_mode(A,B,m_rara,n_rara,Rg,gamma,omega,h2,k2,salto)
elseif modo == 3
    TE_mode(A,B,m_rara,n_rara,Rg,gamma,omega,h2,k2,salto)
    figure
    TM_mode(A,B,m_rara,n_rara,Rg,gamma,omega,h2,k2,salto)
end

function TM_mode = TM_mode(a,b,m,n,Rg,gamma,omega,h2,k2,salto)
    syms x y z
    
    E0 = 1;
    epsilon = (8.854.*10^-12);
    mu = 4.*pi.*10^-7;
    
    %cara frontal (x,y)
    Z = .00001;
    
    [xplot, yplot] = meshgrid(0:a/salto:a, 0:a/salto:b);
    Exs = real(-(gamma./h2).*((m.*pi)/a).*E0.*cos((((m.*pi)/a).*x)).*sin(((n.*pi)/b).*y).*exp(-gamma.*Z));
    Eys = real(-(gamma./h2).*((n.*pi)/b).*E0.*sin((((m.*pi)/a).*x)).*cos(((n.*pi)/b).*y).*exp(-gamma.*Z));
    Hxs = real(((1i.*omega.*epsilon)/h2).*((n.*pi)/b).*E0.*sin((((m.*pi)/a).*x)).*cos(((n.*pi)/b).*y).*exp(-gamma.*Z));
    Hys = real(-((1i.*omega.*epsilon)/h2).*((m.*pi)/a).*E0.*cos((((m.*pi)/a).*x)).*sin(((n.*pi)/b).*y).*exp(-gamma.*Z));
    
    Exs_evaluada= subs(Exs, {x,y}, {xplot,yplot});
    Eys_evaluada= subs(Eys, {x,y}, {xplot,yplot});
    Hxs_evaluada= subs(Hxs, {x,y}, {xplot,yplot});
    Hys_evaluada= subs(Hys, {x,y}, {xplot,yplot});
    
    E_mag = sqrt(Exs_evaluada.^2+Eys_evaluada.^2);
    H_mag = sqrt(Hxs_evaluada.^2+Hys_evaluada.^2);
    
    tiledlayout(2,2)
    nexttile
    contourf(xplot,yplot,E_mag,30)
    hold on
    quiver(xplot,yplot,Exs_evaluada,Eys_evaluada,"r")
    title("TM-electric field Front view")
    xlabel("x axis")
    ylabel("y axis")
    nexttile
    
    contourf(xplot,yplot,H_mag,30)
    hold on
    quiver(xplot,yplot,Hxs_evaluada,Hys_evaluada,"r")
    title("TM-magnetic field Front view")
    xlabel("x axis")
    ylabel("y axis")
    nexttile([1 2])
    
    quiver(xplot,yplot,Hxs_evaluada,Hys_evaluada)
    hold on
    quiver(xplot,yplot,Exs_evaluada,Eys_evaluada)
    title("TM-MAgnetic & Electric vector field front view")
    legend("Magnetic field", "Electric field")
    figure
    
    %cara lateral (y,z)
    X = .0001;
    
    [yplot, zplot] = meshgrid(0:a/salto:b, 0:a/salto:Rg);
    Eys = real(-(gamma/h2).*((n.*pi)/b).*E0.*sin((((m.*pi)/a).*X)).*cos(((n.*pi)/b).*y).*exp(-gamma.*z));
    Ezs = real(E0*sin(((m*pi*X)/a))*sin(((n*pi*y)/a)).*exp(-gamma.*z));
    Hys = real(-((1i.*omega.*epsilon)/h2).*((m.*pi)/a).*E0.*cos((((m.*pi)/a).*X)).*sin(((n.*pi)/b).*y).*exp(-gamma.*z));
    
    Eys_evaluada= subs(Eys, {y,z}, {yplot,zplot});
    Ezs_evaluada= subs(Ezs, {y,z}, {yplot,zplot});
    Hys_evaluada= subs(Hys, {y,z}, {yplot,zplot});
    
    E_mag = sqrt(Eys_evaluada.^2+Ezs_evaluada.^2);
    H_mag = sqrt(Hys_evaluada.^2);
    [q,p] = size(Hys_evaluada);
    Hzs = zeros(q,p);
    
    tiledlayout(2,2)
    nexttile
    contourf(zplot,yplot,E_mag,15)
    hold on
    quiver(zplot,yplot,Ezs_evaluada,Eys_evaluada,"r")
    title("TM-electric field Side view")
    xlabel("z axis")
    ylabel("y axis")
    nexttile
    
    contourf(zplot,yplot,H_mag,30)
    hold on
    quiver(zplot,yplot,Hzs,Hys_evaluada,"r")
    title("TM-magnetic field side view")
    xlabel("z axis")
    ylabel("y axis")
    nexttile([1 2])
    
    quiver(zplot,yplot,Hzs,Hys_evaluada)
    hold on
    quiver(zplot,yplot,Ezs_evaluada,Eys_evaluada)
    title("TM-Magnetic & Electric vector field side view")
    legend("Magnetic field", "Electric field")
    hold off
    figure
    
    %cara top (x,z)
    Y = .00001;
    [zplot, xplot] = meshgrid(0:a/salto:Rg, 0:a/salto:a);
    Exs = real((-gamma./h2).*((m.*pi)/a).*E0.*cos((((m.*pi)/a).*x)).*sin(((n.*pi)/b).*Y).*exp(-gamma.*z));
    Ezs = real(E0*sin(((m*pi*x)/a))*sin(((n*pi*Y)/a)).*exp(-gamma.*z));
    Hxs = real(((1i.*omega.*epsilon)/h2).*((n.*pi)/b).*E0.*sin((((m.*pi)/a).*x)).*cos(((n.*pi)/b).*Y).*exp(-gamma.*z));
    
    Exs_evaluada= subs(Exs, {x,z}, {xplot,zplot});
    Ezs_evaluada= subs(Ezs, {x,z}, {xplot,zplot});
    Hxs_evaluada= subs(Hxs, {x,z}, {xplot,zplot});
    [q,p] = size(Hxs_evaluada);
    Hzs = zeros(q,p);
    
    E_mag = sqrt(Exs_evaluada.^2+Ezs_evaluada.^2);
    H_mag = sqrt(Hxs_evaluada.^2+Hzs.^2);
    
    tiledlayout(2,2)
    nexttile
    contourf(xplot,zplot,E_mag,30)
    hold on
    quiver(xplot,zplot,Exs_evaluada,Ezs_evaluada, "r")
    title("TM-electric field top view")
    xlabel("x axis")
    ylabel("z axis")
    nexttile
    
    contourf(xplot,zplot,H_mag,30)
    hold on
    quiver(xplot,zplot,Hxs_evaluada,Hzs, "r")
    title("TM-magnetic field top view")
    xlabel("x axis")
    ylabel("z axis")
    nexttile([1 2])
    
    quiver(xplot,zplot,Hxs_evaluada,Hzs)
    hold on
    quiver(xplot,zplot,Exs_evaluada,Ezs_evaluada)
    title("TM-Magnetic & Electric vector field top view")
    legend("Magnetic field", "Electric field")
    hold off

end

function TE_mode = TE_mode(a,b,m,n,Rg,gamma,omega,h2,k2,salto)
    syms x y z

    %z = 0:.1:Rg;
    H0 = 1;
    epsilon = (8.854.*10^-12);
    mu = 4.*pi.*10^-7;
    
    %cara frontal (x,y)
    Z = .00001;
    
    [xplot, yplot] = meshgrid(0:a/salto:a, 0:b/salto:b);
    Hxs = real((gamma./h2).*((m.*pi)/a).*H0.*sin((((m.*pi)/a).*x)).*cos(((n.*pi)/b).*y).*exp(-gamma.*Z));
    Hys = real((gamma./h2).*((n.*pi)/b).*H0.*cos((((m.*pi)/a).*x)).*sin(((n.*pi)/b).*y).*exp(-gamma.*Z));
    Exs = real(((1i.*omega.*mu)/h2).*((n.*pi)/b).*H0.*cos((((m.*pi)/a).*x)).*sin(((n.*pi)/b).*y).*exp(-gamma.*Z));
    Eys = real(((-1i.*omega.*mu)/h2).*((m.*pi)/a).*H0.*sin((((m.*pi)/a).*x)).*cos(((n.*pi)/b).*y).*exp(-gamma.*Z));
    
    Exs_evaluada= subs(Exs, {x,y}, {xplot,yplot});
    Eys_evaluada= subs(Eys, {x,y}, {xplot,yplot});
    Hxs_evaluada= subs(Hxs, {x,y}, {xplot,yplot});
    Hys_evaluada= subs(Hys, {x,y}, {xplot,yplot});
    
    E_mag = sqrt(Exs_evaluada.^2+Eys_evaluada.^2);
    H_mag = sqrt(Hxs_evaluada.^2+Hys_evaluada.^2);
    
    tiledlayout(2,2)
    nexttile
    contourf(xplot,yplot,E_mag,30)
    hold on
    quiver(xplot,yplot,Exs_evaluada,Eys_evaluada,"r")
    title("TE-Electric field Front view")
    xlabel("x axis")
    ylabel("y axis")
    nexttile
    
    contourf(xplot,yplot,H_mag,30)
    hold on
    quiver(xplot,yplot,Hxs_evaluada,Hys_evaluada, "r")
    title("TE-Magnetic field Front view")
    xlabel("x axis")
    ylabel("y axis")
    nexttile([1 2])
    
    quiver(xplot,yplot,Hxs_evaluada,Hys_evaluada)
    hold on
    quiver(xplot,yplot,Exs_evaluada,Eys_evaluada)
    title("TE-Magnetic & Electric vector field front view")
    legend("Magnetic field", "Electric field")
    figure
    
    %cara lateral (y,z)
    X = .0001;
    
    [yplot, zplot] = meshgrid(0:a/salto:b, 0:a/salto:Rg);
    Hys = real((gamma./h2).*((n.*pi)/b).*H0.*cos((((m.*pi)/a).*X)).*sin(((n.*pi)/b).*y).*exp(-gamma.*z));
    Hzs = real(H0.*cos((m*pi*X)/a).*cos((n*pi*y)/b).*exp(-gamma.*z));
    Eys = real(((-1i.*omega.*mu)/h2).*((m.*pi)/a).*H0.*sin((((m.*pi)/a).*X)).*cos(((n.*pi)/b).*y).*exp(-gamma.*z));
    
    Hys_evaluada= subs(Hys, {y,z}, {yplot,zplot});
    Hzs_evaluada= subs(Hzs, {y,z}, {yplot,zplot});
    Eys_evaluada= subs(Eys, {y,z}, {yplot,zplot});
    [q,p] = size(Hys_evaluada);
    Ezs = zeros(q,p);
    
    H_mag = sqrt(Hys_evaluada.^2+Hzs_evaluada.^2);
    E_mag = sqrt(Eys_evaluada.^2);
    
    tiledlayout(2,2)
    nexttile
    contourf(zplot,yplot,E_mag,15)
    hold on
    quiver(zplot,yplot,Ezs,Eys_evaluada,"r")
    title("TE-electric field Side(left view) view")
    xlabel("z axis")
    ylabel("y axis")
    nexttile
    
    contourf(zplot,yplot,H_mag,30)
    hold on
    quiver(zplot,yplot,Hzs_evaluada,Hys_evaluada,"r")
    title("TE-magnetic field side(left view) view")
    xlabel("z axis")
    ylabel("y axis")
    nexttile([1 2])
    
    quiver(zplot,yplot,Hzs_evaluada,Hys_evaluada)
    hold on
    quiver(zplot,yplot,Ezs,Eys_evaluada)
    title("TE-MAgnetic & Electric vector field side view")
    legend("Magnetic field", "Electric field")
    hold off
    figure
    
    %cara top (x,z)
    Y = .00001;
    [zplot, xplot] = meshgrid(0:a/salto:Rg, 0:a/salto:a);
    Hxs = real((gamma./h2).*((m.*pi)/a).*H0.*sin((((m.*pi)/a).*x)).*cos(((n.*pi)/b).*Y).*exp(-gamma.*z));
    Hzs = real(H0.*cos((m*pi*x)/a).*cos((n*pi*Y)/b).*exp(-gamma.*z));
    Exs = real(((1i.*omega.*mu)/h2).*((n.*pi)/b).*H0.*cos((((m.*pi)/a).*x)).*sin(((n.*pi)/b).*Y).*exp(-gamma.*z));
    
    Hxs_evaluada= subs(Hxs, {x,z}, {xplot,zplot});
    Hzs_evaluada= subs(Hzs, {x,z}, {xplot,zplot});
    Exs_evaluada= subs(Exs, {x,z}, {xplot,zplot});
    [q,p] = size(Hxs_evaluada);
    Ezs = zeros(q,p);
    
    H_mag = sqrt(Hxs_evaluada.^2+Hzs_evaluada.^2);
    E_mag = sqrt(Exs_evaluada.^2+Ezs.^2);
    
    tiledlayout(2,2)
    nexttile
    contourf(xplot,zplot,E_mag,30)
    hold on
    quiver(xplot,zplot,Exs_evaluada,Ezs,"r")
    title("TE-electric field bottom view")
    xlabel("x axis")
    ylabel("z axis")
    nexttile
    
    contourf(xplot,zplot,H_mag,30)
    hold on
    quiver(xplot,zplot,Hxs_evaluada,Hzs_evaluada,"r")
    title("TE-magnetic field bottom view")
    xlabel("x axis")
    ylabel("z axis")
    nexttile([1 2])
    
    quiver(xplot,zplot,Hxs_evaluada,Hzs_evaluada)
    hold on
    quiver(xplot,zplot,Exs_evaluada,Ezs)
    title("TE-MAgnetic & Electric vector field bottom view")
    legend("Magnetic field", "Electric field")
    hold off

end
