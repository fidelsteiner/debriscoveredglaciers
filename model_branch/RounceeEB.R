for(t in 1:datasize){
n_iterations[t] <- 0;
Ts_past[t] <- 0;
Td[N,t] = 273.15;
# Initially assume Ts = Tair, for all other time steps assume it's equal to previous Ts
if{t == 1
Td[1,t] = Tair[t]}
else{Td[1,t] = Td[1,t-1]}
                                
# Calculate temperature profile in the debris
#For t = 0, which is i = 1, assume initial condition of linear temperature profile in the debris
                                if(t == 1){
                                    Td_gradient <- (Td[1,1] - Td[N,1])/d}
                                    for(j in 2:(N-1)){
                                        Td[j,1] <- Td(1,1) - (j*h)*Td_gradient}
                                    end
                                else
                                    % Perform Crank-Nicholson Scheme
                                        for j = 2:1:(N-1)
                                            % Equations A8 in Reid and Brock (2010) 
                                                a_Crank(j,i) = C;
                                                b_Crank(j,i) = 2*C+1;
                                                c_Crank(j,i) = C;

                                            % Equations A9 in Reid and Brock (2010) 
                                                if j == 2
                                                    d_Crank(j,i) = C*Td(1,i) + C*Td(1,i-1) + (1-2*C)*Td(j,i-1) + C*Td(j+1,i-1);
                                                elseif j < (N-1)
                                                    d_Crank(j,i) = C*Td(j-1,i-1) + (1-2*C)*Td(j,i-1) + C*Td(j+1,i-1);
                                                elseif j == (N-1)
                                                    d_Crank(j,i) = 2*C*Td(N,i) + C*Td(N-2,i-1) + (1-2*C)*Td(N-1,i-1);
                                                end
                                                    % note notation:
                                                        % "i-1" refers to the past
                                                        % "j-1" refers to the cell above it
                                                        % "j+1" refers to the cell below it          

                                            % Equations A10 and A11 in Reid and Brock (2010)
                                                if j == 2
                                                    A_Crank(j,i) = b_Crank(j,i);
                                                    S_Crank(j,i) = d_Crank(j,i);
                                                else
                                                    A_Crank(j,i) = b_Crank(j,i) - a_Crank(j,i)/A_Crank(j-1,i)*c_Crank(j-1,i);
                                                    S_Crank(j,i) = d_Crank(j,i) + a_Crank(j,i)/A_Crank(j-1,i)*S_Crank(j-1,i);
                                                end
                                        end

                                    % Equations A12 in Reid and Brock (2010)
                                        for j = N-1:-1:2
                                            if j == (N-1)
                                                Td(j,i) = S_Crank(j,i)/A_Crank(j,i);
                                            else
                                                Td(j,i) = 1/A_Crank(j,i)*(S_Crank(j,i)+c_Crank(j,i)*Td(j+1,i));
                                            end
                                        end
                                end
                                % Assume snow-free surface and compute fluxes normally
                                    % Calculate Surface Energy Fluxes
                                        if Rain_PyrStat(i) > 0
                                            eS_Saturated(i) = 611*exp(-Lv/R*(1/Td(1,i)-1/273.15));
                                            eS(i) = eS_Saturated(i);
                                            eZ(i) = RH_PyrStat(i)*eZ_Saturated(i);
                                            LE_Benn(i) = 0.622*density_air_0/P0*Lv*A_Benn*u_PyrStat(i)*(eZ(i)-eS(i));
                                        else
                                            eS_Saturated(i) = 0;
                                            LE_Benn(i) = 0;
                                        end
                                        Rn(i) = Sin(i)*(1-Albedo) + emissivity*(Lin_PyrStat(i) - (5.67*10^-8*Td(1,i)^4));
                                        H_Benn(i) = density_air_0*(P/P0)*cA*A_Benn*u_PyrStat(i)*(Tair(i)-Td(1,i));
                                        P_Flux(i) = density_water*cW*Rain_PyrStat(i)/(delta_t)*(Tair(i)-Td(1,i));
                                        Qc(i) = k*(Td(2,i) - Td(1,i))/h;
                                        F_Ts_Benn(i) = Rn(i) + H_Benn(i) + LE_Benn(i) + Qc(i) + P_Flux(i);
                                        dRn(i) = -4*emissivity*(5.67*10^-8)*Td(1,i)^3;
                                        dH_Benn(i) = -1*density_air_0*P/P0*cA*A_Benn*u_PyrStat(i);
                                        if Rain_PyrStat(i) > 0
                                            dLE_Benn(i) = -0.622*density_air_0/P0*Lv*A_Benn*u_PyrStat(i)*611*exp(-Lv/R*(1/Td(1,i)-1/273.15))*(Lv/R*Td(1,i)^-2);
                                        else
                                            dLE_Benn(i) = 0;
                                        end
                                        dP_Flux(i) = -density_water*cW*Rain_PyrStat(i)/(delta_t);
                                        dQc(i) = -k/h;
                                        dF_Ts_Benn(i) = dRn(i) + dH_Benn(i) + dLE_Benn(i) + dQc(i) + dP_Flux(i);
                                % Newton-Raphson method to solve for surface temperature
                                        while abs(Td(1,i)-Ts_past(i)) > 0.01 & n_iterations < 100
                                            n_iterations(i) = n_iterations(i) + 1;
                                            Ts_past(i) = Td(1,i);
                                            % max step size is 1 degree C
                                                Td(1,i) = Ts_past(i) - F_Ts_Benn(i)/dF_Ts_Benn(i);
                                                    if (Td(1,i) - Ts_past(i)) > 1
                                                        Td(1,i) = Ts_past(i) + 1;
                                                    elseif (Td(1,i) - Ts_past(i)) < -1
                                                        Td(1,i) = Ts_past(i) - 1;
                                                    end
                                     % Calculate temperature profile in the debris
                                            % For t = 0, which is i = 1, assume initial condition of linear temperature profile in the debris
                                                if i == 1
                                                    Td_gradient = (Td(1,1) - Td(N,1))/debris_thickness;
                                                    for j = 2:1:(N-1)
                                                        Td(j,1) = Td(1,1) - (j*h)*Td_gradient;
                                                    end
                                                else
                                                    % Perform Crank-Nicholson Scheme
                                                        for j = 2:1:(N-1)
                                                            % Equations A8 in Reid and Brock (2010) 
                                                                    a_Crank(j,i) = C;
                                                                    b_Crank(j,i) = 2*C+1;
                                                                    c_Crank(j,i) = C;

                                                            % Equations A9 in Reid and Brock (2010) 
                                                                    if j == 2
                                                                        d_Crank(j,i) = C*Td(1,i) + C*Td(1,i-1) + (1-2*C)*Td(j,i-1) + C*Td(j+1,i-1);
                                                                    elseif j < (N-1)
                                                                        d_Crank(j,i) = C*Td(j-1,i-1) + (1-2*C)*Td(j,i-1) + C*Td(j+1,i-1);
                                                                    elseif j == (N-1)
                                                                        d_Crank(j,i) = 2*C*Td(N,i) + C*Td(N-2,i-1) + (1-2*C)*Td(N-1,i-1);
                                                                    end
                                                                        % note notation:
                                                                            % "i-1" refers to the past
                                                                            % "j-1" refers to the cell above it
                                                                            % "j+1" refers to the cell below it          

                                                            % Equations A10 and A11 in Reid and Brock (2010)
                                                                    if j == 2
                                                                        A_Crank(j,i) = b_Crank(j,i);
                                                                        S_Crank(j,i) = d_Crank(j,i);
                                                                    else
                                                                        A_Crank(j,i) = b_Crank(j,i) - a_Crank(j,i)/A_Crank(j-1,i)*c_Crank(j-1,i);
                                                                        S_Crank(j,i) = d_Crank(j,i) + a_Crank(j,i)/A_Crank(j-1,i)*S_Crank(j-1,i);
                                                                    end
                                                        end

                                                            % Equations A12 in Reid and Brock (2010)
                                                                    for j = N-1:-1:2
                                                                        if j == (N-1)
                                                                            Td(j,i) = S_Crank(j,i)/A_Crank(j,i);
                                                                        else
                                                                            Td(j,i) = 1/A_Crank(j,i)*(S_Crank(j,i)+c_Crank(j,i)*Td(j+1,i));
                                                                        end
                                                                    end
                                                end
                                            % Assume snow-free surface and compute fluxes normally
                                                    % Calculate Surface Energy Fluxes
                                                    if Rain_PyrStat(i) > 0
                                                        eS_Saturated(i) = 611*exp(-Lv/R*(1/Td(1,i)-1/273.15));
                                                        eS(i) = eS_Saturated(i);
                                                        eZ(i) = RH_PyrStat(i)*eZ_Saturated(i);
                                                        LE_Benn(i) = 0.622*density_air_0/P0*Lv*A_Benn*u_PyrStat(i)*(eZ(i)-eS(i));
                                                    else
                                                        eS_Saturated(i) = 0;
                                                        LE_Benn(i) = 0;
                                                    end
                                                    Rn(i) = Sin(i)*(1-Albedo) + emissivity*(Lin_PyrStat(i) - (5.67*10^-8*Td(1,i)^4));
                                                    H_Benn(i) = density_air_0*(P/P0)*cA*A_Benn*u_PyrStat(i)*(Tair(i)-Td(1,i));
                                                    P_Flux(i) = density_water*cW*Rain_PyrStat(i)/(delta_t)*(Tair(i)-Td(1,i));
                                                    Qc(i) = k*(Td(2,i) - Td(1,i))/h;
                                                    F_Ts_Benn(i) = Rn(i) + H_Benn(i) + LE_Benn(i) + Qc(i) + P_Flux(i);
                                                    dRn(i) = -4*emissivity*(5.67*10^-8)*Td(1,i)^3;
                                                    dH_Benn(i) = -1*density_air_0*P/P0*cA*A_Benn*u_PyrStat(i);
                                                    if Rain_PyrStat(i) > 0
                                                        dLE_Benn(i) = -0.622*density_air_0/P0*Lv*A_Benn*u_PyrStat(i)*611*exp(-Lv/R*(1/Td(1,i)-1/273.15))*(Lv/R*Td(1,i)^-2);
                                                    else
                                                        dLE_Benn(i) = 0;
                                                    end
                                                    dP_Flux(i) = -density_water*cW*Rain_PyrStat(i)/(delta_t);
                                                    dQc(i) = -k/h;

                                                    dF_Ts_Benn(i) = dRn(i) + dH_Benn(i) + dLE_Benn(i) + dQc(i) + dP_Flux(i);
                                                if n_iterations == 100
                                                    Td(1,i) = (Td(1,i) + Ts_past(i)) / 2;
                                                end
                                        end
                                        Qc_ice(i) = k*(Td(N-1,i) - Td(N,i))/h;
                                        if Qc_ice(i) < 0
                                            Qc_ice(i) = 0;
                                        end
                                        Melt(i) = Qc_ice(i)*delta_t / (density_ice*Lf);
                                            % meters of ice melt
                        end % end Crank-Nicholson Newton-Raphson method for loop