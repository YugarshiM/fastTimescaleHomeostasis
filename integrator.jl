function integrator(hash,tgt1,tgt2,tgt3,gf,gs,ICs,tEndSecs,reltolVal,chunkSizeSecs,tauz)

    print("Starting: ",hash, "\n")

    cSize = chunkSizeSecs * 1000.0;
    tEnd = tEndSecs * 1000.0;

    #####################################################################
    ### Parameters

    ## Sensor Parameters
    alphaOffset = .3;
    #alphaOffset = -1.0; # alpha always off
    #alphaOffset = 10.0; # alpha always on
    GF =  gf;
    GS =  gs;
    GD =  1.0;
    gamma = 10e-7;
    delta = 10e-5;
    sensorFuseDelta = .001;

    ## Targets
    FBar = tgt1
    SBar = tgt2
    DBar = tgt3

    ## Homeostatic Time Scales
    tauG = tauz[1];
    tauHalfac = tauz[2];
    tauS = tauz[3];
    tauAlpha = tauz[4];

    ## Reversal Potentials
    EK(t) = -80.0;
    EH = -20.0;
    ENa = 30.0;
    EL = -50.0;

    #####################################################################

    ### Define generic functions
    ## Generic Activation/Inactication Curves
    sigmd(x,translate,center,reactance) = 1/(1+exp(((x-translate) + center)/reactance));
    ## Define the time constant function for ion channels
    tauX(Volt, CT, DT, AT, BT, halfAc) = CT - DT/(1 + exp(((Volt-halfAc) + AT)/BT));
    ## Some intrinsic currents require a special time constant function
    spectau(Volt, CT, DT, AT, BT, AT2, BT2, halfAc) = CT + DT/(exp(((Volt-halfAc) + AT)/BT) + exp(((Volt-halfAc) + AT2)/BT2));
    ## Define the ionic currents; q is the exponent of the activation variable m
    iIonic(g, m, h, q, Volt, Erev) = g*(m^q)*h*(Volt - Erev);

    #####################################################################

    # Ion Channel Acitvation/Inactivation Curves
    NaMinf(V,trns) = sigmd(V,trns, 25.5, -5.29);  # m^3
    NaHinf(V,trns) = sigmd(V,trns, 48.9, 5.18);  # h
    CaTMinf(V,trns) = sigmd(V,trns, 27.1, -7.20);  # m^3
    CaTHinf(V,trns) = sigmd(V,trns, 32.1, 5.50);  # h
    CaSMinf(V,trns) = sigmd(V,trns, 33.0, -8.1);  # m^3
    CaSHinf(V,trns) = sigmd(V,trns, 60.0, 6.20);  # h
    HMinf(V,trns) = sigmd(V,trns, 70.0, 6.0);  # m
    KdMinf(V,trns) = sigmd(V,trns, 12.3, -11.8);  # m^4
    KCaMinf(V,trns,IntCa) = (IntCa/(IntCa + 3.0))*sigmd(V,trns, 28.3, -12.6); # m^4
    AMinf(V,trns) = sigmd(V,trns, 27.2, -8.70); # m^3
    AHinf(V,trns) = sigmd(V,trns, 56.9, 4.90);  # h

    # Time Constants (ms)
    tauNaM(V,trns) = tauX(V, 1.32, 1.26, 120.0, -25.0, trns);
    tauNaH(V,trns) = tauX(V, 0.00, -0.67, 62.9, -10.0, trns)*tauX(V, 1.50, -1.00, 34.9, 3.60, trns);
    tauCaTM(V,trns) = tauX(V, 21.7, 21.3, 68.1, -20.5, trns);
    tauCaTH(V,trns) = tauX(V, 105.0, 89.8, 55.0, -16.9, trns);
    tauCaSM(V,trns) = spectau(V, 1.40, 7.00, 27.0, 10.0, 70.0, -13.0, trns);
    tauCaSH(V,trns) = spectau(V, 60.0, 150.0, 55.0, 9.00, 65.0, -16.0, trns);
    tauHM(V,trns) = tauX(V, 272.0, -1499.0, 42.2, -8.73, trns);
    tauKdM(V,trns) = tauX(V, 7.20, 6.40, 28.3, -19.2, trns);
    tauKCaM(V,trns) = tauX(V, 90.3, 75.1, 46.0, -22.7, trns);
    tauAM(V,trns) = tauX(V, 11.6, 10.4, 32.9, -15.2, trns);
    tauAH(V,trns) = tauX(V, 38.6, 29.2, 38.9, -26.5, trns);

    # Calcium Reverse Potential
    R = 8.314*1000.0;  # Ideal Gas Constant (*10^3 to put into mV)
    temp = 10.0;  # Temperature in Celcius; Temperature of the Sea lawl
    z = 2.0;  # Valence of Calcium Ions
    Far = 96485.33;  # Faraday's Constant
    CaOut = 3000.0;  # Outer Ca Concentration (uM)
    ECa(CaIn) = ((R*(273.15 + temp))/(z*Far))*NaNMath.log(CaOut/CaIn);

    # Ionic Currents (mV / ms)
    iNa(g,m,h,V) = iIonic(g, m, h, 3.0, V, ENa);
    iCaT(g,m,h,V,CaIn) = iIonic(g, m, h, 3.0, V, ECa(CaIn));
    iCaS(g,m,h,V,CaIn) = iIonic(g, m, h, 3.0, V, ECa(CaIn));
    iH(g,m,V) = iIonic(g,  m, 1.0, 1.0, V, EH);
    iKd(g,m,V,t) = iIonic(g, m, 1.0, 4.0, V, EK(t));
    iKCa(g,m,V,t) = iIonic(g, m, 1.0, 4.0, V, EK(t));
    iA(g,m,h,V,t) = iIonic(g, m, h, 3.0, V, EK(t));
    iL(V) = iIonic(.01, 1.0, 1.0, 1.0, V, EL);
    iApp(t) = 0;

    # Sensor Equations
    F(FM,FH) = GF*FM^2*FH;
    S(SM,SH) = GS*SM^2*SH;
    D(DM)    = GD*DM^2;

    sensorFuse(EF,ES,ED) = exp(-(((EF/(sensorFuseDelta*100))^8 + (ES/(sensorFuseDelta*8))^8 + (ED/(sensorFuseDelta*15))^8)^(1/8)));

    ############################################################################################
    #=
    # Capacitance Implictly Set to 1
    u[1] - V
    u[2 thru 12] - current gating // (M - activation, D - Inactivation)
      2 - NaM
      3 - NaH
      4 - CaTM
      5 - CaTH
      6 - CaSM
      7 - CaSH
      8 - HM
      9 - KdM
      10 - KCaM
      11 - AM
      12 - AH
    u[13] - IntCa
    u[14 thru 20] - dynamic maximal conductances
      14 - gNa
      15 - gCaT
      16 - gCaS
      17 - gH
      18 - gKd
      19 - gKCa
      20 - gA
    u[21 thru  31] - dynamic centers of activation/inactivation curves
      21 - NaM
      22 - NaH
      23 - CaTM
      24 - CaTH
      25 - CaSM
      26 - CaSH
      27 - HM
      28 - KdM
      29 - KCaM
      30 - AM
      31 - AH
    u[32 thru 36] - Activation/Inactivation of Sensors
      32 - FM
      33 - FH
      34 - SM
      35 - SH
      36 - DM
    u[37 thru 40] - Combining Calcium Sensor Information / Homeostatic Gating
      37 - EF
      38 - ES
      39 - ED
      40 - Alpha
    =#

    function f(du,u,p,t)
        du[1] = -iL(u[1])-iNa(u[14],u[2],u[3],u[1])-iCaT(u[15],u[4],u[5],u[1],u[13])-iCaS(u[16],u[6],u[7],u[1],u[13])-iH(u[17],u[8],u[1])-iKd(u[18],u[9],u[1],t)-iKCa(u[19],u[10],u[1],t)-iA(u[20],u[11],u[12],u[1],t)
        du[2]  = (NaMinf(u[1],u[21]) - u[2])/tauNaM(u[1],u[21])
        du[3]  = (NaHinf(u[1],u[22]) - u[3])/tauNaH(u[1],u[22])
        du[4]  = (CaTMinf(u[1],u[23]) - u[4])/tauCaTM(u[1],u[23])
        du[5]  = (CaTHinf(u[1],u[24]) - u[5])/tauCaTH(u[1],u[24])
        du[6]  = (CaSMinf(u[1],u[25]) - u[6])/tauCaSM(u[1],u[25])
        du[7]  = (CaSHinf(u[1],u[26]) - u[7])/tauCaSH(u[1],u[26])
        du[8]  = (HMinf(u[1],u[27]) - u[8])/tauHM(u[1],u[27])
        du[9]  = (KdMinf(u[1],u[28]) - u[9])/tauKdM(u[1],u[28])
        du[10] = (KCaMinf(u[1],u[29],u[13]) - u[10])/tauKCaM(u[1],u[29])
        du[11] = (AMinf(u[1],u[30]) - u[11])/tauAM(u[1],u[30])
        du[12] = (AHinf(u[1],u[31]) - u[12])/tauAH(u[1],u[31])
        du[13] = (-.94*(iCaT(u[15],u[4],u[5],u[1],u[13])+iCaS(u[16],u[6],u[7],u[1],u[13]))-u[13]+.05)/20
        du[14] = ((1*(FBar-F(u[32],u[33])) + 0*(SBar-S(u[34],u[35])) + 0*(DBar-D(u[36])))*u[14]-gamma*u[14]^3)*u[40]/tauG
        du[15] = ((0*(FBar-F(u[32],u[33])) + 1*(SBar-S(u[34],u[35])) + 0*(DBar-D(u[36])))*u[15]-gamma*u[15]^3)*u[40]/tauG
        du[16] = ((0*(FBar-F(u[32],u[33])) + 1*(SBar-S(u[34],u[35])) + 0*(DBar-D(u[36])))*u[16]-gamma*u[16]^3)*u[40]/tauG
        du[17] = ((0*(FBar-F(u[32],u[33])) + 1*(SBar-S(u[34],u[35])) + 1*(DBar-D(u[36])))*u[17]-gamma*u[17]^3)*u[40]/tauG
        du[18] = ((1*(FBar-F(u[32],u[33])) + -1*(SBar-S(u[34],u[35])) + 0*(DBar-D(u[36])))*u[18]-gamma*u[18]^3)*u[40]/tauG
        du[19] = ((0*(FBar-F(u[32],u[33])) + -1*(SBar-S(u[34],u[35])) + -1*(DBar-D(u[36])))*u[19]-gamma*u[19]^3)*u[40]/tauG
        du[20] = ((0*(FBar-F(u[32],u[33])) + -1*(SBar-S(u[34],u[35])) + -1*(DBar-D(u[36])))*u[20]-gamma*u[20]^3)*u[40]/tauG
        du[21] = ((-1*(FBar-F(u[32],u[33])) + 0*(SBar-S(u[34],u[35])) + 0*(DBar-D(u[36])))-delta*u[21]^3)*u[40]/tauHalfac
        du[22] = ((1*(FBar-F(u[32],u[33])) + 0*(SBar-S(u[34],u[35])) + 0*(DBar-D(u[36])))-delta*u[22]^3)*u[40]/tauHalfac
        du[23] = ((0*(FBar-F(u[32],u[33])) + -1*(SBar-S(u[34],u[35])) + 0*(DBar-D(u[36])))-delta*u[23]^3)*u[40]/tauHalfac
        du[24] = ((0*(FBar-F(u[32],u[33])) + 1*(SBar-S(u[34],u[35])) + 0*(DBar-D(u[36])))-delta*u[24]^3)*u[40]/tauHalfac
        du[25] = ((0*(FBar-F(u[32],u[33])) + -1*(SBar-S(u[34],u[35])) + 0*(DBar-D(u[36])))-delta*u[25]^3)*u[40]/tauHalfac
        du[26] = ((0*(FBar-F(u[32],u[33])) + 1*(SBar-S(u[34],u[35])) + 0*(DBar-D(u[36])))-delta*u[26]^3)*u[40]/tauHalfac
        du[27] = ((0*(FBar-F(u[32],u[33])) + -1*(SBar-S(u[34],u[35])) + -1*(DBar-D(u[36])))-delta*u[27]^3)*u[40]/tauHalfac
        du[28] = ((-1*(FBar-F(u[32],u[33])) + 1*(SBar-S(u[34],u[35])) + 0*(DBar-D(u[36])))-delta*u[28]^3)*u[40]/tauHalfac
        du[29] = ((0*(FBar-F(u[32],u[33])) + 1*(SBar-S(u[34],u[35])) + 1*(DBar-D(u[36])))-delta*u[29]^3)*u[40]/tauHalfac
        du[30] = ((0*(FBar-F(u[32],u[33])) + 1*(SBar-S(u[34],u[35])) + 1*(DBar-D(u[36])))-delta*u[30]^3)*u[40]/tauHalfac
        du[31] = ((0*(FBar-F(u[32],u[33])) + -1*(SBar-S(u[34],u[35])) + -1*(DBar-D(u[36])))-delta*u[31]^3)*u[40]/tauHalfac
        du[32] = (sigmd(iCaT(u[15],u[4],u[5],u[1],u[13])+iCaS(u[16],u[6],u[7],u[1],u[13]),0,14.8,1) - u[32])/.5
        du[33] = (sigmd(-1*(iCaT(u[15],u[4],u[5],u[1],u[13])+iCaS(u[16],u[6],u[7],u[1],u[13])),0,-9.8,1) - u[33])/1.5
        du[34] = (sigmd(iCaT(u[15],u[4],u[5],u[1],u[13])+iCaS(u[16],u[6],u[7],u[1],u[13]),0,7.2,1) - u[34])/50
        du[35] = (sigmd(-1*(iCaT(u[15],u[4],u[5],u[1],u[13])+iCaS(u[16],u[6],u[7],u[1],u[13])),0,-2.8,1) - u[35])/60
        du[36] = (sigmd(iCaT(u[15],u[4],u[5],u[1],u[13])+iCaS(u[16],u[6],u[7],u[1],u[13]),0,3,1) - u[36])/500
        du[37] = ((FBar - F(u[32],u[33])) - u[37])/tauS
        du[38] = ((SBar - S(u[34],u[35])) - u[38])/tauS
        du[39] = ((DBar - D(u[36])) - u[39])/tauS
        du[40] = (sigmd(-1*sensorFuse(u[37],u[38],u[39]),0,alphaOffset,-.01)-u[40])/tauAlpha
    end

    ############################################################################################

    print("ICs: ", ICs, "\n");

    # Savetimes just enforces when to save
    # Actually, time steps of integrator are adaptive!
    tspan = (0.0,tEnd);
    SaveTimes = (0.0:.5:tEnd); #in milliseconds
    prob = ODEProblem(f,ICs,tspan);

    function save_func(u,t,integ)
      return copy(u)     
    end
    saved_values = SavedValues(Float64, Matrix{Float64});
    collection_cb = SavingCallback(save_func, saved_values, saveat = SaveTimes);

    diskTimes = collect.([cSize:cSize:tEnd]);
    condition(u, t, integrator) = t ∈ diskTimes[1]
    function save_to_disk_and_clear_cb!(integ)

        fileNum = lpad(Int(integ.t/cSize),3,"0");

        V        = [eachVal[1] for eachVal in saved_values.saveval];
        V        = hcat(V...);
        conds    = [eachVal[14:20] for eachVal in saved_values.saveval];
        conds    = hcat(conds...);
        halfAcs  = [eachVal[21:31] for eachVal in saved_values.saveval];
        halfAcs  = hcat(halfAcs...);
        alpha    = [eachVal[40] for eachVal in saved_values.saveval];
        alpha    = hcat(alpha...);
        erSens   = [eachVal[37:39] for eachVal in saved_values.saveval];
        erSens   = hcat(erSens...);
        #sens    = [eachVal[32:36] for eachVal in saved_values.saveval];
        #sens    = hcat(sens...);
        EKvals   = EK.(saved_values.t);

        outFile = matopen(hash*"/xxxNEWxxx_$(fileNum).mat","w");
          write(outFile,"tme",saved_values.t);
          #write(outFile,"sol",saved_values.saveval);
          write(outFile,"conds",conds);
          write(outFile,"V",V)
          write(outFile,"halfAcs",halfAcs)
          write(outFile,"alfa",alpha)
          write(outFile,"ICs",ICs)
          write(outFile,"taus",[tauG tauHalfac tauS tauAlpha])
          write(outFile,"errorSensors",erSens)
          write(outFile,"EKvals",EKvals)
          write(outFile,"lastValue",saved_values.saveval[end])
        close(outFile);
        empty!(saved_values.t)
        empty!(saved_values.saveval)
    end
    saveToDisk_cb = DiscreteCallback(condition, save_to_disk_and_clear_cb!)

    cbset = CallbackSet(collection_cb,saveToDisk_cb);

    @time sim = solve(prob, Tsit5(), maxiters=1e10, save_on=false, reltol = reltolVal, tstops = diskTimes[1], callback = cbset);

    # rename files
    files = readdir(hash)
    mat_files = filter(f -> startswith(f, "xxxNEWxxx"), files)

    for file in mat_files
        new_name = hash[end-5:end]*file[10:end]
        a = hash*"/"*file;
        b = hash*"/"*new_name;
        mv(a,b)
    end
end

