import numpy as np

def computeAirIndexCiddor(wlen, P, T, H, xc):

    P     = float(P);
    T     = float(T);
    H     = float(H);
    xc    = float(xc);
    
    # h is the fractionnal humidity value
    h = H / 100.0;

    # Sigma is in microns^-1
    sigma = np.divide(1.0e-6, wlen);
    sigma2 = sigma * sigma;
    sigma4 = sigma2 * sigma2;
    sigma6 = sigma4 * sigma2;
    print("wlen");
    print(wlen);

    print("sigma");
    print(sigma);

    # Temperatures in Kelvin !
    T2 = T * T;

    # Following the steps from Ciddor 1996 in Appendix B
    #************ Point 1 ************
    # Find svp (see below eq. 4)
    A =  1.2378847e-5;
    B = -1.9121316e-2;
    C = 33.93711047;
    D = -6.3431645e3;
    svp = np.exp(A * T2 + B * T + C + D / T);
    print("svp")
    print(svp)

    # find pw
    pw = h * svp;

    # find f (see below eq. 4)
    t     = T - 273.15;
    t2    = t**2;
    alpha = 1.00062;
    beta  = 3.14e-8;
    gamma = 5.6e-7;
    f     = alpha + beta * P + gamma * t2;
    print("f")
    print(f)

    # find xw
    xw     = f * h * svp / P;

    #************ Point 2 ************
    # Ciddor 1996, eq. 1
    k0 = 238.0185;
    k1 = 5792105.0;
    k2 = 57.362;
    k3 = 167917.0;
    n_as = 1 + 1e-8 * (k1 / (k0 - sigma2) + k3 / (k2 - sigma2));

    # Ciddor 1996, eq. 2
    n_axs = 1 + (n_as - 1)*(1 + 0.534e-6 * (xc - 450));
    n_gasx = 1 + 1e-8 * (k1 * (k0 + sigma2) / (k0 - sigma2)**2 + k3 * (k2 + sigma2) / (k2 - sigma2)**2) * (1 + 0.534e-6 * (xc - 450));
    print("n_gasx");
    print(n_gasx);

    #************ Point 3 ************
    # Find Ma (done in the function computeBIPM_Density)

    #************ Point 4 ************
    # Find Za (done in point 6)
    # Ciddor 1996, eq. 12
    Pa  = 101325;
    Ta  = 15 + 273.15;
    xwa = 0;

    #************ Point 5 ************
    # Find Zw (done in point 7)
    Pw  = 1333;
    Tw  = 20 + 273.15;
    xww = 1;

    #************ Point 6 ************
    # compute rho_axs
    rho_axs = computeBIPM_Density(Pa, Ta, xwa, xc);
    print("rho_axs");
    print(rho_axs);

    # compute rho_ws
    rho_ws = computeBIPM_Density(Pw, Tw, xww, xc);

    #************ Point 7 ************
    # Compressibility of moist air under experimental conditions
    # (eq. 12) Done in point 8 through computeBIPM_Density

    #************ Point 8 ************
    # Density of the dry component of the moist air
    # rho_a = P * Ma * (1 - xw) / (Z * R * T);
    rho_a = computeBIPM_Density(P, T, xw, xc);

    #************ Point 9 ************
    # Density of the water vapour component of the moist air
    #rho_w = P * Mw * xw / (Z * R * T);
    rho_w = computeBIPM_Density(pw, T, 1, xc);

    # Ciddor 1996, eq. 3
    #**************** vapeur d'eau standard:
    cf =   1.022;
    w0 = 295.235;
    w1 =   2.6422;
    w2 =  -0.032380;
    w3 =   0.004028;
    n_ws  = 1 + 1e-8 * cf * (w0 + w1 * sigma2 + w2 * sigma4 + w3 * sigma6);
    n_gws = 1 + 1e-8 * cf * (w0 + 3 * w1 * sigma2 + 5 * w2 * sigma4 + 7 * w3 * sigma6);
    print("n_ws");
    print(n_ws);

    #************ Point 10 ************
    # Ciddor 1996, eq. 5
    n_prop = 1 + (rho_a / rho_axs) * (n_axs - 1) + (rho_w / rho_ws) * (n_ws - 1);

    n_gprop = 1 + (rho_a / rho_axs) * (n_gasx - 1) + (rho_w / rho_ws) * (n_gws - 1);
    print("n_prop");
    print(n_prop);
        
    return n_prop

#******************************************************************************

def computeCompressibility(P, T, xw):
    # P in pascals
    # T in Kelvins
    # xw is the molar fraction of water vapor in moist air
    P     = float(P);
    T     = float(T);
    xw    = float(xw);

    # Coefficients to compute Z
    a0 =  1.58123e-6;
    a1 = -2.9331e-8;
    a2 =  1.1043e-10;
    b0 =  5.707e-6;
    b1 = -2.051e-8;
    c0 =  1.9898e-4;
    c1 = -2.376e-6;
    d  =  1.83e-11;
    e  = -0.765e-8;

    # t in Celsius
    t = T - 273.15;

    # Compute compressibility of gas, from Ciddor 1996, eq. 12
    Z = 1 - (P / T) * (a0 + a1 * t + a2 * t**2 + (b0 + b1 * t) * xw + (c0 + c1 * t) * xw**2) + (P / T) ** 2 * (d + e * xw**2);

    return Z;

#******************************************************************************

def computeBIPM_Density(P, T, xw, xc):
    # P in pascals
    # T in Kelvins
    # xw is the molar fraction of water vapor in moist air
    # xc ppm of CO2
    P     = float(P);
    T     = float(T);
    xw    = float(xw);
    
    # gas constant,
    R = 8.314510;

    # Mw is the molar mass of water vapor
    Mw = 0.018015;

    # molar mass of dry air containing xc ppm of CO2
    Ma = 1e-3 * (28.9635 + 12.011e-6 * (xc - 400));

    # Compressibility of gas
    Z = computeCompressibility(P, T, xw);

    # BIPM density, Ciddor 1996 eq. 4
    rho = (P * Ma / (Z * R * T)) * (1 - xw * (1 - Mw / Ma));

    return rho;

#******************************************************************************

def computeAirIndexMathar(wlen, P, T, H):
     #  DESCRIPTION
     #  Compute the chromatic air index, using the formulae from
     #  Mathar (2007, journal of optics A, 9, 470)
    
    P     = float(P);
    T     = float(T);
    H     = float(H);
    
    # wlen in m, convert it to cm
    sigma = 1.0 / wlen * 1e-2 / 1.0;

    # T in Kelvins
    Tref = 273.15 + 17.5;

    # Pref in Pascals
    Pref = 75000.;

    # Href in percentage of humidity
    Href = 10.;

    if any(wlen >=1.3e-6) and any(wlen <= 2.5e-6):
        print("JHK bands")
        
        sigmaref = 1e4 / 2.25;
        
        Cjref = [  0.200192e-3, 0.113474e-9, -0.424595e-14, 0.100957e-16, -0.293315e-20, 0.307228e-24];

        CjT = [  0.588625e-1, -0.385766e-7, 0.888019e-10, -0.567650e-13, 0.166615e-16, -0.174845e-20];

        CjTT = [-3.01513, 0.406167e-3, -0.514544e-6, 0.343161e-9, -0.101189e-12, 0.106749e-16];

        CjH = [-0.103945e-7, 0.136858e-11, -0.171039e-14, 0.112908e-17, -0.329925e-21, 0.344747e-25];

        CjHH = [  0.573256e-12, 0.186367e-16, -0.228150e-19, 0.150947e-22, -0.441214e-26, 0.461209e-30];

        CjP = [  0.267085e-8, 0.135941e-14, 0.135295e-18, 0.818218e-23, -0.222957e-26, 0.249964e-30];

        CjPP = [  0.609186e-17, 0.519024e-23, -0.419477e-27, 0.434120e-30, -0.122445e-33, 0.134816e-37];

        CjTH = [ 0.497859e-4, -0.661752e-8, 0.832034e-11, -0.551793e-14, 0.161899e-17, -0.169901e-21];

        CjTP =  [0.779176e-6, 0.396499e-12, 0.395114e-16, 0.233587e-20, -0.636441e-24, 0.716868e-28];

        CjHP = [-0.206567e-15, 0.106141e-20, -0.149982e-23, 0.984046e-27, -0.288266e-30, 0.299105e-34];
        
    elif any(wlen >=2.8e-6) and any(wlen <= 4.2e-6):
        print("L band")
        sigmaref = 1e4 / 3.4;
        
        Cjref = [  0.200049e-3, 0.145221e-9, 0.250951e-12, -0.745834e-15, -0.161432e-17, 0.352780e-20];

        CjT = [  0.588432e-1, -0.825182e-7, 0.137982e-9, 0.352420e-13, -0.730651e-15, -0.167911e-18];

        CjTT = [-3.13579, 0.694124e-3, -0.500604e-6, -0.116668e-8, 0.209644e-11, 0.591037e-14];

        CjH = [-0.108142e-7, 0.230102e-11, -0.154652e-14, -0.323014e-17, 0.630616e-20, 0.173880e-22];

        CjHH = [  0.586812e-12, 0.312198e-16, -0.197792e-19, -0.461945e-22, 0.788398e-25, 0.245580e-27];

        CjP = [ 0.266900e-8, 0.168162e-14, 0.353075e-17, -0.963455e-20, -0.223079e-22, 0.453166e-25];

        CjPP = [ 0.608860e-17, 0.461560e-22, 0.184282e-24, -0.524471e-27, -0.121299e-29, 0.246512e-32];

        CjTH = [ 0.517962e-4, -0.112149e-7, 0.776507e-11, 0.172569e-13, -0.320582e-16, -0.899435e-19];

        CjTP =  [0.778638e-6, 0.446396e-12, 0.784600e-15, -0.195151e-17, -0.542083e-20, 0.103530e-22];

        CjHP = [-0.217243e-15, 0.104747e-20, -0.523689e-23, 0.817386e-26, 0.309913e-28, -0.363491e-31];

    elif any(wlen >=4.35e-6) and any(wlen <= 5.2e-6):
        print("M band")
        
        sigmaref = 1e4 / 4.8;
        
        Cjref = [ 0.200020e-3, 0.275346e-9, 0.325702e-12, -0.693603e-14, 0.285610e-17, 0.338758e-18];

        CjT = [ 0.590035e-1, -0.375764e-6, 0.134585e-9, 0.124316e-11, 0.508510e-13, -0.189245e-15];

        CjTT = [-4.09830, 0.250037e-2, 0.275187e-6, -0.653398e-8, -0.310589e-9, 0.127747e-11];

        CjH = [-0.140463e-7, 0.839350e-11, -0.190929e-14, -0.121399e-16, -0.898863e-18, 0.364662e-20];

        CjHH = [ 0.543605e-12, 0.112802e-15, -0.229979e-19, -0.191450e-21, -0.120352e-22, 0.500955e-25];

        CjP = [ 0.266898e-8, 0.273629e-14, 0.463466e-17, -0.916894e-19, 0.136685e-21, 0.413687e-23];

        CjPP = [ 0.610706e-17, 0.116620e-21, 0.244736e-24, -0.497682e-26, 0.742024e-29, 0.224625e-30];

        CjTH = [ 0.674488e-4, -0.406775e-7, 0.289063e-11, 0.819898e-13, 0.468386e-14, -0.191182e-16];

        CjTP =  [0.778627e-6, 0.593296e-12, 0.145042e-14, 0.489815e-17, 0.327941e-19, 0.128020e-21];

        CjHP = [-0.211676e-15, 0.487921e-20, -0.682545e-23, 0.942802e-25, -0.946422e-27, -0.153682e-29];
        
    elif any(wlen >=7.5e-6) and any(wlen <= 14.1e-6):
        
        print("N band")
        sigmaref = 1e4 / 10.1;
        
        Cjref = [ 0.199885e-3, 0.344739e-9, -0.273714e-12, 0.393383e-15, -0.569488e-17, 0.164556e-19];
        
        CjT = [ 0.593900e-1, -0.172226e-5, 0.237654e-8, -0.381812e-11, 0.305050e-14, -0.157464e-16];

        CjTT = [-6.50355, 0.103830e-1, -0.139464e-4, 0.220077e-7, -0.272412e-10, 0.126364e-12];

        CjH = [-0.221938e-7, 0.347377e-10, -0.465991e-13, 0.735848e-16, -0.897119e-19, 0.380817e-21];

        CjHH = [ 0.393524e-12, 0.464083e-15, -0.621764e-18, 0.981126e-21, -0.121384e-23, 0.515111e-26];

        CjP = [ 0.266809e-8, 0.695247e-15, 0.159070e-17, -0.303451e-20, -0.661489e-22, 0.178226e-24];

        CjPP = [ 0.610508e-17, 0.227694e-22, 0.786323e-25, -0.174448e-27, -0.359791e-29, 0.978307e-32];

        CjTH = [ 0.106776e-3, -0.168516e-6, 0.226201e-9, -0.356457e-12, 0.437980e-15, -0.194545e-17];

        CjTP =  [0.778368e-6, 0.216404e-12, 0.581805e-15, -0.189618e-17, -0.198869e-19, 0.589381e-22];

        CjHP = [-0.206365e-15, 0.300234e-19, -0.426519e-22, 0.684306e-25, -0.467320e-29, 0.126117e-30];
    else:
        print("Error!!!");

    nm1 = np.zeros(len(wlen));
    for j in range(6):
        Cj = Cjref[j] + CjT[j] * (1/T - 1/Tref) + CjTT[j] * (1/T - 1/Tref)**2 + CjH[j]  * (H - Href) + CjHH[j] * (H - Href)**2 + CjP[j]  * (P - Pref) + CjPP[j] * (P - Pref)**2 + CjTH[j] * (1/T - 1/Tref) * (H - Href) + CjTP[j] * (1/T - 1/Tref) * (P - Pref) + CjHP[j] * (H - Href) * (H - Href);

        print(Cj)
        nm1 += Cj * (sigma - sigmaref) ** float(j);
    
    n_index = 1 + nm1;

    return n_index;


