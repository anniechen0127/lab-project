mat=1.8e20 #moles in atm
alpha=1e6/mat/12*1e15 #conversion factor from Pg to atmospheric carbon concentration
lamb=3.8/4.8 #Climate sensitivity(Wm-2K-1)
kappa_l=10e2 #I think in the table this is 9.3*100, Land surface heat capacity Ka-1(Wm-2)-1
c_amp=1.1 #carbon feedback amplification factor
beta_l=3.5 #Bioshpere CO2 fertilization parameter (Pgppm-1)
beta_o=2.4 #Ocean carbon diffusion parameter(Pgppm-1)

param_clim_ocn_init=600
gamma_l=-0.13 #Biosphere temperature response(PgK-1)
gamma_o=-0.2 #Ocean carbon solubility response(PgK-1)


beta_od=.5 #Deep shallow ocean carbon diffusion coefficient(Pgppm-1)
docn_init=1.0000e+05
kappa_o=12.6 #Atmosphere ocean diffusion coefficient Ka-1(Wm-2)-1
Do=0.21 #I change this from kdeep to Do


aco2c = 290  # atmospheric co2 concentration in ppm
cina = aco2c/alpha  # carbon in atm ; ppm -> Pg

oco2c = 290  # the ocean is in balance with atm
cino = 600
rho_o=oco2c/cino  # 

odco2c = 290 # the deep ocean is
cinod = 1e6
rho_od= odco2c/cinod
