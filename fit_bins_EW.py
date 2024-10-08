import argparse
from datetime import timedelta
import time
start_time = time.time()
def fit_emission_absorption_lines_w(cube_name, survey, path_cubes, z, initfwhm, weights = True, rest = True ,EW = True):

        '''
        Fit emission and absorption lines across a broad wavelength range
        for PMAS/MUSE IFU cubes. Ignore a transition if either the
        wavelength coverage or its combination with the redshift does
        not include it. MPFIT is used for the fitting processes.

        In addition, calculate cross-correlation coefficients between 
        each pixel of a PMAS/MUSE IFU cube and the optical spectra 
        of two SNRs: Cas A and SNR 4449. Focus on five regions (rest frame):

        4700-5300 A (O III)
        6100-6500 A (O I)
        6500-6620 A (Halpha + N II)
        6620-6900 A (S II)
        7000-7600 A (Ca II + Ni II) | [Ca II] 7291, 7323 A, [Ni II] 7378 A (McCully+14)

        IMPORTANT: Assume that the stellar contribution to the continuum 
        has not been subtracted.

        AUTHORS: Héctor Martínez Rodríguez & Lluís Galbany
        '''

        import numpy as np
        import matplotlib
        # gui_env = ['TKAgg','GTKAgg','Qt4Agg','WXAgg']
#        matplotlib.use('TKAgg', force=True)#warn=False,
        import matplotlib.pyplot as plt
        from astropy.table import Table,vstack,Column
        from astropy.io import ascii
        import astropy.io.fits as fits
        from scipy import interpolate
        import os, sys, time
        from astropy.io import fits 

        # MPFIT
        from lib.mpfit import mpfit
        # My packages
        from lib.math_functions import gaussian, polynomial # Functions needed my multipeak
        from lib.multipeak import multipeak, multipeak_fit # Multipeak (Gaussian functions)
        from lib.read_IFU_cube import read_IFU_cube # Read IFU cube (wavelength, flux)
        from lib.spectral_line_norm import cut_and_normalize_spectrum_cc # Cross-correlation coefficients

        
        start_time = time.time()
        if survey in ('PMAS'): # initfwhm depends on the resolution: 5.0 for PMAS, 3.0 for MUSE

            initfwhm = 5.0

        elif survey in ('MUSE'):

            initfwhm = 3.0

        else:

            raise NameError('Please, introduce a valid survey: PMAS/MUSE')



        # Assume rest frame unless otherwise specified

        if rest:

            z = 0.0




        ## Calculate Halpha and Hbeta EWs

        #if 'starsub' not in (cube_name):

        #    EW = True 

        #else:

        #    EW = False




        pathw = os.getcwd() + '/'
        path_write = os.getcwd() + '/'
        path_data = os.getcwd() + '/data/'

        #cube_name = 'NGC2906.resam_PMAS'
        #cube_name = 'NGC2906.starsub_PMAS'
        #cube_name = 'NGC2906.resam_MUSE'
        #cube_name = 'NGC2906.starsub_MUSE'
        cube = read_IFU_cube(cube_name, survey, path_cubes)
        
  #      hi_data = fits.open(cube_name+'.fits')  
  #      cube = SpectralCube.read(hi_data) 



        l = cube[0]
        Flux = cube[1]
        Errimage = cube[2]
        hdr = cube[3]
    
    #   cube_flux = fits.open('gas.fits')[0]
    #   Flux = SpectralCube.read(cube_flux) 
    #   
    #   cube_err = fits.open('gas.fits')[1]
    #   Errimage = SpectralCube.read(cube_err) 
    #   
    #   hdr = cube_flux.header
    #   
    #   l = np.array([i for i in np.arange(hdr['CRVAL3'],hdr['CRVAL3']+hdr['NAXIS3'],hdr['CRPIX3'])])

        lnorm = l/(1+z) # Cross-correlation coefficients. Compare with rest-frame spectra




        # CAUTION: ORIGINAL DATA ARE INVERTED. THE EFFECT WILL BE COMPENSATED IN THE FINAL DATA CUBE
        nrows = np.shape(Flux)[0]
        rows = np.arange(nrows)
        ncols = np.shape(Flux)[1]
        cols = np.arange(ncols)




        #### All the transitions (not including Halpha + NII_a,b + Halpha_SNR; OIII_SNR; OIII, OI, CaII blendings when there is a SNR ) ####

        ### N.B.: [Ca II] 7291, 7323 A, [Ni II] 7378 A (McCully+14) -> CaII blending at 7323 A does NOT belong to the galactic emissions ###
        ### Consequently, Ca II is included at the end, after all the transitions ###

        ### N.B.: The previous version of this code fitted the Hbeta_SNR component. However, this is not too interesting, because ###
        ### Halpha_SNR is more prominent. Now, this version looks for 1) a possible OIII_SNR component around 5000 A and 2) many SNR ###
        ### transitions blended around OIII (5000 A), O I (6300 A) and CaII (7323 A), as it happens with CasA and SNR4449 ### 


        ids = ['[O II] 3727',\
               '[Ne III] 3868',\
               'Heps 3970',\
               'Hdelt 4101',\
               'Hgamma 4340',\
               '[O III] 4363',\
               'Hbeta 4861',\
               '[Ne I] 4931',\
               '[O III] 4959',\
               '[O III] 5007',\
               '[Mg I] 5175',\
               '[???] 5200',\
               '[He II] 5413',\
               '[O I] 5577',\
               '[N II] 5756',\
               'DIB 5780',\
               '[He I] 5875',\
               '[Na I] 5889',\
               '[Na I] 5895',\
               'DIB 6283',\
               '[O I] 6300',\
               '[S III] 6312',\
               '[O I] 6364',\
               '[S II] 6716',\
               '[S II] 6731',\
               '[???] 6675',\
               '[???] 6825',\
               '[Ar III] 7138',\
               '[K I] 7665',\
               '[K I] 7699',\
               '[Ca II] 8498',\
               '[Ca II] 8542',\
               '[Ca II] 8662',\
               '[S III] 9069',\
               '[Ca II] 7323']    

        lam = [3726.04,\
               3868.76,\
               3970.1,\
               4101.7,\
               4340.47,\
               4363.209,\
               4861.33,\
               4930.944,\
               4958.911,\
               5006.843,\
               5175.0,\
               5200.0,\
               5412.0,\
               5577.34,\
               5756.0,\
               5780.0,\
               5875.62,\
               5889.0,\
               5895.0,\
               6283.0,\
               6300.304,\
               6312.06,\
               6363.776,\
               6716.440,\
               6730.815,\
               6675.0,\
               6825.0,\
               7137.0,\
               7665.0,\
               7699.0,\
               8498.0,\
               8542.0,\
               8662.0,\
               9069.0,\
               7323.0] # 34 transitions (NOT including CaII)




        ids = np.asarray(ids)
        lam = np.asarray(lam)



        # S II is a doublet and thus will be fit separately. Same for O II and the combination Hgamma + O III (4340 & 4363)
        wh_HgOIII_em = np.where((lam == 4340.47) | (lam == 4363.209)) [0]
        wh_SII_em = np.where((lam == 6716.440) | (lam == 6730.815)) [0]
        wh_OII_em = np.where((lam == 3726.04)) [0] # [3726.04 & 3728.80]


        # 5175.0, 5780.0, 6283.0, 7665.0, 7699.0, 8498.0, 8542.0, 8662.0 -> absorptions
        wh_abs = np.where((lam == 5175.0) | (lam == 5780.0) | (lam == 6283.0) | (lam == 7665.0) | (lam == 7699.0) | (lam == 8498.0)\
                          | (lam == 8542.0) | (lam == 8662.0)) [0]


        # Na I presents a double absorption
        wh_NaI_abs = np.where((lam == 5889.0) | (lam == 5895.0)) [0]


        # Finally, the three possible blendings around OIII, OI and CaII
        # N.B.: OIII will be fitted separately (OIII + OII_SNR), as well as OI (single emission). 
        # These are huge Gaussians that will be used just to detect SNRs!
        wh_SNR_bl = np.where((lam == 5006.843) | (lam == 6300.304) | (lam == 7323.0)) [0]



        # Index location of the single emission lines. Note that CaII_blending is removed (wh_SNR_bl[-1])
        # Opposite set
        indwh = np.concatenate((wh_HgOIII_em, wh_SII_em, wh_OII_em, wh_abs, wh_NaI_abs, np.asarray([wh_SNR_bl[-1]])), axis=0) 
        nlam = len(lam)

        wh_single_em = range(nlam)

        for iidx in sorted(indwh, reverse = True): # Remove S II, O II, Hgamma + O III, Na I and single absorptions
            #print(indwh,iidx,wh_single_em,wh_single_em[iidx])
            wh_single_em=list(wh_single_em)
            del wh_single_em[iidx]




        wh_HgOIII_em = list(wh_HgOIII_em)
        wh_SII_em = list(wh_SII_em)
        wh_OII_em = list(wh_OII_em)
        wh_abs = list(wh_abs)
        wh_NaI_abs = list(wh_NaI_abs)
        wh_SNR_bl = list(wh_SNR_bl)



            
        ########################## Important definitions ##########################

        lres = [6548.05, 6562.80, 6583.45] # S II, Halpha, N II
        lim0 = [-43. + lres[1], 47. + lres[1]] # Fit limits for N II, Halpha, N II



        # Flux for each spectral line: lam (- CaII) + lres + Halpha_SNR + OIII, OI, CaII blendings = 41 transitions
        flux = np.zeros( ( nrows, ncols, nlam + 6 ) ) 
        eflux = np.zeros( ( nrows, ncols, nlam + 6 ) ) # Flux errors

        vel = np.zeros( ( nrows, ncols, nlam + 6 ) ) # Velocity for each spectral line
        evel = np.zeros( ( nrows, ncols, nlam + 6 ) ) # Velocity errors

        linepars = np.zeros( ( nrows, ncols, 3*(nlam + 6) ) ) # Flux, lambda, sigma for each line
        elinepars = np.zeros( ( nrows, ncols, 4*(nlam + 6) ) ) # Errors for each line + std(spectrum - fit) * np.max(spectrum) [continuum scatter]

        fluxint = np.zeros( ( nrows, ncols, nlam + 6 ) ) # Integrated flux for each spectral line




        ## EWs: Only calculated in resam cubes, not in starsub
        ew = np.zeros( ( nrows, ncols, 2 ) ) # EW for Halpha and Hbeta
        dew = np.zeros( ( nrows, ncols, 2 ) ) # EW error for Halpha and Hbeta




        clight = 2.99792458E5 # Speed of light in vacuum [km/s]




        # Structure of linepars: [flux, lambda, sigma] for

        # '[NII] 6548', 'Halpha 6563', '[NII] 6583',..., 'Hbeta SNR', 'Halpha SNR'

        # elinepars has the same structure and includes max(spectrum) * std(spectrum-fit) [continuum scatter] in the last component




        #####################################################################################################################################

        # CROSS-CORRELATION COEFFICIENTS



        # Spectra for known SNRs
 #####       CasA = np.loadtxt(path_data + 'casa_integratedspec_01mar11_deredden.dat')
#####        SNR4449 = np.loadtxt(path_data + 'spec_snr4449_mmt2008_z.flm')


        # IMPORTANT READ THIS!!!!!!: 

        # PMAS's wavelength bins are given in multiples of 10, e.g. 6500, 6510, 6520, 6530,...
        # Therefore focusing on a region to find the cross-correlation coefficient (CCC) could retrieve
        # several discrete CCCs instead of just one (as it happens up there).
        # Hence, I will interpolate CasA and SNR4449 to prevent this from happening

 ########       CasA_interp = interpolate.interp1d(CasA[:,0], CasA[:,1], kind='linear');
 ########       SNR4449_interp = interpolate.interp1d(SNR4449[:,0], SNR4449[:,1], kind='linear');


 ########       # Wavelength limits for the cross-correlation
        lammin_cc = [4700, 6100, 6500, 6620, 7000]
        lammax_cc = [5300, 6500, 6620, 6900, 7600]


        nreg_cc = len(lammin_cc)


 ########       # Cross-correlation data cubes

 ########       cc_CasA = np.zeros( ( nrows, ncols, nreg_cc ) )
 ########       cc_SNR4449 = np.zeros( ( nrows, ncols, nreg_cc ) )




        #####################################################################################################################################

        # MPFIT EXAMPLE IN PYTHON


        #def myfunct(p, fjac=None, x=None, y=None, err=None):
        #    # Parameter values are passed in "p"
        #    # If fjac==None then partial derivatives should not be
        #    # computed.  It will always be None if MPFIT is called with default
        #    # flag.
        #    model = x**4
        #    # Non-negative status value means MPFIT should continue, negative means
        #    # stop the calculation.
        #    status = 0
        #    return [status, (y-model)/err]



        #x = np.arange(100)
        #p = [5.7, 2.2, 500., 1.5, 2000.]
        #y = ( p[0] + p[1]*x + p[2]*x**2 + p[3]*sqrt(x) + p[4]*np.exp(-x))
        #err = np.std(y)
        #fa = {'x':x, 'y':y, 'err':err}
        #m = mpfit(myfunct , p, functkw=fa)
        #print 'status = ', m.status
        #if (m.status <= 0): print 'error message = ', m.errmsg
        #print 'parameters = ', m.params

        #####################################################################################################################################




        # Loop over rows and columns
        #t0 = time.time()


        for i in rows:
            
            for j in cols:
                
                #print 'Row #' + str(i), 'Column #' + str(j)
                
                sp = Flux[i][j]
                #esp = np.ones(len(sp)) # Default errors
                ins = np.where( (l * (1+z) > 3800) & (l * (1+z) <= 9000) )[0] # Index limit
                #print(sp)
                maxsp = np.max(sp[ins]) 


                # IMPORTANT: the final errors cannot be trusted if esp = np.ones (see the MPFIT documentation):
                # [If the fit is unweighted (i.e. no errors were given, or the weights were uniformly set to unity),...]
                if weights: 

                    esp = Errimage[i][j]
                    # Correct both bad values and bad pixels with 
                    wh_badesp = np.where( (esp <=0.) | (esp > 50 * np.median(esp) ) )[0]
                    esp[wh_badesp] = np.median(esp)

                else:

                    esp = np.ones(len(sp)) # Default errors


                # Halpha wavelength (maximum flux). Beware of edge/signal-to-noise effects, reason why these wavelength limits are included
                ins = np.where( ( l >= lres[1] * (1+z) - 5. ) & ( l <= lres[1] * (1+z) + 5. ) )[0] 
                cn = len(ins)
                l1 = l[ins]
                sp1 = sp[ins]


                if np.sum(sp) > 0.: # Remove bad spectra

                
                    lha = l1[np.argmax(sp1)] # Halpha wavelength
                    # N.B.: not necessarily the maximum flux if there is a remnant; e.g., CasA and SNR4449 have a stronger OIII emission at 5007 angstroms
                    lim = [lha - 43., lha + 47.] # Fit limits in wavelength (- 43 angstroms, + 47 angstroms)

                    ins = np.where( (l > lim[0]) & (l < lim[1]) )[0] # Index limit
                    l1 = l[ins] # Lambda values within the range
                    sp1 = sp[ins] # Flux values within the range
                    # This line is commented out due to OIII (see above)
                    #maxsp = np.max(sp1)
                    sp1 = sp1/maxsp # Normalized flux values withing the range
                    esp1 = esp[ins]


                    if weights:

                        esp1 = esp[ins]/maxsp


                    z0 = (lha-lres[1])/lres[1] # First approximation for redshift

                        
                    # FIT A "FOURTH" PEAK CORRESPONDING TO A SNR (Halpha_SNR from shock lies on top of Halpha, NII_ab, NII_b from Galaxy)


                    # This array is going to have all the relevant information from the fits corresponding to the three transitions:
                    # par[0, 3, 6] = [flux_NII_A, flux_Halpha, flux_NII_B]         # Normalized flux
                    # par[1, 4, 7] = [lambda_NII_A, lambda_Halpha, lambda_NII_B]   # Wavelength
                    # par[2, 5, 8] = [sigma_NII_A, sigma_Halpha, sigma_NII_B]      # Sigma
                    # par[9, 10, 11] = [flux_Halpha_SNR, lambda_Halpha_SNR, sigma_Halpha_SNR]


                    par = np.zeros(12)
                    par[[x*3 + 2 for x in [0,1,2]]] = initfwhm/2.35482  # Sigma
                    par[[x*3 + 1 for x in [0,1,2]]] = [(1+z0)*x for x in lres] # Wavelength

                    par[11] = 3*initfwhm/2.35482  # 3-sigma for fourth peak (roughly the width of a single "peak" encompassing the three lines)
                    par[10] = par[4] # Halpha wavelength for fourth peak


                    # Explanation on par[3] and par[9]:
                    # When normalizing the flux, if there is noise whose flux is < 0, if par[3] = np.max(sp1[ins]) and par[9] = np.max(sp1[ins])
                    # or par[9] = 0, subtracting np.min(sp1) afterwards will result in a peak height > 1. Thus, an error will arise in the
                    # calculation (ERROR: parameters are not within the PARINFO limits) with status = 0
                    ins = np.where( (l1 > par[1] - 5.) & (l1 < par[1] + 5.) )[0] # Index limit for NII_A. Fit within +- 5 angstroms
                    par[0] = np.max(sp1[ins])
                    ins = np.where( (l1 > par[4] - 5.) & (l1 < par[4] + 5.) )[0] # Index limit for Halpha. Fit within +- 5 angstroms
                    #par[3] = np.max(sp1[ins])
                    par[3] = np.min(sp1[ins]) + 0.97
                    ins = np.where( (l1 > par[7] - 5.) & (l1 < par[7] + 5.) )[0] # Index limit for NII_B. Fit within +- 5 angstroms
                    par[6] = np.max(sp1[ins])
                    
                    ins = np.where( (l1 > par[1] - 15.) & (l1 < par[7] + 15.) )[0] # Index limit for fourth peak
                    par[9] = np.min(sp1[ins]) + 0.03 #np.max(sp1[ins]) # Fourth peak. Should end up being zero if there is no SNR!

                    # These correspond to the continuum (assumes a starting flat slope)
                    par = np.concatenate((par,[np.min(sp1),0.]), axis=0) # Add new elements to the list
                    par[[x*3 for x in [0,1,2,3]]] = par[[x*3 for x in [0,1,2,3]]] - np.min(sp1) # Normalized fluxes minus smallest value
                    #par = list(par) # Convert from np.array to list for MPFIT

                    parinfo = [{'value' : 0., 'fixed' : 0, 'limited' : [0,0], 'limits' : [0.,0.]} for y in range(14)]
                        

                    # Set values for parinfo. These will be used by MPFIT
                    for ii in range(3):

                        parinfo[3*ii+2]['value'] = initfwhm/2.35482
                        parinfo[3*ii]['value'] = par[3*ii]
                        parinfo[3*ii+1]['value'] = par[3*ii+1]

                        parinfo[3*ii+1]['limited'] = [1,1]
                        parinfo[3*ii+1]['limits'][0]  = par[3*ii+1]-5 # Typically, around 200-300 km/s
                        parinfo[3*ii+1]['limits'][1]  = par[3*ii+1]+5
                        parinfo[3*ii+2]['limited'] = [1,1]
                        parinfo[3*ii+2]['limits'][0]  = 0.5*par[3*ii+2]
                        parinfo[3*ii+2]['limits'][1]  = 2.0*par[3*ii+2]

                        parinfo[3*ii]['limited'][0] = 1
                        parinfo[3*ii]['limits'][0]  = 0.

                        
                    parinfo[12]['value'] = par[12]
                    parinfo[13]['value'] = par[13]


                    # Fourth peak

                    parinfo[11]['value'] = 3*initfwhm/2.35482 # Extra broad sigma
                    parinfo[9]['value'] = par[9] # Flux
                    parinfo[10]['value'] = par[10] # Wavelength

                    parinfo[10]['limited'] = [1,1]
                    #parinfo[10]['limits'][0]  = par[10]-5
                    #parinfo[10]['limits'][1]  = par[10]+5
                    parinfo[10]['limits'][0]  = par[10]-25 # 1200 km/s
                    parinfo[10]['limits'][1]  = par[10]+25
                    parinfo[11]['limited'] = [1,1]
                    parinfo[11]['limits'][0]  = 0.8*par[11] # Avoid overlapping with upper sigma limit for regular peaks
                    parinfo[11]['limits'][1]  = 3.5*par[11] # ~10-sigma upper limit for the fit

                    parinfo[9]['limited'][0] = 1
                    parinfo[9]['limits'][0]  = 0.


                    # MPFIT structure: a = {'covar':, 'niter', 'params', 'nfev', 'debug', 'perror', 'damp', 'status', 'errmsg'}
                    # Set plot = True to see spectrum and fit
                    a = multipeak_fit(l1, sp1, par, npeak = 4, parinfo = parinfo, fit = True, err = esp1, plot = False, silent = True);

                    #plt.xlabel(r'$\rm{\lambda \, [\AA]}$', fontsize=13);
                    #plt.ylabel(r'$\rm{Normalized \ flux}$', fontsize=13);
                    #plt.title(r'[Halpha + NII] 6548 & 6563 & 6583', fontsize=13);
                    #plt.tight_layout();
                    #plt.savefig(path_pictures + 'Halpha_NII_' + str(bool(i)) + 'SNR.pdf', bbox_to_inches = 'tight');



                    # Need to use fit later, otherwise change "return p, ft" to "return p" in the function
                    # Python doesn’t support call by reference (equivalent to intent(out) or intent(inout) in Fortran)
                    fit = a[1] 
                    a = a[0]

                    #print 'status = ',a.status, ',\b', 'DOF = ', len(l1) - len(a.params), ',\b', 'chi2 = ', a.fnorm, ',\b', 'niter = ', a.niter


                    #The formal 1-sigma errors in each parameter, computed from the
                    #covariance matrix.  If a parameter is held fixed, or if it touches a
                    #boundary, then the error is reported as zero.

                    sigma = a.perror


                    if a.status > 0:
                        
                        ins = -1
                        
                        ins1 = np.where( l1 < a.params[1] - 3*a.params[2] ) [0] # 3-sigma from NII_A
                        cn1 = len(ins1)
                        if cn1 > 0:
                            ins = [ins, ins1]
                            
                        ins2 = np.where( l1 > a.params[7] + 3*a.params[8] ) [0] # 3-sigma from NII_B
                        cn2 = len(ins2)
                        if cn2 > 0:
                            ins = [ins, ins2]
                            
                        if cn1 + cn2 > 5:
                            ins = ins[1:]
                            std = np.std(sp1[ins] - fit[ins])
                            ins0 = np.where( np.abs( l1 - a.params[4] ) <= 2*a.params[5] ) [0]
                            num = len(ins0)
                            
                           
                        flt = maxsp*(sp1 - polynomial(l1, a.params[12:])) # Recall the fourth peak
                        

                        for k in [0, 1, 2]: # Three lines

                            #print k,',',
                            #print 3*(k),':',3*(k)+3,',',
                            #print 4*(k),':',4*(k)+4
                        
                            flux[i, j, k] = np.sqrt(2.0*np.pi) * a.params[3*k] * a.params[3*k+2] * maxsp
                            eflux[i, j, k] = np.sqrt(2.0*np.pi) * np.sqrt( (sigma[3*k] * a.params[3*k+2])**2 + (a.params[3*k] * sigma[3*k+2])**2 ) * maxsp 
                            # S/N = flux[1]/eflux[1] + flux[1]/(maxsp*np.sqrt(num)*std)
                        
                            vel[i, j, k] = (a.params[3*k+1] - lres[k]) / lres[k] * clight
                            evel[i, j, k] = sigma[3*k+1] / lres[k] * clight
                        
                            linepars[i, j, 3*k : 3*k+3] = [a.params[3*k] * maxsp, a.params[3*k+1], a.params[3*k+2]]
                            elinepars[i, j, 4*k : 4*k+4] = [sigma[3*k] * maxsp, sigma[3*k+1], sigma[3*k+2], std*maxsp]
                        
                        
                            ins = np.where( np.abs( l1 - a.params[3*k+1] ) <= 3*a.params[3*k+2] ) [0] # 3-sigma limit
                            fluxint[i, j, k] = np.trapz(y = flt[ins], x = l1[ins])
                        
                        
                        
                        # Parameters for the fourth peak. They are placed in the end

                        flux[i, j, -1] = np.sqrt(2.0*np.pi) * a.params[9] * a.params[11] * maxsp
                        eflux[i, j, -1] = np.sqrt(2.0*np.pi) * np.sqrt( (sigma[9] * a.params[11])**2 + (a.params[9] * sigma[11])**2 ) * maxsp 
                        
                        vel[i, j, -1] = (a.params[10] - lres[1]) / lres[1] * clight
                        evel[i, j, -1] = sigma[10] / lres[1] * clight
                        
                        linepars[i, j, -3:] = [a.params[9] * maxsp, a.params[10], a.params[11]]
                        elinepars[i, j, -4:] = [sigma[9] * maxsp, sigma[10], sigma[11], std*maxsp]
                        
                        
                        ins = np.where( np.abs( l1 - a.params[10] ) <= 3*a.params[11] ) [0] # 3-sigma limit
                        fluxint[i, j, -1] = np.trapz(y = flt[ins], x = l1[ins])




                        ## EW for galactic Halpha

                        if EW:
    
                            lcont = polynomial(l1, a.params[12:]) # Recall the fourth peak

                            a = multipeak_fit(l1, sp1/lcont -1, par, npeak = 4, parinfo = parinfo, fit = True, err = esp1/lcont, plot = False, silent = True);
                            
                        #    #plt.xlabel(r'$\rm{\lambda \, [\AA]}$', fontsize=13);
                        #    #plt.ylabel(r'$\rm{Normalized \ flux}$', fontsize=13);
                        #    #plt.title(r'[Halpha + NII] 6548 & 6563 & 6583', fontsize=13);
                        #    #plt.tight_layout();
                        #    #plt.savefig(path_pictures + 'Halpha_NII_EW_' + str(bool(i)) + 'SNR.pdf', bbox_to_inches = 'tight');
                        
                            
                            fit = a[1] 
                            a = a[0]
                            
                            
                            ins = -1
                            
                            ins1 = np.where( l1 < a.params[1] - 3*a.params[2] ) [0] # 3-sigma from NII_A
                            cn1 = len(ins1)
                            if cn1 > 0:
                                ins = [ins, ins1]
                                
                            ins2 = np.where( l1 > a.params[7] + 3*a.params[8] ) [0] # 3-sigma from NII_B
                            cn2 = len(ins2)
                            if cn2 > 0:
                                ins = [ins, ins2]
                                
                            if cn1 + cn2 > 5:
                                ins = ins[1:]
                                std = np.std(sp1[ins] - fit[ins])
                                ins0 = np.where( np.abs( l1 - a.params[4] ) <= 2*a.params[5] ) [0]
                                num = len(ins0)
                                
                                
                        #    # (1+z) added by Lluis because flux has not been multiplied im map_emi_line    
                            ew[i, j, 0] = np.sqrt(2.0*np.pi) * a.params[3] * a.params[5] #* (1+z) 
                            dew[i, j, 0] = np.sqrt(2.0*np.pi) * np.sqrt( (sigma[3] * a.params[5])**2 + (a.params[3] * sigma[5])**2 ) #* (1+z)




                        z0 = (a.params[4]-lres[1])/lres[1] # Fitted redshift




                    # Record velocity and sigma to fix lambda and sigma for Hbeta when there is a SNR
                    # Note that it is outside of the "if a.status > 0" loop to be overwritten in each step.
                    # Unlike z0, these parameters are not defined before the main loop
                    #vel_Balmer_remnant = (a.params[10] - lres[1]) / lres[1] * clight
                    #sigma_Balmer_remnant = a.params[11]
                    #lambda_Hb_remnant = lam[wh_Hb_em] * (vel_Balmer_remnant / clight + 1)
                    #print vel_Balmer_remnant, lambda_Hb_remnant, sigma_Balmer_remnant




                    # ['[Ne III] 3868', 'Heps 3970', 'Hdelt 4101','Hbeta 4861', '[Ne I] 4931', '[O III] 4959',\
                    # '[O III] 4959', '[???] 5200', '[He II] 5413', '[O I] 5577', '[N II] 5756', '[He I] 5875', '[O I] 6300',\
                    #  '[S III] 6312', '[O I] 6364', '[???] 6675', '[???] 6825', '[Ar III] 7138', '[S III] 9069']



                    for ii in wh_single_em: # [1, 2, 3, 6, 7, 8, 9, 11, 12, 13, 14, 16, 20, 21, 22, 25, 26, 27, 33]

                        #print ii+3,',',
                        #print 3*(ii+3),':',3*(ii+3)+3,',',
                        #print 4*(ii+3),':',4*(ii+3)+4
                        
                        lam0 = lam[ii] * (1+z0)
                        #print lam[ii], ids[ii], lam0
                        ins = np.where( (l >= lam0 - 5) & (l <= lam0 + 5) )[0]   # Within +- 5 angstroms
                        cn = len(ins) # Equivalent no np.count_nonzero(ins)
                        
                        if cn > 0:
                            
                            l1 = l[ins]
                            sp1 = sp[ins]
                            ll = l1[np.argmax(sp1)]
                            
                            ins = np.where( (l >= ll - 25) & (l <= ll + 25) & ( l >= 3760. / (1+z) ) )[0]
                            cn = len(ins)
                            inns = np.where( (l >= ll - 50) & (l <= ll + 50) ) [0]
                            cn = len(inns)

                            if cn > 15:
                                
                                l1 = l[ins]
                                sp1 = sp[ins]
                                esp1 = esp[ins]
                                sp1 = sp1/maxsp
                    

                                if weights:

                                    esp1 = esp[ins]/maxsp

                                
                                parinfo = [{'value' : 0., 'fixed' : 0, 'limited' : [0,0], 'limits' : [0.,0.]} for y in range(5)]
                                par = [max(sp1), ll, initfwhm/2.35482, min(sp1), 0.]
                                
                                # Analogous to first parinfo structure
                                parinfo[2]['value'] = initfwhm/2.35482
                                parinfo[0]['value'] = par[0]
                                parinfo[1]['value'] = ll

                                parinfo[1]['limited'] = [1,1]
                                parinfo[1]['limits'][0]  = ll-5
                                parinfo[1]['limits'][1]  = ll+5
                                parinfo[2]['limited'] = [1,1]
                                parinfo[2]['limits'][0]  = 0.5*initfwhm/2.35482
                                parinfo[2]['limits'][1]  = 2.0*initfwhm/2.35482

                                parinfo[0]['limited'][0] = 1
                                parinfo[0]['limits'][0]  = 0.0
                                parinfo[3]['value'] = par[3]
                                parinfo[4]['value'] = 0.0
                                
                                
                                a = multipeak_fit(l1, sp1, par, npeak = 1, parinfo = parinfo, fit = True, err = esp1, plot = False, silent = True)
                                #plt.xlabel(r'$\rm{\lambda \, [\AA]}$', fontsize=13);
                                #plt.ylabel(r'$\rm{Normalized \ flux}$', fontsize=13);
                                #plt.title(ids[ii], fontsize=13);
                                #plt.tight_layout();
                                #if ('Hbeta' in ids[ii]) | ('Hdelt' in ids[ii]) | ('Heps' in ids[ii]):
                                    #plt.savefig(path_pictures + ids[ii].split()[0] + '_' + str(bool(i)) + 'SNR.pdf', bbox_to_inches = 'tight');
                                    
                                #else:
                                    #plt.savefig(path_pictures + ids[ii].split()[0][1:] + ids[ii].split()[1][:-1] + '_' + str(bool(i)) + 'SNR.pdf',\
                                #                bbox_to_inches = 'tight');
                                
                                fit = a[1] 
                                a = a[0]       
                                sigma = a.perror


                                if a.status > 0:
                        
                                    ind = -1
                        
                                    ins1 = np.where( l1 < a.params[1] - 3*a.params[2] ) [0] # Lower 3-sigma from line
                                    cn1 = len(ins1)
                                    if cn1 > 0:
                                        ind = [ind, ins1]
                            
                                    ins2 = np.where( l1 > a.params[1] + 3*a.params[2] ) [0] # Upper 3-sigma from line
                                    cn2 = len(ins2)
                                    if cn2 > 0:
                                        ind = [ind, ins2]
                            
                                    if cn1 + cn2 > 5:
                                        ind = ind[1:]
                                        std = np.std(sp1[ind] - fit[ind])
                                        ind = np.where( np.abs( l1 - a.params[1] ) <= 2*a.params[2] ) [0]
                                        num = len(ind)
                            
                            
                           
                                    flux[i, j, ii + 3] = np.sqrt(2.0*np.pi) * a.params[0] * a.params[2] * maxsp
                                    eflux[i, j, ii + 3] = np.sqrt(2.0*np.pi) * np.sqrt( (sigma[0] * a.params[2])**2 + (a.params[0] * sigma[2])**2 ) * maxsp

                                    vel[i, j, ii + 3] = (a.params[1] - lam[ii]) / lam[ii] * clight
                                    evel[i, j, ii + 3] = sigma[1] / lam[ii] * clight
                        
                                    linepars[i, j, 3*(3+ii) : 3*(3+ii)+3] = [a.params[0] * maxsp, a.params[1], a.params[2]]
                                    elinepars[i, j, 4*(3+ii) : 4*(3+ii)+4] = [sigma[0] * maxsp, sigma[1], sigma[2], std*maxsp]
                        
                                    flt = maxsp*(sp1 - polynomial(l1, a.params[3:5]))
                            
                            
                                    ins = np.where( np.abs( l1 - a.params[1] ) <= 3*a.params[2] ) [0] # 3-sigma limit
                                    if len(l1[ins]) == 1:
                                    
                                        fluxint[i, j, ii + 3] = flt[ins]
                                    
                                    else:
                                    
                                        fluxint[i, j, ii + 3] = np.trapz(y = flt[ins], x = l1[ins])

                    #----------------------------------------------------------------------------------------------------------------------------------#    




                                    ## EW for galactic Hbeta

                                    if EW:

                                        if 'Hbeta' in ids[ii]:
                        
                                            lcont = polynomial(l1, a.params[3:])

                                            a = multipeak_fit(l1, sp1/lcont -1, par, npeak = 1, parinfo = parinfo, fit = True, err = esp1/lcont, plot =False, silent = True);
                                        
                                           # plt.xlabel(r'$\rm{\lambda \, [\AA]}$', fontsize=13);
                                           # plt.ylabel(r'$\rm{Normalized \ flux}$', fontsize=13);
                                           # plt.title(ids[ii], fontsize=13);
                                           # plt.tight_layout();
                                           # plt.savefig(path_pictures + 'Hbeta_EW_' + str(bool(i)) + 'SNR.pdf', bbox_to_inches = 'tight');
                                    
                                        
                                            fit = a[1] 
                                            a = a[0]
                                        
                                        
                                            ins = -1
                                        
                                            ins1 = np.where( l1 < a.params[1] - 3*a.params[2] ) [0]
                                            cn1 = len(ins1)
                                            if cn1 > 0:
                                                ins = [ins, ins1]
                                            
                                                ins2 = np.where( l1 > a.params[1] + 3*a.params[2] ) [0]
                                                cn2 = len(ins2)
                                            if cn2 > 0:
                                                ins = [ins, ins2]
                                            
                                            if cn1 + cn2 > 5:
                                                ins = ins[1:]
                                                std = np.std(sp1[ins] - fit[ins])
                                                ind = np.where( np.abs( l1 - a.params[1] ) <= 2*a.params[2] ) [0]
                                                num = len(ind)
                                            
                                               
                                            ew[i, j, 1] = np.sqrt(2.0*np.pi) * a.params[0] * a.params[2] #* (1+z) 
                                            dew[i, j, 1] = np.sqrt(2.0*np.pi) * np.sqrt( (sigma[0] * a.params[2])**2 + (a.params[0] * sigma[2])**2 ) #* (1+z)

                                    
                                
                                
                    #----------------------------------------------------------------------------------------------------------------------------------#

                


         


                    #Cross-correlation coefficients




                    sp = Flux[i][j]
                
                    # Limit the spectrum to the cc regions
                    ins_cc = np.where((lnorm > lammin_cc[0]) & (lnorm < lammax_cc[-1]))
                    l_cc = lnorm[ins_cc]
                    sp_cc = sp[ins_cc] 

                    for kk in range(nreg_cc):

                        # Save only the spectral flux (see the [1] component in the end)
                       # CasA_norm = cut_and_normalize_spectrum_cc(l_cc, CasA_interp(l_cc), [lammin_cc[kk], lammax_cc[kk]]) [1]
                        #SNR4449_norm = cut_and_normalize_spectrum_cc(l_cc, SNR4449_interp(l_cc), [lammin_cc[kk], lammax_cc[kk]]) [1]

                        sp_norm = cut_and_normalize_spectrum_cc(l_cc, sp_cc, [lammin_cc[kk], lammax_cc[kk]]) [1]


                        # Cross-correlation coefficients. Equivalent to
                        #cc_CasA[i, j, kk] = np.correlate(CasA_norm, sp_norm)
                        #cc_SNR4449[i, j, kk] = np.correlate(SNR4449_norm, sp_norm)
                     #   cc_CasA[i, j, kk] = np.sum(CasA_norm * sp_norm)
                      #  cc_SNR4449[i, j, kk] = np.sum(SNR4449_norm * sp_norm)


                        

        ############################################################################################################################################
                
        #t1 = time.time()




        # Record final FITS data cube (structure identical to original data's) with all the fitted parameters


        prihdr = fits.Header()
        prihdr = hdr
        prihdr['COMMENT'] = '= Results from MPFIT applied to %s data cube' %(survey)
        prihdr['DATE_MPFIT'] = time.strftime("%m/%d/%Y")
        prihdu = fits.PrimaryHDU(header=prihdr)
        
        print("--- %s seconds ---" % (time.time() - start_time))

        #tbhdu = [0]*(12) # header, cols
        #clmns = [0]*(11) # flux, eflux, vel, evel, linepars, elinepars, fluxint, cc_CasA, cc_SNR4449, ew, dew
        tbhdu = [0]*(4) # header, cols
        clmns = [0]*(3) # flux, eflux, vel, evel, linepars, elinepars, fluxint, cc_CasA, cc_SNR4449

        hdu = fits.PrimaryHDU(header=prihdr)
        tbhdu[0] = hdu


        clmns[0] = fits.ImageHDU(flux.T, name = 'flux')
        clmns[0].header['CONTENT'] = 'LINE FLUX: GAUSSIAN FIT [erg/s/A/cm^2]'


#        clmns[7] = fits.ImageHDU(cc_CasA.T, name = 'CasA')
#        clmns[7].header['CONTENT'] = 'CROSS-CORRELATION COEFFICIENTS FOR CAS A'
#        clmns[7].header['PLANE01'] = '4700-5300 A (O III)'
#        clmns[7].header['PLANE02'] = '6100-6500 A (O I)'
#        clmns[7].header['PLANE03'] = '6500-6620 A (Halpha + N II)'
#        clmns[7].header['PLANE04'] = '6620-6900 A (S II)'
#        clmns[7].header['PLANE05'] = '7000-7600 A (Ca II + Ni II)'#
#

#        clmns[8] = fits.ImageHDU(cc_SNR4449.T, name = 'SNR4449')
#        clmns[8].header['CONTENT'] = 'CROSS-CORRELATION COEFFICIENTS FOR SNR 4449'
#        clmns[8].header['PLANE01'] = '4700-5300 A (O III)'
#        clmns[8].header['PLANE02'] = '6100-6500 A (O I)'
#        clmns[8].header['PLANE03'] = '6500-6620 A (Halpha + N II)'
#        clmns[8].header['PLANE04'] = '6620-6900 A (S II)'
#        clmns[8].header['PLANE05'] = '7000-7600 A (Ca II + Ni II)'#


        clmns[1] = fits.ImageHDU(ew.T, name = 'ew')
        clmns[1].header['CONTENT'] = 'EQUIVALENT WIDTH','A]'
        clmns[1].header['PLANE01'] = 'Halpha 6562.77', 'A'



        clmns[2] = fits.ImageHDU(dew.T, name = 'dew')
        clmns[2].header['CONTENT'] = 'EQUIVALENT WIDTH ERRORS','A]'
        clmns[2].header['PLANE01'] = 'Halpha 6562.77', 'A'




        #for c in range(11):
        for c in range(3):

            # Final table encompassing header and values
            tbhdu[c+1] = clmns[c]

            
        hdulist = fits.HDUList(tbhdu)
        hdulist.writeto(path_write + cube_name + '_%s.EWs.binned.fits' % (survey), overwrite = True)

def MakeGasCube(cube_name):
        import numpy as np
        import glob as glob
        from astropy.io import ascii
        from astropy.io import fits
        import time, os
        from datetime import timedelta

        cubo=fits.open('Final_Sample/'+cube_name+'/'+cube_name+'.resam.fits')[0]
        hdr=cubo.header

        outs = glob.glob('Final_Sample/'+cube_name+'/adabin_out_SII/'+cube_name+'_*.BN')

        inputs = glob.glob('Final_Sample/'+cube_name+'/adabin_out_SII/'+cube_name+'_*.txt')
        bin_map = fits.open('Final_Sample/'+cube_name+'/seg_adabin_SNSII_10_'+cube_name+'.fits')[0].data
	
        file_1 = max(outs, key=os.path.getsize)
 
        ini_model = int(list(ascii.read(file_1, data_start=235)['col1'])[0])
        end_model = int(list(ascii.read(file_1, data_start=235)['col1'])[-1])
        
        cubo = np.zeros((1 + end_model - ini_model, 1, int(np.max(bin_map))))
        errores = np.zeros((1 + end_model - ini_model, 1, int(np.max(bin_map))))
        
        for bin_index in range(1, int(np.max(bin_map))):
            try:
                try:
                    file_ = [s for s in outs if "_bin_"+str(bin_index)+"." in s][0]
                    file_in = [s for s in inputs if "_bin_"+str(bin_index)+"." in s][0]
                except IndexError: 'list index out of range'
                
                data_out = list(ascii.read(file_, data_start=235)['col3'])

                data_in = list(ascii.read(file_, data_start=235)['col2'])
          
               
                match_indices = [i for i, num in enumerate(data_out) if "*" in str(num)]

                if len(match_indices) >0:
    
                    for ind in match_indices:
                      
                        data_out[ind] = (float(data_out[ind-1]))
                    data_out = [float(num) for num in data_out]

                match_indices = [i for i, num in enumerate(data_in) if "*" in str(num)]

                if len(match_indices) >0:

                    for ind in match_indices:
                       
                       # data_in[ind] = (float(data_in[ind-1]))
                        data_in[ind] = 0.
                    data_in = [float(num) for num in data_in]

                error = list(ascii.read(file_in)['col3'])
                start_error = int(ascii.read(file_in)['col1'][0])
                end_error = int(ascii.read(file_in)['col1'][-1])

                end = 1

                if start_error >= ini_model:
                    diff = start_error-ini_model
                    error = [0] * diff + error
                    start = 0
                elif start_error < ini_model:
                    start = ini_model - start_error

                if end_error >= end_model:
                    end = end_error - end_model
                    errores[:, 0, bin_index] = 10**20 * np.array(error[start:-end])
                elif end_error < end_model:
                    diff = end_model - end_error
                    error =  error + [0] * diff
                    errores[:, 0, bin_index] = 10**20 * np.array(error[start:])

                cubo[:, 0, bin_index] = (np.array(data_in) / np.array(data_out))
        
            except KeyError: 'col3' 
     
        hdr['CDELT3'] = 1.0
        hdr['CD3_3'] = 1.0
        hdr['NAXIS3'] = len(data_out)
        hdr['CRVAL3'] = ini_model


        hdu_flux=fits.PrimaryHDU(cubo, header=hdr)
        hdu_error=fits.ImageHDU(errores, header=hdr,name='errores')
        hdu_new=fits.HDUList([hdu_flux, hdu_error])
        hdu_new.writeto('binned_stellar_specs_'+cube_name+'.fits', overwrite=True)

PARSER = argparse.ArgumentParser()
PARSER.add_argument('-N', '--Cubename', type=str, default=None)

if __name__ == '__main__':
    args = PARSER.parse_args()
    MakeGasCube(args.Cubename)
    elapsed = (time.time() - start_time)
    print(str(timedelta(seconds=elapsed)))
    fit_emission_absorption_lines_w('binned_stellar_specs_'+args.Cubename,'MUSE','',0,3)
    elapsed = (time.time() - start_time)
    print(str(timedelta(seconds=elapsed)))
