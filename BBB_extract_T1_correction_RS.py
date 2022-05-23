# BBB DCE-MRI T1 TIME COURSE CALCULATION - BBB_t1TimeCourse_noBSL_140917.py
#===============================================================================
# AUTHOR: Donnie Cameron
# START DATE: 07/09/2017
# LAST UPDATED: 13/10/2017
# REQUIRED PACKAGES: numpy, scipy, matplotlib, nibabel
#===============================================================================
'''DESCRIPTION: This script takes DCE-MRI data, DESPOT-HIFI T1w data, and brain
tissue segmentation masks to calculate regional T1 values across a time course.
Time course data are then plotted together for observation.
   Currently the first point in the DCE-MRI dataset is used as the baseline, due
to timing issues between DCE and DESPOT1 acquisitions.

Method from HEYE ET AL. (2016). Tracer Kinetic Modelling for DCE-MRI Quantifi-
cation of Subtle Blood-Brain Barrier Permeability. NeuroImage 125:446–55.

Blood T1 values from LU ET AL. (2004). Determining the longitudinal relaxation
time (T1) of blood at 3.0 Tesla. Magn Reson Med 52:679-82.

Gadovist relaxivities from PINTASKE ET AL. (2008). Relaxivity of Gadopentetate
Dimeglumine (Magnevist), Gadobutrol (Gadovist), and Gadobenate Dimeglumine
(MultiHance) in human blood plasma at 0.2, 1.5, and 3 Tesla. [ERRATUM in Invest
Radiol 41:859.] Invest Radiol 41:213-21'''

## IMPORT PACKAGES

import os
os.chdir('C:\\Users\\Rashed & Mahbuba\\Desktop\\New Laptop lenovo 25 july\\Rashed RA work\\BBB research Donnie\\BBB_sample_data_fromDonnie\\data\\20_N2144')

# Change the above once the script is finalised


import math
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
plt.style.use( 'ggplot' )
plt.ion()
from scipy.optimize import minimize, newton, basinhopping, least_squares

from datetime import datetime
startTime = datetime.now()

    ## FUNCTION DEFINITIONS
    ## function to calculate rho_calc and T1_calc

def t1ROI_fun_initial( x, alpha_pRad, alpha_PRad, TR, TR_IR, TI_IR, \
sig_SPGR, sig_IR, ir_eff, rho_scale ):
    '''This function performs linear fits to ROI-based SPGR signal / sin (flip *
    b1) vs ROI-based SPGR signal / tan (flip * b1), for a set of flip angles,
    and estimates T1 and rho (M0) from the slope and intercept of the fit.
    DCam - list input variable requirements here'''

    kappa = x # Estimated B1

    # Flip and other terms
    sina = np.sin( alpha_pRad * kappa )
    cosa = np.cos( alpha_pRad * kappa )
    tana = np.tan( alpha_pRad * kappa )

    #== Linear fit to SPGR data =======#
    X = sig_SPGR / tana
    Y = sig_SPGR / sina

    linfit = np.polyfit( X, Y, 1 )
    slope = linfit[ 0 ]
    intercept = linfit[ 1 ]

    global rho_calc  # Need to make these variables
    global T1_calc   # global to allow export

    rho_calc = intercept / ( 1.0 - slope )
    T1_calc = - TR / np.log( slope )

    #==================================#
    # Calculate spoiled GRE residual

    E1 = np.exp( - TR / T1_calc ) # Can now calculate this factor from T1_Calc

    SPGR_res = abs(rho_calc * ( 1. - E1 ) * sina / ( 1. - E1 * cosa ))

    # Calculate IR spoiled GRE residual

    # TI = TI_IR .* TI_SCALE;      # Re-scale TI.
    # full_TR_IR = TI + READOUT_PULSES * TR_IR
    IR_res = rho_scale * rho_calc * np.sin( alpha_PRad * kappa ) * \
    ( 1 + ir_eff * np.exp( - TI_IR / T1_calc ) + \
    np.exp( - TR_IR / T1_calc ) )
    IR_res = abs( IR_res )

    res = np.linalg.norm( SPGR_res - sig_SPGR ) \
    + np.linalg.norm( IR_res - sig_IR )

    # Use np.ma to mask NaNs and infs from calculated residuals
    res = np.linalg.norm( np.ma.masked_invalid( SPGR_res ) - sig_SPGR ) \
    + np.linalg.norm( np.ma.masked_invalid( IR_res ) - sig_IR )

    return res

def t1ROI_fun(x, alpha_pRad, alpha_PRad, TR, TR_IR, TI_IR, TD_IR, N, PECentre, sig_SPGR, sig_IR):
    '''This function performs linear fits to ROI-based SPGR signal / sin (flip *
    b1) vs ROI-based SPGR signal / tan (flip * b1), for a set of flip angles,
    and estimates T1 and rho (M0) from the slope and intercept of the fit.
    DCam - list input variable requirements here'''

    rho_s = x[1]
    kappa = x[2]
    T1_s  = x[0]

    # Flip and other terms

    sina = np.sin(np.asarray(alpha_pRad) * kappa)
    cosa = np.cos(np.asarray(alpha_pRad) * kappa)
    tana = np.tan(np.asarray(alpha_pRad) * kappa)



    #==================================#
    # Calculate spoiled GRE residual

    #E1 = np.exp( - TR / T1_calc ) # Can now calculate this factor from T1_Calc
    # s = S0. * (((1 - exp(-TR_s. / T1_s)). * sin(b_rad)). / (1 - exp(-TR_s. / T1_s). * cos(b_rad)))
    #SPGR_res = rho_s * (((1. - np.exp(-TR / T1_s)) * np.sina) / (1. - np.exp(-TR / T1_s) * np.cosa))
    SPGR_res = abs(rho_s * (((1. - np.exp(-TR / T1_s)) * np.sin(np.asarray(alpha_pRad) * kappa)) / (
            1. - np.exp(-TR / T1_s) * np.cos(np.asarray(alpha_pRad) * kappa))))


    # Calculate IR spoiled GRE residual

    # TI = TI_IR .* TI_SCALE;      # Re-scale TI.
    # full_TR_IR = TI + READOUT_PULSES * TR_IR

    tau = TR_IR * N  # % duration of the acquisition block

    #T1Star = ((1 / T1_s) - (1. / TR_s). * log(cos(b_rad))). ^ (-1);
    T1_s = 1./ ((1. / T1_calc) + (1. / TR_IR * np.log(np.cos(alpha_PRad * kappa))))
    #M0Star = S0 * ((1 - exp(-TR_s / T1_s)). / (1 - exp(-TR_s. / T1Star)));
    #S0=rho_s; M0Star=rho_calc; T1star=T1_calc;
    rho_s = rho_calc / ((1. - np.exp(-TR_IR / T1_s)) / (1. - np.exp(-TR_IR / T1_calc)))

    #A1 = M0Star. * (1 - exp(-tau. / T1Star));
    A1 = rho_calc * (1. - np.exp(-tau / T1_calc))
    #A2 = S0 * (1 - exp(-TD_s / T1_s));

    #A3 = S0 * (1 - exp(-TI_s / T1_s));
    #B1 = exp(-tau. / T1Star);
    #B2 = exp(-TD_s / T1_s);
    #B3 = -exp(-TI_s / T1_s);
    A2 = rho_s * (1. - np.exp(-TD_IR / T1_s))
    A3 = rho_s * (1. - np.exp(-TI_IR / T1_s))
    B1 = np.exp(-tau / T1_calc)
    B2 = np.exp(-TD_IR / T1_s)
    B3 = -np.exp(-TI_IR / T1_s)

    #A = A3 + A2. * B3 + A1. * B2. * B3;
    #B = B1. * B2. * B3;

    #M1 = A. / (1 - B);
    A = A3 + A2 * B3 + A1 * B2 * B3
    B = B1 * B2 * B3

    M1 = A / (1. - B)

    #s = (M0Star + (M1 - M0Star). * exp(-(PECentre. * tau). / T1Star)). * sin(b_rad)

    IR_res = (rho_calc + (M1 - rho_calc) * np.exp(-(PECentre * tau) / T1_calc)) * np.sin(alpha_PRad*kappa)
    IR_res = abs(IR_res)



    res = np.linalg.norm( SPGR_res - sig_SPGR ) \
    + np.linalg.norm( IR_res - sig_IR )

    # Use np.ma to mask NaNs and infs from calculated residuals
    res = np.linalg.norm( np.ma.masked_invalid( SPGR_res ) - sig_SPGR ) \
    + np.linalg.norm( np.ma.masked_invalid( IR_res ) - sig_IR )

    return res

# to calculate T1_calc, rho_calc, kappa_calc
def minT1_initial( alpha_pRad, alpha_PRad, TR, TR_IR, TI_IR, sigSPGR, sigIR, \
ir_eff, rho_scale ):
    '''This function allows simple T1 calculation (via t1ROI_fun minimisation)
    for a selection of ROIs or masks.
    DCam - list input variable requirements here'''

    x_start = ( 1.0 ) # Minimisation starting value for kappa

    param = minimize( t1ROI_fun_initial, x_start, \
    method = 'Nelder-Mead', \
    args = ( alpha_pRad, alpha_PRad, TR, TR_IR, TI_IR, \
    sigSPGR, sigIR, ir_eff, rho_scale ) )

    return param

# this function will give actual rho_s, kappa and T1_s
def minT1( alpha_pRad, alpha_PRad, TR, TR_IR, TI_IR, TD_IR, N, PECentre, sig_SPGR, sig_IR ):
    '''This function allows simple T1 calculation (via t1ROI_fun minimisation)
    for a selection of ROIs or masks.
    DCam - list input variable requirements here'''

    minT1_bounds=( ( 0.0 , None ), ( 0.0 , None) , (0.0 , None) )
    minT1_param_init = basinhopping(t1ROI_fun, (T1_calc, rho_calc , 1.0), stepsize= 2.0 ,  minimizer_kwargs={'args': (alpha_pRad, alpha_PRad, TR, TR_IR, TI_IR, TD_IR, N, PECentre, sig_SPGR, sig_IR)})

    #minT1_param = minimize(t1ROI_fun, (1.0, 1.0, 1.0), \
                     #method='Nelder-Mead', \
                     #args=(alpha_pRad, alpha_PRad, TR, TR_IR, TI_IR, TD_IR, N, PECentre, sig_SPGR, sig_IR))
    #minT1_param = minimize(t1ROI_fun, (T1_calc, rho_calc , 1.0), method='L-BFGS-B', \
                           #args=(alpha_pRad, alpha_PRad, TR, TR_IR, TI_IR, TD_IR, N, PECentre, sig_SPGR, sig_IR),
                           #bounds=minT1_bounds, options={'maxiter': 100})

    # minT1_param = minimize(t1ROI_fun, minT1_param_init.x, method = 'L-BFGS-B', \
    # args = ( alpha_pRad, alpha_PRad, TR, TR_IR, TI_IR, TD_IR, N, PECentre, sig_SPGR, sig_IR), bounds = minT1_bounds, options = {'maxiter':100} )



    minT1_param = least_squares(t1ROI_fun, [T1_calc, rho_calc , 1.0] , jac='2-point',
                                bounds=([0.0, 0.0, 0.0], [np.Inf, np.Inf, np.Inf]), method='trf', ftol=1e-08,
                                xtol=1e-08, gtol=1e-08,
                                x_scale=1.0, loss='linear', f_scale=1.0, diff_step=None, tr_solver=None, tr_options={},
                                jac_sparsity=None, max_nfev=None, verbose=0,
                                args=(alpha_pRad, alpha_PRad, TR, TR_IR, TI_IR, TD_IR, N, PECentre, sig_SPGR, sig_IR),
                                kwargs={})
    return minT1_param


def GdConcFun( Ci, sigEn, r1, r2, T1, TR, TE, alpha_b ):
    '''This function solved via SciPy's 'newton' function (Newton-Raphson or
    secant method) to determine [Gadovist] via equations from Armitage et al.
    DCam - list input variable requirements here'''

    P = TR / T1
    Q = r1 * Ci * TR
    cosab = np.cos( alpha_b )
    expP = np.exp( - P ) # -P in Armitage (2005, 2011), +P in Heyes (2016)
    # Armitage parameters give realistic (non-neg.) concentration-time curves.

    # Subtract RHS of Armitage equation from LHS to solve ...
    return sigEn - np.exp( - r2 * Ci * TE ) * \
    ( ( 1. - np.exp( - P - Q ) - cosab * ( expP - np.exp( - 2 * P - Q ) ) ) / \
    ( 1. - expP - cosab * ( np.exp( - P - Q ) - np.exp( - 2 * P - Q ) ) ) ) + 1.

def patlakFun( param, tVec, conc_plasma, conc_time, excPts ):
    '''This function is used to compute Patlak model contrast uptake parameters
    via global and local non-linear minimisation.
    DCam - list input variable requirements here'''

    v_p = param[ 0 ]      # fractional interstitial volume
    K_trans = param[ 1 ]  # volume transfer constant

    patlak = v_p * conc_plasma[ excPts : -1 ] + \
    K_trans * np.trapz( conc_plasma, x = tVec )  # Use integral of full VIF


    res = np.linalg.norm( patlak - conc_time[ excPts : -1 ] )

    return res

def patlakFit( param, tVec, conc_plasma, conc_time ):
    '''This function is used for plotting the Patlak model fit.
    DCam - list input variable requirements here'''

    v_p = param[ 0 ]      # fractional interstial volume
    K_trans = param[ 1 ]  # volume transfer constant

    # return v_p * conc_plasma + K_trans * np.trapz( conc_plasma, x = tVec )
    return v_p * conc_plasma + K_trans * np.trapz(conc_plasma, x=tVec)

## MAIN CODE LOOP
# Load spoiled GRE images (in nifti format) from the working directory
T1img1 = nib.load('despot1_5_BET_MSK.nii.gz')
T1img1_data = T1img1.get_data( )
hdr_T1img1 = T1img1.header

T1img2 = nib.load('despot1_10_BET_MSK.nii.gz')
T1img2_data = T1img2.get_data( )

T1img3 = nib.load('despot1_15_BET_MSK.nii.gz')
T1img3_data = T1img3.get_data( )

B1img = nib.load('despot1_IR5_BET_MSK.nii.gz')
B1img_data = B1img.get_data( )

# Load DCE time series data from working directory
DCEimg = nib.load('dce_BET_MSK.nii.gz')
DCEimg_data = DCEimg.get_data( )

# Load parenchyma/CSF segmentation masks from working directory
WM = nib.load('t1S_seg_2.nii.gz') # White matter
WM_mask = WM.get_data( )
GM = nib.load('t1S_seg_1.nii.gz') # Grey matter
GM_mask = GM.get_data( )
CSF = nib.load('t1S_seg_0.nii.gz') # Cerebrospinal fluid
CSF_mask = CSF.get_data( )

# TIME SERIES SIGNAL ENHANCEMENT CALCULATIONS
#==============================================
# Calculate timings of each acquisition time point. ** USER INPUT **
acq_len = 591.5  # Acquisition duration (s) (subtracted 31.5 so pt 0 = time 0)
nPts = len( DCEimg_data[ 1, 1, 1, : ] )
tVec = np.arange( 0, nPts ) * acq_len / nPts

# The vascular input function is calculated from a voxel in the superior
# sagittal sinus. Give co-ordinates for that voxel here.

ss_vox = [ 43, 34, 90 ] # Voxel positioned in superior sagittal sinus

# Calculate median signal intensities for each tissue type, for each time point
# Initialise arrays to hold median tissue SIs per image
med_WM = np.zeros( len( DCEimg_data[ 1, 1, 1, : ] ) )  # No baseline
med_GM = np.zeros( len( DCEimg_data[ 1, 1, 1, : ] ) )
med_CSF = np.zeros( len( DCEimg_data[ 1, 1, 1, : ] ) )
# No median measurement for sagittal sinus - we use one voxel only
SS = np.zeros( len( DCEimg_data[ 1, 1, 1, : ] ) )

# No median signal intensity at baseline from DESPOT1_15deg ( different TR/TE )

# Calculate median signal per time point by looping through DCE images
for i in range( 0, len( DCEimg_data[ 1, 1, 1, : ] ) ):
    DCEtPt = DCEimg_data[ :, :, :, i ]  # Split DCE data into time points
    med_WM[ i ] = np.median( DCEtPt[ np.where( WM_mask > 0 ) ] )
    med_GM[ i ] = np.median( DCEtPt[ np.where( GM_mask > 0 ) ] )
    med_CSF[ i ] = np.median( DCEtPt[ np.where( CSF_mask > 0 ) ] )
    SS[ i ] = DCEtPt[ ss_vox[ 0 ], ss_vox[ 1 ], ss_vox[ 2 ] ]

# Convert to signal enhancement from baseline (ratio, use % for plotting only)
sigEn_WM = ( med_WM - med_WM[ 0 ] ) / med_WM[ 0 ]
sigEn_GM = ( med_GM - med_GM[ 0 ] ) / med_GM[ 0 ]
sigEn_CSF = ( med_CSF - med_CSF[ 0 ] ) / med_CSF[ 0 ]
sigEn_SS = ( SS - SS[ 0 ] ) / SS[ 0 ]

'''DCam. The first few participants have different TR/TE for DCE-MRI. Will use
t = 0 DCE timepoint as baseline instead ...'''

# Plotting code - shows signal intensity time-curves for WM, GM, CSF, SS
fig = plt.figure( )
fig.suptitle( 'Signal Enhancement Curves', fontsize = 14 )
ax1 = plt.gca( )
plt1, = plt.plot( tVec, sigEn_WM * 100, '#56B4E9', label = 'White Matter', \
linewidth = 2 )
plt2, = plt.plot( tVec, sigEn_GM * 100, '#0072B2', label = 'Grey Matter' , \
linewidth = 2 )
#plt3, = plt.plot( tVec, sigEn_CSF * 100, '#CC79A7', label = 'CSF' , \
#linewidth = 2 )
ax1.set_xlabel( 'Time (s)' )
ax1.set_ylabel( 'Signal Enhancement (%)' )
ax1.grid(None)
# Put SS on similar scale
ax2 = ax1.twinx( )
plt4, = plt.plot( tVec, sigEn_SS * 100, '#D55E00', label = 'Sag. Sinus' , \
linewidth = 2 )
ax2.set_ylabel( 'SSS Signal Enhancement (%)', color = '#D55E00' )
ax2.tick_params( 'y', colors = '#D55E00' )
ax2.set_xlim( [ 0, np.max( tVec ) + 20 ] )
ax2.grid(None)
plt.show( )
plt.legend( handles = [ plt1, plt2, plt4 ], loc = 4, frameon = False )

plt.savefig('sigEn_plot.png', format = 'png', bbox_inches = 'tight', dpi = 500)

# BASELINE T1 CALCULATIONS
#===========================
# Adapted from pixel-by-pixel DESPOT1 mapping code 'BBB_t1maps_MAIN_070917'
ir_eff = - 1 + np.cos( 0.97 * math.pi) # Inversion effiency (constant)
rho_scale = 0.975;                     # Scale of IR-SPGR M0 relative to SPGR

alpha_p = [ 5, 10, 15 ] # List of prescribed flip angles for spGRE (deg).
alpha_P = 5             # List of prescribed flip angles for IR-spGRE (deg).
TI_IR = 0.45            # Inversion time (s)
TR = 0.005568           # Repetition time of spoiled GRE sequence (s)
TE = 0.001508           # Echo time of spoiled GRE sequence (s)
TR_IR = 0.8             # Time between successive inversion pulses (s)
TD_IR=0
N= 96
PECentre=0

alpha_pRad = [ 0 ] * len( alpha_p ) # Initialise list of flips in rads

for fa in range( 0, len( alpha_p ) ):
    alpha_pRad[ fa ] = alpha_p[ fa ] * math.pi / 180


alpha_PRad = alpha_P * math.pi / 180    # IR flip angle in rads

medWM_sigSPGR = [ np.median( T1img1_data[ np.where( WM_mask > 0 ) ] ), \
np.median( T1img2_data[ np.where( WM_mask > 0 ) ] ), \
med_WM[ 0 ] ] # Arrange WM SPGR data and redefine med_WM[ 0 ] for clarity
medWM_sigIR = np.median( B1img_data[ np.where( WM_mask > 0 ) ] )

medGM_sigSPGR = [ np.median( T1img1_data[ np.where( GM_mask > 0 ) ] ), \
np.median( T1img2_data[ np.where( GM_mask > 0 ) ] ), \
med_GM[ 0 ] ] # Arrange GM SPGR data and redefine med_GM[ 0 ] for clarity
medGM_sigIR = np.median( B1img_data[ np.where( GM_mask > 0 ) ] )

medCSF_sigSPGR = [ np.median( T1img1_data[ np.where( CSF_mask > 0 ) ] ), \
np.median( T1img2_data[ np.where( CSF_mask > 0 ) ] ), \
med_CSF[ 0 ] ] # Arrange CSF SPGR data and redefine med_CSF[ 0 ] for clarity
medCSF_sigIR = np.median( B1img_data[ np.where( CSF_mask > 0 ) ] )

''' DCam, 110917. Need to make code more modular - it should accept a variable
number of arguments (tissue masks).'''

# Minimisation - white matter
HIFIparam_initial = minT1_initial( alpha_pRad, alpha_PRad, TR, TR_IR, TI_IR, \
medWM_sigSPGR, medWM_sigIR, ir_eff, rho_scale )

RHO_WM_calc = rho_calc
KAPPA_WM_calc = HIFIparam_initial.x[ 0 ] * 100
T1_WM_calc = T1_calc

# Minimisation - grey matter
HIFIparam_initial = minT1_initial( alpha_pRad, alpha_PRad, TR, TR_IR, TI_IR, \
medGM_sigSPGR, medGM_sigIR, ir_eff, rho_scale )

RHO_GM_calc = rho_calc
KAPPA_GM_calc = HIFIparam_initial.x[ 0 ] * 100
T1_GM_calc = T1_calc

# Minimisation - CSF
HIFIparam_initial = minT1_initial( alpha_pRad, alpha_PRad, TR, TR_IR, TI_IR, \
medCSF_sigSPGR, medCSF_sigIR, ir_eff, rho_scale )

RHO_CSF_calc = rho_calc
KAPPA_CSF_calc = HIFIparam_initial.x[ 0 ] * 100
T1_CSF_calc = T1_calc


# Minimisation - white matter
HIFIparam = minT1( alpha_pRad, alpha_PRad, TR, TR_IR, TI_IR, \
TD_IR, N, PECentre, medWM_sigSPGR, medWM_sigIR )

T1_WM = HIFIparam.x[ 0 ]
RHO_WM = HIFIparam.x[ 1 ]
KAPPA_WM = (HIFIparam.x[ 2 ]) * 100

# Minimisation - grey matter
HIFIparam = minT1( alpha_pRad, alpha_PRad, TR, TR_IR, TI_IR, \
TD_IR, N, PECentre, medGM_sigSPGR, medGM_sigIR)

T1_GM = HIFIparam.x[ 0 ]
RHO_GM = HIFIparam.x[ 1 ]
KAPPA_GM = (HIFIparam.x[ 2 ]) * 100


# Minimisation - CSF
HIFIparam = minT1( alpha_pRad, alpha_PRad, TR, TR_IR, TI_IR, \
TD_IR, N, PECentre, medCSF_sigSPGR, medCSF_sigIR )


T1_CSF = HIFIparam.x[ 0 ]
RHO_CSF = HIFIparam.x[ 1 ]
KAPPA_CSF = (HIFIparam.x[ 2 ]) * 100

'''DCam, 110917. Try to exclude more non-brain voxels by eroding masks. TEST'''

# CONTRAST AGENT CONCENTRATION CALCULATION
#==========================================
'''This is where the ROI-based contrast agent concentration is calculated.
PINTASKE 2006 (Erratum) gives r1 = 4.5 and r2 = 6.3 ... This r1 is similar
to the whole-blood r1 given by SHEN (2015). Assume similar r2'''
# Constant variables:
r1 = 4.5  # Gadovist r2 in human blood plasma, from Pintaske (err - 2008)
r2 = 6.3  # Gadovist r2 in human blood plasma, from Pintaske (err - 2008)
T1_Blood = 1.584  # Venous blood T1, taken  from Lu et al.

excPts = 5 # Number of points to exclude from Patlak model fitting

# Initialise Gd concentration arrays for time course
WM_GdConc = np.zeros( len( DCEimg_data[ 1, 1, 1, : ] ) )  # No baseline
GM_GdConc = np.zeros( len( DCEimg_data[ 1, 1, 1, : ] ) )
CSF_GdConc = np.zeros( len( DCEimg_data[ 1, 1, 1, : ] ) )
SS_GdConc = np.zeros( len( DCEimg_data[ 1, 1, 1, : ] ) )

C = 0.02   # Initial guess for contrast conc. (mM). (Should be ~ 0.02 - 2)
#this was 0.8 in DC code, armitage suggests 0.02 for GM and 0.009 for WM

# TR and TE for DCE-MRI sequence (different from DESPOT1)
#TR_CE = 0.006449; Value in DC code
#TE_CE = 0.002412; value in DC code

TR_CE= 0.006379
TE_CE= 0.002468

# Correct for B1 errors in flip angles using KAPPA from DESPOT-HIFI
corrAlpha_WM = ( KAPPA_WM / 100 ) * alpha_pRad[ -1 ]
corrAlpha_GM = ( KAPPA_GM / 100 ) * alpha_pRad[ -1 ]
corrAlpha_CSF = ( KAPPA_CSF / 100 ) * alpha_pRad[ -1 ]

# Loop over DCE-MRI time points
for i in range( 1 , len( DCEimg_data[ 1, 1, 1, : ] ) ):

    conc = newton( GdConcFun, C, \
    args = ( sigEn_WM[ i ], r1, r2, T1_WM, TR_CE, TE_CE, corrAlpha_WM ), maxiter=5000 )
    WM_GdConc[ i ] = conc

    conc = newton( GdConcFun, C, \
    args = ( sigEn_GM[ i ], r1, r2, T1_GM, TR_CE, TE_CE, corrAlpha_GM ) , maxiter=5000)
    GM_GdConc[ i ] = conc

    conc = newton( GdConcFun, C, \
    args = ( sigEn_CSF[ i ], r1, r2, T1_CSF, TR_CE, TE_CE, corrAlpha_CSF ), maxiter=5000)
    CSF_GdConc[ i ] = conc

    # Assume optimal flip angle in SS
    conc = newton( GdConcFun, C, \
    args = ( sigEn_SS[ i ], r1, r2, T1_Blood, TR_CE, TE_CE, alpha_pRad[ -1 ] ), maxiter=2000 )
    SS_GdConc[ i ] = conc

# FITTING PATLAK FUNCTION
#========================
'''Here we fit the Patlak function to the concentration time curves. The first
step in curve-fitting is to apply the basin-hopping minimisation approach to 
find an approximate global minimum, then L-BFGS-B for precise estimation.'''
# Constant variables:
HCt = 0.423                 # Mean haematocrit measurement from all subjects
conc_plasma = SS_GdConc / ( 1 - HCt ) # Calculate [Gd_plasma] from [Gd_blood]
'''DCam, 131017. Use per-participant haematocrit to calculate plasma [Gd]'''
vp_est = 0.2     # v_p starting value for minimisation. (Heye gives ~ 0.008)
Kt_est = 0.0004  # K_trans starting value for minimisation.
bounds = ( ( 0, 1 ), ( 0, None ) ) # Restrict to positive values

# Basin-Hopping global minimisation
WM_param_init = basinhopping( patlakFun, [ vp_est, Kt_est ], \
minimizer_kwargs = { 'args': ( tVec, conc_plasma, WM_GdConc, excPts )} )

GM_param_init = basinhopping( patlakFun, [ vp_est, Kt_est ], \
minimizer_kwargs = { 'args': ( tVec, conc_plasma, GM_GdConc, excPts )} )

CSF_param_init = basinhopping( patlakFun, [ vp_est, Kt_est ], \
minimizer_kwargs = { 'args': ( tVec, conc_plasma, CSF_GdConc, excPts )} )

# Local minimisation
WM_param = minimize( patlakFun, WM_param_init.x, method = 'L-BFGS-B', \
args = ( tVec, conc_plasma, WM_GdConc, excPts ), bounds = bounds )
print( WM_param.x )

GM_param = minimize( patlakFun, GM_param_init.x, method = 'L-BFGS-B',\
args = ( tVec, conc_plasma, GM_GdConc, excPts ), bounds = bounds )
print( GM_param.x )

CSF_param  = minimize( patlakFun, CSF_param_init.x, method = 'L-BFGS-B',\
args = ( tVec, conc_plasma, CSF_GdConc, excPts ), bounds = bounds )
print( CSF_param.x )

# Plotting - shows fitted Gad concentration time-curves for WM, GM, CSF, SS
fig2 = plt.figure( )
fig2.suptitle( 'Concentration-Time Curves', fontsize = 14 )
ax1 = plt.gca( )
plt1, = plt.plot( tVec, patlakFit( WM_param.x, tVec, conc_plasma, WM_GdConc ), \
'#56B4E9', label = 'WM Patlak Fit',linewidth = 2 )
plt2 = plt.scatter( tVec, WM_GdConc, s = 20, c = '#56B4E9', \
label = 'White Matter', linewidths = 0 )
plt3, = plt.plot( tVec, patlakFit( GM_param.x, tVec, conc_plasma, GM_GdConc ), \
'#0072B2', label = 'GM Patlak Fit', linewidth = 2 )
plt4 = plt.scatter( tVec, GM_GdConc, s = 20, c = '#0072B2', \
label = 'Grey Matter', linewidths = 0 )
#plt5, = plt.plot( tVec, patlakFit( CSF_param.x, tVec, conc_plasma, CSF_GdConc ),\
#'#CC79A7', linewidth = 2 )
#plt6 = plt.scatter( tVec, CSF_GdConc, s = 20, c = '#CC79A7', \
#label = 'CSF', linewidths = 0 )
ax1.set_xlabel( 'Time (s)' )
ax1.set_ylabel( 'Contrast Agent Concentration (mM)' )
ylim1 = ax1.get_ylim( )
#ax1.set_ylim( [ 0, ylim1[ 1 ] ] )
ax1.grid(None)
# Put SS on similar scale
ax2 = ax1.twinx( )
plt7, = plt.plot( tVec, SS_GdConc,'#D55E00', ls = '--', label = 'Sag. Sinus', \
linewidth = 2 )
ax2.set_ylabel( 'SSS Contrast Agent Concentration (mM)', color = '#D55E00' )
ylim2 = ax2.get_ylim( )
#ax2.set_ylim( [ 0, ylim2[ 1 ] ] )
ax2.tick_params( 'y', colors = '#D55E00' )
ax2.set_xlim( [ 0, np.max( tVec ) + 20 ] )
ax2.grid(None)
plt.show( )
plt.legend( handles = [ plt1, plt2, plt3, plt4, plt7 ], loc = 1, ncol = 2, \
frameon = False )

plt.savefig('GdConc_plot.png', format = 'png', bbox_inches = 'tight', dpi = 500)

print( datetime.now() - startTime )

GM_fit=(patlakFit( GM_param.x, tVec, conc_plasma, GM_GdConc ))
WM_fit=(patlakFit( WM_param.x, tVec, conc_plasma, GM_GdConc ))


















































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































