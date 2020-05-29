# -*- coding: utf-8 -*-
"""
Created on Wed May 13 15:20:10 2020
@author: Nick

"Phs_4". It is the subroutine calculating the cross section for a hard sphere ("hs") 
carrying an electric charge ("4"). 
"P" is for "potential": the interaction potential between a hard sphere carrying a 
charge and a polarizable particle. The radius dependence of a charge interaction with an 
induced dipole is to the power of -4. Hence "4". 

P16, 12, 8 ... are the collision integrals from:

Transport Properties of Ions in Gases by Edward A. Mason Earl W. McDaniel 

These are all just interpolating around the tabulated cross sections and 
a bit of a parameterization of values outside tabulated data are applied.

The tabulated data are available open access from the publisher.

"""

import math
from scipy import interpolate

def phs4(reducedTemperature):
#   Reduced collision integral for a hard sphere - 4 potential
    Xt = [0.01,0.04,0.09,0.16,0.25,0.36,0.49,0.64,0.81,
          1.00,1.21,1.44,1.69,1.96,2.25,2.56,2.89,3.24,3.61,4.00,4.84,
          5.76,6.76,7.84,9.00,10.24,11.56,12.96,14.44,16.00,19.36,
          23.04,27.04,31.36,36.00,40.96,46.24,51.84,57.76,64.00,
          100.00,1000.00]    
    YtG = [13.67,6.640,4.343,3.213,2.548,2.117,1.823,1.617,
           1.472,1.368,1.292,1.236,1.194,1.161,1.136,1.116,1.100,
           1.087,1.076,1.067,1.054,1.044,1.036,1.031,1.026,1.023,
           1.020,1.017,1.015,1.014,1.010,1.005,1.000,1.000,1.000,
           1.000,1.000,1.000,1.000,1.000,1.000,1.000]
#sets up to see if the temperature is even in range It might not be neccesary with pchip interpolation,
#but it will at least catch some odd errors coming through.
    if reducedTemperature > Xt[41]: 
        q=1.00
        return float(q)
    if reducedTemperature < Xt[0]:
        q=1.4691/math.sqrt(reducedTemperature)
        return float(q)
#interpolates using scipy "pchip" tested to give best results with this data set.  
    q = interpolate.pchip_interpolate(Xt, YtG, reducedTemperature)    
#this is a numpy object so you have to convert it to something useful to the rest of the program.
    return float(q)


#potential for 16,6,4
def P16_6_4(t,g):
    
    Xt = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,
         0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.50,0.6,
         0.70,0.8,0.90,1.,1.2,1.4,1.6,1.8,2.,2.5,3.,3.5,
         4.,4.5,5.0,6.,8.,10.,20.,40.,60.,80.,100.,200.,
         400.,600.,800.,1000.]

    Xg = [0.,0.2,0.4,0.6,0.8,1.0]

    #an Xt x Xg matrix 
    Ygt = [[16.093,14.284,12.383,10.418,8.3866,6.3605],
           [11.255,10.036,8.7789,7.5146,6.2526,5.0363],
           [9.1204,8.1620,7.1910,6.2285,5.2871,4.3946],
           [7.8507,7.0484,6.2472,5.4609,7.7023,3.9901],
           [6.9858,6.2902,5.6042,4.9360,4.2978,3.7025],
           [6.3486,5.7316,5.1301,4.5475,3.9955,3.4830],
           [5.8539,5.2985,4.7616,4.2447,3.7578,3.3076],
           [5.4547,4.9494,4.4650,3.9999,3.5642,3.1626],
           [5.1267,4.6601,4.2187,3.7966,3.4024,3.0398],
           [4.8486,4.4168,4.0099,3.6235,3.2640,2.9339],
           [3.8978,3.5865,3.2998,3.0293,2.7820,2.5561],
           [3.3136,3.0799,2.8643,2.6633,2.4782,2.3108],
           [2.9011,2.7214,2.5556,2.4006,2.2582,2.1284],
           [2.5888,2.4485,2.3188,2.1973,2.0848,1.9822],
           [2.3430,2.2323,2.1292,2.0324,1.9428,1.8608],
           [2.1452,2.0565,1.9739,1.8958,1.8233,1.7569],
           [1.9834,1.9113,1.8442,1.7807,1.7214,1.6671],
           [1.8489,1.7892,1.7345,1.6823,1.6333,1.5888],
           [1.6362,1.5973,1.5606,1.5239,1.4903,1.4592],
           [1.4820,1.4546,1.4271,1.4034,1.3795,1.3571],
           [1.3644,1.3446,1.3262,1.3068,1.2908,1.2756],
           [1.2717,1.2584,1.2448,1.2324,1.2198,1.2081],
           [1.1978,1.1888,1.1796,1.1701,1.1616,1.1535],
           [1.0881,1.0851,1.0801,1.0762,1.0720,1.0680],
           [1.0106,1.0102,1.0096,1.0081,1.0067,1.0054],
           [.95331,.95523,.95629,.95723,.95777,.95760],
           [.90908,.91240,.91536,.91776,.91894,.92035],
           [.87415,.87906,.88236,.88553,.88824,.89007],
           [.81243,.81891,.82396,.82794,.83175,.83508],
           [.77127,.77879,.78492,.78992,.79423,.79779],
           [.74175,.74981,.75645,.76205,.76679,.77098],
           [.71939,.72776,.73468,.74057,.74564,.75010],
           [.70159,.71027,.71738,.72343,.72869,.73332],
           [.68696,.69584,.70318,.70935,.71470,.71944],
           [.66395,.67322,.68079,.68722,.69273,.69762],
           [.63286,.64214,.64997,.65659,.66232,.66738],
           [.61159,.62102,.62887,.63557,.64137,.64649],
           [.55536,.56478,.57266,.57940,.58530,.59042],
           [.50763,.51680,.52448,.53105,.53682,.54183],
           [.48217,.49110,.49860,.50502,.51062,.51556],
           [.46499,.47374,.48108,.48738,.49283,.49771],
           [.45214,.46073,.46794,.47412,.47946,.48429],
           [.41459,.42263,.42939,.43520,.44019,.44476],
           [.38027,.38773,.39402,.39942,.40417,.40833],
           [.36152,.36866,.37466,.37983,.38435,.38836],
           [.34877,.35568,.36149,.36650,.37089,.37476],
           [.33918,.34592,.35159,.35647,.36074,.36453]]


    if (t.real > Xt[46]):
        q=0.9229*(16.*(1.+g)/(16.*(3.+g)-12*(1.+g)))**0.125/(t**0.125)
        return q      
    if (t.real < Xt[0]):
        q=1.4714*(48.*(1.-g)/(16.*(3.+g)-12.*(1.+g)))**0.5/(t**0.5)
        return q      
    if ((g.real > Xg[5]) or (g.real < Xg[0])):
        print(' Gamma out of Range for the 16-6-4 Potential ! ')
        q=0
        return q

    function_16_6_4 = interpolate.interp2d(Xg, Xt, Ygt, kind='cubic')   
    q = function_16_6_4(g, t)

    return float(q)

def P12_6_4(t,g):

    Xt = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,
          0.1,0.15,0.2,0.25,
          0.3,0.35,0.4,0.45,0.50,0.6,0.70,0.8,0.90,1.,1.2,1.4,1.6,1.8,
          2.,2.5,3.,3.5,4.,4.5,5.0,6.,8.,10.,20.,40.,
          60.,80.,100.,200.,400.,600.,800.,1000.]

    Xg = [0.,0.2,0.4,0.6,0.8,1.0]

    Ygt = [[17.189,15.412,13.501,11.437,9.2149,6.8973],
           [12.040,10.842,9.5800,8.2499,6.8636,5.4699],
           [9.7638,8.8259,7.8493,6.8367,5.8002,4.7765],
           [8.4086,7.6262,6.8189,5.9916,5.1555,4.3383],
           [7.7849,6.8084,6.1159,5.4126,4.7089,4.0257],
           [6.8043,6.2056,5.5967,4.9834,4.3746,3.7864],
           [6.2735,5.7368,5.1935,4.6484,4.1114,3.5945],
           [5.8512,5.3580,4.8672,4.3772,3.8969,3.4355],
           [5.4971,5.0480,4.5961,4.1508,3.7169,3.3007],
           [5.1997,4.7837,4.3691,3.9583,3.5623,3.1837],
           [4.1716,3.8736,3.5805,3.2966,3.0223,2.7634],
           [3.5272,3.3069,3.0899,2.8786,2.6774,2.4849],
           [3.0654,2.8996,2.7355,2.5748,2.4211,2.2749],
           [2.7135,2.5878,2.4612,2.3370,2.2176,2.1034],
           [2.4367,2.3405,2.2419,2.1442,2.0500,1.9598],
           [2.2167,2.1400,2.0617,1.9843,1.9094,1.8373],
           [2.0347,1.9763,1.9123,1.8503,1.7897,1.7319],
           [1.8814,1.8381,1.7875,1.7365,1.6875,1.6399],
           [1.6512,1.6222,1.5873,1.5557,1.5222,1.4895],
           [1.4823,1.4626,1.4404,1.4155,1.3948,1.3734],
           [1.3528,1.3420,1.3258,1.3107,1.2948,1.2786],
           [1.2533,1.2461,1.2365,1.2254,1.2149,1.2046],
           [1.1753,1.1707,1.1636,1.1572,1.1498,1.1426],
           [1.0557,1.0567,1.0555,1.0530,1.0504,1.0482],
           [.97301,.97667,.97769,.97887,.97922,.97860],
           [.91159,.91732,.92070,.92330,.92519,.92717],
           [.86473,.87137,.87630,.87980,.88342,.88613],
           [.82778,.83471,.84082,.84515,.84949,.85350],
           [.76182,.77066,.77785,.78383,.78877,.79368],
           [.71782,.72750,.73548,.74248,.74860,.75400],
           [.68645,.69642,.70485,.71222,.71884,.72478],
           [.66248,.67274,.68145,.68907,.69590,.70209],
           [.64327,.65380,.66278,.67062,.67756,.68389],
           [.62742,.63819,.64733,.65535,.66247,.66886],
           [.60267,.61363,.62301,.63120,.63854,.64515],
           [.56877,.57986,.58934,.59776,.60526,.61203],
           [.54553,.55666,.56626,.57471,.58225,.58908],
           [.48364,.49471,.50427,.51271,.52025,.52710],
           [.43117,.44187,.45114,.45932,.46665,.47330],
           [.40339,.41379,.42279,.43075,.43788,.44435],
           [.38480,.39494,.40372,.41149,.41846,.42478],
           [.37097,.38090,.38950,.39710,.40394,.41012],
           [.33109,.34028,.34852,.35531,.36166,.36741],
           [.29540,.30382,.31113,.31761,.32345,.32873],
           [.27627,.28425,.29118,.29731,.30284,.30786],
           [.26342,.27109,.27775,.28365,.28897,.29379],
           [.25386,.26129,.26775,.27347,.27862,.28329]]

    if (t > Xt[46]):
        q=0.9022*(12*(1.+g)/(12.*(3.+g)-12*(1.+g)))**0.1666/(t**0.1666)
        return q
    if (t < Xt[0]):
        q=1.4714*(36.*(1.-g)/(12.*(3.+g)-12.*(1.+g)))**0.5/(t**0.5)
        return q
    if ((g > Xg[5]) or (g < Xg[0])):
        print(' Gamma out of Range for the 12-6-4 Potential ! ')
        q=0.
        return q
    
    function_12_6_4 = interpolate.interp2d(Xg, Xt, Ygt, kind='cubic')   
    q = function_12_6_4(g, t)   
    
    return float(q)

def P8_6_4(t,g):

    Xt = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,
          0.1,0.15,0.2,0.25,
          0.3,0.35,0.4,0.45,0.50,0.6,0.70,0.8,0.90,1.,1.2,1.4,1.6,1.8,
          2.,2.5,3.,3.5,4.,4.5,5.0,6.,8.,10.,20.,40.,
          60.,80.,100.,200.,400.,600.,800.,1000.]

    Xg = [0.,0.2,0.4,0.6,0.8,1.0]

    Ygt = [[20.157,18.726,16.992,14.844,12.123,8.6628],
           [14.143,13.186,12.043,10.652,8.9385,6.8378],
           [11.475,10.729,9.8441,8.7817,7.4943,5.9445],
           [9.8838,9.2611,8.5293,7.6590,6.6162,5.3759],
           [8.7970,8.2595,7.6303,6.8876,6.0065,4.9681],
           [7.9968,7.5167,6.9631,6.3150,5.5501,4.6545],
           [7.3733,6.9457,6.4463,5.8649,5.1897,4.4027],
           [6.8748,6.4783,6.0272,5.5055,4.8943,4.1926],
           [6.4555,6.0940,5.6796,5.2017,4.6511,4.0134],
           [6.1001,5.7656,5.3825,4.9438,4.4403,3.8601],
           [4.8459,4.6108,4.3437,4.0397,3.6925,3.2971],
           [4.0340,3.8649,3.6714,3.4523,3.2011,2.9147],
           [3.4475,3.3240,3.1791,3.0191,2.8324,2.6193],
           [3.0039,2.9114,2.8024,2.6820,2.5412,2.3791],
           [2.6589,2.5910,2.5046,2.4132,2.3049,2.1805],
           [2.3832,2.3307,2.2663,2.1966,2.1104,2.0129],
           [2.1639,2.1219,2.0691,2.0159,1.9499,1.8715],
           [1.9832,1.9494,1.9085,1.8620,1.8124,1.7509],
           [1.7055,1.6856,1.6579,1.6309,1.5971,1.5576],
           [1.5058,1.4929,1.4747,1.4592,1.4359,1.4105],
           [1.3563,1.3493,1.3366,1.3270,1.3136,1.2958],
           [1.2407,1.2369,1.2307,1.2254,1.2158,1.2058],
           [1.1501,1.1488,1.1449,1.1435,1.1387,1.1321],
           [1.0158,1.0187,1.0195,1.0212,1.0215,1.0217],
           [.92242,.92714,.93104,.93579,.93906,.94186],
           [.85449,.86064,.86578,.87192,.87767,.88272],
           [.80239,.80971,.81617,.82266,.82993,.83710],
           [.76122,.76932,.77682,.78482,.79174,.80015],
           [.68849,.69757,.70637,.71570,.72495,.73444],
           [.64058,.65035,.65979,.66960,.67960,.68996],
           [.60601,.61622,.62625,.63647,.64681,.65760],
           [.57962,.59010,.60043,.61100,.62172,.63277],
           [.55861,.56928,.57978,.59060,.60156,.61291],
           [.54134,.55214,.56281,.57375,.58488,.59640],
           [.51426,.52518,.53608,.54708,.55851,.57026],
           [.47700,.48815,.49930,.51054,.52200,.53381],
           [.45141,.46263,.47392,.48520,.49674,.50868],
           [.38377,.39493,.40615,.41734,.42877,.44061],
           [.32743,.33822,.34893,.35972,.37073,.38213],
           [.29820,.30867,.31898,.32942,.34007,.35110],
           [.27890,.28912,.29912,.30928,.31964,.33036],
           [.26472,.27471,.28448,.29439,.30450,.31497],
           [.22469,.23392,.24295,.25204,.26131,.27090],
           [.19026,.19865,.20691,.21514,.22351,.23216],
           [.17246,.18036,.18813,.19585,.20369,.21180],
           [.16080,.16834,.17577,.18312,.19060,.19832],
           [.15228,.15954,.16669,.17377,.18096,.18840]] 

    if (t > Xt[46]):
        q=0.8661*(8.*(1.+g)/(8.*(3.+g)-12*(1.+g)))**0.25/(t**0.25)
        return q

    if (t < Xt[0]):
        q=1.4714*(24.*(1.-g)/(8.*(3.+g)-12.*(1.+g)))**0.5/(t**0.5)
        return q

    if ((g > Xg[5]) or (g < Xg[0])):
        print(' Gamma out of Range for the 8-6-4 Potential ! ')
        q=0.
        return q

    function_8_6_4 = interpolate.interp2d(Xg, Xt, Ygt, kind='cubic')   
    q = function_8_6_4(g, t)   
    
    return float(q)
