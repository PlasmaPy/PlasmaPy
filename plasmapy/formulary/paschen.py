#@title
"""r

Paschen curve


The Paschen law calculates the breakdown voltage V in function of pressure times distance(pd) product in a self-sustained small gap DC Townsend discharges

V=Bpd/ln(Apd/ln(1+1/γ)) 

respectively the breakdown field

(E/p)=(V/d)∗(1/p)=B/ln(Apd/ln(1+1/γ)) 

Where:

V Breakdown voltage (V)

p Gas pressure (Pa)

d Electrode distance (m)

A Material constant (ionisations/(Pa m)

B Material constant (V/(Pa m)

γ  Second Townsend coefficient


The breakdown voltage and breakdown field depend on the gas, the electrode material, the pressure and discharge gap used.

The coefficients A and B depend on the nature of the gas used and they are only valid in some specific regions of the reduced field (E/p)(see later).



Major assumptions:

1D
Electrons produced at cathode due to ion impact only.
E field constant as well as the coefficients A and B


Deviations from the Paschen curve can originate from:

Electron attachment

Electrode material

Electrical field distribution (inhomogeneous field)

Space charge effects

Impurities in the gas

Small electrode gap

……
IMPORTANT

The given Paschen law is applicable only for small gap Townsend DC discharges.

Transition from the Townsend to Streamer-leader regime! If the E field due to space charge is about the applied field then transition to streamer discharge is observed (Pd>200 Torr cm =266 Pa m ; Raizer)

Streamers occur when the product of the ionization coefficient  α , and the gap distance d, exceeds the Meek’s criterion. Rule of thumb  α d = 18-21 for air according to Reather and Meek(ref) ).

For high frequency discharges (microwave, RF) deviation from the above DC Paschen law are observed and the DC Paschen law is stricely not applicable for RF. Calculation of breakdown voltages for high frequency discharges will be discussed in a future notebook.


"""
import numpy as np



def get_paschen_constants (gas,electrode): 

# Function to get the necessary constants A anb B to calculate the Paschen breakdown voltage

# Supported gases 

    gases=["Air","N2","H2","He","Ne","Ar","Kr","Xe"]

    
#  paschen_constants contains the coefficents A and B  for the estimation of the
#       First Townsend Ionization Coefficent
#       (exponential fit to the First Townsend Ionization coefficient)
#       as adapted from
#       E.Nasser, Fundamentals of Gaseous Ionization and Plasma Electronics,
#       Wiley-Interscience, New York 1971
#    
#  format: paschen_constants dir {"gas":[A,B]}
#  units: A in [Ionisation/(Pa m)] and B in [V/(Pa m)]

    paschen_constants={"Air":[11,274],
         "N2":[9.0, 257],
         "H2":[3.8,104],
         "He":[2.3,26],
         "Ne":[3.0, 75],
         "Ar":[11,135],
         "Kr":[13,180],
         "Xe":[20,263]}


# Supported electrode materials

    material=["Al","Cu","Ni","Pt","C","W","Fe"]
    
#  Townsend_gamma is the Second Townsend Ionization coefficient as given by
#       A.Beroual and I. Fonfana, Discharge in Long Air Gap Modeling and Application
#       IOP Publishing Ltd 2016
#       ISBN 978-0-7503-1236-3 (ebook)
#       ISBN 978-0-7503-1237-0 (print)
#
#
    
    
    townsend_gamma={"Air":{"Al":0.035,"Cu":0.025,"Ni":0.036,"Pt":0.017,"C":None,"W":None,"Fe":0.02},
        "N2":{"Al":0.1,"Cu":0.066,"Ni":0.077,"Pt":0.59,"C":None,"W":None,"Fe":0.059},
        "H2":{"Al":0.095,"Cu":0.05,"Ni":0.053,"Pt":0.02,"C":0.014,"W":None,"Fe":0.061},
        "He":{"Al":0.021,"Cu":None,"Ni":0.015,"Pt":0.01,"C":None,"W":None,"Fe":0.015},
        "Ne":{"Al":0.053,"Cu":0.02,"Ni":0.031,"Pt":0.023,"C":None,"W":0.045,"Fe":0.022},
        "Ar":{"Al":0.12,"Cu":0.058,"Ni":0.058,"Pt":0.058,"C":None,"W":None,"Fe":0.058},
        "Kr":{"Al":None,"Cu":None,"Ni":None,"Pt":None,"C":None,"W":None,"Fe":None},
        "Xe":{"Al":None,"Cu":None,"Ni":None,"Pt":None,"C":None,"W":None,"Fe":None}}
 
 # Check if the asked gas and electrode material is supported
    resg=gas in gases
    rese=electrode in material
 # If the gas is supported get the constants A and B   
    print(resg,rese)
    if resg==True :
        print(gas)  
        A=paschen_constants[gas][0]
        B=paschen_constants[gas][1]
        print(A,B)
    
# get the townsend_gamma coefficient for the the gas/electrode combination
        if rese==True:
            gam=townsend_gamma[gas]
            print(gam)
            gn=gam[electrode]
            print (gn)
  # test if townsend_gamma exists (None) for the demanded gas/electrode configuration 
  # if not a default townsend_gamma value of 0.01 is taken
        
            if gn is None:
                gn=0.01           
                print("default")
                print(gn)
        else:
# if the electrode material is not supportes set townsend_gamma to default = 0.01
            gn=0.01
            print("default")
# create output dir {const)        
        const={"A":A,"B":B,"gam":gn}
        print(const)
        return const
# if gas is not supported set const=None
    else :
        const=None
        print("No constants for this gas available",const)
        return const



def breakdown_voltage(distance,pressure,A,B,gam):
# calculate the breakdown voltage according to the Paschen law
# pressure input can be a list (since in th Paschen law pd is the important paramter
# if several distance at const pressure are asked for just inverted the two values
    x=list(distance *p for p in pressure)
    g=np.log(1+(1/gam))
    vb=list((pd,B*pd/(np.log((A*pd)/g))) for pd in x)
    return vb


def minimum_breakdown_voltage(A,B,gam):
# calculate the minimum breakdown voltage [V] and the corresponding pd [Pa m] value
    g=np.log(1+(1/gam))
    vmin=2.718*(B/A)*g
    pdmin=2.718*(g/A)
    return (vmin,pdmin)
    
  
# Some example of Paschen law calculations
# get Paschen constants
c=get_paschen_constants(gas="N2",electrode="Ni")
print("Paschen constants", c)

# calculate Paschen law
# distance [m] and pressure [Pa]
d=1
p=[1,2]
if c==None:
# calculate breakdown voltage for the given Paschen constants     
    vb1=breakdown_voltage(0.1,[100,200],11,135,0.058)
    print(vb1)
else:
# calculate Paschen law using in-built constants
    vb2=breakdown_voltage(d,p,**c)
    print("Breakdown voltage",vb2)
    min=minimum_breakdown_voltage(**c)
    print ("minimum breakdown voltage",min[0],"p*d at minimum breakdown voltage",min[1])

