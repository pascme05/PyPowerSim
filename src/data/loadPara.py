#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         loadPara
# Date:         08.02.2026
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.1
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function loads the parameters for the switching devices, capacitors, and transformers. This includes experimental,
data, topology, and electrical as well as thermal parameter information. The parameters are summarized in one common
dataPara variable.
Inputs:     1) setup:   includes all simulation variables
            2) path:    includes all path variables
Outputs:    1) setup:   extended setup variable
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================
from src.data.loadParaSwi import loadParaSwi
from src.data.loadParaCap import loadParaCap
from src.data.loadParaTra import loadParaTra

# ==============================================================================
# External
# ==============================================================================


#######################################################################################################################
# Function
#######################################################################################################################
def loadPara(setup, path):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("START: Loading Parameter")
    print("------------------------------------------")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    dataPara = {'Swi': {}, 'Cap': {}}

    ###################################################################################################################
    # Loading
    ###################################################################################################################
    # ==============================================================================
    # Switches
    # ==============================================================================
    if setup['Top']['sourceType'] == 'DAB':
        pri_name = setup['Exp'].get('SwiPri', setup['Exp']['Swi'])
        sec_name = setup['Exp'].get('SwiSec', setup['Exp']['Swi'])
        paraSwiPri = loadParaSwi(pri_name, path['parPath'], setup)
        paraSwiSec = loadParaSwi(sec_name, path['parPath'], setup)
        paraSwi = paraSwiPri
    else:
        paraSwi = loadParaSwi(setup['Exp']['Swi'], path['parPath'], setup)
        paraSwiPri = paraSwi
        paraSwiSec = paraSwi

    # ==============================================================================
    # Capacitor
    # ==============================================================================
    if setup['Top']['sourceType'] == 'DAB':
        try:
            pri_name = setup['Exp'].get('CapPri', setup['Exp']['Cap'])
            sec_name = setup['Exp'].get('CapSec', setup['Exp']['Cap'])
            paraCapPri = loadParaCap(pri_name, path['parPath'], setup)
            paraCapSec = loadParaCap(sec_name, path['parPath'], setup)
        except:
            paraCapPri = []
            paraCapSec = []
        paraCap = paraCapPri
    else:
        paraCap = loadParaCap(setup['Exp']['Cap'], path['parPath'], setup)
        paraCapPri = paraCap
        paraCapSec = paraCap

    # ==============================================================================
    # Transformers
    # ==============================================================================
    # DCDC Topologies
    if setup['Top']['sourceType'] == 'DAB':
        paraTra = loadParaTra(setup['Exp']['Trafo'], path['parPath'], setup)
    else:
        paraTra = []

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    dataPara['Swi'] = paraSwi
    dataPara['SwiPri'] = paraSwiPri
    dataPara['SwiSec'] = paraSwiSec
    dataPara['Cap'] = paraCap
    dataPara['CapPri'] = paraCapPri
    dataPara['CapSec'] = paraCapSec
    dataPara['Tra'] = paraTra

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Loading Parameter")
    print("------------------------------------------")

    return dataPara
