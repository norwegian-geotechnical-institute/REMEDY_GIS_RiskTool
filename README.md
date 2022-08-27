# REMEDY_GIS_RiskTool
An open-source GIS-based tool based on the [GIBV method](https://www.sciencedirect.com/science/article/pii/S0886779820306271) to quantify building damage risk due to deep excavation works from both wall deformation and drawdown-induced consolidation settlements. This tool is a result of the research project [REMEDY / Begrens Skade 2](https://www.sciencedirect.com/science/article/pii/S0886779820306271) (2017-2022) funded by several industry partners and the Norwegian Research Council (NFR). 

## Introduction
One of the main targets for the REMEDY project is to provide research and development that will enable the industry to perform more risk-informed planning, design and execution. Deep excavation and foundation work in soft clays may cause large settlements, resulting in damage to neighbouring buildings and structures. The costs related to these types of damage can be substantial. There is significant potential for reducing such costs if the causes are better understood and the risks are assessed during the project phases. 

**Functionality**
For all cornerpoints on building polygons within the zone of influence surrounding the footprint of a deep excavation, the following is calculated:
1. Analysis of vertical greenfield settlement due to retaining wall deformation, with distance from wall, based on experience database
2. Analysis of groundwater drawdown-induced consolidation settlements in clay, with distance from wall, based on experience database 

## Installation
Copy the source files to a local drive or server drive, keep all the python files, the toolbox (tbx) file and the json file in the same folder. Create a folder for the execution log and the results, and unzip the folder containing the layer files. Make sure the directories specified on the top in the ArcGIS wrappers point to the log and lyr folders. 

### For use with ArcGIS Pro
Download the latest version of ESRI ArcGIS Pro Standard. Create a new project according to the user manual, and connect to the folder where the code is located, using "Add folder connection". Then add the toolbox file "BegrensSkade_GIBV.tbx" to the ArcGIS project by using "Add toolbox". All necessary python packages should be available with the latest ArcGIS installation. Expand the toolbox and run any of the three programs, follow the user manual if necessary.

### For use with QGIS
The python file "BegrensSkade.py" contains three functions for the three programs "Excavation, ImpactMap and Tunnel". These functions are independent of arcpy and can be run in a open-source environment, for example QGIS, without the ArcGIS wrappers. The function arguments are described in the code and in the user manual. Required python libraries must be installed manually by the user. Input feature shapefiles and rasters have to be in the same coordinate system in this option. 

## Feedback and Contribute
Have you identified a problem with the code? Have a feature request? We want to hear about it!
Submit an issue or start a discussion!

## Test
N/A

## Useful notes
1. The NGI project website is located here: https://www.ngi.no/Prosjekter/BegrensSkade-II-REMEDY-Risk-Reduction-of-Groundwork-Damage
2. The journal paper describing the methodology is located here: https://www.sciencedirect.com/science/article/pii/S0886779820306271

## License
Licensed under the MIT license.
