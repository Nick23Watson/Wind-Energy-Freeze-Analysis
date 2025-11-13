# Wind-Energy-Freeze-Analysis
A research project exploring the comparison of freezing events to freezing predicted by reanalysis data, specifically to compare turbine performance between winterized and non-winterized turbines (Texas vs. Iowa).

In 2021, Texas experienced a severe winter storm known as Winter Storm Uri (February 13-17). This storm caused severe power loss to most of the state, along with dangerous conditions for inhabitants of the state. When this storm hit, the wind turbines in Texas were also impacted, contributing to the power loss. Some of this performance loss can be attributed to the icing that occurred on the turbine blades, impeding performance. (Insert various stats and references). Other states, such as Iowa, "winterize" their turbines - add mitigation to prevent the build up of ice on the turbine blades. This helps reduce the negative effects of icing on turbine blades. The purpose of this project is as follows:
1. Determine if ERA5 reanalysis can be used to accurately identify icing events by comparing to observational data
2. Compare icing events in Texas and Iowa, and compare turbine power during these events

Goal 1 is important because it could be used to identify and understand icing events in other parts of the world. Goal 2 is important because it shows the importance of "winterizing" turbines even in areas that icing is not expected. 

## Navigating the GitHub
There are several folders pertaining to different parts of the project. 
- Data: contains both raw (ERA5 and ASOS data) and processed data (working on saving new files with processed data such as total wind speed instead of separate components, etc.)
- Scripts: Contains each of the necessary files to perform the analysis on the raw data
- Results: Contains relevant plots and figures for each of the analyses performed
