# Wind-Energy-Freeze-Analysis
A research project exploring the comparison of freezing events to freezing predicted by reanalysis data, specifically to compare turbine performance between winterized and non-winterized turbines (Texas vs. Iowa).

In 2021, Texas experienced a severe winter storm known as Winter Storm Uri (February 13-17). This storm caused severe power loss to most of the state, along with dangerous conditions for inhabitants of the state. When this storm hit, the wind turbines in Texas were also impacted, contributing to the power loss. Some of this performance loss can be attributed to the icing that occurred on the turbine blades, impeding performance. (Insert various stats and references). Other states, such as Iowa, "winterize" their turbines - add mitigation to prevent the build up of ice on the turbine blades. This helps reduce the negative effects of icing on turbine blades. The purpose of this project is as follows:
1. Determine if ERA5 reanalysis can be used to accurately identify icing events by comparing results to observational data (ASOS).
2. Determine if ERA5 reanalysis can be a useful tool in identifying icing events to supplement the work done using ASOS data.

Goal 1 is important because it shows that ERA5 could be used to identify and understand icing events in other parts of the world. Goal 2 is important because it determine if ERA5 is capable of filling in gaps in the ASOS results due to ASOS acting as singular data points depending on where the respective weather stations are.

## Navigating the GitHub
There are several folders pertaining to different parts of the project. 
- Data: contains raw (ERA5 and ASOS) data. Processing was done entirely within the code and output the results, so there is no processed data.
- Scripts: Contains each of the necessary files to perform the analysis on the raw data
- Results: Contains relevant plots and figures for each of the analyses performed
