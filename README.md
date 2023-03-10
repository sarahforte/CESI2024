# CESI2024
Scripts to calculate the Canadian Environmental Sustainability Indicators for water quantity

Current R scripts are in two sets:
•	Functions scripts build on the work started by Joe. They are modular functions that are used for more than one indicator.
o	CESI data sorting functions
o	CESI indicator functions
o	CESI trend functions
o	CESI mapping functions

•	Indicator scripts that call on the various functions to calculate indicator and make products used in narrative. Each includes a list of the CESI function files called on.
o	CESI étiage
o	CESI crue
o	CESI rendement
o	CESI statistiques > this uses output from previous three to calculate statistics and make summary figures included in the narrative
 
All CESI indicators use the same colour palette, as described in the document CESI_Colour_Palette.pdf. The functions above use these colours.

Starting in the 2022-2023 work plan, HSU agreed to provide output maps that look like the other CESI indicators, so their typical north arrow and scale bar are added to the maps.
