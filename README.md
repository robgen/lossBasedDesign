lossBasedDesign
============
[![GitHub Stars](https://img.shields.io/github/stars/robgen/lossBasedDesign.svg)](https://github.com/robgen/lossBasedDesign/stargazers) [![GitHub Issues](https://img.shields.io/github/issues/robgen/lossBasedDesign.svg)](https://github.com/robgen/lossBasedDesign/issues) [![Current Version](https://img.shields.io/badge/version-1.0.0-green.svg)](https://github.com/robgen/lossBasedDesign)

This is a Matlab class that allows performing Direct Loss Based seismic Design (DLBD) for reinforced concrete frames and cantilever walls.

Based on the paper
Gentile R, Calvi GM, 2022. Direct loss-based seismic design of concrete structures. Earthquake Engineering and Structural Dynamics. Under Review

## Buy me a coffee
Whether you use this project, have learned something from it, or just like it, please consider supporting it with a small donation, so I can dedicate more time on open-source projects like this :)

<a href="http://paypal.me/robgen" target="_blank"><img src="https://www.paypalobjects.com/webstatic/mktg/logo/pp_cc_mark_74x46.jpg" alt="Paypal" style="height: auto !important;width: auto !important;" ></a>

## Dependencies
surrogatePSDM (https://github.com/robgen/surrogatedPSDM)
VULNERABILITYbuilding, EALcalculator, CSMnbs, intersections (https://github.com/robgen/robSeismicAnalyses)

## Setup
- From the repo robSeismicAnalyses (https://github.com/robgen/robSeismicAnalyses), copy the functions VULNERABILITYbuilding.m, EALcalculator.m, CSMnbs.m, and intersections.m. Add them to a folder on the Matlab path
- Clone the surrogatePSDM repo (https://github.com/robgen/surrogatedPSDM) to any folder in your computer. Add this folder to your Matlab path
- Clone the lossBasedDesign repo to any folder in your computer. Add this folder to the Matlab path and you are ready to go (or just cd to this folder)

## Usage
A full demo of this class is given in the file exampleLBD.m

## License
This project is licensed under the terms of the **Creative Commons “Attribution-No Derivatives 4.0 International”** license. This software is supplied "AS IS" without any warranties and support. The Author assumes no responsibility or liability for the use of the software. The Author reserves the right to make changes in the software without notification.
