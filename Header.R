# Name: Header.R
# Auth: u.niazi@imperial.ac.uk
# Date: 26 Aug 2015
# Desc: header file with functions and libraries to load

######################### libraries 
library(annotate)
library(lumiHumanAll.db)
library(lumiHumanIDMapping)
library(lumi)
library(limma)
library(downloader)
source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')
library(NMF)
source('../CGraphClust/CGraphClust.R')
######################### end libraries


######################## global variables
p.old = par()

########################


######################### functions
# Name: f_Plot3DPCA
# Args: mComp = n X 3 matrix with first 3 components as the 3 column vectors
#       color = colours for the points
#       ... additional arguments to the plot function
# Rets: none
# Desc: takes the first 3 components and plots the data in a 3d plot
f_Plot3DPCA = function(mComp, color, ...) {
  x = mComp[,1]
  y = mComp[,2]
  z = mComp[,3]
  if (!require(scatterplot3d)) stop('scatterplot3d library required')
  scatterplot3d(x, y, z, color, ...)
}

######################### end functions