##########################
# Abstract data measured #
##########################

abstract type AbstractDataMeasured{T} end

"Gets the grid where the data is measured"
grid(datam::AbstractDataMeasured) = datam.grid

"Gets the region where the data is measured"
roi(datam::AbstractDataMeasured) = datam.roi

"Gets the time where the data is measured"
time_measured(datam::AbstractDataMeasured) = datam.mtime


################################
# measured data implementation #
################################
# include("../Images/data_imgs.jl")
