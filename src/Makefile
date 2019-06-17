# Start of the Makefile
# Defining variables
sources = $(wildcard *.f90)
objects = $(sources:%.f90=%.o)
modules = $(wildcard *.mod)
outputs := $(wildcard bLakeJob.sub.*)

ifeq ($(DEBUG), 1)
FFLAGS = -traceback -CB -fpe0 -g
else
FFLAGS = -O3 -traceback -fpe0
endif
FFLAGS += -I$(NETCDF_HOME)/include $(LAPACK_INCLUDE_F95) -I/usr/include
#FFLAGS += -I$(NETCDF_HOME)/include -I/usr/include
LDFLAGS =
LIBS = -L$(NETCDF_HOME)/lib -lpnetcdf $(LINK_LAPACK95) -L/usr/lib64
#LIBS = -L$(NETCDF_HOME)/lib -lpnetcdf -L/usr/lib64

f90comp = mpif90 $(FFLAGS) -c
f90link = mpif90 $(LDFLAGS) -o
# source file search variable VPATH = : : :
# Makefile
ALBM:	$(objects)
		$(f90link) ALBM $(objects) $(LIBS)
bLake.o :	bLake.f90				simulation_mod.o				\
				bayesian_mod.o			read_data_mod.o				\
				sensitivity_mod.o
				$(f90comp) bLake.f90
bayesian_mod.o :		bayesian_mod.f90     	sim_coupler_mod.o		\
							io_utilities_mod.o		read_data_mod.o		\
							shr_ctrl_mod.o				shr_typedef_mod.o		\
							shr_param_mod.o			costfunc_mod.o
							$(f90comp) bayesian_mod.f90
sensitivity_mod.o :	sensitivity_mod.f90		sim_coupler_mod.o		\
							io_utilities_mod.o		read_data_mod.o		\
							shr_ctrl_mod.o				shr_typedef_mod.o		\
							shr_param_mod.o
							$(f90comp) sensitivity_mod.f90
simulation_mod.o :	simulation_mod.f90      sim_coupler_mod.o		\
							shr_param_mod.o			shr_typedef_mod.o		\
							read_data_mod.o
							$(f90comp) simulation_mod.f90
costfunc_mod.o	:		costfunc_mod.f90			math_utilities_mod.o	\
							read_data_mod.o			data_buffer_mod.o
							$(f90comp) costfunc_mod.f90
sim_coupler_mod.o :	sim_coupler_mod.f90		thermal_mod.o			\
							soil_thermal_mod.o		shr_ctrl_mod.o			\
							bubble_mod.o				carbon_cycle_mod.o	\
							diagenesis_mod.o			read_data_mod.o		\
							shr_kind_mod.o				shr_typedef_mod.o		\
							boundary_mod.o				math_utilities_mod.o
							$(f90comp) sim_coupler_mod.f90
diagenesis_mod.o :	diagenesis_mod.f90		data_buffer_mod.o			\
							phy_utilities_mod.o		bg_utilities_mod.o		\
							shr_ctrl_mod.o				shr_param_mod.o
							$(f90comp) diagenesis_mod.f90
carbon_cycle_mod.o :	carbon_cycle_mod.f90		phy_utilities_mod.o		\
							data_buffer_mod.o			bg_utilities_mod.o		\
							shr_ctrl_mod.o				shr_param_mod.o			\
							shr_typedef_mod.o
							$(f90comp) carbon_cycle_mod.f90
bubble_mod.o :		bubble_mod.f90		data_buffer_mod.o				\
						phy_utilities_mod.o	shr_ctrl_mod.o				\
						shr_kind_mod.o			shr_typedef_mod.o
						$(f90comp) bubble_mod.f90
soil_thermal_mod.o :	soil_thermal_mod.f90	phy_utilities_mod.o	\
							data_buffer_mod.o		shr_kind_mod.o			\
							shr_ctrl_mod.o			shr_param_mod.o
							$(f90comp) soil_thermal_mod.f90			
thermal_mod.o :	thermal_mod.f90		phy_utilities_mod.o		\
						bg_utilities_mod.o	data_buffer_mod.o			\
						shr_ctrl_mod.o			shr_param_mod.o			\
						shr_typedef_mod.o
						$(f90comp) thermal_mod.f90
boundary_mod.o :	boundary_mod.f90		data_buffer_mod.o			\
						bg_utilities_mod.o	phy_utilities_mod.o		\
						radiation_mod.o
						$(f90comp) boundary_mod.f90
radiation_mod.o :	radiation_mod.f90		shr_kind_mod.o				\
						phy_const_mod.o		math_utilities_mod.o		\
						radiation_io_mod.o	data_buffer_mod.o			\
						bg_utilities_mod.o
						$(f90comp) radiation_mod.f90
data_buffer_mod.o :	data_buffer_mod.f90	shr_ctrl_mod.o				\
							shr_typedef_mod.o		phy_utilities_mod.o  	\
							read_data_mod.o
							$(f90comp) data_buffer_mod.f90
radiation_io_mod.o :	radiation_io_mod.f90		shr_kind_mod.o			\
							shr_ctrl_mod.o				io_utilities_mod.o	\
							math_utilities_mod.o
							$(f90comp) radiation_io_mod.f90
read_data_mod.o	:	read_data_mod.f90 		shr_ctrl_mod.o 			\
							shr_typedef_mod.o			shr_param_mod.o			\
							io_utilities_mod.o		phy_utilities_mod.o		\
							math_utilities_mod.o		bg_const_mod.o
							$(f90comp) read_data_mod.f90
io_utilities_mod.o :	io_utilities_mod.f90		shr_kind_mod.o			\
							shr_typedef_mod.o			shr_ctrl_mod.o			\
							math_utilities_mod.o
							$(f90comp) io_utilities_mod.f90
bg_utilities_mod.o :		bg_utilities_mod.f90		bg_const_mod.o		\
								phy_const_mod.o			shr_param_mod.o	\
								shr_ctrl_mod.o
								$(f90comp) bg_utilities_mod.f90
phy_utilities_mod.o :   phy_utilities_mod.f90   phy_const_mod.o			\
                        shr_ctrl_mod.o				shr_typedef_mod.o			\
								shr_kind_mod.o
								$(f90comp) phy_utilities_mod.f90
math_utilities_mod.o :	math_utilities_mod.f90	shr_kind_mod.o		\
								shr_ctrl_mod.o
								$(f90comp) math_utilities_mod.f90
bg_const_mod.o :	bg_const_mod.f90		shr_kind_mod.o			\
						shr_ctrl_mod.o
						$(f90comp) bg_const_mod.f90
phy_const_mod.o :	phy_const_mod.f90		shr_kind_mod.o
						$(f90comp) phy_const_mod.f90
shr_param_mod.o :	shr_param_mod.f90		shr_kind_mod.o
						$(f90comp) shr_param_mod.f90
shr_ctrl_mod.o :	shr_ctrl_mod.f90		shr_kind_mod.o			\
						shr_typedef_mod.o
						$(f90comp) shr_ctrl_mod.f90
shr_typedef_mod.o :	shr_typedef_mod.f90	shr_kind_mod.o
							$(f90comp) shr_typedef_mod.f90
shr_kind_mod.o :	shr_kind_mod.f90
						$(f90comp) shr_kind_mod.f90
# Cleaning everything
.PHONY :	clean
clean :		
		-rm -f ALBM $(objects) $(modules)
# End of the makefile
