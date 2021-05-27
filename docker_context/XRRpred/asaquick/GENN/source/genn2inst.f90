!       GENN V.3 (train)
!
!       GEneral Neural Network, V2.2 2013 (c) -h for help
!       Eshel Faraggi, efaraggi@gmail.com
!	You Must Obtain License from above address to use this program
!	
!       g77 -ffixed-line-length-none -Wall -W -Wsurprising  -O2 -o intgrt.e intgrt.f
!       ifort -132 -O2 -o intgrt.e intgrt.f
!	-check bound -g -traceback
!	ifort -openmp -check bound -g -traceback  -o seder.e seder.f90
!	gfortran -g -ffree-line-length-none -Wall -fbounds-check -o seder.e seder.f90 
! 	gfortran -ffree-line-length-none -O2 -o seder.e seder.f90
!       Eshel Faraggi 2011, 2012, 2013 (c)

include "genn2inst.h"

      program genn2inst
      use all_sizes
      use NN_input
      use data_files
      use randseed
      use NN_wei
      use NN_wei_sv
      use NN_owei
      use NN_nodes
      use DB_output
      implicit none
      
      call readcmd
      call init_all
      call seder2train()
      
      stop
      end
!_______________________________________________________________________

