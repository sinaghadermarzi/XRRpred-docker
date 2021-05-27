!       GENN V.2 -- SEDER MODULE
!       Eshel Faraggi 04/2011 (c)
!!! PLEASE DO NOT USE THIS PROGRAM/CODE YET!

module all_sizes
  integer nte,naax,nrx,nip,nvo,nw1,nw2,nrxo,nfilt,mpout,ngi,ngo,gmmn
  integer ntotp,nisf,nosf,ntrap,ntstp,ntraa, ntsta,icv, jdf, isv, dbv
  integer iug,itbmax,ifilt,ioffl
  double precision fcv,alpha,ratel,rmoment
end module

module NN_input
  double precision, allocatable :: protin(:,:), protut(:,:) ! sequence inputs / outputs
  double precision, allocatable :: globin(:,:), globut(:,:) ! global inputs / outputs
  integer, allocatable :: ibi(:,:), iprotdg(:)	  ! window boundaries
end module

module data_files
  character*3600 datdir, fylein, arg, weifl, oweifl, weiadd, prdfl, prdlist, weiav, foldpart
  character*120 sinfyles(100), sutfyles(100)
  character, allocatable :: protid(:)*3600, aares(:)*12, aacid(:)*12
  integer, allocatable :: lstt(:), lsglo(:)	! map of 10 fold cv order, residue/protein
  integer nsinp,nsutp
end module

module randseed
   integer itsd
end module

module NN_wei
  double precision, allocatable :: wi(:,:,:), wh(:,:), wo(:,:,:)	! Main Sequence weights (ordered)
  double precision, allocatable :: wgi(:,:), wgo(:,:)
end module

module NN_owei					! Old weights (momentum calculation)
  double precision, allocatable :: oldwi(:,:,:), oldwh(:,:), oldwo(:,:,:)
  double precision, allocatable :: oldwgi(:,:), oldwgo(:,:)
end module

module NN_decays
  double precision, allocatable :: dij1(:,:), dij2(:,:), dijh(:,:)
  double precision dfac		! Guiding (decay) parameter
end module
  
module NN_nodes ! Main nodes
  double precision, allocatable :: resv(:,:), hnu1(:), hnu2(:), phi2(:,:), phi2p(:,:)
  double precision, allocatable :: gval(:), gphi2(:), gphi2p(:)
end module

module DB_output ! database indexed predictions
  double precision, allocatable :: phiut(:,:), phiutf(:,:)
  double precision, allocatable :: gphiut(:,:), gphiutf(:,:)
  integer, allocatable :: icphi(:)
end module

      subroutine seder()
      use all_sizes
      use NN_input
      use data_files
      use randseed
      use NN_wei
      use NN_owei
      use NN_decays
      use NN_nodes
      use DB_output
      implicit none
      integer itb,ne,ix,ip,iglo,iq,icond,j,imax,itbopt, nav
      double precision trer, ofer,teer,cormax
      double precision, allocatable :: sphiut(:,:,:)
      character tmpfl2*3600

      if (weiav.ne."") then
	if (prdfl.eq."") then
	  write(0,*) "Must have -pr1/l with -aw, aborting."
	  stop
	endif
        open(unit=29,file=weiav,status='old')
        allocate( sphiut(nvo,naax,2) )
        sphiut = 0.d0
        nav = 0
      endif
      if (prdfl.ne."") then
        write(0,*) "GEneral Neural Network (predict), V2.2 2011 (c) -h for help"
        write(0,*) "Eshel Faraggi, efaraggi@gmail.com"
        write(0,*) "You Must Obtain a License From Above Address to Use GENN"
        write(0,*) "Weights file: ", trim(weifl)
        write(0,'(200A15)') "INS/OUTS: ", (trim(sinfyles(j)),j=1,nsinp), "/", (trim(sutfyles(j)),j=1,nsutp)
      endif
107   continue
      if (weiav.ne."") then
        read(29,*,end=117) weifl
        if (dbv.gt.0) call init_all
        nav = nav + 1
      endif

      write(tmpfl2,'(a,i6.6,a,i11.11,a)') "genn.", ntotp, "_", itsd, "."
      call netinit(nrx,nip,nw1,nw2,nvo,nrxo,ngi,ngo,iug)
      tmpfl2 = trim(tmpfl2)//trim(weiadd)
      if (weifl.eq."") then
        weifl = trim(tmpfl2)//".wei"
      else
	call readwei(weifl,nrx,nip,nw1,nw2,nvo,nrxo,ngi,ngo,iug,alpha)
        if (prdfl.eq."") weifl = trim(tmpfl2)//".wei"
      endif

      ! Predictions
      if (prdfl.eq."") then
        itb = 0
      else
        call initphi(naax,nvo)
        do ix=1,naax
          iglo = lsglo(ix)
          call readwin(ix,iglo,nip,nvo,nisf,nosf,ngi,ngo)
          call calcnet(ix,nrx,nosf,nip,nw1,nw2,nvo,ngi,ngo,alpha,ibi(1,ix))
        enddo
        call normphi(naax,nvo)
        call calcgerr(1,naax,trer)
        write(0,'(A10,100(F11.4))') "MAIN: ",trer
        if (ifilt.eq.1) then
          call savinp(naax,nvo)
          call netinit(nrx,nvo,nfilt,nfilt,nvo,nrxo,ngi,ngo,iug)
          weifl = trim(weifl)//".filter"
          call initphi(naax,nvo)
          call readwei(weifl,nrx,nvo,nfilt,nfilt,nvo,nrxo,ngi,ngo,iug,alpha)
          do ix=1,naax
            iglo = lsglo(ix)
            call readwinf(ix,iglo,nvo,nvo,nisf,nosf,ngi,ngo)
            call calcnet(ix,nrx,nosf,nvo,nfilt,nfilt,nvo,ngi,ngo,alpha,ibi(1,ix))
          enddo
          call normphi(naax,nvo)
          call calcgerr(1,naax,trer)
          write(0,'(A10,100(F11.4))') "Filter: ",trer
        endif

        if (weiav.ne."") then
          do ix=1,naax
            do j=1,nvo
              sphiut(j,ix,1) = sphiut(j,ix,1) + phiut(j,ix)
              sphiut(j,ix,2) = sphiut(j,ix,2) + phiut(j,ix)*phiut(j,ix)
            enddo
          enddo
        endif
      endif
      if (weiav.ne."") goto 107
117   continue
      if (prdfl.ne."") then
        if (weiav.ne."") then
          close(29)
          do ix=1,naax
            do j=1,nvo
              phiut(j,ix) = sphiut(j,ix,1) / nav
              sphiut(j,ix,2) = dsqrt(sphiut(j,ix,2)/nav - phiut(j,ix)*phiut(j,ix))
            enddo
          enddo
        endif

        if (mpout.eq.0) then
          do ix=1,naax
            iglo = lsglo(ix)
            if (weiav.ne."") then
              write(*,'(A,1x,A7,1x,A10,1x,200G12.4)') trim(protid(iglo)), trim(aacid(ix)), trim(aares(ix)), (phiut(j,ix),j=1,nvo),(sphiut(j,ix,2),j=1,nvo)
            else
              write(*,'(A,1x,A7,1x,A10,1x,200G12.4)') trim(protid(iglo)), trim(aacid(ix)), trim(aares(ix)), (phiut(j,ix),j=1,nvo)
            endif
          enddo
        elseif (mpout.eq.2) then
          do ix=1,naax
            iglo = lsglo(ix)
            trer = nvo ; teer = nvo
            do j=1,nvo
              trer = trer + protut(j,ix)
              teer = teer + phiut(j,ix)
            enddo
            if (trer.le.0.d0) trer = 1.d0
            if (teer.le.0.d0) teer = 1.d0
            do j=1,nvo
              protut(j,ix) = (1.d0+protut(j,ix)) / trer
              phiut(j,ix) = (1.d0+phiut(j,ix)) / teer
            enddo
            write(*,'(A,1x,A7,1x,A10,1x,I5,I5,200G12.4)') trim(protid(iglo)), trim(aacid(ix)), trim(aares(ix)), imax(nvo,protut(1,ix)), imax(nvo,phiut(1,ix)), (phiut(j,ix),j=1,nvo), (protut(j,ix),j=1,nvo)
          enddo
        else          ! -po 1 for comparison print prediction, natives
          do ix=1,naax
            iglo = lsglo(ix)
            write(*,'(A,1x,A7,1x,A10,1x,200G12.4)') trim(protid(iglo)), trim(aacid(ix)), trim(aares(ix)), (phiut(j,ix),j=1,nvo) , (protut(j,ix),j=1,nvo)
          enddo
        endif
        stop
      endif

      if (oweifl.ne."") weifl = oweifl
      
      if (mpout.lt.1) then
        cormax = 20.d0
      else
        cormax = -20.d0
      endif

      tmpfl2 = trim(weifl)//".fort.15"
      open(unit=ioffl,file=tmpfl2,status='new')
      write(ioffl,*) "GEneral Neural Network, V2.2 2011 (c) -h for help"
      write(ioffl,*) "Eshel Faraggi, efaraggi@gmail.com"
      write(ioffl,*) "You Must Obtain a License From Above Address to Use GENN"
      write(ioffl,*) "Options: ", trim(weiadd)
      write(ioffl,*) "Weights file: ", trim(weifl)
      write(ioffl,'(a,i15)') "Random Seed: ", itsd
      write(ioffl,'(200A15)') "INS/OUTS: ", (trim(sinfyles(j)),j=1,nsinp), "/", (trim(sutfyles(j)),j=1,nsutp)
      write(ioffl,*) trim(foldpart)

      ne = 1
!!!!! training epochs (main network):
      do while ((ne.le.nte).and.(itb.lt.itbmax))
        call initphi(naax,nvo)
        do ix=1,ntraa
          ip = lstt(ix)
          iglo = lsglo(ip)
          do iq=1,iprotdg(ip)
            call readwin(ip,iglo,nip,nvo,nisf,nosf,ngi,ngo)
            call calcnet(ip,nrx,nosf,nip,nw1,nw2,nvo,ngi,ngo,alpha,ibi(1,ip))
            call backprop(nrx,nip,nw1,nw2,nrxo,nvo,ngi,ngo,alpha,ratel,rmoment,ibi(1,ip))
          enddo
        enddo
        do ix=ntraa+1,naax
          ip = lstt(ix)
          iglo = lsglo(ip)
          call readwin(ip,iglo,nip,nvo,nisf,nosf,ngi,ngo)
          call calcnet(ip,nrx,nosf,nip,nw1,nw2,nvo,ngi,ngo,alpha,ibi(1,ip))
        enddo
        call normphi(naax,nvo)

        call calcierr(trer,ofer,teer,icond,cormax)
        if ((icond.eq.1).or.(ne.le.1)) then
          call outwei(weifl,nrx,nip,nw1,nw2,nvo,nrxo,ngi,ngo,iug,alpha)
          cormax = ofer
          itbopt = ne
          itb = 0
          write(ioffl,'(I10,100(F11.4))') ne,trer,ofer,teer
        else
          itb = itb + 1
        endif
        ne = ne + 1
      enddo

      call readwei(weifl,nrx,nip,nw1,nw2,nvo,nrxo,ngi,ngo,iug,alpha)
      call initphi(naax,nvo)
      do ix=1,naax
        iglo = lsglo(ix)
        call readwin(ix,iglo,nip,nvo,nisf,nosf,ngi,ngo)
        call calcnet(ix,nrx,nosf,nip,nw1,nw2,nvo,ngi,ngo,alpha,ibi(1,ix))
      enddo
      call normphi(naax,nvo)
      call calcierr(trer,ofer,teer,icond,cormax)
      write(ioffl,'(A10,I10,1x,I10,1x,3(F11.4))') "MAIN: ",ne-1,itbopt,trer,ofer,teer

      if (ifilt.eq.1) then
        write(ioffl,*)
        call savinp(naax,nvo)
        write(ioffl,*) "Filter:"
        weifl = trim(weifl)//".filter"
        call netinit(nrx,nvo,nfilt,nfilt,nvo,nrxo,ngi,ngo,iug)
        itb = 0
        ne = 1
        do while ((ne.le.nte).and.(itb.lt.itbmax))
          call initphi(naax,nvo)
          do ix=1,ntraa
            ip = lstt(ix)
            iglo = lsglo(ip)
            do iq=1,iprotdg(ip)
              call readwinf(ip,iglo,nvo,nvo,nisf,nosf,ngi,ngo)
              call calcnet(ip,nrx,nosf,nvo,nfilt,nfilt,nvo,ngi,ngo,alpha,ibi(1,ip))
              call backprop(nrx,nvo,nfilt,nfilt,nrxo,nvo,ngi,ngo,alpha,ratel,rmoment,ibi(1,ip))
            enddo
          enddo
          do ix=ntraa+1,naax
            ip = lstt(ix)
            iglo = lsglo(ip)
            call readwinf(ip,iglo,nvo,nvo,nisf,nosf,ngi,ngo)
            call calcnet(ip,nrx,nosf,nvo,nfilt,nfilt,nvo,ngi,ngo,alpha,ibi(1,ip))
          enddo
          call normphi(naax,nvo)

          call calcierr(trer,ofer,teer,icond,cormax)
          if ((icond.eq.1).or.(ne.le.1)) then
            call outwei(weifl,nrx,nvo,nfilt,nfilt,nvo,nrxo,ngi,ngo,iug,alpha)
            cormax = ofer
            itbopt = ne
            itb = 0
            write(ioffl,'(I10,100(F11.4))') ne,trer,ofer,teer
          else
            itb = itb + 1
          endif
        ne = ne + 1
        enddo

        call readwei(weifl,nrx,nvo,nfilt,nfilt,nvo,nrxo,ngi,ngo,iug,alpha)
        call initphi(naax,nvo)
        do ix=1,naax
          iglo = lsglo(ix)
          call readwinf(ix,iglo,nvo,nvo,nisf,nosf,ngi,ngo)
          call calcnet(ix,nrx,nosf,nvo,nfilt,nfilt,nvo,ngi,ngo,alpha,ibi(1,ix))
        enddo
        call normphi(naax,nvo)
        call calcierr(trer,ofer,teer,icond,cormax)
        write(ioffl,'(A10,I10,1x,I10,1x,100(F11.4))') "Filter: ",ne-1,itbopt,trer,ofer,teer
      endif
      close(ioffl)
      
      end
!_______________________________________________________________________
!--
      subroutine readcmd
      use all_sizes
      use NN_input
      use data_files
      use randseed
      use NN_wei
      use NN_owei
      use NN_decays
      use NN_nodes
      use DB_output
      implicit none
      integer inttget,nargs,icrd,K,j,i1,i2,i3
      integer, allocatable :: seeds(:)
      character tmpfl*3600

!       Defaults
      fylein = "genn.in"	 !name of input file
      datdir = ""		 !directory where input files are
      weifl = ""                 !weight file
      oweifl = ""                !weights output file
      prdfl = ""                 !id for input files for prediction
      prdlist = ""               !file of ids for prediction
      nte = 1000		 !maximum number of epochs
      ntotp = 0			 !total number of input files to use (if not given all files are used)
      itsd = inttget(7)	 	 !seed for random number generator (0 - constant seed (debugging))
      icv = 10			 !set to cross validate on
      fcv = 0.1d0		 !fraction of dataset to cross validate on (fraction of test fold from dataset)
      nrx = 21			 !size of input window
      nrxo = 1			 !size of output window (not implemented yet)
      nw1 = 21			 !nodes in hidden layer 1 (hl1)
      nw2 = 21			 !nodes in hidden layer 2 (hl2)
      nfilt = 11		 !nodes in filter hidden layer
      iug = 1			 !use guiding?
      dfac = 2.d0		 !guiding factor
      alpha = 0.2		 !activation parameter
      itbmax = 100		 !maximum number of non-improving training epochs
      ratel = 0.001		 !learning rate
      rmoment = 0.4		 !momentum
      ifilt = 0			 !run filter network? (0-no/1-yeah)
      mpout = 0			 !probability outputs?
      ngi = 0			 !number of global inputs
      ngo = 0			 !number of global outputs
      jdf = 0			 !use degeneracy file
      ioffl = 15		 !file number for first network weights
      isv = 0 			 !no server version
      gmmn = 0			 !normalize with default z-score method (otherwise normalize with max/min)?
      weiadd = ""
      weiav = ""
      dbv = 0 			 !different db between different weights (prediction)

      nargs=iargc()
      if (nargs.ge.1) then
        icrd = 1
!       Switches:
        do while (icrd.le.nargs)
          call getarg(icrd,arg)
          if (arg.eq."-l") then			 !db_list_file
            call getarg(icrd+1,tmpfl)
            fylein = trim(tmpfl)
	    i1 = 1 + SCAN(fylein, "/", .TRUE.)
	    tmpfl = fylein(i1:)
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-d") then		 !data_dir
            call getarg(icrd+1,tmpfl)
            datdir = trim(tmpfl)
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-wf") then		 !weight file
            call getarg(icrd+1,tmpfl)
            weifl = trim(tmpfl)
            if (tmpfl(25:25).eq.".") then
              weiadd = trim(weiadd)//trim(arg)//tmpfl(13:24)
            else
              weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
            endif
          elseif (arg.eq."-owf") then		 !weight file
            call getarg(icrd+1,tmpfl)
            oweifl = trim(tmpfl)
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-pr1") then		 !predict one file
            call getarg(icrd+1,tmpfl)
            prdfl = trim(tmpfl)
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-prl") then		 !predict one file
            call getarg(icrd+1,tmpfl)
            prdfl = trim(tmpfl)
            prdlist = trim(tmpfl)
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-m") then		 !max_num_epochs
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) nte
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-np") then		 !number_db_files_to_use
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) ntotp
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-r") then		 !random_seed
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) itsd
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-cv") then		 !test_fold_index
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) icv
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-f") then		 !fraction_of_test_fold_of_db
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) fcv
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-wi") then		 !input_window
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) nrx
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-wo") then		 !output_window
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) nrxo
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-h1") then		 !num_hid_layer1_nodes
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) nw1
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-h2") then		 !num_hid_layer2_nodes
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) nw2
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-hf") then		 !num_filter_hidden_layer_nodes
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) nfilt
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-nf") then		 !run_filter_network
            ifilt = 1
            weiadd = trim(weiadd)//trim(arg)
          elseif (arg.eq."-mi") then		 !max_num_bad_iterations
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) itbmax
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-gf") then		 !guiding_factor
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) dfac
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-ng") then		 !dont_use_distance_guiding
            iug = 0
            weiadd = trim(weiadd)//trim(arg)
          elseif (arg.eq."-a") then		 !activation_parameter
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) alpha
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-u") then		 !learning_rate
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) ratel
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-p") then		 !momentum
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) rmoment
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-gi") then		 !num_global_inputs
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) ngi
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-go") then		 !num_global_outputs
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) ngo
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-po") then		 !probability_output
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) mpout
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-gm") then		 !normalize_globin_w_minmax
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) gmmn
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-df") then		 !Use_degeneracy_files
            jdf = 1
            weiadd = trim(weiadd)//trim(arg)
          elseif (arg.eq."-sv") then		 !Generate_server(only_overfitset)
            isv = 1
            weiadd = trim(weiadd)//trim(arg)
          elseif (arg.eq."-aw") then		 !average over weights in file
            call getarg(icrd+1,tmpfl)
            weiav = trim(tmpfl)
            open(unit=33,file=weiav,status='old')
            read(33,*) weifl
            close(33)
            fylein = ""
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
!            write(0,'(1x,130(A),1x)',advance='no') trim(arg), " ", trim(tmpfl)
          elseif (arg.eq."-dw") then		 !different db b/w different weight files
            dbv = 1
            weiadd = trim(weiadd)//trim(arg)
          elseif (arg.eq."-h") then		 !print_this_help
            call printhelp()
            stop
          elseif (arg(1:1).eq."-") then		 !No negatives?
            write(0,*) "Aborting, didnt recognize ", trim(arg)
            stop
!          else
!            write(0,*) "(aborting) Unrecognized param.: ", arg
!            stop
          endif
          icrd = icrd + 1
        enddo
      endif
      call random_seed(size=K)
      allocate(seeds(K))
      seeds = (/ (1+itsd/j, j=1, K) /)
      call random_seed(put=seeds)
      deallocate(seeds)

      end
!_______________________________________________________________________
!--
      subroutine init_all
      use all_sizes
      use NN_input
      use data_files
      use randseed
      use NN_wei
      use NN_owei
      use NN_decays
      use NN_nodes
      use DB_output
      implicit none
      integer nipmax,nw1max,nw2max

      if ((fylein.ne."").or.(weifl.ne."")) then
        call initdb()
        nipmax = max(nip,nvo)
        nw1max = max(nw1,nfilt)
        nw2max = max(nw2,nfilt)
        if (allocated(wi)) then
          deallocate( wi , wh , wo )
          deallocate( oldwi , oldwh , oldwo )
          deallocate( dij1 , dij2 , dijh )
          deallocate( resv , hnu1 , hnu2 , phi2 , phi2p )
          deallocate( phiut, phiutf, icphi )
          if (ngi.gt.0) deallocate( wgi , oldwgi , gval )
          if (ngo.gt.0) deallocate( wgo , oldwgo , gphi2 , gphi2p , gphiut )
        endif
        allocate( wi(nrx+1,nipmax,nw1max) , wh(nw1max+1,nw2max) , wo(nw2max+1,nrxo,nvo) )
        allocate( oldwi(nrx+1,nipmax,nw1max) , oldwh(nw1max+1,nw2max) , oldwo(nw2max+1,nrxo,nvo) )
        allocate( dij1(nrx,nw1max) , dij2(nw2max,nrxo) , dijh(nw1max,nw2max) )
        allocate( resv(nrx,nipmax) , hnu1(nw1max) , hnu2(nw2max) , phi2(nrxo,nvo) , phi2p(nrxo,nvo) )
        allocate( phiut(nvo,naax), phiutf(nvo,naax), icphi(naax) )
        if (ngi.gt.0) allocate( wgi(ngi,nw1max) , oldwgi(ngi,nw1max) , gval(ngi) )
        if (ngo.gt.0) allocate( wgo(nw2max+1,ngo) , oldwgo(nw2max+1,ngo) , gphi2(ngo) , &
                              gphi2p(ngo) , gphiut(ngo,naax) )
      endif
      end
!_______________________________________________________________________
!--
      subroutine printhelp()
      use all_sizes
      use data_files
      implicit none
      
      write(*,*) "Options:"

          write(*,'(A10,A40,A20)') "Switch","Description                      ", "Default"
          write(*,'(A10,A40,A20)') "-l ","db_list_file                     ", trim(fylein)
          write(*,'(A10,A40,A20)') "-d ","data_dir                         ", trim(datdir)
          write(*,'(A10,A40,A20)') "-wf","weight input file                ", trim(weifl)
          write(*,'(A10,A40,A20)') "-owf","weights output file             ", trim(oweifl)
          write(*,'(A10,A40,A20)') "-pr1","predict single file (give id)   ", trim(prdfl)
          write(*,'(A10,A40,A20)') "-prl","file of ids to predict          ", trim(prdlist)
          write(*,'(A10,A40,I20)') "-m ","max_num_epochs                   ", nte
          write(*,'(A10,A40,A20)') "-np","number_db_files_to_use           ","whole_db"
          write(*,'(A10,A40,A20)') "-r ","random_seed                      ", "time_dep"
          write(*,'(A10,A40,I20)') "-cv","test_fold_index                  ",icv
          write(*,'(A10,A40,D20.6)') "-f ","fraction_of_test_fold_of_db      ", fcv
          write(*,'(A10,A40,I20)') "-wi","input_window                     ",nrx
          write(*,'(A10,A40,I20)') "-wo","output_window                    ",nrxo
          write(*,'(A10,A40,I20)') "-h1","num_hid_layer1_nodes             ", nw1
          write(*,'(A10,A40,I20)') "-h2","num_hid_layer2_nodes             ", nw2
          write(*,'(A10,A40,I20)') "-mi","max_num_bad_iterations           ", itbmax
          write(*,'(A10,A40,D20.6)') "-a ","activation_parameter             ", alpha
          write(*,'(A10,A40,D20.6)') "-u ","learning_rate                    ", ratel
          write(*,'(A10,A40,D20.6)') "-p ","momentum                         ", rmoment

          write(*,'(A10,A40,A20)') "-hf","num_filter_hidden_layer_nodes    ","11"
          write(*,'(A10,A40,A20)') "-nf","run_network_filter               "
          write(*,'(A10,A40,A20)') "-mi","max_num_bad_iterations           ","100"      
          write(*,'(A10,A40,A20)') "-gf","guiding_factor                   ","2.d0"
          write(*,'(A10,A40,A20)') "-ng","dont_use_distance_guiding        "
          write(*,'(A10,A40,A20)') "-gi","num_global_inputs                ","0"
          write(*,'(A10,A40,A20)') "-go","num_global_outputs               ","0"
          write(*,'(A10,A40,A20)') "-so","average_over_output_(nvo)        "
          write(*,'(A10,A40,I20)') "-po","print output 0-reg,1-compr,2-prob",mpout
          write(*,'(A10,A40,A20)') "-df","Use_degeneracy_files             "
          write(*,'(A10,A40,A20)') "-sv","Generate_server(only_overfitset) "
          write(*,'(A10,A40,A20)') "-aw","average over weights in file     "
          write(*,'(A10,A40,A20)') "-dw","diff nn struc for diff wei files "
          write(*,'(A10,A40,A20)') "-gm","glob.norm.? 0-no,1-maxmin,2-zscor",gmmn
          write(*,'(A10,A40,A20)') "-h ","print_this_help                  "
      write(*,*) 
      write(*,*) "Train NN prediction. In genn.in first line is database directory,"
      write(*,*) "second line contains input types (files in ids dirs bellow) third line"
      write(*,*) "contains output types. Rest of lines are ids for"
      write(*,*) "subdirs in database directory (ids dirs). "
      write(*,*) "Global input/output files should be stored in genn.gin/genn.gut in ids dirs"
      write(*,*) "-np controls how many files out of db_list_file"
      write(*,*) "to use. Filter is useful for probability output. Use global"
      write(*,*) "inputs/outputs only in addition to sequence based values. For"
      write(*,*) "only global inputs use genn2inst. Degeneracy files should be"
      write(*,*) "stored in id-dirs under file name genn.dgn"

      end
!_______________________________________________________________________
!--
      subroutine initdb()
      use all_sizes
      use NN_input
      use data_files
      implicit none
      character arbig*3600, cchk(50000)*12, tempc*12, tempc2*12
      integer i,j,k,nfold,numval,ife,ifb,icountel,iof,nf2,j2s,j2r,ntao, numlns, jj, jj2, kk, j3, numres, inta, ib
      integer ioerr, nips
      double precision atmp(20000)
      character tmpfl*3600
      
      if (weifl.ne."") then
        open(unit=22,file=weifl,status='old')
        read(22,*) sinfyles(1)
        if (datdir.eq."") datdir = sinfyles(1)
        read(22,'(A)') arbig
        nsinp = icountel(arbig)
        backspace(22)
        read(22,*) (sinfyles(j),j=1,nsinp)
        read(22,'(A)') arbig
        nsutp = icountel(arbig)
        backspace(22)
        read(22,*) (sutfyles(j),j=1,nsutp)
        read(22,*) nrx,nip,nw1,nw2,nvo,nrxo,ngi,ngo,alpha
        if (prdfl.ne."") then
          if (prdlist.eq."") then
            ntotp = 1
            arg = prdfl
          else
            call numpro(i,prdlist)                !get number of database files
            if ((ntotp.le.0).or.(ntotp.gt.i+3)) then
              ntotp = i + 3
            endif
            close(22)
            open(unit=22,file=prdlist,status='old')
          endif
          j = 0
          do i=1,ntotp
            if (prdlist.ne."") then
              read(22,'(A)')arg
            else
              arg = prdfl
            endif
            tmpfl = trim(datdir)//"/"//trim(arg)//"/"//trim(sinfyles(1))
            open(unit=33,file=tmpfl,status='old')
            do while (j.ge.0)
              read(33,*,end=308)
              j = j + 1
            enddo
308         continue
            close(33)
          enddo
          naax = j
          if (allocated(lstt)) then
            deallocate( lstt , lsglo , protin , iprotdg , protut , ibi , protid,  aares, aacid )
	    if (allocated(globin)) deallocate(globin)
	    if (allocated(globut)) deallocate(globut)
          endif
          allocate( lstt(naax) , lsglo(naax) , protin(nip+1,naax+1) , iprotdg(naax) , protut(nvo,naax) , ibi(4,naax) , protid(ntotp),  aares(naax), aacid(naax) )
          if (ngi.gt.0) allocate( globin(ngi,ntotp) )
          if (ngo.gt.0) allocate( globut(ngo,ntotp) )
          nisf = (nrx+1)/2
          nosf = (nrxo+1)/2
          rewind(22)
          j = 1
          do i=1,ntotp
            ntao = j - 1
            nips = 0
            if (prdlist.ne."") then
              read(22,'(A)')arg
            else
              arg = prdfl
            endif
            protid(i) = trim(arg)
            do k=1,nsinp
              tmpfl = trim(datdir)//"/"//trim(arg)//"/"//trim(sinfyles(k))
              open(unit=33,file=tmpfl,status='old')
              read(33,'(A)') arbig
              backspace(33)
              numlns = icountel(arbig) - 2
              jj = 0
              do while (jj.ge.0)
                jj = jj + 1
                read(33,*,iostat=ioerr,end=317) tempc,tempc2,(protin(nips+kk,ntao+jj),kk=1,numlns)
                if (ioerr.ne.0) then
                  write(0,*) "Reading problem (assuming 0): ", jj, trim(tmpfl)
                  do kk=1,numlns
                    protin(nips+kk,ntao+jj) = 0.d0
                  enddo
                endif
                aares(ntao+jj) = tempc2
                aacid(ntao+jj) = tempc
                if (k.eq.1) then
                  cchk(jj) = tempc
                else
                  if (tempc.ne.cchk(jj)) then
                    write(0,*) trim(tmpfl), " MIS-ALIGN(in), aborting: ", tempc, jj, k, cchk(jj)
                    stop
                  endif
                endif
              enddo
317           continue
              close(33)
              jj2 = jj - 1
              nips = nips + numlns    !number of inputs (per residue)
            enddo
            if (nips.ne.nip) then
              write(0,*) "nip.ne.nips, aborting.", nip, nips
              stop
            endif
            nips = 0
            do k=1,nsutp
              tmpfl = trim(datdir)//"/"//trim(arg)//"/"//trim(sutfyles(k))
              open(unit=33,file=tmpfl,status='old')
              read(33,'(A)') arbig
              backspace(33)
              numlns = icountel(arbig) - 2
              do jj=1,jj2
                read(33,*,end=319) tempc,tempc2,(protut(nips+kk,ntao+jj),kk=1,numlns)
                if (tempc.ne.cchk(jj)) then
                  write(0,*) trim(tmpfl), " MIS-ALIGN(out), aborting: ", tempc, jj, k, cchk(jj)
                  stop
                endif
              enddo
319           continue
              close(33)
              nips = nips + numlns    !number of inputs (per residue)
            enddo
            if (nips.ne.nvo) then
              write(0,*) "nvo.ne.nips, aborting.", nvo, nips
              stop
            endif
            if (ngi.gt.0) then
              tmpfl = trim(datdir)//"/"//trim(arg)//"/genn.gin"
              open(unit=33,file=tmpfl,status='old')
              read(33,*) (globin(j3,i),j3=1,ngi)
              close(33)
            endif
            if (ngo.gt.0) then
              tmpfl = trim(datdir)//"/"//trim(arg)//"/genn.gut"
              open(unit=34,file=tmpfl,status='old')
              read(34,*) (globut(j3,i),j3=1,ngo)
              close(34)
            endif
            numres = 0
            do jj=1,jj2
              lsglo(j) = i
              j = j + 1
              numres = numres + 1
            enddo
!* Setetetup the edges for each amino acid
            do k=1,numres
              inta = ntao + k
              ib = 1 + nisf - k
              if(ib.lt.1) ib = 1
              ibi(1,inta) = ib
              ib = numres + nisf - k
              if(ib.gt.nrx) ib = nrx
              ibi(2,inta) = ib
              ib = 1 + nosf - k
              if(ib.lt.1) ib = 1
              ibi(3,inta) = ib
              ib = numres + nosf - k
              if(ib.gt.nrxo) ib = nrxo
              ibi(4,inta) = ib
            enddo
          enddo
          close(22)
        endif
      endif

      if (prdfl.ne."") then
        if (weifl.eq."") then
          write(0,*) "Can't predict w/o -wf, aborting."
          stop
        endif
      else
        call numpro(i,fylein)          !get number of database files
        if ((ntotp.le.0).or.(ntotp.gt.i)) then
          ntotp = i
        endif
        nfold = nint(fcv*ntotp)
        ntstp = ntotp-nfold
        numval = nint(0.05d0 * ntotp)
        ntrap = ntstp - numval
        ife = nint(fcv * icv * ntotp)
        ifb = ife - nfold

        ntraa = 0
        open(unit=22,file=fylein,status='old')
        read(22,*) sinfyles(1)
        if (datdir.eq."") datdir = sinfyles(1)
        read(22,'(A)') arbig
        nsinp = icountel(arbig)
        backspace(22)
        read(22,*) (sinfyles(j),j=1,nsinp)
        read(22,'(A)') arbig
        nsutp = icountel(arbig)
        backspace(22)
        read(22,*) (sutfyles(j),j=1,nsutp)

        naax = 0
        j = 0; ntsta = 0; iof = 0
        do i=1,ntotp
          if ((i.le.ifb).or.(i.gt.ife)) iof = iof + 1
          if ((iof.eq.(ntrap+1)).and.(ntraa.eq.0)) ntraa = ntsta
          read(22,'(A)',end=346)arg
          k = 1
          tmpfl = trim(datdir)//"/"//trim(arg)//"/"//sutfyles(k)
          open(unit=33,file=tmpfl,status='old')
          do 335 while (j.ge.0)
            read(33,*,end=336)
            if ((i.le.ifb).or.(i.gt.ife)) ntsta = ntsta + 1 
            j = j + 1
335       continue
336       continue
          close(33)
        enddo
346     continue
        close(33)
        naax = j                  !number of residues
        
        nip = 0
        do k=1,nsinp
          tmpfl = trim(datdir)//"/"//trim(arg)//"/"//sinfyles(k)
          open(unit=33,file=tmpfl,status='old')
          read(33,'(A)') arbig
          close(33)
          nip = nip + icountel(arbig) - 2 !number of inputs (per residue)
        enddo

        nvo = 0
        do k=1,nsutp
          tmpfl = trim(datdir)//"/"//trim(arg)//"/"//sutfyles(k)
          open(unit=33,file=tmpfl,status='old')
          read(33,'(A)') arbig
          close(33)
          nvo = nvo + icountel(arbig) - 2 !number of output targets (per residue)
        enddo
        
        if ((naax.lt.20).or.(nip.lt.1).or.(nvo.lt.1).or.(ntotp.lt.1)) then
          write(0,*) "Something wrong with database, aborting."
          write(0,*) naax, nip, nvo, ntotp
          stop
        endif

        allocate( lstt(naax) , lsglo(naax) , protin(nip+1,naax+1) , iprotdg(naax) , protut(nvo,naax) , ibi(4,naax) , protid(ntotp),  aares(naax), aacid(naax) )
        if (ngi.gt.0) allocate( globin(ngi,ntotp) )
        if (ngo.gt.0) allocate( globut(ngo,ntotp) )
        nisf = (nrx+1)/2
        nosf = (nrxo+1)/2
        
        rewind(22)
        read(22,*); read(22,*); read(22,*)

        if ((numval.lt.1).or.(nfold.lt.1).or.(ntsta.eq.naax).or.(isv.eq.1)) then
          nf2 = nint(fcv*naax)
          ntsta = naax - nf2
          numval = nint(0.05d0 * naax)
          ntraa = ntsta - numval
          ife = nint(fcv * icv * naax)
          ifb = ife - nf2
          if (isv.eq.1) then
            ntrap = nint(0.95d0 * ntotp)
            ntstp = ntotp
            ntraa = nint(0.95d0 * naax)
            ntsta = naax
          endif
          write(foldpart,'(A24,9(I9),F8.3,2(I9))') "Residue fold partition: ", nip, nvo, ntrap, ntstp, ntotp, ntraa, ntsta, naax, icv,fcv,ngi,ngo
        else
          write(foldpart,'(A24,9(I9),F8.3,2(I9))') "Protein fold partition: ", nip, nvo, ntrap, ntstp, ntotp, ntraa, ntsta, naax, icv,fcv,ngi,ngo
        endif
        j2s = ntsta + 1
        j2r = 1
        j = 1
        do i=1,ntotp
          ntao = j - 1
          nip = 0
          nvo = 0
          if ((prdfl.eq."").or.(prdlist.ne."")) then
            read(22,'(A)')arg
          endif
          protid(i) = trim(arg)
          do k=1,nsinp
            tmpfl = trim(datdir)//"/"//trim(arg)//"/"//trim(sinfyles(k))
            open(unit=33,file=tmpfl,status='old')
            read(33,'(A)') arbig
            backspace(33)
            numlns = icountel(arbig) - 2
            jj = 0
            do while (j.gt.0)
              jj = jj + 1
              read(33,*,iostat=ioerr,end=297) tempc,tempc2,(protin(nip+kk,ntao+jj),kk=1,numlns)
              if (ioerr.ne.0) then
                write(0,*) "Reading problem(T) (assuming 0): ", jj, trim(tmpfl)
                do kk=1,numlns
                  protin(nip+kk,ntao+jj) = 0.d0
                enddo
              endif
              aares(ntao+jj) = tempc2
              aacid(ntao+jj) = tempc
              if (k.eq.1) then
                cchk(jj) = tempc
              else
                if (tempc.ne.cchk(jj)) then
                  write(0,*) trim(tmpfl), " MIS-ALIGN(in), aborting: ", tempc, jj, k, cchk(jj)
                  stop
                endif
              endif
            enddo
297         continue
            jj2 = jj - 1
            close(33)
            nip = nip + numlns    !number of inputs (per residue)
          enddo
          
          do k=1,nsutp
            tmpfl = trim(datdir)//"/"//trim(arg)//"/"//trim(sutfyles(k))
            open(unit=33,file=tmpfl,status='old')
            read(33,'(A)') arbig
            backspace(33)
            numlns = icountel(arbig) - 2
            do jj=1,jj2
              read(33,*,end=301) tempc,tempc2,(protut(nvo+kk,ntao+jj),kk=1,numlns)
              if (tempc.ne.cchk(jj)) then
                write(0,*) trim(tmpfl), " MIS-ALIGN(out), aborting: ", tempc, jj, k, cchk(jj)
                stop
              endif
            enddo
            close(33)
301         continue
            nvo = nvo + numlns    !number of inputs (per residue)
          enddo

          if (ngi.gt.0) then
            tmpfl = trim(datdir)//"/"//trim(arg)//"/genn.gin"
            open(unit=33,file=tmpfl,status='old')
            read(33,*) (globin(j3,i),j3=1,ngi)
            close(33)
          endif
          if (ngo.gt.0) then
            tmpfl = trim(datdir)//"/"//trim(arg)//"/genn.gut"
            open(unit=34,file=tmpfl,status='old')
            read(34,*) (globut(j3,i),j3=1,ngo)
            close(34)
          endif
          if (jdf.eq.1) then
            tmpfl = trim(datdir)//"/"//trim(arg)//"/genn.dgn"
            open(unit=35,file=tmpfl,status='old')
          endif
          numres = 0
          do jj=1,jj2
            if (jdf.eq.1) then
              read(35,*,end=376) iprotdg(j)
            else
              iprotdg(j) = 1
            endif
            lsglo(j) = i
            if ((ifb.lt.i).and.(i.le.ife)) then
              lstt(j2s) = j
              j2s = j2s + 1
            else
              lstt(j2r) = j
              j2r = j2r + 1
            endif
            j = j + 1
            numres = numres + 1
          enddo
376       continue
          close(35)

!* Setetup the edges for each amino acid
          do k=1,numres
            inta = ntao + k
            ib = 1 + nisf - k
            if(ib.lt.1) ib = 1
            ibi(1,inta) = ib
            ib = numres + nisf - k
            if(ib.gt.nrx) ib = nrx
            ibi(2,inta) = ib
            ib = 1 + nosf - k
            if(ib.lt.1) ib = 1
            ibi(3,inta) = ib
            ib = numres + nosf - k
            if(ib.gt.nrxo) ib = nrxo
            ibi(4,inta) = ib
          enddo

        enddo
        close(22)
      endif

      do j=1,ngi
        do i=1,ntotp
          atmp(i) = globin(j,i)
        enddo
        if (gmmn.eq.1) then
          call maxminnorm(ntotp,atmp(1),ioffl)
        elseif (gmmn.eq.2) then
          call stdnorm(ntotp,atmp(1),ioffl)
        endif
        do i=1,ntotp
          globin(j,i) = atmp(i)
        enddo
      enddo

      end
!_______________________________________________________________________
!--
      integer function icountel(arbig)
      implicit none
      character arbig*3600, curc
      integer lenarg, i, ib, ic
      lenarg = len_trim(arbig)
      icountel = 0
      curc = " "
      do i=1,lenarg
        ib = ichar(arbig(i:i))
        ic = ichar(curc)
        if ((ib.ne.9).and.(ib.ne.32)) then
          if ((ic.eq.9).or.(ic.eq.32)) then
            icountel = icountel + 1
          endif
        endif
        curc = arbig(i:i)
      enddo
      end
!_______________________________________________________________________
!--
      subroutine numpro(itmp,fylein)
      implicit none
      integer itmp, i
      character*3600 fylein
      
      open(unit=71,file=fylein,status='old')	!count number of database files 
      i = 0
      do while (i.ge.0)
        read(71,*,end=226)
        i = i + 1
      enddo
226   continue
      close(71)
      itmp = i-3
      end
!_______________________________________________________________________
!--
      subroutine permute(ib,ie)
      use data_files
      use randseed
      implicit none
      integer i, j, ib, ie, itmp
      double precision atmp
      
      do i=ib,ie
	call RANDOM_NUMBER(atmp)
        j = ib + nint((ie - ib) * atmp)
        itmp = lstt(i)
        lstt(i) = lstt(j)
        lstt(j) = itmp
      enddo

      end
!_______________________________________________________________________
!-- Initialize the neural network (the weights are randomized)
      subroutine netinit(nrx,nip,nw1,nw2,nvo,nrxo,ngi,ngo,iug)
      use randseed
      use NN_wei
      use NN_owei
      use NN_decays
      implicit none
      integer i,j,k, iug,jc,jc2,mc
      double precision atmp, delhi
      integer nrx,nip,nw1,nw2,nvo,nrxo,ngi,ngo
      
!**     randomizing weights
      do i=1,ngi
        do j=1,nw1
          call RANDOM_NUMBER(atmp)
          atmp = atmp - 0.5
          wgi(i,j) = atmp
          oldwgi(i,j) = atmp / 10.d0
        enddo
      enddo

      do i=1,nrx+1
        do j=1,nip
          do k=1,nw1
            call RANDOM_NUMBER(atmp)
            atmp = atmp - 0.5
            wi(i,j,k) = atmp
            oldwi(i,j,k) = atmp / 10.d0
          enddo
        enddo
      enddo

      do i=1,nw1+1
        do j=1,nw2
          call RANDOM_NUMBER(atmp)
          atmp = atmp - 0.5
          wh(i,j) = atmp
          oldwh(i,j) = atmp / 10.d0
        enddo
      enddo

      do i=1,nw2+1
        do j=1,nrxo
          do k=1,nvo
            call RANDOM_NUMBER(atmp)
            atmp = atmp - 0.5
            wo(i,j,k) = atmp
            oldwo(i,j,k) = atmp / 10.d0
          enddo
        enddo
      enddo

      do j=1,ngo
        do i=1,nw2+1
          call RANDOM_NUMBER(atmp)
          atmp = atmp - 0.5
          wgo(i,j) = atmp
          oldwgo(i,j) = atmp / 10.d0
        enddo
      enddo

!**     initialize distance matrix
      if (nw1.gt.1) then
        delhi = (nrx-1.d0)/(nw1-1.d0)
      else
        delhi = 0.d0
        dfac = 1.d0
      endif
      do i=1,nrx
        do k=1,nw1
          if (iug.eq.0) then
            dij1(i,k) = 1.d0
          else
            dij1(i,k) = dfac / dsqrt(1.d0+((k-1)*delhi-i+1)*((k-1)*delhi-i+1))
          endif
        enddo
      enddo

      jc = (nw1+1)/2
      jc2 = (nw2+1)/2
      mc = (nrxo+1)/2
      do i=1,nw1
        do j=1,nw2
          if (iug.eq.0) then
            dijh(i,j) = 1.d0
          else
            dijh(i,j) = dfac / dsqrt(1.d0+delhi*delhi*(i-jc-j+jc2)*(i-jc-j+jc2))
          endif
        enddo
      enddo

      do i=1,nw2
        do j=1,nrxo
          if (iug.eq.0) then
            dij2(i,j) = 1.d0
          else
            dij2(i,j) = dfac / dsqrt(1.d0+((i-jc)*delhi-j+mc)*((i-jc)*delhi-j+mc))
          endif
        enddo
      enddo

      end
!___________________________________________________________________
!-- Initialize dbindex angle array
      subroutine initphi(naax,nvo)
      use DB_output
      implicit none
      integer i,j
      integer naax,nvo

      do i=1,naax
        do j=1,nvo
          phiut(j,i) = 0.d0
        enddo
        icphi(i) = 0
      enddo

      end
!___________________________________________________________________
!-- Read the current window
      subroutine readwin(iwin,iglo,nip,nvo,nisf,nosf,ngi,ngo)
      use NN_input
      use NN_nodes
      implicit none
      integer i,j,irw
      integer iwin,iglo,nip,nvo,nisf,nosf,ngi,ngo

      do i=1,ngi
        gval(i) = globin(i,iglo)
      enddo
      
      do i=ibi(1,iwin),ibi(2,iwin)
        irw = iwin-nisf+i
        do j=1,nip
          resv(i,j) = protin(j,irw)
        enddo
      enddo

      do i=ibi(3,iwin),ibi(4,iwin)
        irw = iwin-nosf+i
        do j=1,nvo
          phi2(i,j) = protut(j,irw)
        enddo
      enddo
      
      do i=1,ngo
        gphi2(i) = globut(i,iglo)
      enddo

      end
!___________________________________________________________________
!-- Calculate the neural network's output
      subroutine calcnet(iwin,nrx,nosf,nip,nw1,nw2,nvo,ngi,ngo,alpha,ibij)
      use NN_wei
      use NN_decays
      use NN_nodes
      use DB_output
      implicit none
      integer i,j,k, ikan
      integer iwin,nrx,nosf,nip,nw1,nw2,nvo,ngi,ngo,ibij(4)
      double precision alpha, ztemp
      
!** Calculate nueron values of the first hidden layer
      do i=1,nw1
        ztemp=0.d0
	do j=ibij(1),ibij(2)
	  do k=1,nip
	    ztemp = ztemp + resv(j,k) * wi(j,k,i) * dij1(j,i)
	  enddo
	enddo
        ztemp = ztemp + wi(nrx+1,1,i)	!bias
        do j=1,ngi
          ztemp = ztemp + gval(j) * wgi(j,i)
        enddo
	hnu1(i) = tanh(alpha*ztemp)
      enddo

      do i=1,nw2
        ztemp=0.d0
	do j=1,nw1
          ztemp = ztemp + hnu1(j) * wh(j,i) * dijh(j,i)
	enddo
        ztemp = ztemp + wh(nw1+1,i)
	hnu2(i) = tanh(alpha*ztemp)
      enddo

!** Calculate output nueron values
      do i=ibij(3),ibij(4)
        ikan = iwin+i-nosf
        do k=1,nvo
          ztemp=0.d0
          do j=1,nw2
            ztemp = ztemp + hnu2(j) * wo(j,i,k) * dij2(j,i)
          enddo
          ztemp = ztemp + wo(nw2+1,i,k)
          phi2p(i,k) = tanh(alpha*ztemp)
          phiut(k,ikan) = phiut(k,ikan) + phi2p(i,k)
        enddo
        icphi(ikan) = icphi(ikan) + 1
      enddo
      
      do i=1,ngo
        ztemp=0.d0
        do j=1,nw2
          ztemp = ztemp + hnu2(j) * wgo(j,i)
        enddo
        ztemp = ztemp + wgo(nw2+1,i)
        gphi2p(i) = tanh(alpha*ztemp)
        gphiut(i,iwin) = gphi2p(i)
      enddo

      end
!___________________________________________________________________
!-- This was the tanh transfer function.  f(x) = 1/(1+exp(alpha*x)) is another choice
!-- Note: the derivative f'(x) = alpha*f(x)*(1-f(x)) should be used then in calcnet/f
!      double precision function tfunc(xoutp)
!      implicit none
!      double precision xoutp
!
!      tfunc = tanh(xoutp)
!
!      return
!      end function tfunc
!___________________________________________________________________
!-- Calculate error over DB
      subroutine calcierr(trer,ofer,teer,icond,cormax)
      use all_sizes
      use NN_input
      use data_files
      use DB_output
      implicit none
      integer i,j,ip,iglo, imax
      double precision trer,ofer,teer,cormax, gtrer, gofer, gteer
      integer icond

      trer = 0.d0; ofer = 0.d0; teer = 0.d0
      gtrer = 0.d0; gofer = 0.d0; gteer = 0.d0
      if (mpout.lt.2) then		!GENERAL OUTPUT
        do i=1,ntraa
          ip = lstt(i)
          iglo = lsglo(ip)
          do j=1,nvo
            trer = trer + dabs(protut(j,ip)-phiut(j,ip))
          enddo
          do j=1,ngo
            gtrer = gtrer + dabs(globut(j,iglo) - gphiut(j,ip))
          enddo
        enddo
        trer = trer / (nvo * ntraa)
        if (ngo.gt.0) then
          gtrer = gtrer / (ngo * ntraa)
        else
          gtrer = 0.d0
        endif
        trer = trer + gtrer
        do i=ntraa+1,ntsta
          ip = lstt(i)
          iglo = lsglo(ip)
          do j=1,nvo
            ofer = ofer + dabs(protut(j,ip)-phiut(j,ip))
          enddo
          do j=1,ngo
            gofer = gofer + dabs(globut(j,iglo) - gphiut(j,ip))
          enddo
        enddo
        ofer = ofer / (nvo * (ntsta - ntraa))
        if (ngo.gt.0) then
          gofer = gofer / (ngo * (ntsta - ntraa))
        else
          gofer = 0.d0
        endif
        ofer = ofer + gofer
        if (ofer.lt.cormax) then
          icond = 1
        else
          icond = 0
        endif
        do i=ntsta+1,naax
          ip = lstt(i)
          iglo = lsglo(ip)
          do j=1,nvo
            teer = teer + dabs(protut(j,ip)-phiut(j,ip))
          enddo
          do j=1,ngo
            gteer = gteer + dabs(globut(j,iglo) - gphiut(j,ip))
          enddo
        enddo
        if (naax.gt.ntsta) then
          teer = teer / (nvo * (naax - ntsta))
        else
          teer = 0.d0
        endif
        if (ngo.gt.0) then
          gteer = gteer / (ngo * (naax - ntsta))
        else
          gteer = 0.d0
        endif
        teer = teer + gteer
      else				!PROBABILITY OUTPUT
        do i=1,ntraa
          ip = lstt(i)
          iglo = lsglo(ip)
          if (imax(nvo,protut(1,ip)).eq.imax(nvo,phiut(1,ip))) trer = trer + 1.d0
          do j=1,ngo
            gtrer = gtrer + dabs(globut(j,iglo) - gphiut(j,ip))
          enddo
        enddo
        trer = trer / ntraa
        if (ngo.gt.0) then
          gtrer = gtrer / (ngo * ntraa)
        else
          gtrer = 0.d0
        endif
        trer = trer - gtrer
        do i=ntraa+1,ntsta
          ip = lstt(i)
          iglo = lsglo(ip)
          if (imax(nvo,protut(1,ip)).eq.imax(nvo,phiut(1,ip))) ofer = ofer + 1.d0
          do j=1,ngo
            gofer = gofer + dabs(globut(j,iglo) - gphiut(j,ip))
          enddo
        enddo
        ofer = ofer / (ntsta - ntraa)
        if (ngo.gt.0) then
          gofer = gofer / (ngo * (ntsta - ntraa))
        else
          gofer = 0.d0
        endif
        ofer = ofer - gofer
        if (ofer.gt.cormax) then
          icond = 1
        else
          icond = 0
        endif
        do i=ntsta+1,naax
          ip = lstt(i)
          iglo = lsglo(ip)
          if (imax(nvo,protut(1,ip)).eq.imax(nvo,phiut(1,ip))) teer = teer + 1.d0
          do j=1,ngo
            gteer = gteer + dabs(globut(j,iglo) - gphiut(j,ip))
          enddo
        enddo 
        if (naax.gt.ntsta) then
          teer = teer / (naax - ntsta)
        else
          teer = 0.d0
        endif
        if (ngo.gt.0) then
          gteer = gteer / (ngo * (naax - ntsta))
        else
          gteer = 0.d0
        endif
        teer = teer - gteer
      endif

      end
!___________________________________________________________________
!-- general error calculator
      subroutine calcgerr(istrt,istop,trer)
      use all_sizes
      use NN_input
      use data_files
      use DB_output
      implicit none
      integer i,ip,iglo,j,imax
      double precision trer, gtrer
      integer istrt, istop
      trer = 0.d0;
      gtrer = 0.d0;
      if (mpout.lt.2) then		!GENERAL OUTPUT
        do i=istrt,istop
          ip = i
          iglo = lsglo(ip)
          do j=1,nvo
            trer = trer + dabs(protut(j,ip)-phiut(j,ip))
          enddo
          do j=1,ngo
            gtrer = gtrer + dabs(globut(j,iglo) - gphiut(j,ip))
          enddo
        enddo
        if (istop.ge.istrt) then
          trer = trer / (nvo * (istop - istrt + 1))
        else
          trer = 0.d0
        endif
        if (ngo.gt.0) then
          gtrer = gtrer / (ngo * (istop - istrt + 1))
        else
          gtrer = 0.d0
        endif
        trer = trer + gtrer
      else				!PROBABILITY OUTPUT
        do i=istrt,istop
          ip = i
          iglo = lsglo(ip)
          if (imax(nvo,protut(1,ip)).eq.imax(nvo,phiut(1,ip))) trer = trer + 1.d0
          do j=1,ngo
            gtrer = gtrer + dabs(globut(j,iglo) - gphiut(j,ip))
          enddo
        enddo
        if (istop.ge.istrt) then
          trer = trer / (istop - istrt + 1)
        else
          trer = 0.d0
        endif
        if (ngo.gt.0) then
          gtrer = gtrer / (ngo * (istop - istrt + 1))
        else
          gtrer = 0.d0
        endif
        trer = trer - gtrer
      endif

      end
!___________________________________________________________________
!-- Back propogation
      subroutine backprop(nrx,nip,nw1,nw2,nrxo,nvo,ngi,ngo,alpha,ratel,rmoment,ibij)
      use NN_wei
      use NN_owei
      use NN_decays
      use NN_nodes
      implicit none
      integer ibegi,iendi,ibego,iendo, i,j,k, kk
      integer nrx,nip,nw1,nw2,nrxo,nvo,ngi,ngo,ibij(4)
      double precision alpha,ratel,rmoment, errotmp, ztemp
      
      double precision errh1(nrx,nw1),birh1(nw1),errh2(nw1,nw2),birh2(nw2),biro2(nrxo,nvo),erro2(nw2,nrxo,nvo)
      double precision, allocatable :: gerro2(:)
      
      if (ngo.gt.0) allocate( gerro2(ngo) )
      
      ibegi = ibij(1)
      iendi = ibij(2)
      ibego = ibij(3)
      iendo = ibij(4)
!** Calculate errors at the output layer and the RMSE
      do i=ibego,iendo
        do k=1,nvo
          errotmp = phi2(i,k) - phi2p(i,k)
          biro2(i,k) = errotmp * alpha * (1.d0-phi2p(i,k)*phi2p(i,k))
        enddo

        do j=1,nw2
          do k=1,nvo
            erro2(j,i,k) = biro2(i,k) * dij2(j,i)
          enddo
        enddo
      enddo
      
      do i=1,ngo
        errotmp = gphi2(i) - gphi2p(i)
        gerro2(i) = errotmp * alpha * (1.d0-gphi2p(i)*gphi2p(i))
      enddo
!** Calculate the errors at the second hidden layer
      do k=1,nw2
        ztemp=0.d0
        do j=ibego,iendo
          do kk=1,nvo
  	    ztemp = ztemp + wo(k,j,kk) * erro2(k,j,kk)
          enddo
        enddo
        do j=1,ngo
          ztemp = ztemp + wgo(k,j) * gerro2(j)
        enddo
        birh2(k) = ztemp * alpha  * (1 - hnu2(k) * hnu2(k))
        do j=1,nw1
	  errh2(j,k) = birh2(k) * dijh(j,k)
        enddo
      enddo

!** Calculate the errors at the first hidden layer
      do k=1,nw1
        ztemp=0.d0
        do j=1,nw2
	  ztemp = ztemp + wh(k,j)*errh2(k,j)
        enddo
        birh1(k) = ztemp * alpha  * (1 - hnu1(k) * hnu1(k))
        do j=ibegi,iendi
	  errh1(j,k) = birh1(k) * dij1(j,k)
        enddo

!** Update weights:
!** Update input/hidden layer weights (wi)
        do i=ibegi,iendi
          do j=1,nip
            ztemp = ratel*errh1(i,k)*resv(i,j) + rmoment*oldwi(i,j,k)
            wi(i,j,k) = wi(i,j,k) + ztemp
            oldwi(i,j,k) = ztemp
          enddo
        enddo
        
        do i=1,ngi
          ztemp = ratel*birh1(k)*gval(i) + rmoment*oldwgi(i,k)
          wgi(i,k) = wgi(i,k) + ztemp
          oldwgi(i,k) = ztemp
        enddo

!** Update input bias weights (wi(nrx+1,1,1-nw1))
        ztemp = ratel*birh1(k) + rmoment*oldwi(nrx+1,1,k)
        wi(nrx+1,1,k) = wi(nrx+1,1,k) + ztemp
        oldwi(nrx+1,1,k) = ztemp
      enddo

!** Update hidden layer weights (wh)
      do k=1,nw2
        do j=1,nw1
          ztemp = ratel*errh2(j,k)*hnu1(j) + rmoment*oldwh(j,k)
          wh(j,k) = wh(j,k) + ztemp
          oldwh(j,k) = ztemp
        enddo

!** Update hidden bias weights (wh(nw1+1,1-nw2))
        ztemp = ratel*birh2(k) + rmoment*oldwh(nw1+1,k)
        wh(nw1+1,k) = wh(nw1+1,k) + ztemp
        oldwh(nw1+1,k) = ztemp

!** Update hidden/output layer weights (wo)
        do i=ibego,iendo
          do kk=1,nvo
            ztemp = ratel*erro2(k,i,kk)*hnu2(k) + rmoment*oldwo(k,i,kk)
            wo(k,i,kk) = wo(k,i,kk) + ztemp
            oldwo(k,i,kk) = ztemp
          enddo
        enddo
        
        do i=1,ngo
          ztemp = ratel*gerro2(i)*hnu2(k) + rmoment*oldwgo(k,i)
          wgo(k,i) = wgo(k,i) + ztemp
          oldwgo(k,i) = ztemp
        enddo
      enddo

!** Update hidden/output layer bias weights (wo)
        k=nw2+1
        do i=ibego,iendo
          do kk=1,nvo
            ztemp = ratel*biro2(i,kk) + rmoment*oldwo(k,i,kk)
            wo(k,i,kk) = wo(k,i,kk) + ztemp
            oldwo(k,i,kk) = ztemp
          enddo
        enddo
        do i=1,ngo
          ztemp = ratel*gerro2(i) + rmoment*oldwgo(k,i)
          wgo(k,i) = wgo(k,i) + ztemp
          oldwgo(k,i) = ztemp
        enddo

      end
!___________________________________________________________________
!-- Write weights
      subroutine outwei(tmpfl,nrx,nip,nw1,nw2,nvo,nrxo,ngi,ngo,iug,alpha)
      use data_files
      use NN_wei
      implicit none
      integer i,j,k 
      integer ifl,nrx,nip,nw1,nw2,nrxo,nvo,ngi,ngo,iug
      double precision alpha
      character tmpfl*3600

      ifl = 77
      open(unit=ifl,file=tmpfl,status='unknown')
      write(ifl,'(a,a,a)') '"',trim(datdir),'"'
      write(ifl,'(100(a,1x))') (trim(sinfyles(j)),j=1,nsinp)
      write(ifl,'(100(a,1x))') (trim(sutfyles(j)),j=1,nsutp)
      write(ifl,'(9(I12),f14.5)') nrx,nip,nw1,nw2,nvo,nrxo,ngi,ngo,iug,alpha
      do i=1,nrx+1
	do j=1,nip
	  do k=1,nw1
	    write(ifl,*) wi(i,j,k)
	  enddo
	enddo
      enddo
      do i=1,nw1+1
        do j=1,nw2
	  write(ifl,*) wh(i,j)
        enddo
      enddo
      do i=1,nw2+1
        do j=1,nrxo
          do k=1,nvo
	    write(ifl,*) wo(i,j,k)
          enddo
        enddo
      enddo
      do i=1,ngi
        do j=1,nw1
          write(ifl,*) wgi(i,j)
        enddo
      enddo
      do i=1,ngo
        do j=1,nw2+1
          write(ifl,*) wgo(j,i)
        enddo
      enddo
      close(ifl)

      end
!___________________________________________________________________
!-- Read weights
      subroutine readwei(tmpfl,nrx,nip,nw1,nw2,nvo,nrxo,ngi,ngo,iug,alpha)
      use NN_wei
      implicit none
      integer i,j,k 
      integer ifl,nrx,nip,nw1,nw2,nrxo,nvo,ngi,ngo,iug,nipc,nvoc
      double precision alpha
      character tmpfl*3600

      ifl = 77
      open(unit=ifl,file=tmpfl,status='unknown')
      read(ifl,*)
      read(ifl,*)
      read(ifl,*)
      read(ifl,'(9(I12),f14.5)') nipc,nip,nw1,nw2,nvo,nvoc,ngi,ngo,iug,alpha
      if ((nipc.ne.nrx).or.(nvoc.ne.nrxo)) then
        write(0,*) "Different number of inputs/output b/w data and weight files, aborting.",nrx,nipc,nrxo,nvoc
        stop
      endif
      do i=1,nrx+1
	do j=1,nip
	  do k=1,nw1
	    read(ifl,*) wi(i,j,k)
	  enddo
	enddo
      enddo
      do i=1,nw1+1
        do j=1,nw2
	  read(ifl,*) wh(i,j)
        enddo
      enddo
      do i=1,nw2+1
        do j=1,nrxo
          do k=1,nvo
	    read(ifl,*) wo(i,j,k)
          enddo
        enddo
      enddo
      do i=1,ngi
        do j=1,nw1
          read(ifl,*) wgi(i,j)
        enddo
      enddo
      do i=1,ngo
        do j=1,nw2+1
          read(ifl,*) wgo(j,i)
        enddo
      enddo
      close(ifl)

      end
!___________________________________________________________________
!-- Initialize input for filter NN
      subroutine savinp(naax,nvo)
      use DB_output
      implicit none
      integer i,j
      integer naax,nvo

      do i=1,naax
        do j=1,nvo
          phiutf(j,i) = phiut(j,i)
        enddo
      enddo

      end
!___________________________________________________________________
!-- Read the filter window
      subroutine readwinf(iwin,iglo,nip,nvo,nisf,nosf,ngi,ngo)
      use NN_input
      use NN_nodes
      use DB_output
      implicit none
      integer i,j,irw
      integer iwin,iglo,nip,nvo,nisf,nosf,ngi,ngo

      do i=1,ngi
        gval(i) = globin(i,iglo)
      enddo
      
      do i=ibi(1,iwin),ibi(2,iwin)
        irw = iwin-nisf+i
        do j=1,nip
          resv(i,j) = phiutf(j,irw)
        enddo
      enddo

      do i=ibi(3,iwin),ibi(4,iwin)
        irw = iwin-nosf+i
        do j=1,nvo
          phi2(i,j) = protut(j,irw)
        enddo
      enddo

      do i=1,ngo
        gphi2(i) = globut(i,iglo)
      enddo

      end
!___________________________________________________________________
!-- get integer time
      integer function inttget(iflg)
      implicit none
      integer ita3,jta3,ita4log,ita4,krndm
      double precision ta1,ta2,ta3,ta4,dta4
      integer iflg
      integer*4 timeArray(3)    ! Holds the hour, minute, and second

      krndm = 22
      if (iflg.eq.0) then
        inttget = 7
      else
        call itime(timeArray)     ! Get the current time for seed
        ta1 = timeArray(1) / 1.74
        ta2 = timeArray(2) / 4.35
        ta3 = timeArray(3) / 2.8
        jta3 = timeArray(3)
        ita3 = exp(ta3)
        ta3 = (mod(ita3,61)/2.8)
        ta4 = exp(ta1) + exp(ta2) + exp(ta3) + krndm
        ita4log = nint(log10(ta4))
        ita4 = ta4
        dta4 = 10**ita4log * (ta4 - ita4)
        inttget = abs(dta4 - ta4) + ta3 + krndm
!        inttget = 64021
      endif
        

      return
      end
!___________________________________________________________________
!-- argmax index
      integer function imax(isov,vector)
      implicit none
      integer i
      integer isov
      double precision vector(isov)
      
      imax = 1;
      do i=2,isov
        if (vector(i).gt.vector(imax)) imax = i
      enddo

      return
      end
!___________________________________________________________________
!-- Normalize based on std (z-score)
      subroutine stdnorm(isov,vector,iounit)
      implicit none
      integer i, isov, iounit
      double precision vector(isov), davr, dstd, mavrg, mstdv
      
      davr = mavrg(isov,vector)
      dstd = mstdv(isov,vector,davr)
      write(*,*) "globnorm(avr/std): ", davr, dstd
      if (dstd.gt.0.d0) then
        do i=1,isov
          vector(i) = (vector(i) - davr) / dstd
        enddo
      else
        do i=1,isov
          vector(i) = (vector(i) - davr)
        enddo
      endif

      return
      end
!___________________________________________________________________
!-- Normalize based on max/min
      subroutine maxminnorm(isov,vector,iounit)
      implicit none
      integer i, isov, iounit
      double precision vector(isov), dmax, dmin
      
      dmax = vector(1)
      dmin = vector(1)
      do i=2,isov
        if (vector(i).gt.dmax) dmax = vector(i)
        if (vector(i).lt.dmin) dmin = vector(i)
      enddo
      do i=1,isov
        vector(i) = 2.d0 * ((vector(i) - dmin) / (dmax - dmin) - 0.5d0)
      enddo
      write(*,*) "globnorm(min/max): ", dmin, dmax

      return
      end
!___________________________________________________________________
!-- Average function
      double precision function mavrg(isov,vector)
      implicit none
      integer i, isov
      double precision vector(isov)
      
      mavrg = vector(1);
      do i=2,isov
        mavrg = mavrg + vector(i)
      enddo
      
      if (isov.gt.0) then
        mavrg = mavrg / isov
      else
        mavrg = 0.d0
      endif

      return
      end
!___________________________________________________________________
!-- STDEV function
      double precision function mstdv(isov,vector,dvar)
      implicit none
      integer i, isov
      double precision vector(isov), dvar
      
      mstdv = (vector(1)-dvar)*(vector(1)-dvar)
      do i=2,isov
        mstdv = mstdv + (vector(i)-dvar)*(vector(i)-dvar)
      enddo
      
      if (isov.gt.0) then
        mstdv = mstdv / isov
        mstdv = sqrt(mstdv)
      else
        mstdv = 1.d0
      endif

      return
      end
!___________________________________________________________________
!-- Normalize over nrxo window (averaging) can add weights based on distance from central residue
      subroutine normphi(naax,nvo)
      use DB_output
      implicit none
      integer i,j
      integer naax,nvo

      do i=1,naax
        do j=1,nvo
          phiut(j,i) = phiut(j,i) / icphi(i)
        enddo
      enddo

      end
!___________________________________________________________________


