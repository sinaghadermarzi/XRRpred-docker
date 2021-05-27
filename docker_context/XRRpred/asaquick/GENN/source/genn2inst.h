!       GENN V.2 -- SEDER MODULE
!       Eshel Faraggi 04/2011 (c)
!!! PLEASE DO NOT USE THIS PROGRAM/CODE YET!

module all_sizes
  integer nte,nip,nop,nw1,nw2
  integer ntotp,ntrap,jdf, itbmax, dbv
  double precision fcv,alpha,ratel,rmoment
end module

module NN_input
  integer, allocatable :: iprotdg(:)
  double precision, allocatable :: protin(:), protut(:) ! sequence inputs / outputs
end module

module data_files
  integer ioffl
  character*3600 datdir, fylein, arg, tmpfl, weifl, oweifl, weiadd, prdfl, prdlist, weiav, foldpart
  character*3600 sinfyles(100), sutfyles(100)
  character, allocatable :: protid(:)*3600
  integer, allocatable :: lstt(:), lstt2(:)	! map of training testing
  integer nsinp,nsutp
end module

module randseed
   integer itsd
end module

module NN_wei
  double precision, allocatable :: wi(:,:), wh(:,:), wo(:,:)	! Main Sequence weights
end module

module NN_wei_sv
  double precision, allocatable :: wisv(:,:), whsv(:,:), wosv(:,:) ! save weights
end module

module NN_owei					! Old weights (momentum calculation)
  double precision, allocatable :: oldwi(:,:), oldwh(:,:), oldwo(:,:)
end module

module NN_nodes ! Main nodes
  double precision, allocatable :: hnu1(:), hnu2(:), phi2p(:)
end module

module DB_output ! database indexed predictions
  double precision, allocatable :: phiut(:)
end module

      subroutine seder2train()
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
      integer j,itb,ne,ix,ip,iq,icond,ipi,ipo, nopti, nav
      double precision trer, ofer,cormax, tramax
      double precision, allocatable :: sphiut(:,:)
      character tmpfl2*3600,outfl*3600

      if (weiav.ne."") then
	if (prdfl.eq."") then
	  write(0,*) "Must have -pr1/l with -aw, aborting."
	  stop
	endif
        open(unit=29,file=weiav,status='old')
        allocate( sphiut(nop*ntotp,2) )
        sphiut = 0.d0
        nav = 0
      endif
      if (prdfl.ne."") then
        write(0,*) "GEneral Neural Network (predict), V2.2 2011 (c) -h for help"
        write(0,*) "Eshel Faraggi, efaraggi@mamiris.com"
        write(0,*) "Research and Information Systems"
        write(0,*) "You Must Obtain a License From Above Address to use this program"
        write(0,*) "Weights file: ", trim(weifl)
        write(0,'(200A15)') "INS/OUTS: ", (trim(sinfyles(j)),j=1,nsinp), "/", (trim(sutfyles(j)),j=1,nsutp)
      endif
107   continue
      if (weiav.ne."") then
        read(29,*,end=117) weifl
        if (dbv.gt.0) call init_all
        nav = nav + 1
      endif

      write(tmpfl2,'(a,i6.6,a,i11.11,a)') "train.", ntotp, "_", itsd, "."
      tmpfl2 = trim(tmpfl2)//trim(weiadd)
      call netinit(nip,nw1,nw2,nop)
      if (weifl.eq."") then
        weifl = trim(tmpfl2)//".wei"
      else
	call readwei(weifl,nip,nw1,nw2,nop,alpha)
        if (prdfl.eq."") weifl = trim(tmpfl2)//".wei"
      endif
!      call permute(ntotp,1,ntotp,lstt(1))

      if (prdfl.eq."") then
        itb = 0
      else
        do ip=1,ntotp
          ipi = 1 + nip * (ip - 1)
          ipo = 1 + nop * (ip - 1)
          call calcnet(ip,nip,nw1,nw2,nop,alpha,protin(ipi))
          if (weiav.ne."") then
            do ix=ipo,ipo+nop-1
              sphiut(ix,1) = sphiut(ix,1) + phiut(ix) !write(*,'(a,1000(f14.5,1x))') trim(protid(ip)), (phiut(ix), protut(ix),ix=ipo,ipo+nop-1)
              sphiut(ix,2) = sphiut(ix,2) + phiut(ix)*phiut(ix)
            enddo
          endif
        enddo
      endif
      if (weiav.ne."") goto 107
117   continue
      if (prdfl.ne."") then
        if (weiav.ne."") then
          close(29)
          do ip=1,ntotp
            ipo = 1 + nop * (ip - 1)
            do ix=ipo,ipo+nop-1
              phiut(ix) = sphiut(ix,1) / nav
              sphiut(ix,2) = dsqrt(sphiut(ix,2)/nav - phiut(ix)*phiut(ix))
            enddo
          enddo
        endif
        do ip=1,ntotp
          ipo = 1 + nop * (ip - 1)
          do ix=ipo,ipo+nop-1
            if (weiav.ne."") then
              write(*,'(A,1x,3G12.4)') trim(protid(ip)), phiut(ix), sphiut(ix,2), protut(ix)
            else
              write(*,'(A,1x,2G12.4)') trim(protid(ip)), phiut(ix), protut(ix)
            endif
          enddo
        enddo
        stop
      endif

      if (oweifl.ne."") weifl = oweifl

      cormax = 20.d0
      outfl = trim(weifl)//".fort.15"
      open(unit=ioffl,file=outfl,status='new')
      write(ioffl,*) "GEneral Neural Network, V3.1 (train) 2011 (c) -h for help"
      write(ioffl,*) "Eshel Faraggi, efaraggi@gmail.com"
      write(ioffl,*) "Research and Information Systems"
      write(ioffl,*) "You Must Obtain a License From Above Address to Use GENN"
      write(ioffl,*) "Options: ", trim(weiadd)
      write(ioffl,*) "Weights file: ", trim(weifl)
      write(ioffl,'(a2,a,i15)') "# ", "Random Seed: ", itsd
      write(ioffl,'(200A15)') "INS/OUTS: ", (trim(sinfyles(j)),j=1,nsinp), "/", (trim(sutfyles(j)),j=1,nsutp)
      write(ioffl,*) trim(foldpart)

      ne = 1
!!!!! training epochs (main network):
      do while ((ne.le.nte).and.(itb.lt.itbmax))
        call permute(ntrap,1,ntrap,lstt2(1))
        do ix=1,ntrap
          ip = lstt2(ix)
          ip = lstt(ip)
          ipi = 1 + nip * (ip - 1)
          ipo = 1 + nop * (ip - 1)
          do iq=1,iprotdg(ip)
            call calcnet(ip,nip,nw1,nw2,nop,alpha,protin(ipi))
            call backprop(nip,nw1,nw2,nop,alpha,ratel,rmoment,protin(ipi),protut(ipo))
          enddo
        enddo
        do ix=ntrap+1,ntotp
          ip = lstt(ix)
          ipi = 1 + nip * (ip - 1)
          call calcnet(ip,nip,nw1,nw2,nop,alpha,protin(ipi))
        enddo

        call calcierr(trer,ofer,icond,cormax)
        if ((icond.eq.1).or.(ne.le.1)) then
          call savewei()
          tramax = trer
          cormax = ofer
          nopti = ne
          itb = 0
        else
          itb = itb + 1
        endif
        write(ioffl,'(I10,100(F11.4))') ne,trer,ofer

        ne = ne + 1
      enddo
      write(ioffl,'(a,i7,a,d12.4,a,d12.4)') "Min error at ",nopti, ": training, ",tramax,"; overfit, ",cormax
      close(ioffl)

      call writewei(weifl,nip,nw1,nw2,nop,alpha)

      end
!_______________________________________________________________________
!--
      subroutine readcmd
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
      integer inttget,nargs,icrd,K,j
      integer, allocatable :: seeds(:)

!       Defaults
      fylein = "genn.in"	 !names of input files
      datdir = ""		 !directory where input files are
      weifl = ""                 !weight file
      oweifl = ""                !weights output file
      prdfl = ""                 !id for input files for prediction
      prdlist = ""               !file of ids for prediction
      nte = 1000		 !maximum number of epochs
      ntotp = 0			 !total number of input files to use (else use all)
      itsd = inttget(7)	 !seed for random number generator (0 - constant seed (debugging))
      fcv = 0.05d0		 !fraction of dataset to for overfitprotection set
      nw1 = 21			 !nodes in hidden layer 1 (hl1)
      nw2 = 21			 !nodes in hidden layer 2 (hl2)
      alpha = 0.2		 !activation parameter
      itbmax = 100		 !maximum number of non-improving training epochs
      ratel = 0.001		 !learning rate
      rmoment = 0.4		 !momentum
      jdf = 0			 !use degeneracy file
      ioffl = 15		 !file number for progress file
      weiadd = ""
      weiav = ""
      dbv = 0 			 !different db between different weights (prediction). Can save time not reading db everytime

      nargs=iargc()
      if (nargs.ge.1) then
        icrd = 1
!       Switches:
        do while (icrd.le.nargs)
          call getarg(icrd,arg)
          if (arg.eq."-l") then			 !db_list_file
            call getarg(icrd+1,tmpfl)
            fylein = trim(tmpfl)
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-d") then		 !data_dir
            call getarg(icrd+1,tmpfl)
            datdir = trim(tmpfl)
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-wf") then		 !weight file
            call getarg(icrd+1,tmpfl)
            weifl = trim(tmpfl)
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
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
          elseif (arg.eq."-f") then		 !fraction_of_overfit_set_of_db
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) fcv
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-h1") then		 !num_hid_layer1_nodes
            call getarg(icrd+1,tmpfl)
            weiadd = trim(weiadd)//trim(arg)//"_"//trim(tmpfl)
            read(tmpfl,*) nw1
          elseif (arg.eq."-h2") then		 !num_hid_layer2_nodes
            call getarg(icrd+1,tmpfl)
            weiadd = trim(weiadd)//trim(arg)//"_"//trim(tmpfl)
            read(tmpfl,*) nw2
          elseif (arg.eq."-mi") then		 !max_num_bad_iterations
            call getarg(icrd+1,tmpfl)
            read(tmpfl,*) itbmax
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
          elseif (arg.eq."-a") then		 !activation_parameter
            call getarg(icrd+1,tmpfl)
            weiadd = trim(weiadd)//trim(arg)//"_"//trim(tmpfl)
            read(tmpfl,*) alpha
          elseif (arg.eq."-u") then		 !learning_rate
            call getarg(icrd+1,tmpfl)
            weiadd = trim(weiadd)//trim(arg)//"_"//trim(tmpfl)
            read(tmpfl,*) ratel
          elseif (arg.eq."-p") then		 !momentum
            call getarg(icrd+1,tmpfl)
            weiadd = trim(weiadd)//trim(arg)//"_"//trim(tmpfl)
            read(tmpfl,*) rmoment
          elseif (arg.eq."-df") then		 !Use_degeneracy_files
            jdf = 1
            weiadd = trim(weiadd)//trim(arg)
          elseif (arg.eq."-dw") then		 !Different input structure for weight files
            dbv = 1
            weiadd = trim(weiadd)//trim(arg)
          elseif (arg.eq."-aw") then		 !average over weights in file
            call getarg(icrd+1,tmpfl)
            weiav = trim(tmpfl)
            open(unit=33,file=weiav,status='old')
            read(33,*) weifl
            close(33)
            fylein = ""
            weiadd = trim(weiadd)//trim(arg)//trim(tmpfl)
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
      use NN_wei_sv
      use NN_owei
      use NN_nodes
      use DB_output
      implicit none

      if ((fylein.ne."").or.(weifl.ne."")) then
        call initdb()
        if (allocated(wi)) then
          deallocate( wi , wh , wo )
          deallocate( wisv , whsv , wosv )
          deallocate( oldwi , oldwh , oldwo )
          deallocate( hnu1 , hnu2 , phi2p )
          deallocate( phiut )
        endif
        allocate( wi(nip+1,nw1) , wh(nw1+1,nw2) , wo(nw2+1,nop) )
        allocate( wisv(nip+1,nw1) , whsv(nw1+1,nw2) , wosv(nw2+1,nop) )
        allocate( oldwi(nip+1,nw1) , oldwh(nw1+1,nw2) , oldwo(nw2+1,nop) )
        allocate( hnu1(nw1) , hnu2(nw2) , phi2p(nop) )
        allocate( phiut(nop*ntotp) )
      endif

      end
!_______________________________________________________________________
!--
      subroutine printhelp()
      use all_sizes
      use data_files
      implicit none
      
      print*, "Options:"

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
          write(*,'(A10,A40,D20.6)') "-f ","fraction_of_overfit_set_of_db  ", fcv
          write(*,'(A10,A40,I20)') "-h1","num_hid_layer1_nodes             ", nw1
          write(*,'(A10,A40,I20)') "-h2","num_hid_layer2_nodes             ", nw2
          write(*,'(A10,A40,I20)') "-mi","max_num_bad_iterations           ", itbmax
          write(*,'(A10,A40,D20.6)') "-a ","activation_parameter             ", alpha
          write(*,'(A10,A40,D20.6)') "-u ","learning_rate                    ", ratel
          write(*,'(A10,A40,D20.6)') "-p ","momentum                         ", rmoment
          write(*,'(A10,A40,A20)') "-df","Use_degeneracy_files             "
          write(*,'(A10,A40,A20)') "-aw","average over weights in file     "
          write(*,'(A10,A40,A20)') "-dw","diff nn struc for diff wei files "
          write(*,'(A10,A40,A20)') "-h ","print_this_help                  "
      write(*,*) 
      write(*,*) "Train NN prediction. db_list_file first line is database directory,"
      write(*,*) "second line contains input types (files in ids dirs bellow) third line"
      write(*,*) "contains output types. Rest of lines are protein ids (7 char max) for"
      write(*,*) "subdirs in database directory (ids dirs). "
      write(*,*) "Degeneracy files should be stored in id-dirs under file name genn.dgn"

      end
!_______________________________________________________________________
!--
      subroutine initdb()
      use all_sizes
      use NN_input
      use data_files
      implicit none
      character arbig*3600
      integer i,j,k,numval,icountel,jj,kk, nipln, nopln
      double precision dtmp(1000)
      
      if (prdfl.ne."") then
        if (prdlist.eq."") then
          ntotp = 1
        else
          fylein = trim(prdlist)
          call numpro(i)		!get number of database files
          if ((ntotp.le.0).or.(ntotp.gt.i)) then
            ntotp = i + 3
          endif
        endif
        fylein = trim(weifl)
      else
        call numpro(i)		!get number of database files
        if ((ntotp.le.0).or.(ntotp.gt.i)) then
          ntotp = i
        endif
      endif
      numval = nint(fcv*ntotp)
      ntrap = ntotp-numval

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

      if (prdfl.ne."") then
        read(22,*) nip,nw1,nw2,nop,alpha
        if (prdlist.eq."") then
          arg = prdfl
        else
          close(22)
          open(unit=22,file=prdlist,status='old')
          read(22,'(A)') arg
        endif
      else
        read(22,'(A)') arg
      endif

      j = 0;
      nip = 0
      do k=1,nsinp
        tmpfl = trim(datdir)//"/"//trim(arg)//"/"//sinfyles(k)
        open(unit=33,file=tmpfl,status='old')
        read(33,'(A)',end=336) arbig
        nipln = icountel(arbig)
        backspace(33)
        do 335 while (j.ge.0)
          read(33,'(A)',end=336) arbig
          nip = nip + nipln	!number of inputs
335     continue
336     continue
        close(33)
      enddo

      nop = 0
      do k=1,nsutp
        tmpfl = trim(datdir)//"/"//trim(arg)//"/"//sutfyles(k)
        open(unit=33,file=tmpfl,status='old')
        read(33,'(A)',end=339) arbig
        nopln = icountel(arbig)
        backspace(33)
        do 338 while (j.ge.0)
          read(33,'(A)',end=339) arbig
          nop = nop + nopln	!number of outputs
338     continue
339     continue
        close(33)
      enddo
      
      if ((nip.lt.1).or.(nop.lt.1)) then
        write(0,*) "Something wrong with database (nip/nop), aborting.",nip,nop
        stop
      endif
      if ((ntotp.lt.20).and.(prdfl.eq."")) then
        write(0,*) "Something wrong with database (not enough proteins for training), aborting.", ntotp
        stop
      endif
      write(foldpart,'(a2,a,i10,i10,i10)') "# ", "DB: nip,nop,ntotp = ",nip, nop, ntotp

      if (allocated(lstt)) then
        deallocate( lstt,lstt2,protin,iprotdg,protut,protid )
      endif
      allocate( lstt(ntotp), lstt2(ntrap) )
      allocate( protin(nip*ntotp), iprotdg(ntotp), protut(nop*ntotp), protid(ntotp) )
      do k=1,ntotp
        lstt(k) = k
      enddo
      lstt2 = (/ (k, k=1,ntrap) /)
      
      rewind(22)
      if (prdfl.eq."") then
        read(22,*); read(22,*); read(22,*)
      endif

      jj = 0
      kk = 0
      i = 1
      do while(i.le.ntotp)
        if ((prdfl.eq."").or.(prdlist.ne."")) then
          read(22,'(A)') arg
        else
          arg = prdfl
        endif
        do k=1,nsinp
          tmpfl = trim(datdir)//"/"//trim(arg)//"/"//sinfyles(k)
          open(unit=33,file=tmpfl,status='old')
          read(33,'(A)') arbig
          do while (arbig(1:1).eq."#") 
            read(33,'(A)') arbig
          enddo
          backspace(33)
          read(33,'(A)') arbig
          nipln = icountel(arbig)
          backspace(33)

	  do while (jj.ge.0)
            read(33,*,end=297) (dtmp(j),j=1,nipln)
            do j=1,nipln
              jj = jj + 1
              protin(jj) = dtmp(j)
            enddo
          enddo
297       continue
          close(33)
        enddo
        
        do k=1,nsutp
          tmpfl = trim(datdir)//"/"//trim(arg)//"/"//sutfyles(k)
          open(unit=33,file=tmpfl,status='old')
          read(33,'(A)') arbig
          do while (arbig(1:1).eq."#") 
            read(33,'(A)') arbig
          enddo
          backspace(33)
          read(33,'(A)') arbig
          nopln = icountel(arbig)
          backspace(33)

          do while(kk.ge.0)
            read(33,*,end=301) (dtmp(j),j=1,nopln)
            do j=1,nopln
              kk = kk + 1
              protut(kk) = dtmp(j)
            enddo
          enddo
          close(33)
301       continue
        enddo

        if (jdf.eq.1) then
          tmpfl = trim(datdir)//"/"//trim(arg)//"/genn.dgn"
          open(unit=35,file=tmpfl,status='old')
        endif
        if (jdf.eq.1) then
          read(35,*,end=376) iprotdg(i)
        else
          iprotdg(i) = 1
        endif
        protid(i) = trim(arg)
376     continue
        close(35)
        i = i + 1
      enddo

      close(22)

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
      subroutine numpro(itmp)
      use data_files
      implicit none
      integer itmp, i
      
      open(unit=1,file=fylein,status='old')	!count number of database files 
      i = 0
      do while (i.ge.0)
        read(1,*,end=226)
        i = i + 1
      enddo
226   continue
      close(1)
      itmp = i-3
      end
!_______________________________________________________________________
!-- Initialize the neural network (the weights are randomized)
      subroutine netinit(nip,nw1,nw2,nop)
      use randseed
      use NN_wei
      use NN_owei
      implicit none
      integer i,j,k
      double precision atmp
      integer nip,nw1,nw2,nop
      
!***    randomizing weights
        do j=1,nip+1
          do k=1,nw1
            call RANDOM_NUMBER(atmp)
            atmp = atmp - 0.5
            wi(j,k) = atmp
            oldwi(j,k) = atmp / 10.d0
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
        do k=1,nop
          call RANDOM_NUMBER(atmp)
          atmp = atmp - 0.5
          wo(i,k) = atmp
          oldwo(i,k) = atmp / 10.d0
        enddo
      enddo

      end
!___________________________________________________________________
!-- Calculate the neural network's output
      subroutine calcnet(iwin,nip,nw1,nw2,nop,alpha,resv)
      use NN_wei
      use NN_nodes
      use DB_output
      implicit none
      integer i,j,k
      integer iwin,nip,nw1,nw2,nop
      double precision alpha, ztemp, resv(nip)
      double precision tfunc !transfer function

!** Calculate nueron values of the first hidden layer
      do i=1,nw1
        ztemp=0.d0
        do k=1,nip
	  ztemp = ztemp + resv(k) * wi(k,i)
	enddo
        ztemp = ztemp + wi(nip+1,i)	!bias
	hnu1(i) = tfunc(alpha*ztemp)
      enddo

      do i=1,nw2
        ztemp=0.d0
	do j=1,nw1
          ztemp = ztemp + hnu1(j) * wh(j,i)
	enddo
        ztemp = ztemp + wh(nw1+1,i)
	hnu2(i) = tfunc(alpha*ztemp)
      enddo

!** Calculate output nueron values
      do k=1,nop
        ztemp=0.d0
        do j=1,nw2
          ztemp = ztemp + hnu2(j) * wo(j,k)
        enddo
        ztemp = ztemp + wo(nw2+1,k)
        phi2p(k) = tfunc(alpha*ztemp)
        phiut((iwin-1)*nop+k) = phi2p(k)
      enddo
      
      end
!___________________________________________________________________
!-- This is the tanh transfer function.  f(x) = 1/(1+exp(alpha*x)) is another choice
!-- Note: the derivative f'(x) = alpha*f(x)*(1-f(x)) should be used then in calcnet/f
      double precision function tfunc(xoutp)
      implicit none
      double precision xoutp

      tfunc = tanh(xoutp)

      return
      end function tfunc
!___________________________________________________________________
!-- Calculate error over DB
      subroutine calcierr(trer,ofer,icond,cormax)
      use all_sizes
      use NN_input
      use data_files
      use DB_output
      implicit none
      integer i,j,ipr
      double precision trer,ofer,cormax
      integer icond

        trer = 0.d0; ofer = 0.d0
        do i=1,ntrap
          ipr = lstt(i)
          ipr = nop * (ipr - 1)
          do j=1+ipr,ipr+nop
            trer = trer + dabs(protut(j)-phiut(j))
          enddo
        enddo
        trer = trer / (nop * ntrap)
        do i=ntrap+1,ntotp
          ipr = lstt(i)
          ipr = nop * (ipr - 1)
          do j=1+ipr,ipr+nop
            ofer = ofer + dabs(protut(j)-phiut(j))
          enddo
        enddo
        ofer = ofer / (nop * (ntotp - ntrap))
        if (ofer.lt.cormax) then
          icond = 1
        else
          icond = 0
        endif

      end
!___________________________________________________________________
!-- Back propogation
      subroutine backprop(nip,nw1,nw2,nop,alpha,ratel,rmoment,resv,phi2)
      use NN_wei
      use NN_owei
      use NN_nodes
      implicit none
      integer j,k
      integer nip,nw1,nw2,nop
      double precision alpha,ratel,rmoment, errotmp, ztemp, resv(nip), phi2(nop)
      double precision errh1(nw1),errh2(nw2),erro(nop)

!** Calculate errors at the output/hidden2 layers
      do j=1,nw2
	ztemp=0.d0
	do k=1,nop
          errotmp = phi2(k) - phi2p(k)
	  erro(k) = errotmp * alpha * (1.d0-phi2p(k)*phi2p(k))
	  ztemp = ztemp + wo(j,k) * erro(k)
	enddo
        errh2(j) = ztemp * alpha * (1 - hnu2(j) * hnu2(j))
      enddo
      
!** Calculate the errors at the first hidden layer
      do k=1,nw1
        ztemp=0.d0
        do j=1,nw2
	  ztemp = ztemp + wh(k,j)*errh2(j)
        enddo
        errh1(k) = ztemp * alpha  * (1 - hnu1(k) * hnu1(k))

!** Update weights:
!** Update input/hidden layer weights (wi)
        do j=1,nip
          ztemp = ratel*errh1(k)*resv(j) + rmoment*oldwi(j,k)
          wi(j,k) = wi(j,k) + ztemp
          oldwi(j,k) = ztemp
        enddo
!** Update input bias weights (wi(nrx+1,1,1-nw1))
        ztemp = ratel*errh1(k) + rmoment*oldwi(nip+1,k)
        wi(nip+1,k) = wi(nip+1,k) + ztemp
        oldwi(nip+1,k) = ztemp
      enddo

!** Update hidden layer weights (wh)
      do k=1,nw2
        do j=1,nw1
          ztemp = ratel*errh2(k)*hnu1(j) + rmoment*oldwh(j,k)
          wh(j,k) = wh(j,k) + ztemp
          oldwh(j,k) = ztemp
        enddo
!** Update hidden bias weights (wh(nw1+1,1-nw2))
        ztemp = ratel*errh2(k) + rmoment*oldwh(nw1+1,k)
        wh(nw1+1,k) = wh(nw1+1,k) + ztemp
        oldwh(nw1+1,k) = ztemp

!** Update hidden/output layer weights (wo)
        do j=1,nop
          ztemp = ratel*erro(j)*hnu2(k) + rmoment*oldwo(k,j)
          wo(k,j) = wo(k,j) + ztemp
          oldwo(k,j) = ztemp
        enddo
      enddo
!** Update hidden/output layer bias weights (wo)
      k=nw2+1
      do j=1,nop
	ztemp = ratel*erro(j) + rmoment*oldwo(k,j)
	wo(k,j) = wo(k,j) + ztemp
	oldwo(k,j) = ztemp
      enddo

      end
!___________________________________________________________________
!-- save weights
      subroutine savewei()
      use NN_wei
      use NN_wei_sv
      implicit none

      wisv = wi
      whsv = wh
      wosv = wo

      end
!___________________________________________________________________
!-- Write weights
      subroutine writewei(tmpfl2,nip,nw1,nw2,nop,alpha)
      use data_files
      use NN_wei_sv 
      implicit none
      integer i,j,iflo
      integer nip,nw1,nw2,nop
      double precision alpha
      character tmpfl2*3600

      iflo = 48
      open(unit=iflo,file=tmpfl2,status='new')
      write(iflo,'(a,a,a)') '"',trim(datdir),'"'
      write(iflo,'(100(a,1x))') (trim(sinfyles(j)),j=1,nsinp)
      write(iflo,'(100(a,1x))') (trim(sutfyles(j)),j=1,nsutp)
      write(iflo,*) nip,nw1,nw2,nop, alpha
      do i=1,nip+1
        do j=1,nw1
          write(iflo,*) wisv(i,j)
        enddo
      enddo
      do i=1,nw1+1
        do j=1,nw2
	  write(iflo,*) whsv(i,j)
        enddo
      enddo
      do i=1,nw2+1
        do j=1,nop
          write(iflo,*) wosv(i,j)
        enddo
      enddo
      close(iflo)

      end
!___________________________________________________________________
!-- Read weights
      subroutine readwei(tmpfl,nip,nw1,nw2,nop,alpha)
      use NN_wei
      use NN_owei
      implicit none
      integer i,j,iflo
      integer nip,nw1,nw2,nop,nipc,nopc
      double precision alpha
      character tmpfl*3600

      iflo = 48
      open(unit=iflo,file=tmpfl,status='old')
      read(iflo,*)
      read(iflo,*)
      read(iflo,*)
      read(iflo,*) nipc,nw1,nw2,nopc,alpha
      if ((nipc.ne.nip).or.(nopc.ne.nop)) then
        write(0,*) "Different number of inputs/output b/w data and weight files, aborting.",nip,nipc,nop,nopc
        stop
      endif

      do i=1,nip+1
        do j=1,nw1
          read(iflo,*) wi(i,j)
	  oldwi(i,j) = wi(i,j) / 10.d0
        enddo
      enddo
      do i=1,nw1+1
        do j=1,nw2
	  read(iflo,*) wh(i,j)
	  oldwh(i,j) = wh(i,j) / 10.d0
        enddo
      enddo
      do i=1,nw2+1
        do j=1,nop
          read(iflo,*) wo(i,j)
	  oldwo(i,j) = wo(i,j) / 10.d0
        enddo
      enddo
      close(iflo)

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
        ta4 = exp(ta1) + exp(ta2) + exp(ta3) + krndm + iflg
        ita4log = nint(log10(ta4))
        ita4 = ta4
        dta4 = 10**ita4log * (ta4 - ita4)
        inttget = abs(dta4 - ta4) + ta3 + krndm
!        inttget = 64021
      endif
        

      return
      end
!___________________________________________________________________
!--
      subroutine permute(isz,ib,ie,iv)
      use randseed
      implicit none
      integer isz, i, j, ib, ie, itmp, iv(isz)
      double precision atmp
      
      do i=ib,ie
	call RANDOM_NUMBER(atmp)
        j = ib + nint((ie - ib) * atmp)
        itmp = iv(i)
        iv(i) = iv(j)
        iv(j) = itmp
      enddo

      end
!_______________________________________________________________________


