        Program GLOC
        !-------------------------HISTORY-----------------------------
        ! As of 15 Mar 2023 this program is CGEM version 1.1.0.1
        ! Versioning is as follows
        ! 4th number, any change to comments etc.
        ! 3rd number, any change to user interface, e.g., new features
        ! 2nd number, any bug fix affecting results
        ! 1st number, changes to basic model design and assumptions
        !
        !
        ! Estimates time to G-LOC from G0 based on
        ! calculated blood flow through the brain
        ! v 0.0.0.0 Created by Kyle Copeland, FAA, CAMI AAM-610 (Radiobiolgy Research), OCT-NOV 2012
        ! Modified by Kyle Copeland, FAA, CAMI AAM-631 (Numerical   
        ! Sciences Research) OCT 2109 - JAN 2020
        !   October-Dec 2019
        !      Corrected HMAP
        !      Added vision effects: lower limit of vision loss and 
        !         blackout
        !      Added ability to allow user-input dGzdt profile to  
        !         test for symptoms using for any flight segment 
        !         acceleration profile against 6 standard pilots or 
        !         a user-defined pilot
        !      Added PBG option
        !      Revised AGSM to be based on pressure
        !      Revised g-suit to be based on pressure and coverage
        ! v 1 (publication)  January-April 2020
        !      Revised output headers to match CGEM User Guide
        !      appendices
        !      Fix Lifebank tracking error in sub bankbal
        ! Modified by Kyle Copeland FAA, CAMI AAM-631 (Health Safety
        ! Information) FEB 2023
        ! v 1.0.0.1  February 2023 
        !      Revised comments for DOT publishing of source code
        ! v 1.1.0.1  March 2023
        !      Revised comments for DOT publishing of source code
        !      Converted sub eyebankbal to sub eyebbal, now based on 
        !      april 2020 sub bankbal (was from 2019 v of bankbal)
        !-------------------------NOTICES---------------------------            
        ! If this program is incorporated into other software, a  
        ! statement identifying it may be required under 17 U.S.C.  
        ! 403 to appear with any copyright notice.
        !
        ! This document is disseminated under the sponsorship of the
        ! U.S. Department of Transportation in the interest of
        ! information exchange. The United States Government assumes
        ! no liability for the contents thereof. 
        !
        !-------------------------PURPOSE---------------------------
        ! The CGEM code calculates occurrences of medical symptoms   
        ! (unconsciousness, greyout, and blackout) as well as symptom   
        ! abatements during user-defined exposures to accelerations. 
        ! A resource flow based model of cell metabolism is used
        ! along with modeling of changes in blood flow to important 
        ! regions. Outputs include report of cell performance and 
        ! symptoms likely to be observed in experimental participants
        ! due to the accelerations. The user controls all input from 
        ! file gloc_inp.dat, including other input and output file 
        ! names.      
        !                                                            
        !-------------------------AUTHOR-----------------------------                                                 
        ! Kyle Copeland, FAA, CAMI, Aerospace Medical Research 
        ! Division, AAM-600
        !                                                            
        !-----------------------REFERENCES---------------------------                                              
        ! Copeland, K. and J. E. Whinnery, Cerebral Blood Flow Based 
        !  Computer Modeling of Gz-Induced Effects, DOT/FAA/AM-23/6,  
        !  Office of Aerospace Medicine, Washington, DC, 2023.
        !------------------------------------------------------------

        IMPLICIT NONE
        CHARACTER(12)::outname, egpname, egpoutname
        CHARACTER(6)::sex
        CHARACTER(52)::comments
        Real(4):: Gnorm, Ghigh, fnorm, fmax, fcon, flife, gtm, beta, f0
        Real(4):: reserve, deltah, G0, gmax, flowmin, flowmax, flowlive
        Real(4):: HCON, HBD, smp, smh, gsh, gsp, seattilt, st, agsm, pgb
        Real(4):: bankcon, banklife, otherstrain, agtpcorr, sbc, tenlim
        Integer(4)::male,who,igfile,j,k,n1,n2,ne1,ne2,non1,non2
        Integer(4)::death,edeath,odeath 
        Real(4):: BSP,BDP,MSP,MDP,Drugdelay,hsuit,coverage,agsmp
        Real(4):: B,bcs,becs,bels,bloodO2,blmin,bls,boc,boec,boel,bol
        Real(4):: boncs,bonls,bonoc,cover,d1,dgdt,dgdtdown,dgdtup
        Real(4):: bonol,f,fog,fon,g,gc,Gf,gg,gmaxtime,gtol,GU
        Real(4):: hlap,howtall,hreduce,osp,pbg,pbgp,rp,sfp
        Real(4):: t,tblack,te,tgrey,th,tnt,ts,tsn,tu
        Real(4):: Gb,Geff,gz,Heye,HED,CBF,EBF,EGF
        Common /CARDIO/BSP,BDP,MSP,MDP,Drugdelay
        Common /SUBJECTS/Gnorm,Ghigh,fnorm,fmax,fcon,flife,gtm,beta,f0, &
     & reserve,flowmin,flowmax,flowlive,howtall,seattilt,male
        Common /ANTIG/agsmp,Hreduce,osp,smp,cover,sfp,pbgp,tenlim
        Common /banks/bankcon, banklife

! load subject and experiment data

        Open(unit=2,file='gloc_inp.dat', status='OLD')
        Print*, ' CGEM v.1.1.0.1, 15 MAR 2023 '
        Read(unit=2,FMT=10) Gnorm, comments
        Print*, Gnorm, comments
        Read(unit=2,FMT=10) Ghigh, comments
        Print*, Ghigh, comments
        Read(unit=2,FMT=10) fnorm, comments
        Print*, fnorm, comments
        Read(unit=2,FMT=10) fmax, comments
        Print*, fmax, comments
        Read(unit=2,FMT=10) fcon, comments
        Print*, fcon, comments
        Read(unit=2,FMT=10) flife, comments
        Print*, flife, comments
        Read(unit=2,FMT=10) gtm, comments
        Print*, gtm, comments
        Read(unit=2,FMT=10) beta, comments
        Print*, beta, comments
        Read(unit=2,FMT=10) BSP, comments
        Print*, BSP, comments
        Read(unit=2,FMT=10) BDP, comments
        Print*, BDP, comments
        Read(unit=2,FMT=10) MSP, comments
        Print*, MSP, comments
        Read(unit=2,FMT=10) MDP, comments
        Print*, MDP, comments
        Read(unit=2,FMT=10) bankcon, comments
        Print*, bankcon, comments
        Read(unit=2,FMT=10) banklife, comments
        Print*, banklife, comments
        Read(unit=2,FMT=10) gmaxtime, comments
        Print*, gmaxtime, comments
        Read(unit=2,FMT=10) dgdtup, comments
        Print*, dgdtup, comments
        Read(unit=2,FMT=10) dgdtdown, comments
        Print*, dgdtdown, comments
        Read(unit=2,FMT=15) male, comments
        Print*, male, comments
        Read(unit=2,FMT=10) howtall, comments
        Print*, howtall, comments
        Read(unit=2,FMT=10) smp, comments
        Print*, smp, comments
          IF (smp.gt.15.) smp = 15.
          IF (smp.LT.0.) smp = 0.
        Print*, 'Suit max presure revised to ', smp
        Read(unit=2,FMT=10) sbc, comments
        Print*, sbc, comments
        Read(unit=2,FMT=10) agsm, comments
        Print*, agsm, comments
        Read(unit=2,FMT=10) pbg, comments
        Print*, pbg, comments
        Read(unit=2,FMT=10) otherstrain, comments
        Print*, otherstrain, comments
        Read(unit=2,FMT=10) tenlim, comments
        Print*, tenlim, comments
        Read(unit=2,FMT=10) seattilt, comments
        Print*, seattilt, comments
        Read(unit=2,FMT=10) Drugdelay, comments
        Print*, Drugdelay, comments
        READ(Unit=2,FMT=20) outname
        Print*, 'Results to ', outname
        Read(unit=2,FMT=15) who , comments
        Print*, who, comments
        Read(unit=2,FMT=15) igfile, comments
        Print*, igfile, comments
        IF ((who.le.6).and.(who.ge.1)) then ! replace pilot dataset with a std set
           CALL Subject(who,fnorm,fmax,fcon,flife,beta,howtall,male)
        ENDIF
        if(igfile.eq.1) then !a custom defined experiment
           READ(Unit=2,FMT=20) egpname
           Print*, 'Experimental test Gz Profile from ' , egpname
           READ(Unit=2,FMT=20) egpoutname
           Print*, 'Experimental results to ', egpoutname
        endif
!reads
10      FORMAT (F10.4,A52)
15      FORMAT (I2,A52)
20      FORMAT (A12)
!writes
25      FORMAT (10F8.3)
        Rewind 2
           if (male.eq.1) then
              sex='male  '
           else
              sex='female'
           endif

        CALL AGCME(smp,sbc,agsm,pbg,otherstrain,sfp,pbgp,cover,agsmp)

        IF (sbc.gt.0.) then ! a suit is present
           Hreduce=hsuit(smp) !suit inflation raises heart, reducing heart-brain distance
        Else
           Hreduce=0.0 !no suit
        ENDIF
        Open(unit=1,file=outname,status='UNKNOWN')

        IF (igfile.eq.1) then
           CALL custom (egpoutname,egpname)
           STOP
        ENDIF
! dynamic variables
        GU=1.0
! parameters for run
!        drudelay=0.0
        st=seattilt ! angle from vertical of seat (can reduce brain alt above heart)
!        gs=gsuit ! added g tolerance from gear worn by subject
        gtol=gtm ! gtolerance multiplier, 1 for normal people
        G0=Gnorm ! ~1.0 g,  min Gz for experiment
        gmax=Ghigh ! e.g. 9.4 g, max Gz for some experiments
        flowmin=fcon ! 17-20 dl/min, flow needed to maintain consiousness indefinitely
        flowmax=fmax ! 110.0 dl/min, max physiological blood flow through brain
        flowlive=flife ! 9.0 dl/min, minimum flow needed to indefinitely sustain life
        f0=fnorm      !~52.0 dl/min, normal resting flow through brain
!       gtol = tolerance multiplier for experimental population vs normal people
!       reserve=5.0 seconds of consiousness if blood flow is stopped
        B=beta*1000. ! 2-3000., heart rate response time constant in ms
        rp=0.3  ! dgdt of maxed out physiological response
        bcs=bankcon
        bls=banklife
        becs=5.04
        bels=180.0
        boncs=5.04
        bonls=180.0

        HCON = HBD(howtall,male,st)
        HEYE = HED(howtall,male,st)
        write(1,*) ' CGEM v.1.1.0.1, 15 MAR 2023 '
        write(1,*) 'Equipment Related Parameters:'
        write(1,*) ' Seat tilt', st,' degrees'!KC 20200402 added to drugdelay
        write(1,*) ' Suit coverage fraction', sbc, '(0.0-0.7) '
        write(1,*) ' Suit max pressure', smp, 'psi, (0.0-15.0) '
        write(1,*) ' AGSM effectiveness',agsm, '(0-1) '
        write(1,*) ' Pressure breathing gear max effect',pbg, 'mmHg'
        Write(1,*) 'Subject Physiology Parameters'
        write(1,*) ' Muscle tensing HLAP at start',otherstrain,' mmHg'
        write(1,*) ' Tensing HLAP limit ',tenlim, ' mmHg'
        write(1,*) ' Pharma heart response delay ',Drugdelay,' s'
        write(1,*) ' Brain-heart dist, G0, min flow, no-flow awake(s)'
        write(1,*) HCON, G0, flowmin, bankcon
        write(1,*) ' max G,         F0,    max flow, Eye-heart dist'
        write(1,*) gmax, f0, flowmax, HEYE
        write(1,*) ' Subject height, sex, time@Gmax, no-flow alive(s)'
        write(1,*) howtall, male, sex, gmaxtime, banklife
        write(1,*) ' Units for distances are cm'
        write(1,*) ' Units for flows are dl/min'
        write(1,*) ' Units for acceleration are G'
        write(1,*) ' Units for times are seconds'
        write(1,*) ' Flow needed to avoid eventual death', flowlive,' s'
        write(1,*) ' Heart response time constant (ms):', B
        Write(1,*) 'Experimental Information '
        write(1,*) ' Rampdown rate for thus run:', dgdtdown
        write(1,*) ' T-GLOC T-RECON T-GREY  T-BLACK  dGz/dt   Gz@C      &
     &  Gz@U    Gz@B    Gz@G  minLife'
        Print*, 'Parameters:'
        write(*,*) ' Max suit pressure', smp,' psi'!KC typo fixed 20200402
        write(*,*) ' Seat tilt', st,' degrees'
        write(*,*) ' Suit coverage fraction', sbc, '(0.0-0.7) '
        write(*,*) ' Suit max pressure', smp, 'psi, (0.0-15.0) '
        write(*,*) ' AGSM effectiveness',agsm, '(0-1) '
        write(*,*) ' Pressure breathing gear max effect',pbg, 'mmHg'
        write(*,*) ' Muscle tensing HLAP at start',otherstrain,' mmHg'
        write(*,*) ' Tensing HLAP limit ',tenlim, ' mmHg'
        write(*,*) ' Pharma heart response delay ',Drugdelay,' s'
        Print*, ' Brain-heart dist, G0, min flow, no-flow awake(s)'
        write(*,*) HCON, G0, flowmin, bankcon
        write(*,*) ' max G,         F0,    max flow, Eye-heart dist'
        write(*,*) gmax, f0, flowmax, HEYE
        Print*, ' Subject height, sex, time@Gmax, no-flow alive(s)'
        Print*, howtall, male,sex,gmaxtime, banklife
        Print*, ' Units for distances are cm'
        Print*, ' Units for flows are dl/min'
        Print*, ' Units for acceleration are G'
        Print*, ' Units for times are seconds'
        Print*, ' Flow needed to avoid eventual death', flowlive,' s'
        Print*, ' Heart response time constant (ms):', B
        Print*, ' Rampdown rate for thus run:', dgdtdown
        write(*,*) ' T-GLOC T-RECON T-GREY  T-BLACK  dGz/dt   Gz@C      &
     &  Gz@U    Gz@B    Gz@G  minLife'
! std conciousness to uncon and back experiments
        If (dgdtup.ne.0.0 .and. dgdtdown.ne.0.0) then
            goto 300 !skip loops, do it just once
        endif
        do k = 1,1!3 !KC 20200402 from 1,3 to 1,1 !adjust grid for finer resolution
         do j = 1,1000!9 !KC 20200402 from 1,9 to 1,1000
          dgdt = j*10.**k/1000. !onset rate for Gz, Gz/sec
          n1=0 !starting state is conscious 0
          n2=0 ! unconcious = 1
               ! returned to conscious 2
          boc = bankcon
          bol = banklife
          ne1=0 !starting state is good vision 0
          ne2=0 ! greyout = 1
                ! returned to full vision 2
          boec = 5.04
          boel = banklife
          non1=0 !starting state is optic nerves work 0
          non2=0 ! optic nerves shut down = 1
                 ! returned to function 2
          bonoc = 5.04
          bonol = banklife
          t=0.0
          te=0.0
          tu=0.0
          th=0.0
          tgrey=0.
          tblack=0.
          Gg=0.0
          Gb=0.0
          Gc=0.0
          tsn=0.0
          tnt=0.0
!          Print *, k,j,dgdt !diagnostic
          do while (n1.eq.n2) !loops until subject is unconscious
           ! instantaneous value of Gz
           G=GZ(t,G0,dgdt,gmax,GU)
           Gf=Geff(G,G0) !effective Gz
! effective cerebral blood flow adjusted for blood oxygen content
! gsuit helps flow, but no effect on lung compression
! (the most modern suits incorporate a variable pressure cpap (PBG) to help with
!  lung compression AGSM straining and alleviate lung compression)
           CALL HLBP(t,B,G,d1,HLAP,tnt,tsn)
           F=CBF(fnorm,G0,HCON,Gf,HLAP)*BloodO2(G,st)
           FON=EBF(fnorm,G0,HEYE,Gf,HLAP)*BloodO2(G,st)
           FOG=EGF(fnorm,G0,HEYE,Gf,HLAP)*BloodO2(G,st)
! metabolic results of flow
           call bankbal(F,flowmin,flowlive,boc,bol,n1,n2,bcs,bls,death)
           call eyebbal(FON,flowmin,flowlive,bonoc,bonol,non1,non2,     &
     &                    boncs,bonls,odeath)
           call eyebbal(FOG,flowmin,flowlive,boec,boel,ne1,ne2,         &
     &                    becs,bels,edeath)
           if (non1.ne.non2) then
               Gb=G
               tblack=t/1000. !blackout in seconds
           endif
           if (ne1.ne.ne2) then
               Gg=G
               tgrey=t/1000. !greyout in seconds
           endif
           if (n1.eq.n2) then
              t=t+1.
           else
               GU=G
               ts=t/1000. !LOCINTI
!              write(*,25) ts*gtol, dgdt, GLOC, F, GCON, bol
!              write(1,25) ts*gtol, dgdt, GLOC, F, GCON, bol
           endif
          end do
! Point of GLOC for experiment j,k
! End of time-to-unconsciousness loop

100       do while (t.lt.((gmaxtime+ts)*1000.))
! Begin hold for gmaxtime seconds
! Must advance the bank values during hold down time
!         No need to recalculate F, boc, bol
!           if ( t.eq.0 .or. (t.gt.0.0 .and. MOD(t,1000.).eq.0)) then
!              print*, t, G, F, boc, bol
!           endif
! metabolic results of flow needed to check for death
!           CALL HLBP(t,B,G,d1,HLAP,tnt,tsn)
!           F=CBF(fnorm,G0,HCON,Gf,HLAP)*BloodO2(G,st)
!           FON=EBF(fnorm,G0,HEYE,Gf,HLAP)*BloodO2(G,st)
!           FOG=EGF(fnorm,G0,HEYE,Gf,HLAP)*BloodO2(G,st)
           call bankbal(F,flowmin,flowlive,boc,bol,n1,n2,bcs,bls,death)
           call eyebbal(FON,flowmin,flowlive,bonoc,bonol,non1,non2,     &
     &                    boncs,bonls,odeath)
           call eyebbal(FOG,flowmin,flowlive,boec,boel,ne1,ne2,         &
     &                    becs,bels,edeath)
           if (non1.ne.non2) then
               Gb=G
               tblack=t/1000. !LOCT or RTC
           endif
           if (ne1.ne.ne2) then
               Gg=G
               tgrey=t/1000. !LOCT or RTC
           endif

             IF (death.eq.1) then
                Print*, 'Subject death has occurred!', t, bol, death
!                Write(1,110) G
                t=(gmaxtime+ts)*1000.
110             FORMAT ('For G=',F9.3,' Subject death has occurred!')
             endif
           t=t+1.
          end do
          blmin=bol
          th=t ! time at end of hold
! End of time held unconscious at high Gz
200    CONTINUE !Begin return to conciousness by lowering Gz
          IF (dgdtdown.eq.0.) then
            dgdt = -dgdt ! assume symmetry with onset rate for Gz, Gz/sec
          ELSE
            dgdt=-1.*dgdtdown
          ENDIF
          n1=1 ! starting state is conscious if n1= 0
          n2=1 ! unconscious if n1 = 1
               ! subject conscious when n2=2
          do while (n2.eq.1)
            t=t+1.0 ! time in ms
! instantaneous value of Gz
            G=GZ(t-th,G0,dgdt,gmax,GU)
            Gf=Geff(G,G0) !effective Gz
! effective cerebral blood flow adjusted for blood oxygen content
! gsuit helps flow, but no effect on lung compression
           CALL HLBP(t-th,B,G,d1,HLAP,tnt,tsn)
           F=CBF(fnorm,G0,HCON,Gf,HLAP)*BloodO2(G,st)
           FON=EBF(fnorm,G0,HEYE,Gf,HLAP)*BloodO2(G,st)
           FOG=EGF(fnorm,G0,HEYE,Gf,HLAP)*BloodO2(G,st)
! metabolic results of flow
           call bankbal(F,flowmin,flowlive,boc,bol,n1,n2,bcs,bls,death)
           call eyebbal(FON,flowmin,flowlive,bonoc,bonol,non1,non2,     &
     &                    boncs,bonls,odeath)
           call eyebbal(FOG,flowmin,flowlive,boec,boel,ne1,ne2,         &
     &                    becs,bels,edeath)
           IF (bol.lt.blmin) then
              blmin=bol
           ENDIF
            if (n1.ne.n2) then
              te=t/1000. !total experimental time to return to con
              tu=te-ts-gmaxtime
              GC=G
!gtm affects only time to unconciousness
!        write(1,*) ' T-GLOC T-RECON T-GREY T-BLACK dGz/dt Gz@C Gz@U '
!     &  ,'  Gz@G  Gz@B  lifebankmin'
               write(*,25) ts*gtm, tu, tgrey, tblack, -dgdt, GC, GU, Gb, &
     &          Gg, blmin
               write(1,25) ts*gtm, tu, tgrey, tblack, -dgdt, GC, GU, Gb, &
     &          Gg, blmin
            endif
          end do
        end do
      end do
      STOP
300   CONTINUE
          dgdt = dgdtup !onset rate for Gz, Gz/sec
          n1=0 !starting state is conscious 0
          n2=0 ! unconcious = 1
               ! returned to conscious 2
          boc = bankcon
          bol = banklife
          ne1=0 !starting state is good vision 0
          ne2=0 ! greyout = 1
                ! returned to full vision 2
          boec = 5.04
          boel = banklife
          non1=0 !starting state is optic nerves work 0
          non2=0 ! optic nerves shut down = 1
                 ! returned to function 2
          bonoc = 5.04
          bonol = banklife
          t=0.0
!          Print *, k,j,dgdt !diagnostic
          do while (n1.eq.n2) !loops until subject is unconscious
           ! instantaneous value of Gz
           G=GZ(t,G0,dgdt,gmax,GU)
           Gf=Geff(G,G0) !effective Gz
! effective cerebral blood flow adjusted for blood oxygen content
! gsuit helps flow, but no effect on lung compression
           CALL HLBP(t,B,G,d1,HLAP,tnt,tsn)
           F=CBF(fnorm,G0,HCON,Gf,HLAP)*BloodO2(G,st)
           FON=EBF(fnorm,G0,HEYE,Gf,HLAP)*BloodO2(G,st)
           FOG=EGF(fnorm,G0,HEYE,Gf,HLAP)*BloodO2(G,st)
! metabolic results of flow
           call bankbal(F,flowmin,flowlive,boc,bol,n1,n2,               &
     &                    bcs,bls,death)
           call eyebbal(FON,flowmin,flowlive,bonoc,bonol,non1,non2,     &
     &                    boncs,bonls,odeath)
           call eyebbal(FOG,flowmin,flowlive,boec,boel,ne1,ne2,         &
     &                    becs,bels,edeath)
           if (non1.ne.non2) then
               Gb=G
               tblack=t/1000. !LOCT or RTC
           endif
           if (ne1.ne.ne2) then
               Gg=G
               tgrey=t/1000. !LOCT or RTC
           endif
           if (n1.eq.n2) then
              t=t+1.
           else
               GU=G
               ts=t/1000.
           endif
          end do
! Point of GLOC for experiment j,k
! End of time-to-unconsciousness loop
400       do while (t.lt.((gmaxtime+ts)*1000.))
! Begin hold for gmaxtime seconds
! Must advance the bank values during hold down time
!         No need to recalculate F, boc, bol
!            but do check for death
!           if ( t.eq.0 .or. (t.gt.0.0 .and. MOD(t,1000.).eq.0)) then
!              print*, t, G, F, boc, bol
!           endif
           call bankbal(F,flowmin,flowlive,boc,bol,n1,n2,               &
     &                    bcs,bls,death)
           call eyebbal(FON,flowmin,flowlive,bonoc,bonol,non1,non2,     &
     &                    boncs,bonls,odeath)
           call eyebbal(FOG,flowmin,flowlive,boec,boel,ne1,ne2,         &
     &                    becs,bels,edeath)
           if (non1.ne.non2) then
               Gb=G
               tblack=t/1000. !LOCT or RTC
           endif
           if (ne1.ne.ne2) then
               Gg=G
               tgrey=t/1000. !LOCT or RTC
           endif
             IF (death.eq.1) then
                Print*, 'Subject death has occurred!', t, bol, death
                t=(gmaxtime+ts)*1000.
             endif
               t=t+1.
          end do
          blmin=bol
          th=t ! time at end of hold
! End of time held unconscious at high Gz
500    CONTINUE !Begin return to conciousness by lowering Gz
          dgdt=-1*dgdtdown
          n1=1 ! starting state is conscious if n1= 0
          n2=1 ! unconcious if n1 = 1
               ! subject conscious when n2=2
          do while (n2.eq.1)
            t=t+1.0 ! time in ms
! instantaneous value of Gz
            G=GZ(t-th,G0,dgdt,gmax,GU)
            Gf=Geff(G,G0) !effective Gz
! effective cerebral blood flow adjusted for blood oxygen content
! BTW gsuit helps flow, but no effect on lung compression
           CALL HLBP(t-th,B,G,d1,HLAP,tnt,tsn)
           F=CBF(fnorm,G0,HCON,Gf,HLAP)*BloodO2(G,st)
           FON=EBF(fnorm,G0,HEYE,Gf,HLAP)*BloodO2(G,st)
           FOG=EGF(fnorm,G0,HEYE,Gf,HLAP)*BloodO2(G,st)
! metabolic results of flow
           call bankbal(F,flowmin,flowlive,boc,bol,n1,n2,               &
     &                    bcs,bls,death)
           call eyebbal(FON,flowmin,flowlive,bonoc,bonol,non1,non2,     &
     &                    boncs,bonls,odeath)
           call eyebbal(FOG,flowmin,flowlive,boec,boel,ne1,ne2,         &
     &                    becs,bels,edeath)
            IF (bol.lt.blmin) then
              blmin=bol
            ENDIF
            if (n1.ne.n2) then
              te=t/1000. !total experimental time to return to con
              tu=te-ts-gmaxtime !
              GC=G
!gtm affects only time to unconciousness
!        write(1,*) ' T-GLOC T-RECON T-GREY T-BLACK dGz/dt Gz@C Gz@U '
!     &  ,'  Gz@G  Gz@B  lifebankmin'
              write(*,25) ts*gtm, tu, tgrey, tblack, -dgdt, GC, GU, Gb, &
     &          Gg, blmin
              write(1,25) ts*gtm, tu, tgrey, tblack, -dgdt, GC, GU, Gb, &
     &          Gg, blmin
            endif
          end do
600   STOP
      end program GLOC

       SUBROUTINE bankbal(f, fc, fl, bc, bl, n1, n2, bcs, bls, d)
       ! What happens to bank accounts as flow changes?
       ! Does subject change state-of-consciousness?
       IMPLICIT none
       INTEGER(4)::n1,n2,d
       Real(4)::f,fc,fl,bc,bl,bcs,bls
       Real(4)::cfactor,dpc,dpl,extra,pca
       LOGICAL::Odd

       IF ((Odd(n2)).and.(n1.eq.0)) then
         n1=1 !flip starting state of conciousness
       ELSEIF ((.NOT.(Odd(n2))).and.(n1.eq.1)) then
         n1=0 !flip starting state of conciousness
       ELSE

!      Cfactor converts life points to con points drawn from flow
       cfactor = fl/(fc-fl)
       extra=0.0
!       cfactor = 0.0!

       endif
!      A. conciousness-to-unconciousness path
       d=0
       if (n1.eq.0) then
!        1. flow not yet reduce below threshold
         if (f.ge.fc.and.bc.eq.bcs.and.bl.eq.bls) then
            return ! no changes needed to anything
!        2. flow reduce below fc and temporarily returned by heart action
!           to >= fc so bank can be refilled some at this time
         elseif (f.ge.fc.and.bc.lt.bcs.and.bl.eq.bls) then
            bc = bc + dpc(f,fc,fl)
            if (bc.gt.bcs) then !Can refill bank, but not overfill
               bc=bcs
            endif
            return
!        3. flow reduce below fc and fl and temporarily returned by heart action
!           to >= fc so bank(s) can be refilled some at this time
         elseif (f.ge.fc.and.bl.lt.bls) then
            bl = bl + dpl(f,fl)
            !a.) !Can refill, but not overfill a bank
            if (bl.gt.bls) then
               extra = bl-bls
                bl=bls
            endif
            !b.) add any extra to bc since bl is full
            ! i.) normalize point values
               pca=extra*cfactor
            ! ii.) add to bc
               bc=bc+pca-.001
               if (bc.gt.bcs) then !Can refill bank, but not overfill
                  bc=bcs
               endif

            return
!        4. Flow is below threshold for c but not l,
!           i.e., f.lt.fc, f.gt.fl, bc.le.bcs, bl.eq.bls
         elseif ((f.lt.fc).and.(f.ge.fl).and.(bl.eq.bls)) then
            bc = bc + dpc(f,fc,fl)
            if (bc.le.0) then
             ! bank can't be overdrawn!
               bc = 0
               n2=n2+1 ! subject becomes unconscious, but life not in danger.
            endif
            return
!        5.Flow is below c threshold but not l, lifebank can replenish:
!          f.lt.fc, f.gt.fl, bl.lt.bls
         elseif ((f.lt.fc).and.(f.ge.fl).and.(bl.lt.bls)) then
!         this represents an opportunity to add to bank l at expense of bank c
          ! a) Try to refil bl first
            bl = bl + dpl(f,fl)
          ! i) !Can refill, but not overfill a bank
            if (bl.gt.bls) then
               extra = bl-bls
               bl=bls
            endif
          ! b) add any extra to bc since bl is full
            ! i.) normalize point values
               pca=extra*cfactor
            ! ii.) add to bc
               bc=bc+pca
          ! c) then add losses as if all flow was to bl
               bc = bc - .001
               if (bc.gt.bcs) then !Can refill bank, but not overfill
                  bc = bcs
               endif
               if (bc.le.0) then
               ! bank can't be overdrawn!
                  bc = 0
                  n2 = n2 + 1
               endif
               return
!        6.Flow is below both thresholds: f.lt.fl, thus withdraw from both banks
         else
             bc=bc+dpc(f,fc,fl)
             bl=bl+dpl(f,fl)
             if (bc.le.0.) then
               ! bank can't be overdrawn!
                bc = 0
               ! subject becomes unconscious
                n2 = n2 + 1
             endif
             if (bl.le.0.) then
              ! subject could be BRAIN DEAD !!!!
                d=1 !fixed to integer 20230310 by KC
             endif
             return
         endif
       endif
!
!      B. unconsciousness-to-consciousness path
       if (n1.eq.1) then ! confirm unconscious state
          if (f.lt.fl) then
!        1.Flow is below both thresholds: f.lt.fl,
!          thus, if possible, withdraw from both banks
             bc=bc+dpc(f,fc,fl)
             bl=bl+dpl(f,fl)
             if (bc.le.0.) then
               ! bank can't be overdrawn!
                bc = 0
               ! subject remains unconscious
             endif
             if (bl.le.0.) then
              ! subject could be BRAIN DEAD !!!!
                d=1
             endif
             return
          endif
!
!        2.Flow is below consciousness threshold: f.lt.fc,f.ge.fl,bl.lt.bls
!          thus withdraw from bc only, try to fill bl
          if (f.lt.fc.and.f.ge.fl.and.bl.lt.bls) then
             bl=bl+dpl(f,fl)
             ! no change in state
             if (bl.ge.bls) then
                bl=bls
             endif
             return
          endif
!
!        3.Flow is below consciousness threshold,
!          bl is full: f.lt.fc,f.ge.fl,bl.ge.bls
!          thus a stable state, g induced coma with no threat of death!
          if (f.lt.fc.and.f.ge.fl.and.bl.eq.bls) then
               ! no change in banks or states
               return
          endif
!
!        4.Flow is above all thresholds, bl not full: f.ge.fc,bl.lt.bls
!          thus try to fill bl, use any leftover to fill bc
          if (f.ge.fc.and.bl.lt.bls) then
          ! a) add as much as possible to bl
             bl=bl+dpl(f,fl)
!             Print*, 'bls, bl, dpl', bls, bl, dpl(f,fl)
             extra=0.0
             if (bl.ge.bls) then
                extra=bl-bls
                bl=bls
             endif
          ! b) then add any extra to bc since bl is full
!             bc=bc+dpc(f,fc,fl) !This line treats as if extra is maxed.
!             Print*, 'dpc, bc', dpc(f,fc,fl), bc
             bc=bc+extra*cfactor-.001
!             Print*, 'extra, cfactor, bc', extra, cfactor, bc

             if (bc.ge.bcs) then !Can refill bank, but not overfill
                  bc = bcs
                  n2 = n2 + 1 ! change from unc to con! Wake up when bank refilled!
             elseif (bc.le.0) then ! bank can't be overdrawn!
                bc = 0
             else
             ! no adjustments to bc and remain unconscious
             endif
             return
          endif
!        5.Flow is above all thresholds, bl is full: f.ge.fc,bl.eq.bls
!          thus try to fill bc
          if (f.ge.fc.and.bl.eq.bls) then
          ! Commit all extra to refill bc since bl is full
            bc=bc+dpc(f,fc,fl)
             if (bc.ge.bcs) then !Can refill bank, but not overfill
                  bc = bcs
                  n2 = n2 + 1 ! change from unc to con! Wake up when bank refilled!
             else
             ! no adjustments to bc and remain unconscious
             endif
             return
          endif
        ENDIF
        END SUBROUTINE bankbal

       SUBROUTINE eyebbal(f, fc, fl, bc, bl, n1, n2, bcs, bls, d)
       ! This is bankbal with an added fast recovery constant
       ! Vision effect recovery is fast 
       ! KC Notes on using bankbal for vision effects
       ! What happens to retinal cells/vision bank account as flow changes?
       ! Does subject experience greyout or permanent damage?
       ! Assumptions and data from T Whinnery and EM Forster, 2015,
       ! Neurologic state transitions in the eye and brain: kinetics of loss and
       ! recovery of vision and consciousness. Vis Neurosci. 32:E008.
       ! doi: 10.1017/S095252381500005X.
       ! ~2.74 seconds to recover when stress is removed for both blackout and greyout
       ! (E Harvey, Physiological effects of positive G forces, IAC#750)
       ! 5.04 seconds of reserves that maintain full function (average from Kydd
       ! and Ashley 1970, Physiologic responses to short duration Gz.
       ! Warminster, PA: Naval Air Development Center, 1970;
       ! Report No. NADC MR-7012.
       IMPLICIT none
       INTEGER(4)::n1,n2,d
       Real(4)::f,fc,fl,bc,bl,bcs,bls
       Real(4)::cfactor,dpc,dpl,extra,pca,refill
       LOGICAL::Odd

       refill = 5.04/2.74 !fast recovery rate from references

       IF ((Odd(n2)).and.(n1.eq.0)) then
         n1=1 !flip starting state of conciousness
       ELSEIF ((.NOT.(Odd(n2))).and.(n1.eq.1)) then
         n1=0 !flip starting state of conciousness
       ELSE

!      Cfactor converts life points to con points drawn from flow
       cfactor = fl/(fc-fl)
       extra=0.0
!       cfactor = 0.0!

       endif
!      A. conciousness-to-unconciousness path
       d=0
       if (n1.eq.0) then
!        1. flow not yet reduce below threshold
         if (f.ge.fc.and.bc.eq.bcs.and.bl.eq.bls) then
            return ! no changes needed to anything
!        2. flow reduce below fc and temporarily returned by heart action
!           to >= fc so bank can be refilled some at this time
         elseif (f.ge.fc.and.bc.lt.bcs.and.bl.eq.bls) then
            bc = bc + dpc(f,fc,fl)*refill
            if (bc.gt.bcs) then !Can refill bank, but not overfill
               bc=bcs
            endif
            return
!        3. flow reduce below fc and fl and temporarily returned by heart action
!           to >= fc so bank(s) can be refilled some at this time
         elseif (f.ge.fc.and.bl.lt.bls) then
            bl = bl + dpl(f,fl)
            !a.) !Can refill, but not overfill a bank
            if (bl.gt.bls) then
               extra = bl-bls
                bl=bls
            endif
            !b.) add any extra to bc since bl is full
            ! i.) normalize point values
               pca=extra*cfactor*refill
            ! ii.) add to bc
               bc=bc+pca-.001
               if (bc.gt.bcs) then !Can refill bank, but not overfill
                  bc=bcs
               endif

            return
!        4. Flow is below threshold for c but not l,
!           i.e., f.lt.fc, f.gt.fl, bc.le.bcs, bl.eq.bls
         elseif ((f.lt.fc).and.(f.ge.fl).and.(bl.eq.bls)) then
            bc = bc + dpc(f,fc,fl)
            if (bc.le.0) then
             ! bank can't be overdrawn!
               bc = 0
               n2=n2+1 ! subject becomes unconscious, but life not in danger.
            endif
            return
!        5.Flow is below c threshold but not l, lifebank can replenish:
!          f.lt.fc, f.gt.fl, bl.lt.bls
         elseif ((f.lt.fc).and.(f.ge.fl).and.(bl.lt.bls)) then
!         this represents an opportunity to add to bank l at expense of bank c
          ! a) Try to refil bl first
            bl = bl + dpl(f,fl)
          ! i) !Can refill, but not overfill a bank
            if (bl.gt.bls) then
               extra = bl-bls
               bl=bls
            endif
          ! b) add any extra to bc since bl is full
            ! i.) normalize point values
               pca=extra*cfactor*refill
            ! ii.) add to bc
               bc=bc+pca
          ! c) then add losses as if all flow was to bl
               bc = bc - .001
               if (bc.gt.bcs) then !Can refill bank, but not overfill
                  bc = bcs
               endif
               if (bc.le.0) then
               ! bank can't be overdrawn!
                  bc = 0
                  n2 = n2 + 1
               endif
               return
!        6.Flow is below both thresholds: f.lt.fl, thus withdraw from both banks
         else
             bc=bc+dpc(f,fc,fl)
             bl=bl+dpl(f,fl)
             if (bc.le.0.) then
               ! bank can't be overdrawn!
                bc = 0
               ! subject becomes unconscious
                n2 = n2 + 1
             endif
             if (bl.le.0.) then
              ! subject could be BRAIN DEAD !!!!
                d=1 !fixed to integer 20230310 by KC
             endif
             return
         endif
       endif
!
!      B. unconsciousness-to-consciousness path
       if (n1.eq.1) then ! confirm unconscious state
          if (f.lt.fl) then
!        1.Flow is below both thresholds: f.lt.fl,
!          thus, if possible, withdraw from both banks
             bc=bc+dpc(f,fc,fl)
             bl=bl+dpl(f,fl)
             if (bc.le.0.) then
               ! bank can't be overdrawn!
                bc = 0
               ! subject remains unconscious
             endif
             if (bl.le.0.) then
              ! subject could be BRAIN DEAD !!!!
                d=1
             endif
             return
          endif
!
!        2.Flow is below consciousness threshold: f.lt.fc,f.ge.fl,bl.lt.bls
!          thus withdraw from bc only, try to fill bl
          if (f.lt.fc.and.f.ge.fl.and.bl.lt.bls) then
             bl=bl+dpl(f,fl)
             ! no change in state
             if (bl.ge.bls) then
                bl=bls
             endif
             return
          endif
!
!        3.Flow is below consciousness threshold,
!          bl is full: f.lt.fc,f.ge.fl,bl.ge.bls
!          thus a stable state, g induced coma with no threat of death!
          if (f.lt.fc.and.f.ge.fl.and.bl.eq.bls) then
               ! no change in banks or states
               return
          endif
!
!        4.Flow is above all thresholds, bl not full: f.ge.fc,bl.lt.bls
!          thus try to fill bl, use any leftover to fill bc
          if (f.ge.fc.and.bl.lt.bls) then
          ! a) add as much as possible to bl
             bl=bl+dpl(f,fl)
!             Print*, 'bls, bl, dpl', bls, bl, dpl(f,fl)
             extra=0.0
             if (bl.ge.bls) then
                extra=bl-bls
                bl=bls
             endif
          ! b) then add any extra to bc since bl is full
!             bc=bc+dpc(f,fc,fl)*refill !This line treats as if extra is maxed.
!             Print*, 'dpc, bc', dpc(f,fc,fl), bc
             bc=bc+extra*refill*cfactor-.001
!             Print*, 'extra, cfactor, bc', extra, cfactor, bc

             if (bc.ge.bcs) then !Can refill bank, but not overfill
                  bc = bcs
                  n2 = n2 + 1 ! change from unc to con! Wake up when bank refilled!
             elseif (bc.le.0) then ! bank can't be overdrawn!
                bc = 0
             else
             ! no adjustments to bc and remain unconscious
             endif
             return
          endif
!        5.Flow is above all thresholds, bl is full: f.ge.fc,bl.eq.bls
!          thus try to fill bc
          if (f.ge.fc.and.bl.eq.bls) then
          ! Commit all extra to refill bc since bl is full
            bc=bc+dpc(f,fc,fl)*refill
             if (bc.ge.bcs) then !Can refill bank, but not overfill
                  bc = bcs
                  n2 = n2 + 1 ! change from unc to con! Wake up when bank refilled!
             else
             ! no adjustments to bc and remain unconscious
             endif
             return
          endif
        ENDIF
        END SUBROUTINE eyebbal

        FUNCTION dpl(f, fl)
         IMPLICIT NONE
         Real(4)::f,fl,plu,ple,dpl
           plu = .001 ! a withdrawl for current life operations
            ple = .001*(f)/(fl) !fraction of needed flow available
                  !to restore resources
            dpl = ple-plu ! delta pl, ie, cell life time lost or gained
        end FUNCTION dpl

        FUNCTION dpc(f, fc, fl)
        IMPLICIT NONE
        Real(4)::dpc,f,fc,fl,fcl,pcu,pce
            pcu = .001 ! a withdrawl for current con operations
            fcl=fc-fl
            if (f.ge.fl) then
              pce = .001*(f-fl)/(fcl) !fraction of needed flow available
                    !to restore conbank account
            else
              pce = 0.0
            endif
            dpc = pce-pcu ! delta pc, ie, con time lost or gained
        end FUNCTION dpc

        FUNCTION HBD(h,m,s)
        IMPLICIT NONE
        !      KC  Modified from head-heart to include head-neck, neck-heart segments
        !          20191029

        ! h = subject height
        ! m = female or male (0 or 1)
        ! s = seatback angle from vertical in degrees (e.g., 30 for F-16)

        REAL(4)::BRAIN,normaliz,dZl,dZu,h,HEART,ans,HBD,s,NECK
        integer(4)::m
        Real(4)::agtpcorr,Hreduce
        Common /ANTIG/agtpcorr,Hreduce

        ! Heart-brain delta Z distance
        ! Based on ICRP Pub 110 phantom data hight and center of mass distances
        ! Note: Recruiting site lists allowable range as 162.56 to 195.58 cm
        ! Note: for male, eye c.o.m. is lower than brain c.o.m. by 2.7 cm
        !        and it is well known that the brain center needed to maintain
        !        conciousness is below eye level. Thus I instead choose the effective
        !       brain height distance as EBH = (eye-c.o.m + base of brain)/2
        ! NOTE: This must be adjusted for seat incline that does not include the head
        !        pilot facing forwards = below base of brain to center of conciousness.
        ! NECK = distance that is vertical regardless of seat incline to keep,
        !       based on c.o.m. of the cortical cervical spine.
        ! Note: heart movement with gsuit pressurization is normalized to person height
         If (m.eq.0) then
        ! Female data from ICRP Pub 110
           BRAIN = (156.4 + (163.*316./348.))/2
           HEART = 128.2
           NECK = 149.0
           normaliz = h/163.0
         else !use male data
           BRAIN = 164.3 ! (166.9 + 161.7)/2.
           HEART = 134.1! heart blood c.o.m. height in std man
           NECK = 156.2
           normaliz = h/176.0
         endif
         dZu = BRAIN-NECK
         dZl = NECK-HEART
        ans = normaliz*(dZu+dZl*cos(s*2.*3.141592/360.))
        ans = ans-Hreduce*cos(s*2.*3.141592/360.) !adjust for suit push-up effect
        HBD=ans
        END FUNCTION HBD

        FUNCTION HED(h,m,s)
        IMPLICIT NONE
        ! h = subject height
        ! m = female or male (0 or 1)
        ! s = seatback angle from vertical in degrees (e.g., 30 for F-16)

        REAL(4)::EYE,normaliz,dZl,dZu,h,HEART,ans,HED,s,NECK
        integer(4)::m
        Real(4)::agtpcorr,Hreduce
        Common /ANTIG/agtpcorr,Hreduce

        ! Heart-eye delta Z distance
        ! Based on ICRP Pub 110 ct scan based phantom data height and c.o.mass dist
        ! NOTE: This must be adjusted for seat incline that does not include the head
        !        pilot facing forwards = below base of brain to center of conciousness.
        ! NECK = distance that is vertical regardless of seat incline to keep,
        !       based on c.o.m. of the cortical cervical spine.
        ! Note: heart movement with gsuit pressurization is normalized to person height
         If (m.eq.0) then
        ! Female data from ICRP Pub 110
           EYE = 156.4
           HEART = 128.2
           NECK = 149.0
           normaliz = h/163.0
         else !use male data
           EYE = 166.9
           HEART = 134.1! heart blood c.o.m. height in std man
           NECK = 156.2
           normaliz = h/176.0
         endif
         dZu = EYE-NECK
         dZl = NECK-HEART
        ans = normaliz*(dZu+dZl*cos(s*2.*3.141592/360.))
        ans = ans-Hreduce*cos(s*2.*3.141592/360.)!adjust for suit push-up effect
        HED=ans
        END FUNCTION HED

        ! G force at time t
        FUNCTION GZ(t,G0,dgdt,gmax,GU)
        IMPLICIT NONE
        ! dgdt in G/s, t in ms
        ! assume t=0 at start of ramping, whether up or down

        Real(4)::gmax,Z,GZ,G0,GU,t,dgdt

        if (dgdt.ge.0.) then !ramping up to gmax
           Z = G0+t*dgdt/1000.
           if (Z.ge.gmax) then
               Z=gmax
           endif
        else
           Z = GU+t*dgdt/1000.
           if (Z.lt.G0) then !ramping back down to G0
               Z=G0
           endif
        endif
        GZ = Z
        end function GZ

        Function Geff(G,g0) !effective g
        IMPLICIT NONE
        REAL(4)::G,g0,ge,Geff
           ge=G
           IF (ge.lt.0.) then
              ge=-ge !treat negative as positive for flows effects
           endif
           Geff=ge
        END FUNCTION Geff

        ! flow rate through brain
        Function CBF(fnorm,G0,H,G,HMAP)
        ! Fun fact, normal blood volume in cerebral matter is 100-130 ml
        IMPLICIT NONE
        Real(4)::H,GZ,G0,G,t,fnorm,dPdGz,dPdH,Rmin
        Real(4)::ICP,CPP,BMAP,HMAP,CF,VR,CVR,CBF,CFT
        COMMON /VR/VR
        ICP = 9.0 ! average intracranial pressure from
        dPdH=-0.7335  ! mmHg/cm, from relative densities of blood and Hg
                      ! H= Height of vital brain area above heart in cm.
                      ! A constant value of H neglects any spinal compression.
!        dPdH=-1.066  ! mmHg/cm, 32/30cm from Harvey
        dPdGz = H*dPdH !BMAP drop per Gz. ~20 for H =27.27
        BMAP=HMAP+dPdGz*G   ! Mean Arterial Pressure, Brain
        CPP = BMAP-ICP      ! Cerebral Perfusion Pressure
        if (CPP.LT.0.) then ! No negative (draining) pressures allowed
           CPP=0.           ! Equivalent to cells holding most of the reserve
        endif               ! in themselves.
                            ! Blood pressure <0 means no resupply possible.
        CFT=CPP/1.6 !trial CF with normal vascular resistance of 1.6
!        CF=CPP/CVR(CFT,G,fnorm,CPP)
        VR=CVR(CFT,G,fnorm,CPP)
        CF=CPP/VR
        if (CF.gt.fnorm) then
           CF=fnorm! !Autonomic goal is to maintain normal flows
           !  and concious effort to increase flows do no good
           !  above the normal flows (all banks should be full if f>fc)
           !  the flow is limited here to fnorm
           !  for accuracy during the unc-to-con return
        endif
        CBF=CF
        end function CBF

        ! flow rate through eye 
        Function EBF(fnorm,G0,H,G,HMAP)
        IMPLICIT NONE
        ! Fun fact, ocular (eye) pressure is higher by 2 mmHg than surrounding brain
        ! Optic nerve function is not subject to this extra pressure, but retina is!
        Real(4)::H,GZ,G0,G,t,fnorm,dPdGz,dPdH,VR
        Real(4)::IOP,OPP,BMAP,HMAP,COF,COFT,EBF
        COMMON /VR/VR
        IOP = 11.0 !minimal introcular pressure to overcome for any eye function.
        ! Blackout occurs when none of the light receptors can function 
        ! (expected at ~ 5 SECONDS AFTER IOP = ICP+2, where ICP = 9)
        ! (Tipton, 1983, AFAMRL-TH-83-047 p.19)
        dPdH=-0.7335  ! mmHg/cm, from relative densities of blood and Hg
                      ! H= Height of vital brain area above heart in cm.
                      ! A constant value of H neglects any spinal compression.
 !       dPdH=-1.066  ! mmHg/cm, 32/30cm from Harvey
        dPdGz = H*dPdH !BMAP drop per Gz. ~20 for H =27.27
        BMAP=HMAP+dPdGz*G   ! Mean Arterial Pressure, Brain
        OPP = BMAP-IOP      ! Ocular Perfusion Pressure
        if (OPP.LT.0.) then ! No negative (draining) pressures allowed
           OPP=0.           ! Equivalent to cells holding most of the reserve
        endif               ! in themselves.
                            ! Blood pressure <0 means no resupply possible.
        COF=OPP/VR
        if (COF.gt.fnorm) then
           COF=fnorm
        endif
        EBF=COF
        end function EBF

        ! flow rate through eye

        Function EGF(fnorm,G0,H,G,HMAP)
        IMPLICIT NONE 
        ! Same as EBF except H->H+1.6, KC 20230310
        ! Fun fact, ocular (eye) pressure is higher by 2 mmHg than surrounding brain
        ! Optic nerve function is not subject to this extra pressure, but retina is!
        Real(4)::H,GZ,G0,G,t,fnorm,dPdGz,dPdH,VR
        Real(4)::IOP,OPP,BMAP,HMAP,COF,EGF
        COMMON /VR/VR
        IOP = 22.0 ! introcular pressure that must be overcome for full eye function.
        ! less than this results in grey out and eventually black out, when none of
        ! the light receptors can function (B.O expected ~5 SECONDS AFTER IOP = ICP+2)
        ! (Tipton, 1983, AFAMRL-TH-83-047 p.19)
        dPdH=-0.7335  ! mmHg/cm, from relative densities of blood and Hg
                       ! H= Height of vital brain area above heart in cm.
                       ! A constant value of H neglects any spinal compression.
                       ! Added height 1.6 cm for top of eye
!        dPdH=-1.066  ! mmHg/cm, 32/30cm from Harvey

        dPdGz = (H+1.6)*dPdH !BMAP drop per Gz. ~20 for H =27.27 
                            !Add 1.6 for top of the eye relative to center
        BMAP=HMAP+dPdGz*G   ! Mean Arterial Pressure, Brain
        OPP = BMAP-IOP      ! Ocular Perfusion Pressure
        if (OPP.LT.0.) then ! No negative (draining) pressures allowed
           OPP=0.           ! Equivalent to cells holding most of the reserve
        endif               ! in themselves.
                            ! Blood pressure <0 means no resupply possible.
        COF=OPP/VR
        if (COF.gt.fnorm) then
           COF=fnorm
        endif
        EGF=COF
        end function EGF

        FUNCTION CVR(Flow,G,fnorm,P) !cerebral vascular resistance
        IMPLICIT NONE
        Real(4)::G,Flow,Rmin,Rnorm,Rtest,fnorm,P,C1,CVR
        Integer(4)::i

        Rmin=.4 ! Cerebral Vascular Resistance: a value of
                ! 0.4 assumes subject response is maxed out.
                ! Normal resting resistance = 1.6
                ! The time response to regulate flow rate is
                ! less than 1 sec, so no time dependence is
                ! needed here unless dGz/dt is well above 9 g.
        Rnorm=1.6

        C1 = 1.134 !if P=9 and Rnorm=.16, c1 corrects result to a normal
                   !brain flow of 49.6 from 56.25

        IF (G.LE.1.4) then ! normal conditions prevail, no percieved stresses
           CVR=Rnorm
           Return
        ELSE ! under stress, automatic attempt to correct flow to at least fnorm
           ! max correction to fnorm/Rmin
           DO i=160,40,-1
              Rtest=Real(i)/100.
              Flow= P/Rtest
!              Flow= P/(C1 * Rtest) !Pressure and resistance are not precise enough
!                   to justify C1 inclusion
              IF (Flow.ge.fnorm) then
                 exit
              Else
                 CONTINUE
              End IF
           END DO
        ENDIF
        CVR=Rtest
        END FUNCTION CVR

        Function BloodO2(G,st)
                IMPLICIT NONE
        Real(4)::G,st,egz,BloodO2,percent,pi
        ! Blood O2 is reduced by high Gz lost lung capacity.
        ! At 9 G this is about 50%
        ! It is about 10% at 5 Gz, and seems to be fairly linear

         pi=3.12415926
         egz = G*cos(st*pi/180.)
         If (egz.gt.13.5) then
         percent = 0.0
         elseif (egz.lt.4.5) then
         percent = 100.0
         else
         percent = 100.0*(13.5-egz)/9.0
         endif
         BloodO2 = percent/100.
        end function BloodO2
!
        SUBROUTINE HLBP(t,B,G,Delay1,HLAP,tnt,tsn)
        IMPLICIT NONE
        ! Maximum Heart modulated rise in MAP due to outside stress
        ! Exponential rise from F0+c1 to c2 in about 14 seconds
        ! (D.A. Tipton, 1983, AFAMRL-TR-83-047)
        ! C1 = minumum additional BP rise
        ! t = time from start in milliseconds
        ! tnt = local total negative Gz expsure time
        ! tsn = time stamp of last negative Gz exposure
        ! B = exponential rise constant
        ! C2 = exponential rise tweaker
        Real(4)::t,B,Bt,C1,C2,G,Delay1,Delay2,Delay3,tadj,tnt,tsn,thg
        Real(4)::BSP,BDP,MSP,MDP,Drugdelay,osp,smp,ospgain,altosp

        Real(4)::agtpcorr,Hreduce,pcorr,altpcorr,cover
        Real(4)::agsmp,pbgp,sva,pbgpg,sfp,sfpg,tenlim
        Real(4)::HMAP_max,HMAP_rest,HMAP_rise,HLAP,dpdg
        Common /ANTIG/agsmp,Hreduce,osp,smp,cover,sfp,pbgp,tenlim
        !agtpcorr = max added pressure at heart level
        !Hreduce = upwards movement of heart due to gsuit inflation
        Common /CARDIO/BSP,BDP,MSP,MDP,Drugdelay
        ! BSP = baseline systolic blood pressure
        ! BDP = baseline diastolic blood pressure
        ! MSP = max systolic blood pressure
        ! MDP = max diastolic blood pressure
        ! DrugDelay = possible BP increase delay from
        !             pharmacological sources

!       Normal max excited bp is 220/100, thus HMAP cannot exceed 160.
!       Individual will depend on aerobic fitness, etc.
!       Some drugs will change these values, particularly systolic pressure.
!
        HMAP_rest = (BSP+BDP)/2
        HMAP_max = (MSP+MDP)/2

        IF (G.lt.1.) then
           agtpcorr=0 !gear and AGSM no help with negative g
        ELSEIf (G.lt.9.and.G.ge.1) then
           pbgpg=pbgp*(G-1.)/8. ! variable breathing assist maxed at 9 G
           sfpg=sfp*(G-1.)/8. ! variable suit inflation maxed at 9 G
           IF (agsmp.ge.osp) osp = 0.0 !AGSM and other straining do not add
           ! since AGSM is more effective than gripping the seat.
           ! Use the greater of agsmp or sfpg
               sva=0.
           IF (agsmp.ge.sfpg) then
               sva = agsmp ! AGSM supercedes the g suit effectiveness
         ! at raising ITP (heart level change will still matter, of course)
           ELSE
         !(agsmp.lt.sfpg)
               sva = sfpg ! AGSM less than the g suit effectiveness
         ! at raising ITP (heart level change will still matter, of course)
           ENDIF
           !Add effects to get a final max correction to pressure
           agtpcorr = sva+pbgpg+osp !adding pbg pressure here assumes it is balanced,
           ! if pbg is unbalanced it cannot add to agsm.
        ELSE !(G.ge.9.) then
           pbgpg=pbgp ! variable breathing assist maxed at 9 G
           sfpg=sfp ! variable suit inflation maxed at 9 G
           IF (agsmp.ge.osp) osp = 0.0 !AGSM and other straining do not add
           ! since AGSM is more effective than gripping the seat
           ! Use the greater of agsmp or sfpg
               sva=0.
           IF (agsmp.ge.sfp) then
               sva = agsmp ! AGSM supercedes the g suit effectiveness
         ! at raising ITP (heart level change will still matter, of course)
           ELSE
         !(agsmp.lt.sfp)
               sva = sfp ! AGSM does not supercede the g suit effectiveness
         ! at raising ITP (heart level change will still matter, of course)
           ENDIF
           !Add effects to get a final max correction to pressure
           agtpcorr = sva+pbgpg+osp !adding pbg pressure here assumes it is balanced,
           ! if pbg is unbalanced it cannot add to agsm.
        endif
        pcorr=agtpcorr

! Next correction bit is from Buick.
        IF (smp.eq.0.) then ! blood can pool in lower extremities
           dpdg=-1.5 ! 3mmHg/g loss in HL Systolic pressure
        ELSE ! blood force up by suit, pressure tends to build slightly
           dpdg=2.5*cover ! 5mmHg/G gain in HL Systolic pressure from std. suit
        endif

        HMAP_rise = HMAP_max-HMAP_rest+pcorr+dpdg*G
!        Print*, t,B,G,delay1,HLAP
! Examples of possible heart actions for normal people
      ! initially very nervous/excited test subject
!        HMAP_rise=20.0 !Max rise in heart level Mean Arterial Pressure
!        HMAP_rest=140.0 !resting heart level Mean Arterial Pressure 180/100
      ! initially somewhat nervous/excited test subject
!        HMAP_rise=45.0 !Max rise in heart level Mean Arterial Pressure
!        HMAP_rest=115.0 !resting heart level Mean Arterial Pressure 140/90
      ! initially calm test subject
!        HMAP_rise=60.0 !Max rise in heart level Mean Arterial Pressure
!        HMAP_rest=100.0 !resting heart level Mean Arterial Pressure 120/80
! A swami who is always calm as a test subject
!        HMAP_rise=0.0 !Max rise in heart level Mean Arterial Pressure
!        HMAP_rest=100.0 !resting heart level Mean Arterial Pressure

! Heart rate does not start to rise in relaxed subjects until =>1.4Gz
           tadj = t-Delay1 ! start at baseline heart action

           thg=tadj !time above 1.4 Gz
!           thg=0 ! comment this line to include muscle stress buildup after high Gz begins

           ospgain = 50.*thg/30000. ! to mimic effect on bp of building muscle tension
           ! up to 50 mmHg even in previously relaxed subjects over ~30 seconds
           !(except swamis and zen masters, of course;-), who remain relaxed.
           IF (ospgain.gt.50.) ospgain=50.
           !else cannot exceed 60 mmHg, physical limit of simple muscle tensing
           altosp=osp+ospgain
           IF (altosp.gt.tenlim) altosp=tenlim !max allowed for this run
           ! pre-exposure build-up cannot decrease

           IF (pcorr.gt.altosp) altosp = 0.0
           !no need to adjust if subject tensing is already accounted for

           delay2 = Drugdelay*1000. ! Different for beta blockers like metprolol,
                     ! which will slow the inital reaction by a few seconds
                     ! Time scale here is ms
           tadj=tadj-delay2
           IF (tsn.ne.0.0) THEN
              delay3 = tnt
           ELSE
              delay3=0.
           ENDIF
           tadj=tadj-delay3
           IF (tadj.lt.0) then
              tadj=0
           endif
        C2 = 1.0

        Bt = 1/B ! B is in milliseconds
        C1= (1.0 - exp(-Bt * tadj))
!        If (C1.gt.1.)
!        Print*, C1, t, tadj
        HLAP = HMAP_rest + HMAP_rise * C1 + altosp
        !assume test subjects transition from relaxed state to maxed bp
      end subroutine HLBP

      SUBROUTINE custom(outname,egpname)
      ! test the dgdt profile in egpname for gloc and rtc
        IMPLICIT NONE
        CHARACTER(12)::outname, egpname
        CHARACTER(29)::lhp1,lhp2,lhp3
        CHARACTER(87)::lastheader
        Real(4):: Gnorm, Ghigh, fnorm, fmax, fcon, flife, gtm, beta, f0
        Real(4):: reserve, deltah, G0, gmax, flowmin, flowmax, flowlive
        Real(4):: HCON,HBD,HEYE,HED,HLAP,seattilt,st,tm,tm1,tnt
        Real(4):: bankcon,banklife,BSP,BDP,MSP,MDP,tsn,f,bonol
        Real(4):: b,bcs,becs,bels,boc,boec,boel,bol,boncs,bonls,bonoc
        Real(4):: bls,dgdt,fog,fon,G,Geff,Gf,GU,gtol,howtall,rp,t
        Real(4):: t1,BLOODO2,CBF,EBF,EGF,Drugdelay
        Real(4), Dimension(1000,1000)::profile
        Integer(4)::male,i,ii,leo,n1,n2,ne1,ne2,neold,nold,non1,non2
        Integer(4)::nonold,death,edeath,odeath,nolines,totalt,ttot
        Common /CARDIO/BSP,BDP,MSP,MDP,Drugdelay
        Common /SUBJECTS/Gnorm,Ghigh,fnorm,fmax,fcon,flife,gtm,beta,f0, &
     & reserve,flowmin,flowmax,flowlive,howtall,seattilt,male
        Common /banks/bankcon, banklife
! dynamic variables
        GU=1.0
! parameters for run
        st=seattilt ! angle from vertical of seat (can reduce brain dist
                    ! above heart and change lung capacity)
        gtol=gtm ! gtolerance multiplier, 1 for normal people
        G0=Gnorm ! ~1.0 g,  min Gz for experiment
        gmax=Ghigh ! e.g. 9.4 g, max Gz for some experiments
        flowmin=fcon ! 17-20 dl/min, flow needed to maintain consiousness indefinitely
        flowmax=fmax ! 110.0 dl/min, max physiological blood flow through brain
        flowlive=flife ! 9.0 dl/min, minimum flow needed to indefinitely sustain life
        f0=fnorm      !~52.0 dl/min, normal resting flow through brain
        B=beta*1000. ! 2-3000., heart rate response time constant in ms
        rp=0.3  ! dgdt of maxed out physiological response
        bcs=bankcon
        bls=banklife
        becs=5.04
        bels=180.0
        boncs=5.04
        bonls=180.0
        HCON = HBD(howtall,male,st)
        HEYE = HED(howtall,male,st)
        lhp1='        Time      G   Geff  '
        lhp2='C_Bank F_Con   F_Vis   F_BO '
        lhp3='  BO_Bank  HMAP    C   V   B'
        lastheader=lhp1//lhp2//lhp3
1211    FORMAT(A87)
        OPEN(unit=4,file=outname,status='unknown')
        write(4,*) ' CGEM v.1.1.0.1, 15 MAR 2023 '
        write(4,*) ' Subject Parameters:'
        write(4,*) ' Brain-heart dist, G0, min flow, no-flow awake(sec)'
        write(4,*) HCON, G0, flowmin, bankcon
        write(4,*) ' max G,         F0,    max flow, Eye-heart dist'
        write(4,*) gmax, f0, flowmax, HEYE
        write(4,*) ' Subject height, male, no-flow alive(sec)'
        write(4,*) howtall, male, banklife
        write(4,*) ' Units for distances are cm'
        write(4,*) ' Units for flows are dl/min'
        write(4,*) ' Units for acceleration are G'
        write(4,*) ' Units for times are seconds'
!        write(4,*) ' Added tolerance from G suit ', gs,' g' !obsolete code
        write(4,*) ' Flow needed to avoid eventual death', flowlive,' s'
        write(4,*) ' USING these BPs',BSP,BDP,MSP,MDP
        write(4,*) ' Heart response time constant (ms):', B
        write(4,*) ' -------------------- Results -------------------'
        write(4,1211) lastheader

        write(*,*) ' CGEM v.1.1.0.1, 15 MAR 2023 '
        Print*, ' Subject Parameters:'
        Print*, ' Brain-heart dist, G0, min flow, no-flow awake(sec)'
        Print*, HCON, G0, flowmin, bankcon
        write(*,*) ' max G,         F0,    max flow, Eye-heart dist'
        write(*,*) gmax, f0, flowmax, HEYE
        Print*, ' Subject height, male, no-flow alive(sec)'
        Print*, howtall, male, banklife
        Print*, ' Units for distances are cm'
        Print*, ' Units for flows are dl/min'
        Print*, ' Units for acceleration are G'
        Print*, ' Units for times are seconds'
!        write(*,*) ' Added tolerance from G suit ', gs,' g'
        Print*, ' Using these BPs',BSP,BDP,MSP,MDP
        Print*, ' Flow needed to avoid eventual death', flowlive,' s'
        Print*, ' Heart response time constant (ms):', B
        Print*, ' -------------------- Results -------------------'
        WRITE(*,1211) lastheader
!        Print*, ' totalt,i,G,boc,F,FOG,FON,bonoc,HLAP,n2,ne2,non2'

          n1=0 !starting state is conscious 0
          n2=0 ! unconscious = 1
               ! returned to conscious 2
          death=0 ! subject is not dead
          boc = bankcon
          bol = banklife
          ne1=0 !starting state is good vision 0
          ne2=0 ! greyout = 1
                ! returned to full vision 2
          boec = 5.04
          boel = 180.0
          non1=0 !starting state is optic nerves work 0
          non2=0 ! optic nerves shut down = 1
                 ! returned to function 2
          bonoc = 5.04
          bonol = 180.0
          t=0.0
          tnt=0.0
          tsn=0.0
          G=G0
          F=fnorm
          FON=fnorm
          FOG=fnorm
          totalt=0
        Open(unit=3,file=egpname, status='OLD')
        read(3,*) nolines
         Do ii=1, nolines
          read(3,*) dgdt, leo
          DO i= 0, leo-1
           totalt=totalt+1
           tm=Real(totalt)
           CALL HLBP(tm,B,G,tm1,HLAP,tnt,tsn)
           IF (mod(tm,1000.).eq.0.) then
           ttot=totalt/1000
           Write(*,700)ttot,G,Gf,boc,F,FOG,FON,bonoc,HLAP,n2,ne2,non2
           Write(4,700)ttot,G,Gf,boc,F,FOG,FON,bonoc,HLAP,n2,ne2,non2
           endif
!           Print*, t,Gf,tr
!           WRITE(*,701)FON,bonol,bonoc,flowmin,flowlive,boncs,bonls
!           read (3,*) dgdt
           G=G+dgdt*.001
           IF (G.le.0.) then
              ! negative Gz exposure induces a delayed response to
              ! an immediately following exposure to +Gz
              tnt=tnt+1. !total time of Gz exposure
              tsn=Real(totalt) !time of last negative Gz exposure
              IF (tsn.GT.5000.) tsn=5000. !max interference is 5 seconds
           endif
           IF (G.le.1.4) then
              t1=Real(totalt)
           endif
           Gf=Geff(G,G0) !effective Gz, corrected for gear
! effective cerebral blood flow adjusted for blood oxygen content
! gsuit helps flow, but does not stop lung compression
           t=real(i)
           tm= Real(totalt)
           tm1=t1
           CALL HLBP(tm,B,G,tm1,HLAP,tnt,tsn)
           F=CBF(fnorm,G0,HCON,Gf,HLAP)*BloodO2(G,st)
           FON=EBF(fnorm,G0,HEYE,Gf,HLAP)*BloodO2(G,st)
           FOG=EGF(fnorm,G0,HEYE,Gf,HLAP)*BloodO2(G,st)
! metabolic results of flow
           nold=n2
           neold=ne2
           nonold=non2
           call bankbal(F,flowmin,flowlive,boc,bol,n1,n2,bcs,bls,death)
           call eyebbal(FON,flowmin,flowlive,bonoc,bonol,non1,non2,     &
     &                    boncs,bonls,odeath)
           call eyebbal(FOG,flowmin,flowlive,boec,boel,ne1,ne2,         &
     &                    becs,bels,edeath)

           IF ((neold.ne.ne2).or.(nonold.ne.non2).or.(nold.ne.n2)) then
           Write(*,700)totalt,G,Gf,boc,F,FOG,FON,bonoc,HLAP,n2,ne2,non2
           Write(4,700)totalt,G,Gf,boc,F,FOG,FON,bonoc,HLAP,n2,ne2,non2
           endif
          end do
         end do
         close(3)
         close(4)
700   FORMAT(1i12,2f7.3,6f8.3,3i4)
! 701   FORMAT(3x,7f8.3)
      END SUBROUTINE

      SUBROUTINE Subject(I,fnorm,fmax,fcon,flife,beta,howtall,male)
      ! There are six standard subjects: best male, mediam male, worst male
      ! best female, med. female, worst female.
      ! Or a user defined subject from gloc_inp.dat
        Real(4):: fnorm, fmax, fcon, flife, gtm, beta, howtall
        Real(4):: bankcon, banklife
        Real(4):: BSP,BDP,MSP,MDP,Drugdelay
        Integer(4)::male,i

        Common /banks/bankcon,banklife
        Common /CARDIO/BSP,BDP,MSP,MDP,Drugdelay

       IF (I.eq.1) Then
      ! 1. best male

        fnorm = 54. !normal blood flow rate through brain in dl/min (45-54)
        fmax = fnorm/0.4  !maximum blood flow rate through brain (about 110)
        fcon = 18. !blood flow needed through brain to maintain consciousness (18-20)
        flife = 8. !blood flow needed through brain to maintain brain cell life (8-10)
        gtm = 1.0 ! G tolerance scale factor, usually = 1.0, but can be as great as 1.53 using breathing techniques
        beta = 2.0 !time constant in seconds for heart rate response function
 !      (This should be about 1/7 the time to ramp up to full response.)
        bankcon = 15 !seconds of consciousness if oxygen/blood flow stops (5-15).
        banklife = 180!seconds of life if oxygen/blood flow stops.
        male = 1 ! 0=female 1=male
        howtall = 162.5 ! in cm 162.5-195.6 cm
        BSP=130 ! male at top of allowable uncontrolled range
        BDP=90
        MSP=213 !95th percentile male from FRIEND Cohort Age 20-40
        MDP=98

      Elseif (I.eq.2) then
            ! 2. midrange male

        fnorm = 49.5 !normal blood flow rate through brain in dl/min (45-54)
        fmax = fnorm/0.4  !maximum blood flow rate through brain (about 110)
        fcon = 19. !blood flow needed through brain to maintain consciousness (18-20)
        flife = 9. !blood flow needed through brain to maintain brain cell life (8-10)
        gtm = 1.0 ! G tolerance scale factor, usually = 1.0, but can be as great as 1.53 using breathing techniques
        beta = 2.5 !time constant in seconds for heart rate response function
 !      (This should be about 1/7 the time to ramp up to full response.)
        bankcon = 7.1 !seconds of consciousness if oxygen/blood flow stops (5-15).
        banklife = 180!seconds of life if oxygen/blood flow stops.
        male = 1 ! 0=female 1=male
        howtall = 179.0 ! in cm 162.5-195.6 cm Age 20-40
        BSP=120
        BDP=80
        MSP=177 !50th percentile
        MDP=80

      Elseif (I.eq.3) then
      ! 3. worst male

        fnorm = 45. !normal blood flow rate through brain in dl/min (45-54)
        fmax = fnorm/0.4  !maximum blood flow rate through brain (about 110)
        fcon = 20. !blood flow needed through brain to maintain consciousness (18-20)
        flife = 10. !blood flow needed through brain to maintain brain cell life (8-10)
        gtm = 1.0 ! G tolerance scale factor, usually = 1.0, but can be as great as 1.53 using breathing techniques
        beta = 3.0 !time constant in seconds for heart rate response function
 !      (This should be about 1/7 the time to ramp up to full response.)
        bankcon = 5 !seconds of consciousness if oxygen/blood flow stops (5-15).
        banklife = 180!seconds of life if oxygen/blood flow stops.
        male = 1 ! 0=female 1=male
        howtall = 195.6 ! in cm 162.5-195.6 cm
        BSP=100
        BDP=60
        MSP=147 !5th percentile from FRIEND Cohort Age 20-40
        MDP=59

! There are no obvious physiological differences to declare below.
! It is based on being able to fit in the seat. There will be differences in GTM,
! in median height of the population, etc.

      Elseif (I.eq.4) then
      ! 1. best female

        fnorm = 54. !normal blood flow rate through brain in dl/min (45-54)
        fmax = fnorm/0.4  !maximum blood flow rate through brain (about 110)
        fcon = 18. !blood flow needed through brain to maintain consciousness (18-20)
        flife = 8. !blood flow needed through brain to maintain brain cell life (8-10)
        gtm = 1.0 ! G tolerance scale factor, usually = 1.0, but can be as great as 1.53 using breathing techniques
        beta = 2.0 !time constant in seconds for heart rate response function
 !      (This should be about 1/7 the time to ramp up to full response.)
        bankcon = 15 !seconds of consciousness if oxygen/blood flow stops (5-15).
        banklife = 180!seconds of life if oxygen/blood flow stops.
        male = 0 ! 0=female 1=male
        howtall = 162.5 ! in cm 162.5-195.6 cm
        BSP=130
        BDP=90
        MSP=187 !95th percentile from FRIEND Cohort Age 20-40
        MDP=93

      Elseif (I.eq.5) then
         ! 5. midrange female

        fnorm = 49.5 !normal blood flow rate through brain in dl/min (45-54)
        fmax = fnorm/0.4  !maximum blood flow rate through brain (about 110)
        fcon = 19. !blood flow needed through brain to maintain consciousness (18-20)
        flife = 9. !blood flow needed through brain to maintain brain cell life (8-10)
        gtm = 1.0 ! G tolerance scale factor, usually = 1.0, but can be as great as 1.53 using breathing techniques
        beta = 2.5 !time constant in seconds for heart rate response function
 !      (This should be about 1/7 the time to ramp up to full response.)
        bankcon = 7.1 !seconds of consciousness if oxygen/blood flow stops (5-15).
        bank2life = 180!seconds of life if oxygen/blood flow stops.
        male = 0 ! 0=female 1=male
        howtall = 179.0 ! in cm 162.5-195.6 cm
        BSP=120
        BDP=80
        MSP=157 !50th percentile from FRIEND Cohort Age 20-40
        MDP=76

      Else !if (I.eq.6) then
      ! 6. worst female

        fnorm = 45. !normal blood flow rate through brain in dl/min (45-54)
        fmax = fnorm/0.4  !maximum blood flow rate through brain (about 110)
        fcon = 20. !blood flow needed through brain to maintain consciousness (18-20)
        flife = 10. !blood flow needed through brain to maintain brain cell life (8-10)
        gtm = 1.0 ! G tolerance scale factor, usually = 1.0, but can be as great as 1.53 using breathing techniques
        beta = 3.0 !time constant in seconds for heart rate response function
 !      (This should be about 1/7 the time to ramp up to full response.)
        bankcon = 5 !seconds of consciousness if oxygen/blood flow stops (5-15).
        banklife = 180!seconds of life if oxygen/blood flow stops.
        male = 0 ! 0=female 1=male
        howtall = 195.6 ! in cm 162.5-195.6 cm
        BSP=100
        BDP=60
        MSP=131 !5th percentile from FRIEND Cohort Age 20-40
        MDP=60

      END IF
      END SUBROUTINE

      SUBROUTINE AGCME(smpsi,sbc,agsm,pbg,ostrain,sfp,pbgp,cover,agsmp)
      !anti-g countermeasure effectiveness
      Real(4), INTENT(IN)::smpsi,sbc,agsm,pbg,ostrain
      Real(4), INTENT(OUT)::pbgp,sfp,cover,agsmp

      Real(4)::Csmp,Cpbg,CagsmPiPa,Csbc
      Real(4)::agsmITPmax,osp
      !useful constants
      ! HMAP heart level Mean arterial pressure
      ! ITP intrathoracic pressure
      Csmp = 3 !mmHg Hlap per PSI of suit inflation for normal mil
      Csbc = 0.35 ! fraction of body coverage expected for Csmp=3
      cover = sbc/Csbc
      IF (cover.gt.2.) cover=2.0 !effect maxes out at ~6 mmHg
      IF (cover.LT.0.) cover=0.0 ! no negative inflations!
      sfp=smpsi*Csmp*cover

      CagsmPiPa = 0.75 ! conversion of deltaITP to deltaHMAP from agsm
      agsmITPmax = 130 ! maximally effective agsm itp
      agsmp = agsm * agsmITPmax * CagsmPiPa

      Cpbg = 1.0 !conversion of positive breathing presure to heart level MAP
      pbgp = Cpbg*pbg
      IF (pbgp.GT.60.) pbgp = 60.
      IF (pbgp.LT.0.) pbgp = 0.

      osp = ostrain

      END SUBROUTINE

      FUNCTION hsuit(psi)
      Real(4)::hsuit, psi
      ! heart level change from gsuit inflation pressing heart upwards
      hsuit = psi*0.6 !expect 6mm/psi
      END FUNCTION

      function Odd(n)
      logical::Odd
      integer(4)::n,test

      Odd = .FALSE.
      test = MOD(n,2) ! get remainder of division of n by 2

      IF (test.eq.1) then
         Odd = .TRUE.
      endif
      return
      end function Odd
