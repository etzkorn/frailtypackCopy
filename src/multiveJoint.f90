    
!-------------------------------------
!             Joint multive
!-------------------------------------

! Organizational Overview of:
! (0) Variable definitions
! (1) Unpack model controls
! (2) Parameter initialization (if specified by user)
! (3) Copy data
! (4) Sort Unique Dates for each event type
! (5) Construct Vector of Spline Nodes (Z_i)
! (6) Map dates to sequence of unique dates
! (7) Calculate penalties
! (8) Copy initialized parameters to "b"
! (9) Model Optimization
! (10) Organize Hessian Matrix
! (11) Evaluate hazard and survival estimates over a grid
! (12) Calculate AIC
! (13) Calculate Fitted Values
! (14) Deallocation and return

! Supporting Functions
! (a) vecspli
! (b) vecpenP
! (c) susp
! (d) cosp ! calcul les points pour les fonctions ! et leur bandes de confiance
! (e) conf1
! (f) conf
! (g) isp
! (h) mmsp
! (i) multi ! matrix multiply A by B with the result in C
! (j) func30 ! calcul de l integrant, pour un effet aleatoire donne frail et un groupe donne auxig (cf funcpa)

! istop codes:
!    1: convergence
!    2: maximum number of iterations
!    4: error (any calculation error)

!------------------------------------------------------------------------------
    subroutine joint_multiv(controls,nobsEvent,k0,&
    tt00,tt10,tt0meta0,tt1meta0,ic0,icmeta0, &
    groupe0,groupe0meta,groupe0dc,tt0dc0,tt1dc0,icdc0,icdc20,&
    nbvar,&
    vax0,vaxdc0,vaxmeta0,vaxdc20,&
    noVarEvent, &
    np,b,H_hessOut,HIHOut,resOut,LCV,critCV,mtEvent,mt1Event,x1Out,lamOut,xSu1,suOut,x2Out, & 
    lam2Out,xSu2,su2Out,x3Out,lam3Out,xSu3,su3Out,x4Out,lam4Out,xSu4,su4Out, &
    ni,cptEvent,ResMartingaleEvent,frailtyEstimates, &
    linearpred,linearpreddc,linearpredM,linearpreddc2,ziOut1,ziOutdc,ziOutmeta, &
    time,timedc,timeM,ghNodes0,ghWeights0,tolerance0)

    use parametersmultiv
    use residusMmultiv
    use comonmultiv
    use taillesmultiv
    use splines
    use optim
    use sortiemultive

    IMPLICIT NONE  
!-------------------------------------------------
! Inputs (in order of entry)
    integer,dimension(12),intent(in)::controls
    !controls = c(maxit = maxit[1], # [1] 
    !     initialize = initialize, # [2] 
    !     typeof = typeof, # [3]
    !     equidistant = 1, # [4]
    !     irep = !crossVal, # [5] irep
    !     recurrentAG = recurrentAG, # [6] ag0
    !     nbIntervEvent = 0, # [7] nbIntervEvent
    !     n.knots = n.knots, # [8]
    !     event2.ind = event2.ind, # [9]
    !     terminal2.ind = terminal2.ind, # [10]
    !     GHpoints = GHpoints, # [11]
    !     jointGeneral = as.integer(jointGeneral)) # [12] typeJoint0    
    integer,dimension(3),intent(in)::nobsEvent
    !    nobsEvent=c(nsujet0,ng0,nsujetmeta0)
    !        length/number of rows of data for recurrent1, terminal1, recurrent2
    double precision,dimension(4)::k0,kappaCV
    double precision,dimension(nobsEvent(1))::tt00,tt10
    double precision,dimension(nobsEvent(3))::tt0meta0,tt1meta0
    integer,dimension(nobsEvent(1))::ic0
    integer,dimension(nobsEvent(3))::icmeta0 
    integer,dimension(nobsEvent(1))::groupe0
    integer,dimension(nobsEvent(3))::groupe0meta 
    integer,dimension(nobsEvent(2))::groupe0dc
    double precision,dimension(nobsEvent(2))::tt0dc0,tt1dc0
    integer,dimension(nobsEvent(2))::icdc0,icdc20
    integer,dimension(4),intent(in)::nbvar
    !    nbvar=c(nva10,nva20,nva30,nva40)
    !        length/number of columns of data for recurrent1, terminal1, recurrent2
    double precision,dimension(nobsEvent(1),nbvar(1)),intent(in):: vax0
    double precision,dimension(nobsEvent(2),nbvar(2)),intent(in):: vaxdc0
    double precision,dimension(nobsEvent(3),nbvar(3)),intent(in):: vaxmeta0
    double precision,dimension(nobsEvent(2),nbvar(4)),intent(in):: vaxdc20
    integer,dimension(4),intent(in)::noVarEvent
    !    noVarEvent=c(noVar1,noVar2,noVar3,noVar4)
    !        indicator for no variables for a given event
    integer::np
    double precision,dimension(np)::b
    double precision,dimension(np,np)::H_hessOut,HIHOut
    double precision::resOut
    double precision,dimension(2)::LCV
    integer,dimension(6),intent(out)::critCV
    !    critCV=c(ier,istop,istopshared(1),istopshared(2),istopshared(3),istopshared(4))
    !        indicators for errors
    !         istop codes:
    !            1: convergence
    !            2: maximum number of iterations
    !            4: error (any calculation error)
    !        ier is often -1 when hessian matrix fails to invert
    
    integer,dimension(4),intent(in)::mtEvent
    !    mtEvent=c(mt1,mt2,mt3,mt4)
    integer,dimension(4),intent(in)::mt1Event
    !    mt1Event=c(mt11,mt12,mt13,mt14)

    double precision,dimension(mtEvent(1))::x1Out
    double precision,dimension(mtEvent(1),3)::lamOut
    double precision,dimension(mt1Event(1))::xSu1
    double precision,dimension(mt1Event(1),3)::suOut

    double precision,dimension(mtEvent(2))::x2Out
    double precision,dimension(mtEvent(2),3)::lam2Out
    double precision,dimension(mt1Event(2))::xSu2
    double precision,dimension(mt1Event(2),3)::su2Out

    double precision,dimension(mtEvent(3))::x3Out, x4Out
    double precision,dimension(mtEvent(3),3)::lam3Out, lam4Out
    double precision,dimension(mt1Event(3))::xSu3, xSu4
    double precision,dimension(mt1Event(3),3)::su3Out, su4Out

    integer::ni
    integer,dimension(4),intent(out)::cptEvent
    !    cptEvent=c(cpt,cpt_dc,cptmeta,cpt_dc2)
    !        counts of events of each type
    double precision,dimension(nobsEvent(2),3)::ResMartingaleEvent
    !    ResMartingaleEvent=c(Res_martingale,Res_martingaledc,Res_martingale2,Res_martingaledc2)
    double precision,dimension(nobsEvent(2),5)::frailtyEstimates
    !    frailtyEstimates=c(frailtypred,frailtypred2,frailtyvar,frailtyvar2,frailtyCorr)
    
    double precision,dimension(nobsEvent(1))::linearpred
    double precision,dimension(nobsEvent(2))::linearpreddc,linearpreddc2
    double precision,dimension(nobsEvent(3))::linearpredM

    double precision,dimension(controls(8)+6)::ziOut1,ziOutdc,ziOutmeta
    double precision,dimension(controls(7)+1)::time,timedc,timeM

    double precision,dimension(controls(11))::ghNodes0, ghWeights0
    double precision,dimension(3)::tolerance0

! End Inputs    
!-------------------------------------------------
! Other Variables
    integer::maxit0,maxit00,mt11,mt12,mt13,mt14,nsujetmeta0,ag0,initialize,&
    event2_ind0,terminal2_ind0
    integer::nbIntervEvent,nz0
    integer::nsujet0,ng0,nva10,nva20,nva30,nva40,mt1,mt2,mt3,mt4,irep
    integer::equidistant,equidistant0
    integer::ss,sss
    double precision,dimension(:),allocatable::b01,b02,b03,b04!,b02joint
    integer::np1,np2,np3,np4,np2joint
    double precision,dimension(4)::shape_weib,scale_weib    
    integer::noVar1,noVar2,noVar3,noVar4
    integer::cpt,cptmeta,cpt_dc,cpt_dc2,ier
    integer::groupe,groupemeta,ij,kk,j,k,n,ii,iii,iii2,cptstr1,cptstr2   &
    ,i,ic,icmeta,icdc,icdc2,istop,cptni,cptni1,cptni2,nb_echec,nb_echecor,id,cptbiais &
    ,cptauxdc   
    double precision::tt0,tt0meta,tt0dc,tt1,tt1meta,tt1dc,h,res,min,mindc,max,pord, &
    maxdc,maxt,maxtmeta,maxtdc,moy_peh0,moy_peh1,lrs,BIAIS_moy,mint,mintdc,mintmeta
    double precision,dimension(2)::res01
!AD: add for new marq
    double precision::ca,cb,dd
    double precision,external::funcpaMultivSplines,funcpaMultivCpm,funcpaMultivWeib
    !external::funcpaMultivSplines,funcpaMultivCpm,funcpaMultivWeib
!cpm
    integer::indd,ent,entdc,entmeta,typeof0,nzsha
    double precision::temp,kappa1    
    integer::nbinterv0
!predictor
    double precision,dimension(nobsEvent(2))::Res_martingale,Res_martingaledc,Res_martingale2,&
    frailtypred,frailtypred2,frailtyvar,frailtyvar2,frailtyCorr
    double precision,external::funcpamultires
    double precision,dimension(1,nbvar(1))::coefBeta
    double precision,dimension(1,nbvar(2))::coefBetadc    
    double precision,dimension(1,nbvar(3))::coefBetaM        
    !double precision,dimension(1,nbvar(4))::coefBetadc2        
    double precision,dimension(1,nobsEvent(1))::XBeta
    double precision,dimension(1,nobsEvent(2))::XBetadc!, XBetadc2
    double precision,dimension(1,nobsEvent(3))::XBetaM
    double precision,dimension(:,:),allocatable::ve1,ve2,ve3,ve4
!-------- Parametres shared
    double precision::ddls
    double precision,dimension(:),allocatable::str00
    integer::cpts
    double precision,dimension(:,:),allocatable::H_hess0,HIH0
    double precision,dimension(2)::LCVs,shapeweibs,scaleweibs,k0s
!     double precision,dimension(:),allocatable::x1Outs,x2Outs
    double precision,dimension(:,:),allocatable::xTOuts !en plus
!     double precision,dimension(:,:),allocatable::lamOuts,lam2Outs
    double precision,dimension(:,:,:),allocatable::lamTOuts !en plus
!     double precision,dimension(100,3)::suOuts,su2Outs
    double precision,dimension(100,3,1)::suTOuts !en plus
!     double precision,dimension(100)::xSu1s,xSu2s
    double precision,dimension(100,1)::xSuTs !en plus
    double precision,dimension(:),allocatable::zis
    double precision,dimension(nobsEvent(2))::Resmartingales,frailtypreds,frailtysds,frailtyvars
    double precision,dimension(:),allocatable::linearpreds,martingaleCoxs,times
    integer::timedepMul,nbinnerknots0Mul,qorder0Mul
    integer,dimension(:),allocatable::filtretps0Mul
    double precision,dimension(0:100,1)::BetaTpsMatMul
    double precision,dimension(3)::EPS
! End other variables
!----------------------------------------------------------------------------
! (1) Unpack Model Controls
   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'\nmultiveJoint.f90:: Unpacking model controls.'
   !debug: write(1,*)'multiveJoint.f90:: controls = ', controls
   !debug: write(1,*)'multiveJoint.f90:: nobsEvent = ', nobsEvent
   !debug: write(1,*)'multiveJoint.f90:: k0 = ', k0
    !write(1,*)'multiveJoint.f90:: tt00 = ', tt00
    !write(1,*)'multiveJoint.f90:: tt10 = ', tt10
    !write(1,*)'multiveJoint.f90:: tt0meta0 = ', tt0meta0
    !write(1,*)'multiveJoint.f90:: tt1meta0 = ', tt1meta0
    !write(1,*)'multiveJoint.f90:: ic0 = ', ic0
    !write(1,*)'multiveJoint.f90:: icmeta0 = ', icmeta0
    !write(1,*)'multiveJoint.f90:: groupe0 = ', groupe0
    !write(1,*)'multiveJoint.f90:: groupe0meta = ', groupe0meta
    !write(1,*)'multiveJoint.f90:: groupe0dc = ', groupe0dc
    !write(1,*)'multiveJoint.f90:: tt0dc0 = ', tt0dc0
    !write(1,*)'multiveJoint.f90:: tt1dc0 = ', tt1dc0
    !write(1,*)'multiveJoint.f90:: icdc0 = ', icdc0
    !write(1,*)'multiveJoint.f90:: icdc20 = ', icdc20
   !debug: write(1,*)'multiveJoint.f90:: nbvar = ', nbvar
    !write(1,*)'multiveJoint.f90:: vax0 = ', vax0
    !write(1,*)'multiveJoint.f90:: vaxdc0 = ', vaxdc0
    !write(1,*)'multiveJoint.f90:: vaxmeta0 = ', vaxmeta0
    !write(1,*)'multiveJoint.f90:: vaxdc20 = ', vaxdc20
   !debug: write(1,*)'multiveJoint.f90:: noVarEvent = ', noVarEvent
   !debug: write(1,*)'multiveJoint.f90:: np = ', np
   !debug: write(1,*)'multiveJoint.f90:: b = ', b
   !debug: write(1,*)'multiveJoint.f90:: nz0 = ', nz0
   !debug: write(1,*)'multiveJoint.f90:: ghNodes0 = ', ghNodes0
   !debug: write(1,*)'multiveJoint.f90:: ghWeights0 = ', ghWeights0
   !debug: close(1) 

    maxit0=controls(1)
    maxit00=controls(1)
    initialize = controls(2)
    typeof0 = controls(3)
    equidistant0 = controls(4)
    irep = controls(5)
    ag0 = controls(6)
    nbinterv0=controls(7)
    nz0 = controls(8)
    event2_ind0 = controls(9)
    terminal2_ind0 = controls(10)
    ghPoints = controls(11)
    typeJoint = controls(12)

    nzloco=nz0
    nzdc=nz0
    nzmeta=nz0
    nz1=nz0
    nz2=nz0
    nz3=nz0

    nsujet0=nobsEvent(1)
    ng0=nobsEvent(2) 
    nsujetmeta0=nobsEvent(3)
    nva10=nbvar(1)
    nva20=nbvar(2)
    nva30=nbvar(3) ! LE: change this?
    nva40=nbvar(4)
    noVar1=noVarEvent(1)
    noVar2=noVarEvent(2)
    noVar3=noVarEvent(3)
    noVar4=noVarEvent(4)

    mt1=mtEvent(1)
    mt2=mtEvent(2)
    mt3=mtEvent(3)
    mt4=mtEvent(4)
    mt11=mt1Event(1)
    mt12=mt1Event(2)
    mt13=mt1Event(3)
    mt14=mt1Event(4)
    kappaCV = 0.d0
    equidistant=equidistant0
    maxiter = maxit0 !maxiter joint
    typeof = typeof0
    ier = 0
    istop = 0

!convergence criteria
    epsa = tolerance0(1) ! change in parameter values
    epsb = tolerance0(2) ! change in likelihood
    epsd = tolerance0(3) ! change in likelihood derivative. most difficult to achieve? derivative of likelihood

!----------------------------------------------------------------------
! (2) Parameter Initialization

   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: (2) Initializing parameters.'
   !debug: close(1) 


    select case(typeof)
        case(0) ! splines
            np1 = nz0 + 2 + nva10 + 1 
            np2 = nz0 + 2 + nva20 + 1 
            np3 = nz0 + 2 + nva30 + 1
            np4 = nz0 + 2 + nva40 + 1
        case(1) ! peicewise
            k0=0.d0
            np1 = nbinterv0 + nva10 + 1 
            np2 = nbinterv0 + nva20 + 1 
            np3 = nbinterv0 + nva30 + 1
            np4 = nbinterv0 + nva40 + 1
        case(2) ! weibull
            k0=0.d0
            np1 = 2 + nva10 +1 
            np2 = 2 + nva20 +1
            np3 = 2 + nva30 +1
            np4 = 2 + nva40 +1
            !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
            !write(1,*)'multiveJoint.f90:: nva10...nva40 = ', nva10, nva20, nva30, nva40
            !close(1) 
            !np2joint = 4 + nva10 + nva20 + 2 ! four baseline hazard parameters, frailty variance, and link (alpha)
    end select
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: np1...np4 = ', np1, np2, np3, np4
    !close(1) 

    if(event2_ind0==1 .AND. terminal2_ind0==1)then
       !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
       !debug: write(1,*) 'Two recurrent and two terminal events not yet implemented.'
       !debug: close(1) 
        goto 1000
    endif

    if((event2_ind0.eq.1).AND.(terminal2_ind0.eq.0))then
       !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
       !debug: write(1,*) 'Two recurrent and one terminal event not yet implemented.'
       !debug: close(1) 
        goto 1000
    endif

    allocate(b01(np1),b02(np2),b03(np3),b04(np4))!,b02joint(np2joint))
    !if(event2_ind0.eq.1)then
    !    allocate(b03(np3))
    !endif

! Without Initializaiton
       !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
       !debug: write(1,*)'multiveJoint.f90:: Parameter initialization skipped.'
       !debug: close(1) 
        b01 = 0.25d0 ! recurrent
        !if(event2_ind0.eq.1)then
        b03 = 0.25d0 ! meta
        !end if
        b02 = 0.25d0 ! loco
        b04 = 0.25d0 ! loco

    if(typeof.ne.0) kappaCV=1.d0 

    !k0 = kappaCV
    
   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: End parameter initialization.'
   !debug: close(1) 
! End (1):  Initializing parameters.
!------------------------------------------------------------------

    Res_martingale=0.d0
    Res_martingaledc=0.d0
    Res_martingale2=0.d0
    !Res_martingaledc2 = 0.d0
    frailtypred=0.d0
    frailtypred2=0.d0
    frailtyvar=0.d0
    frailtyvar2=0.d0
    frailtyCorr=0.d0
    linearpred=0.d0
    linearpreddc=0.d0
    linearpreddc2 = 0.d0
    linearpredM=0.d0

    if(typeof.eq.0)then !splines
        ziOut1=0.d0
        ziOutdc=0.d0
        ziOutmeta=0.d0
    endif

    
    if(typeof.ne.0)then
        nbinterv = nbinterv0
    endif

    vectn=nz0+2
	
! This sets the delta used to approximate the derivatives
! 1 gives differences of th=1.d-3
! 3 gives differences of th=1.d-5
    model = 3
    
                
    lrs = 0.d0
    moy_peh0 = 0.d0
    moy_peh1 = 0.d0
    
    nb_echec = 0
    nb_echecor = 0
    nb0recu = 0
    moyrecu =0.d0             
        
    ngmax=ng0
    ng=ng0
    nsujetmax=nsujet0
    nsujet=nsujet0
    
   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: Allocation (1).'
   !debug: close(1) 

    allocate(ResidusRec(ngmax),Residusdc(ngmax),&
    !ResidusRec2(ngmax),&
    Rrec(ngmax),Nrec(ngmax),&
    Rdc(ngmax),Ndc(ngmax),&
    !Rrec2(ngmax),Nrec2(ngmax),&
    !Rdc2(ngmax), Ndc2(ngmax),&
    vuu(2),nig(ngmax),&!nigmeta(ngmax),
    cdc(ngmax),cdc2(ngmax),t0dc(ngmax),t1dc(ngmax),&
    t0(nsujet),t1(nsujet), g(nsujet), c(nsujet),aux(2*nsujet),&
    aux1(ngmax),aux2(ngmax),res1(ngmax),res4(ngmax),res3(ngmax),mi(ngmax))!,res1meta(ngmax),res3meta(ngmax))
    if(event2_ind.eq.1)then
        allocate(nigmeta(ngmax),res1meta(ngmax),res3meta(ngmax),ResidusRec2(ngmax),&
        Rrec2(ngmax),Nrec2(ngmax))
    endif

    shape_weib = 0.d0
    scale_weib = 0.d0

    if(event2_ind.eq.1)then
        nsujetmeta=nsujetmeta0
        nsujetmetamax=nsujetmeta0
        allocate(t0meta(nsujetmeta),t1meta(nsujetmeta),cmeta(nsujetmeta),gmeta(nsujetmeta),auxmeta(2*nsujetmeta))
        ndatemeta=2*nsujetmeta 
    endif
    ndatemaxdc=2*ng0     
    
   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: Allocation (2).'
   !debug: close(1) 
    if(typeof.eq.0)then
        allocate(nt0dc(ngmax),nt1dc(ngmax),nt0(nsujetmax),nt1(nsujetmax),mm3dc(ndatemaxdc),mm2dc(ndatemaxdc),&
        mm1dc(ndatemaxdc),mmdc(ndatemaxdc),im3dc(ndatemaxdc),im2dc(ndatemaxdc),im1dc(ndatemaxdc),imdc(ndatemaxdc))
        if(event2_ind.eq.1)then
            allocate(nt0meta(nsujetmeta),nt1meta(nsujetmeta),mm3meta(ndatemeta),immeta(ndatemeta),mm2meta(ndatemeta),&
            mm1meta(ndatemeta),mmmeta(ndatemeta),im3meta(ndatemeta),im2meta(ndatemeta),im1meta(ndatemeta))
        endif
    endif  

    nst=3
    ni=0
!---start of iteration simulations       
    id=1
    cptni=0
    cptni1=0
    cptni2=0
    biais_moy=0.d0
    cptbiais=0
    cptaux=0
    cptauxdc=0
    ij=0
    kk=0
    groupe=0
    groupemeta=0
    n=0
    !nzloco=0
    !nzdc=0
    !nzmeta=0
!LE: not sure why these were reset to 0, but it is causing problems

    effet=1
    res01(1)=0.d0
    res01(2)=0.d0
!------------  between no file and subject number -----        
    nvarmax=ver
    
   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: Allocation (3).'
   !debug: close(1) 
    allocate(vax(nva10),vaxdc(nva20),vaxdc2(nva40))
    if(event2_ind0 .eq.1)then
        allocate(vaxmeta(nva30))
    endif      
    nva1=nva10
    nva2=nva20
    nva3=nva30
    nva4=nva40
    nva = nva1+nva2+nva3+nva4
    nvarmax=nva
    event2_ind = event2_ind0
    terminal2_ind = terminal2_ind0
    
   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: Allocation (4).'
   !debug: close(1) 
    allocate(ve(nsujetmax,nvarmax),vedc(ngmax,nvarmax),&!vemeta(nsujetmeta,nvarmax),&
    vedc2(ngmax,nvarmax), &
    ve1(nsujetmax,nva1),ve2(ngmax,nva2),&!ve3(nsujetmeta,nva3),
    ve4(ngmax,nva4),&
    filtre(nva10),filtre2(nva20),&!filtre3(nva30),
    filtre4(nva40))
    nig=0

    if(event2_ind0 .eq. 1)then
        allocate(vemeta(nsujetmeta,nvarmax),ve3(nsujetmeta,nva3),filtre3(nva30))
        nigmeta=0 
    endif
    
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: noVar1, noVar2, noVar3, noVar4 = ', &
    !  noVar1, noVar2, noVar3, noVar4
    !close(1) 
! Recurrent 1
    if(noVar1.eq.1)then 
        filtre=0
        nva1=0
    else
        filtre=1
    endif
  
! Terminal 1
    if(noVar2.eq.1)then 
        filtre2=0
        nva2=0
    else
        filtre2=1
    endif  
    
    if((noVar1.eq.1).or.(noVar2.eq.1))then
        nva = nva1+nva2 
    endif

! Recurrent 2
    if(event2_ind0.eq.1)then
        if(noVar3.eq.1)then 
            filtre3=0
            nva3=0
        else
            filtre3=1
        endif
    endif

! Terminal 2
    if(noVar4.eq.1)then 
        filtre4=0
        nva4=0
    else
        filtre4=1
    endif
    
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: filtre1, filtre2, filtre3, filtre4 = ', &
    !  filtre, filtre2, filtre3, filtre4
    !close(1) 

!----------------------------------------------------------
! (3) Begin Copy Data
    maxt = 0.d0
    mint = 0.d0

!    maxtmeta = 0.d0
!    mintmeta = 0.d0

    maxtdc = 0.d0
    mintdc = 0.d0

    cpt = 0
    cptmeta =0 
    cptcens = 0
    cptcensmeta = 0
    cpt_dc = 0
    cpt_dc2 = 0
    k = 0
    cptstr1 = 0
    cptstr2 = 0


!----- Copy Terminal event data
   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: (3a) Copying data for terminal events.'
    !write(1,*)'multiveJoint.f90:: vaxdc0(1:100,1:nva2)=',vaxdc0(1:100,1:nva2)
   !debug: close(1) 

    do k = 1,ngmax
        if(k.eq.1)then
            mintdc = tt0dc0(k) ! min assignment just once 
        endif

        tt0dc=tt0dc0(k) ! copy to temporary single values 
        tt1dc=tt1dc0(k) ! (numbers cannot be copied directly from one vector to another)
        icdc=icdc0(k)
        icdc2=icdc20(k)
        groupe=groupe0dc(k)
        
        t0dc(k) =  tt0dc 
        t1dc(k) = tt1dc  
        cdc(k) = icdc
        cdc2(k) = icdc2
        cpt_dc = cpt_dc + icdc ! cpt_dc = count of observed terminal events
        cpt_dc2 = cpt_dc2 + icdc2 
        if(tt0dc.gt.0.d0)then
            cptauxdc=cptauxdc+1
        endif
        
        ! Copy Model Matrices
        do j=1,nva20        
            vaxdc(j)=vaxdc0(k,j)
        enddo                
        do j=1,nva40        
            vaxdc2(j)=vaxdc20(k,j)
        enddo                               

        iii = 0
        iii2 = 0   
        do ii = 1,nva20
            if(filtre2(ii).eq.1)then
                iii2 = iii2 + 1
                vedc(k,iii2) = dble(vaxdc(ii))
            endif
        end do   

        iii = 0
        iii2 = 0
        do ii = 1,nva40
            if(filtre4(ii).eq.1)then
                iii2 = iii2 + 1
                vedc2(k,iii2) = dble(vaxdc2(ii))
            endif
        end do

        ! Calculate maximum and minimum observed times
        if(maxtdc.lt.t1dc(k))then
            maxtdc = t1dc(k)
        endif
        if(mintdc.gt.t0dc(k))then
            mintdc = t0dc(k)
        endif
    end do

    if(typeof .ne. 0)then 
        cens = maxtdc
    endif

    k = 0
    cptstr1 = 0
    cptstr2 = 0
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: vedc(1:100,1:nva2)=',vedc(1:100,1:nva2)
    !close(1) 

! ----------------------------------------
! Copy Recurrent Event 1 Data
   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: (3b) Copying data for recurrent event 1.'
    !write(1,*)'multiveJoint.f90:: nsujet = ', nsujet
    !write(1,*)'multiveJoint.f90:: nva10 = ', nva10
   !debug: close(1) 

    do i = 1,nsujet     !sur les observations

        !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
        !write(1,*)'multiveJoint.f90:: groupe0(i)=', groupe0(i)
        !write(1,*)'multiveJoint.f90:: tt00(i)=', tt00(i)
        !write(1,*)'multiveJoint.f90:: tt10(i)=', tt10(i)
        !write(1,*)'multiveJoint.f90:: ic0(i)=', ic0(i)
        !close(1) 

        if(i.eq.1)then
            !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
            !write(1,*)'multiveJoint.f90:: setting mint'
            !close(1) 
            mint = tt00(i) ! affectation du min juste une fois
            !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
            !write(1,*)'multiveJoint.f90::  mint = ', mint
            !close(1) 
        endif

        tt0=tt00(i)
        !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
        !write(1,*)'multiveJoint.f90:: tt0=', tt0
        !close(1) 
        tt1=tt10(i)
        !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
        !write(1,*)'multiveJoint.f90:: tt1=', tt1
        !close(1) 
        ic=ic0(i)
        !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
        !write(1,*)'multiveJoint.f90:: ic=', ic
        !close(1) 
        groupe=groupe0(i)
        !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
        !write(1,*)'multiveJoint.f90:: groupe=', groupe
        !close(1) 
        t0(i) = tt0
        !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
        !write(1,*)'multiveJoint.f90:: t0(i)=', t0(i)
        !close(1) 
        t1(i) = tt1
        !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
        !write(1,*)'multiveJoint.f90:: t1(i)=', t1(i)
        !close(1) 
        g(i) = groupe
        !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
        !write(1,*)'multiveJoint.f90:: g(i)=', g(i)
        !close(1) 
        iii = 0
        iii2 = 0

        do j=1,nva10
            vax(j)=vax0(i,j)
        end do
        
        !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
        !write(1,*)'multiveJoint.f90:: vax copied. vax = ', vax
        !close(1) 

        if(tt0.gt.0.d0)then
            cptaux=cptaux+1
        endif
        
        if(ic.eq.1)then
            cpt = cpt + 1
            c(i)=1
            nig(groupe) = nig(groupe)+1 ! nb d event recurr dans un groupe
  
            do ii = 1,nva10
                if(filtre(ii).eq.1)then
                    iii = iii + 1
                    ve(i,iii) = dble(vax(ii)) !ici sur les observations
                endif
            end do
        else 
            if(ic.eq.0)then
                cptcens=cptcens+1
                c(i) = 0 

                do ii = 1,nva10
                    if(filtre(ii).eq.1)then
                    iii = iii + 1
                    ve(i,iii) = dble(vax(ii))
                    endif
                end do 
            endif
        endif

        !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
        !write(1,*)'multiveJoint.f90:: cptcens=', cptcens
        !write(1,*)'multiveJoint.f90:: nig(groupe)=', nig(groupe)
        !close(1) 

        if(maxt.lt.t1(i))then
            maxt = t1(i)
        endif
        if(mint.gt.t0(i))then
            mint = t0(i)
        endif
    end do 

    nsujet=i-1

! Copy Recurrent Event 2 Data
    if(event2_ind0.eq.1)then
       !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
       !debug: write(1,*)'multiveJoint.f90:: (3c) Copying data for recurrent event 2.'
       !debug: close(1) 
        do i = 1,nsujetmeta     !sur les observations

            if(i.eq.1)then
                mintmeta = tt0meta0(i) ! affectation du min juste une fois
            endif

            tt0meta=tt0meta0(i)
            tt1meta=tt1meta0(i)
            icmeta=icmeta0(i)
            groupemeta=groupe0meta(i)

            do j=1,nva30
                vaxmeta(j)=vaxmeta0(i,j)
            enddo

            if(tt0meta.gt.0.d0)then
                cptauxmeta=cptauxmeta+1
            endif

            if(icmeta.eq.1)then
                cptmeta = cptmeta + 1
                cmeta(i)=1
                t0meta(i) = tt0meta 
                t1meta(i) = tt1meta  
                t1meta(i) = t1meta(i)
                gmeta(i) = groupemeta
                nigmeta(groupemeta) = nigmeta(groupemeta)+1 ! nb d event recurr dans un groupe
                iii = 0
                iii2 = 0

                do ii = 1,nva30
                    if(filtre3(ii).eq.1)then
                        iii = iii + 1
                        vemeta(i,iii) = dble(vaxmeta(ii)) !ici sur les observations
                    endif
                end do
            else 
                if(icmeta.eq.0)then
                    cptcens=cptcens+1
                    cmeta(i) = 0 
                    iii = 0
                    iii2 = 0
                    do ii = 1,nva30
                        if(filtre3(ii).eq.1)then
                            iii = iii + 1
                            vemeta(i,iii) = dble(vaxmeta(ii))
                        endif
                    end do 
                    t0meta(i) =  tt0meta
                    t1meta(i) = tt1meta
                    t1meta(i) = t1meta(i)
                    gmeta(i) = groupemeta
                endif
            endif
            if(maxtmeta.lt.t1meta(i))then
                maxtmeta = t1meta(i)
            endif
            if(mintmeta.gt.t0meta(i))then
                mintmeta = t0meta(i)
            endif
        end do 
    endif
    !write(*,*) 'cpt = ', cpt
    !write(*,*) 'cptmeta = ', cptmeta
    !write(*,*) 'cpt_dc =', cpt_dc
    !write(*,*) 'cpt_dc2 =', cpt_dc2
! End Copy Data
!------------------------------------------------------------------

    ndatemax=2*nsujet

    allocate(date(ndatemax),datedc(ndatemaxdc))!,datemeta(ndatemeta))

    ! LE: changed datedc(ndatemax) to datedc(ndatemaxdc)
    if(event2_ind0.eq.1)then
        allocate(datemeta(ndatemeta))
    endif
    
    if(typeof.eq.0)then
        allocate(mm3(ndatemax),mm2(ndatemax) &
        ,mm1(ndatemax),mm(ndatemax),im3(ndatemax),im2(ndatemax),im1(ndatemax),im(ndatemax))
    endif

!--------------------------------------------------
! (4) Sort Unique Dates for each event type
   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: (4) Calculating grid of dates for spline evaluation.'
    !write(1,*)'multiveJoint.f90:: maxtdc = ', maxtdc
    !write(1,*)'multiveJoint.f90:: ngmax = ', ngmax
    !write(1,*)'multiveJoint.f90:: size(aux) = ', size(aux)
   !debug: close(1) 

! Terminal Event 1 (Same grid for Terminal Event 2)
! aux: the ordered unique values of c(t0dc, t1dc).
! the last value is repeated if there are fewer unique values than 2*ng.

! datedc : the ordered unique values of c(t0dc, t1dc).
! ndatedc: the number of ordered unique values in c(t0dc, t1dc).
    mindc = 0.d0
    maxdc = maxtdc 
    do i = 1,2*ngmax  
        do k = 1,ngmax 
            if((t0dc(k).ge.mindc))then
                if(t0dc(k).lt.maxdc)then
                    maxdc = t0dc(k) 
                endif
            endif
            if((t1dc(k).ge.mindc))then
                if(t1dc(k).lt.maxdc)then 
                    maxdc = t1dc(k)
                endif
            endif
        end do   
        aux(i) = maxdc
        mindc = maxdc + 1.d-12
        maxdc = maxtdc
    end do

    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: mindc = ', mindc
    !write(1,*)'multiveJoint.f90:: maxdc = ', maxdc
    !close(1) 

    datedc(1) = aux(1)
    k = 1
    do i=2,2*ngmax
        if(aux(i).gt.aux(i-1))then
            k = k+1
            datedc(k) = aux(i) 
        endif
    end do 
    
    if(typeof == 0)then
        ndatedc = k   
    endif 
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Size of date grid for terminal 1 = '
    !write(1,*) ndatedc
    !write(1,*)'multiveJoint.f90:: Dates for recurrent 1 = '
    !write(1,*) datedc
    !close(1) 

! Recurrent Event 1
    min = 0.d0
    aux =0.d0
    max = maxt
    
    do i = 1,2*nsujet
        do k = 1,nsujet
            if((t0(k).ge.min))then
                if(t0(k).lt.max)then
                    max = t0(k)
                endif
            endif
            if((t1(k).ge.min))then
                if(t1(k).lt.max)then
                    max = t1(k)
                endif
            endif
        end do   
        aux(i) = max
        min = max + 1.d-12
        max = maxt
    end do

    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: min = ', min
    !write(1,*)'multiveJoint.f90:: max = ', max
    !close(1) 

    date(1) = aux(1)
    k = 1
    do i=2,2*nsujet
        if(aux(i).gt.aux(i-1))then
            k = k+1
            date(k) = aux(i)
        endif
    end do 

    if(typeof==0)then
        ndate = k
    endif
    
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Size of date grid for recurrent 1 = '
    !write(1,*) ndate
    !close(1) 

    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Size of date grid for recurrent 1 = '
    !write(1,*) ndate
    !write(1,*)'multiveJoint.f90:: Dates for recurrent 1 = '
    !write(1,*) date
    !close(1) 

! Recurrent Event 2 (Second recurrent event not currently supported)

    !min = 0.d0
    !auxmeta =0.d0
    !max = maxtmeta

    !do i = 1,2*nsujetmeta
    !    do k = 1,nsujetmeta
    !        if((t0meta(k).ge.min))then
    !            if(t0meta(k).lt.max)then
    !                max = t0meta(k)
    !            endif
    !        endif
    !        if((t1meta(k).ge.min))then
    !            if(t1meta(k).lt.max)then
    !                max = t1meta(k)
    !            endif
    !        endif
    !    end do   
    !    auxmeta(i) = max
    !    min = max + 1.d-12
    !    max = maxtmeta
    !end do

    !datemeta(1) = auxmeta(1)
    !k = 1
    !do i=2,2*nsujetmeta
    !    if(auxmeta(i).gt.auxmeta(i-1))then
    !        k = k+1
    !        datemeta(k) = auxmeta(i)
    !    endif
    !end do 

    !if(typeof==0)then
    !    ndatemeta = k
    !endif
! End (4) Sort unique dates for each type of event
!--------------------------------------------

!--------------------------------------------
! (5) Construction of Nodes for Splines (only equidistant allowed at this point)
    if(typeof.eq.0)then
       !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
       !debug: write(1,*)'multiveJoint.f90:: (5) Constructing nodes for splines.'
        !write(1,*)'multiveJoint.f90:: nzloco = ', nzloco
        !write(1,*)'multiveJoint.f90:: nzdc = ', nzdc
       !debug: close(1) 

        if(equidistant.eq.0)then ! spline nodes are placed on percentiles (not currently allowed)
            ! recurrent
            !i=0
            !j=0
            !do i=1,nsujet
            !    if(t1(i).ne.(0.d0).and.c(i).eq.1)then
            !        j=j+1
            !    endif
            !end do
            !nbrecu=j

            !allocate(t2(nbrecu))
            !j=0
            !do i=1,nsujet
            !    if(t1(i).ne.(0.d0).and.c(i).eq.1)then
            !        j=j+1
            !        t2(j)=t1(i)
            !    endif
            !end do

            !allocate(zi(-2:(nzloco+3)))
            !ndate = k
            !zi(-2) = mint
            !zi(-1) = mint
            !zi(0) = mint
            !zi(1) = mint
            !j=0
            !do j=1,nzloco-2
            !    pord = dble(j)/(dble(nzloco)-1.d0)
            !    call percentile3(t2,nbrecu,pord,zi(j+1))
            !end do
            !zi(nzloco) = maxt
            !zi(nzloco+1) = maxt
            !zi(nzloco+2) = maxt
            !zi(nzloco+3) = maxt
            !ziOut1 = zi
            !deallocate(t2)

            ! death
            !i=0
            !j=0
            !do i=1,ng
            !    if(t1dc(i).ne.(0.d0).and.cdc(i).eq.1)then
            !        j=j+1
            !    endif
            !end do
            !nbdeces=j

            !allocate(t3(nbdeces))
            !j=0
            !do i=1,ng
            !    if(t1dc(i).ne.(0.d0).and.cdc(i).eq.1)then
            !        j=j+1
            !        t3(j)=t1dc(i)
            !    endif
            !end do

            !allocate(zidc(-2:(nzdc+3)))
            !zidc(-2) = mintdc
            !zidc(-1) = mintdc
            !zidc(0) = mintdc
            !zidc(1) = mintdc
            !j=0
            !do j=1,nzdc-2
            !    pord = dble(j)/(dble(nzdc)-1.d0)
            !    call percentile3(t3,nbdeces,pord,zidc(j+1))
            !end do
            !zi(nzdc) = maxtdc
            !zi(nzdc+1) = maxtdc
            !zi(nzdc+2) = maxtdc
            !zi(nzdc+3) = maxtdc
            !ziOutdc = zidc
            !deallocate(t3)

            ! meta
            !i=0
            !j=0
            !do i=1,nsujetmeta
            !    if(t1meta(i).ne.(0.d0).and.cmeta(i).eq.1)then
            !        j=j+1
            !    endif
            !end do
            !nbrecumeta=j

            !allocate(t4(nbrecumeta))
            !j=0
            !do i=1,nsujetmeta
            !    if(t1meta(i).ne.(0.d0).and.cmeta(i).eq.1)then
            !        j=j+1
            !        t4(j)=t1meta(i)
            !    endif
            !end do

            !allocate(zimeta(-2:(nzmeta+3)))
            !zimeta(-2) = mintmeta
            !zimeta(-1) = mintmeta
            !zimeta(0) = mintmeta
            !zimeta(1) = mintmeta
            !j=0
            !do j=1,nzmeta-2
            !    pord = dble(j)/(dble(nzmeta)-1.d0)
            !    call percentile3(t4,nbrecumeta,pord,zimeta(j+1))
            !end do
            !zimeta(nzmeta) = maxtmeta
            !zimeta(nzmeta+1) = maxtmeta
            !zimeta(nzmeta+2) = maxtmeta
            !zimeta(nzmeta+3) = maxtmeta
            !ziOutmeta = zimeta
            !deallocate(t4)

        else ! nodes are equidistant
            ! h: duration between two nodes

            ! recurrent
            nzmax=nzloco+3
            allocate(zi(-2:nzmax))
            zi(-2) = date(1)
            zi(-1) = date(1)
            zi(0) = date(1)
            zi(1) = date(1)
            h = (date(ndate)-date(1))/dble(nzloco-1)
            do i=2,nzloco-1
                zi(i) =zi(i-1) + h
            end do
            zi(nzloco) = date(ndate)
            zi(nzloco+1)=zi(nzloco)
            zi(nzloco+2)=zi(nzloco)
            zi(nzloco+3)=zi(nzloco)
            ziOut1 = zi

            ! terminal 1
            allocate(zidc(-2:(nzdc+3)))
            zidc(-2) = datedc(1)
            zidc(-1)= datedc(1)
            zidc(0)= datedc(1)
            zidc(1)= datedc(1)
            h =(datedc(ndatedc)-datedc(1))/dble(nzdc-1)
            do i=2,nzdc-1
                zidc(i) =zidc(i-1)+h
            end do
            zidc(nzdc) = datedc(ndatedc)
            zidc(nzdc+1)=zidc(nzdc)
            zidc(nzdc+2)=zidc(nzdc)
            zidc(nzdc+3)=zidc(nzdc)
            ziOutdc = zidc

            ! second recurrent event
            if(event2_ind0.eq.1)then
                allocate(zimeta(-2:(nzmeta+3)))
                zimeta(-2) = datemeta(1)
                zimeta(-1)= datemeta(1)
                zimeta(0)= datemeta(1)
                zimeta(1)= datemeta(1)
                h =(datemeta(ndatemeta)-datemeta(1))/dble(nzmeta-1)
                do i=2,nzmeta-1
                    zimeta(i) =zimeta(i-1)+h
                end do
                zimeta(nzmeta) = datemeta(ndatemeta)
                zimeta(nzmeta+1)=zimeta(nzmeta)
                zimeta(nzmeta+2)=zimeta(nzmeta)
                zimeta(nzmeta+3)=zimeta(nzmeta)
                ziOutmeta = zimeta
            endif
            
            ! terminal 2
            ! we do not need another section for terminal 2. 
            ! it will be the same as terminal 1.

         endif
    endif
! End (5) construction of nodes for splines
!--------------------------------------------------------------

!---------------------------------------------------------------------- 
! (6) Map each observed time to an index in the ordered date vector (splines only)
    if(typeof == 0)then
       !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
       !debug: write(1,*)'multiveJoint.f90:: (6) Mapping observations to date grid (splines only).'
       !debug: close(1) 

! nb0recu: number of subjects without terminal event
! moyrecu: number of recurrent events observed
! indictronqdc: indicator for truncation
! nt0dc: the date index for the start time for person k
! nt1dc: the date index for the terminal time for person k
! nt0: the date index for the start at risk time for a given event
! nt1: the date index for the terminal at risk time for a given event 

! Terminal 1: nt0dc, nt1dc (will also be used for Terminal 2)
        indictronqdc=0
        do k=1,ngmax
            if(nig(k).eq.0.00d0)then
                nb0recu = nb0recu + 1 
            endif
            moyrecu =  moyrecu + dble(nig(k))

            if(t0dc(k).eq.0.00d0)then
                nt0dc(k) = 0
            endif

            if(t0dc(k).ne.0.00d0)then
                indictronqdc=1
            endif

            do j=1,ndatedc ! LE, check for approximate equality
                if((t0dc(k).le.(datedc(j)+0.0001d0)).and.(t0dc(k).ge.(datedc(j)-0.0001d0)))then
                    nt0dc(k)=j
                endif
                if((t1dc(k).le.(datedc(j)+0.0001d0)).and.(t1dc(k).ge.(datedc(j)-0.0001d0)))then
                    nt1dc(k)=j
                endif
            end do
        end do
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Date indexes for terminal 1 = '
    !write(1,*) nt0dc
    !write(1,*)'multiveJoint.f90:: Dates for terminal 1 = '
    !write(1,*) datedc
    !write(1,*)'multiveJoint.f90:: Date indexes (2) for terminal 1 = '
    !write(1,*) nt1dc
    !close(1) 

! Recurrent 1: nt0,nt1
        indictronq=0
        do i=1,nsujet 
            if(t0(i).eq.0.00d0)then
                nt0(i) = 0
            endif
            if(t0(i).ne.0.00d0)then
                indictronq=1
            endif
            do j=1,ndate ! check ifeach date is approximately equal to one of the index dates
                if((t0(i).le.(date(j)+0.0001d0)).and.(t0(i).ge.(date(j)-0.0001d0)))then
                    nt0(i)=j
                endif
                if((t1(i).le.(date(j)+0.0001d0)).and.(t1(i).ge.(date(j)-0.0001d0)))then
                    nt1(i)=j
                endif
            end do
        end do 
        !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
        !write(1,*)'multiveJoint.f90:: Date indexes for recurrent 1 = '
        !write(1,*) nt0
        !write(1,*)'multiveJoint.f90:: Date indexes (2) for recurrent 1 = '
        !write(1,*) nt1
        !close(1) 

! Recurrent 2: nt0meta, nt1meta
!    if(event2_ind0.eq.1)then
!        indictronqmeta=0
!        do i=1,nsujetmeta

!            if(t0meta(i).eq.0.d0)then
!                    nt0meta(i) = 0
!            endif
!            if(t0meta(i).ne.0.d0)then
!                indictronqmeta=1
!            endif
!            do j=1,ndatemeta
!                if(datemeta(j).eq.t0meta(i))then
!                    nt0meta(i)=j
!                endif
!                if(datemeta(j).eq.t1meta(i))then
!                    nt1meta(i)=j
!                endif
!            end do
!        end do
!    endif
    endif

! End (6): Map each observed time to an index in the ordered date vector
!----------------------------------------------------------------------

!------------------------------------------------------
! (7) Calculating Penalties (Splines only)
     if(typeof.eq.0)then
       !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
       !debug: write(1,*)'multiveJoint.f90:: (7) Calculating spline values on date grid. (splines only)'
       !debug: close(1) 

        call vecspli(ndate,ndatedc,ndatemeta)

        !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
        !write(1,*)'multiveJoint.f90:: Spline values calculated.'
        !close(1) 

        ! First Recurrent Event
        allocate(m3m3(nzmax),m2m2(nzmax),m1m1(nzmax),mmm(nzmax),m3m2(nzmax),m3m1(nzmax),&
        m3m(nzmax),m2m1(nzmax),m2m(nzmax),m1m(nzmax))

        ! Terminal Events
        allocate(m3m3b(nzdc),m2m2b(nzdc),m1m1b(nzdc), &
        mmmb(nzdc),m3m2b(nzdc),m3m1b(nzdc),m3mb(nzdc),m2m1b(nzdc),m2mb(nzdc),m1mb(nzdc),&
        m3m3d(nzdc),m2m2d(nzdc),m1m1d(nzdc), &
        mmmd(nzdc),m3m2d(nzdc),m3m1d(nzdc),m3md(nzdc),m2m1d(nzdc),m2md(nzdc),m1md(nzdc))

        ! Second recurrent event
        if(event2_ind0.eq.1)then
            allocate(m3m3c(nzmeta),m2m2c(nzmeta),m1m1c(nzmeta),mmmc(nzmeta),m3m2c(nzmeta),m3m1c(nzmeta),&
            m3mc(nzmeta),m2m1c(nzmeta),m2mc(nzmeta),m1mc(nzmeta))
        endif

        !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
        !write(1,*)'multiveJoint.f90:: Calculating penalties for splines.'
        !close(1) 

        call vecpenP(nzloco+2,zi,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m)
        call vecpenP(nzdc+2,zidc,m3m3b,m2m2b,m1m1b,mmmb,m3m2b,m3m1b,m3mb,m2m1b,m2mb,m1mb)
        if(event2_ind0.eq.1)then
            call vecpenP(nzmeta+2,zimeta,m3m3c,m2m2c,m1m1c,mmmc,m3m2c,m3m1c,m3mc,m2m1c,m2mc,m1mc)
        endif
        call vecpenP(nzdc+2,zidc,m3m3d,m2m2d,m1m1d,mmmd,m3m2d,m3m1d,m3md,m2m1d,m2md,m1md)

        !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
        !write(1,*)'multiveJoint.f90:: Spline penalties calculated.'
        !close(1) 
    endif
! End (7) calculating penalties (splines only)
!------------------------------------------------------
    npmax=np

    allocate(I_hess(npmax,npmax),H_hess(npmax,npmax),Hspl_hess(npmax,npmax) &
    ,PEN_deri(npmax,1),hess(npmax,npmax),v((npmax*(npmax+3)/2)),I1_hess(npmax,npmax) &
    ,H1_hess(npmax,npmax),I2_hess(npmax,npmax),H2_hess(npmax,npmax),HI2(npmax,npmax) & 
    ,HIH(npmax,npmax),IH(npmax,npmax),HI(npmax,npmax),BIAIS(npmax,1))

    if(typeof.ne.0)then
        allocate(vvv((npmax*(npmax+1)/2)))
    endif


!------------------------------------------------------
! (7b) calculating piecewise constant stuff

! This section removed because piecewise constant hazards will not be supported.
    !if(typeof == 1)then

!----->   LOCO  
    !    j=0
    !    do i=1,nsujet
    !        if(t1(i).ne.(0.d0).and.c(i).eq.1)then
    !            j=j+1
    !        endif
    !    end do
    !    nbrecu=j
    !    allocate(t2(nbrecu))
!
    !    j=0
    !    do i=1,nsujet
    !        if(t1(i).ne.(0.d0).and.c(i).eq.1)then
    !            j=j+1
    !            t2(j)=t1(i)
    !        endif
    !    end do
    !    indd=1
    !    do while (indd.eq.1)
    !        indd=0
    !        do i=1,nbrecu-1
    !            if(t2(i).gt.t2(i+1))then
    !                temp=t2(i)
    !                t2(i)=t2(i+1)
    !                t2(i+1)=temp
    !                indd=1
    !            endif
    !        end do
    !    end do
    !    ent=int(nbrecu/nbintervR)
    !    allocate(ttt(0:nbintervR))

    !    ttt(0)=0.d0
    !    ttt(nbintervR)=cens
    !    j=0
    !    do j=1,nbintervR-1
    !        if(equidistant.eq.0)then
    !            ttt(j)=(t2(ent*j)+t2(ent*j+1))/(2.d0)
    !        else
    !            ttt(j)=(cens/nbintervR)*j
    !        endif
    !    end do
    !    time = ttt

!----->   DECES  
! nbdeces: number of deaths

    !    j=0
    !    do i=1,ngmax
    !        if(t1dc(i).ne.(0.d0).and.cdc(i).eq.1)then
    !            j=j+1
    !        endif
    !    end do
    !    nbdeces=j

    !    allocate(t3(nbdeces))
    !    j=0
    !    do i=1,ngmax
    !        if(t1dc(i).ne.(0.d0).and.cdc(i).eq.1)then
    !            j=j+1
    !            t3(j)=t1dc(i)
    !        endif
    !    end do

    !    indd=1
    !    do while (indd.eq.1)
    !        indd=0
    !        do i=1,nbdeces-1
    !            if(t3(i).gt.t3(i+1))then
    !                temp=t3(i)
    !                t3(i)=t3(i+1)
    !                t3(i+1)=temp
    !                indd=1
    !            endif
    !        end do
    !    end do

    !    entdc=int(nbdeces/nbintervDC)
    !    allocate(tttdc(0:nbintervDC))
    !    tttdc(0)=0.d0
    !    tttdc(nbintervDC)=cens

    !    j=0
    !    do j=1,nbintervDC-1
    !        if(equidistant.eq.0)then
    !            tttdc(j)=(t3(entdc*j)+t3(entdc*j+1))/(2.d0)
    !        else
    !            tttdc(j)=(cens/nbintervDC)*j
    !        endif
    !    end do     
    !    timedc = tttdc

!----->   META
    !    j=0

    !    do i=1,nsujetmeta
    !        if(t1meta(i).ne.(0.d0).and.cmeta(i).eq.1)then
    !            j=j+1
    !        endif
    !    end do
    !    nbrecumeta=j
    !    allocate(t4(nbrecumeta))

    !    j=0
    !    do i=1,nsujetmeta
    !        if(t1meta(i).ne.(0.d0).and.cmeta(i).eq.1)then
    !            j=j+1
    !            t4(j)=t1meta(i)
    !        endif
    !    end do

    !    indd=1
    !    do while (indd.eq.1)
    !        indd=0
    !        do i=1,nbrecumeta-1
    !            if(t4(i).gt.t4(i+1))then
    !                temp=t4(i)
    !                t4(i)=t4(i+1)
    !                t4(i+1)=temp
    !                indd=1
    !            endif
    !        end do
    !    end do

    !    entmeta=int(nbrecumeta/nbintervM)

    !    allocate(tttmeta(0:nbintervM))

    !    tttmeta(0)=0.d0
    !    tttmeta(nbintervM)=cens

    !    j=0
    !    do j=1,nbintervM-1
    !        if(equidistant.eq.0)then
    !            tttmeta(j)=(t4(entmeta*j)+t4(entmeta*j+1))/(2.d0)
    !        else
    !            tttmeta(j)=(cens/nbintervM)*j
    !        endif
    !    end do
    !    timeM = tttmeta
    !    deallocate(t2,t3,t4)
    !endif
! End (7b) Piecewise constant stuff
!------------------------------------------------------

    ca=0.d0
    cb=0.d0
    dd=0.d0
    if(typeof.ne.0)then
        allocate(kkapa(4))
    endif
!------------------------------------------------------------------
! (8) Copy Initialized Parameters to b

!     write(*,*)''
!     write(*,*)'=============== Modele final ================='
!     write(*,*)'========= typeof = (0:Splines, 1:Piecewise, 2:Weibull) ',typeof    
!     write(*,*)'========= nombre de parametres ',np
!     write(*,*)'========= nombre de sujet, d obs meta, d obs loco ',ng,nsujetmeta,nsujet

! LE: I am not sure that these are necessary
    indic_alpha=4
    indic_eta=3
    indic_rho=0
    indic_a1=1
    indic_a2=2

    !write(*,*)'Number of variables', nva1, nva2, nva3, nva4, 'Total:', nva
    !write(*,*)''
    !write(*,*)'Number of parameters',np1, np2, np3, np4, 'Total:',np
    !write(*,*)''
    !write(*,*)'Initialized Parameter Vector: b = ', b

! END (8) Copy Initialized Parameters to b
!------------------------------------------------------------------

!------------------------------------------------------------------
! (9) Optimization


    allocate(ghNodes(ghPoints),ghWeights(ghPoints))
    ghNodes=ghNodes0
    ghWeights=ghWeights0
    res=0.d0

   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: (9) Optimizing Parameters.'
   !debug: write(1,*)'multiveJoint.f90:: starting b = ', b
   !debug: close(1) 

    if((event2_ind.eq.1).and.(terminal2_ind.eq.0))then
        !write(*,*) 'Model for two recurrent events and one terminal event.'
        select case(typeof)
            case(0)
                !call marq98(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaMultivSplines)
            case(2)
                !call marq98(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaMultivWeib)
        end select
    endif

    if((event2_ind.eq.0).and.(terminal2_ind.eq.1))then
       !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
       !debug: write(1,*) 'multiveJoint.f90:: Model for one recurrent event and two terminal events.'
       !debug: close(1) 
        select case(typeof)
            case(0)
                call marq98(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaMultivSplines)
            case(2)
                call marq98(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaMultivWeib)
        end select
    endif

    if((event2_ind.eq.1).and.(terminal2_ind.eq.1))then
        !write(*,*) 'Model for two recurrent events and two terminal events.'        
        select case(typeof)
            case(0)
                !call marq98(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaMultivSplines)
            case(2)
                !call marq98(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaMultivWeib)
        end select
    endif

   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: Optimization complete.'
   !debug: close(1) 

! END (9) Optimizaion
!-----------------------------------------------------------

    deallocate(ghNodes,ghWeights)

    if(typeof .ne. 0)then
        deallocate(kkapa)
    endif

    resOut=res
    critCV(1)=ier
    critCV(2)=istop

    if(istop.ne.1)then
        goto 1000
    endif

!-----------------------------------------------------------
! (10) Organize Hessian Matrix
   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: (10) Organizing hessian matricies...'
   !debug: close(1) 

    call multi(I_hess,H_hess,np,np,np,IH)
    call multi(H_hess,IH,np,np,np,HIH)

    do ss=1,npmax
        do sss=1,npmax
            HIHOut(ss,sss) = HIH(ss,sss)
            H_hessOut(ss,sss)= H_hess(ss,sss)
        end do  
    end do
! (10) Organize Hessian Matrix
!-----------------------------------------------------------

    !if((effet.eq.1).and.(ier.eq.-1))then
    !    v((np-nva-indic_alpha)*(np-nva-indic_alpha+1)/2)=10.d10
        ! what does this mean?
    !endif

    res01(effet+1)=res
!-----------------------------------------------------------
! (11) Calculate Lambda and survival estimates
   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: (11) Calculating hazard and survival functions...'
   !debug: close(1) 
    !select case(typeof)
        !case(0)
            !call distanceJ_splines(nzloco,nzdc,nzmeta,b,mt1,mt2,mt3,x1Out,lamOut,suOut,x2Out,lam2Out,su2Out,&
            !x3Out,lam3Out,su3Out)
        !case(1)
            !Call distanceJ_cpm(b,nbintervR+nbintervDC+nbintervM,mt1,mt2,mt3,x1Out,lamOut,xSu1,suOut,x2Out, &
            !lam2Out,xSu2,su2Out,x3Out,lam3Out,xSu3,su3Out)
        !case(2)
            !Call distanceJ_weib(b,np,mt1,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out,x3Out,lam3Out,xSu3,su3Out)
            !scale_weib(1) = etaR
            !shape_weib(1) = betaR
            !scale_weib(2) = etaD
            !shape_weib(2) = betaD
            !if(event2_ind0.eq.1)then
            !    scale_weib(3) = etaM
            !    shape_weib(3) = betaM
            !endif
            !scale_weib(4) = etaD2
            !shape_weib(4) = betaD2
    !end select
! END (11) Calculate Lambda and survival estimates
!-----------------------------------------------------------

!-----------------------------------------------------------
! (12) Calculate Likelihood Cross-Validation Criterion
! LCV(1) = The approximate like cross-validation Criterion
! LCV(2) = Akaike information Criterion 
!     calcul de la trace, pour le LCV (likelihood cross validation)
    !LCV=0.d0
    !if(typeof == 0)then
!        write(*,*)'The approximate like cross-validation Criterion in the non parametric case'
        !call multi(H_hess,I_hess,np,np,np,HI)    
        !do i =1,np
        !    LCV(1) = LCV(1) + HI(i,i)
        !end do
        !LCV(1) = (LCV(1) - resnonpen) / nsujet
        !else
!        write(*,*)'=========> Akaike information Criterion <========='
        !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
        !write(1,*)'multiveJoint.f90:: Calculating AIC...'
        !close(1) 

        !LCV(2) = (1.d0 / nsujet) *(np - resOut)
!        write(*,*)'======== AIC :',LCV(2)
    !endif
! END (12) Likelihood Cross-Validation Criterion
!-----------------------------------------------------------

1000 continue ! If an error is produced earlier, we jump to here


!-----------------------------------------------------------
! (13) Calculate Fitted Values
   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: (13) Calculating fitted values...'
   !debug: close(1) 

    !write(*,*)'=========== coefBeta loco =========='
    !coefBeta(1,:) = b((np-nva+1):(np-nva+nva1))
    !print*,coefBeta

    !write(*,*)'=========== coefBeta dc =========='
    !coefBetadc(1,:) = b((np-nva+nva1+1):(np-nva+nva1+nva2))
    !print*,coefBetadc

    !write(*,*)'=========== coefBeta meta =========='
    !if(Event2_ind0.eq.1)then
    !    coefBetaM(1,:) = b((np-nva+nva1+nva2+1):(np-nva+nva1+nva2+nva3))
        !print*,coefBetaM(1,:)
    !endif

    !write(*,*)'=========== coefBeta dc2 =========='
    !coefBetadc2(1,:) = b((np-nva4+1) : np)
    !print*,coefBetadc

    !do i=1,nsujet
    !   do j=1,nva1
    !        ve1(i,j)=ve(i,j)
    !    end do
    !end do

    !do i=1,ngmax
    !    do j=1,nva2
    !        ve2(i,j)=vedc(i,j)
    !    end do
    !end do

    !if(event2_ind0.eq.1)then
    !    do i=1,nsujetmeta
    !        do j=1,nva3
    !            ve3(i,j)=vemeta(i,j)
    !        end do
    !    end do
    !endif

    !do i=1,ngmax
    !    do j=1,nva4
    !        ve4(i,j)=vedc2(i,j)
    !    end do
    !end do

   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: Multiplying design matricies by beta....'
   !debug: close(1) 

    !Xbeta = matmul(coefBeta,transpose(ve1))
    !Xbetadc = matmul(coefBetadc,transpose(ve2))
!    Xbetadc2 = matmul(coefBetadc2,transpose(ve4))
    !if(event2_ind0.eq.1)then
    !    XbetaM = matmul(coefBetaM,transpose(ve3))
    !endif

    if((istop.eq.1).and.(effet.eq.1))then
!        print*,'======== Call Residus Martingale ==========='
        deallocate(I_hess,H_hess)

        !allocate(vres((2*(2+3)/2)),I_hess(2,2),H_hess(2,2))

        !effetres = effet

        !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
        !write(1,*)'multiveJoint.f90:: Calculating martingale residuals....'
        !close(1) 

        !Call Residus_Martingale_multive(b,np,funcpamultires,Res_martingale,Res_martingaledc,Res_martingale2,&
        !frailtypred,frailtypred2,frailtyvar,frailtyvar2,frailtyCorr)

        !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
        !write(1,*)'multiveJoint.f90:: Calculating individual specific log hazards....'
        !close(1) 

        !do i=1,nsujet
        !    linearpred(i)=Xbeta(1,i)+frailtypred(g(i))
        !end do

        !do i=1,ng
        !    linearpreddc(i)=Xbetadc(1,i)+alpha1*frailtypred(g(i))+alpha2*frailtypred2(gmeta(i))
        !    linearpreddc2(i)=Xbetadc2(1,i)+alpha1*frailtypred(g(i))+alpha2*frailtypred2(gmeta(i))
        !end do

        !do i=1,nsujetmeta
        !    linearpredM(i)=XbetaM(1,i)+frailtypred2(gmeta(i))
        !end do

        !deallocate(I_hess,H_hess,vres)
    else
        deallocate(I_hess,H_hess)
    endif
! END (13) Calculate Fitted Values
!-----------------------------------------------------------

!-----------------------------------------------------------
! (14) Deallocating All variables
   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: (14) Deallocating variables....'
   !debug: close(1) 
    deallocate(b01,b02,b04)

    deallocate(ResidusRec,Residusdc,Rrec,Nrec,Rdc,Ndc)!,Rdc2,Ndc2)
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Deallocated (1)'
    !close(1) 
    deallocate(vuu)
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Deallocated (2)'
    !write(1,*)'multiveJoint.f90:: nig(1:3) = ', nig(1:3)
    !close(1) 
    deallocate(nig)
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: cdc(1:3) = ', cdc(1:3)
    !close(1) 
    deallocate(cdc)
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: cdc2(1:3) = ', cdc2(1:3)
    !close(1) 
    deallocate(cdc2)
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Deallocated (3)'
    !close(1) 
    deallocate(t0dc,t1dc)
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Deallocated (4)'
    !close(1) 
    deallocate(aux1,aux2)
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Deallocated (5)'
    !close(1) 
    deallocate(res1,res4,res3)
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Deallocated (6)'
    !close(1) 
    deallocate(mi)
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Deallocated (7)'
    !close(1) 
    deallocate(t0,t1,c)
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Deallocated (8)'
    !close(1) 
    deallocate(g,aux)
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Deallocated (9)'
    !close(1) 
    deallocate(vax,vaxdc,vaxdc2)
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Deallocated (10)'
    !close(1) 
    deallocate(ve,vedc,vedc2)
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Deallocated (11)'
    !close(1) 
    deallocate(ve1,ve2,ve4)
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Deallocated (12)'
    !close(1) 
    deallocate(filtre,filtre2,filtre4)
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Deallocated (13)'
    !close(1) 
    deallocate(Hspl_hess,PEN_deri,hess,v,I1_hess,H1_hess,&
    I2_hess,H2_hess,HI2,HIH,IH,HI,BIAIS)
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Deallocated (14)'
    !close(1) 
    deallocate(date,datedc)
    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Deallocating 2nd recurrent event stuff'
    !close(1) 
    if(event2_ind0.eq.1)then
        deallocate(t0meta,t1meta,ResidusRec2,cmeta,gmeta,auxmeta,vemeta,&
        res1meta,res3meta,ve3,nigmeta,vaxmeta,filtre3,datemeta,b03,Rrec2,Nrec2)
    endif

    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*)'multiveJoint.f90:: Deallocating spline variables.'
    !close(1) 

    if(typeof.eq.0)then
        deallocate(nt0dc,nt1dc,nt0,nt1,mm3dc,mm2dc,mm1dc,mmdc,im3dc,im2dc,im1dc,imdc,mm3,mm2,&
        mm1,mm,im3,im2,im1,im,zi,zidc,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,&
        m2m,m1m)
        deallocate(m3m3b,m2m2b,m1m1b,mmmb,m3m2b,m3m1b,m3mb,m2m1b,m2mb,m1mb,m3m3d,m2m2d,m1m1d,&
        mmmd,m3m2d,m3m1d,m3md,m2m1d,m2md,m1md)
        if(event2_ind0.eq.1)then
            deallocate(zimeta,nt0meta,nt1meta,mm3meta,immeta,mm2meta,mm1meta,mmmeta,im3meta,im2meta,im1meta,&
            m3m3c,m2m2c,m1m1c,mmmc,m3m2c,m3m1c,m3mc,m2m1c,m2mc,m1mc)
        endif
    endif

    if(typeof.ne.0)then
        deallocate(vvv)
    endif

    !Only relevant for piecewise constant hazards
    !if(typeof.eq.1)then
        !deallocate(ttt,tttdc,betacoef)
        !if(event2_ind0.eq.1)then
        !    deallocate(tttmeta)
        !endif
    !endif

   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: All Variables Deallocated.'
   !debug: close(1) 

    !cptEvent(1)=cpt
    !cptEvent(2)=cpt_dc
    !if(event2_ind0.eq.1)then
    !    cptEvent(3)=cptmeta
    !endif
    !cptEvent(4)=cpt_dc2

    !ResMartingaleEvent(1:nobsEvent(3),1)=Res_martingale(1:nobsEvent(3))
    !ResMartingaleEvent(1:nobsEvent(3),2)=Res_martingaledc(1:nobsEvent(3))
    !ResMartingaleEvent(1:nobsEvent(3),3)=Res_martingale2(1:nobsEvent(3))

    !frailtyEstimates(1:nobsEvent(3),1)=frailtypred(1:nobsEvent(3))
    !frailtyEstimates(1:nobsEvent(3),2)=frailtypred2(1:nobsEvent(3))
    !frailtyEstimates(1:nobsEvent(3),3)=frailtyvar(1:nobsEvent(3))
    !frailtyEstimates(1:nobsEvent(3),4)=frailtyvar2(1:nobsEvent(3))
    !frailtyEstimates(1:nobsEvent(3),5)=frailtyCorr(1:nobsEvent(3))


   !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
   !debug: write(1,*)'multiveJoint.f90:: complete....'
   !debug: close(1) 

    return

    end subroutine joint_multiv
!-----------------------------------------------------------------------










!-----------------------------------------------------------------------
! Supporting Functions

!========================== VECSPLI ==============================
! LE: This needed to be changed for the competing joint models
    subroutine vecspli(ndate,ndatedc)
    use taillesmultiv
    use comonmultiv,only:date,datedc,vectn,zidc,zi,mm3,mm2,mm1,mm,im3,im2,im1,im &
    ,mm3dc,mm2dc,mm1dc,mmdc,im3dc,im2dc,im1dc,imdc

    IMPLICIT NONE

    integer,intent(in)::ndate,ndatedc
    integer::i,j,k
    double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2

! Calculate hazard u(ti) for recurrent event 1
!----------  calcul de u(ti)  ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)

    j=0
    do i=1,ndate-1
        do k = 2,vectn-2
            if((date(i).ge.zi(k-1)).and.(date(i).lt.zi(k)))then
                j = k-1
            endif
        end do 
        ht = date(i)-zi(j)
        htm= date(i)-zi(j-1)
        h2t= date(i)-zi(j+2)
        ht2 = zi(j+1)-date(i)
        ht3 = zi(j+3)-date(i)
        hht = date(i)-zi(j-2)
        h = zi(j+1)-zi(j)
        hh= zi(j+1)-zi(j-1)
        h2= zi(j+2)-zi(j)
        h3= zi(j+3)-zi(j)
        h4= zi(j+4)-zi(j)
        h3m= zi(j+3)-zi(j-1)
        h2n=zi(j+2)-zi(j-1)
        hn= zi(j+1)-zi(j-2)
        hh3 = zi(j+1)-zi(j-3)
        hh2 = zi(j+2)-zi(j-2)
        mm3(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
        mm2(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
              *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
        mm1(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
              h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
        mm(i)  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
        im3(i) = (0.25d0*(date(i)-zi(j-3))*mm3(i))+(0.25d0*hh2 &
            *mm2(i))+(0.25d0*h3m*mm1(i))+(0.25d0*h4*mm(i))
        im2(i) = (0.25d0*hht*mm2(i))+(h3m*mm1(i)*0.25d0) &
            +(h4*mm(i)*0.25d0)
        im1(i) = (htm*mm1(i)*0.25d0)+(h4*mm(i)*0.25d0)
        im(i)  = ht*mm(i)*0.25d0

    end do
!AD: add for death 
!----------  calcul de u(ti) :  STRATE2 ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)
    j=0
    do i=1,ndatedc-1
        do k = 2,vectn-2
            ! Which nodes does the date fall between
            ! Added -0.0001 to account for numerical error
            if((datedc(i).ge.(zidc(k-1)-0.0001d0)).and.(datedc(i).lt.zidc(k)))then
                j = k-1
            endif
        end do 
        ht = datedc(i)-zidc(j)
        htm= datedc(i)-zidc(j-1)
        h2t= datedc(i)-zidc(j+2)
        ht2 = zidc(j+1)-datedc(i)
        ht3 = zidc(j+3)-datedc(i)
        hht = datedc(i)-zidc(j-2)
        h = zidc(j+1)-zidc(j)
        hh= zidc(j+1)-zidc(j-1)
        h2= zidc(j+2)-zidc(j)
        h3= zidc(j+3)-zidc(j)
        h4= zidc(j+4)-zidc(j)
        h3m= zidc(j+3)-zidc(j-1)
        h2n=zidc(j+2)-zidc(j-1)
        hn= zidc(j+1)-zidc(j-2)
        hh3 = zidc(j+1)-zidc(j-3)
        hh2 = zidc(j+2)-zidc(j-2)

        mm3dc(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
        mm2dc(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
              *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
        mm1dc(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
              h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
        mmdc(i)  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
        im3dc(i) = (0.25d0*(datedc(i)-zi(j-3))*mm3dc(i))+(0.25d0*hh2 &
              *mm2dc(i))+(0.25d0*h3m*mm1dc(i))+(0.25d0*h4*mmdc(i))
        im2dc(i) = (0.25d0*hht*mm2dc(i))+(h3m*mm1dc(i)*0.25d0) &
              +(h4*mmdc(i)*0.25d0)
        im1dc(i) = (htm*mm1dc(i)*0.25d0)+(h4*mmdc(i)*0.25d0)
        imdc(i)  = ht*mmdc(i)*0.25d0

!        open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
!        write(1,*) "vecspli:: ndatedc = ", ndatedc
!        write(1,*) "vecspli:: datedc(i) = ", datedc(i)
!        write(1,*) "vecspli:: j = ", j
!        write(1,*) "vecspli:: vectn(2) = ", vectn(2)
!        write(1,*) "vecspli:: zidc(j) = ", zidc(j)
!        write(1,*) "vecspli:: ht = ", ht
!        write(1,*) "vecspli:: htm = ", htm
!        write(1,*) "vecspli:: h2t = ", h2t
!        write(1,*) "vecspli:: ht2 = ", ht2
!        write(1,*) "vecspli:: ht3 = ", ht3
!        write(1,*) "vecspli:: hht = ", hht
!        write(1,*) "vecspli:: h = ", h
!        write(1,*) "vecspli:: hh = ", hh
!        write(1,*) "vecspli:: h2 = ", h2
!        write(1,*) "vecspli:: h3 = ", h3
!        write(1,*) "vecspli:: h4 = ", h4
!        write(1,*) "vecspli:: h3m = ", h3m
!        write(1,*) "vecspli:: h2n = ", h2n
!        write(1,*) "vecspli:: hn = ", hn
!        write(1,*) "vecspli:: hh3 = ", hh3
!        write(1,*) "vecspli:: hh2 = ", hh2
!        write(1,*) "vecspli:: mm3dc(i) = ", mm3dc(i)
!        write(1,*) "vecspli:: mm2dc(i) = ", mm2dc(i)
!        write(1,*) "vecspli:: mm1dc(i) = ", mm1dc(i)
!        write(1,*) "vecspli:: mmdc(i) = ", mmdc(i)
!        write(1,*) "vecspli:: im3dc(i) = ", im3dc(i)
!        write(1,*) "vecspli:: im2dc(i) = ", im2dc(i)
!        write(1,*) "vecspli:: im1dc(i) = ", im1dc(i)
!        write(1,*) "vecspli:: imdc(i) = ", imdc(i)
!        close(1)
    end do


!AD:end
!    j=0
!    do i=1,ndatemeta-1
!        do k = 2,vectn(3)-2
!            if((datemeta(i).ge.zimeta(k-1)).and.(datemeta(i).lt.zimeta(k)))then
!                j = k-1
!            endif
!        end do 
!        ht = datemeta(i)-zimeta(j)
!        htm= datemeta(i)-zimeta(j-1)
!        h2t= datemeta(i)-zimeta(j+2)
!        ht2 = zimeta(j+1)-datemeta(i)
!        ht3 = zimeta(j+3)-datemeta(i)
!        hht = datemeta(i)-zimeta(j-2)
!        h = zimeta(j+1)-zimeta(j)
!        hh= zimeta(j+1)-zimeta(j-1)
!        h2= zimeta(j+2)-zimeta(j)
!        h3= zimeta(j+3)-zimeta(j)
!        h4= zimeta(j+4)-zimeta(j)
!        h3m= zimeta(j+3)-zimeta(j-1)
!        h2n=zimeta(j+2)-zimeta(j-1)
!        hn= zimeta(j+1)-zimeta(j-2)
!        hh3 = zimeta(j+1)-zimeta(j-3)
!        hh2 = zimeta(j+2)-zimeta(j-2)
!        mm3meta(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
!        mm2meta(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
!        *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n)) 
!        mm1meta(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
!        h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
!        mmmeta(i)  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
!        im3meta(i) = (0.25d0*(datemeta(i)-zimeta(j-3))*mm3meta(i))+(0.25d0*hh2  &
!        *mm2meta(i))+(0.25d0*h3m*mm1meta(i))+(0.25d0*h4*mmmeta(i))
!        im2meta(i) = (0.25d0*hht*mm2meta(i))+(h3m*mm1meta(i)*0.25d0)+(h4*mmmeta(i)*0.25d0)
!        im1meta(i) = (htm*mm1meta(i)*0.25d0)+(h4*mmmeta(i)*0.25d0)
!        immeta(i)  = ht*mmmeta(i)*0.25d0
!    end do

    end subroutine vecspli

!--------------------------------------------------------------------------

    subroutine vecpenP(n,zi,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m) 

    IMPLICIT NONE

    integer,intent(in)::n
    integer::i
    double precision::h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2,a3,a2,b2 &
    ,c2,a1,b1,c1,a0,x3,x2,x
    double precision,dimension(-2:(n+1))::zi
    double precision,dimension(n-2),intent(inout)::m3m3,m2m2,m1m1,mmm,m3m2, &
    m3m1,m3m,m2m1,m2m,m1m

!*********************************************************************
        m3m3=0.d0
    m2m2=0.d0
    m1m1=0.d0
    mmm=0.d0
    m3m2=0.d0
    m3m1=0.d0
    m3m=0.d0
    m2m1=0.d0
    m2m=0.d0
    m1m=0.d0
    do i=1,n-3
        h = zi(i+1)-zi(i)

        hh= zi(i+1)-zi(i-1)
        h2= zi(i+2)-zi(i)
        h3= zi(i+3)-zi(i)
        h4= zi(i+4)-zi(i)
        h3m= zi(i+3)-zi(i-1)
        h2n=zi(i+2)-zi(i-1)
        hn= zi(i+1)-zi(i-2)
        hh3 = zi(i+1)-zi(i-3)
        hh2 = zi(i+2)-zi(i-2)
        a3 = h*hh*hn*hh3
        a2 = hh2*hh*h*hn
        b2 = hh2*h2n*hh*h
        c2 = hh2*h2*h*h2n
        a1 = h3m*h2n*hh*h
        b1 = h3m*h2*h*h2n
        c1 = h3m*h3*h2*h
        a0 = h4*h3*h2*h
        x3 = zi(i+1)*zi(i+1)*zi(i+1)-zi(i)*zi(i)*zi(i)
        x2 = zi(i+1)*zi(i+1)-zi(i)*zi(i)
        x  = zi(i+1)-zi(i)

        m3m3(i) = (192.d0*h/(hh*hn*hh3*hh*hn*hh3))
        m2m2(i) = 64.d0*(((3.d0*x3-(3.d0*x2*(2.d0*zi(i+1)+zi(i-2) &
        ))+x*(4.d0*zi(i+1)*zi(i+1)+zi(i-2)*zi(i-2)+4.d0*zi(i+1) &
        *zi(i-2)))/(a2*a2)))
        m2m2(i) = m2m2(i) + 64.d0*(((3.d0*x3-(3.d0*x2*(zi(i+2)  &
        +zi(i-1)+zi(i+1)))+x*(zi(i+2)*zi(i+2)+zi(i-1)*zi(i-1) &
        +zi(i+1)*zi(i+1)+2.d0*zi(i+2)*zi(i-1)+2.d0*zi(i+2) &
        *zi(i+1)+2.d0*zi(i-1)*zi(i+1)))/(b2*b2)))
        m2m2(i) = m2m2(i) +64.d0*((3.d0*x3-(3.d0*x2*(2.d0*zi(i+2) &
        +zi(i)))+x*(4.d0*zi(i+2)*zi(i+2)+zi(i)*zi(i)+4.d0*zi(i+2) &
        *zi(i)))/(c2*c2))

        m2m2(i) = m2m2(i) +128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+2) &
        +zi(i-1)+3.d0*zi(i+1)+zi(i-2)))+x*(2.d0*zi(i+1)*zi(i+2) &
        +2.d0*zi(i+1)*zi(i-1)+2.d0*zi(i+1)*zi(i+1)+zi(i-2)*zi(i+2) &
        +zi(i-2)*zi(i-1)+zi(i-2)*zi(i+1)))/(a2*b2))
        m2m2(i) = m2m2(i) + 128.d0*((3.d0*x3-(1.5d0* & 
        x2*(2.d0*zi(i+2)+zi(i)+2.d0*zi(i+1)+zi(i-2)))+x* &
        (4.d0*zi(i+1)*zi(i+2)+2.d0*zi(i+1)*zi(i)+2.d0*zi(i-2) &
        *zi(i+2)+zi(i-2)*zi(i)))/(a2*c2))
        m2m2(i) = m2m2(i) + 128.d0*((3.d0*x3-(1.5d0*x2 &
        *(3.d0*zi(i+2)+zi(i)+zi(i-1)+zi(i+1)))+x*(zi(i+2)*zi(i)+ &
        2.d0*zi(i-1)*zi(i+2)+zi(i)*zi(i-1)+2.d0*zi(i+1)*zi(i+2) &
        +zi(i+1)*zi(i)+2.d0*zi(i+2)*zi(i+2)))/(b2*c2))
        m1m1(i) = 64.d0*((3.d0*x3-(3.d0*x2*(2.d0*zi(i-1)+zi(i+1))) &
        +x*(4.d0*zi(i-1)*zi(i-1)+zi(i+1)*zi(i+1)+4.d0*zi(i-1) &
        *zi(i+1)))/(a1*a1))
        m1m1(i) = m1m1(i) + 64.d0*((3.d0*x3-(3.d0*x2*(zi(i-1)+zi(i) &     
        +zi(i+2)))+x*(zi(i-1)*zi(i-1)+zi(i)*zi(i)+zi(i+2)* &
        zi(i+2)+2.d0*zi(i-1)*zi(i)+2.d0*zi(i-1)*zi(i+2)+2.d0* &
        zi(i)*zi(i+2)))/(b1*b1))
        m1m1(i) = m1m1(i) + 64.d0*((3.d0*x3-(3.d0*x2*(zi(i+3) &
        +2.d0*zi(i)))+x*(zi(i+3)*zi(i+3)+4.d0*zi(i)*zi(i) &
        +4.d0*zi(i+3)*zi(i)))/(c1*c1)) 
        m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(3.d0 &
        *zi(i-1)+zi(i)+zi(i+2)+zi(i+1)))+x*(2.d0*zi(i-1)*zi(i-1) &
        +2.d0*zi(i-1)*zi(i)+2.d0*zi(i-1)*zi(i+2)+zi(i+1)*zi(i-1) &
        +zi(i+1)*zi(i)+zi(i+1)*zi(i+2)))/(a1*b1))
        m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+3)+ &
        2.d0*zi(i)+2.d0*zi(i-1)+zi(i+1)))+x*(2.d0*zi(i-1)*zi(i+3) &
        +4.d0*zi(i-1)*zi(i)+zi(i+1)*zi(i+3)+2.d0*zi(i+1)*zi(i))) &
        /(a1*c1))    
        m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+3)+3.d0 &
        *zi(i)+zi(i-1)+zi(i+2)))+x*(zi(i-1)*zi(i+3)+2.d0*zi(i-1) &    
        *zi(i)+zi(i+3)*zi(i)+2.d0*zi(i)*zi(i)+zi(i+2)*zi(i+3) &
        +2.d0*zi(i+2)*zi(i)))/(b1*c1))
        mmm(i) = (192.d0*h/(h4*h3*h2*h4*h3*h2))
        m3m2(i) = 192.d0*(((-x3+(0.5d0*x2*(5.d0*zi(i+1)+zi(i-2) &
        ))-x*(2.d0*zi(i+1)*zi(i+1)+zi(i+1)*zi(i-2)))/(a3*a2)) &
        +((-x3+(0.5d0*x2*(4.d0*zi(i+1)+zi(i-1)+zi(i+2)))-x* &
        (zi(i+1)*zi(i+2)+zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(a3*b2)) &
        +((-x3+(0.5d0*x2*(3.d0*zi(i+1)+2.d0*zi(i+2)+zi(i)))-x* &
        (2.d0*zi(i+1)*zi(i+2)+zi(i+1)*zi(i)))/(a3*c2)) )
        m3m1(i) = 192.d0*(((x3-(0.5d0*x2*(4.d0*zi(i+1)+2.d0*zi(i-1) &
        ))+x*(2.d0*zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(a3*a1)) &
        +((x3-(0.5d0*x2*(3.d0*zi(i+1)+zi(i+2)+zi(i-1)+zi(i))) &
        +x*(zi(i+1)*zi(i-1)+zi(i+1)*zi(i)+zi(i+1)*zi(i+2)))/(b1*a3)) &
        +((x3-(0.5d0*x2*(3.d0*zi(i+1)+zi(i+3)+2.d0*zi(i)))+x*(zi(i+1) &
        *zi(i+3)+2.d0*zi(i+1)*zi(i)))/(c1*a3)) )
        m3m(i) = 576.d0*((-(x3/3.d0)+(0.5d0*x2*(zi(i+1)+zi(i))) &
        -x*zi(i+1)*zi(i))/(a3*a0))
        m2m1(i) = 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0*zi(i-1)+3.d0* &
        zi(i+1)+zi(i-2)))-x*(4.d0*zi(i+1)*zi(i-1)+2.d0*zi(i+1) &
        *zi(i+1)+2.d0*zi(i-2)*zi(i-1)+zi(i-2)*zi(i+1)))/(a2*a1))
        m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i-1)+ &
        zi(i)+zi(i+2)+2.d0*zi(i+1)+zi(i-2)))-x*(2.d0*zi(i+1)*zi(i-1) &
        +2.d0*zi(i+1)*zi(i)+2.d0*zi(i+1)*zi(i+2)+zi(i-2)*zi(i-1)+ &
        zi(i-2)*zi(i)+zi(i-2)*zi(i+2)))/(a2*b1))
        m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3)+2.d0 &
        *zi(i)+2.d0*zi(i+1)+zi(i-2)))-x*(2.d0*zi(i+1)*zi(i+3)+4.d0 &
        *zi(i+1)*zi(i)+zi(i-2)*zi(i+3)+2.d0*zi(i-2)*zi(i)))/(a2*c1))
        m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2* &
        (3.d0*zi(i-1)+2.d0*zi(i+1)+zi(i+2)))-x*(2.d0*zi(i+2)*zi(i-1) &
        +zi(i+2)*zi(i+1)+2.d0*zi(i-1)*zi(i-1)+3.d0 &
        *zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(b2*a1))
        m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0 &
        *zi(i-1)+zi(i)+2.d0*zi(i+2)+zi(i+1)))-x*(zi(i+2)*zi(i-1) &
        +zi(i+2)*zi(i)+zi(i+2)*zi(i+2)+zi(i-1)*zi(i-1)+zi(i-1) &
        *zi(i)+zi(i-1)*zi(i+2)+zi(i+1)*zi(i-1)+zi(i+1)*zi(i) &
        +zi(i+1)*zi(i+2)))/(b2*b1))
        m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3) &
        +2.d0*zi(i)+zi(i+2)+zi(i-1)+zi(i+1)))-x*(zi(i+2)*zi(i+3) &
        +2.d0*zi(i+2)*zi(i)+zi(i-1)*zi(i+3)+2.d0*zi(i-1)*zi(i) &
        +zi(i+1)*zi(i+3)+2.d0*zi(i+1)*zi(i)))/(b2*c1))
        m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0*zi(i-1) &
        +zi(i+1)+2.d0*zi(i+2)+zi(i)))-x*(4.d0*zi(i+2)*zi(i-1)+2.d0* &
        zi(i+2)*zi(i+1)+2.d0*zi(i)*zi(i-1)+zi(i)*zi(i+1)))/(c2*a1))
        m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i-1) &
        +2.d0*zi(i)+3.d0*zi(i+2)))-x*(2.d0*zi(i+2)*zi(i-1)+2.d0 &
        *zi(i+2)*zi(i)+2.d0*zi(i+2)*zi(i+2)+zi(i)*zi(i-1)+zi(i) &
        *zi(i)+zi(i)*zi(i+2)))/(c2*b1))
        m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3) &
        +3.d0*zi(i)+2.d0*zi(i+2)))-x*(2.d0*zi(i+2)*zi(i+3)+4.d0 &
        *zi(i+2)*zi(i)+zi(i)*zi(i+3)+2.d0*zi(i)*zi(i)))/(c2*c1))
        m2m(i) = 192.d0*(((x3-(0.5d0*x2*(3.d0*zi(i)+2.d0*zi(i+1) &
        +zi(i-2)))+x*(2.d0*zi(i+1)*zi(i)+zi(i-2)*zi(i)))/(a2*a0)) &
        +((x3-(0.5d0*x2*(3.d0*zi(i)+zi(i+2)+zi(i-1)+zi(i+1))) &
        +x*(zi(i+2)*zi(i)+zi(i-1)*zi(i)+zi(i+1)*zi(i)))/(b2*a0)) &
        +((x3-(0.5d0*x2*(4.d0*zi(i)+2.d0*zi(i+2)))+x*(2.d0*zi(i+2) &
        *zi(i)+zi(i)*zi(i)))/(c2*a0)) )
        m1m(i) = 192.d0*(((-x3+(0.5d0*x2*(3.d0*zi(i)+2.d0*zi(i-1) &
        +zi(i+1)))-x*(2.d0*zi(i-1)*zi(i)+zi(i+1)*zi(i)))/(a1*a0)) &
        +((-x3+(0.5d0*x2*(4.d0*zi(i)+zi(i-1)+zi(i+2))) &
        -x*(zi(i-1)*zi(i)+zi(i)*zi(i)+zi(i+2)*zi(i)))/(b1*a0)) &
        +((-x3+(0.5d0*x2*(5.d0*zi(i)+zi(i+3)))-x*(zi(i+3)*zi(i) &
        +2.d0*zi(i)*zi(i)))/(c1*a0)) )

    end do

    end subroutine vecpenP




!==========================  SUSP  ====================================
    subroutine susp(x,the,n,su,lam,zi)

    use taillesmultiv

    IMPLICIT NONE 

    integer,intent(in)::n
    double precision,intent(out)::lam,su
    double precision,dimension(-2:npmax),intent(in)::zi,the
    double precision,intent(in)::x
    integer::j,k,i
    double precision::ht,ht2,h2,som,htm,h2t,h3,h2n,hn, &
    im,im1,im2,mm1,mm3,ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2 &
    ,h,gl,hh

    gl=0.d0
    som = 0.d0
    do k = 2,n+1
        if((x.ge.zi(k-1)).and.(x.lt.zi(k)))then
            j = k-1
            if(j.gt.1)then
                do i=2,j
                    som = som+the(i-4)
                end do
            endif
            ht = x-zi(j)
            htm= x-zi(j-1)
            h2t= x-zi(j+2)
            ht2 = zi(j+1)-x
            ht3 = zi(j+3)-x
            hht = x-zi(j-2)
            h = zi(j+1)-zi(j)
            hh= zi(j+1)-zi(j-1)
            h2= zi(j+2)-zi(j)
            h3= zi(j+3)-zi(j)
            h4= zi(j+4)-zi(j)
            h3m= zi(j+3)-zi(j-1)
            h2n=zi(j+2)-zi(j-1)
            hn= zi(j+1)-zi(j-2)
            hh3 = zi(j+1)-zi(j-3)
            hh2 = zi(j+2)-zi(j-2)
            mm3 = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
            mm2 = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
            *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
            mm1 = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
            h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
            mm  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
            im3 = (0.25d0*(x-zi(j-3))*mm3)+(0.25d0*hh2*mm2) &
            +(0.25d0*h3m*mm1)+(0.25d0*h4*mm)
            im2 = (0.25d0*hht*mm2)+(h3m*mm1*0.25d0)+(h4*mm*0.25d0)
            im1 = (htm*mm1*0.25d0)+(h4*mm*0.25d0)
            im  = ht*mm*0.25d0
            gl = som +(the(j-3)*im3)+(the(j-2)*im2)+(the(j-1)*im1)+(the(j)*im)
            lam = (the(j-3)*mm3)+(the(j-2)*mm2)+(the(j-1)*mm1)+(the(j)*mm)
        endif
    end do

    if(x.ge.zi(n))then
        som = 0.d0
        do i=1,n+1
            som = som+the(i-3)
        end do
        gl = som
    endif

    su  = dexp(-gl)

    return

    end subroutine susp

!==========================  COSP  ====================================
! calcul les points pour les fonctions
! et leur bandes de confiance

    subroutine cosp(x,the,n,y,zi,binf,su,bsup,lbinf,lam,lbsup)

    use taillesmultiv

    IMPLICIT NONE

    integer,intent(in)::n
    double precision,intent(in)::x
    double precision,intent(out)::lam,su
    double precision,intent(out)::binf,bsup,lbinf,lbsup
    double precision,dimension(npmax,npmax),intent(in)::y
    double precision,dimension(-2:npmax),intent(in)::the,zi
    integer::j,k,i
    double precision::ht,ht2,h2,som,pm,htm,h2t,h3,h2n,hn, &
    im,im1,im2,mm1,mm3,ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2, &
    h,gl,hh

    j=0
    gl=0.d0
    som = 0.d0
    do k = 2,n-1
        if((x.ge.zi(k-1)).and.(x.lt.zi(k)))then
            j = k-1
            if(j.gt.1)then
                do i=2,j
                som = som+the(i-4)
                end do  
            endif  
            ht = x-zi(j)
            htm= x-zi(j-1)
            h2t= x-zi(j+2)
            ht2 = zi(j+1)-x
            ht3 = zi(j+3)-x
            hht = x-zi(j-2)
            h = zi(j+1)-zi(j)
            hh= zi(j+1)-zi(j-1)
            h2= zi(j+2)-zi(j)
            h3= zi(j+3)-zi(j)
            h4= zi(j+4)-zi(j)
            h3m= zi(j+3)-zi(j-1)
            h2n=zi(j+2)-zi(j-1)
            hn= zi(j+1)-zi(j-2)
            hh3 = zi(j+1)-zi(j-3)
            hh2 = zi(j+2)-zi(j-2)
            mm3 = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
            mm2 = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
            *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
            mm1 = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
            h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
            mm  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
            im3 = (0.25d0*(x-zi(j-3))*mm3)+(0.25d0*hh2*mm2) &
            +(0.25d0*h3m*mm1)+(0.25d0*h4*mm)
            im2 = (0.25d0*hht*mm2)+(h3m*mm1*0.25d0)+(h4*mm*0.25d0)
            im1 = (htm*mm1*0.25d0)+(h4*mm*0.25d0)
            im  = ht*mm*0.25d0
            gl = som +(the(j-3)*im3)+(the(j-2)*im2)+(the(j-1)*im1)+(the(j)*im)
            lam = (the(j-3)*mm3)+(the(j-2)*mm2)+(the(j-1)*mm1)+(the(j)*mm)
        endif
    end do

    if(x.ge.zi(n))then
        som = 0.d0
        do i=1,n
            som = som+the(i-3)
        end do
        gl = som
    endif

    call conf(x,j,n,y,pm,zi)

    binf = dexp(-gl + 1.96d0*pm)
    su  = dexp(-gl)
    bsup = dexp(-gl - 1.96d0*pm)

    call conf1(x,j,n,y,pm,zi)
    lbinf = lam - 1.96d0*pm
    lbsup = lam + 1.96d0*pm
!         write(*,*)'lbinf apres conf1',lbinf,lam,pm

    return

    end subroutine cosp


!=====================  CONF1  =============================


    subroutine  conf1(x,ni,n,y,pm,zi)

    use taillesmultiv

    IMPLICIT NONE

    integer,intent(in)::ni,n
    double precision,intent(in)::x
    double precision,dimension(-2:npmax),intent(in)::zi
    double precision,dimension(npmax,npmax),intent(in)::y
    double precision,intent(out)::pm
    integer::i,j
    double precision::res,mmsp
    double precision,dimension(npmax)::vecti,aux

    do i=1,n
        vecti(i) = mmsp(x,ni,i,zi)
    end do

    do i=1,n
        aux(i) = 0.d0
        do j=1,n
            aux(i) = aux(i) - y(i,j)*vecti(j)
        end do
    end do

    res = 0.d0
    do i=1,n
        res = res + aux(i)*vecti(i)
    end do

    res = -res
    pm = dsqrt(res)

    end subroutine  conf1

!=====================  CONF  =============================

    subroutine  conf(x,ni,n,y,pm,zi)

    use taillesmultiv

    IMPLICIT NONE

    integer,intent(in)::ni,n
    double precision,intent(in)::x
    double precision,dimension(-2:npmax),intent(in)::zi
    double precision,dimension(npmax,npmax),intent(in)::y
    double precision,intent(out)::pm
    integer::i,j
    double precision::res,isp
    double precision,dimension(52)::vecti,aux

    do i=1,n
        vecti(i) = isp(x,ni,i,zi)
    end do

    do i=1,n
    aux(i) = 0.d0
    do j=1,n
        aux(i) = aux(i) - y(i,j)*vecti(j)
    end do
    end do

    res = 0.d0
    do i=1,n
    res = res + aux(i)*vecti(i)
    end do
    res=-res
    pm = dsqrt(res)

    end subroutine  conf


!==========================   ISP   ==================================

    double precision function isp(x,ni,ns,zi)

    use taillesmultiv

    IMPLICIT NONE

    integer,intent(in)::ni,ns
    double precision,intent(in)::x
    double precision,dimension(-2:npmax),intent(in)::zi
    double precision::val,mmsp



    if(x.eq.zi(ni))then
        if(ni.le.ns-3)then
            val = 0.d0
            else
                if(ni.le.ns-2)then
                    val = ((zi(ni)-zi(ni-1))*mmsp(x,ni,ns,zi))*0.25d0
                else
                    if(ni.eq.ns-1)then
                        val = ((zi(ni)-zi(ni-2))*mmsp(x,ni,ns,zi)+ &
                        (zi(ni+3)-zi(ni-1))*mmsp(x,ni,ns+1,zi))*0.25d0
                    else
                        if(ni.eq.ns)then
                            val = ((zi(ni)-zi(ni-3))*mmsp(x,ni,ns,zi)+ &
                            (zi(ni+2)-zi(ni-2))*mmsp(x,ni,ns+1,zi) &
                            +(zi(ni+3)-zi(ni-1))*mmsp(x,ni,ns+2,zi))*0.25d0
                            else
                                val = 1.d0
                            endif
                        endif
                endif
        endif
    else
        if(ni.lt.ns-3)then
            val = 0.d0
        else
            if(ni.eq.ns-3)then
                val = (x-zi(ni))*mmsp(x,ni,ns,zi)*0.25d0
            else
                if(ni.eq.ns-2)then
                    val = ((x-zi(ni-1))*mmsp(x,ni,ns,zi)+ &
                    (zi(ni+4)-zi(ni))*mmsp(x,ni,ns+1,zi))*0.25d0
                else
                    if(ni.eq.ns-1)then
                        val =((x-zi(ni-2))*mmsp(x,ni,ns,zi)+ &
                        (zi(ni+3)-zi(ni-1))*mmsp(x,ni,ns+1,zi) &
                        +(zi(ni+4)-zi(ni))*mmsp(x,ni,ns+2,zi))*0.25d0
                    else
                        if(ni.eq.ns)then
                            val =((x-zi(ni-3))*mmsp(x,ni,ns,zi)+ &
                            (zi(ni+2)-zi(ni-2))*mmsp(x,ni,ns+1,zi) &
                            +(zi(ni+3)-zi(ni-1))*mmsp(x,ni,ns+2,zi) &
                            +(zi(ni+4)-zi(ni))*mmsp(x,ni,ns+3,zi))*0.25d0
                        else
                            val = 1.d0
                        endif
                    endif
                endif
            endif
        endif
    endif

    isp = val

    return

    end function isp

!==========================  MMSP   ==================================

    double precision function mmsp(x,ni,ns,zi)

    use taillesmultiv

    IMPLICIT NONE 

    integer,intent(in)::ni,ns
    double precision,intent(in)::x
    double precision,dimension(-2:npmax),intent(in)::zi
    double precision::val

    if(ni.lt.ns-3)then
        val = 0.d0
    else
        if(ns-3.eq.ni)then
            if(x.eq.zi(ni))then
                val = 0.d0
            else  
                val = (4.d0*(x-zi(ni))*(x-zi(ni)) &
                *(x-zi(ni)))/((zi(ni+4)-zi(ni))*(zi(ni+3) &
                -zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
            endif
        else 
            if(ns-2.eq.ni)then
                if(x.eq.zi(ni))then
                    val = (4.d0*(zi(ni)-zi(ni-1))*(zi(ni)-zi(ni-1))) &
                    /((zi(ni+3)-zi(ni-1))*(zi(ni+2)-zi(ni-1)) &
                    *(zi(ni+1)-zi(ni-1)))
                else  
                    val = (4.d0*(x-zi(ni-1))*(x-zi(ni-1)) &
                    *(zi(ni+1)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+2) &
                    -zi(ni-1))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)-zi(ni))) &
                    +   (4.d0*(x-zi(ni-1))*(x-zi(ni)) &
                    *(zi(ni+2)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+2) &
                    -zi(ni))*(zi(ni+1)-zi(ni))*(zi(ni+2)-zi(ni-1))) &
                    +   (4.d0*(x-zi(ni))*(x-zi(ni)) &
                    *(zi(ni+3)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+3) &
                    -zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
                endif
            else
                if(ns-1.eq.ni)then
                    if(x.eq.zi(ni))then
                        val = (4.d0*((zi(ni)-zi(ni-2))*(zi(ni+1) &
                        -zi(ni)))/((zi(ni+2)-zi(ni-2))*(zi(ni+1) &
                        -zi(ni-1))*(zi(ni+1)-zi(ni-2)))) &
                        +((4.d0*((zi(ni)-zi(ni-1))*(zi(ni+2)-zi(ni))) &
                        /((zi(ni+2)-zi(ni-2))*(zi(ni+2)-zi(ni-1)) &
                        *(zi(ni+1)-zi(ni-1)))))
                    else
                        val = (4.d0*((x-zi(ni-2))*(zi(ni+1) &
                        -x)*(zi(ni+1)-x))/((zi(ni+2) &
                        -zi(ni-2))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)- &
                        zi(ni))*(zi(ni+1)-zi(ni-2)))) &
                        +((4.d0*((x-zi(ni-1))*(zi(ni+2)-x)  &
                        *(zi(ni+1)-x))/((zi(ni+2)-zi(ni-2)) &
                        *(zi(ni+2)-zi(ni-1))*(zi(ni+1)-zi(ni-1))* &
                        (zi(ni+1)-zi(ni))))) &
                        +((4.d0*((zi(ni+2)-x)*(zi(ni+2)-x) &
                        *(x-zi(ni)))/((zi(ni+2)-zi(ni-2)) &
                        *(zi(ni+2)-zi(ni))*(zi(ni+2)-zi(ni-1))* &
                        (zi(ni+1)-zi(ni)))))
                    endif
                else
                    if(ni.eq.ns)then
                            if(x.eq.zi(ni))then
                            val =(4.d0*(x-zi(ni+1))*(x &
                            -zi(ni+1))/((zi(ni+1)-zi(ni-1))*(zi(ni+1) &
                            -zi(ni-2))*(zi(ni+1)-zi(ni-3))))
                        else   
                            val =(4.d0*(x-zi(ni+1))*(x &
                            -zi(ni+1))*(zi(ni+1)-x)/((zi(ni+1) &
                            -zi(ni-1))*(zi(ni+1)-zi(ni-2))*(zi(ni+1) &
                            -zi(ni))*(zi(ni+1)-zi(ni-3))))
                        endif
                    else
                        val = 0.d0
                    endif
                endif
            endif
        endif
    endif

    mmsp = val

    return

    end function mmsp

!-------------------------------------------------------
! Matrix Multiplecation 
! multiply A by B with the result in C

    subroutine multi(A,B,IrowA,JcolA,JcolB,C)
!     remarque :  jcolA=IrowB
    use taillesmultiv

    IMPLICIT NONE

    integer,intent(in)::IrowA,JcolA,JcolB
    double precision,dimension(npmax,npmax),intent(in):: A,B
    double precision,dimension(npmax,npmax),intent(out)::C
    integer::i,j,k
    double precision::sum

    do I=1,IrowA
        do J=1,JcolB
            sum=0
            do K=1,JcolA
                sum=sum+A(I,K)*B(K,J)
            end do
            C(I,J)=sum
        end do
    end do

    return

    end subroutine multi
!-------------------------------------------------------

    double precision  function func30(frail1,frail2)
! calcul de l integrant, pour un effet aleatoire donne frail et un groupe donne auxig (cf funcpa)

    use taillesmultiv
    use comonmultiv,only:nigmeta,alpha1,alpha2,res1meta,res3meta,&
    nig,auxig,alpha,theta,eta,aux1,res1,res3,cdc 

    implicit none

    double precision::frail1,frail2

    func30=0.d0
    func30 = frail1*(cdc(auxig)*alpha1+nig(auxig)) &
    +frail2*(cdc(auxig)*alpha2+nigmeta(auxig)) &
    -dexp(frail1)*(res1(auxig)-res3(auxig))-dexp(frail2)*(res1meta(auxig)-res3meta(auxig)) &
    -dexp(frail1*alpha1+frail2*alpha2)*aux1(auxig) &
    +(2.d0*((2.d0*dexp(alpha)/(dexp(alpha)+1.d0))-1.d0) &
    *frail1*frail2/sqrt(theta*eta)  &
    -(frail1**2.d0)/theta -(frail2**2.d0)/eta) &
    /(2.d0*(1.d0-((2.d0*dexp(alpha)/(dexp(alpha)+1.d0))-1.d0)**2.d0))
    func30 = dexp(func30)

    return

    end function func30

