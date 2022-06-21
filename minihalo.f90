program main
  implicit none
  integer i,j,k
  double precision z,dVdz,z_f,a_c,dz_f,dz,Trad
  double precision dndz
  ! double precision dnda_BBKS,n_BBKS,dndz_PS,n_PS,dndz
  double precision rho_over_M,R_UCMH,g_ff,f_cool,epsilon
  double precision T_vir,T,t_comp,t_ff,t_cool,t_life,x_e,n_b_sq,epsilon_nu,E_th
  double precision nu_obs(60),nIhalo(60),dIdz(60),zmin(60)
  double precision I_obs,TR,nu,B_obs,I_ind,f_sky

  integer,parameter::Nz_k=100000
  integer,parameter::Nz_j=10000
  integer,parameter::nz_switch=0
  ! 0 is Press-Schechter, 1 is Bardeen, Bond, Kaiser and Szalay
  integer,parameter::obs_switch=0
  ! 0 is 21-cm line, 1 is CMB foreground
  integer,parameter::gas_profile=1
  ! 0 is uniform gas density, 1 is Makino model.
  double precision,parameter::zmax=1000.0
  double precision,parameter::Mpc_over_m=3.08567758d22
  double precision,parameter::Mpc_over_cm=3.08567758d24
  double precision,parameter::eV_over_J=1.60218d-19
  double precision,parameter::Omega_m=0.3062
  double precision,parameter::Omega_b=0.0493
  double precision,parameter::Omega_mhsq=0.1408
  double precision,parameter::Omega_bhsq=0.0224
  double precision,parameter::Pi=3.141592653589d0
  double precision,parameter::c=2.99792d10        ! [cm/s]
  double precision,parameter::h_P=6.62607004d-34  ! [J s]
  double precision,parameter::kB=1.38064852d-23   ! [J/K]
  double precision,parameter::solarmass=1.9891d30 ! [kg]
  double precision,parameter::Jtoerg=1.0d7
  double precision,parameter::del_c=1.69d0        ! critical overdensity for spherical collapse
  ! double precision,parameter::del_c=1.686d0       ! critical overdensity for spherical collapse
  double precision,parameter::rho_c=2.775d11      ! rho_{crit,0} h^{-2} [h^2 M_sun / Mpc^3]
  double precision,parameter::m_H=1.672621898d-27 ! [kg]
  double precision,parameter::z_eq=3384.0         ! redshift for matter-radiation equality,
  double precision rho_crit_z

  ! model parameters:
  ! double precision,parameter::A_mat=5.0d2**2.0
  double precision,parameter::A_mat=(61.0d0/0.7817397046698262)**2.0
  ! double precision,parameter::A_mat=15421.430445460266
  ! double precision,parameter::A_mat=63.7d0**2.0
  ! Amplitude of (matter) power spectrum, A_mat = sigma_0^2 in Abe et al. (2021)
  double precision,parameter::M_minihalo=1.0d10   ! [M_sun]
  double precision k_spike, P_zeta, z_cut
  double precision,parameter::sigma_0(9)=(/0.55722,0.8382,1.17039,1.54698,1.96149,2.40840,2.88310,3.3818,3.90135/)
  double precision,parameter::A_s=2.10d-9
  double precision,parameter::n_s=0.965
  double precision A_mat_LCDM
  double precision outtest

  rho_over_M=Omega_mhsq*rho_c/M_minihalo

  ! ! minoda test:
  ! open(27,status='replace',form='formatted',file="concen_param.txt")
  ! ! do j=Nz_j,1,-1
  ! !    z=zmax*dble(j)/dble(Nz_j)
  ! !    outtest=rho_boc_sq_int(A_mat,z)
  ! ! end do
  ! outtest=rho_boc_sq_int(A_mat,4.977023564332112d0)
  ! outtest=rho_boc_sq_int(A_mat,10.476157527896651d0)
  ! outtest=rho_boc_sq_int(A_mat,20.09233002565048d0)
  ! outtest=rho_boc_sq_int(A_mat,50.9413801481638d0)
  ! outtest=rho_boc_sq_int(A_mat,107.22672220103243d0)
  ! outtest=rho_boc_sq_int(A_mat,205.65123083486534d0)

  ! ! do j=1,9
  ! !    outtest=rho_boc_sq_int(sigma_0(j)*sigma_0(j),7.0d0)
  ! ! end do
  ! close(27)
  ! stop


  ! print*,rho_boc_sq_int(A_mat,1000.0d0),(200.0*Omega_b/Omega_m)**2.0
  ! print*,rho_boc_sq_int(A_mat,9.0d0),(200.0*Omega_b/Omega_m)**2.0
  ! print*,rho_boc_sq_int(A_mat,1000.0d0),(200.0*Omega_b/Omega_m)
  ! print*,rho_boc_sq_int(A_mat,9.0d0),(200.0*Omega_b/Omega_m)
  ! stop

  open(5,status='replace',form='formatted',file="TRad.txt")
  open(15,status='replace',form='formatted',file="halo.txt")
  open(25,status='replace',form='formatted',file="Iobs_M12.txt")
  open(26,status='replace',form='formatted',file="Inu_int.txt")

  if(obs_switch.eq.0)then
     print*,"Calculation at 21-cm frequencies."
     do i=1,60
        ! input frequency
        nu_obs(i)=1.0d6*(2.0*i)
        zmin(i)=1.42d9/nu_obs(i)-1.0
        dIdz(i)=0.0
        print*,1.0d-6*nu_obs(i),"MHz",zmin(i)
     end do
  else if(obs_switch.eq.1)then
     print*,"Calculation at CMB frequencies."
     ! do j=Nz_j,1,-1
     !    z=zmax*dble(j)/dble(Nz_j)
     !    if(sqrt(A_mat)*D_plus(z)>del_c)exit
     ! end do
     k_spike=2.0*Pi/((M_minihalo/1.30d11)/(Omega_mhsq/0.112))**(1.0/3.0)
     print*,"k_s=",k_spike,"Mpc^-1"
     A_mat_LCDM=(log(4.0d-7*k_spike*c/sqrt(3.0*(1.0+z_eq)*Omega_mhsq))+0.577-3.5)**2.0
     !     A_mat_LCDM=81.0*(1.0+z_eq)*(1.0+z_eq)*A_s*(k_spike/0.05)**(n_s-1.0)*A_mat_LCDM
     A_mat_LCDM=81.0*(1.0+z_eq)*(1.0+z_eq)*A_s*A_mat_LCDM
     print*,"A_mat(LCDM)=",A_mat_LCDM
     if(sqrt(2.0*A_mat_LCDM)/del_c>1.0)then
        z_cut=dble(int(sqrt(2.0*A_mat_LCDM)/del_c))+1.0
        print*,"z_cut in double =",sqrt(2.0*A_mat_LCDM)/del_c
     else
        z_cut=0.0
     end if
     ! z_cut=z
     ! z_cut=10.0
     print*,"z_cut=",z_cut
     do i=1,60
        ! input frequency
        nu_obs(i)=7.0d6*10.0**(0.1*i)
        zmin(i)=z_cut
        ! zmin(i)=0.0
        ! minoda check
        dIdz(i)=0.0
     end do
     write(26,*)"### M_halo=",M_minihalo,"Msun,  sigma0=",sqrt(A_mat)
     write(26,*)"### nu_obs=",1.0e-9*nu_obs(50),"GHz,  z_cut=",z_cut
     write(26,*)"### z, fsky*n*Iind, dVdz, Iobs, r_com(z)"
  else
     print*,"Error: using undefined obs_switch=",obs_switch
     stop
  end if

  do j=Nz_j,1,-1
     z=zmax*dble(j)/dble(Nz_j)
     dz=z*(1.0-dble(j-1)/dble(j))

     if(mod(j,10).eq.1)then
        write(15,*)z,dndz_PS(1.0/(1.0+z),A_mat,M_minihalo),dndz_BBKS(1.0/(1.0+z),A_mat,M_minihalo)
     end if

     do i=1,60
        nIhalo(i)=0.0
     end do

     do k=Nz_k,1,-1
        z_f=z*exp(dble(k)/dble(Nz_k)*log(zmax/z))
        dz_f=z_f-z*exp(dble(k-1)/dble(Nz_k)*log(zmax/z))
        a_c=1.0/(1.0+z_f)
        R_UCMH = a_c*(3.0/800.0/Pi/rho_over_M)**(1.0/3.0)*Mpc_over_cm
        T_vir=1.98d4*(1.0+z_f)/10.0*(200.0/18.0/Pi/Pi*Omega_mhsq)**(1.0/3.0)
        T_vir=T_vir*(M_minihalo*1.0d-8)**(2.0/3.0)
        f_cool=0.5*(1.0+ionization(T_vir))/ionization(T_vir)
        t_life=2.0/3.0*Mpc_over_cm*(1.0/Hz(z)-1.0/Hz(z_f))
        t_comp=f_cool*4.2d14*(20.0/(1.0+z_f))**4.0 ! [sec]
        ! minoda added 2021/12/16

        ! cooling judgment
        if(t_life/t_comp>20.0)then
           do i=1,60
              nIhalo(i)=nIhalo(i)+0.0
           end do
        else
           x_e = ionization(T_vir)
           rho_crit_z=0.67*0.67*rho_c*(Omega_m/a_c/a_c/a_c+(1.0-Omega_m))
           n_b_sq=(rho_crit_z*solarmass/(Mpc_over_cm)**3.0/m_H)**2.0
           ! [warning] minoda changed for debug:
           n_b_sq=rho_boc_sq_int(A_mat,z_f)*n_b_sq
           ! n_b_sq=40000.0*n_b_sq


           epsilon=1.4d-27*n_b_sq*x_e*x_e*sqrt(T_vir)*exp(-h_P*nu/kB/T_vir)*1.2 ! velocity-averaged Gaunt-factor
           E_th = 1.5*rho_crit_z*solarmass/(Mpc_over_cm)**3.0/m_H*kB*T_vir*Jtoerg ! erg/cm^3
           t_ff = E_th/epsilon
           t_cool=min(t_comp,t_ff)


           if(t_life/t_cool>20.0)then
              do i=1,60
                 nIhalo(i)=nIhalo(i)+0.0
              end do
           else
              T=T_vir*exp(-1.0*t_life/t_cool)
              x_e = ionization(T)

              ! if(mod(k,100).eq.0)then
              !    print*,"minoda wrote: z_f=",z_f,"t_life=",t_life,"t_comp=",t_comp
              !    print*,"T_vir=",T_vir,"T=",T,"x=",x_e
              ! end if

              ! if(j.eq.1000.and.k<5)then
              !    print*,"z=",z,"z_f=",z_f
              !    print*,"M_minihalo=",M_minihalo
              !    print*,"T_vir=",T_vir,"T=",T,"x_e=",x_e
              ! end if
              if(nz_switch.eq.0)then
                 dndz=dndz_PS(a_c,A_mat,M_minihalo)
              else if(nz_switch.eq.1)then
                 dndz=dndz_BBKS(a_c,A_mat,M_minihalo)
              else
                 print*,"Error: using undefined nz_switch=",nz_switch
                 stop
              end if

              ! minoda: moved before cooling judgment
              ! rho_crit_z=0.67*0.67*rho_c*(Omega_m/a_c/a_c/a_c+(1.0-Omega_m))
              ! n_b_sq=(rho_crit_z*solarmass/(Mpc_over_cm)**3.0/m_H)**2.0
              ! n_b_sq=rho_boc_sq_int(A_mat,z_f)*n_b_sq


              do i=1,60
                 if (z<zmin(i))then
                    nIhalo(i)=nIhalo(i)+0.0
                 else
                    nu = (1.0+z_f)*nu_obs(i)
                    epsilon_nu=6.8d-38*n_b_sq*x_e*x_e/sqrt(T)*exp(-h_P*nu/kB/T)
                    g_ff=log(exp(5.96-sqrt(3.0)/Pi*log(1.0d-9*nu/(1.0e-4*T)**1.5))+exp(1.0))
                    ! g_ff=sqrt(3.0)/Pi*log(4.0/1.781*(kB*T)**1.5/h_P/nu/sqrt(13.6*eV_over_J))
                    epsilon_nu=epsilon_nu*g_ff
                    I_ind=min(4.0/3.0*R_UCMH*epsilon_nu,2.0*kB*T*(nu/c)**2.0*1.0d7)
                    f_sky=Pi*(R_UCMH/r(z)*(1+z)/Mpc_over_cm)**2.0/4.0/Pi
                    nIhalo(i)=nIhalo(i)+dz_f*dndz*I_ind*f_sky
                    if(i.eq.50.and.j.eq.100.and.k.eq.1)then
                    ! if(i.eq.50.and.j.eq.100.and.k.eq.1)then
                       print*,"M_minihalo=",M_minihalo
                       print*,"A_mat=",A_mat,"sigma0=",sqrt(A_mat)
                       print*,"zmin=",zmin(i)
                       print*,"z=",z,"z_f=",z_f
                       print*,"nu_obs=",nu_obs(i)

                       print*,"T_vir=",T_vir
                       print*,"t_life=",t_life,"t_cool=",t_cool
                       print*,"exp(-t_l/t_c)=",exp(-1.0*t_life/t_cool)
                       print*,"T=",T,"sqrt(T)=",sqrt(T)
                       print*,"x_e=",x_e
                       print*,"exp(-h_P*nu/kB/T)=",exp(-h_P*nu/kB/T)
                       print*,"g=",g_ff

                       print*,"n_b_sq=",n_b_sq
                       print*,"epsilon=",epsilon
                       print*,"R_UCMH=",R_UCMH
                       print*,"I_ind=",I_ind

                       print*,"r_z=",r(z)
                       print*,"Hz=",Hz(z)
                       print*,"f_sky=",f_sky
                       print*,"dndz_PS=",dndz
                       print*,"black=",2.0*kB*T*(nu/c)**2.0*1.0d7
                       print*,"nIhalo=",nIhalo(i)
                       ! stop
                    end if
                 end if
              end do

              ! if(mod(k,1000).eq.0)print*,"minoda wrote (7): z_f=",z_f,"dz_f=",dz_f
           end if
        end if

     end do
     dVdz=4.0*Pi*r(z)*r(z)*c/Hz(z)
     !  dVdz=4.0*Pi*c/Hz(z)
     do i=1,60
        dIdz(i)=dIdz(i)+nIhalo(i)*dVdz*dz/(1.0+z)/(1.0+z)/(1.0+z)
        if(i.eq.50.and.j.le.1000.and.mod(j,10).eq.0)then
           write(26,*)z,nIhalo(i),dVdz,dIdz(i),r(z)
        end if
     end do

     if(z>200.0.and.mod(j,10).eq.0)then
        print*,"A",z,j,dz
     else if(z<200.1.and.z>30.0.and.mod(j,1).eq.0)then
        print*,"B",z,j,dz
     else if(z<30.1.and.z.ge.1.0.and.mod(j,1).eq.0)then
        print*,"C",z,j,dz
        ! print*,dIdz(10),dIdz(20),dIdz(30),dIdz(40),dIdz(50),dIdz(60)
     else if(z<1.0)then
        print*,"D",z,j,dz
     end if

  end do


  print*,"M_minihalo=",M_minihalo,"A_mat=",A_mat,"sigma0=",sqrt(A_mat)
  k_spike=2.0*Pi/((M_minihalo/1.30e11)/(Omega_mhsq/0.112))**(1.0/3.0)
  ! P_zeta=k_spike*c/dsqrt(3.0*(1.0+z_eq)*Omega_mhsq)
  P_zeta=log(4.0d-7*k_spike*c/sqrt(3.0*(1.0+z_eq)*Omega_mhsq))+0.577-3.5
  P_zeta=A_mat/81.0/(1.0+z_eq)/(1.0+z_eq)/P_zeta/P_zeta
  ! P_zeta=A_mat/81.0/(1.0+z_eq)/(1.0+z_eq)/(log(4.0d-7*k_spike*c/sqrt(3.0*(1.0+z_eq)*Omega_mhsq))+0.577-3.5)**2.0
  print*,"k_spike=",k_spike,"P_zeta=",P_zeta

  do i=1,60
     Trad=1.0d-7*dIdz(i)*0.5*c*c/nu_obs(i)/nu_obs(i)/kB
     !     write(5,*)nu_obs(i),dIdz(i)
     write(5,*)1.42040575d9/nu_obs(i)-1.0,Trad
     write(25,*)1.0d-9*nu_obs(i),1.0d23*dIdz(i)
     ! write(5,*)1.0d-9*nu_obs(i),1.0d23*dIdz(i)
     if(nu_obs(i)<200.0d6)then
        print*,1.0d-6*nu_obs(i),"MHz",Trad
     end if
     if(nu_obs(i)>10.0d9)then
        print*,1.0d-9*nu_obs(i),"GHz",1.0d23*dIdz(i),"Jy/str"
     end if
  end do
  close(5)
  close(15)
  close(25)

  ! do i=0,1000
  !    z=1000.0-i
  ! open(3,status='replace',form='formatted',file="dndz_PS.txt")
  ! open(4,status='replace',form='formatted',file="nz_PS.txt")
  ! n_PS=0.0
  ! do i=1000,0,-1
  !    z=dble(i)
  !    a_c=1.0/(1.0+z)
  !    dndz_PS = dexp(-0.5*del_c*del_c/A_mat/a_c/a_c)
  !    dndz_PS = 2.0*rho_over_M*del_c/dsqrt(2.0*Pi*A_mat)*dndz_PS
  !    write(3,*)z,dndz_PS
  !    n_PS = n_PS + dndz_PS
  !    write(4,*)z,n_PS
  ! end do
  ! close(3)
  ! close(4)
  !    ! print*,z,n_BBKS,derfc(del_c*(1.d0+z)/dsqrt(A_mat))/(4.0*Pi*Pi*dsqrt(27.d0))
  !    write(3,*)z,n_BBKS,n_PS
  ! end do

contains
  double precision function ionization(Temp)
    double precision Temp
    double precision A_rec,C_coll,U_T
    double precision,parameter::kB=1.38064852d-23   ! [J/K]
    A_rec=1.14d-13*4.309*(Temp*1.0d-4)**(-0.6166)
    A_rec=A_rec/(1.0+0.6703*(Temp*1.0d-4)**0.53)
    U_T=13.6*1.60218d-19/kB/Temp
    C_coll=0.291d-7*U_T**0.39/dexp(U_T)/(0.232+U_T)
    ionization=C_coll/(A_rec+C_coll)
  end function ionization

  double precision function G_of_x(x,n)
    double precision x,n
    ! G_of_x=G_tilde_true*G_tilde_true
    G_of_x=x/((log(1.0+x)-x/(1.0+x))**((5.0+n)/6.0))
  end function G_of_x


  double precision function rho_boc_sq_int(A_inp,z_inp) ! 3.0 * int^1_0 r^2 (rho_gas(r,z)/rho_crit(z))^2 dr
    double precision A_inp,z_inp,dr_scale,r_scale,rho_boc
    double precision y,F_y,Int_t,t,dt
    double precision n_eff,alpha,A_con,B_con,C_con,G_tilde,G_x
    double precision G_tmin,G_tmax,G_true,G_error,NewGx,nu_mat
    double precision aofz,O_mz
    double precision,parameter::del_c=1.69d0
    ! double precision,parameter::del_c=1.686d0
    double precision,parameter::Omega_m=0.3062
    double precision,parameter::Omega_b=0.0493
    double precision,parameter::Omega_mhsq=0.14
    double precision,parameter::Omega_bhsq=0.0224
    integer ix_radius,i_Gt
    integer,parameter::max_radius=1000
    integer,parameter::i_Gtend=100


    dr_scale=1.0/dble(max_radius)
    r_scale=0.0
    rho_boc_sq_int=0.0

    ! Determining concentration parameter
    ! n_eff=2.0*0.8-3.0           ! 2.0*dln\sigma/dlnR-3.0,
    n_eff=-3.0                  ! dlnP(k)/dlnk
    ! https://www.uio.no/studier/emner/matnat/astro/AST4320/h12/undervisningsmateriale/psmassfunction.pdf
    ! minoda added 2011.11.10
    ! -dlnD/dln(1+z), alpha is 1 in matter-dominated era.
    aofz=1.0/(1.0+z_inp)
    O_mz = Omega_m/(Omega_m+(1.0-Omega_m)*aofz*aofz*aofz)
    alpha=O_mz**0.55
    ! alpha=1.0
    A_con=2.67*(1.0+1.23*(n_eff+3.0))
    B_con=3.92*(1.0+1.30*(n_eff+3.0))
    C_con=1.0+0.19*(1.0-alpha)

    open(17,status='replace',form='formatted',file="Gx_debug.txt")
    ! !!!!!!!!! for debug !!!!!!!!!!
    ! do i_Gt=0,100
    !    G_tilde=10.0**(-2.0+0.05*i_Gt)
    !    G_x=G_tilde/((log(1.0+G_tilde)-G_tilde/(1.0+G_tilde))**((5.0+n_eff)/6.0))
    !    write(17,*)G_tilde,G_x
    ! end do
    ! stop
    ! to determine G_tilde.
    ! nu_mat=del_c/dsqrt(A_inp)/D_plus(z_inp)
    nu_mat=del_c/dsqrt(A_inp)*(1.0+z_inp)
    G_true=A_con/nu_mat*(1.0+nu_mat*nu_mat/B_con)
    G_tmin=1.0d-2
    G_tmax=1.0d4
    if(G_true<G_of_x(G_tmin,n_eff).or.G_of_x(G_tmax,n_eff)<G_true)then
       write(17,*)"Minoda Error: z=",z_inp,", G_true =",G_true,","
       write(17,*)"there is no G_tilde in range of (",G_tmin,",",G_tmax,")."
       stop
    else
       write(17,*)"trial, G_tilde, G(G_tilde), G_true, error[%]"
       do i_Gt=1,i_Gtend
          G_tilde=dsqrt(G_tmin*G_tmax)
          NewGx=G_of_x(G_tilde,n_eff)
          G_error=abs(NewGx-G_true)/G_true
          if(G_error<1.0d-10)exit
          if(NewGx>G_true)then
             G_tmax=G_tilde
          else
             G_tmin=G_tilde
          end if
          ! write(17,*)i_Gt-1,G_tilde,NewGx,G_true,G_error*100.0
       end do
       ! write(17,*)i_Gt-1,G_tilde,NewGx,G_true,G_error*100.0
       if(i_Gt>i_Gtend)then
          write(17,*)"Minoda Error: relative error is too large."
          write(17,*)"Please increase i_Gtend or limit range of (G_tmin,G_tmax)."
          stop
       else
          ! write(17,*)"Successfully determined G_tilde."
       end if
    end if

    y=C_con*G_tilde
    ! minoda: must be commented!
    ! write(27,*)y
    ! write(27,*)z_inp,y
    ! write(17,*)"concentraton parameter is",y
    close(17)

    F_y=log(1.0+y)-y/(1.0+y)
    Int_t=0.0
    t=0.0
    dt=y*dr_scale
    do ix_radius=1,max_radius
       t=t+dt
       Int_t=Int_t+(1.0+t)**(2.0*y/F_y/t)*t*t*dt
    end do


    ! open(16,status='replace',form='formatted',file="reduced_density.txt")
    do ix_radius=1,max_radius
       r_scale=r_scale+dr_scale
       if(gas_profile.eq.0)then
          ! minoda: uniform density profile
          rho_boc=200.0*Omega_bhsq/Omega_mhsq
       else if(gas_profile.eq.1)then
          ! minoda: hydrostatic profile
          rho_boc=200.0/3.0*y**3.0/Int_t*Omega_bhsq/Omega_mhsq*(1.0+r_scale*y)**(2.0/r_scale/F_y)
       else
          print*,"Error: Check the gas profile parameter!!!"
          stop
       end if
       ! write(16,*)r_scale,rho_boc
       rho_boc_sq_int=rho_boc_sq_int+rho_boc*rho_boc*r_scale*r_scale*dr_scale
       ! minoda: warning!!! for debug!!!
       ! rho_boc_sq_int=rho_boc_sq_int+rho_boc*r_scale*r_scale*dr_scale
    end do
    rho_boc_sq_int=3.0*rho_boc_sq_int
    ! close(16)
  end function rho_boc_sq_int

  double precision function Hz(z_inp) ! unit is [cm/s/Mpc]
    double precision z_inp
    double precision,parameter::h_red=0.7
    double precision,parameter::Omega_m=0.3062
    double precision,parameter::Omega_mhsq=0.14
    Hz=1.0d7*h_red*sqrt(1.0-Omega_m+Omega_m*(1.0+z_inp)*(1.0+z_inp)*(1.0+z_inp))
    !  Hz=1.0d7*sqrt(1.0-Omega_mhsq+Omega_mhsq*(1.0+z_inp)*(1.0+z_inp)*(1.0+z_inp))
  end function Hz

  double precision function r(z_inp)
    double precision z_inp,z_int,dz
    integer ii
    dz=0.001*z_inp
    r=0.0
    do ii=1,1000
       z_int=ii*dz
       r=r+c/Hz(z_int)*dz
    end do
  end function r

  double precision function D_plus(z_input)
    double precision z_input,aofz,D0,O_mz
    double precision,parameter::Omega_m0=0.3062
    aofz=1.0/(1.0+z_input)
    O_mz = Omega_m0/(Omega_m0+(1.0-Omega_m0)*aofz*aofz*aofz)
    D0=2.5*Omega_m0/(Omega_m0**(4.0/7.0)-(1.0-Omega_m0)+(1.0+0.5*Omega_m0)*(1.0+(1.0-Omega_m0)/70.0))
    D_plus=2.5*aofz*O_mz/(O_mz**(4.0/7.0)-(1.0-O_mz)+(1.0+0.5*O_mz)*(1.0+(1.0-O_mz)/70.0))/D0
  end function D_plus

  double precision function dndz_PS(a_c,A_mat,M_halo)
    double precision a_c,A_mat,M_halo,rho_over_M
    double precision,parameter::Pi=3.141592653589d0
    double precision,parameter::del_c=1.69d0 ! critical overdensity for spherical collapse
    ! double precision,parameter::del_c=1.686d0 ! critical overdensity for spherical collapse
    double precision,parameter::Omega_mhsq=0.1408
    double precision,parameter::rho_c=2.775d11      ! [h^2 M_sun / Mpc^3]
    double precision,parameter::Omega_m0=0.3062
    double precision Dplus, O_mz, D0
    rho_over_M=Omega_mhsq*rho_c/M_halo
    dndz_PS = 2.0*rho_over_M*del_c/dsqrt(2.0*Pi*A_mat) ! minoda modified 21.09.15
    ! ! O_mz = Omega_m0/(Omega_m0+(1.0-Omega_m0)*a_c*a_c*a_c)
    ! ! D0=2.5*Omega_m0/(Omega_m0**(4.0/7.0)-(1.0-Omega_m0)+(1.0+0.5*Omega_m0)*(1.0+(1.0-Omega_m0)/70.0))
    ! ! Dplus=2.5*a_c*O_mz/(O_mz**(4.0/7.0)-(1.0-O_mz)+(1.0+0.5*O_mz)*(1.0+(1.0-O_mz)/70.0))/D0
    ! ! ! dndz_PS = 2.0*rho_over_M*del_c/dsqrt(2.0*Pi*A_mat*a_c*a_c)
    ! ! dndz_PS = dexp(-0.5*del_c*del_c/A_mat/Dplus/Dplus)*dndz_PS
    ! dndz_PS = dexp(-0.5*del_c*del_c/A_mat/D_plus(1.0/a_c-1.0)/D_plus(1.0/a_c-1.0))*dndz_PS
    dndz_PS = dexp(-0.5*del_c*del_c/A_mat/a_c/a_c)*dndz_PS
  end function dndz_PS

  double precision function dndz_BBKS(a_c,A_mat,M_halo)
    double precision a_c,A_mat,M_halo
    double precision nu_a
    double precision,parameter::Pi=3.141592653589d0
    double precision,parameter::del_c=1.69d0 ! critical overdensity for spherical collapse
    ! double precision,parameter::del_c=1.686d0 ! critical overdensity for spherical collapse
    double precision f_a,f_b,f_c,f_d,f_e,h_nu

    nu_a = del_c/dsqrt(A_mat)/a_c
    f_a = derf(dsqrt(2.5d0)*nu_a)+derf(0.5d0*dsqrt(2.5d0)*nu_a)
    f_b = 0.5d0*(nu_a*nu_a*nu_a-3.d0*nu_a)
    f_c = dsqrt(0.4d0/Pi)
    f_d = (31.d0/4.d0*nu_a*nu_a+8.d0/5.d0)*dexp(-5.d0/8.d0*nu_a*nu_a)
    f_e = (0.5d0*nu_a*nu_a-8.d0/5.d0)*dexp(-2.5d0*nu_a*nu_a)
    ! h_nu = nu_a/(4.0*Pi*Pi*dsqrt(27.d0))
    h_nu = nu_a/(4.0*Pi*Pi*dsqrt(27.d0))*dexp(-0.5*nu_a*nu_a)
    h_nu = h_nu*(f_a*f_b + f_c*(f_d+f_e))
    dndz_BBKS = h_nu*a_c*4.0d13/M_halo
    ! print*,1./a_c-1.0,nu_a
  end function dndz_BBKS

end program
