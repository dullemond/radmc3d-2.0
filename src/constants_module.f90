module constants_module

  doubleprecision pi,pihalf,twopi,fourpi,sqrtpi
  parameter(pi=3.14159265358979323846264338328d0)
  parameter(pihalf=1.57079632679489661923132169164d0)
  parameter(twopi=6.28318530717958647692528676656d0)
  parameter(fourpi=12.5663706143591729538505735331d0)
  parameter(sqrtpi=1.772453850905515881919d0)

  doubleprecision :: parsec,au
  parameter(parsec=3.08572d18,au=1.49598d13)

  doubleprecision :: cc,kk,mp,hh
  parameter(cc  = 2.9979245800000d10)     ! Light speed             [cm/s]
  parameter(kk  = 1.3807d-16)             ! Bolzmann's constant     [erg/K]
  parameter(mp  = 1.6726d-24)             ! Mass of proton          [g]
  parameter(hh  = 6.6262000d-27)          ! Planck constant 

  integer :: stdi, stdo, fflo, ffli
  logical :: radmc_as_child

end module constants_module
