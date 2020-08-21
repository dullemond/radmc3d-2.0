!========================================================================
!                          POLARIZATION MODULE
!
!                Based on subroutines kindly provided by 
!                Michiel Min, Aug 2009 (Thanks Michiel!!)
!
!------------------------------------------------------------------------
!
! DEFINITIONS OF THE STOKES VECTOR COMPONENTS
!
! The Stokes vector is (I,Q,U,V) or, dependent on how you look at it,
! (E,Q,U,V). This vector has meaning in a plane perpendicular to the
! direction of motion of the photon or light ray. There is a degree of
! freedom of chosing the angle of orientation of the (x',y') coordinate
! system in that plane and its right/left-handedness. It is important
! to define exactly how the coordinate system is set up, and how the
! Q, U and V Stokes components are defined with respect to these.
!
! CONVENTIONS: 
! There are different conventions for how to set up the coordinate system
! and define the Stokes vectors. Our definition follows the IAU 1974 
! definition as described in:
!
!   Hamaker & Bregman (1996) A&AS 117, pp.161
!
! In this convention the x' axis points to the north on the sky, while
! the y' axis points to the east on the sky. The z' axis points to the
! observer. This coordinate system is positively right-handed. The 
! radiation moves toward positive z'. Angles in the (x',y') plane
! are measured counter-clockwise (angle=0 means positive x' direction,
! angle=pi/2 means positive y' direction).
!
! In the following we will (completely consistent still with the IAU
! definitions above) define "up" to be positive y' (=east) and "right" to be
! positive x' (=north). So, the (x',y') coordinates are in a plane
! perpendicular to the photon propagation, and oriented as seen by the
! observer of that photon. So the direction of propagation is toward you,
! while y' points up and x' points to the right, just as one would normally
! orient it. Compared to the IAU definition we will point our head in east
! direction, so to speak. 
!
! Note that this is fully equivalent to adjusting the IAU definition to have
! x' pointing west and y' pointing north, which is perhaps a bit more
! intuitive, since most images in the literature have this orientation.  So
! for convenience of communication, let us simply adjust the IAU 1974
! definition to have positive $x'$ ("right") pointing west and positive $y'$
! ("up") pointing north. It will have no further consequences for the
! definitions and internal workings of RADMC-3D because RADMC-3D does not
! know what ``north'' and ``east'' are.
!
! The (Q,U) definition (linear polarization) is such that a linearly
! polarized ray with Q=+I,U=V=0 has the electric field in the (x',y')=(1,0)
! direction, while Q=-I,U=V=0 has the electric field in the (x',y')=(0,1)
! direction. If we have Q=0, U=+I, V=0 then the E-field points in the x'=y'
! direction, while Q=0, U=-I, V=0 the E-field points in the x'=-y'
! direction.
!
! The (V) definition (circular polarization) is such that (quoting directly
! from the Hamaker & Bregman paper): "For *right*-handed circularly
! polarized radiation, the position angle of the electric vector at any
! point increases with time; this implies that the y' component of the field
! lags the x' component. Also the electric vectors along the line of sight
! form a *left*-handed screw. The Stokes V is positive for right-handed
! circular polarization."
!
! We can put these definitions into the standard formulae:
!
!   Q = I*cos(2*beta)*cos(2*chi)
!   U = I*cos(2*beta)*sin(2*chi)
!   V = I*sin(2*beta)
!
! The angle chi is the angle of the E-field in the (x',y') coordinates,
! measured counter-clockwise from x' (consistent with our definition
! of angles). Example: chi = 45 deg = pi/4, then cos(2*chi)=0 and
! sin(2*chi)=1, meaning that Q=0 and U/I=+1. Indeed this is consistent 
! with the above definition that U/I=+1 is E_x'=E_y'. 
!
! The angle 2*beta is the phase difference between the y'-component of the
! E-field and the x'-component of the E-field such that for 0<beta<pi/2 the
! E-field rotates in a counter-clockwise sense. In other words: the y'-wave
! lags 2*beta behind the x' wave. Example: if we have beta=pi/4, i.e.
! 2*beta=pi/2, then cos(2*beta)=0 and sin(2*beta)=1, so we have Q=U=0 and
! V/I=+1. This corresponds to the y' wave being lagged pi/2 behind the x'
! wave, meaning that we have a counter-clockwise rotation. If we use the
! right-hand-rule and point the thumb into the direction of propagation
! (toward us) then the fingers indeed point in counter-rotating direction,
! meaning that V/I=+1 is righthanded polarized radiation.
!
! In terms of the electric fields:
!
!   E_x'(t) = E_h cos(wt-Delta_h)
!   E_y'(t) = E_v cos(wt-Delta_v)
!
! (with E_h>0 and E_v>0 and Delta_h,v are the phase lags of the components
! with respect to some arbitrary phase) we can write the Stokes components
! as:
!
!   I = E_h^2 + E_v^2
!   Q = E_h^2 - E_v^2
!   U = 2 E_h E_v cos(Delta)
!   V = 2 E_h E_v sin(Delta)
!
! with 
!
!   Delta = Delta_v - Delta_h = 2*beta
!
! NOTE: The IAU 1974 definition is different from the definitions used
!       in the Planck mission, for instance. So be careful. There is
!       something said about this on the website of the healpix software:
!       http://healpix.jpl.nasa.gov/html/intronode12.htm
!
!       The IAU 1974 definition is also different from the Mishchenko
!       book and papers (see below).
!
! If we want to use these definitions for a photon moving somewhere in
! space (where "north" and "east" are not defined), all these definition
! are therefore only meaningful if we define the orientation of the (x',y')
! basis. To do this we introduce a unit vector (Sx,Sy,Sz) that is
! perpendicular to the direction of motion and is, by definition, pointing
! in the (x',y')=(0,1) direction. In other words: the S-vector points into
! the negative-Q direction.
!
! We can transform from a '-basis to a ''-basis by rotating the S-vector
! counter-clockwise (as seen by the observer watching the radiation) by an
! angle ang. Any vector (x',y') in the '-basis will become a vector
! (x'',y'') in a ''-basis, given by the transformation:
!
!  (x'') = ( cos ang   sin ang ) (x')
!  (y'') = ( -sin ang  cos ang ) (y')
!
! NOTE: We choose (x',y') to be the usual counter-clockwise basis for the
! observer seeing the radiation. Rotating the basis in counter-clockwise
! direction means rotating the vector in that basis in clockwise direction,
! hence the sign convention in the matrix. 
!
! If we have (I,Q,U,V) in the '-basis (which we might have written as
! (I',Q',U',V') but by convention we drop the '), the (I'',Q'',U'',V'') in
! the ''-basis becomes
!
!  (I'')   (   1        0            0         0   ) (I)
!  (Q'') = (   0    cos(2*ang)  sin(2*ang)     0   ) (Q)
!  (U'')   (   0   -sin(2*ang)  cos(2*ang)     0   ) (U)
!  (V'')   (   0        0            0         1   ) (V)
!
! This form of matrix and its sign convention can be understood in the
! following way: The Q and U can be seen as I*cos(2*beta)*UnitVector, where
! UnitVector can be expressed as a complex number exp(2*i*chi). However, the
! actual vector pointing in the major-axis direction is exp(i*chi) in the
! complex plane. The transformed complex number (rotating clockwise! by ang)
! means multiplying that with exp(-i*ang).  So exp(i*chi'') =
! exp(-i*ang)*exp(i*chi) = exp(i*(chi-ang)). So chi''=chi-ang (not
! surprisingly!). For UnitVector'' ==exp(i*chi'') this means UnitVector'' =
! exp(-2*i*ang)*UnitVector.  In other words: it is like rotating with angle
! 2*ang.
!
!------------------------------------------------------------------------
!
! SCATTERING OFF RANDOMLY ORIENTED NON-SPHERICAL PARTICLES
!
! Now let's look at the scattering. For this we need the scattering Mueller
! matrix Z(theta,phi). Here theta is the degree of directional change in
! radian, with theta=0 meaning forward scattering, theta=pi/2 is
! perpendicular scattering and theta=pi is back-scattering. The phi angle is
! the scattering direction in the (x',y') plane. It defines the "scattering
! plane" spanned by the incoming and outgoing photon direction. For phi=0
! this scattering plane is in the (x',y')=(1,0) direction. In that case,
! the S-vector acts as the rotation axis of the scattering. 
!
! NOTE: Our convention is that the local 2-D (x',y') basis is such that it
! defines a 2-D plane perpendicular to the INCOMING photon. This photon is
! pointing toward the observer, i.e. if we make 3-D vectors of the x'- and
! y'-basis: e_x' and e_y', then e_z'=cross(e_x',e_y') points toward the
! observer and gives the direction of propagation of the incoming photon.
! The S-vector is by definition e_y'. If phi=0 and 0<theta<pi, the
! scattering takes place in the x'-direction, i.e. in 3-D this is the plane
! spanned by (e_x',e_z'), i.e. the plane perpendicular to S==e_y'. For phi=0
! and 0<theta<pi the new photon propagation direction after scattering will
! be in positive x' direction, i.e. n_new = a_z*e_z' + a_x*e_x', with
! a_x>=0.
!
! The Stokes parameters of the scattered flux off a single grain can be
! written in terms of the Stokes parameters of the incoming flux using the
! scattering matrix Z:
! 
!   ( I_scat )            ( Z_11   Z_12   Z_13   Z_14 ) ( I_inc )
!   ( Q_scat )    m_grain ( Z_21   Z_22   Z_23   Z_24 ) ( Q_inc )
!   ( U_scat )  = ------- ( Z_31   Z_32   Z_33   Z_34 ) ( U_inc ) 
!   ( V_scat )      r^2   ( Z_41   Z_42   Z_43   Z_44 ) ( V_inc )
!
! where (I,Q,U,V)_inc = flux_inc and (I,Q,U,V)_scat = flux_scat. The
! scattering matrix is Z is the matrix-valued differential cross section for
! the scattering per unit mass of dust (hence the factor m_grain in the
! equation above). Note that the book of van de Hulst calls this matrix F
! instead of Z. The book by Mishchenko calls it Z for aligned particles, but
! F for randomly oriented particles, because in these two cases they
! define the (x',y')_inc and (x',y')_scat basis vectors differently
! (in the fixed alignment case the particle orientation is fixed in
! the lab frame while in the randomly oriented particle case the 
! x',y' for incoming radiation is oriented w.r.t. the outgoing scattering
! angle). For exact normalizations of Z, see below in the section
! comparing our conventions to standard books.
!
! For randomly oriented particles we need to specify Z only as a function of
! theta. This is, however, only true for randomly oriented particles for
! which each helicity has equal amount of counter-helicity. If we have
! magnetic or gravitational alignment of particles then the full
! Z(theta,phi) must be specified.
!
! The scattering matrix Z(theta) is thus given for a '-basis such that the
! S-vector is perpendicular to the scattering plane. Typically your photon
! initially has some arbitrary perpendicular S-vector. Scattering may then
! take place in any arbitrary scattering plane. If we want to use the
! Z(theta), for that scattering direction, you first have to rotate the
! '-basis to the ''-basis for which the NEW S-vector (i.e. the
! x''-direction) is again perpendicular to the scattering plane. Then we can
! apply the Z(theta), giving us the new (I,Q,U,V). This also automatically
! defines the new orientation of the (x',y') plane which, after the
! scattering event, is of course also rotated in 3-D space. The 3-D version
! of the y' axis stays unaltered during scattering and the 3-D version of
! the x'-axis is rotated around the y'-axis (=S-vector) by an angle theta
! in right-handed direction (if you put the S-vector pointing toward you,
! the x' axis rotates counter-clockwise).
!
! If we want to calculate the scattering source function in a given
! direction of the observer, for a given incident angle of radiation, we
! know precisely the theta and phi angles in advance. We can then indeed
! simply rotate the '-basis, multiply the scattering matrix, and rotate to
! the user-specified S-vector. This is relatively simple, and is done
! by the subroutine polarization_randomorient_scatsource().
! 
! However, if we want to choose a new photon direction and polarization for
! a photon package in a Monte Carlo code, then we do not know the theta and
! phi in advance. In fact, we run into the problem that we cannot a-priori
! rotate the '-basis into the scattering direction.  On the contrary, we
! want to know the scattering probability as a function of theta and phi, so
! we NEED the full Z(theta,phi) rather than just Z(theta). Fortunately we
! can still do this with rotation, because multiplying with Z(theta,phi) is
! the same as (1) first rotating the basis by angle ang=phi (using the
! formulae above), (2) multiplying with the Z(theta) matrix, and (3)
! rotating back to the original basis by choosing ang=-phi. 
!
! So let us assume that we have the full scattering matrix Z(theta,phi). Then
! if we integrate the outgoing I_out(theta,phi) over 4*pi for a given
! unpolarized input intensity I, we should get kappa_scat*I:
!
!      /+1             /2pi                             
!      |   dcos(theta) |   dphi  Z(theta,phi)[1,1] * I   
!      /-1             /0                                
! 
!                               = kappa_scat * I
!
! That is how we define the normalization of the scattering Mueller
! matrix. The Z(theta,phi)[1,1] is thus the angular differential cross
! section of scattering per unit mass of dust. For isotropic scattering it
! would thus be Z11=kappa_scat/(4*pi).
!
! For polarized input intensity this should hold as well, as long as the
! grains are randomly oriented and have no netto helicity (i.e. equal
! left- and right-helicity). 
!
! We can verify this by doing the integrals for the q, u and v components:
!
!                      /+1     /2pi      4                      
! kappa_scat(q,u,v) =  |   dmu |   dphi Sum Z(theta,phi)[1,i] * (1,q,u,v)[i]
!                      /-1     /0       i=1                         
!
! where mu=cos(theta0), q=Q/I, u=U/I, v=V/I (i.e. all values between -1 and 1).
! We can work out this integral a bit deeper, because the integration over
! dphi can be done as follows:
! 
!                      /+1      4                  /2pi         
! kappa_scat(q,u,v) =  |   dmu Sum Z(theta)[1,i] * |   dphi Rotate_phi{(1,q,u,v)}[i]
!                      /-1     i=1                 /0               
!
! i.e. we now have only Z(mu) instead of Z(mu,phi). We have
!
!  Rotate_phi{(1,q,u,v)} = 
!     (1,cos(2*phi)*q+sin(2*phi)*u,-sin(2*phi)*q+cos(2*phi)*u,v)
!
! The integral over 2*pi of the middle two components is 0, so we obtain
!
!                          /+1      4 
! kappa_scat(q,u,v) = 2 pi |   dmu Sum Z(theta)[1,i] * (1,0,0,v)[i]
!                          /-1     i=1
!
!                          /+1      
!                   = 2 pi |   dmu ( Z[1,1] + Z[1,4]*v )
!                          /-1     
!
! As one can see, linearly polarized light has the same scattering cross
! section as unpolarized light. This is because of the random orientation
! of the dust particles, and will NOT be the case when the particles are
! somehow aligned. But note that circularly polarized light might have
! a different scattering cross-section if Z[1,4].ne.0. This means that if
! we have netto helicity of our particles, then we also have dependence
! on the input polarization state. 
!
! For scattering off randomly oriented particles with a plane of symmetry
! (i.e. non-helical) the scattering matrix takes the form (e.g. Mishchenko et
! al. 1996 JQRST 55, 535):
!
!           ( Z_11   Z_12    0      0   )
!           ( Z_12   Z_22    0      0   )
!  Z(mu) =  (  0      0     Z_33   Z_34 )
!           (  0      0    -Z_34   Z_44 )
!
! This has 6 independent matrix elements.
!
! For the somewhat more special case of Mie scattering off spherical
! particles the scattering matrix is even somewhat simpler:
!
!           ( Z_11   Z_12    0      0   )
!           ( Z_12   Z_11    0      0   )
!  Z(mu) =  (  0      0     Z_33   Z_34 )
!           (  0      0    -Z_34   Z_33 )
!
! This has only 4 independent matrix elements! A famous code that can
! compute these, and is freely available on the internet, is the code
! by Bohren & Huffman from the appendix of their book. A version that
! has been improved by Bruce Draine in the early 90s is the bhmie code
! from the scatterlib library:
!
!   http://code.google.com/p/scatterlib/
!
! The callbhmie() program from that bhmie library computes not only
! the scattering and absorption cross sections, but also the Z_11,
! Z_12, Z_33 and Z_34 (which they call S_11, S_12, S_33 and S_34,
! but please see below for the different normalization). 
!
! In this code we allow for the 6-element matrix form, i.e. more general
! than the Mie scattering case, i.e. allowing for prolate or oblate
! particles. 
!
! In any case, for these randomly-oriented or spherical grains we have
! Z_14=0, so that
!
!                           /+1      
!  kappa_scat(q,u,v) = 2 pi |   dmu Z(mu)[1,1]
!                           /-1     
!
! in other words:
!
!  kappa_scat(q,u,v) = kappa_scat
!
! This means that to find out the scattering cross section, we do not need
! any information about the polarization state. This is good, because
! otherwise even the unpolarized Monte Carlo scattering model would be
! wrong, because the polarization state would then affect the further
! scattering behavior. This is, in the case of randomly oriented and/or
! spherical grains not the case. For now we shall use this as a working
! hypothesis. NOTE: It is still not perfectly correct, because while the
! total cross section for scattering is the same, the direction in which the
! photons scatter may be different at the second and later scatterings. But
! it is usually not that problematic.
!
! For Monte Carlo we are concerned with photons (or photon packages). We use
! tau_scat = rho*kappa_scat*pathlength to find out IF a photon scatters.
! 
! ONCE we know that it scatters, the angular probability distribution
! function for the direction of scattering is:
!
!               Z(mu,phi)[1,1] + Z(mu,phi)[1,2] * q + Z(mu,phi)[1,3] * u + Z(mu,phi)[1,4] * v
!  P(mu,phi) = -------------------------------------------------------------------------------
!                                            kappa_scat
!
! where the q, u and v are those of the incoming radiation.
!
! The P(mu,phi) is normalized in (mu,phi) space to unity:
! 
!   /+1     /2pi
!   |   dmu |   dphi P(mu,phi) = 1
!   /-1     /0
!
! so P(mu,phi) is the probability per dmu and per dphi. 
!
! So to find the new direction of scattering we use P(mu,phi) as the
! random-number distribution function and pick a direction. Once we know the
! direction, we still need to figure out what the q,u,v are. For that we
! apply the Z(mu,phi) in the 3-stage way we discussed above.
!
! We can work out P(mu,phi) in terms of Z(mu) matrix elements by performing,
! at a given value of phi, a basis rotation of (x',y') to (x'',y'') such
! that the x''-axis points into the direction of scattering, i.e.
! the e_x'' basis vector equals the e_x' vector rotated by phi in 
! counter-clockwise direction. According to the above rules for
! rotation (with ang=phi) the new q'', u'' and v'' then become:
!
!   q'' =  cos(2*phi)*q + sin(2*phi)*u
!   u'' = -sin(2*phi)*q + cos(2*phi)*u
!
! Since in the (x'',y'') coordinates the scattering is now in the positive
! x''-direction, we have that the probability for scattering is
! P=(Z(mu)[1,1]+Z(mu)[1,2]*q'')/kappa_scat because the other matrix elements
! are zero in this coordinate system.  Inserting the above expression for
! q'' yields:
!
!               Z(mu)[1,1] + Z(mu)[1,2] * ( cos(2*phi)*q + sin(2*phi)*u )
!  P(mu,phi) = -----------------------------------------------------------
!                                    kappa_scat
!
! We can integrate this over phi to get the P(mu):
!
!                   Z(mu)[1,1]
!   P(mu) =  2 pi  ------------
!                   kappa_scat
!
! This shows that, for the probability of scattering with a certain
! angle theta=acos(mu), the input polarization state is irrelevant.
! In other words: the mu-phase function for scattering is entirely
! given by Z(mu)[1,1]. Only the phi-phase function depends on the
! input polarization state.
!
! We can also integrate P(mu,phi) over mu:
!
!
!             Zint[1,1] + Zint[1,2] * ( cos(2*phi)*q + sin(2*phi)*u )
!   P(phi) = ---------------------------------------------------------
!                           2*pi * kappa_scat
!
! with 
!
!                    /+1
!   Zint[1,i] = 2 pi |  dmu Z(mu)[1,i]
!                    /-1
!
! Note that Zint[1,1] == kappa_scat.
!
! Now let's return to the question how to find the new direction of the
! photon using P(mu,phi) as probability distribution function. We do this in
! two stages.
!
!  Stage 1: Find phi using P(phi)
!
!  Stage 2: Compute P(mu|phi), which is the probability function for
!           scattering at an angle mu, given angle phi. P(mu|phi) is:
!
!                                    /+1
!            P(mu|phi) = P(mu,phi) / |  dmu' P(mu',phi)
!                                    /-1 
!
!                        P(mu,phi)
!                      = ---------
!                         P(phi)
!
!            Find mu from P(mu|phi).
!
! Once we found both mu and phi, we can rotate (1,q,u,v) to the new phi, and
! multiply with the scattering matrix to get (i'',q'',u'',v''). However, since
! by definition i'' should become 1, because we have dealt with the
! amplitude of the scattering by using kappa_scat and the P(mu,phi), we
! obtain the new q,u,v as:
!
!   q_new = q''/i'',   u_new = u''/i''  and v_new = v''/i''
!   
! That should be it for the Monte Carlo scattering event. 
!
! NOTE: For randomly-oriented and mirror-symmetric particle ensembles
!       the extinction matrix is scalar, i.e. the extinction coefficient
!       is the same for I,Q,U,V and there is no mixing. In other words:
!       K = diag(1,1,1,1)*kappa. This will NOT be the case for oriented
!       particle ensembles.
! 
!------------------------------------------------------------------------
!
! SCATTERING OFF ALIGNED NON-SPHERICAL PARTICLES
!
! Suppose we have particles which have a plane of symmetry (i.e. no
! helicity) and can be aligned. Take, for instance, prolate or oblate
! ellipsoids. If we have gravity or a magnetic field, then we can align
! these particles with their axis of symmetry along the gravity or 
! B-field direction. Then the scattering Mueller matrix depends on
! three angles instead of just one. 
!
! We rotate the (x',y') plane such that the symmetry axis of the dust
! particle is in the y'-z'-plane. The angle of the particle can be 
! written as xi, such that for xi=0 the symmetry axis of the
! particle is pointing in the y'-direction and for xi=pi/2 it is
! pointing in the z'-direction, and for xi=pi/4 it is in the y'=z'
! direction. The scattering angles remain theta and phi. In this 
! case we can, unfortunately, not reduce angles through symmetries.
! We will have to account for the full Z(xi,theta,phi)-dependence.
! Also: We cannot assume that the matrix has the upper-right and
! lower-left matrix elements zero.
!
! If we wish to include this at some point, we must embed the 
! cross-section routine (e.g. the T-matrix code of Mishchenko)
! into the code, because it would make no sense to precalculate
! and store a 4-dimensional array (3 angles and 1 freq), I think.
!
! It is going to be difficult to implement this into the Monte
! Carlo code, because that involves integrals over the cross
! section. For a "last scattering method" this would, however,
! be no problem.
!
! FOR NOW WE DO NOT INCLUDE THIS POSSIBILITY IN RADMC-3D
!
!------------------------------------------------------------------------
!
! ABSORPTION AND EMISSION BY ALIGNED NON-SPHERICAL NON-HELICAL PARTICLES
! 
! If particles are randomly oriented and non-helical, then the thermal
! emission and absorption are non-polarized. However, if the grains are
! aligned (still non-helical) then the thermal emission and absorption are
! going to be linearly polarized. Like the scattering, the absorption
! is handled through a Mueller matrix. 
!
! Let us assume grains with a rotation axis and no helicity. Let us 
! rotate the (x',y') coordinates such that the symmetry axis of the
! grain lies in the y'-z'-plane. Let us define an angle xi such that
! the symmetry axis of the grain is described by a vector 
! n=(0,cos(xi),sin(xi)). The absorption-emission process is then 
! described by:
!
!     (I)    ( A_1 B_nu(T) )   ( A_1    A_2     0      0   ) (I)
!   d (Q)    ( A_2 B_nu(T) )   ( A_2    A_1     0      0   ) (Q)
!  -- (U)  = (      0      ) - (  0      0     A_1     0   ) (U)
!  dt (V)    (      0      )   (  0      0      0     A_1  ) (V)
! 
! where 
!
!  A_1 = 0.5 * ( A_hori + A_vert )
!  A_2 = 0.5 * ( A_hori - A_vert )
!
! with 
!
!  A_hori = A_hori(xi,nu) = rho * kappa_hori(xi,nu)
!  A_vert = A_vert(xi,nu) = rho * kappa_vert(xi,nu)
!
! The A_hori and A_vert are the absorption coefficients for horizontally
! polarized light (along the x'-axis) and vertically polarized light
! (along the y'-axis). 
!
! We need a table of kappa_hori and kappa_vert for a set of angles
! xi (the angle of the grain symmetry axis with the y'-axis) and
! frequencies nu. 
!
! As of version 0.41 this is now included (see below).
!
!------------------------------------------------------------------------
!
! THE "LAST SCATTERING METHOD" FOR POLARIZATION MAPS
!
! Suppose we are interested in the polarized images only. If multiple
! scattering is relevant, in principle we would need to follow also the
! history of polarization of each photon package before it is scattered into
! the line of sight. This would require us to implement polarization into
! the Monte Carlo method. This is not a problem, but it may be extra work
! and not so easy to do right. A simple approximation would be to compute,
! when ray-tracing for the final image, the scattering source function
! simply on the basis of the assumption that the input photon is
! non-polarized. If the scattering source function is dominated by
! e.g. direct stellar radiation then this approximation is presumably very
! good. If, however, the second scatterings also provide a substantial
! contribution to the source function, then this might be not so
! good. Qualitatively, however, I do not think that this is going to make a
! big difference, because second scatterers are much more isotropically
! distributed, and hence their contributions to the polarization will be
! minor. So all in all I believe that the "last scattering method", in which
! only the polarization to the last scattering is taken into account, is a
! very good approximation.
!
!------------------------------------------------------------------------
!
! REFERENCES:
!   Book by Mishchenko, Travis & Lacis "Scattering, absorption
!    and emission of light by small particles" Cambridge Univ Press:
!    http://www.giss.nasa.gov/staff/mmishchenko/books.html
!   Book by Bohren & Huffmann "Absorption and scattering of light
!    by small particles", Wiley-VCH
!
!   NOTE: These books use slightly different definition of the Stokes
!   parameters. The I and Q are defined the same (their hat theta = our x' ;
!   their hat phi = our y'). They point the hat theta axis downward and the
!   hat phi to the right, so that means that while our definitions of I and
!   Q are the same, what we would call horizontal, they call vertical.
!
!   But their U and V signs are exactly opposite to ours (we use the 
!   IAU 1974 standard definition, see above). For them U=I means
!   E_x' = -E_y' and V=I means clockwise polarization. So the transformation
!   between their and our Stokes is:
!
!     I_ours  =  I_mishch   = I_bohrenhuffman
!     Q_ours  =  Q_mishch   = Q_bohrenhuffman
!     U_ours  =  -U_mishch  = -U_bohrenhuffman
!     V_ours  =  -V_mishch  = -V_bohrenhuffman
!
!   For the Mueller matrices the Z_11, Z_12, Z_21, Z_22, as well as the
!   Z_33, Z_34, Z_43, Z_44 stay the same, while any possible non-zero
!   Z_13, Z_14, Z_23, Z_24, and Z_31, Z_32, Z_41, Z_42 components would
!   become negative.
!
!   Our normalization of the scattering matrix is the same as Mishchenko,
!   but different from Bohren & Hufmann:
!
!    m_grain * Z_ours = Z_mishch = S_bohrenhuffman / k^2
!
!   with k the wave number such that i*k*r is the wave phase. In other
!   words: Bohren & Huffman's S matrix (e.g. Eq. 4.77 of their book, the
!   2004 edition) is dimensionless (or: 1/ster), while ours has 
!   cm^2/ster/gram-of-dust.
! 
!   BH's and Mishch's definition of the scattering plane is the same as we
!   use. If you look at Mishchenko's Fig. 4.1 you see that the scattering
!   plane (the plane in which both the incoming and the outgoing direction
!   lie) is the x-z plane, i.e. the hat theta-n-plane. Their hat phi
!   direction remains untouched. Since we have hat theta = x', this implies
!   in our language that the x'-z' plane is the scattering plane and
!   positive Theta is that we scatter to the right (when looking into the
!   beam of radiation). This is also our definition.
!
! NOTE: Some additions for aligned grains are at the very end.
!
!========================================================================

module polarization_module
use mathroutines_module
use rtglobal_module
use constants_module

type photon
   double precision :: E,Q,U,V
   double precision :: n(1:3),s(1:3)
end type photon

integer, parameter :: pol_nphi=128
double precision :: pol_lin(1:pol_nphi+1)
double precision :: pol_intsin2phi(1:pol_nphi+1),pol_intcos2phi(1:pol_nphi+1)
logical :: polarization_initialized=.false.

contains

!------------------------------------------------------------------------
!                    INITIALIZE POLARIZATION 
!------------------------------------------------------------------------
subroutine polarization_init()
  implicit none
  integer :: i
  double precision :: phi
  pol_lin(1)           = 0.d0
  pol_intcos2phi(1)    = 0.d0
  pol_intsin2phi(1)    = 0.d0
  do i=2,pol_nphi
     phi               = (i-1.0d0)*twopi/(1.d0*pol_nphi)
     pol_lin(i)        = phi
     pol_intcos2phi(i) = 0.5d0*sin(2*phi)
     pol_intsin2phi(i) = 0.5d0*(1.d0-cos(2*phi))
  enddo
  pol_lin(pol_nphi+1)        = twopi
  pol_intcos2phi(pol_nphi+1) = 0.d0
  pol_intsin2phi(pol_nphi+1) = 0.d0
  polarization_initialized = .true.
end subroutine polarization_init


!------------------------------------------------------------------------
!               MONTE CARLO SCATTERING WITH POLARIZATION
!                    (RANDOMLY ORIENTED PARTICLES)
!
! We assume that a scattering event *is* happening. We assume that it
! has already been determined with which dust grain the scattering 
! event happens. Now we want to find with an appropriate random 
! procedure which direction the scattered photon goes, and what its
! polarization state is.
!
! Before calling this routine you must first determine when/where a
! scattering event takes place. Since we are dealing with randomly
! oriented particles without helicity the total scattering cross
! section is independent of the polarization state of the incident
! light. So we can simply use \int Z11 d\Omega as the kappa_scat.
! Once we determined that a scattering event takes place, this
! subroutine can tell you where the photon moves next and which
! polarization state it has. 
!
! ARGUMENTS:
!   phot%n(1:3)     : direction of propagation (changed on output)
!   phot%s(1:3)     : direction of the reference S vector (changed on output)
!   phot%E,Q,U,V    : Stokes vector of the photon (changed on output), 
!                     in units of erg/(sec*Hz)
!   nmu             : The number of mu=cos(theta) angles at which the
!                     scattering Mueller matrix elements are given.
!   mui             : Values of mu. Has nmu elements. We must have 
!                     mui(1) = -1.d0 and mui(nmu) = +1.d0.
!   thetai          : Values of theta. Has nmu elements.
!                     The reason why both mui and thetai should be given
!                     eventhough mui==cos(thetai), is for speed.
!   Z11,...,Z44     : Scattering matrix coefficient arrays (they are
!                     the Z(mu) values). Each of thes arrays has nmu
!                     elements. The values are specified at the mui
!                     gridpoints. They are angular differential cross
!                     sections and have the dimension of
!                     cm^2/steradian/gram-of-dust.
!   CZ11,CZ12       : These are the cumulative integrals of Z11(mu)
!                     and Z12(mu), e.g. CZ12(1)=0, CZ12(nmu)=Zint12.
!                     See compcm below for more information.
!   compcm          : If you set the flag compcm to true,
!                     then the CZ11 and CZ12 arrays will be calculated
!                     on-the-fly (though the arrays must still be
!                     provided, even if their values are irrelevant). 
!                     But this is of course time-consuming. So better is
!                     to compute these arrays a-priori.
!
! RETURNS:
!   phot%n(1:3)     : New photon propagation direction
!   phot%s(1:3)     : New photon S-vector
!   phot%E,Q,U,V    : New photon Stokes vector
!
!------------------------------------------------------------------------
subroutine polarization_randomorient_mc_scatter(phot,nmu,mui, &
                            thetai,Z11,Z12,Z22,Z33,Z34,Z44,   &
                            CZ11,CZ12,compcm)
  implicit none
  type(photon) :: phot
  integer :: nmu
  logical :: compcm
  double precision :: mui(1:nmu),thetai(1:nmu)
  double precision :: Zint11,Zint12
  double precision :: Z11(1:nmu),Z12(1:nmu),Z22(1:nmu)
  double precision :: Z33(1:nmu),Z34(1:nmu),Z44(1:nmu)
  double precision :: CZ11(1:nmu),CZ12(1:nmu)
  double precision :: theta,x,y,z
  double precision :: E,Q,U,V,phi,r,eps,eps1,cos2a,sin2a
  double precision :: zz11,zz12,zz22,zz33,zz34,zz44
  double precision :: a1,a2,a3,val1,val2,rnumber,dummy
  integer :: i,jlo
  logical :: pol
  !
  ! Make sure the polarization module has been initialized
  !
  if(.not.polarization_initialized) call polarization_init()
  !
  ! Check mu-grid
  !
  if(nmu.lt.2) then
     write(stdo,*) 'ERROR: If you use the full phase function for dust scattering,'
     write(stdo,*) '       you must have at least two mu=cos(theta) points, one for'
     write(stdo,*) '       mu=-1 and one for mu=+1.'
     stop
  endif
  if(mui(1).gt.mui(nmu)) then
     if((mui(1).ne.1.d0).or.(mui(nmu).ne.-1.d0)) then
        write(stdo,*) 'ERROR in polarization module: Mu-grid must go from -1 to 1 or 1 to -1.'
        stop
     endif
  else
     if((mui(nmu).ne.1.d0).or.(mui(1).ne.-1.d0)) then
        write(stdo,*) 'ERROR in polarization module: Mu-grid must go from -1 to 1 or 1 to -1.'
        stop
     endif
  endif
  ! Bugfix 2017.03.15: Now check also if the theta grid is ok
  if(thetai(1).lt.thetai(nmu)) then
     if((thetai(1).ne.0.d0).or.(abs(thetai(nmu)-pi).gt.1.d-6)) then
        write(stdo,*) 'ERROR in polarization module: theta-grid must go from 0 to pi or pi to 0.'
        stop
     endif
  else
     if((thetai(nmu).ne.0.d0).or.(abs(thetai(1)-pi).gt.1.d-6)) then
        write(stdo,*) 'ERROR in polarization module: theta-grid must go from 0 to pi or pi to 0.'
        stop
     endif
  endif
  !
  !##########################################################################
  ! Check the self-consistency of the s and n vectors
  ! (this can be removed for speed, once the code is validated)
  dummy = phot%n(1)**2+phot%n(2)**2+phot%n(3)**2
  if(abs(dummy-1.d0).gt.1d-10) then
     write(stdo,*) 'ERROR in polarization module: n-vector not unit vector'
     stop 41
  endif
  dummy = phot%s(1)**2+phot%s(2)**2+phot%s(3)**2
  if(abs(dummy-1.d0).gt.1d-10) then
     write(stdo,*) 'ERROR in polarization module: s-vector not unit vector'
     stop 42
  endif
  dummy = phot%n(1)*phot%s(1)+phot%n(2)*phot%s(2)+phot%n(3)*phot%s(3)
  if(abs(dummy).gt.1d-10) then
     write(stdo,*) 'ERROR in polarization module: s-vector not perpendicular to n-vector'
     stop 43
  endif
  !##########################################################################
  !
  ! Determine if input photon is polarized
  !
  if((phot%Q.eq.0.d0).and.(phot%U.eq.0.d0).and.(phot%V.eq.0.d0)) then
     pol = .false.
  else
     pol = .true.
  endif
  !
  ! If the CZ11 and CZ12 arrays are not precomputed by the user,
  ! we compute them here on-the-fly. We use Simpson's integration rule.
  !
  if(compcm) then
!###############################################
!    Switch this off for now
     stop
!###############################################
     CZ11(1) = 0.d0
     CZ12(1) = 0.d0
     do i=2,nmu
        CZ11(i) = CZ11(i-1) + twopi*abs(mui(i)-mui(i-1))*0.5d0*(Z11(i-1)+Z11(i))
        CZ12(i) = CZ12(i-1) + twopi*abs(mui(i)-mui(i-1))*0.5d0*(Z12(i-1)+Z12(i))
     enddo
  endif
  !
  ! From this, we can extract the Zint11 and Zint12, which are the
  ! integrals of Z over 4*pi steradian.
  !
  Zint11 = CZ11(nmu)
  Zint12 = CZ12(nmu)
  !
  ! Check for consistency
  !
  if(abs(Zint12).gt.Zint11) then
     write(stdo,*) 'ERROR in polarization module: Zint12>Zint11'
     stop
  endif
  !
  ! Choose angle of scattering compared to the current (x',y') coordinate 
  ! system of the polarization.
  !
  if(pol) then
     !
     ! Incoming photon is polarized
     !
     ! Choose the direction phi (counterclockwise from the x'-axis) in which
     ! the scattering takes place. In the header of this module we explain
     ! how this is done and what the sign conventions are (in particular:
     ! phi is measured counter-clockwise from the x-axis).
     !
     ! The cumulative probability distribution function int P(phi) dphi is:
     !
     !  /phi                /phi  Zint[1,1] + Zint[1,2] * ( cos(2*phi')*q + sin(2*phi')*u )
     !  |   dphi' P(phi') = |    ----------------------------------------------------------- dphi'
     !  /0                  /0                      2*pi * kappa_scat
     !
     ! This can be written as the sum of three precalculated/pretabulated
     ! functions multiplied by three numbers: Zint[1,1]/(2*pi*kappa),
     ! Zint[1,2]*q/(2*pi*kappa) and Z[1,2]*u/(2*pi*kappa). In the method below
     ! we dont care about the normalization (since we simply upscale the
     ! random number accordingly), the three numbers become Z[1,1], Z[1,2]*q
     ! and Z[1,2]*u. Note that the full integrals of the cos and sin parts
     ! integrate to 0 at phi=pi, so that is why we need to normalize the
     ! rnumber only against Z[1,1]. The fact that we can thus use
     ! pretabulated functions can save us a lot of computation time.
     !
     a1      = Zint11
     a2      = Zint12*phot%Q/phot%E
     a3      = Zint12*phot%U/phot%E
     rnumber = a1 * pol_lin(pol_nphi+1) * ran2(iseed)
     call hunt3(a1,a2,a3,pol_lin,pol_intcos2phi,pol_intsin2phi, &
                pol_nphi+1,rnumber,jlo)
     if((jlo.lt.1).or.(jlo.ge.pol_nphi+1)) then
        write(stdo,*) 'ERROR in polarization module: hunt3 failed...'
        stop
     endif
     val1    = (a1*pol_lin(jlo)  +a2*pol_intcos2phi(jlo)  +a3*pol_intsin2phi(jlo))
     val2    = (a1*pol_lin(jlo+1)+a2*pol_intcos2phi(jlo+1)+a3*pol_intsin2phi(jlo+1))
     eps     = (rnumber-val1) / (val2-val1)
     if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 3901
     phi     = (1.d0-eps)*pol_lin(jlo) + eps*pol_lin(jlo+1)
     !
  else
     !
     ! If incoming photon is unpolarized, phi angle is random.
     ! This goes much faster, of course!
     !
     phi=ran2(iseed)*twopi
  endif
  !
  ! The scattering plane is now chosen.
  !
  ! Now rotate the Stokes vector to the new plane, i.e. such that the new
  ! reference plane for the polatization (determined by the unit vector
  ! (s(1),s(2),s(3)) that is perpendicular to (n(1),n(2),n(3)) becomes
  ! perpendicular to the plane in which the scattering deflection takes
  ! place. Note that the s-vector points in the y' direction while the
  ! scattering takes place in the positive x' direction.
  !
  call polarization_rotateunitvector(phot%s,phot%n,phi)
  !
  ! Then rotate the Stokes vector to the new reference system. 
  !
  if(pol) then
     !
     ! Compute the cos(2*phi) and sin(2*phi)
     !
     cos2a = cos(2*phi)
     sin2a = sin(2*phi)
     !
     ! Do the rotation
     !
     Q =  phot%Q*cos2a + phot%U*sin2a
     U = -phot%Q*sin2a + phot%U*cos2a
     !
     ! Copy U and Q back to the photon structure
     ! 
     phot%Q = Q  
     phot%U = U
     !
  endif
  !
  ! From now on the scattering happens in the plane defined by the
  ! S-vector (being, by definition, perpendicular to that plane).
  ! So the scattering now happens in the positive x' direction.
  !
  ! Now pick a random mu=cos(theta) from the distribution
  !
  rnumber = ran2(iseed)*(CZ11(nmu)*phot%E+CZ12(nmu)*phot%Q)
  a1      = phot%E
  a2      = phot%Q
  call hunt2(a1,a2,CZ11,CZ12,nmu,rnumber,jlo)
  if((jlo.lt.1).or.(jlo.ge.nmu)) then
     write(stdo,*) 'ERROR in polarization module: hunt2 failed...'
     stop
  endif
  val1    = (a1*CZ11(jlo)  +a2*CZ12(jlo))
  val2    = (a1*CZ11(jlo+1)+a2*CZ12(jlo+1))
  eps     = (rnumber-val1) / (val2-val1)
  eps1    = 1.d0-eps
  if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 3902
  theta   = eps1*thetai(jlo) + eps*thetai(jlo+1)
  !
  ! Now determine the scattered stokes vector by multiplying 
  ! it with the scattering Mueller matrix
  !
  ! First interpolate the scattering Mueller matrix in mu
  !
  zz11 = eps1*Z11(jlo) + eps*Z11(jlo+1)
  zz12 = eps1*Z12(jlo) + eps*Z12(jlo+1)
  zz22 = eps1*Z22(jlo) + eps*Z22(jlo+1)
  zz33 = eps1*Z33(jlo) + eps*Z33(jlo+1)
  zz34 = eps1*Z34(jlo) + eps*Z34(jlo+1)
  zz44 = eps1*Z44(jlo) + eps*Z44(jlo+1)
  !
  ! Then multiply with Stokes vector
  !
  E    =  zz11 * phot%E + zz12 * phot%Q
  Q    =  zz12 * phot%E + zz22 * phot%Q
  U    =  zz33 * phot%U + zz34 * phot%V
  V    = -zz34 * phot%U + zz44 * phot%V
  !
  ! Now renormalize the photon package to conserve energy.
  ! The reason why we do this is because we assume that
  ! a scattering *is* happening, and we only need to 
  ! know the polarized state. The energy of the package
  ! remains - by definition - the way it is. 
  !
  if(E.ne.0d0) then
     phot%Q = phot%E*Q/E
     phot%U = phot%E*U/E
     phot%V = phot%E*V/E
  else
     phot%E = phot%E
     phot%Q = 0.d0
     phot%U = 0.d0
     phot%V = 0.d0
  endif
  !
  ! Rotate the propagation vector around the rotation axis given by the
  ! S-vector in a positive corkscrew direction by angle theta. Since the
  ! S-vector points upward in the (x',y')-plane, and since the photon moves
  ! toward us (when we see the x',y' plane), positive corkscrew rotation is
  ! indeed in positive x'-direction, as we defined it to be.
  !
  call polarization_rotateunitvector(phot%n,phot%s,theta)
  !
  ! Renormalize propagation direction and S-vector
  !
  r         = sqrt(phot%n(1)**2+phot%n(2)**2+phot%n(3)**2)
  phot%n(1) = phot%n(1)/r
  phot%n(2) = phot%n(2)/r
  phot%n(3) = phot%n(3)/r
  r         = sqrt(phot%s(1)**2+phot%s(2)**2+phot%s(3)**2)
  phot%s(1) = phot%s(1)/r
  phot%s(2) = phot%s(2)/r
  phot%s(3) = phot%s(3)/r
  !
  ! Done...
  !
  return
end subroutine polarization_randomorient_mc_scatter


!------------------------------------------------------------------------
!              COMPUTE THE TOTAL SCATTERING CROSS SECTION
!                    AND CUMULATIVE CROSS SECTION
!                    (RANDOMLY ORIENTED PARTICLES)
!
! This subroutine computes the cumulative scattering cross section:
!
!               /mu                                 /mu                
!   CZ11 = 2 pi |   Z11(mu') dmu'       CZ12 = 2 pi |   Z12(mu') dmu'  
!               /-1                                 /-1                
!
! and 
!                     /+1               
!   kappa_scat = 2 pi |   Z11(mu) dmu  = CZ11(mu=1)
!                     /-1                  
!
! and
!                     /+1               
!   g = < mu > = 2 pi |   Z11(mu) mu dmu / kappa_scat
!                     /-1                  
! 
!------------------------------------------------------------------------
subroutine polarization_total_scattering_opacity(nmu,mui, &
                            Z11,Z12,Z22,Z33,Z34,Z44,      &
                            kappa_scat,g,CZ11,CZ12)
  implicit none
  integer :: nmu
  double precision :: mui(1:nmu)
  double precision :: Z11(1:nmu),Z12(1:nmu),Z22(1:nmu)
  double precision :: Z33(1:nmu),Z34(1:nmu),Z44(1:nmu)
  double precision,optional :: kappa_scat,g,CZ11(1:nmu),CZ12(1:nmu)
  double precision :: CumZ11,CumZ12,CumZ11mu,dummy
  integer :: i
  !
  CumZ11 = 0.d0
  CumZ12 = 0.d0
  CumZ11mu = 0.d0
  if(present(CZ11)) CZ11(1)=CumZ11
  if(present(CZ12)) CZ12(1)=CumZ12
  do i=2,nmu
     dummy    = twopi*abs(mui(i)-mui(i-1))*0.5d0
     CumZ11   = CumZ11 + dummy*(Z11(i-1)+Z11(i))
     CumZ12   = CumZ12 + dummy*(Z12(i-1)+Z12(i))
     CumZ11mu = CumZ11mu + dummy*(Z11(i-1)+Z11(i)) * &
                0.5d0*(mui(i-1)+mui(i))
     if(present(CZ11)) CZ11(i)=CumZ11
     if(present(CZ12)) CZ12(i)=CumZ12
  enddo
  kappa_scat = CumZ11
  if(present(g)) g = CumZ11mu/kappa_scat
end subroutine polarization_total_scattering_opacity


!------------------------------------------------------------------------
!            COMPUTE SCATTERING SOURCE FUNCTION CONTRIBUTION
!                    (RANDOMLY ORIENTED PARTICLES)
!
! In the Monte Carlo code we need to continuously add contributions to
! the scattering source function in each cell visited by the photon.
! This subroutine computes the contributions for all four Stokes
! components. We calculate the j_scat,nu source function for scattering,
! which has the dimension erg/(s*cm^3*Hz*ster). Let us assume that
! the cell our photons pass through has volume V, our photon package
! has luminosity E_phot in units of erg/(s*Hz), the length of the 
! path through the cell that the photon package follows is called L,
! then we need to add the following contribution to j_scat
! (where j_scat is a Stokes vector):
!
!   / j_scat_E \   / j_scat_E \   L                  / E_phot \
!   | j_scat_Q | = | j_scat_Q | + - rho Z(ninc.nobs) | Q_phot |
!   | j_scat_U |   | j_scat_U |   V                  | U_phot |
!   \ j_scat_V /   \ j_scat_V /                      \ V_phot /
!
! here Z is the scattering matrix (see below) with dimension of
! cm^2/(gram*ster). Because the Z matrix is already "per steradian"
! (in contrast to the opacity kappa_scat), we do not need to 
! divide by 4*pi anymore.
!
! This subroutine computes and returns:
!
!   / src(1) \                / E_phot \
!   | src(2) | = Z(ninc.nobs) | Q_phot |
!   | src(3) |                | U_phot |
!   \ src(4) /                \ V_phot /
!
! ARGUMENTS:
!   phot%n(1:3)       : direction of propagation
!   phot%s(1:3)       : direction of the reference S vector
!   phot%E,Q,U,V      : Stokes vector of the photon package, in units of
!                       erg/(sec*Hz). It is, so to speak, the "energy"
!                       (or better: "luminosity") of the photon package.
!   nobs(1:3)         : Unit vector pointing toward observer
!   sobs(1:3)         : Unit vector defining the reference direction
!                       used by the observer. NOTE: sobs must be 
!                       perpendicular to nobs!
!   nmu               : The number of mu=cos(theta) angles at which the
!                       scattering Mueller matrix elements are given.
!   mui               : Values of mu. Has nmu elements.
!   Z11,...,Z44       : Scattering matrix coefficient arrays (they are
!                       the Z(mu) values). Each of these arrays has nmu
!                       elements. The normalization is such that for
!                       isotropic unpolarized scattering Z11=kapscat/(4*pi).
!
! RETURNS:
!   src(1:4)          : The scattering source function for E,Q,U,V, 
!                       (see definition above). This must still be 
!                       multiplied by the rho_dust, divided by the cell
!                       volume and multiplied by ds (path length). 
!------------------------------------------------------------------------
subroutine polarization_randomorient_scatsource(phot,nobs,sobs,nmu,mui, &
                                   Z11,Z12,Z22,Z33,Z34,Z44,src)
  implicit none
  type(photon) :: phot
  integer :: nmu
  double precision :: nobs(1:3)
  double precision :: sobs(1:3)
  double precision :: src(1:4)
  double precision :: mui(1:nmu)
  double precision :: Z11(1:nmu),Z12(1:nmu),Z22(1:nmu)
  double precision :: Z33(1:nmu),Z34(1:nmu),Z44(1:nmu)
  double precision :: sscat(1:3),slen,cosa,sina,cos2a,sin2a
  double precision :: zz11,zz12,zz22,zz33,zz34,zz44
  double precision :: eps,eps1,E,Q,U,V,sq,su,mu,inp,dummy
  integer :: jlo
  logical :: pol
  !
  !##########################################################################
  ! Check the self-consistency of the s and n vectors
  ! (this can be removed for speed, once the code is validated)
  dummy = phot%n(1)**2+phot%n(2)**2+phot%n(3)**2
  if(abs(dummy-1.d0).gt.1d-10) then
     write(stdo,*) 'ERROR in polarization module: n-vector not unit vector'
     stop 41
  endif
  dummy = phot%s(1)**2+phot%s(2)**2+phot%s(3)**2
  if(abs(dummy-1.d0).gt.1d-10) then
     write(stdo,*) 'ERROR in polarization module: s-vector not unit vector'
     stop 42
  endif
  dummy = phot%n(1)*phot%s(1)+phot%n(2)*phot%s(2)+phot%n(3)*phot%s(3)
  if(abs(dummy).gt.1d-10) then
     write(stdo,*) 'ERROR in polarization module: s-vector not perpendicular to n-vector'
     stop 43
  endif
  dummy = nobs(1)**2+nobs(2)**2+nobs(3)**2
  if(abs(dummy-1.d0).gt.1d-10) then
     write(stdo,*) 'ERROR in polarization module: nobs-vector not unit vector'
     stop 51
  endif
  dummy = sobs(1)**2+sobs(2)**2+sobs(3)**2
  if(abs(dummy-1.d0).gt.1d-10) then
     write(stdo,*) 'ERROR in polarization module: sobs-vector not unit vector'
     stop 52
  endif
  dummy = nobs(1)*sobs(1)+nobs(2)*sobs(2)+nobs(3)*sobs(3)
  if(abs(dummy).gt.1d-10) then
     write(stdo,*) 'ERROR in polarization module: sobs-vector not perpendicular to nobs-vector'
     stop 53
  endif
  !##########################################################################
  !
  ! Determine if input photon is polarized
  !
  if((phot%Q.eq.0.d0).and.(phot%U.eq.0.d0).and.(phot%V.eq.0.d0)) then
     pol = .false.
  else
     pol = .true.
  endif
  !
  ! Determine scattering plane
  !
  ! Make the cross product between photon propagation vector and the
  ! direction vector toward the observer.  The sign convention is such that
  ! the scattering will be toward positive x' in the (x',y')-plane, when the
  ! photon moves toward the observer. Example: For an observer at
  ! z=+infinity, i.e. nobs=(0,0,1), a photon coming (as seen by the
  ! observer) from the right and scattering into the line-of-sight, the
  ! photon propagation would be e.g. phot%n=(-1,0,0), the S-vector will be
  ! sscat=(0,1,0). Indeed, the photon will then scatter toward positive x'.
  !
  sscat(1) = phot%n(2) * nobs(3) - phot%n(3) * nobs(2)
  sscat(2) = phot%n(3) * nobs(1) - phot%n(1) * nobs(3)
  sscat(3) = phot%n(1) * nobs(2) - phot%n(2) * nobs(1)
  !
  ! Normalize this vector
  !
  slen    = sscat(1)**2 + sscat(2)**2 + sscat(3)**2
  if(slen.lt.1.d-14) then
     !
     ! If the length of the cross product is so small, then
     ! evidently the n and nobs are virtually parallel to
     ! each other. We can then choose an arbitrary S vector
     ! as long as it is perpendicular to n.
     !
     call polarization_make_s_vector(phot%n,sscat)
  else
     !
     ! Normalize to unit
     !
     slen     = sqrt(slen)
     sscat(1) = sscat(1) / slen
     sscat(2) = sscat(2) / slen
     sscat(3) = sscat(3) / slen
  endif
  !
  !#####################################################################
  ! The following is only for self-consistency check; can be removed later
  inp     = sscat(1)*phot%n(1) + sscat(2)*phot%n(2) + sscat(3)*phot%n(3) 
  if(abs(inp).gt.1d-6) then
     write(stdo,*) 'ERROR in polarization module: Something went wrong with sscat in '
     write(stdo,*) '      polarization_scatsource().'
     stop
  endif
  !#####################################################################
  !
  ! If input photon is polarized then we must rotate the Stokes
  ! vector to a new basis such that the x'-direction is in the
  ! scattering plane, i.e. the new S-vector is equal to the
  ! S-vector that determines the scattering plane (=sscat).
  !
  if(pol) then
     !
     ! Determine the cos(ang) between this sscat and the S-vector 
     ! of the photon
     !
     cosa = sscat(1)*phot%s(1) + sscat(2)*phot%s(2) + sscat(3)*phot%s(3) 
     !
     ! Now the cross- and inner product formula for finding sin(a). The sign
     ! convention is such that if the scattering-plane S-vector is
     ! counter-clockwise from the photon's own S-vector when the photon
     ! propagation direction is pointing toward the observer, then the angle
     ! is positive.
     !
     sina = ( phot%s(2)*sscat(3) - phot%s(3)*sscat(2) ) * phot%n(1) + &
            ( phot%s(3)*sscat(1) - phot%s(1)*sscat(3) ) * phot%n(2) + &
            ( phot%s(1)*sscat(2) - phot%s(2)*sscat(1) ) * phot%n(3)
     !
     ! Since we need cos(2*ang) and sin(2*ang) for the rotation of the
     ! Stokes vector, we compute them here.
     !
     cos2a = cosa**2 - sina**2
     sin2a = 2d0*sina*cosa
     !
     ! Now rotate the Stokes vector of the photon to the new S-vector
     !
     E =  phot%E
     Q =  phot%Q*cos2a + phot%U*sin2a
     U = -phot%Q*sin2a + phot%U*cos2a
     V =  phot%V
     !
  else
     !
     ! Unpolarized, so just copy
     !
     E =  phot%E
     Q =  phot%Q
     U =  phot%U
     V =  phot%V
     !
  endif
  !
  ! Find the scattering angle mu = cos(theta)
  !
  mu = phot%n(1)*nobs(1) + phot%n(2)*nobs(2) + phot%n(3)*nobs(3) 
  !
  ! Find this in the mui table
  !
  call hunt(mui,nmu,mu,jlo)
  eps  = (mu-mui(jlo))/(mui(jlo+1)-mui(jlo))
  eps1 = 1.d0-eps
  if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 445
  !
  ! Check whether the input photon is polarized
  !
  if(pol) then
     !
     ! Yes, input photon is polarized. Do the full thing.
     !
     ! Find the scattering Mueller matrix elements
     !
     zz11 = eps1*Z11(jlo) + eps*Z11(jlo+1)
     zz12 = eps1*Z12(jlo) + eps*Z12(jlo+1)
     zz22 = eps1*Z22(jlo) + eps*Z22(jlo+1)
     zz33 = eps1*Z33(jlo) + eps*Z33(jlo+1)
     zz34 = eps1*Z34(jlo) + eps*Z34(jlo+1)
     zz44 = eps1*Z44(jlo) + eps*Z44(jlo+1)
     !
     ! Then multiply with Stokes vector
     !
     src(1) =  zz11 * E + zz12 * Q
     src(2) =  zz12 * E + zz22 * Q
     src(3) =  zz33 * U + zz34 * V
     src(4) = -zz34 * U + zz44 * V
  else
     !
     ! Input photon is non-polarized, so we can do things
     ! with less effort
     !
     ! Find the scattering Mueller matrix elements
     !
     zz11 = eps1*Z11(jlo) + eps*Z11(jlo+1)
     zz12 = eps1*Z12(jlo) + eps*Z12(jlo+1)
     !
     ! Then multiply with Stokes vector
     !
     src(1) = zz11 * E
     src(2) = zz12 * E
     src(3) = 0.d0
     src(4) = 0.d0
     !
  endif
  !
  ! Finally we must rotate the Stokes source vector to the S-vector given by
  ! the observer.
  !
  ! Determine the cos(ang) between S-vector of the photon and the S-vector
  ! given by the observer.
  !
  cosa = sscat(1)*sobs(1) + sscat(2)*sobs(2) + sscat(3)*sobs(3) 
  !
  ! Now the cross- and inner product formula for sin(a). The sign convention
  ! is such that if the observer-given S-vector is counter-clockwise from
  ! the scattering-plane S-vector, as seen from the observer's perspective,
  ! then the angle is positive.
  !
  sina = ( sscat(2)*sobs(3) - sscat(3)*sobs(2) ) * nobs(1) + &
         ( sscat(3)*sobs(1) - sscat(1)*sobs(3) ) * nobs(2) + &
         ( sscat(1)*sobs(2) - sscat(2)*sobs(1) ) * nobs(3)
  !
  ! Since we need cos(2*ang) and sin(2*ang) for the rotation of the
  ! Stokes vector, we compute them here.
  !
  cos2a  = cosa**2 - sina**2
  sin2a  = 2d0*sina*cosa
  !
  ! Now rotate the Stokes vector of the photon to the new S-vector
  !
  sq     = src(2)
  su     = src(3)
  src(2) =  sq*cos2a + su*sin2a
  src(3) = -sq*sin2a + su*cos2a
  !
  ! Done... 
  ! 
  ! NOTE: This src(1...4) vector still has to be multiplied by the length of
  ! the photon path element and the dust density to get the actual source
  ! term. You can check whether this normalization is correct, because the
  ! first component (src(1)) should then give precisely the same as normal,
  ! for the case when you do not take account of polarization. The scattering
  ! matrix is the angular differential cross section instead of the full
  ! cross-section; therefore we do not need to divide by 4*pi anymore 
  ! (in the "normal" scattering in the monte carlo module we still divided
  ! by 4*pi because we used the kappa_scat, which is 4*pi times larger
  ! than our differential cross section here).
  !
end subroutine polarization_randomorient_scatsource


!-----------------------------------------------------------------------
!              ROTATE VECTOR (V1,V2,V3) AROUND VECTOR (S1,S2,S3)
!                NOTE: BOTH VECTORS MUST BE UNIT VECTORS!
!
! This subroutine rotates the input unit vector V around the rotation axis
! defined by the unit vector S by an angle ang. If ang is positive, then the
! rotation is such that if (s1,s2,s3)=(0,0,1), i.e. the z-axis, then the
! vector rotates counter-clockwise in the x,y plane:
!
!    y    (z-axis is pointing toward you) 
!    ^   
!    |    
!    |    
!    +------> x
!
!-----------------------------------------------------------------------
subroutine polarization_rotateunitvector(v,s,ang)
  IMPLICIT NONE
  double precision :: v(1:3),s(1:3),yy(1:3),ang,inp
  double precision :: cosa,sina
  !
  ! Make the cos(ang) and sin(ang)
  !
  cosa=cos(ang)
  sina=sin(ang)
  !
  ! Compute the inner product of the two unit vectors
  !
  inp   = v(1)*s(1) + v(2)*s(2) + v(3)*s(3)
  !
  ! Do the rotation
  !
  yy(1) = s(1)*inp + (v(1)-inp*s(1))*cosa + (s(2)*v(3)-s(3)*v(2))*sina
  yy(2) = s(2)*inp + (v(2)-inp*s(2))*cosa + (s(3)*v(1)-s(1)*v(3))*sina
  yy(3) = s(3)*inp + (v(3)-inp*s(3))*cosa + (s(1)*v(2)-s(2)*v(1))*sina
  !
  ! Now copy back
  !
  v(1)  = yy(1)
  v(2)  = yy(2)
  v(3)  = yy(3)
  return
end subroutine polarization_rotateunitvector
	
	
!-----------------------------------------------------------------------
!      MAKE ARBITRARY BUT CONSISTENT S-VECTOR FOR GIVEN N-VECTOR
!-----------------------------------------------------------------------
subroutine polarization_make_s_vector(n,s)
  implicit none
  double precision :: n(1:3),s(1:3),dum
  if(abs(n(1)).gt.0.5d0) then
     s(1) = -n(3)
     s(2) = 0.d0
     s(3) = n(1)
  else
     s(1) = 0.d0
     s(2) = n(3)
     s(3) = -n(2)
  endif
  dum = sqrt( s(1)*s(1) + s(2)*s(2) + s(3)*s(3) )
  if(dum.lt.1.d-5) stop 6027
  s(:) = s(:) / dum
end subroutine polarization_make_s_vector


!-----------------------------------------------------------------------
!            VERSION OF THE NUMERICAL RECIPES HUNT() ROUTINE
!-----------------------------------------------------------------------
subroutine hunt2(a1,a2,x1,x2,n,x,jlo)
  integer :: jlo,n
  doubleprecision :: x,x1(n),x2(n),a1,a2
  integer :: inc,jhi,jm
  logical :: ascnd
  ascnd=(a1*x1(n)+a2*x2(n).gt.a1*x1(1)+a2*x2(1))
  if(jlo.le.0.or.jlo.gt.n)then
     jlo=0
     jhi=n+1
     goto 3
  endif
  inc=1
  if(x.ge.a1*x1(jlo)+a2*x2(jlo).eqv.ascnd)then
1    jhi=jlo+inc
     if(jhi.gt.n)then
        jhi=n+1
     else if(x.ge.a1*x1(jhi)+a2*x2(jhi).eqv.ascnd)then
        jlo=jhi
        inc=inc+inc
        goto 1
     endif
  else
     jhi=jlo
2    jlo=jhi-inc
     if(jlo.lt.1)then
        jlo=0
     else if(x.lt.a1*x1(jlo)+a2*x2(jlo).eqv.ascnd)then
        jhi=jlo
        inc=inc+inc
        goto 2
     endif
  endif
3 if(jhi-jlo.eq.1)return
  jm=(jhi+jlo)/2
  if(x.gt.a1*x1(jm)+a2*x2(jm).eqv.ascnd)then
     jlo=jm
  else
     jhi=jm
  endif
  goto 3
END SUBROUTINE hunt2


!-----------------------------------------------------------------------
!            VERSION OF THE NUMERICAL RECIPES HUNT() ROUTINE
!-----------------------------------------------------------------------
subroutine hunt3(a1,a2,a3,x1,x2,x3,n,x,jlo)
  integer :: jlo,n
  doubleprecision :: x,x1(n),x2(n),x3(n),a1,a2,a3
  integer :: inc,jhi,jm
  logical :: ascnd
  ascnd=(a1*x1(n)+a2*x2(n)+a3*x3(n).gt.a1*x1(1)+a2*x2(1)+a3*x3(1))
  if(jlo.le.0.or.jlo.gt.n)then
     jlo=0
     jhi=n+1
     goto 3
  endif
  inc=1
  if(x.ge.a1*x1(jlo)+a2*x2(jlo)+a3*x3(jlo).eqv.ascnd)then
1    jhi=jlo+inc
     if(jhi.gt.n)then
        jhi=n+1
     else if(x.ge.a1*x1(jhi)+a2*x2(jhi)+a3*x3(jhi).eqv.ascnd)then
        jlo=jhi
        inc=inc+inc
        goto 1
     endif
  else
     jhi=jlo
2    jlo=jhi-inc
     if(jlo.lt.1)then
        jlo=0
     else if(x.lt.a1*x1(jlo)+a2*x2(jlo)+a3*x3(jlo).eqv.ascnd)then
        jhi=jlo
        inc=inc+inc
        goto 2
     endif
  endif
3 if(jhi-jlo.eq.1)return
  jm=(jhi+jlo)/2
  if(x.gt.a1*x1(jm)+a2*x2(jm)+a3*x3(jm).eqv.ascnd)then
     jlo=jm
  else
     jhi=jm
  endif
  goto 3
END SUBROUTINE hunt3



!=======================================================================
!                    ROUTINES FOR ALIGNED GRAINS
! 
! The code above assumes randomly oriented grains or spherical grains.
! Below are routines for implementing aligned grains. This will be 
! implemented step-by-step. First only the polarized emission by aligned
! grains. We will not include the physics of how grains are aligned
! and by how much. We will simply assume grains to be aligned in a
! certain direction.
!=======================================================================


!-------------------------------------------------------------------
!                INITIALIZE THE ALIGNED GRAINS MODE
!-------------------------------------------------------------------
subroutine aligned_grains_init(action)
  implicit none
  integer :: action
  logical :: fex
  !
  ! Read the directional field
  !
  call read_grainaligndir_field(action)
  !
  ! Make sure that we have a global cos(eta) grid for the angle
  ! between the line-of-sight and the directional field. This is
  ! necessary because different opacity files can have different
  ! angular grids. 
  !
  inquire(file='alignment_angular_grid.inp',exist=fex)
  if(fex) then
     !
     ! Read the alignment angle grid from the file
     ! alignment_angular_grid.inp
     !
     call read_align_angular_grid(action)
  else
     !
     ! Generate the global alignment angle grid
     ! Let us take the number of points to be 20
     !
     call create_align_angular_grid(20)
  endif
  !
end subroutine aligned_grains_init


!-------------------------------------------------------------------
!                READ THE GLOBAL ALIGNMENT ANGLE ARRAY 
!
! The eta-grid is used for the treatment of aligned grains. 
!-------------------------------------------------------------------
subroutine read_align_angular_grid(action)
  implicit none
  integer :: imu,ierr,action,iformat
  logical :: fex
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(align_mui_grid)) return
  elseif(action.eq.2) then
     if(allocated(align_mui_grid)) deallocate(align_mui_grid)
     if(allocated(align_etai_grid)) deallocate(align_etai_grid)
  endif
  !
  ! Check which file is present
  !
  inquire(file='alignment_angular_grid.inp',exist=fex)
  if(.not.fex) then
     write(stdo,*) 'ERROR: When using aligned grains you must define a global'
     write(stdo,*) '       angular grid with an input file called:'
     write(stdo,*) '       alignment_angular_grid.inp'
     write(stdo,*) '       (see manual for its contents).'
     stop
  endif
  !
  ! Message
  !
  write(stdo,*) 'Reading global angular grid for alignment of grains...'
  call flush(stdo)
  !
  ! Read this file
  !
  open(unit=1,file='alignment_angular_grid.inp',status='old',err=701)
  read(1,*) iformat
  if(iformat.ne.1) then
     write(stdo,*) 'ERROR: File format of alignment_angular_grid.inp must be 1'
     stop
  endif
  read(1,*) align_munr
  !
  ! Allocate the mu angular array. 
  !
  if(allocated(align_mui_grid)) deallocate(align_mui_grid)
  allocate(align_mui_grid(1:align_munr),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in RTGLobal Module: Could not allocate align_mui_grid'
     stop 
  endif
  if(allocated(align_etai_grid)) deallocate(align_etai_grid)
  allocate(align_etai_grid(1:align_munr),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in RTGLobal Module: Could not allocate align_etai_grid'
     stop 
  endif
  !
  ! Now read the angular grid from 0 all the way to 180
  !
  do imu=1,align_munr
     read(1,*,end=702,err=702) align_etai_grid(imu)
     if(imu.gt.1) then
        if(align_etai_grid(imu).le.align_etai_grid(imu-1)) then
           write(stdo,*) 'ERROR in file alignment_angular_grid.inp:'
           write(stdo,*) '   theta grid must be monotonically increasing'
           write(stdo,*) '   from 0 to 90 (the rest is a mirror copy from 90 to 180)'
           stop
        endif
     endif
     align_mui_grid(imu) = cos(align_etai_grid(imu)*pi/180.d0)
  enddo
  close(1)
  !
  ! Do some elementary checks
  !
  if(align_etai_grid(1).ne.0.d0) then
     write(stdo,*) 'ERROR while reading alignment_angular_grid.inp file:'
     write(stdo,*) '      The first theta value MUST be 0'
     stop
  endif
  if(align_etai_grid(align_munr).ne.90.d0) then
     write(stdo,*) 'ERROR while reading alignment_angular_grid.inp file:'
     write(stdo,*) '      The last theta value MUST be 90'
     stop
  endif
  align_mui_grid(1)          = 1.d0
  align_mui_grid(align_munr) = 0.d0
  !
  goto 710
701 continue
  write(stdo,*) 'Could not open file alignment_angular_grid.inp'
  stop 13
702 continue
  write(stdo,*) 'alignment_angular_grid.inp: either wrong format or other reading error'
  stop 
710 continue
  return
end subroutine read_align_angular_grid


!-------------------------------------------------------------------
!                READ THE GLOBAL ALIGNMENT ANGLE ARRAY 
!
! The eta-grid is used for the treatment of aligned grains. 
!-------------------------------------------------------------------
subroutine create_align_angular_grid(nang)
  implicit none
  integer :: nang,imu,ierr
  !
  ! Store the nang
  !
  align_munr = nang
  !
  ! Allocate the mu angular array. 
  !
  if(allocated(align_mui_grid)) deallocate(align_mui_grid)
  allocate(align_mui_grid(1:align_munr),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in RTGLobal Module: Could not allocate align_mui_grid'
     stop 
  endif
  if(allocated(align_etai_grid)) deallocate(align_etai_grid)
  allocate(align_etai_grid(1:align_munr),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in RTGLobal Module: Could not allocate align_etai_grid'
     stop 
  endif
  !
  ! Make the grid between eta
  !
  do imu=1,align_munr
     align_etai_grid(imu) = 0.5d0 * pi * (imu-1.d0) / (align_munr-1.d0)
     align_mui_grid(imu)    = cos(align_etai_grid(imu))
  enddo
  if(abs(align_mui_grid(1)-1.d0).gt.1d-12) stop 538
  if(abs(align_mui_grid(align_munr)).gt.1d-12) stop 539
  align_mui_grid(1)          = 1.d0
  align_mui_grid(align_munr) = 0.d0
  !
end subroutine create_align_angular_grid


!-------------------------------------------------------------------
!                 READ THE ALIGNMENT DIRECTION FIELD
!-------------------------------------------------------------------
subroutine read_grainaligndir_field(action)
  implicit none
  character*80 :: filename1,filename2,filename3
  integer :: ispec,ierr,index,action,icell,i,idum,precis,style,reclen
  integer(kind=8) :: iiformat,nn,kk,reclen8
  logical :: fex1,fex2,fex3
  double precision :: dummy
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(grainalign_dir)) then
        !
        ! The alignment direction is already determined.
        ! Let us make sure that it is normalized to unity.
        !
        do icell=1,nrcells
           index = cellindex(icell)
           dummy = grainalign_dir(1,index)**2 + grainalign_dir(2,index)**2 + grainalign_dir(3,index)**2 
           if(dummy.eq.0.d0) then
              write(stdo,*) 'ERROR: Grain alignment direction vector has length 0'
              stop 73328
           endif
           grainalign_dir(1,index) = grainalign_dir(1,index) / dummy
           grainalign_dir(2,index) = grainalign_dir(2,index) / dummy
           grainalign_dir(3,index) = grainalign_dir(3,index) / dummy
        enddo
        return
     endif
  endif
  !
  ! Default
  !
  precis = 8
  !
  ! Message
  !
  write(stdo,*) 'Reading grain alignment direction field...'
  call flush(stdo)
  !
  ! Create the file name and find if the input file is present and if
  ! it is in text format (.inp) or unformatted (.uinp)
  !
  filename1 = 'grainalign_dir.inp'
  filename2 = 'grainalign_dir.uinp'
  filename3 = 'grainalign_dir.binp'
  inquire(file=filename1,exist=fex1)
  inquire(file=filename2,exist=fex2)
  inquire(file=filename3,exist=fex3)
  idum=0
  if(fex1) idum=idum+1
  if(fex2) idum=idum+1
  if(fex3) idum=idum+1
  if(idum.gt.1) then
     write(stdo,*) 'ERROR: Found more than one file grainalign_dir.*inp'
     stop
  endif
  if(idum.eq.0) then
     write(stdo,*) 'Could not find any grainalign_dir.*inp. For aligned grain mode, we need this. Aborting.'
     stop 3903
  endif
  !
  ! Switch formatted/unformatted
  !
  if(fex1) then
     !
     ! Open formatted ascii style
     !
     style = 1
     open(unit=1,file=filename1)
     !
     ! Read format number
     !
     read(1,*) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of '//TRIM(filename1)//' is invalid/unknown.'
        write(stdo,*) 'Format number = ',iiformat
        stop
     endif
     !
     ! Read number of grid points
     !
     read(1,*) nn
  elseif(fex2) then
     !
     ! Open f77-style unformatted (with records)
     !
     style = 2
     open(unit=1,file=filename2,form='unformatted')
     !
     ! Read format number
     !
     read(1) iiformat,reclen8
     reclen = reclen8
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of '//TRIM(filename2)//' is invalid/unknown.'
        write(stdo,*) 'Format number = ',iiformat
        stop
     endif
     read(1) nn
  else
     !
     ! Open C-compliant binary
     !
     style = 3
     open(unit=1,file=filename3,status='old',access='stream')
     !
     ! Read format number
     !
     read(1) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of '//TRIM(filename3)//' is invalid/unknown.'
        write(stdo,*) 'Format number = ',iiformat
        stop
     endif
     read(1) nn
     precis = nn
     read(1) nn
  endif
  !
  ! Do some checks
  !
  if(nn.ne.nrcellsinp) then
     write(stdo,*) 'ERROR: grainalign_dir.*inp does not have same number'
     write(stdo,*) '       of cells as the grid.'
     write(stdo,*) nn,nrcellsinp
     stop
  endif
  !
  ! Create the grainalign_dir arrays
  !
  if(allocated(grainalign_dir)) deallocate(grainalign_dir)
  if(allocated(grainalign_eff)) deallocate(grainalign_eff)
  allocate(grainalign_dir(1:3,1:nrcells),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate the grainalign_dir() array'
     stop
  endif
  allocate(grainalign_eff(1:nrcells),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate the grainalign_eff() array'
     stop
  endif
  !
  ! Now read the grainalign_dir field
  !
  call read_vectorfield(1,style,precis,3,3,nrcellsinp,1,1,reclen=reclen, &
                        vector0=grainalign_dir)
  !
  ! Close file
  !
  close(1)
  !
  ! Normalize this vector field
  !
  do icell=1,nrcells
     index = cellindex(icell)
     dummy = grainalign_dir(1,index)**2 + grainalign_dir(2,index)**2 + grainalign_dir(3,index)**2 
     dummy = sqrt(dummy)
     if(dummy.gt.0.d0) then
        !
        ! The dir vector must be of length unity
        !
        grainalign_dir(1,index) = grainalign_dir(1,index) / dummy
        grainalign_dir(2,index) = grainalign_dir(2,index) / dummy
        grainalign_dir(3,index) = grainalign_dir(3,index) / dummy
     else
        !
        ! When vector has length 0, then it does not matter in which 
        ! direction this unit vector points. Just take something,
        ! for consistency.
        !
        grainalign_dir(1,index) = 1.0
        grainalign_dir(2,index) = 0.0
        grainalign_dir(3,index) = 0.0
     endif
     !
     ! The length of the original vector gives the efficiency
     ! of the alignment. It cannot be >1
     !
     if(dummy.gt.1.01d0) then
        write(stdo,*) 'ERROR: Length of alignment vector > 1. This would mean ', &
             'more than 100% efficient alignment. Aborting.'
        stop
     endif
     if(dummy.gt.1.d0) then
        dummy = 1.d0
     endif
     grainalign_eff(index) = dummy
  enddo
  !
end subroutine read_grainaligndir_field


!-----------------------------------------------------------------------
!             "RANDOM" THERMAL EMISSION FOR ALIGNED GRAINS
!
! This routine will create the starting photon package for thermal
! emission from aligned grains.
!-----------------------------------------------------------------------
subroutine polarization_random_aligned_thermemis(phot,nmu,mui, &
                         orth,para,opcumul,aligndir,aligneff)
  implicit none
  type(photon) :: phot
  integer :: nmu,imu
  double precision :: mui(1:nmu),aligndir(1:3),aligneff
  double precision :: orth(1:nmu),para(1:nmu),opcumul(1:nmu)
  double precision :: mu,eps,phi,perp(1:3),or,pa,rn,dummy
  !###################################
  dummy = aligndir(1)**2+aligndir(2)**2+aligndir(3)**2
  if(abs(dummy-1.d0).gt.1d-5) stop 4007
  !###################################
  !
  ! Find out at which theta this photon is emitted
  !
  rn  = ran2(iseed)
  call hunt(opcumul,nmu,rn,imu)
  if((imu.lt.1).or.(imu.ge.nmu)) stop 3726
  eps = (rn-opcumul(imu)) / (opcumul(imu+1)-opcumul(imu))
  mu  = (1.d0-eps)*mui(imu) + eps*mui(imu+1)
  !
  ! Since the table only goes from mu=0 to mu=1, we must
  ! also randomly flip sign
  !
  rn  = ran2(iseed)
  if(rn.ge.0.5d0) then
     mu = -mu
  endif
  !
  ! Find out at which phi this photon is emitted
  !
  phi = ran2(iseed) * twopi
  !
  ! Create the direction
  !
  ! ...first create an arbitrary unit vector perpendicular 
  !    to the alignment direction. We can use the subroutine
  !    for creating an S-vector.
  !
  call polarization_make_s_vector(aligndir,perp)
  !
  ! ...now rotate randomly (0..2*pi) around alignment direction
  !
  call polarization_rotateunitvector(perp,aligndir,phi)
  !
  ! ...then create the direction vector using these two unit vectors
  !
  phot%n(:) = mu * aligndir(:) + sqrt(1.d0-mu*mu) * perp(:)
  !
  !#####################################
  dummy = sqrt(phot%n(1)**2+phot%n(2)**2+phot%n(3)**2)
  if(abs(dummy-1.d0).gt.1d-5) stop 8291
  !#####################################
  !
  ! Create an appropriate S-vector: this should lie in the
  ! plane spanned between the phot%n direction vector and
  ! the aligndir vector. Use a kind of Gram-Schmidt process
  ! for this.
  !
  dummy = phot%n(1)*aligndir(1) + phot%n(2)*aligndir(2) + phot%n(3)*aligndir(3) 
  phot%s(:) = aligndir(:) - dummy*phot%n(:)
  dummy = sqrt(phot%s(1)*phot%s(1) + phot%s(2)*phot%s(2) + phot%s(3)*phot%s(3))
  if(dummy.gt.1d-5) then
     !
     ! Normalize the S-vector
     !
     phot%s(:) = phot%s(:) / dummy
  else
     !
     ! The photon emission direction lies so close to the alignment
     ! direction that we assume it to be unpolarized
     !
     call polarization_make_s_vector(phot%n,phot%s)
  endif
  !
  ! Now create the polarization state, normalized to 1.
  ! Since the mu-dependence is already taken care of,
  ! the orth and para must be normelized to 1.
  !
  or     = (1.d0-eps) * orth(imu) + eps * orth(imu+1)
  pa     = (1.d0-eps) * para(imu) + eps * para(imu+1)
  dummy  = or + pa
  if(dummy.le.0.d0) stop 3003
  or     = or / dummy
  pa     = pa / dummy
  phot%E = 1.d0
  phot%Q = aligneff * ( or - pa )
  phot%U = 0.d0
  phot%V = 0.d0
  !
end subroutine polarization_random_aligned_thermemis



end module polarization_module
