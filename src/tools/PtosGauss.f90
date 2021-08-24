!     ------------------------------------------------------------------
Subroutine Jacobian (Jac,Det,NDim,NodG,XL,dF)
    !     ------------------------------------------------------------------
    Implicit real*8 (a-h,o-z)
    Dimension dF(NDim,NodG),XL(1)
    !     ojo que en algun lugar dF(NDimE,NodEL)
    !     Hay que hacer compatible esta cuestion de diemensionamientos
    !
    Real*8 JacInv(NDim,NDim)
    Real*8 Jac(NDim,NDim)
    !
    !     Do iGC=1,NDim
    !     Do iLC=1,NDim
    !     JacInv(iGC,iLC)= 0.0d0
    !     Enddo
    !     Enddo
    JacInv = 0.0d0
    !
    Do Nod=1,NodG
        iPNod = (Nod-1)*NDim
        Do iGC=1,NDim
            Do iLC=1,NDim
                JacInv(iLC,iGC) = JacInv(iLC,iGC) + &
                    &              XL( iPNod + iGC ) * dF(iLC,Nod)
            Enddo
        Enddo
    Enddo
    !     Ya no (Uso de rutina IMSL para invertir JacInv y obtener Jac)
    !
    Select Case (NDim)
        Case (1)
            Determinant = JacInv(1,1)
            if (Determinant .NE.0.0d0 ) then
                Jac(1,1)    = 1.0d0 / Determinant
            else
                Write (*,*)'WARNING, Element Jacobian = 0'
                Write (*,*)'WARNING, Element Jacobian = 0'
                Write (*,*)'WARNING, Element Jacobian = 0'
                Stop
            endif
        !
        Case (2)
            Determinant = JacInv(1,1)*JacInv(2,2)-JacInv(1,2)*JacInv(2,1)
            if (Determinant .NE.0.0d0 ) then
                Jac(1,1)    = JacInv(2,2) / Determinant
                Jac(1,2)    =-JacInv(1,2) / Determinant
                Jac(2,1)    =-JacInv(2,1) / Determinant
                Jac(2,2)    = JacInv(1,1) / Determinant
            else
                Write (*,*)'WARNING, Element Jacobian = 0 , dim 2'
                Write (*,*)'WARNING, Element Jacobian = 0 , dim 2'
                Write (*,*)'WARNING, Element Jacobian = 0 , dim 2'
                Stop
            endif
        !
        Case (3)
            A1 = JacInv(3,3)*JacInv(2,2) - JacInv(3,2)*JacInv(2,3)
            B1 = JacInv(3,3)*JacInv(2,1) - JacInv(3,1)*JacInv(2,3)
            C1 = JacInv(3,2)*JacInv(2,1) - JacInv(3,1)*JacInv(2,2)
            !
            A2 = JacInv(3,3)*JacInv(1,2) - JacInv(3,2)*JacInv(1,3)
            B2 = JacInv(3,3)*JacInv(1,1) - JacInv(3,1)*JacInv(1,3)
            C2 = JacInv(3,2)*JacInv(1,1) - JacInv(3,1)*JacInv(1,2)
            !
            A3 = JacInv(2,3)*JacInv(1,2) - JacInv(2,2)*JacInv(1,3)
            B3 = JacInv(2,3)*JacInv(1,1) - JacInv(2,1)*JacInv(1,3)
            C3 = JacInv(2,2)*JacInv(1,1) - JacInv(2,1)*JacInv(1,2)
            !
            Determinant= JacInv(1,1)*A1 - JacInv(1,2)*B1 + JacInv(1,3)*C1
            !
            if (Determinant .NE.0.0d0 ) then
                Jac(1,1)    = A1 / Determinant
                Jac(1,2)    =-A2 / Determinant
                Jac(1,3)    = A3 / Determinant
                !
                Jac(2,1)    =-B1 / Determinant
                Jac(2,2)    = B2 / Determinant
                Jac(2,3)    =-B3 / Determinant
                !
                Jac(3,1)    = C1 / Determinant
                Jac(3,2)    =-C2 / Determinant
                Jac(3,3)    = C3 / Determinant
            else
                Write (*,*)'WARNING, Element Jacobian = 0 , dim 3' , Determinant
                Write (*,*)'WARNING, Element Jacobian = 0 , dim 3' , Determinant
                Write (*,*)'WARNING, Element Jacobian = 0 , dim 3' , Determinant
                Stop
            endif
        !
        Case default
            Write (*,*)'Invalid Case for Determinant of element"s Jacobian'
            Write (*,*)'Dimension of Simplex not contained in 1-3',NDim
            Write (*,*)'WARNING, Jacobian of an element = 0'
            Stop
    !     Determinant=0.0d0
    End Select
    !
    !     Det = Determinant(NDim,JacInv)
    Det = Determinant
    !
    if (Det < 0.0d0) then
        Write (*,*)'!!!!! WARNING, Jacobian of an element < 0'
    End if
    !
    !     CALL DLINRG (NDim, JacInv, NDim, Jac, NDim)
    !
    Return
End Subroutine Jacobian
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Real*8 function Determinant(NDim,Jac)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Implicit real*8 (a-h,o-z)
    Real*8 Jac(NDim,NDim)
    !
    Select Case (NDim)
        Case (1)
            Determinant = Jac(1,1)
        Case (2)
            Determinant = Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1)
        Case (3)
            Determinant = Jac(1,1)*(Jac(3,3)*Jac(2,2)-Jac(3,2)*Jac(2,3)) &
                &        - Jac(1,2)*(Jac(3,3)*Jac(2,1)-Jac(3,1)*Jac(2,3)) &
                &        + Jac(1,3)*(Jac(3,2)*Jac(2,1)-Jac(3,1)*Jac(2,2))
        Case default
            !     Write (*,*)'Invalid Case for Determinant of element"s Jacobian'
            !     Write (*,*)'Dimension of Simplex not contained in 1-3',NDim
            Determinant=0.0d0
    !
    End Select
    !
    Return
End
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine LocalShapeDer (NodEL,NDim,dF,GaussPoint,Iorder,iBu)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Implicit real*8 (a-h,o-z)
    Dimension dF(NDim,NodEL),GaussPoint(NDim)
    Integer d(NDim)
    !     Ejemplo: ndim=3, d(3) vector indicador
    !     de la derivada de primer orden = d/dy-> [0,1,0]
    !     iBu = 1 Bubble element, iBu=0 Linear_Element
    !     NLSF=NDim+1+iBu : Number of Linear Shape Functions
    !     Deriv=1
    NLSF=NDim+1+iBu
    d=0
    !
    Do Nod=1,NodEL
        Do iLC=1,NDim
            !
            d(iLC)=1
            if (Iorder .eq. 1 .and. Nod .GT. NLSF) then
                dF(iLC,Nod)= 0.0d0
            else
                dF(iLC,Nod)= SFSimplex(GaussPoint,d,Iorder,Nod,NDim)
            endif
            d(iLC)=0
        !
        Enddo
    Enddo
    !
    Return
End Subroutine LocalShapeDer
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine GlobalShapeDer (NodEL,NDim,dF_L,dF_G,Jac)
    Implicit real*8 (a-h,o-z)
    Dimension dF_L(NDim,NodEL),dF_G(NDim,NodEL)
    Real*8 Jac(NDim,NDim)
    !
    !     Do Nod=1,NodEL
    !     Do iGC=1,NDim
    !     dF_G(iGC,Nod)=0.0d0
    !     Enddo
    !     Enddo
    dF_G = 0.0d0
    !    C
    Do Nod=1,NodEL
        Do iLC=1,NDim
            Do iGC=1,NDim
                dF_G(iGC,Nod)= dF_G (iGC,Nod)+ dF_L(iLC,Nod)*Jac(iGC,iLC)
            Enddo
        Enddo
    Enddo
    !
    Return
End Subroutine GlobalShapeDer
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine ShapeF (NodEL,NDim,F,GaussPoint,iOrder,iBu)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Implicit real*8 (a-h,o-z)
    Dimension F(NodEL),GaussPoint(NDim)
    !     iBu = 1 Bubble element, iBu=0 Linear_Element
    !     NLSF=NDim+1+iBu : Number of Linear Shape Functions
    !     iOrder = Interpolation order: 1=Linear; 2=Cuadratic
    Integer d(NDim)
    !
    d=0
    NLSF=NDim+1+iBu
    !
    Do Nod=1,NodEL
        if (Iorder .eq. 1 .and. Nod .GT. NLSF) then
            F(Nod)= 0.0d0
        else
            F(Nod)= SFSimplex(GaussPoint,d,iOrder,Nod,NDim)
        endif
    Enddo
    !
    Return
End Subroutine ShapeF
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine GaussRule (Psi,W,NDim,NGP,iSimplex)
    Implicit real*8 (a-h,o-z)
    Dimension PSI(NDim,NGP),W(NGP)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    CaseSimplex: Select case (iSimplex)
        !     iSimplex = 0 -> The element is a Simplex
        !     iSimplex = 1 -> The element is a HiperCube
        !
        Case (0)
            !     The element is a simplex
            !
            CaseNDimSimplex: Select case (NDim)
                !
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                !     gauss points for  tetrahedral elements
                !     last modifief 19/9/01
                !
                Case (3)
                    !
                    CaseNGP3d: Select case (NGP)
                        !
                        Case (1)
                                         !     tetrahedral 1-points formula, degree of precision 1
                                         !     Gauss Points Intrinsic Coordinates
                            PSI(1,1) = 0.25D0
                            PSI(2,1) = 0.25D0
                            PSI(3,1) = 0.25D0
                                         !  Weights/6
                            W(1) = 1.0D0*0.166666666666667D0
                        !
                        Case (4)
                            !     tetrahedral 4-points formula, degree of precision 2
                            !     Gauss Points Intrinsic Coordinates
                            PSI(1,1) = 0.1381966011250105D0
                            PSI(2,1) = 0.1381966011250105D0
                            PSI(3,1) = 0.1381966011250105D0
                            PSI(1,2) = 0.5854101966249685D0
                            PSI(2,2) = 0.1381966011250105D0
                            PSI(3,2) = 0.1381966011250105D0
                            PSI(1,3) = 0.1381966011250105D0
                            PSI(2,3) = 0.5854101966249685D0
                            PSI(3,3) = 0.1381966011250105D0
                            PSI(1,4) = 0.1381966011250105D0
                            PSI(2,4) = 0.1381966011250105D0
                            PSI(3,4) = 0.5854101966249685D0

                                         !     Weights/6
                            WT   = 0.25D0*0.1666666666666667D0
                            W(1) = WT
                            W(2) = WT
                            W(3) = WT
                            W(4) = WT
                        !
                        Case (5)
                                         !     tetrahedral 5-points formula, degree of precision 3
                                         !     Gauss Points Intrinsic Coordinates
                            PSI(1,1) = 0.25D0
                            PSI(2,1) = 0.25D0
                            PSI(3,1) = 0.25D0
                            PSI(1,2) = 0.1666666666666667D0
                            PSI(2,2) = 0.1666666666666667D0
                            PSI(3,2) = 0.1666666666666667D0
                            PSI(1,3) = 0.3333333333333333D0
                            PSI(2,3) = 0.1666666666666667D0
                            PSI(3,3) = 0.1666666666666667D0
                            PSI(1,4) = 0.1666666666666667D0
                            PSI(2,4) = 0.3333333333333333D0
                            PSI(3,4) = 0.1666666666666667D0
                            PSI(1,5) = 0.1666666666666667D0
                            PSI(2,5) = 0.1666666666666667D0
                            PSI(3,5) = 0.3333333333333333D0
                                         !     Weights/6
                            W(1) = - 0.80D0*0.1666666666666667D0
                            WT   =   0.45D0*0.1666666666666667D0
                            W(2) = WT
                            W(3) = WT
                            W(4) = WT
                            W(5) = WT
                        !
                        Case (11)
                              !    Tetrahedral 11-points formula, degree of precision 4
                              !    Gauss Points Intrinsic Coordinates
                              ! See reference: Yen Liu and Marcel Vinokur,
                              !   "Exact Integrations of Polynomials and
                              !    Symmetric Quadrature Formulas over Arbitrary Polyhedral Grids",
                              !    NASA Report > NASA TM-112202; Ames Research Center A-976805
                            PSI(1,1) = 0.25D0
                            PSI(2,1) = 0.25D0
                            PSI(3,1) = 0.25D0
                            !
                                                   ! Vertex based points
                            Alpha = 5.0d0/7.0d0
                            Psi1  = 0.25d0 + 0.75d0 * Alpha
                            Psi2  = 0.25d0 * ( 1.0d0 - Alpha )
                            PSI(1,2) = Psi2
                            PSI(2,2) = Psi2
                            PSI(3,2) = Psi2
                            !
                            PSI(1,3) = Psi1
                            PSI(2,3) = Psi2
                            PSI(3,3) = Psi2
                            !
                            PSI(1,4) = Psi2
                            PSI(2,4) = Psi1
                            PSI(3,4) = Psi2
                            !
                            PSI(1,5) = Psi2
                            PSI(2,5) = Psi2
                            PSI(3,5) = Psi1
                            !
                                                   ! Side based points
                            Beta  = Dsqrt( 70.0d0 ) / 28.0d0
                            Psi1  = 0.25d0 + 0.5d0 * Beta
                            Psi2  = 0.25d0 - 0.5d0 * Beta
                            PSI(1,6) = Psi1
                            PSI(2,6) = Psi1
                            PSI(3,6) = Psi2
                            !
                            PSI(1,7) = Psi1
                            PSI(2,7) = Psi2
                            PSI(3,7) = Psi1
                            !
                            PSI(1,8) = Psi2
                            PSI(2,8) = Psi1
                            PSI(3,8) = Psi1
                            !
                            PSI(1,9) = Psi2
                            PSI(2,9) = Psi2
                            PSI(3,9) = Psi1
                            !
                            PSI(1,10)= Psi2
                            PSI(2,10)= Psi1
                            PSI(3,10)= Psi2
                            !
                            PSI(1,11)= Psi1
                            PSI(2,11)= Psi2
                            PSI(3,11)= Psi2
                            !
                                                   !     Weights/6
                            W(1) = - ( 1.48D2 / 1.875D3 ) * 0.1666666666666667D0
                            WT   =   ( 3.43D2 / 7.5D3   ) * 0.1666666666666667D0
                            W(2) = WT
                            W(3) = WT
                            W(4) = WT
                            W(5) = WT
                            !
                            WTB   =   ( 5.6D1 / 3.75D2   ) * 0.1666666666666667D0
                            W(6 ) = WTB
                            W(7 ) = WTB
                            W(8 ) = WTB
                            W(9 ) = WTB
                            W(10) = WTB
                            W(11) = WTB
                    !
                    EndSelect CaseNGP3d
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                !     gauss points for triangular elements
                Case (2)
                    !
                    CaseNGP2d: Select case (NGP)
                        !
                        Case (1)
                            !     triangle 1-points formula, degree of precision 1
                            !     Gauss Points Intrinsic Coordinates
                            PSI(1,1) = 0.333333333333333D0
                            PSI(2,1) = 0.333333333333333D0
                            !     Weights/2
                            W(1) = 1.0D0 *0.5d0
                        !
                        Case (3)
                            !     triangle 3-points formula, degree of precision 2
                            !     Gauss Points Intrinsic Coordinates
                            PSI(1,1) = 0.166666666666667D0
                            PSI(2,1) = 0.166666666666667D0
                            PSI(1,2) = 0.666666666666667D0
                            PSI(2,2) = 0.166666666666667D0
                            PSI(1,3) = 0.166666666666667D0
                            PSI(2,3) = 0.666666666666667D0
                            !     Weights/2
                            WT   = 0.333333333333333D0 *0.5d0
                            W(1) = WT
                            W(2) = WT
                            W(3) = WT
                        !
                        Case (4)
                            !     triangle 4-points formula, degree of precision 3
                            !     Gauss Points Intrinsic Coordinates
                            PSI(1,1) = 0.2D0
                            PSI(2,1) = 0.2D0
                            PSI(1,2) = 0.6D0
                            PSI(2,2) = 0.2D0
                            PSI(1,3) = 0.2D0
                            PSI(2,3) = 0.6D0
                            PSI(1,4) = 0.333333333333333D0
                            PSI(2,4) = 0.333333333333333D0
                            !     Weights/2
                            WT   = 0.520833333333333D0 *0.5d0
                            W(1) = WT
                            W(2) = WT
                            W(3) = WT
                            W(4) = - 0.56250D0 *0.5d0
                        !
                        Case (6)
                            !     triangle 6-points formula, degree of precision 4
                            !     Gauss Points Intrinsic Coordinates
                            PSI(1,1) = 0.091576213509771D0
                            PSI(2,1) = 0.091576213509771D0
                            PSI(1,2) = 0.445948490915965D0
                            PSI(2,2) = 0.445948490915965D0
                            PSI(1,3) = 0.816847572980459D0
                            PSI(2,3) = 0.091576213509771D0
                            PSI(1,4) = 0.108103018168070D0
                            PSI(2,4) = 0.445948490915965D0
                            PSI(1,5) = 0.091576213509771D0
                            PSI(2,5) = 0.816847572980459D0
                            PSI(1,6) = 0.445948490915965D0
                            PSI(2,6) = 0.108103018168070D0
                            !     Weights/2
                            WT1  = 0.109951743655322D0 *0.5d0
                            WT2  = 0.223381589678011D0 *0.5d0
                            W(1) = WT1
                            W(2) = WT2
                            W(3) = WT1
                            W(4) = WT2
                            W(5) = WT1
                            W(6) = WT2
                        !
                        Case (7)
                            !     triangle 7-points formula, degree of precision 5
                            !     Gauss Points Intrinsic Coordinates
                            PSI(1,1) = 0.333333333333333D0
                            PSI(2,1) = 0.333333333333333D0
                            PSI(1,2) = 0.101286507323456D0
                            PSI(2,2) = 0.101286507323456D0
                            PSI(1,3) = 0.470142064105115D0
                            PSI(2,3) = 0.470142064105115D0
                            PSI(1,4) = 0.797426985353087D0
                            PSI(2,4) = 0.101286507323456D0
                            PSI(1,5) = 0.059715871789770D0
                            PSI(2,5) = 0.470142064105115D0
                            PSI(1,6) = 0.101286507323456D0
                            PSI(2,6) = 0.797426985353087D0
                            PSI(1,7) = 0.470142064105115D0
                            PSI(2,7) = 0.059715871789770D0
                            !     Weights/2
                            W(1) = 0.225000000000000d0 *0.5d0
                            WT1  = 0.125939180544827d0 *0.5d0
                            WT2  = 0.132394152788506d0 *0.5d0
                            W(2) = WT1
                            W(3) = WT2
                            W(4) = WT1
                            W(5) = WT2
                            W(6) = WT1
                            W(7) = WT2
                    !
                    EndSelect CaseNGP2d
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                !Gauss Points for 1D segments as simplexs(gauss coordinates between[0-1])
                !
                Case (1)
                    !     NDim=1
                    !     integrates exactly a polynomial of order 2*NGP-1
                    CaseNGP1D: Select case (NGP)
                        Case (1)
                            PSI(1,1) = 0.5D0
                            W(1) = 1.0d0
                        !
                        Case (2)
                            GaussPoint1 = 1.0d0/Dsqrt(3.0d0)
                            PSI(1,1) = 0.5d0*( 1.0d0 - GaussPoint1 )
                            PSI(1,2) = 0.5d0*( 1.0d0 + GaussPoint1 )
                            W(1) = 0.5d0
                            W(2) = 0.5d0
                        !
                        Case (3)
                            GaussPoint1 = Dsqrt(0.6d0)
                            PSI(1,1) = 0.5d0*( 1.0d0 - GaussPoint1 )
                            PSI(1,2) = 0.5d0*( 1.0d0 + GaussPoint1 )
                            PSI(1,3) = 0.5d0
                            W1 = 2.5d0/9.0d0
                            W(1) = W1
                            W(2) = W1
                            W(3) = 4.0d0/9.0d0
                        !
                        Case (4)
                            R = Dsqrt(1.2d0)
                            GaussPoint1 = Dsqrt( (3.0d0+2.0d0*R) / 7.0d0 )
                            GaussPoint2 = Dsqrt( (3.0d0-2.0d0*R) / 7.0d0 )
                            PSI(1,1) = 0.5d0*( 1.0d0 - GaussPoint1 )
                            PSI(1,2) = 0.5d0*( 1.0d0 + GaussPoint1 )
                            PSI(1,3) = 0.5d0*( 1.0d0 - GaussPoint2 )
                            PSI(1,4) = 0.5d0*( 1.0d0 + GaussPoint2 )
                            W1 = 1.0d0 / (12.0d0*R)
                            W(1) = (0.25d0-W1)
                            W(2) = (0.25d0-W1)
                            W(3) = (0.25d0+W1)
                            W(4) = (0.25d0+W1)
                    !
                    EndSelect CaseNGP1D
            !
            EndSelect CaseNDimSimplex
        !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Case (1)
            !     The element is a Hipercube
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !
            CaseNDimCube: Select case (NDim)
                Case (1)
                    CaseNGPCube: Select case (NGP)
                        !     integrates exactly a polynomial of order 2*NGP-1
                        Case (1)
                            PSI(1,1) = 0.0D0
                            !     Weights
                            W(1) = 2.0d0
                        !
                        Case (2)
                            GaussPoint1 = 1.0d0/Dsqrt(3.0d0)
                            PSI(1,1) = -GaussPoint1
                            PSI(1,2) =  GaussPoint1
                            !     Weights
                            W(1) = 1.0d0
                            W(2) = 1.0d0
                        !
                        Case (3)
                            GaussPoint1 = Dsqrt(0.6d0)
                            PSI(1,1) = -GaussPoint1
                            PSI(1,2) =  GaussPoint1
                            PSI(1,3) =  0.0d0
                            !     Weights
                            W1 = 5.0d0/9.0d0
                            W(1) = W1
                            W(2) = W1
                            W(3) = 8.0d0/9.0d0
                        !
                        Case (4)
                            R = Dsqrt(1.2d0)
                            GaussPoint1 = Dsqrt( (3.0d0+2.0d0*R) / 7.0d0 )
                            GaussPoint2 = Dsqrt( (3.0d0-2.0d0*R) / 7.0d0 )
                            PSI(1,1) = -GaussPoint1
                            PSI(1,2) =  GaussPoint1
                            PSI(1,3) = -GaussPoint2
                            PSI(1,4) =  GaussPoint2
                            !     Weights
                            W1 = 1.0d0 / (6.0d0*R)
                            W(1) = (0.5d0-W1)
                            W(2) = (0.5d0-W1)
                            W(3) = (0.5d0+W1)
                            W(4) = (0.5d0+W1)
                    !
                    EndSelect CaseNGPCube

                case (2)
                    CaseNGPCubeDim2: Select case (NGP)

                        case (1)

                            PSI(1,1) = 0.0d0
                            PSI(2,1) = 0.0d0

                            W(1) = 4.0d0

                        case (4)

                            PSI(1,1) =    5.773502691896257d-01
                            PSI(1,2) =    5.773502691896257d-01
                            PSI(1,3) =   -5.773502691896257d-01
                            PSI(1,4) =   -5.773502691896257d-01
                            PSI(2,1) =    5.773502691896257d-01
                            PSI(2,2) =   -5.773502691896257d-01
                            PSI(2,3) =    5.773502691896257d-01
                            PSI(2,4) =   -5.773502691896257d-01
                            W(1) =    1.0d0
                            W(2) =    1.0d0
                            W(3) =    1.0d0
                            W(4) =    1.0d0



                        case (9)

                            PSI(1,1) =    7.745966692414834d-01
                            PSI(1,2) =    7.745966692414834d-01
                            PSI(1,3) =    7.745966692414834d-01
                            PSI(1,4) =    0.000000000000000d+00
                            PSI(1,5) =    0.000000000000000d+00
                            PSI(1,6) =    0.000000000000000d+00
                            PSI(1,7) =   -7.745966692414834d-01
                            PSI(1,8) =   -7.745966692414834d-01
                            PSI(1,9) =   -7.745966692414834d-01
                            PSI(2,1) =    7.745966692414834d-01
                            PSI(2,2) =    0.000000000000000d+00
                            PSI(2,3) =   -7.745966692414834d-01
                            PSI(2,4) =    7.745966692414834d-01
                            PSI(2,5) =    0.000000000000000d+00
                            PSI(2,6) =   -7.745966692414834d-01
                            PSI(2,7) =    7.745966692414834d-01
                            PSI(2,8) =    0.000000000000000d+00
                            PSI(2,9) =   -7.745966692414834d-01
                            W(1) =    3.086419753086406d-01
                            W(2) =    4.938271604938261d-01
                            W(3) =    3.086419753086406d-01
                            W(4) =    4.938271604938261d-01
                            W(5) =    7.901234567901234d-01
                            W(6) =    4.938271604938261d-01
                            W(7) =    3.086419753086406d-01
                            W(8) =    4.938271604938261d-01
                            W(9) =    3.086419753086406d-01


                    EndSelect CaseNGPCubeDim2




                case (3)
                    CaseNGPCubeDim3: Select case (NGP)

                        case (1)

                            PSI(1,1) = 0.0d0
                            PSI(2,1) = 0.0d0
                            PSI(3,1) = 0.0d0
                            W(1) = 8.0d0

                        case (8)
                            PSI(1,1) =    5.773502691896257d-01
                            PSI(1,2) =    5.773502691896257d-01
                            PSI(1,3) =    5.773502691896257d-01
                            PSI(1,4) =    5.773502691896257d-01
                            PSI(1,5) =   -5.773502691896257d-01
                            PSI(1,6) =   -5.773502691896257d-01
                            PSI(1,7) =   -5.773502691896257d-01
                            PSI(1,8) =   -5.773502691896257d-01
                            PSI(2,1) =    5.773502691896257d-01
                            PSI(2,2) =    5.773502691896257d-01
                            PSI(2,3) =   -5.773502691896257d-01
                            PSI(2,4) =   -5.773502691896257d-01
                            PSI(2,5) =    5.773502691896257d-01
                            PSI(2,6) =    5.773502691896257d-01
                            PSI(2,7) =   -5.773502691896257d-01
                            PSI(2,8) =   -5.773502691896257d-01
                            PSI(3,1) =    5.773502691896257d-01
                            PSI(3,2) =   -5.773502691896257d-01
                            PSI(3,3) =    5.773502691896257d-01
                            PSI(3,4) =   -5.773502691896257d-01
                            PSI(3,5) =    5.773502691896257d-01
                            PSI(3,6) =   -5.773502691896257d-01
                            PSI(3,7) =    5.773502691896257d-01
                            PSI(3,8) =   -5.773502691896257d-01
                            W(1) =    1.0d0
                            W(2) =    1.0d0
                            W(3) =    1.0d0
                            W(4) =    1.0d0
                            W(5) =    1.0d0
                            W(6) =    1.0d0
                            W(7) =    1.0d0
                            W(8) =    1.0d0


                        case (27)

                            PSI(1,1) =    7.745966692414834d-01
                            PSI(1,2) =    7.745966692414834d-01
                            PSI(1,3) =    7.745966692414834d-01
                            PSI(1,4) =    7.745966692414834d-01
                            PSI(1,5) =    7.745966692414834d-01
                            PSI(1,6) =    7.745966692414834d-01
                            PSI(1,7) =    7.745966692414834d-01
                            PSI(1,8) =    7.745966692414834d-01
                            PSI(1,9) =    7.745966692414834d-01
                            PSI(1,10) =    0.000000000000000d+00
                            PSI(1,11) =    0.000000000000000d+00
                            PSI(1,12) =    0.000000000000000d+00
                            PSI(1,13) =    0.000000000000000d+00
                            PSI(1,14) =    0.000000000000000d+00
                            PSI(1,15) =    0.000000000000000d+00
                            PSI(1,16) =    0.000000000000000d+00
                            PSI(1,17) =    0.000000000000000d+00
                            PSI(1,18) =    0.000000000000000d+00
                            PSI(1,19) =   -7.745966692414834d-01
                            PSI(1,20) =   -7.745966692414834d-01
                            PSI(1,21) =   -7.745966692414834d-01
                            PSI(1,22) =   -7.745966692414834d-01
                            PSI(1,23) =   -7.745966692414834d-01
                            PSI(1,24) =   -7.745966692414834d-01
                            PSI(1,25) =   -7.745966692414834d-01
                            PSI(1,26) =   -7.745966692414834d-01
                            PSI(1,27) =   -7.745966692414834d-01
                            PSI(2,1) =    7.745966692414834d-01
                            PSI(2,2) =    7.745966692414834d-01
                            PSI(2,3) =    7.745966692414834d-01
                            PSI(2,4) =    0.000000000000000d+00
                            PSI(2,5) =    0.000000000000000d+00
                            PSI(2,6) =    0.000000000000000d+00
                            PSI(2,7) =   -7.745966692414834d-01
                            PSI(2,8) =   -7.745966692414834d-01
                            PSI(2,9) =   -7.745966692414834d-01
                            PSI(2,10) =    7.745966692414834d-01
                            PSI(2,11) =    7.745966692414834d-01
                            PSI(2,12) =    7.745966692414834d-01
                            PSI(2,13) =    0.000000000000000d+00
                            PSI(2,14) =    0.000000000000000d+00
                            PSI(2,15) =    0.000000000000000d+00
                            PSI(2,16) =   -7.745966692414834d-01
                            PSI(2,17) =   -7.745966692414834d-01
                            PSI(2,18) =   -7.745966692414834d-01
                            PSI(2,19) =    7.745966692414834d-01
                            PSI(2,20) =    7.745966692414834d-01
                            PSI(2,21) =    7.745966692414834d-01
                            PSI(2,22) =    0.000000000000000d+00
                            PSI(2,23) =    0.000000000000000d+00
                            PSI(2,24) =    0.000000000000000d+00
                            PSI(2,25) =   -7.745966692414834d-01
                            PSI(2,26) =   -7.745966692414834d-01
                            PSI(2,27) =   -7.745966692414834d-01
                            PSI(3,1) =    7.745966692414834d-01
                            PSI(3,2) =    0.000000000000000d+00
                            PSI(3,3) =   -7.745966692414834d-01
                            PSI(3,4) =    7.745966692414834d-01
                            PSI(3,5) =    0.000000000000000d+00
                            PSI(3,6) =   -7.745966692414834d-01
                            PSI(3,7) =    7.745966692414834d-01
                            PSI(3,8) =    0.000000000000000d+00
                            PSI(3,9) =   -7.745966692414834d-01
                            PSI(3,10) =    7.745966692414834d-01
                            PSI(3,11) =    0.000000000000000d+00
                            PSI(3,12) =   -7.745966692414834d-01
                            PSI(3,13) =    7.745966692414834d-01
                            PSI(3,14) =    0.000000000000000d+00
                            PSI(3,15) =   -7.745966692414834d-01
                            PSI(3,16) =    7.745966692414834d-01
                            PSI(3,17) =    0.000000000000000d+00
                            PSI(3,18) =   -7.745966692414834d-01
                            PSI(3,19) =    7.745966692414834d-01
                            PSI(3,20) =    0.000000000000000d+00
                            PSI(3,21) =   -7.745966692414834d-01
                            PSI(3,22) =    7.745966692414834d-01
                            PSI(3,23) =    0.000000000000000d+00
                            PSI(3,24) =   -7.745966692414834d-01
                            PSI(3,25) =    7.745966692414834d-01
                            PSI(3,26) =    0.000000000000000d+00
                            PSI(3,27) =   -7.745966692414834d-01
                            W(1) =    1.714677640603555d-01
                            W(2) =    2.743484224965694d-01
                            W(3) =    1.714677640603555d-01
                            W(4) =    2.743484224965694d-01
                            W(5) =    4.389574759945121d-01
                            W(6) =    2.743484224965694d-01
                            W(7) =    1.714677640603555d-01
                            W(8) =    2.743484224965694d-01
                            W(9) =    1.714677640603555d-01
                            W(10) =    2.743484224965694d-01
                            W(11) =    4.389574759945121d-01
                            W(12) =    2.743484224965694d-01
                            W(13) =    4.389574759945121d-01
                            W(14) =    7.023319615912208d-01
                            W(15) =    4.389574759945121d-01
                            W(16) =    2.743484224965694d-01
                            W(17) =    4.389574759945121d-01
                            W(18) =    2.743484224965694d-01
                            W(19) =    1.714677640603555d-01
                            W(20) =    2.743484224965694d-01
                            W(21) =    1.714677640603555d-01
                            W(22) =    2.743484224965694d-01
                            W(23) =    4.389574759945121d-01
                            W(24) =    2.743484224965694d-01
                            W(25) =    1.714677640603555d-01
                            W(26) =    2.743484224965694d-01
                            W(27) =    1.714677640603555d-01



                    EndSelect CaseNGPCubeDim3

            !
            EndSelect CaseNDimCube
    !
    EndSelect CaseSimplex
    !
    Return
end Subroutine GaussRule

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine LocalShapeDerQ (NodEL,NDim,dF,GaussPoint,Iorder,iS)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Implicit real*8 (a-h,o-z)
    Dimension dF(NDim,NodEL),GaussPoint(NDim)
    Integer d(NDim)
    ! Ejemplo: ndim=3, d(3) vector indicador de la derivada de primer orden = d/dy-> [0,1,0]
    !    iS = 1 Serendipity element, iS=0 Lagrangian element
    d=0

    Do Nod=1,NodEL
        Do iLC=1,NDim
            d(iLC)=1
            dF(iLC,Nod)= &
                &            SFQuadrilateral(GaussPoint,d,Iorder,Nod,NDim,iS)
            d(iLC)=0
        Enddo
    Enddo

    Return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine ShapeFQ (NodEL,NDim,F,GaussPoint,iOrder,iS)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Implicit real*8 (a-h,o-z)
    Dimension F(NodEL),GaussPoint(NDim)
    !    iS = 1 Serendipity element, iS=0 Lagrangian element
    !    iOrder = Interpolation order: 1=Linear; 2=Cuadratic
    Integer d(NDim)

    d=0
    Do Nod=1,NodEL
        F(Nod)= SFQuadrilateral(GaussPoint,d,iOrder,Nod,NDim,iS)
    Enddo
   
    Return
End
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
