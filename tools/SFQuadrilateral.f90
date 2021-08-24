!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!       A function that calculates the values of shape functions
!       at specified points for Quadrilateral Elements
!
!
!_______________________________________________________________________
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Real*8 Function SFQuadrilateral (PsiCor,d,Poly_Order,node,NDim,iS)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !-----------------------------------------------------------------------
    Implicit real*8 (a-h,o-z)
    Integer Poly_Order,d(NDim),dX(1),dY(1),dZ(1),Deriv
    Dimension Sfn(ndim+1),PsiCor(ndim),PsiCorX(1),PsiCorY(1),PsiCorZ(1)
    !     NDim = dimension of the quadrilateral
    !     NDim= 1 Not considered, resort to Simplex Elements
    !     NDim= 2 Square
    !     NDim= 3 Cube
    !     Poly_Order: complete Polynomial degree of the shape function
    !        iS=0: family of Lagrangian elements (valid for Poly_Order>1)
    !        NDim=2 ->  9-node element
    !        NDim=3 -> 27-node element
    !        iS=1: family of Serendipity elements (valid for Poly_Order>1)
    !        NDim=2 ->  8-node element
    !        NDim=3 -> 20-node element
    !     node: number (ordinal) of the shape function)
    !     PsiCor(ndim)= Simplex intrinsec coordinates
    !     Deriv: order  of Partial DERIVATIVE (Function=:0)
    !     d(1): order of derivation with respect to Psi1
    !     d(2): order of derivation with respect to Psi2
    !     d(3): order of derivation with respect to Psi3
    !     Sfn(i): Linear Shape Function  corresponding to Node i

    SFQuadrilateral = 0.0D0
    Xi = 0.0D0
    Yi = 0.0D0
    Zi = 0.0D0
    Deriv = sum(d)
    Select Case (NDim)
           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Case (3)
            ! Cube Element
            X  = PsiCor(1)
            Y  = PsiCor(2)
            Z  = PsiCor(3)
            XX = X*X
            YY = Y*Y
            ZZ = Z*Z

            ! Xi: Coordenada X del nodo correspondiente a la funcion de forma Ni
            ! Yi: Coordenada Y del nodo correspondiente a la funcion de forma Ni
            ! Zi: Coordenada Z del nodo correspondiente a la funcion de forma Ni

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Select Case (Poly_Order)
                   ! Linear:1, Quadratic:2, Cubic:3, Etc.
                Case (1)
                    ! Poly_Order Linear
                    Select Case (node)
                        Case(1)
                            Xi = -1.0d0
                            Yi = -1.0d0
                            Zi = -1.0d0
                        Case(2)
                            Xi =  1.0d0
                            Yi = -1.0d0
                            Zi = -1.0d0
                        Case(3)
                            Xi =  1.0d0
                            Yi =  1.0d0
                            Zi = -1.0d0
                        Case(4)
                            Xi = -1.0d0
                            Yi =  1.0d0
                            Zi = -1.0d0
                        Case(5)
                            Xi = -1.0d0
                            Yi = -1.0d0
                            Zi =  1.0d0
                        Case(6)
                            Xi =  1.0d0
                            Yi = -1.0d0
                            Zi =  1.0d0
                        Case(7)
                            Xi =  1.0d0
                            Yi =  1.0d0
                            Zi =  1.0d0
                        Case(8)
                            Xi = -1.0d0
                            Yi =  1.0d0
                            Zi =  1.0d0
                    End Select
                    Select Case (Deriv)
                        Case (0)
                            ! Function Values
                            SFQuadrilateral = 0.125d0*(1.0d0+Xi*X)*(1.0d0+Yi*Y)*(1.0d0+Zi*Z)
                        Case (1)
                            ! First Order Derivative
                            If(d(1).eq.1)Then                                          ! d/dx
                                SFQuadrilateral = 0.125d0*Xi*(1.0d0+Yi*Y)*(1.0d0+Zi*Z)
                            ElseIf(d(2).eq.1)Then                                      ! d/dy
                                SFQuadrilateral = 0.125d0*(1.0d0+Xi*X)*Yi*(1.0d0+Zi*Z)
                            ElseIf(d(3).eq.1)Then                                      ! d/dz
                                SFQuadrilateral = 0.125d0*(1.0d0+Xi*X)*(1.0d0+Yi*Y)*Zi
                            EndIf
                        Case (2)
                            ! Second Order Derivative
                            If(d(1).eq.2)Then                                          ! d2/dx2
                                SFQuadrilateral = 0.0d0
                            ElseIf(d(2).eq.2)Then                                      ! d2/dy2
                                SFQuadrilateral = 0.0d0
                            ElseIf(d(3).eq.2)Then                                      ! d2/dz2
                                SFQuadrilateral = 0.0d0
                            ElseIf(d(1)+d(2).eq.2.and.d(1).ne.2.and.d(2).ne.2)Then     ! d2/dxdy
                                SFQuadrilateral = 0.125d0*Xi*Yi*(1.0d0+Zi*Z)
                            ElseIf(d(1)+d(3).eq.2.and.d(1).ne.2.and.d(3).ne.2)Then     ! d2/dxdz
                                SFQuadrilateral = 0.125d0*Xi*(1.0d0+Yi*Y)*Zi
                            ElseIf(d(2)+d(3).eq.2.and.d(2).ne.2.and.d(3).ne.2)Then     ! d2/dydz
                                SFQuadrilateral = 0.125d0*(1.0d0+Xi*X)*Yi*Zi
                            EndIf
                    End Select

                Case (2)
                    ! Poly_Order Cuadratic
                    If(iS.eq.0)Then
                        !----------------------------------------------------------------------
                        ! Lagrangian Element (27-node element)
                        !----------------------------------------------------------------------

                        ! Construido por multiplicacion de las funciones de forma del 1D
                        ! fX: funcion de variable X
                        ! gY: funcion de variable Y
                        ! hZ: funcion de variable Z

                        NDim1D = 1
                        iOrder = 2
                        PsiCorX(1) = X
                        PsiCorY(1) = Y
                        PsiCorZ(1) = Z
                        Select Case (node)
                            Case(1)
                                nodeX = 1
                                nodeY = 1
                                nodeZ = 1
                            Case(2)
                                nodeX = 2
                                nodeY = 1
                                nodeZ = 1
                            Case(3)
                                nodeX = 2
                                nodeY = 2
                                nodeZ = 1
                            Case(4)
                                nodeX = 1
                                nodeY = 2
                                nodeZ = 1
                            Case(5)
                                nodeX = 1
                                nodeY = 1
                                nodeZ = 2
                            Case(6)
                                nodeX = 2
                                nodeY = 1
                                nodeZ = 2
                            Case(7)
                                nodeX = 2
                                nodeY = 2
                                nodeZ = 2
                            Case(8)
                                nodeX = 1
                                nodeY = 2
                                nodeZ = 2
                            Case(9)
                                nodeX = 3
                                nodeY = 1
                                nodeZ = 1
                            Case(10)
                                nodeX = 2
                                nodeY = 3
                                nodeZ = 1
                            Case(11)
                                nodeX = 3
                                nodeY = 2
                                nodeZ = 1
                            Case(12)
                                nodeX = 1
                                nodeY = 3
                                nodeZ = 1
                            Case(13)
                                nodeX = 1
                                nodeY = 1
                                nodeZ = 3
                            Case(14)
                                nodeX = 2
                                nodeY = 1
                                nodeZ = 3
                            Case(15)
                                nodeX = 2
                                nodeY = 2
                                nodeZ = 3
                            Case(16)
                                nodeX = 1
                                nodeY = 2
                                nodeZ = 3
                            Case(17)
                                nodeX = 3
                                nodeY = 1
                                nodeZ = 2
                            Case(18)
                                nodeX = 2
                                nodeY = 3
                                nodeZ = 2
                            Case(19)
                                nodeX = 3
                                nodeY = 2
                                nodeZ = 2
                            Case(20)
                                nodeX = 1
                                nodeY = 3
                                nodeZ = 2
                            Case(21)
                                nodeX = 3
                                nodeY = 3
                                nodeZ = 1
                            Case(22)
                                nodeX = 3
                                nodeY = 1
                                nodeZ = 3
                            Case(23)
                                nodeX = 2
                                nodeY = 3
                                nodeZ = 3
                            Case(24)
                                nodeX = 3
                                nodeY = 2
                                nodeZ = 3
                            Case(25)
                                nodeX = 1
                                nodeY = 3
                                nodeZ = 3
                            Case(26)
                                nodeX = 3
                                nodeY = 3
                                nodeZ = 2
                            Case(27)
                                nodeX = 3
                                nodeY = 3
                                nodeZ = 3
                        End Select
                        Select Case (Deriv)
                            Case(0)
                                ! Function Values
                                dX(1) = 0
                                dY(1) = 0
                                dZ(1) = 0
                            Case(1)
                                ! First Order Derivative
                                If(d(1).eq.1)Then                ! d/dx
                                    dX(1) = 1
                                    dY(1) = 0
                                    dZ(1) = 0
                                ElseIf(d(2).eq.1)Then            ! d/dy
                                    dX(1) = 0
                                    dY(1) = 1
                                    dZ(1) = 0
                                Else                             ! d/dz
                                    dX(1) = 0
                                    dY(1) = 0
                                    dZ(1) = 1
                                EndIf
                            Case(2)
                                ! Second Order Derivative
                                If(d(1).eq.2)Then                ! d2/dx2
                                    dX(1) = 2
                                    dY(1) = 0
                                    dZ(1) = 0
                                ElseIf(d(2).eq.2)Then            ! d2/dy2
                                    dX(1) = 0
                                    dY(1) = 2
                                    dZ(1) = 0
                                ElseIf(d(3).eq.2)Then            ! d2/dz2
                                    dX(1) = 0
                                    dY(1) = 0
                                    dZ(1) = 2
                                ElseIf(d(1)+d(2).eq.2)Then       ! d2/dxdy
                                    dX(1) = 1
                                    dY(1) = 1
                                    dZ(1) = 0
                                ElseIf(d(1)+d(3).eq.2)Then       ! d2/dxdz
                                    dX(1) = 1
                                    dY(1) = 0
                                    dZ(1) = 1
                                Else                             ! d2/dydz
                                    dX(1) = 0
                                    dY(1) = 1
                                    dZ(1) = 1
                                EndIf
                        End Select
                        fX = SF1D(PsiCorX,dX,iOrder,nodeX,NDim1D,0)
                        gY = SF1D(PsiCorY,dY,iOrder,nodeY,NDim1D,0)
                        hZ = SF1D(PsiCorZ,dZ,iOrder,nodeZ,NDim1D,0)
                        SFQuadrilateral = fX * gY * hZ
                       !----------------------------------------------------------------------
                    ElseIf(iS.eq.1)Then
                        !----------------------------------------------------------------------
                        ! Serendipity Element (20-node element)
                        !----------------------------------------------------------------------
                        Select Case (node)
                            Case(1)
                                Xi = -1.0d0
                                Yi = -1.0d0
                                Zi = -1.0d0
                            Case(2)
                                Xi =  1.0d0
                                Yi = -1.0d0
                                Zi = -1.0d0
                            Case(3)
                                Xi =  1.0d0
                                Yi =  1.0d0
                                Zi = -1.0d0
                            Case(4)
                                Xi = -1.0d0
                                Yi =  1.0d0
                                Zi = -1.0d0
                            Case(5)
                                Xi = -1.0d0
                                Yi = -1.0d0
                                Zi =  1.0d0
                            Case(6)
                                Xi =  1.0d0
                                Yi = -1.0d0
                                Zi =  1.0d0
                            Case(7)
                                Xi =  1.0d0
                                Yi =  1.0d0
                                Zi =  1.0d0
                            Case(8)
                                Xi = -1.0d0
                                Yi =  1.0d0
                                Zi =  1.0d0
                            Case(9)
                                Xi =  0.0d0
                                Yi = -1.0d0
                                Zi = -1.0d0
                            Case(10)
                                Xi =  1.0d0
                                Yi =  0.0d0
                                Zi = -1.0d0
                            Case(11)
                                Xi =  0.0d0
                                Yi =  1.0d0
                                Zi = -1.0d0
                            Case(12)
                                Xi = -1.0d0
                                Yi =  0.0d0
                                Zi = -1.0d0
                            Case(13)
                                Xi = -1.0d0
                                Yi = -1.0d0
                                Zi =  0.0d0
                            Case(14)
                                Xi =  1.0d0
                                Yi = -1.0d0
                                Zi =  0.0d0
                            Case(15)
                                Xi =  1.0d0
                                Yi =  1.0d0
                                Zi =  0.0d0
                            Case(16)
                                Xi = -1.0d0
                                Yi =  1.0d0
                                Zi =  0.0d0
                            Case(17)
                                Xi =  0.0d0
                                Yi = -1.0d0
                                Zi =  1.0d0
                            Case(18)
                                Xi =  1.0d0
                                Yi =  0.0d0
                                Zi =  1.0d0
                            Case(19)
                                Xi =  0.0d0
                                Yi =  1.0d0
                                Zi =  1.0d0
                            Case(20)
                                Xi = -1.0d0
                                Yi =  0.0d0
                                Zi =  1.0d0
                        End Select
                        Select Case (Deriv)
                            Case (0)
                                ! Function Values
                                Select Case (node)
                                    Case (1:8)
                                        SFQuadrilateral = 0.125d0*(1.0d0+Xi*X)*(1.0d0+Yi*Y)*(1.0d0+Zi*Z)*(Xi*X+Yi*Y+Zi*Z-2)
                                    Case (9,11,17,19)
                                        SFQuadrilateral = 0.25d0*(1.0d0-XX)*(1.0d0+Yi*Y)*(1.0d0+Zi*Z)
                                    Case (10,12,18,20)
                                        SFQuadrilateral = 0.25d0*(1.0d0+Xi*X)*(1.0d0-YY)*(1.0d0+Zi*Z)
                                    Case (13:16)
                                        SFQuadrilateral = 0.25d0*(1.0d0+Xi*X)*(1.0d0+Yi*Y)*(1.0d0-ZZ)
                                End Select
                            Case (1)
                                ! First Order Derivative
                                Select Case (node)
                                    Case (1:8)
                                        If(d(1).eq.1)Then                                 ! d/dx
                                            SFQuadrilateral = 0.125d0*Xi*(1.0d0+Yi*Y)*(1.0d0+Zi*Z)*(2.0d0*Xi*X+Yi*Y+Zi*Z-1.0D0)
                                        ElseIf(d(2).eq.1)Then                             ! d/dy
                                            SFQuadrilateral = 0.125d0*(1.0d0+Xi*X)*Yi*(1.0d0+Zi*Z)*(Xi*X+2.0d0*Yi*Y+Zi*Z-1.0D0)
                                        Else                                              ! d/dz
                                            SFQuadrilateral = 0.125d0*(1.0d0+Xi*X)*(1.0d0+Yi*Y)*Zi*(Xi*X+Yi*Y+2.0d0*Zi*Z-1.0D0)
                                        EndIf
                                    Case (9,11,17,19)
                                        If(d(1).eq.1)Then                                 ! d/dx
                                            SFQuadrilateral = -0.50d0*X*(1.0d0+Yi*Y)*(1.0d0+Zi*Z)
                                        ElseIf(d(2).eq.1)Then                             ! d/dy
                                            SFQuadrilateral = 0.25d0*(1.0d0-XX)*Yi*(1.0d0+Zi*Z)
                                        Else                                              ! d/dz
                                            SFQuadrilateral = 0.25d0*(1.0d0-XX)*(1.0d0+Yi*Y)*Zi
                                        EndIf
                                    Case (10,12,18,20)
                                        If(d(1).eq.1)Then                                 ! d/dx
                                            SFQuadrilateral = 0.25d0*Xi*(1.0d0-YY)*(1.0d0+Zi*Z)
                                        ElseIf(d(2).eq.1)Then                             ! d/dy
                                            SFQuadrilateral = -0.50d0*(1.0d0+Xi*X)*Y*(1.0d0+Zi*Z)
                                        Else                                              ! d/dz
                                            SFQuadrilateral = 0.25d0*(1.0d0+Xi*X)*(1.0d0-YY)*Zi
                                        EndIf
                                    Case (13:16)
                                        If(d(1).eq.1)Then                                 ! d/dx
                                            SFQuadrilateral = 0.25d0*Xi*(1.0d0+Yi*Y)*(1.0d0-ZZ)
                                        ElseIf(d(2).eq.1)Then                             ! d/dy
                                            SFQuadrilateral = 0.25d0*(1.0d0+Xi*X)*Yi*(1.0d0-ZZ)
                                        Else                                              ! d/dz
                                            SFQuadrilateral = -0.50d0*(1.0d0+Xi*X)*(1.0d0+Yi*Y)*Z
                                        EndIf
                                End Select
                            Case (2)
                                ! Second Order Derivative
                                Select Case (node)
                                    Case (1:8)
                                        If(d(1).eq.2)Then                                 ! d2/dx2
                                            SFQuadrilateral = 0.25d0*(1.0d0+Yi*Y)*(1.0d0+Zi*Z)
                                        ElseIf(d(2).eq.2)Then                             ! d2/dy2
                                            SFQuadrilateral = 0.25d0*(1.0d0+Xi*X)*(1.0d0+Zi*Z)
                                        ElseIf(d(3).eq.2)Then                             ! d2/dz2
                                            SFQuadrilateral = 0.25d0*(1.0d0+Xi*X)*(1.0d0+Yi*Y)
                                        ElseIf(d(1)+d(2).eq.2)Then                        ! d2/dxdy
                                            SFQuadrilateral = 0.125d0*Xi*Yi*(1.0d0+Zi*Z)*(2.0d0*Xi*X+2.0d0*Yi*Y+Zi*Z)
                                        ElseIf(d(1)+d(3).eq.2)Then                        ! d2/dxdz
                                            SFQuadrilateral = 0.125d0*Xi*Zi*(1.0d0+Yi*Y)*(2.0d0*Xi*X+Yi*Y+2.0d0*Zi*Z)
                                        Else                                              ! d2/dydz
                                            SFQuadrilateral = 0.125d0*Yi*Zi*(1.0d0+Xi*X)*(Xi*X+2.0d0*Yi*Y+2.0d0*Zi*Z)
                                        EndIf
                                    Case (9,11,17,19)
                                        If(d(1).eq.2)Then                                 ! d2/dx2
                                            SFQuadrilateral = -0.50d0*(1.0d0+Yi*Y)*(1.0d0+Zi*Z)
                                        ElseIf(d(2).eq.2)Then                             ! d2/dy2
                                            SFQuadrilateral = 0.0d0
                                        ElseIf(d(3).eq.2)Then                             ! d2/dz2
                                            SFQuadrilateral = 0.0d0
                                        ElseIf(d(1)+d(2).eq.2)Then                        ! d2/dxdy
                                            SFQuadrilateral = -0.50d0*X*Yi*(1.0d0+Zi*Z)
                                        ElseIf(d(1)+d(3).eq.2)Then                        ! d2/dxdz
                                            SFQuadrilateral = -0.50d0*X*(1.0d0+Yi*Y)*Zi
                                        Else                                              ! d2/dydz
                                            SFQuadrilateral = 0.25d0*(1.0d0-XX)*Yi*Zi
                                        EndIf
                                    Case (10,12,18,20)
                                        If(d(1).eq.2)Then                                 ! d2/dx2
                                            SFQuadrilateral = 0.0d0
                                        ElseIf(d(2).eq.2)Then                             ! d2/dy2
                                            SFQuadrilateral = -0.50d0*(1.0d0+Xi*X)*(1.0d0+Zi*Z)
                                        ElseIf(d(3).eq.2)Then                             ! d2/dz2
                                            SFQuadrilateral = 0.0d0
                                        ElseIf(d(1)+d(2).eq.2)Then                        ! d2/dxdy
                                            SFQuadrilateral = -0.50d0*Xi*Y*(1.0d0+Zi*Z)
                                        ElseIf(d(1)+d(3).eq.2)Then                        ! d2/dxdz
                                            SFQuadrilateral = 0.25d0*Xi*(1.0d0-YY)*Zi
                                        Else                                              ! d2/dydz
                                            SFQuadrilateral = -0.50d0*(1.0d0+Xi*X)*Y*Zi
                                        EndIf
                                    Case (13:16)
                                        If(d(1).eq.2)Then                                 ! d2/dx2
                                            SFQuadrilateral = 0.0d0
                                        ElseIf(d(2).eq.2)Then                             ! d2/dy2
                                            SFQuadrilateral = 0.0d0
                                        ElseIf(d(3).eq.2)Then                             ! d2/dz2
                                            SFQuadrilateral = -0.50d0*(1.0d0+Xi*X)*(1.0d0+Yi*Y)
                                        ElseIf(d(1)+d(2).eq.2)Then                        ! d2/dxdy
                                            SFQuadrilateral = 0.25d0*Xi*Yi*(1.0d0-ZZ)
                                        ElseIf(d(1)+d(3).eq.2)Then                        ! d2/dxdz
                                            SFQuadrilateral = -0.50d0*Xi*(1.0d0+Yi*Y)*Z
                                        Else                                              ! d2/dydz
                                            SFQuadrilateral = -0.50d0*(1.0d0+Xi*X)*Yi*Z
                                        EndIf
                                End Select
                        End Select
                    EndIf

            End Select

           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Case (2)
            ! Square Element
            X  = PsiCor(1)
            Y  = PsiCor(2)
            XX = X*X
            YY = Y*Y

            ! Xi: Coordenada X del nodo correspondiente a la funcion de forma Ni
            ! Yi: Coordenada Y del nodo correspondiente a la funcion de forma Ni

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Select Case (Poly_Order)
                   ! Linear:1, Quadratic:2, Cubic:3, Etc.
                Case (1)
                    ! Poly_Order Linear
                    Select Case (node)
                        Case(1)
                            Xi = -1.0d0
                            Yi = -1.0d0
                        Case(2)
                            Xi =  1.0d0
                            Yi = -1.0d0
                        Case(3)
                            Xi =  1.0d0
                            Yi =  1.0d0
                        Case(4)
                            Xi = -1.0d0
                            Yi =  1.0d0
                    End Select
                    Select Case (Deriv)
                        Case (0)
                            ! Function Values
                            SFQuadrilateral = 0.25d0*(1.0d0+Xi*X)*(1.0d0+Yi*Y)
                        Case (1)
                            ! First Order Derivative
                            If(d(1).eq.1)Then                                     ! d/dx
                                SFQuadrilateral = 0.25d0*Xi*(1.0d0+Yi*Y)
                            Else                                                  ! d/dy
                                SFQuadrilateral = 0.25d0*(1.0d0+Xi*X)*Yi
                            EndIf
                        Case (2)
                            ! Second Order Derivative
                            If(d(1).eq.2)Then                                     ! d2/dx2
                                SFQuadrilateral = 0.0d0
                            ElseIf(d(2).eq.2)Then                                 ! d2/dy2
                                SFQuadrilateral = 0.0d0
                            Else                                                  ! d2/dxdy
                                SFQuadrilateral = 0.25d0*Xi*Yi
                            EndIf
                    End Select

                Case (2)
                    ! Poly_Order Cuadratic
                    If(iS.eq.0)Then
                        !----------------------------------------------------------------------
                        ! Lagrangian Element
                        !----------------------------------------------------------------------

                        ! Construido por multiplicacion de las funciones de forma del 1D
                        ! fX: funcion de variable X
                        ! gY: funcion de variable Y

                        NDim1D = 1
                        iOrder = 2
                        PsiCorX(1) = X
                        PsiCorY(1) = Y
                        Select Case (node)
                            Case(1)
                                nodeX = 1
                                nodeY = 1
                            Case(2)
                                nodeX = 2
                                nodeY = 1
                            Case(3)
                                nodeX = 2
                                nodeY = 2
                            Case(4)
                                nodeX = 1
                                nodeY = 2
                            Case(5)
                                nodeX = 3
                                nodeY = 1
                            Case(6)
                                nodeX = 2
                                nodeY = 3
                            Case(7)
                                nodeX = 3
                                nodeY = 2
                            Case(8)
                                nodeX = 1
                                nodeY = 3
                            Case(9)
                                nodeX = 3
                                nodeY = 3
                        End Select
                        Select Case (Deriv)
                            Case(0)
                                ! Function Values
                                dX(1) = 0
                                dY(1) = 0
                            Case(1)
                                ! First Order Derivative
                                If(d(1).eq.1)Then                           ! d/dx
                                    dX(1) = 1
                                    dY(1) = 0
                                Else                                        ! d/dy
                                    dX(1) = 0
                                    dY(1) = 1
                                EndIf
                            Case(2)
                                ! Second Order Derivative
                                If(d(1).eq.2)Then                           ! d2/dx2
                                    dX(1) = 2
                                    dY(1) = 0
                                ElseIf(d(2).eq.2)Then                       ! d2/dy2
                                    dX(1) = 0
                                    dY(1) = 2
                                Else                                        ! d2/dxdy
                                    dX(1) = 1
                                    dY(1) = 1
                                EndIf
                        End Select
                        fX = SF1D(PsiCorX,dX,iOrder,nodeX,NDim1D,0)
                        gY = SF1D(PsiCorY,dY,iOrder,nodeY,NDim1D,0)
                        SFQuadrilateral = fX * gY
                       !----------------------------------------------------------------------
                    ElseIf(iS.eq.1.or.iS.eq.2)Then
                        !----------------------------------------------------------------------
                        ! Serendipity or less than Serendipity Element (coincide for NDim=2)
                        !----------------------------------------------------------------------
                        Select Case (node)
                            Case(1)
                                Xi = -1.0d0
                                Yi = -1.0d0
                            Case(2)
                                Xi =  1.0d0
                                Yi = -1.0d0
                            Case(3)
                                Xi =  1.0d0
                                Yi =  1.0d0
                            Case(4)
                                Xi = -1.0d0
                                Yi =  1.0d0
                            Case(5)
                                Xi =  0.0d0
                                Yi = -1.0d0
                            Case(6)
                                Xi =  1.0d0
                                Yi =  0.0d0
                            Case(7)
                                Xi =  0.0d0
                                Yi =  1.0d0
                            Case(8)
                                Xi = -1.0d0
                                Yi =  0.0d0
                        End Select
                        Select Case (Deriv)
                            Case (0)
                                ! Function Values
                                Select Case (node)
                                    Case (1:4)
                                        SFQuadrilateral = 0.25d0*(1.0d0+Xi*X)*(1.0d0+Yi*Y)*(Xi*X+Yi*Y-1)
                                    Case (5,7)
                                        SFQuadrilateral = 0.50d0*(1.0d0-XX)*(1.0d0+Yi*Y)
                                    Case (6,8)
                                        SFQuadrilateral = 0.50d0*(1.0d0+Xi*X)*(1.0d0-YY)
                                End Select
                            Case (1)
                                ! First Order Derivative
                                Select Case (node)
                                    Case (1:4)
                                        If(d(1).eq.1)Then                             ! d/dx
                                            SFQuadrilateral = 0.25d0*Xi*(1.0d0+Yi*Y)*(2.0d0*Xi*X+Yi*Y)
                                        Else                                          ! d/dy
                                            SFQuadrilateral = 0.25d0*Yi*(1.0d0+Xi*X)*(2.0d0*Yi*Y+Xi*X)
                                        EndIf
                                    Case (5,7)
                                        If(d(1).eq.1)Then                             ! d/dx
                                            SFQuadrilateral = -X*(1.0d0+Yi*Y)
                                        Else                                          ! d/dy
                                            SFQuadrilateral = 0.50d0*(1-XX)*Yi
                                        EndIf
                                    Case (6,8)
                                        If(d(1).eq.1)Then                             ! d/dx
                                            SFQuadrilateral = 0.50d0*Xi*(1-YY)
                                        Else                                          ! d/dy
                                            SFQuadrilateral = -(1.0d0+Xi*X)*Y
                                        EndIf
                                End Select
                            Case (2)
                                ! Second Order Derivative
                                Select Case (node)
                                    Case (1,4)
                                        If(d(1).eq.2)Then                              ! d2/dx2
                                            SFQuadrilateral = 0.50d0*(1.0D0+Yi*Y)
                                        ElseIf(d(2).eq.2)Then                          ! d2/dy2
                                            SFQuadrilateral = 0.50d0*(1.0D0+Xi*X)
                                        Else                                           ! d2/dxdy
                                            SFQuadrilateral = 0.25d0*(1.0d0+2.0d0*Xi*X+2.0d0*Yi*Y)
                                        EndIf
                                    Case (5,7)
                                        If(d(1).eq.2)Then                              ! d2/dx2
                                            SFQuadrilateral = -(1.0d0+Yi*Y)
                                        ElseIf(d(2).eq.2)Then                          ! d2/dy2
                                            SFQuadrilateral = 0.0d0
                                        Else                                           ! d2/dxdy
                                            SFQuadrilateral = -X*Yi
                                        EndIf
                                    Case (6,8)
                                        If(d(1).eq.2)Then                              ! d2/dx2
                                            SFQuadrilateral = 0.0d0
                                        ElseIf(d(2).eq.2)Then                          ! d2/dy2
                                            SFQuadrilateral = -(1.0d0+Xi*X)
                                        Else                                           ! d2/dxdy
                                            SFQuadrilateral = -Xi*Y
                                        EndIf
                                End Select
                        End Select
                       !----------------------------------------------------------------------
                    EndIf

            End Select

           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Case (1)
            ! Line Element
            X  = PsiCor(1)
            XX = X*X

            ! Xi: Coordenada X del nodo correspondiente a la funcion de forma Ni

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Select Case (Poly_Order)
                   ! Linear:1, Quadratic:2, Cubic:3, Etc.
                Case (1)
                    ! Poly_Order Linear
                    Select Case (Deriv)
                        Case (0)
                            ! Function Values
                            Select Case (node)
                                Case(1)
                                    SFQuadrilateral = 0.50d0*(1.0d0-X)
                                Case(2)
                                    SFQuadrilateral = 0.50d0*(1.0d0+X)
                            End Select
                        Case (1)
                            ! First Order Derivative
                            Select Case (node)
                                Case(1)
                                    SFQuadrilateral = -0.50d0
                                Case(2)
                                    SFQuadrilateral = 0.50d0
                            End Select
                        Case (2)
                            ! Second Order Derivative
                            Select Case (node)
                                Case(1)
                                    SFQuadrilateral = 0.0d0
                                Case(2)
                                    SFQuadrilateral = 0.0d0
                            End Select
                    End Select

                Case (2)
                    ! Poly_Order Cuadratic
                    Select Case (Deriv)
                        Case (0)
                            ! Function Values
                            Select Case (node)
                                Case(1)
                                    SFQuadrilateral = -0.50d0*X*(1.0d0-X)
                                Case(2)
                                    SFQuadrilateral = 0.50d0*X*(1.0d0+X)
                                Case(3)
                                    SFQuadrilateral = 1.0d0-XX
                            End Select
                        Case (1)
                            ! First Order Derivative
                            Select Case (node)
                                Case(1)
                                    SFQuadrilateral = -0.50d0*(1.0d0-2.0d0*X)
                                Case(2)
                                    SFQuadrilateral = 0.50d0*(1.0d0+2.0d0*X)
                                Case(3)
                                    SFQuadrilateral = -2.0d0*X
                            End Select
                        Case (2)
                            ! Second Order Derivative
                            Select Case (node)
                                Case(1)
                                    SFQuadrilateral = 1.0d0
                                Case(2)
                                    SFQuadrilateral = 1.0d0
                                Case(3)
                                    SFQuadrilateral = -2.0d0
                            End Select
                    End Select

            End Select
    End Select

    Return
End Function SFQuadrilateral
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!       A function that calculates the values of shape functions
!       at specified points for line Elements in (-1,1) for building
!                quadrilateral Lagrangian elements
!
!_______________________________________________________________________
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Real*8 Function SF1D (PsiCor,d,Poly_Order,node,NDim,iS)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !-----------------------------------------------------------------------
    Implicit real*8 (a-h,o-z)
    Integer Poly_Order,d(NDim),dX(1),dY(1),dZ(1),Deriv
    Dimension Sfn(ndim+1),PsiCor(ndim),PsiCorX(1),PsiCorY(1),PsiCorZ(1)
    !     NDim= 1
    !     Poly_Order: complete Polynomial degree of the shape function
    !     iS=0: family of Lagrangian elements (valid for Poly_Order>1)
    !     node: number (ordinal) of the shape function)
    !     PsiCor(ndim)= Simplex intrinsec coordinates
    !     Deriv: order  of Partial DERIVATIVE (Function=:0)

    Deriv = sum(d)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Line Element

    X  = PsiCor(1)
    XX = X*X
    SF1D = 0.0D0

    ! Xi: Coordenada X del nodo correspondiente a la funcion de forma Ni

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Select Case (Poly_Order)
           ! Linear:1, Quadratic:2, Cubic:3, Etc.
        Case (1)
            ! Poly_Order Linear
            Select Case (Deriv)
                Case (0)
                    ! Function Values
                    Select Case (node)
                        Case(1)
                            SF1D = 0.50d0*(1.0d0-X)
                        Case(2)
                            SF1D = 0.50d0*(1.0d0+X)
                    End Select
                Case (1)
                    ! First Order Derivative
                    Select Case (node)
                        Case(1)
                            SF1D = -0.50d0
                        Case(2)
                            SF1D = 0.50d0
                    End Select
                Case (2)
                    ! Second Order Derivative
                    Select Case (node)
                        Case(1)
                            SF1D = 0.0d0
                        Case(2)
                            SF1D = 0.0d0
                    End Select
            End Select

        Case (2)
            ! Poly_Order Cuadratic
            Select Case (Deriv)
                Case (0)
                    ! Function Values
                    Select Case (node)
                        Case(1)
                            SF1D = -0.50d0*X*(1.0d0-X)
                        Case(2)
                            SF1D = 0.50d0*X*(1.0d0+X)
                        Case(3)
                            SF1D = 1.0d0-XX
                    End Select
                Case (1)
                    ! First Order Derivative
                    Select Case (node)
                        Case(1)
                            SF1D = -0.50d0*(1.0d0-2.0d0*X)
                        Case(2)
                            SF1D = 0.50d0*(1.0d0+2.0d0*X)
                        Case(3)
                            SF1D = -2.0d0*X
                    End Select
                Case (2)
                    ! Second Order Derivative
                    Select Case (node)
                        Case(1)
                            SF1D = 1.0d0
                        Case(2)
                            SF1D = 1.0d0
                        Case(3)
                            SF1D = -2.0d0
                    End Select
            End Select

    End Select

    Return
End Function SF1D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


! keywords
! Derivadas de funciones de forma lineares para Brics / Quads / Cubos
! Shape derivative functions linear 
! Joaquin Aranciaga
! Santiago 
! Souza Neto
! Fbar
! 2015
! Residual Stress
! ndim=2 y ndim=3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
SubRoutine localShapeFunctionHyperCube (Ndim, lsf, Gauss)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Implicit Real*8 (a-h, o-z)

    Real*8 lsf

          ! http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch11.d/AFEM.Ch11.pdf

    Dimension lsf(8), Gauss(Ndim)

    if (Ndim == 3) then
        psi_1 = Gauss(1)
        psi_2 = Gauss(2)
        psi_3 = Gauss(3)

        !     Funciones de forma del paralelepopedo lineales en los puntos
        !     de Gauss.
      
              !lsf(1) = (1.0d0 - psi_2) * (1.0d0 - psi_3) * psi_1 !lsf(1) = (1.0d0 - psi_1) * (1.0d0 - psi_2) * (1.0d0 - psi_3)
        !lsf(2) = (1.0d0 - psi_3) * psi_1 * psi_2 !lsf(2) = (1.0d0 - psi_2) * (1.0d0 - psi_3) * psi_1
        !lsf(3) = (1.0d0 - psi_1) * (1.0d0 - psi_3) * psi_2 !lsf(3) = (1 - psi_3) * psi_1 * psi_2
        !lsf(4) = (1.0d0 - psi_1) * (1.0d0 - psi_2) * (1.0d0 - psi_3) !lsf(4) = (1 - psi_1) * (1 - psi_3) * psi_2
        !lsf(5) = (1.0d0 - psi_2) * psi_1 * psi_3 !lsf(5) = (1.0d0 - psi_1) * (1.0d0 - psi_2) * psi_3
        !lsf(6) = psi_1 * psi_2 * psi_3 !lsf(6) = (1.0d0 - psi_2) * psi_1 * psi_3
        !lsf(7) = (1.0d0 - psi_1) * psi_2 * psi_3 !lsf(7) = psi_1 * psi_2 * psi_3
        !lsf(8) = (1.0d0 - psi_1) * (1.0d0 - psi_2) * psi_3 !lsf(8) = (1.0d0 - psi_1) * psi_2 * psi_3
	   
        lsf(1) = (1.0d0 - psi_1) * (1.0d0 - psi_2) * (1.0d0 - psi_3)
        lsf(2) = (1.0d0 - psi_2) * (1.0d0 - psi_3) * (psi_1 + 1.0d0)
        lsf(3) = (1.0d0 - psi_3) * (psi_1 + 1.0d0) * (psi_2 + 1.0d0)
        lsf(4) = (1.0d0 - psi_1) * (1.0d0 - psi_3) * (psi_2 + 1.0d0)
        lsf(5) = (1.0d0 - psi_1) * (1.0d0 - psi_2) * (psi_3 + 1.0d0)
        lsf(6) = (1.0d0 - psi_2) * (psi_1 + 1.0d0) * (psi_3 + 1.0d0)
        lsf(7) = (psi_1 + 1.0d0) * (psi_2 + 1.0d0) * (psi_3 + 1.0d0)
        lsf(8) = (1.0d0 - psi_1) * (psi_2 + 1.0d0) * (psi_3 + 1.0d0)
        lsf = lsf /8.0d0
    end if
	
    if (Ndim == 2) then
        psi_1 = Gauss(1)
        psi_2 = Gauss(2)
       
        lsf(1) = (1.0d0 - psi_1) * (1.0d0 - psi_2)
        lsf(2) = (1.0d0 - psi_2) * (psi_1 + 1.0d0)
        lsf(3) = (psi_1 + 1.0d0) * (psi_2 + 1.0d0)
        lsf(4) = (1.0d0 - psi_1) * (psi_2 + 1.0d0)
        lsf = lsf /4.0d0
    end if

    Return

End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
SubRoutine localShapeFunctionDerivativesHyperCube &
    (NodG, Ndim, dlsf, Gauss)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Implicit Real*8 (a-h, o-z)

    Dimension dlsf(Ndim, NodG), Gauss(Ndim)
	
    if (Ndim == 3) then
        psi_1 = Gauss(1)
        psi_2 = Gauss(2)
        psi_3 = Gauss(3)

        dlsf(1,1) = - (1.0d0 - psi_2) * (1.0d0 - psi_3)
        dlsf(1,2) =   (1.0d0 - psi_2) * (1.0d0 - psi_3)
        dlsf(1,3) =   (1.0d0 - psi_3) * (psi_2 + 1.0d0)
        dlsf(1,4) = - (1.0d0 - psi_3) * (psi_2 + 1.0d0)
        dlsf(1,5) = - (1.0d0 - psi_2) * (psi_3 + 1.0d0)
        dlsf(1,6) =   (1.0d0 - psi_2) * (psi_3 + 1.0d0)
        dlsf(1,7) =   (psi_2 + 1.0d0) * (psi_3 + 1.0d0)
        dlsf(1,8) = - (psi_2 + 1.0d0) * (psi_3 + 1.0d0)

        dlsf(2,1) = - (1.0d0 - psi_1) * (1.0d0 - psi_3)
        dlsf(2,2) = - (1.0d0 - psi_3) * (psi_1 + 1.0d0)
        dlsf(2,3) =   (1.0d0 - psi_3) * (psi_1 + 1.0d0)
        dlsf(2,4) =   (1.0d0 - psi_1) * (1.0d0 - psi_3)
        dlsf(2,5) = - (1.0d0 - psi_1) * (psi_3 + 1.0d0)
        dlsf(2,6) = - (psi_1 + 1.0d0) * (psi_3 + 1.0d0)
        dlsf(2,7) =   (psi_1 + 1.0d0) * (psi_3 + 1.0d0)
        dlsf(2,8) =   (1.0d0 - psi_1) * (psi_3 + 1.0d0)

        dlsf(3,1) = - (1.0d0 - psi_1) * (1.0d0 - psi_2)
        dlsf(3,2) = - (1.0d0 - psi_2) * (psi_1 + 1.0d0)
        dlsf(3,3) = - (psi_1 + 1.0d0) * (psi_2 + 1.0d0)
        dlsf(3,4) = - (1.0d0 - psi_1) * (psi_2 + 1.0d0)
        dlsf(3,5) =   (1.0d0 - psi_1) * (1.0d0 - psi_2)
        dlsf(3,6) =   (1.0d0 - psi_2) * (psi_1 + 1.0d0)
        dlsf(3,7) =   (psi_1 + 1.0d0) * (psi_2 + 1.0d0)
        dlsf(3,8) =   (1.0d0 - psi_1) * (psi_2 + 1.0d0)

        dlsf = dlsf / 8.0d0

    end if
	
    if (Ndim == 2) then
        psi_1 = Gauss(1)
        psi_2 = Gauss(2)

        dlsf(1,1) = - (1.0d0 - psi_2)
        dlsf(1,2) =   (1.0d0 - psi_2)
        dlsf(1,3) =   (psi_2 + 1.0d0)
        dlsf(1,4) = - (psi_2 + 1.0d0)
        dlsf(2,1) = - (1.0d0 - psi_1)
        dlsf(2,2) = - (psi_1 + 1.0d0)
        dlsf(2,3) =   (psi_1 + 1.0d0)
        dlsf(2,4) =   (1.0d0 - psi_1)

        dlsf = dlsf/4.0d0

    end if
    Return

End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
