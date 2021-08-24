C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C     A function that calculates the values of shape functions
C     at specified points for Simplex Elements (Triangles,TetraHedra)
C
C
C_______________________________________________________________________
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      Real*8 Function SFSimplex (PsiCor, d, Poly_Order, node, NDim)
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C-----------------------------------------------------------------------
      Implicit real*8 (a-h, o-z)
      INTEGER Poly_Order, d(NDim), Deriv
      Dimension Sfn(ndim+1), PsiCor(ndim)
C     NDim = dimension of the simplex
C     NDim = 1 Line
C     NDim = 2 Triangle
C     NDim = 3 Tetrahedra
C     Poly_Order: complete Polynomial degree of the shape function
C     node: number (ordinal) of the shape function)
C     PsiCor(ndim)= Simplex intrinsec coordinates
C     Deriv: order  of Partial DERIVATIVE (Function=:0)
C     d(1): order of derivation with respect to Psi1
C     d(2): order of derivation with respect to Psi2
C     d(3): order of derivation with respect to Psi3
C     Sfn(i): Linear Shape Function  corresponding to Node i
C     (Sfn(2) = PsiCor(1) coordinate)
C     sfN3: Linear Shape Function  corresponding to Node 3
C     (Psi2 coordinate)
C     sfN4: Linear Shape Function  corresponding to Node 4
C     (Psi3 coordinate)
C     sfN(1)= 1-Psi1-psi2-psi3
C
C     sfN1= 1.0d0 - sfN2 - sfN3 - sfN4
C
      SFSimplex = 0.0D0
      Deriv = sum(d)
      Select Case (NDim)
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Case (3)
C     TetraHedral Element
         SfN(1) = 1.0d0 - PsiCor(1) - PsiCor(2) - PsiCor(3)
         SfN(2) = PsiCor(1)
         SfN(3) = PsiCor(2)
         SfN(4) = PsiCor(3)
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         Select Case (Poly_Order)
c     linear:1, cuadratic:2, cubic:3, etc
         Case (1)
c     Poly_Order Linear
            Select Case (Deriv)
            Case (0)
C     Function Values
               Select Case (node)
               Case (:4)
                  SFSimplex = SfN(node)
               Case Default
                  SFSimplex = 256.d0*SfN(1)*SfN(2)*sfN(3)*sfN(4)
               End Select
            Case (1)
c     First order derivative
               Select Case (node)
               Case (1)
                  SFSimplex = -1.0d0
               Case (2:4)
                  SFSimplex = 1.0d0 * d(node-1)
               Case (5)
c     SFSimplex =256.d0*sfN1*sfN2*sfN3*sfN4
                  Select Case (d(1)+2*d(2))
                  Case (1)
                     SFSimplex =256.d0*sfN(3)*sfN(4)*( sfN(1) - sfN(2) )
C
                  Case (2)
                     SFSimplex =256.d0*sfN(2)*sfN(4)*( sfN(1) - sfN(3) )
C
                  Case Default
                     SFSimplex =256.d0*sfN(2)*sfN(3)*( sfN(1) - sfN(4) )
                  End Select
C
               End Select
            Case (2)
c     Second order derivative
               Select Case (node)
               Case (:4)
                  SFSimplex = 0.0d0
               Case (5)
c     SFSimplex =256.d0*sfN1*sfN2*sfN3*sfN4
                  Select Case (d(1)+3*d(2)+7*d(3))
                  Case (2)
                     SFSimplex=-512.0d0*sfN(3)*sfN(4)
                  Case (6)
                     SFSimplex=-512.0d0*sfN(2)*sfN(4)
                  Case (14)
                     SFSimplex=-512.0d0*sfN(3)*sfN(2)
                  Case (4)
                     SFSimplex=256.0d0*sfN(4)*(sfN(1)-sfN(2)-sfN(3))
                  Case (8)
                     SFSimplex=256.0d0*sfN(3)*(sfN(1)-sfN(2)-sfN(4))
                  Case (10)
                     SFSimplex=256.0d0*sfN(2)*(sfN(1)-sfN(3)-sfN(4))
                  End Select
               End Select
C
            Case (3)
C     Third order derivative
               Select Case (node)
               Case (:4)
                  SFSimplex = 0.0d0
               Case (5)
                  if (d(1).eq.3 .or. d(2).eq.3 .or. d(3).eq.3) then
C     d111-d222-d333
                     SFSimplex = 0.0d0
                  else if (d(1).eq.2) then
C     d112-d113-d121-d211-d131-d311
C     d112-d113
                     SFSimplex = - 512.0d0*sfN(4-d(3))
C
                  else if (d(2).eq.2) then
C     d221-d223
                     SFSimplex = - 512.0d0*sfN(3+d(1)-d(3))
                  else if (d(3).eq.2) then
C     d331-d332
                     SFSimplex = - 512.0d0*sfN(2+d(1))
                  else
C     if(d(1).eq.1 .and. d(2).eq.1. and. d(3).eq.1) then
                     SFSimplex =256.d0*( sfN(1)-sfN(2)-sfN(3)-sfN(4))
                  Endif
C
               End Select
            Case (4)
C     Fourth order derivative
               Select Case (node)
               Case (:4)
                  SFSimplex = 0.0d0
               Case (5)
                  SFSimplex = -256.0d0*d(1)*d(2)*d(3)
               End Select
            Case (5:)
C     fifth order derivative
               SFSimplex = 0.0d0
C     Deriv
            End Select
         Case (2)
c     Poly_Order Cuadratic
C
            Select Case (Deriv)
            Case (0)
               Select Case (node)
               Case (1:4)
                  SFSimplex = sfN(node)*(2.0d0*sfN(node)-1.0d0)
C
               Case (5)
                  SFSimplex = 4.0d0*sfN(1)*sfN(2)
               Case (6)
                  SFSimplex = 4.0d0*sfN(3)*sfN(2)
               Case (7)
                  SFSimplex = 4.0d0*sfN(1)*sfN(3)
               Case (8)
                  SFSimplex = 4.0d0*sfN(1)*sfN(4)
               Case (9)
                  SFSimplex = 4.0d0*sfN(2)*sfN(4)
               Case (10)
                  SFSimplex = 4.0d0*sfN(3)*sfN(4)
               Case (11)
                  SFSimplex = 256.d0*sfN(1)*sfN(2)*sfN(3)*sfN(4)
               End Select
            Case (1)
c     First order derivative
               Select Case (node)
               Case (1)
                  SFSimplex = 1.0d0-4.0d0*sfN(1)
               Case (2:4)
                  SFSimplex = (4.0d0*sfN(node)-1.0d0)*d(node-1)
C
               Case (5)
                  SFSimplex = 4.0d0*(sfN(1)*d(1)-sfN(2))
               Case (6)
                  SFSimplex = 4.0d0*(sfN(3)*d(1)+sfN(2)*d(2))
               Case (7)
                  SFSimplex = 4.0d0*(sfN(1)*d(2)-sfN(3))
               Case (8)
                  SFSimplex = 4.0d0*(sfN(1)*d(3)-sfN(4))
               Case (9)
                  SFSimplex = 4.0d0*(sfN(4)*d(1)+sfN(2)*d(3))
               Case (10)
                  SFSimplex = 4.0d0*(sfN(4)*d(2)+sfN(3)*d(3))
               Case (11)
c     SFSimplex =256.d0*sfN1*sfN2*sfN3*sfN4
                  Select Case (d(1)+2*d(2))
                  Case (1)
                     SFSimplex =256.d0*sfN(3)*sfN(4)*(sfN(1)-sfN(2))
C
                  Case (2)
                     SFSimplex =256.d0*sfN(2)*sfN(4)*(sfN(1)-sfN(3))
C
                  Case Default
                     SFSimplex =256.d0*sfN(2)*sfN(3)*(sfN(1)-sfN(4))
                  End Select
               End Select
            Case (2)
c     Second order derivative
               Select Case (node)
               Case (1)
                  SFSimplex = 4.0d0
               Case (2:4)
                  Select Case (d(node-1))
                  Case (2)
                     SFSimplex = 4.0d0
                  Case Default
                     SFSimplex = 0.0d0
                  End Select
               Case (5)
                  SFSimplex = -4.0d0*d(1)
               Case (6)
                  SFSimplex =  4.0d0*d(1)*d(2)
               Case (7)
                  SFSimplex = -4.0d0*d(2)
               Case (8)
                  SFSimplex = -4.0d0*d(3)
               Case (9)
                  SFSimplex =  4.0d0*d(1)*d(3)
               Case (10)
                  SFSimplex =  4.0d0*d(2)*d(3)
               Case (11)
c     SFSimplex =256.d0*sfN1*sfN2*sfN3*sfN4
                  Select Case (d(1)+3*d(2)+7*d(3))
                  Case (2)
                     SFSimplex=-512.0d0*sfN(3)*sfN(4)
                  Case (6)
                     SFSimplex=-512.0d0*sfN(2)*sfN(4)
                  Case (14)
                     SFSimplex=-512.0d0*sfN(3)*sfN(2)
                  Case (4)
                     SFSimplex=256.0d0*sfN(4)*(sfN(1)-sfN(2)-sfN(3))
                  Case (8)
                     SFSimplex=256.0d0*sfN(3)*(sfN(1)-sfN(2)-sfN(4))
                  Case (10)
                     SFSimplex=256.0d0*sfN(2)*(sfN(1)-sfN(3)-sfN(4))
                  End Select
c     Node
               End Select
            Case (3)
C     Third order derivative
               Select Case (node)
               Case (:10)
                  SFSimplex = 0.0d0
               Case (11)
                  if (d(1).eq.3 .or. d(2).eq.3 .or. d(3).eq.3) then
C     d111-d222-d333
                     SFSimplex = 0.0d0
                  else if (d(1).eq.2) then
C     d112-d113-d121-d211-d131-d311
C     d112-d113
                     SFSimplex = - 512.0d0*sfN(4-d(3))
C
                  else if (d(2).eq.2) then
C     d221-d223
                     SFSimplex = - 512.0d0*sfN(3+d(1)-d(3))
                  else if (d(3).eq.2) then
C     d331-d332
                     SFSimplex = - 512.0d0*sfN(2+d(1))
                  else
C     if(d(1).eq.1 .and. d(2).eq.1. and. d(3).eq.1) then
                     SFSimplex =256.d0*( sfN(1)-sfN(2)-sfN(3)-sfN(4))
                  Endif
c     Node
               End Select
C
            Case (4)
C     Fourth order derivative
               Select Case (node)
               Case (:10)
                  SFSimplex = 0.0d0
               Case (11)
                  icoef = d(1)*d(2)*d(3)
                  SFSimplex = -256.0d0*icoef
               End Select
            Case (5:)
C     Fifth order derivative & so
               SFSimplex = 0.0d0
c     Deriv
            End Select
c     Poly_Order
         End Select
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Case (2)
C     Triangle Element
         SfN(1) = 1.0d0 - PsiCor(1) - PsiCor(2)
         SfN(2) = PsiCor(1)
         SfN(3) = PsiCor(2)
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         Select Case (Poly_Order)
c     linear:1, cuadratic:2, cubic:3, etc
         Case (1)
c     Poly_Order Linear
            Select Case (Deriv)
            Case (0)
C     Function Values
               Select Case (node)
               Case (:4)
                  SFSimplex = SfN(node)
               Case Default
                  SFSimplex = 27.d0*SfN(1)*SfN(2)*sfN(3)
               End Select
            Case (1)
c     First order derivative
               Select Case (node)
               Case (1)
                  SFSimplex = -1.0d0
C
               Case (2:3)
                  SFSimplex = 1.0d0 * d(node-1)
               Case (4)
c     SFSimplex =27.d0*sfN1*sfN2*sfN3
                  Select Case (d(1))
                  Case (1)
                     SFSimplex =27.d0*sfN(3)*( sfN(1) - sfN(2) )
                  Case Default
                     SFSimplex =27.d0*sfN(2)*( sfN(1) - sfN(3) )
                  End Select
C
               End Select
            Case (2)
c     Second order derivative
               Select Case (node)
               Case (:3)
                  SFSimplex = 0.0d0
               Case (4)
c     SFSimplex =27.d0*sfN1*sfN2*sfN3
                  Select Case (d(1)+3*d(2))
                  Case (2)
                     SFSimplex=-54.0d0*sfN(3)
                  Case (6)
                     SFSimplex=-54.0d0*sfN(2)
                  Case (4)
                     SFSimplex=27.0d0*(sfN(1)-sfN(2)-sfN(3))
                  End Select
               End Select
C
            Case (3)
C     Third order derivative
               Select Case (node)
               Case (:3)
                  SFSimplex = 0.0d0
               Case (4)
                  if (d(1).eq.3 .or. d(2).eq.3) then
C     d111-d222
                     SFSimplex = 0.0d0
                  else
C     d112-d221
                     SFSimplex = - 54.0d0
                  Endif
C
               End Select
            Case (4:)
C     Fourth order derivative and so forth
C
               SFSimplex = 0.0d0
C     Deriv
            End Select
         Case (2)
c     Poly_Order Cuadratic
C
            Select Case (Deriv)
            Case (0)
               Select Case (node)
               Case (1:3)
                  SFSimplex = sfN(node)*(2.0d0*sfN(node)-1.0d0)
C
               Case (4)
                  SFSimplex = 4.0d0*sfN(1)*sfN(2)
               Case (5)
                  SFSimplex = 4.0d0*sfN(3)*sfN(2)
               Case (6)
                  SFSimplex = 4.0d0*sfN(1)*sfN(3)
               Case (7)
                  SFSimplex = 27.d0*sfN(1)*sfN(2)*sfN(3)
               End Select
            Case (1)
c     First order derivative
               Select Case (node)
               Case (1)
                  SFSimplex = 1.0d0-4.0d0*sfN(1)
               Case (2:3)
                  SFSimplex = (4.0d0*sfN(node)-1.0d0)*d(node-1)
C
               Case (4)
                  SFSimplex = 4.0d0*(sfN(1)*d(1)-sfN(2))
               Case (5)
                  SFSimplex = 4.0d0*(sfN(3)*d(1)+sfN(2)*d(2))
               Case (6)
                  SFSimplex = 4.0d0*(sfN(1)*d(2)-sfN(3))
C
               Case (7)
c     SFSimplex =27.d0*sfN1*sfN2*sfN3
                  Select Case (d(1))
                  Case (1)
                     SFSimplex =27.d0*sfN(3)*( sfN(1) - sfN(2) )
                  Case Default
                     SFSimplex =27.d0*sfN(2)*( sfN(1) - sfN(3) )
                  End Select
               End Select
            Case (2)
c     Second order derivative
               Select Case (node)
               Case (1)
                  SFSimplex = 4.0d0
               Case (2:3)
                  Select Case (d(node-1))
                  Case (2)
                     SFSimplex = 4.0d0
                  Case Default
                     SFSimplex = 0.0d0
                  End Select
               Case (4)
                  SFSimplex = -4.0d0*d(1)
               Case (5)
                  SFSimplex =  4.0d0*d(1)*d(2)
               Case (6)
                  SFSimplex = -4.0d0*d(2)
               Case (7)
c     SFSimplex =27.d0*sfN1*sfN2*sfN3
                  Select Case (d(1)+3*d(2))
                  Case (2)
                     SFSimplex=-54.0d0*sfN(3)
                  Case (6)
                     SFSimplex=-54.0d0*sfN(2)
                  Case (4)
                     SFSimplex=27.0d0*(sfN(1)-sfN(2)-sfN(3))
                  End Select
c     Node
               End Select
            Case (3)
C     Third order derivative
               Select Case (node)
               Case (:6)
                  SFSimplex = 0.0d0
               Case Default
                  if (d(1).eq.3 .or. d(2).eq.3) then
C     d111-d222
                     SFSimplex = 0.0d0
                  else
C     d112-d221
                     SFSimplex = - 54.0d0
                  Endif
               End Select
            Case (4:)
C     Fourth order derivative and so forth
               SFSimplex = 0.0d0
c     Deriv
            End Select
c     Poly_Order
         End Select
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Case (1)
C     Line Element
         SfN(1) = 1.0d0 - PsiCor(1)
         SfN(2) = PsiCor(1)
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         Select Case (Poly_Order)
c     linear:1, cuadratic:2, cubic:3, etc
         Case (1)
c     Poly_Order Linear
            Select Case (Deriv)
            Case (0)
C     Function Values
               Select Case (node)
               Case (:2)
                  SFSimplex = SfN(node)
               Case Default
                  SFSimplex = 4.d0*SfN(1)*SfN(2)
               End Select
            Case (1)
c     First order derivative
               Select Case (node)
               Case (1)
                  SFSimplex = -1.0d0
               Case (2)
                  SFSimplex =  1.0d0
               Case default
                  SFSimplex =  4.0d0*( SfN(1)-SfN(2) )
               End Select
            Case (2)
c     Second order derivative and so forth
               Select Case (node)
               Case (3)
                  SFSimplex = -8.0d0
               Case default
                  SFSimplex = 0.0d0
               End Select
            Case (3:)
               SFSimplex = 0.0d0
C     Deriv
            End Select
         Case (2)
c     Poly_Order Cuadratic
C
            Select Case (Deriv)
            Case (0)
               Select Case (node)
               Case (1:2)
                  SFSimplex = sfN(node)*(2.0d0*sfN(node)-1.0d0)
C
               Case (3)
                  SFSimplex = 4.0d0*sfN(1)*sfN(2)
               End Select
            Case (1)
c     First order derivative
               Select Case (node)
               Case (1)
                  SFSimplex = 1.0d0 - 4.0d0*sfN(1)
               Case (2)
                  SFSimplex = (4.0d0*sfN(2)-1.0d0)
C
               Case (3)
                  SFSimplex = 4.0d0*(sfN(1)-sfN(2))
               End Select
            Case (2)
c     Second order derivative
               Select Case (node)
               Case (3)
                  SFSimplex = -8.0d0
               Case default
                  SFSimplex = 4.0d0
c     Node
               End Select
            Case (3:)
C     Third order derivative and so forth
               SFSimplex = 0.0d0
c     Deriv
            End Select
c     Poly_Order
         End Select
C
C
C
C     NDim
      End Select
C
C
C
      RETURN
      End
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
