subroutine executerElement(id_Elem_Family, AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, & 
							Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
	implicit none
	
	integer :: id_Elem_Family,MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)
	

	
	Select Case (Id_Elem_Family)
		Case (100)
			Call  trivialElement(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (101)
			Call  PolyElement(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (102)
			Call  NullElement(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (110) 
			Call DirichletNode(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (111) 
			Call DirichletNodeLag(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (112) 
			Call zeroVariable(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (120) 
			Call NodalForce(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (121) 
			Call NodalForceLine(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (130)
			Call setMaterialParam(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (131)
			Call setDamageParam(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (132)
			Call  IncrementVariables(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (133)
			Call insertDeformation(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (134)
			Call insertDeformationGeneral(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (140)
			Call posProcElem(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
!~ 		Case (141)
!~ 			Call posProcElemOld(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
!~ 								Sol0,Sol1,CommonPar,&
!~ 								Param,JParam,DelT,DTm,Time)
		Case (142)
			Call  Cell2Point(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (143)
			Call posProcElemNew(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (150)
			Call Viscosity(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (200)
			Call FSgen(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (201)
			Call FSgenMixedStab(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (203)
			Call FSgen_newDamage(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (210)
		!> 210 - NeumannFiniteStrain.f90 (old 533)
			Call NeumannFiniteStrain(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (211)
			Call NeumannFiniteStrainSpatial(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (212)
			Call NeumannRefTraction(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
!~ 		Case (300)
!~ 			Call damageEvolution(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
!~ 								Sol0,Sol1,CommonPar,&
!~ 								Param,JParam,DelT,DTm,Time) # deprecated
		Case (400)
			Call StressHom(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)		
		Case (410)
			Call TotalDisp(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (4410)
			Call TotalDisp3D(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (420)
			Call minRestrictionBC2DExact(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (421)
			Call minRestrictionRVE(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (422)
			Call enforcePeriodic2D(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (423)
			Call MRmeanVolRVE(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (424)
			Call enforcePeriodic2D_inc(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (425)
			Call enforcePeriodic2D_inc_copy(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (430)
			Call canonicalProblem(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)  
!~ 		Case (431)
!~ 			Call computeTangent(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
!~ 								Sol0,Sol1,CommonPar,&
!~ 								Param,JParam,DelT,DTm,Time)
		Case (432)
			Call TangentHom(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (433)
			Call computeMinDetQ(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (434)
			Call computeMinDetQNew(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (440)
			Call DecisionBifurcation(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (441)
			Call LocDomainBC(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (442)
			Call MarkLocPointS(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (443)
			Call LocPointsBC(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (500)
			Call nonlinearFibres(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (501)
			Call nonlinearFibresDamage(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (502)
			Call nonlinearFibresDamage_localised(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (503)
			Call NonLinFibresGen(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (504)
			Call NonlinearFibresQuang(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (505)
			Call nonlinearFibresDamage_viscous(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (5505)
			Call nonlinearFibresDamage_viscous3D(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (506)
			Call torsionalSpring(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (507) ! Newest one
			Call networkConstraintGen_pureInc(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (5507) ! Newest one
			Call networkConstraintGen_pureInc3D(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (508) ! Newest one
			Call networkConstraintLinear(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (5508) ! Newest one
			Call networkConstraintLinear3D(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (509) ! Newest one
			Call networkConstraintGen(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (5509) ! Newest one
			Call networkConstraintGen3D(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (510)
			Call networkConstraint(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (511)
			Call networkConstraint_noVolume(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (512)
			Call affineBoundary(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (513)
			Call networkConstraint_delta(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (514)
			Call networkConstraint_RVEnormal_new_gen(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (515)
			Call networkConstraintTheta(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (516)
			Call affineBoundaryTheta(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (517)
			Call networkConstraint_RVEnormal(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (518)
			Call networkConstraint_RVEnormal_noIncrement(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (519)
			Call networkConstraint_RVEnormal_new(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (520)
			Call computeAnisotropyTensor(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (521)
			Call computeAnisotropyTensorInv(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)	
		Case (530)
			Call damageEvolutionFibres(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)	
		Case (550)
			Call PosProcFibres(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)		
		Case (5550)
			Call PosProcFibres3D(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (560)
			Call CanonicalProblemFibres(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)		
		Case (561)
			Call TangentHomFibres(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)		
		Case (600)
			Call arcLengthInconsistent(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (601)
			Call zeroVariable_arcLength(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (602)
			Call computeCoeficientsCilindrical_arcLength(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)	
		Case (603)
			Call computeIncrementEstimation(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)		
		Case (604)
			Call chooseIncrementLambda(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)	
		Case (605)
			Call incrementDisplacementsArcLength(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)	
		Case (606)
			Call enforcePeriodic2D_arcLength(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (607)
			Call incrementBoundaryMultipliersArcLength(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)		
		Case (608)
			Call IncrementVariables_Conditional(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (609)
			Call IncrementVariablesArcLength(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)	
		Case (610)
			Call arcLengthConsistent(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)		
		Case (611)
			Call incrementLambda(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)		
		Case (612)
			Call incrementLambdaArcLength(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (613)
			Call computeCoeficientsSpherical_arcLength(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)					
		Case (614)
			Call posProcElem_arcLength(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)		
		Case (615)
			Call zeroVariable_arcLength2(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)			
		Case (616)
			Call zeroVariable_arcLength_generic(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)			
		Case (617)
			Call computeCoeficientsArcLength_generic(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)			
		Case (700)
			Call arcLength_simple(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (701)
			Call zeroVariable_arcLength_simple(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)		
		Case (702)
			Call computeCoeficientsCilindrical_arcLength_simple(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)	
		Case (703)
			Call computeIncrementEstimation_simple(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (704)
			Call chooseIncrementLambda_simple(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)					
		Case (705)
			Call incrementDisplacementsArcLength_simple(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)	
		Case (706)
			Call computeCoeficientsSpherical_arcLength_simple(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)	
		Case Default
          Write(*,*) '!!!!!!! WARNING, ELEMENT NOT IDENTIFIED', Id_Elem_Family

       End Select 

end subroutine

subroutine executerSymbolic(id_Elem_Family, Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
	implicit none
	
	Integer :: id_Elem_Family, MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
	
	Select Case (Id_Elem_Family)
		Case (100)
		!> 100 - TrivialElement.f90
			Call  trivialElementS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (101)
		!> 101 - PolyElement.f90
			Call  PolyElementS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (102)
		!> 102 - NullElement.f90
			Call  NullElementS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (110) 
		!> 110 - DirichletNode.f90 (old 518)
			Call DirichletNodeS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (111) 
		!> 111 - DirichletNodeLag.f90 (old 519)
			Call DirichletNodeLagS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (112) 
		!> 112 - zeroVariable.f90 (old 519)
			Call zeroVariableS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (120) 
		!> 120 - NodalForce.f90 (old 1008)
			Call NodalForceS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (121) 
		!> 121 - NodalForceLine.f90 !!! (old 539) For Cook Problem load
			Call NodalForceLineS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (130)
		!> 130 - setMaterialParam.f90 (old 901)
			Call setMaterialParamS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) 
		Case (131)
		!> 131 - setDamageParam.f90 (old 902)
			Call setDamageParamS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (132)
		!> 132 - IncrementVariables.f90 (old 104)
			Call  IncrementVariablesS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (133)
		!> 133 - insertDeformation.f90  (old 900)
			Call insertDeformationS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) 
		Case (134)
		!> 134 - insertDeformationGeneral.f90 (old 903)
			Call insertDeformationGeneralS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) 
		Case (140)
		!> 140 - posProcElem.f90 (old 797)
			Call posProcElemS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!~ 		Case (141)
!~ 		!> 141 - posProcElemOld.f90 (old 806)
!~ 			Call PosProcElemOldS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (142)
		!> 142 - Cell2Point.f90 (old 103)
			Call  Cell2PointS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)	
		Case (143)
		!> 143 - PosProcElemNew.f90 (old 103)
			Call  PosProcElemNewS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (150)
		!> 150 - Viscosity.f90 (0ld 514) 
			Call ViscosityS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (200)
		!> 200 - FSgen.f90 (old 800)
			Call FSgenS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) 
		Case (201)
		!> 201 - FSgenMixedStab.f90 (old 798)
			Call FSgenMixedStabS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (203)
		!> 203 - FSgen_newDamage.f90 (old 798)
			Call FSgen_newDamageS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (210)
		!> 210 - NeumannFiniteStrain.f90 (old 533)
			Call NeumannFiniteStrainS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (211)
		!> 211 - NeumannFiniteStrainSpatial.f90 (old 534)
			Call NeumannFiniteStrainSpatialS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (212)
		!> 212 - NeumannRefTraction.f90 (old 535)
			Call NeumannRefTractionS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!~ 		Case (300)
!~ 		!> 300 - damageEvolution.f90 (old 300)
!~ 			Call damageEvolutionS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) # deprecated
		Case (400)
		!> 400 - StressHom.f90 (old 1024)
			Call StressHomS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)		
		Case (410)
		!> 410 - TotalDisp.f90 (old 802)
			Call TotalDispS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) 
		Case (4410)
		!> 4410 - TotalDisp3D.f90 (old 802)
			Call TotalDisp3DS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) 
		Case (420)
		!> 420 - minRestrictionBC2DExact.f90 (old 512)
			Call minRestrictionBC2DExactS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (421)
		!> 421 - minRestrictionRVE.f90 (old 807)
			Call minRestrictionRVES(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (422)
		!> 422 - enforcePeriodic2D.f90 (old 502)
			Call enforcePeriodic2DS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (423)
		!> 423 - MRmeanVolRVE.f90 (old 807)
			Call MRmeanVolRVES(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (424)
		!> 424 - enforcePeriodic2D_inc.f90 (old ..)
			Call enforcePeriodic2D_incS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (425)
		!> 425 - enforcePeriodic2D_inc_copy.f90 (old ..)
			Call enforcePeriodic2D_inc_copyS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (430)
		!> 430 - canonicalProblem.f90 (old 804)
			Call canonicalProblemS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)  
!~ 		Case (431)
!~ 		!> 431 - computeTangent.f90 (old 541)
!~ 			Call computeTangentS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (432)
		!> 432 - TangentHom.f90 (old 542)
			Call TangentHomS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (433)
		!> 433 - computeMinDetQ.f90 (old 700)
			Call computeMinDetQS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (434)
		!> 434 - computeMinDetQNew.f90 (old 700)
			Call computeMinDetQNewS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (440)
		!> 440 - DecisionBifurcation.f90 (old 600)
			Call DecisionBifurcationS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (441)
		!> 441 - LocDomainBC.f90 (old 601)
			Call LocDomainBCS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (442)
		!> 442 - MarkLocPoint.f90 (old 701)
			Call MarkLocPointS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (443)
		!> 443 - LocPointsBC.f90 (old 702)
			Call LocPointsBCS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (500)
		!> 500 - nonlinearFibres.f90 (old 1033)
			Call nonlinearFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (501)
		!> 501 - nonlinearFibresDamage.f90 (old 1033)
			Call nonlinearFibresDamageS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (502)
		!> 502 - nonlinearFibresDamage_localised.f90 (old 1033)
			Call nonlinearFibresDamage_localisedS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (503)
		!> 503 - NonLinFibresGen.f90 (old 1033)
			Call NonLinFibresGenS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (504)
		!> 504 - NonLinearFibresQuang.f90 (old 1033)
			Call NonlinearFibresQuangS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (505)
		!> 505 - nonlinearFibresDamage_viscous.f90
			Call nonlinearFibresDamage_viscousS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (5505)
		!> 5505 - nonlinearFibresDamage_viscous3D.f90
			Call nonlinearFibresDamage_viscous3DS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (506)
		!> 506 - torsionalSpring.f90
			Call torsionalSpringS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (507) ! Newest
		!> 507 - networkConstraintGen_pureInc.f90 
			Call networkConstraintGen_pureIncS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (5507) ! Newest
		!> 5507 - networkConstraintGen_pureInc.f90 
			Call networkConstraintGen_pureInc3DS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (508) ! Newest
		!> 508 - networkConstraintLinear.f90 
			Call networkConstraintLinearS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (5508) ! Newest
		!> 5508 - networkConstraintLinear3D.f90 
			Call networkConstraintLinear3DS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (509) ! Newest
		!> 509 - networkConstraintGen.f90 
			Call networkConstraintGenS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (5509) ! Newest
		!> 5509 - networkConstraintGen3D.f90 
			Call networkConstraintGen3DS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (510)
		!> 510 - networkConstraint.f90 (old 1034)
			Call networkConstraintS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (511)
		!> 511 - networkConstraint_noVolume.f90 (old 1034)
			Call networkConstraint_noVolumeS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (512)
		!> 512 - affineBoundary.f90
			Call affineBoundaryS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (513)
		!> 513 - networkConstraint_delta.f90 (old 1034)
			Call networkConstraint_deltaS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (514)
		!> 514 - networkConstraint_RVEnormal_new_gen.f90 
			Call networkConstraint_RVEnormal_new_genS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (515)
		!> 515 - networkConstraintTheta.f90 (old 1034)
			Call networkConstraintThetaS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (516)
		!> 516 - affineBoundaryTheta.f90
			Call affineBoundaryThetaS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (517)
		!> 517 - networkConstraint_RVEnormal.f90
			Call networkConstraint_RVEnormalS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (518)
		!> 518 - networkConstraint_RVEnormal_.f90
			Call networkConstraint_RVEnormal_noIncrementS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (519)
		!> 519 - networkConstraint_RVEnormal_new.f90
			Call networkConstraint_RVEnormal_newS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (520)
		!> 520 - computeAnisotropyTensor.f90 (old 1015)
			Call computeAnisotropyTensorS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (521)
		!> 521 - computeAnisotropyTensorInv.f90 (old 1016)
			Call computeAnisotropyTensorInvS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)	
		Case (530)
		!> 530 - damageEvolutionFibres.f90 (old 1013)
			Call damageEvolutionFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)	
		Case (550)
		!> 550 - posProcFibres.f90 (old 1022)
			Call PosProcFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (5550)
		!> 5550 - posProcFibres3D.f90 (old 1022)
			Call PosProcFibres3DS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (560)
		!> 560 - CanonicalProblemFibres.f90 (old 1022)
			Call CanonicalProblemFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (561)
		!> 561 - TangentHomFibres.f90 (old 1022)
			Call TangentHomFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (600)
		!> 600 - arcLengthInconsistent.f90 (old 798)
			Call arcLengthInconsistentS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (601)
		!> 601 - zeroVariable_arcLength.f90 (old 798)
			Call zeroVariable_arcLengthS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (602)
		!> 602 - computeCoeficientsCilindrical_arcLength.f90 (old 1022)
			Call computeCoeficientsCilindrical_arcLengthS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (603)
		!> 603 - incrementIncrementEstimation.f90 (old ..)
			Call computeIncrementEstimationS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (604)
		!> 604 - chooseIncrementLambda.f90 (old ..)
			Call chooseIncrementLambdaS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (605)
		!> 605 - incrementDisplacementsArcLength.f90 (old 1022)
			Call incrementDisplacementsArcLengthS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (606)
		!> 606 - enforcePeriodic2D_arcLength.f90 (old ..)
			Call enforcePeriodic2D_arcLengthS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (607)
		!> 607 - incrementBoundaryMultipliersArcLength.f90 (old ..)
			Call incrementBoundaryMultipliersArcLengthS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (608)
		!> 608 - IncrementVariables_Conditional.f90 (old 104)
			Call IncrementVariables_ConditionalS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)			
		Case (609)
		!> 609 - IncrementVariablesArcLength.f90 (old ..)
			Call IncrementVariablesArcLengthS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (610)
		!> 610 - arcLengthConsistent.f90 (old 1022)
			Call arcLengthConsistentS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (611)
		!> 611 - incrementLambda.f90 (old 1022)
			Call incrementLambdaS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (612)
		!> 612 - incrementLambdaArcLength.f90 (old 1022)
			Call incrementLambdaArcLengthS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (613)
		!> 613 - computeCoeficientsSpherical_arcLength.f90 (old 1022)
			Call computeCoeficientsSpherical_arcLengthS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (614)
		!> 614 - PosProcElem_arcLength.f90 (old 103)
			Call  PosProcElem_arcLengthS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (615)
		!> 615 - zeroVariable_arcLength2.f90 (old 798)
			Call zeroVariable_arcLength2S(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (616)
		!> 616 - zeroVariable_arcLength_generic.f90 (old 798)
			Call zeroVariable_arcLength_genericS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (617)
		!> 617 - computeCoeficientsArcLength_generic.f90 (old 798)
			Call computeCoeficientsArcLength_genericS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (700)
		!> 700 - arcLength_simple.f90 (old 798)
			Call arcLength_simpleS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (701)
		!> 701 - zeroVariable_arcLength_simple.f90 (old 1022)
			Call zeroVariable_arcLength_simpleS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (702)
		!> 702 - computeCoeficientsCilindrical_arcLength_simple.f90 (old 1022)
			Call computeCoeficientsCilindrical_arcLength_simpleS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (703)
		!> 703 - computeIncrementEstimation_simple.f90 (old ..)
			Call computeIncrementEstimation_simpleS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (704)
		!> 704 - chooseIncrementLambda_simple.f90 (old ..)
			Call chooseIncrementLambda_simpleS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (705)
		!> 705 - incrementDisplacementsArcLength_simple.f90 (old 1022)
			Call incrementDisplacementsArcLength_simpleS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (706)
		!> 706 - computeCoeficientsSpherical_arcLength_simple.f90 (old 1022)
			Call computeCoeficientsSpherical_arcLength_simpleS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
    Case Default
           Write(*,*) '!!!!!!! WARNING, ELEMENT NOT IDENTIFIED'	, Id_Elem_Family
           Write(*,*) 'while trying to determine the Coupling Matrix'	


	End Select

end subroutine  



