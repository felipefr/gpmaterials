cmake_minimum_required(VERSION 2.6)

enable_language (Fortran)
 
# definição do nome do project. Não precisa ser obrigatoriamente o nome do diretório
project(gpmaterials)

ADD_DEFINITIONS(-fPIC -g -fbacktrace)

SET(Directoires ../src/tools ../src/FiniteStrain  ../src/Multiscale  ../src/damage  ../src/fibres  
		../src/otherMaterials ../src/arcLength ../src/arcLength_simple)
SET(fileExtensions .f .for .f90 .F .F90 .FOR) 
SET(src)
foreach(d ${Directoires})
	foreach(e ${fileExtensions})
		file(GLOB dSources ${d}/*${e})
		LIST(APPEND src ${dSources})
	endforeach(e)
endforeach(d)

add_library(objlib OBJECT ${src})
add_library(gpmaterials SHARED $<TARGET_OBJECTS:objlib>)
add_library(gpmaterials_static STATIC $<TARGET_OBJECTS:objlib>)

#add_library(gpmaterials SHARED ${src})

#~ add_custom_target(
	#~ multiscalenewlib
    #~ DEPENDS ../src/Multiscale/multiscaleNewLib.f90
    #~ )

#~ add_dependencies(globalvariables multiscalenewlib)
#~ add_dependencies(multiscalenewlib globalvariables)
