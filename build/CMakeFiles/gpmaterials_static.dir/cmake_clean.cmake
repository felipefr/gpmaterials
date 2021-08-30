file(REMOVE_RECURSE
  "libgpmaterials_static.a"
  "libgpmaterials_static.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/gpmaterials_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
