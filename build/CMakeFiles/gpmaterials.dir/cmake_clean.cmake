file(REMOVE_RECURSE
  "libgpmaterials.pdb"
  "libgpmaterials.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/gpmaterials.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
