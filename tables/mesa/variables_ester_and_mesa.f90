   module var_mesa_ester
   
   
   implicit none
   
   
      integer, save :: species
      logical, save :: init=.false.
      integer, pointer, dimension(:), save :: net_iso, chem_id
      integer, save ::  kap_handle
      real(8),allocatable,  dimension(:), save :: xa
!    contains
   end module var_mesa_ester
