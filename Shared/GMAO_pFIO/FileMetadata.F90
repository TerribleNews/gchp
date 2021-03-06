#include "pFIO_ErrLog.h"
#include "unused_dummy.H"

module pFIO_FileMetadataMod
   use pFIO_KeywordEnforcerMod
   use pFIO_StringIntegerMapMod
   use pFIO_StringIntegerMapUtilMod
   use pFIO_ThrowMod
   use pFIO_ConstantsMod
   use pFIO_UtilitiesMod
   use pFIO_ErrorHandlingMod
   use pFIO_UnlimitedEntityMod
   use pFIO_AttributeMod
   use pFIO_VariableMod
   use pFIO_CoordinateVariableMod
   use pFIO_StringVariableMapMod
   use pFIO_StringVariableMapUtilMod
   use pFIO_StringAttributeMapMod
   use pFIO_StringVectorMod
   implicit none
   private

   public :: FileMetadata 

   type :: FileMetadata
      private
      type (StringIntegerMap) :: dimensions
      type (Variable) :: global
      type (StringVariableMap) :: variables
   contains

      procedure :: get_dimensions
      procedure :: add_dimension
      procedure :: get_dimension
      procedure :: modify_dimension

      procedure :: get_attributes
      generic :: add_attribute => add_attribute_0d, add_attribute_1d
      procedure :: add_attribute_0d
      procedure :: add_attribute_1d
      procedure :: get_attribute
      procedure :: has_attribute

      procedure :: get_variable
      procedure :: get_coordinate_variable
      procedure :: add_variable
      procedure :: get_variables
      procedure :: modify_variable
      procedure :: has_dimension

      generic :: operator(==) => equal
      generic :: operator(/=) => not_equal
      procedure :: equal
      procedure :: not_equal

      procedure :: serialize
      procedure :: deserialize
      procedure :: is_coordinate_variable

   end type FileMetadata


contains


   function get_dimensions(this) result(dimensions)
      type (StringIntegerMap), pointer :: dimensions
      class (FileMetadata), target, intent(in) :: this

      dimensions => this%dimensions

   end function get_dimensions


   subroutine add_dimension(this, dim_name, extent, unusable, rc)
      class (FileMetadata), target, intent(inout) :: this
      character(len=*), intent(in) :: dim_name
      integer, intent(in) :: extent

      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      _UNUSED_DUMMY(unusable)
      call this%dimensions%insert(dim_name, extent)
      _RETURN(_SUCCESS)      
   end subroutine add_dimension

   subroutine modify_dimension(this, dim_name, extent, unusable, rc)
      class (FileMetadata), target, intent(inout) :: this
      character(len=*), intent(in) :: dim_name
      integer, intent(in) :: extent

      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      _UNUSED_DUMMY(unusable)

      call this%dimensions%set(dim_name, extent)
      _RETURN(_SUCCESS)      
      
   end subroutine modify_dimension

   function has_dimension(this, dim_name, unusable, rc) result(isPresent)
      class (FileMetadata), target, intent(in) :: this
      character(len=*), intent(in) :: dim_name
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      type (StringIntegerMapIterator) :: iter
      logical :: isPresent

      _UNUSED_DUMMY(unusable)
      iter = this%dimensions%find(dim_name)
      isPresent = (iter /=this%dimensions%end())
      _RETURN(_SUCCESS)
      
   end function has_dimension

   integer function get_dimension(this, dim_name, unusable, rc) result(extent)
      class (FileMetadata), target, intent(in) :: this
      character(len=*), intent(in) :: dim_name
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      type (StringIntegerMapIterator) :: iter
      character(len=255) :: msg

      _UNUSED_DUMMY(unusable)
      iter = this%dimensions%find(dim_name)

      if (iter /= this%dimensions%end()) then
         extent = this%dimensions%at(dim_name)
         _RETURN(_SUCCESS)
      else
         extent = 0
         msg = 'FileMetadata::get_dimension() - no such dimension <'//dim_name//'>.'
         _ASSERT(.false., trim(msg))
      end if
      
   end function get_dimension


   subroutine add_attribute_0d(this, attr_name, attr_value, unusable, rc)
      class (FileMetadata), target, intent(inout) :: this
      character(len=*), intent(in) :: attr_name
      class (*), intent(in) :: attr_value
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      _UNUSED_DUMMY(unusable)
      call this%global%add_attribute(attr_name, attr_value)
      _RETURN(_SUCCESS)
   end subroutine add_attribute_0d

   subroutine add_attribute_1d(this, attr_name, values, unusable, rc)
      class (FileMetadata), target, intent(inout) :: this
      character(len=*), intent(in) :: attr_name
      class (*), intent(in) :: values(:)

      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      _UNUSED_DUMMY(unusable)
      call this%global%add_attribute(attr_name, values)
      _RETURN(_SUCCESS)
   end subroutine add_attribute_1d


   function get_attribute(this, attr_name, unusable, rc) result(ref)
      type (Attribute), pointer :: ref
      class (FileMetadata), target, intent(in) :: this
      character(len=*), intent(in) :: attr_name
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc
      character(len=255) :: msg

      _UNUSED_DUMMY(unusable)

      ref => this%global%get_attribute(attr_name)
      msg = 'FileMetadata::get_attribute() - no such attribute <'//attr_name//'>.'
      _ASSERT(associated(ref), trim(msg))
      _RETURN(_SUCCESS)
   end function get_attribute


   ! No RC is necessary - no failure mode.
   logical function has_attribute(this, attr_name)
      class (FileMetadata), target, intent(in) :: this
      character(len=*), intent(in) :: attr_name

      has_attribute = this%global%is_attribute_present(attr_name)
      
   end function has_attribute


   function get_attributes(this, rc ) result(attributes)
      type (StringAttributeMap), pointer :: attributes
      class (FileMetadata), target, intent(in) :: this
      integer, optional, intent(out) :: rc

      attributes => this%global%get_attributes()
      _RETURN(_SUCCESS)
   end function get_attributes


   function get_variable(this, var_name, unusable, rc) result(var)
      class (Variable), pointer :: var
      class (FileMetadata), target, intent(in) :: this
      character(len=*), intent(in) :: var_name
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      _UNUSED_DUMMY(unusable)
      var => this%variables%at(var_name)
      _RETURN(_SUCCESS)
   end function get_variable

   ! Returns null pointer unless var_name is a key corresponding to
   ! a CoordinateVariable value.
   ! rc returns _SUCCESS unless the var_name is not found at all.

   function get_coordinate_variable(this, var_name, unusable, rc) result(var)
      class (CoordinateVariable), pointer :: var
      class (FileMetadata), target, intent(in) :: this
      character(len=*), intent(in) :: var_name
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      class (Variable), pointer :: tmp
      character(len=255) :: msg      

      _UNUSED_DUMMY(unusable)

      tmp => this%variables%at(var_name)

      msg = 'can not find '//trim(var_name)
      _ASSERT(associated(tmp), trim(msg))

      select type (tmp)
      class is (CoordinateVariable)
         var => tmp
      class default
         var => null()
      end select

      _RETURN(_SUCCESS)
      
   end function get_coordinate_variable

   logical function is_coordinate_variable(this, var_name, unusable, rc) 
      class (FileMetadata),target, intent(in) :: this
      character(*), intent(in) :: var_name

      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      class (Variable), pointer :: tmp
      character(len=255) :: msg      

      _UNUSED_DUMMY(unusable)

      tmp => this%variables%at(var_name)

      msg = 'can not find the varaible '//trim(var_name)
      _ASSERT(associated(tmp), trim(msg))
      select type (tmp)
      class is (CoordinateVariable)
         is_coordinate_variable = .true.
      class default
         is_coordinate_variable = .false.
      end select
      
      _RETURN(_SUCCESS)
   end function is_coordinate_variable


   function get_variables(this, rc ) result(variables)
      type (StringVariableMap), pointer :: variables
      class (FileMetadata), target, intent(in) :: this
      integer, optional, intent(out) :: rc

      variables => this%variables
      _RETURN(_SUCCESS)
   end function get_variables


   subroutine add_variable(this, var_name, var, unusable, rc)
      class (FileMetadata), target, intent(inout) :: this
      character(len=*), intent(in) :: var_name
      class (Variable), intent(in) :: var
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      type (StringVector), pointer :: dims
      type (StringVectorIterator) :: iter
      integer, pointer :: dim_this
      character(len=:), pointer :: dim_name
      type (UnlimitedEntity), pointer :: const_value_ptr
      integer, allocatable :: shp(:), shp_const(:)
      integer :: empty(0)
      character(len=255) :: msg

      _UNUSED_DUMMY(unusable)

      ! ensure all of var's dimensions are defined
      shp = empty
      dims => var%get_dimensions()
      iter = dims%begin()
      do while (iter /= dims%end())

         dim_name => iter%get()
         dim_this => this%dimensions%at(dim_name)
         if ( .not. associated(dim_this) ) then
            write(*,'(4a)') 'FileMetadata::add_variable() - undefined dimension: ', &
                  Trim(dim_name), ' for var: ', Trim(var_name)
            _ASSERT( associated(dim_this), trim(msg))
         end if
         shp =[shp,dim_this]
         call iter%next()
      end do

      const_value_ptr => var%get_const_value()
      if ( .not. const_value_ptr%is_empty() ) then
         shp_const = const_value_ptr%get_shape()
         _ASSERT( all(shp == shp_const), "const_value shape does not match dims")
      endif

      call this%variables%insert(var_name, var)
      _RETURN(_SUCCESS)
      
   end subroutine add_variable

   subroutine modify_variable(this, var_name, var, unusable, rc)
      class (FileMetadata), target, intent(inout) :: this
      character(len=*), intent(in) :: var_name
      class (Variable), intent(in) :: var
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      type (StringVector), pointer :: dims
      type (StringVectorIterator) :: iter
      integer, pointer :: dim_this
      character(len=:), pointer :: dim_name
      character(len=255) :: msg

      _UNUSED_DUMMY(unusable)

      ! ensure all of var's dimensions are defined
      dims => var%get_dimensions()
      iter = dims%begin()
      do while (iter /= dims%end())
         dim_name => iter%get()
         dim_this => this%dimensions%at(dim_name)
         if ( .not. associated(dim_this) ) then
            msg = "FileMetadata:: modify_variable() - undefined dimension " &
                  // dim_name
            _ASSERT( associated(dim_this), trim(msg) )
         end if
         call iter%next()
      end do

      call this%variables%set(var_name, var)

      _RETURN(_SUCCESS)
      
   end subroutine modify_variable

   subroutine add_var_attribute_0d(this, var_name, attr_name, value, unusable, rc)
      class (FileMetadata), target, intent(inout) :: this
      character(len=*), intent(in) :: var_name
      character(len=*), intent(in) :: attr_name
      class (*), intent(in) :: value
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      class (Variable), pointer :: var

      _UNUSED_DUMMY(unusable)

      var => this%get_variable(var_name)
      call var%add_attribute(attr_name, value)

      _RETURN(_SUCCESS)
      
   end subroutine add_var_attribute_0d

   subroutine add_var_attribute_1d(this, var_name, attr_name, values, unusable, rc)
      class (FileMetadata), target, intent(inout) :: this
      character(len=*), intent(in) :: var_name
      character(len=*), intent(in) :: attr_name
      class (*), intent(in) :: values(:)

      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      class (Variable), pointer :: var

      _UNUSED_DUMMY(unusable)

      var => this%get_variable(var_name)
      call var%add_attribute(attr_name, values)
      _RETURN(_SUCCESS)
   end subroutine add_var_attribute_1d


   logical function equal(a, b)
      class (FileMetadata), target, intent(in) :: a
      class (FileMetadata), target, intent(in) :: b

      equal = same_dimensions(a,b)
      if (.not. equal) return

      equal = same_attributes(a,b)
      if (.not. equal) return

      equal = same_variables(a,b)

   contains

      logical function same_dimensions(a, b) result(equal)
         class (FileMetadata), target, intent(in) :: a
         class (FileMetadata), target, intent(in) :: b
         type (StringIntegerMapIterator) :: iter
         integer, pointer :: dim_a, dim_b
         character(len=:), pointer :: dim_name

         equal = a%dimensions%size() == b%dimensions%size()
         if (.not. equal) return
         
         iter = a%dimensions%begin()
         do while (iter /= a%dimensions%end())

            dim_name => iter%key()
            dim_b => b%dimensions%at(dim_name)

            equal = (associated(dim_b))
            if (.not. equal) return

            dim_a => iter%value()
            equal = (dim_a == dim_b)
            if (.not. equal) return

            call iter%next()
         end do
      end function same_dimensions

      logical function same_attributes(a, b) result(equal)
         class (FileMetadata), target, intent(in) :: a
         class (FileMetadata), target, intent(in) :: b

         type (StringAttributeMapIterator) :: iter
         type (Attribute), pointer :: attr_a, attr_b
         character(len=:), pointer :: attr_name

         equal = (a%global == b%global)

      end function same_attributes
      
      logical function same_variables(a, b) result(equal)
         class (FileMetadata), target, intent(in) :: a
         class (FileMetadata), target, intent(in) :: b

         type (StringVariableMapIterator) :: iter
         class (Variable), pointer :: var_a, var_b
         character(len=:), pointer :: var_name

         equal = a%variables%size() == b%variables%size()
         if (.not. equal) return

         iter = a%variables%begin()
         do while (iter /= a%variables%end())
            
            var_name => iter%key()
            var_b => b%variables%at(var_name)
            
            equal = (associated(var_b))
            if (.not. equal) return
            
            var_a => iter%value()
            equal = (var_a == var_b)
            if (.not. equal) return
            
            call iter%next()
         end do

      end function same_variables
      
   end function equal

   logical function not_equal(a, b)
      class (FileMetadata), target, intent(in) :: a
      class (FileMetadata), target, intent(in) :: b

      not_equal = .not. (a == b)
   end function not_equal

   subroutine serialize(this, buffer, rc)
      class (FileMetadata), intent(in) :: this
      integer,allocatable,intent(inout) :: buffer(:)
      integer, optional, intent(out) :: rc
      integer :: length 
      integer, allocatable :: tmp_buffer(:)

      if(allocated(buffer)) deallocate(buffer)
            
      call StringIntegerMap_serialize(this%dimensions, tmp_buffer)
      buffer = [tmp_buffer]
      call this%global%serialize(tmp_buffer)
      buffer = [buffer,tmp_buffer]
      call StringVariableMap_serialize(this%variables, tmp_buffer)
      buffer = [buffer,tmp_buffer]

      length = serialize_buffer_length(length) + size(buffer)
      buffer = [serialize_intrinsic(length),buffer]
      _RETURN(_SUCCESS)
   end subroutine

   subroutine deserialize(this, buffer, rc)
      class (FileMetadata), intent(inout) :: this
      integer, intent(in) :: buffer(:)
      integer, optional, intent(out) :: rc
      integer :: n,length
      integer :: status 
      n = 1
      call deserialize_intrinsic(buffer(n:),length)
      _ASSERT(length <= size(buffer), "length does not match")

      length = serialize_buffer_length(length)
      n = n+length
      this%dimensions = StringIntegerMap_deserialize(buffer(n:))
      call deserialize_intrinsic(buffer(n:),length)
      n = n + length
      call deserialize_intrinsic(buffer(n:),length)
      call this%global%deserialize(buffer(n:n+length-1), status)
      _VERIFY(status)
      n = n + length
      this%variables = StringVariableMap_deserialize(buffer(n:))
      _RETURN(_SUCCESS)
   end subroutine deserialize

end module pFIO_FileMetadataMod
