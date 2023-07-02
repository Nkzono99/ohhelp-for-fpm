module ohhelp
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, ohhelp-for-fpm!"
  end subroutine say_hello
end module ohhelp
