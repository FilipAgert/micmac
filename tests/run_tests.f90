program run_tests
  use test_micmac
  implicit none

  print *, "==============================="
  print *, "  FORTRAN TEST SUITE RUNNER"
  print *, "==============================="

  call run_tests_micmac()

  print *, "==============================="
  print *, "         TESTING DONE"
  print *, "==============================="

end program run_tests
