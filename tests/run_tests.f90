program run_tests
  use test_micmac
  use test_fitting
  implicit none

  print *, "==============================="
  print *, "  FORTRAN TEST SUITE RUNNER"
  print *, "==============================="

  call run_tests_micmac()
  call run_tests_fitting()

  print *, "==============================="
  print *, "         TESTING DONE"
  print *, "==============================="

end program run_tests
