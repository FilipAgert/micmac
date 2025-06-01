program run_tests
  use test_micmac
  use test_fitting
  use test_minimise
  implicit none

  print *, "==============================="
  print *, "  FORTRAN TEST SUITE RUNNER"
  print *, "==============================="

  call run_tests_micmac()
  call run_tests_fitting()
  call run_tests_minimise()

  print *, "==============================="
  print *, "         TESTING DONE"
  print *, "==============================="

end program run_tests
