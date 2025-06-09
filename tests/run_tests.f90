program run_tests
  use test_micmac
  use test_fitting
  use test_optimise
  use test_ho
  implicit none

  print *, "==============================="
  print *, "  FORTRAN TEST SUITE RUNNER"
  print *, "==============================="

  !call run_tests_micmac()
  call run_tests_fitting()
  call run_tests_optimise()
  call run_tests_ho()

  print *, "==============================="
  print *, "         TESTING DONE"
  print *, "==============================="

end program run_tests
