//-------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the TopOptVox authors
//
// This file is part of the TopOptVox.
//
// The TopOptVox is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of TopOptVox.
//
// Author: Zhidong Brian Zhang
// Created on: April 12, 2020
//
// ---------------------------------------------------------------------

/*
 * timer.cc
 */

#include "timer.h"

Timer::Timer () {
  setup();
  timestampYMDHMS ();
}

Timer::~Timer(){}

PetscErrorCode Timer::setup(){
  PetscErrorCode ierr = 0;
  // Get the number of processes
  ierr = MPI_Comm_size (PETSC_COMM_WORLD, &num_proc);
  // Get the individual process ID
  ierr = MPI_Comm_rank (PETSC_COMM_WORLD, &rank);
  CHKERRQ (ierr);
  return ierr;
}

void Timer::timer_start () {
  if (rank == 0) t1 = std::chrono::high_resolution_clock::now ();
}

void Timer::timer_stop () {
  if (rank == 0) t2 = std::chrono::high_resolution_clock::now ();
}

void Timer::timer_print () {
  if (rank == 0)
    {
      double dur = std::chrono::duration<double> (t2 - t1).count ();
      std::cout << "\nXXX consumed time: " << std::setprecision (3) << dur << "s.\n" << std::endl;
    }
}

void Timer::timer_print (const std::string &oname) {
  if (rank == 0)
    {
      double dur = std::chrono::duration<double> (t2 - t1).count ();
      std::cout << '\n' << oname << " consumed time: " << std::setprecision (3) << dur << "s.\n" << std::endl;
    }
}

void Timer::timestampYMDHMS () {
  if (rank == 0)
    {
      const int TIME_SIZE = 80;
      static char time_buffer[TIME_SIZE];
      std::tm *tm_ptr;
      time_t now;
      time (&now);
      tm_ptr = localtime (&now);
      std::strftime (time_buffer, TIME_SIZE, "%b %d %Y %T", tm_ptr);
      std::cout << '\n';
      std::cout << '\n';
      std::cout << time_buffer << std::endl;
      std::cout << '\n';
      std::cout << '\n';
    }
}
