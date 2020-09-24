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
 * timer.h
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <petsc.h>


class Timer {
  private:

    /*
     * time variables
     */
    std::chrono::high_resolution_clock::time_point t1, t2;

  public:
    /*
     * rank information
     */
    PetscInt rank;
    PetscInt num_proc;

    /*
     * Constructor
     */
    Timer ();

    /*
     * Destructor
     */
    ~Timer ();

    /*
     * Setup
     */
    PetscErrorCode setup();

    /*
     * timer starter
     */
    void timer_start ();

    /*
     * timer stopper
     */
    void timer_stop ();

    /*
     * print time
     */
    void timer_print ();
    void timer_print (const std::string &name);

    /*
     * Print the current Year-Month-Day-Hour-Minute-Second
     */
    void timestampYMDHMS ();
};

#endif /* OPT_TIMER_H_ */
