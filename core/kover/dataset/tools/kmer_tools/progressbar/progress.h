/*
*	Kover: Learn interpretable computational phenotyping models from k-merized genomic data
*	Copyright (C) 2015  Alexandre Drouin & Gael Letarte St-Pierre
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PROGRESS_H
#define PROGRESS_H

#include <iostream>
#include "callable.h"
    
class ProgressBar: public Callable
{
public:
    ProgressBar(unsigned int goal, bool visible, std::string toolname);
    virtual  void operator()();
    void display();
private:
    unsigned int progress;
    unsigned int goal;
    bool visible;
    std::string toolname;
};
    

#endif
