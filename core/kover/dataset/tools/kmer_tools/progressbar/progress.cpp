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

#include "progress.h"

ProgressBar::ProgressBar(unsigned int goal, bool visible, std::string toolname):progress(0), goal(goal), visible(visible), toolname(toolname)
{if (visible)
{
    std::cout << std::endl;
    display();
}
}

void ProgressBar::operator ()()
{
    ++progress;
    if (visible)
    {
        display();
    }
}

void ProgressBar::display()
{
    float percentage = (1.0 * progress / goal);
    unsigned int width = 50;
    
    std::cout << toolname ;
    std::cout << " [";
    unsigned int position = width * percentage;
    for (unsigned int i = 0; i < width; ++i) 
    {
        if (i < position) std::cout << "=";
        else if (i == position) std::cout << ">";
        else std::cout << " ";
    }
    if (progress == goal)
    {
        std::cout << "] " << int(percentage * 100.0) << "%\n";
    }
    else
    {
        std::cout << "] " << int(percentage * 100.0) << " %\r";
        std::cout.flush();
    }
    
}
