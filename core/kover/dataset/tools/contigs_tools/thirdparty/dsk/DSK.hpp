/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#ifndef _DSK_HPP_
#define _DSK_HPP_

/********************************************************************************/

#include <gatb/gatb_core.hpp>

/********************************************************************************/

/** \brief Kmer counting class
 *
 * This is the high level class for running DSK, the counting kmer algorithm. It is an
 * implementation of the Tool interface.
 *
 * The real job is delegated to another class: DSKAlgorithm.
 * Actually, DSK analyzes the kmer size chosen by the user and
 * then decides what kind of DSKAlgorithm to use according to the
 * kmer size.
 */
class DSK : public gatb::core::tools::misc::impl::Tool
{
public:

    /** Constructor. */
    DSK ();

    /** Get the default Storage mode.
     * \return the storage mode (likely HDF5).  */
    static StorageMode_e getStorageMode() { return STORAGE_HDF5; }

private:

    /** \copydoc Tool::execute. */
    void  execute ();
};

/********************************************************************************/

#endif /* _DSK_HPP_ */

