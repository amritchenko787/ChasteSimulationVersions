/*

Copyright (c) 2005-2018, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef STRAIGHTVERTEXMESHGENERATOR_AB3_HPP_
#define STRAIGHTVERTEXMESHGENERATOR_AB3_HPP_

#include <cmath>
#include <vector>

#include "MutableVertexMesh.hpp"

/**
 * Honeycomb mesh generator that creates a 2D honeycomb mesh (with equal distance
 * between nodes) for use in vertex simulations.
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 */
class StraightVertexMeshGenerator_ab3
{
protected:

    /** A pointer to the mesh this class creates */
    MutableVertexMesh<2,2>* mpMesh;

public:

    /**
     * Constructor.
     *
     * @param numElementsAcross  The number of columns of elements in the mesh
     * @param numElementsUp  The number of rows of elements in the mesh
     * @param isFlatBottom  Whether to enforce a flat bottom to the mesh (defaults to false)
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangement (defaults to 0.01)
     * @param t2Threshold the maximum threshold distance for Type 2 swaps (defaults to 0.001)
     * @param elementArea the element area, which has default value 0.5*sqrt(3.0)
     */
    StraightVertexMeshGenerator_ab3(unsigned numElementsAcross,
                                 unsigned numElementsUp,
				 double ycoordinate,     //modified 7-12-21
				 bool isSubstrate=false,     //added 9-16-21
				 bool isCollagen=false,     //added 7-12-21
				 bool isLongCollagen=false,     //added 8-12-21
				 bool isHaptotaxis=false,	//added 8-14-23
				 double numNodesRemoved=0,	//added  8-13-21
				 unsigned shortFiberConnectivity=0,	//added  12-25-21	
				 unsigned numSubNode = 3,	//added 3-21-22
				 bool isFlatBottom=false,
                                 double cellRearrangementThreshold=0.01,
                                 double t2Threshold=0.001,
                                 double elementArea=0.5*sqrt(3.0));

    /**
     * Null constructor for derived classes to call.
     */
    StraightVertexMeshGenerator_ab3()
    {
    }

    /**
     * Destructor - deletes the mesh object and pointer.
     */
    virtual ~StraightVertexMeshGenerator_ab3();

    /**
     * @return a 2D honeycomb mesh
     */
    virtual MutableVertexMesh<2,2>* GetMesh();
};

#endif /*STRAIGHTVERTEXMESHGENERATOR_AB3_HPP_*/
