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

#include <cassert>

#include "MeshAttributes_ab.hpp"
#include "Exception.hpp"

template<unsigned SPACE_DIM>
MeshAttributes_ab<SPACE_DIM>::MeshAttributes_ab()
    :   mAttributes(std::vector<double>()),
        mRegion(0u),
        mAppliedForce(zero_vector<double>(SPACE_DIM)),
        mRadius(0.0),
        mNeighbourIndices(std::vector<unsigned>()),
        mNeighboursSetUp(false),
        mIsParticle(false)
{
}

template<unsigned SPACE_DIM>
std::vector<double>& MeshAttributes_ab<SPACE_DIM>::rGetAttributes()
{
    return mAttributes;
}

template<unsigned SPACE_DIM>
void MeshAttributes_ab<SPACE_DIM>::AddAttribute(double attribute)
{
    mAttributes.push_back(attribute);
}

template<unsigned SPACE_DIM>
unsigned MeshAttributes_ab<SPACE_DIM>::GetRegion()
{
    return mRegion;
}

template<unsigned SPACE_DIM>
void MeshAttributes_ab<SPACE_DIM>::SetRegion(unsigned region)
{
    mRegion = region;
}

template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM>& MeshAttributes_ab<SPACE_DIM>::rGetAppliedForce()
{
    return mAppliedForce;
}

//added 7-5-20
template<unsigned SPACE_DIM>
std::vector<double>& MeshAttributes_ab<SPACE_DIM>::rGetSubstrateNodeLocation()
{
    return mSubstrateNodeLocation;
}


//added 7-14-20
template<unsigned SPACE_DIM>
double& MeshAttributes_ab<SPACE_DIM>::rGetSubstrateNodeDistance()
{
    return mSubstrateNodeDistance;
}


//added 7-16-20
template<unsigned SPACE_DIM>
unsigned& MeshAttributes_ab<SPACE_DIM>::rGetSubstrateNodeIndex()
{
    return mSubstrateNodeIndex;
}


//added 7-14-20
template<unsigned SPACE_DIM>
void MeshAttributes_ab<SPACE_DIM>::HasFocal(bool flag1)
{
    mflag1 = flag1;
}

//added 7-14-20
template<unsigned SPACE_DIM>
bool MeshAttributes_ab<SPACE_DIM>::IsFocal()
{
    return mflag1;
}


//added 7-5-20
template<unsigned SPACE_DIM>
void MeshAttributes_ab<SPACE_DIM>::AddSubstrateNodeLocation(const c_vector<double, SPACE_DIM>& rSubstrateNodeLocation)
{
    //mSubstrateNodeLocation = rSubstrateNodeLocation;
	
	for (unsigned i=0; i < SPACE_DIM; i++)
        {
		mSubstrateNodeLocation.push_back(rSubstrateNodeLocation[i]);              //added 7-5-20
        } 
}


//added 7-14-20
template<unsigned SPACE_DIM>
void MeshAttributes_ab<SPACE_DIM>::AddSubstrateNodeDistance(double& rSubstrateNodeDistance)
{
    mSubstrateNodeDistance = rSubstrateNodeDistance;	
}


//added 7-16-20
template<unsigned SPACE_DIM>
void MeshAttributes_ab<SPACE_DIM>::AddSubstrateNodeIndex(unsigned& rSubstrateNodeIndex)
{
    mSubstrateNodeIndex = rSubstrateNodeIndex;	
}



template<unsigned SPACE_DIM>
void MeshAttributes_ab<SPACE_DIM>::AddAppliedForceContribution(const c_vector<double, SPACE_DIM>& rForceContribution)
{
    mAppliedForce += rForceContribution;
}

//added 7-7-20
template<unsigned SPACE_DIM>
void MeshAttributes_ab<SPACE_DIM>::ClearNodeLocation()
{
    //std::cout<<"loveandhate222255555"<<std::endl;
    mSubstrateNodeLocation.clear();
}


template<unsigned SPACE_DIM>
void MeshAttributes_ab<SPACE_DIM>::ClearAppliedForce()
{
    //std::cout<<"loveandhate222255555"<<std::endl;
    for (unsigned d = 0; d < SPACE_DIM; d++)
    {
        mAppliedForce[d] = 0.0;
    }
}

template<unsigned SPACE_DIM>
void MeshAttributes_ab<SPACE_DIM>::AddNeighbour(unsigned index)
{
    mNeighbourIndices.push_back(index);
}

template<unsigned SPACE_DIM>
void MeshAttributes_ab<SPACE_DIM>::ClearNeighbours()
{
    mNeighbourIndices.clear();
}

template<unsigned SPACE_DIM>
void MeshAttributes_ab<SPACE_DIM>::RemoveDuplicateNeighbours()
{
    sort( mNeighbourIndices.begin(), mNeighbourIndices.end() );
    mNeighbourIndices.erase( unique( mNeighbourIndices.begin(), mNeighbourIndices.end() ), mNeighbourIndices.end() );
}

template<unsigned SPACE_DIM>
bool MeshAttributes_ab<SPACE_DIM>::NeighboursIsEmpty()
{
    return mNeighbourIndices.empty();
}

template<unsigned SPACE_DIM>
void MeshAttributes_ab<SPACE_DIM>::SetNeighboursSetUp(bool flag)
{
    mNeighboursSetUp = flag;
};

template<unsigned SPACE_DIM>
bool MeshAttributes_ab<SPACE_DIM>::GetNeighboursSetUp()
{
    return mNeighboursSetUp;
};

template<unsigned SPACE_DIM>
std::vector<unsigned>& MeshAttributes_ab<SPACE_DIM>::rGetNeighbours()
{
    return mNeighbourIndices;
};


template<unsigned SPACE_DIM>
bool MeshAttributes_ab<SPACE_DIM>::IsParticle()
{
    return mIsParticle;
}

template<unsigned SPACE_DIM>
void MeshAttributes_ab<SPACE_DIM>::SetIsParticle(bool isParticle)
{
    mIsParticle = isParticle;
}

template<unsigned SPACE_DIM>
double MeshAttributes_ab<SPACE_DIM>::GetRadius()
{
    return mRadius;
}

template<unsigned SPACE_DIM>
void MeshAttributes_ab<SPACE_DIM>::SetRadius(double radius)
{
    if (radius < 0.0)
    {
        EXCEPTION("Trying to set node attributes mRadius to a negative value.");
    }

    mRadius = radius;
}

// Explicit instantiation
template class MeshAttributes_ab<1>;
template class MeshAttributes_ab<2>;
template class MeshAttributes_ab<3>;
